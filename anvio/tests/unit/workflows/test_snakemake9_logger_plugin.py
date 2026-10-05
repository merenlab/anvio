"""Tests for the Snakemake 9 logger plugin used by anvi'o workflows."""

import logging
import os
import sys
import types
import tempfile
import unittest
from unittest.mock import Mock, patch

from snakemake_interface_logger_plugins.common import LogEvent
from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry

import anvio.workflows as workflows
from anvio.errors import ConfigError
from anvio.workflows.logger import register_workflow_logger
from anvio.workflows.scripts.manifest import initialize_manifest
from snakemake_interface_common.exceptions import InvalidPluginException
from snakemake_logger_plugin_anvio import LogHandler


class WorkflowLoggerRegistrationTestCase(unittest.TestCase):
    def test_registers_importable_plugin_once(self):
        registry = LoggerPluginRegistry()
        with patch.dict(registry.plugins, {}, clear=True):
            register_workflow_logger()
            plugin = registry.get_plugin('anvio')
            self.assertIs(plugin.log_handler, LogHandler)
            with patch.object(registry, 'register_plugin') as register:
                register_workflow_logger()
            register.assert_not_called()
            self.assertIs(registry.get_plugin('anvio'), plugin)

    def test_preserves_existing_registration(self):
        registry = LoggerPluginRegistry()
        existing = object()
        with patch.dict(registry.plugins, {'anvio': existing}, clear=True):
            # Even an unavailable bundled module must not disturb an existing plugin.
            with patch.dict(sys.modules, {'snakemake_logger_plugin_anvio': None}):
                register_workflow_logger()
            self.assertIs(registry.get_plugin('anvio'), existing)

    def test_interface_import_failure_preserves_cause(self):
        with patch.dict(sys.modules, {'snakemake_interface_logger_plugins.registry': None}):
            with self.assertRaises(ConfigError) as raised:
                register_workflow_logger()
        self.assertIsInstance(raised.exception.__cause__, ImportError)
        self.assertIn('Snakemake interfaces', str(raised.exception))

    def test_non_class_handler_validation_preserves_cause(self):
        registry = LoggerPluginRegistry()
        malformed = types.ModuleType('snakemake_logger_plugin_anvio')
        malformed.LogHandler = object()
        with patch.dict(registry.plugins, {}, clear=True):
            with patch.dict(sys.modules, {'snakemake_logger_plugin_anvio': malformed}):
                with self.assertRaises(ConfigError) as raised:
                    register_workflow_logger()
        self.assertIsInstance(raised.exception.__cause__, TypeError)
        self.assertIn("workflow logger", str(raised.exception))


    def make_workflow(self, directory):
        workflow = workflows.WorkflowSuperClass.__new__(workflows.WorkflowSuperClass)
        workflow.args = types.SimpleNamespace(workflow='contigs', config_file='config.json')
        workflow.name = 'contigs'
        workflow.dirs_dict = {'LOGS_DIR': directory}
        workflow.run = Mock()
        workflow.rules = []
        workflow.config = {}
        workflow.additional_params = []
        workflow.list_dependencies = False
        workflow.dry_run_only = False
        workflow.save_workflow_graph = False
        workflow.get_max_num_cpus_requested_by_the_workflow = Mock(return_value=1)
        workflow.sanity_checks = Mock()
        workflow.dry_run = Mock()
        workflow.pre_execution_checks = Mock()
        return workflow

    def test_registration_failures_precede_process_and_manifest_changes(self):
        registry = LoggerPluginRegistry()
        invalid = InvalidPluginException('anvio', 'invalid handler')
        with tempfile.TemporaryDirectory() as directory:
            manifest = os.path.join(directory, 'contigs-workflow-manifest.tsv')
            for old_manifest in (None, 'existing-manifest.tsv'):
                for failure in ('import', 'validation'):
                    with self.subTest(old_manifest=old_manifest, failure=failure):
                        workflow = self.make_workflow(directory)
                        original_argv = sys.argv
                        with patch.dict(os.environ):
                            if old_manifest is None:
                                os.environ.pop('ANVIO_WORKFLOW_MANIFEST_PATH', None)
                            else:
                                os.environ['ANVIO_WORKFLOW_MANIFEST_PATH'] = old_manifest
                            original_env = dict(os.environ)
                            with patch.dict(registry.plugins, {}, clear=True):
                                context = (patch.dict(sys.modules, {'snakemake_logger_plugin_anvio': None})
                                           if failure == 'import' else
                                           patch.object(registry, 'register_plugin', side_effect=invalid))
                                with context, patch.object(workflows, 'snakemake_main') as main:
                                    with self.assertRaises(ConfigError) as raised:
                                        workflow.go(skip_dry_run=True)
                            main.assert_not_called()
                            self.assertIs(sys.argv, original_argv)
                            self.assertEqual(dict(os.environ), original_env)
                            self.assertFalse(os.path.exists(manifest))
                            self.assertIsInstance(raised.exception.__cause__,
                                                  ImportError if failure == 'import' else InvalidPluginException)

    def test_dry_run_and_dependency_listing_do_not_register(self):
        with tempfile.TemporaryDirectory() as directory:
            for mode in ('dry_run', 'dependencies'):
                with self.subTest(mode=mode):
                    workflow = self.make_workflow(directory)
                    workflow.dry_run_only = mode == 'dry_run'
                    workflow.list_dependencies = mode == 'dependencies'
                    original_argv = sys.argv
                    original_env = dict(os.environ)
                    with patch('anvio.workflows.logger.register_workflow_logger') as register:
                        captured_argv = []
                        with patch.object(workflows, 'snakemake_main', side_effect=lambda: captured_argv.extend(sys.argv)) as main:
                            if mode == 'dependencies':
                                with self.assertRaises(SystemExit):
                                    workflow.go()
                                self.assertNotIn('--logger', captured_argv)
                            else:
                                workflow.go()
                                main.assert_not_called()
                    register.assert_not_called()
                    workflow.dry_run.assert_called_once()
                    self.assertEqual(sys.argv, original_argv)
                    self.assertEqual(dict(os.environ), original_env)
                    self.assertFalse(os.path.exists(os.path.join(directory, 'contigs-workflow-manifest.tsv')))


class Snakemake9LoggerPluginTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.manifest_path = os.path.join(self.temp_dir.name, 'workflow-manifest.tsv')
        initialize_manifest(self.manifest_path)
        self.old_manifest_path = os.environ.get('ANVIO_WORKFLOW_MANIFEST_PATH')
        os.environ['ANVIO_WORKFLOW_MANIFEST_PATH'] = self.manifest_path
        self.handler = LogHandler(common_settings=None, settings=None)

    def tearDown(self):
        if self.old_manifest_path is None:
            os.environ.pop('ANVIO_WORKFLOW_MANIFEST_PATH', None)
        else:
            os.environ['ANVIO_WORKFLOW_MANIFEST_PATH'] = self.old_manifest_path
        self.temp_dir.cleanup()

    def test_job_info_and_completion_write_manifest_row(self):
        job_info = logging.LogRecord('snakemake', logging.INFO, '', 0, 'Job info', (), None)
        job_info.event = LogEvent.JOB_INFO
        job_info.jobid = 4
        job_info.rule_name = 'anvi_profile'
        job_info.wildcards = {'group': 'G01', 'readset': 'S01'}
        job_info.log = ['00_LOGS/metagenomics/anvi_profile/G01-S01.log']
        self.handler.emit(job_info)

        job_finished = logging.LogRecord('snakemake', logging.INFO, '', 0, 'Finished', (), None)
        job_finished.event = LogEvent.JOB_FINISHED
        job_finished.job_id = 4
        self.handler.emit(job_finished)

        log_path = logging.LogRecord('snakemake', logging.INFO, '', 0,
                                     'Complete log(s): .snakemake/log/run.snakemake.log', (), None)
        self.handler.emit(log_path)

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'succeeded\tanvi_profile\tG01\tS01\t'
                                 '00_LOGS/metagenomics/anvi_profile/G01-S01.log\t'
                                 '.snakemake/log/run.snakemake.log')


if __name__ == '__main__':
    unittest.main()
