"""Tests for the Snakemake 9 logger plugin used by anvi'o workflows."""

import logging
import os
import tempfile
import unittest

from snakemake_interface_logger_plugins.common import LogEvent
from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry

from anvio.workflows.scripts.manifest import initialize_manifest
from snakemake_logger_plugin_anvio import LogHandler


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

    def test_plugin_is_discoverable_by_snakemake(self):
        self.assertTrue(LoggerPluginRegistry().is_installed('anvio'))

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
