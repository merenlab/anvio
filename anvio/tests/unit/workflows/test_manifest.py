"""
    Unit tests for workflow manifest helpers.
"""

import os
import sys
import types
import tempfile
import unittest
import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = REPO_ROOT / 'workflows' / 'manifest.py'
LOG_HANDLER_PATH = REPO_ROOT / 'workflows' / 'snakemake_log_handler.py'


def load_module_from_path(module_name, path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def load_log_handler_module(manifest_module):
    anvio_module = types.ModuleType('anvio')
    workflows_module = types.ModuleType('anvio.workflows')

    old_modules = {name: sys.modules.get(name) for name in ('anvio',
                                                            'anvio.workflows',
                                                            'anvio.workflows.manifest')}
    sys.modules['anvio'] = anvio_module
    sys.modules['anvio.workflows'] = workflows_module
    sys.modules['anvio.workflows.manifest'] = manifest_module

    try:
        return load_module_from_path('test_snakemake_log_handler', LOG_HANDLER_PATH)
    finally:
        for name, module in old_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


class WorkflowManifestTestCase(unittest.TestCase):
    def setUp(self):
        self.manifest = load_module_from_path('test_workflow_manifest', MANIFEST_PATH)
        self.temp_dir = tempfile.TemporaryDirectory()
        self.manifest_path = os.path.join(self.temp_dir.name, 'metagenomics-workflow-manifest.tsv')


    def tearDown(self):
        self.temp_dir.cleanup()


    def test_manifest_row_contains_requested_columns(self):
        self.manifest.initialize_manifest(self.manifest_path)
        self.manifest.append_manifest_row(self.manifest_path,
                                          'succeeded',
                                          'anvi_profile',
                                          group='G01',
                                          read='S01',
                                          log_path='00_LOGS/G01-S01-anvi_profile.log',
                                          snakemake_log_path='.snakemake/log/2026-01-21T005403.112748.snakemake.log')

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[0], 'status\trule\tgroup\tread\tlog_path\tsnakemake_log_path')
        self.assertEqual(rows[1], 'succeeded\tanvi_profile\tG01\tS01\t00_LOGS/G01-S01-anvi_profile.log\t.snakemake/log/2026-01-21T005403.112748.snakemake.log')


    def test_update_snakemake_log_path_backfills_existing_rows(self):
        self.manifest.initialize_manifest(self.manifest_path)
        self.manifest.append_manifest_row(self.manifest_path,
                                          'failed',
                                          'bowtie',
                                          group='G01',
                                          read='S02',
                                          log_path='00_LOGS/G01-S02-bowtie.log')

        self.manifest.update_snakemake_log_path(self.manifest_path,
                                                '.snakemake/log/2026-01-21T005403.112748.snakemake.log')

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'failed\tbowtie\tG01\tS02\t00_LOGS/G01-S02-bowtie.log\t.snakemake/log/2026-01-21T005403.112748.snakemake.log')


class SnakemakeLogHandlerTestCase(unittest.TestCase):
    def setUp(self):
        self.manifest = load_module_from_path('test_workflow_manifest_for_handler', MANIFEST_PATH)
        self.log_handler = load_log_handler_module(self.manifest)
        self.temp_dir = tempfile.TemporaryDirectory()
        self.manifest_path = os.path.join(self.temp_dir.name, 'metagenomics-workflow-manifest.tsv')
        self.manifest.initialize_manifest(self.manifest_path)

        self.old_manifest_path = os.environ.get('ANVIO_WORKFLOW_MANIFEST_PATH')
        os.environ['ANVIO_WORKFLOW_MANIFEST_PATH'] = self.manifest_path


    def tearDown(self):
        if self.old_manifest_path is None:
            os.environ.pop('ANVIO_WORKFLOW_MANIFEST_PATH', None)
        else:
            os.environ['ANVIO_WORKFLOW_MANIFEST_PATH'] = self.old_manifest_path

        self.temp_dir.cleanup()


    def test_handler_records_rule_status_and_wildcards(self):
        self.log_handler.log_handler({'level': 'job_info',
                                      'jobid': 7,
                                      'name': 'anvi_profile',
                                      'wildcards': {'group': 'G01', 'readset': 'S01'},
                                      'log': ['00_LOGS/G01-S01-anvi_profile.log']})
        self.log_handler.log_handler({'level': 'job_finished',
                                      'jobid': 7})

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'succeeded\tanvi_profile\tG01\tS01\t00_LOGS/G01-S01-anvi_profile.log\t')


    def test_handler_backfills_complete_snakemake_log_path(self):
        self.log_handler.log_handler({'level': 'job_info',
                                      'jobid': 9,
                                      'name': 'anvi_merge',
                                      'wildcards': {'group': 'G01'},
                                      'log': ['00_LOGS/G01-anvi_merge.log']})
        self.log_handler.log_handler({'level': 'job_error',
                                      'jobid': 9})
        self.log_handler.log_handler({'level': 'error',
                                      'msg': 'Exiting because a job execution failed. Look above for error message\nComplete log: .snakemake/log/2026-01-21T005403.112748.snakemake.log'})

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'failed\tanvi_merge\tG01\t\t00_LOGS/G01-anvi_merge.log\t.snakemake/log/2026-01-21T005403.112748.snakemake.log')


if __name__ == '__main__':
    unittest.main()
