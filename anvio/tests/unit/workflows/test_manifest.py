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
MANIFEST_PATH = REPO_ROOT / 'workflows' / 'scripts' / 'manifest.py'
LOG_HANDLER_PATH = REPO_ROOT / 'workflows' / 'scripts' / 'snakemake_log_handler.py'

WORKFLOW_MANIFEST_CASES = [
    {'workflow': 'contigs',
     'status': 'succeeded',
     'rule': 'anvi_gen_contigs_database',
     'wildcards': {'group': 'G01'},
     'log_path': '00_LOGS/anvi_gen_contigs_database/G01.log',
     'expected_row': 'succeeded\tanvi_gen_contigs_database\tG01\t\t00_LOGS/anvi_gen_contigs_database/G01.log\t'},
    {'workflow': 'metagenomics',
     'status': 'failed',
     'rule': 'bowtie',
     'wildcards': {'group': 'G01', 'readset': 'S01'},
     'log_path': '00_LOGS/bowtie/G01-S01.log',
     'expected_row': 'failed\tbowtie\tG01\tS01\t00_LOGS/bowtie/G01-S01.log\t'},
    {'workflow': 'pangenomics',
     'status': 'succeeded',
     'rule': 'anvi_pan_genome',
     'wildcards': {},
     'log_path': '00_LOGS/anvi_pan_genome/PROJECT.log',
     'expected_row': 'succeeded\tanvi_pan_genome\t\t\t00_LOGS/anvi_pan_genome/PROJECT.log\t'},
    {'workflow': 'phylogenomics',
     'status': 'succeeded',
     'rule': 'iqtree',
     'wildcards': {},
     'log_path': '00_LOGS/iqtree/PROJECT.log',
     'expected_row': 'succeeded\tiqtree\t\t\t00_LOGS/iqtree/PROJECT.log\t'},
    {'workflow': 'trnaseq',
     'status': 'succeeded',
     'rule': 'anvi_trnaseq',
     'wildcards': {'sample_name': 'S01'},
     'log_path': '00_LOGS/anvi_trnaseq/S01.log',
     'expected_row': 'succeeded\tanvi_trnaseq\t\t\t00_LOGS/anvi_trnaseq/S01.log\t'},
    {'workflow': 'ecophylo',
     'status': 'succeeded',
     'rule': 'combine_sequence_data',
     'wildcards': {'group': 'Ribosomal_S3'},
     'log_path': '00_LOGS/combine_sequence_data/Ribosomal_S3.log',
     'expected_row': 'succeeded\tcombine_sequence_data\tRibosomal_S3\t\t00_LOGS/combine_sequence_data/Ribosomal_S3.log\t'},
    {'workflow': 'sra_download',
     'status': 'failed',
     'rule': 'prefetch',
     'wildcards': {'accession': 'SRR000001'},
     'log_path': '00_LOGS/prefetch/SRR000001.log',
     'expected_row': 'failed\tprefetch\t\t\t00_LOGS/prefetch/SRR000001.log\t'},
]


def load_module_from_path(module_name, path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def load_log_handler_module(manifest_module):
    anvio_module = types.ModuleType('anvio')
    workflows_module = types.ModuleType('anvio.workflows')

    scripts_module = types.ModuleType('anvio.workflows.scripts')

    old_modules = {name: sys.modules.get(name) for name in ('anvio',
                                                            'anvio.workflows',
                                                            'anvio.workflows.scripts',
                                                            'anvio.workflows.scripts.manifest')}
    sys.modules['anvio'] = anvio_module
    sys.modules['anvio.workflows'] = workflows_module
    sys.modules['anvio.workflows.scripts'] = scripts_module
    sys.modules['anvio.workflows.scripts.manifest'] = manifest_module

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
                                          log_path='00_LOGS/anvi_profile/G01-S01.log',
                                          snakemake_log_path='.snakemake/log/2026-01-21T005403.112748.snakemake.log')

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[0], 'status\trule\tgroup\tread\tlog_path\tsnakemake_log_path')
        self.assertEqual(rows[1], 'succeeded\tanvi_profile\tG01\tS01\t00_LOGS/anvi_profile/G01-S01.log\t.snakemake/log/2026-01-21T005403.112748.snakemake.log')


    def test_update_snakemake_log_path_backfills_existing_rows(self):
        self.manifest.initialize_manifest(self.manifest_path)
        self.manifest.append_manifest_row(self.manifest_path,
                                          'failed',
                                          'bowtie',
                                          group='G01',
                                          read='S02',
                                          log_path='00_LOGS/bowtie/G01-S02.log')

        self.manifest.update_snakemake_log_path(self.manifest_path,
                                                '.snakemake/log/2026-01-21T005403.112748.snakemake.log')

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'failed\tbowtie\tG01\tS02\t00_LOGS/bowtie/G01-S02.log\t.snakemake/log/2026-01-21T005403.112748.snakemake.log')


class WorkflowManifestAllWorkflowsTestCase(unittest.TestCase):
    def setUp(self):
        self.manifest = load_module_from_path('test_workflow_manifest_all_workflows', MANIFEST_PATH)
        self.temp_dir = tempfile.TemporaryDirectory()


    def tearDown(self):
        self.temp_dir.cleanup()


    def test_manifest_file_can_be_initialized_for_each_runnable_workflow(self):
        for workflow_case in WORKFLOW_MANIFEST_CASES:
            workflow = workflow_case['workflow']
            manifest_path = os.path.join(self.temp_dir.name, '00_LOGS', f'{workflow}-workflow-manifest.tsv')

            with self.subTest(workflow=workflow):
                self.manifest.initialize_manifest(manifest_path)

                with open(manifest_path) as manifest_file:
                    rows = manifest_file.read().splitlines()

                self.assertEqual(rows, ['status\trule\tgroup\tread\tlog_path\tsnakemake_log_path'])


    def test_manifest_rows_can_be_written_for_each_runnable_workflow(self):
        for workflow_case in WORKFLOW_MANIFEST_CASES:
            workflow = workflow_case['workflow']
            manifest_path = os.path.join(self.temp_dir.name, f'{workflow}-workflow-manifest.tsv')

            with self.subTest(workflow=workflow):
                self.manifest.initialize_manifest(manifest_path)
                self.manifest.append_manifest_row(manifest_path,
                                                  workflow_case['status'],
                                                  workflow_case['rule'],
                                                  group=workflow_case['wildcards'].get('group', ''),
                                                  read=workflow_case['wildcards'].get('readset', ''),
                                                  log_path=workflow_case['log_path'])

                with open(manifest_path) as manifest_file:
                    rows = manifest_file.read().splitlines()

                self.assertEqual(rows[1], workflow_case['expected_row'])


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
                                      'log': ['00_LOGS/anvi_profile/G01-S01.log']})
        self.log_handler.log_handler({'level': 'job_finished',
                                      'jobid': 7})

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'succeeded\tanvi_profile\tG01\tS01\t00_LOGS/anvi_profile/G01-S01.log\t')


    def test_handler_backfills_complete_snakemake_log_path(self):
        self.log_handler.log_handler({'level': 'job_info',
                                      'jobid': 9,
                                      'name': 'anvi_merge',
                                      'wildcards': {'group': 'G01'},
                                      'log': ['00_LOGS/anvi_merge/G01.log']})
        self.log_handler.log_handler({'level': 'job_error',
                                      'jobid': 9})
        self.log_handler.log_handler({'level': 'error',
                                      'msg': 'Exiting because a job execution failed. Look above for error message\nComplete log: .snakemake/log/2026-01-21T005403.112748.snakemake.log'})

        with open(self.manifest_path) as manifest_file:
            rows = manifest_file.read().splitlines()

        self.assertEqual(rows[1], 'failed\tanvi_merge\tG01\t\t00_LOGS/anvi_merge/G01.log\t.snakemake/log/2026-01-21T005403.112748.snakemake.log')


    def test_handler_records_representative_jobs_for_each_runnable_workflow(self):
        for jobid, workflow_case in enumerate(WORKFLOW_MANIFEST_CASES, start=100):
            workflow = workflow_case['workflow']
            manifest_path = os.path.join(self.temp_dir.name, f'{workflow}-workflow-manifest.tsv')

            with self.subTest(workflow=workflow):
                self.manifest.initialize_manifest(manifest_path)
                os.environ['ANVIO_WORKFLOW_MANIFEST_PATH'] = manifest_path
                self.log_handler.JOBS.clear()
                self.log_handler.SNAKEMAKE_LOG_PATH = ''

                self.log_handler.log_handler({'level': 'job_info',
                                              'jobid': jobid,
                                              'name': workflow_case['rule'],
                                              'wildcards': workflow_case['wildcards'],
                                              'log': [workflow_case['log_path']]})

                if workflow_case['status'] == 'succeeded':
                    self.log_handler.log_handler({'level': 'job_finished',
                                                  'jobid': jobid})
                else:
                    self.log_handler.log_handler({'level': 'job_error',
                                                  'jobid': jobid})

                with open(manifest_path) as manifest_file:
                    rows = manifest_file.read().splitlines()

                self.assertEqual(rows[1], workflow_case['expected_row'])


if __name__ == '__main__':
    unittest.main()
