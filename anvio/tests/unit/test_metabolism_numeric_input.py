"""Regression tests for numeric missing values in metabolism inputs."""

import os
import tempfile
import unittest

import pandas as pd

import anvio.terminal as terminal
from anvio.errors import ConfigError
from anvio.metabolism.dbaccess import KeggEstimatorArgs
from anvio.metabolism.enrichment import KeggModuleEnrichment
from anvio.metabolism.input_validation import PANDAS_DEFAULT_MISSING_MARKERS
from anvio.metabolism.input_validation import parse_numeric_column


class TestNumericColumnDiagnostics(unittest.TestCase):
    def test_duplicate_index_reports_each_invalid_row_position(self):
        data = pd.DataFrame(
            {'sample': ['first-sample', 'second-sample'], 'score': ['bad', 'inf']},
            index=[7, 7],
        )

        with self.assertRaises(ConfigError) as error:
            parse_numeric_column(data, 'score', 'modules.tsv', identifier_column='sample')

        message = str(error.exception)
        for expected in ['score', 'row 2', 'row 3', 'first-sample', 'second-sample', "'bad'", 'inf']:
            self.assertIn(expected, message)


class TestEnzymesNumericInputs(unittest.TestCase):
    def load_enzymes_file(self, contents):
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        enzymes_path = os.path.join(temp_dir.name, 'enzymes.tsv')
        with open(enzymes_path, 'w') as enzymes_file:
            enzymes_file.write(contents)

        loader = KeggEstimatorArgs.__new__(KeggEstimatorArgs)
        loader.enzymes_txt = enzymes_path
        loader.add_coverage = True
        loader.progress = terminal.Progress(verbose=False)
        loader.run = terminal.Run(verbose=False, log_file_path=None)
        return loader.get_enzymes_of_interest_df()

    def test_preserves_missing_numeric_markers_and_literal_na_gene_id(self):
        enzymes = self.load_enzymes_file(
            'gene_id\tenzyme_accession\tsource\tcoverage\tdetection\n'
            'NA\tK00001\tKOfam\t\tNA\n'
            'gene-2\tK00002\tKOfam\tN/A\tNaN\n'
        )

        self.assertEqual(enzymes['gene_id'].tolist(), ['NA', 'gene-2'])
        self.assertTrue(enzymes['coverage'].isna().all())
        self.assertTrue(enzymes['detection'].isna().all())

    def test_accepts_pandas_default_missing_markers_in_numeric_columns(self):
        rows = [f'gene-{index}\tK00001\tKOfam\t{marker}\t{marker}'
                for index, marker in enumerate(sorted(PANDAS_DEFAULT_MISSING_MARKERS))]
        enzymes = self.load_enzymes_file(
            'gene_id\tenzyme_accession\tsource\tcoverage\tdetection\n' + '\n'.join(rows) + '\n'
        )

        self.assertTrue(enzymes['coverage'].isna().all())
        self.assertTrue(enzymes['detection'].isna().all())

    def test_rejects_malformed_coverage_with_field_and_gene_context(self):
        with self.assertRaises(ConfigError) as error:
            self.load_enzymes_file(
                'gene_id\tenzyme_accession\tsource\tcoverage\tdetection\n'
                'bad-coverage\tK00001\tKOfam\tbogus\t0.5\n'
            )
        message = str(error.exception)
        self.assertIn('coverage', message)
        self.assertIn('bad-coverage', message)
        self.assertIn('bogus', message)

    def test_rejects_infinite_detection(self):
        with self.assertRaises(ConfigError) as error:
            self.load_enzymes_file(
                'gene_id\tenzyme_accession\tsource\tcoverage\tdetection\n'
                'bad-detection\tK00001\tKOfam\t2.5\tinf\n'
            )
        message = str(error.exception)
        self.assertIn('detection', message)
        self.assertIn('bad-detection', message)
        self.assertIn('inf', message)


class TestEnrichmentNumericInputs(unittest.TestCase):
    def run_enrichment_input(self, modules_text):
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        modules_path = os.path.join(temp_dir.name, 'modules.tsv')
        groups_path = os.path.join(temp_dir.name, 'groups.tsv')
        output_path = os.path.join(temp_dir.name, 'enrichment-input.tsv')

        with open(modules_path, 'w') as modules_file:
            modules_file.write(modules_text)
        with open(groups_path, 'w') as groups_file:
            groups_file.write('sample\tgroup\nNA\tA\nblank\tA\nmissing\tB\n')

        enrichment = KeggModuleEnrichment.__new__(KeggModuleEnrichment)
        enrichment.modules_txt = modules_path
        enrichment.groups_txt = groups_path
        enrichment.sample_header_in_modules_txt = 'sample'
        enrichment.module_completion_threshold = 0.75
        enrichment.use_stepwise_completeness = False
        enrichment.include_missing = False
        enrichment.just_do_it = False
        enrichment.progress = terminal.Progress(verbose=False)
        enrichment.run = terminal.Run(verbose=False, log_file_path=None)
        enrichment.get_enrichment_input(output_path)
        return pd.read_csv(output_path, sep='\t', keep_default_na=False)

    def test_preserves_missing_completeness_and_literal_na_sample(self):
        output = self.run_enrichment_input(
            'sample\tmodule\tmodule_name\tpathwise_module_completeness\n'
            'NA\tM00001\tTest module\t0.8\n'
            'blank\tM00001\tTest module\t\n'
            'missing\tM00001\tTest module\tN/A\n'
        )

        self.assertEqual(output['sample_ids'].tolist(), ['NA'])
        self.assertEqual(output['accession'].tolist(), ['M00001'])

    def test_rejects_malformed_completeness(self):
        with self.assertRaises(ConfigError) as error:
            self.run_enrichment_input(
                'sample\tmodule\tmodule_name\tpathwise_module_completeness\n'
                'NA\tM00001\tTest module\t0.8\n'
                'blank\tM00001\tTest module\t\n'
                'bad-sample\tM00001\tTest module\toops\n'
            )
        message = str(error.exception)
        self.assertIn('pathwise_module_completeness', message)
        self.assertIn('bad-sample', message)
        self.assertIn('oops', message)

    def test_rejects_infinite_completeness(self):
        with self.assertRaises(ConfigError) as error:
            self.run_enrichment_input(
                'sample\tmodule\tmodule_name\tpathwise_module_completeness\n'
                'NA\tM00001\tTest module\t0.8\n'
                'blank\tM00001\tTest module\t\n'
                'bad-sample\tM00001\tTest module\t-inf\n'
            )
        message = str(error.exception)
        self.assertIn('pathwise_module_completeness', message)
        self.assertIn('bad-sample', message)
        self.assertIn('-inf', message)
