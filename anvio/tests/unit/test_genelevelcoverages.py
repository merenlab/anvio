"""Tests for gene-level coverage serialization and recovery."""

import gzip
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

import anvio
import anvio.db as db
import anvio.tables as t
from anvio.tables.genelevelcoverages import TableForGeneLevelCoverages
from anvio.errors import ConfigError


class FakeDatabase:
    def __init__(self, table_data=None):
        self.table_data = table_data or {}
        self._exec_many = mock.Mock()
        self.disconnect = mock.Mock()

    def get_meta_value(self, key):
        return True

    def get_table_as_dict(self, table_name):
        return self.table_data

    def disconnect(self):
        pass

    def remove_meta_key_value_pair(self, key):
        pass

    def set_meta_value(self, key, value):
        pass

    def update_meta_value(self, key, value):
        pass


class GeneLevelCoverageSerializationTestCase(unittest.TestCase):
    def setUp(self):
        self.table = object.__new__(TableForGeneLevelCoverages)
        self.table.db_path = 'unused.db'
        self.table.table_name = 'gene_level_stats'
        self.table.parameters = {}
        self.table.mode = 'STANDARD'
        self.table.progress = mock.Mock()
        self.table.run = mock.Mock()
        self.table.print_info = mock.Mock()
        self.table.check_split_names = mock.Mock()
        self.table.check_params = mock.Mock()

    def test_read_recovers_non_outlier_positions_as_bool(self):
        coverage = np.array([5, 65535, 2], dtype=np.uint16)
        non_outliers = np.array([True, False, True], dtype=bool)
        entry = {
            'gene_callers_id': 1,
            'sample_name': 'S1',
            'gene_coverage_values_per_nt': gzip.compress(coverage.tobytes()),
            'non_outlier_positions': gzip.compress(non_outliers.tobytes()),
        }
        database = FakeDatabase({'entry': entry})

        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=database), \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            result = self.table.read()

        recovered = result[1]['S1']
        np.testing.assert_array_equal(recovered['gene_coverage_values_per_nt'], coverage)
        np.testing.assert_array_equal(recovered['non_outlier_positions'], non_outliers)
        self.assertEqual(recovered['non_outlier_positions'].dtype, np.dtype(bool))

    def test_read_accepts_historical_and_current_mask_encodings(self):
        for dtype in (np.uint16, bool):
            for values in ([False, False, False], [True, True, True], [True, False, True]):
                with self.subTest(dtype=dtype, values=values):
                    entry = {
                        'gene_callers_id': 1, 'sample_name': 'S1',
                        'gene_coverage_values_per_nt': gzip.compress(np.array([5, 8, 2], dtype=np.uint16).tobytes()),
                        'non_outlier_positions': gzip.compress(np.array(values, dtype=dtype).tobytes()),
                    }
                    with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=FakeDatabase({'entry': entry})), \
                            mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
                        mask = self.table.read()[1]['S1']['non_outlier_positions']
                    np.testing.assert_array_equal(mask, values)
                    self.assertEqual(mask.dtype, np.dtype(bool))
                    self.assertEqual(len(mask), 3)

    def test_read_rejects_malformed_mask_length(self):
        entry = {
            'gene_callers_id': 1, 'sample_name': 'S1',
            'gene_coverage_values_per_nt': gzip.compress(np.array([5, 8, 2], dtype=np.uint16).tobytes()),
            'non_outlier_positions': gzip.compress(bytes([1, 0, 1, 0])),
        }
        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=FakeDatabase({'entry': entry})), \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            with self.assertRaisesRegex(ConfigError, r'coverage has 3\s+positions'):
                self.table.read()

    def test_read_rejects_invalid_historical_mask_values(self):
        entry = {
            'gene_callers_id': 1, 'sample_name': 'S1',
            'gene_coverage_values_per_nt': gzip.compress(np.array([5, 8, 2], dtype=np.uint16).tobytes()),
            'non_outlier_positions': gzip.compress(np.array([1, 2, 0], dtype=np.uint16).tobytes()),
        }
        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=FakeDatabase({'entry': entry})), \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            with self.assertRaisesRegex(ConfigError, 'zero and one'):
                self.table.read()

    def test_standard_read_rejects_missing_non_outlier_positions(self):
        entry = {
            'gene_callers_id': 1, 'sample_name': 'S1',
            'gene_coverage_values_per_nt': gzip.compress(np.array([5, 8, 2], dtype=np.uint16).tobytes()),
        }
        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=FakeDatabase({'entry': entry})), \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            with self.assertRaises(KeyError):
                self.table.read()

    def test_inseq_sqlite_roundtrip_preserves_coverage_and_insertion_metrics(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = os.path.join(temp_dir, 'genes.db')
            database = db.DB(db_path, anvio.__genes__version__, new_database=True)
            database.set_meta_value('db_type', 'genes')
            database.set_meta_value('collection_name', 'collection')
            database.set_meta_value('bin_name', 'bin')
            database.set_meta_value('gene_level_coverages_stored', False)
            database.create_table(
                t.gene_level_inseq_stats_table_name,
                t.gene_level_inseq_stats_table_structure,
                t.gene_level_inseq_stats_table_types,
            )
            database.disconnect()

            coverage = np.array([0, 2, 65535, 4], dtype=np.uint16)
            inseq_entry = {
                'gene_callers_id': 27,
                'sample_name': 'sample-A',
                'mean_coverage': 102.75,
                'insertions': 3,
                'insertions_normalized': 1.25,
                'mean_disruption': 0.75,
                'below_disruption': 2,
                'gene_coverage_values_per_nt': coverage,
            }
            writer = TableForGeneLevelCoverages(db_path, {}, 'INSEQ', ignore_splits_name_check=True)
            writer.run = mock.Mock()
            writer.progress = mock.Mock()
            writer.print_info = mock.Mock()
            writer.store({27: {'sample-A': inseq_entry}})

            reader = TableForGeneLevelCoverages(db_path, {}, 'INSEQ', ignore_splits_name_check=True)
            reader.run = mock.Mock()
            reader.progress = mock.Mock()
            reader.print_info = mock.Mock()
            recovered = reader.read()[27]['sample-A']

        np.testing.assert_array_equal(recovered['gene_coverage_values_per_nt'], coverage)
        self.assertIsNone(recovered['non_outlier_positions'])
        for field in ('mean_coverage', 'insertions', 'insertions_normalized', 'mean_disruption', 'below_disruption'):
            self.assertEqual(recovered[field], inseq_entry[field])

    def test_read_cleans_up_database_and_progress_after_invalid_mask(self):
        entry = {
            'gene_callers_id': 1, 'sample_name': 'S1',
            'gene_coverage_values_per_nt': gzip.compress(np.array([5, 8, 2], dtype=np.uint16).tobytes()),
            'non_outlier_positions': gzip.compress(bytes([1, 0, 1, 0])),
        }
        database = FakeDatabase({'entry': entry})
        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=database) as constructor, \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            with self.assertRaises(ConfigError):
                self.table.read()
        constructor.assert_called_once()
        database.disconnect.assert_called_once()
        self.table.progress.end.assert_called_once()

    def test_read_cleans_up_database_when_preflight_fails(self):
        database = FakeDatabase()
        self.table.check_params.side_effect = ConfigError('Invalid cached parameters')
        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=database), \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            with self.assertRaises(ConfigError):
                self.table.read()
        database.disconnect.assert_called_once()
        self.table.progress.new.assert_not_called()
        self.table.progress.end.assert_not_called()

    def test_store_writes_non_outlier_positions_as_bool(self):
        coverage = np.array([5, 65535, 2], dtype=np.uint16)
        non_outliers = np.array([True, False, True], dtype=bool)
        self.table.table_structure = ['gene_coverage_values_per_nt', 'non_outlier_positions']
        database = FakeDatabase()
        data = {1: {'S1': {
            'gene_coverage_values_per_nt': coverage,
            'non_outlier_positions': non_outliers,
        }}}

        with mock.patch('anvio.tables.genelevelcoverages.db.DB', return_value=database), \
                mock.patch('anvio.tables.genelevelcoverages.utils.get_required_version_for_db', return_value=1):
            self.table.store(data)

        stored = database._exec_many.call_args.args[1][0]
        recovered_coverage = np.frombuffer(gzip.decompress(stored[0]), dtype=np.uint16)
        recovered_non_outliers = np.frombuffer(gzip.decompress(stored[1]), dtype=bool)
        np.testing.assert_array_equal(recovered_coverage, coverage)
        np.testing.assert_array_equal(recovered_non_outliers, non_outliers)


if __name__ == '__main__':
    unittest.main()
