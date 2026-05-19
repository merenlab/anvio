# pylint: disable=line-too-long
"""
    Unit tests for gene call table helpers.
"""

import unittest

import numpy

import anvio

from anvio.errors import ConfigError
from anvio.tables.genecalls import TablesForGeneCalls


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "The Anvi'o Project"
__email__ = "anvio@anvio.org"


class GeneCallsIntegerNormalizationTestCase(unittest.TestCase):
    def setUp(self):
        self.gene_calls_table = TablesForGeneCalls.__new__(TablesForGeneCalls)


    def test_numpy_integer_values_are_normalized_to_python_ints(self):
        gene_calls_dict = {
            numpy.int64(1): {
                'contig': 'contig_01',
                'start': numpy.int32(0),
                'stop': numpy.uint64(300),
                'direction': 'f',
                'partial': numpy.int8(0),
                'call_type': numpy.int64(1),
                'source': 'external',
                'version': 'unknown',
            }
        }

        normalized = self.gene_calls_table.normalize_gene_calls_dict_integer_types(gene_calls_dict)

        self.assertEqual(list(normalized.keys()), [1])
        self.assertIs(type(list(normalized.keys())[0]), int)

        for field_name in ['start', 'stop', 'partial', 'call_type']:
            self.assertIs(type(normalized[1][field_name]), int)

        self.gene_calls_table.check_gene_calls_dict(normalized)


    def test_boolean_values_are_not_accepted_as_integers(self):
        with self.assertRaises(ConfigError):
            self.gene_calls_table.normalize_gene_call_integer_value(True, 'partial')


    def test_integer_normalization_does_not_silently_collapse_duplicate_keys(self):
        gene_calls_dict = {
            '1': {
                'contig': 'contig_01',
                'start': 0,
                'stop': 300,
                'direction': 'f',
                'partial': 0,
                'call_type': 1,
                'source': 'external',
                'version': 'unknown',
            },
            numpy.int64(1): {
                'contig': 'contig_02',
                'start': 0,
                'stop': 300,
                'direction': 'f',
                'partial': 0,
                'call_type': 1,
                'source': 'external',
                'version': 'unknown',
            }
        }

        with self.assertRaises(ConfigError):
            self.gene_calls_table.normalize_gene_calls_dict_integer_types(gene_calls_dict)


if __name__ == '__main__':
    unittest.main()
