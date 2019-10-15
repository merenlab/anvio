# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Unit tests for the HomogeneityCalculator class.
"""

import unittest

import anvio
import anvio.homogeneityindex as homogeneityindex

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ryan Moore"
__email__ = "moorer@udel.edu"


class HomogeneityCalculatorTestCase(unittest.TestCase):
    def setUp(self):
        self.match_score = 3
        self.fully_functionally_conserved_score = 2 / 3
        self.calculator = homogeneityindex.HomogeneityCalculator(quick_homogeneity=False)

    def test_single_sequence_cluster(self):
        actual = self.calculator.compute_functional_index(['AAAAA'])
        expected = 1.0

        self.assertEqual(actual, expected)

    def test_empty_gene_clusters(self):
        actual = self.calculator.compute_functional_index([])
        expected = 0.0

        self.assertEqual(actual, expected)

    def test_identical_sequences_get_a_perfect_score(self):
        seq1 = 'AAAAA'
        seq2 = 'AAAAA'
        seq3 = 'AAAAA'

        actual = self.calculator.compute_functional_index([seq1, seq2, seq3])
        expected = 1.0

        self.assertEqual(actual, expected)

    def test_functionally_conserved_sequences_increase_the_score(self):
        """Functionally conserved residue pairs have a score of 2/3 each."""
        seq1 = 'LVI'
        seq2 = 'CHA'
        seq3 = 'VIL'

        actual = self.calculator.compute_functional_index([seq1, seq2, seq3])
        expected = self.fully_functionally_conserved_score

        self.assertEqual(actual, expected)

    def test_non_conserved_sequnces_get_the_lowest_score(self):
        """Two sequences with no functionally conserved residues get the lowest score."""
        seq1 = 'AAAAA'  # Nonpolar
        seq2 = 'FFFFF'  # Aromatic

        actual = self.calculator.compute_functional_index([seq1, seq2])
        expected = 0.0

        self.assertEqual(actual, expected)

    def test_two_gaps_dont_reduce_the_score(self):
        """Sequences should not be penalized for having the same gap position."""
        seq1 = 'AA-AA'
        seq2 = 'AA-AA'

        actual = self.calculator.compute_functional_index([seq1, seq2])
        expected = 1.0

        self.assertEqual(actual, expected)

    def test_a_gap_with_a_residue_is_a_mismatch(self):
        """If one sequence has a residue where another has a gap, that is a mismatch."""
        seq1 = 'AAAAA'
        seq2 = 'AA-AA'

        actual = self.calculator.compute_functional_index([seq1, seq2])
        expected = 0.8

        self.assertEqual(actual, expected)

    def test_two_gaps_count_as_a_match_and_increase_the_score(self):
        """If a pair of residues are both gaps, they count as an exact match and increase the similarity score."""
        seq1 = 'AAAAA'
        seq1_gapped = 'AA-AA'

        two_seq_score = self.calculator.compute_functional_index([seq1, seq1_gapped])
        three_seq_score = self.calculator.compute_functional_index([seq1, seq1_gapped, seq1_gapped])

        self.assertLess(two_seq_score, three_seq_score)

    def test_weird_residues_never_get_full_score(self):
        """Even when they match exactly, these residues will get the functionally conserved score rather than the
        exact match score. """
        for residue in ['J', 'B', 'Z']:
            seq = residue * 5

            actual = self.calculator.compute_functional_index([seq, seq])
            expected = self.fully_functionally_conserved_score

            self.assertEqual(actual, expected)

    def test_two_Xs_count_as_a_mismatch(self):
        """Two Xs count as a mismatch even when they are in the same alignment column in two sequences."""
        seq = 'X' * 5

        actual = self.calculator.compute_functional_index([seq, seq])
        expected = 0.0

        self.assertEqual(actual, expected)

    # This is currently broken.
    @unittest.expectedFailure
    def test_empty_sequences(self):
        seq1 = ''
        seq2 = ''

        actual = self.calculator.compute_functional_index([seq1, seq2])
        expected = 1.0

        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
