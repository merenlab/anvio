# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Unit tests for the HomogeneityCalculator class.
"""

import unittest

import anvio
import anvio.homogeneityindex as homogeneityindex

__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ryan Moore"
__email__ = "moorer@udel.edu"


class ComputeFunctionalIndexTestCase(unittest.TestCase):
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

    def test_a_gap_with_a_residue_does_not_reduce_the_score(self):
        """Sequences should not be penalized if one sequence has a residue where another has a gap."""
        seq1 = 'AAAAA'
        seq2 = 'AA-AA'

        actual = self.calculator.compute_functional_index([seq1, seq2])
        expected = 1.0

        self.assertEqual(actual, expected)

    def test_adding_in_more_seqs_with_gaps_does_not_increase_the_score(self):
        """If you add in more sequences with a gap in the same position, it does not increase the score.  See
        https://github.com/merenlab/anvio/issues/1265 for more info. """
        seq1 = 'AAAAA'
        seq1_gapped = 'AA-AA'

        two_seq_score = self.calculator.compute_functional_index([seq1, seq1_gapped])
        three_seq_score = self.calculator.compute_functional_index([seq1, seq1_gapped, seq1_gapped])

        self.assertEqual(two_seq_score, three_seq_score)

    def test_that_it_handles_gaps_properly(self):
        """For more info, see https://github.com/merenlab/anvio/issues/1265."""
        expected = 1.0

        self.assertEqual(self.calculator.compute_functional_index(['A-A', 'A-A', 'A-A']), expected)
        self.assertEqual(self.calculator.compute_functional_index(['AAA', 'AAA', 'A-A']), expected)
        self.assertEqual(self.calculator.compute_functional_index(['-AA', 'A-A', 'AA-']), expected)

    def test_weird_residues_never_get_full_score(self):
        """Even when they match exactly, these residues will get the functionally conserved score rather than the
        exact match score. """
        for residue in ['J', 'B', 'Z']:
            seq = residue * 5

            actual = self.calculator.compute_functional_index([seq, seq])
            expected = self.fully_functionally_conserved_score

            self.assertEqual(actual, expected)

    def test_ambiguous_bases_are_not_dependent_on_comparison_order(self):
        """Ambiguous residues are treated as functionally conserved with the residues they represent.  The order of
        comparison should not matter. """
        # B group
        self.assertEqual(self.calculator.compute_functional_index(['B', 'N']),
                         self.calculator.compute_functional_index(['N', 'B']))
        self.assertEqual(self.calculator.compute_functional_index(['B', 'D']),
                         self.calculator.compute_functional_index(['D', 'B']))

        # Z group
        self.assertEqual(self.calculator.compute_functional_index(['Z', 'Q']),
                         self.calculator.compute_functional_index(['Q', 'Z']))
        self.assertEqual(self.calculator.compute_functional_index(['Z', 'E']),
                         self.calculator.compute_functional_index(['E', 'Z']))

        # J group
        self.assertEqual(self.calculator.compute_functional_index(['J', 'L']),
                         self.calculator.compute_functional_index(['L', 'J']))
        self.assertEqual(self.calculator.compute_functional_index(['J', 'I']),
                         self.calculator.compute_functional_index(['I', 'J']))

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


class ConvertSequencesToBinaryArrayTestCase(unittest.TestCase):
    def setUp(self):
        s1 = 'A-N-'  # 0101
        s2 = 'AR-D'  # 0010
        s3 = 'AR-D'  # 0010
        s4 = 'A--D'  # 0110

        self.gene_sequences = [s1, s2, s3, s4]

        self.calculator = homogeneityindex.HomogeneityCalculator(quick_homogeneity=False)

    def test_by_residue(self):
        actual = self.calculator.convert_sequences_to_binary_array(self.gene_sequences,
                                                                   bygene=False)

        # Each entry represents a column in the alignment.
        expected = [0b0000, 0b1001, 0b0111, 0b1000]

        self.assertEqual(actual, expected)

    def test_by_gene(self):
        actual = self.calculator.convert_sequences_to_binary_array(self.gene_sequences,
                                                                   bygene=True)

        # Each entry represents a sequence in the alignment.
        expected = [0b0101, 0b0010, 0b0010, 0b0110]

        self.assertEqual(actual, expected)


class ComputeGeometrixIndexTestCase(unittest.TestCase):
    def setUp(self):
        s1 = 'A-N--'  # 01011
        s2 = 'AR-D-'  # 00101
        s3 = 'AR-D-'  # 00101
        s4 = 'A--D-'  # 01101

        self.gene_sequences = [s1, s2, s3, s4]

        self.calculator = homogeneityindex.HomogeneityCalculator(quick_homogeneity=False)
        
        # For the geometric homogeneity index.  This looks like a lot, but it is good to double check that the 
        # by-hand calculation matches the code, and to list it out explicitly so it is clearer how the algorithm 
        # works. 
        num_seqs = len(self.gene_sequences)
        num_aln_cols = len(self.gene_sequences[0])

        max_similarities_per_aln_col = num_seqs
        max_similarities_per_seq = num_aln_cols

        num_comparisons_per_aln_col = num_aln_cols - 1
        num_comparisons_per_seq = num_seqs - 1

        # First do the pairwise comparisons....

        # 1v2, 1v3, 1v4, (1v5 for by residue)
        aln_col1_similarities = [2, 1, 3, 0]
        seq1_similarities = [2, 2, 3]

        # 2v1, 2v3, 2v4, (2v5 for by residue)
        aln_col2_similarities = [2, 1, 3, 2]
        seq2_similarities = [2, 5, 4]

        # 3v1, 3v2, 3v4, (3v5 for by residue)
        aln_col3_similarities = [1, 1, 0, 3]
        seq3_similarities = [2, 5, 4]

        # 4v1, 4v2, 4v3, (4v5 for by residue)
        aln_col4_similarities = [3, 3, 0, 1]
        seq4_similarities = [3, 4, 4]

        # 5v1, 5v2, 5v3, 5v4
        aln_col5_similarities = [0, 2, 3, 1]

        # Then each column has a similarity score w.r.t. the other columns.
        aln_col1_similarity_score = sum(aln_col1_similarities) / max_similarities_per_aln_col / num_comparisons_per_aln_col
        aln_col2_similarity_score = sum(aln_col2_similarities) / max_similarities_per_aln_col / num_comparisons_per_aln_col
        aln_col3_similarity_score = sum(aln_col3_similarities) / max_similarities_per_aln_col / num_comparisons_per_aln_col
        aln_col4_similarity_score = sum(aln_col4_similarities) / max_similarities_per_aln_col / num_comparisons_per_aln_col
        aln_col5_similarity_score = sum(aln_col5_similarities) / max_similarities_per_aln_col / num_comparisons_per_aln_col

        # Also, each seq has a similarity score w.r.t. the other sequences.
        seq1_similarity_score = sum(seq1_similarities) / max_similarities_per_seq / num_comparisons_per_seq
        seq2_similarity_score = sum(seq2_similarities) / max_similarities_per_seq / num_comparisons_per_seq
        seq3_similarity_score = sum(seq3_similarities) / max_similarities_per_seq / num_comparisons_per_seq
        seq4_similarity_score = sum(seq4_similarities) / max_similarities_per_seq / num_comparisons_per_seq

        # The quick geo score is the mean of all alignment column similarity scores.
        self.quick_geometric_similarity = sum([aln_col1_similarity_score,
                                               aln_col2_similarity_score,
                                               aln_col3_similarity_score,
                                               aln_col4_similarity_score,
                                               aln_col5_similarity_score]) / num_aln_cols

        # The by sequence similarity score is the mean of all sequence similarity scores.
        by_seq_similarity = sum([seq1_similarity_score,
                                 seq2_similarity_score,
                                 seq3_similarity_score,
                                 seq4_similarity_score]) / num_seqs

        # Finally the full geometric similarity score is the average of the quick score (by residue) and the by
        # sequence score.
        self.full_geometric_similarity = (self.quick_geometric_similarity + by_seq_similarity) / 2

    def test_quick_geometric_homogeneity_index(self):
        actual = self.calculator.compute_geometric_index(self.gene_sequences,
                                                         quick_homogeneity=True)

        self.assertEqual(actual, self.quick_geometric_similarity)

    def test_geometric_homogeneity_index(self):
        actual = self.calculator.compute_geometric_index(self.gene_sequences,
                                                         quick_homogeneity=False)

        self.assertEqual(actual, self.full_geometric_similarity)

    def test_that_quick_geo_index_is_invariant_to_sequence_order(self):
        aln1 = ['A-A', 'AA-', '-AA', '-A-']
        aln2 = ['AA-', 'A-A', '-AA', '-A-']

        self.assertEqual(self.calculator.compute_geometric_index(aln1, quick_homogeneity=True),
                         self.calculator.compute_geometric_index(aln2, quick_homogeneity=True))

    def test_that_full_geo_index_is_invariant_to_sequence_order(self):
        aln1 = ['A-A', 'AA-', '-AA', '-A-']
        aln2 = ['AA-', 'A-A', '-AA', '-A-']

        self.assertEqual(self.calculator.compute_geometric_index(aln1, quick_homogeneity=False),
                         self.calculator.compute_geometric_index(aln2, quick_homogeneity=False))


if __name__ == '__main__':
    unittest.main()
