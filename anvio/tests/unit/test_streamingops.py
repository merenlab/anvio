# -*- coding: utf-8
# pylint: disable=line-too-long
"""
Unit tests for streaming accumulators in streamingops.py.
"""

import unittest
import numpy as np
from collections import OrderedDict

import anvio
import anvio.streamingops as streamingops

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Florian Trigodet"
__email__ = ""


class TestCoverageAccumulator(unittest.TestCase):
    """Test the CoverageAccumulator class."""

    def test_initialization(self):
        """Test that accumulator initializes with correct size."""
        acc = streamingops.CoverageAccumulator(1000)
        self.assertEqual(len(acc.coverage), 1000)
        self.assertEqual(acc.num_reads, 0)
        self.assertTrue(np.all(acc.coverage == 0))

    def test_single_read_coverage(self):
        """Test coverage update from a single read."""
        acc = streamingops.CoverageAccumulator(100)
        acc.update([(10, 30)])  # 20bp read at position 10-30

        # Positions 10-29 should have coverage 1
        self.assertTrue(np.all(acc.coverage[10:30] == 1))
        # Other positions should be 0
        self.assertTrue(np.all(acc.coverage[:10] == 0))
        self.assertTrue(np.all(acc.coverage[30:] == 0))
        self.assertEqual(acc.num_reads, 1)

    def test_overlapping_reads(self):
        """Test that overlapping reads accumulate coverage."""
        acc = streamingops.CoverageAccumulator(100)
        acc.update([(10, 30)])
        acc.update([(20, 40)])

        # Positions 10-19: coverage 1
        self.assertTrue(np.all(acc.coverage[10:20] == 1))
        # Positions 20-29: coverage 2 (overlap)
        self.assertTrue(np.all(acc.coverage[20:30] == 2))
        # Positions 30-39: coverage 1
        self.assertTrue(np.all(acc.coverage[30:40] == 1))
        self.assertEqual(acc.num_reads, 2)

    def test_multiple_blocks_per_read(self):
        """Test read with multiple coverage blocks (e.g., from deletions)."""
        acc = streamingops.CoverageAccumulator(100)
        # Read with two blocks: 10-20 and 30-40 (gap from 20-30, e.g., deletion)
        acc.update([(10, 20), (30, 40)])

        self.assertTrue(np.all(acc.coverage[10:20] == 1))
        self.assertTrue(np.all(acc.coverage[20:30] == 0))  # Gap
        self.assertTrue(np.all(acc.coverage[30:40] == 1))

    def test_boundary_clipping(self):
        """Test that reads extending beyond split boundaries are clipped."""
        acc = streamingops.CoverageAccumulator(50)
        # Read that extends beyond both boundaries
        acc.update([(-10, 60)])

        # Should be clipped to [0, 50)
        self.assertTrue(np.all(acc.coverage == 1))

    def test_finalize(self):
        """Test finalize returns correct statistics."""
        acc = streamingops.CoverageAccumulator(10)
        acc.update([(0, 5)])   # First 5 positions covered
        acc.update([(0, 10)])  # All 10 positions covered

        result = acc.finalize()

        self.assertEqual(result['num_reads'], 2)
        self.assertEqual(result['min'], 1)
        self.assertEqual(result['max'], 2)
        self.assertEqual(result['detection'], 1.0)  # All positions covered
        self.assertAlmostEqual(result['mean'], 1.5)  # (2+2+2+2+2+1+1+1+1+1)/10


class TestSNVAccumulator(unittest.TestCase):
    """Test the SNVAccumulator class."""

    def test_initialization(self):
        """Test that accumulator initializes with correct shape."""
        acc = streamingops.SNVAccumulator(100)
        self.assertEqual(acc.allele_counts.shape, (5, 100))
        self.assertTrue(np.all(acc.allele_counts == 0))

    def test_update_from_vectorized(self):
        """Test SNV accumulation from a vectorized read."""
        acc = streamingops.SNVAccumulator(100)

        # Create a mock vectorized array
        # Columns: [ref_pos, query_seq, mapping_type, ref_seq]
        # 5 mapped bases at positions 10-14
        vectorized = np.array([
            [10, ord('A'), 0, ord('A')],
            [11, ord('C'), 0, ord('C')],
            [12, ord('G'), 0, ord('T')],  # Mismatch: G vs T
            [13, ord('T'), 0, ord('T')],
            [14, ord('A'), 0, ord('A')],
        ], dtype=np.int32)

        acc.update_from_vectorized(vectorized, split_start=0)

        # Check that allele counts are updated
        # Position 10: A (index 0)
        self.assertEqual(acc.allele_counts[0, 10], 1)
        # Position 11: C (index 1)
        self.assertEqual(acc.allele_counts[1, 11], 1)
        # Position 12: G (index 2) - the query base, not reference
        self.assertEqual(acc.allele_counts[2, 12], 1)
        # Position 13: T (index 3)
        self.assertEqual(acc.allele_counts[3, 13], 1)
        # Position 14: A (index 0)
        self.assertEqual(acc.allele_counts[0, 14], 1)

    def test_update_filters_non_mapped(self):
        """Test that only mapped positions (type 0) are counted."""
        acc = streamingops.SNVAccumulator(100)

        # Create a vectorized array with insertion (type 1)
        vectorized = np.array([
            [10, ord('A'), 0, ord('A')],  # Mapped
            [10, ord('T'), 1, -1],        # Insertion (should be ignored)
            [10, ord('G'), 1, -1],        # Insertion (should be ignored)
            [11, ord('C'), 0, ord('C')],  # Mapped
        ], dtype=np.int32)

        acc.update_from_vectorized(vectorized, split_start=0)

        # Position 10 should only have 1 A (the insertion bases ignored)
        self.assertEqual(acc.allele_counts[0, 10], 1)  # A
        self.assertEqual(acc.allele_counts[3, 10], 0)  # T (from insertion, should be 0)

    def test_split_start_offset(self):
        """Test that split_start correctly offsets positions."""
        acc = streamingops.SNVAccumulator(100)

        vectorized = np.array([
            [110, ord('A'), 0, ord('A')],
            [111, ord('C'), 0, ord('C')],
        ], dtype=np.int32)

        acc.update_from_vectorized(vectorized, split_start=100)

        # Positions should be 10 and 11 (110-100 and 111-100)
        self.assertEqual(acc.allele_counts[0, 10], 1)  # A at position 10
        self.assertEqual(acc.allele_counts[1, 11], 1)  # C at position 11

    def test_get_coverage(self):
        """Test that get_coverage returns sum across alleles."""
        acc = streamingops.SNVAccumulator(100)

        # Add two reads at position 10
        vectorized1 = np.array([[10, ord('A'), 0, ord('A')]], dtype=np.int32)
        vectorized2 = np.array([[10, ord('C'), 0, ord('A')]], dtype=np.int32)

        acc.update_from_vectorized(vectorized1, split_start=0)
        acc.update_from_vectorized(vectorized2, split_start=0)

        coverage = acc.get_coverage()
        self.assertEqual(coverage[10], 2)  # Two reads at position 10


class TestINDELAccumulator(unittest.TestCase):
    """Test the INDELAccumulator class."""

    def test_initialization(self):
        """Test that accumulator initializes correctly."""
        acc = streamingops.INDELAccumulator(100)
        self.assertEqual(len(acc.indels), 0)

    def test_insertion_accumulation(self):
        """Test that insertions are accumulated correctly."""
        acc = streamingops.INDELAccumulator(100)

        # Add an insertion at position 10
        insertions = [(10, np.array([ord('A'), ord('T'), ord('G')]))]
        deletions = []

        acc.update(insertions, deletions, 'test_split', 'A' * 100)

        self.assertEqual(len(acc.indels), 1)
        indel = list(acc.indels.values())[0]
        self.assertEqual(indel['type'], 'INS')
        self.assertEqual(indel['pos'], 10)
        self.assertEqual(indel['sequence'], 'ATG')
        self.assertEqual(indel['length'], 3)
        self.assertEqual(indel['count'], 1)

    def test_deletion_accumulation(self):
        """Test that deletions are accumulated correctly."""
        acc = streamingops.INDELAccumulator(100)

        insertions = []
        deletions = [(20, 5)]  # 5bp deletion at position 20

        acc.update(insertions, deletions, 'test_split', 'A' * 100)

        self.assertEqual(len(acc.indels), 1)
        indel = list(acc.indels.values())[0]
        self.assertEqual(indel['type'], 'DEL')
        self.assertEqual(indel['pos'], 20)
        self.assertEqual(indel['length'], 5)
        self.assertEqual(indel['count'], 1)

    def test_duplicate_indel_counting(self):
        """Test that identical indels increment count."""
        acc = streamingops.INDELAccumulator(100)

        insertions = [(10, np.array([ord('A'), ord('T')]))]
        deletions = []

        # Add same insertion 3 times
        acc.update(insertions, deletions, 'test_split', 'A' * 100)
        acc.update(insertions, deletions, 'test_split', 'A' * 100)
        acc.update(insertions, deletions, 'test_split', 'A' * 100)

        self.assertEqual(len(acc.indels), 1)
        indel = list(acc.indels.values())[0]
        self.assertEqual(indel['count'], 3)

    def test_different_indels_tracked_separately(self):
        """Test that different indels are tracked separately."""
        acc = streamingops.INDELAccumulator(100)

        # Two different insertions
        acc.update([(10, np.array([ord('A')]))], [], 'test_split', 'A' * 100)
        acc.update([(10, np.array([ord('T')]))], [], 'test_split', 'A' * 100)
        acc.update([(20, np.array([ord('A')]))], [], 'test_split', 'A' * 100)

        self.assertEqual(len(acc.indels), 3)

    def test_out_of_bounds_ignored(self):
        """Test that indels outside split boundaries are ignored."""
        acc = streamingops.INDELAccumulator(50)

        insertions = [(-5, np.array([ord('A')])), (60, np.array([ord('T')]))]
        deletions = []

        acc.update(insertions, deletions, 'test_split', 'A' * 50)

        self.assertEqual(len(acc.indels), 0)


if __name__ == '__main__':
    unittest.main()
