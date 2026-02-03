# -*- coding: utf-8
# pylint: disable=line-too-long
"""
Unit tests for the Read class in bamops.py, specifically the new
extract_variant_evidence() method for single-pass variant extraction.
"""

import unittest
from unittest.mock import MagicMock
import numpy as np

import anvio
import anvio.bamops as bamops

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Florian Trigodet"
__email__ = ""


def create_mock_pysam_read(cigartuples, query_sequence, reference_start, reference_sequence=None):
    """Create a mock pysam.AlignedSegment for testing.

    Parameters
    ==========
    cigartuples : list of tuples
        CIGAR operations as (operation, length) tuples.
        Operations: 0=M (match/mismatch), 1=I (insertion), 2=D (deletion),
                    3=N (skip), 4=S (soft clip), 5=H (hard clip)

    query_sequence : str
        The read sequence (what was sequenced)

    reference_start : int
        0-based reference position where alignment starts

    reference_sequence : str, optional
        The reference sequence aligned to. If None, will be derived from
        query_sequence (treating all as matches).
    """
    mock_read = MagicMock()
    mock_read.cigartuples = cigartuples
    mock_read.query_sequence = query_sequence
    mock_read.reference_start = reference_start

    # Calculate reference_end from cigar
    ref_consumed = 0
    for op, length in cigartuples:
        if op in (0, 2, 3):  # M, D, N consume reference
            ref_consumed += length
    mock_read.reference_end = reference_start + ref_consumed

    # Handle reference sequence
    if reference_sequence is not None:
        mock_read.has_tag = lambda x: x == 'MD'
        mock_read.get_reference_sequence = lambda: reference_sequence
    else:
        mock_read.has_tag = lambda x: False

    return mock_read


class TestExtractVariantEvidenceBasic(unittest.TestCase):
    """Test basic functionality of extract_variant_evidence()"""

    def test_perfect_match_no_variants(self):
        """A read that perfectly matches the reference should have no SNV evidence."""
        # 10bp read, perfect match
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 10)],  # 10M
            query_sequence='ACGTACGTAC',
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=0)

        # Should have one coverage block
        self.assertEqual(len(result['coverage_blocks']), 1)
        self.assertEqual(result['coverage_blocks'][0], (100, 110))

        # No SNVs (perfect match)
        self.assertEqual(len(result['snv_evidence']), 0)

        # No indels
        self.assertEqual(len(result['insertions']), 0)
        self.assertEqual(len(result['deletions']), 0)

    def test_snv_detection(self):
        """Test that SNVs (mismatches) are correctly detected."""
        # 10bp read with 2 mismatches at positions 102 and 108
        # query:     A C T T A C G T G C
        # reference: A C G T A C G T A C
        # position:  0 1 2 3 4 5 6 7 8 9  (relative to ref_start)
        # Mismatches at position 2 (T vs G) and position 8 (G vs A)
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 10)],  # 10M
            query_sequence='ACTTACGTGC',
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=100)

        # Two SNVs
        self.assertEqual(len(result['snv_evidence']), 2)

        # Check positions (relative to split_start=100)
        snv_positions = [snv[0] for snv in result['snv_evidence']]
        self.assertIn(2, snv_positions)  # 102 - 100 = 2
        self.assertIn(8, snv_positions)  # 108 - 100 = 8

    def test_insertion_detection(self):
        """Test that insertions are correctly detected."""
        # Read with a 3bp insertion at position 105
        # Reference: ACGTA---CGTAC (positions 100-109, 10bp)
        # Read:      ACGTATTACGTAC (13bp)
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 5), (1, 3), (0, 5)],  # 5M3I5M
            query_sequence='ACGTATTACGTAC',
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=100)

        # One insertion
        self.assertEqual(len(result['insertions']), 1)
        ins_pos, ins_seq = result['insertions'][0]
        self.assertEqual(ins_pos, 4)  # Position 104 relative to split_start, insertion reported at last match position
        self.assertEqual(len(ins_seq), 3)  # 3bp insertion
        # Check the inserted sequence is 'TTA' (as ord values)
        self.assertEqual(chr(ins_seq[0]), 'T')
        self.assertEqual(chr(ins_seq[1]), 'T')
        self.assertEqual(chr(ins_seq[2]), 'A')

    def test_deletion_detection(self):
        """Test that deletions are correctly detected."""
        # Read with a 3bp deletion at position 105-107
        # Reference: ACGTACGTAC (positions 100-109, 10bp)
        # Read:      ACGTA--TAC (7bp aligned)
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 5), (2, 3), (0, 2)],  # 5M3D2M
            query_sequence='ACGTATAC',
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=100)

        # One deletion
        self.assertEqual(len(result['deletions']), 1)
        del_pos, del_len = result['deletions'][0]
        self.assertEqual(del_pos, 5)  # Position 105 relative to split_start
        self.assertEqual(del_len, 3)  # 3bp deletion

    def test_complex_cigar(self):
        """Test a complex CIGAR with matches, mismatches, insertion, and deletion."""
        # CIGAR: 3M1I2M2D3M = 3 match, 1 insert, 2 match, 2 delete, 3 match
        # Reference positions: 100-102 (3M), 103-104 (2M), 105-106 (2D), 107-109 (3M)
        # Total ref consumed: 3 + 2 + 2 + 3 = 10
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 3), (1, 1), (0, 2), (2, 2), (0, 3)],
            query_sequence='ACGXTAGTA',  # 3 + 1 + 2 + 3 = 9bp (X is the insertion)
            reference_start=100,
            reference_sequence='ACGTAGGGTA'  # Reference at aligned positions
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=100)

        # Should have coverage blocks (non-contiguous due to deletion)
        self.assertEqual(len(result['coverage_blocks']), 3)

        # One insertion
        self.assertEqual(len(result['insertions']), 1)

        # One deletion
        self.assertEqual(len(result['deletions']), 1)
        del_pos, del_len = result['deletions'][0]
        self.assertEqual(del_len, 2)


class TestExtractVariantEvidenceSplitRelative(unittest.TestCase):
    """Test that positions are correctly relative to split_start."""

    def test_split_start_offset(self):
        """Positions should be relative to split_start."""
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 10)],
            query_sequence='ACTTACGTAC',  # Mismatch at position 102 (T vs G)
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)

        # With split_start=0, positions are absolute
        result0 = read.extract_variant_evidence(split_start=0)
        self.assertEqual(result0['coverage_blocks'][0], (100, 110))
        self.assertEqual(result0['snv_evidence'][0][0], 102)  # Absolute position

        # With split_start=100, positions are relative to split
        result100 = read.extract_variant_evidence(split_start=100)
        self.assertEqual(result100['coverage_blocks'][0], (0, 10))
        self.assertEqual(result100['snv_evidence'][0][0], 2)  # Relative position


class TestExtractVariantEvidenceSkipFlags(unittest.TestCase):
    """Test the skip_SNV and skip_INDEL flags."""

    def test_skip_snv(self):
        """When skip_SNV=True, snv_evidence should be empty."""
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 10)],
            query_sequence='ACTTACGTAC',  # Has a mismatch
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=0, skip_SNV=True)

        # SNV evidence should be empty
        self.assertEqual(len(result['snv_evidence']), 0)

        # Coverage should still be computed
        self.assertEqual(len(result['coverage_blocks']), 1)

    def test_skip_indel(self):
        """When skip_INDEL=True, insertions and deletions should be empty."""
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 5), (1, 2), (0, 5)],  # Has an insertion
            query_sequence='ACGTATTCGTAC',
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=0, skip_INDEL=True)

        # INDEL evidence should be empty
        self.assertEqual(len(result['insertions']), 0)
        self.assertEqual(len(result['deletions']), 0)

        # Coverage and SNVs should still be computed
        self.assertGreater(len(result['coverage_blocks']), 0)


class TestExtractVariantEvidenceVectorized(unittest.TestCase):
    """Test that the vectorized array is included for SCV processing."""

    def test_vectorized_array_included(self):
        """The result should include the vectorized array."""
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 10)],
            query_sequence='ACGTACGTAC',
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)
        result = read.extract_variant_evidence(split_start=0)

        # Vectorized array should be present
        self.assertIn('vectorized', result)
        self.assertIsNotNone(result['vectorized'])

        # Should be a 2D numpy array with 4 columns
        self.assertEqual(result['vectorized'].ndim, 2)
        self.assertEqual(result['vectorized'].shape[1], 4)


class TestExtractVariantEvidenceCompareToOldMethods(unittest.TestCase):
    """Test that extract_variant_evidence produces equivalent results to old methods."""

    def test_coverage_blocks_match_get_blocks(self):
        """Coverage blocks should match what get_blocks() returns."""
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 5), (1, 2), (0, 3), (2, 2), (0, 4)],  # Complex CIGAR
            query_sequence='ACGTATTCGTACGT',
            reference_start=100,
            reference_sequence='ACGTACGTACGTAC'
        )

        read = bamops.Read(mock_read)

        # Get coverage blocks from new method
        result = read.extract_variant_evidence(split_start=0)
        new_blocks = result['coverage_blocks']

        # Get coverage blocks from old method
        old_blocks = read.get_blocks()

        # Should be identical
        self.assertEqual(new_blocks, old_blocks)

    def test_snv_positions_match_aligned_sequence_method(self):
        """SNV positions should match those found via get_aligned_sequence_and_reference_positions."""
        mock_read = create_mock_pysam_read(
            cigartuples=[(0, 10)],
            query_sequence='ACTTACGTGC',  # Mismatches at positions 2 and 7
            reference_start=100,
            reference_sequence='ACGTACGTAC'
        )

        read = bamops.Read(mock_read)

        # Get SNVs from new method
        result = read.extract_variant_evidence(split_start=0)
        new_snv_positions = set(snv[0] for snv in result['snv_evidence'])

        # Get SNVs using old approach
        aligned_seq, ref_positions = read.get_aligned_sequence_and_reference_positions()
        ref_seq_ord = read.reference_sequence

        # Find mismatches manually (mimicking old SNV detection)
        old_snv_positions = set()
        for i, (pos, query_ord) in enumerate(zip(ref_positions, aligned_seq)):
            # The reference sequence index needs to be calculated relative to alignment
            ref_idx = pos - read.reference_start
            if ref_idx < len(ref_seq_ord) and query_ord != ref_seq_ord[ref_idx]:
                old_snv_positions.add(int(pos))

        # Should match
        self.assertEqual(new_snv_positions, old_snv_positions)


if __name__ == '__main__':
    unittest.main()
