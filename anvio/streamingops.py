# -*- coding: utf-8
# pylint: disable=line-too-long
"""
Streaming accumulators for the profiler redesign.

These classes accumulate variant evidence from reads processed via
Read.extract_variant_evidence() and produce output compatible with
the existing profile database tables.
"""

import numpy as np
from collections import OrderedDict

import anvio
import anvio.constants as constants
import anvio.terminal as terminal

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Florian Trigodet"
__email__ = ""


run = terminal.Run()
progress = terminal.Progress()


class CoverageAccumulator:
    """Accumulates coverage data for a split.

    This accumulator tracks per-position coverage counts within a split.
    Coverage is incremented for each aligned base from reads.

    Parameters
    ==========
    split_length : int
        Length of the split in base pairs.

    Attributes
    ==========
    coverage : numpy.ndarray
        Array of coverage counts, one per position in the split.
    num_reads : int
        Number of reads processed.

    Notes
    =====
    - Memory usage: O(split_length) - bounded by split size (~20KB default)
    - Call update() for each read's coverage blocks
    - Call finalize() to get the final coverage array and statistics
    """

    def __init__(self, split_length):
        self.split_length = split_length
        self.coverage = np.zeros(split_length, dtype=np.int32)
        self.num_reads = 0

    def update(self, coverage_blocks):
        """Update coverage with blocks from a read.

        Parameters
        ==========
        coverage_blocks : list of tuples
            List of (start, end) tuples representing covered regions.
            Positions should be relative to split start (0-indexed within split).
        """
        for start, end in coverage_blocks:
            # Clip to split boundaries (reads may extend beyond split)
            start = max(0, start)
            end = min(self.split_length, end)
            if start < end:
                self.coverage[start:end] += 1
        self.num_reads += 1

    def finalize(self):
        """Finalize and return coverage data.

        Returns
        =======
        dict with keys:
            'coverage': numpy array of per-position coverage
            'num_reads': number of reads processed
            'mean': mean coverage
            'std': standard deviation of coverage
            'min': minimum coverage
            'max': maximum coverage
            'median': median coverage
            'detection': fraction of positions with coverage > 0
        """
        c = self.coverage
        return {
            'coverage': c,
            'num_reads': self.num_reads,
            'mean': np.mean(c) if len(c) > 0 else 0.0,
            'std': np.std(c) if len(c) > 0 else 0.0,
            'min': np.min(c) if len(c) > 0 else 0,
            'max': np.max(c) if len(c) > 0 else 0,
            'median': np.median(c) if len(c) > 0 else 0.0,
            'detection': np.sum(c > 0) / len(c) if len(c) > 0 else 0.0,
        }


class SNVAccumulator:
    """Accumulates SNV (single nucleotide variant) evidence for a split.

    This accumulator tracks per-position allele counts within a split.
    It uses a 2D array with dimensions (5 nucleotides, split_length).

    Parameters
    ==========
    split_length : int
        Length of the split in base pairs.

    Attributes
    ==========
    allele_counts : numpy.ndarray
        2D array of shape (5, split_length) tracking counts for A, C, G, T, N.

    Notes
    =====
    - Memory usage: O(5 * split_length) - bounded by split size
    - Call update() for each read's SNV evidence
    - Call finalize() with additional context to get SNV profiles
    """

    # Nucleotide to array index mapping (matches constants.nucleotides order)
    NT_TO_INDEX = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3, ord('N'): 4,
                   ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3, ord('n'): 4}

    def __init__(self, split_length):
        self.split_length = split_length
        # 5 rows: A, C, G, T, N (matching constants.nucleotides)
        self.allele_counts = np.zeros((5, split_length), dtype=np.float64)

    def update(self, snv_evidence, coverage_blocks):
        """Update allele counts with evidence from a read.

        Parameters
        ==========
        snv_evidence : list of tuples
            List of (pos, ref_ord, query_ord) tuples for mismatches.
            Positions should be relative to split start.

        coverage_blocks : list of tuples
            List of (start, end) tuples for coverage regions.
            Used to count reference alleles at non-mismatch positions.
        """
        # First, increment reference alleles for all covered positions
        # This is implicit - we track all observed alleles
        # Actually, for efficiency, we should track all aligned bases, not just mismatches
        # But the current design only gives us mismatches...

        # For now, we'll need to update this to work with the full aligned sequence
        # This accumulator expects the calling code to pass ALL aligned bases, not just mismatches
        pass

    def update_from_aligned_pairs(self, aligned_positions, aligned_bases_ord):
        """Update allele counts from aligned position/base pairs.

        Parameters
        ==========
        aligned_positions : array-like
            Reference positions (relative to split start) for each aligned base.

        aligned_bases_ord : array-like
            Ordinal values of the aligned bases (query sequence).
        """
        for pos, base_ord in zip(aligned_positions, aligned_bases_ord):
            if 0 <= pos < self.split_length:
                idx = self.NT_TO_INDEX.get(base_ord)
                if idx is not None:
                    self.allele_counts[idx, pos] += 1

    def update_from_vectorized(self, vectorized, split_start):
        """Update allele counts from a vectorized read array.

        Parameters
        ==========
        vectorized : numpy.ndarray
            The vectorized read array from Read.vectorize().
            Columns: [ref_pos, query_seq, mapping_type, ref_seq]

        split_start : int
            The reference position of the split start.
        """
        # Filter to only mapped positions (mapping_type == 0)
        mapped_mask = vectorized[:, 2] == 0
        mapped = vectorized[mapped_mask]

        for row in mapped:
            pos = int(row[0]) - split_start
            base_ord = int(row[1])
            if 0 <= pos < self.split_length:
                idx = self.NT_TO_INDEX.get(base_ord)
                if idx is not None:
                    self.allele_counts[idx, pos] += 1

    def get_coverage(self):
        """Get per-position coverage from allele counts."""
        return self.allele_counts.sum(axis=0)

    def finalize(self):
        """Finalize and return the allele counts array.

        Returns
        =======
        numpy.ndarray
            The allele counts array of shape (5, split_length).

        Notes
        =====
        The returned array is used by ProcessNucleotideCounts to compute
        the final SNV profiles with filtering and statistics.
        """
        return self.allele_counts


class INDELAccumulator:
    """Accumulates INDEL (insertion/deletion) evidence for a split.

    This accumulator tracks indel events by position and sequence/length.

    Parameters
    ==========
    split_length : int
        Length of the split in base pairs.

    Attributes
    ==========
    indels : dict
        Dictionary mapping indel hash to indel entry (OrderedDict).

    Notes
    =====
    - Memory usage: O(number of unique indels) - typically small
    - Call update() for each read's insertions and deletions
    - Call finalize() to get the indel profiles dict
    """

    def __init__(self, split_length):
        self.split_length = split_length
        self.indels = {}

    def update(self, insertions, deletions, split_name, split_sequence,
               per_position_info=None):
        """Update indel counts with evidence from a read.

        Parameters
        ==========
        insertions : list of tuples
            List of (pos, sequence_as_ords) tuples.

        deletions : list of tuples
            List of (pos, length) tuples.

        split_name : str
            Name of the split (for output).

        split_sequence : str
            The split's reference sequence (for getting reference base).

        per_position_info : dict, optional
            Per-position gene annotation info. If provided, indel entries
            will include gene-related fields.
        """
        # Process insertions
        for ins_pos, ins_seq_ord in insertions:
            if ins_pos < 0 or ins_pos >= self.split_length:
                continue

            ins_seq = ''.join(chr(x) for x in ins_seq_ord)
            indel_hash = hash((ins_pos, ins_seq))

            if indel_hash in self.indels:
                self.indels[indel_hash]['count'] += 1
            else:
                entry = self._create_indel_entry(
                    'INS', ins_pos, ins_seq, len(ins_seq),
                    split_name, split_sequence, per_position_info
                )
                self.indels[indel_hash] = entry

        # Process deletions
        for del_pos, del_len in deletions:
            if del_pos < 0 or del_pos >= self.split_length:
                continue

            indel_hash = hash((del_pos, del_len))

            if indel_hash in self.indels:
                self.indels[indel_hash]['count'] += 1
            else:
                entry = self._create_indel_entry(
                    'DEL', del_pos, '', del_len,
                    split_name, split_sequence, per_position_info
                )
                self.indels[indel_hash] = entry

    def _create_indel_entry(self, indel_type, pos, sequence, length,
                            split_name, split_sequence, per_position_info):
        """Create an indel entry OrderedDict."""
        entry = OrderedDict([
            ('split_name', split_name),
            ('pos', pos),
            ('pos_in_contig', None),  # Will be set later with split.start
            ('corresponding_gene_call', -1),
            ('in_noncoding_gene_call', 0),
            ('in_coding_gene_call', 0),
            ('base_pos_in_codon', 0),
            ('codon_order_in_gene', -1),
            ('cov_outlier_in_split', 0),
            ('cov_outlier_in_contig', 0),
            ('reference', split_sequence[pos] if pos < len(split_sequence) else 'N'),
            ('type', indel_type),
            ('sequence', sequence),
            ('length', length),
            ('count', 1),
        ])

        # Add per-position info if available
        if per_position_info is not None and pos < self.split_length:
            for key in ['corresponding_gene_call', 'in_noncoding_gene_call',
                        'in_coding_gene_call', 'base_pos_in_codon', 'codon_order_in_gene']:
                if key in per_position_info:
                    entry[key] = per_position_info[key][pos]

        return entry

    def finalize(self):
        """Finalize and return the indels dictionary.

        Returns
        =======
        dict
            Dictionary mapping indel hash to indel entry (OrderedDict).
        """
        return self.indels


class SCVAccumulator:
    """Accumulates SCV (single codon variant) evidence for genes in a split.

    This accumulator tracks per-codon allele counts for each gene that
    overlaps the split and has SNVs.

    Parameters
    ==========
    split_length : int
        Length of the split in base pairs.

    Attributes
    ==========
    gene_allele_counts : dict
        Dictionary mapping gene_id to codon allele counts array.
    gene_calls : dict
        Dictionary mapping gene_id to gene call info.
    reference_codon_sequences : dict
        Dictionary mapping gene_id to reference codon sequence.

    Notes
    =====
    - Memory usage: O(number of active genes * codons per gene)
    - This is the most complex accumulator due to gene-aware processing
    - Call update_from_vectorized() for each read
    - Call finalize() to get SCV profiles per gene
    """

    # Codon to array index (64 codons)
    CDN_TO_INDEX = {codon: i for i, codon in enumerate(constants.codons)}

    def __init__(self, split_length, per_position_info=None, genes_with_snvs=None):
        """Initialize SCV accumulator.

        Parameters
        ==========
        split_length : int
            Length of the split.

        per_position_info : dict, optional
            Per-position gene annotation info with keys like
            'corresponding_gene_call', 'codon_order_in_gene', etc.

        genes_with_snvs : set, optional
            Set of gene IDs that have SNVs. Only these genes will be tracked.
        """
        self.split_length = split_length
        self.per_position_info = per_position_info or {}
        self.genes_with_snvs = genes_with_snvs or set()

        self.gene_allele_counts = {}
        self.gene_calls = {}
        self.reference_codon_sequences = {}

    def set_genes_with_snvs(self, genes_with_snvs):
        """Set the genes that have SNVs (for filtering)."""
        self.genes_with_snvs = genes_with_snvs

    def update_from_vectorized(self, vectorized, split_start):
        """Update codon counts from a vectorized read array.

        This method extracts codons from the read and updates the
        per-gene codon allele counts.

        Parameters
        ==========
        vectorized : numpy.ndarray
            The vectorized read array from Read.vectorize().

        split_start : int
            The reference position of the split start.

        Notes
        =====
        This is a complex operation that handles:
        - Determining which genes the read overlaps
        - Extracting gapless segments
        - Converting to codon sequences
        - Handling forward/reverse genes
        - Filtering to codons that fully covered
        """
        if not self.per_position_info or 'corresponding_gene_call' not in self.per_position_info:
            return

        # Get the reference positions covered by this read
        mapped_mask = vectorized[:, 2] == 0
        if not np.any(mapped_mask):
            return

        ref_start = int(vectorized[mapped_mask, 0].min())
        ref_end = int(vectorized[mapped_mask, 0].max()) + 1

        # Convert to split-relative positions
        split_rel_start = ref_start - split_start
        split_rel_end = ref_end - split_start

        # Clip to split boundaries
        split_rel_start = max(0, split_rel_start)
        split_rel_end = min(self.split_length, split_rel_end)

        if split_rel_start >= split_rel_end:
            return

        # Find genes overlapped by this read
        gene_ids_in_read = self.per_position_info['corresponding_gene_call'][split_rel_start:split_rel_end]
        genes_in_read = set(gene_ids_in_read)

        for gene_id in genes_in_read:
            if gene_id == -1:
                continue
            if gene_id not in self.genes_with_snvs:
                continue

            # Complex codon extraction logic would go here
            # This is a simplified placeholder - full implementation would
            # mirror the logic in contigops.Auxiliary.run_SCVs
            pass

    def finalize(self):
        """Finalize and return SCV profiles per gene.

        Returns
        =======
        dict
            Dictionary mapping gene_id to SCV profile dict.
        """
        return {
            'gene_allele_counts': self.gene_allele_counts,
            'gene_calls': self.gene_calls,
            'reference_codon_sequences': self.reference_codon_sequences,
        }


class StreamingSplitProcessor:
    """Processes a split using streaming accumulators.

    This class coordinates the accumulators to process reads for a single
    split in a streaming fashion.

    Parameters
    ==========
    split : contigops.Split
        The split object to process.

    skip_SNV : bool, False
        Skip SNV profiling.

    skip_INDEL : bool, False
        Skip INDEL profiling.

    profile_SCVs : bool, False
        Profile single codon variants.

    Notes
    =====
    This class is the main entry point for the streaming profiler.
    It replaces the multi-pass approach in contigops.Auxiliary.
    """

    def __init__(self, split, skip_SNV=False, skip_INDEL=False, profile_SCVs=False):
        self.split = split
        self.skip_SNV = skip_SNV
        self.skip_INDEL = skip_INDEL
        self.profile_SCVs = profile_SCVs

        # Initialize accumulators
        self.coverage_acc = CoverageAccumulator(split.length)

        if not skip_SNV:
            self.snv_acc = SNVAccumulator(split.length)
        else:
            self.snv_acc = None

        if not skip_INDEL:
            self.indel_acc = INDELAccumulator(split.length)
        else:
            self.indel_acc = None

        if profile_SCVs and not skip_SNV:
            self.scv_acc = SCVAccumulator(
                split.length,
                per_position_info=getattr(split, 'per_position_info', None)
            )
        else:
            self.scv_acc = None

    def process_read(self, read, split_start):
        """Process a single read and update all accumulators.

        Parameters
        ==========
        read : bamops.Read
            The read to process.

        split_start : int
            The reference position of the split start.
        """
        # Extract all variant evidence in a single pass
        evidence = read.extract_variant_evidence(
            split_start=split_start,
            skip_SNV=self.skip_SNV,
            skip_INDEL=self.skip_INDEL
        )

        # Update coverage
        self.coverage_acc.update(evidence['coverage_blocks'])

        # Update SNV accumulator
        if self.snv_acc is not None:
            self.snv_acc.update_from_vectorized(evidence['vectorized'], split_start)

        # Update INDEL accumulator
        if self.indel_acc is not None:
            self.indel_acc.update(
                evidence['insertions'],
                evidence['deletions'],
                self.split.name,
                self.split.sequence,
                getattr(self.split, 'per_position_info', None)
            )

        # Update SCV accumulator
        if self.scv_acc is not None:
            self.scv_acc.update_from_vectorized(evidence['vectorized'], split_start)

    def finalize(self):
        """Finalize all accumulators and return results.

        Returns
        =======
        dict with keys:
            'coverage': coverage accumulator results
            'snv_allele_counts': SNV allele counts array (or None)
            'indels': INDEL profiles dict (or None)
            'scv': SCV profiles per gene (or None)
        """
        results = {
            'coverage': self.coverage_acc.finalize(),
            'snv_allele_counts': self.snv_acc.finalize() if self.snv_acc else None,
            'indels': self.indel_acc.finalize() if self.indel_acc else None,
            'scv': self.scv_acc.finalize() if self.scv_acc else None,
        }
        return results
