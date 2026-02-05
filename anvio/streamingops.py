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

    def update_from_vectorized(self, vectorized, split_start, position_limit=None):
        """Update allele counts from a vectorized read array.

        Parameters
        ==========
        vectorized : numpy.ndarray
            The vectorized read array from Read.vectorize().
            Columns: [ref_pos, query_seq, mapping_type, ref_seq]

        split_start : int
            The reference position of the split start.

        position_limit : int, optional
            If provided, only count positions in range [0, position_limit).
            Used for chunk-based processing where the accumulator is smaller
            than the full split.
        """
        limit = position_limit if position_limit is not None else self.split_length

        # Filter to only mapped positions (mapping_type == 0)
        mapped_mask = vectorized[:, 2] == 0
        mapped = vectorized[mapped_mask]

        if len(mapped) == 0:
            return

        # Compute positions relative to split start
        positions = mapped[:, 0] - split_start
        bases = mapped[:, 1]

        # Filter to valid position range [0, limit)
        valid_mask = (positions >= 0) & (positions < limit)
        positions = positions[valid_mask]
        bases = bases[valid_mask]

        if len(positions) == 0:
            return

        # Map nucleotide ord values to array indices (vectorized)
        # A=65, C=67, G=71, T=84, N=78 (uppercase)
        # a=97, c=99, g=103, t=116, n=110 (lowercase)
        base_indices = np.full(len(bases), -1, dtype=np.int32)
        base_indices[(bases == 65) | (bases == 97)] = 0   # A/a
        base_indices[(bases == 67) | (bases == 99)] = 1   # C/c
        base_indices[(bases == 71) | (bases == 103)] = 2  # G/g
        base_indices[(bases == 84) | (bases == 116)] = 3  # T/t
        base_indices[(bases == 78) | (bases == 110)] = 4  # N/n

        # Filter out unknown bases
        valid_bases_mask = base_indices >= 0
        positions = positions[valid_bases_mask]
        base_indices = base_indices[valid_bases_mask]

        # Scatter-add: increment allele_counts at (base_index, position) pairs
        np.add.at(self.allele_counts, (base_indices, positions), 1)

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

            # Convert ordinal array to string efficiently using bytes conversion
            ins_seq = bytes(ins_seq_ord.astype(np.uint8)).decode('ascii')
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


class SharedSequenceStore:
    """Stores contig sequences in shared memory for multi-process access.

    This class packs all contig sequences into a single shared memory buffer,
    allowing multiple worker processes to access sequences without duplicating
    the data in each process's memory space.

    Parameters
    ==========
    contig_sequences : dict
        Dictionary mapping contig names to {'sequence': str} dicts.
        This is the format returned by init_contig_sequences().

    Attributes
    ==========
    shm : multiprocessing.shared_memory.SharedMemory
        The shared memory object containing all sequences.
    index : dict
        Maps contig_name -> (offset, length) in the shared buffer.
    shm_name : str
        Name of the shared memory block (needed for workers to attach).

    Notes
    =====
    - The owner (main process) must call close() and unlink() when done
    - Workers should call close() but NOT unlink()
    - Sequences are stored as ASCII bytes in the shared buffer

    Example
    =======
    # In main process:
    store = SharedSequenceStore(contig_sequences)
    # Pass store.shm_name and store.index to workers
    ...
    store.close()
    store.unlink()

    # In worker process:
    store = SharedSequenceStore.from_existing(shm_name, index)
    seq = store.get_sequence('contig_name', start, end)
    ...
    store.close()  # Don't unlink!
    """

    def __init__(self, contig_sequences):
        from multiprocessing import shared_memory

        # Calculate total size needed
        total_size = sum(len(v['sequence']) for v in contig_sequences.values())

        # Create shared memory
        self.shm = shared_memory.SharedMemory(create=True, size=total_size)
        self.shm_name = self.shm.name
        self._is_owner = True

        # Build index and copy sequences to shared memory
        self.index = {}
        offset = 0
        buffer = self.shm.buf

        for contig_name, data in contig_sequences.items():
            seq = data['sequence']
            seq_bytes = seq.encode('ascii')
            length = len(seq_bytes)

            # Copy to shared memory
            buffer[offset:offset + length] = seq_bytes

            # Record in index
            self.index[contig_name] = (offset, length)
            offset += length

    @classmethod
    def from_existing(cls, shm_name, index):
        """Attach to existing shared memory (for worker processes).

        Parameters
        ==========
        shm_name : str
            Name of the shared memory block.
        index : dict
            The index mapping contig_name -> (offset, length).

        Returns
        =======
        SharedSequenceStore
            A store instance attached to existing shared memory.
        """
        from multiprocessing import shared_memory

        instance = cls.__new__(cls)
        instance.shm = shared_memory.SharedMemory(name=shm_name)
        instance.shm_name = shm_name
        instance.index = index
        instance._is_owner = False
        return instance

    def get_sequence(self, contig_name, start=None, end=None):
        """Get a contig sequence or subsequence.

        Parameters
        ==========
        contig_name : str
            Name of the contig.
        start : int, optional
            Start position (0-indexed). If None, starts from beginning.
        end : int, optional
            End position (exclusive). If None, goes to end.

        Returns
        =======
        str
            The sequence or subsequence.
        """
        offset, length = self.index[contig_name]

        if start is None:
            start = 0
        if end is None:
            end = length

        # Adjust for the contig's position in the buffer
        buf_start = offset + start
        buf_end = offset + end

        return bytes(self.shm.buf[buf_start:buf_end]).decode('ascii')

    def close(self):
        """Close access to shared memory. Call from all processes."""
        self.shm.close()

    def unlink(self):
        """Remove shared memory. Only call from the owner (main process)."""
        if self._is_owner:
            self.shm.unlink()

    def as_dict_proxy(self):
        """Return a dict-like proxy for compatibility with contig_sequences interface.

        Returns a SharedSequenceDict that provides the same interface as the
        original contig_sequences dict: d[contig_name]['sequence']
        """
        return SharedSequenceDict(self)


class SharedSequenceDict:
    """Dict-like proxy for SharedSequenceStore.

    Provides the same interface as the contig_sequences dict, allowing code
    that accesses contig_sequences[contig_name]['sequence'] to work
    transparently with shared memory.

    This class implements __len__, __getitem__, __contains__, and keys()
    to be compatible with typical usage patterns.
    """

    def __init__(self, store):
        """
        Parameters
        ==========
        store : SharedSequenceStore
            The underlying shared memory store.
        """
        self._store = store

    def __len__(self):
        """Return number of contigs."""
        return len(self._store.index)

    def __contains__(self, contig_name):
        """Check if contig exists."""
        return contig_name in self._store.index

    def __getitem__(self, contig_name):
        """Return a dict-like object with 'sequence' key for the contig."""
        if contig_name not in self._store.index:
            raise KeyError(contig_name)
        return _SharedSequenceEntry(self._store, contig_name)

    def keys(self):
        """Return contig names."""
        return self._store.index.keys()

    def __iter__(self):
        """Iterate over contig names."""
        return iter(self._store.index.keys())


class _SharedSequenceEntry:
    """Dict-like object representing a single contig's entry.

    Provides access to the 'sequence' key, fetching from shared memory.
    """

    def __init__(self, store, contig_name):
        self._store = store
        self._contig_name = contig_name

    def __getitem__(self, key):
        if key == 'sequence':
            return self._store.get_sequence(self._contig_name)
        raise KeyError(key)

    def __contains__(self, key):
        return key == 'sequence'

    def get(self, key, default=None):
        if key == 'sequence':
            return self._store.get_sequence(self._contig_name)
        return default


class SharedNtPositionsStore:
    """Stores nt_positions_info arrays in shared memory for multi-process access.

    This class packs all per-contig numpy arrays (uint8) into a single shared
    memory buffer, allowing multiple worker processes to access the data without
    duplicating it in each process's memory space.

    Parameters
    ==========
    nt_positions_info : dict
        Dictionary mapping contig names to numpy arrays (uint8).
        This is the format returned by the nt_positions_info property.

    Attributes
    ==========
    shm : multiprocessing.shared_memory.SharedMemory
        The shared memory object containing all arrays.
    index : dict
        Maps contig_name -> (offset, length) in the shared buffer.
    shm_name : str
        Name of the shared memory block (needed for workers to attach).

    Notes
    =====
    - The owner (main process) must call close() and unlink() when done
    - Workers should call close() but NOT unlink()
    - Arrays are stored as raw uint8 bytes in the shared buffer
    """

    def __init__(self, nt_positions_info):
        from multiprocessing import shared_memory

        # Calculate total size needed
        total_size = sum(arr.nbytes for arr in nt_positions_info.values())

        if total_size == 0:
            # Handle empty case
            self.shm = None
            self.shm_name = None
            self.index = {}
            self._is_owner = True
            return

        # Create shared memory
        self.shm = shared_memory.SharedMemory(create=True, size=total_size)
        self.shm_name = self.shm.name
        self._is_owner = True

        # Build index and copy arrays to shared memory
        self.index = {}
        offset = 0
        buffer = self.shm.buf

        for contig_name, arr in nt_positions_info.items():
            length = arr.nbytes

            # Copy to shared memory
            buffer[offset:offset + length] = arr.tobytes()

            # Record in index
            self.index[contig_name] = (offset, length)
            offset += length

    @classmethod
    def from_existing(cls, shm_name, index):
        """Attach to existing shared memory (for worker processes).

        Parameters
        ==========
        shm_name : str
            Name of the shared memory block.
        index : dict
            The index mapping contig_name -> (offset, length).

        Returns
        =======
        SharedNtPositionsStore
            A store instance attached to existing shared memory.
        """
        from multiprocessing import shared_memory

        instance = cls.__new__(cls)
        if shm_name is None:
            instance.shm = None
            instance.shm_name = None
            instance.index = {}
        else:
            instance.shm = shared_memory.SharedMemory(name=shm_name)
            instance.shm_name = shm_name
            instance.index = index
        instance._is_owner = False
        return instance

    def get_array(self, contig_name):
        """Get a contig's nt_positions array.

        Parameters
        ==========
        contig_name : str
            Name of the contig.

        Returns
        =======
        numpy.ndarray
            The uint8 array for this contig.
        """
        if contig_name not in self.index:
            raise KeyError(contig_name)

        offset, length = self.index[contig_name]

        # Create numpy array view of shared memory
        return np.frombuffer(self.shm.buf, dtype=np.uint8, count=length, offset=offset)

    def close(self):
        """Close access to shared memory. Call from all processes."""
        if self.shm is not None:
            self.shm.close()

    def unlink(self):
        """Remove shared memory. Only call from the owner (main process)."""
        if self._is_owner and self.shm is not None:
            self.shm.unlink()

    def as_dict_proxy(self):
        """Return a dict-like proxy for compatibility with nt_positions_info interface.

        Returns a SharedNtPositionsDict that provides the same interface as the
        original nt_positions_info dict: d[contig_name] returns numpy array
        """
        return SharedNtPositionsDict(self)


class SharedNtPositionsDict:
    """Dict-like proxy for SharedNtPositionsStore.

    Provides the same interface as the nt_positions_info dict, allowing code
    that accesses nt_positions_info[contig_name] to work transparently with
    shared memory.
    """

    def __init__(self, store):
        self._store = store

    def __len__(self):
        """Return number of contigs."""
        return len(self._store.index)

    def __bool__(self):
        """Return True if there are any contigs."""
        return len(self._store.index) > 0

    def __contains__(self, contig_name):
        """Check if contig exists."""
        return contig_name in self._store.index

    def __getitem__(self, contig_name):
        """Return the numpy array for the contig."""
        return self._store.get_array(contig_name)

    def keys(self):
        """Return contig names."""
        return self._store.index.keys()

    def __iter__(self):
        """Iterate over contig names."""
        return iter(self._store.index.keys())

    def get(self, contig_name, default=None):
        """Get array for contig, or default if not found."""
        if contig_name in self._store.index:
            return self._store.get_array(contig_name)
        return default
