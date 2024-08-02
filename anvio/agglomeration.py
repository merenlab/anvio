#!/usr/bin/env python
# -*- coding: utf-8

import gc
import os
import math
import argparse

from collections import deque

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.vmatch import Vmatch
from anvio.sequence import Alignment, AlignedQuery, AlignedTarget


__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print


class Agglomerator:
    def __init__(self, seq_names, seq_strings, num_threads=1):
        """This class agglomerates sequences into clusters.

        Parameters
        ==========
        seq_names : list
            List of name strings corresponding to the sequences in `seq_strings`

        seq_strings : list
            List of sequence strings corresponding to the names in `seq_names`

        num_threads : int, 1
            Number of threads available for alignment
        """
        self.seq_names = seq_names
        self.seq_strings = seq_strings
        self.num_threads = num_threads
        self.agglom_aligned_query_dict = None
        self.agglom_aligned_ref_dict = None


    def agglomerate(self, max_mismatch_freq=0, priority_function=None):
        """Agglomerate sequences by aligning all to all and then remapping alignments to seed references.

        Sets the attributes, `agglom_aligned_query_dict` and `agglom_aligned_ref_dict`:
        the former maps each input sequence name to a `sequence.AlignedQuery` object
        and the latter maps input sequence names corresponding to cluster seed sequences to `sequence.AlignedReference` objects.
        AlignedQuery and AlignedReference objects each contain a list of Alignment objects which relate queries to references.

        Parameters
        ==========
        max_mismatch_freq : float, 0
            The maximum mismatch frequency, lying in the interval [0, 1), allowed in alignments.
            The higher the value, the larger the agglomerated clusters.

        priority_function : function reference, None
            The priority function should map the input sequences to a numeric rank of how clusters will be seeded.
            By default, the priority function ranks in descending order of sequence length,
            then descending order of number of alignments, then in ascending order of sequence name.
            (The longest sequence with the most alignments is the first to seed a cluster.)
        """
        progress = terminal.Progress()
        progress.new("Agglomerating")
        progress.update("Writing FASTA file of sequences")
        temp_dir_path = filesnpaths.get_temp_directory_path()
        fasta_path = os.path.join(temp_dir_path, 'seqs.fa')
        seq_dict = {}
        with open(fasta_path, 'w') as fasta_file:
            for name, seq_string in zip(self.seq_names, self.seq_strings):
                fasta_file.write(f">{name}\n{seq_string}\n")
                seq_dict[name] = seq_string
        progress.end()

        align_df = Vmatch(argparse.Namespace(match_mode='query_substring_with_mismatches',
                                             fasta_db_file=fasta_path,
                                             fasta_query_file=fasta_path,
                                             num_threads=self.num_threads,
                                             max_hamming_dist=math.ceil(max(map(len, self.seq_strings)) * max_mismatch_freq),
                                             min_ident=int(100 - 100 * max_mismatch_freq),
                                             align_output_length=10,
                                             temp_dir=temp_dir_path)).search_queries()

        pid = "Parsing alignments"
        progress.new(pid)
        # The dictionary of aligned queries is named `agglom_aligned_query_dict` to indicate that
        # its AlignedQuery objects are modified during agglomeration.
        agglom_aligned_query_dict = {}
        aligned_ref_dict = {}
        num_processed_aligns = -1
        parsing_progress_interval = 10000
        total_align_count = len(align_df)
        pp_total_align_count = pp(total_align_count)
        for query_name, query_align_df in align_df.groupby('query_name'):
            query_seq_string = seq_dict[query_name]
            query_length = len(query_seq_string)
            aligned_query = AlignedQuery(query_seq_string, query_name)
            agglom_aligned_query_dict[query_name] = aligned_query
            for target_name, query_start_in_target, mismatch_positions in zip(query_align_df['target_name'],
                                                                              query_align_df['query_start_in_target'],
                                                                              query_align_df['mismatch_positions']):
                num_processed_aligns += 1
                if num_processed_aligns % parsing_progress_interval == 0:
                    pp_progress_interval_end = pp(total_align_count if num_processed_aligns + parsing_progress_interval > total_align_count else num_processed_aligns + parsing_progress_interval)
                    progress.update_pid(pid)
                    progress.update(f"{pp(num_processed_aligns + 1)}-{pp_progress_interval_end}/{pp_total_align_count}")

                try:
                    aligned_target = aligned_ref_dict[target_name]
                except KeyError:
                    aligned_target = AlignedTarget(seq_dict[target_name], target_name)
                    aligned_ref_dict[target_name] = aligned_target

                # Convert the positions of mismatches into a cigar tuple for the alignment. The
                # search method ensured that each alignment contains at least one mismatch.
                cigartuples = []
                prev_mismatch_pos = -2
                for mismatch_num, mismatch_pos in enumerate(map(int, mismatch_positions.split(','))):
                    if prev_mismatch_pos == -2:
                        # This is the first mismatch in the alignment.
                        if mismatch_pos > 0:
                            # There is not a mismatch at the first position of the query.
                            cigartuples.append((7, mismatch_pos))
                        cigartuples.append((8, 1))
                    elif mismatch_pos == prev_mismatch_pos + 1:
                        # This mismatch follows another mismatch.
                        cigartuples[-1] = (8, cigartuples[-1][1] + 1)
                    else:
                        cigartuples.append((7, mismatch_pos - prev_mismatch_pos - 1))
                        cigartuples.append((8, 1))
                    prev_mismatch_pos = mismatch_pos
                if query_length - prev_mismatch_pos > 1:
                    cigartuples.append((7, query_length - prev_mismatch_pos - 1))

                alignment = Alignment(0, query_start_in_target, cigartuples, aligned_query, aligned_target)
                # The Alignment doesn't need to be added to the AlignedQuery object, as these are
                # changed later when queries are remapped.
                aligned_target.alignments.append(alignment)
        del seq_dict
        gc.collect()
        progress.end()

        pid = "Agglomerating aligned reference seqs"
        progress.new(pid)
        progress.update("...")

        if priority_function is None:
            priority_function = lambda aligned_ref: (-len(aligned_ref.seq_string),
                                                     -len(aligned_ref.alignments),
                                                     aligned_ref.name)

        for agglom_aligned_query in agglom_aligned_query_dict.values():
            agglom_aligned_query.alignments = []

        # Agglomerated clusters should preferentially be seeded
        # by the longest reference sequences with the most alignments.
        ordered_ref_names = [aligned_ref.name for aligned_ref
                             in sorted(aligned_ref_dict.values(), key=priority_function)]
        ordered_ref_inputs = [(name, i) for i, name in enumerate(ordered_ref_names)]

        # This dict is used to track which sequences have been agglomerated.
        # Keys are sequence names; values are priority rank.
        # When a sequence is agglomerated, either as the reference seed of a cluster or a member,
        # this dict is updated with the priority of its reference seed
        # if the priority is lower (stronger) than the existing priority for the sequence in the dict.
        processed_ref_dict = {name: len(ordered_ref_names) for name in ordered_ref_names}

        agglom_aligned_refs = []
        agglom_progress_interval = 1000
        total_ref_count = len(ordered_ref_inputs)
        pp_total_ref_count = pp(total_ref_count)
        num_processed_refs = -1
        for agglom_ref_priority, name in enumerate(ordered_ref_names):
            num_processed_refs += 1
            if num_processed_refs % agglom_progress_interval == 0:
                pp_progress_interval_end = pp(total_ref_count if num_processed_refs + agglom_progress_interval > total_ref_count else num_processed_refs + agglom_progress_interval)
                progress.update_pid(pid)
                progress.update(f"{pp(num_processed_refs + 1)}-{pp_progress_interval_end}/{pp_total_ref_count}")

            if agglom_ref_priority >= processed_ref_dict[name]:
                # The reference sequence has already been processed,
                # as it mapped to another reference sequence that had been processed.
                continue

            processed_ref_dict[name] = agglom_ref_priority

            aligned_ref = aligned_ref_dict[name]
            agglom_aligned_ref = AlignedTarget(aligned_ref.seq_string, name=name)

            # Track the references agglomerated with this seed.
            presently_processed_ref_names = [name]

            remapping_stack = deque()
            remapping_stack.append((name, aligned_ref, 0, {}, []))

            while remapping_stack:
                remapping_item = remapping_stack.pop()
                current_ref_name = remapping_item[0]
                current_aligned_ref = remapping_item[1]
                # Record mismatches between query sequences and the agglomerated reference sequence,
                # with the coordinate system being nucleotide positions in the agglomerated reference.
                current_ref_start_in_agglom_ref = remapping_item[2]
                agglom_ref_mismatch_dict = remapping_item[3]
                current_ref_mismatches_to_agglom_ref = remapping_item[4]

                next_remapping_items = []
                for alignment in current_aligned_ref.alignments:
                    alignment_length = alignment.alignment_length
                    aligned_query = alignment.aligned_query
                    query_name = aligned_query.name
                    query_seq_string = aligned_query.seq_string

                    if query_name in presently_processed_ref_names:
                        # The query has already been agglomerated with this seed,
                        # as it mapped to another agglomerated sequence.
                        continue

                    try:
                        prev_ref_priority = processed_ref_dict[query_name]
                    except KeyError:
                        # No sequences aligned to the query sequence.
                        continue

                    # In the next iteration, the current query sequence will be investigated as a reference.
                    presently_processed_ref_names.append(query_name)
                    if agglom_ref_priority < prev_ref_priority:
                        processed_ref_dict[query_name] = agglom_ref_priority

                    # Get the mismatches between the query and the current reference,
                    # with the coordinate system being nucleotide positions in the current reference.
                    query_mismatches_to_current_ref_in_alignment_frame = []
                    current_ref_seq_string = alignment.aligned_target.seq_string
                    alignment_start_in_current_ref = alignment.target_start
                    current_ref_pos = alignment_start_in_current_ref
                    for cigartuple in alignment.cigartuples:
                        if cigartuple[0] == 8:
                            for incremental_pos in range(cigartuple[1]):
                                mismatch_pos = current_ref_pos + incremental_pos
                                current_ref_nt = current_ref_seq_string[mismatch_pos]
                                query_mismatches_to_current_ref_in_alignment_frame.append(
                                    (mismatch_pos - alignment_start_in_current_ref, current_ref_nt)
                                )
                        current_ref_pos += cigartuple[1]

                    # Position of the alignment in the coordinate system of the agglomerated reference sequence
                    alignment_start_in_agglom_ref = current_ref_start_in_agglom_ref + alignment_start_in_current_ref
                    alignment_end_in_agglom_ref = alignment_start_in_agglom_ref + alignment_length

                    # Record mismatches between the query and current reference
                    # that are also mismatches between the query and agglomerated reference.
                    # When a new mismatch with the agglomerated reference is encountered,
                    # it is added to the dict of all agglomerated reference mismatches.
                    query_mismatches_to_agglom_ref = []
                    query_mismatches_to_agglom_ref_in_alignment_frame = []
                    for alignment_pos, current_ref_nt in query_mismatches_to_current_ref_in_alignment_frame:
                        agglom_ref_pos = alignment_pos + alignment_start_in_agglom_ref
                        agglom_ref_nt = agglom_ref_mismatch_dict.get(agglom_ref_pos)
                        if agglom_ref_nt:
                            query_nt = query_seq_string[alignment_pos]
                            if agglom_ref_nt == query_nt:
                                continue
                        else:
                            # This nucleotide in the agglomerated reference sequence has matched all other aligned sequences thus far.
                            agglom_ref_nt = current_ref_nt
                            agglom_ref_mismatch_dict[agglom_ref_pos] = agglom_ref_nt
                        query_mismatches_to_agglom_ref.append((agglom_ref_pos, agglom_ref_nt))
                        query_mismatches_to_agglom_ref_in_alignment_frame.append((alignment_pos, agglom_ref_nt))

                    # Record mismatches between the query and agglomerated reference
                    # at positions where the query matches the current reference.
                    query_mismatch_to_current_ref_in_agglom_ref_frame_positions = [
                        alignment_pos + alignment_start_in_agglom_ref
                        for alignment_pos, _ in query_mismatches_to_current_ref_in_alignment_frame
                    ]
                    for agglom_ref_pos, agglom_ref_nt in current_ref_mismatches_to_agglom_ref:
                        if agglom_ref_pos in query_mismatch_to_current_ref_in_agglom_ref_frame_positions:
                            # The mismatch position has already been considered,
                            # as there is a mismatch between the query and current reference at this position as well.
                            continue
                        if alignment_start_in_agglom_ref <= agglom_ref_pos < alignment_end_in_agglom_ref:
                            # Only consider mismatches within the bounds of the alignment between the query and current reference.
                            query_mismatches_to_agglom_ref.append((agglom_ref_pos, agglom_ref_nt))
                            query_mismatches_to_agglom_ref_in_alignment_frame.append(
                                (agglom_ref_pos - alignment_start_in_agglom_ref, agglom_ref_nt)
                            )

                    # Change the properties of the alignment to reflect remapping to the agglomerated reference.
                    cigartuples = []
                    # Sort mismatches by position.
                    query_mismatches_to_agglom_ref_in_alignment_frame.sort(key=lambda query_mismatch_item: query_mismatch_item[0])
                    prev_alignment_pos = -1
                    prev_agglom_ref_nt = ''
                    for alignment_pos, agglom_ref_nt in query_mismatches_to_agglom_ref_in_alignment_frame:
                        if alignment_pos > prev_alignment_pos + 1:
                            cigartuples.append((7, alignment_pos - prev_alignment_pos - 1))
                        if cigartuples:
                            if cigartuples[-1][0] == 8:
                                cigartuples[-1] = (8, cigartuples[-1][1] + 1)
                            else:
                                cigartuples.append((8, 1))
                        else:
                            cigartuples.append((8, 1))
                        prev_alignment_pos = alignment_pos
                        prev_agglom_ref_nt = agglom_ref_nt
                    if alignment_length > prev_alignment_pos + 1:
                        cigartuples.append((7, alignment_length - prev_alignment_pos - 1))

                    agglom_aligned_query = agglom_aligned_query_dict[query_name]
                    agglom_alignment = Alignment(alignment.query_start,
                                                 alignment_start_in_agglom_ref,
                                                 cigartuples,
                                                 aligned_query=agglom_aligned_query,
                                                 aligned_target=agglom_aligned_ref)
                    agglom_aligned_query.alignments.append(agglom_alignment)
                    agglom_aligned_ref.alignments.append(agglom_alignment)

                    next_remapping_items.append(
                        (query_name,
                         aligned_ref_dict[query_name],
                         alignment_start_in_agglom_ref,
                         dict(agglom_ref_mismatch_dict.items()),
                         [mismatch_tuple for mismatch_tuple in query_mismatches_to_agglom_ref])
                    )

                for next_remapping_item in next_remapping_items[::-1]:
                    remapping_stack.append(next_remapping_item)

            agglom_aligned_refs.append(agglom_aligned_ref)
        agglom_aligned_ref_dict = {ref.name: ref for ref in agglom_aligned_refs}

        self.agglom_aligned_query_dict = agglom_aligned_query_dict
        self.agglom_aligned_ref_dict = agglom_aligned_ref_dict
        progress.end()
