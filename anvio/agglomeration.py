#!/usr/bin/env python
# -*- coding: utf-8

import anvio
import anvio.terminal as terminal

from anvio.sequence import Aligner, AlignedTarget, Alignment

from collections import deque


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print


class Agglomerator:
    def __init__(self, seq_names, seq_strings, num_threads=1, progress=None):
        """This class agglomerates sequences into clusters.

        Parameters
        ==========
        seq_names : list
            List of name strings corresponding to the sequences in `seq_strings`

        seq_strings : list
            List of sequence strings corresponding to the names in `seq_names`

        num_threads : int, 1
            Number of threads available for alignment

        progress : terminal.Progress object, None
        """
        self.seq_names = seq_names
        self.seq_strings = seq_strings
        self.num_threads = num_threads
        if not progress:
            progress = terminal.Progress()
            progress.new("Agglomerating")
        self.progress = progress
        self.agglom_aligned_query_dict = None
        self.agglom_aligned_ref_dict = None


    def agglomerate(self,
                    max_mismatch_freq=0,
                    priority_function=None,
                    alignment_target_chunk_size=20000,
                    alignment_progress_interval=100000,
                    agglom_progress_interval=10000):
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

        alignment_target_chunk_size : int, 20000
            The chunk size for k-mer dict formation from alignment target sequences:
            all queries are aligned to the first chunk of targets, then the second chunk, etc.
            The default balances speed and memory consumption:
            the larger the chunks, the faster the mapping and the more memory consumed.

        alignment_progress_interval : int, 100000
            The number of queries aligned to a chunk of targets between alignment progress statements

        agglomeration_progress_interval : int, 10000
            The number of alignment references remapped between progress statements
        """
        self.progress.update("Aligning sequences to themselves")

        if priority_function is None:
            priority_function = lambda aligned_ref: (-len(aligned_ref.seq_string),
                                                     -len(aligned_ref.alignments),
                                                     aligned_ref.name)

        # The `aligned_query_dict` output of `Aligner.align`
        # is named `agglomerated_aligned_query_dict` and modified during agglomeration.
        (agglom_aligned_query_dict,
         aligned_ref_dict) = Aligner(self.seq_names,
                                     self.seq_strings,
                                     self.seq_names,
                                     self.seq_strings,
                                     num_threads=self.num_threads,
                                     progress=self.progress).align(max_mismatch_freq=max_mismatch_freq,
                                                                   target_chunk_size=alignment_target_chunk_size,
                                                                   query_progress_interval=alignment_progress_interval)

        for agglom_aligned_query in agglom_aligned_query_dict.values():
            agglom_aligned_query.alignments = []

        self.progress.update("Agglomerating alignments")

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
        processed_input_count = 0
        total_input_count = len(ordered_ref_inputs)
        for agglom_ref_priority, name in enumerate(ordered_ref_names):

            if agglom_ref_priority >= processed_ref_dict[name]:
                # The reference sequence has already been processed,
                # as it mapped to another reference sequence that had been processed.
                processed_input_count += 1
                if processed_input_count % agglom_progress_interval == 0:
                    self.progress.update("%s/%s seqs processed in agglomerative remapping"
                                         % (pp(processed_input_count), pp(total_input_count)))
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

                    if query_name == current_ref_name:
                        # Ignore a sequence mapping to itself.
                        continue

                    if query_name in presently_processed_ref_names:
                        # The query has already been agglomerated with this seed,
                        # as it mapped to another agglomerated sequence.
                        continue

                    # In the next iteration, the current query sequence will be investigated as a reference.
                    presently_processed_ref_names.append(query_name)
                    if agglom_ref_priority < processed_ref_dict[query_name]:
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
            processed_input_count += 1
            if processed_input_count % agglom_progress_interval == 0:
                self.progress.update("%s/%s seqs processed in agglomerative remapping"
                                     % (pp(processed_input_count), pp(total_input_count)))
        if processed_input_count % agglom_progress_interval != 0:
            self.progress.update("%s/%s seqs processed in agglomerative remapping"
                                 % (pp(total_input_count), pp(total_input_count)))

        agglom_aligned_ref_dict = {ref.name: ref for ref in agglom_aligned_refs}

        self.agglom_aligned_query_dict = agglom_aligned_query_dict
        self.agglom_aligned_ref_dict = agglom_aligned_ref_dict
