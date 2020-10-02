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
            An existing Progress object to be updated
            With the default of None, no progress will be printed.
        """

        self.seq_names = seq_names
        self.seq_strings = seq_strings
        self.num_threads = num_threads
        self.progress = progress
        self.agglomerated_aligned_query_dict = None
        self.agglomerated_aligned_reference_dict = None


    def agglomerate(self,
                    max_mismatch_freq=0,
                    priority_function=None,
                    alignment_target_chunk_size=20000,
                    alignment_progress_interval=100000,
                    agglomeration_progress_interval=10000):
        """Agglomerate sequences by aligning all to all and then remapping alignments to seed references.

        Sets the attributes, `agglomerated_aligned_query_dict` and `agglomerated_aligned_reference_dict`:
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

        if self.progress:
            self.progress.update("Aligning sequences to themselves")

        if priority_function is None:
            priority_function = lambda aligned_reference: (-len(aligned_reference.seq_string),
                                                           -len(aligned_reference.alignments),
                                                           aligned_reference.name)

        # The `aligned_query_dict` output of `Aligner.align`
        # is named `agglomerated_aligned_query_dict` and modified during agglomeration.
        (agglomerated_aligned_query_dict,
         aligned_reference_dict) = Aligner(self.seq_names,
                                           self.seq_strings,
                                           self.seq_names,
                                           self.seq_strings,
                                           num_threads=self.num_threads,
                                           progress=self.progress).align(max_mismatch_freq=max_mismatch_freq,
                                                                         target_chunk_size=alignment_target_chunk_size,
                                                                         query_progress_interval=alignment_progress_interval)

        for agglomerated_aligned_query in agglomerated_aligned_query_dict.values():
            agglomerated_aligned_query.alignments = []

        if self.progress:
            self.progress.update("Agglomerating alignments")

        # Agglomerated clusters should preferentially be seeded
        # by the longest reference sequences with the most alignments.
        ordered_reference_names = [aligned_reference.name for aligned_reference
                                   in sorted(aligned_reference_dict.values(), key=priority_function)]
        ordered_reference_inputs = [(name, i) for i, name in enumerate(ordered_reference_names)]

        # This dict is used to track which sequences have been agglomerated.
        # Keys are sequence names; values are priority rank.
        # When a sequence is agglomerated, either as the reference seed of a cluster or a member,
        # this dict is updated with the priority of its reference seed
        # if the priority is lower (stronger) than the existing priority for the sequence in the dict.
        processed_reference_dict = {name: len(ordered_reference_names) for name in ordered_reference_names}

        agglomerated_aligned_references = []
        processed_input_count = 0
        total_input_count = len(ordered_reference_inputs)
        for agglomerated_reference_priority, name in enumerate(ordered_reference_names):

            if agglomerated_reference_priority >= processed_reference_dict[name]:
                # The reference sequence has already been processed,
                # as it mapped to another reference sequence that had been processed.
                processed_input_count += 1
                if self.progress:
                    if processed_input_count % agglomeration_progress_interval == 0:
                        self.progress.update("%d/%d sequences processed in agglomerative remapping"
                                             % (processed_input_count, total_input_count))
                continue

            processed_reference_dict[name] = agglomerated_reference_priority

            aligned_reference = aligned_reference_dict[name]
            agglomerated_aligned_reference = AlignedTarget(aligned_reference.seq_string, name=name)

            # Track the references agglomerated with this seed.
            presently_processed_reference_names = [name]

            remapping_stack = deque()
            remapping_stack.append((name, aligned_reference, 0, {}, []))

            while remapping_stack:
                remapping_item = remapping_stack.pop()
                current_reference_name = remapping_item[0]
                current_aligned_reference = remapping_item[1]
                # Record mismatches between query sequences and the agglomerated reference sequence,
                # with the coordinate system being nucleotide positions in the agglomerated reference.
                current_reference_start_in_agglomerated_reference = remapping_item[2]
                agglomerated_reference_mismatch_dict = remapping_item[3]
                current_reference_mismatches_to_agglomerated_reference = remapping_item[4]

                next_remapping_items = []
                for alignment in current_aligned_reference.alignments:
                    alignment_length = alignment.alignment_length
                    aligned_query = alignment.aligned_query
                    query_name = aligned_query.name
                    query_seq_string = aligned_query.seq_string

                    if query_name == current_reference_name:
                        # Ignore a sequence mapping to itself.
                        continue

                    if query_name in presently_processed_reference_names:
                        # The query has already been agglomerated with this seed,
                        # as it mapped to another agglomerated sequence.
                        continue

                    # In the next iteration, the current query sequence will be investigated as a reference.
                    presently_processed_reference_names.append(query_name)
                    if agglomerated_reference_priority < processed_reference_dict[query_name]:
                        processed_reference_dict[query_name] = agglomerated_reference_priority

                    # Get the mismatches between the query and the current reference,
                    # with the coordinate system being nucleotide positions in the current reference.
                    query_mismatches_to_current_reference_in_alignment_frame = []
                    current_reference_seq_string = alignment.aligned_target.seq_string
                    alignment_start_in_current_reference = alignment.target_start
                    current_reference_index = alignment_start_in_current_reference
                    for cigartuple in alignment.cigartuples:
                        if cigartuple[0] == 8:
                            for incremental_index in range(cigartuple[1]):
                                mismatch_index = current_reference_index + incremental_index
                                current_reference_nucleotide = current_reference_seq_string[mismatch_index]
                                query_mismatches_to_current_reference_in_alignment_frame.append((mismatch_index - alignment_start_in_current_reference,
                                                                                                current_reference_nucleotide))
                        current_reference_index += cigartuple[1]

                    # Position of the alignment in the coordinate system of the agglomerated reference sequence
                    alignment_start_in_agglomerated_reference = current_reference_start_in_agglomerated_reference + alignment_start_in_current_reference
                    alignment_end_in_agglomerated_reference = alignment_start_in_agglomerated_reference + alignment_length

                    # Record mismatches between the query and current reference
                    # that are also mismatches between the query and agglomerated reference.
                    # When a new mismatch with the agglomerated reference is encountered,
                    # it is added to the dict of all agglomerated reference mismatches.
                    query_mismatches_to_agglomerated_reference = []
                    query_mismatches_to_agglomerated_reference_in_alignment_frame = []
                    for alignment_index, current_reference_nucleotide in query_mismatches_to_current_reference_in_alignment_frame:
                        agglomerated_reference_index = alignment_index + alignment_start_in_agglomerated_reference
                        agglomerated_reference_nucleotide = agglomerated_reference_mismatch_dict.get(agglomerated_reference_index)
                        if agglomerated_reference_nucleotide:
                            query_nucleotide = query_seq_string[alignment_index]
                            if agglomerated_reference_nucleotide == query_nucleotide:
                                continue
                        else:
                            # This nucleotide in the agglomerated reference sequence has matched all other aligned sequences thus far.
                            agglomerated_reference_nucleotide = current_reference_nucleotide
                            agglomerated_reference_mismatch_dict[agglomerated_reference_index] = agglomerated_reference_nucleotide
                        query_mismatches_to_agglomerated_reference.append((agglomerated_reference_index,
                                                                        agglomerated_reference_nucleotide))
                        query_mismatches_to_agglomerated_reference_in_alignment_frame.append((alignment_index,
                                                                                            agglomerated_reference_nucleotide))

                    # Record mismatches between the query and agglomerated reference
                    # at positions where the query matches the current reference.
                    query_mismatch_to_current_reference_in_agglomerated_reference_frame_indices = [alignment_index + alignment_start_in_agglomerated_reference
                                                                                                for alignment_index, _
                                                                                                in query_mismatches_to_current_reference_in_alignment_frame]
                    for agglomerated_reference_index, agglomerated_reference_nucleotide in current_reference_mismatches_to_agglomerated_reference:
                        if agglomerated_reference_index in query_mismatch_to_current_reference_in_agglomerated_reference_frame_indices:
                            # The mismatch position has already been considered,
                            # as there is a mismatch between the query and current reference at this position as well.
                            continue
                        if alignment_start_in_agglomerated_reference <= agglomerated_reference_index < alignment_end_in_agglomerated_reference:
                            # Only consider mismatches within the bounds of the alignment between the query and current reference.
                            query_mismatches_to_agglomerated_reference.append((agglomerated_reference_index,
                                                                            agglomerated_reference_nucleotide))
                            query_mismatches_to_agglomerated_reference_in_alignment_frame.append((agglomerated_reference_index - alignment_start_in_agglomerated_reference,
                                                                                                agglomerated_reference_nucleotide))

                    # Change the properties of the alignment to reflect remapping to the agglomerated reference.
                    cigartuples = []
                    # Sort mismatches by position.
                    query_mismatches_to_agglomerated_reference_in_alignment_frame.sort(key=lambda query_mismatch_item: query_mismatch_item[0])
                    prev_alignment_index = -1
                    prev_agglomerated_reference_nucleotide = ''
                    for alignment_index, agglomerated_reference_nucleotide in query_mismatches_to_agglomerated_reference_in_alignment_frame:
                        if alignment_index > prev_alignment_index + 1:
                            cigartuples.append((7, alignment_index - prev_alignment_index - 1))
                        if cigartuples:
                            if cigartuples[-1][0] == 8:
                                cigartuples[-1] = (8, cigartuples[-1][1] + 1)
                            else:
                                cigartuples.append((8, 1))
                        else:
                            cigartuples.append((8, 1))
                        prev_alignment_index = alignment_index
                        prev_agglomerated_reference_nucleotide = agglomerated_reference_nucleotide
                    if alignment_length > prev_alignment_index + 1:
                        cigartuples.append((7, alignment_length - prev_alignment_index - 1))

                    agglomerated_aligned_query = agglomerated_aligned_query_dict[query_name]
                    agglomerated_alignment = Alignment(alignment.query_start,
                                                    alignment_start_in_agglomerated_reference,
                                                    cigartuples,
                                                    aligned_query=agglomerated_aligned_query,
                                                    aligned_target=agglomerated_aligned_reference)
                    agglomerated_aligned_query.alignments.append(agglomerated_alignment)
                    agglomerated_aligned_reference.alignments.append(agglomerated_alignment)

                    next_remapping_items.append((query_name,
                                                aligned_reference_dict[query_name],
                                                alignment_start_in_agglomerated_reference,
                                                dict(agglomerated_reference_mismatch_dict.items()),
                                                [mismatch_tuple for mismatch_tuple in query_mismatches_to_agglomerated_reference]))

                for next_remapping_item in next_remapping_items[::-1]:
                    remapping_stack.append(next_remapping_item)

            agglomerated_aligned_references.append(agglomerated_aligned_reference)
            processed_input_count += 1
            if self.progress:
                if processed_input_count % agglomeration_progress_interval == 0:
                    self.progress.update("%d/%d sequences processed in agglomerative remapping"
                                         % (processed_input_count, total_input_count))
        if self.progress:
            if processed_input_count % agglomeration_progress_interval != 0:
                self.progress.update("%d/%d sequences processed in agglomerative remapping"
                                     % (total_input_count, total_input_count))

        agglomerated_aligned_reference_dict = {ref.name: ref for ref in agglomerated_aligned_references}

        self.agglomerated_aligned_query_dict = agglomerated_aligned_query_dict
        self.agglomerated_aligned_reference_dict = agglomerated_aligned_reference_dict
