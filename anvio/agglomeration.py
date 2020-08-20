#!/usr/bin/env python
# -*- coding: utf-8

import anvio.terminal as terminal

from anvio.sequence import Aligner, AlignedTarget, Alignment

import functools
import multiprocessing
import sys

from copy import deepcopy
from collections import deque


# Increase the recursion limit from the default of 1,000.
# Multiprocessing can raise an exception when returning results when pickle hits the recursion limit.
# The multiprocessing chunk size is set to the relatively small value (for tRNA-seq datasets) of 1,000,
# limiting the number of items returned, but the size of each item can still be large, so play it safe.
sys.setrecursionlimit(10000)
chunksize=1000


class Agglomerator:
    def __init__(self, seq_names, seq_strings, max_mismatch_freq=0, num_threads=1, progress=None):
        self.seq_names = seq_names
        self.seq_strings = seq_strings
        self.max_mismatch_freq = max_mismatch_freq

        self.num_threads = num_threads
        self.progress = progress


    def agglomerate(self):
        if self.progress:
            self.progress.update("Aligning sequences to themselves")
        # The `aligned_query_dict` output of `Aligner.align`
        # is named `agglomerated_aligned_query_dict` and modified during agglomeration.
        (agglomerated_aligned_query_dict,
         aligned_reference_dict) = Aligner(self.seq_names,
                                           self.seq_strings,
                                           self.seq_names,
                                           self.seq_strings,
                                           num_threads=self.num_threads,
                                           progress=self.progress).align(max_mismatch_freq=self.max_mismatch_freq)
        for agglomerated_aligned_query in agglomerated_aligned_query_dict.values():
            agglomerated_aligned_query.alignments = []

        if self.progress:
            self.progress.update("Agglomerating alignments")
        # Agglomerated clusters should preferentially be seeded
        # by the longest reference sequences with the most alignments.
        ordered_reference_names = [aligned_reference.name for aligned_reference
                                   in sorted(aligned_reference_dict.values(),
                                             key=lambda aligned_reference: (-len(aligned_reference.seq_string),
                                                                            -len(aligned_reference.alignments),
                                                                            aligned_reference.name))]
        ordered_reference_inputs = [(ordered_reference_name, i)
                                    for i, ordered_reference_name in enumerate(ordered_reference_names)]

        # This dict is used to track which sequences have been agglomerated.
        # Keys are sequence names; values are agglomerated reference index from `ordered_reference_names`.
        # When a sequence is agglomerated, either as the seed of a cluster or a member,
        # this dict is updated with the index of its agglomerated reference sequence
        # IF the index is lower than the existing index for the sequence in the dict.
        # After multiprocessing, this dict indicates which clusters are spurious
        # due to their seeds also being queries aligned to higher-priority seeds (lower index in `ordered_reference_names`).
        # Spurious clusters that need to be removed post facto are an inevitability of parallelization
        # due to lower-priority references occasionally being processed before higher-priority references in different processes.
        manager = multiprocessing.Manager()
        processed_reference_dict = manager.dict([(ordered_reference_name, len(ordered_reference_names))
                                                 for ordered_reference_name in ordered_reference_names])
        lock = multiprocessing.Lock()

        # Set static parameters in the multiprocessing target function.
        target = functools.partial(remap_queries,
                                   processed_reference_dict=processed_reference_dict,
                                   aligned_reference_dict=aligned_reference_dict,
                                   agglomerated_aligned_query_dict=agglomerated_aligned_query_dict)

        # The lock used with the shared dict, `processed_reference_dict`,
        # cannot be serialized and so cannot be passed to spawned processes.
        # This is circumvented by including the lock as a global variable in the forked parent process.
        pool = multiprocessing.Pool(self.num_threads, initializer=initialize, initargs=(lock, ))

        agglomerated_references = []
        processed_input_count = 0
        total_input_count = len(ordered_reference_inputs)

        for agglomerated_reference in pool.imap_unordered(target,
                                                          ordered_reference_inputs,
                                                          chunksize=chunksize):
            if agglomerated_reference:
                agglomerated_references.append(agglomerated_reference)
            processed_input_count += 1

            if processed_input_count % 10000 == 0:
                if self.progress:
                    self.progress.update("%d/%d sequences processed in agglomerative remapping"
                                        % (processed_input_count, total_input_count))
        pool.close()
        pool.join() # allow processes to terminate
        if self.progress:
            self.progress.update("%d/%d sequences processed in agglomerative remapping"
                                % (total_input_count, total_input_count))

        removal_indices = []
        candidate_agglomerated_reference_names = [agglomerated_reference.name
                                                  for agglomerated_reference in agglomerated_references]
        for reference_name, lowest_reference_index in processed_reference_dict.items():
            if reference_name in candidate_agglomerated_reference_names:
                if lowest_reference_index < ordered_reference_names.index(reference_name):
                    removal_indices.append(candidate_agglomerated_reference_names.index(reference_name))
        removal_indices.sort(reverse=True)
        for removal_index in removal_indices:
            agglomerated_references.pop(removal_index)

        agglomerated_aligned_reference_dict = {agglomerated_reference.name: agglomerated_reference
                                               for agglomerated_reference in agglomerated_references}

        return agglomerated_aligned_query_dict, agglomerated_aligned_reference_dict


def initialize(_lock):
    global lock
    lock = _lock


class RemappingItem:
    def __init__(self,
                 reference_name,
                 aligned_reference,
                 reference_start_in_agglomerated_reference,
                 agglomerated_reference_mismatch_dict,
                 reference_mismatches_to_agglomerated_reference):
        self.reference_name = reference_name
        self.aligned_reference = aligned_reference
        self.reference_start_in_agglomerated_reference = reference_start_in_agglomerated_reference
        self.agglomerated_reference_mismatch_dict = agglomerated_reference_mismatch_dict
        self.reference_mismatches_to_agglomerated_reference = reference_mismatches_to_agglomerated_reference


def remap_queries(reference_input_item,
                  processed_reference_dict,
                  aligned_reference_dict,
                  agglomerated_aligned_query_dict):
    input_reference_name, agglomerated_reference_priority = reference_input_item

    with lock:
        if agglomerated_reference_priority >= processed_reference_dict[input_reference_name]:
            # The reference sequence has already been processed,
            # as it mapped to another reference sequence that had been processed.
            return None

        processed_reference_dict[input_reference_name] = agglomerated_reference_priority

    aligned_reference = aligned_reference_dict[input_reference_name]
    agglomerated_aligned_reference = AlignedTarget(aligned_reference.seq_string, name=input_reference_name)

    # Track the references processed within this function call.
    presently_processed_reference_names = [input_reference_name]

    remapping_stack = deque()
    remapping_stack.append(RemappingItem(input_reference_name, aligned_reference, 0, {}, []))

    while remapping_stack:
        remapping_item = remapping_stack.pop()
        # Record mismatches between query sequences and the agglomerated reference sequence,
        # with the coordinate system being nucleotide positions in the agglomerated reference.
        current_reference_start_in_agglomerated_reference = remapping_item.reference_start_in_agglomerated_reference
        agglomerated_reference_mismatch_dict = remapping_item.agglomerated_reference_mismatch_dict
        current_reference_mismatches_to_agglomerated_reference = remapping_item.reference_mismatches_to_agglomerated_reference

        next_remapping_items = []
        for alignment in remapping_item.aligned_reference.alignments:
            alignment_length = alignment.alignment_length
            aligned_query = alignment.aligned_query
            query_name = aligned_query.name
            query_seq_string = aligned_query.seq_string

            if query_name == remapping_item.reference_name:
                # Ignore a sequence mapping to itself.
                continue

            if query_name in presently_processed_reference_names:
                # The query was already processed in this function call,
                # as it mapped to another agglomerated query.
                continue

            # In the next iteration, the current query sequence will be investigated as a reference.
            presently_processed_reference_names.append(query_name)
            with lock:
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

            # Position of the alignment in the coordinate system of the agglomerated reference sequence.
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
                if agglomerated_reference_index in agglomerated_reference_mismatch_dict:
                    agglomerated_reference_nucleotide = agglomerated_reference_mismatch_dict[agglomerated_reference_index]
                    query_nucleotide = query_seq_string[alignment_index]
                    if agglomerated_reference_nucleotide == query_nucleotide:
                        continue
                else:
                    agglomerated_reference_nucleotide = current_reference_nucleotide
                    agglomerated_reference_mismatch_dict[agglomerated_reference_index] = agglomerated_reference_nucleotide
                query_mismatches_to_agglomerated_reference.append((agglomerated_reference_index,
                                                                   agglomerated_reference_nucleotide))
                query_mismatches_to_agglomerated_reference_in_alignment_frame.append((alignment_index,
                                                                                      agglomerated_reference_nucleotide))

            # In addition to mismatches between the query and current reference
            # that are also mismatches to the agglomerated reference,
            # record mismatches between the query and agglomerated reference
            # at positions where the current reference matches the agglomerated reference.
            query_mismatch_to_current_reference_in_agglomerated_reference_frame_indices = [alignment_index + alignment_start_in_agglomerated_reference
                                                                                           for alignment_index, _
                                                                                           in query_mismatches_to_current_reference_in_alignment_frame]
            for agglomerated_reference_index, agglomerated_reference_nucleotide in current_reference_mismatches_to_agglomerated_reference:
                if agglomerated_reference_index in query_mismatch_to_current_reference_in_agglomerated_reference_frame_indices:
                    # The mismatch position has already been considered,
                    # as there is a mismatch between the query and current reference at this position as well.
                    continue
                if alignment_start_in_agglomerated_reference <= agglomerated_reference_index < alignment_end_in_agglomerated_reference:
                    # Ignore mismatches that lie outside the bounds of the alignment between the query and current reference.
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

            next_remapping_items.append(RemappingItem(query_name,
                                                      aligned_reference_dict[query_name],
                                                      alignment_start_in_agglomerated_reference,
                                                      dict(agglomerated_reference_mismatch_dict.items()),
                                                      [mismatch_tuple for mismatch_tuple in current_reference_mismatches_to_agglomerated_reference]))

        for next_remapping_item in next_remapping_items[::-1]:
            remapping_stack.append(next_remapping_item)

    return agglomerated_aligned_reference
