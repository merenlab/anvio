#!/usr/bin/env python
# -*- coding: utf-8

import anvio.terminal as terminal

import sys

from copy import deepcopy


# For huge datasets, the default recursion limit of 1,000 can be exceeded.
sys.setrecursionlimit(10000)

class Agglomerator:
    def __init__(self, aligned_query_dict, aligned_target_dict, progress=terminal.Progress()):

        self.aligned_query_dict = aligned_query_dict
        self.aligned_target_dict = aligned_target_dict
        self.agglomerated_aligned_query_dict = deepcopy(self.aligned_query_dict)
        for agglomerated_aligned_query in self.agglomerated_aligned_query_dict.values():
            agglomerated_aligned_query.alignments = []
        self.agglomerated_aligned_reference_dict = {}

        self.progress = progress

        self.processed_reference_names = []
        self.agglomerated_reference_count = 0


    def agglomerate(self, query_can_occur_in_multiple_agglomerations=False):

        self.progress.new("Agglomerating sequence alignments")

        processed_reference_count = 0
        input_reference_count = len(self.aligned_target_dict)
        for reference_name, aligned_target in sorted(self.aligned_target_dict.items(),
                                                     key=lambda aligned_target_item: (-len(aligned_target_item[1].seq),
                                                                                      aligned_target_item[0])):
            self.agglomerated_aligned_reference = deepcopy(aligned_target)
            self.agglomerated_aligned_reference.alignments = []
            self.remap_queries(reference_name, query_can_occur_in_multiple_agglomerations=query_can_occur_in_multiple_agglomerations)
            processed_reference_count += 1

            self.progress.update("%d/%d sequences processed: %d agglomerated references found"
                                 % (processed_reference_count, input_reference_count, self.agglomerated_reference_count))

        for reference_name in self.agglomerated_aligned_reference_dict:
            self.agglomerated_aligned_query_dict.pop(reference_name)

        self.progress.end()


    def remap_queries(self,
                      current_reference_name,
                      recursion_count=0,
                      agglomerated_reference_name=None,
                      processed_reference_names_in_recursion=None,
                      current_reference_start_in_agglomerated_reference=0,
                      agglomerated_reference_mismatch_dict=None,
                      current_reference_mismatches_to_agglomerated_reference=None,
                      query_can_occur_in_multiple_agglomerations=False):

        if recursion_count == 0:
            if current_reference_name in self.processed_reference_names:
                # The reference sequence has already been processed,
                # as it mapped to another reference sequence that had been processed.
                return
            self.processed_reference_names.append(current_reference_name)

            agglomerated_reference_name = current_reference_name
            # Track the references processed in the initial and recursive function calls.
            processed_reference_names_in_recursion = [current_reference_name]
            agglomerated_reference_mismatch_dict = {}
            # Record mismatches between query sequences and the agglomerated reference sequence,
            # with the coordinate system being nucleotide positions in the agglomerated reference.
            current_reference_mismatches_to_agglomerated_reference = []

        aligned_reference = self.aligned_target_dict[current_reference_name]
        for alignment in aligned_reference.alignments:
            alignment_length = alignment.alignment_length
            aligned_query = alignment.aligned_query
            query_name = aligned_query.name
            query_seq = aligned_query.seq

            if query_name == current_reference_name:
                # Ignore a sequence mapping to itself.
                continue

            if query_name in processed_reference_names_in_recursion:
                # The query was already processed as a reference in a previous layer of recursion.
                continue

            if not query_can_occur_in_multiple_agglomerations:
                if query_name in self.processed_reference_names:
                    # The query was already processed as a reference,
                    # but not within the current recursive call tree.
                    continue

            # Get the mismatches between the query and the current reference,
            # with the coordinate system being nucleotide positions in the current reference.
            query_mismatches_to_current_reference_in_alignment_frame = []
            current_reference_seq = alignment.aligned_target.seq
            alignment_start_in_current_reference = alignment.target_start
            current_reference_index = alignment_start_in_current_reference
            for cigartuple in alignment.cigartuples:
                if cigartuple[0] == 8:
                    for incremental_index in range(cigartuple[1]):
                        mismatch_index = current_reference_index + incremental_index
                        current_reference_nucleotide = current_reference_seq[mismatch_index]
                        query_mismatches_to_current_reference_in_alignment_frame.append((mismatch_index - alignment_start_in_current_reference, current_reference_nucleotide))
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
                    query_nucleotide = query_seq[alignment_index]
                    if agglomerated_reference_nucleotide == query_nucleotide:
                        continue
                else:
                    agglomerated_reference_nucleotide = current_reference_nucleotide
                    agglomerated_reference_mismatch_dict[agglomerated_reference_index] = agglomerated_reference_nucleotide
                query_mismatches_to_agglomerated_reference.append((agglomerated_reference_index, agglomerated_reference_nucleotide))
                query_mismatches_to_agglomerated_reference_in_alignment_frame.append((alignment_index, agglomerated_reference_nucleotide))

            # In addition to mismatches between the query and current reference
            # that are also mismatches to the agglomerated reference,
            # record mismatches between the query and agglomerated reference
            # at positions where the current reference matches the agglomerated reference.
            query_mismatch_to_current_reference_in_agglomerated_reference_frame_indices = [
                alignment_index + alignment_start_in_agglomerated_reference for alignment_index, _ in query_mismatches_to_current_reference_in_alignment_frame]
            for agglomerated_reference_index, agglomerated_reference_nucleotide in current_reference_mismatches_to_agglomerated_reference:
                if agglomerated_reference_index in query_mismatch_to_current_reference_in_agglomerated_reference_frame_indices:
                    # The mismatch position has already been considered,
                    # as there is a mismatch between the query and current reference at this position as well.
                    continue
                if alignment_start_in_agglomerated_reference <= agglomerated_reference_index < alignment_end_in_agglomerated_reference:
                    # Ignore mismatches that lie outside the bounds of the alignment between the query and current reference.
                    query_mismatches_to_agglomerated_reference.append(
                        (agglomerated_reference_index, agglomerated_reference_nucleotide))
                    query_mismatches_to_agglomerated_reference_in_alignment_frame.append(
                        (agglomerated_reference_index - alignment_start_in_agglomerated_reference, agglomerated_reference_nucleotide))

            # In the next layer of recursion, the current query sequence is investigated as a reference.
            processed_reference_names_in_recursion.append(query_name)
            self.processed_reference_names.append(query_name)

            self.remap_queries(query_name,
                               recursion_count=recursion_count + 1,
                               agglomerated_reference_name=agglomerated_reference_name,
                               processed_reference_names_in_recursion=processed_reference_names_in_recursion,
                               current_reference_start_in_agglomerated_reference=alignment_start_in_agglomerated_reference,
                               agglomerated_reference_mismatch_dict=deepcopy(agglomerated_reference_mismatch_dict),
                               current_reference_mismatches_to_agglomerated_reference=deepcopy(query_mismatches_to_agglomerated_reference),
                               query_can_occur_in_multiple_agglomerations=query_can_occur_in_multiple_agglomerations)

            # Change the properties of the alignment to reflect remapping to the agglomerated reference.
            agglomerated_aligned_query = self.agglomerated_aligned_query_dict[query_name]
            agglomerated_alignment = deepcopy(alignment)
            agglomerated_aligned_query.alignments.append(agglomerated_alignment)
            self.agglomerated_aligned_reference.alignments.append(agglomerated_alignment)

            agglomerated_alignment.aligned_query = agglomerated_aligned_query
            agglomerated_alignment.aligned_target = self.agglomerated_aligned_reference
            agglomerated_alignment.target_start = alignment_start_in_agglomerated_reference
            agglomerated_alignment.target_end = alignment_start_in_agglomerated_reference + alignment_length

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
            agglomerated_alignment.cigartuples = cigartuples

        if recursion_count == 0:
            self.agglomerated_aligned_reference_dict[agglomerated_reference_name] = self.agglomerated_aligned_reference
            self.agglomerated_reference_count += 1
