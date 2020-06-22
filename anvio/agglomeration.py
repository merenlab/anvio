#!/usr/bin/env python
# -*- coding: utf-8

import anvio.fastalib
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal

import argparse
import os
import pandas as pd
import pysam
import re
import sys
import tempfile

from anvio.errors import ConfigError, FilesNPathsError

from copy import deepcopy


sys.setrecursionlimit(10000)

class Agglomerator:
    def __init__(
        self,
        input_fasta_path=None,
        input_bam_path=None,
        output_fasta_path=None,
        output_bam_path=None,
        replicates_path=None,
        order_by_replicate_abundance=True,
        max_possible_alignments=2,
        sort_index_bam_output=True,
        run=None,
        progress=None,
        verbose=False):

        self.input_fasta_path = input_fasta_path
        self.input_bam_path = input_bam_path
        self.output_bam_path = output_bam_path
        self.replicates_path = replicates_path
        self.order_by_replicate_abundance = order_by_replicate_abundance
        self.max_possible_alignments = max_possible_alignments
        self.sort_index_bam_output = sort_index_bam_output
        self.output_fasta_path = output_fasta_path
        if run:
            self.run = run
        else:
            self.run = terminal.Run()
        if progress:
            self.progress = progress
        else:
            self.progress = terminal.Progress()
            self.progress.new("Agglomerating sequence alignments")
        self.verbose = verbose

        self.make_entries_for_replicate_seqs = False
        self.replicate_name_dict = None

        self.processed_ref_names = []
        self.ref_seq = None

        self.sanity_check()

        self.init()


    def sanity_check(self):
        if not filesnpaths.is_file_exists(self.input_fasta_path, dont_raise=True):
            raise FilesNPathsError(
                "Your `input_fasta_path`, %s, does not exist." % self.input_fasta_path)

        if not filesnpaths.is_file_exists(self.input_bam_path, dont_raise=True):
            raise FilesNPathsError(
                "Your `input_bam_path`, %s, does not exist." % self.input_bam_path)

        if filesnpaths.is_file_exists(self.output_bam_path, dont_raise=True):
            raise FilesNPathsError(
                "Your `output_bam_path`, %s, already exists." % self.output_bam_path)

        if filesnpaths.is_file_exists(self.output_fasta_path, dont_raise=True):
            raise FilesNPathsError(
                "Your `output_fasta_path`, %s, already exists." % self.output_fasta_path)

        if self.replicates_path:
            if not filesnpaths.is_file_exists(self.replicates_path, dont_raise=True):
                raise FilesNPathsError(
                    "Your `replicates_path`, %s, does not exist." % self.replicates_path)

        if self.max_possible_alignments < 2:
            raise ConfigError(
                "Your value for `max_possible_alignments`, %d, "
                "was less than the required minimum value of 2." % self.max_possible_alignments)


    def init(self):
        if self.replicates_path is None and self.order_by_replicate_abundance is True:
            self.order_by_replicate_abundance = False

        self.input_fasta = anvio.fastalib.SequenceSource(self.input_fasta_path, lazy_init=False)
        self.ref_seq_count = self.input_fasta.total_seq

        self.input_bam = pysam.AlignmentFile(self.input_bam_path, 'rb')

        if self.sort_index_bam_output:
            self.raw_output_bam_path = os.path.splitext(self.output_bam_path)[0] + "-RAW.bam"
        else:
            self.raw_output_bam_path = self.output_bam_path
        self.raw_output_bam = pysam.AlignmentFile(
            self.raw_output_bam_path, 'wb', header=self.input_bam.header.as_dict())

        self.output_fasta = anvio.fastalib.FastaOutput(self.output_fasta_path)

        if not self.replicates_path:
            return

        replicates_df = pd.read_csv(self.replicates_path, sep='\t', header=0)
        if len(replicates_df.columns) == 2:
            replicates_df.columns = ['ref_name', 'replicate_count']
        elif len(replicates_df.columns) == 3:
            replicates_df.columns = ['ref_name', 'replicate_count', 'replicate_names']
            self.make_entries_for_replicate_seqs = True
            self.replicate_name_dict = {}
            for ref_name, replicate_names in zip(
                replicates_df['ref_name'], replicates_df['replicate_names']):
                self.replicate_name_dict[ref_name] = replicate_names.split(',')
        else:
            raise ConfigError(
                "Your `replicates_txt` file, %s, "
                "should have two or three tab-separated columns without headers. "
                "The first two columns should always be sequence name and replicate count. "
                "The third column should be provided "
                "if you need individual entries for each replicate sequence in the output BAM file. "
                "The third column is replicate sequence names. "
                "Each entry is replicate sequence names separated by commas, "
                "e.g., seq1,seq2,seq3" % self.replicates_path)


    def agglomerate(self):

        num_processed_ref_seqs = 0
        while next(self.input_fasta):
            self.ref_seq = self.input_fasta.seq
            self.remap_queries(self.input_fasta.id)

            num_processed_ref_seqs += 1
            self.progress.increment(num_processed_ref_seqs)
            if self.verbose:
                self.progress.update(
                    "%d/%d sequences processed" % (num_processed_ref_seqs, self.ref_seq_count))

        self.raw_output_bam.close()
        self.output_fasta.close()

        if self.sort_index_bam_output:
            pysam.sort('-o', self.output_bam_path, self.raw_output_bam_path)
            pysam.index(self.output_bam_path)
            os.remove(self.raw_output_bam_path)


    def remap_queries(
        self,
        rn_name,
        recursion_count=0,
        r0_name=None,
        r_name_in_recursion_list=None,
        rn_start_in_r0=0,
        m0_dict=None,
        rn_m0_in_r0_list=None):

        """
        ABBREVIATIONS
        =============
        Root reference = r0
        Current reference = rn
        Query = q
        Mismatch relative to current reference = mn
        Mismatch relative to root reference = m0
        Mismatched nucleotide position index in current reference = mni
        Mismatched root nucleotide position index in root reference = m0i
        Mismatched root nucleotide position index in alignment = m0ai
        Nucleotide = nt
        Alignment = a
        Start position index = start
        End position index = end
        """

        if recursion_count == 0:
            if rn_name in self.processed_ref_names:
                # The reference sequence has already been processed,
                # because it mapped to another reference sequence which was processed.
                return
            self.processed_ref_names.append(rn_name)

            r0_name = rn_name
            # Track the references processed in the root and recursive function calls.
            r_name_in_recursion_list = [rn_name]
            m0_dict = {}
            # Record mismatches between query sequences and the root reference sequence,
            # with the coordinate system being nucleotide positions in the root reference.
            rn_m0_in_r0_list = []

        # These are AlignedSegment objects,
        # each recording a mapping of a query to the current reference.
        alis = [ali for ali in self.input_bam.fetch(rn_name)]
        for ali in alis:
            q_name = ali.query_name

            if q_name == rn_name:
                # Ignore a reference sequence mapping to itself.
                continue

            if q_name in r_name_in_recursion_list:
                # The query was already processed as a reference in a previous layer of recursion.
                continue

            if q_name in self.processed_ref_names:
                # The query was already processed as a reference,
                # but not within the current recursive call tree.
                continue

            # Get the mismatches between the query and the current reference,
            # with the coordinate system being nucleotide positions in the current reference.
            q_mn_in_rn_list = self.get_q_mn_in_rn_list(ali.get_tag('MD'))

            # Position of the alignment in the coordinate system of the root reference sequence.
            a_start_in_r0 = rn_start_in_r0 + ali.reference_start

            # Add newly found mismatch positions to lists
            # relating the positions of mismatches in the root reference
            # to the nucleotide in the root reference.
            # The first list uses positions in the root reference.
            q_m0_in_r0_list = []
            # The second dictionary uses positions in the current reference.
            q_m0_in_a_list = []
            for mni, mn_nt in q_mn_in_rn_list:
                m0i = mni + rn_start_in_r0
                if m0i in m0_dict:
                    m0_nt = m0_dict[m0i]
                else:
                    m0_nt = mn_nt
                    m0_dict[m0i] = m0_nt
                q_m0_in_r0_list.append((m0i, m0_nt))
                q_m0_in_a_list.append((mni, m0_nt))

            # Find all mismatches between the query and root reference.
            # Mismatches between the query and the current reference were just found,
            # but in previous recursive layers,
            # other mismatches between references and the root reference were found.
            # Ignore mismatches that lie outside the bounds of the current alignment.
            q_m0i_list = [t[0] for t in q_m0_in_r0_list]
            for rn_m0i, rn_m0_nt in rn_m0_in_r0_list:
                if rn_m0i in q_m0i_list:
                    # The mismatch position has already been considered.
                    # Multiple mutations in different recursive layers
                    # separate the query from the root reference at this position.
                    continue
                if rn_m0i >= a_start_in_r0:
                    q_m0_in_r0_list.append((rn_m0i, rn_m0_nt))
                    q_m0_in_a_list.append((rn_m0i - a_start_in_r0, rn_m0_nt))

            # In the next recursion, the query sequence is investigated as a reference.
            r_name_in_recursion_list.append(q_name)
            self.processed_ref_names.append(q_name)

            self.remap_queries(
                q_name,
                recursion_count=recursion_count + 1,
                r0_name=r0_name,
                r_name_in_recursion_list=r_name_in_recursion_list,
                rn_start_in_r0=a_start_in_r0,
                m0_dict=deepcopy(m0_dict),
                rn_m0_in_r0_list=deepcopy(q_m0_in_r0_list))

            # Change the alignment object to reflect the remapping to the root reference.
            ali.reference_name = r0_name
            ali.reference_start = a_start_in_r0
            ali.set_tag('NM', len(q_m0_in_r0_list))
            ali.set_tag('XM', len(q_m0_in_r0_list))

            md_tag = ''
            # Sort the mismatches by position.
            q_m0_in_a_list.sort(key=lambda t: t[0])
            prev_m0ai = -1
            prev_m0_nt = ''
            for m0ai, m0_nt in q_m0_in_a_list:
                md_tag += str(m0ai - prev_m0ai - 1) + m0_nt
                prev_m0ai = m0ai
                prev_m0_nt = m0_nt
            md_tag += str(ali.alen - prev_m0ai - 1)
            ali.set_tag('MD', md_tag)

            if self.make_entries_for_replicate_seqs:
                # The user wanted an identical entry in the BAM file for each replicate sequence.
                for replicate_name in self.replicate_name_dict[q_name]:
                    ali.query_name = replicate_name
                    self.raw_output_bam.write(ali)
            else:
                self.raw_output_bam.write(ali)

        if recursion_count == 0:
            # Write a FASTA file of root, or "seed", sequences.
            self.output_fasta.write_id(rn_name)
            self.output_fasta.write_seq(self.ref_seq, split=False)

        return


    def get_q_mn_in_rn_list(self, md_tag):
        # The format of an MD tag is [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
        # If two mismatches or deletions are adjacent
        # without a run of identical bases between them,
        # a ‘0’ (indicating a 0-length run) separates them.
        md_parts = re.split('(\D+)', md_tag)
        q_mn_in_rn_list = []

        mi = int(md_parts[0])

        iterator = iter(md_parts[1:])
        for s in iterator:
            q_mn_in_rn_list.append((mi, s))
            mi += int(next(iterator)) + 1
        return q_mn_in_rn_list