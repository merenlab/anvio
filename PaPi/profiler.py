#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys
import numpy
import pysam
import random
import operator
from PaPi.utils import Progress
from PaPi.utils import Run 
from PaPi.utils import pretty_print as pp 
import PaPi.utils as utils
from PaPi.contig_stats_essential import Essential
from PaPi.contig_stats_auxiliary import Auxiliary
from PaPi.contig_composition import Composition


class BAMProfiler:
    """Creates an über class for BAM file operations"""
    def __init__(self, args):
        self.bam = None
        self.input_file_path = args.input_file
        self.output_directory = args.output_directory
        self.references_dict = {}
        self.mapped = None
        self.contigs_of_interest = None
        self.no_trehading = False
        self.number_of_threads = 4
        self.min_contig_length = args.min_contig_length
        self.min_mean_coverage = args.min_mean_coverage

        self.list_contigs_and_exit = args.list_contigs

        if args.contigs:
            if os.path.exists(args.contigs):
                self.contigs_of_interest = [c.strip() for c in open(args.contigs).readlines() if c.strip() and not c.startswith('#')]
            else:
                self.contigs_of_interest = [c.strip() for c in args.contigs.split(',')] if args.contigs else None

        self.progress = Progress()
        self.comm = Run()

        self.init()
        self.profile()
        #self.report()


    def init(self):
        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = pysam.Samfile(self.input_file_path, 'rb')
        self.progress.end()
        self.comm.info('Input BAM file', self.input_file_path)

        self.progress.new('Init')
        self.progress.update('Initializing the output directory ...')
        self.init_output_directory()
        self.progress.end()
        self.comm.info('Output directory', self.output_directory)

        self.references = self.bam.references
        self.raw_lengths = self.bam.lengths

        try:
            self.num_reads_mapped = self.bam.mapped
        except ValueError:
            raise utils.ConfigError, "It seems the BAM file is not indexed. See 'papi-init-bam' script."

        self.comm.info('Total reads mapped', pp(int(self.num_reads_mapped)))
        self.comm.info('Total number of contigs in file', pp(len(self.references)))

        if self.list_contigs_and_exit:
            print "\nContigs in the file:\n"
            for (reference, length) in zip(self.references, self.raw_lengths):
                print "\t- %s (%s)" % (reference, pp(int(length)))
            print
            sys.exit()

        if self.contigs_of_interest:
            indexes = [self.references.index(r) for r in self.contigs_of_interest if r in self.references]
            self.references = [self.references[i] for i in indexes]
            self.raw_lengths = [self.raw_lengths[i] for i in indexes]
            self.comm.info('Total num contigs selected for analysis', pp(len(self.references)))

        contigs_longer_than_M = set()
        for i in range(0, len(self.references)):
            if self.raw_lengths[i] > self.min_contig_length:
                contigs_longer_than_M.add(i)
        if not len(contigs_longer_than_M):
            raise utils.ConfigError, "0 contigs larger than %s nts." % pp(self.min_contig_length)
        else:
            self.references = [self.references[i] for i in contigs_longer_than_M]
            self.raw_lengths = [self.raw_lengths[i] for i in contigs_longer_than_M]
            self.comm.info('Contigs with raw length longer than M', len(self.references))


    def init_output_directory(self):
        if not self.output_directory:
            self.output_directory = self.input_file_path + '-PaPi-OUTPUT'
        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except:
                self.progress.end()
                raise utils.ConfigError, "Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                                (self.output_directory)
        if not os.access(self.output_directory, os.W_OK):
            self.progress.end()
            raise utils.ConfigError, "You do not have write permission for the output directory: '%s'" % self.output_directory


    def generate_output_destination(self, postfix, directory = False):
        return_path = os.path.join(self.output_dir, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def profile(self):
        """Big deal function"""

        self.progress.new('Profiling the BAM file for Essential Stats')
        # So we start with essential stats. In the section below, we will simply go through each reference (contig),
        # in the BAM file and populate the references dictionary for the first time. There are two major sections,
        # one for no_threading option, and the other with multiple threads.
        if self.no_trehading:
            for i in range(0, len(self.references)):
                reference = self.references[i]
                self.references_dict[reference] = {}

                #fill in basics
                self.progress.update('Essential stats for "%s" (%d of %d) ...' % (reference, i + 1, len(self.references)))
                self.references_dict[reference]['essential'] = Essential(reference, self.bam.pileup(reference)).report()

        else:
            def worker(reference, shared_references_dict):
                shared_references_dict[reference] = Essential(reference, self.bam.pileup(reference)).report()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_references_dict = mp.get_empty_shared_dict()

            # arrange processes
            processes_to_run = []
            for reference in self.references:
                processes_to_run.append((reference, shared_references_dict),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)
            for reference in self.references:
                self.references_dict[reference] = {}
                self.references_dict[reference]['essential'] = shared_references_dict.pop(reference)
        self.progress.end()


        # breath in, breath out. filtering based on M and C starts.
        self.progress.new('Filtering contigs based on min-length')
        # this is important:
        # paired-end libraries with large inserts can cover long areas with large empty areas in between. after
        # analyzing the coverage across each contig, we know the actual lenght of real nucleotides (this information
        # is held by ['essential']['length']). so we will further eliminate contigs that are kinda useless:
        self.progress.update('Screening actual contig lengths ...')
        references_to_discard = set()
        for reference in self.references_dict:
            if self.references_dict[reference]['essential']['length'] < self.min_contig_length:
                references_to_discard.add(reference)

        if len(references_to_discard):
            for reference in references_to_discard:
                self.references_dict.pop(reference)
            self.references = self.references_dict.keys()
            self.progress.end()
            self.comm.info('Total number of contigs after precise M elimination', pp(len(self.references)))
        else:
            self.progress.end()


        self.progress.new('Filtering contigs based on mean coverage')
        # this is also important. here we are going to remove any contig with a mean coverage less than C; mean
        # coverage info is stored in ['essential']['mean_coverage']. the mean coverage does not include areas
        # where zero reads mapped. 
        self.progress.update('Screening coverage for each contig ...')
        references_to_discard = set()
        for reference in self.references_dict:
            if self.references_dict[reference]['essential']['mean_coverage'] < self.min_mean_coverage:
                references_to_discard.add(reference)

        if len(references_to_discard):
            for reference in references_to_discard:
                self.references_dict.pop(reference)
            self.references = self.references_dict.keys()
            self.progress.end()
            self.comm.info('Total number of contigs after C', pp(len(self.references)))
        else:
            self.progress.end()


        # QA/QC is done. Now we go into Auxiliary analyses.
        self.progress.new('Computing auxiliary stats')
        if self.no_trehading:
            for i in range(0, len(self.references)):
                reference = self.references[i]
                # fill in entropy and representatives
                self.progress.update('Auxiliary stats for "%s" (%d of %d) ...' % (reference, i + 1, len(self.references)))
                self.references_dict[reference]['auxiliary'] = Auxiliary(reference,
                                                                         self.bam.pileup(reference),
                                                                         self.references_dict[reference]['essential']\
                                                                         ).report()
        else:
            def worker(reference, shared_references_dict):
                shared_references_dict[reference] = Auxiliary(reference,
                                                              self.bam.pileup(reference),
                                                              self.references_dict[reference]['essential']\
                                                              ).report()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_references_dict = mp.get_empty_shared_dict()

            # arrange processes
            processes_to_run = []
            for reference in self.references:
                processes_to_run.append((reference, shared_references_dict),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)
            for reference in self.references:
                self.references_dict[reference]['auxiliary'] = shared_references_dict.pop(reference)
        self.progress.end()


        # it is time to fill in the tetranucleotide frequency info per contig
        self.progress.new('TNF Stats')
        if self.no_trehading:
            for i in range(0, len(self.references)):
                reference = self.references[i]
                self.progress.update('Computing TNF for "%s" (%d of %d) ...' % (reference, i + 1, len(self.references)))
                self.references_dict[reference]['composition'] = Composition(reference,
                                                                             self.references_dict[reference]['auxiliary']\
                                                                             ).report()
        else:
            def worker(reference, shared_references_dict):
                shared_references_dict[reference] = Composition(reference,
                                                                self.references_dict[reference]['auxiliary']\
                                                                ).report()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_references_dict = mp.get_empty_shared_dict()

            # arrange processes
            processes_to_run = []
            for reference in self.references:
                processes_to_run.append((reference, shared_references_dict),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)
            for reference in self.references:
                self.references_dict[reference]['composition'] = shared_references_dict.pop(reference)
        self.progress.end()

        # Profiling is done.


    def report(self):
        # time to report stuff:
        #
        # TNF matrix. basic stats. representative sequences file (phymbll?). 