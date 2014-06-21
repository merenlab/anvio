#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
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
        self.references_dict = {}
        self.mapped = None
        self.contigs_of_interest = None

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

        self.references = self.bam.references
        self.lengths = self.bam.lengths

        if self.list_contigs_and_exit:
            print "\nContigs in the file:\n"
            for (reference, length) in zip(self.references, self.lengths):
                print "\t- %s (%s)" % (reference, pp(int(length)))
            print
            sys.exit()

        if self.contigs_of_interest:
            indexes = [self.references.index(r) for r in self.contigs_of_interest if r in self.references]
            self.references = [self.references[i] for i in indexes]
            self.lengths = [self.lengths[i] for i in indexes]
            l = len(self.references)
            self.comm.info('Contigs selected for analysis', ', '.join(self.references[0:3]) \
                                                                + ' ... (%d) more' % (l-3) if l > 3 else '')
 
        try:
            self.num_reads_mapped = self.bam.mapped
        except ValueError:
            raise utils.ConfigError, "It seems the BAM file is not indexed. See 'papi-init-bam' script."

        self.comm.info('Total num references', pp(len(self.references))) 
        self.comm.info('Total reads mapped', pp(int(self.num_reads_mapped)))


    def profile(self):
        self.progress.new('Profiling the BAM file')

        # Go through each reference, populate references dictionary...
        for i in range(0, len(self.references)):
            reference = self.references[i]
            self.references_dict[reference] = {}

            #fill in basics
            self.progress.update('Essential stats for "%s" (%d of %d) ...' % (reference, i + 1, len(self.references)))
            self.references_dict[reference]['essential'] = Essential(reference, self.bam.pileup(reference)).report()
            self.references_dict[reference]['essential']['length'] = self.lengths[i]

            # fill in entropy and representatives
            self.progress.update('Auxiliary stats for "%s" (%d of %d) ...' % (reference, i + 1, len(self.references)))
            self.references_dict[reference]['auxiliary'] = Auxiliary(reference, self.bam.pileup(reference),
                                                                     self.references_dict[reference]['essential']).report()

            # fill in tetranucleotide frequency
            self.progress.update('Computing TFN for "%s" (%d of %d) ...' % (reference, i + 1, len(self.references)))
            self.references_dict[reference]['composition'] = Composition(reference, self.references_dict[reference]['auxiliary']).report()


        self.progress.end()

        print self.references_dict
        print len(self.references_dict['contig_8']['essential']['coverage'])
        print len(self.references_dict['contig_8']['auxiliary']['rep_seq'])
        print self.references_dict['contig_8']['essential']['length']

    def report(self):
        # BASIC CONTIG STATS
        refs_sorted_by_occurence = sorted([(self.references[ref]['basics']['length'], ref) for ref in self.references], reverse=True)

        fields = ["length", "mean_coverage", "median_coverage", "min_coverage", "max_coverage"]
        text_output = ""
        text_output += '\t'.join(["%-30s" % field for field in ["contig"] + fields]) + '\n'
        for tpl in refs_sorted_by_occurence:
            reference = tpl[1]
            ref_obj = self.references[reference]
            text_output += "\t".join(["%-30s" % reference] + ["%-30s" % pp(int(ref_obj[field])) for field in fields]) + "\n"

        self.contigs_basic_stats = self.input_file_path + '-CONTIGS.txt'
        if not (os.path.exists(self.contigs_basic_stats) and self.contigs_of_interest):
            f = open(self.contigs_basic_stats, 'w').write(text_output)
            self.comm.info('Contigs basic stats', self.contigs_basic_stats)

        if self.list_contigs_and_exit:
            print
            print "%d contigs found in the BAM file:" % (len(self.references))
            print
            print text_output
            sys.exit()
           
