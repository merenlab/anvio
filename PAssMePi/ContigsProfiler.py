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
from Oligotyping.lib.entropy import entropy
from Oligotyping.utils.utils import Progress
from Oligotyping.utils.utils import Run 
from Oligotyping.utils.utils import pretty_print as pp 


class BamProfiler:
    def __init__(self, sam_file_path):
        self.bam = None
        self.sam_file_path = args.input_file
        self.Sorted = args.sorted # we're masking sorted
        self.indexed = args.indexed
        self.references = {}
        self.mapped = None
        self.contigs_of_interest = None

        self.list_contigs_and_exit = args.list_contigs

        if args.contigs:
            if os.path.exists(args.contigs):
                self.contigs_of_interest = [c.strip() for c in open(args.contigs).readlines() if c.strip() and not c.startswith('#')]
            else:
                self.contigs_of_interest = [c.strip() for c in args.contigs.split(',')] if args.contigs else None

        if self.Sorted:
            self.sam_file_sorted_path = self.sam_file_path
        else: 
            self.sam_file_sorted_path = None

        self.progress = Progress()
        self.comm = Run()

        self.init()

        if self.contigs_of_interest:
            self.comm.info('Contigs selected for analysis', ','.join(self.references.keys()))


        self.column_entropy_profiles = {}
        self.contig_entropy_profiles = {}
        self.generate_column_entropy_profile()
        self.generate_report()


    def init(self):
        if not self.Sorted:
            self.progress.new('SORT')
            self.progress.update('Sorting BAM File... May take a while depending on the size.')
            pysam.sort(self.sam_file_path, self.sam_file_path + '-sorted')
            self.sam_file_sorted_path = self.sam_file_path + '-sorted.bam'
            self.indexed = False
            self.progress.end()
            self.comm.info('Sorted BAM File', self.sam_file_sorted_path)

        if not self.indexed:
            self.progress.new('INDEX')
            self.progress.update('Analyzing BAM File.')
            pysam.index(self.sam_file_sorted_path)
            self.progress.end()
            self.comm.info('Indexed BAM File', self.sam_file_sorted_path)

        self.progress.new('READ')
        self.progress.update('Reading BAM File')
        self.bam = pysam.Samfile(self.sam_file_sorted_path, 'rb')
        self.progress.end()

        self.comm.info('Input BAM file', self.sam_file_sorted_path)

        references = self.bam.references
        lengths = self.bam.lengths

        if self.contigs_of_interest:
            indexes = [references.index(r) for r in self.contigs_of_interest if r in references]
            references = [references[i] for i in indexes]
            lengths = [lengths[i] for i in indexes]
 
        for i in range(0, len(references)):
            ref = references[i]
            length = lengths[i]
            self.references[ref] = {"length": length}

        self.num_reads_mapped = self.bam.mapped

        self.comm.info('References', ','.join(['"%s (%s)"' % (ref, pp(int(self.references[ref]['length']))) \
                                    for ref in self.references.keys()[0:5]])) 
        self.comm.info('Total reads mapped', pp(int(self.num_reads_mapped)))

        self.progress.new('Analyzing Coverage Stats')
        bad_refs = []
        for reference in self.references:
            self.progress.update('"%s" ...' % (reference))
            coverage = []

            for pileupcolumn in self.bam.pileup(reference):
                coverage.append(pileupcolumn.n)

            if not coverage:
                bad_refs.append(reference)
            else:
                self.references[reference]['min_coverage'] = numpy.median(coverage) / 2
                self.references[reference]['max_coverage'] = numpy.max(coverage)
                self.references[reference]['median_coverage'] = numpy.median(coverage)
                self.references[reference]['mean_coverage'] = numpy.mean(coverage)
        self.progress.end()

        # FIXME CLEARING BAD CONTIGS
        for reference in bad_refs:
            self.references.pop(reference)

        # BASIC CONTIG STATS
        refs_sorted_by_occurence = sorted([(self.references[ref]['length'], ref) for ref in self.references], reverse=True)

        fields = ["length", "mean_coverage", "median_coverage", "min_coverage", "max_coverage"]
        text_output = ""
        text_output += '\t'.join(["%-30s" % field for field in ["contig"] + fields]) + '\n'
        for tpl in refs_sorted_by_occurence:
            reference = tpl[1]
            ref_obj = self.references[reference]
            text_output += "\t".join(["%-30s" % reference] + ["%-30s" % pp(int(ref_obj[field])) for field in fields]) + "\n"

        self.contigs_basic_stats = self.sam_file_path + '-CONTIGS.txt'
        if not (os.path.exists(self.contigs_basic_stats) and self.contigs_of_interest):
            f = open(self.contigs_basic_stats, 'w').write(text_output)
            self.comm.info('Contigs basic stats', self.contigs_basic_stats)

        if self.list_contigs_and_exit:
            print
            print "%d contigs found in the BAM file:" % (len(self.references))
            print
            print text_output
            sys.exit()
           
