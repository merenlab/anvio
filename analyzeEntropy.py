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


class SamProfiler:
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
           

    def generate_column_entropy_profile(self):
        ref_keys = self.references.keys()
        for i in range(0, len(ref_keys)):
            reference = ref_keys[i]
            ref_obj = self.references[ref_keys[i]]
            self.progress.new('Analyzing Contig "%s" (%d of %d)' % (reference, i + 1, len(ref_keys)))

            self.column_entropy_profiles[reference] = {}
            
            ep = self.column_entropy_profiles[reference]

            self.progress.update('working on %s (%d of %d) // AC: %d :: MC: %d' % (reference,
                                                                                   i + 1,
                                                                                   len(ref_keys),
                                                                                   ref_obj['mean_coverage'],
                                                                                   ref_obj['min_coverage']))

            for pileupcolumn in self.bam.pileup(reference):
                if pileupcolumn.pos % 500 == 0:
                    self.progress.update('Generating Entropy Profiles: %s (%.2f%%)' % (pp(pileupcolumn.pos),
                                                                                       pileupcolumn.pos * 100.0 / ref_obj['length']))

                column = ''.join([pileupread.alignment.seq[pileupread.qpos] for pileupread in pileupcolumn.pileups])

                ep[pileupcolumn.pos] = ColumnEntropyProfile(column,
                                                            pileupcolumn.pos,
                                                            ref_obj['max_coverage'],
                                                            ref_obj['median_coverage'])

            self.contig_entropy_profiles[reference] = ContigEntropy(ref_obj, ep).summarize()

            self.progress.end()

        self.comm.info('Stuff is processed', True)


    def generate_report(self):
        self.columns_report_file_path = self.sam_file_path + '-COLUMNS-REPORT.txt'
        self.contigs_report_file_path = self.sam_file_path + '-CONTIGS-REPORT.txt'

        columns_report = open(self.columns_report_file_path, 'w')
        contigs_report = open(self.contigs_report_file_path, 'w')

        columns_report.write('%s\n' % ('\t'.join(["contig", "pos", "coverage", "entropy", "entropy_n", "competing_nt", "freq_A", "freq_T", "freq_C", "freq_G", "stars"])))
        for reference in self.column_entropy_profiles:
            self.progress.new('Reporting: "%s"' % (reference))
            ep = self.column_entropy_profiles[reference]

            positions = sorted(ep.keys())

            max_coverage = max([x.coverage for x in ep.values()])

            for i in range(0, len(positions)):
                position = positions[i]
                if position % 500 == 0:
                    self.progress.update('%.2f%%' % (i * 100.0 / len(positions)))
                column = ep[position]
                if column.entropy == 0:
                    continue
                info_line = '%s\t%.7d\t%d\t%.3f\t%.3f\t%s\t%s\t%s' % (reference,
                                                 column.pos,
                                                 column.coverage,
                                                 column.entropy,
                                                 column.normalized_entropy,
                                                 column.competing_nucleotides,
                                                 '\t'.join(['%.2f' % (0.0 if not column.nucleotide_frequencies.has_key(n) else column.nucleotide_frequencies[n]) for n in 'ATCG']),
                                                 '*' * int((column.entropy * 50) * (column.coverage * 1.0 / max_coverage))) 

                columns_report.write(info_line + '\n')
            self.progress.end()
        self.comm.info("Columns Report", self.columns_report_file_path)
        columns_report.close()

        
        contigs_report.write('%s\n' % ('\t'.join(["contig", "average_entropy", "average_normalized_entropy"])))
        for reference in self.references:
            cp = self.contig_entropy_profiles[reference]
            info_line = '%s\t%.4f\t%.4f' % (reference,
                                            cp['average_entropy'],
                                            cp['average_normalized_entropy']) 
            contigs_report.write(info_line + '\n')
        self.comm.info("Contigs Report", self.contigs_report_file_path)
        contigs_report.close()


class ContigEntropy:
    def __init__(self, reference, entropy_profile):
        self.entropy_profile = entropy_profile
        self.reference = reference

    def summarize(self):
        average_entropy = sum(e.entropy for e in self.entropy_profile.values()) * 1000.0 / self.reference['length']
        average_normalized_entropy = sum(e.normalized_entropy for e in self.entropy_profile.values()) * 1000.0 / self.reference['length']
        return {'average_entropy': average_entropy,
                'average_normalized_entropy': average_normalized_entropy}


class ColumnEntropyProfile:
    def __init__(self, column, pos, contig_max_coverage, contig_median_coverage):
        self.pos = pos
        self.contig_max_coverage = contig_max_coverage
        self.contig_median_coverage = contig_median_coverage
        self.coverage = len(column)
        self.nucleotide_counts = {}
        self.nucleotide_frequencies = {}
        self.entropy = 0.0
        self.normalized_entropy = 0.0
        self.competing_nucleotides = ''

        nucleotides = list(set(column))
        self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])
        nucleotides_sorted_by_occurence = [x[0] for x in sorted(self.nucleotide_counts.iteritems(), key=operator.itemgetter(1), reverse=True)]

        if len(nucleotides_sorted_by_occurence) == 1:
            self.entropy = 0
            self.normalized_entropy = 0
            return

        if len(nucleotides_sorted_by_occurence) > 2:
            self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])

            n1 = nucleotides_sorted_by_occurence[0]
            n2 = nucleotides_sorted_by_occurence[1]

            column = ''.join([n for n in column if n in [n1, n2]])

        nucleotides = list(set(column))
        denominator = len(column)
        self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])
        nucleotides_sorted_by_occurence = [x[0] for x in sorted(self.nucleotide_counts.iteritems(), key=operator.itemgetter(1), reverse=True)]
        self.nucleotide_frequencies = dict([(n, self.nucleotide_counts[n] * 1.0 / denominator) for n in nucleotides])
        n1 = nucleotides_sorted_by_occurence[0]
        n2 = nucleotides_sorted_by_occurence[1]
        self.competing_nucleotides = ''.join(sorted([n1, n2]))

        if self.nucleotide_counts[n1] * 1.0 / self.nucleotide_counts[n2] > 10:
            self.entropy = 0
            self.normalized_entropy = 0
            self.competing_nucleotides = ''
            return

        if len(column) < 100:
            self.entropy = entropy(column)
            self.normalized_entropy = 0
            return

        self.competing_nucleotides = ''.join(sorted([n1, n2]))
        self.entropy = entropy(column)
        self.normalized_entropy = self.entropy * (self.coverage * 1.0 / contig_max_coverage)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='BAM file profiler for entropy analysis')
    parser.add_argument('input_file', metavar = 'FILE PATH',
                        help = 'SAM file to analyze')
    parser.add_argument('--sorted', action = 'store_true', default = False,
                        help = 'Flag to define whether BAM file is already sorted')
    parser.add_argument('--indexed', action = 'store_true', default = False,
                        help = 'Flag to define whether BAM file is already indexed')
    parser.add_argument('--list-contigs', action = 'store_true', default = False,
                        help = 'Whend declared, lists contigs in the BAM file and\
                                exits without any further analysis.')
    parser.add_argument('--contigs', default = None,
                        help = 'It is possible to analyze only a group of contigs from\
                                a given BAM file. Contigs of interest can be specified\
                                using a comma separated list, or in a text file where\
                                each line contains a contig name.')

    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print 'No such file: "%s"' % args.input_file
        sys.exit()

    rofiler = SamProfiler(args)
