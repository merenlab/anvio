# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import operator
from scipy import log2 as log
from numpy import sqrt


class ColumnEntropyProfile:
    def __init__(self, column, pos, contig_max_coverage, contig_median_coverage):
        self.pos = pos
        self.contig_max_coverage = contig_max_coverage
        self.coverage = len(column)
        self.nucleotide_counts = {}
        self.nucleotide_frequencies = {}
        self.consensus_nucleotide = None
        self.entropy = 0.0
        self.n2n1ratio = 0
        self.normalized_entropy = 0.0
        self.competing_nucleotides = ''
        self.competing_nucleotide_ratio = 0

        nucleotides = list(set(column))
        self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])
        nucleotides_sorted_by_occurence = [x[0] for x in sorted(self.nucleotide_counts.iteritems(), key=operator.itemgetter(1), reverse=True)]

        if len(nucleotides_sorted_by_occurence) == 1:
            # no variation.
            self.entropy = 0
            self.normalized_entropy = 0
            self.consensus_nucleotide = column[0]
            return

        if len(nucleotides_sorted_by_occurence) > 2:
            self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])

            n1 = nucleotides_sorted_by_occurence[0]
            n2 = nucleotides_sorted_by_occurence[1]
            self.n2n1ratio = self.nucleotide_counts[n2] / self.nucleotide_counts[n1] * 1.0

            column = ''.join([n for n in column if n in [n1, n2]])

        nucleotides = list(set(column))
        denominator = float(len(column))
        self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])
        nucleotides_sorted_by_occurence = [x[0] for x in sorted(self.nucleotide_counts.iteritems(), key=operator.itemgetter(1), reverse=True)]
        self.nucleotide_frequencies = dict([(n, self.nucleotide_counts[n] * 1.0 / denominator) for n in nucleotides])
        n1 = nucleotides_sorted_by_occurence[0]
        n2 = nucleotides_sorted_by_occurence[1]
        self.consensus_nucleotide = n1
        self.competing_nucleotides = ''.join(sorted([n1, n2]))

        if 'N' in self.competing_nucleotides or 'n' in self.competing_nucleotides:
            self.entropy = 0
            self.normalized_entropy = 0
            self.competing_nucleotides = ''
            return

        if self.nucleotide_counts[n1] * 1.0 / self.nucleotide_counts[n2] > 10:
            self.entropy = 0
            self.normalized_entropy = 0
            self.competing_nucleotides = ''
            return

        if len(column) < 4:
            self.entropy = 0
            self.normalized_entropy = 0
            self.competing_nucleotides = ''
            return

        self.competing_nucleotides = ''.join(sorted([n1, n2]))
        self.entropy = entropy(column)
        self.competing_nucleotide_ratio = self.nucleotide_frequencies[self.competing_nucleotides[0]]
        self.normalized_entropy = self.entropy * (self.coverage * 1.0 / contig_median_coverage)
        if self.normalized_entropy > 1:
            self.normalized_entropy = 1


    def generate_report(self):
        self.columns_report_file_path = self.sam_file_path + '-COLUMNS-REPORT.txt'
        self.contigs_report_file_path = self.sam_file_path + '-CONTIGS-REPORT.txt'

        columns_report = open(self.columns_report_file_path, 'w')
        contigs_report = open(self.contigs_report_file_path, 'w')

        columns_report.write('%s\n' % ('\t'.join(["contig", "pos", "coverage", "entropy", "entropy_n", "competing_nt", "competing_nucleotide_ratio", "freq_A", "freq_T", "freq_C", "freq_G", "stars"])))

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

                #if column.entropy == 0:
                #    continue

                info_line = '%s\t%.7d\t%d\t%.3f\t%.3f\t%s\t%.2f\t%s\t%s' % (reference,
                                                 column.pos,
                                                 column.coverage,
                                                 column.entropy,
                                                 column.normalized_entropy,
                                                 column.competing_nucleotides,
                                                 column.competing_nucleotide_ratio,
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


def entropy(l, l_qual = None, expected_qual_score = 40, sqrt_norm = False):
    l = l.upper() 
    
    valid_chars = set(['A', 'T', 'C', 'G', '-'])

    if sqrt_norm:
        l_normalized = ''
        for char in valid_chars:
            l_normalized += char * int(round(sqrt(l.count(char))))
        l = l_normalized


    E_Cs = []
    for char in valid_chars:
        P_C = (l.count(char) * 1.0 / len(l)) + 0.0000000000000000001
        E_Cs.append(P_C * log(P_C))
   
    if l_qual:
        # return weighted entropy
        return -(sum(E_Cs) * (l_qual['mean'] / expected_qual_score))
    else:
        # return un-weighted entropy
        return -(sum(E_Cs))

