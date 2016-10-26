# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

import os
import sys
import glob
import string

from collections import Counter

import anvio

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


clustering_configs_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/clusterconfigs')
clustering_configs = {}

default_port_number = 8080

blank_default = "tnf-splits"
single_default = "tnf"
merged_default = "tnf-cov"
pan_default="presence-absence"

max_num_splits_for_hierarchical_clustering = 20000

# default methods for hierarchical cluster analyses
distance_metric_default = 'euclidean'
linkage_method_default = 'ward'

# this is to have a common language across multiple modules when genomes (whether they are MAGs,
# SAGs, or isolate genomes):
essential_genome_info = ['gc_content', 'num_contigs', 'num_splits', 'total_length', 'num_genes', 'percent_complete', 'percent_redundancy',
                         'genes_are_called', 'avg_gene_length', 'num_genes_per_kb', ]


for run_type_and_default_config_tuples in [('single', single_default), ('merged', merged_default), ('blank', blank_default)]:
    run_type, default_config = run_type_and_default_config_tuples
    if not os.path.exists(os.path.join(clustering_configs_dir, run_type, default_config)):
        print "Error: The default clustering configuration file for %s runs, '%s',\n\
       is missing from data/clusterconfigs dir! I don't know how this happened,\n\
       but I can't fix this! Anvi'o needs an adult :(" % (run_type, default_config)
        sys.exit()

for dir in [d.strip('/').split('/')[-1] for d in glob.glob(os.path.join(clustering_configs_dir, '*/'))]:
    clustering_configs[dir] = {}
    for config in glob.glob(os.path.join(clustering_configs_dir, dir, '*')):
        clustering_configs[dir][os.path.basename(config)] = config

IS_ESSENTIAL_FIELD = lambda f: (not f.startswith('__')) and (f not in ["contig", "GC_content", "length"])
IS_AUXILIARY_FIELD = lambda f: f.startswith('__')
allowed_chars = string.ascii_letters + string.digits + '_' + '-' + '.'
digits = string.digits
complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV',\
                               'tgcayrkmvhdbTGCAYRKMVHDB')

nucleotides = 'ATCGN'

AA_to_single_letter_code = Counter({'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
                                    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
                                    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                                    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'STP': '*',
                                    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y',
                                    'Val': 'V'})

codon_to_AA = Counter({'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
                       'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
                       'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                       'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
                       'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
                       'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
                       'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                       'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
                       'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
                       'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
                       'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                       'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
                       'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
                       'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
                       'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'STP', 'TAG': 'STP',
                       'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'STP', 'TGG': 'Trp'})

codon_to_AA_RC = Counter({'AAA': 'Phe', 'AAC': 'Val', 'AAG': 'Leu', 'AAT': 'Ile',
                          'ACA': 'Cys', 'ACC': 'Gly', 'ACG': 'Arg', 'ACT': 'Ser',
                          'AGA': 'Ser', 'AGC': 'Ala', 'AGG': 'Pro', 'AGT': 'Thr',
                          'ATA': 'Tyr', 'ATC': 'Asp', 'ATG': 'His', 'ATT': 'Asn',
                          'CAA': 'Leu', 'CAC': 'Val', 'CAG': 'Leu', 'CAT': 'Met',
                          'CCA': 'Trp', 'CCC': 'Gly', 'CCG': 'Arg', 'CCT': 'Arg',
                          'CGA': 'Ser', 'CGC': 'Ala', 'CGG': 'Pro', 'CGT': 'Thr',
                          'CTA': 'STP', 'CTC': 'Glu', 'CTG': 'Gln', 'CTT': 'Lys',
                          'GAA': 'Phe', 'GAC': 'Val', 'GAG': 'Leu', 'GAT': 'Ile',
                          'GCA': 'Cys', 'GCC': 'Gly', 'GCG': 'Arg', 'GCT': 'Ser',
                          'GGA': 'Ser', 'GGC': 'Ala', 'GGG': 'Pro', 'GGT': 'Thr',
                          'GTA': 'Tyr', 'GTC': 'Asp', 'GTG': 'His', 'GTT': 'Asn',
                          'TAA': 'Leu', 'TAC': 'Val', 'TAG': 'Leu', 'TAT': 'Ile',
                          'TCA': 'STP', 'TCC': 'Gly', 'TCG': 'Arg', 'TCT': 'Arg',
                          'TGA': 'Ser', 'TGC': 'Ala', 'TGG': 'Pro', 'TGT': 'Thr',
                          'TTA': 'STP', 'TTC': 'Glu', 'TTG': 'Gln', 'TTT': 'Lys'})

pretty_names = {}

def get_pretty_name(key):
    if key in pretty_names:
        return pretty_names[key]
    else:
        return key
