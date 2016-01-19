# -*- coding: utf-8 -*-

import os
import sys
import glob
import string

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

single_default = "tnf"
merged_default = "tnf-cov"

if not os.path.exists(os.path.join(clustering_configs_dir, 'single', single_default)):
    print "Error: The default clustering configuration file for single runs, '%s',\n\
       is missing from data/clusterconfigs dir! I can't fix this!." % (single_default)
    sys.exit()

if not os.path.exists(os.path.join(clustering_configs_dir, 'merged', merged_default)):
    print "Error: The default clustering configuration file for merged runs, '%s',\n\
       is missing from data/clusterconfigs dir! I can't fix this!." % (merged_default)
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

codon_to_AA = {'ATA':'Ile', 'ATC':'Ile', 'ATT':'Ile', 'ATG':'Met',
               'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACT':'Thr',
               'AAC':'Asn', 'AAT':'Asn', 'AAA':'Lys', 'AAG':'Lys',
               'AGC':'Ser', 'AGT':'Ser', 'AGA':'Arg', 'AGG':'Arg',
               'CTA':'Leu', 'CTC':'Leu', 'CTG':'Leu', 'CTT':'Leu',
               'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCT':'Pro',
               'CAC':'His', 'CAT':'His', 'CAA':'Gln', 'CAG':'Gln',
               'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGT':'Arg',
               'GTA':'Val', 'GTC':'Val', 'GTG':'Val', 'GTT':'Val',
               'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCT':'Ala',
               'GAC':'Asp', 'GAT':'Asp', 'GAA':'Glu', 'GAG':'Glu',
               'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGT':'Gly',
               'TCA':'Ser', 'TCC':'Ser', 'TCG':'Ser', 'TCT':'Ser',
               'TTC':'Phe', 'TTT':'Phe', 'TTA':'Leu', 'TTG':'Leu',
               'TAC':'Tyr', 'TAT':'Tyr', 'TAA':'STP', 'TAG':'STP',
               'TGC':'Cys', 'TGT':'Cys', 'TGA':'STP', 'TGG':'Trp'}

pretty_names = {}

def get_pretty_name(key):
    if pretty_names.has_key(key):
        return pretty_names[key]
    else:
        return key
