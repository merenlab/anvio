# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

import os
import sys
import glob
import numpy
import string

from collections import Counter

import anvio

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


# these are the atomic data that are generated for each contig profiled
# based on read recruitment results. anvio/contigops.py has the details:
essential_data_fields_for_anvio_profiles = ['std_coverage',
                                            'mean_coverage',
                                            'mean_coverage_Q2Q3',
                                            'detection',
                                            'abundance',
                                            'variability']

# this is to distinguish fields that are often useless for clustering ops
# and other purposes
IS_ESSENTIAL_FIELD = lambda f: (not f.startswith('__')) and (f not in ["contig", "GC_content", "length"])

default_pdb_database_path = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PDB.db')
default_modeller_database_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/db')
default_modeller_scripts_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/scripts')

default_interacdome_data_path = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Interacdome')

clustering_configs_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/clusterconfigs')
clustering_configs = {}

default_scgs_taxonomy_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY')
default_scgs_for_taxonomy = ['Ribosomal_S2',
                             'Ribosomal_S3_C',
                             'Ribosomal_S6',
                             'Ribosomal_S7',
                             'Ribosomal_S8',
                             'Ribosomal_S9',
                             'Ribosomal_S11',
                             'Ribosomal_S20p',
                             'Ribosomal_L1',
                             'Ribosomal_L2',
                             'Ribosomal_L3',
                             'Ribosomal_L4',
                             'Ribosomal_L6',
                             'Ribosomal_L9_C',
                             'Ribosomal_L13',
                             'Ribosomal_L16',
                             'Ribosomal_L17',
                             'Ribosomal_L20',
                             'Ribosomal_L21p',
                             'Ribosomal_L22',
                             'ribosomal_L24',
                             'Ribosomal_L27A']
default_hmm_source_for_scg_taxonomy = set(["Bacteria_71"])

default_trna_taxonomy_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/TRNA_TAXONOMY')
default_anticodons_for_taxonomy = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
                                   'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT',
                                   'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTC', 'CTG',
                                   'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA',
                                   'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG',
                                   'TAT', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTC', 'TTG',
                                   'TTT']
default_hmm_source_for_trna_genes = set(["Transfer_RNAs"])

# The following block of constants are used in the tRNA-seq workflow.
TRNA_FEATURE_NAMES = ['trna_his_position_0',
                      'acceptor_stem',
                      'fiveprime_acceptor_stem_sequence',
                      'position_8',
                      'position_9',
                      'd_arm',
                      'd_stem',
                      'fiveprime_d_stem_sequence',
                      'd_loop',
                      'threeprime_d_stem_sequence',
                      'position_26',
                      'anticodon_arm',
                      'anticodon_stem',
                      'fiveprime_anticodon_stem_sequence',
                      'anticodon_loop',
                      'threeprime_anticodon_stem_sequence',
                      'v_loop',
                      't_arm',
                      't_stem',
                      'fiveprime_t_stem_sequence',
                      't_loop',
                      'threeprime_t_stem_sequence',
                      'threeprime_acceptor_stem_sequence',
                      'discriminator',
                      'threeprime_terminus']
TRNA_SEED_FEATURE_THRESHOLD_CHOICES = TRNA_FEATURE_NAMES[TRNA_FEATURE_NAMES.index('acceptor_stem'): TRNA_FEATURE_NAMES.index('anticodon_loop') + 1]
TRNASEQ_CHECKPOINTS = ('profile', 'normalize', 'map_fragments', 'substitutions', 'indels')

default_port_number = int(os.environ['ANVIO_PORT']) if 'ANVIO_PORT' in os.environ else 8080

blank_default = "tnf"
single_default = "tnf"
merged_default = "tnf-cov"
pan_default = "presence-absence"
trnaseq_default = "cov"

default_gene_caller = "prodigal"

# see https://github.com/merenlab/anvio/issues/1358
gene_call_types = {'CODING': 1,
                   'NONCODING': 2,
                   'UNKNOWN': 3}

max_num_items_for_hierarchical_clustering = 20000

# max coverage depth to read from BAM files using pysam.
# this parameter also can be set later using command line parameters
# we use uint16 as dtype for numpy arrays when we work on & store coverages
# which has limit of 65536, so this constant needs to be smaller than that.
# If you change this value please change all dtypes.
max_depth_for_coverage = 60000

# default methods for hierarchical cluster analyses
distance_metric_default = 'euclidean'
linkage_method_default = 'ward'


# Whether a cigarstring operation consumes the read, reference, or both
#
#Here are the possible bam operations.
#
#    M       BAM_CMATCH      0
#    I       BAM_CINS        1
#    D       BAM_CDEL        2
#    N       BAM_CREF_SKIP   3
#    S       BAM_CSOFT_CLIP  4
#    H       BAM_CHARD_CLIP  5
#    P       BAM_CPAD        6
#    =       BAM_CEQUAL      7
#    X       BAM_CDIFF       8
#
#Notes
#=====
#- A description of what possible cigar operations are possible, see
#  https://imgur.com/a/fiQZXNg, which comes from here:
#  https://samtools.github.io/hts-specs/SAMv1.pdf
cigar_consumption = numpy.array([
    (1, 1),
    (1, 0),
    (0, 1),
    (0, 1),
    (1, 0),
    (0, 0),
    (0, 0),
    (1, 1),
    (1, 1),
])

# this is to have a common language across multiple modules when genomes (whether they are MAGs,
# SAGs, or isolate genomes):
essential_genome_info = ['gc_content', 'num_contigs', 'num_splits', 'total_length', 'num_genes', 'percent_completion', 'percent_redundancy',
                         'genes_are_called', 'avg_gene_length', 'num_genes_per_kb', ]

levels_of_taxonomy = ["t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
levels_of_taxonomy_unknown = {"t_domain": 'Unknown_domains',
                              "t_phylum": 'Unknown_phyla',
                              "t_class": 'Unknown_classes',
                              "t_order": 'Unknown_orders',
                              "t_family": 'Unknown_families',
                              "t_genus": 'Unknown_genera',
                              "t_species": 'Unknown_species'}

for run_type, default_config in [('single', single_default),
                                 ('merged', merged_default),
                                 ('trnaseq', trnaseq_default),
                                 ('blank', blank_default)]:
    if not os.path.exists(os.path.join(clustering_configs_dir, run_type, default_config)):
        print()
        print(f"Error: Although there is a run type defined in the anvi'o constants for \n"
              f"       '{run_type}', the default clustering configuration file for it, namely \n"
              f"       '{default_config}', is missing from the 'anvio/data/clusterconfigs' dir. \n"
              f"       If you are a developer and getting this error, please make sure the file \n"
              f"       is in anvi'o distribution. If you are a user and getting this error, it \n"
              f"       something went terribly wrong with your installation :(\n")
        sys.exit()

for dir in [d.strip('/').split('/')[-1] for d in glob.glob(os.path.join(clustering_configs_dir, '*/'))]:
    clustering_configs[dir] = {}
    for config in glob.glob(os.path.join(clustering_configs_dir, dir, '*')):
        clustering_configs[dir][os.path.basename(config)] = config

allowed_chars = string.ascii_letters + string.digits + '_' + '-' + '.'
digits = string.digits
complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')

unambiguous_nucleotides = set(list('ATCG'))
nucleotides = sorted(list(unambiguous_nucleotides)) + ['N']

WC_BASE_PAIRS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}
# In tRNA, wobble base pairing, including G/U, is common
WC_PLUS_WOBBLE_BASE_PAIRS = {
    'A': ('T', ),
    'T': ('A', 'G'),
    'C': ('G', ),
    'G': ('C', 'T')
}

AA_atomic_composition = {'Ala': Counter({"C":3,  "H":7,  "N":1, "O":2, "S":0}),
                         'Arg': Counter({"C":6,  "H":14, "N":4, "O":2, "S":0}),
                         'Asn': Counter({"C":4,  "H":8,  "N":2, "O":3, "S":0}),
                         'Asp': Counter({"C":4,  "H":7,  "N":1, "O":4, "S":0}),
                         'Cys': Counter({"C":3,  "H":7,  "N":1, "O":2, "S":1}),
                         'Gln': Counter({"C":5,  "H":10, "N":2, "O":3, "S":0}),
                         'Glu': Counter({"C":5,  "H":9,  "N":1, "O":4, "S":0}),
                         'Gly': Counter({"C":2,  "H":5,  "N":1, "O":2, "S":0}),
                         'His': Counter({"C":6,  "H":9,  "N":3, "O":2, "S":0}),
                         'Ile': Counter({"C":6,  "H":13, "N":1, "O":2, "S":0}),
                         'Leu': Counter({"C":6,  "H":13, "N":1, "O":2, "S":0}),
                         'Lys': Counter({"C":6,  "H":14, "N":2, "O":2, "S":0}),
                         'Met': Counter({"C":5,  "H":11, "N":1, "O":2, "S":1}),
                         'Phe': Counter({"C":9,  "H":11, "N":1, "O":2, "S":0}),
                         'Pro': Counter({"C":5,  "H":9,  "N":1, "O":2, "S":0}),
                         'Ser': Counter({"C":3,  "H":7,  "N":1, "O":3, "S":0}),
                         'Thr': Counter({"C":4,  "H":9,  "N":1, "O":3, "S":0}),
                         'Trp': Counter({"C":11, "H":12, "N":2, "O":2, "S":0}),
                         'Tyr': Counter({"C":9,  "H":11, "N":1, "O":3, "S":0}),
                         'Val': Counter({"C":5,  "H":11, "N":1, "O":2, "S":0})}

# taken from http://prowl.rockefeller.edu/aainfo/volume.htm
# volume reference: A.A. Zamyatin, Protein Volume in Solution, Prog. Biophys. Mol. Biol. 24(1972)107-123.
# surface area reference: C. Chotia, The Nature of the Accessible and Buried Surfaces in Proteins, J. Mol. Biol., 105(1975)1-14.
AA_geometry = Counter({'Ala': {"volume":88.6,  "area":115},
                       'Arg': {"volume":173.4, "area":225},
                       'Asn': {"volume":111.1, "area":150},
                       'Asp': {"volume":114.1, "area":160},
                       'Cys': {"volume":108.5, "area":135},
                       'Gln': {"volume":138.4, "area":190},
                       'Glu': {"volume":143.8, "area":180},
                       'Gly': {"volume":60.1,  "area":75},
                       'His': {"volume":153.2, "area":195},
                       'Ile': {"volume":166.7, "area":175},
                       'Leu': {"volume":166.7, "area":170},
                       'Lys': {"volume":168.6, "area":200},
                       'Met': {"volume":162.9, "area":185},
                       'Phe': {"volume":189.9, "area":210},
                       'Pro': {"volume":112.7, "area":145},
                       'Ser': {"volume":89.0,  "area":115},
                       'Thr': {"volume":116.1, "area":140},
                       'Trp': {"volume":227.8, "area":255},
                       'Tyr': {"volume":193.6, "area":230},
                       'Val': {"volume":140.0, "area":155}})

AA_to_codons = Counter({'Ala': ['GCA', 'GCC', 'GCG', 'GCT'],
                        'Arg': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
                        'Asn': ['AAC', 'AAT'],
                        'Asp': ['GAC', 'GAT'],
                        'Cys': ['TGC', 'TGT'],
                        'Gln': ['CAA', 'CAG'],
                        'Glu': ['GAA', 'GAG'],
                        'Gly': ['GGA', 'GGC', 'GGG', 'GGT'],
                        'His': ['CAC', 'CAT'],
                        'Ile': ['ATA', 'ATC', 'ATT'],
                        'Leu': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
                        'Lys': ['AAA', 'AAG'],
                        'Met': ['ATG'],
                        'Phe': ['TTC', 'TTT'],
                        'Pro': ['CCA', 'CCC', 'CCG', 'CCT'],
                        'STP': ['TAA', 'TAG', 'TGA'],
                        'Ser': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
                        'Thr': ['ACA', 'ACC', 'ACG', 'ACT'],
                        'Trp': ['TGG'],
                        'Tyr': ['TAC', 'TAT'],
                        'Val': ['GTA', 'GTC', 'GTG', 'GTT']})

AA_to_single_letter_code = Counter({'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
                                    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
                                    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                                    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'STP': '*',
                                    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y',
                                    'Val': 'V'})

amino_acids = sorted(list(AA_to_single_letter_code.keys()))

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

anticodon_to_AA = Counter({'AAA': 'Phe', 'AAC': 'Val', 'AAG': 'Leu', 'AAT': 'Ile',
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

codon_to_codon_RC = Counter({'AAA': 'TTT', 'AAC': 'GTT', 'AAG': 'CTT', 'AAT': 'ATT',
                             'ACA': 'TGT', 'ACC': 'GGT', 'ACG': 'CGT', 'ACT': 'AGT',
                             'AGA': 'TCT', 'AGC': 'GCT', 'AGG': 'CCT', 'AGT': 'ACT',
                             'ATA': 'TAT', 'ATC': 'GAT', 'ATG': 'CAT', 'ATT': 'AAT',
                             'CAA': 'TTG', 'CAC': 'GTG', 'CAG': 'CTG', 'CAT': 'ATG',
                             'CCA': 'TGG', 'CCC': 'GGG', 'CCG': 'CGG', 'CCT': 'AGG',
                             'CGA': 'TCG', 'CGC': 'GCG', 'CGG': 'CCG', 'CGT': 'ACG',
                             'CTA': 'TAG', 'CTC': 'GAG', 'CTG': 'CAG', 'CTT': 'AAG',
                             'GAA': 'TTC', 'GAC': 'GTC', 'GAG': 'CTC', 'GAT': 'ATC',
                             'GCA': 'TGC', 'GCC': 'GGC', 'GCG': 'CGC', 'GCT': 'AGC',
                             'GGA': 'TCC', 'GGC': 'GCC', 'GGG': 'CCC', 'GGT': 'ACC',
                             'GTA': 'TAC', 'GTC': 'GAC', 'GTG': 'CAC', 'GTT': 'AAC',
                             'TAA': 'TTA', 'TAC': 'GTA', 'TAG': 'CTA', 'TAT': 'ATA',
                             'TCA': 'TGA', 'TCC': 'GGA', 'TCG': 'CGA', 'TCT': 'AGA',
                             'TGA': 'TCA', 'TGC': 'GCA', 'TGG': 'CCA', 'TGT': 'ACA',
                             'TTA': 'TAA', 'TTC': 'GAA', 'TTG': 'CAA', 'TTT': 'AAA'})

conserved_amino_acid_groups = {
    'Nonpolar': ['L','V','I','M','C','H','A'],
    'Aromatic': ['F','W','Y'],
    'Bases': ['K','R','H'],
    'Neutral Amines': ['Q', 'N'],
    'Acids': ['D','E'],
    'Polar and Nonpolar': ['H','Y'],
    'Mostly nonpolar': ['S','T'],
    'B': ['B','N','D'],
    'Z': ['Z','Q','E'],
    'J': ['J','L','I'],
    'None': []
}

conserved_amino_acid_groups['N'] = conserved_amino_acid_groups['Neutral Amines'] + ['B']
conserved_amino_acid_groups['D'] = conserved_amino_acid_groups['Acids'] + ['B']
conserved_amino_acid_groups['Q'] = conserved_amino_acid_groups['Neutral Amines'] + ['Z']
conserved_amino_acid_groups['E'] = conserved_amino_acid_groups['Acids'] + ['Z']
conserved_amino_acid_groups['LI'] = conserved_amino_acid_groups['Nonpolar'] + ['J']


amino_acid_property_group = {}
for key in ['A','V','M','C']:
    amino_acid_property_group[key] = 'Nonpolar'
for key in ['F','W']:
    amino_acid_property_group[key] = 'Aromatic'
for key in ['K','R']:
    amino_acid_property_group[key] = 'Bases'
for key in ['H','Y']:
    amino_acid_property_group[key] = 'Polar and Nonpolar'
for key in ['S','T']:
    amino_acid_property_group[key] = 'Mostly nonpolar'
for key in ['G','P','X']:
    amino_acid_property_group[key] = 'None'
amino_acid_property_group['B'] = 'B'
amino_acid_property_group['Z'] = 'Z'
amino_acid_property_group['J'] = 'J'
amino_acid_property_group['N'] = 'N'
amino_acid_property_group['D'] = 'D'
amino_acid_property_group['Q'] = 'Q'
amino_acid_property_group['E'] = 'E'
amino_acid_property_group['L'] = 'LI'
amino_acid_property_group['I'] = 'LI'

codons = sorted(list(set(codon_to_AA.keys())))
coding_codons = [x for x in codons if codon_to_AA[x] != "STP"]

is_synonymous = {}
for i in coding_codons:
    is_synonymous[i] = {}
    for j in coding_codons:
        if codon_to_AA[i] == codon_to_AA[j]:
            is_synonymous[i][j] = True
        else:
            is_synonymous[i][j] = False

pretty_names = {}

def get_pretty_name(key):
    if key in pretty_names:
        return pretty_names[key]
    else:
        return key


def get_nt_to_num_lookup(d):
    D = {order: ord(nt) for nt, order in d.items()}
    lookup = 5 * numpy.ones(max(D.values()) + 1, dtype=numpy.uint8)

    for order, num in D.items():
        lookup[num] = order

    return lookup


def get_codon_to_num_lookup(reverse_complement=False):
    nts = sorted(list(unambiguous_nucleotides))
    as_ints = [ord(nt) for nt in nts]

    size = max(as_ints) + 1
    lookup = 64 * numpy.ones((size, size, size), dtype=numpy.uint8)

    num_to_codon = dict(enumerate(codons))
    if reverse_complement:
        num_to_codon = {k: codon_to_codon_RC[codon] for k, codon in num_to_codon.items()}

    D = {tuple([ord(nt) for nt in codon]): k for k, codon in num_to_codon.items()}

    for a in as_ints:
        for b in as_ints:
            for c in as_ints:
                lookup[a, b, c] = D[(a, b, c)]

    return lookup


# See utils.nt_seq_to_codon_num_array etc. for utilization of these lookup arrays
nt_to_num_lookup = get_nt_to_num_lookup({'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4})
nt_to_RC_num_lookup = get_nt_to_num_lookup({'A': 3, 'C': 2, 'G': 1, 'T': 0, 'N': 4})
codon_to_num_lookup = get_codon_to_num_lookup(reverse_complement=False)
codon_to_RC_num_lookup = get_codon_to_num_lookup(reverse_complement=True)


# anvi'o news stuff
anvio_news_url = "https://raw.githubusercontent.com/merenlab/anvio/master/NEWS.md"
