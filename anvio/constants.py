# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

import os
import sys
import glob
import numpy
import string

from collections import Counter

import anvio

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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

default_prostt5_weight_path = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PROSTT5/weights')
choice_of_pangenome = ["structure", "sequence"]

clustering_configs_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/clusterconfigs')
clustering_configs = {}

default_scgs_taxonomy_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY')
default_scgs_for_taxonomy = ['Ribosomal_L1',
                             'Ribosomal_L13',
                             'Ribosomal_L14',
                             'Ribosomal_L16',
                             'Ribosomal_L17',
                             'Ribosomal_L19',
                             'Ribosomal_L2',
                             'Ribosomal_L20',
                             'Ribosomal_L21p',
                             'Ribosomal_L22',
                             'Ribosomal_L27A',
                             'Ribosomal_L3',
                             'Ribosomal_L4',
                             'Ribosomal_L5',
                             'Ribosomal_S11',
                             'Ribosomal_S15',
                             'Ribosomal_S16',
                             'Ribosomal_S2',
                             'Ribosomal_S6',
                             'Ribosomal_S7',
                             'Ribosomal_S8',
                             'Ribosomal_S9']

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

# The following block of constants is used in the tRNA-seq workflow.
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
# (This does not apply to the tRNA-seq workflow, which stores coverages as uint32.)
max_depth_for_coverage = 60000

# default methods for hierarchical cluster analyses
distance_metric_default = 'euclidean'
linkage_method_default = 'ward'

# The purpose of the `fetch_filters` dictionary below is to filter reads as they are
# read from BAM files especially during anvi'o profiling (the primary client of this
# dictionary is `anvio/bamops.py`). Essentially, any combination of the following
# properties defined in the `read` object returned by the `fetch` function of pysam
# can be added to this dictionary to create new filters that are then globally applied
# to 'fetched' reads during profiling to exclude those that return `false`:
#
#     >>> 'aend', 'alen', 'aligned_pairs', 'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples',
#         'compare', 'flag', 'from_dict', 'fromstring', 'get_aligned_pairs', 'get_blocks', 'get_cigar_stats',
#         'get_forward_qualities', 'get_forward_sequence', 'get_overlap', 'get_reference_positions',
#         'get_reference_sequence', 'get_tag', 'get_tags', 'has_tag', 'header', 'infer_query_length',
#         'infer_read_length', 'inferred_length', 'is_duplicate', 'is_paired', 'is_proper_pair',
#         'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary',
#         'is_unmapped', 'isize', 'mapping_quality', 'mapq', 'mate_is_reverse', 'mate_is_unmapped', 'mpos',
#         'mrnm', 'next_reference_id', 'next_reference_name', 'next_reference_start', 'opt', 'overlap', 'pnext',
#         'pos', 'positions', 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end',
#         'query_alignment_length', 'query_alignment_qualities', 'query_alignment_sequence',
#         'query_alignment_start', 'query_length', 'query_name', 'query_qualities', 'query_sequence',
#         'reference_end', 'reference_id', 'reference_length', 'reference_name', 'reference_start', 'rlen',
#         'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags', 'template_length', 'tid', 'tlen',
#         'to_dict', 'to_string', 'tostring'
#
# Please note that these variable names may change across versions of pysam. See anvio/bamops.py for most
# up-to-date usage of these filters since we are terrible at updating comments elsewhere in the code after
# making significant changes to our modules :/
fetch_filters = {None                 : None,
                 'proper-pairs'       : lambda x: not x.mate_is_unmapped,
                 'double-forwards'    : lambda x: x.is_paired and not x.is_reverse and not x.mate_is_reverse and not x.mate_is_unmapped and x.reference_name == x.next_reference_name,
                 'double-reverses'    : lambda x: x.is_paired and x.is_reverse and x.mate_is_reverse and not x.mate_is_unmapped and x.reference_name == x.next_reference_name,
                 'inversions'         : lambda x: (x.is_paired and not x.is_reverse and not x.mate_is_reverse and not x.mate_is_unmapped and x.reference_name == x.next_reference_name and (abs(x.tlen) < 2000)) or \
                                                  (x.is_paired and x.is_reverse and x.mate_is_reverse and not x.mate_is_unmapped and x.reference_name == x.next_reference_name and (abs(x.tlen) < 2000)),
                 'distant-inversions'  : lambda x: (x.is_paired and not x.is_reverse and not x.mate_is_reverse and not x.mate_is_unmapped and x.reference_name == x.next_reference_name) or \
                                                  (x.is_paired and x.is_reverse and x.mate_is_reverse and not x.mate_is_unmapped and x.reference_name == x.next_reference_name),
                 'single-mapped-reads': lambda x: x.mate_is_unmapped,
                 'distant-pairs-1K'   : lambda x: x.is_paired and not x.mate_is_unmapped and abs(x.tlen) > 1000}

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

AA_to_anticodons = Counter({'Ala': ['AGC', 'CGC', 'GGC', 'TGC'],
                            'Arg': ['ACG', 'CCG', 'CCT', 'GCG', 'TCG', 'TCT'],
                            'Asn': ['ATT', 'GTT'],
                            'Asp': ['ATC', 'GTC'],
                            'Cys': ['ACA', 'GCA'],
                            'Gln': ['CTG', 'TTG'],
                            'Glu': ['CTC', 'TTC'],
                            'Gly': ['ACC', 'CCC', 'GCC', 'TCC'],
                            'His': ['ATG', 'GTG'],
                            'Ile': ['AAT', 'GAT', 'TAT'],
                            'Leu': ['AAG', 'CAA', 'CAG', 'GAG', 'TAA', 'TAG'],
                            'Lys': ['CTT', 'TTT'],
                            'Met': ['CAT'],
                            'Phe': ['AAA', 'GAA'],
                            'Pro': ['AGG', 'CGG', 'GGG', 'TGG'],
                            'STP': ['CTA', 'TCA', 'TTA'],
                            'Ser': ['ACT', 'AGA', 'CGA', 'GCT', 'GGA', 'TGA'],
                            'Thr': ['AGT', 'CGT', 'GGT', 'TGT'],
                            'Trp': ['CCA'],
                            'Tyr': ['ATA', 'GTA'],
                            'Val': ['AAC', 'CAC', 'GAC', 'TAC']})

AA_to_single_letter_code = Counter({'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
                                    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
                                    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                                    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'STP': '*',
                                    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y',
                                    'Val': 'V'})

amino_acids = sorted(list(AA_to_single_letter_code.keys()))

# Standard genetic code (translation table 1 at the following link)
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
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


# The KEGG KOfams of multidomain ribosomal proteins are used as a reference
# in the calculation of codon usage bias metrics by anvi-get-codon-frequencies.
# This list is taken from the BRITE hierarchy for the ribosome.
# 09180 Brite Hierarchies ->
#  09182 Protein families: genetic information processing ->
#   03011 Ribosome [BR:ko03011]
# The list can be seen here: https://www.genome.jp/brite/ko00001+K02981
ribosomal_protein_kofams = [
    "K02981", # RP-S2e, RPS2; small subunit ribosomal protein S2e
    "K02985", # RP-S3e, RPS3; small subunit ribosomal protein S3e
    "K02984", # RP-S3Ae, RPS3A; small subunit ribosomal protein S3Ae
    "K02987", # RP-S4e, RPS4; small subunit ribosomal protein S4e
    "K02989", # RP-S5e, RPS5; small subunit ribosomal protein S5e
    "K02991", # RP-S6e, RPS6; small subunit ribosomal protein S6e
    "K02993", # RP-S7e, RPS7; small subunit ribosomal protein S7e
    "K02995", # RP-S8e, RPS8; small subunit ribosomal protein S8e
    "K02997", # RP-S9e, RPS9; small subunit ribosomal protein S9e
    "K02947", # RP-S10e, RPS10; small subunit ribosomal protein S10e
    "K02949", # RP-S11e, RPS11; small subunit ribosomal protein S11e
    "K02951", # RP-S12e, RPS12; small subunit ribosomal protein S12e
    "K02953", # RP-S13e, RPS13; small subunit ribosomal protein S13e
    "K02955", # RP-S14e, RPS14; small subunit ribosomal protein S14e
    "K02958", # RP-S15e, RPS15; small subunit ribosomal protein S15e
    "K02957", # RP-S15Ae, RPS15A; small subunit ribosomal protein S15Ae
    "K02960", # RP-S16e, RPS16; small subunit ribosomal protein S16e
    "K02962", # RP-S17e, RPS17; small subunit ribosomal protein S17e
    "K02964", # RP-S18e, RPS18; small subunit ribosomal protein S18e
    "K02966", # RP-S19e, RPS19; small subunit ribosomal protein S19e
    "K02969", # RP-S20e, RPS20; small subunit ribosomal protein S20e
    "K02971", # RP-S21e, RPS21; small subunit ribosomal protein S21e
    "K02973", # RP-S23e, RPS23; small subunit ribosomal protein S23e
    "K02974", # RP-S24e, RPS24; small subunit ribosomal protein S24e
    "K02975", # RP-S25e, RPS25; small subunit ribosomal protein S25e
    "K02976", # RP-S26e, RPS26; small subunit ribosomal protein S26e
    "K02978", # RP-S27e, RPS27; small subunit ribosomal protein S27e
    "K02977", # RP-S27Ae, RPS27A, UBA80; ubiquitin-small subunit ribosomal protein S27Ae
    "K02979", # RP-S28e, RPS28; small subunit ribosomal protein S28e
    "K02980", # RP-S29e, RPS29; small subunit ribosomal protein S29e
    "K02983", # RP-S30e, RPS30; small subunit ribosomal protein S30e
    "K02998", # RP-SAe, RPSA; small subunit ribosomal protein SAe
    "K02925", # RP-L3e, RPL3; large subunit ribosomal protein L3e
    "K02930", # RP-L4e, RPL4; large subunit ribosomal protein L4e
    "K02932", # RP-L5e, RPL5; large subunit ribosomal protein L5e
    "K02934", # RP-L6e, RPL6; large subunit ribosomal protein L6e
    "K02937", # RP-L7e, RPL7; large subunit ribosomal protein L7e
    "K02936", # RP-L7Ae, RPL7A; large subunit ribosomal protein L7Ae
    "K02938", # RP-L8e, RPL8; large subunit ribosomal protein L8e
    "K02940", # RP-L9e, RPL9; large subunit ribosomal protein L9e
    "K02866", # RP-L10e, RPL10; large subunit ribosomal protein L10e
    "K02865", # RP-L10Ae, RPL10A; large subunit ribosomal protein L10Ae
    "K02868", # RP-L11e, RPL11; large subunit ribosomal protein L11e
    "K02870", # RP-L12e, RPL12; large subunit ribosomal protein L12e
    "K02873", # RP-L13e, RPL13; large subunit ribosomal protein L13e
    "K02872", # RP-L13Ae, RPL13A; large subunit ribosomal protein L13Ae
    "K02875", # RP-L14e, RPL14; large subunit ribosomal protein L14e
    "K02877", # RP-L15e, RPL15; large subunit ribosomal protein L15e
    "K02880", # RP-L17e, RPL17; large subunit ribosomal protein L17e
    "K02883", # RP-L18e, RPL18; large subunit ribosomal protein L18e
    "K02882", # RP-L18Ae, RPL18A; large subunit ribosomal protein L18Ae
    "K02885", # RP-L19e, RPL19; large subunit ribosomal protein L19e
    "K02889", # RP-L21e, RPL21; large subunit ribosomal protein L21e
    "K02891", # RP-L22e, RPL22; large subunit ribosomal protein L22e
    "K02894", # RP-L23e, RPL23; large subunit ribosomal protein L23e
    "K02893", # RP-L23Ae, RPL23A; large subunit ribosomal protein L23Ae
    "K02896", # RP-L24e, RPL24; large subunit ribosomal protein L24e
    "K02898", # RP-L26e, RPL26; large subunit ribosomal protein L26e
    "K02901", # RP-L27e, RPL27; large subunit ribosomal protein L27e
    "K02900", # RP-L27Ae, RPL27A; large subunit ribosomal protein L27Ae
    "K02903", # RP-L28e, RPL28; large subunit ribosomal protein L28e
    "K02905", # RP-L29e, RPL29; large subunit ribosomal protein L29e
    "K02908", # RP-L30e, RPL30; large subunit ribosomal protein L30e
    "K02910", # RP-L31e, RPL31; large subunit ribosomal protein L31e
    "K02912", # RP-L32e, RPL32; large subunit ribosomal protein L32e
    "K02915", # RP-L34e, RPL34; large subunit ribosomal protein L34e
    "K02918", # RP-L35e, RPL35; large subunit ribosomal protein L35e
    "K02917", # RP-L35Ae, RPL35A; large subunit ribosomal protein L35Ae
    "K02920", # RP-L36e, RPL36; large subunit ribosomal protein L36e
    "K02922", # RP-L37e, RPL37; large subunit ribosomal protein L37e
    "K02921", # RP-L37Ae, RPL37A; large subunit ribosomal protein L37Ae
    "K02923", # RP-L38e, RPL38; large subunit ribosomal protein L38e
    "K02924", # RP-L39e, RPL39; large subunit ribosomal protein L39e
    "K02927", # RP-L40e, RPL40, UBA52; ubiquitin-large subunit ribosomal protein L40e
    "K02928", # RP-L41e, RPL41; large subunit ribosomal protein L41e
    "K02929", # RP-L44e, RPL44; large subunit ribosomal protein L44e
    "K02941", # RP-LP0, RPLP0; large subunit ribosomal protein LP0
    "K02942", # RP-LP1, RPLP1; large subunit ribosomal protein LP1
    "K02943", # RP-LP2, RPLP2; large subunit ribosomal protein LP2
    "K02967", # RP-S2, MRPS2, rpsB; small subunit ribosomal protein S2
    "K02988", # RP-S5, MRPS5, rpsE; small subunit ribosomal protein S5
    "K02990", # RP-S6, MRPS6, rpsF; small subunit ribosomal protein S6
    "K02992", # RP-S7, MRPS7, rpsG; small subunit ribosomal protein S7
    "K02996", # RP-S9, MRPS9, rpsI; small subunit ribosomal protein S9
    "K02946", # RP-S10, MRPS10, rpsJ; small subunit ribosomal protein S10
    "K02948", # RP-S11, MRPS11, rpsK; small subunit ribosomal protein S11
    "K02950", # RP-S12, MRPS12, rpsL; small subunit ribosomal protein S12
    "K02954", # RP-S14, MRPS14, rpsN; small subunit ribosomal protein S14
    "K02956", # RP-S15, MRPS15, rpsO; small subunit ribosomal protein S15
    "K02959", # RP-S16, MRPS16, rpsP; small subunit ribosomal protein S16
    "K02961", # RP-S17, MRPS17, rpsQ; small subunit ribosomal protein S17
    "K02963", # RP-S18, MRPS18, rpsR; small subunit ribosomal protein S18
    "K02970", # RP-S21, MRPS21, rpsU; small subunit ribosomal protein S21
    "K16174", # MRPS18B, MRPS18-2; small subunit ribosomal protein S18b, mitochondrial
    "K17401", # MRPS22; small subunit ribosomal protein S22
    "K17402", # MRPS23; small subunit ribosomal protein S23
    "K17403", # MRPS24; small subunit ribosomal protein S24
    "K17404", # MRPS25; small subunit ribosomal protein S25
    "K17405", # MRPS26; small subunit ribosomal protein S26
    "K17406", # MRPS27; small subunit ribosomal protein S27
    "K17407", # MRPS28; small subunit ribosomal protein S28
    "K17408", # DAP3, MRPS29; small subunit ribosomal protein S29
    "K17409", # MRPS30; small subunit ribosomal protein S30
    "K17410", # MRPS31; small subunit ribosomal protein S31
    "K17411", # MRPS33; small subunit ribosomal protein S33
    "K17412", # MRPS34; small subunit ribosomal protein S34
    "K17413", # MRPS35; small subunit ribosomal protein S35
    "K17414", # MRPS36; small subunit ribosomal protein S36
    "K17415", # MRP21; small subunit ribosomal protein MRP21
    "K17417", # YMR31; small subunit ribosomal protein YMR-31
    "K19032", # PSRP3; 30S ribosomal protein 3
    "K19033", # PSRP4, RPS31; 30S ribosomal protein S31
    "K02863", # RP-L1, MRPL1, rplA; large subunit ribosomal protein L1
    "K02886", # RP-L2, MRPL2, rplB; large subunit ribosomal protein L2
    "K02906", # RP-L3, MRPL3, rplC; large subunit ribosomal protein L3
    "K02926", # RP-L4, MRPL4, rplD; large subunit ribosomal protein L4
    "K02931", # RP-L5, MRPL5, rplE; large subunit ribosomal protein L5
    "K02933", # RP-L6, MRPL6, rplF; large subunit ribosomal protein L6
    "K02939", # RP-L9, MRPL9, rplI; large subunit ribosomal protein L9
    "K02864", # RP-L10, MRPL10, rplJ; large subunit ribosomal protein L10
    "K02867", # RP-L11, MRPL11, rplK; large subunit ribosomal protein L11
    "K02935", # RP-L7, MRPL12, rplL; large subunit ribosomal protein L7/L12
    "K02871", # RP-L13, MRPL13, rplM; large subunit ribosomal protein L13
    "K02874", # RP-L14, MRPL14, rplN; large subunit ribosomal protein L14
    "K02876", # RP-L15, MRPL15, rplO; large subunit ribosomal protein L15
    "K02878", # RP-L16, MRPL16, rplP; large subunit ribosomal protein L16
    "K02879", # RP-L17, MRPL17, rplQ; large subunit ribosomal protein L17
    "K02881", # RP-L18, MRPL18, rplR; large subunit ribosomal protein L18
    "K02884", # RP-L19, MRPL19, rplS; large subunit ribosomal protein L19
    "K02887", # RP-L20, MRPL20, rplT; large subunit ribosomal protein L20
    "K02888", # RP-L21, MRPL21, rplU; large subunit ribosomal protein L21
    "K02890", # RP-L22, MRPL22, rplV; large subunit ribosomal protein L22
    "K02892", # RP-L23, MRPL23, rplW; large subunit ribosomal protein L23
    "K02895", # RP-L24, MRPL24, rplX; large subunit ribosomal protein L24
    "K02899", # RP-L27, MRPL27, rpmA; large subunit ribosomal protein L27
    "K02902", # RP-L28, MRPL28, rpmB; large subunit ribosomal protein L28
    "K02907", # RP-L30, MRPL30, rpmD; large subunit ribosomal protein L30
    "K02911", # RP-L32, MRPL32, rpmF; large subunit ribosomal protein L32
    "K02913", # RP-L33, MRPL33, rpmG; large subunit ribosomal protein L33
    "K02914", # RP-L34, MRPL34, rpmH; large subunit ribosomal protein L34
    "K02916", # RP-L35, MRPL35, rpmI; large subunit ribosomal protein L35
    "K02919", # RP-L36, MRPL36, rpmJ; large subunit ribosomal protein L36
    "K17418", # MRPL37; large subunit ribosomal protein L37
    "K17419", # MRPL38; large subunit ribosomal protein L38
    "K17420", # MRPL39; large subunit ribosomal protein L39
    "K17421", # MRPL40; large subunit ribosomal protein L40
    "K17422", # MRPL41; large subunit ribosomal protein L41
    "K17423", # MRPL42; large subunit ribosomal protein L42
    "K17424", # MRPL43; large subunit ribosomal protein L43
    "K17425", # MRPL44; large subunit ribosomal protein L44 [EC:3.1.26.-]
    "K17426", # MRPL45; large subunit ribosomal protein L45
    "K17427", # MRPL46; large subunit ribosomal protein L46
    "K17428", # MRPL47, NCM1; large subunit ribosomal protein L47
    "K17429", # MRPL48; large subunit ribosomal protein L48
    "K17430", # MRPL49, NOF1; large subunit ribosomal protein L49
    "K17431", # MRPL50; large subunit ribosomal protein L50
    "K17432", # MRPL51; large subunit ribosomal protein L51
    "K17433", # MRPL52; large subunit ribosomal protein L52
    "K17434", # MRPL53; large subunit ribosomal protein L53
    "K17435", # MRPL54; large subunit ribosomal protein L54
    "K17436", # MRPL55; large subunit ribosomal protein L55
    "K17437", # MRPL15; large subunit ribosomal protein L15
    "K17438", # MRPL25; large subunit ribosomal protein L25
    "K17439", # MRPL35; large subunit ribosomal protein L35
    "K17440", # MRP49; large subunit ribosomal protein MRP49
    "K19034", # PSRP5; 50S ribosomal protein 5
    "K19035", # PSRP6; 50S ribosomal protein 6
    "K02945", # RP-S1, rpsA; small subunit ribosomal protein S1
    "K02982", # RP-S3, rpsC; small subunit ribosomal protein S3
    "K02986", # RP-S4, rpsD; small subunit ribosomal protein S4
    "K02994", # RP-S8, rpsH; small subunit ribosomal protein S8
    "K02952", # RP-S13, rpsM; small subunit ribosomal protein S13
    "K02965", # RP-S19, rpsS; small subunit ribosomal protein S19
    "K02968", # RP-S20, rpsT; small subunit ribosomal protein S20
    "K02897", # RP-L25, rplY; large subunit ribosomal protein L25
    "K02904", # RP-L29, rpmC; large subunit ribosomal protein L29
    "K02909", # RP-L31, rpmE; large subunit ribosomal protein L31
    "K07590", # RP-L7A, rplGB; large subunit ribosomal protein L7A
    "K02869", # RP-L12, rpl12; large subunit ribosomal protein L12
    "K02944", # RP-LX, rplX; large subunit ribosomal protein LX
]


# anvi'o news stuff
anvio_news_url = "https://raw.githubusercontent.com/merenlab/anvio/master/NEWS.md"
