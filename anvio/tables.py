# -*- coding: utf-8
# pylint: disable=line-too-long
""" Table schemas for databases."""

from anvio.constants import codon_to_AA

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


contigs_db_version = "8"
profile_db_version = "17"
pan_db_version = "4"
samples_info_db_version = "2"
auxiliary_hdf5_db_version = "1"
genomes_storage_hdf5_db_vesion = "3"
users_db_version = "1"

versions_for_db_types = {'contigs': contigs_db_version,
                         'profile': profile_db_version,
                         'pan': pan_db_version}


####################################################################################################
#
#     TABLE DESCRIPTIONS SPECIFIC FOR THE PAN DATABASE (THE REST COMES FROM THE PROFILE DATABASE)
#
####################################################################################################

pan_protein_clusters_table_name        = 'protein_clusters'
pan_protein_clusters_table_structure   = ['entry_id', 'gene_caller_id', 'protein_cluster_id', 'genome_name', 'alignment_summary']
pan_protein_clusters_table_types       = ['numeric' ,     'numeric'   ,         'str'       ,      'str'   ,        'str'       ]


####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE CONTIGS DATABASE
#
####################################################################################################


contig_sequences_table_name            = 'contig_sequences'
contig_sequences_table_structure       = ['contig', 'sequence']
contig_sequences_table_types           = [  'str' ,   'str'   ]

contigs_info_table_name                = 'contigs_basic_info'
contigs_info_table_structure           = ['contig', 'length' , 'gc_content', 'num_splits']
contigs_info_table_types               = [  'str' , 'numeric',   'numeric' ,   'numeric' ]

splits_info_table_name                 = 'splits_basic_info'
splits_info_table_structure            = ['split', 'order_in_parent' , 'start' ,  'end'  , 'length' , 'gc_content', 'gc_content_parent', 'parent' ]
splits_info_table_types                = ['text' ,     'numeric     ','numeric','numeric', 'numeric',   'numeric' ,      'numeric'     ,  'text'  ]


# following tables deal with open reading frames found in contis by a gene caller (such as prodigal), and their functional annotations and stuff.

genes_in_contigs_table_name             = 'genes_in_contigs'
genes_in_contigs_table_structure        = ['gene_callers_id', 'contig', 'start' , 'stop'  , 'direction', 'partial', 'source', 'version']
genes_in_contigs_table_types            = [     'numeric'   ,  'text' ,'numeric','numeric',   'text'   , 'numeric',  'text' ,   'text' ]

genes_in_splits_table_name             = 'genes_in_splits'
genes_in_splits_table_structure        = ['entry_id', 'split', 'gene_callers_id', 'start_in_split', 'stop_in_split', 'percentage_in_split']
genes_in_splits_table_types            = [ 'numeric',  'text',      'numeric'   ,    'numeric'    ,    'numeric'   ,       'numeric'      ]

genes_in_splits_summary_table_name     = 'genes_in_splits_summary'
genes_in_splits_summary_table_structure = ['split', 'num_genes', 'avg_gene_length', 'ratio_coding']
genes_in_splits_summary_table_types     = [ 'text',  'numeric' ,     'numeric'    ,   'numeric'   ]

gene_protein_sequences_table_name      = 'gene_protein_sequences'
gene_protein_sequences_table_structure = ['gene_callers_id', 'sequence']
gene_protein_sequences_table_types     = [     'numeric'   ,   'text'  ]

gene_function_calls_table_name         = 'gene_functions'
gene_function_calls_table_structure    = ['entry_id', 'gene_callers_id', 'source', 'accession', 'function', 'e_value']
gene_function_calls_table_types        = [ 'numeric',     'numeric'    ,  'text' ,    'text'   ,   'text'  , 'numeric']

# tables for taxonomy

taxon_names_table_name                 = 'taxon_names'
taxon_names_table_structure            = ['taxon_id', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxon_names_table_types                = [ 'numeric',   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

splits_taxonomy_table_name             = 'splits_taxonomy'
splits_taxonomy_table_structure        = ['split', 'taxon_id',]
splits_taxonomy_table_types            = [ 'text',  'numeric',]

genes_taxonomy_table_name              = 'genes_taxonomy'
genes_taxonomy_table_structure         = ['gene_callers_id', 'taxon_id',]
genes_taxonomy_table_types             = [    'numeric'    ,  'numeric',]

# the followitn three tables keep hmm hits. they require the gene calls to be made.

hmm_hits_info_table_name               = 'hmm_hits_info'
hmm_hits_info_table_structure          = ['source', 'ref' , 'search_type', 'domain', 'genes']
hmm_hits_info_table_types              = [ 'text' , 'text',    'text'    ,  'text' , 'text' ]       # This one here is the id that apperas in gene_calls table
                                                                                         #         /
hmm_hits_table_name                    = 'hmm_hits'                                      # _______|_______
hmm_hits_table_structure               = ['entry_id', 'source', 'gene_unique_identifier', 'gene_callers_id', 'gene_name', 'gene_hmm_id', 'e_value']
hmm_hits_table_types                   = [ 'numeric',  'text' ,          'text'         ,      'numeric'   ,   'text'   ,     'text'   , 'numeric']

hmm_hits_splits_table_name             = 'hmm_hits_in_splits'
hmm_hits_splits_table_structure        = ['entry_id', 'hmm_hit_entry_id', 'split', 'percentage_in_split', 'source']
hmm_hits_splits_table_types            = [ 'numeric',      'numeric'    ,  'text',       'numeric'      ,  'text' ]


####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE PROFILE DATABASE
#
####################################################################################################

clusterings_table_name               = 'clusterings'
clusterings_table_structure          = ['clustering', 'newick']
clusterings_table_types              = [   'str'    ,  'str'  ]

states_table_name                    = 'states'
states_table_structure               = ['name', 'content', 'last_modified']
states_table_types                   = ['text',  'text'  ,      'text'    ]

variable_aas_table_name              = 'variable_amino_acid_frequencies'
variable_aas_table_structure         = ['entry_id', 'sample_id', 'corresponding_gene_call', 'codon_order_in_gene', 'reference', 'departure_from_reference', 'coverage'] + sorted(list(set(codon_to_AA.values())))
variable_aas_table_types             = [ 'numeric',    'text'  ,        'numeric'         ,       'numeric'      ,    'text'  ,          'numeric'        , 'numeric' ] + ['numeric'] * len(list(set(codon_to_AA.values())))

variable_nts_table_name              = 'variable_nucleotide_positions'
variable_nts_table_structure         = ['entry_id', 'sample_id', 'split_name',   'pos'  , 'pos_in_contig', 'corresponding_gene_call', 'in_partial_gene_call', 'in_complete_gene_call', 'base_pos_in_codon', 'codon_order_in_gene', 'coverage', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'departure_from_reference', 'competing_nts', 'reference',    'A'   ,    'T'   ,    'C'   ,    'G'   ,    'N'   ]
variable_nts_table_types             = [ 'numeric',    'text'  ,    'text'   , 'numeric',    'numeric'   ,        'numeric'         ,       'numeric'       ,       'numeric'        ,       'numeric'    ,       'numeric'      , 'numeric' ,          'bool'       ,          'bool'        ,          'numeric'        ,      'text'    ,    'text'  , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']

gene_coverages_table_name            = 'gene_coverages'
gene_coverages_table_structure       = ['entry_id', 'gene_callers_id', 'sample_id', 'mean_coverage']
gene_coverages_table_types           = [ 'numeric',     'numeric'    ,   'text'   ,    'numeric'   ]

views_table_name                     = 'views'
views_table_structure                = ['view_id', 'target_table']
views_table_types                    = [  'str'  ,      'str'    ]

# notice that atomic data table is the only table that doesn't have a name. because how we use this table is a bit tricky.
# for single profiles, contents of this table is stored as "atomic data", however, for merged profiles,
# each column of the atomic data table becomes its own table, where the row names remain identical, yet columns
# become sample names. 
atomic_data_table_structure          = ['contig', 'std_coverage', 'mean_coverage', 'mean_coverage_Q2Q3', 'max_normalized_ratio', 'relative_abundance', 'detection', 'abundance', 'variability', '__parent__']
atomic_data_table_types              = [ 'text' ,   'numeric'   ,    'numeric'   ,      'numeric'      ,        'numeric'      ,      'numeric'     ,     'numeric'    ,  'numeric' ,   'numeric'  ,    'text'   ]


####################################################################################################
#
#     DESCRIPTIONS FOR TABLES THAT ARE BOTH IN CONTIGS and PROFILE DATABASES
#
####################################################################################################

collections_info_table_name           = 'collections_info'
collections_info_table_structure      = ['collection_name', 'num_splits', 'num_bins', 'bin_names']
collections_info_table_types          = [      'text'     ,  'numeric'  ,  'numeric',    'text'  ]

collections_bins_info_table_name      = 'collections_bins_info'
collections_bins_info_table_structure = ['entry_id', 'collection_name', 'bin_name', 'source', 'html_color']
collections_bins_info_table_types     = [ 'numeric',       'text'     ,   'text'  ,  'text' ,    'text'   ]

collections_contigs_table_name        = 'collections_of_contigs'
collections_contigs_table_structure   = ['entry_id', 'collection_name', 'contig', 'bin_name']
collections_contigs_table_types       = [ 'numeric',       'text'     ,  'text' ,   'text'  ]

collections_splits_table_name         = 'collections_of_splits'
collections_splits_table_structure    = ['entry_id', 'collection_name', 'split', 'bin_name']
collections_splits_table_types        = [ 'numeric',       'text'     , 'text' ,   'text'  ]


####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE SAMPLES INFORMATION DATABASE
#
####################################################################################################

# samples information table is where what we know about samples will be stored. different experiments
# will have different attributes (i.e., gestational age will be an attribute for samples in a human gut
# dataset, and in contrast pH or temperature will be more relevant to a marine sample). these info will
# be generated from the input files, and the columns of samples_information table will be determined
# on the fly (see anvio/dbops/SamplesInfoDatabase/create for details)
samples_information_table_name             = 'samples_information'

samples_order_table_name                   = 'samples_order'
samples_order_table_structure              = ['attributes', 'basic', 'newick']
samples_order_table_types                  = [    'str'   ,  'str' ,  'str' ]

# this table will map user defined attributes to ascii-only, reliable aliases. those aliases will be used
# in the sampels_information table, and every time that information is read, the dictionary should be
# generated by converting aliases into real attributes defined in this table. easy peasy.
samples_attribute_aliases_table_name       = 'samples_attribute_aliases'
samples_attribute_aliases_table_structure  = ['alias', 'attribute']
samples_attribute_aliases_table_types      = [ 'str' ,     'str'  ]

