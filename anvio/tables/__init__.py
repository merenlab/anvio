# -*- coding: utf-8
# pylint: disable=line-too-long
""" Table schemas for databases."""

from anvio.constants import codons, nucleotides, essential_genome_info


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


contigs_db_version = "12"
profile_db_version = "30"
pan_db_version = "12"
auxiliary_data_version = "2"
structure_db_version = "1"
genomes_storage_vesion = "6"

versions_for_db_types = {'contigs': contigs_db_version,
                         'profile': profile_db_version,
                         'structure': structure_db_version,
                         'pan': pan_db_version,
                         'genomestorage': genomes_storage_vesion,
                         'auxiliary data for coverages': auxiliary_data_version}

####################################################################################################
#
#     TABLE DESCRIPTIONS SPECIFIC FOR THE PAN DATABASE (THE REST COMES FROM THE PROFILE DATABASE)
#
####################################################################################################

pan_gene_clusters_table_name           = 'gene_clusters'
pan_gene_clusters_table_structure      = ['entry_id', 'gene_caller_id', 'gene_cluster_id', 'genome_name', 'alignment_summary']
pan_gene_clusters_table_types          = ['numeric' ,     'numeric'   ,      'str'       ,     'str'    ,        'str'       ]


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
genes_in_contigs_table_types            = [    'numeric'    ,  'text' ,'numeric','numeric',   'text'   , 'numeric',  'text' ,   'text' ]

genes_in_splits_table_name             = 'genes_in_splits'
genes_in_splits_table_structure        = ['entry_id', 'split', 'gene_callers_id', 'start_in_split', 'stop_in_split', 'percentage_in_split']
genes_in_splits_table_types            = [ 'numeric',  'text',      'numeric'   ,    'numeric'    ,    'numeric'   ,       'numeric'      ]

gene_amino_acid_sequences_table_name      = 'gene_amino_acid_sequences'
gene_amino_acid_sequences_table_structure = ['gene_callers_id', 'sequence']
gene_amino_acid_sequences_table_types     = [     'numeric'   ,   'text'  ]

gene_function_calls_table_name         = 'gene_functions'
gene_function_calls_table_structure    = ['entry_id', 'gene_callers_id', 'source', 'accession', 'function', 'e_value']
gene_function_calls_table_types        = [ 'numeric',     'numeric'    ,  'text' ,    'text'  ,   'text'  , 'numeric']

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
hmm_hits_info_table_types              = [ 'text' , 'text',    'text'    ,  'text' , 'text' ]       # This one here is the id that apper as in gene_calls table
                                                                                         #         /
hmm_hits_table_name                    = 'hmm_hits'                                      # _______|_______
hmm_hits_table_structure               = ['entry_id', 'source', 'gene_unique_identifier', 'gene_callers_id', 'gene_name', 'gene_hmm_id', 'e_value']
hmm_hits_table_types                   = [ 'numeric',  'text' ,          'text'         ,      'numeric'   ,   'text'   ,     'text'   , 'numeric']

hmm_hits_splits_table_name             = 'hmm_hits_in_splits'
hmm_hits_splits_table_structure        = ['entry_id', 'hmm_hit_entry_id', 'split', 'percentage_in_split', 'source']
hmm_hits_splits_table_types            = [ 'numeric',      'numeric'    ,  'text',       'numeric'      ,  'text' ]

# following table keeps nt poisition info

nt_position_info_table_name       = 'nt_position_info'
nt_position_info_table_structure  = ['contig_name', 'position_info']
nt_position_info_table_types      = [    'str'    ,      'blob'    ]


####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE PROFILE DATABASE
#
####################################################################################################

item_orders_table_name               = 'item_orders'
item_orders_table_structure         = ['name', 'type', 'data', 'additional']
item_orders_table_types             = ['text', 'text', 'text',    'text'   ]

item_additional_data_table_name      = 'item_additional_data'
item_additional_data_table_structure = ['entry_id', 'item_name', 'data_key', 'data_value', 'data_type', 'data_group']
item_additional_data_table_types     = [ 'numeric',    'text'  ,   'text'  ,    'text'   ,    'text'  ,    'text'   ]

layer_orders_table_name              = 'layer_orders'
layer_orders_table_structure         = ['data_key', 'data_type', 'data_value']
layer_orders_table_types             = [  'text'  ,    'text'  ,    'text'   ]

layer_additional_data_table_name      = 'layer_additional_data'
layer_additional_data_table_structure = ['entry_id', 'item_name', 'data_key', 'data_value', 'data_type', 'data_group']
layer_additional_data_table_types     = [ 'numeric',    'text'  ,   'text'  ,    'text'   ,    'text'  ,    'text'   ]

states_table_name                    = 'states'
states_table_structure               = ['name', 'content', 'last_modified']
states_table_types                   = ['text',  'text'  ,      'text'    ]

variable_codons_table_name           = 'variable_codons'
variable_codons_table_structure      = ['entry_id', 'sample_id', 'corresponding_gene_call', 'codon_order_in_gene', 'reference', 'departure_from_reference', 'coverage'] + codons
variable_codons_table_types          = [ 'numeric',    'text'  ,        'numeric'         ,       'numeric'      ,    'text'  ,          'numeric'        , 'numeric' ] + ['numeric'] * len(codons)

variable_nts_table_name              = 'variable_nucleotides'
variable_nts_table_structure         = ['entry_id', 'sample_id', 'split_name',   'pos'  , 'pos_in_contig', 'corresponding_gene_call', 'in_partial_gene_call', 'in_complete_gene_call', 'base_pos_in_codon', 'codon_order_in_gene', 'coverage', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'departure_from_reference', 'competing_nts', 'reference'] + nucleotides
variable_nts_table_types             = [ 'numeric',    'text'  ,    'text'   , 'numeric',    'numeric'   ,        'numeric'         ,       'numeric'       ,       'numeric'        ,       'numeric'    ,       'numeric'      , 'numeric' ,          'bool'       ,          'bool'        ,          'numeric'        ,      'text'    ,    'text'  ] + ['numeric'] * len(nucleotides)

views_table_name                     = 'views'
views_table_structure                = ['view_id', 'target_table']
views_table_types                    = [  'str'  ,      'str'    ]

# notice that atomic data table is the only table that doesn't have a name. because how we use this table is a bit tricky.
# for single profiles, contents of this table is stored as "atomic data", however, for merged profiles,
# each column of the atomic data table becomes its own table, where the row names remain identical, yet columns
# become sample names.
atomic_data_table_structure          = ['contig', 'std_coverage', 'mean_coverage', 'mean_coverage_Q2Q3', 'max_normalized_ratio', 'relative_abundance', 'detection', 'abundance', 'variability', '__parent__']
atomic_data_table_types              = [ 'text' ,   'numeric'   ,    'numeric'   ,      'numeric'      ,        'numeric'      ,      'numeric'     ,   'numeric' ,  'numeric' ,   'numeric'  ,    'text'   ]


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
#     TABLE DESCRIPTIONS FOR THE PROFILE AUXILIARY COVERAGE DATABASE
#
####################################################################################################

split_coverages_table_name       = 'split_coverages'
split_coverages_table_structure  = ['split_name', 'sample_name', 'coverages']
split_coverages_table_types      = [    'str'   ,     'str'    ,   'blob'   ]


####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE GENOME STORAGE
#
####################################################################################################

genome_info_table_name       = 'genome_info'
genome_info_table_structure  = ['genome_name', 'genome_hash', 'external_genome'] + essential_genome_info
genome_info_table_types      = [    'str'    ,     'text'   ,     'numeric'    ] + ['numeric'] * len(essential_genome_info)

gene_info_table_name       = 'gene_info'
gene_info_table_structure  = ['genome_name', 'gene_caller_id', 'aa_sequence', 'dna_sequence', 'partial', 'length' ]
gene_info_table_types      = [    'str'    ,     'numeric'   ,    'text'    ,     'text'    , 'numeric', 'numeric']

genome_gene_function_calls_table_name      = 'gene_function_calls'
genome_gene_function_calls_table_structure = ['genome_name', ] + gene_function_calls_table_structure[:]
genome_gene_function_calls_table_types     = [    'str'    , ] + gene_function_calls_table_types[:]

tables_without_unique_entry_ids = [genome_gene_function_calls_table_name]

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE STRUCTURE DB
#
####################################################################################################

structure_pdb_data_table_name       = 'structures'
structure_pdb_data_table_structure  = ['corresponding_gene_call', 'pdb_content']
structure_pdb_data_table_types      = [         'integer'       ,    'blob'    ]

structure_templates_table_name       = 'templates'
structure_templates_table_structure  = ['entry_id' , 'corresponding_gene_call' , 'pdb_id' , 'chain_id' , 'ppi']
structure_templates_table_types      = ['integer'  , 'integer'                 , 'text'   , 'text'     , 'real']

structure_models_table_name       = 'models'
structure_models_table_structure  = ['entry_id' , 'corresponding_gene_call' , 'molpdf' , 'GA341_score' , 'DOPE_score' , 'picked_as_best']
structure_models_table_types      = ['integer'  , 'integer'                 , 'real'   , 'real'        , 'real'       , 'integer']

# The FULL table structure is defined in the StructureDatabase class based on what annotation
# sources are found and/or requested.
structure_residue_info_table_name       = 'residue_info'
structure_residue_info_table_structure  = ['entry_id', 'corresponding_gene_call', 'codon_order_in_gene', 'contact_numbers', 'codon', 'amino_acid', 'codon_number']
structure_residue_info_table_types      = ['integer',         'integer'        ,        'integer'     ,   'text'          , 'text',  'text',       'integer']

residue_info_sources = {"DSSP":        {"structure": ['codon_order_in_gene' , 'aa'   , 'sec_struct' , 'rel_solvent_acc' , 'phi'  , 'psi'  , 'NH_O_1_index' , 'NH_O_1_energy' , 'O_NH_1_index' , 'O_NH_1_energy' , 'NH_O_2_index' , 'NH_O_2_energy' , 'O_NH_2_index' , 'O_NH_2_energy'],
                                        "types":     ['integer'             , 'text' , 'text'       , 'real'            , 'real' , 'real' , 'integer'      , 'real'          , 'integer'      , 'real'          , 'integer'      , 'real'          , 'integer'      , 'real']},
                       }
