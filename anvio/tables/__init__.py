# -*- coding: utf-8
# pylint: disable=line-too-long
""" Table schemas for databases."""

from anvio.constants import codons, nucleotides, essential_genome_info, TRNA_FEATURE_NAMES, THREEPRIME_VARIANTS

import itertools

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


contigs_db_version = "19"
profile_db_version = "35"
genes_db_version = "6"
pan_db_version = "14"
auxiliary_data_version = "2"
structure_db_version = "2"
genomes_storage_vesion = "7"
trnaseq_db_version = "1"
workflow_config_version = "1"
metabolic_modules_db_version = "2"

versions_for_db_types = {'contigs': contigs_db_version,
                         'profile': profile_db_version,
                         'genes': genes_db_version,
                         'structure': structure_db_version,
                         'pan': pan_db_version,
                         'genomestorage': genomes_storage_vesion,
                         'auxiliary data for coverages': auxiliary_data_version,
                         'trnaseq': trnaseq_db_version,
                         'config': workflow_config_version,
                         'modules': metabolic_modules_db_version}


####################################################################################################
#
#     TABLE DESCRIPTIONS SPECIFIC FOR THE PAN DATABASE (THE REST COMES FROM THE PROFILE DATABASE)
#
####################################################################################################

pan_gene_clusters_table_name           = 'gene_clusters'
pan_gene_clusters_table_structure      = ['gene_caller_id', 'gene_cluster_id', 'genome_name', 'alignment_summary']
pan_gene_clusters_table_types          = [    'numeric'   ,      'str'       ,     'str'    ,        'str'       ]


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


# following tables deal with open reading frames found in contigs by a gene caller (such as prodigal), and their functional annotations and stuff.

genes_in_contigs_table_name             = 'genes_in_contigs'
genes_in_contigs_table_structure        = ['gene_callers_id', 'contig', 'start' , 'stop'  , 'direction', 'partial', 'call_type', 'source', 'version']
genes_in_contigs_table_types            = [    'numeric'    ,  'text' ,'numeric','numeric',   'text'   , 'numeric',  'numeric' ,  'text' ,   'text' ]

genes_in_splits_table_name             = 'genes_in_splits'
genes_in_splits_table_structure        = ['split', 'gene_callers_id', 'start_in_split', 'stop_in_split', 'percentage_in_split']
genes_in_splits_table_types            = [ 'text',      'numeric'   ,    'numeric'    ,    'numeric'   ,       'numeric'      ]

gene_amino_acid_sequences_table_name      = 'gene_amino_acid_sequences'
gene_amino_acid_sequences_table_structure = ['gene_callers_id', 'sequence']
gene_amino_acid_sequences_table_types     = [     'numeric'   ,   'text'  ]

gene_function_calls_table_name         = 'gene_functions'
gene_function_calls_table_structure    = ['gene_callers_id', 'source', 'accession', 'function', 'e_value']
gene_function_calls_table_types        = [    'numeric'    ,  'text' ,    'text'  ,   'text'  , 'numeric']

taxon_names_table_name                 = 'taxon_names'
taxon_names_table_structure            = ['taxon_id', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxon_names_table_types                = [ 'numeric',   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

splits_taxonomy_table_name             = 'splits_taxonomy'
splits_taxonomy_table_structure        = ['split', 'taxon_id',]
splits_taxonomy_table_types            = [ 'text',  'numeric',]

genes_taxonomy_table_name              = 'genes_taxonomy'
genes_taxonomy_table_structure         = ['gene_callers_id', 'taxon_id',]
genes_taxonomy_table_types             = [    'numeric'    ,  'numeric',]

hmm_hits_info_table_name               = 'hmm_hits_info'
hmm_hits_info_table_structure          = ['source', 'ref' , 'search_type', 'domain', 'genes']
hmm_hits_info_table_types              = [ 'text' , 'text',    'text'    ,  'text' , 'text' ]       # This one here is the id that apper as in gene_calls table
                                                                                         #         /
hmm_hits_table_name                    = 'hmm_hits'                                      # _______|_______
hmm_hits_table_structure               = ['entry_id', 'source', 'gene_unique_identifier', 'gene_callers_id', 'gene_name', 'gene_hmm_id', 'e_value']
hmm_hits_table_types                   = [ 'numeric',  'text' ,          'text'         ,      'numeric'   ,   'text'   ,     'text'   , 'numeric']

hmm_hits_splits_table_name             = 'hmm_hits_in_splits'
hmm_hits_splits_table_structure        = ['hmm_hit_entry_id', 'split', 'percentage_in_split', 'source']
hmm_hits_splits_table_types            = [     'numeric'    ,  'text',       'numeric'      ,  'text' ]

scg_taxonomy_table_name                = 'scg_taxonomy'
scg_taxonomy_table_structure           = ['gene_callers_id', 'gene_name', 'accession', 'percent_identity', 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
scg_taxonomy_table_types               = [    'numeric'    ,    'text'  ,    'text'  ,       'text'      ,   'text'  ,   'text'  ,   'text' ,  'text'  ,   'text'  ,   'text' ,   'text'   ]

trna_taxonomy_table_name                = 'trna_taxonomy'
trna_taxonomy_table_structure           = ['gene_callers_id', 'amino_acid', 'anticodon', 'accession', 'percent_identity', 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
trna_taxonomy_table_types               = [    'numeric'    ,    'text'   ,    'text'  ,    'text'  ,       'text'      ,   'text'  ,   'text'  ,   'text' ,  'text'  ,   'text'  ,   'text' ,   'text'   ]

nt_position_info_table_name            = 'nt_position_info'
nt_position_info_table_structure       = ['contig_name', 'position_info']
nt_position_info_table_types           = [    'str'    ,      'blob'    ]

nucleotide_additional_data_table_name      = 'nucleotide_additional_data'
nucleotide_additional_data_table_structure = ['item_name', 'data_key', 'data_value', 'data_type', 'data_group']
nucleotide_additional_data_table_types     = [   'text'  ,   'text'  ,    'text'   ,    'text'  ,    'text'   ]

amino_acid_additional_data_table_name      = 'amino_acid_additional_data'
amino_acid_additional_data_table_structure = ['item_name', 'data_key', 'data_value', 'data_type', 'data_group']
amino_acid_additional_data_table_types     = [   'text'  ,   'text'  ,    'text'   ,    'text'  ,    'text'   ]

gene_level_coverage_stats_table_name      = 'gene_level_coverage_stats'
gene_level_coverage_stats_table_structure = ['gene_callers_id', 'sample_name', 'mean_coverage', 'detection', 'non_outlier_mean_coverage', 'non_outlier_coverage_std', 'gene_coverage_values_per_nt', 'non_outlier_positions']
gene_level_coverage_stats_table_types     = [    'numeric'    ,     'text'   ,    'numeric'   ,  'numeric' ,         'numeric'          ,          'numeric'        ,             'blob'           ,          'blob'        ]

gene_level_inseq_stats_table_name      = 'gene_level_inseq_stats'
gene_level_inseq_stats_table_structure = ['gene_callers_id', 'sample_name', 'mean_coverage', 'insertions', 'insertions_normalized', 'mean_disruption', 'below_disruption', 'gene_coverage_values_per_nt']
gene_level_inseq_stats_table_types     = [    'numeric'    ,     'text'   ,    'numeric'   ,  'numeric'  ,        'numeric'       ,      'numeric'   ,      'numeric'    ,              'blob'          ]

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE PROFILE DATABASE
#
####################################################################################################

item_orders_table_name               = 'item_orders'
item_orders_table_structure          = ['name', 'type', 'data', 'additional']
item_orders_table_types              = ['text', 'text', 'text',    'text'   ]

item_additional_data_table_name      = 'item_additional_data'
item_additional_data_table_structure = ['item_name', 'data_key', 'data_value', 'data_type', 'data_group']
item_additional_data_table_types     = [   'text'  ,   'text'  ,    'text'   ,    'text'  ,    'text'   ]

layer_orders_table_name              = 'layer_orders'
layer_orders_table_structure         = ['data_key', 'data_type', 'data_value']
layer_orders_table_types             = [  'text'  ,    'text'  ,    'text'   ]

layer_additional_data_table_name      = 'layer_additional_data'
layer_additional_data_table_structure = ['item_name', 'data_key', 'data_value', 'data_type', 'data_group']
layer_additional_data_table_types     = [   'text'  ,   'text'  ,    'text'   ,    'text'  ,    'text'   ]

states_table_name                    = 'states'
states_table_structure               = ['name', 'content', 'last_modified']
states_table_types                   = ['text',  'text'  ,      'text'    ]

variable_codons_table_name           = 'variable_codons'
variable_codons_table_structure      = ['sample_id', 'corresponding_gene_call', 'codon_order_in_gene', 'reference', 'departure_from_reference', 'coverage'] + codons
variable_codons_table_types          = [   'text'  ,        'numeric'         ,       'numeric'      ,    'text'  ,          'numeric'        , 'numeric' ] + ['numeric'] * len(codons)

variable_nts_table_name              = 'variable_nucleotides'
variable_nts_table_structure         = ['sample_id', 'split_name',   'pos'  , 'pos_in_contig', 'corresponding_gene_call', 'in_noncoding_gene_call', 'in_coding_gene_call', 'base_pos_in_codon', 'codon_order_in_gene', 'coverage', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'departure_from_reference', 'competing_nts', 'reference'] + nucleotides
variable_nts_table_types             = [   'text'  ,    'text'   , 'numeric',    'numeric'   ,        'numeric'         ,       'numeric'         ,     'numeric'        ,       'numeric'    ,       'numeric'      , 'numeric' ,          'bool'       ,          'bool'        ,          'numeric'        ,      'text'    ,    'text'  ] + ['numeric'] * len(nucleotides)

indels_table_name                    = 'indels'
indels_table_structure               = ['sample_id', 'split_name', 'pos'    , 'pos_in_contig', 'corresponding_gene_call', 'in_noncoding_gene_call', 'in_coding_gene_call' , 'base_pos_in_codon', 'codon_order_in_gene', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'reference', 'type', 'sequence', 'length' , 'count'  , 'coverage']
indels_table_types                   = ['text'     , 'text'      , 'integer', 'integer'      , 'integer'                , 'integer'               , 'integer'             , 'integer'          , 'integer'            , 'integer'             , 'integer'              , 'text'     , 'text', 'text'    , 'integer', 'integer', 'integer']

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
collections_bins_info_table_structure = ['collection_name', 'bin_name', 'source', 'html_color']
collections_bins_info_table_types     = [      'text'     ,   'text'  ,  'text' ,    'text'   ]

collections_contigs_table_name        = 'collections_of_contigs'
collections_contigs_table_structure   = ['collection_name', 'contig', 'bin_name']
collections_contigs_table_types       = [      'text'     ,  'text' ,   'text'  ]

collections_splits_table_name         = 'collections_of_splits'
collections_splits_table_structure    = ['collection_name', 'split', 'bin_name']
collections_splits_table_types        = [      'text'     , 'text' ,   'text'  ]


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

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE STRUCTURE DB
#
####################################################################################################

pdb_data_table_name       = 'structures'
pdb_data_table_structure  = ['corresponding_gene_call', 'pdb_content']
pdb_data_table_types      = [         'integer'       ,    'blob'    ]

templates_table_name       = 'templates'
templates_table_structure  = ['corresponding_gene_call' , 'pdb_id' , 'chain_id' , 'proper_percent_similarity', 'percent_similarity', 'align_fraction']
templates_table_types      = ['integer'                 , 'text'   , 'text'     , 'real',                      'real',               'real']

models_table_name       = 'models'
models_table_structure  = ['corresponding_gene_call' , 'molpdf' , 'GA341_score' , 'DOPE_score' , 'picked_as_best']
models_table_types      = ['integer'                 , 'real'   , 'real'        , 'real'       , 'integer']

residue_info_table_name       = 'residue_info'
residue_info_table_structure  = ['corresponding_gene_call', 'codon_order_in_gene', 'contact_numbers', 'codon', 'amino_acid', 'codon_number', 'sec_struct' , 'rel_solvent_acc' , 'phi'  , 'psi'  , 'NH_O_1_index' , 'NH_O_1_energy' , 'O_NH_1_index' , 'O_NH_1_energy' , 'NH_O_2_index' , 'NH_O_2_energy' , 'O_NH_2_index' , 'O_NH_2_energy']
residue_info_table_types      = [        'integer'        ,        'integer'     ,   'text'          , 'text',  'text',       'integer',      'text'       , 'real'            , 'real' , 'real' , 'integer'      , 'real'          , 'integer'      , 'real'          , 'integer'      , 'real'          , 'integer'      , 'real']

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE KEGG MODULES DB
#
####################################################################################################

module_table_name = "kegg_modules"
module_table_structure = ['module', 'data_name', 'data_value', 'data_definition', 'line']
module_table_types     = [ 'str'  ,   'str'    ,     'str'   ,       'str'      ,'numeric' ]

pathway_table_name = "kegg_pathway_maps"
pathway_table_structure = ['pathway_map', 'data_name', 'data_value', 'data_definition', 'line']
pathway_table_types     = [ 'str'  ,   'str'    ,     'str'   ,       'str'      ,'numeric' ]

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE TRNASEQ DB
#
####################################################################################################

trnaseq_sequences_table_name            = 'sequences'
trnaseq_sequences_table_structure       = ['name', 'read_count', 'profiled_or_mapped', 'sequence']
trnaseq_sequences_table_types           = ['str' , 'numeric'   , 'numeric'           , 'str']

trnaseq_feature_table_name              = 'feature'
trnaseq_feature_table_structure         = ['name', 'has_complete_feature_set', 'anticodon_sequence', 'amino_acid', 'sequence_length', 'features_start', 'features_stop', 'num_conserved', 'num_unconserved', 'num_paired', 'num_unpaired', 'num_in_extrapolated_fiveprime_feature', 'num_extra_fiveprime' , 'num_extra_threeprime'] + list(itertools.chain(*zip([f + '_start' for f in TRNA_FEATURE_NAMES], [f + '_stop' for f in TRNA_FEATURE_NAMES]))) + ['alpha_start', 'alpha_stop', 'beta_start', 'beta_stop']
trnaseq_feature_table_types             = ['str' , 'bool'                    , 'str'               , 'str'       , 'numeric'        , 'numeric'       , 'numeric'      , 'numeric'      , 'numeric'        , 'numeric'   , 'numeric'     , 'numeric'                              , 'numeric'             , 'numeric'             ] + ['str'] * len(TRNA_FEATURE_NAMES) * 2                                                                              + ['numeric'    , 'numeric'   , 'numeric'   , 'numeric']

trnaseq_unconserved_table_name          = 'feature_unconserved_nucleotides'
trnaseq_unconserved_table_structure     = ['name', 'pos'    , 'observed_nucleotide', 'expected_nucleotides']
trnaseq_unconserved_table_types         = ['str' , 'numeric', 'str'                , 'str']

trnaseq_unpaired_table_name             = 'feature_unpaired_nucleotides'
trnaseq_unpaired_table_structure        = ['name', 'fiveprime_pos', 'threeprime_pos', 'observed_fiveprime_nucleotide', 'observed_threeprime_nucleotide']
trnaseq_unpaired_table_types            = ['str' , 'numeric'      , 'numeric'       , 'str'                          , 'str']

trnaseq_trimmed_table_name              = 'trimmed'
trnaseq_trimmed_table_structure         = ['name', 'unique_seq_count', 'read_count', 'profiled_or_mapped', 'sequence', 'normalized_seq_representation', 'fiveprime_unique_seq_count', 'fiveprime_read_count'] + [threeprime_variant + '_read_count' for threeprime_variant in THREEPRIME_VARIANTS]
trnaseq_trimmed_table_types             = ['str' , 'numeric'         , 'numeric'   , 'numeric'           , 'str'     , 'numeric'                      , 'numeric'                   , 'numeric'             ] + ['numeric' for _ in THREEPRIME_VARIANTS]

trnaseq_normalized_table_name           = 'normalized'
trnaseq_normalized_table_structure      = ['name', 'trimmed_seq_count', 'mean_specific_coverage', 'mean_nonspecific_coverage', 'specific_coverages', 'nonspecific_coverages', 'modified_seq_representation', 'specific_read_count', 'nonspecific_read_count', 'count_of_specific_reads_with_extra_fiveprime', 'count_of_nonspecific_reads_with_extra_fiveprime', 'specific_mapped_read_count', 'nonspecific_mapped_read_count', 'specific_longer_fiveprime_extensions', 'specific_longer_fiveprime_extension_read_counts', 'nonspecific_longer_fiveprime_extensions', 'nonspecific_longer_fiveprime_extension_read_counts'] + [threeprime_variant + '_specific_read_count' for threeprime_variant in THREEPRIME_VARIANTS] + [threeprime_variant + '_nonspecific_read_count' for threeprime_variant in THREEPRIME_VARIANTS]
trnaseq_normalized_table_types          = ['str' , 'numeric'          , 'numeric'               , 'numeric'                  , 'str'               , 'str'                  , 'numeric'                    , 'numeric'            , 'numeric'               , 'numeric'                                     , 'numeric'                                        , 'numeric'                   , 'numeric'                      , 'str'                                 , 'str'                                            , 'str'                                    , 'str'                                               ] + ['numeric' for _ in THREEPRIME_VARIANTS                                                  ] + ['numeric' for _ in THREEPRIME_VARIANTS]

trnaseq_modified_table_name             = 'modified'
trnaseq_modified_table_structure        = ['name', 'mean_specific_coverage', 'mean_nonspecific_coverage', 'specific_coverages', 'nonspecific_coverages', 'substitution_positions', 'substitution_A_specific_coverage', 'substitution_C_specific_coverage', 'substitution_G_specific_coverage', 'substitution_T_specific_coverage', 'substitution_A_nonspecific_coverage', 'substitution_C_nonspecific_coverage', 'substitution_G_nonspecific_coverage', 'substitution_T_nonspecific_coverage', 'deletion_positions', 'deletion_specific_coverage', 'deletion_nonspecific_coverage', 'consensus_sequence', 'count_of_normalized_seqs_without_dels', 'names_of_normalized_seqs_without_dels', 'count_of_normalized_seqs_with_dels', 'names_of_normalized_seqs_with_dels', 'specific_read_count', 'nonspecific_read_count', 'count_of_specific_reads_without_extra_fiveprime', 'count_of_specific_reads_with_extra_fiveprime', 'specific_mapped_read_count', 'nonspecific_mapped_read_count', 'specific_longer_fiveprime_extensions', 'specific_longer_fiveprime_extension_read_counts', 'nonspecific_longer_fiveprime_extensions', 'nonspecific_longer_fiveprime_extension_read_counts'] + [threeprime_variant + '_specific_read_count' for threeprime_variant in THREEPRIME_VARIANTS] + [threeprime_variant + '_nonspecific_read_count' for threeprime_variant in THREEPRIME_VARIANTS]
trnaseq_modified_table_types            = ['str' , 'numeric'               , 'numeric'                  , 'str'               , 'str'                  , 'str'                   , 'str'                             , 'str'                             , 'str'                             , 'str'                             , 'str'                                , 'str'                                , 'str'                                , 'str'                                , 'str'               , 'str'                       , 'str'                          , 'str'               , 'numeric'                              , 'str'                                  , 'numeric'                           , 'str'                               , 'numeric'            , 'numeric'               , 'numeric'                                        , 'numeric'                                     , 'numeric'                   , 'numeric'                      , 'str'                                 , 'str'                                            , 'str'                                    , 'str'                                               ] + ['numeric' for _ in THREEPRIME_VARIANTS                                                  ] + ['numeric' for _ in THREEPRIME_VARIANTS]

####################################################################################################
#
#     META META META
#
####################################################################################################

tables_without_unique_entry_ids = [genome_gene_function_calls_table_name]

requires_unique_entry_id = {
    'self': False,
    'kmer_contigs': False,
    'kmer_splits': False,
    'atomic_data_splits': False,
    'atomic_data_contigs': False,
    'std_coverage_splits': False,
    'std_coverage_contigs': False,
    'mean_coverage_splits': False,
    'mean_coverage_contigs': False,
    'mean_coverage_Q2Q3_splits': False,
    'mean_coverage_Q2Q3_contigs': False,
    'max_normalized_ratio_splits': False,
    'max_normalized_ratio_contigs': False,
    'relative_abundance_splits': False,
    'relative_abundance_contigs': False,
    'detection_splits': False,
    'detection_contigs': False,
    'abundance_splits': False,
    'abundance_contigs': False,
    'variability_splits': False,
    'variability_contigs': False,
    'gene_cluster_frequencies': False,
    'gene_cluster_presence_absence': False,
    'clusters': False,
    'protein_clusters': False,
    'clusterings': False,
    'PC_frequencies': False,
    'PC_presence_absence': False,
    'additional_data': False,
    'samples_order': False,
    'samples_attribute_aliases': False,
    'samples_information': False,
    'variable_nucleotide_positions': False,
    'variable_amino_acid_frequencies': False,
    'gene_protein_sequences': False,
    'genes_in_splits_summary': False,
    pan_gene_clusters_table_name: True,
    genes_in_splits_table_name: True,
    gene_function_calls_table_name: True,
    hmm_hits_splits_table_name: True,
    scg_taxonomy_table_name: True,
    trna_taxonomy_table_name: True,
    nucleotide_additional_data_table_name: True,
    amino_acid_additional_data_table_name: True,
    gene_level_coverage_stats_table_name: True,
    gene_level_inseq_stats_table_name: True,
    item_additional_data_table_name: True,
    layer_additional_data_table_name: True,
    variable_codons_table_name: True,
    variable_nts_table_name: True,
    indels_table_name: True,
    collections_bins_info_table_name: True,
    collections_contigs_table_name: True,
    collections_splits_table_name: True,
    templates_table_name: True,
    models_table_name: True,
    residue_info_table_name: True,
    contig_sequences_table_name: False,
    contigs_info_table_name: False,
    splits_info_table_name: False,
    genes_in_contigs_table_name: False,
    gene_amino_acid_sequences_table_name: False,
    taxon_names_table_name: False,
    splits_taxonomy_table_name: False,
    genes_taxonomy_table_name: False,
    hmm_hits_info_table_name: False,
    hmm_hits_table_name: False,
    nt_position_info_table_name: False,
    item_orders_table_name: False,
    layer_orders_table_name: False,
    states_table_name: False,
    views_table_name: False,
    collections_info_table_name: False,
    split_coverages_table_name: False,
    genome_info_table_name: False,
    gene_info_table_name: False,
    genome_gene_function_calls_table_name: True,
    pdb_data_table_name: False,
    module_table_name: False,
    pathway_table_name: False,
    trnaseq_sequences_table_name: False,
    trnaseq_feature_table_name: False,
    trnaseq_unconserved_table_name: False,
    trnaseq_unpaired_table_name: False,
    trnaseq_trimmed_table_name: False,
    trnaseq_normalized_table_name: False,
    trnaseq_modified_table_name: False
}
