# -*- coding: utf-8
# pylint: disable=line-too-long
""" Table schemas for databases."""

from anvio.constants import codons, nucleotides, essential_genome_info, TRNA_FEATURE_NAMES

import itertools

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


contigs_db_version = "24"
profile_db_version = "40"
genes_db_version = "6"
pan_db_version = "22"
auxiliary_data_version = "2"
structure_db_version = "2"
genomes_storage_vesion = "7"
trnaseq_db_version = "2"
workflow_config_version = "3"
metabolic_modules_db_version = "4"

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

pan_gc_tracker_table_name           = 'gc_tracker' # the purpose of this table is to keep track of the sequence-based GCs in structure mode
pan_gc_tracker_table_structure      = ['gene_caller_id', 'gene_cluster_id', 'genome_name', 'alignment_summary']
pan_gc_tracker_table_types          = [    'numeric'   ,      'str'       ,     'str'    ,        'str'       ]

pan_gc_psgc_associations_table_name      = 'gc_psgc_associations' # which sequence-based GCs ended up in which PSGCs
pan_gc_psgc_associations_table_structure = ['gene_cluster_id', 'protein_structure_informed_gene_cluster_id']
pan_gc_psgc_associations_table_types     = [     'str'       ,                    'str'                    ]

pan_reaction_network_reactions_table_name        = 'pan_reaction_network_reactions'
pan_reaction_network_reactions_table_structure   = ['modelseed_reaction_id', 'modelseed_reaction_name', 'ko_kegg_reaction_source', 'ko_ec_number_source', 'other_kegg_reaction_ids', 'other_ec_numbers', 'metabolite_modelseed_ids', 'stoichiometry', 'compartments', 'reversibility']
pan_reaction_network_reactions_table_types       = [         'text'        ,            'text'        ,          'text'          ,         'text'       ,         'text'           ,       'text'      ,           'text'          ,      'text'    ,     'text'    ,      'bool'    ]

pan_reaction_network_metabolites_table_name      = 'pan_reaction_network_metabolites'
pan_reaction_network_metabolites_table_structure = ['modelseed_compound_id', 'modelseed_compound_name', 'kegg_aliases', 'formula', 'charge' , 'smiles']
pan_reaction_network_metabolites_table_types     = [         'text'        ,           'text'         ,     'text'    ,   'text' , 'numeric',  'text' ]

pan_reaction_network_kegg_table_name             = 'pan_reaction_network_kegg'
pan_reaction_network_kegg_table_structure        = ['kegg_id', 'name', 'modules', 'pathways', 'brite_categorization']
pan_reaction_network_kegg_table_types            = [ 'text'  , 'text',  'text'  ,   'text'  ,         'text'        ]

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

reaction_network_reactions_table_name        = 'reaction_network_reactions'
reaction_network_reactions_table_structure   = ['modelseed_reaction_id', 'modelseed_reaction_name', 'ko_kegg_reaction_source', 'ko_ec_number_source', 'other_kegg_reaction_ids', 'other_ec_numbers', 'metabolite_modelseed_ids', 'stoichiometry', 'compartments', 'reversibility']
reaction_network_reactions_table_types       = [         'text'        ,            'text'        ,           'text'         ,         'text'       ,           'text'         ,       'text'      ,          'text'           ,      'text'    ,     'text'    ,      'bool'    ]

reaction_network_metabolites_table_name      = 'reaction_network_metabolites'
reaction_network_metabolites_table_structure = ['modelseed_compound_id', 'modelseed_compound_name', 'kegg_aliases', 'formula', 'charge' , 'smiles']
reaction_network_metabolites_table_types     = [         'text'        ,           'text'         ,     'text'    ,   'text' , 'numeric',  'text' ]

reaction_network_kegg_table_name             = 'reaction_network_kegg'
reaction_network_kegg_table_structure        = ['kegg_id', 'name', 'modules', 'pathways', 'brite_categorization']
reaction_network_kegg_table_types            = [ 'text'  , 'text',  'text'  ,   'text'  ,         'text'        ]

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
#     ADDITIONAL TABLE DESCRIPTIONS FOR THE TRNASEQ VARIANT OF THE CONTIGS DATABASE
#
####################################################################################################

trna_seed_feature_table_name            = 'trna_feature'
trna_seed_feature_table_structure       = ['gene_callers_id'] + list(itertools.chain(*zip([f + '_start' for f in TRNA_FEATURE_NAMES[: -1]], [f + '_stop' for f in TRNA_FEATURE_NAMES[: -1]]))) + ['alpha_start', 'alpha_stop', 'beta_start', 'beta_stop']
trna_seed_feature_table_types           = ['gene_callers_id'] + ['str'] * len(TRNA_FEATURE_NAMES[: -1]) * 2                                                                                    + ['numeric'    , 'numeric'   , 'numeric'   , 'numeric']


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

# names of view tables are automatically generated as a function of the atomic data
# generated by anvi-profiler. a list of current data items turned into view talbes
# can be found in `contsants.essential_data_fields_for_anvio_profiles`
view_table_structure = ['item', 'layer',  'value' ]
view_table_types     = ['text', 'text' , 'numeric']

protein_abundances_table_name        = 'protein_abundances'
protein_abundances_table_structure   = ['protein_id', 'reference_source', 'reference_id', 'sample_name', 'abundance_value']
protein_abundances_table_types       = [  'numeric' ,       'text'      ,     'text'    ,     'text'   ,     'numeric'    ]

metabolite_abundances_table_name      = 'metabolite_abundances'
metabolite_abundances_table_structure = ['reference_source', 'reference_id', 'sample_name', 'abundance_value']
metabolite_abundances_table_types     = [      'text'      ,     'text'    ,     'text'   ,     'numeric'    ]

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

module_table_name       = "modules"
module_table_structure  = ['module', 'data_name', 'data_value', 'data_definition', 'line']
module_table_types      = [ 'str'  ,   'str'    ,     'str'   ,       'str'      ,'numeric' ]

pathway_table_name      = "kegg_pathway_maps"
pathway_table_structure = ['pathway_map', 'data_name', 'data_value', 'data_definition', 'line']
pathway_table_types     = [ 'str'       ,   'str'    ,     'str'   ,       'str'      ,'numeric' ]

brite_table_name        = "brite_hierarchies"
brite_table_structure   = ['hierarchy_accession', 'hierarchy_name', 'ortholog_accession', 'ortholog_name', 'categorization']
brite_table_types       = [        'str'        ,       'str'     ,         'str'       ,      'str'     ,      'str'      ]

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE TRNASEQ DB
#
####################################################################################################

trnaseq_sequences_table_name            = 'sequences'
trnaseq_sequences_table_structure       = ['name', 'read_count', 'id_info', 'sequence']
trnaseq_sequences_table_types           = ['str' , 'numeric'   , 'str'    , 'str']

trnaseq_feature_table_name              = 'feature'
trnaseq_feature_table_structure         = ['name', 'has_complete_feature_set', 'anticodon_sequence', 'amino_acid', 'sequence_length', 'features_start', 'num_conserved', 'num_unconserved', 'num_paired', 'num_unpaired', 'num_in_extrapolated_fiveprime_feature', 'num_extra_fiveprime' , 'acceptor_length'] + list(itertools.chain(*zip([f + '_start' for f in TRNA_FEATURE_NAMES], [f + '_stop' for f in TRNA_FEATURE_NAMES]))) + ['alpha_start', 'alpha_stop', 'beta_start', 'beta_stop']
trnaseq_feature_table_types             = ['str' , 'bool'                    , 'str'               , 'str'       , 'numeric'        , 'numeric'       , 'numeric'      , 'numeric'        , 'numeric'   , 'numeric'     , 'numeric'                              , 'numeric'             , 'numeric'        ] + ['str'] * len(TRNA_FEATURE_NAMES) * 2                                                                              + ['numeric'    , 'numeric'   , 'numeric'   , 'numeric']

trnaseq_unconserved_table_name          = 'unconserved'
trnaseq_unconserved_table_structure     = ['name', 'pos'    , 'observed_nucleotide', 'expected_nucleotides']
trnaseq_unconserved_table_types         = ['str' , 'numeric', 'str'                , 'str']

trnaseq_unpaired_table_name             = 'unpaired'
trnaseq_unpaired_table_structure        = ['name', 'fiveprime_pos', 'threeprime_pos', 'observed_fiveprime_nucleotide', 'observed_threeprime_nucleotide']
trnaseq_unpaired_table_types            = ['str' , 'numeric'      , 'numeric'       , 'str'                          , 'str']

trnaseq_trimmed_table_name              = 'trimmed'
trnaseq_trimmed_table_structure         = ['name', 'unique_seq_count', 'read_count', 'id_info', 'sequence', 'normalized_seq_representation', 'fiveprime_unique_seq_count', 'fiveprime_read_count', 'threeprime_termini', 'threeprime_terminus_read_counts']
trnaseq_trimmed_table_types             = ['str' , 'numeric'         , 'numeric'   , 'str'    , 'str'     , 'numeric'                      , 'numeric'                   , 'numeric'             , 'str'               , 'str']

trnaseq_normalized_table_name           = 'normalized'
trnaseq_normalized_table_structure      = ['name', 'trimmed_seq_count', 'id_info', 'mean_specific_coverage', 'mean_nonspecific_coverage', 'specific_coverages', 'nonspecific_coverages', 'specific_read_count', 'nonspecific_read_count', 'count_of_specific_reads_with_extra_fiveprime', 'count_of_nonspecific_reads_with_extra_fiveprime', 'specific_mapped_read_count', 'nonspecific_mapped_read_count', 'specific_long_fiveprime_extensions', 'specific_long_fiveprime_extension_read_counts', 'nonspecific_long_fiveprime_extensions', 'nonspecific_long_fiveprime_extension_read_counts', 'specific_threeprime_termini', 'specific_threeprime_terminus_read_counts', 'nonspecific_threeprime_termini', 'nonspecific_threeprime_terminus_read_counts']
trnaseq_normalized_table_types          = ['str' , 'numeric'          , 'str'    , 'numeric'               , 'numeric'                  , 'str'               , 'str'                  , 'numeric'            , 'numeric'               , 'numeric'                                     , 'numeric'                                        , 'numeric'                   , 'numeric'                      , 'str'                               , 'str'                                          , 'str'                                  , 'str'                                             , 'str'                        , 'str'                                     , 'str'                           , 'str']

trnaseq_modified_table_name             = 'modified'
trnaseq_modified_table_structure        = ['name', 'mean_specific_coverage', 'mean_nonspecific_coverage', 'specific_coverages', 'nonspecific_coverages', 'substitution_positions', 'substitution_A_specific_coverage', 'substitution_C_specific_coverage', 'substitution_G_specific_coverage', 'substitution_T_specific_coverage', 'substitution_A_nonspecific_coverage', 'substitution_C_nonspecific_coverage', 'substitution_G_nonspecific_coverage', 'substitution_T_nonspecific_coverage', 'insertion_starts', 'insertion_seqs', 'insertion_specific_coverages', 'insertion_nonspecific_coverages', 'deletion_starts', 'deletion_lengths', 'deletion_specific_coverages', 'deletion_nonspecific_coverages', 'consensus_sequence', 'count_of_normalized_seqs_without_indels', 'names_of_normalized_seqs_without_indels', 'count_of_normalized_seqs_with_indels', 'names_of_normalized_seqs_with_indels', 'specific_read_count', 'nonspecific_read_count', 'count_of_specific_reads_with_extra_fiveprime', 'count_of_nonspecific_reads_with_extra_fiveprime', 'specific_mapped_read_count', 'nonspecific_mapped_read_count', 'specific_long_fiveprime_extensions', 'specific_long_fiveprime_extension_read_counts', 'nonspecific_long_fiveprime_extensions', 'nonspecific_long_fiveprime_extension_read_counts', 'specific_threeprime_termini', 'specific_threeprime_terminus_read_counts', 'nonspecific_threeprime_termini', 'nonspecific_threeprime_terminus_read_counts']
trnaseq_modified_table_types            = ['str' , 'numeric'               , 'numeric'                  , 'str'               , 'str'                  , 'str'                   , 'str'                             , 'str'                             , 'str'                             , 'str'                             , 'str'                                , 'str'                                , 'str'                                , 'str'                                , 'str'             , 'str'           , 'str'                         , 'str'                            , 'str'            , 'str'             , 'str'                        , 'str'                           , 'str'               , 'numeric'                                , 'str'                                    , 'numeric'                             , 'str'                                 , 'numeric'            , 'numeric'               , 'numeric'                                     , 'numeric'                                        , 'numeric'                   , 'numeric'                      , 'str'                               , 'str'                                          , 'str'                                  , 'str'                                             , 'str'                        , 'str'                                     , 'str'                           , 'str']

####################################################################################################
#
#     META META META
#
####################################################################################################

# DO NOT remove any table name here, even if a table is no longer in use. any table that has ever used
# in any anvi'o database must remain here for migration scripts to work properly:
table_requires_unique_entry_id = {'self': False,
                                  'kmer_contigs': False,
                                  'kmer_splits': False,
                                  'std_coverage_splits': True,
                                  'std_coverage_contigs': True,
                                  'mean_coverage_splits': True,
                                  'mean_coverage_contigs': True,
                                  'mean_coverage_Q2Q3_splits': True,
                                  'mean_coverage_Q2Q3_contigs': True,
                                  'detection_splits': True,
                                  'detection_contigs': True,
                                  'abundance_splits': True,
                                  'abundance_contigs': True,
                                  'variability_splits': True,
                                  'variability_contigs': True,
                                  'gene_cluster_frequencies': True,
                                  'gene_cluster_presence_absence': True,
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
                                  'gene_coverages': False,
                                  'mean_coverage_Q1Q3_splits': False,
                                  'mean_coverage_Q1Q3_contigs': False,
                                  'portion_covered_contigs': False,
                                  'portion_covered_splits': False,
                                  'frequency_view': False,           #
                                  'presence_absence_view': False,    # These two tables are for pangeomes
                                  'functions_frequency_view': True,        #
                                  'functions_presence_absence_view': True, # And these two are for anvi-display-functions stuff
                                  'atomic_data_splits': False,
                                  'atomic_data_contigs': False,
                                  'max_normalized_ratio_contigs': False,
                                  'relative_abundance_contigs': False,
                                  'max_normalized_ratio_splits': False,
                                  'relative_abundance_splits': False,
                                  pan_gene_clusters_table_name: True,
                                  pan_gc_tracker_table_name: True,
                                  pan_gc_psgc_associations_table_name: False, # the first entry is always the de novo gene cluster name, which should be unique
                                  'gene_cluster_function_reactions': False, # renamed to 'pan_reaction_network_reactions'
                                  pan_reaction_network_reactions_table_name: False,
                                  'gene_cluster_function_metabolites': False, # renamed to 'pan_reaction_network_metabolites'
                                  pan_reaction_network_metabolites_table_name: False,
                                  pan_reaction_network_kegg_table_name: False,
                                  genes_in_splits_table_name: True,
                                  gene_function_calls_table_name: True,
                                  'gene_function_reactions': False, # renamed to 'reaction_network_reactions'
                                  reaction_network_reactions_table_name: False,
                                  'gene_function_metabolites': False, # renamed to 'reaction_network_metabolites'
                                  reaction_network_metabolites_table_name: False,
                                  reaction_network_kegg_table_name: False,
                                  hmm_hits_splits_table_name: True,
                                  scg_taxonomy_table_name: True,
                                  trna_taxonomy_table_name: True,
                                  nucleotide_additional_data_table_name: True,
                                  amino_acid_additional_data_table_name: True,
                                  gene_level_coverage_stats_table_name: True,
                                  gene_level_inseq_stats_table_name: True,
                                  trna_seed_feature_table_name: False,
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
                                  protein_abundances_table_name: False,
                                  metabolite_abundances_table_name: False,
                                  collections_info_table_name: False,
                                  split_coverages_table_name: False,
                                  genome_info_table_name: False,
                                  gene_info_table_name: False,
                                  genome_gene_function_calls_table_name: True,
                                  pdb_data_table_name: False,
                                  'kegg_modules': False,        # no longer in use as of modules db v3
                                  module_table_name: False,
                                  pathway_table_name: False,
                                  brite_table_name: False,
                                  trnaseq_sequences_table_name: False,
                                  trnaseq_feature_table_name: False,
                                  trnaseq_unconserved_table_name: False,
                                  trnaseq_unpaired_table_name: False,
                                  trnaseq_trimmed_table_name: False,
                                  trnaseq_normalized_table_name: False,
                                  trnaseq_modified_table_name: False,
        }


def is_known_table(table_name):
    # these are transitionary tables programmers often use in migration scripts
    if table_name.endswith('_TEMP'):
        return True

    if table_name in table_requires_unique_entry_id:
        return True
    else:
        return False


def is_table_requires_unique_entry_id(table_name):
    if table_name.endswith('_TEMP'):
        return False

    if table_name not in table_requires_unique_entry_id:
        raise Exception(f"You and your table '{table_name}' are lost :(")

    return table_requires_unique_entry_id[table_name]
