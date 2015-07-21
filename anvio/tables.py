# -*- coding: utf-8
""" Table schemas for databases."""


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


annotation_db_version = "3"
profile_db_version = "5"


####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE ANNOTATION DATABASE
#
####################################################################################################


contig_sequences_table_name          = 'contig_sequences'
contig_sequences_table_structure     = ['contig', 'sequence']
contig_sequences_table_types         = [  'str' ,   'str'   ]

contigs_info_table_name              = 'contigs_basic_info'
contigs_info_table_structure         = ['contig', 'length' , 'gc_content', 'num_splits']
contigs_info_table_types             = [  'str' , 'numeric',   'numeric' ,   'numeric' ]

splits_info_table_name               = 'splits_basic_info'
splits_info_table_structure          = ['split', 'order_in_parent' , 'start' ,  'end'  , 'length' , 'gc_content', 'gc_content_parent', 'parent' ]
splits_info_table_types              = ['text' ,     'numeric     ','numeric','numeric', 'numeric',   'numeric' ,      'numeric'     ,  'text'  ]

genes_contigs_table_name             = 'genes_in_contigs'
genes_contigs_table_structure        = ['prot', 'contig', 'start', 'stop'   , 'direction', 'figfam', 'function', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
genes_contigs_table_types            = ['text',  'text' ,'numeric','numeric',   'text'   ,  'text' ,   'text'  ,   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

genes_splits_summary_table_name      = 'genes_in_splits_summary'
genes_splits_summary_table_structure = ['split', 'taxonomy', 'num_genes', 'avg_gene_length', 'ratio_coding', 'ratio_hypothetical', 'ratio_with_tax', 'tax_accuracy']
genes_splits_summary_table_types     = [ 'text',   'text'  ,  'numeric' ,     'numeric'    ,   'numeric'   ,      'numeric'      ,     'numeric'   ,   'numeric'   ]

genes_splits_table_name              = 'genes_in_splits'
genes_splits_table_structure         = ['entry_id', 'split', 'prot', 'start_in_split', 'stop_in_split', 'percentage_in_split']
genes_splits_table_types             = [ 'numeric',  'text', 'text',    'numeric'    ,    'numeric'   ,       'numeric'      ]

hmm_hits_info_table_name             = 'hmm_hits_info'
hmm_hits_info_table_structure        = ['source', 'ref' , 'search_type', 'genes']
hmm_hits_info_table_types            = [ 'text' , 'text',    'text'    , 'text' ]

hmm_hits_contigs_table_name          = 'hmm_hits_in_contigs'
hmm_hits_contigs_table_structure     = ['entry_id', 'source', 'gene_unique_identifier', 'contig', 'start' , 'stop'  , 'gene_name', 'gene_id', 'e_value']
hmm_hits_contigs_table_types         = [ 'numeric',  'text' ,          'text'         ,  'text' ,'numeric','numeric',   'text'   ,  'text'  , 'numeric']

hmm_hits_splits_table_name           = 'hmm_hits_in_splits'
hmm_hits_splits_table_structure      = ['entry_id', 'source', 'gene_unique_identifier', 'gene_name', 'split', 'percentage_in_split', 'e_value']
hmm_hits_splits_table_types          = [ 'numeric',  'text' ,          'text'         ,   'text'   ,  'text',       'numeric'      , 'numeric']

collections_info_table_name          = 'collections_info'
collections_info_table_structure     = ['source', 'num_splits', 'num_clusters']
collections_info_table_types         = [ 'text' ,  'numeric'  ,   'numeric'   ]

collections_colors_table_name        = 'collections_colors'
collections_colors_table_structure   = ['entry_id', 'source', 'cluster_id', 'htmlcolor']
collections_colors_table_types       = [ 'numeric',  'text' ,    'text'   ,    'text'  ]

collections_contigs_table_name       = 'collections_of_contigs'
collections_contigs_table_structure  = ['entry_id', 'source', 'contig', 'cluster_id']
collections_contigs_table_types      = [ 'numeric',  'text' ,  'text' ,    'text'   ]

collections_splits_table_name        = 'collections_of_splits'
collections_splits_table_structure   = ['entry_id', 'source', 'split', 'cluster_id']
collections_splits_table_types       = [ 'numeric',  'text' , 'text' ,    'text'   ]

states_table_name                    = 'states'
states_table_structure               = ['name', 'content', 'last_modified']
states_table_types                   = ['text',  'text'  ,      'text'    ]

####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE ANNOTATION DATABASE
#
####################################################################################################

# notice that metadata table is the only table that doesn't have a name. because it is a bit tricky.
# for single profiles, metadata_table is stored as "metadata", however, for merged profiles,
# each metadata item in the metadata_table_structure becomes its own table, where the values of splits
# for a given table name class is stored accross samples.
metadata_table_structure             = ['contig', 'std_coverage', 'mean_coverage', 'normalized_coverage', 'max_normalized_ratio', 'relative_abundance', 'portion_covered', 'abundance', 'variability', '__parent__']
metadata_table_types                 = [ 'text' ,   'numeric'   ,    'numeric'   ,       'numeric'      ,        'numeric'      ,      'numeric'     ,     'numeric'    ,  'numeric' ,   'numeric'  ,    'text'   ]

clusterings_table_name               = 'clusterings'
clusterings_table_structure          = ['clustering', 'newick' ]
clusterings_table_types              = [   'str'    ,  'str'   ]

variable_positions_table_name        = 'variable_positions'
variable_positions_table_structure   = ['entry_id', 'sample_id', 'split_name',   'pos'  , 'coverage', 'n2n1ratio', 'competing_nts', 'consensus',    'A'   ,    'T'   ,    'C'   ,    'G'   ,    'N'   ]
variable_positions_table_types       = [ 'numeric',    'text'  ,    'text'   , 'numeric',  'numeric',  'numeric' ,      'text'    ,    'text'  , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']

gene_coverages_table_name            = 'gene_coverages'
gene_coverages_table_structure       = ['entry_id', 'prot', 'sample_id', 'mean_coverage']
gene_coverages_table_types           = [ 'numeric', 'text',   'text'   ,    'numeric'   ]

views_table_name                     = 'views'
views_table_structure                = ['view_id', 'target_table']
views_table_types                    = [  'str'  ,      'str'    ]
