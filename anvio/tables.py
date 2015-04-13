# -*- coding: utf-8

""" Table schemas for databases."""

import os

import anvio.db as db
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.utils import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()

# ANNOTATION DATABASE TABLES:

annotation_db_version = "1.0.0"

contig_sequences_table_name          = 'contig_sequences'
contig_sequences_table_structure     = ['contig', 'sequence']
contig_sequences_table_types         = [  'str' ,   'str'   ]

contigs_info_table_name              = 'contigs_basic_info'
contigs_info_table_structure         = ['contig', 'length' , 'gc_content']
contigs_info_table_types             = [  'str' , 'numeric',   'numeric' ]

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
hmm_hits_contigs_table_structure     = ['entry_id', 'source', 'contig', 'start' , 'stop'  , 'gene_name', 'gene_id', 'e_value']
hmm_hits_contigs_table_types         = [ 'numeric',  'text' ,  'text' ,'numeric','numeric',   'text'   ,  'text'  , 'numeric']

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


# PRFOFILE DATABASE TABLES

profile_db_version = "1.0.0"

clusterings_table_name               = 'clusterings'
clusterings_table_structure          = ['clustering', 'newick' ]
clusterings_table_types              = [   'str'    ,  'str'   ]

variable_positions_table_name       = 'variable_positions'
variable_positions_table_structure  = ['entry_id', 'sample_id', 'split_name',   'pos'  , 'coverage', 'n2n1ratio', 'competing_nts', 'consensus',    'A'   ,    'T'   ,    'C'   ,    'G'   ,    'N'   ]
variable_positions_table_types      = [ 'numeric',    'text'  ,    'text'   , 'numeric',  'numeric',  'numeric' ,      'text'    ,    'text'  , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']

gene_coverages_table_name            = 'gene_coverages'
gene_coverages_table_structure       = ['entry_id', 'prot', 'sample_id', 'mean_coverage']
gene_coverages_table_types           = [ 'numeric', 'text',   'text'   ,    'numeric'   ]


class Table(object):
    """Superclass for rudimentary needs and operations for annotation db tables"""
    def __init__(self, db_path, version, run=run, progress=progress, quiet = False):
        if not os.path.exists(db_path):
            raise ConfigError, "Database ('%s') does not exist. You must create one first." % db_path

        self.quiet = quiet
        self.db_type = None
        self.db_path = db_path
        self.version = version
        self.next_available_id = {}

        self.splits_info = None
        self.contigs_info = None
        self.split_length = None

        self.run = run
        self.progress = progress

        database = db.DB(self.db_path, version)
        self.db_type = database.get_meta_value('db_type')

        if self.db_type == 'annotation': 
            # FIXME: a better design is required. the salient point is, "Table" must serve for both profile db
            # and annotation db calls.
            self.split_length = database.get_meta_value('split_length')
            self.contigs_info = database.get_table_as_dict(contigs_info_table_name)
            self.splits_info  = database.get_table_as_dict(splits_info_table_name)

            self.contig_name_to_splits = {}
            for split_name in self.splits_info:
                parent = self.splits_info[split_name]['parent']
                if self.contig_name_to_splits.has_key(parent):
                    self.contig_name_to_splits[parent].append(split_name)
                else:
                    self.contig_name_to_splits[parent] = [split_name]

        database.disconnect()


    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError, "If you need unique ids, you must call 'set_next_available_id' first"

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1


    def set_next_available_id(self, table):
        database = db.DB(self.db_path, self.version)
        table_content = database.get_table_as_dict(table)
        if table_content:
            self.next_available_id[table] = max(table_content.keys()) + 1
        else:
            self.next_available_id[table] = 0

        database.disconnect()


    def export_contigs_in_db_into_FASTA_file(self):
        if self.db_type != 'annotation':
            return None

        database = db.DB(self.db_path, self.version)
        contig_sequences_table = database.get_table_as_dict(contig_sequences_table_name)
        database.disconnect()

        self.progress.new('Exporting contigs into a FASTA file')
        self.progress.update('...')
        contigs_fasta_path = os.path.join(filesnpaths.get_temp_directory_path(), 'contigs.fa')
        contigs_fasta = u.FastaOutput(contigs_fasta_path)
        for contig in contig_sequences_table:
            contigs_fasta.write_id(contig)
            contigs_fasta.write_seq(contig_sequences_table[contig]['sequence'], split=False)

        self.progress.end()
        self.run.info('FASTA for contigs', contigs_fasta_path)

        return contigs_fasta_path


    def delete_entries_for_key(self, table_column, key, tables_to_clear = []):
        # FIXME: this shoudl be in db.py
        # removes rows from each table in 'tables_to_remove' where 'table_column' equals 'value'
        database = db.DB(self.db_path, self.version)

        table_content = database.get_table_as_dict(tables_to_clear[0])
        if key in table_content:
            self.run.warning('Previous entries for "%s" is being removed from "%s"' % (key, ', '.join(tables_to_clear)))
            for table_name in tables_to_clear:
                database._exec('''DELETE FROM %s WHERE %s = "%s"''' % (table_name, table_column, key))

        database.disconnect()
