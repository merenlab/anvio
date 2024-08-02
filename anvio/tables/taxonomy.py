# -*- coding: utf-8
# pylint: disable=line-too-long

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table
from anvio.errors import ConfigError
from anvio.tables.genecalls import GenesInSplits


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TaxonNamesTable(object):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

 
    def populate_taxon_names_table(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        db_entries = [tuple([t_name_id] + [self.taxon_names_dict[t_name_id][t_level] for t_level in t.taxon_names_table_structure[1:]]) for t_name_id in self.taxon_names_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.taxon_names_table_name, db_entries)

        database.disconnect()
        self.run.info('Taxon names table', 'Updated with %d unique taxon names' % len(db_entries))


class TablesForGeneLevelTaxonomy(Table, TaxonNamesTable):
    """Populate all tables with taxonomy information.

       Essentially takes in a dictionary of genes and taxon calls, populates three tables: taxon_names,
       splits_taxonomy, and gene_taxonomy when 'create' is called."""

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        utils.is_contigs_db(self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress)
        TaxonNamesTable.__init__(self, self.db_path, self.run, self.progress)

        # this class keeps track of genes that occur in splits, and responsible
        # for generating the necessary table in the contigs database
        self.genes_in_splits = GenesInSplits()


    def create(self, genes_taxonomy_dict, taxon_names_dict, source='unkown source'):
        self.source = source

        if not self.genes_are_called:
            raise ConfigError("Something is wrong. The contigs database says that genes were now called, and here "
                               "you are trying to populate taxonomy tables for genes. No, thanks.")

        self.init_gene_calls_dict()

        self.genes_taxonomy_dict = genes_taxonomy_dict
        self.taxon_names_dict = taxon_names_dict

        self.sanity_check()

        # oepn connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        self.splits_info = database.get_table_as_dict(t.splits_info_table_name)

        # test whether there are already genes tables populated
        taxonomy_source = database.get_meta_value('gene_level_taxonomy_source')
        if taxonomy_source:
            self.run.warning('Previous taxonomy information from "%s" is being replaced with the incoming data '
                             'through "%s".' % (taxonomy_source, self.source))
            database._exec('''DELETE FROM %s''' % (t.splits_taxonomy_table_name))
            database._exec('''DELETE FROM %s''' % (t.taxon_names_table_name))
            database._exec('''DELETE FROM %s''' % (t.genes_taxonomy_table_name))

        # populate taxon mames table
        self.populate_taxon_names_table()

        # populate genes taxonomy table
        self.populate_genes_taxonomy_table()

        # compute and push split taxonomy information.
        self.populate_splits_taxonomy_table()

        # set the source
        database.remove_meta_key_value_pair('gene_level_taxonomy_source')
        database.set_meta_value('gene_level_taxonomy_source', self.source)

        # disconnect like a pro.
        database.disconnect()


    def populate_genes_taxonomy_table(self):
        # open connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # push taxonomy data
        db_entries = [(gene_callers_id, self.genes_taxonomy_dict[gene_callers_id]) for gene_callers_id in self.genes_taxonomy_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.genes_taxonomy_table_name, db_entries)

        # disconnect
        database.disconnect()

        self.run.info('Genes taxonomy table', 'Taxonomy stored for %d gene calls' % len(db_entries))


    def sanity_check(self):
        """Basic checks to make sure things are in order at least to a minimum extent."""

        self.progress.new('Sanity checking ...')
        self.progress.update('Comparing gene caller ids in the incoming data and in the contigs database ..')
        gene_caller_ids_in_database = set(self.gene_calls_dict.keys())
        gene_caller_ids_in_taxonomy_dict = set(self.genes_taxonomy_dict.keys())
        gene_caller_ids_missing_in_db = gene_caller_ids_in_taxonomy_dict.difference(gene_caller_ids_in_database)
        self.progress.end()

        run.info("Num gene caller ids in the db", pp(len(gene_caller_ids_in_database)))
        run.info("Num gene caller ids in the incoming data", pp(len(gene_caller_ids_in_taxonomy_dict)))

        if gene_caller_ids_missing_in_db:
            raise ConfigError("Taxonomy information for genes you are trying to import into the database contains "
                               "%s gene caller ids that do not appear to be in the database. This is a step you must "
                               "be very careful to make sure you are not importing annotations for genes that have "
                               "nothing to do with your contigs database. To make sure of that, you should always work "
                               "with `anvi-get-dna-sequences-for-gene-calls` or `anvi-get-aa-sequences-for-gene-calls` programs "
                               "to get the data to annotate. For instance one of the gene caller ids you have in your "
                               "input data that does not appear in the database is this one: '%s'. Anvi'o hopes it makes "
                               "sense to you, because it definitely does not make any sense to anvi'o :("\
                                                        % (len(gene_caller_ids_missing_in_db), str(gene_caller_ids_missing_in_db.pop())))

        # check whether input matrix dict
        keys_found =  list(self.taxon_names_dict.values())[0].keys()
        missing_keys = [key for key in t.taxon_names_table_structure[1:] if key not in keys_found]
        if len(missing_keys):
            raise ConfigError("Anvi'o is trying to get ready to create tables for taxonomy, but there is something "
                               "wrong :( The taxonomy names dict (one of the required input dictionaries to the class "
                               "seems to be missing a one or more keys that are necessary to finish the job. Here is "
                               "a list of missing keys: %s. The complete list of input keys should contain these: %s."\
                                        % (', '.join(missing_keys), ', '.join(t.taxon_names_table_structure[1:])))

        if not len(self.taxon_names_dict):
            raise ConfigError("Anvi'o is trying to get ready to create tables for taxonomy, but taxonomy names dict "
                               "(one of the required input dictionaries to the class responsible for this task) seems "
                               "to be empty.")


    def populate_splits_taxonomy_table(self):
        """Populate the taxonomy information per split"""

        # build a dictionary for fast access to all genes identified within a contig
        gene_caller_ids_in_contigs = {}
        for gene_callers_id in self.genes_taxonomy_dict:
            contig = self.gene_calls_dict[gene_callers_id]['contig']
            if contig in gene_caller_ids_in_contigs:
                gene_caller_ids_in_contigs[contig].add(gene_callers_id)
            else:
                gene_caller_ids_in_contigs[contig] = set([gene_callers_id])

        contigs_without_annotation = list(set(self.contigs_info.keys()) - set(gene_caller_ids_in_contigs.keys()))

        for contig in contigs_without_annotation:
            gene_caller_ids_in_contigs[contig] = set([])

        splits_dict = {}

        num_splits_processed = 0
        num_splits_with_taxonomy = 0

        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                num_splits_processed += 1

                splits_dict[split_name] = None
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                taxon_name_ids = []
                for gene_callers_id in gene_caller_ids_in_contigs[contig]:
                    if self.gene_calls_dict[gene_callers_id]['stop'] > start and self.gene_calls_dict[gene_callers_id]['start'] < stop:
                        taxon_name_ids.append(self.genes_taxonomy_dict[gene_callers_id])

                if not taxon_name_ids:
                    continue

                if len(set(taxon_name_ids)) == 1:
                    splits_dict[split_name] = taxon_name_ids[0]
                else:
                    d = Counter()
                    for taxon_name_id in taxon_name_ids:
                        d[taxon_name_id] += 1

                    most_frequent_taxon_name_id, occurrence = d.most_common()[0]
                    splits_dict[split_name] = most_frequent_taxon_name_id

                num_splits_with_taxonomy += 1

        # open connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # push taxonomy data
        db_entries = [(split, splits_dict[split]) for split in splits_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.splits_taxonomy_table_name, db_entries)

        # disconnect
        database.disconnect()

        self.run.info('Splits taxonomy', 'Input data from "%s" annotated %d of %d splits (%.1f%%) with taxonomy.'\
                                            % (self.source, num_splits_with_taxonomy, num_splits_processed,
                                               num_splits_with_taxonomy * 100.0 / num_splits_processed))




        
