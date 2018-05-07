# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TableForGeneFunctions(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        self.run = run
        self.progress = progress

        utils.is_contigs_db(self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        self.set_next_available_id(t.gene_function_calls_table_name)


    def create(self, functions_dict, drop_previous_annotations_first = False):
        self.sanity_check()

        # incoming stuff:
        gene_function_sources = set([v['source'] for v in list(functions_dict.values())])
        unique_num_genes = len(set([v['gene_callers_id'] for v in list(functions_dict.values())]))

        # oepn connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # are there any previous annotations in the db:
        gene_function_sources_in_db = database.get_meta_value('gene_function_sources')
        gene_function_sources_in_db = set(gene_function_sources_in_db.split(',') if gene_function_sources_in_db else [])

        # difference between sources in the db, and incoming sources:
        gene_function_sources_both_in_db_and_incoming_dict = gene_function_sources.intersection(gene_function_sources_in_db)

        # here we will do some magic. there are mulitple scenarios to consider here based on whether there
        # are functions already in the database, whether some of them matches to the incoming functions, etc.
        # let's go case-by-case:
        if not gene_function_sources_in_db:
            # set the sources and continue
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources)))

        elif gene_function_sources_in_db and drop_previous_annotations_first:
            # there are gene calls, but the user wants everything to be dropeped.
            self.run.warning("As per your request, anvi'o is DROPPING all previous function calls from %d sources\
                              before adding the incoming data, which contains %d entries originating from %d sources: %s" \
                                    % (len(gene_function_sources_in_db), len(functions_dict),
                                       len(gene_function_sources), ', '.join(gene_function_sources)))

            # clean the table and reset the next available ids
            database._exec('''DELETE FROM %s''' % (t.gene_function_calls_table_name))
            self.reset_next_available_id_for_table(t.gene_function_calls_table_name)

            # set the sources
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(gene_function_sources))

        elif gene_function_sources_in_db and gene_function_sources_both_in_db_and_incoming_dict:
            # some of the functions in the incoming dict match to what is already in the db. remove
            self.run.warning("Some of the annotation sources you want to add into the database are already in the db. So\
                              anvi'o will REPLACE those with the incoming data from these sources: %s" % \
                                            ', '.join(gene_function_sources_both_in_db_and_incoming_dict))

            # remove those entries for matching sources:
            for source in gene_function_sources_both_in_db_and_incoming_dict:
                database._exec('''DELETE FROM %s WHERE source = "%s"''' % (t.gene_function_calls_table_name, source))

            # set the sources
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources_in_db.union(gene_function_sources))))

        else:
            # fuctions in the db, but none of them match with the incoming annotation sources. totally new stuff.
            # good then. update sources
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources_in_db.union(gene_function_sources))))

        # push the data
        db_entries = [tuple([self.next_id(t.gene_function_calls_table_name)] + [functions_dict[v][h] for h in t.gene_function_calls_table_structure[1:]]) for v in functions_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.gene_function_calls_table_name, db_entries)

        # disconnect like a pro.
        database.disconnect()

        self.run.info('Gene functions', '%d function calls from %d sources for %d unique gene calls has\
                                        been added to the contigs database.' % \
                                            (len(functions_dict), len(gene_function_sources), unique_num_genes))


    def sanity_check(self):
        pass



