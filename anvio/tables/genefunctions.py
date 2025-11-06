# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.tables.tableops import Table


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
P = terminal.pluralize


class TableForGeneFunctions(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        self.run = run
        self.progress = progress

        utils.is_contigs_db(self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)


    def add_empty_sources_to_functional_sources(self, gene_function_sources):
        if type(gene_function_sources) is not set:
            raise ConfigError('The programmer who called this function forgot that gene_function_sources must be of '
                              'type %s. If this is not your falut, please contact an anvi\'o developer.' % set)
        # open connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        self.add_new_sources_to_functional_sources(gene_function_sources, database)

        # disconnect like a pro.
        database.disconnect()


    def get_function_sources_in_db(self, database):
        gene_function_sources_in_db = database.get_meta_value('gene_function_sources')
        return set(gene_function_sources_in_db.split(',') if gene_function_sources_in_db else [])


    def drop_functions(self, database, sources_to_drop=[]):
        if not len(sources_to_drop):
            self.run.warning("Someone called 'drop functions' function with an empty list of sources. Cray. "
                             "pretending that it didn't happen, and returning gracefully.")
            return

        sources_in_db = self.get_function_sources_in_db(database)

        if not len(sources_in_db):
            self.run.warning(f"Someone called 'drop functions' with {P('source', len(sources_to_drop))} to drop, but then "
                             f"there is are no such functions in the database ANYWAY. Anvi'o is confused, but will ignore this and "
                             f"return gracefully. Things may explode downstream yo.")
            return

        sources_not_known = [s for s in sources_to_drop if s not in sources_in_db]
        sources_to_remain = [s for s in sources_in_db if s not in sources_to_drop]
        sources_to_drop = [s for s in sources_to_drop if s in sources_in_db]

        if len(sources_not_known) and len(sources_to_drop):
            self.run.warning(f"Some sources the 'drop functions' was asked to drop are not in the database ({', '.join(sources_not_known)}). "
                             f"But others are ({', '.join(sources_to_drop)}). So anvi'o will drop the ones that are actually in the "
                             f"database, and will take credit also for those that were not in the database in the first place.")
        elif not len(sources_to_drop):
            self.run.warning(f"Someone is trying to drop functions ({', '.join(sources_not_known)}) that do not occur in their "
                             f"contigs database :/ Anvi'o will pretend that this didn't happen.")
            return

        self.run.warning(f"Dropping {P('functional annotation source', len(sources_to_drop))} yo: {', '.join(sources_to_drop)}.")

        # remove the functions
        sources = ','.join([f'"{s}"' for s in sources_to_drop])
        database._exec(f'''DELETE FROM {t.gene_function_calls_table_name} WHERE source IN ({sources})''')

        # update the self table
        database.remove_meta_key_value_pair('gene_function_sources')
        database.set_meta_value('gene_function_sources', ','.join(sources_to_remain))


    def add_new_sources_to_functional_sources(self, gene_function_sources, database, drop_previous_annotations_first=False):
        """Update the self table with new functional annotation sources."""

        # are there any previous annotations in the db:
        gene_function_sources_in_db = self.get_function_sources_in_db(database)

        # difference between sources in the db, and incoming sources:
        gene_function_sources_both_in_db_and_incoming_dict = gene_function_sources.intersection(gene_function_sources_in_db)

        # here we will do some magic. there are mulitple scenarios to consider here based on whether there
        # are functions already in the database, whether some of them matches to the incoming functions, etc.
        # let's go case-by-case:
        if not gene_function_sources_in_db:
            # there are no funcions in this database. set the self table value to those that are
            # coming in
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources)))

        elif gene_function_sources_in_db and drop_previous_annotations_first:
            # there are gene calls, but the user wants everything to be dropeped first
            # FIXME: this is an artifact from times where we didn't have an `anvi-delete-functions`
            # program. there should be no reason in which the user should want to drop all functions
            # from a contigs database. the parameter should be removed, and this clause should be
            # removed..
            self.run.warning("As per your request, anvi'o is going to drop ALL FUNCTIONS IN THE DATABASE")

            # drop all previous annotations
            self.drop_functions(database, sources_to_drop=gene_function_sources_in_db)

            # everything is gone. set the self value to incoming functions.
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources)))

        elif gene_function_sources_in_db and gene_function_sources_both_in_db_and_incoming_dict:
            # some of the functions in the incoming dict match to what is already in the db. remove
            if len(gene_function_sources_both_in_db_and_incoming_dict) == 1:
                self.run.warning(f"A functional annotation source you wish to add to the database ({list(gene_function_sources_both_in_db_and_incoming_dict)[0]}) "
                                 f"is already in the database. Anvi'o will first drop the existing one so the incoming annotation could REPLACE it.")
            else:
                self.run.warning(f"Some functional annotation sources you wish to add to the database ({', '.join(gene_function_sources_both_in_db_and_incoming_dict)}) "
                                 f"are already in the database. Anvi'o will first drop them so the incoming annotations could REPLACE them.")

            # remove those entries for matching sources:
            self.drop_functions(database, gene_function_sources_both_in_db_and_incoming_dict)

            # set the sources
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources_in_db.union(gene_function_sources))))

        else:
            # fuctions in the db, but none of them match with the incoming annotation sources. totally new stuff.
            # good then. update sources
            database.remove_meta_key_value_pair('gene_function_sources')
            database.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources_in_db.union(gene_function_sources))))


    def create(self, functions_dict, drop_previous_annotations_first=False):
        self.sanity_check()

        # open connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # Add the new sources to existing sources
        gene_function_sources = set([v['source'] for v in list(functions_dict.values())])
        self.add_new_sources_to_functional_sources(gene_function_sources, database, drop_previous_annotations_first=drop_previous_annotations_first)

        unique_num_genes = len(set([v['gene_callers_id'] for v in list(functions_dict.values())]))

        # push the data
        db_entries = [tuple([functions_dict[v][h] for h in t.gene_function_calls_table_structure]) for v in functions_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.gene_function_calls_table_name, db_entries)

        # disconnect like a pro.
        database.disconnect()

        sources_string = ", ".join(gene_function_sources)
        self.run.info('Gene functions', f"{P('function call', len(functions_dict))} from {P('source', len(gene_function_sources))} ({sources_string}) "
                                        f"for {P('unique gene call', unique_num_genes)} have been added to the contigs database.")


    def sanity_check(self):
        pass
