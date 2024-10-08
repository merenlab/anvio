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
__maintainer__ = "Quentin Clayssen"
__email__ = "quentin.clayssen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TableForSCGTaxonomy(Table):
    def __init__(self, db_path, run=run, progress=progress, profile_db_path=False):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        utils.is_contigs_db(self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress)


    def add(self, blastp_search_output):
        """Incrementally adds new hits to a contigs database.

           It is essential to run the member function `update_db_self_table_values` once adding new hits are complete.
           At the time of writing this class w couldn't find a better way to do it.
        """

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        entries=[]
        for gene_callers_id, scg_name, scg_hits in blastp_search_output:
            # go back if there is nothing to do
            if not len(scg_hits):
                continue

            for scg_hit in scg_hits:
                entries.append([gene_callers_id, scg_name] + [scg_hit[f] for f in t.scg_taxonomy_table_structure[2:]])

        self.database.insert_many(t.scg_taxonomy_table_name, entries)
        self.database.disconnect()


    def update_db_self_table_values(self, taxonomy_was_run=False, database_version=None):
        """Updates the self table in contigs db.

        The purpose of this function is to clarify whether scg taxonomy was run for a contigs
        database, and if yes, which version of the local database was used to keep track of
        versions.

        Paremeters
        ==========
        taxonomy_was_run: bool, False
            Set True if taxonomy was run successfuly.
        database_version: str, None
            This sould be read from the ctx.target_database_release in taxonomyops.
        """

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        self.database.update_meta_value("scg_taxonomy_was_run", taxonomy_was_run)
        self.database.update_meta_value("scg_taxonomy_database_version", database_version)
        self.database.disconnect()


    def get_accession(self,taxonomy):
        for level, taxon in reversed(list(taxonomy.items())):
            if taxon == "NA" :
                continue
            code=abs(hash(level+taxon)) % (10 ** 8)
            accession=taxon+"_"+str(code)
            break

        return(accession)


    def get_data_for_taxonomy_estimation(self):
        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        #FIXME Argument for select return genes

        dictonnary_taxonomy_by_index=self.database.get_table_as_dict(t.scg_taxonomy_table_name)
        self.database.disconnect()
        if not len(dictonnary_taxonomy_by_index):
            raise ConfigError("Your contigs database does not seem to contain any information anvi'o can use to "
                              "estimate taxonomy of anything. Please try running the program 'anvi-run-scg-taxonomy' "
                              "first.")
        else:
            return(dictonnary_taxonomy_by_index)
