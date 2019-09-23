# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import hashlib
import traceback
import random
import string
from collections import OrderedDict

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

import anvio.hmmops as hmmops

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.tables.tableops import Table
from anvio.parsers import parser_modules
from anvio.dbops import ContigsSuperclass
from anvio.tables.genecalls import TablesForGeneCalls


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
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
        self.profile_db_path = profile_db_path

        utils.is_contigs_db(self.db_path)

        if profile_db_path:
            utils.is_profile_db(self.profile_db_path)
            self.profile_db_path = profile_db_path

            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress)

        self.set_next_available_id(t.scg_taxonomy_table_name)


    def add(self, blastp_search_output):
        """Incrementally adds new hits to a contigs database.

           It is essential to run the member functio `update_self_value` once adding new hits are complete.
           At the time of writing this class w couldn't find a better way to do it.
        """

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        entries=[]
        for gene_callers_id, scg_name, scg_hits in blastp_search_output:
            # go back if there is nothing to do
            if not len(scg_hits):
                continue

            for scg_hit in scg_hits:
                entries.append([self.next_id(t.scg_taxonomy_table_name), gene_callers_id, scg_name] + [scg_hit[f] for f in t.scg_taxonomy_structure[3:]])

        self.database.insert_many(t.scg_taxonomy_table_name, entries)
        self.database.disconnect()


    def update_self_value(self):
        """Updates the self table in contigs db to clarify that scg taxonomy were run"""

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        self.database.update_meta_value("scg_taxonomy_was_run", True)
        self.database.disconnect()


    def get_accession(self,taxonomy):
        for level, taxon in reversed(list(taxonomy.items())):
            if taxon == "NA" :
                continue
            code=abs(hash(level+taxon)) % (10 ** 8)
            accession=taxon+"_"+str(code)
            break

        return(accession)


    def taxonomy_estimation_to_congis(self,possibles_taxonomy):
        try:
            self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
            self.database.insert_many(t.scg_taxonomy_table_name, possibles_taxonomy)
        except:
            self.run.warning(traceback.print_exc(), header='Anvi\'o fail the enter the result in %s' % self.db_pat, lc="red")
        finally:
            self.database.disconnect()


    def taxonomy_estimation_to_profile(self,possibles_taxonomy):
        try:
            self.bin_database = db.DB(self.profile_db_path, utils.get_required_version_for_db(self.profile_db_path))
            self.bin_database.insert_many(t.collection_taxonomy_estimation_name, possibles_taxonomy)
        except:
            self.run.warning(traceback.print_exc(), header='Anvi\'o fail the enter the result in %s' % self.profile_db_path, lc="red")
        finally:
            self.bin_database.disconnect()


    def get_dic_id_bin(self, args):
        self.bin_database = db.DB(self.profile_db_path, utils.get_required_version_for_db(self.profile_db_path))
        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        dic_id_bin={}

        splits_dict = ccollections.GetSplitNamesInBins(args).get_dict()
        run.info('Init', '%d splits in %d bin(s)' % (
            sum([len(v) for v in list(splits_dict.values())]), len(splits_dict)))

        self.progress.new('Load HMM resulst')

        s = hmmops.SequencesForHMMHits(self.db_path)

        hits_in_splits, split_name_to_bin_id = s.get_hmm_hits_in_splits(splits_dict)

        dic_genes_in_splits=self.database.get_table_as_dict("genes_in_splits")

        self.progress.end()

        self.progress.new('Aligment result by Bin', progress_total_items=len(dic_genes_in_splits))

        for split in dic_genes_in_splits.values():
            if split['split'] in split_name_to_bin_id:
                if split_name_to_bin_id[split['split']] not in dic_id_bin:
                    dic_id_bin[split_name_to_bin_id[split['split']]]=[split['gene_callers_id']]
                else:
                    dic_id_bin[split_name_to_bin_id[split['split']]]+=[split['gene_callers_id']]
            self.progress.increment()
        self.progress.end()


        self.bin_database.disconnect()
        self.database.disconnect()

        return(dic_id_bin)


    def get_data_for_taxonomy_estimation(self):
        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        #FIXME Argument for select return genes

        dictonnary_taxonomy_by_index=self.database.get_table_as_dict(t.scg_taxonomy_table_name)
        self.database.disconnect()
        if not len(dictonnary_taxonomy_by_index):
            traceback.print_exc()
            raise ConfigError("Your contigs database does not seem to contain any information anvi'o can use to\
                               estimate taxonomy of anything. Please try running the program 'anvi-run-scg-taxonomy'\
                               first.")
        else:
            return(dictonnary_taxonomy_by_index)
