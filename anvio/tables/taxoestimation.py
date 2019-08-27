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
__status__ = "Development"




run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

def timer(function):
    import time

    def timed_function(*args, **kwargs):
        n = 1
        start = time.time()
        for i in range(n):
            x = function(*args, **kwargs)
        end = time.time()
        print('Average time per call over {} calls for function \'{}\': {:6f} seconds'.format(
            n, function.__name__, (end - start) / n))
        return x
    return timed_function



class TablesForTaxoestimation(Table):
    def __init__(self, db_path, run=run, progress=progress,profile_db_path=False):

        self.db_path = db_path
        self.run = run
        self.progress = progress
        self.profile_db_path = profile_db_path

        utils.is_contigs_db(self.db_path)

        if profile_db_path:
            utils.is_profile_db(self.profile_db_path)
            self.profile_db_path = profile_db_path

            utils.is_profile_db_and_contigs_db_compatible(
                self.profile_db_path, self.db_path)


        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress)


    def alignment_result_to_congigs(self, table_index, diamond_output, source="GTDB"):
        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        entries=[]
        #FIXME source problem if to much hit
        for result in diamond_output :
            table_index+=1
            SCG=result[1]
            gene_callers_id=result[0]
            if not len(result[3]):
                continue
            if len(result[3]) < 5:
                accession=self.get_accession(result[2])
                entries+=[tuple([table_index, gene_callers_id, SCG, "Consensus", accession, result[3][0]["bestident"]]+list(result[2].values()))]
                for consider_taxonomy in result[3]:
                    table_index+=1
                    entries+=[tuple([table_index, gene_callers_id, SCG, source, consider_taxonomy["accession"], consider_taxonomy["bestident"]]+ list(consider_taxonomy["taxonomy"].values()))]
            else:
                accession=self.get_accession(result[2])
                entries+=[tuple([table_index, gene_callers_id, SCG, source+"_simplified", accession, result[3][0]["bestident"]]+list(result[2].values()))]
        
        self.database.insert_many(t.scg_taxonomy_table_name, entries)
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

    def get_dic_id_bin(self,args):

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
            raise ConfigError("Anvi'o could not find the data for the taxonomic estimation,\
                               you should try to run 'anvi-diamond-for-taxonomy'")
        else:
            return(dictonnary_taxonomy_by_index)
