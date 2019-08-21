# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import hashlib
import traceback
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
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


blast_hits_table_name                    = 'blast_hits'
blast_hits_table_structure               = ['match_id' , 'gene_callers_id', 'gene_name', 'taxon_id', 'pourcentage_identity', 'bitscore']
blast_hits_table_types                   = ['text'     ,       'text'    ,      'text'   ,   'text'   ,     'text'   ,         'text']


collection_taxonomy_estimation_name             = 'collection_taxonomy_estimation'
collection_taxonomy_estimation_structure        = ['entry_id', 'collection_name', 'bin_name', 'source'  , 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
collection_taxonomy_estimation_types            = [ 'numeric',   'text'   ,        'text'  ,  'text',      'text',   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

scg_taxonomy_estimation_name      = 'scg_taxonomy_estimation'
scg_taxonomy_estimation_structure = ['gene_caller_id',      'gene_name',  'source' , 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
scg_taxonomy_estimation_types     = [ 'numeric',             'text'    ,  'text',      'text'  ,   'text'  ,  'text'  ,  'text'   ,  'text'  ,  'text'  ,   'text'   ]


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


    def alignment_result_to_congigs(self,diamond_output):

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        entries=[]


        for result in diamond_output :
            estimation_id=0
            SCG=result[1]
            gene_callers_id=result[0]
            if not len(result[3]):
                continue
            entries+=[tuple([gene_callers_id, SCG, "Anvio", estimation_id, " "]+list(result[2].values()))]
            for consider_taxonomy in result[3]:
                estimation_id+=1
                if len(list(consider_taxonomy["taxonomy"].values())):
                    entries+=[tuple([gene_callers_id, SCG, "GTDB", consider_taxonomy["code"], consider_taxonomy["bestident"]]+ list(consider_taxonomy["taxonomy"].values()))]

        self.database.insert_many(t.scg_taxonomy_estimation_name, entries)
        self.database.disconnect()

    def taxonomy_estimation_to_congis(self,possibles_taxonomy):

        try:
            self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
            self.database.insert_many(t.scg_taxonomy_estimation_name, possibles_taxonomy)
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

    @timer
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

        try:
            dic_blast_hits=self.database.get_table_as_list_of_tuples(t.scg_taxonomy_estimation_name)
        except:
            traceback.print_exc()
            raise ConfigError("Anvi'o could not find the data for the taxonomic estimation,\
                               you should try to run 'anvi-diamond-for-taxonomy'")

        print(dic_blast_hits)
        self.database.disconnect()
        return(dic_blast_hits)
