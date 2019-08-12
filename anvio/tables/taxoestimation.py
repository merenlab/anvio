# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import hashlib
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


taxonomy_estimation_bin_name             = 'taxonomy_estimation_bin'
taxonomy_estimation_bin_structure        = ['entry_id', 'collection_name', 'bin_name', 'source'  , 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxonomy_estimation_bin_types            = [ 'numeric',   'text'   ,        'text'  ,  'text',      'text',   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

taxonomy_estimation_metagenome_name      = 'taxonomy_estimation_metagenome'
taxonomy_estimation_metagenome_structure = ['gene_caller_id',      'gene_name',  'source' , 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxonomy_estimation_metagenome_types     = [ 'numeric',             'text'    ,  'text',      'text'  ,   'text'  ,  'text'  ,  'text'   ,  'text'  ,  'text'  ,   'text'   ]


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

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


    def alignment_result_to_congigs(self,diamond_output,taxonomy_dict,match_id):

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        taxonomy_dictonnary=self.database.get_table_as_dict(t.taxon_names_table_name)

        list_taxo=[]
        entries=[]

        for result in diamond_output :
            SCG=result[0]
            for line_hit_to_split in result[1].split('\n')[1:-2]:
                if not line_hit_to_split.startswith('Query') or len(line_hit_to_split):
                    line_hit=line_hit_to_split.split('\t')

                    try:
                        entries+=[tuple([match_id,line_hit[0],SCG,line_hit[1],line_hit[2],line_hit[11]])]
                    except:
                        print("error parsing output aligment: %s" % (' '.join(line_hit)))
                        continue



                match_id+=1

                if line_hit[1] not in list_taxo and line_hit[1] not in taxonomy_dictonnary:
                    list_taxo+=[line_hit[1]]

        taxo_entries=[tuple([t_name_id]+list(taxonomy_dict[t_name_id].values())) for t_name_id in list_taxo]
        self.database.insert_many(t.blast_hits_table_name, entries)
        self.database.insert_many(t.taxon_names_table_name, taxo_entries)

        self.database.disconnect()
        return match_id

    def taxonomy_estimation_to_congis(self,possibles_taxonomy):

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        self.database.insert_many(t.taxonomy_estimation_metagenome_name, possibles_taxonomy)
        self.database.disconnect()

    def taxonomy_estimation_to_profile(self,possibles_taxonomy):

        self.bin_database = db.DB(self.profile_db_path, utils.get_required_version_for_db(self.profile_db_path))
        self.bin_database.insert_many(t.taxonomy_estimation_bin_name, possibles_taxonomy)
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

        try:
            dic_blast_hits=self.database.get_table_as_dict(t.blast_hits_table_name)
            taxonomy_dict =OrderedDict(self.database.get_table_as_dict(t.taxon_names_table_name))
        except:
            raise ConfigError("Anvi'o could not find the data for the taxonomic estimation,\
                               you should try to run 'anvi-diamond-for-taxonomy'")

        taxon_id_missing_taxonomy=[]
        for blast_hit in dic_blast_hits.values():
            if blast_hit['taxon_id'] not in taxonomy_dict:
                taxon_id_missing_taxonomy+=blast_hit['taxon_id']
        if len(taxon_id_missing_taxonomy):
            print(taxon_id_missing_taxonomy)
            raise ConfigError("it seams , some taxon id of blast hit are not linked to any phylogenie,\
                               you should run 'anvi-diamond-for-taxonomy'.%d" % (','.join(taxon_id_missing_taxonomy)))

        self.database.disconnect()
        return(dic_blast_hits,taxonomy_dict)
