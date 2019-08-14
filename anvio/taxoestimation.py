#!/usr/bin/env python3
# -*- coding: utf-8

import argparse
import os
import re
import sys
import time
import copy
import traceback


import multiprocessing
from collections import Counter

import anvio
import pickle
from anvio.drivers.diamond import Diamond
import pandas as pd
from collections import OrderedDict

from tabulate import tabulate
import anvio.db as db
import pandas as pd
import numpy as np

import anvio.ccollections as ccollections
import anvio.fastalib as f
import anvio.filesnpaths as filesnpaths
import anvio.hmmops as hmmops
import anvio.hmmopswrapper as hmmopswrapper
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.tables.tableops import Table
from anvio.tables.taxoestimation import TablesForTaxoestimation
import anvio.tables as t
from anvio.dbops import ContigsSuperclass
from anvio.drivers import Aligners, driver_modules
from anvio.errors import ConfigError, FilesNPathsError


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
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)
aligners = Aligners()

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


class SCGsdiamond:
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.taxonomy_file_path = A('taxonomy_file')
        self.taxonomy_database_path = A('taxonomy_database')
        self.write_buffer_size = int(A('write_buffer_size') if A('write_buffer_size') is not None else 500)
        self.db_path=A('contigs_db')
        self.core=A('num_threads')
        self.num_process=A('contigs_db')

        self.max_target_seqs=20
        self.evalue=1e-05
        self.min_pct_id=90


        self.num_process=args.num_process

        if not self.core:
            self.core="1"

        if not args.num_process:
            self.num_process="2"

        self.initialized = False

        self.metagenome=False

        if args.metagenome:
            self.metagenome=True

        self.source="unknow"




        if not self.taxonomy_file_path:
            self.taxonomy_file_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG/mergedb/matching_taxonomy.tsv')
            self.source="https://gtdb.ecogenomic.org/"

        if not self.taxonomy_database_path:
            self.taxonomy_database_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG/mergedb/species/')



        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))


        self.SCGs = [db for db in os.listdir(
            self.taxonomy_database_path) if db.endswith(".dmnd")]

        self.taxonomic_levels_parser = {'d': 't_domain',
                                        'p': "t_phylum",
                                        'c': "t_class",
                                        'o': "t_order",
                                        'f': "t_family",
                                        'g': "t_genus",
                                        's': "t_species"}

        self.taxonomic_levels = [
            't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]

        self.SCG_DB_PATH = lambda SCG: os.path.join(
            self.taxonomy_database_path, SCG)

        self.SCG_FASTA_DB_PATH = lambda SCG: os.path.join(self.taxonomy_database_path,
                                                          [db for db in os.listdir(self.taxonomy_database_path) if db.endwith(".dmnd")])

        self.sanity_check()

        self.taxonomy_dict=OrderedDict()


    def sanity_check(self):
        if not filesnpaths.is_file_exists(self.taxonomy_file_path, dont_raise=True):
            raise ConfigError("Anvi'o could not find taxonomy file '%s'. You must declare one before continue."\
                               % self.taxonomy_file_path)
        filesnpaths.is_file_exists(self.taxonomy_database_path)

        if not len(self.SCGs):
            raise ConfigError(
                "This class can't be used with out a list of single-copy core genes.")

        if not len(self.SCGs) == len(set(self.SCGs)):
            raise ConfigError("Each member of the list of SCGs you wish to use with this class must\
                               be unique and yes, you guessed right. You have some repeated gene\
                               names.")

        SCGs_missing_databases = [
            SCG for SCG in self.SCGs if not filesnpaths.is_file_exists(self.SCG_DB_PATH(SCG))]
        if len(SCGs_missing_databases):
            raise ConfigError("Even though anvi'o found the directory for databases for taxonomy stuff,\
                               your setup seems to be missing %d databases required for everything to work\
                               with the current genes configuration of this class. Here are the list of\
                               genes for which we are missing databases: '%s'." % (', '.join(missing_databases)))


    def init(self):
        if self.initialized:
            return

        # initialize taxonomy dict. we should do this through a database in the long run.
        with open(self.taxonomy_file_path, 'r') as taxonomy_file:
            self.progress.new("Loading taxonomy file")
            for accession, taxonomy_text in [l.strip('\n').split('\t') for l in taxonomy_file.readlines() if not l.startswith('#') and l]:
                # taxonomy_text kinda looks like this:
                #
                #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter sp1
                #
                d = {}
                for token, taxon in [e.split('__', 1) for e in taxonomy_text.split(';')]:
                    if token in self.taxonomic_levels_parser:
                        d[self.taxonomic_levels_parser[token]] = taxon

                self.taxonomy_dict[accession] = d

                self.progress.increment()

        self.progress.end()
        self.initialized = True


    def get_hmm_sequences_dict_into_type_multi(self, hmm_sequences_dict):
        hmm_sequences_dict_per_type = {}

        for entry_id in hmm_sequences_dict:
            entry = hmm_sequences_dict[entry_id]

            name = entry['gene_name']

            if name in hmm_sequences_dict_per_type:
                hmm_sequences_dict_per_type[name][entry_id] = entry
            else:
                hmm_sequences_dict_per_type[name] = {entry_id: entry}

        return hmm_sequences_dict_per_type



    def predict_from_SCGs_dict_multiseq(self, hmm_sequences_dict):
        """Takes an HMMs dictionary, and yields predictions"""

        self.init()

        hmm_sequences_dict_per_type = self.get_hmm_sequences_dict_into_type_multi(
            hmm_sequences_dict)

        num_listeprocess = len(hmm_sequences_dict_per_type)



        aligners="Diamond"

        #self.run.warning('', header='Aligment for %s ' % self.db_path, lc='green')
        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Taxonomy', self.taxonomy_file_path)
        self.run.info('Database reference', self.taxonomy_database_path)
        self.run.info('Source', self.source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('Number of SCGs', len(hmm_sequences_dict_per_type))
        self.run.info('SCGs', ','.join(list(hmm_sequences_dict_per_type.keys())))


        self.run.warning('', header='Parameters for aligments with %s for taxonomy' % aligners, lc='green')
        self.run.info('Blast type', "Blastp")
        self.run.info('Maximum number of target sequences', self.max_target_seqs)
        self.run.info('Minimum bit score to report alignments', self.min_pct_id)
        self.run.info('Number aligment running in same time', self.num_process)
        self.run.info('Number of CPUs will be used for each aligment',self.core)

        """tmp_dir = filesnpaths.get_temp_directory_path()

        self.hmm_scan_output = os.path.join(tmp_dir, 'hmm.output')
        self.hmm_scan_hits = os.path.join(tmp_dir, 'hmm.hits')"""



        self.tables_for_taxonomy = TablesForTaxoestimation(self.db_path, run, progress)
        self.tables_for_taxonomy.delete_contents_of_table(t.blast_hits_table_name)

        self.progress.new('Computing SCGs aligments', progress_total_items=num_listeprocess)
        self.progress.update('Initializing %d process...' % int(self.num_process))

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        diamond_output=[]
        match_id=0

        for SCG in hmm_sequences_dict_per_type:

            sequence=""
            for entry in hmm_sequences_dict_per_type[SCG].values():
                if 'sequence' not in entry or 'gene_name' not in entry:
                    raise ConfigError("The `get_filtered_dict` function got a parameter that\
                                       does not look like the way we expected it. This function\
                                       expects a dictionary that contains keys `gene_name` and `sequence`.")

                sequence = sequence+">"+str(entry['gene_callers_id'])+"\n"+entry['sequence']+"\n"
                entry['hits']=[]
            input_queue.put([SCG,sequence])

        workers = []
        for i in range(0, int(self.num_process)):

            worker = multiprocessing.Process(target=self.get_raw_blast_hits_multi,
                args=(input_queue,output_queue))

            workers.append(worker)
            worker.start()

        finish_process = 0
        while finish_process < num_listeprocess:
            try:
                diamond_output += [output_queue.get()]


                if self.write_buffer_size > 0 and len(diamond_output) % self.write_buffer_size == 0:

                    match_id=self.tables_for_taxonomy.alignment_result_to_congigs(diamond_output,self.taxonomy_dict,match_id)
                    diamond_output=[]

                finish_process += 1

                self.progress.increment(increment_to=finish_process)
                progress.update("Processed %s of %s SGCs aligment in %s threads with %s cores." % (finish_process, num_listeprocess,int(self.num_process),self.core))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()


        match_id=self.tables_for_taxonomy.alignment_result_to_congigs(diamond_output,self.taxonomy_dict,match_id)
        progress.end()
        self.run.info('Number of match selected',match_id)


    def get_raw_blast_hits_multi(self, input_queue,output_queue):

        while True:
            d = input_queue.get(True)
            db_path = self.SCG_DB_PATH(d[0])
            diamond = Diamond(db_path,run=run_quiet, progress= progress_quiet)
            diamond.max_target_seqs = self.max_target_seqs
            diamond.evalue = self.evalue
            diamond.min_pct_id = self.min_pct_id
            diamond.num_process = self.core

            diamond_output = diamond.blastp_stdin_multi(d[1])

            output_queue.put([d[0],diamond_output])


class SCGsTaxomy:

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.db_path=A('contigs_db')
        self.profile_db_path=A('profile_db')
        self.output_file_path=A('output_file')

        self.collection_name=args.collection_name

        self.initialized = False

        self.bin_id=args.bin_id

        self.metagenome=False




        if args.metagenome:
            self.metagenome=True

        self.cut_off_methode = args.cut_off_methode

        self.methode = args.methode

        self.hits_per_gene={}

        self.pident_level_path=os.path.join(os.path.dirname(
            anvio.__file__), 'data/misc/SCG/mergedb/dico_low_ident.pickle')

        with open(self.pident_level_path, 'rb') as handle:
            self.dicolevel = pickle.load(handle)

        self.taxonomic_levels_parser = {'t_domain' : 'd__',
                                        "t_phylum" : 'p__',
                                        "t_class" : 'c__',
                                        "t_order" : 'o__',
                                        "t_family" : 'f__',
                                        "t_genus" : 'g__' ,
                                        "t_species" : 's__' }


    def sanity_check(self):
        filesnpaths.is_file_exists(self.db_path)
        filesnpaths.is_file_exists(self.taxonomy_database_path)

        utils.is_contigs_db(self.db_path)

        if self.profile_db_path:
            utils.is_profile_db(self.profile_db_path)

        if not len(self.SCGs):
            raise ConfigError(
                "This class can't be used with out a list of single-copy core genes.")

        if not len(self.SCGs) == len(set(self.SCGs)):
            raise ConfigError("Each member of the list of SCGs you wish to use with this class must\
                               be unique and yes, you guessed right. You have some repeated gene\
                               names.")

        SCGs_missing_databases = [
            SCG for SCG in self.SCGs if not filesnpaths.is_file_exists(self.SCG_DB_PATH(SCG))]
        if len(SCGs_missing_databases):
            raise ConfigError("Even though anvi'o found the directory for databases for taxonomy stuff,\
                               your setup seems to be missing %d databases required for everything to work\
                               with the current genes configuration of this class. Here are the list of\
                               genes for which we are missing databases: '%s'." % (', '.join(missing_databases)))


    def init(self):

        filesnpaths.is_output_file_writable(self.output_file_path)


        self.tables_for_taxonomy = TablesForTaxoestimation(self.db_path, run, progress, self.profile_db_path)

        self.dic_blast_hits,self.taxonomy_dict=self.tables_for_taxonomy.get_data_for_taxonomy_estimation()

        if not (self.dic_blast_hits):
            raise ConfigError("Anvi'o can't make a taxonomy estimation because aligment didn't return any match or you forgot to run 'anvi-diamond-for-taxonomy'.")
            self.run.info('taxonomy estimation not possible for:', self.db_path)


        if self.profile_db_path and not self.metagenome:
            self.identifier="Bin_id"
            self.entries_db_profile=[]


            self.dic_id_bin = self.tables_for_taxonomy.get_dic_id_bin(self.args)


        else:
            self.entries_db_contigs=[]
            self.identifier="Gene_id"


        self.progress.new('Association alimgent result by %s' % self.identifier, progress_total_items=len(self.dic_blast_hits))

        for query in self.dic_blast_hits.values():

            if self.profile_db_path and not self.metagenome:
                for bin_id,bin_gene_callers_id in self.dic_id_bin.items():
                    if int(query['gene_callers_id']) in bin_gene_callers_id:
                        var=bin_id
                        break
                    else:
                        var=None

                if var==None:
                    continue

            else:
                var=query['gene_callers_id']

            hit=[{'accession':query['taxon_id'], 'pident':float(query['pourcentage_identity']), 'bitscore': float(query['bitscore'])}]

            if var not in self.hits_per_gene:
                self.hits_per_gene[var]={}
            if query['gene_name'] not in self.hits_per_gene[var]:
                self.hits_per_gene[var][query['gene_name']]=[]

            self.hits_per_gene[var][query['gene_name']] = self.hits_per_gene[var][query['gene_name']] + hit
            self.progress.increment()
        self.progress.end()



        self.initialized = True


    def estimate_taxonomy(self,source="GTDB"):

        self.init()


        self.run.warning('', header='Taxonomy estimation for %s' % self.db_path, lc='green')
        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Source', source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('Taxonomy assignment for', self.identifier)
        self.run.info('output file for taxonomy', self.output_file_path)


        possibles_taxonomy=[]
        entry_id=0
        stdout_taxonomy=[]
        possibles_taxonomy.append([self.identifier,'domain','phylum','class','order','family','genus','species'])

        for name, SCGs_hit_per_gene in self.hits_per_gene.items():
            #print(SCGs_hit_per_gene)

            taxonomy = self.get_consensus_taxonomy(
                SCGs_hit_per_gene, name)

            if not taxonomy:
                self.run.info('taxonomy estimation not possible for:', name)
                possibles_taxonomy.append([name]+["NA"]*7)
                continue

            if self.metagenome or not self.profile_db_path:
                if str(list(taxonomy.values())[-1]) not in stdout_taxonomy and str(list(taxonomy.values())[-1]) :
                    stdout_taxonomy.append(str(list(taxonomy.values())[-1]))

                possibles_taxonomy.append([name]+list(taxonomy.values()))

                self.entries_db_contigs+=[(tuple([name,list(SCGs_hit_per_gene.keys())[0],source]+list(taxonomy.values())))]
                #output_taxonomy+=name+'\t'+'\t'.join(list(taxonomy.values()))+'\n'

            if self.profile_db_path and not self.metagenome:

                possibles_taxonomy.append([name]+list(taxonomy.values()))


                self.run.info('Bin name',
                              name, nl_before=1)
                self.run.info('estimate taxonomy',
                              '/'.join(list(taxonomy.values())))

                self.entries_db_profile+=[(tuple([entry_id,self.collection_name,name,source]+list(taxonomy.values())))]
                entry_id+=1
                #output_taxonomy+=name+'\t'+'\t'.join(list(taxonomy.values()))+'\n'


        if self.metagenome or not self.profile_db_path:
            if len(possibles_taxonomy):
                self.run.info('Possible presence ','|'.join(list(stdout_taxonomy)))
                self.tables_for_taxonomy.taxonomy_estimation_to_congis(self.entries_db_contigs)


        if self.profile_db_path and not self.metagenome:
                self.tables_for_taxonomy.taxonomy_estimation_to_profile(self.entries_db_profile)



        try:
            with open(self.output_file_path,"w") as output:
                output_taxonomy=['\t'.join(line) for line in possibles_taxonomy]
                output.write('\n'.join(output_taxonomy))
        except:
            self.run.warning(traceback.print_exc(), header='Anvi\'o to creat the file for taxonomy result %s' % self.output_file_path, lc="red")


        self.show_taxonomy_estimation(possibles_taxonomy)


    def show_taxonomy_estimation(self, possibles_taxonomy):
        self.run.warning(None, header='Taxonomy estimation' , lc="yellow")
        possibles_taxonomy_dataframe=pd.DataFrame(possibles_taxonomy,columns=possibles_taxonomy[0])
        possibles_taxonomy_dataframe.set_index(self.identifier, inplace=True)
        possibles_taxonomy_dataframe=possibles_taxonomy_dataframe.sort_values(by=['domain','phylum','class','order','family','genus','species'],ascending=False)
        print(tabulate(possibles_taxonomy_dataframe, headers="firstrow",
                       tablefmt="fancy_grid", numalign="right"))


    def show_hits(self, name, gene_name, hits):
        self.run.warning(None, header='%s / %s' %
                         (name, gene_name), lc="green")
        header = ['%id', 'bitscore', 'taxonomy']
        table = []

        for hit in hits:
            table.append([hit['pident'], hit['bitscore'],
                          ' / '.join(self.taxonomy_dict[hit['accession']].values())])

        print(tabulate(table, headers=header,
                       tablefmt="fancy_grid", numalign="right"))


    def show_table_score(self, name, selected_entrys_by_score):
        self.run.warning(None, header='%s' % (name), lc="yellow")
        header = ['Average bitscore', 'taxonomy']
        table = []

        for code, score in sorted(selected_entrys_by_score.items(), key=lambda x: (-x[1], x[0])):
            table.append(
                [score, ' / '.join(self.taxonomy_dict[code].values())])

        print(tabulate(table, headers=header,
                       tablefmt="fancy_grid", numalign="right"))


    def show_matrix_rank_orginal(self, name, matrix, list_position_entry, list_position_ribosomal):
        show_matrix = [sublist[:6] for sublist in matrix]
        show_list_position_ribosomal = list_position_ribosomal[:6]
        header = show_list_position_ribosomal
        table = []
        i = 0

        for individue in show_matrix:
            taxonomyindividue = list(
                self.taxonomy_dict[list_position_entry[i]].values())
            line = [taxonomyindividue[-1]] + individue
            table.append(line)
            i += 1

        self.run.warning(None, header='%s' % (name), lc="blue")
        print(tabulate(table, headers=header,
                       tablefmt="fancy_grid", numalign="right"))
        if len(show_list_position_ribosomal[6:]):
            show_matrix = [sublist[6:] for sublist in show_matrix]
            show_list_position_ribosomal = show_list_position_ribosomal[6:]
            self.show_matrix_rank(name, matrix, show_matrix,
                                  show_list_position_ribosomal)

    def show_matrix_rank(self, name, matrix, list_position_entry, list_position_ribosomal):
        headers =list(matrix.keys())

        self.run.warning(None, header='%s' % (name), lc="blue")
        print(tabulate(matrix, headers=headers,
                       tablefmt="fancy_grid", numalign="right"))


    def get_consensus_taxonomy(self, SCGs_hit_per_gene, name):
        """Different methode for assignation"""

        consensus_taxonomy=self.solo_hits(SCGs_hit_per_gene)
        if consensus_taxonomy:

            return(consensus_taxonomy)

        else:

            if self.methode == "friedman":
                consensus_taxonomy = self.rank_assignement(SCGs_hit_per_gene, name)
            return(consensus_taxonomy)

    def get_matching_gene(self, SCGs_hit_per_gene):
        matching_genes = [
            gene for gene in SCGs_hit_per_gene if len(SCGs_hit_per_gene[gene])]
        number_of_matching_genes = len(matching_genes)
        return matching_genes

    def get_consensus_taxonomy_with_score_by_entry(self, score_by_entry, name, cut_off_methode):
        try:
            maxscore = max(score_by_entry.values())
        except:
            self.run.info("Estimate Taxonomy of " + str(name), "N/A")
            return
        selected_entrys_by_score = {
            code: score for code, score in score_by_entry.items() if score > (float(maxscore) * float(cut_off_methode))}

        self.run.info("Number of taxonomy use for the consensus",
                      len(selected_entrys_by_score))

        if anvio.DEBUG:
            self.show_table_score(name, selected_entrys_by_score)

        taxonomy = []
        for code, score in sorted(selected_entrys_by_score.items()):
            taxonomy.append(self.taxonomy_dict[code])

        self.assign_taxonomie_solo_hit(taxonomy)



    def fill_matrix(self, name, emptymatrix_ident, SCGs_hit_per_gene, list_position_entry,
                    list_position_ribosomal, matchinggenes):

        for SCG in list_position_ribosomal:
            for entry in SCGs_hit_per_gene[SCG]:
                if entry['accession'] in list_position_entry:
                    emptymatrix_ident.at[entry['accession'], SCG]= float(entry['pident'])

        return(emptymatrix_ident)


    def make_liste_individue(self, SCGs_hit_per_gene, matchinggenes):
        liste_individue = []
        liste_ribo= []
        for SCG in matchinggenes:
            if SCG not in liste_ribo:
                liste_ribo+=[SCG]
            for entry in SCGs_hit_per_gene[SCG]:
                if entry['accession'] not in liste_individue:
                    liste_individue+=[entry['accession']]
        return(liste_individue,liste_ribo)


    def make_rank_matrix(self, name, SCGs_hit_per_gene):
        matchinggenes = self.get_matching_gene(SCGs_hit_per_gene)

        list_position_entry, list_position_ribosomal=self.make_liste_individue( SCGs_hit_per_gene, matchinggenes)


        emptymatrix = pd.DataFrame(columns=list_position_ribosomal, index=list_position_entry)

        matrix_pident= self.fill_matrix(name, emptymatrix, SCGs_hit_per_gene, list_position_entry,
                                                      list_position_ribosomal, matchinggenes)

        if anvio.DEBUG:
            self.show_matrix_rank(
                name, matrix_pident, list_position_entry, list_position_ribosomal)

        return(matrix_pident)


    def rank_assignement(self, SCGs_hit_per_gene, name):
        try:
            matrix_pident = self.make_rank_matrix(name, SCGs_hit_per_gene)

            taxonomy= self.make_list_taxonomy(matrix_pident)

            assignation_reduce = self.reduce_assignation_level(taxonomy)
            assignation = self.assign_taxonomie_solo_hit(assignation_reduce)

        except:
            traceback.print_exc()
            self.run.warning(SCGs_hit_per_gene, header='Fail matrix')
            assignation=[]

        finally:
            return(assignation)

    def assign_taxonomie_solo_hit(self, taxonomy):
        if not taxonomy or not taxonomy[0]:
            return 0
        assignation = taxonomy[0]
        for taxon in taxonomy[1:]:
            for level in taxon:
                if taxon[level] not in assignation.values():
                    assignation[level]=''
        return(assignation)


    def make_list_taxonomy(self, matrix_pident):
        taxonomy = []
        matrix_rank = matrix_pident.rank(method='min', ascending=False, na_option='bottom')
        series_sum_rank=matrix_rank.sum(axis=1)
        minimum_rank=series_sum_rank.min()
        topseries=series_sum_rank.loc[series_sum_rank[:] == minimum_rank]
        for individue, row in topseries.items():
            bestSCG = pd.to_numeric(matrix_pident.loc[individue,:]).idxmax()
            bestident = matrix_pident.loc[individue, bestSCG]
            taxonomy.append({"bestSCG" : bestSCG,"bestident" : bestident,"taxo" : OrderedDict(self.taxonomy_dict[individue])})
        return(taxonomy)


    def reduce_assignation_level(self,taxonomy):
        reduce_taxonomy=[]
        for possibilitie in taxonomy:
            for level,value in reversed(possibilitie["taxo"].items()):
                if possibilitie["bestident"] < float(self.dicolevel[possibilitie["bestSCG"]][self.taxonomic_levels_parser[level]+value]):
                    possibilitie["taxo"][level]=''
                else:
                    break
            reduce_taxonomy.append(possibilitie["taxo"])

        return(reduce_taxonomy)


    def solo_hits(self, SCGs_hit_per_gene,cutoff_solo_hit=0.90):
        taxonomy = []
        for entry in SCGs_hit_per_gene.values():
            if len(entry) == 1 and entry[0]['pident'] > cutoff_solo_hit:
                taxonomy.append(self.taxonomy_dict[entry[0]['accession']])
        if taxonomy:
            assignation=self.assign_taxonomie_solo_hit(taxonomy)
            return assignation
        else:
            return False
