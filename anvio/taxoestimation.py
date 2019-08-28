#!/usr/bin/env python3
# -*- coding: utf-8

import argparse
import copy
import multiprocessing
import os
import pickle
import re
import sys
import time
import traceback
import tarfile
import shutil

from collections import Counter, OrderedDict

import anvio
import anvio.ccollections as ccollections
import anvio.db as db
import anvio.fastalib as f
import anvio.filesnpaths as filesnpaths
import anvio.hmmops as hmmops
import anvio.hmmopswrapper as hmmopswrapper
import anvio.tables as t
import anvio.terminal as terminal
import anvio.utils as utils
import numpy as np
import pandas as pd
from copy import deepcopy
from anvio.dbops import ContigsSuperclass
from anvio.drivers import Aligners, driver_modules
from anvio.drivers.diamond import Diamond
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.tableops import Table
from anvio.tables.taxoestimation import TablesForTaxoestimation
from tabulate import tabulate

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


class TaxonomyEstimation:

    def __init__(self, taxonomy_dict, dicolevel=None, reduction=False, run=run, progress=progress):

        
        if not dicolevel:
            self.pident_level_path = os.path.join(
                    os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY/GTDB/MIN_PCT_ID_PER_TAXONOMIC_LEVEL.pickle')

            with open(self.pident_level_path, 'rb') as handle:
                self.dicolevel = pickle.load(handle)

        if reduction:
            self.reduction=True
        else:
            self.reduction=False

        self.taxonomy_dict = taxonomy_dict

        self.taxonomic_levels_parser = {'t_domain': 'd__',
                                        "t_phylum": 'p__',
                                        "t_class": 'c__',
                                        "t_order": 'o__',
                                        "t_family": 'f__',
                                        "t_genus": 'g__',
                                        "t_species": 's__'}

        self.taxonomic_levels = [
            't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]

    def show_table_score(self, name, selected_entrys_by_score):
        self.run.warning(None, header='%s' % (name), lc="yellow")
        header = ['Average bitscore', 'taxonomy']
        table = []

        for code, score in sorted(selected_entrys_by_score.items(), key=lambda x: (-x[1], x[0])):
            table.append(
                [score, ' / '.join(self.taxonomy_dict[code].values())])

        print(tabulate(table, headers=header,
                       tablefmt="fancy_grid", numalign="right"))

    def show_matrix_rank(self, name, matrix, list_position_entry, list_position_ribosomal):
        headers = list(matrix.keys())

        self.run.warning(None, header='%s' % (name), lc="blue")
        print(tabulate(matrix, headers=headers,
                       tablefmt="fancy_grid", numalign="right"))

    def get_consensus_taxonomy(self, SCGs_hit_per_gene, name):
        """Different methode for assignation"""
        consensus_taxonomy = self.solo_hits(SCGs_hit_per_gene)
        if consensus_taxonomy:
            for SCG, entry in SCGs_hit_per_gene.items():
                taxonomy = [{"bestSCG": SCG, "bestident": entry[0]['pident'],
                            "accession": entry[0]['accession'], "taxonomy":OrderedDict(self.taxonomy_dict[entry[0]['accession']])}]

            return(consensus_taxonomy, taxonomy)

        else:
            consensus_taxonomy, taxonomy = self.rank_assignement(
                SCGs_hit_per_gene, name)
            return(consensus_taxonomy, taxonomy)

    def get_matching_gene(self, SCGs_hit_per_gene):
        matching_genes = [
            gene for gene in SCGs_hit_per_gene if len(SCGs_hit_per_gene[gene])]
        number_of_matching_genes = len(matching_genes)
        return matching_genes

    def fill_matrix(self, name, emptymatrix_ident, SCGs_hit_per_gene, list_position_entry,
                    list_position_ribosomal, matchinggenes):

        for SCG in list_position_ribosomal:
            for entry in SCGs_hit_per_gene[SCG]:
                if entry['accession'] in list_position_entry:
                    emptymatrix_ident.at[entry['accession'], SCG] = float(
                        entry['pident'])

        return(emptymatrix_ident)

    def make_liste_individue(self, SCGs_hit_per_gene, matchinggenes):
        liste_individue = []
        liste_ribo = []
        for SCG in matchinggenes:
            if SCG not in liste_ribo:
                liste_ribo += [SCG]
            for entry in SCGs_hit_per_gene[SCG]:
                if entry['accession'] not in liste_individue:
                    liste_individue += [entry['accession']]
        return(liste_individue, liste_ribo)

    def make_rank_matrix(self, name, SCGs_hit_per_gene):
        matchinggenes = self.get_matching_gene(SCGs_hit_per_gene)

        list_position_entry, list_position_ribosomal = self.make_liste_individue(
            SCGs_hit_per_gene, matchinggenes)

        emptymatrix = pd.DataFrame(
            columns=list_position_ribosomal, index=list_position_entry)

        matrix_pident = self.fill_matrix(name, emptymatrix, SCGs_hit_per_gene, list_position_entry,
                                         list_position_ribosomal, matchinggenes)

        if anvio.DEBUG:
            self.show_matrix_rank(
                name, matrix_pident, list_position_entry, list_position_ribosomal)

        return(matrix_pident)

    def rank_assignement(self, SCGs_hit_per_gene, name, reduction=False):
        try:
            matrix_pident = self.make_rank_matrix(name, SCGs_hit_per_gene)

            taxonomy = self.make_list_taxonomy(matrix_pident)
            taxonomy_to_reduction = deepcopy(taxonomy)
            assignation_reduce = self.reduce_assignation_level(taxonomy_to_reduction)
            assignation = self.assign_taxonomie_solo_hit(assignation_reduce)
        except:
            traceback.print_exc()
            self.run.warning(SCGs_hit_per_gene, header='Fail matrix')
            assignation = []

        finally:
            return(assignation, taxonomy)

    def assign_taxonomie_solo_hit(self, taxonomy):
        if not taxonomy or not taxonomy[0]:
            return 0
        assignation = taxonomy[0]
        if taxonomy[1:]:
            for taxon in taxonomy[1:]:
                for level in taxon:
                    #print(taxon[level],assignation.values())
                    if taxon[level] not in list(assignation.values()) and taxon[level]!='NA':
                        assignation[level] = 'NA'
                       

        return(assignation)

    def make_list_taxonomy(self, matrix_pident):
        taxonomy = []
        matrix_rank = matrix_pident.rank(
            method='min', ascending=False, na_option='bottom')
        series_sum_rank = matrix_rank.sum(axis=1)
        minimum_rank = series_sum_rank.min()
        top_series = series_sum_rank.loc[series_sum_rank[:] == minimum_rank]
        for individue, row in top_series.items():
            bestSCG = pd.to_numeric(matrix_pident.loc[individue, :]).idxmax()
            bestident = matrix_pident.loc[individue, bestSCG]
            taxonomy.append({"bestSCG": bestSCG, "bestident": bestident,
                             "accession": individue, "taxonomy":OrderedDict(self.taxonomy_dict[individue])})
        return(taxonomy)

    def reduce_assignation_level(self, taxonomy):
        reduce_taxonomy = []
        for possibilitie in taxonomy:
            for level, value in reversed(possibilitie["taxonomy"].items()):
                if possibilitie["taxonomy"][level] == 'NA':
                    continue
                if possibilitie["bestident"] < float(self.dicolevel[possibilitie["bestSCG"]]
                                                     [self.taxonomic_levels_parser[level] + value]):
                    possibilitie["taxonomy"][level] = 'NA'
                else:
                    break
            reduce_taxonomy.append(possibilitie["taxonomy"])

        return(reduce_taxonomy)

    def solo_hits(self, SCGs_hit_per_gene, cutoff_solo_hit=0.90):
        consensus_taxonomy = []
        for SCG, entry in SCGs_hit_per_gene.items():
            if len(entry) == 1 and entry[0]['pident'] > cutoff_solo_hit:
                consensus_taxonomy.append(
                    self.taxonomy_dict[entry[0]['accession']])
        if consensus_taxonomy:
            assignation = self.assign_taxonomie_solo_hit(consensus_taxonomy)
            return assignation
        else:
            return False


class SCGsdiamond(TaxonomyEstimation):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.taxonomy_file_path = A('taxonomy_file')
        self.taxonomy_database_path = A('taxonomy_database')
        self.write_buffer_size = int(A('write_buffer_size') if A(
            'write_buffer_size') is not None else 1000)
        self.db_path = A('contigs_db')
        self.core = A('num_threads')
        self.num_process = A('contigs_db')

        self.max_target_seqs = 20
        self.evalue = 1e-05
        self.min_pct_id = 90

        self.num_process = args.num_process

        if not self.core:
            self.core = "1"

        if not args.num_process:
            self.num_process = "1"

        self.initialized = False

        self.metagenome = False

        if args.metagenome:
            self.metagenome = True

        self.source = "unknow"

        if not self.taxonomy_file_path or not self.taxonomy_database_path:
            scg_ref_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG_TAXONOMY/GTDB/')
            scg_compres_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG.tar.gz')
            if not filesnpaths.is_file_exists(scg_ref_path, dont_raise=True) and \
             filesnpaths.is_file_exists(scg_compres_path, dont_raise=True) and filesnpaths.is_output_dir_writable(scg_ref_path):
                shutil.unpack_archive(scg_compres_path, scg_compres_path[:-10])

        if not self.taxonomy_file_path:
            self.taxonomy_file_path = os.path.join(scg_ref_path, 'ACCESSION_TO_TAXONOMY.txt')

        self.source = "https://gtdb.ecogenomic.org/"

        if not self.taxonomy_database_path:
            self.taxonomy_database_path = os.path.join(scg_ref_path, 'SCGs')

        self.database = db.DB(
            self.db_path, utils.get_required_version_for_db(self.db_path))

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

        self.taxonomy_dict = OrderedDict()

    def sanity_check(self):
        if not filesnpaths.is_file_exists(self.taxonomy_file_path, dont_raise=True):
            raise ConfigError("Anvi'o could not find taxonomy file '%s'. You must declare one before continue."
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

        TaxonomyEstimation.__init__(self, self.taxonomy_dict)
        self.initialized = True

    def get_hmm_sequences_dict_into_type_multi(self, hmm_sequences_dict,hmm_sequences_best_hit_dict):
        hmm_sequences_dict_per_type = {}
        
        
        for entry_id in hmm_sequences_dict:
            entry = hmm_sequences_dict[entry_id]

            SCG = entry['gene_name']
            if entry_id in hmm_sequences_best_hit_dict:
                entry['gene_callers_id'] = str(entry['gene_callers_id'])+"_best_hit"

            if SCG in hmm_sequences_dict_per_type:
                hmm_sequences_dict_per_type[SCG][entry_id] = entry
            else:
                hmm_sequences_dict_per_type[SCG] = {entry_id: entry}

        return hmm_sequences_dict_per_type

    def predict_from_SCGs_dict_multiseq(self, hmm_sequences_dict,hmm_sequences_best_hit_dict):
        """Takes an HMMs dictionary, and yields predictions"""

        self.init()

        hmm_sequences_dict_per_type = self.get_hmm_sequences_dict_into_type_multi(
            hmm_sequences_dict,hmm_sequences_best_hit_dict)

        num_listeprocess = len(hmm_sequences_dict_per_type)

        aligners = "Diamond"

        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Taxonomy', self.taxonomy_file_path)
        self.run.info('Database reference', self.taxonomy_database_path)
        self.run.info('Source', self.source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('Number of SCGs', len(hmm_sequences_dict_per_type))
        self.run.info('SCGs', ','.join(
            list(hmm_sequences_dict_per_type.keys())))

        self.run.warning(
            '', header='Parameters for aligments with %s for taxonomy' % aligners, lc='green')
        self.run.info('Blast type', "Blastp")
        self.run.info('Maximum number of target sequences',
                      self.max_target_seqs)
        self.run.info('Minimum bit score to report alignments',
                      self.min_pct_id)
        self.run.info('Number aligment running in same time', self.num_process)
        self.run.info(
            'Number of CPUs will be used for each aligment', self.core)

        self.tables_for_taxonomy = TablesForTaxoestimation(
            self.db_path, run, progress)
        self.tables_for_taxonomy.delete_contents_of_table(
            t.scg_taxonomy_table_name)

        self.progress.new('Computing SCGs aligments',

                          progress_total_items=num_listeprocess)
        self.progress.update('Initializing %d process...' %
                             int(self.num_process))

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        diamond_output = []
        table_index=0
        self.genes_estimation_output=[]

        for SCG in hmm_sequences_dict_per_type:

            sequence = ""
            for entry in hmm_sequences_dict_per_type[SCG].values():
                if 'sequence' not in entry or 'gene_name' not in entry:
                    raise ConfigError("The `get_filtered_dict` function got a parameter that\
                                       does not look like the way we expected it. This function\
                                       expects a dictionary that contains keys `gene_name` and `sequence`.")

                sequence = sequence + ">" + \
                    str(entry['gene_callers_id']) + \
                    "\n" + entry['sequence'] + "\n"
                entry['hits'] = []
            input_queue.put([SCG, sequence])

        workers = []
        for i in range(0, int(self.num_process)):

            worker = multiprocessing.Process(target=self.get_raw_blast_hits_multi,
                                             args=(input_queue, output_queue))

            workers.append(worker)
            worker.start()

        finish_process = 0
        while finish_process < num_listeprocess:
            try:
                diamond_output+=output_queue.get()

                if self.write_buffer_size > 0 and len(diamond_output) % self.write_buffer_size == 0:

                    table_index = self.tables_for_taxonomy.alignment_result_to_congigs(
                        table_index,diamond_output)
                    diamond_output = []

                finish_process += 1

                self.progress.increment(increment_to=finish_process)
                progress.update("Processed %s of %s SGCs aligment in %s processus with %s cores." % (
                    finish_process, num_listeprocess, int(self.num_process), self.core))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        table_index = self.tables_for_taxonomy.alignment_result_to_congigs(table_index,diamond_output)
        progress.end()

    def show_hits(self, hits, name):
        self.run.warning(None, header='%s' %
                         name, lc="green")
        header = ['%id', 'bitscore', 'taxonomy']
        table = []

        for hit in hits:
            table.append([hit['pident'], hit['bitscore'],
                          ' / '.join(self.taxonomy_dict[hit['accession']].values())])

        print(tabulate(table, headers=header,
                       tablefmt="fancy_grid", numalign="right"))

    def get_raw_blast_hits_multi(self, input_queue, output_queue):
        while True:
            d = input_queue.get(True)
            db_path = self.SCG_DB_PATH(d[0])
            diamond = Diamond(db_path, run=run_quiet, progress=progress_quiet)
            diamond.max_target_seqs = self.max_target_seqs
            diamond.evalue = self.evalue
            diamond.min_pct_id = self.min_pct_id
            diamond.num_process = self.core

            SCG = d[0]
            diamond_output = diamond.blastp_stdin_multi(d[1])
            hits_per_gene = {}
            genes_estimation_output=[]

            for line_hit_to_split in diamond_output.split('\n'):
                if len(line_hit_to_split) and not line_hit_to_split.startswith('Query'):
                    line_hit = line_hit_to_split.split('\t')
                    gene_callers_id = line_hit[0]
                    hit = [dict(zip(['accession', 'pident', 'bitscore'], [
                        line_hit[1], float(line_hit[2]), float(line_hit[11])]))]

                    if gene_callers_id not in hits_per_gene:
                        hits_per_gene[gene_callers_id] = {}
                    if SCG not in hits_per_gene[gene_callers_id]:
                        hits_per_gene[gene_callers_id][SCG] = []

                    hits_per_gene[gene_callers_id][SCG] += hit


            if anvio.DEBUG:
                self.show_hits(hits_per_gene[SCG], SCG)
            for gene_callers_id, SCGs_hit_per_gene in hits_per_gene.items():
                consensus_taxonomy, taxonomy = TaxonomyEstimation.get_consensus_taxonomy(self,
                    SCGs_hit_per_gene, gene_callers_id)
                genes_estimation_output.append([gene_callers_id, SCG, consensus_taxonomy, taxonomy])
            output_queue.put(genes_estimation_output)


class SCGsTaxonomy(TaxonomyEstimation):

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.output_file_path = A('output_file')

        self.collection_name = args.collection_name

        self.initialized = False

        self.bin_id = args.bin_id

        self.metagenome = False

        if args.metagenome:
            self.metagenome = True

        self.cut_off_methode = args.cut_off_methode

        self.hits_per_gene = {}


        self.taxonomic_levels_parser = {"t_domain": 'd__',
                                        "t_phylum": 'p__',
                                        "t_class": 'c__',
                                        "t_order": 'o__',
                                        "t_family": 'f__',
                                        "t_genus": 'g__',
                                        "t_species": 's__'}
        self.taxonomy_dict=dict()

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
        self.tables_for_taxonomy = TablesForTaxoestimation(
            self.db_path, run, progress, self.profile_db_path)

        self.dictonnary_taxonomy_by_index= self.tables_for_taxonomy.get_data_for_taxonomy_estimation()

        if not (self.dictonnary_taxonomy_by_index):
            raise ConfigError(
                "Anvi'o can't make a taxonomy estimation because aligment didn't return any match or you forgot to run 'anvi-diamond-for-taxonomy'.")

        if self.metagenome:
             self.initialized = True
             return



        if self.profile_db_path:
            self.hits_per_gene={}


            contigs_db = db.DB(self.db_path, anvio.__contigs__version__)

            collection_to_split = ccollections.GetSplitNamesInBins(self.args).get_dict()

            split_to_gene_callers_id = dict()
            bin_to_gene_callers_id = dict()



            for row in contigs_db.get_all_rows_from_table('genes_in_splits'):
                split_name, gene_callers_id = row[1], row[2]

                if split_name not in split_to_gene_callers_id:
                    split_to_gene_callers_id[split_name] = set()

                split_to_gene_callers_id[split_name].add(gene_callers_id)

            for bin_name in collection_to_split:
                for split in collection_to_split[bin_name]:
                    if bin_name not in bin_to_gene_callers_id:
                        bin_to_gene_callers_id[bin_name] = set()

                    if split in split_to_gene_callers_id:
                        bin_to_gene_callers_id[bin_name].update(split_to_gene_callers_id[split])

            
            for gene_estimation in self.dictonnary_taxonomy_by_index.values():
                if gene_estimation["source"]=="Consensus":
                    continue


                if gene_estimation['accession'] not in self.taxonomy_dict:
                    self.taxonomy_dict[gene_estimation['accession']]= {"t_domain": gene_estimation['t_domain'],
                                                                "t_phylum": gene_estimation['t_phylum'],
                                                                "t_class": gene_estimation['t_class'],
                                                                "t_order": gene_estimation['t_order'],
                                                                "t_family": gene_estimation['t_family'],
                                                                "t_genus": gene_estimation['t_genus'],
                                                                "t_species": gene_estimation['t_species']}

                for bin_id, gene_callers_id in bin_to_gene_callers_id.items():
                    #gene_id=int(gene_estimation['gene_caller_id'].replace("_best_hit",""))
                    if gene_estimation['gene_caller_id'] in gene_callers_id:

                        hit = [{'accession': gene_estimation['accession'], 'pident':float(
                            gene_estimation['pourcentage_identity'])}]

                        if bin_id not in self.hits_per_gene:
                            self.hits_per_gene[bin_id] = {}
                        if gene_estimation['gene_name'] not in self.hits_per_gene[bin_id]:
                            self.hits_per_gene[bin_id][gene_estimation['gene_name']] = []

                        self.hits_per_gene[bin_id][gene_estimation['gene_name']] += hit


            taxonomyestimation=TaxonomyEstimation.__init__(self,self.taxonomy_dict)

            self.initialized = True
            return

        else:
            self.SCGs_hit_per_gene={}
            for gene_estimation in self.dictonnary_taxonomy_by_index.values():
                self.taxonomy_dict[gene_estimation['accession']]= {"t_domain": gene_estimation['t_domain'],
                                                            "t_phylum": gene_estimation['t_phylum'],
                                                            "t_class": gene_estimation['t_class'],
                                                            "t_order": gene_estimation['t_order'],
                                                            "t_family": gene_estimation['t_family'],
                                                            "t_genus": gene_estimation['t_genus'],
                                                            "t_species": gene_estimation['t_species']}

                hit = [{'accession': gene_estimation['accession'], 'pident':float(gene_estimation['pourcentage_identity'])}]
                if gene_estimation['gene_name'] not in self.SCGs_hit_per_gene:
                    self.SCGs_hit_per_gene[gene_estimation['gene_name']] = []
                self.SCGs_hit_per_gene[gene_estimation['gene_name']] += hit
            taxonomyestimation=TaxonomyEstimation.__init__(self,self.taxonomy_dict)
            consensus_taxonomy, taxonomy = TaxonomyEstimation.get_consensus_taxonomy(self,
                    self.SCGs_hit_per_gene, self.db_path)
            output_full_genome=[['Genome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'],[self.db_path.replace(".db","")]+list(consensus_taxonomy.values())]
            self.show_taxonomy_estimation(output_full_genome)
            self.generate_outpu_file(output_full_genome)


        self.initialized = True
        

    def estimate_taxonomy(self, source="GTDB",number_scg=21):

        

        self.run.warning('', header='Taxonomy estimation for %s' %
                         self.db_path, lc='green')
        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Source', source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('output file for taxonomy', self.output_file_path)

        self.init()


        if self.metagenome:
            self.estimate_taxonomy_for_metagenome()


        if self.profile_db_path:
            possibles_taxonomy = []
            entry_id = 0
            entries_db_profile = []
            dictionary_bin_taxonomy_estimation=dict()
            possibles_taxonomy.append(
                ['Genome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

            for bin_id, SCGs_hit_per_gene in self.hits_per_gene.items():
                consensus_taxonomy, taxonomy = TaxonomyEstimation.get_consensus_taxonomy(self,
                    SCGs_hit_per_gene, bin_id)

                dictionary_bin_taxonomy_estimation[bin_id]={"consensus_taxonomy":consensus_taxonomy, "taxonomy_use_for_consensus":taxonomy}

                possibles_taxonomy.append([bin_id] + list(consensus_taxonomy.values()))


                entries_db_profile += [
                    (tuple([entry_id, self.collection_name, bin_id, source] + list(consensus_taxonomy.values())))]
                entry_id += 1

                self.tables_for_taxonomy.taxonomy_estimation_to_profile(
                    entries_db_profile)

            self.show_taxonomy_estimation(possibles_taxonomy)
            self.generate_outpu_file(possibles_taxonomy)

    def estimate_taxonomy_for_metagenome(self):
        output_genes_estimation = []
        estimate_taxonomy_presences=[]
        dictonarry_presence={}
        dictonnary_number_appear=dict()
        estimate_taxonomy_presences = [{"t_domain": "Unknow",
                                        "t_phylum": "NA",
                                        "t_class": "NA",
                                        "t_order": "NA",
                                        "t_family": "NA",
                                        "t_genus": "NA",
                                        "t_species": "NA"}]

        output_genes_estimation.append(
            ['genes_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

        for gene_estimation in self.dictonnary_taxonomy_by_index.values():
            if gene_estimation["source"] == "GTDB" :
                continue
            
            if gene_estimation['accession'] not in self.taxonomy_dict:
                taxonomy = {"t_domain": gene_estimation['t_domain'],
                                                            "t_phylum": gene_estimation['t_phylum'],
                                                            "t_class": gene_estimation['t_class'],
                                                            "t_order": gene_estimation['t_order'],
                                                            "t_family": gene_estimation['t_family'],
                                                            "t_genus": gene_estimation['t_genus'],
                                                            "t_species": gene_estimation['t_species']}
                self.taxonomy_dict[gene_estimation['accession']]=taxonomy

            output_genes_estimation.append([gene_estimation['gene_caller_id']] + list(taxonomy.values()))

            new=False
            for taxon_presente in taxonomy.values():
                if taxon_presente=="NA":
                    continue
                if taxon_presente not in dictonarry_presence and taxon_presente!="NA":
                    dictonarry_presence[taxon_presente]={gene_estimation['gene_name']: [gene_estimation['accession']]}
                    new=True
                else:
                    if gene_estimation['gene_name'] not in dictonarry_presence[taxon_presente]:
                        dictonarry_presence[taxon_presente][gene_estimation['gene_name']]=[gene_estimation['accession']]
                    else:
                        dictonarry_presence[taxon_presente][gene_estimation['gene_name']]+=[gene_estimation['accession']]
                        continue
            if new:
                for estimate_taxonomy_presence in estimate_taxonomy_presences:
                    share = { level : taxonomy[level] for level, taxon  in estimate_taxonomy_presence.items() & taxonomy.items() if taxon!="Bacteria" and taxon!="Archea" and taxon!="NA"}
                    difference = { level : taxonomy[level] for level, taxon  in estimate_taxonomy_presence.items() - taxonomy.items() if taxon=="NA"}
                    if difference and share:
                        estimate_taxonomy_presence.clear()
                        estimate_taxonomy_presence.update(taxonomy)
                        new=False
                        break
                if new:
                    estimate_taxonomy_presences+=[taxonomy]


        estimate_taxonomy_presences.pop(0)
        output=[['metagenome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species','number of SCG']]
        num_metagenome=1

        for estimate_taxonomy_presence in estimate_taxonomy_presences:
            for level in estimate_taxonomy_presence.values():
                if level not in dictonnary_number_appear:
                    dictonnary_number_appear[level]=1
                else:
                    dictonnary_number_appear[level]+=1

        if len(estimate_taxonomy_presences)>1:
            for estimate_taxonomy_presence in estimate_taxonomy_presences:
                output+=[[self.db_path.replace(".db","")+"_genome_"+str(num_metagenome)]+list(estimate_taxonomy_presence.values())+[dictonnary_number_appear[list(estimate_taxonomy_presence.values())[-1]]]]
                num_metagenome+=1
        else:
            output+=[[self.db_path.replace(".db","")]+list(estimate_taxonomy_presence.values())]
        



        self.show_taxonomy_estimation(output)
        

        outpu_appear=[["taxon","number of scg"]]
        for level, appear in dictonnary_number_appear.items():
            outpu_appear+=[[level],[appear]]

        self.show_taxonomy_estimation(outpu_appear)

        self.generate_outpu_file(output,append=True)
        for taxon, list_scgs in dictonarry_presence.items():
            if taxon=="NA":
                continue
            len_max_scg=("NA",0)
            for SCG, list_appear_scg in list_scgs.items():
                if len(list_appear_scg) > len_max_scg[1]:
                    len_max_scg=(SCG,len(list_appear_scg))
            if dictonnary_number_appear[taxon] < len_max_scg[1]:
                continue
                self.run.warning("%s is estimate %d time but it seams that the SCG %s have %d appear for this taxons, it could mean you have "\
                 % (taxon , dictonnary_number_appear[taxon], len_max_scg[0], len_max_scg[1]))


    def generate_outpu_file(self,output_data,header=False,append=False):
        if filesnpaths.is_output_file_writable(self.output_file_path,ok_if_exists=append):
            with open(self.output_file_path, "a") as output_file:
                output_data = ['\t'.join(line)
                                   for line in output_data]
                output_data='\n'.join(output_data)
                output_file.write(output_data)

    def show_taxonomy_estimation(self, possibles_taxonomy):
        self.run.warning(None, header='Taxonomy estimation', lc="yellow")
        print(tabulate(possibles_taxonomy, headers="firstrow",
                       tablefmt="fancy_grid", numalign="right"))



    def show_taxonomy_estimation_genes(self, possibles_taxonomy):
        self.run.warning(None, header='Taxonomy estimation', lc="yellow")
        possibles_taxonomy_dataframe = pd.DataFrame(
            possibles_taxonomy, columns=possibles_taxonomy[0])
        possibles_taxonomy_dataframe.set_index("bin_id", inplace=True)
        possibles_taxonomy_dataframe = possibles_taxonomy_dataframe.sort_values(
            by=['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'], ascending=False)
        print(tabulate(possibles_taxonomy_dataframe, headers="firstrow",
                       tablefmt="fancy_grid", numalign="right"))