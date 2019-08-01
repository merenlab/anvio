#!/usr/bin/env python3
# -*- coding: utf-8

import argparse
import os
import re
import sys
import time
import copy

import multiprocessing
from collections import Counter

import anvio
import pickle
from collections import OrderedDict

from tabulate import tabulate
import anvio.db as db


import anvio.ccollections as ccollections
import anvio.fastalib as f
import anvio.filesnpaths as filesnpaths
import anvio.hmmops as hmmops
import anvio.hmmopswrapper as hmmopswrapper
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.tables.tableops import Table
#from anvio.tables.taxoestimation import Tabletaxo
import anvio.tables as t
from anvio.dbops import ContigsSuperclass
from anvio.drivers import Aligners, driver_modules
from anvio.drivers.diamond import Diamond
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
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.taxonomy_file_path = A('taxonomy_file')
        self.taxonomy_database_path = A('taxonomy_database')
        self.db_path=A('contigs_db')
        self.core=int(A('num_core_by_threads'))


        self.num_threads=args.num_threads

        if not self.core:
            self.num_threads="2"

        if not args.num_threads:
            self.num_threads="2"

        self.initialized = False

        self.metagenome=False

        if args.metagenome:
            self.metagenome=True


        if not self.taxonomy_file_path:
            self.taxonomy_file_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG/mergedb/matching_taxonomy.tsv')

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
        sequence_by_SCG = {}


        match_id=0
        num_listeprocess = len(hmm_sequences_dict_per_type)


        Sequence_queu=[]
        sequence_by_SCG = []
        self.run.info('SCGs', ','.join(list(hmm_sequences_dict_per_type.keys())))

        self.progress.new('Computing SCGs aligments', progress_total_items=num_listeprocess)
        #self.progress.update('Initializing %d threads...' % self.num_threads)


        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        dico_taxo={}
        entries=[]

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
        for i in range(0, self.num_threads):
            #self.progress.update('Initializing %d core...' % self.core)
            worker = multiprocessing.Process(target=self.get_raw_blast_hits_multi,
                args=(match_id,input_queue,output_queue,self.core))

            workers.append(worker)
            worker.start()


        finish_process = 0
        while finish_process < num_listeprocess:
            try:
                diamond_output = output_queue.get()

                for line_hit in [line.split('\t') for line in diamond_output[1].split('\n')[1:-2]]:

                    entries+=[tuple([match_id,int(line_hit[0]),diamond_output[0],line_hit[1],line_hit[2],line_hit[11]])]


                    match_id+=1

                    dico_taxo[line_hit[1]]=list(self.taxonomy_dict[line_hit[1]].values())



                #self.database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % "taxon_names", taxo)
                finish_process += 1
                self.progress.increment(increment_to=finish_process)
                progress.update("ok")
                #progress.update("Processed %d of %d SGCs aligment in %d threads with %d cores." % (finish_process, num_listeprocess,self.num_threads,self.core))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                self.database.disconnect()
                break

        for worker in workers:
            worker.terminate()

        self.database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.blast_hits_table_name, entries)
        progress.end()
        self.run.info('Number of hit', match_id)

        self.database.disconnect()


    def get_raw_blast_hits_multi(self, match_id, input_queue,output_queue,core, max_target_seqs=20, evalue=1e-05, min_pct_id=90):

        while True:
            d = input_queue.get(True)
            db_path = self.SCG_DB_PATH(d[0])
            diamond = Diamond(db_path,run=run_quiet, progress= progress_quiet)
            diamond.max_target_seqs = max_target_seqs
            diamond.evalue = evalue
            diamond.min_pct_id = min_pct_id
            diamond.num_threads = core

            diamond_output = diamond.blastp_stdin_multi(d[1])


            output_queue.put([d[0],diamond_output])


    def get_raw_blast_hits_multi_original(self, d,match_id, max_target_seqs=20, evalue=1e-05, min_pct_id=90):

        sequence=""
        bin_dict_id={}
        for id, entry in d.items():
            if 'sequence' not in entry or 'gene_name' not in entry:
                raise ConfigError("The `get_filtered_dict` function got a parameter that\
                                   does not look like the way we expected it. This function\
                                   expects a dictionary that contains keys `gene_name` and `sequence`.")

            sequence = sequence+">"+id+"\n"+entry['sequence']+"\n"
            entry['hits']=[]
            bin_dict_id[id]=entry['bin_id']

        db_path = self.SCG_DB_PATH(entry['gene_name'])
        diamond = Diamond(db_path,run=run_quiet, progress= progress_quiet)
        diamond.max_target_seqs = max_target_seqs
        diamond.evalue = evalue
        diamond.min_pct_id = min_pct_id
        diamond.num_threads = self.num_threads

        diamond_output = diamond.blastp_stdin_multi(sequence)


        for line_hit in [line.split('\t') for line in diamond_output.split('\n') if line.startswith('Bacteria')]:


            entries=[tuple([match_id,bin_dict_id[line_hit[0]],line_hit[0],entry['gene_name'],line_hit[1],line_hit[2],line_hit[11]])]

            self.database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % self.blast_hits_table_name, entries)

            taxo=[tuple([line_hit[1]]+list(self.taxonomy_dict[line_hit[1]].values()))]

            try:
                self.database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % "taxon_names", taxo)
            except:
                pass

            hit=dict(zip(['accession', 'pident', 'bitscore'], [
                       float(line_hit[i]) if i > 1 else line_hit[i] for i in [1, 2, 11]]))
            taxo=self.taxonomy_dict[line_hit[1]]
            d[line_hit[0]]['hits']=d[line_hit[0]]['hits']+[hit]

            match_id+=1

        return d, match_id

class SCGsTaxomy:

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.db_path=A('contigs_db')
        self.profile_database_path=A('profile_db')

        self.collection_name=args.collection_name

        self.initialized = False

        self.bin_id=args.bin_id

        self.metagenome=False

        if args.metagenome:
            self.metagenome=True

        self.cut_off_methode = args.cut_off_methode

        self.methode = args.methode

        self.taxonomy_dict = {}

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

        self.contigs_database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))


    def sanity_check(self):
        filesnpaths.is_file_exists(self.db_path)
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

        self.dic_blast_hits=self.contigs_database.get_table_as_dict('blast_hits')
        self.taxonomy_dict =OrderedDict(self.contigs_database.get_table_as_dict('taxon_names'))


        self.dic_id_bin={}
        if self.profile_database_path and self.metagenome:
            self.run.info('Assignment level', "Bin")
            utils.is_profile_db_and_contigs_db_compatible(
                self.profile_database_path, self.db_path)

            self.bin_database = db.DB(self.profile_database_path, utils.get_required_version_for_db(self.profile_database_path))

            splits_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()

            run.info('Init', '%d splits in %d bin(s)' % (
                sum([len(v) for v in list(splits_dict.values())]), len(splits_dict)))

            s = hmmops.SequencesForHMMHits(self.db_path)

            hits_in_splits, split_name_to_bin_id = s.get_hmm_hits_in_splits(splits_dict)

            dic_genes_in_splits=self.contigs_database.get_table_as_dict("genes_in_splits")

            for split in dic_genes_in_splits.values():
                if split['split'] in split_name_to_bin_id:
                    if split_name_to_bin_id[split['split']] not in self.dic_id_bin:
                        self.dic_id_bin[split_name_to_bin_id[split['split']]]=[split['gene_callers_id']]
                    else:
                        self.dic_id_bin[split_name_to_bin_id[split['split']]]+=[split['gene_callers_id']]
        else:
            self.run.info('Assignment level', "Gene")
            self.bin_database=None

        self.hits_per_gene={}
        for query in self.dic_blast_hits.values():

            if self.bin_database and not self.metagenome:
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

        self.initialized = True


    def estimate_taxonomy(self):

        self.init()

        possibles_taxonomy=[]
        entry_id=0
        self.run.info('Assignment level', "Gene")
        for name, SCGs_hit_per_gene in self.hits_per_gene.items():

            taxonomy = self.get_consensus_taxonomy(
                SCGs_hit_per_gene, name)

            if not taxonomy:
                self.run.info('taxonomy estimation not possible for:', name)
                continue


            if not self.metagenome:
                self.run.info('estimate taxonomy',
                              '/'.join(list(taxonomy.values())))

                entries=[tuple([entry_id,self.collection_name,name,"GTDB"]+list(taxonomy.values()))]

                if self.bin_database:
                    self.bin_database.insert_many(t.taxonomy_estimation_bin_name, entries)


            if self.metagenome or not self.bin_database:

                entries=[tuple([name,list(SCGs_hit_per_gene.keys())[0],"GTDB"]+list(taxonomy.values()))]

                self.contigs_database.insert_many(t.taxonomy_estimation_metagenome_name, entries)

                if str(list(taxonomy.values())[-1]) not in possibles_taxonomy and taxonomy:
                    possibles_taxonomy.append(str(list(taxonomy.values())[-1]))

            entry_id+=1

        if self.metagenome or not self.bin_database:
            if len(possibles_taxonomy):
                self.run.info('Possible presence ','|'.join(list(possibles_taxonomy)))

        if self.bin_database:
            self.bin_database.disconnect()
        self.contigs_database.disconnect()



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


    def show_matrix_rank(self, name, matrix, list_position_entry, list_position_ribosomal):
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


    def get_consensus_taxonomy(self, SCGs_hit_per_gene, name):
        """Different methode for assignation"""

        consensus_taxonomy=self.solo_hits(SCGs_hit_per_gene)

        if consensus_taxonomy:

            return(consensus_taxonomy)

        else:

            if self.methode == "friedman":
                consensus_taxonomy = self.rank_assignement(SCGs_hit_per_gene, name)

            """if self.methode == "tree":
                self.make_tree_with_hit(
                    hits_per_gene, hmm_sequences_dict_per_type, name)
                consensus_taxonomy = "tree"
                return(consensus_taxonomy)"""

            if self.methode == "bitscore":
                score_by_entry = self.get_cumul_hit_per_gene(SCGs_hit_per_gene)
                if anvio.DEBUG:
                    self.show_table_score(name, score_by_entry)
                consensus_taxonomy = self.get_consensus_taxonomy_with_score_by_entry(
                    score_by_entry, name, self.cut_off_methode)


            return(consensus_taxonomy)


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


    def get_cumul_hit_per_gene(self, SCGs_hit_per_gene):
        """add bitscore of eatch  per query"""
        matching_genes = [
            gene for gene in SCGs_hit_per_gene if len(SCGs_hit_per_gene[gene])]
        cumul_hit_per_gene = {}
        number_of_matching_genes = len(matching_genes)
        for matching_gene in matching_genes:
            for entry in SCGs_hit_per_gene[matching_gene]:
                if entry["accession"] not in cumul_hit_per_gene:
                    cumul_hit_per_gene[entry['accession']] = (
                        float(entry['bitscore']) / number_of_matching_genes)
                else:
                    cumul_hit_per_gene[entry['accession']] =\
                        float(cumul_hit_per_gene[entry['accession']])\
                        + (float(entry['bitscore']) / number_of_matching_genes)
        return cumul_hit_per_gene


    def get_matching_gene(self, SCGs_hit_per_gene):
        matching_genes = [
            gene for gene in SCGs_hit_per_gene if len(SCGs_hit_per_gene[gene])]
        number_of_matching_genes = len(matching_genes)
        return matching_genes


    def align_hit_sequence(self, entry, SCGs_hit_per_gene):

        sequences_match = []
        sequences_match.append((entry['bin_id'], entry['sequence']))
        db_path = self.SCG_FASTA_DB_PATH(entry['gene_name'])
        fasta_db = f.ReadFasta(db_path)
        maxscore = 0

        for hit in SCGs_hit_per_gene[entry['gene_name']]:
            index = fasta_db.ids.index(hit['accession'])
            species_name = hit['taxonomy']['t_species'].split(" ")
            species_name = species_name[1]
            if float(hit['bitscore']) > float(maxscore):
                maxscore = hit['bitscore']
                sequences_match.append(
                    (species_name, fasta_db.sequences[index]))
            if float(hit['bitscore']) < (float(maxscore) * float(self.cut_off_methode)):
                continue
            sequences_match.append((species_name, fasta_db.sequences[index]))
        return(sequences_match)


    def align_with_muscle(self, entry, sequences_match, name):

        self.align_with = "muscle"
        aligners.select(self.align_with)
        r = terminal.Run()
        program = driver_modules['phylogeny']['default']
        aligner = aligners.select("muscle", quiet=True)
        alignments = aligner(run=r).run_stdin(sequences_match)

        temp_file_path_muscle = filesnpaths.get_temp_file_path(prefix='ANVIO_muscle_%s'
                                                               % str(name) + str(entry['gene_name']) + ".fa")

        with open(temp_file_path_muscle, 'w') as muscle:
            for id, sequence in alignments.items():
                muscle.write('>%s\n%s\n' % (id, sequence))
        output_file_path = str(name) + "_" + str(entry['gene_name']) + ".tree"

        filesnpaths.is_output_file_writable(output_file_path)

        program().run_command(temp_file_path_muscle, output_file_path)


    def make_tree_with_hit(self, SCGs_hit_per_gene, hmm_sequences_dict_per_type, name):
        matching_genes = self.get_matching_gene(SCGs_hit_per_gene)
        for entry in hmm_sequences_dict_per_type[name].values():
            if entry['gene_name'] in matching_genes:
                sequences_match = self.align_hit_sequence(entry, SCGs_hit_per_gene)
            if len(sequences_match) > 3:
                self.align_with_muscle(entry, sequences_match, name)
            else:
                print("for " + str(name) + " the protein: " + str(entry['gene_name']) + "\t match " +
                      SCGs_hit_per_gene[entry['gene_name']][0]['taxonomy']['t_species'])


    def make_dico_position_entry(self, SCGs_hit_per_gene, matchinggenes):
        dico_position_entry = {}
        for hit in matchinggenes:
            i = 0
            for entry in SCGs_hit_per_gene[hit]:
                if entry['accession'] not in dico_position_entry:
                    dico_position_entry[entry['accession']] = 1
                else:
                    dico_position_entry[entry['accession']] += 1
                i += 1
        return(dico_position_entry)


    def make_list_position_entry(self, SCGs_hit_per_gene, dico_position_entry, list_position_ribosomal, matchinggenes):
        list_position_entry = []
        j = 0
        len_position = len(list_position_ribosomal)
        key_value = sorted(dico_position_entry,
                           key=lambda x: dico_position_entry[x], reverse=True)
        for key in key_value:
            # parameter number of individue in the first matrix "j <= 5 " and appear in minimun 50% of blast
            if j< 10 and dico_position_entry[key] > (len_position * 0.5):
                list_position_entry.append(key)
                j += 1

        list_position_ribosomal = self.make_list_position_ribosomal(
            matchinggenes, SCGs_hit_per_gene, list_position_entry)
        if not len(list_position_entry) or not len(list_position_ribosomal) :
            sys.exit(status=None)
        return(list_position_entry, list_position_ribosomal)


    def make_list_position_ribosomal(self, matchinggenes, SCGs_hit_per_gene, list_position_entry):
        list_position_ribosomal = []
        for hit in matchinggenes:
            k = 0
            if not list_position_entry:
                list_position_ribosomal.append(hit)
            else:
                for entry in SCGs_hit_per_gene[hit]:
                    if entry['accession'] in list_position_entry:
                        k += 1
                # parameter, k number of
                if k >= len(list_position_entry):
                    list_position_ribosomal.append(hit)

        return(list_position_ribosomal)


    def creat_list_position(self, SCGs_hit_per_gene, matchinggenes):
        dico_position_entry = self.make_dico_position_entry(
            SCGs_hit_per_gene, matchinggenes)
        list_position_entry = []
        list_position_ribosomal = self.make_list_position_ribosomal(
            matchinggenes, SCGs_hit_per_gene, list_position_entry)
        list_position_entry, last_list_position_ribosomal = self.make_list_position_entry(
            SCGs_hit_per_gene, dico_position_entry, list_position_ribosomal, matchinggenes)

        return(list_position_entry, last_list_position_ribosomal)


    def emptylistmaker(self, numbr_genes, default_value='NA'):
        listofempty = [default_value] * numbr_genes
        return listofempty


    def make_emptymatrix(self, list_position_entry, list_position_ribosomal):
        matrix = []
        numbr_genes = len(list_position_ribosomal)
        maxposition = len(list_position_entry)
        i = 0
        # do 'matrix=[(['NA'] *len(list_position_ribosomal))]*(len(dico_position_entry))' don't work.\
        # subsitution will impact all items on the same index in lists
        while i < maxposition:
            emptylist = self.emptylistmaker(numbr_genes)
            matrix.append(emptylist)
            i += 1
        return(matrix)


    def fill_matrix(self, name, emptymatrix, SCGs_hit_per_gene, list_position_entry,
                    list_position_ribosomal, matchinggenes):
        maxrank = 1

        emptymatrix_ident= copy.deepcopy(emptymatrix)

        for hit in list_position_ribosomal:
            rank = 0
            lastscore= 1000
            for entry in SCGs_hit_per_gene[hit]:
                if entry['accession'] in list_position_entry:

                    # parameter for considere 2 hit have same rank rank
                    if float(entry['pident']) < lastscore:
                        rank += 1
                        lastscore=float(entry['pident'])

                    emptymatrix, emptymatrix_ident = self.fill_position_matrix(
                        emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, entry['pident'], rank, entry, hit)
                    lastscore=entry['pident']


        return(emptymatrix, emptymatrix_ident, maxrank)


    def fill_position_matrix(self, emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, pident, rank, entry, hit):
        emptymatrix[list_position_entry.index(
            entry['accession'])][list_position_ribosomal.index(hit)] = rank

        emptymatrix_ident[list_position_entry.index(
            entry['accession'])][list_position_ribosomal.index(hit)] = float(pident)

        return(emptymatrix,emptymatrix_ident)

    def fill_NA_matrix(self, matrix, matrix_pident, penality,max_target_seqs=20):
        j=0
        for liste in matrix:

            for n, i in enumerate(liste):
                if str(i) == 'NA':
                    liste[n] = (int(penality) + int(max_target_seqs))
                    matrix_pident[j][n]=0
            j+=1
        return(matrix,matrix_pident)


    def make_rank_matrix(self, name, SCGs_hit_per_gene):
        matchinggenes = self.get_matching_gene(SCGs_hit_per_gene)

        list_position_entry, list_position_ribosomal = self.creat_list_position(
            SCGs_hit_per_gene, matchinggenes)

        emptymatrix = self.make_emptymatrix(
            list_position_entry, list_position_ribosomal)

        matrix, matrix_pident, maxrank= self.fill_matrix(name, emptymatrix, SCGs_hit_per_gene, list_position_entry,
                                                      list_position_ribosomal, matchinggenes)

        final_matrix,matrix_pident_final = self.fill_NA_matrix(matrix, matrix_pident, maxrank)

        if anvio.DEBUG:
            self.show_matrix_rank(
                name, final_matrix, list_position_entry, list_position_ribosomal)

            self.show_matrix_rank(
                name, matrix_pident_final, list_position_entry, list_position_ribosomal)

        return(final_matrix, matrix_pident_final, list_position_entry, list_position_ribosomal)


    def statistic_test_friedman(self, name, matrix, matrixlist_position_entry, list_position_ribosomal):
        pvalue = 0
        newmatrix = matrix
        new_matrixlist_position_entry = matrixlist_position_entry
        while len(newmatrix) > 2:
            statistic, pvalue = ss.friedmanchisquare(*newmatrix)
            print("statistic:%s pvalue:%f" % (statistic, pvalue))
            if pvalue >= 0.05 or pvalue == "nan":
                return(newmatrix, new_matrixlist_position_entry)
                break
            new_matrixlist_position_entry, newmatrix = self.del_higher_ranked_individu_matrix(
                newmatrix, new_matrixlist_position_entry)
        return(newmatrix, new_matrixlist_position_entry)


    def rank_assignement(self, SCGs_hit_per_gene, name):

        matrix, matrix_pident, matrixlist_position_entry, list_position_ribosomal = self.make_rank_matrix(
            name, SCGs_hit_per_gene)
        start_make_list_taxonomy = time.perf_counter()

        taxonomy= self.make_list_taxonomy(
            matrix_pident, matrix, matrixlist_position_entry,list_position_ribosomal)

        assignation_reduce = self.reduce_assignation_level(taxonomy)
        assignation = self.assign_taxonomie_solo_hit(assignation_reduce)

        return(assignation)


    def del_higher_ranked_individu_matrix(self, matrix, matrixlist_position_entry):
        maxrank = 0
        i = 0
        for liste in matrix:
            summist = sum(list(liste))
            if summist > maxrank:
                maxrank = summist
                dellist = i
            i += 1
        del matrix[dellist]
        del matrixlist_position_entry[dellist]
        return matrixlist_position_entry, matrix


    def assign_taxonomie_solo_hit(self, taxonomy):
        if not taxonomy or not taxonomy[0]:
            return 0
        assignation = taxonomy[0]

        for s in taxonomy[1:]:
            for key in s:
                if s[key] not in assignation.values():
                    if anvio.DEBUG:
                        self.run.info('CONSIDER taxonomy','/'.join(s.values()))
                    assignation[key]=''
        return(assignation)


    def make_list_taxonomy(self, matrix_pident, matrix, list_position_entry,list_position_ribosomal):
        taxonomy = []
        miniscore = 10000
        i = 0
        while i < len(matrix):
            summist = sum(matrix[i])
            bestident = float(max(matrix_pident[i]))
            bestSCG=list_position_ribosomal[matrix_pident[i].index(bestident)]
            if summist == miniscore:
                taxonomy.append({"bestSCG" : bestSCG,"bestident" : bestident,"taxo" : OrderedDict(self.taxonomy_dict[list_position_entry[i]])})
            if summist < miniscore:
                taxonomy=[]
                taxonomy.append({"bestSCG" : bestSCG,"bestident" : bestident,"taxo" : OrderedDict(self.taxonomy_dict[list_position_entry[i]])})
                miniscore = summist
            i += 1
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


    def solo_hits(self, SCGs_hit_per_gene):
        taxonomy = []
        for entry in SCGs_hit_per_gene.values():
            if len(entry) == 1 and entry[0]['pident'] > 0.90:
                taxonomy.append(self.taxonomy_dict[entry[0]['accession']])
        if taxonomy:
            assignation=self.assign_taxonomie_solo_hit(taxonomy)
            return assignation
        else:
            return False
