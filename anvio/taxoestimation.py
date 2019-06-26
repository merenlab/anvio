
#!/usr/bin/env python3
# -*- coding: utf-8

import argparse
import os
import re
import sys
import time
from collections import Counter

import anvio
import anvio.ccollections as ccollections
import anvio.fastalib as f
import anvio.filesnpaths as filesnpaths
import anvio.hmmops as hmmops
import anvio.hmmopswrapper as hmmopswrapper
import anvio.terminal as terminal
import anvio.utils as utils
from anvio.dbops import ContigsSuperclass
from anvio.drivers import Aligners, driver_modules
from anvio.drivers.diamond import Diamond
from anvio.errors import ConfigError, FilesNPathsError
from tabulate import tabulate

run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)
aligners = Aligners()


class SCGTaxonomy:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None
        self.taxonomy_file_path = A('taxonomy_file')
        self.taxonomy_database_path = A('taxonomy_database')

        self.num_threads=args.num_threads

        if not args.num_threads:
            self.num_threads="2"

        self.initialized = False

        if not self.taxonomy_file_path:
            self.taxonomy_file_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG/mergedb/matching_taxonomy.tsv')

        if not self.taxonomy_database_path:
            self.taxonomy_database_path = os.path.join(os.path.dirname(
                anvio.__file__), 'data/misc/SCG/mergedb/species/')

        self.cut_off_methode = args.cut_off_methode

        self.methode = args.methode

        self.profile_db = args.profile_db

        """self.SCGs = ["Ribosomal_L1", "Ribosomal_L13", "Ribosomal_L16", "Ribosomal_L17",
                     "Ribosomal_L2", "Ribosomal_L20", "Ribosomal_L21p", "Ribosomal_L22",
                     "Ribosomal_L27A", "Ribosomal_L3", "Ribosomal_L4", "Ribosomal_L6",
                     "Ribosomal_L9_C", "Ribosomal_S11", "Ribosomal_S2", "Ribosomal_S20p",
                     "Ribosomal_S3_C", "Ribosomal_S6", "Ribosomal_S7", "Ribosomal_S9",
                     "ribosomal_L24"]"""

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

        print(self.SCG_DB_PATH)

        """self.SCG_DB_PATH = lambda db: os.path.join(
            self.taxonomy_database_path, [db for db in os.listdir(self.taxonomy_database_path) if db.endswith(".dmnd")])

        self.SCG_DB_PATH = lambda SCG: os.path.join(
            self.taxonomy_database_path, "%s.dmnd" % SCG)

        self.SCG_FASTA_DB_PATH = lambda SCG: os.path.join(
            self.taxonomy_database_path, SCG)"""

        self.SCG_FASTA_DB_PATH = lambda SCG: os.path.join(self.taxonomy_database_path,
                                                          [db for db in os.listdir(self.taxonomy_database_path) if db.endwith(".dmnd")])

        self.taxonomy_dict = {}

        self.sanity_check()

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

    def sanity_check(self):
        filesnpaths.is_file_exists(self.taxonomy_file_path)
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

    @timer
    def init(self):
        if self.initialized:
            return
        # initialize taxonomy dict. we should do this through a database in the long run.
        with open(self.taxonomy_file_path, 'r') as taxonomy_file:
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

        self.initialized = True

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
            print(list_position_entry[i])
            taxonomyindividue = list(
                self.taxonomy_dict[list_position_entry[i]].values())
            print(taxonomyindividue)
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

    @timer
    def get_hmm_sequences_dict_into_type(self, hmm_sequences_dict):
        hmm_sequences_dict_per_type = {}

        for entry_id in hmm_sequences_dict:
            entry = hmm_sequences_dict[entry_id]

            """if self.profile_db:
                name = entry['bin_id']
            else:"""
            name = entry['gene_name']

            if name in hmm_sequences_dict_per_type:
                hmm_sequences_dict_per_type[name][entry_id] = entry
            else:
                hmm_sequences_dict_per_type[name] = {entry_id: entry}

        return hmm_sequences_dict_per_type

    @timer
    def predict_from_SCGs_dict(self, hmm_sequences_dict):
        """Takes an HMMs dictionary, and yields predictions"""

        self.init()

        # split hmm_sequences_dict
        hmm_sequences_dict_per_type = self.get_hmm_sequences_dict_into_type(
            hmm_sequences_dict)

        for name in hmm_sequences_dict_per_type:
            start_predict_from_SCGs_dict = time.perf_counter()
            hits_per_gene = {}

            self.run.info('Bin name', name)
            self.run.info('Num SCGs', len(hmm_sequences_dict_per_type[name]))
            if self.profile_db:
                self.run.info('SCGs', ', '.join(
                    [e['gene_name'] for e in hmm_sequences_dict_per_type[name].values()]))
            else:
                self.run.info('contigs with SCGs', ', '.join(
                    [e['contig'] for e in hmm_sequences_dict_per_type[name].values()]))

            j = 0

            for entry in hmm_sequences_dict_per_type[name].values():
                gene_name = entry['gene_name']

                start_get_raw_blast_hits = time.perf_counter()

                hits = self.get_raw_blast_hits(entry)

                end_get_raw_blast_hits = time.perf_counter()
                print("\n time predict_from_SCGs_dict for a bin : ",
                      end_get_raw_blast_hits - start_get_raw_blast_hits)
                continue

                if not hits:
                    j += 1
                    continue

                # replace accessions with taxonomy
                for hit in hits:
                    hit['taxonomy'] = self.taxonomy_dict[hit['accession']]

                hits_per_gene[gene_name] = hits

                if anvio.DEBUG:
                    self.show_hits(name, gene_name, hits)

                if not self.profile_db and hits:
                    taxonomy = self.get_consensus_taxonomy(
                        hits_per_gene, name, hmm_sequences_dict_per_type)
                    entry['taxonomy'] = taxonomy
                    # print(entry)

            if j >= len(hits_per_gene):
                self.run.info(name, "diamond didn't return any match \n")
                continue
            if self.profile_db:
                taxonomy = self.get_consensus_taxonomy(
                    hits_per_gene, name, hmm_sequences_dict_per_type)

            end_predict_from_SCGs_dict = time.perf_counter()
            print("\n time predict_from_SCGs_dict for a bin : ",
                  end_predict_from_SCGs_dict - start_predict_from_SCGs_dict)
    @timer
    def predict_from_SCGs_dict_multiseq(self, hmm_sequences_dict):
        """Takes an HMMs dictionary, and yields predictions"""

        self.init()

        # split hmm_sequences_dict
        hmm_sequences_dict_per_type = self.get_hmm_sequences_dict_into_type(
            hmm_sequences_dict)
        j=0
        for name in hmm_sequences_dict_per_type:
            start_predict_from_SCGs_dict = time.perf_counter()
            hits_per_gene = {}

            self.run.info('SCGs', name)
            self.run.info('Num SCGs', len(hmm_sequences_dict_per_type[name]))

            """self.run.info('SCGs', ', '.join(
                [e['gene_name'] for e in hmm_sequences_dict_per_type[name].values()]))"""

            start_get_raw_blast_hits = time.perf_counter()

            hmm_sequences_dict_per_type[name] = self.get_raw_blast_hits_multi(hmm_sequences_dict_per_type[name])

            end_get_raw_blast_hits = time.perf_counter()
            print("\n time predict_from_SCGs_dict for a bin : ",
                  end_get_raw_blast_hits - start_get_raw_blast_hits)

            continue
            # replace accessions with taxonomy
            for hit in hmm_sequences_dict_per_type[name]['hits']:
                hit['taxonomy'] = self.taxonomy_dict[hit['accession']]

            ###
            hits_per_gene[gene_name] = hmm_sequences_dict_per_type[name]['hits']

            if anvio.DEBUG:
                self.show_hits(name, gene_name, hmm_sequences_dict_per_type[name]['hits'])

            taxonomy = self.get_consensus_taxonomy(
                hits_per_gene, name, hmm_sequences_dict_per_type)
            entry['taxonomy'] = taxonomy
                # print(entry)

            if j >= len(hits_per_gene):
                self.run.info(name, "diamond didn't return any match \n")
                continue
            if self.profile_db:
                taxonomy = self.get_consensus_taxonomy(
                    hits_per_gene, name, hmm_sequences_dict_per_type)

                end_predict_from_SCGs_dict = time.perf_counter()
                print("\n time predict_from_SCGs_dict for a bin : ",
                      end_predict_from_SCGs_dict - start_predict_from_SCGs_dict)



    @timer
    def get_raw_blast_hits(self, d, max_target_seqs=20, evalue=1e-05, min_pct_id=90):
        """Takes a dictionary that contains `gene_name` and `sequence`, and returns
           filtered BLAST hits against the corresopnding database. I.e.,

            >>> d = {'sequence': 'MVRVKKGVNALKTRRNILKQAKGFRGPRKSKEKLAYEQLVHSYTSAFAHRRDKKGDFRRLWNVRINAALRPLGHTYSKFIGAMNKKGMEVDRKTLSDLAQNAPESFERLVKQVTA',
                     'gene_name': 'Ribosomal_L20'}
            >>> self.get_raw_blast_hits(d)
        """

        if 'sequence' not in d or 'gene_name' not in d:
            raise ConfigError("The `get_filtered_dict` function got a parameter that\
                               does not look like the way we expected it. This function\
                               expects a dictionary that contains keys `gene_name` and `sequence`.")

        db_path = self.SCG_DB_PATH(d['gene_name'])
        if db_path:
            sequence = d['sequence']

            diamond = Diamond(db_path)
            diamond.max_target_seqs = max_target_seqs
            diamond.evalue = evalue
            diamond.min_pct_id = min_pct_id
            diamond.num_threads=self.num_threads

            diamond_output = diamond.blastp_stdin(sequence)

            hits = []
            for entry in [line.split('\t') for line in diamond_output.split('\n') if line.startswith('seq')]:
                accession = entry[1]

                # dict(zip(['accession', 'pident', 'length', 'mismatch', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], [float(entry[i]) if i > 1 else entry[i] for i in range(1, 12)]))
                hit = dict(zip(['accession', 'pident', 'bitscore'], [
                           float(entry[i]) if i > 1 else entry[i] for i in [1, 2, 11]]))

                hits.append(hit)

            return hits
        else:
            raise ConfigError("No matching database for" + d['gene_name'])



    @timer
    def get_raw_blast_hits_multi(self, d, max_target_seqs=20, evalue=1e-05, min_pct_id=90):

        sequence=""
        for id, entry in d.items():
            if 'sequence' not in entry or 'gene_name' not in entry:
                raise ConfigError("The `get_filtered_dict` function got a parameter that\
                                   does not look like the way we expected it. This function\
                                   expects a dictionary that contains keys `gene_name` and `sequence`.")

                #sequence = sequence+">%s\n%s\n"% (id,entry['sequence'])
            sequence = sequence+">"+id+"\n"+entry['sequence']+"\n"
            entry['hits']=[]

        db_path = self.SCG_DB_PATH(entry['gene_name'])
        diamond = Diamond(db_path)
        diamond.max_target_seqs = 999
        diamond.evalue = evalue
        diamond.min_pct_id = min_pct_id

        diamond_output = diamond.blastp_stdin_multi(sequence)


        for entry in [line.split('\t') for line in diamond_output.split('\n') if line.startswith('Bacteria')]:
            accession = entry[1]


            # dict(zip(['accession', 'pident', 'length', 'mismatch', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], [float(entry[i]) if i > 1 else entry[i] for i in range(1, 12)]))
            hit=dict(zip(['accession', 'pident', 'bitscore'], [
                       float(entry[i]) if i > 1 else entry[i] for i in [1, 2, 11]]))
            d[entry[0]]['hits']=d[entry[0]]['hits']+[hit]

        return d





    @timer
    def get_consensus_taxonomy(self, hits_per_gene, name, hmm_sequences_dict_per_type):
        """Different methode for assignation"""

        if not self.solo_hits(hits_per_gene):

            self.run.info("Assignation methode : ", self.methode)

            if self.methode == "friedman":
                consensus_taxonomy = self.rank_assignement(hits_per_gene, name)

            if self.methode == "tree":
                self.make_tree_with_hit(
                    hits_per_gene, hmm_sequences_dict_per_type, name)
                consensus_taxonomy = "tree"
                return(consensus_taxonomy)

            if self.methode == "bitscore":
                score_by_entry = self.get_cumul_hit_per_gene(hits_per_gene)
                if anvio.DEBUG:
                    self.show_table_score(name, score_by_entry)
                consensus_taxonomy = self.get_consensus_taxonomy_with_score_by_entry(
                    score_by_entry, name, self.cut_off_methode)

            if not self.profile_db:
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

    def get_cumul_hit_per_gene(self, hits_per_gene):
        """add bitscore of eatch blast per query"""
        matching_genes = [
            gene for gene in hits_per_gene if len(hits_per_gene[gene])]
        cumul_hit_per_gene = {}
        number_of_matching_genes = len(matching_genes)

        self.run.info('Num SCGs with match(s)', number_of_matching_genes)
        self.run.info('SCGs with match(s)', ', '.join(matching_genes))

        for matching_gene in matching_genes:
            for entry in hits_per_gene[matching_gene]:
                if entry["accession"] not in cumul_hit_per_gene:
                    cumul_hit_per_gene[entry['accession']] = (
                        float(entry['bitscore']) / number_of_matching_genes)
                else:
                    cumul_hit_per_gene[entry['accession']] =\
                        float(cumul_hit_per_gene[entry['accession']])\
                        + (float(entry['bitscore']) / number_of_matching_genes)
        return cumul_hit_per_gene

    def get_matching_gene(self, hits_per_gene):
        matching_genes = [
            gene for gene in hits_per_gene if len(hits_per_gene[gene])]
        number_of_matching_genes = len(matching_genes)
        self.run.info('Num SCGs with match(s)', number_of_matching_genes)
        self.run.info('SCGs with match(s)', ', '.join(matching_genes))
        return matching_genes

    def align_hit_sequence(self, entry, hits_per_gene):

        sequences_match = []
        sequences_match.append((entry['bin_id'], entry['sequence']))
        db_path = self.SCG_FASTA_DB_PATH(entry['gene_name'])
        fasta_db = f.ReadFasta(db_path)
        maxscore = 0

        for hit in hits_per_gene[entry['gene_name']]:
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

    def make_tree_with_hit(self, hits_per_gene, hmm_sequences_dict_per_type, name):
        matching_genes = self.get_matching_gene(hits_per_gene)
        for entry in hmm_sequences_dict_per_type[name].values():
            if entry['gene_name'] in matching_genes:
                sequences_match = self.align_hit_sequence(entry, hits_per_gene)
            if len(sequences_match) > 3:
                self.align_with_muscle(entry, sequences_match, name)
            else:
                print("for " + str(name) + " the protein: " + str(entry['gene_name']) + "\t match " +
                      hits_per_gene[entry['gene_name']][0]['taxonomy']['t_species'])

    def make_dico_position_entry(self, hits_per_gene, matchinggenes):
        dico_position_entry = {}
        for hit in matchinggenes:
            i = 0
            for entry in hits_per_gene[hit]:
                if entry['accession'] not in dico_position_entry:
                    dico_position_entry[entry['accession']] = 1
                else:
                    dico_position_entry[entry['accession']] += 1
                i += 1
        return(dico_position_entry)

    def make_list_position_entry(self, hits_per_gene, dico_position_entry, list_position_ribosomal, matchinggenes):
        list_position_entry = []
        j = 0
        len_position = len(list_position_ribosomal)
        key_value = sorted(dico_position_entry,
                           key=lambda x: dico_position_entry[x], reverse=True)
        for key in key_value:
            # parameter number of individue in the first matrix "j <= 5 " and appear in minimun 50% of blast
            if j <= 5 and dico_position_entry[key] > (len_position * 0.5):
                list_position_entry.append(key)
                j += 1
                #print(self.taxonomy_dict[key]['t_species'] ,key,dico_position_entry[key])

        list_position_ribosomal = self.make_list_position_ribosomal(
            matchinggenes, hits_per_gene, list_position_entry)
        return(list_position_entry, list_position_ribosomal)

    def make_list_position_ribosomal(self, matchinggenes, hits_per_gene, list_position_entry):
        list_position_ribosomal = []
        for hit in matchinggenes:
            k = 0
            if not list_position_entry:
                list_position_ribosomal.append(hit)
            else:
                for entry in hits_per_gene[hit]:
                    if entry['accession'] in list_position_entry:
                        k += 1
                # parameter, k number of
                if k >= len(list_position_entry):
                    list_position_ribosomal.append(hit)

        return(list_position_ribosomal)

    def creat_list_position(self, hits_per_gene, matchinggenes):
        dico_position_entry = self.make_dico_position_entry(
            hits_per_gene, matchinggenes)
        list_position_entry = []
        list_position_ribosomal = self.make_list_position_ribosomal(
            matchinggenes, hits_per_gene, list_position_entry)
        list_position_entry, last_list_position_ribosomal = self.make_list_position_entry(
            hits_per_gene, dico_position_entry, list_position_ribosomal, matchinggenes)

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
        # subsitution will impact all items on the index in lists
        while i < maxposition:
            emptylist = self.emptylistmaker(numbr_genes)
            matrix.append(emptylist)
            i += 1
        return(matrix)

    def fill_matrix(self, name, emptymatrix, hits_per_gene, list_position_entry,
                    list_position_ribosomal, matchinggenes):
        maxrank = 1
        bestident = 0
        for hit in list_position_ribosomal:
            rank = 1
            bestscore = 0
            perfectident = False
            for entry in hits_per_gene[hit]:
                if entry['accession'] in list_position_entry:
                    if rank > maxrank:
                        maxrank = rank
                    # parameter
                    if entry['pident'] > bestident:
                        bestident = entry['pident']
                    if perfectident and entry['pident'] != 100:
                        rank += 1
                        emptymatrix = self.fill_position_matrix(
                            emptymatrix, list_position_entry, list_position_ribosomal, rank, entry, hit)
                        continue
                    if bestscore < entry['bitscore']:
                        bestscore = entry['bitscore']
                        if entry['pident'] == 100:
                            perfectident = True
                        emptymatrix = self.fill_position_matrix(
                            emptymatrix, list_position_entry, list_position_ribosomal, rank, entry, hit)
                        continue
                    # parameter for considere 2 hit have same rank rank
                    if entry['pident'] == 100 or float(entry['bitscore']) >= (bestscore * 1):
                        emptymatrix = self.fill_position_matrix(
                            emptymatrix, list_position_entry, list_position_ribosomal, rank, entry, hit)
                        continue
                    rank += 1
                    emptymatrix = self.fill_position_matrix(
                        emptymatrix, list_position_entry, list_position_ribosomal, rank, entry, hit)

        return(emptymatrix, maxrank, bestident)

    def fill_position_matrix(self, emptymatrix, list_position_entry, list_position_ribosomal, rank, entry, hit):
        emptymatrix[list_position_entry.index(
            entry['accession'])][list_position_ribosomal.index(hit)] = rank
        return(emptymatrix)

    def fill_NA_matrix(self, matrix, penality, max_target_seqs=20):
        for liste in matrix:
            for n, i in enumerate(liste):
                if str(i) == 'NA':
                    liste[n] = (int(penality) + int(max_target_seqs))
        return(matrix)

    @timer
    def make_rank_matrix(self, name, hits_per_gene):
        matchinggenes = self.get_matching_gene(hits_per_gene)
        list_position_entry, list_position_ribosomal = self.creat_list_position(
            hits_per_gene, matchinggenes)
        emptymatrix = self.make_emptymatrix(
            list_position_entry, list_position_ribosomal)
        matrix, maxrank, bestident = self.fill_matrix(name, emptymatrix, hits_per_gene, list_position_entry,
                                                      list_position_ribosomal, matchinggenes)
        final_matrix = self.fill_NA_matrix(matrix, maxrank)

        return(final_matrix, list_position_entry, list_position_ribosomal, bestident)

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

    def rank_assignement(self, hits_per_gene, name):

        matrix, matrixlist_position_entry, list_position_ribosomal, bestident = self.make_rank_matrix(
            name, hits_per_gene)
        if anvio.DEBUG:
            self.show_matrix_rank(
                name, matrix, matrixlist_position_entry, list_position_ribosomal)
        start_make_list_taxonomy = time.perf_counter()

        taxonomy = self.make_list_taxonomy(
            matrix, matrixlist_position_entry, bestident)
        taxonomy_reduce = self.reduce_assignation_level(taxonomy, bestident)
        assignation = self.assign_taxonomie_solo_hit(taxonomy_reduce)
        return(assignation)
        end_make_list_taxonomy = time.perf_counter()
        print("\n time make_list_taxonomy : ",
              end_make_list_taxonomy - start_make_list_taxonomy)

        """posthoc_conover=sp.posthoc_conover_friedman(matrix)
        print("posthoc_conover")
        print(posthoc_conover)
        posthoc_miller_friedman=sp.posthoc_miller_friedman(matrix)
        print("posthoc_miller_friedman")
        print(posthoc_miller_friedman)
        posthoc_siegel_friedman=sp.posthoc_siegel_friedman(matrix)
        print("posthoc_siegel_friedman")
        print(posthoc_siegel_friedman)
        posthoc_durbin=sp.posthoc_durbin(matrix)
        print("posthoc_durbin")
        print(posthoc_durbin)

        posthoc_nemenyi_friedman=sp.posthoc_nemenyi_friedman(matrix)
        print("posthoc_nemenyi_friedman")
        print(posthoc_nemenyi_friedman)"""

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
        # print(taxonomy)
        assignation = taxonomy[0]
        for s in taxonomy[1:]:
            for key in s:
                if s[key] not in assignation.values():
                    assignation.pop(key, None)
        self.run.info('estimate taxonomy',
                      '/'.join(list(assignation.values())))
        return(assignation)

    def make_dicoidentitielevel(self):
        dicoidentitielevel = {}
        dicoidentitielevel[95] = 't_species'
        dicoidentitielevel[90] = 't_genus'
        dicoidentitielevel[85] = 't_family'
        dicoidentitielevel[80] = 't_order'
        return(dicoidentitielevel)

    def make_list_taxonomy(self, matrix, list_position_entry, bestident):
        taxonomy = []
        miniscore = 10000
        i = 0
        for liste in matrix:
            summist = sum(list(liste))
            if summist == miniscore:
                taxonomy.append(self.taxonomy_dict[list_position_entry[i]])
            if summist < miniscore:
                taxonomy = []
                taxonomy.append(self.taxonomy_dict[list_position_entry[i]])
                miniscore = summist
            i += 1
        return(taxonomy)

    def reduce_assignation_level(self, taxonomy, bestident):
        dicoidentitielevel = self.make_dicoidentitielevel()
        for key, value in dicoidentitielevel.items():
            if bestident < int(key):
                for possibilitie in taxonomy:
                    possibilitie.pop(value, None)
        return(taxonomy)

    def assignation_list(self, taxonomy):
        if len(taxonomy) == 1:
            self.run.info('estimate taxonomy rank', taxonomy[0])
        if len(taxonomy) >= 2:
            for possibilities in taxonomy:
                self.run.info('Relevant match', possibilities)

            self.assign_taxonomie_solo_hit(taxonomy)

    def solo_hits(self, hits_per_gene):
        taxonomy = []
        for hit in hits_per_gene:
            for entry in hits_per_gene[hit]:
                if len(hits_per_gene[hit]) == 1 and entry['pident'] > 0.90:
                    taxonomy.append(self.taxonomy_dict[entry['accession']])
        if taxonomy:
            self.assign_taxonomie_solo_hit(taxonomy)
            return True
        else:
            return False
