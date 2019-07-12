
#!/usr/bin/env python3
# -*- coding: utf-8

import argparse
import os
import re
import sys
import time
from collections import Counter

import anvio

import anvio.db as db
import anvio.ccollections as ccollections
import anvio.fastalib as f
import anvio.filesnpaths as filesnpaths
import anvio.hmmops as hmmops
import anvio.hmmopswrapper as hmmopswrapper
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.tables.tableops import Table
from anvio.dbops import ContigsSuperclass
from anvio.drivers import Aligners, driver_modules
from anvio.drivers.diamond import Diamond
from anvio.errors import ConfigError, FilesNPathsError
from tabulate import tabulate


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
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


        self.num_threads=args.num_threads

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


        self.blast_hits_table_name                    = 'blast_hits'
        self.blast_hits_table_structure               = ['match_id','bin_id',  'gene_callers_id', 'gene_name', 'taxon_id', 'pourcentage_identity', 'bitscore']
        self.blast_hits_table_types                   = ['text'     , 'text'     ,   'text'    ,      'text'   ,   'text'   ,     'text'   ,         'text']

        self.taxon_names_table_name                 = 'taxon_names'
        self.taxon_names_table_structure            = ['taxon_id', 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
        self.taxon_names_table_types                = [ 'numeric UNIQUE',   'text',   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        self.database.drop_table(self.blast_hits_table_name)
        self.database.drop_table(self.taxon_names_table_name)

        self.database.create_table(self.blast_hits_table_name, self.blast_hits_table_structure, self.blast_hits_table_types)
        self.database.create_table(self.taxon_names_table_name, self.taxon_names_table_structure, self.taxon_names_table_types)


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

        self.taxonomy_dict={}


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


    @timer
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
        hmm_sequences_dict_per_type = self.get_hmm_sequences_dict_into_type_multi(
            hmm_sequences_dict)
        j=0
        hits_per_gene = {}
        start_blast_sort = time.perf_counter()
        match_id=0
        for SCG in hmm_sequences_dict_per_type:



            self.run.info('SCGs', SCG)
            self.run.info('Num sequence', len(hmm_sequences_dict_per_type[SCG]))

            """self.run.info('SCGs', ', '.join(
                [e['gene_SCG'] for e in hmm_sequences_dict_per_type[SCG].values()]))"""



            hmm_sequences_dict_per_type[SCG],match_id = self.get_raw_blast_hits_multi(hmm_sequences_dict_per_type[SCG],match_id)


            start_get_raw_blast_hits = time.perf_counter()


            if self.metagenome:
                var='gene_callers_id'
                possibles_taxonomy=[]
            else:
                var='bin_id'
            for query in hmm_sequences_dict_per_type[SCG].values():

                if anvio.DEBUG:
                    self.show_hits(query[var], SCG, query['hits'])
                if query[var] not in hits_per_gene:
                    hits_per_gene[query[var]]={}
                if SCG not in hits_per_gene[query[var]]:
                    hits_per_gene[query[var]][SCG]=[]
                if len(query['hits']):
                    hits_per_gene[query[var]][SCG] = query['hits']
                else:
                    hits_per_gene[query[var]][SCG] = hits_per_gene[query[var]][SCG] + query['hits']


        self.database.disconnect()








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
            diamond.num_threads = self.num_threads

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
    def get_raw_blast_hits_multi(self, d,match_id, max_target_seqs=20, evalue=1e-05, min_pct_id=90):

        sequence=""
        bin_dict_id={}
        for id, entry in d.items():
            if 'sequence' not in entry or 'gene_name' not in entry:
                raise ConfigError("The `get_filtered_dict` function got a parameter that\
                                   does not look like the way we expected it. This function\
                                   expects a dictionary that contains keys `gene_name` and `sequence`.")

                #sequence = sequence+">%s\n%s\n"% (id,entry['sequence'])
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
            #print(taxo)

            try:
                self.database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % "taxon_names", taxo)
            except:
                pass

            #dict(zip(['accession', 'pident', 'length', 'mismatch', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], [float(entry[i]) if i > 1 else entry[i] for i in range(1, 12)]))
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

        self.num_threads=args.num_threads

        if not args.num_threads:
            self.num_threads="2"

        self.initialized = False

        self.metagenome=False

        if args.metagenome:
            self.metagenome=True

        self.cut_off_methode = args.cut_off_methode

        self.methode = args.methode

        self.taxonomy_dict = {}

        #self.profile_db = args.profile_db

    def estimate_taxonomy(self):

        self.database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))



        """A=self.database.get_table_as_dataframe('blast_hits')
        B=self.database.get_table_as_dataframe('taxon_names')
        print(A)
        print(B)
        """

        A=self.database.get_table_as_dict('blast_hits')
        self.taxonomy_dict =self.database.get_table_as_dict('taxon_names')





        #print(A)
        if self.metagenome:
            var='gene_callers_id'
            possibles_taxonomy=[]
        else:
            var='bin_id'
        hits_per_gene={}
        for query in A.values():
            hit=[{'accession':query['taxon_id'], 'pident':float(query['pourcentage_identity']), 'bitscore': float(query['bitscore'])}]



            if query[var] not in hits_per_gene:
                hits_per_gene[query[var]]={}
            if query['gene_name'] not in hits_per_gene[query[var]]:
                hits_per_gene[query[var]][query['gene_name']]=[]

            hits_per_gene[query[var]][query['gene_name']] = hits_per_gene[query[var]][query['gene_name']] + hit








        for name, SCGs_hit_per_gene in hits_per_gene.items():

            if not self.metagenome:
                self.run.info('Taxo for ', name)


            taxonomy = self.get_consensus_taxonomy(
                SCGs_hit_per_gene, name)

            if not taxonomy:
                self.run.info('taxonomy estimation not possible for:', name)
                continue


            if not self.metagenome:
                self.run.info('estimate taxonomy',
                              '/'.join(list(taxonomy.values())))


            if self.metagenome and str(list(taxonomy.values())[-1]) not in possibles_taxonomy:
                possibles_taxonomy.append(str(list(taxonomy.values())[-1]))

        if self.metagenome:
            self.run.info('Possible presence ','|'.join(list(possibles_taxonomy)))



        """start_predict_from_SCGs_dict = time.perf_counter()

        end_predict_from_SCGs_dict = time.perf_counter()
        print("\n time predict_from_SCGs_dict for a bin : ",
              end_predict_from_SCGs_dict - start_predict_from_SCGs_dict)"""



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

    @timer
    def get_consensus_taxonomy(self, hits_per_gene, name):
        """Different methode for assignation"""

        consensus_taxonomy=self.solo_hits(hits_per_gene)

        if consensus_taxonomy:

            return(consensus_taxonomy)

        else:

            if self.methode == "friedman":
                consensus_taxonomy = self.rank_assignement(hits_per_gene, name)

            """if self.methode == "tree":
                self.make_tree_with_hit(
                    hits_per_gene, hmm_sequences_dict_per_type, name)
                consensus_taxonomy = "tree"
                return(consensus_taxonomy)"""

            if self.methode == "bitscore":
                score_by_entry = self.get_cumul_hit_per_gene(hits_per_gene)
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

    def get_cumul_hit_per_gene(self, hits_per_gene):
        """add bitscore of eatch  per query"""
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
            if j< 10 and dico_position_entry[key] > (len_position * 0.5):
                list_position_entry.append(key)
                j += 1
                #print(self.taxonomy_dict[key]['t_species'] ,key,dico_position_entry[key])

        list_position_ribosomal = self.make_list_position_ribosomal(
            matchinggenes, hits_per_gene, list_position_entry)
        if not len(list_position_entry) or not len(list_position_ribosomal) :
            sys.exit(status=None)
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
        bestSCG=None
        emptymatrix_ident=emptymatrix
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
                        bestSCG=hit
                    if perfectident and entry['pident'] != 100:

                        rank += 1

                        emptymatrix, emptymatrix_ident = self.fill_position_matrix(
                            emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, entry['pident'], rank, entry, hit)

                        continue

                    if bestscore < entry['bitscore']:
                        bestscore = entry['bitscore']

                        if entry['pident'] == 100:
                            perfectident = True

                        emptymatrix, emptymatrix_ident = self.fill_position_matrix(
                            emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, entry['pident'], rank, entry, hit)

                        continue
                    # parameter for considere 2 hit have same rank rank
                    if entry['pident'] == 100 or float(entry['bitscore']) >= (bestscore * 1):
                        emptymatrix, emptymatrix_ident = self.fill_position_matrix(
                            emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, entry['pident'], rank, entry, hit)

                        continue
                    rank += 1
                    emptymatrix, emptymatrix_ident = self.fill_position_matrix(
                        emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, entry['pident'], rank, entry, hit)


        return(emptymatrix, emptymatrix_ident, maxrank, bestident,bestSCG)


    def fill_position_matrix(self, emptymatrix, emptymatrix_ident, list_position_entry, list_position_ribosomal, pident, rank, entry, hit):
        emptymatrix[list_position_entry.index(
            entry['accession'])][list_position_ribosomal.index(hit)] = rank
        emptymatrix_ident[list_position_entry.index(
            entry['accession'])][list_position_ribosomal.index(hit)] = pident
        return(emptymatrix,emptymatrix_ident)

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

        matrix, matrix_pident, maxrank, bestident, bestSCG = self.fill_matrix(name, emptymatrix, hits_per_gene, list_position_entry,
                                                      list_position_ribosomal, matchinggenes)
        if anvio.DEBUG:
            self.show_matrix_rank(
                name, matrix, list_position_entry, list_position_ribosomal)

            self.show_matrix_rank(
                name, matrix_pident, list_position_entry, list_position_ribosomal)

        final_matrix = self.fill_NA_matrix(matrix, maxrank)

        return(final_matrix, matrix_pident, list_position_entry, list_position_ribosomal, bestident, bestSCG)

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

        matrix, matrix_pident, matrixlist_position_entry, list_position_ribosomal, bestident, bestSCG = self.make_rank_matrix(
            name, hits_per_gene)
        """if anvio.DEBUG:
            self.show_matrix_rank(
                name, matrix, matrixlist_position_entry, list_position_ribosomal)"""
        start_make_list_taxonomy = time.perf_counter()

        taxonomy,bestlist = self.make_list_taxonomy(
            matrix, matrixlist_position_entry,list_position_ribosomal)

        print(bestlist)
        assignation = self.assign_taxonomie_solo_hit(taxonomy)
        #assignation_reduce = self.reduce_assignation_level(matrix_pident, assignation, bestlist)
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

    def del_higher_ranked_individu_matrix(self, matrix, matrixlist_position_entry,):
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
                        self.run.info('CONSIDER taxonomy','/'.join(assignation.values()))
                    assignation.pop(key, None)
        return(assignation)

    def make_dicoidentitielevel(self):
        dicoidentitielevel = {}
        dicoidentitielevel[95] = 't_species'
        dicoidentitielevel[90] = 't_genus'
        dicoidentitielevel[85] = 't_family'
        dicoidentitielevel[80] = 't_order'
        return(dicoidentitielevel)

    def make_list_taxonomy(self, matrix, list_position_entry,list_position_ribosomal):
        taxonomy = []
        bestlist=[]
        miniscore = 10000
        i = 0
        for liste in matrix:
            summist = sum(list(liste))
            bestident = max(list(liste))
            bestSCG=list_position_ribosomal[liste.index(bestident)]
            if summist == miniscore:
                taxonomy.append(self.taxonomy_dict[list_position_entry[i]])
                bestlist.append([bestSCG,bestident])
            if summist < miniscore:
                taxonomy = []

                taxonomy.append(self.taxonomy_dict[list_position_entry[i]])
                miniscore = summist
                bestlist=[[bestSCG,bestident]]
            i += 1
        return(taxonomy,bestlist)

    def reduce_assignation_level(self, taxonomy, bestident, bestSCG):
        dicoidentitielevel = self.make_dicoidentitielevel()
        for key, value in dicoidentitielevel.items():
            if bestident < int(key):
                for possibilitie in taxonomy:
                    possibilitie.pop(value, None)
        return(taxonomy)

    def low_ident(self,matrix_pident,taxonomy):
        "lol"


    def reduce_assignation_level(self,matrix_pident, taxonomy, bestident, bestSCG):
        for possibilitie in taxonomy:
            while bestident < dico[bestSCG][list(taxonomy.values())[-1]]:
                possibilitie.pop(list(taxonomy.key())[-1], None)
        return(taxonomy)


    def solo_hits(self, hits_per_gene):
        taxonomy = []
        for hit in hits_per_gene:
            for entry in hits_per_gene[hit]:
                if len(hits_per_gene[hit]) == 1 and entry['pident'] > 0.90:
                    taxonomy.append(self.taxonomy_dict[entry['accession']])
        if taxonomy:
            assignation=self.assign_taxonomie_solo_hit(taxonomy)
            return assignation
        else:
            return False
