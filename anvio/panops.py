# -*- coding: utf-8
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import copy
import shutil
import hashlib

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.ccollections as ccollections
import anvio.clustering as clustering
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.drivers.blast import BLAST
from anvio.drivers.diamond import Diamond
from anvio.drivers.mcl import MCL

from anvio.errors import ConfigError, FilesNPathsError

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Pangenome:
    def __init__(self, args = None, run = run, progress = progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        input_file_for_internal_genomes = A('internal_genomes')
        input_file_for_external_genomes = A('external_genomes')
        self.num_threads = A('num_threads')
        self.output_dir = A('output_dir')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')
        self.min_percent_identity = A('min_percent_identity')
        self.PC_min_occurrence = A('min_occurrence')
        self.mcl_inflation = A('mcl_inflation') 
        self.sensitive = A('sensitive') 
        self.maxbit = A('maxbit') 
        self.use_ncbi_blast = A('use_ncbi_blast') 

        self.genomes = {}

        fields_for_internal_genomes_input = ['name', 'bin_id', 'collection_id', 'profile_db_path', 'contigs_db_path']
        fields_for_external_genomes_input = ['name', 'contigs_db_path']

        self.log_file_path = None

        internal_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_for_internal_genomes, expected_fields = fields_for_internal_genomes_input) if input_file_for_internal_genomes else {}
        external_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_for_external_genomes, expected_fields = fields_for_external_genomes_input) if input_file_for_external_genomes else {}

        self.internal_genome_names = internal_genomes_dict.keys()
        self.external_genome_names = external_genomes_dict.keys()

        if len(self.internal_genome_names) + len(self.external_genome_names) != len(set(self.internal_genome_names + self.external_genome_names)):
            raise ConfigError, "Each entry both in internal and external genome descriptions should have a unique 'name'. This does not\
                                seem to be the case with your input :/"

        # convert relative paths to absolute paths and MERGE internal and external genomes into self.genomes:
        for source, input_file in [(external_genomes_dict, input_file_for_external_genomes), (internal_genomes_dict, input_file_for_internal_genomes)]:
            for genome_name in source:
                self.genomes[genome_name] = source[genome_name]
                for db_path_var in ['contigs_db_path', 'profile_db_path']:
                    if db_path_var not in self.genomes[genome_name]:
                        continue
                    path = self.genomes[genome_name][db_path_var]
                    if not path.startswith('/'):
                        self.genomes[genome_name][db_path_var] = os.path.abspath(os.path.join(os.path.dirname(input_file), path))

        # to be filled during init:
        self.hash_to_genome_name = {}
        self.view_data = {}
        self.view_data_presence_absence = {}
        self.additional_view_data = {}


    def get_output_file_path(self, file_name):
        output_file_path = os.path.join(self.output_dir, file_name)

        return output_file_path


    def check_programs(self):
        utils.is_program_exists('diamond')
        utils.is_program_exists('mcl')


    def check_params(self):
        # deal with the output directory:
        try:
            filesnpaths.is_file_exists(self.output_dir)
        except FilesNPathsError:
            filesnpaths.gen_output_directory(self.output_dir, delete_if_exists = self.overwrite_output_destinations)

        filesnpaths.is_output_dir_writable(self.output_dir)
        self.output_dir = os.path.abspath(self.output_dir)

        if not self.log_file_path:
            self.log_file_path = self.get_output_file_path('log.txt')
        
        filesnpaths.is_output_file_writable(self.log_file_path)
        os.remove(self.log_file_path) if os.path.exists(self.log_file_path) else None

        if type(self.maxbit) != float:
            raise ConfigError, "maxbit value must be of type float :("

        if self.maxbit < 0 or self.maxbit > 1:
            raise ConfigError, "Well. maxbit must be between 0 and 1. Yes. Very boring."

        if type(self.min_percent_identity) != float:
            raise ConfigError, "Minimum percent identity value must be of type float :("

        if self.min_percent_identity < 0 or self.min_percent_identity > 100:
            raise ConfigError, "Minimum percent identity must be between 0%% and 100%%. Although your %.2f%% is\
                                pretty cute, too." % self.min_percent_identity


        if len([c for c in self.genomes.values() if 'contigs_db_path' not in c]):
            raise ConfigError, "self.genomes does not seem to be a properly formatted dictionary for\
                                the anvi'o class Pangenome."

        for genome_name in self.genomes:
            if not os.path.exists(self.genomes[genome_name]['contigs_db_path']):
                raise ConfigError, "The contigs database for genome %s is not where the input data suggested where\
                                    it would be.." % genome_name
            if genome_name in self.internal_genome_names and not os.path.exists(self.genomes[genome_name]['profile_db_path']):
                raise ConfigError, "The profile database for genome %s is not where the input data suggested where\
                                    it would be.." % genome_name


    def init_external_genomes(self):
        self.progress.new('Initializing external genomes')
        for genome_name in self.external_genome_names:
            c = self.genomes[genome_name]
            c['name'] = genome_name

            self.progress.update('working on %s' % (genome_name))

            contigs_db_summary = summarizer.get_contigs_db_info_dict(c['contigs_db_path'])

            for key in contigs_db_summary:
                c[key] = contigs_db_summary[key]

            c['genome_entry_hash'] = c['contigs_db_hash']

            self.hash_to_genome_name[c['genome_entry_hash']] = genome_name
        self.progress.end()

        # if two contigs db has the same hash, we are kinda f'd:
        if len(set([self.genomes[genome_name]['genome_entry_hash'] for genome_name in self.external_genome_names])) != len(self.external_genome_names):
            raise ConfigError, 'Not all hash values are unique across all contig databases you provided. Something\
                                very fishy is going on :/'

        # make sure genes are called in every contigs db:
        genomes_missing_gene_calls = [g for g in self.external_genome_names if not self.genomes[genome_name]['genes_are_called']]
        if len(genomes_missing_gene_calls):
            raise ConfigError, 'Genes must have been called during the generation of contigs database for this workflow to work. However,\
                                these external genomes do not have gene calls: %s' % (', '.join(genomes_missing_gene_calls))

        self.run.info('External genomes', '%d have been initialized.' % len(self.external_genome_names))


    def init_internal_genomes(self):
        self.progress.new('Initializing internal genomes')

        # to not initialize things over and over again:
        unique_profile_db_path_to_internal_genome_name = {}
        for profile_path in set([self.genomes[g]['profile_db_path'] for g in self.internal_genome_names]):
            unique_profile_db_path_to_internal_genome_name[profile_path] = [g for g in self.internal_genome_names if self.genomes[g]['profile_db_path'] == profile_path]

        for profile_db_path in unique_profile_db_path_to_internal_genome_name:
            self.collections = ccollections.Collections()
            self.collections.populate_collections_dict(profile_db_path, anvio.__profile__version__)

            for genome_name in unique_profile_db_path_to_internal_genome_name[profile_db_path]:
                self.progress.update('working on %s' % (genome_name))
                c = self.genomes[genome_name]

                dbops.is_profile_db_and_contigs_db_compatible(c['profile_db_path'], c['contigs_db_path'])

                # set name
                c['name'] = genome_name

                collection_dict = self.collections.get_collection_dict(c['collection_id'])
                bins_info_dict = self.collections.get_bins_info_dict(c['collection_id'])

                if c['bin_id'] not in bins_info_dict:
                    self.progress.end()
                    raise ConfigError, "You betrayed us :( Genome %s does not appear to be a valid bin in collection %s in %s"\
                                % (c['bin_id'], c['collection_id'], c['profile_db_path'])


                split_names_of_interest = collection_dict[c['bin_id']]
                if not len(split_names_of_interest):
                    raise ConfigError, "There are 0 splits defined for bin id %s in collection %s..." % (c['bin_id'], c['collection_id'])


                contigs_db_summary = summarizer.get_contigs_db_info_dict(c['contigs_db_path'], split_names = split_names_of_interest, include_AA_counts = True)
                for key in contigs_db_summary:
                    c[key] = contigs_db_summary[key]

                # set hash
                c['genome_entry_hash'] = hashlib.sha224('_'.join([split_names_of_interest[0], split_names_of_interest[-1], c['contigs_db_hash']])).hexdigest()
                self.hash_to_genome_name[c['genome_entry_hash']] = genome_name

        self.progress.end()

        if len(set([self.genomes[genome_name]['genome_entry_hash'] for genome_name in self.internal_genome_names])) != len(self.internal_genome_names):
            raise ConfigError, "Not all hash values are unique across internal genomes. This is almost impossible to happen unless something very\
                                wrong with your workflow :/ Please let the developers know if you can't figure this one out"

        # make sure genes are called in every contigs db:
        genomes_missing_gene_calls = [g for g in self.internal_genome_names if not self.genomes[genome_name]['genes_are_called']]
        if len(genomes_missing_gene_calls):
            raise ConfigError, 'Genes must have been called during the generation of contigs database for this workflow to work. However,\
                                these external genomes do not have gene calls: %s' % (', '.join(genomes_missing_gene_calls))

        self.run.info('Internal genomes', '%d have been initialized.' % len(self.internal_genome_names))


    def init_genomes(self):
        self.init_external_genomes()
        self.init_internal_genomes()


    def gen_combined_proteins_unique_FASTA(self):
        self.progress.new('Storing combined protein sequences')
        combined_proteins_FASTA_path = self.get_output_file_path('combined-proteins.fa')
        output_file = open(combined_proteins_FASTA_path, 'w')

        for genome_name in self.genomes:
            g = self.genomes[genome_name]
            self.progress.update('Working on %s ...' % genome_name)
            contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
            protein_sequences_dict = contigs_db.db.get_table_as_dict(t.gene_protein_sequences_table_name)
            for gene_caller_id in g['gene_caller_ids']:
                output_file.write('>%s_%d\n' % (g['genome_entry_hash'], gene_caller_id))
                output_file.write('%s\n' % protein_sequences_dict[gene_caller_id]['sequence'])
            contigs_db.disconnect()

        output_file.close()
        self.progress.end()

        # unique the FASTA file
        unique_proteins_FASTA_path, unique_proteins_names_file_path, unique_proteins_names_dict = utils.unique_FASTA_file(combined_proteins_FASTA_path, store_frequencies_in_deflines = False)

        self.run.info('Num protein sequences', '%s' % pp(sum([g['num_genes'] for g in self.genomes.values()])))
        self.run.info('Num unique protein sequences', '%s' % pp(len(unique_proteins_names_dict)))
        self.run.info('Combined protein sequences FASTA', combined_proteins_FASTA_path)
        self.run.info('Unique protein sequences FASTA', unique_proteins_FASTA_path)

        return unique_proteins_FASTA_path, unique_proteins_names_dict


    def run_diamond(self, unique_proteins_fasta_path, unique_proteins_names_dict):
        diamond = Diamond(unique_proteins_fasta_path, run = self.run, progress = self.progress,
                          num_threads = self.num_threads, overwrite_output_destinations = self.overwrite_output_destinations)

        diamond.names_dict = unique_proteins_names_dict
        diamond.target_db_path = self.get_output_file_path('.'.join(unique_proteins_fasta_path.split('.')[:-1]))
        diamond.search_output_path = self.get_output_file_path('diamond-search-results')
        diamond.tabular_output_path = self.get_output_file_path('diamond-search-results.txt')

        diamond.sensitive = self.sensitive

        return diamond.get_blastall_results()


    def run_blast(self, unique_proteins_fasta_path, unique_proteins_names_dict):
        self.run.warning("You elected to use NCBI's blastp for protein search. Running blastp will be significantly\
                          slower than DIAMOND (although, anvi'o developers are convinced that you *are*\
                          doing the right thing, so, kudos to you).")
        blast = BLAST(unique_proteins_fasta_path, run = self.run, progress = self.progress,
                          num_threads = self.num_threads, overwrite_output_destinations = self.overwrite_output_destinations)

        blast.names_dict = unique_proteins_names_dict
        blast.log_file_path = self.log_file_path
        blast.target_db_path = self.get_output_file_path('.'.join(unique_proteins_fasta_path.split('.')[:-1]))
        blast.search_output_path = self.get_output_file_path('blast-search-results.txt')


        return blast.get_blastall_results()


    def run_search(self, unique_proteins_fasta_path, unique_proteins_names_dict):
        if self.use_ncbi_blast:
            return self.run_blast(unique_proteins_fasta_path, unique_proteins_names_dict)
        else:
            return self.run_diamond(unique_proteins_fasta_path, unique_proteins_names_dict)


    def run_mcl(self, mcl_input_file_path):
        mcl = MCL(mcl_input_file_path, run = self.run, progress = self.progress, num_threads = self.num_threads)

        mcl.inflation = self.mcl_inflation
        mcl.clusters_file_path = self.get_output_file_path('mcl-clusters.txt')
        mcl.log_file_path = self.log_file_path

        return mcl.get_clusters_dict()


    def gen_mcl_input(self, blastall_results):
        self.progress.new('Processing search results')
        self.progress.update('...')

        all_ids = set([])

        # mapping for the fields in the blast output
        mapping = [str, str, float, int, int, int, int, int, int, int, float, float]

        # here we perform an initial pass on the blast results to fill the dict that will hold
        # the bit score for each gene when it was blasted against itself. this dictionary
        # will then be used to calculate the 'maxbit' value between two genes, which I learned
        # from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8). ITEP defines maxbit as 
        # 'bit score between target and query / min(selfbit for query, selbit for target)'. This
        # heuristic approach provides a mean to set a cutoff to eliminate weak matches between
        # two genes. maxbit value reaches to 1 for hits between two genes that are almost identical.
        self_bit_scores = {}
        line_no = 1
        self.progress.update('(initial pass of the serach results to set the self bit scores ...)')
        for line in open(blastall_results):
            fields = line.strip().split('\t')

            try:
                query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                    [mapping[i](fields[i]) for i in range(0, len(mapping))]
            except Exception, e:
                self.progress.end()
                raise ConfigError, "Something went wrong while processing the blastall output file in line %d.\
                                    Here is the error from the uppoer management: '''%s'''" % (line_no, e)
            line_no += 1
            all_ids.add(query_id)
            all_ids.add(subject_id)

            if query_id == subject_id:
                self_bit_scores[query_id] = bit_score

        self.progress.end()

        ids_without_self_search = all_ids - set(self_bit_scores.keys())
        if len(ids_without_self_search):
            search_tool = 'BLAST' if self.use_ncbi_blast else 'DIAMOND'
            self.run.warning("%s did not retun search results for %d of %d the protein sequences in your input FASTA file.\
                              Anvi'o will do some heuristic magic to complete the missing data in the search output to recover\
                              from this. But since you are a scientist, here are the protein sequence IDs for which %s\
                              failed to report self search results: %s." \
                                                    % (search_tool, len(ids_without_self_search), len(all_ids), \
                                                       search_tool, ', '.join(ids_without_self_search)))

        # HEURISTICS TO ADD MISSING SELF SEARCH RESULTS
        # we are here, because protein sequences in ids_without_self_search did not have any hits in the search output
        # although they were in the FASTA file the target database were built from. so we will make sure they are not
        # missing from self_bit_scores dict, or mcl_input (additional mcl inputs will be stored in the following dict)
        additional_mcl_input_lines = {}
        for id_without_self_search in ids_without_self_search:
            entry_hash, gene_caller_id = id_without_self_search.split('_')

            try:
                genome_name = self.hash_to_genome_name[entry_hash]
            except KeyError:
                raise ConfigError, "Something horrible happened. This can only happend if you started a new analysis with\
                                    additional genomes without cleaning the previous work directory. Sounds familiar?"

            # divide the DNA length of the gene by three to get the AA length, and multiply that by two to get an approximate
            # bit score that would have recovered from a perfect match
            self_bit_scores[id_without_self_search] = (self.genomes[genome_name]['gene_lengths'][int(gene_caller_id)] / 3.0) * 2

            # add this SOB into additional_mcl_input_lines dict.
            additional_mcl_input_lines[id_without_self_search] = '%s\t%s\t1.0\n' % (id_without_self_search, id_without_self_search)


        # CONTINUE AS IF NOTHING HAPPENED
        self.run.info('Min percent identity', self.min_percent_identity)
        self.run.info('Maxbit', self.maxbit)
        self.progress.new('Processing search results')

        mcl_input_file_path = self.get_output_file_path('mcl-input.txt')
        mcl_input = open(mcl_input_file_path, 'w')

        line_no = 1
        num_edges_stored = 0
        for line in open(blastall_results):
            fields = line.strip().split('\t')

            query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                [mapping[i](fields[i]) for i in range(0, len(mapping))]

            line_no += 1

            if line_no % 5000 == 0:
                self.progress.update('Lines processed %s ...' % pp(line_no))

            #
            # FILTERING BASED ON PERCENT IDENTITY
            #
            if perc_id < self.min_percent_identity:
                continue

            #
            # FILTERING BASED ON MAXBIT
            #
            maxbit = bit_score / min(self_bit_scores[query_id], self_bit_scores[subject_id])
            if maxbit < self.maxbit:
                continue

            mcl_input.write('%s\t%s\t%f\n' % (query_id, subject_id, perc_id / 100.0))
            num_edges_stored += 1

        # add additional lines if there are any:
        for line in additional_mcl_input_lines.values():
            mcl_input.write(line)
            num_edges_stored += 1

        mcl_input.close()

        self.progress.end()
        self.run.info('Filtered search results', '%s edges stored' % pp(num_edges_stored))
        self.run.info('MCL input', '%s' % mcl_input_file_path)

        return mcl_input_file_path


    def gen_data_from_protein_clusters(self, protein_clustering_dict):
        self.progress.new('Generating view data')
        self.progress.update('...')

        PCs = protein_clustering_dict.keys()

        for PC in PCs:
            self.view_data[PC] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.view_data_presence_absence[PC] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.additional_view_data[PC] = {'num_genes_in_pc': 0, 'num_genomes_pc_has_hits': 0}
            for entry_hash, gene_caller_id in [e.split('_') for e in protein_clustering_dict[PC]]:
                try:
                    genome_name = self.hash_to_genome_name[entry_hash]
                except KeyError:
                    raise ConfigError, "Something horrible happened. This can only happend if you started a new analysis with\
                                        additional genomes without cleaning the previous work directory. Sounds familiar?"
                self.view_data[PC][genome_name] += 1
                self.view_data_presence_absence[PC][genome_name] = 1
                self.additional_view_data[PC]['num_genes_in_pc'] += 1
            self.additional_view_data[PC]['num_genomes_pc_has_hits'] = len([True for genome in self.view_data[PC] if self.view_data[PC][genome] > 0])

        self.progress.end()

        #
        # FILTERING BASED ON OCCURRENCE
        #
        PCs_of_interest = set([])
        for PC in PCs:
            if self.additional_view_data[PC]['num_genomes_pc_has_hits'] >= self.PC_min_occurrence:
                PCs_of_interest.add(PC)

        for PC in PCs:
            if PC not in PCs_of_interest:
                self.view_data.pop(PC)
                self.view_data_presence_absence.pop(PC)
                self.additional_view_data.pop(PC)

        if self.PC_min_occurrence > 1:
            self.run.info('PCs min occurrence', '%d (the filter removed %s PCs)' % (self.PC_min_occurrence, (len(protein_clustering_dict) - len(PCs_of_interest))))

        view_data_file_path = self.get_output_file_path('anvio-view-data.txt')
        additional_view_data_file_path = self.get_output_file_path('anvio-additional-view-data.txt')
        view_data_presence_absence_file_path = self.get_output_file_path('anvio-view-data-prsence-absence.txt')

        utils.store_dict_as_TAB_delimited_file(self.view_data, view_data_file_path, headers = ['contig'] + sorted(self.genomes.keys()))
        utils.store_dict_as_TAB_delimited_file(self.additional_view_data, additional_view_data_file_path, headers = ['contig'] + sorted(self.additional_view_data.values()[0].keys()))
        utils.store_dict_as_TAB_delimited_file(self.view_data_presence_absence, view_data_presence_absence_file_path, headers = ['contig'] + sorted(self.view_data_presence_absence.values()[0].keys()))

        # here's where we finalize experimental data for clustering
        experimental_data = copy.deepcopy(self.view_data_presence_absence)
        for PC in self.additional_view_data:
            for i in range(0, int(len(self.genomes) / 2)):
                experimental_data[PC]['num_genomes_pc_has_hits_%d' % i] = self.additional_view_data[PC]['num_genomes_pc_has_hits']
        experimental_data_file_path = self.get_output_file_path('anvio-experimental-data-for-clustering.txt')
        utils.store_dict_as_TAB_delimited_file(experimental_data, experimental_data_file_path, headers = ['contig'] + sorted(experimental_data.values()[0].keys()))

        self.run.info("Anvi'o view data for protein clusters", view_data_file_path)
        self.run.info("Anvi'o additional view data", additional_view_data_file_path)

        return view_data_file_path, view_data_presence_absence_file_path, additional_view_data_file_path, experimental_data_file_path


    def gen_samples_info_file(self):
        samples_info_dict = {}
        samples_info_file_path = self.get_output_file_path('anvio-samples-information.txt')

        # set headers
        headers = ['total_length']

        for h in ['percent_complete', 'percent_redundancy']:
            if self.genomes.values()[0].has_key(h):
                headers.append(h)

        headers.extend(['gc_content', 'num_genes', 'avg_gene_length', 'num_genes_per_kb'])

        for c in self.genomes.values():
            new_dict = {}
            for header in headers:
                new_dict[header] = c[header]

            samples_info_dict[c['name']] = new_dict

        utils.store_dict_as_TAB_delimited_file(samples_info_dict, samples_info_file_path, headers = ['samples'] + headers)

        self.run.info("Anvi'o samples information", samples_info_file_path)

        return samples_info_file_path


    def gen_samples_order_file(self, view_data_file_path):
        self.progress.new('Hierarchical clustering of the (transposed) view data')
        self.progress.update('..')


        newick = clustering.get_newick_tree_data(view_data_file_path, transpose = True)

        samples_order_file_path = self.get_output_file_path('anvio-samples-order.txt')
        samples_order = open(samples_order_file_path, 'w')
        samples_order.write('attributes\tbasic\tnewick\n')
        samples_order.write('protein_clusters\t\t%s\n' % newick)
        samples_order.close()

        self.progress.end()

        self.run.info("Anvi'o samples order", samples_order_file_path)

        return samples_order_file_path


    def gen_ad_hoc_anvio_run(self, view_data_file_path, experimental_data_file_path, additional_view_data_file_path, samples_info_file_path, samples_order_file_path):
        ad_hoc_run = AdHocRunGenerator(view_data_file_path, run = self.run, progress = self.progress)

        ad_hoc_run.matrix_data_for_clustering = experimental_data_file_path
        ad_hoc_run.additional_view_data_file_path = additional_view_data_file_path
        ad_hoc_run.samples_info_file_path = samples_info_file_path
        ad_hoc_run.samples_order_file_path = samples_order_file_path

        ad_hoc_run.output_directory = self.get_output_file_path('anvio-run-files')
        ad_hoc_run.delete_output_directory_if_exists = True

        ad_hoc_run.generate()


    def sanity_check(self):
        self.check_programs()

        if type(self.mcl_inflation) != float:
            raise ConfigError, "Well, MCL likes its inflation parameter in 'float' form..."

        if self.mcl_inflation > 100 or self.mcl_inflation < 0.1:
            raise ConfigError, "MCL inflation parameter should have a reasonable value :/ Like between 0.1 and 100.0."

        if type(self.genomes) != type({}):
            raise ConfigError, "self.genomes must be a dict. Anvi'o needs an adult :("

        if len(self.genomes) < 2:
            raise ConfigError, "There must be at least two genomes for this workflow to work. You have like '%d' of them :/" \
                    % len(self.genomes)

        self.check_params()

        self.run.log_file_path = self.log_file_path
        self.run.info('Args', (str(self.args)), quiet = True)


    def process(self):
        self.sanity_check()

        # initialize genomes
        self.init_genomes()

        # first we will export all proteins 
        unique_proteins_FASTA_path, unique_proteins_names_dict = self.gen_combined_proteins_unique_FASTA()

        # run search
        blastall_results = self.run_search(unique_proteins_FASTA_path, unique_proteins_names_dict)

        # generate MCL input from filtered blastall_results
        mcl_input_file_path = self.gen_mcl_input(blastall_results)

        # get clusters from MCL
        protein_clustering_dict = self.run_mcl(mcl_input_file_path)

        # create view data from protein clusters
        view_data_file_path, view_data_presence_absence_file_path, additional_view_data_file_path, experimental_data_file_path = self.gen_data_from_protein_clusters(protein_clustering_dict)

        # gen samples info and order files
        samples_info_file_path = self.gen_samples_info_file()
        samples_order_file_path = self.gen_samples_order_file(view_data_file_path)

        # gen ad hoc anvi'o run
        self.gen_ad_hoc_anvio_run(view_data_presence_absence_file_path, experimental_data_file_path, additional_view_data_file_path, samples_info_file_path, samples_order_file_path)

        # done
        self.run.info('log file', self.run.log_file_path)
        self.run.quit()


class AdHocRunGenerator:
    """From a matrix file to full-blown anvi'o interface.
    
       This is a class to take in a view data matrix at minimum, and create all
       necessary files for an anvi'o interactive interface call in manual mode."""

    def __init__(self, view_data_path, run = run, progress = progress):
        self.run = run
        self.progress = progress

        self.view_data_path = view_data_path

        self.tree_file_path = None
        self.matrix_data_for_clustering = None
        self.additional_view_data_file_path = None

        self.samples_info_file_path = None
        self.samples_order_file_path = None

        self.output_directory = os.path.abspath('./ad-hoc-anvio-run-directory')
        self.delete_output_directory_if_exists = False

        self.sanity_checked = False


    def sanity_check(self):
        filesnpaths.is_file_tab_delimited(self.view_data_path)
        if self.tree_file_path:
            filesnpaths.is_proper_newick(self.tree_file_path)

        self.check_output_directory()

        new_view_data_path = self.get_output_file_path('view_data.txt')
        shutil.copyfile(self.view_data_path, new_view_data_path)
        self.view_data_path = new_view_data_path

        if self.tree_file_path:
            new_tree_path = self.get_output_file_path('tree.txt')
            shutil.copyfile(self.tree_file_path, new_tree_path)
            self.tree_file_path = new_tree_path

        if self.additional_view_data_file_path:
            new_additional_view_data_file_path = self.get_output_file_path('additional_view_data.txt')
            shutil.copyfile(self.additional_view_data_file_path, new_additional_view_data_file_path)
            self.additional_view_data_file_path = new_additional_view_data_file_path

        self.sanity_checked = True


    def is_good_to_go(self):
        if not self.sanity_checked:
            raise ConfigError, "AdHocRunGenerator :: You gotta be nice, and run sanity check first :/"


    def get_output_file_path(self, file_name):
        return os.path.join(self.output_directory, file_name)


    def check_output_directory(self):
        if os.path.exists(self.output_directory) and not self.delete_output_directory_if_exists:
            raise ConfigError, "AdHocRunGenerator will not work with an existing directory. Please provide a new\
                                path, or use the bool member 'delete_output_directory_if_exists' to overwrite\
                                any existing directory."

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists = self.delete_output_directory_if_exists)


    def generate(self):
        self.sanity_check()

        if not self.tree_file_path:
            self.gen_clustering_of_view_data()

        self.gen_samples_db()


        self.run.info("Ad hoc anvi'o run files", self.output_directory)


    def gen_clustering_of_view_data(self):
        self.is_good_to_go()

        self.progress.new('Hierarchical clustering')
        self.progress.update('..')

        self.tree_file_path = self.get_output_file_path('tree.txt')

        if self.matrix_data_for_clustering:
            clustering.get_newick_tree_data(self.matrix_data_for_clustering, self.tree_file_path)
        else:
            clustering.get_newick_tree_data(self.view_data_path, self.tree_file_path)

        self.progress.end()

        self.run.info('Tree', self.tree_file_path)


    def gen_samples_db(self):
        if not self.samples_info_file_path and not self.samples_order_file_path:
            return

        samples_db_output_path = self.get_output_file_path('samples.db')
        s = dbops.SamplesInformationDatabase(samples_db_output_path, run = self.run, progress = self.progress, quiet = True)
        s.create(self.samples_info_file_path, self.samples_order_file_path)

