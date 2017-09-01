# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import math

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.summarizer as summarizer
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths

from anvio.drivers.blast import BLAST
from anvio.drivers.diamond import Diamond
from anvio.drivers.mcl import MCL
from anvio.drivers import Aligners

from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomestorage import GenomeStorage

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
aligners = Aligners()


class Pangenome(GenomeStorage):
    def __init__(self, args=None, run=run, progress=progress):
        GenomeStorage.__init__(self, args, run, progress)
        self.init_genomes_data_storage()

        self.args = args
        self.run = run
        self.progress = progress

        self.max_num_PCs_for_hierarchical_clustering = constants.max_num_items_for_hierarchical_clustering

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.project_name = A('project_name')
        self.output_dir = A('output_dir')
        self.num_threads = A('num_threads')
        self.skip_alignments = A('skip_alignments')
        self.align_with = A('align_with')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')
        self.min_percent_identity = A('min_percent_identity')
        self.PC_min_occurrence = A('min_occurrence')
        self.mcl_inflation = A('mcl_inflation')
        self.sensitive = A('sensitive')
        self.minbit = A('minbit')
        self.use_ncbi_blast = A('use_ncbi_blast')
        self.exclude_partial_gene_calls = A('exclude_partial_gene_calls')
        self.description_file_path = A('description')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')

        # when it is time to organize PCs
        self.linkage = A('linkage') or constants.linkage_method_default
        self.distance = A('distance') or constants.distance_metric_default

        self.log_file_path = None

        # to be filled during init:
        self.protein_sequences_dict = {}
        self.view_data = {}
        self.view_data_presence_absence = {}
        self.additional_view_data = {}
        self.aligner = None

        # we don't know what we are about
        self.description = None


    def generate_pan_db(self):
        meta_values = {'internal_genome_names': ','.join(self.internal_genome_names),
                       'external_genome_names': ','.join(self.external_genome_names),
                       'num_genomes': len(self.genomes),
                       'min_percent_identity': self.min_percent_identity,
                       'pc_min_occurrence': self.PC_min_occurrence,
                       'mcl_inflation': self.mcl_inflation,
                       'default_view': 'PC_presence_absence',
                       'use_ncbi_blast': self.use_ncbi_blast,
                       'diamond_sensitive': self.sensitive,
                       'minbit': self.minbit,
                       'exclude_partial_gene_calls': self.exclude_partial_gene_calls,
                       'gene_alignments_computed': False if self.skip_alignments else True,
                       'genomes_storage_hash': self.genomes_storage_hash,
                       'project_name': self.project_name,
                       'PCs_clustered': False,
                       'description': self.description if self.description else '_No description is provided_',
                      }

        dbops.PanDatabase(self.pan_db_path, quiet=False).create(meta_values)


    def get_output_file_path(self, file_name, delete_if_exists=False):
        output_file_path = os.path.join(self.output_dir, file_name)

        if delete_if_exists:
            if os.path.exists(output_file_path):
                os.remove(output_file_path)

        return output_file_path


    def check_programs(self):
        if self.use_ncbi_blast:
            utils.is_program_exists('blastp')
        else:
            utils.is_program_exists('diamond')

        utils.is_program_exists('mcl')


    def check_params(self):
        # check the project name:
        if not self.project_name:
            raise ConfigError("Please set a project name, and be prepared to see it around as (1) anvi'o will use\
                                that name to set the output directory and to name various output files such as the\
                                databases that will be generated at the end of the process. If you set your own output\
                                directory name, you can have multiple projects in it and all of those projects can use\
                                the same intermediate files whenever possible.")

        utils.is_this_name_OK_for_database('pan project name', self.project_name, stringent=False)

        # if the user did not set a specific output directory name, use the project name
        # for it:
        self.output_dir = self.output_dir if self.output_dir else self.project_name

        # deal with the output directory:
        try:
            filesnpaths.is_file_exists(self.output_dir)
        except FilesNPathsError:
            filesnpaths.gen_output_directory(self.output_dir, delete_if_exists=self.overwrite_output_destinations)

        filesnpaths.is_output_dir_writable(self.output_dir)
        self.output_dir = os.path.abspath(self.output_dir)

        if not self.log_file_path:
            self.log_file_path = self.get_output_file_path('log.txt')

        filesnpaths.is_output_file_writable(self.log_file_path)
        os.remove(self.log_file_path) if os.path.exists(self.log_file_path) else None

        if not isinstance(self.minbit, float):
            raise ConfigError("minbit value must be of type float :(")

        if self.minbit < 0 or self.minbit > 1:
            raise ConfigError("Well. minbit must be between 0 and 1. Yes. Very boring.")

        if not isinstance(self.min_percent_identity, float):
            raise ConfigError("Minimum percent identity value must be of type float :(")

        if self.min_percent_identity < 0 or self.min_percent_identity > 100:
            raise ConfigError("Minimum percent identity must be between 0%% and 100%%. Although your %.2f%% is\
                               pretty cute, too." % self.min_percent_identity)


        if len([c for c in list(self.genomes.values()) if 'genome_hash' not in c]):
            raise ConfigError("self.genomes does not seem to be a properly formatted dictionary for\
                               the anvi'o class Pangenome.")

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering\
                               while also asking it to enforce it.")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'rU').read()

        if not self.skip_alignments:
            self.aligner = aligners.select(self.align_with)

        self.pan_db_path = self.get_output_file_path(self.project_name + '-PAN.db')


    def run_diamond(self, unique_proteins_fasta_path, unique_proteins_names_dict):
        diamond = Diamond(unique_proteins_fasta_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads, overwrite_output_destinations=self.overwrite_output_destinations)

        diamond.names_dict = unique_proteins_names_dict
        diamond.target_db_path = self.get_output_file_path(filesnpaths.get_name_from_file_path(unique_proteins_fasta_path))
        diamond.search_output_path = self.get_output_file_path('diamond-search-results')
        diamond.tabular_output_path = self.get_output_file_path('diamond-search-results.txt')

        diamond.sensitive = self.sensitive

        return diamond.get_blastall_results()


    def run_blast(self, unique_proteins_fasta_path, unique_proteins_names_dict):
        self.run.warning("You elected to use NCBI's blastp for protein search. Running blastp will be significantly\
                          slower than DIAMOND (although, anvi'o developers are convinced that you *are*\
                          doing the right thing, so, kudos to you).")
        blast = BLAST(unique_proteins_fasta_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads, overwrite_output_destinations=self.overwrite_output_destinations)

        blast.names_dict = unique_proteins_names_dict
        blast.log_file_path = self.log_file_path
        blast.target_db_path = self.get_output_file_path(filesnpaths.get_name_from_file_path(unique_proteins_fasta_path))
        blast.search_output_path = self.get_output_file_path('blast-search-results.txt')

        return blast.get_blastall_results()


    def run_search(self, unique_proteins_fasta_path, unique_proteins_names_dict):
        if self.use_ncbi_blast:
            return self.run_blast(unique_proteins_fasta_path, unique_proteins_names_dict)
        else:
            return self.run_diamond(unique_proteins_fasta_path, unique_proteins_names_dict)


    def run_mcl(self, mcl_input_file_path):
        mcl = MCL(mcl_input_file_path, run=self.run, progress=self.progress, num_threads=self.num_threads)

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
        # will then be used to calculate the 'minbit' value between two genes, which I learned
        # from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8). ITEP defines minbit as
        # 'bit score between target and query / min(selfbit for query, selbit for target)'. This
        # heuristic approach provides a mean to set a cutoff to eliminate weak matches between
        # two genes. minbit value reaches to 1 for hits between two genes that are almost identical.
        self_bit_scores = {}
        line_no = 1
        self.progress.update('(initial pass of the serach results to set the self bit scores ...)')
        for line in open(blastall_results):
            fields = line.strip().split('\t')

            try:
                query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                    [mapping[i](fields[i]) for i in range(0, len(mapping))]
            except Exception as e:
                self.progress.end()
                raise ConfigError("Something went wrong while processing the blastall output file in line %d.\
                                    Here is the error from the uppoer management: '''%s'''" % (line_no, e))
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
                raise ConfigError("Something horrible happened. This can only happend if you started a new analysis with\
                                    additional genomes without cleaning the previous work directory. Sounds familiar?")

            # divide the DNA length of the gene by three to get the AA length, and multiply that by two to get an approximate
            # bit score that would have recovered from a perfect match
            self_bit_scores[id_without_self_search] = (self.genomes[genome_name]['gene_lengths'][int(gene_caller_id)] / 3.0) * 2

            # add this SOB into additional_mcl_input_lines dict.
            additional_mcl_input_lines[id_without_self_search] = '%s\t%s\t1.0\n' % (id_without_self_search, id_without_self_search)


        # CONTINUE AS IF NOTHING HAPPENED
        self.run.info('Min percent identity', self.min_percent_identity)
        self.run.info('Minbit', self.minbit)
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
            # FILTERING BASED ON MINBIT
            #
            minbit = bit_score / min(self_bit_scores[query_id], self_bit_scores[subject_id])
            if minbit < self.minbit:
                continue

            mcl_input.write('%s\t%s\t%f\n' % (query_id, subject_id, perc_id / 100.0))
            num_edges_stored += 1

        # add additional lines if there are any:
        for line in list(additional_mcl_input_lines.values()):
            mcl_input.write(line)
            num_edges_stored += 1

        mcl_input.close()

        self.progress.end()
        self.run.info('Filtered search results', '%s edges stored' % pp(num_edges_stored))
        self.run.info('MCL input', '%s' % mcl_input_file_path)

        return mcl_input_file_path


    def process_PCs(self, PCs_dict):
        self.progress.new('Generating view data')
        self.progress.update('...')

        PCs = list(PCs_dict.keys())

        for PC in PCs:
            self.view_data[PC] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.view_data_presence_absence[PC] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.additional_view_data[PC] = {'num_genes_in_pc': 0, 'num_genomes_pc_has_hits': 0, 'SCG': 0}

            for gene_entry in PCs_dict[PC]:
                genome_name = gene_entry['genome_name']

                self.view_data[PC][genome_name] += 1
                self.view_data_presence_absence[PC][genome_name] = 1
                self.additional_view_data[PC]['num_genes_in_pc'] += 1

            self.additional_view_data[PC]['SCG'] = 1 if set(self.view_data[PC].values()) == set([1]) else 0

            self.additional_view_data[PC]['num_genomes_pc_has_hits'] = len([True for genome in self.view_data[PC] if self.view_data[PC][genome] > 0])

        self.progress.end()

        ########################################################################################
        #                           FILTERING BASED ON OCCURRENCE
        ########################################################################################
        PCs_of_interest = set([])
        for PC in PCs:
            if self.additional_view_data[PC]['num_genomes_pc_has_hits'] >= self.PC_min_occurrence:
                PCs_of_interest.add(PC)

        removed_PCs = 0
        for PC in PCs:
            if PC not in PCs_of_interest:
                self.view_data.pop(PC)
                self.view_data_presence_absence.pop(PC)
                self.additional_view_data.pop(PC)
                PCs_dict.pop(PC)
                removed_PCs += 1

        if self.PC_min_occurrence > 1:
            self.run.info('PCs min occurrence', '%d (the filter removed %d PCs)' % (self.PC_min_occurrence, removed_PCs))

        ########################################################################################
        #            CAN WE CLUSTER THIS STUFF? DOES THE USER WANT US TO TRY REGARDLESS?
        ########################################################################################
        if len(PCs_dict) > self.max_num_PCs_for_hierarchical_clustering:
            if self.enforce_hierarchical_clustering:
                self.run.warning("You have %s PCs, which exceeds the number of PCs anvi'o is comfortable to cluster. But\
                                  since you have used the flag `--enforce-hierarchical-clustering`, anvi'o will attempt\
                                  to create a hierarchical clustering of your PCs anyway. It may take a bit of \
                                  time. Pour yourself a coffee. Or go to a nice vacation. See you in 10 mins, or next year \
                                  or never." % pp(len(PCs_dict)))
            else:
                self.run.warning("It seems you have %s protein clusters in your pangenome. This exceeds the soft limit\
                                  of %s for anvi'o to attempt to create a hierarchical clustering of your protein clusters\
                                  (which becomes the center tree in all anvi'o displays). If you want a hierarchical\
                                  clustering to be done anyway, please see the flag `--enforce-hierarchical-clustering`." \
                                            % (pp(len(PCs_dict)), pp(self.max_num_PCs_for_hierarchical_clustering)))
                self.skip_hierarchical_clustering = True

        ########################################################################################
        #                           STORING FILTERED DATA IN THE DB
        ########################################################################################
        table_structure=['PC'] + sorted(self.genomes.keys())
        table_types=['text'] + ['numeric'] * len(self.genomes)
        dbops.TablesForViews(self.pan_db_path).create_new_view(
                                        data_dict=self.view_data,
                                        table_name='PC_frequencies',
                                        table_structure=table_structure,
                                        table_types=table_types,
                                        view_name = 'PC_frequencies')

        dbops.TablesForViews(self.pan_db_path).create_new_view(
                                        data_dict=self.view_data_presence_absence,
                                        table_name='PC_presence_absence',
                                        table_structure=table_structure,
                                        table_types=table_types,
                                        view_name = 'PC_presence_absence')

        additional_data_structure = ['PC', 'num_genomes_pc_has_hits', 'num_genes_in_pc', 'SCG']
        dbops.TablesForViews(self.pan_db_path).create_new_view(
                                        data_dict=self.additional_view_data,
                                        table_name='additional_data',
                                        table_structure=additional_data_structure,
                                        table_types=['text', 'numeric', 'numeric', 'numeric'],
                                        view_name = None)

        # add additional data structure to the self table, so we can have them initially ordered
        # in the interface the way additional_data_structure suggests:
        pan_db = dbops.PanDatabase(self.pan_db_path, quiet=True)
        pan_db.db.set_meta_value('additional_data_headers', ','.join(additional_data_structure[1:]))
        pan_db.disconnect()

        ########################################################################################
        #                   RETURN THE -LIKELY- UPDATED PROTEIN CLUSTERS DICT
        ########################################################################################
        return PCs_dict


    def cluster_PCs(self):
        """Uses a clustering configuration to add hierarchical clustering of protein clusters into the pan db

        Note how this function cheats the system to create an enchanced clustering configuration:
        We want to use the clustering configurations for pan genomomic analyses to order
        protein clusters. however, we want to add something into the clustering configuraiton
        file, which depends on the number of genomes we have. this addition is 'num_genomes_pc_has_hits'
        data, which pulls together protein clusters that are distributed across genomes similarly based
        on this extra bit of inofrmation. becasue the clustering configurations framework in anvi'o
        does not allow us to have variable information in these recipes, we are going to generate one
        on the fly to have a more capable one."""

        if self.skip_hierarchical_clustering:
            return

        updated_clustering_configs = {}

        for config_name in constants.clustering_configs['pan']:
            config_path = constants.clustering_configs['pan'][config_name]

            # now we have the config path. we first get a temporary file path:
            enhanced_config_path = filesnpaths.get_temp_file_path()

            # setup the additional section based on the number of genomes we have:
            if config_name == 'presence-absence':
                additional_config_section="""\n[AdditionalData !PAN.db::additional_data]\ncolumns_to_use = %s\nnormalize = False\n""" \
                                        % ','.join(['num_genomes_pc_has_hits'] * (int(round(len(self.genomes) / 2))))
            elif config_name == 'frequency':
                additional_config_section="""\n[AdditionalData !PAN.db::additional_data]\ncolumns_to_use = %s\nnormalize = False\nlog=True\n""" \
                                        % ','.join(['num_genes_in_pc'] * (int(round(math.sqrt(len(self.genomes))))))

            # write the content down in to file at the new path:
            open(enhanced_config_path, 'w').write(open(config_path).read() + additional_config_section)

            # update the clustering configs:
            updated_clustering_configs[config_name] = enhanced_config_path

            dbops.do_hierarchical_clusterings(self.pan_db_path, updated_clustering_configs, database_paths={'PAN.db': self.pan_db_path},\
                                              input_directory=self.output_dir, default_clustering_config=constants.pan_default,\
                                              distance=self.distance, linkage=self.linkage, run=self.run, progress=self.progress)


    def gen_samples_db(self):
        samples_info_file_path = self.gen_samples_info_file()
        samples_order_file_path = self.gen_samples_order_file()

        samples_db_output_path = self.get_output_file_path(self.project_name + '-SAMPLES.db', delete_if_exists=True)

        s = dbops.SamplesInformationDatabase(samples_db_output_path, run=self.run, progress=self.progress, quiet=True)
        s.create(samples_info_file_path, samples_order_file_path)


    def gen_samples_order_file(self):
        self.progress.new('Samples DB')
        self.progress.update('Copmputing the hierarchical clustering of the (transposed) view data')

        samples_order_file_path = self.get_output_file_path(self.project_name + '-samples-order.txt')
        samples_order = open(samples_order_file_path, 'w')
        samples_order.write('attributes\tbasic\tnewick\n')

        for clustering_tuple in [('PC presence absence', self.view_data), ('PC frequencies', self.view_data_presence_absence)]:
            v, d = clustering_tuple
            newick = clustering.get_newick_tree_data_for_dict(d, transpose=True, distance = self.distance, linkage=self.linkage)
            samples_order.write('%s\t\t%s\n' % (v, newick))

        samples_order.close()

        self.progress.end()

        self.run.info("Anvi'o samples order", samples_order_file_path)

        return samples_order_file_path


    def gen_samples_info_file(self):
        self.progress.new('Samples DB')
        self.progress.update('Generating the samples information file ..')

        samples_info_dict = {}
        samples_info_file_path = self.get_output_file_path(self.project_name + '-samples-information.txt')

        # set headers
        headers = ['total_length']

        for h in ['percent_completion', 'percent_redundancy']:
            if h in list(self.genomes.values())[0]:
                headers.append(h)

        headers.extend(['gc_content', 'num_genes', 'avg_gene_length', 'num_genes_per_kb'])

        for c in list(self.genomes.values()):
            new_dict = {}
            for header in headers:
                new_dict[header] = c[header]

            samples_info_dict[c['name']] = new_dict

        utils.store_dict_as_TAB_delimited_file(samples_info_dict, samples_info_file_path, headers=['samples'] + headers)

        self.progress.end()
        self.run.info("Anvi'o samples information", samples_info_file_path)

        return samples_info_file_path


    def gen_ad_hoc_anvio_run(self, view_data_file_path, experimental_data_file_path, additional_view_data_file_path, samples_info_file_path):
        ad_hoc_run = summarizer.AdHocRunGenerator(view_data_file_path, run=self.run, progress=self.progress)

        ad_hoc_run.matrix_data_for_clustering = experimental_data_file_path
        ad_hoc_run.additional_view_data_file_path = additional_view_data_file_path
        ad_hoc_run.samples_info_file_path = samples_info_file_path

        ad_hoc_run.output_directory = self.get_output_file_path(os.path.basename(self.output_dir))
        ad_hoc_run.delete_output_directory_if_exists = True

        ad_hoc_run.generate()


    def sanity_check(self):
        self.check_programs()

        if not isinstance(self.mcl_inflation, float):
            raise ConfigError("Well, MCL likes its inflation parameter in 'float' form...")

        if self.mcl_inflation > 100 or self.mcl_inflation < 0.1:
            raise ConfigError("MCL inflation parameter should have a reasonable value :/ Like between 0.1 and 100.0.")

        if not isinstance(self.genomes, type({})):
            raise ConfigError("self.genomes must be a dict. Anvi'o needs an adult :(")

        if len(self.genomes) < 2:
            raise ConfigError("There must be at least two genomes for this workflow to work. You have like '%d' of them :/" \
                    % len(self.genomes))

        if self.skip_alignments and self.align_with:
            raise ConfigError("You are asking anvi'o to skip aligning sequences within your protein clusters, and then you \
                               are also asking it to use '%s' for aligning sequences within your protein clusters. It is easy \
                               to ignore this and skip the alignment, but anvi'o gets nervous when it realizes her users are \
                               being inconsistent. Please make up your mind, and come back as the explicit person you are" \
                                                                            % self.align_with)

        self.check_params()

        self.run.log_file_path = self.log_file_path
        self.run.info('Args', (str(self.args)), quiet=True)


    def store_PCs(self, PCs_dict):
        self.progress.new('Storing protein clusters in the database')
        self.progress.update('...')

        table_for_PCs = dbops.TableForProteinClusters(self.pan_db_path, run=self.run, progress=self.progress)

        num_genes_in_PCs = 0
        for pc_name in PCs_dict:
            for gene_entry in PCs_dict[pc_name]:
                table_for_PCs.add(gene_entry)
                num_genes_in_PCs += 1

        self.progress.end()

        table_for_PCs.store()

        pan_db = dbops.PanDatabase(self.pan_db_path, quiet=True)
        pan_db.db.set_meta_value('num_PCs', len(PCs_dict))
        pan_db.db.set_meta_value('num_genes_in_PCs', num_genes_in_PCs)
        pan_db.disconnect()

        self.run.info('protein clusters info', '%d PCs stored in the database' % len(PCs_dict))


    def gen_PCs_dict_from_mcl_clusters(self, mcl_clusters):
        self.progress.new('Generating the protein clusters dictionary from raw MCL clusters')
        self.progress.update('...')

        PCs_dict = {}

        for PC in mcl_clusters:
            PCs_dict[PC] = []

            for entry_hash, gene_caller_id in [e.split('_') for e in mcl_clusters[PC]]:
                try:
                    genome_name = self.hash_to_genome_name[entry_hash]
                except KeyError:
                    self.progress.end()
                    raise ConfigError("Something horrible happened. This can only happen if you started a new analysis with\
                                        additional genomes without cleaning the previous work directory. Sounds familiar?")

                PCs_dict[PC].append({'gene_caller_id': int(gene_caller_id), 'protein_cluster_id': PC, 'genome_name': genome_name, 'alignment_summary': ''})

        self.progress.end()

        return PCs_dict


    def compute_alignments_for_PCs(self, PCs_dict):
        if self.skip_alignments:
            self.run.warning('Skipping gene alignments.')
            return PCs_dict

        r = terminal.Run()
        r.verbose = False

        self.progress.new('Aligning genes in protein sequences')
        self.progress.update('...')
        pc_names = list(PCs_dict.keys())
        num_pcs = len(pc_names)
        for i in range(0, num_pcs):
            self.progress.update('%d of %d' % (i, num_pcs)) if i % 10 == 0 else None
            pc_name = pc_names[i]

            if len(PCs_dict[pc_name]) == 1:
                # this sequence is a singleton and does not need alignment
                continue

            gene_sequences_in_pc = []
            for gene_entry in PCs_dict[pc_name]:
                sequence = self.genomes_storage.get_gene_sequence(gene_entry['genome_name'], gene_entry['gene_caller_id'])
                gene_sequences_in_pc.append(('%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id']), sequence),)

            # alignment
            alignments = self.aligner(run=r).run_stdin(gene_sequences_in_pc)

            for gene_entry in PCs_dict[pc_name]:
                gene_entry['alignment_summary'] = utils.summarize_alignment(alignments['%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id'])])

        self.progress.end()

        return PCs_dict


    def process(self):
        # check sanity
        self.sanity_check()

        # gen pan_db
        self.generate_pan_db()

        # get all protein sequences:
        combined_aas_FASTA_path = self.get_output_file_path('combined-aas.fa')
        unique_aas_FASTA_path, unique_aas_names_dict = self.genomes_storage.gen_combined_aa_sequences_FASTA(combined_aas_FASTA_path, exclude_partial_gene_calls=self.exclude_partial_gene_calls)

        # run search
        blastall_results = self.run_search(unique_aas_FASTA_path, unique_aas_names_dict)

        # generate MCL input from filtered blastall_results
        mcl_input_file_path = self.gen_mcl_input(blastall_results)

        # get clusters from MCL
        mcl_clusters = self.run_mcl(mcl_input_file_path)

        # we have the raw protein clusters dict, but we need to re-format it for following steps
        PCs_dict = self.gen_PCs_dict_from_mcl_clusters(mcl_clusters)
        del mcl_clusters

        # compute alignments for genes within each PC (or don't)
        PCs_dict = self.compute_alignments_for_PCs(PCs_dict)

        # populate the pan db with results
        PCs_dict = self.process_PCs(PCs_dict)

        # store protein clusters dict into the db
        self.store_PCs(PCs_dict)

        # generate a hierarchical clustering of protein clusters (or don't)
        self.cluster_PCs()

        # gen samples info and order files
        self.gen_samples_db()

        # done
        self.run.info('log file', self.run.log_file_path)
        self.run.quit()
