# -*- coding: utf-8
"""
    Classes for pan operations.

    anvi-pan-genome is the default client of this module
"""

import os
import shutil
import tempfile

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.clustering as clustering
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

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
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        input_file_for_contig_dbs = A('input_contig_dbs')
        self.num_threads = A('num_threads')
        self.output_dir = A('output_dir')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')
        self.min_percent_identity = A('min_percent_identity')

        self.temp_files_to_remove_later = []

        self.contig_dbs = utils.get_TAB_delimited_file_as_dictionary(input_file_for_contig_dbs, expected_fields = ['name', 'path']) if input_file_for_contig_dbs else {}

        # convert relative paths to absolute paths
        for contigs_db in self.contig_dbs:
            path = self.contig_dbs[contigs_db]['path']
            if not path.startswith('/'):
                self.contig_dbs[contigs_db]['path'] = os.path.abspath(os.path.join(os.path.dirname(input_file_for_contig_dbs), path))

        # to be filled during init:
        self.hash_to_contigs_db_name = {}
        self.view_data = {}


    def get_output_file_path(self, file_name, temp_file = False):
        output_file_path = os.path.join(self.output_dir, file_name)

        if temp_file:
            self.temp_files_to_remove_later.append(output_file_path)

        return output_file_path


    def remove_temp_files(self):
        if self.debug:
            return

        for file_path in self.temp_files_to_remove_later:
            os.remove(file_path)


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

        if type(self.min_percent_identity) != float:
            raise ConfigError, "Minimum percent identity value must be of type float :("

        if self.min_percent_identity < 20 or self.min_percent_identity > 100:
            raise ConfigError, "Minimum percent identity must be between 20%% and 100%%. Although your %.2f%% is\
                                pretty cute, too." % self.min_percent_identity


    def init_contig_dbs(self):
        if type(self.contig_dbs) != type({}):
            raise ConfigError, "self.contig_dbs must be of type dict. Anvi'o needs an adult :("

        if not len(self.contig_dbs):
            raise ConfigError, "There is no contig databases to process..."

        if len(self.contig_dbs) < 2:
            raise ConfigError, "There must be at least two contigs databases for this to work :/"

        if len([c for c in self.contig_dbs.values() if 'path' not in c]):
            raise ConfigError, "self.contig_dbs does not seem to be a properly formatted dictionary for\
                                the anvi'o class Pangenome. You did something very wrong."

        missing_dbs = [c['path'] for c in self.contig_dbs.values() if not os.path.exists(c['path'])]
        if len(missing_dbs):
            raise ConfigError, "%d of %d of your contigs databases are not found where they were supposed to be \
                                based on the description you provided :( Here is one that is missing: '%s'" \
                                                % (len(missing_dbs), len(self.contig_dbs), missing_dbs[0])

        # just go over the contig dbs to make sure they all are OK, AAAAAND set some stuff for later use.
        self.progress.new('Initializing')
        contig_db_names = self.contig_dbs.keys()
        for i in range(0, len(contig_db_names)):
            contigs_db_name = contig_db_names[i]
            c = self.contig_dbs[contigs_db_name]
            c['name'] = contigs_db_name

            self.progress.update('%d of %d ... %s' % (i + 1, len(self.contig_dbs), contigs_db_name))

            contigs_db_summary = summarizer.get_contigs_db_info_dict(c['path'])

            for key in contigs_db_summary:
                c[key] = contigs_db_summary[key]

            self.hash_to_contigs_db_name[c['contigs_db_hash']] = contigs_db_name
        self.progress.end()

        # if two contigs db has the same hash, we are kinda f'd:
        if len(set([c['contigs_db_hash'] for c in self.contig_dbs.values()])) != len(self.contig_dbs):
            raise ConfigError, 'Not all hash values are unique across all contig databases you provided. Something\
                                very fishy is going on :/'

        # make sure genes are called in every contigs db:
        if len([c['genes_are_called'] for c in self.contig_dbs.values()]) != len(self.contig_dbs):
            raise ConfigError, 'Genes are not called in every contigs db in the collection :/'

        self.run.info('Contig DBs', '%d contig databases have been found.' % len(self.contig_dbs))


    def gen_combined_proteins_fasta(self):
        self.progress.new('Storing combined protein sequences')
        output_file_path = self.get_output_file_path('combined_proteins.fa', temp_file = True)
        output_file = open(output_file_path, 'w')

        for c in self.contig_dbs.values():
            self.progress.update('Working on %s ...' % c['name'])
            num_genes = 0
            contigs_db = dbops.ContigsDatabase(c['path'])
            protein_sequences = contigs_db.db.get_table_as_dict(t.gene_protein_sequences_table_name)
            for protein_id in protein_sequences:
                num_genes += 1
                output_file.write('>%s_%d\n' % (c['contigs_db_hash'], protein_id))
                output_file.write('%s\n' % protein_sequences[protein_id]['sequence'])
            contigs_db.disconnect()
            c['num_genes'] = num_genes

        output_file.close()
        self.progress.end()

        self.run.info('ORFs', '%s protein sequences are stored for analysis.' % pp(sum([c['num_genes'] for c in self.contig_dbs.values()])))

        return output_file_path


    def run_diamond(self, combined_proteins_fasta_path):
        diamond = Diamond(combined_proteins_fasta_path, run = self.run, progress = self.progress,
                          num_threads = self.num_threads, overwrite_output_destinations = self.overwrite_output_destinations)

        diamond.log_file_path = self.get_output_file_path('log.txt')
        diamond.target_db_path = self.get_output_file_path('.'.join(combined_proteins_fasta_path.split('.')[:-1]))
        diamond.search_output_path = self.get_output_file_path('diamond-search-results')
        diamond.tabular_output_path = self.get_output_file_path('diamond-search-results.txt')

        return diamond.get_blastall_results()


    def run_mcl(self, mcl_input_file_path):
        mcl = MCL(mcl_input_file_path, run = self.run, progress = self.progress, num_threads = self.num_threads)

        mcl.clusters_file_path = self.get_output_file_path('mcl-clusters.txt')
        mcl.log_file_path = self.get_output_file_path('log.txt')

        return mcl.get_clusters_dict()


    def gen_mcl_input(self, blastall_results):
        self.progress.new('Filtering blastall results')
        self.progress.update('...')

        mcl_input_file_path = self.get_output_file_path('mcl-input.txt')
        mcl_input = open(mcl_input_file_path, 'w')

        mapping = [str, str, float, int, int, int, int, int, int, int, float, float]

        line_no = 1
        num_edges_stored = 0
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

            if line_no % 5000 == 0:
                self.progress.update('Lines processed %s ...' % pp(line_no))

            #
            # FILTERS
            #

            if perc_id < self.min_percent_identity:
                continue

            mcl_input.write('%s\t%s\t%f\n' % (query_id, subject_id, perc_id / 100.0))
            num_edges_stored += 1


        mcl_input.close()

        self.progress.end()
        self.run.info('Filtered diamond results', '%s edges stored' % pp(num_edges_stored))
        self.run.info('MCL input', '%s' % mcl_input_file_path)

        return mcl_input_file_path


    def gen_view_data_from_protein_clusters(self, protein_clustering_dict):
        self.progress.new('Generating view data')
        self.progress.update('...')

        for pc in protein_clustering_dict:
            self.view_data[pc] = dict([(contigs_db_name, 0) for contigs_db_name in self.contig_dbs])
            for contigs_db_hash, gene_callers_id in [e.split('_') for e in protein_clustering_dict[pc]]:
                contigs_db_name = self.hash_to_contigs_db_name[contigs_db_hash]
                self.view_data[pc][contigs_db_name] += 1

        view_data_file_path = self.get_output_file_path('anvio-view-data.txt')
        utils.store_dict_as_TAB_delimited_file(self.view_data, view_data_file_path, headers = ['contig'] + sorted(self.contig_dbs.keys()))

        self.progress.end()
        self.run.info("Anvi'o view data for protein clusters", view_data_file_path)

        return view_data_file_path


    def gen_samples_info_file(self):
        samples_info_dict = {}
        samples_info_file_path = self.get_output_file_path('anvio-samples-information.txt')

        # set headers
        headers = ['total_length']

        for h in ['percent_complete', 'percent_redundancy']:
            if self.contig_dbs.values()[0].has_key(h):
                headers.append(h)

        headers.extend(['gc_content', 'num_genes', 'avg_gene_length', 'num_genes_per_kb'])

        for c in self.contig_dbs.values():
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


    def gen_ad_hoc_anvio_run(self, view_data_file_path, samples_info_file_path, samples_order_file_path):
        ad_hoc_run = AdHocRunGenerator(view_data_file_path, run = self.run, progress = self.progress)

        ad_hoc_run.tree_file_path = '/Users/meren/github/anvio/tests/sandbox/anvi_pangenome_files/test-output/tree.txt'
        ad_hoc_run.samples_info_file_path = samples_info_file_path
        ad_hoc_run.samples_order_file_path = samples_order_file_path

        ad_hoc_run.output_directory = self.get_output_file_path('anvio-run-files')
        ad_hoc_run.delete_output_directory_if_exists = True

        ad_hoc_run.generate()


    def sanity_check(self):
        self.check_programs()
        self.check_params()
        self.init_contig_dbs()


    def process(self):
        self.sanity_check()

        # first we will export all proteins 
        combined_proteins_fasta_path = self.gen_combined_proteins_fasta()

        # run diamond
        blastall_results = self.run_diamond(combined_proteins_fasta_path)

        # generate MCL input from filtered blastall_results
        mcl_input_file_path = self.gen_mcl_input(blastall_results)

        # get clusters from MCL
        protein_clustering_dict = self.run_mcl(mcl_input_file_path)

        # create view data from protein clusters
        view_data_file_path = self.gen_view_data_from_protein_clusters(protein_clustering_dict)

        # gen samples info and order files
        samples_info_file_path = self.gen_samples_info_file()
        samples_order_file_path = self.gen_samples_order_file(view_data_file_path)

        # gen ad hoc anvi'o run
        self.gen_ad_hoc_anvio_run(view_data_file_path, samples_info_file_path, samples_order_file_path)


class AdHocRunGenerator:
    """From a matrix file to full-blown anvi'o interface.
    
       This is a class to take in a view data matrix at minimum, and create all
       necessary files for an anvi'o interactive interface call in manual mode."""

    def __init__(self, view_data_path, tree_file_path = None, samples_info_file_path = None, samples_order_file_path = None, run = run, progress = progress):
        self.run = run
        self.progress = progress

        self.view_data_path = view_data_path
        self.tree_file_path = tree_file_path

        self.samples_info_file_path = samples_info_file_path
        self.samples_order_file_path = samples_order_file_path


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

        self.progress.new('Hierarchical clustering of the view data')
        self.progress.update('..')

        self.tree_file_path = self.get_output_file_path('tree.txt')
        clustering.get_newick_tree_data(self.view_data_path, self.tree_file_path)

        self.progress.end()

        self.run.info('Tree', self.tree_file_path)


    def gen_samples_db(self):
        if not self.samples_info_file_path and not self.samples_order_file_path:
            return

        samples_db_output_path = self.get_output_file_path('samples.db')
        s = dbops.SamplesInformationDatabase(samples_db_output_path, run = self.run, progress = self.progress, quiet = True)
        s.create(self.samples_info_file_path, self.samples_order_file_path)


class MCL:
    def __init__(self, mcl_input_file_path, run = run, progress = progress, num_threads = 1):
        self.run = run
        self.progress = progress

        self.mcl_input_file_path = mcl_input_file_path
        self.num_threads = num_threads

        utils.is_program_exists('mcl')

        self.inflation = 2.0

        self.clusters_file_path = 'mcl-clusters.txt'
        self.log_file_path = 'mcl-log-file.txt'


    def check_output(self, expected_output, process = 'diamond'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError, "Pfft. Something probably went wrong with MCL's '%s' since one of the expected output files are missing.\
                                Please check the log file here: '%s." % (process, self.log_file_path)


    def get_clusters_dict(self):
        self.cluster()

        clusters_dict = {}

        line_no = 1
        for line in open(self.clusters_file_path).readlines():
            clusters_dict['PC_%08d' % line_no] = line.strip().split('\t')

            line_no += 1

        self.run.info('Clusters', '%s clusters processed' % pp(len(clusters_dict)))

        return clusters_dict


    def cluster(self):
        self.progress.new('MCL')
        self.progress.update('clustering (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('mcl %s --abc -I %f -o %s -te %d >> "%s" 2>&1' % (self.mcl_input_file_path,
                                                                         self.inflation,
                                                                         self.clusters_file_path,
                                                                         self.num_threads,
                                                                         self.log_file_path))

        with open(self.log_file_path, "a") as log: log.write('MCL CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()

        self.check_output(self.clusters_file_path, 'makedb')

        self.run.info('MCL output', self.clusters_file_path)


class Diamond:
    def __init__(self, query_fasta, run = run, progress = progress, num_threads = 1, overwrite_output_destinations = False):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations

        utils.is_program_exists('diamond')

        self.tmp_dir = tempfile.gettempdir()

        self.query_fasta = query_fasta
        self.log_file_path = 'diamond-log-file.txt'
        self.target_db_path = 'diamond-target'
        self.search_output_path = 'diamond-search-resuults'
        self.tabular_output_path = 'diamond-search-results.txt'


    def get_blastall_results(self):
        force_makedb, force_blastp, force_view = False, False, False

        if self.overwrite_output_destinations:
            force_makedb = True

        if os.path.exists(self.target_db_path + '.dmnd') and not force_makedb:
            self.run.warning("Notice: A diamond database is found in the output directory, and will be used!")
        else:
            self.makedb()
            force_blastp, forrce_view = True, True

        if os.path.exists(self.search_output_path + '.daa') and not force_blastp:
            self.run.warning("Notice: A DIAMOND search result is found in the output directory: skipping BLASTP!")
        else:
            self.blastp()
            force_view = True

        if os.path.exists(self.tabular_output_path) and not force_view:
            self.run.warning("Notice: A DIAMOND tabular output is found in the output directory. Anvi'o will not generate another one!")
        else:
            self.view()

        return self.tabular_output_path


    def check_output(self, expected_output, process = 'diamond'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError, "Pfft. Something probably went wrong with Diamond's '%s' since one of the expected output files are missing.\
                                Please check the log file here: '%s." % (process, self.log_file_path)


    def makedb(self):
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('diamond makedb --in %s -d %s -p %d >> "%s" 2>&1' % (self.query_fasta,
                                                                         self.target_db_path,
                                                                         self.num_threads,
                                                                         self.log_file_path))

        with open(self.log_file_path, "a") as log: log.write('CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()

        expected_output = self.target_db_path + '.dmnd'
        self.check_output(expected_output, 'makedb')

        self.run.info('Diamond temp search db', expected_output)


    def blastp(self):
        self.progress.new('DIAMOND')
        self.progress.update('running blastp (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('diamond blastp -q %s -d %s -a %s -t %s -p %d >> "%s" 2>&1' % (self.query_fasta,
                                                                                   self.target_db_path,
                                                                                   self.search_output_path,
                                                                                   self.tmp_dir,
                                                                                   self.num_threads,
                                                                                   self.log_file_path))
        with open(self.log_file_path, "a") as log: log.write('CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()

        expected_output = self.search_output_path + '.daa'
        self.check_output(expected_output, 'blastp')

        self.run.info('Diamond blastp results', expected_output)


    def view(self):
        self.progress.new('DIAMOND')
        self.progress.update('generating tabular output (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('diamond view -a %s -o %s -p %d >> "%s" 2>&1' % (self.search_output_path + '.daa',
                                                                     self.tabular_output_path,
                                                                     self.num_threads,
                                                                     self.log_file_path))
        with open(self.log_file_path, "a") as log: log.write('CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()

        self.check_output(self.tabular_output_path, 'view')

        self.run.info('Diamond tabular output file', self.tabular_output_path)
