# coding: utf-8
# pylint: disable=line-too-long
"""Provides the necessary class to profile BAM files."""

import gc
import os
import sys
import copy
import shutil
import argparse
import numpy as np

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from collections import OrderedDict

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.bamops as bamops
import anvio.terminal as terminal
import anvio.contigops as contigops
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.tables.views import TablesForViews
from anvio.tables.indels import TableForIndels
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.tables.variability import TableForVariability
from anvio.tables.codonfrequencies import TableForCodonFrequencies


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"

null_progress = terminal.Progress(verbose=False)
null_run = terminal.Run(verbose=False)
pp = terminal.pretty_print


class BAMProfilerQuick:
    """A class for rapid profiling of BAM files that produce text files rather than profile dbs"""

    def __init__(self, args, skip_sanity_check=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bam_file_paths = A('bam_files')
        self.contigs_db_path = A('contigs_db')
        self.output_file_path = A('output_file')
        self.gene_level_stats = A('gene_mode')
        self.gene_caller = A('gene_caller')
        self.report_minimal = A('report_minimal')

        if not skip_sanity_check:
            self.sanity_check()

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Num BAM files', len(self.bam_file_paths))
        self.run.info('Reporting', 'MINIMAL' if self.report_minimal else 'EVERYTHING', mc="red" if self.report_minimal else "green")


        # to be filled later if necessary
        self.contigs_basic_info = {}
        self.gene_calls_per_contig = {}


    def sanity_check(self):
        if not self.contigs_db_path:
            raise ConfigError("You need to provide an anvi'o contigs database for this to work :/")

        utils.is_contigs_db(self.contigs_db_path)

        if not len(self.bam_file_paths):
            raise ConfigError("You need to provide at least one BAM file for this to work.")

        if not self.output_file_path:
            raise ConfigError("Please provide an output file path.")

        filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=False)

        # find all the bad BAM files
        self.progress.new("Sanity checking BAM files")
        self.progress.update('...')
        bad_bam_files = [f for f in self.bam_file_paths if not filesnpaths.is_file_bam_file(f, dont_raise=True)]
        if len(bad_bam_files):
            self.progress.reset()
            raise ConfigError(f"Not all of your BAM files look like BAM files. Here is the list that "
                              f"samtools didn't like: {', '.join(bad_bam_files)}")
        self.progress.end()

        if self.gene_level_stats:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
            genes_are_called = contigs_db.meta['genes_are_called']
            gene_callers_list = contigs_db.meta['gene_callers']
            contigs_db.disconnect()

            if not genes_are_called:
                raise ConfigError("There are no gene calls in this contigs database :/ You can't use the flag "
                                  "`--report-gene-level-stats`. Yes.")

            gene_callers = [tpl[0] for tpl in gene_callers_list]
            if self.gene_caller not in gene_callers:
                dbops.ContigsDatabase(self.contigs_db_path).list_gene_caller_sources()
                raise ConfigError(f"The gene caller '{self.gene_caller}' is not among those that are found in "
                                  f"the contigs database (and shown above for your convenience).")



    def process(self):
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)

        self.run.warning(None, header="CONTIGS DB", lc="green")
        self.run.info("Number of contigs", pp(contigs_db.meta['num_contigs']))
        self.run.info("Total num nucleotides", pp(contigs_db.meta['total_length']))
        self.run.info("Genes are called", "Yes" if contigs_db.meta['genes_are_called'] else "No :/")

        self.run.info("Reporting mode", "GENES" if self.gene_level_stats else "CONTIGS", mc='green')
        if self.gene_level_stats:
            self.run.info("Gene caller", self.gene_caller)
            self.run.info("Number of genes", pp([tpl[1] for tpl in contigs_db.meta['gene_callers'] if tpl[0] == self.gene_caller][0]))

        self.progress.new('Reading data into memory')
        self.progress.update('Contigs basic info table ...')
        self.contigs_basic_info = contigs_db.db.get_table_as_dict(t.contigs_info_table_name)
        self.progress.end()

        contigs_db.disconnect()

        if self.gene_level_stats:
            self.recover_gene_data()

        self.report_stats()


    def recover_gene_data(self):
        # if we are working with genes, we need to first read the gene calls table into
        # memory
        self.progress.new('Reading data into memory')
        self.progress.update('Gene calls table ...')
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        self.genes_in_contigs = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=f"source == '{self.gene_caller}'", error_if_no_data=True)

        # and then update contigs basic info to easily track gene calls in a given contig
        # for reporting purposes
        self.progress.update('Updating contigs basic info ...')
        for gene_callers_id in self.genes_in_contigs:
            contig_name = self.genes_in_contigs[gene_callers_id]['contig']

            if 'gene_caller_ids' not in self.contigs_basic_info[contig_name]:
                self.contigs_basic_info[contig_name]['gene_caller_ids'] = set([gene_callers_id])
            else:
                self.contigs_basic_info[contig_name]['gene_caller_ids'].add(gene_callers_id)

        self.progress.end()

        contigs_db.disconnect()


    def report_stats(self):
        """Iterates through bam files, reports contigs stats"""

        # the total number of items we will process is equal to the number of contigs that will
        # have to be processed for each BAM file
        total_num_items = len(self.bam_file_paths) * len(self.contigs_basic_info)
        total_num_bam_files = len(self.bam_file_paths)
        contig_names = list(self.contigs_basic_info.keys())
        num_contigs = len(self.contigs_basic_info)

        self.progress.new("Bleep bloop", progress_total_items=total_num_items)
        self.progress.update('...')

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        with open(self.output_file_path, 'w') as output:
            if self.report_minimal and not self.gene_level_stats:
                header = ['contig', 'sample', 'length', 'gc_content', 'num_mapped_reads', 'detection', 'mean_cov']
            elif not self.report_minimal and not self.gene_level_stats:
                header = ['contig', 'sample', 'length', 'gc_content', 'num_mapped_reads', 'detection', 'mean_cov', 'q2q3_cov', 'median_cov', 'min_cov', 'max_cov', 'std_cov']
            elif self.report_minimal and self.gene_level_stats:
                header = ['gene_callers_id', 'contig', 'sample', 'length', 'detection', 'mean_cov']
            else:
                header = ['gene_callers_id', 'contig', 'sample', 'length', 'num_mapped_reads', 'detection', 'mean_cov', 'q2q3_cov', 'median_cov', 'min_cov', 'max_cov', 'std_cov']

            output.write('\t'.join(header) + '\n')

            for i in range(0, total_num_bam_files):
                bam_file_path = self.bam_file_paths[i]

                bam = bamops.BAMFileObject(bam_file_path, 'rb')
                bam_file_name = os.path.splitext(os.path.basename(bam_file_path))[0]

                for j in range(0, num_contigs):
                    contig_name = contig_names[j]

                    if j - 1 == 0 or j % 100 == 0:

                        if mem_tracker.measure():
                            mem_usage = mem_tracker.get_last()
                            mem_diff = mem_tracker.get_last_diff()

                        self.progress.increment(increment_to=((i * num_contigs) + j))
                        self.progress.update(f"BAM {i+1}/{pp(total_num_bam_files)} :: CONTIG {pp(j)}/{pp(num_contigs)} :: MEMORY 🧠 {mem_usage} ({mem_diff})")

                    c = bamops.Coverage()

                    try:
                        c.run(bam, contig_name, read_iterator='fetch', skip_coverage_stats=True)
                    except:
                        continue

                    if self.report_minimal and not self.gene_level_stats:
                        # we are in contigs mode, and want a minimal report
                        mean = np.mean(c.c)
                        detection = np.sum(c.c > 0) / len(c.c)
                        output.write(f"{contig_name}\t"
                                     f"{bam_file_name}\t"
                                     f"{self.contigs_basic_info[contig_name]['length']}\t"
                                     f"{self.contigs_basic_info[contig_name]['gc_content']:.3}\t"
                                     f"{c.num_reads}\t"
                                     f"{detection:.4}\t"
                                     f"{mean:.4}\n")
                    elif not self.report_minimal and not self.gene_level_stats:
                        # we are in contigs mode, but want an extended report
                        C = utils.CoverageStats(c.c, skip_outliers=True)
                        output.write(f"{contig_name}\t"
                                     f"{bam_file_name}\t"
                                     f"{self.contigs_basic_info[contig_name]['length']}\t"
                                     f"{self.contigs_basic_info[contig_name]['gc_content']:.3}\t"
                                     f"{c.num_reads}\t"
                                     f"{C.detection:.4}\t"
                                     f"{C.mean:.4}\t"
                                     f"{C.mean_Q2Q3:.4}\t"
                                     f"{C.median}\t"
                                     f"{C.min}\t"
                                     f"{C.max}\t"
                                     f"{C.std:.4}\n")
                    elif self.gene_level_stats:
                        # we are in gene mode!
                        if 'gene_caller_ids' not in self.contigs_basic_info[contig_name]:
                            # this means we don't have a gene call in this contig. Long hair don't
                            # care. MOVING ON.
                            continue

                        for gene_callers_id in self.contigs_basic_info[contig_name]['gene_caller_ids']:
                            g = self.genes_in_contigs[gene_callers_id]
                            gc = c.c[g['start']:g['stop']]

                            if self.report_minimal:
                                # we want a minimal report
                                mean = np.mean(gc)
                                detection = np.sum(gc > 0) / len(gc)
                                output.write(f"{gene_callers_id}\t"
                                             f"{contig_name}\t"
                                             f"{bam_file_name}\t"
                                             f"{g['stop'] - g['start']}\t"
                                             f"{detection:.4}\t"
                                             f"{mean:.4}\n")
                            else:
                                # calculate the total number of reads mapping to the gene:
                                num_mapped_reads_to_gene = 0
                                for r in bam.fetch_only(contig_name, start=g['start'], end=g['stop']):
                                    num_mapped_reads_to_gene += 1

                                # get coverage stats:
                                GC = utils.CoverageStats(gc, skip_outliers=True)

                                # write everything out:
                                output.write(f"{gene_callers_id}\t"
                                             f"{contig_name}\t"
                                             f"{bam_file_name}\t"
                                             f"{g['stop'] - g['start']}\t"
                                             f"{num_mapped_reads_to_gene}\t"
                                             f"{GC.detection:.4}\t"
                                             f"{GC.mean:.4}\t"
                                             f"{GC.mean_Q2Q3:.4}\t"
                                             f"{GC.median}\t"
                                             f"{GC.min}\t"
                                             f"{GC.max}\t"
                                             f"{GC.std:.4}\n")
                    else:
                        raise ConfigError("We need an adult :(")

                bam.close()

            self.progress.end()


class BAMProfiler(dbops.ContigsSuperclass):
    """Creates an über class for BAM file operations"""

    def __init__(self, args, r=terminal.Run(width=50), p=terminal.Progress()):
        self.args = args
        self.progress = p
        self.run = r

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_file_path = A('input_file')
        self.contigs_db_path = A('contigs_db')
        self.serialized_profile_path = A('serialized_profile')
        self.output_directory = A('output_dir')
        self.list_contigs_and_exit = A('list_contigs')
        self.min_contig_length = A('min_contig_length') or 0
        self.max_contig_length = A('max_contig_length') or sys.maxsize
        self.min_coverage_for_variability = A('min_coverage_for_variability')
        self.contigs_shall_be_clustered = A('cluster_contigs')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.sample_id = A('sample_name')
        self.report_variability_full = A('report_variability_full')
        self.skip_edges = A('skip_edges')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.skip_SNV_profiling = A('skip_SNV_profiling')
        self.skip_INDEL_profiling = A('skip_INDEL_profiling')
        self.profile_SCVs = A('profile_SCVs')
        self.min_percent_identity = A('min_percent_identity')
        self.fetch_filter = A('fetch_filter')
        self.gen_serialized_profile = A('gen_serialized_profile')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.num_threads = int(A('num_threads') or 1)
        self.queue_size = int(A('queue_size') if A('queue_size') is not None else 0)
        self.write_buffer_size_per_thread = int(A('write_buffer_size_per_thread') if A('write_buffer_size_per_thread') is not None else 500)
        self.write_buffer_size = self.write_buffer_size_per_thread * self.num_threads
        self.total_length_of_all_contigs = 0
        self.total_coverage_values_for_all_contigs = 0
        self.total_reads_kept = 0
        self.description_file_path = A('description')

        if self.fetch_filter:
            valid_fetch_filters = [k for k in constants.fetch_filters.keys() if k]

            if self.fetch_filter not in valid_fetch_filters:
                raise ConfigError(f"Your fetch filter '{self.fetch_filter}' is not among those anvi'o knows about :/ If you "
                                  f"would like to try again with a different one, here is a list for your consideration: "
                                  f"{', '.join(valid_fetch_filters)}.")

        # these are the views this class will be filling in:
        self.essential_data_fields_for_anvio_profiles = copy.deepcopy(constants.essential_data_fields_for_anvio_profiles)

        # but the views should not include the view `variability` if SNVs are not going to be profiled:
        if 'variability' in self.essential_data_fields_for_anvio_profiles and self.skip_SNV_profiling:
            self.essential_data_fields_for_anvio_profiles.pop(self.essential_data_fields_for_anvio_profiles.index('variability'))

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        # whether the profile database is a blank (without any BAM files or reads):
        self.blank = A('blank_profile')

        if not self.blank and self.contigs_shall_be_clustered and self.skip_hierarchical_clustering:
            raise ConfigError("You are confused, and confusing anvi'o, too. You can't as hierarchical clustering "
                              "to be performed with one flag, and try to skip it with another one :(")

        if self.blank and self.contigs_shall_be_clustered and self.skip_hierarchical_clustering:
            raise ConfigError("So you want to generate a blank profile, and you both want hierarchical clustering "
                              "of your contigs to be performed, and skipped. No.")

        if self.blank and self.contigs_shall_be_clustered:
            raise ConfigError("When the blank profile is asked to be generated, there is no need to ask for the "
                              "hierarchical clustering of contigs. It is going to be done by default. If it is "
                              "not changing anything, why is anvi'o upset with you? Because. Let's don't use flags "
                              "we don't need.")

        if self.blank and not self.skip_hierarchical_clustering:
            self.contigs_shall_be_clustered = True

        if A('contigs_of_interest'):
            filesnpaths.is_file_exists(args.contigs_of_interest)
            self.contig_names_of_interest = set([c.strip() for c in open(args.contigs_of_interest).readlines()\
                                                                           if c.strip() and not c.startswith('#')])
        else:
            self.contig_names_of_interest = set([])

        if self.list_contigs_and_exit:
            self.list_contigs()
            sys.exit()

        if not self.contigs_db_path:
            raise ConfigError("No contigs database, no profilin'. Bye.")

        # store our run object. 'why?', you may ask. well, keep reading.
        my_run = self.run

        # Initialize contigs db
        dbops.ContigsSuperclass.__init__(self, self.args, r=null_run, p=self.progress)
        self.init_contig_sequences(contig_names_of_interest=self.contig_names_of_interest)
        self.contig_names_in_contigs_db = set(self.contigs_basic_info.keys())

        # restore our run object. OK. take deep breath. because you are a good programmer, you have
        # this voice in your head tellin gyou that the tinkering with self.run here feels kind of out
        # of place. that voice is right, but the voice doesn't know the struggles of poor souls that
        # had to resort to a solution like this. You see, we don't want to see any run messages from
        # ContigsSuper in profiler output. But when we pass a `run=null_run` to the the class, due to
        # inheritance, it also modifies our own `self.run` with the null one, making the profilesuper
        # go all quiet. so here we basically need to re-engage our `self.run`. But then if the user
        # actually ASKED for ProfileSuper to be quiet, then we can't simply just inherit another Run
        # object and move on with our lives in a single line right here. We in fact need to store our
        # previous run object and then re-engage that one instead of a brand new one. There you have
        # it (yes. funny detail, no one cares, but I do because I have no friends -- IF YOU ARE
        # READING THESE LINES YOU ARE AUTOMATICALLY MY FRIEND, so SORRY):
        self.run = my_run

        self.bam = None
        self.contigs = []

        self.database_paths = {'CONTIGS.db': os.path.abspath(self.contigs_db_path)}

        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['blank' if self.blank else 'single']

        # if genes are not called, yet the user is asking for codon frequencies to be profiled, we give
        # a warning and force-turn that flag off.
        if (not self.a_meta['genes_are_called']) and self.profile_SCVs:
            self.run.warning("You asked the codon frequencies to be profiled, but genes were not called "
                             "for your contigs database. Anvi'o is assigning `False` to the --profile-SCVs "
                             "flag, overruling your request like a boss.")
            self.profile_SCVs = False

        # we don't know what we are about
        self.description = None

        # additional layer data will be filled later
        self.layer_additional_keys = []
        self.layer_additional_data = {}


    def init_dirs_and_dbs(self):
        if not self.contigs_db_path:
            raise ConfigError("You can not run profiling without a contigs database. You can create "
                               "one using 'anvi-gen-contigs-database'. Not sure how? Please see the "
                               "tutorial: http://merenlab.org/2015/05/02/anvio-tutorial/")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'rU').read()

        if self.output_directory:
            self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=self.overwrite_output_destinations)
        else:
            output_dir_path = os.path.dirname(os.path.abspath(self.input_file_path))

            if self.sample_id:
                self.output_directory = filesnpaths.check_output_directory(os.path.join(output_dir_path, self.sample_id), ok_if_exists=self.overwrite_output_destinations)
            else:
                raise ConfigError("There is no `self.sample_id`, there is no `self.output_directory` :/ Anvi'o needs an adult :(")

        self.progress.new('Initializing')

        self.progress.update('Creating the output directory ...')
        filesnpaths.gen_output_directory(self.output_directory, self.progress, delete_if_exists=self.overwrite_output_destinations)

        self.progress.update('Creating a new single profile database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.profile_db_path = self.generate_output_destination('PROFILE.db')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        if self.blank:
            # if we are about to generate a blank profile, there is no
            # SNV profiling
            self.skip_SNV_profiling = True

        if self.skip_SNV_profiling:
            self.profile_SCVs = False
            self.skip_INDEL_profiling = True

        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': self.sample_id,
                       'merged': False,
                       'blank': self.blank,
                       'items_ordered': False,
                       'default_view': 'mean_coverage',
                       'min_contig_length': self.min_contig_length,
                       'max_contig_length': self.max_contig_length,
                       'SNVs_profiled': not self.skip_SNV_profiling,
                       'SCVs_profiled': self.profile_SCVs,
                       'INDELs_profiled': not self.skip_INDEL_profiling,
                       'min_percent_identity': self.min_percent_identity or 0,
                       'fetch_filter': self.fetch_filter,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'report_variability_full': self.report_variability_full,
                       'contigs_db_hash': self.a_meta['contigs_db_hash'],
                       'description': self.description if self.description else '_No description is provided_'}
        profile_db.create(meta_values)

        self.progress.update('Creating a new auxiliary database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.auxiliary_db_path = self.generate_output_destination('AUXILIARY-DATA.db')
        self.auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.auxiliary_db_path,
                                                                            self.a_meta['contigs_db_hash'],
                                                                            create_new=True,
                                                                            run=null_run,
                                                                            progress=null_progress)

        self.progress.end()

        if not self.skip_SNV_profiling:
            self.variable_nts_table = TableForVariability(self.profile_db_path, progress=null_progress)

        if self.profile_SCVs:
            self.variable_codons_table = TableForCodonFrequencies(self.profile_db_path, progress=null_progress)

        if not self.skip_INDEL_profiling:
            self.indels_table = TableForIndels(self.profile_db_path, progress=null_progress)


    def _run(self):
        self.check_args()

        self.set_sample_id()

        self.init_dirs_and_dbs()

        self.run.log_file_path = self.generate_output_destination('RUNLOG.txt')
        self.run.info('Sample name set', self.sample_id)
        self.run.info('Description', 'Found (%d characters)' % len(self.description) if self.description else None)
        self.run.info('Profile DB path', self.profile_db_path, display_only=True)
        self.run.info('Contigs DB path', self.contigs_db_path)
        self.run.info('Contigs DB hash', self.a_meta['contigs_db_hash'])
        self.run.info('Command line', utils.get_cmd_line(), align_long_values=False, nl_after=1)

        self.run.info('Minimum percent identity of reads to be profiled', self.min_percent_identity, mc='green')
        self.run.info('Fetch filter engaged', self.fetch_filter, mc='green', nl_after=1)

        self.run.info('Is merged profile?', False)
        self.run.info('Is blank profile?', self.blank)
        self.run.info('Skip contigs shorter than', self.min_contig_length)
        self.run.info('Skip contigs longer than', self.max_contig_length)
        self.run.info('Perform hierarchical clustering of contigs?', self.contigs_shall_be_clustered, nl_after=1)

        self.run.info('Profile single-nucleotide variants (SNVs)?', not self.skip_SNV_profiling)
        self.run.info('Profile single-codon variants (SCVs/+SAAVs)?', self.profile_SCVs)
        self.run.info('Profile insertion/deletions (INDELs)?', not self.skip_INDEL_profiling)
        self.run.info('Minimum coverage to calculate SNVs', self.min_coverage_for_variability)
        self.run.info('Report FULL variability data?', self.report_variability_full)

        self.run.warning("Your minimum contig length is set to %s base pairs. So anvi'o will not take into "
                         "consideration anything below that. If you need to kill this an restart your "
                         "analysis with another minimum contig length value, feel free to press CTRL+C." \
                                                % (pp(self.min_contig_length)))

        if self.max_contig_length < sys.maxsize:
            self.run.warning("Your maximum contig length is set to %s base pairs. Which means anvi'o will remove "
                             "any contigs that are longer than this value." % pp(self.max_contig_length))

        # this is kinda important. we do not run full-blown profile function if we are dealing with a summarized
        # profile...
        if self.blank:
            self.init_mock_profile()

            # creating a null view_data_splits dict:
            for view in self.essential_data_fields_for_anvio_profiles:
                for table_name in [f"{view}_{target}" for target in ['splits', 'contigs']]:
                    TablesForViews(self.profile_db_path).remove(view, table_names_to_blank=[table_name])
                    TablesForViews(self.profile_db_path,
                                   progress=self.progress) \
                                        .create_new_view(view_data=[],
                                                         table_name=table_name,
                                                         skip_sanity_check=True,
                                                         view_name=None if table_name.endswith('contigs') else view)
        elif self.input_file_path:
            self.init_profile_from_BAM()
            if self.num_threads > 1 or self.args.force_multi:
                self.profile_multi_thread()
            else:
                self.profile_single_thread()
        else:
            raise ConfigError("What are you doing? :( Whatever it is, anvi'o will have none of it.")

        # update layer additional data table content
        if self.layer_additional_data:
            self.progress.new("Additional layer data")
            self.progress.update("Updating the profile db ...")
            layer_additional_data_table = TableForLayerAdditionalData(argparse.Namespace(profile_db=self.profile_db_path), r=null_run, p=null_progress)
            layer_additional_data_table.add({self.sample_id: self.layer_additional_data}, self.layer_additional_keys)
            self.progress.end()

            self.run.info("Additional data added to the new profile DB", f"{', '.join(self.layer_additional_keys)}", nl_before=1)

        if self.contigs_shall_be_clustered:
            self.cluster_contigs()

        if self.bam:
            self.bam.close()

        self.run.info_single('Happy 😇', nl_before=1, nl_after=1)

        self.run.quit()


    def generate_variable_codons_table(self):
        if self.skip_SNV_profiling or not self.profile_SCVs:
            return

        for contig in self.contigs:
            for split in contig.splits:
                for gene_callers_id in split.SCV_profiles:
                    # We reorder to the profiles in the order they will appear in the output table
                    split.SCV_profiles[gene_callers_id] = OrderedDict(
                        [(col, split.SCV_profiles[gene_callers_id][col]) for col in t.variable_codons_table_structure]
                    )

                    entries = zip(*split.SCV_profiles[gene_callers_id].values())
                    for entry in entries:
                        self.variable_codons_table.append(entry)

        self.variable_codons_table.store()


    def generate_variable_nts_table(self):
        if self.skip_SNV_profiling:
            return

        for contig in self.contigs:
            for split in contig.splits:
                if split.num_SNV_entries == 0:
                    continue

                # We reorder to the profiles in the order they will appear in the output table
                split.SNV_profiles = OrderedDict(
                    [(col, split.SNV_profiles[col]) for col in t.variable_nts_table_structure]
                )

                entries = zip(*split.SNV_profiles.values())
                for entry in entries:
                    self.variable_nts_table.append(entry)

        self.variable_nts_table.store()


    def generate_indels_table(self):
        if self.skip_INDEL_profiling:
            return

        for contig in self.contigs:
            for split in contig.splits:
                for entry in split.INDEL_profiles.values():
                    self.indels_table.append([self.sample_id] + list(entry.values()))

        self.indels_table.store()


    def store_split_coverages(self):
        for contig in self.contigs:
            for split in contig.splits:
                self.auxiliary_db.append(split.name, self.sample_id, split.coverage.c)

        self.auxiliary_db.store()


    def set_sample_id(self):
        if self.sample_id:
            utils.check_sample_id(self.sample_id)
        else:
            if self.input_file_path:
                self.input_file_path = os.path.abspath(self.input_file_path)
                self.sample_id = os.path.splitext(os.path.basename(self.input_file_path))[0]
                self.sample_id = self.sample_id.replace('-', '_').replace('.', '_').replace(' ', '_')

                if self.sample_id[0] in constants.digits:
                    self.sample_id = 's' + self.sample_id

                if self.fetch_filter:
                    self.sample_id = f"{self.sample_id}_{self.fetch_filter.upper().replace('-', '_').replace('.', '_').replace(' ', '_')}"

                utils.check_sample_id(self.sample_id)
            if self.serialized_profile_path:
                self.serialized_profile_path = os.path.abspath(self.serialized_profile_path)
                self.sample_id = os.path.basename(os.path.dirname(self.serialized_profile_path))


    def check_contigs_without_any_gene_calls(self, contig_names):
        if not self.a_meta['genes_are_called']:
            self.run.warning("The contigs database '%s' does not contain any gene calls. Which means the profiling step "
                             "will not be able to characterize 'gene coverages'. If you are OK with this, anvi'o will be "
                             "OK with it as well." % (self.contigs_db_path))
            return

        contig_names = set(contig_names)
        contigs_without_any_gene_calls = [c for c in contig_names if c not in self.contig_name_to_genes]

        if len(contigs_without_any_gene_calls):
            import random
            P = lambda x: 'are %d contigs' % (x) if x > 1 else 'there is one contig'
            self.run.warning('According to the data generated in the contigs database, there %s in your BAM file '
                             'with 0 gene calls. Which may not be unusual if (a) some of your contigs are very short, '
                             'or (b) your the gene caller was not capable of dealing with the type of data you had. '
                             'If you would like to take a look yourself, here is one contig that is missing any genes: %s"' %\
                                      (P(len(contigs_without_any_gene_calls)), random.choice(contigs_without_any_gene_calls)))


    def list_contigs(self):
        import signal
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = bamops.BAMFileObject(self.input_file_path, 'rb')
        self.bam.fetch_filter = self.fetch_filter
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lengths = self.bam.lengths

        utils.check_contig_names(self.contig_names)

        for tpl in sorted(zip(self.contig_lengths, self.contig_names), reverse=True):
            print('%-40s %s' % (tpl[1], pp(int(tpl[0]))))


    def remove_contigs_based_on_min_max_contig_length(self):
        """Removes contigs that are shorter or longer than user defined values"""

        contigs_with_good_lengths = set()

        for i in range(0, len(self.contig_names)):
            if self.contig_lengths[i] >= self.min_contig_length and self.contig_lengths[i] <= self.max_contig_length:
                contigs_with_good_lengths.add(i)

        if not len(contigs_with_good_lengths):
            try:
                filesnpaths.shutil.rmtree(self.output_directory)
            except:
                pass

            raise ConfigError("Anvi'o applied your min/max length criteria for contigs to filter out the bad ones "
                              "and has bad news: not a single contig in your contigs database was greater than %s "
                              "and smaller than %s nts :( So this profiling attempt did not really go anywhere. "
                              "Please remove your half-baked output directory if it is still there: '%s'." \
                                        % (pp(self.min_contig_length), pp(self.max_contig_length), self.output_directory))
        else:
            self.contig_names = [self.contig_names[i] for i in contigs_with_good_lengths]
            self.contig_lengths = [self.contig_lengths[i] for i in contigs_with_good_lengths]
            self.num_contigs = len(self.contig_names)    # we will store these two
            self.total_length = sum(self.contig_lengths) # into the db in a second.

        contigs_with_good_lengths = set(self.contig_names) # for fast access
        self.split_names = set([])
        self.contig_name_to_splits = {}
        for split_name in sorted(self.splits_basic_info.keys()):
            parent = self.splits_basic_info[split_name]['parent']

            if parent not in contigs_with_good_lengths:
                continue

            self.split_names.add(split_name)

            if parent in self.contig_name_to_splits:
                self.contig_name_to_splits[parent].append(split_name)
            else:
                self.contig_name_to_splits[parent] = [split_name]

        self.num_splits = len(self.split_names)


    def init_profile_from_BAM(self):
        filesnpaths.is_file_bam_file(self.input_file_path)

        self.progress.new('Init')
        self.progress.update('Reading BAM File')

        try:
            self.bam = bamops.BAMFileObject(self.input_file_path)
        except Exception as e:
            raise ConfigError(f"Sorry, this BAM file does not look like a BAM file :( Here is "
                              f"the complaint coming from the depths of the codebase: '{e}'.")

        self.bam.fetch_filter = self.fetch_filter
        self.num_reads_mapped = self.bam.mapped
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lengths = self.bam.lengths

        if self.min_percent_identity:
            # We will permit as many as 1% of sampled reads to be missing an MD tag, for whatever
            # crazy mapping reason, although it is most probable that if 1 read is missing an MD
            # tag, they all are.
            requirement = lambda x: sum(x) / len(x) >= 0.99

            if not self.bam.reads_have_MD_tags(require=requirement):
                raise ConfigError("Your BAM file has reads that do not have MD "
                                  "tags, which give essential information about the "
                                  "reference bases in alignments. Unfortunately, "
                                  "this is required to use --min-percent-identity. "
                                  "Please remove this flag or redo your mapping.")

        utils.check_contig_names(self.contig_names)

        self.run.info('Input BAM', self.input_file_path)
        self.run.info('Output directory path', self.output_directory, display_only=True, nl_after=1)

        self.run.info('Number of reads in the BAM file', pp(int(self.num_reads_mapped)))
        self.run.info('Number of sequences in the contigs DB', pp(len(self.contig_names)))

        if self.contig_names_of_interest:
            indexes = [self.contig_names.index(r) for r in self.contig_names_of_interest if r in self.contig_names]
            self.contig_names = [self.contig_names[i] for i in indexes]
            self.contig_lengths = [self.contig_lengths[i] for i in indexes]
            self.run.info('Number of contigs selected for analysis', pp(len(self.contig_names)), mc='green')

        # it brings good karma to let the user know what the hell is wrong with their data:
        self.check_contigs_without_any_gene_calls(self.contig_names)

        # check for the -M parameter.
        self.remove_contigs_based_on_min_max_contig_length()

        # let's see whether the user screwed up to follow the simple instructions
        # mentioned here: http://merenlab.org/2015/05/01/anvio-tutorial/#preparation
        for contig_name in self.contig_names:
            if contig_name not in self.contig_names_in_contigs_db:
                raise ConfigError("At least one contig name in your BAM file does not match contig names stored in the "
                                   "contigs database. For instance, this is one contig name found in your BAM file: '%s', "
                                   "and this is another one found in your contigs database: '%s'. You may be using an "
                                   "contigs database for profiling that has nothing to do with the BAM file you are "
                                   "trying to profile, or you may have failed to fix your contig names in your FASTA file "
                                   "prior to mapping, which is described here: %s"\
                                        % (contig_name, self.contig_names_in_contigs_db.pop(), 'http://goo.gl/Q9ChpS'))

        self.run.info('Number of contigs to be conisdered (after -M)', self.num_contigs, display_only=True)
        self.run.info('Number of splits', self.num_splits)
        self.run.info('Number of nucleotides', self.total_length)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('num_splits', self.num_splits)
        profile_db.db.set_meta_value('num_contigs', self.num_contigs)
        profile_db.db.set_meta_value('total_length', self.total_length)
        profile_db.disconnect()

        self.layer_additional_data['total_reads_mapped'] = self.num_reads_mapped
        self.layer_additional_keys.append('total_reads_mapped')


    def init_mock_profile(self):
        self.progress.new('Init')
        self.progress.update('...')
        self.num_reads_mapped = 0
        self.progress.end()

        self.contig_names = list(self.contigs_basic_info.keys())
        self.contig_lengths = [self.contigs_basic_info[contig_name]['length'] for contig_name in self.contigs_basic_info]
        self.total_length = sum(self.contig_lengths)
        self.num_contigs = len(self.contig_names)

        utils.check_contig_names(self.contig_names)

        self.run.info('input_bam', None)
        self.run.info('output_dir', self.output_directory, display_only=True)

        # check for the -M parameter.
        self.remove_contigs_based_on_min_max_contig_length()

        self.run.info('num_contigs_after_M', self.num_contigs, display_only=True)
        self.run.info('num_contigs', self.num_contigs, quiet=True)
        self.run.info('num_splits', self.num_splits)
        self.run.info('total_length', self.total_length)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('num_splits', self.num_splits)
        profile_db.db.set_meta_value('num_contigs', self.num_contigs)
        profile_db.db.set_meta_value('total_length', self.total_length)
        profile_db.disconnect()

        self.layer_additional_data['total_reads_mapped'] = self.num_reads_mapped
        self.layer_additional_keys.append('total_reads_mapped')


    def generate_output_destination(self, postfix, directory=False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def populate_gene_info_for_splits(self, contig):
        """Stores info available to the ContigsSuperClass into contigops.Split objects

        Populates Split.per_position_info as a dictionary of arrays, each with length equal to the
        split.
        """

        if self.skip_SNV_profiling:
            return

        info_of_interest = [
            'corresponding_gene_call',
            'codon_order_in_gene',
            'base_pos_in_codon',
            'in_noncoding_gene_call',
            'in_coding_gene_call',
        ]

        if self.profile_SCVs:
            # NOTE The people want SCVs, and that means Split needs to know about gene calls. Rather
            #      than giving each split a dictionary of gene calls, in an ugly but much faster
            #      fashion, we add these 3 pieces of information which provide sufficient gene call
            #      information to calculate SCV. Previously, we parsed the genes_in_contigs dict
            #      which is fine for small databases, but the cost is enormous for very large
            #      (100,000+ contigs) databases. To give perspective, it was about 1s per contig,
            #      which is 30 hours for 100,000 contigs. That code elegantly looked like this:
            #
            #      for split in contig.splits:
            #          gene_ids_in_split = set(gc['gene_callers_id']
            #                                  for key, gc in self.genes_in_splits.items()
            #                                  if key in self.split_name_to_genes_in_splits_entry_ids[split.name])
            #          split.gene_calls = {gene_id: self.genes_in_contigs_dict[gene_id] for gene_id in gene_ids_in_split}
            info_of_interest.extend([
                'forward',
                'gene_start',
                'gene_stop',
            ])

        nt_info = self.get_gene_info_for_each_position(contig.name, info=info_of_interest)

        for split in contig.splits:
            for info in info_of_interest:
                split.per_position_info[info] = nt_info[info][split.start:split.end]


    @staticmethod
    def profile_contig_worker(self, available_index_queue, output_queue):
        bam_file = bamops.BAMFileObject(self.input_file_path)
        bam_file.fetch_filter = self.fetch_filter

        while True:
            try:
                index = available_index_queue.get(True)
                contig_name = self.contig_names[index]
                contig_length = self.contig_lengths[index]

                contig = self.process_contig(bam_file, contig_name, contig_length)
                output_queue.put(contig)

                if contig is not None:
                    # We mark these for deletion the next time garbage is collected
                    for split in contig.splits:
                        del split.coverage
                        del split.auxiliary
                        del split
                    del contig.splits[:]
                    del contig.coverage
                    del contig
            except Exception as e:
                # This thread encountered an error. We send the error back to the main thread which
                # will terminate the job.
                output_queue.put(e)

        # we are closing this object here for clarity, although we are not really closing it since
        # the code never reaches here and the worker is killed by its parent:
        bam_file.close()
        return


    def process_contig(self, bam_file, contig_name, contig_length):
        timer = terminal.Timer(initial_checkpoint_key='Start')

        contig = contigops.Contig(contig_name)
        contig.length = contig_length
        contig.split_length = self.a_meta['split_length']
        contig.skip_SNV_profiling = self.skip_SNV_profiling
        timer.make_checkpoint('Initialization done')

        # populate contig with empty split objects
        for split_name in self.contig_name_to_splits[contig_name]:
            s = self.splits_basic_info[split_name]
            split_sequence = self.contig_sequences[contig_name]['sequence'][s['start']:s['end']]
            split = contigops.Split(split_name, split_sequence, contig_name, s['order_in_parent'], s['start'], s['end'])
            contig.splits.append(split)

        timer.make_checkpoint('Split objects initialized')

        self.populate_gene_info_for_splits(contig)
        timer.make_checkpoint('Gene info for split added')

        # analyze coverage for each split
        contig.analyze_coverage(bam_file, self.min_percent_identity)
        timer.make_checkpoint('Coverage done')

        if not self.skip_SNV_profiling:
            for split in contig.splits:
                split.auxiliary = contigops.Auxiliary(split,
                                                      profile_SCVs=self.profile_SCVs,
                                                      skip_INDEL_profiling=self.skip_INDEL_profiling,
                                                      skip_SNV_profiling=self.skip_SNV_profiling,
                                                      min_coverage_for_variability=self.min_coverage_for_variability,
                                                      report_variability_full=self.report_variability_full,
                                                      min_percent_identity=self.min_percent_identity,
                                                      skip_edges=self.skip_edges)

                split.auxiliary.process(bam_file)

                if split.num_SNV_entries == 0:
                    continue

                # Add these redundant data ad-hoc
                split.SNV_profiles['split_name'] = [split.name] * split.num_SNV_entries
                split.SNV_profiles['sample_id'] = [self.sample_id] * split.num_SNV_entries
                split.SNV_profiles['pos_in_contig'] = split.SNV_profiles['pos'] + split.start

                for gene_id in split.SCV_profiles:
                    split.SCV_profiles[gene_id]['sample_id'] = [self.sample_id] * split.num_SCV_entries[gene_id]
                    split.SCV_profiles[gene_id]['corresponding_gene_call'] = [gene_id] * split.num_SCV_entries[gene_id]

            timer.make_checkpoint('Auxiliary analyzed')

        # output_queue.put(contig) is an expensive operation that does not handle large data
        # structures well. So we delete everything we can
        del contig.coverage.c # only split coverage array is needed
        for split in contig.splits:
            del split.per_position_info

        if anvio.DEBUG:
            timer.gen_report('%s Time Report' % contig.name)

        return contig


    def profile_single_thread(self):
        """The main method for anvi-profile when num_threads is 1"""

        bam_file = bamops.BAMFileObject(self.input_file_path)
        bam_file.fetch_filter = self.fetch_filter

        received_contigs = 0

        self.progress.new('Profiling w/1 thread', progress_total_items=self.num_contigs)

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        self.progress.update('contigs are being processed ...')
        for index in range(self.num_contigs):
            contig_name = self.contig_names[index]
            contig_length = self.contig_lengths[index]

            contig = self.process_contig(bam_file, contig_name, contig_length)

            self.contigs.append(contig)

            received_contigs += 1

            if mem_tracker.measure():
                mem_usage = mem_tracker.get_last()
                mem_diff = mem_tracker.get_last_diff()

            self.progress.increment(received_contigs)
            msg = '%d/%d contigs ⚙  | MEMORY 🧠  %s (%s)' % (received_contigs, self.num_contigs, mem_usage, mem_diff)
            self.progress.update(msg)

            # Here you're about to witness the poor side of Python (or our use of it). Although
            # we couldn't find any refs to these objects, garbage collecter kept them in the
            # memory. So here we are accessing to the atomic data structures in our split
            # objects to try to relieve the memory by encouraging the garbage collector to
            # realize what's up. Afterwards, we explicitly call the garbage collector
            if self.write_buffer_size > 0 and len(self.contigs) % self.write_buffer_size == 0:
                self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
                self.store_contigs_buffer()
                for c in self.contigs:
                    for split in c.splits:
                        del split.coverage
                        del split.auxiliary
                        del split
                    del c.splits[:]
                    del c.coverage
                    del c
                del self.contigs[:]
                gc.collect()

        self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
        self.store_contigs_buffer()
        self.auxiliary_db.close()

        self.progress.end(timing_filepath='anvio.debug.timing.txt' if anvio.DEBUG else None)

        overall_mean_coverage = 1
        if self.total_length_of_all_contigs != 0:
            overall_mean_coverage = self.total_coverage_values_for_all_contigs / self.total_length_of_all_contigs

        # FIXME: We know this is ugly. You can keep your opinion to yourself.
        if overall_mean_coverage > 0.0:
            # avoid dividing by zero
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_splits SET value = value / " + str(overall_mean_coverage) + " * 1.0;")
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_contigs SET value = value / " + str(overall_mean_coverage) + " * 1.0;")

        if not self.skip_SNV_profiling:
            self.layer_additional_data['num_SNVs_reported'] = self.variable_nts_table.num_entries
            self.layer_additional_keys.append('num_SNVs_reported')
            self.run.info("Num SNVs reported", self.layer_additional_data['num_SNVs_reported'], nl_before=1)

        if not self.skip_INDEL_profiling:
            self.layer_additional_data['num_INDELs_reported'] = self.indels_table.num_entries
            self.layer_additional_keys.append('num_INDELs_reported')
            self.run.info("Num INDELs reported", self.layer_additional_data['num_INDELs_reported'])

        if self.profile_SCVs:
            self.layer_additional_data['num_SCVs_reported'] = self.variable_codons_table.num_entries
            self.layer_additional_keys.append('num_SCVs_reported')
            self.run.info("Num SCVs reported", self.layer_additional_data['num_SCVs_reported'])

        if self.total_reads_kept != self.num_reads_mapped:
            # Num reads in profile do not equal num reads in bam
            diff = self.num_reads_mapped - self.total_reads_kept
            perc = (1 - self.total_reads_kept / self.num_reads_mapped) * 100

            if self.contig_names_of_interest:
                self.run.warning(f"There were {pp(diff)} reads present in the BAM file that did not end up being used "
                                 f"by the profiler, which corresponds to about {perc}% of all reads. This is "
                                 f"most likely due to the parameter `--contigs-of-interest`. Regardless, anvi'o "
                                 f"thought you should know about this.")
            elif self.fetch_filter:
                self.run.warning(f"There were {pp(diff)} reads present in the BAM file that did not end up being used "
                                 f"by the profiler, which corresponds to about {perc}% of all reads. This is "
                                 f"most likely due to your `--fetch-filter`, which, depending on the filter, can "
                                 f"make use of a very tiny fraction of all reads. If you think this is much more or "
                                 f"much less than what you would have expected, you may want to investigate it.")
            else:
                self.run.warning(f"There were {pp(diff)} reads present in the BAM file that did not end up being used "
                                 f"by the profiler, which corresponds to about {perc}% of all reads. Since you "
                                 f"don't seem to have used parameters such as `--contigs-of-interest` or `--fetch-filter` "
                                 f"anvi'o is Jon Snow here. There are other reasons why some of your reads may end up "
                                 f"not being considered by the profiler. For instance, if pysam encounteres reads that "
                                 f"it can't deal with (i.e., mapped reads in the BAM file with no defined sequences, or "
                                 f"read with defined sequences without mapping). Another reason can be that you have used "
                                 f"a new paramter in the profiler that removes contigs or reads from consideration. "
                                 f"Regardless, if you think this is worth your attention, now you know.")

        self.layer_additional_data['total_reads_kept'] = self.total_reads_kept
        self.layer_additional_keys.append('total_reads_kept')

        self.check_contigs(num_contigs=received_contigs)


    def profile_multi_thread(self):
        """The main method for anvi-profile when num_threads is >1"""

        manager = multiprocessing.Manager()
        available_index_queue = manager.Queue()
        output_queue = manager.Queue(self.queue_size)

        # put contig indices into the queue to be read from within the worker
        for i in range(0, self.num_contigs):
            available_index_queue.put(i)

        processes = []
        for i in range(0, self.num_threads):
            processes.append(multiprocessing.Process(target=BAMProfiler.profile_contig_worker, args=(self, available_index_queue, output_queue)))

        for proc in processes:
            proc.start()

        received_contigs = 0

        self.progress.new('Profiling w/%d threads' % self.num_threads, progress_total_items=self.num_contigs)
        self.progress.update('initializing threads ...')

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        self.progress.update('contigs are being processed ...')
        while received_contigs < self.num_contigs:
            try:
                contig = output_queue.get()

                if isinstance(contig, Exception):
                    # If thread returns an exception, we raise it and kill the main thread.
                    raise contig

                self.contigs.append(contig)

                received_contigs += 1

                if mem_tracker.measure():
                    mem_usage = mem_tracker.get_last()
                    mem_diff = mem_tracker.get_last_diff()

                self.progress.increment(received_contigs)
                self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | MEMORY 🧠  {mem_usage} ({mem_diff}) ...")

                # Here you're about to witness the poor side of Python (or our use of it). Although
                # we couldn't find any refs to these objects, garbage collecter kept them in the
                # memory. So here we are accessing to the atomic data structures in our split
                # objects to try to relieve the memory by encouraging the garbage collector to
                # realize what's up. Afterwards, we explicitly call the garbage collector
                if self.write_buffer_size > 0 and len(self.contigs) % self.write_buffer_size == 0:
                    self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
                    self.store_contigs_buffer()
                    for c in self.contigs:
                        for split in c.splits:
                            del split.coverage
                            del split.auxiliary
                            del split
                        del c.splits[:]
                        del c.coverage
                        del c
                    del self.contigs[:]
                    gc.collect()

            except KeyboardInterrupt:
                self.run.info_single("Anvi'o profiler received SIGINT, terminating all processes...", nl_before=2)
                break

            except Exception as worker_error:
                # An exception was thrown in one of the profile workers. We kill all processes in this case
                self.progress.end()
                for proc in processes:
                    proc.terminate()
                raise worker_error

        for proc in processes:
            proc.terminate()

        self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
        self.store_contigs_buffer()
        self.auxiliary_db.close()

        self.progress.end(timing_filepath='anvio.debug.timing.txt' if anvio.DEBUG else None)

        overall_mean_coverage = 1
        if self.total_length_of_all_contigs != 0:
            overall_mean_coverage = self.total_coverage_values_for_all_contigs / self.total_length_of_all_contigs

        # FIXME: We know this is ugly. You can keep your opinion to yourself.
        if overall_mean_coverage > 0.0:
            # avoid dividing by zero
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_splits SET value = value / " + str(overall_mean_coverage) + " * 1.0;")
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_contigs SET value = value / " + str(overall_mean_coverage) + " * 1.0;")

        if not self.skip_SNV_profiling:
            self.layer_additional_data['num_SNVs_reported'] = self.variable_nts_table.num_entries
            self.layer_additional_keys.append('num_SNVs_reported')
            self.run.info("Num SNVs reported", self.layer_additional_data['num_SNVs_reported'], nl_before=1)

        if not self.skip_INDEL_profiling:
            self.layer_additional_data['num_INDELs_reported'] = self.indels_table.num_entries
            self.layer_additional_keys.append('num_INDELs_reported')
            self.run.info("Num INDELs reported", self.layer_additional_data['num_INDELs_reported'])

        if self.profile_SCVs:
            self.layer_additional_data['num_SCVs_reported'] = self.variable_codons_table.num_entries
            self.layer_additional_keys.append('num_SCVs_reported')
            self.run.info("Num SCVs reported", self.layer_additional_data['num_SCVs_reported'])

        if self.total_reads_kept != self.num_reads_mapped:
            # Num reads in profile do not equal num reads in bam
            diff = self.num_reads_mapped - self.total_reads_kept
            perc = (1 - self.total_reads_kept / self.num_reads_mapped) * 100
            self.run.warning("There were %d reads present in the BAM file that did not end up being used "
                             "by anvi'o. That corresponds to about %.2f percent of all reads in the bam file. "
                             "This could be either because you supplied --contigs-of-interest, "
                             "or because pysam encountered reads it could not deal with, e.g. they mapped "
                             "but had no defined sequence, or they had a sequence but did not map. "
                             "Regardless, anvi'o thought you should be aware of this." % (diff, perc))

        self.layer_additional_data['total_reads_kept'] = self.total_reads_kept
        self.layer_additional_keys.append('total_reads_kept')

        self.check_contigs(num_contigs=received_contigs)


    def store_contigs_buffer(self):
        for contig in self.contigs:
            self.total_length_of_all_contigs += contig.length
            self.total_coverage_values_for_all_contigs += contig.coverage.mean * contig.length

            # we will divide every abundance after profiling is done.
            contig.abundance = contig.coverage.mean
            for split in contig.splits:
                split.abundance = split.coverage.mean

            self.total_reads_kept += contig.coverage.num_reads

        self.generate_variable_nts_table()
        self.generate_variable_codons_table()
        self.generate_indels_table()
        self.store_split_coverages()

        # the crux of the profiling
        for atomic_data_field in self.essential_data_fields_for_anvio_profiles:
            view_data_splits, view_data_contigs = contigops.get_atomic_data(self.sample_id, self.contigs, atomic_data_field)

            table_name = '_'.join([atomic_data_field, 'splits'])
            TablesForViews(self.profile_db_path,
                           progress=self.progress) \
                                .create_new_view(view_data=view_data_splits,
                                                 table_name=table_name,
                                                 view_name=atomic_data_field,
                                                 append_mode=True)

            table_name = '_'.join([atomic_data_field, 'contigs'])
            TablesForViews(self.profile_db_path,
                           progress=self.progress) \
                                .create_new_view(view_data=view_data_contigs,
                                                 table_name=table_name,
                                                 view_name=None,
                                                 append_mode=True)


    def check_contigs(self, num_contigs=None):
        if not num_contigs:
            num_contigs = len(self.contigs)

        if not num_contigs:
            raise ConfigError("0 contigs to work with. Bye.")


    def cluster_contigs(self):
        default_clustering_config = constants.blank_default if self.blank else constants.single_default

        dbops.do_hierarchical_clustering_of_items(self.profile_db_path, self.clustering_configs, self.split_names, self.database_paths, \
                                                  input_directory=self.output_directory, default_clustering_config=default_clustering_config, \
                                                  distance=self.distance, linkage=self.linkage, run=self.run, progress=self.progress)


    def check_args(self):
        if self.blank:
            self.run.warning("You are about to generate a blank profile. This is what we do when we have nothing "
                             "but a contigs database to play with. Because anvi'o is lazy, it will not check the "
                             "rest of the parameters you may have declred. Most of them will not matter.")

            if not self.output_directory:
                raise ConfigError("If you want to generate a blank profile, you need to declare an output directory path.")
            if not self.sample_id:
                raise ConfigError("Mock profiles require a sample name to be declared. Because :/")
            return

        if (not self.input_file_path) and (not self.serialized_profile_path):
            raise ConfigError("You didn't declare any input files :/ If you intend to create a blank profile without any,\
                                input file, you should be a bit more explicit about your intention (you know, in the help\
                                there is a flag for it and all). Otherwise you should either provide an input BAM file, or\
                                a serialized anvi'o profile. See '--help' maybe?")
        if self.input_file_path and self.serialized_profile_path:
            raise ConfigError("You can't declare both an input file and a serialized profile.")
        if self.serialized_profile_path and (not self.output_directory):
            raise ConfigError("When loading serialized profiles, you need to declare an output directory.")
        if self.input_file_path and not os.path.exists(self.input_file_path):
            raise ConfigError("No such file: '%s'" % self.input_file_path)
        if self.serialized_profile_path and not os.path.exists(self.serialized_profile_path):
            raise ConfigError("No such file: '%s'" % self.serialized_profile_path)
        if not self.min_coverage_for_variability >= 0:
            raise ConfigError("Minimum coverage for variability must be 0 or larger.")
        if not self.min_contig_length >= 0:
            raise ConfigError("Minimum contig length must be 0 or larger.")
        if not self.max_contig_length >= 100:
            raise ConfigError("Maximum contig length can't be less than 100 base pairs.")
        if self.min_contig_length >= self.max_contig_length:
            raise ConfigError("Maximum contig length (%s) must be larger than the minimum "
                              "contig length (%s). Seriously though." % (pp(self.max_contig_length), pp(self.min_contig_length)))

        if self.num_threads < 1:
            raise ConfigError("Nice try. Obviously, number of threds can not be less than 1.")

        if not self.queue_size:
            self.queue_size = self.num_threads * 2

        if not self.write_buffer_size:
            self.run.warning("You set the write buffer size to 0. Which means, the profiling data will be kept in memory until "
                             "the very end of the processing.")

        if self.write_buffer_size < 0:
            raise ConfigError('No. Write buffer size can not have a negative value.')
