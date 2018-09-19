# coding: utf-8
# pylint: disable=line-too-long
"""Provides the necessary class to profile BAM files."""

import os
import sys
import time
import pysam
import shutil
import argparse
import multiprocessing

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
from anvio.tables.codonfrequencies import TableForCodonFrequencies
from anvio.tables.variability import TableForVariability
from anvio.tables.miscdata import TableForLayerAdditionalData


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


pp = terminal.pretty_print

class BAMProfiler(dbops.ContigsSuperclass):
    """Creates an über class for BAM file operations"""
    def __init__(self, args):
        self.args = args

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_file_path = A('input_file')
        self.contigs_db_path = A('contigs_db')
        self.serialized_profile_path = A('serialized_profile')
        self.output_directory = A('output_dir')
        self.list_contigs_and_exit = A('list_contigs')
        self.min_contig_length = A('min_contig_length')
        self.max_contig_length = A('max_contig_length') or sys.maxsize
        self.min_mean_coverage = A('min_mean_coverage')
        self.min_coverage_for_variability = A('min_coverage_for_variability')
        self.contigs_shall_be_clustered = A('cluster_contigs')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.sample_id = A('sample_name')
        self.report_variability_full = A('report_variability_full')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.skip_SNV_profiling = A('skip_SNV_profiling')
        self.profile_SCVs = A('profile_SCVs')
        self.gen_serialized_profile = A('gen_serialized_profile')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.num_threads = int(A('num_threads'))
        self.queue_size = int(A('queue_size'))
        self.write_buffer_size = int(A('write_buffer_size'))
        self.total_length_of_all_contigs = 0
        self.total_coverage_values_for_all_contigs = 0
        self.description_file_path = A('description')

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        # whehther the profile database is a blank (without any BAM files or reads):
        self.blank = A('blank_profile')

        if not self.blank and self.contigs_shall_be_clustered and self.skip_hierarchical_clustering:
            raise ConfigError("You are confused, and confusing anvi'o, too. You can't as hierarchical clustering\
                               to be performed with one flag, and try to skip it with another one :(")

        if self.blank and self.contigs_shall_be_clustered and self.skip_hierarchical_clustering:
            raise ConfigError("So you want to generate a blank profile, and you both want hierarchical clustering\
                               of your contigs to be performed, and skipped. No.")

        if self.blank and self.contigs_shall_be_clustered:
            raise ConfigError("When the blank profile is asked to be generated, there is no need to ask for the\
                               hierarchical clustering of contigs. It is going to be done by default. If it is\
                               not changing anything, why is anvi'o upset with you? Because. Let's don't use flags\
                               we don't need.")

        if self.blank and not self.skip_hierarchical_clustering:
            self.contigs_shall_be_clustered = True

        if args.contigs_of_interest:
            filesnpaths.is_file_exists(args.contigs_of_interest)
            self.contig_names_of_interest = set([c.strip() for c in open(args.contigs_of_interest).readlines()\
                                                                           if c.strip() and not c.startswith('#')])
        else:
            self.contig_names_of_interest = None

        self.progress = terminal.Progress()
        self.run = terminal.Run(width=35)

        if self.list_contigs_and_exit:
            self.list_contigs()
            sys.exit()

        if not self.contigs_db_path:
            raise ConfigError("No contigs database, no profilin'. Bye.")

        # Initialize contigs db
        dbops.ContigsSuperclass.__init__(self, self.args, r=self.run, p=self.progress)
        self.init_contig_sequences()
        self.contig_names_in_contigs_db = set(self.contigs_basic_info.keys())

        self.bam = None
        self.contigs = []

        self.database_paths = {'CONTIGS.db': os.path.abspath(self.contigs_db_path)}

        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['blank' if self.blank else 'single']

        # following variable will be populated during the profiling, and its content will eventually
        # be stored in t.variable_nts_table_name
        self.variable_nts_table_entries = []

        # if genes are not called, yet the user is asking for codon frequencies to be profiled, we give
        # a warning and force-turn that flag off.
        if (not self.a_meta['genes_are_called']) and self.profile_SCVs:
            self.run.warning("You asked the codon frequencies to be profiled, but genes were not called\
                              for your contigs database. Anvi'o is assigning `False` to the profile-codon-frequncies\
                              flag, overruling your request like a boss.")
            self.profile_SCVs = False

        # following variable will be populated while the variable positions table is computed
        self.codons_in_genes_to_profile_SCVs = set([])

        # we don't know what we are about
        self.description = None

        # additional layer data will be filled later
        self.layer_additional_keys = []
        self.layer_additional_data = {}


    def init_dirs_and_dbs(self):
        if not self.contigs_db_path:
            raise ConfigError("You can not run profiling without a contigs database. You can create\
                                one using 'anvi-gen-contigs-database'. Not sure how? Please see the\
                                tutorial: http://merenlab.org/2015/05/02/anvio-tutorial/")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'rU').read()

        self.output_directory = filesnpaths.check_output_directory(self.output_directory or self.input_file_path + '-ANVIO_PROFILE',\
                                                                   ok_if_exists=self.overwrite_output_destinations)

        self.progress.new('Initializing')

        self.progress.update('Creating the output directory ...')
        filesnpaths.gen_output_directory(self.output_directory, self.progress, delete_if_exists=self.overwrite_output_destinations)

        self.progress.update('Creating a new single profile database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.profile_db_path = self.generate_output_destination('PROFILE.db')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        if self.skip_SNV_profiling:
            self.profile_SCVs = False

        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': self.sample_id,
                       'merged': False,
                       'blank': self.blank,
                       'contigs_ordered': self.contigs_shall_be_clustered,
                       'default_view': 'single',
                       'min_contig_length': self.min_contig_length,
                       'max_contig_length': self.max_contig_length,
                       'SNVs_profiled': not self.skip_SNV_profiling,
                       'SCVs_profiled': self.profile_SCVs,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'report_variability_full': self.report_variability_full,
                       'contigs_db_hash': self.a_meta['contigs_db_hash'],
                       'description': self.description if self.description else '_No description is provided_'}
        profile_db.create(meta_values)

        self.progress.update('Creating a new auxiliary database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.auxiliary_db_path = self.generate_output_destination('AUXILIARY-DATA.db')
        self.auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.auxiliary_db_path, self.a_meta['contigs_db_hash'], create_new=True)

        self.progress.end()

        if self.skip_SNV_profiling:
            self.run.warning('Single-nucleotide variation will not be characterized for this profile.')

        if not self.profile_SCVs:
            self.run.warning('Amino acid linkmer frequencies will not be characterized for this profile.')


    def _run(self):
        self.check_args()

        self.set_sample_id()

        self.init_dirs_and_dbs()

        self.run.log_file_path = self.generate_output_destination('RUNLOG.txt')
        self.run.info('anvio', anvio.__version__)
        self.run.info('profiler_version', anvio.__profile__version__)
        self.run.info('sample_id', self.sample_id)
        self.run.info('description', 'Found (%d characters)' % len(self.description) if self.description else None)
        self.run.info('profile_db', self.profile_db_path, display_only=True)
        self.run.info('contigs_db', True if self.contigs_db_path else False)
        self.run.info('contigs_db_hash', self.a_meta['contigs_db_hash'])
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('merged', False)
        self.run.info('blank', self.blank)
        self.run.info('split_length', self.a_meta['split_length'])
        self.run.info('min_contig_length', self.min_contig_length)
        self.run.info('max_contig_length', self.max_contig_length)
        self.run.info('min_mean_coverage', self.min_mean_coverage)
        self.run.info('clustering_performed', self.contigs_shall_be_clustered)
        self.run.info('min_coverage_for_variability', self.min_coverage_for_variability)
        self.run.info('skip_SNV_profiling', self.skip_SNV_profiling)
        self.run.info('profile_SCVs', self.profile_SCVs)
        self.run.info('report_variability_full', self.report_variability_full)

        self.run.warning("Your minimum contig length is set to %s base pairs. So anvi'o will not take into\
                          consideration anything below that. If you need to kill this an restart your\
                          analysis with another minimum contig length value, feel free to press CTRL+C." \
                                                % (pp(self.min_contig_length)))

        if self.max_contig_length < sys.maxsize:
            self.run.warning("Your maximum contig length is set to %s base pairs. Which means anvi'o will remove\
            any contigs that are longer than this value." % pp(self.max_contig_length))

        # this is kinda important. we do not run full-blown profile function if we are dealing with a summarized
        # profile...
        if self.blank:
            self.init_mock_profile()

            # creating a null view_data_splits dict:
            view_data_splits = dict(list(zip(self.split_names, [dict(list(zip(t.atomic_data_table_structure[1:], [None] * len(t.atomic_data_table_structure[1:]))))] * len(self.split_names))))
            TablesForViews(self.profile_db_path).remove('single', table_names_to_blank=['atomic_data_splits'])
            TablesForViews(self.profile_db_path).create_new_view(
                                           data_dict=view_data_splits,
                                           table_name='atomic_data_splits',
                                           table_structure=t.atomic_data_table_structure,
                                           table_types=t.atomic_data_table_types,
                                           view_name='single')
        elif self.input_file_path:
            self.init_profile_from_BAM()
            self.profile()
        else:
            raise ConfigError("What are you doing? :( Whatever it is, anvi'o will have none of it.")

        # update layer additional data table content
        if self.layer_additional_data:
            layer_additional_data_table = TableForLayerAdditionalData(argparse.Namespace(profile_db=self.profile_db_path))
            layer_additional_data_table.add({self.sample_id: self.layer_additional_data}, self.layer_additional_keys)

        if self.contigs_shall_be_clustered:
            self.cluster_contigs()

        if self.bam:
            self.bam.close()

        self.run.quit()


    def generate_variabile_codons_table(self):
        if self.skip_SNV_profiling or not self.profile_SCVs:
            return

        variable_codons_table = TableForCodonFrequencies(self.profile_db_path, progress=self.progress)

        codon_frequencies = bamops.CodonFrequencies()

        codons_in_genes_to_profile_SCVs_dict = {}
        for gene_callers_id, codon_order in self.codons_in_genes_to_profile_SCVs:
            if gene_callers_id not in codons_in_genes_to_profile_SCVs_dict:
                codons_in_genes_to_profile_SCVs_dict[gene_callers_id] = set([])
            codons_in_genes_to_profile_SCVs_dict[gene_callers_id].add(codon_order)

        gene_caller_ids_to_profile = list(codons_in_genes_to_profile_SCVs_dict.keys())

        for i in range(len(gene_caller_ids_to_profile)):
            gene_callers_id = gene_caller_ids_to_profile[i]
            codons_to_profile = codons_in_genes_to_profile_SCVs_dict[gene_callers_id]

            gene_call = self.genes_in_contigs_dict[gene_callers_id]
            contig_name = gene_call['contig']
            codon_frequencies_dict = codon_frequencies.process_gene_call(self.bam, gene_call, self.contig_sequences[contig_name]['sequence'], codons_to_profile)

            for codon_order in codon_frequencies_dict:
                e = codon_frequencies_dict[codon_order]

                db_entry = {'sample_id': self.sample_id, 'corresponding_gene_call': gene_callers_id}
                db_entry['reference'] = e['reference']
                db_entry['coverage'] = e['coverage']
                db_entry['departure_from_reference'] = e['departure_from_reference']
                db_entry['codon_order_in_gene'] = codon_order
                for codon in list(constants.codon_to_AA.keys()):
                    db_entry[codon] = e['frequencies'][codon]

                variable_codons_table.append(db_entry)

        variable_codons_table.store()

        # clear contents of set
        self.codons_in_genes_to_profile_SCVs.clear()

        if len(codon_frequencies.not_reported_items):
            items = codon_frequencies.not_reported_items
            self.run.warning("The profiler of single-codon variants failed to report anything for a\
                              total of %d items, because they looked weird to anvi'o :( Here is a list\
                              of those that did ended up being ignored: '%s'." % (len(items), ', '.join(items)))


    def generate_variabile_nts_table(self):
        if self.skip_SNV_profiling:
            return

        variable_nts_table = TableForVariability(self.profile_db_path, progress=self.progress)

        for contig in self.contigs:
            for split in contig.splits:
                for column_profile in list(split.column_profiles.values()):
                    # let's figure out more about this particular variable position
                    pos_in_contig = column_profile['pos_in_contig']

                    column_profile['in_partial_gene_call'], \
                    column_profile['in_complete_gene_call'],\
                    column_profile['base_pos_in_codon'] = self.get_nt_position_info(contig.name, pos_in_contig)

                    column_profile['sample_id'] = self.sample_id
                    column_profile['corresponding_gene_call'] = -1 # this means there is no gene call that corresponds to this
                                                                   # nt position, which will be updated in the following lines.
                                                                   # yeah, we use '-1', because genecaller ids start from 0 :/
                    column_profile['codon_order_in_gene'] = -1

                    # if this particular position (`pos_in_contig`) falls within a COMPLETE gene call,
                    # we would like to find out which unique gene caller id(s) match to this position.
                    if column_profile['in_complete_gene_call']:
                        corresponding_gene_caller_ids = self.get_corresponding_gene_caller_ids_for_base_position(contig.name, pos_in_contig)

                        # if there are more than one corresponding gene call, this usually indicates an assembly error
                        # just to be on the safe side, we will not report a corresopnding unique gene callers id for this
                        # position
                        if len(corresponding_gene_caller_ids) == 1:
                            # if we are here, it means this nucleotide position is in a complete gene call. we will do two things here.
                            # first, we will store the gene_callers_id that corresponds to this nt position, and then we will store the
                            # order of the corresponding codon in the gene for this nt position.
                            gene_callers_id = corresponding_gene_caller_ids[0]
                            column_profile['corresponding_gene_call'] = gene_callers_id
                            column_profile['codon_order_in_gene'] = self.get_corresponding_codon_order_in_gene(gene_callers_id, contig.name, pos_in_contig)

                            # save this information for later use
                            self.codons_in_genes_to_profile_SCVs.add((gene_callers_id, column_profile['codon_order_in_gene']),)

                    variable_nts_table.append(column_profile)

        variable_nts_table.store()

        self.layer_additional_data['num_SNVs_reported'] = variable_nts_table.num_entries
        self.layer_additional_keys.append('num_SNVs_reported')


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
                self.sample_id = os.path.basename(self.input_file_path).upper().split('.BAM')[0]
                self.sample_id = self.sample_id.replace('-', '_')
                self.sample_id = self.sample_id.replace('.', '_')
                if self.sample_id[0] in constants.digits:
                    self.sample_id = 's' + self.sample_id
                utils.check_sample_id(self.sample_id)
            if self.serialized_profile_path:
                self.serialized_profile_path = os.path.abspath(self.serialized_profile_path)
                self.sample_id = os.path.basename(os.path.dirname(self.serialized_profile_path))


    def check_contigs_without_any_gene_calls(self, contig_names):
        if not self.a_meta['genes_are_called']:
            self.run.warning("The contigs database '%s' does not contain any gene calls. Which means the profiling step\
                              will not be able to characterize 'gene coverages'. If you are OK with this, anvi'o will be\
                              OK with it as well." % (self.contigs_db_path))
            return

        contig_names = set(contig_names)
        contigs_without_any_gene_calls = [c for c in contig_names if c not in self.contig_name_to_genes]

        if len(contigs_without_any_gene_calls):
            import random
            P = lambda x: 'are %d contigs' % (x) if x > 1 else 'there is one contig'
            self.run.warning('According to the data generated in the contigs database, there %s in your BAM file\
                              with 0 gene calls. Which may not be unusual if (a) some of your contigs are very short,\
                              or (b) your the gene caller was not capable of dealing with the type of data you had.\
                              If you would like to take a look yourself, here is one contig that is missing any genes: %s"' %\
                                      (P(len(contigs_without_any_gene_calls)), random.choice(contigs_without_any_gene_calls)))


    def list_contigs(self):
        import signal
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = pysam.Samfile(self.input_file_path, 'rb')
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
            filesnpaths.shutil.rmtree(self.output_directory)
            raise ConfigError("Anvi'o applied your min/max lenght criteria for contigs to filter out the bad ones\
                               and has bad news: not a single contig in your contigs database was greater than %s\
                               and smaller than %s nts :( While at it, anvi'o removed your half-baked output directory,\
                               too." % (pp(self.min_contig_length), pp(self.max_contig_length)))
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
        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = bamops.BAMFileObject(self.input_file_path, run=self.run, progress=self.progress).get()
        self.num_reads_mapped = self.bam.mapped
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lengths = self.bam.lengths

        utils.check_contig_names(self.contig_names)

        self.run.info('input_bam', self.input_file_path)
        self.run.info('output_dir', self.output_directory, display_only=True)
        self.run.info('total_reads_mapped', pp(int(self.num_reads_mapped)))
        self.run.info('num_contigs', pp(len(self.contig_names)))

        if self.contig_names_of_interest:
            indexes = [self.contig_names.index(r) for r in self.contig_names_of_interest if r in self.contig_names]
            self.contig_names = [self.contig_names[i] for i in indexes]
            self.contig_lengths = [self.contig_lengths[i] for i in indexes]
            self.run.info('num_contigs_selected_for_analysis', pp(len(self.contig_names)))

        # it brings good karma to let the user know what the hell is wrong with their data:
        self.check_contigs_without_any_gene_calls(self.contig_names)

        # check for the -M parameter.
        self.remove_contigs_based_on_min_max_contig_length()

        # let's see whether the user screwed up to follow the simple instructions
        # mentioned here: http://merenlab.org/2015/05/01/anvio-tutorial/#preparation
        for contig_name in self.contig_names:
            if contig_name not in self.contig_names_in_contigs_db:
                raise ConfigError("At least one contig name in your BAM file does not match contig names stored in the\
                                    contigs database. For instance, this is one contig name found in your BAM file: '%s',\
                                    and this is another one found in your contigs database: '%s'. You may be using an\
                                    contigs database for profiling that has nothing to do with the BAM file you are\
                                    trying to profile, or you may have failed to fix your contig names in your FASTA file\
                                    prior to mapping, which is described here: %s"\
                                        % (contig_name, self.contig_names_in_contigs_db.pop(), 'http://goo.gl/Q9ChpS'))

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
        self.run.info('total_reads_mapped', pp(int(self.num_reads_mapped)))
        self.run.info('num_contigs', pp(self.num_contigs))

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

    @staticmethod
    def profile_contig_worker(available_index_queue, output_queue, info_dict):
        bam_file = pysam.Samfile(info_dict['input_file_path'], 'rb')
        while True:
            index = available_index_queue.get(True)
            contig_name = info_dict['contig_names'][index]
            contig = contigops.Contig(contig_name)
            contig.length = info_dict['contig_lengths'][index]
            contig.split_length = info_dict['split_length']
            contig.min_coverage_for_variability = info_dict['min_coverage_for_variability']
            contig.skip_SNV_profiling = info_dict['skip_SNV_profiling']
            contig.report_variability_full = info_dict['report_variability_full']

            # populate contig with empty split objects and
            for split_name in info_dict['contig_name_to_splits'][contig_name]:
                s = info_dict['splits_basic_info'][split_name]
                split_sequence = info_dict['contig_sequences'][contig_name]['sequence'][s['start']:s['end']]
                split = contigops.Split(split_name, split_sequence, contig_name, s['order_in_parent'], s['start'], s['end'])
                contig.splits.append(split)

            # analyze coverage for each split
            contig.analyze_coverage(bam_file)

            # test the mean coverage of the contig.
            if contig.coverage.mean < info_dict['min_mean_coverage']:
                 output_queue.put(None)
                 continue

            if not info_dict['skip_SNV_profiling']:
                contig.analyze_auxiliary(bam_file)

            output_queue.put(contig)

            for split in contig.splits:
                del split.coverage
                del split.auxiliary
                del split
            del contig.splits[:]
            del contig.coverage
            del contig

        # we are closing this object here for clarity, although w
        # are not really closing it since the code never reaches here
        # and the worker is killed by its parent:
        bam_file.close()
        return


    def profile(self):
        manager = multiprocessing.Manager()
        info_dict = manager.dict()
        info_dict = {
            'input_file_path': self.input_file_path,
            'contig_names': self.contig_names,
            'contig_lengths': self.contig_lengths,
            'splits_basic_info': self.splits_basic_info,
            'split_length': self.a_meta['split_length'],
            'min_coverage_for_variability': self.min_coverage_for_variability,
            'skip_SNV_profiling': self.skip_SNV_profiling,
            'report_variability_full': self.report_variability_full,
            'contig_name_to_splits': self.contig_name_to_splits,
            'contig_sequences': self.contig_sequences,
            'min_mean_coverage': self.min_mean_coverage
        }

        available_index_queue = manager.Queue()
        output_queue = manager.Queue(self.queue_size)

        # put contig indices into the queue to be read from within
        # the worker
        for i in range(0, self.num_contigs):
            available_index_queue.put(i)

        processes = []
        for i in range(0, self.num_threads):
            processes.append(multiprocessing.Process(target=BAMProfiler.profile_contig_worker, args=(available_index_queue, output_queue, info_dict)))

        for proc in processes:
            proc.start()

        recieved_contigs = 0
        discarded_contigs = 0
        memory_usage = None

        self.progress.new('Profiling using ' + str(self.num_threads) + ' thread%s' % ('s' if self.num_threads > 1 else ''))
        self.progress.update('initializing threads ...')
        # FIXME: memory usage should be generalized.
        last_memory_update = int(time.time())

        self.progress.update('contigs are being processed ...')
        while recieved_contigs < self.num_contigs:
            try:
                contig = output_queue.get()

                # if we have a contig back, it means we are good to go with it,
                # otherwise it is garbage.
                if contig:
                    self.contigs.append(contig)
                else:
                    discarded_contigs += 1

                recieved_contigs += 1

                if (int(time.time()) - last_memory_update) > 5:
                    memory_usage = utils.get_total_memory_usage()
                    last_memory_update = int(time.time())

                self.progress.update('Processed %d of %d contigs. Current memory usage: %s' % \
                            (recieved_contigs, self.num_contigs, memory_usage or '...'))

                # here you're about to witness the poor side of Python (or our use of it).
                # the problem we run into here was the lack of action from the garbage
                # collector on the processed objects. although we couldn't find any refs to
                # these objects, garbage collecter kept them in the memory, and `del` statement
                # on the `split` object did not yield any improvement either. so here we are
                # accessing to the atomic data structures in our split objects to try to relieve
                # the memory by encouraging the garbage collector to realize what's up
                # explicitly.
                if self.write_buffer_size > 0 and len(self.contigs) % self.write_buffer_size == 0:
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
            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for proc in processes:
            proc.terminate()

        self.store_contigs_buffer()
        self.auxiliary_db.close()
        self.progress.end()

        # FIXME: this needs to be checked:
        if discarded_contigs > 0:
            self.run.info('contigs_after_C', pp(recieved_contigs - discarded_contigs))

        overall_mean_coverage = 1
        if self.total_length_of_all_contigs != 0:
            overall_mean_coverage = self.total_coverage_values_for_all_contigs / self.total_length_of_all_contigs

        # FIXME: We know this is ugly. You can keep your opinion to yourself.
        if overall_mean_coverage > 0.0:
            # avoid dividing by zero
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE atomic_data_splits SET abundance = abundance / " + str(overall_mean_coverage) + " * 1.0;")
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE atomic_data_contigs SET abundance = abundance / " + str(overall_mean_coverage) + " * 1.0;")

        self.check_contigs(num_contigs=recieved_contigs-discarded_contigs)


    def store_contigs_buffer(self):
        for contig in self.contigs:
            self.total_length_of_all_contigs += contig.length
            self.total_coverage_values_for_all_contigs += contig.coverage.mean * contig.length

            # we will divide every abundance after profiling is done.
            contig.abundance = contig.coverage.mean
            for split in contig.splits:
                split.abundance = split.coverage.mean

        self.progress.verbose = False
        self.generate_variabile_nts_table()
        self.generate_variabile_codons_table()
        self.store_split_coverages()
        self.progress.verbose = True

        # creating views in the database for atomic data we gathered during the profiling. Meren, please note
        # that the first entry has a view_id, and the second one does not have one. I know you will look at this
        # and be utterly confused 2 months from now. Please go read the description given in the dbops.py for the
        # function create_new_view defined in the class TablesForViews.
        view_data_splits, view_data_contigs = contigops.get_atomic_data_dicts(self.sample_id, self.contigs)

        TablesForViews(self.profile_db_path).create_new_view(
                                        data_dict=view_data_splits,
                                        table_name='atomic_data_splits',
                                        table_structure=t.atomic_data_table_structure,
                                        table_types=t.atomic_data_table_types,
                                        view_name='single',
                                        append_mode=True)

        TablesForViews(self.profile_db_path).create_new_view(
                                        data_dict=view_data_contigs,
                                        table_name='atomic_data_contigs',
                                        table_structure=t.atomic_data_table_structure,
                                        table_types=t.atomic_data_table_types,
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
            self.run.warning("You are about to generate a blank profile. This is what we do when we have nothing\
                              but a contigs database to play with. Because anvi'o is lazy, it will not check the\
                              rest of the parameters you may have declred. Most of them will not matter.")

            if not self.output_directory:
                raise ConfigError("If you want to generate a blank profile, you need to declare an output diretory path.")
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
        if not self.min_mean_coverage >= 0:
            raise ConfigError("Minimum mean coverage must be 0 or larger.")
        if not self.min_contig_length >= 0:
            raise ConfigError("Minimum contig length must be 0 or larger.")
        if not self.max_contig_length >= 100:
            raise ConfigError("Maximum contig length can't be less than 100 base pairs.")
        if self.min_contig_length >= self.max_contig_length:
            raise ConfigError("Maximum contig length (%s) must be larger than the minimum\
                               contig length (%s). Seriously though." % (pp(self.max_contig_length), pp(self.min_contig_length)))

        if self.num_threads < 1:
            raise ConfigError("Nice try. Obviously, number of threds can not be less than 1.")

        if not self.queue_size:
            self.queue_size = self.num_threads * 2

        if not self.write_buffer_size:
            self.run.warning("You set the write buffer size to 0. Which means, the profiling data will be kept in memory until\
                              the very end of the processing.")

        if self.write_buffer_size < 0:
            raise ConfigError('No. Write buffer size can not have a negative value.')
