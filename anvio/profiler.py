# coding: utf-8
# pylint: disable=line-too-long
"""Provides the necessary class to profile BAM files."""

import os
import sys
import pysam
import shutil

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.dictio as dictio
import anvio.bamops as bamops
import anvio.terminal as terminal
import anvio.contigops as contigops
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.clusteringconfuguration import ClusteringConfiguration


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
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
        self.min_mean_coverage = A('min_mean_coverage')
        self.min_coverage_for_variability = A('min_coverage_for_variability')
        self.contigs_shall_be_clustered = A('cluster_contigs')
        self.sample_id = A('sample_name')
        self.report_variability_full = A('report_variability_full')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.skip_SNV_profiling = A('skip_SNV_profiling')
        self.profile_AA_frequencies = A('profile_AA_frequencies')
        self.gen_serialized_profile = A('gen_serialized_profile')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        # whehther the profile database is a blank (without any BAM files or reads):
        self.blank = A('blank_profile')

        if self.blank:
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
            raise ConfigError, "No contigs database, no profilin'. Bye."

        # Initialize contigs db
        dbops.ContigsSuperclass.__init__(self, self.args, r=self.run, p=self.progress)
        self.init_contig_sequences()
        self.contig_names_in_contigs_db = set(self.contigs_basic_info.keys())

        self.bam = None
        self.contigs = {}

        self.database_paths = {'CONTIGS.db': self.contigs_db_path}

        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['blank' if self.blank else 'single']

        self.atomic_data = contigops.AtomicContigSplitData(self.progress)

        # following variable will be populated during the profiling, and its content will eventually
        # be stored in t.variable_nts_table_name
        self.variable_nts_table_entries = []

        # following variable will be populated while the variable positions table is computed
        self.codons_in_genes_to_profile_AA_frequencies = set([])


    def init_dirs_and_dbs(self):
        if not self.contigs_db_path:
            raise ConfigError, "You can not run profiling without a contigs database. You can create\
                                one using 'anvi-gen-contigs-database'. Not sure how? Please see the\
                                tutorial: http://merenlab.org/2015/05/02/anvio-tutorial/"

        self.output_directory = filesnpaths.check_output_directory(self.output_directory or self.input_file_path + '-ANVIO_PROFILE',\
                                                                   ok_if_exists=self.overwrite_output_destinations)

        self.progress.new('Initializing')

        self.progress.update('Creating the output directory ...')
        filesnpaths.gen_output_directory(self.output_directory, self.progress, delete_if_exists=self.overwrite_output_destinations)

        self.progress.update('Creating a new single profile database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.profile_db_path = self.generate_output_destination('PROFILE.db')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        if self.skip_SNV_profiling:
            self.profile_AA_frequencies = False

        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': self.sample_id,
                       'merged': False,
                       'blank': self.blank,
                       'contigs_clustered': self.contigs_shall_be_clustered,
                       'default_view': 'single',
                       'min_contig_length': self.min_contig_length,
                       'SNVs_profiled': not self.skip_SNV_profiling,
                       'AA_frequencies_profiled': self.profile_AA_frequencies,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'report_variability_full': self.report_variability_full,
                       'contigs_db_hash': self.a_meta['contigs_db_hash'],
                       'gene_coverages_computed': self.a_meta['genes_are_called']}
        profile_db.create(meta_values)

        self.progress.end()

        if self.skip_SNV_profiling:
            self.run.warning('Single-nucleotide variation will not be characterized for this profile.')

        if not self.profile_AA_frequencies:
            self.run.warning('Amino acid linkmer frequencies will not be characterized for this profile.')


    def _run(self):
        self.check_args()

        self.set_sample_id()

        self.init_dirs_and_dbs()

        self.run.log_file_path = self.generate_output_destination('RUNLOG.txt')
        self.run.info('anvio', anvio.__version__)
        self.run.info('profiler_version', anvio.__profile__version__)
        self.run.info('sample_id', self.sample_id)
        self.run.info('profile_db', self.profile_db_path, display_only=True)
        self.run.info('contigs_db', True if self.contigs_db_path else False)
        self.run.info('contigs_db_hash', self.a_meta['contigs_db_hash'])
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('merged', False)
        self.run.info('blank', self.blank)
        self.run.info('split_length', self.a_meta['split_length'])
        self.run.info('min_contig_length', self.min_contig_length)
        self.run.info('min_mean_coverage', self.min_mean_coverage)
        self.run.info('clustering_performed', self.contigs_shall_be_clustered)
        self.run.info('min_coverage_for_variability', self.min_coverage_for_variability)
        self.run.info('skip_SNV_profiling', self.skip_SNV_profiling)
        self.run.info('profile_AA_frequencies', self.profile_AA_frequencies)
        self.run.info('report_variability_full', self.report_variability_full)
        self.run.info('gene_coverages_computed', self.a_meta['genes_are_called'])

        # this is kinda important. we do not run full-blown profile function if we are dealing with a summarized
        # profile...
        if self.blank:
            self.init_mock_profile()
        elif self.input_file_path:
            self.init_profile_from_BAM()
            self.profile()
            if self.gen_serialized_profile:
                self.store_profile()
        elif self.serialized_profile_path:
            self.init_serialized_profile()
        else:
            raise ConfigError, "What are you doing? :( Whatever it is, anvi'o will have none of it."

        self.generate_variabile_nts_table()
        self.generate_variabile_aas_table()
        self.generate_gene_coverages_table()
        self.store_split_coverages()

        # creating views in the database for atomic data we gathered during the profiling. Meren, please note
        # that the first entry has a view_id, and the second one does not have one. I know you will look at this
        # and be utterly confused 2 months from now. Please go read the description given in the dbops.py for the
        # function create_new_view defined in the class TablesForViews.
        view_data_splits, view_data_contigs = self.atomic_data.get_data(self.sample_id, self.contigs)
        dbops.TablesForViews(self.profile_db_path).create_new_view(
                                        data_dict=view_data_splits,
                                        table_name='atomic_data_splits',
                                        table_structure=t.atomic_data_table_structure,
                                        table_types=t.atomic_data_table_types,
                                        view_name='single')

        dbops.TablesForViews(self.profile_db_path).create_new_view(
                                        data_dict=view_data_splits,
                                        table_name='atomic_data_contigs',
                                        table_structure=t.atomic_data_table_structure,
                                        table_types=t.atomic_data_table_types,
                                        view_name=None)

        # OK. if this is a blank profile, atomic_data_* tables will be completely empty. but having no
        # split names in the profile database becomes very limiting for downstream analyses, so here
        # we will generate a null atomic_data_splits table, only purpose of which will be to hold the
        # split names
        if self.blank:
            # creating a null view_data_splits dict:
            view_data_splits = dict(zip(self.split_names, [dict(zip(t.atomic_data_table_structure[1:], [None] * len(t.atomic_data_table_structure[1:])))] * len(self.split_names)))
            dbops.TablesForViews(self.profile_db_path).remove('single', table_names_to_blank=['atomic_data_splits'])
            dbops.TablesForViews(self.profile_db_path).create_new_view(
                                           data_dict=view_data_splits,
                                           table_name='atomic_data_splits',
                                           table_structure=t.atomic_data_table_structure,
                                           table_types=t.atomic_data_table_types,
                                           view_name='single')

        if self.contigs_shall_be_clustered:
            self.cluster_contigs()

        runinfo_serialized = self.generate_output_destination('RUNINFO.cp')
        self.run.info('runinfo', runinfo_serialized)
        self.run.store_info_dict(runinfo_serialized, strip_prefix=self.output_directory)

        if self.bam:
            self.bam.close()

        self.run.quit()


    def generate_variabile_aas_table(self):
        if self.skip_SNV_profiling or not self.profile_AA_frequencies:
            # there is nothing to generate really..
            self.run.info('AA_frequencies_table', False, quiet=True)
            return

        self.progress.new('Computing AA frequencies at variable positions')

        variable_aas_table = dbops.TableForAAFrequencies(self.profile_db_path, progress=self.progress)

        aa_frequencies = bamops.AAFrequencies()

        codons_in_genes_to_profile_AA_frequencies_dict = {}
        for gene_call_id, codon_order in self.codons_in_genes_to_profile_AA_frequencies:
            if gene_call_id not in codons_in_genes_to_profile_AA_frequencies_dict:
                codons_in_genes_to_profile_AA_frequencies_dict[gene_call_id] = set([])
            codons_in_genes_to_profile_AA_frequencies_dict[gene_call_id].add(codon_order)

        gene_caller_ids_to_profile = codons_in_genes_to_profile_AA_frequencies_dict.keys()
        num_gene_caller_ids_to_profile = len(gene_caller_ids_to_profile)

        for i in range(0, len(gene_caller_ids_to_profile)):
            gene_caller_id = gene_caller_ids_to_profile[i]
            codons_to_profile = codons_in_genes_to_profile_AA_frequencies_dict[gene_caller_id]

            self.progress.update("Working on gene caller id '%d' (%d of %d) w/ %d codons of interest" \
                                % (gene_caller_id, i + 1, num_gene_caller_ids_to_profile, len(codons_to_profile)))

            gene_call = self.genes_in_contigs_dict[gene_caller_id]
            contig_name = gene_call['contig']
            aa_frequencies_dict = aa_frequencies.process_gene_call(self.bam, gene_call, self.contig_sequences[contig_name]['sequence'], codons_to_profile)

            for codon_order in aa_frequencies_dict:
                e = aa_frequencies_dict[codon_order]

                db_entry = {'sample_id': self.sample_id, 'corresponding_gene_call': gene_caller_id}
                db_entry['reference'] = e['reference']
                db_entry['coverage'] = e['coverage']
                db_entry['departure_from_reference'] = e['departure_from_reference']
                db_entry['codon_order_in_gene'] = codon_order
                for aa in constants.codon_to_AA.values():
                    db_entry[aa] = e['frequencies'][aa]

                variable_aas_table.append(db_entry)

        variable_aas_table.store()
        self.progress.end()
        self.run.info('AA_frequencies_table', True, quiet=True)


    def generate_variabile_nts_table(self):
        if self.skip_SNV_profiling:
            # there is nothing to generate really..
            self.run.info('variable_nts_table', False, quiet=True)
            return

        self.progress.new('NT Variability')
        variable_nts_table = dbops.TableForVariability(self.profile_db_path, progress=self.progress)

        for contig in self.contigs.values():
            for split in contig.splits:
                for column_profile in split.column_profiles.values():
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
                            # first, we will store the gene_caller_id that corresponds to this nt position, and then we will store the
                            # order of the corresponding codon in the gene for this nt position.
                            gene_caller_id = corresponding_gene_caller_ids[0]
                            column_profile['corresponding_gene_call'] = gene_caller_id
                            column_profile['codon_order_in_gene'] = self.get_corresponding_codon_order_in_gene(gene_caller_id, contig.name, pos_in_contig)

                            # save this information for later use
                            self.codons_in_genes_to_profile_AA_frequencies.add((gene_caller_id, column_profile['codon_order_in_gene']),)

                    variable_nts_table.append(column_profile)

        variable_nts_table.store()
        self.progress.end()
        self.run.info('variable_nts_table', True, quiet=True)


    def generate_gene_coverages_table(self):
        gene_coverages_table = dbops.TableForGeneCoverages(self.profile_db_path, progress=self.progress)

        self.progress.new('Profiling genes')
        num_contigs = len(self.contigs)
        contig_names = self.contigs.keys()
        for i in range(0, num_contigs):
            contig = contig_names[i]
            self.progress.update('Processing contig %d of %d' % (i + 1, num_contigs))
            # if no open reading frames were found in a contig, it wouldn't have an entry in the contigs table,
            # therefore there wouldn't be any record of it in contig_ORFs; so we better check ourselves before
            # we wreck ourselves and the ultimately the analysis of this poor user:
            if contig in self.contig_name_to_genes:
                gene_coverages_table.analyze_contig(self.contigs[contig], self.sample_id, self.contig_name_to_genes[contig])

        gene_coverages_table.store()
        self.progress.end()
        self.run.info('gene_coverages_table', True, quiet=True)


    def store_split_coverages(self):
        output_file = self.generate_output_destination('AUXILIARY-DATA.h5')
        split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(output_file, self.a_meta['contigs_db_hash'], create_new=True)

        self.progress.new('Storing split coverages')

        contigs_counter = 1
        for contig_name in self.contigs:
            self.progress.update('working on contig %s of %s' % (pp(contigs_counter), pp(len(self.contigs))))

            for split in self.contigs[contig_name].splits:
                split_coverage_values.append(split.name, self.sample_id, split.coverage.c)

            contigs_counter += 1

        self.progress.end()

        split_coverage_values.close()

        self.run.info('split_coverage_values', 'stored in %s' % output_file, display_only=True)
        self.run.info('split_coverage_values', True, quiet=True)


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


    def init_serialized_profile(self):
        self.progress.new('Init')
        self.progress.update('Reading serialized profile')
        self.contigs = dictio.read_serialized_object(self.serialized_profile_path)
        self.progress.end()

        self.run.info('profile_loaded_from', self.serialized_profile_path)
        self.run.info('num_contigs', pp(len(self.contigs)))

        if self.contig_names_of_interest:
            contigs_to_discard = set()
            for contig in self.contigs:
                if contig not in self.contig_names_of_interest:
                    contigs_to_discard.add(contig)

            if len(contigs_to_discard):
                for contig in contigs_to_discard:
                    self.contigs.pop(contig)
            self.run.info('num_contigs_selected_for_analysis', pp(len(self.contigs)))

        self.check_contigs()

        # it brings good karma to let the user know what the hell is wrong with their data:
        self.check_contigs_without_any_gene_calls(self.contigs.keys())

        contigs_to_discard = set()
        for contig in self.contigs.values():
            if contig.length < self.min_contig_length:
                contigs_to_discard.add(contig.name)
        if len(contigs_to_discard):
            for contig in contigs_to_discard:
                self.contigs.pop(contig)
            self.run.info('contigs_raw_longer_than_M', len(self.contigs))

        self.check_contigs()


    def list_contigs(self):
        import signal
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

        if self.input_file_path:
            self.progress.new('Init')
            self.progress.update('Reading BAM File')
            self.bam = pysam.Samfile(self.input_file_path, 'rb')
            self.progress.end()

            self.contig_names = self.bam.references
            self.contig_lengths = self.bam.lengths

            utils.check_contig_names(self.contig_names)

            for tpl in sorted(zip(self.contig_lengths, self.contig_names), reverse=True):
                print '%-40s %s' % (tpl[1], pp(int(tpl[0])))

        else:
            self.progress.new('Init')
            self.progress.update('Reading serialized profile')
            self.contigs = dictio.read_serialized_object(self.serialized_profile_path)
            self.progress.end()

            self.run.info('profile_loaded_from', self.serialized_profile_path)
            self.run.info('num_contigs', pp(len(self.contigs)))

            for tpl in sorted([(int(self.contigs[contig].length), contig) for contig in self.contigs]):
                print '%-40s %s' % (tpl[1], pp(int(tpl[0])))


    def remove_contigs_that_are_shorter_than_min_contig_length(self):
        """Removes contigs that are shorter than M"""
        contigs_longer_than_M = set()
        for i in range(0, len(self.contig_names)):
            if self.contig_lengths[i] >= self.min_contig_length:
                contigs_longer_than_M.add(i)

        if not len(contigs_longer_than_M):
            raise ConfigError, "0 contigs larger than %s nts." % pp(self.min_contig_length)
        else:
            self.contig_names = [self.contig_names[i] for i in contigs_longer_than_M]
            self.contig_lengths = [self.contig_lengths[i] for i in contigs_longer_than_M]
            self.num_contigs = len(self.contig_names)    # we will store these two
            self.total_length = sum(self.contig_lengths) # into the db in a second.

        contigs_longer_than_M = set(self.contig_names) # for fast access
        self.split_names = set([])
        self.contig_name_to_splits = {}
        for split_name in sorted(self.splits_basic_info.keys()):
            parent = self.splits_basic_info[split_name]['parent']

            if parent not in contigs_longer_than_M:
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
        self.remove_contigs_that_are_shorter_than_min_contig_length()

        # let's see whether the user screwed up to follow the simple instructions
        # mentioned here: http://merenlab.org/2015/05/01/anvio-tutorial/#preparation
        for contig_name in self.contig_names:
            if contig_name not in self.contig_names_in_contigs_db:
                raise ConfigError, "At least one contig name in your BAM file does not match contig names stored in the\
                                    contigs database. For instance, this is one contig name found in your BAM file: '%s',\
                                    and this is another one found in your contigs database: '%s'. You may be using an\
                                    contigs database for profiling that has nothing to do with the BAM file you are\
                                    trying to profile, or you may have failed to fix your contig names in your FASTA file\
                                    prior to mapping, which is described here: %s"\
                                        % (contig_name, self.contig_names_in_contigs_db.pop(), 'http://goo.gl/Q9ChpS')

        self.run.info('num_contigs_after_M', self.num_contigs, display_only=True)
        self.run.info('num_contigs', self.num_contigs, quiet=True)
        self.run.info('num_splits', self.num_splits)
        self.run.info('total_length', self.total_length)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('num_splits', self.num_splits)
        profile_db.db.set_meta_value('num_contigs', self.num_contigs)
        profile_db.db.set_meta_value('total_length', self.total_length)
        profile_db.db.set_meta_value('total_reads_mapped', int(self.num_reads_mapped))
        profile_db.disconnect()


    def init_mock_profile(self):
        self.progress.new('Init')
        self.progress.update('...')
        self.num_reads_mapped = 0
        self.progress.end()

        self.contig_names = self.contigs_basic_info.keys()
        self.contig_lengths = [self.contigs_basic_info[contig_name]['length'] for contig_name in self.contigs_basic_info]
        self.total_length = sum(self.contig_lengths)
        self.num_contigs = len(self.contig_names)

        utils.check_contig_names(self.contig_names)

        self.run.info('input_bam', None)
        self.run.info('output_dir', self.output_directory, display_only=True)
        self.run.info('total_reads_mapped', pp(int(self.num_reads_mapped)))
        self.run.info('num_contigs', pp(self.num_contigs))

        # check for the -M parameter.
        self.remove_contigs_that_are_shorter_than_min_contig_length()

        self.run.info('num_contigs_after_M', self.num_contigs, display_only=True)
        self.run.info('num_contigs', self.num_contigs, quiet=True)
        self.run.info('num_splits', self.num_splits)
        self.run.info('total_length', self.total_length)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('num_splits', self.num_splits)
        profile_db.db.set_meta_value('num_contigs', self.num_contigs)
        profile_db.db.set_meta_value('total_length', self.total_length)
        profile_db.db.set_meta_value('total_reads_mapped', int(self.num_reads_mapped))
        profile_db.disconnect()


    def generate_output_destination(self, postfix, directory=False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def profile(self):
        """Big deal function"""

        # So we start with essential stats. In the section below, we will simply go through each contig
        # in the BAM file and populate the contigs dictionary for the first time.
        for i in range(0, len(self.contig_names)):

            contig_name = self.contig_names[i]

            contig = contigops.Contig(contig_name)
            contig.length = self.contig_lengths[i]
            contig.split_length = self.a_meta['split_length']
            contig.min_coverage_for_variability = self.min_coverage_for_variability
            contig.skip_SNV_profiling = self.skip_SNV_profiling
            contig.report_variability_full = self.report_variability_full

            self.progress.new('Profiling "%s" (%d of %d) (%s nts)' % (contig.name,
                                                                      i + 1,
                                                                      len(self.contig_names),
                                                                      pp(int(contig.length))))

            # populate contig with empty split objects and
            for split_name in self.contig_name_to_splits[contig_name]:
                s = self.splits_basic_info[split_name]
                split_sequence = self.contig_sequences[contig_name]['sequence'][s['start']:s['end']]
                split = contigops.Split(split_name, split_sequence, contig_name, s['order_in_parent'], s['start'], s['end'])
                contig.splits.append(split)

            # analyze coverage for each split
            contig.analyze_coverage(self.bam, self.progress)

            # test the mean coverage of the contig.
            discarded_contigs_due_to_C = set([])
            if contig.coverage.mean < self.min_mean_coverage:
                # discard this contig and continue
                discarded_contigs_due_to_C.add(contig.name)
                self.progress.end()
                continue

            if not self.skip_SNV_profiling:
                contig.analyze_auxiliary(self.bam, self.progress)

            self.progress.end()

            # add contig to the dict.
            self.contigs[contig_name] = contig


        if discarded_contigs_due_to_C:
            self.run.info('contigs_after_C', pp(len(self.contigs)))

        # set contig abundance
        contigops.set_contigs_abundance(self.contigs)

        self.check_contigs()


    def store_profile(self):
        output_file = self.generate_output_destination('PROFILE.cp')
        self.progress.new('Storing Profile')
        self.progress.update('Serializing information for %s contigs ...' % pp(len(self.contigs)))
        dictio.write_serialized_object(self.contigs, output_file)
        self.progress.end()
        self.run.info('profile_dict', output_file)


    def check_contigs(self):
        if not len(self.contigs):
            raise ConfigError, "0 contigs to work with. Bye."


    def cluster_contigs(self):
        default_clustering = constants.blank_default if self.blank else constants.single_default

        for config_name in self.clustering_configs:
            config_path = self.clustering_configs[config_name]

            config = ClusteringConfiguration(config_path, self.output_directory, db_paths=self.database_paths, row_ids_of_interest=self.split_names)

            try:
                clustering_id, newick = clustering.order_contigs_simple(config, distance=self.distance, linkage=self.linkage, progress=self.progress)
            except Exception as e:
                self.run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
                self.progress.end()
                continue

            dbops.add_hierarchical_clustering_to_db(self.profile_db_path, config_name, newick, distance=self.distance, linkage=self.linkage, \
                                                    make_default=config_name == default_clustering, run=self.run)


    def check_args(self):
        if self.blank:
            self.run.warning("You are about to generate a blank profile. This is what we do when we have nothing\
                              but a contigs database to play with. Because anvi'o is lazy, it will not check the\
                              rest of the parameters you may have declred. Most of them will not matter.")

            if not self.output_directory:
                raise ConfigError, "If you want to generate a blank profile, you need to declare an output diretory path."
            if not self.sample_id:
                raise ConfigError, "Mock profiles require a sample name to be declared. Because :/"
            return

        if (not self.input_file_path) and (not self.serialized_profile_path):
            raise ConfigError, "You didn't declare any input files :/ If you intend to create a blank profile without any,\
                                input file, you should be a bit more explicit about your intention (you know, in the help\
                                there is a flag for it and all). Otherwise you should either provide an input BAM file, or\
                                a serialized anvi'o profile. See '--help' maybe?"
        if self.input_file_path and self.serialized_profile_path:
            raise ConfigError, "You can't declare both an input file and a serialized profile."
        if self.serialized_profile_path and (not self.output_directory):
            raise ConfigError, "When loading serialized profiles, you need to declare an output directory."
        if self.input_file_path and not os.path.exists(self.input_file_path):
            raise ConfigError, "No such file: '%s'" % self.input_file_path
        if self.serialized_profile_path and not os.path.exists(self.serialized_profile_path):
            raise ConfigError, "No such file: '%s'" % self.serialized_profile_path
        if not self.min_coverage_for_variability >= 0:
            raise ConfigError, "Minimum coverage for variability must be 0 or larger."
        if not self.min_mean_coverage >= 0:
            raise ConfigError, "Minimum mean coverage must be 0 or larger."
        if not self.min_contig_length >= 0:
            raise ConfigError, "Minimum contig length must be 0 or larger."
