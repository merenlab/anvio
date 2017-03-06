# -*- coding: utf-8
# pylint: disable=line-too-long
"""The library to merge multiple profiles.

The default client of this library is under bin/anvi-merge"""


import os

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.tables as tables
import anvio.terminal as terminal
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
__status__ = "Development"


pp = terminal.pretty_print

run = terminal.Run()
progress = terminal.Progress()


# Let's make sure we have CONCOCT available
__CONCOCT_IS_AVAILABLE__ = False
try:
    import anvio.concoct as concoct
    __CONCOCT_IS_AVAILABLE__ = True
except ImportError as e:
    run.warning("The CONCOCT module could not be imported :( Anvi'o will still be able to perform\
                 the merging, however, the unsupervised binning results will not be available to\
                 you in the resulting merged profile database. This is what the module was upset about:\
                 '''%s''', in case you would like to fix this problem." % e)


class MultipleRuns:
    def __init__(self, args, run=run, progress=progress):
        self.progress = progress
        self.run = run

        self.max_num_splits_for_hierarchical_clustering = constants.max_num_splits_for_hierarchical_clustering

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.sample_id = A('sample_name')
        self.contigs_db_path = A('contigs_db')
        self.input_profile_db_paths = A('input')
        self.output_directory = A('output_dir')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')
        self.skip_concoct_binning = A('skip_concoct_binning')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.description_file_path = A('description')

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        self.split_names = None
        self.sample_ids_found_in_input_dbs = []
        self.normalization_multiplier = {}

        self.profile_dbs_info_dict = {}
        self.profiles = []

        self.merged_profile_db_path = None

        self.clustering_configs = constants.clustering_configs['merged']

        self.database_paths = {'CONTIGS.db': self.contigs_db_path}

        # we don't know what we are about
        self.description = None


    def populate_profile_dbs_info_dict(self):
        improper = []

        for p in self.input_profile_db_paths:
            dbops.is_profile_db(p)

            profile_db = dbops.ProfileDatabase(p)

            if profile_db.meta['db_type'] != 'profile' or profile_db.meta['blank'] or profile_db.meta['merged']:
                improper.append(p)
            else:
                self.profile_dbs_info_dict[p] = profile_db.meta

        proper = [p for p in self.input_profile_db_paths if p not in improper]

        if len(improper) == len(proper):
            raise ConfigError("None of the databases you asked anvi'o to merge were single, non-blank anvi'o profiles. If you\
                               are not testing anvi'o and yet found yourself here, it is safe to assume that something somewhere\
                               in your workflow is quite wrong :/")

        if not len(proper) > 1:
            raise ConfigError("Anvi'o can only merge single, non-blank anvi'o profiles. You have only one database that fits into that\
                               criterion. So there is nothing really to merge here. Yes?")

        if improper:
            self.run.warning("Pleae read carefuly. You sent %d profile databases to anvi'o merger to be merged. However, not\
                              all of them were single, non-blank anvi'o profiles. Anvi'o removed %d of them, and will merge\
                              only the remaining %d. At the end of this warning you will find a list of paths to those databases\
                              anvi'o excluded from merging. If you are not happy with that, please carefully examine what went wrong.\
                              Here are all the paths for excluded databases: %s." \
                                            % (len(self.input_profile_db_paths), len(improper), len(proper), ', '.join(["'%s'" % p for p in improper])))


    def sanity_check(self):
        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=self.overwrite_output_destinations)

        if not self.contigs_db_path:
            raise ConfigError("You must provide a contigs database for this operation.")

        if not os.path.exists(self.contigs_db_path):
            raise ConfigError("Anvi'o couldn't find the contigs database where you said it would be :/")

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering\
                                while also asking it to enforce it.")

        self.populate_profile_dbs_info_dict()

        self.sample_ids_found_in_input_dbs = sorted([v['sample_id'] for v in list(self.profile_dbs_info_dict.values())])
        if len(self.profile_dbs_info_dict) != len(set(self.sample_ids_found_in_input_dbs)):
            raise ConfigError("Sample ids in each single profile database to be merged must be unique. But it is not the case\
                               with your input :/ Here are the sample names in case you would like to find out which ones occur\
                               more than once: '%s'" % (', '.join(self.sample_ids_found_in_input_dbs)))

        # test open the contigs database (and learn its hash while doing it) to make sure we don't have
        # a deal breaker just yet
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, quiet=True)
        contigs_db_hash = contigs_db.meta['contigs_db_hash']
        contigs_db.disconnect()

        for k, p in [('total_length', 'The number of nucleotides described'),
                     ('num_contigs', 'The number of contigs'),
                     ('version', 'The version number'),
                     ('num_splits', 'The number of splits'),
                     ('min_contig_length', 'The minimum contig length (-M) values'),
                     ('min_coverage_for_variability', 'The minimum coverage values to report variability (-V)'),
                     ('report_variability_full', 'Whether to report full variability (--report-variability-full) flags'),
                     ('AA_frequencies_profiled', 'Profile AA frequencies flags (--profile-AA-frequencies)'),
                     ('SNVs_profiled', 'SNV profiling flags (--skip-SNV-profiling)')]:
            v = set([r[k] for r in list(self.profile_dbs_info_dict.values())])
            if len(v) > 1:
                raise ConfigError("%s are not identical for all profiles to be merged, which is a \
                                    deal breaker. All profiles that are going to be merged must be\
                                    run with identical flags and parameters :/" % p)

        # get split names from one of the profile databases. split names must be identical across all
        self.split_names = sorted(list(dbops.get_split_names_in_profile_db(list(self.profile_dbs_info_dict.keys())[0])))

        # make sure all runs were profiled using the same contigs database (if one used):
        hashes_for_profile_dbs = set([r['contigs_db_hash'] for r in self.profile_dbs_info_dict.values()])
        if len(hashes_for_profile_dbs) != 1:
            if None in hashes_for_profile_dbs:
                raise ConfigError("It seems there is at least one run in the mix that was profiled using an\
                                          contigs database, and at least one other that was profiled without using\
                                          one. This is not good. All runs must be profiled using the same contigs\
                                          database, or all runs must be profiled without a contigs database :/")
            else:
                raise ConfigError("It seems these runs were profiled using different contigs databases (or\
                                          different versions of the same contigs database). All runs must be\
                                          profiled using the same contigs database, or all runs must be profiled\
                                          without a contigs database :/")


        # make sure the hash for contigs db is identical across all profile databases:
        if list(hashes_for_profile_dbs)[0] != contigs_db_hash:
            raise ConfigError("The contigs database you provided, which is identified with hash '%s', does\
                                      not seem to match the run profiles you are trying to merge, which share the\
                                      hash identifier of '%s'. What's up with that?" % (contigs_db_hash, list(hashes_for_profile_dbs)[0]))

        # do we have a description file?
        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'rU').read()


    def set_sample_id(self):
        if self.sample_id:
            utils.check_sample_id(self.sample_id)
        else:
            self.sample_id = os.path.basename(self.output_directory)
            self.sample_id = self.sample_id.replace('-', '_')
            if self.sample_id[0] in constants.digits:
                self.sample_id = 's' + self.sample_id
            utils.check_sample_id(self.sample_id)


    def merge_variable_nts_tables(self):
        variable_nts_table = dbops.TableForVariability(self.merged_profile_db_path, progress=self.progress)

        for input_profile_db_path in self.profile_dbs_info_dict:
            sample_profile_db = dbops.ProfileDatabase(input_profile_db_path, quiet=True)
            sample_variable_nts_table = sample_profile_db.db.get_table_as_list_of_tuples(tables.variable_nts_table_name, tables.variable_nts_table_structure)
            sample_profile_db.disconnect()

            for tpl in sample_variable_nts_table:
                entry = tuple([variable_nts_table.next_id(tables.variable_nts_table_name)] + list(tpl[1:]))
                variable_nts_table.db_entries.append(entry)

        variable_nts_table.store()


    def merge_variable_aas_tables(self):
        variable_aas_table = dbops.TableForAAFrequencies(self.merged_profile_db_path, progress=self.progress)

        for input_profile_db_path in self.profile_dbs_info_dict:
            sample_profile_db = dbops.ProfileDatabase(input_profile_db_path, quiet=True)
            sample_variable_aas_table = sample_profile_db.db.get_table_as_list_of_tuples(tables.variable_aas_table_name, tables.variable_aas_table_structure)
            sample_profile_db.disconnect()

            for tpl in sample_variable_aas_table:
                entry = tuple([variable_aas_table.next_id(tables.variable_aas_table_name)] + list(tpl[1:]))
                variable_aas_table.db_entries.append(entry)

        variable_aas_table.store()


    def merge_split_coverage_data(self):
        output_file_path = os.path.join(self.output_directory, 'AUXILIARY-DATA.h5')
        merged_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(output_file_path, self.contigs_db_hash, create_new=True)

        # fill coverages in from all samples
        for input_profile_db_path in self.profile_dbs_info_dict:
            input_file_path = os.path.join(os.path.dirname(input_profile_db_path), 'AUXILIARY-DATA.h5')
            sample_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(input_file_path, self.contigs_db_hash)

            for split_name in self.split_names:
                coverages_dict = sample_split_coverage_values.get(split_name)
                for sample_name in coverages_dict:
                    merged_split_coverage_values.append(split_name, sample_name, coverages_dict[sample_name])

            sample_split_coverage_values.close()

        merged_split_coverage_values.close()


    def set_normalization_multiplier(self):
        # WARNING: here we set normalization values based on total number of reads mapped to each sample.
        # this is not the best way to do it. a better way probably required all reads obtained from each
        # run, yet even that would wrongly assume equal eukaryotic contamination, etc. normalization is a bitch.

        smallest_sample_size = min(self.total_reads_mapped_per_sample.values())

        if smallest_sample_size == 0:
            raise ConfigError("It seems at least one of the samples you are trying to merge has zero hits. Here is a\
                                list of all samples and number of mapped reads they have: %s." \
                                    % ', '.join(['"%s": %s' % (s, pp(self.total_reads_mapped_per_sample[s])) for s in self.total_reads_mapped_per_sample]))

        for input_profile_db_path in self.profile_dbs_info_dict:
            sample_id = self.profile_dbs_info_dict[input_profile_db_path]['sample_id']
            self.normalization_multiplier[input_profile_db_path] = smallest_sample_size * 1.0 / self.total_reads_mapped_per_sample[sample_id]

        PRETTY = lambda x: ', '.join(['%s: %.2f' % (self.profile_dbs_info_dict[s]['sample_id'], x[s]) for s in sorted(list(x.keys()))])
        self.run.warning("anvio just set the normalization values for each sample based on\
                          how many mapped reads they contained. All normalized coverages\
                          will use this information: %s" % PRETTY(self.normalization_multiplier))


    def merge(self):
        self.sanity_check()
        self.set_sample_id()

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=self.overwrite_output_destinations)

        self.run.log_file_path = os.path.join(self.output_directory, 'RUNLOG.txt')

        # set database paths
        self.merged_profile_db_path = os.path.join(self.output_directory, 'PROFILE.db')
        self.samples_db_path = os.path.join(self.output_directory, 'SAMPLES.db')

        profile_db = dbops.ProfileDatabase(self.merged_profile_db_path)

        C = lambda x: list(self.profile_dbs_info_dict.values())[0][x]
        self.contigs_db_hash = C('contigs_db_hash')
        self.min_contig_length = C('min_contig_length')
        self.num_contigs = C('num_contigs')
        self.num_splits = C('num_splits')
        self.min_coverage_for_variability = C('min_coverage_for_variability')
        self.report_variability_full = C('report_variability_full')
        self.AA_frequencies_profiled = C('AA_frequencies_profiled')
        self.SNVs_profiled = C('SNVs_profiled')
        self.total_length = C('total_length')

        if self.num_splits > self.max_num_splits_for_hierarchical_clustering and not self.enforce_hierarchical_clustering:
            self.run.warning("It seems you have more than %s splits in your samples to be merged. This is the\
                              soft limit for anvi'o to attempt to create a hierarchical clustering of your splits\
                              (which becomes the center tree in all anvi'o displays). If you want a hierarchical\
                              clustering to be done anyway, please see the flag `--enforce-hierarchical-clustering`.\
                              But more importantly, please take a look at the anvi'o tutorial to make sure you know\
                              your better options to analyze large metagenomic datasets with anvi'o." \
                                                                % pp(self.max_num_splits_for_hierarchical_clustering))
            self.skip_hierarchical_clustering = True

        if self.num_splits > self.max_num_splits_for_hierarchical_clustering and self.enforce_hierarchical_clustering:
            self.run.warning("Becasue you have used the flag `--enforce-hierarchical-clustering`, anvi'o will attempt\
                              to create a hierarchical clustering of your %s splits. It may take a bit of time..." \
                                                                % pp(self.max_num_splits_for_hierarchical_clustering))

        self.total_reads_mapped_per_sample =  dict([(v['sample_id'], int(v['total_reads_mapped'][v['sample_id']])) for v in list(self.profile_dbs_info_dict.values())])

        sample_ids_list = ', '.join(sorted(self.sample_ids_found_in_input_dbs)) 
        total_reads_mapped_list = ', '.join([str(self.total_reads_mapped_per_sample[sample_id]) for sample_id in self.sample_ids_found_in_input_dbs])

        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': sample_ids_list,
                       'total_reads_mapped': total_reads_mapped_list,
                       'merged': True,
                       'blank': False,
                       'contigs_clustered': not self.skip_hierarchical_clustering,
                       'default_view': 'mean_coverage',
                       'min_contig_length': self.min_contig_length,
                       'SNVs_profiled': self.SNVs_profiled,
                       'AA_frequencies_profiled': self.AA_frequencies_profiled,
                       'num_contigs': self.num_contigs,
                       'num_splits': self.num_splits,
                       'total_length': self.total_length,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'report_variability_full': self.report_variability_full,
                       'contigs_db_hash': self.contigs_db_hash,
                       'description': self.description if self.description else '_No description is provided_'}
        profile_db.create(meta_values)

        # get view data information for both contigs and splits:
        self.atomic_data_fields, self.atomic_data_for_each_run = self.read_atomic_data_tables()

        self.split_parents = self.get_split_parents()

        self.run.info('profiler_version', anvio.__profile__version__)
        self.run.info('output_dir', self.output_directory)
        self.run.info('sample_id', self.sample_id)
        self.run.info('description', 'Found (%d characters)' % len(self.description) if self.description else None)
        self.run.info('profile_db', self.merged_profile_db_path)
        self.run.info('merged', True)
        self.run.info('contigs_db_hash', self.contigs_db_hash)
        self.run.info('num_runs_processed', len(self.sample_ids_found_in_input_dbs))
        self.run.info('merged_sample_ids', sample_ids_list)
        self.run.info('total_reads_mapped', total_reads_mapped_list)
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('clustering_performed', not self.skip_hierarchical_clustering)

        self.set_normalization_multiplier()

        self.progress.new('Merging split coverage values')
        self.merge_split_coverage_data()
        self.progress.end()

        if self.SNVs_profiled:
            self.progress.new('Merging variable positions tables')
            self.merge_variable_nts_tables()
            self.progress.end()
        else:
            self.run.warning("SNVs were not profiled, variable nt positions tables will be empty in the merged profile database.")

        if self.AA_frequencies_profiled:
            self.progress.new('Merging variable AAs tables')
            self.merge_variable_aas_tables()
            self.progress.end()
        else:
            self.run.warning("AA frequencies were not profiled, these tables will be empty in the merged profile database.")

        # critical part:
        self.gen_view_data_tables_from_atomic_data()

        # We cluster? Note: the check is being done in the function!
        self.cluster_contigs_anvio()

        self.progress.end()

        # run CONCOCT, if otherwise is not requested:
        if not self.skip_concoct_binning and __CONCOCT_IS_AVAILABLE__:
            self.bin_contigs_concoct()

        # generate a samples database while you are at it
        self.gen_samples_db_from_view_data()

        self.run.quit()


    def get_normalized_coverage_of_split(self, target, input_profile_db_path, split_name):
        return self.atomic_data_for_each_run[target][input_profile_db_path][split_name]['mean_coverage_Q2Q3'] * self.normalization_multiplier[input_profile_db_path]


    def get_max_normalized_ratio_of_split(self, target, split_name):
        denominator = float(max(self.normalized_coverages[target][split_name].values()))

        d = {}
        for input_profile_db_path in self.profile_dbs_info_dict:
            d[input_profile_db_path] = (self.normalized_coverages[target][split_name][input_profile_db_path] / denominator) if denominator else 0

        return d


    def get_relative_abundance_of_split(self, target, sample_id, split_name):
        denominator = float(sum(self.normalized_coverages[target][split_name].values()))
        return self.normalized_coverages[target][split_name][sample_id] / denominator if denominator else 0


    def gen_samples_db_from_view_data(self):
        """Geenrate a samples db for the merged profile.

           We use the ProfileSuperclass to load all the views we added into the meged profile,
           and generate clusterings of samples for each view to generate a default samples database."""

        self.run.info_single("SAMPLES.db stuff...", nl_before=1, nl_after=1, mc="blue")

        essential_fields = [f for f in self.atomic_data_fields if constants.IS_ESSENTIAL_FIELD(f)]

        class Args: pass
        args = Args()
        args.profile_db = self.merged_profile_db_path

        # initialize views.
        profile_db_super = dbops.ProfileSuperclass(args)
        profile_db_super.load_views(omit_parent_column=True)

        sample_orders = {}
        failed_attempts = []
        self.progress.new('Generating a default SAMPLES.db')
        for essential_field in essential_fields:
            self.progress.update('recovering samples order for "%s"' % (essential_field))
            try:
                sample_orders[essential_field] = \
                        clustering.get_newick_tree_data_for_dict(profile_db_super.views[essential_field]['dict'],
                                                                 distance=self.distance,
                                                                 linkage=self.linkage,
                                                                 transpose=True)
            except:
                failed_attempts.append(essential_field)
        self.progress.end()

        if not len(sample_orders):
            self.run.warning("This may or may not be important: anvi'o attempted to generate a samples\
                              database for this merged profile, however, all attempts to cluster samples\
                              based on view data available in the merged profile failed. No samples db\
                              for you :/")
            return

        if len(failed_attempts):
            self.run.warning("While anvi'o was trying to generate clusterings of samples based on view data\
                              available in the merged profile, clustering of some of the essential data\
                              failed. It is likely not a very big deal, but you shall be the judge of it.\
                              Anvi'o now proceeds to generate a samples db with clusterings it generated\
                              using the view data that worked. Here is the list of stuff that failed: '%s'"\
                              % (', '.join(failed_attempts)))

        samples_order_file_path = filesnpaths.get_temp_file_path()
        samples_order_file = open(samples_order_file_path, 'w')
        samples_order_file.write('attributes\tbasic\tnewick\n')
        for sample_order in sample_orders:
            samples_order_file.write('%s\t%s\t%s\n' % (sample_order, '', sample_orders[sample_order]))
        samples_order_file.close()

        samples_db = dbops.SamplesInformationDatabase(self.samples_db_path, quiet=False)
        samples_db.create(samples_order_path=samples_order_file_path)

        os.remove(samples_order_file_path)

        self.run.info('Samples database', self.samples_db_path)


    def gen_view_data_tables_from_atomic_data(self):
        essential_fields = [f for f in self.atomic_data_fields if constants.IS_ESSENTIAL_FIELD(f)]
        auxiliary_fields = [f for f in self.atomic_data_fields if constants.IS_AUXILIARY_FIELD(f)]

        # setting standard view table structure and types
        view_table_structure = ['contig'] + self.sample_ids_found_in_input_dbs + auxiliary_fields
        view_table_types = ['text'] + ['numeric'] * len(self.sample_ids_found_in_input_dbs) + ['text']

        # generate a dictionary for normalized coverage of each contig across samples per target
        self.normalized_coverages = {'contigs': {}, 'splits': {}}
        for target in ['contigs', 'splits']:
            for split_name in self.split_names:
                self.normalized_coverages[target][split_name] = {}
                for input_profile_db_path in self.profile_dbs_info_dict:
                    self.normalized_coverages[target][split_name][input_profile_db_path] = self.get_normalized_coverage_of_split(target, input_profile_db_path, split_name)

        # generate a dictionary for max normalized ratio of each contig across samples per target
        self.max_normalized_ratios = {'contigs': {}, 'splits': {}}
        for target in ['contigs', 'splits']:
            for split_name in self.split_names:
                self.max_normalized_ratios[target][split_name] = self.get_max_normalized_ratio_of_split(target, split_name)

        self.progress.new('Generating view data tables')
        for target in ['contigs', 'splits']:
            for essential_field in essential_fields:
                self.progress.update('Processing %s for %s ...' % (essential_field, target))

                data_dict = {}
                for split_name in self.split_names:
                    data_dict[split_name] = {'__parent__': self.split_parents[split_name]}

                    for input_profile_db_path in self.profile_dbs_info_dict:
                        sample_id = self.profile_dbs_info_dict[input_profile_db_path]['sample_id']
                        if essential_field == 'normalized_coverage':
                            data_dict[split_name][sample_id] = self.normalized_coverages[target][split_name][input_profile_db_path]
                        elif essential_field == 'max_normalized_ratio':
                            data_dict[split_name][sample_id] = self.max_normalized_ratios[target][split_name][input_profile_db_path]
                        elif essential_field == 'relative_abundance':
                            data_dict[split_name][sample_id] = self.get_relative_abundance_of_split(target, input_profile_db_path, split_name)
                        else:
                            data_dict[split_name][sample_id] = self.atomic_data_for_each_run[target][input_profile_db_path][split_name][essential_field]

                # time to store the data for this view in the profile database
                table_name = '_'.join([essential_field, target])
                dbops.TablesForViews(self.merged_profile_db_path).create_new_view(
                                                data_dict=data_dict,
                                                table_name=table_name,
                                                table_structure=view_table_structure,
                                                table_types=view_table_types,
                                                view_name=essential_field if target == 'splits' else None)

        # if SNVs were not profiled, remove all entries from variability tables:
        if not self.SNVs_profiled:
            dbops.TablesForViews(self.merged_profile_db_path).remove(view_name='variability', table_names_to_blank=['variability_splits', 'variability_contigs'])

        self.progress.end()


    def bin_contigs_concoct(self):
        self.run.info_single("Automatic binning with CONCOCT...", nl_before=1, nl_after=1, mc="blue")

        if not __CONCOCT_IS_AVAILABLE__:
            self.run.warning('The CONCOCT module is not available. Skipping the unsupervised binning\
                              with CONCOCT.')
            return

        class Args:
            pass

        args = Args()
        args.profile_db = self.merged_profile_db_path
        args.contigs_db = self.contigs_db_path
        args.debug = self.debug

        c = concoct.CONCOCT(args)
        c.cluster()
        c.store_clusters_in_db()


    def cluster_contigs_anvio(self):
        # clustering of contigs is done for each configuration file under static/clusterconfigs/merged directory;
        # at this point we don't care what those recipes really require because we already merged and generated
        # every data file that may be required.

        self.run.info_single("Anvi'o hierarchical clustering of contigs...", nl_before=1, nl_after=1, mc="blue")

        if not self.skip_hierarchical_clustering:
            for config_name in self.clustering_configs:
                config_path = self.clustering_configs[config_name]

                config = ClusteringConfiguration(config_path, self.output_directory, db_paths=self.database_paths, row_ids_of_interest=self.split_names)

                try:
                    clustering_id, newick = clustering.order_contigs_simple(config, distance=self.distance, linkage=self.linkage, progress=self.progress)
                except Exception as e:
                    self.run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
                    self.progress.end()
                    continue

                _, distance, linkage = clustering_id.split(':')

                dbops.add_hierarchical_clustering_to_db(self.merged_profile_db_path, config_name, newick, distance=distance, linkage=linkage, make_default=config_name == constants.merged_default, run=self.run)


    def get_split_parents(self):
        parents = {}
        m = list(self.atomic_data_for_each_run['splits'].values())[0]
        for split in m:
            parents[split] = m[split]["__parent__"]
        return parents


    def read_atomic_data_tables(self):
        """reads atomic data for contigs and splits from the database into a dict"""
        atomic_data_table_for_each_run = {}

        for target in ['contigs', 'splits']:
            atomic_data_table_for_each_run[target] = {}

            target_table = 'atomic_data_%s' % target

            for input_profile_db_path in self.profile_dbs_info_dict:
                db = anvio.db.DB(input_profile_db_path, dbops.get_required_version_for_db(input_profile_db_path))
                atomic_data_table_for_each_run[target][input_profile_db_path] = db.get_table_as_dict(target_table)

        atomic_data_table_fields = db.get_table_structure('atomic_data_splits')
        db.disconnect()

        return atomic_data_table_fields, atomic_data_table_for_each_run
