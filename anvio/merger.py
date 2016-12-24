# -*- coding: utf-8
# pylint: disable=line-too-long
"""The library to merge multiple profiles.

The default client of this library is under bin/anvi-merge"""


import os

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.tables as tables
import anvio.dictio as dictio
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
        self.input_runinfo_paths = A('input')
        self.output_directory = A('output_dir')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')
        self.skip_concoct_binning = A('skip_concoct_binning')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        self.split_names = None
        self.merged_sample_ids = []
        self.input_runinfo_dicts = {}
        self.normalization_multiplier = {}
        self.profiles = []

        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['merged']

        self.database_paths = {'CONTIGS.db': self.contigs_db_path}


    def read_runinfo_dict(self, path):
        runinfo = dictio.read_serialized_object(path)
        sample_id = runinfo['sample_id']

        if not sample_id in self.merged_sample_ids:
            self.merged_sample_ids.append(sample_id)
        self.input_runinfo_dicts[sample_id] = runinfo

        input_dir = os.path.dirname(os.path.abspath(path))
        runinfo['input_dir'] = input_dir
        runinfo['profile_db'] = os.path.join(input_dir, 'PROFILE.db')

        return sample_id, runinfo


    def read_runinfo_dicts(self):
        improper = []
        missing_path = []

        for p in self.input_runinfo_paths:
            sample_id, runinfo = self.read_runinfo_dict(p)
            try:
                sample_id, runinfo = self.read_runinfo_dict(p)
            except:
                improper.append(p)
                continue

            # if things are not where they should be, we attempt to reset the directory paths.
            if not os.path.exists(runinfo['profile_db']):
                missing_path.append(p)

        if improper:
            raise ConfigError, "%s seem to be properly formatted anvio object: %s. Are you\
                                           sure these are anvio RUNINFO.cp files?" % \
                                           ('Some RUNINFO files do not' if len(improper) > 1 else "RUNINFO file does not",
                                            ', '.join(improper))

        if missing_path:
            raise ConfigError, "Anvi'o couldn't find any profile databases for %d of %d single runs you provided for merging.\
                                Are you sure your profile databases were generated during profiling? :(" % \
                                            (len(missing_path), len(self.input_runinfo_dicts))


    def sanity_check(self):
        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=self.overwrite_output_destinations)

        if not len(self.input_runinfo_paths) > 1:
            raise ConfigError, "You need to provide at least 2 RUNINFO.cp files for this program\
                                           to be useful."

        if not self.contigs_db_path:
            raise ConfigError, "You must provide a contigs database for this operation."
        if not os.path.exists(self.contigs_db_path):
            raise ConfigError, "Anvi'o couldn't find the contigs database where you said it would be :/"

        missing = [p for p in self.input_runinfo_paths if not os.path.exists(p)]
        if missing:
            raise ConfigError, "%s not found: %s." % ('Some files are' if len(missing) > 1 else "File is",
                                                                 ', '.join(missing))

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError, "You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering\
                                while also asking it to enforce it."

        self.read_runinfo_dicts()

        # test open the contigs database (and learn its hash while doing it) to make sure we don't have
        # a deal breaker just yet
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, quiet=True)
        contigs_db_hash = contigs_db.meta['contigs_db_hash']
        contigs_db.disconnect()

        # test open all profile databases to make sure you are golden with versions
        for runinfo in self.input_runinfo_dicts.values():
            sample_profile_db = anvio.db.DB(runinfo['profile_db'], anvio.__profile__version__)
            sample_profile_db.disconnect()

        if [True for v in self.input_runinfo_dicts.values() if v['merged']]:
            raise ConfigError, "This is very cute, but you can't merge already merged runs. anvio can only merge\
                                      individual profiles (which are generated through anvi-profile program). Sorry."

        if [True for v in self.input_runinfo_dicts.values() if v['blank']]:
            raise ConfigError, "Do you have a blank profile in there? Because it seems you do :/ Well, here is the problem:\
                                blank profiles are merely useful to play with a contigs database when no mapping data is\
                                available, and they are not supposed to be merged."

        for k, p in [('total_length', 'Number of nucleotides described'),
                     ('num_contigs', 'Number of contigs'),
                     ('num_splits', 'Number of splits'),
                     ('split_length', 'Split length (-L)'),
                     ('min_contig_length', 'Minimum contig length (-M)'),
                     ('min_mean_coverage', 'Minimum mean coverage (-C)'),
                     ('min_coverage_for_variability', 'Minimum coverage to report variability (-V)'),
                     ('report_variability_full', 'Report full variability (--report-variability-full)'),
                     ('profile_AA_frequencies', 'Profile AA frequencies parameter (--profile-AA-frequencies)'),
                     ('skip_SNV_profiling', 'Skip SNV profiling parameter (--skip-SNV-profiling)')]:
            v = set([r[k] for r in self.input_runinfo_dicts.values()])
            if len(v) > 1:
                raise ConfigError, "%s is not identical for all profiles to be merged, which is a \
                                    deal breaker. All profiles that are going to be merged must be\
                                    run with identical flags and parameters :/" % p

            # so we carry over this information into the runinfo dict for merged runs:
            self.run.info(k, v.pop())

        # get split names from one of the profile databases. split names must be identical across all
        self.split_names = sorted(list(dbops.get_split_names_in_profile_db(self.input_runinfo_dicts.values()[0]['profile_db'])))

        # make sure all runs were profiled using the same contigs database (if one used):
        sample_runinfos = self.input_runinfo_dicts.values()
        hashes_for_profile_dbs = set([r['contigs_db_hash'] for r in sample_runinfos])
        if len(hashes_for_profile_dbs) != 1:
            if None in hashes_for_profile_dbs:
                raise ConfigError, "It seems there is at least one run in the mix that was profiled using an\
                                          contigs database, and at least one other that was profiled without using\
                                          one. This is not good. All runs must be profiled using the same contigs\
                                          database, or all runs must be profiled without a contigs database :/"
            else:
                raise ConfigError, "It seems these runs were profiled using different contigs databases (or\
                                          different versions of the same contigs database). All runs must be\
                                          profiled using the same contigs database, or all runs must be profiled\
                                          without a contigs database :/"


        # make sure the hash for contigs db is identical across all profile databases:
        if list(hashes_for_profile_dbs)[0] != contigs_db_hash:
            raise ConfigError, "The contigs database you provided, which is identified with hash '%s', does\
                                      not seem to match the run profiles you are trying to merge, which share the\
                                      hash identifier of '%s'. What's up with that?" % (contigs_db_hash, list(hashes_for_profile_dbs)[0])


    def set_sample_id(self):
        if self.sample_id:
            utils.check_sample_id(self.sample_id)
        else:
            self.sample_id = os.path.basename(self.output_directory)
            self.sample_id = self.sample_id.replace('-', '_')
            if self.sample_id[0] in constants.digits:
                self.sample_id = 's' + self.sample_id
            utils.check_sample_id(self.sample_id)


    def is_all_samples_have_it(self, runinfo_variable):
        presence = [runinfo[runinfo_variable] for runinfo in self.input_runinfo_dicts.values() if runinfo_variable in runinfo]

        if not len(presence):
            raise ConfigError, "Something is wrong. Your profiles are not compatible with the merger. Please\
                                      re-profile everything you are trying to merge, and run merger again."

        if presence.count(True) == len(self.input_runinfo_dicts.values()):
            self.run.info(runinfo_variable, True, quiet=True)
        elif presence.count(False) == len(self.input_runinfo_dicts.values()):
            # none has it
            self.run.info(runinfo_variable, False, quiet=True)
            return
        else:
            # some weird shit must have happened.
            raise ConfigError, "Anvi'o is confused. While merging multiple runs, it seem some of the runs have the\
                                      '%s' in their PROFILE.db's, and others do not. This should never happen. Probably\
                                      your best bet is to profile everything (with of course using the same parameters)\
                                      from scratch :/" % runinfo_variable

    def merge_variable_nts_tables(self):
        self.is_all_samples_have_it('variable_nts_table')

        variable_nts_table = dbops.TableForVariability(self.profile_db_path, progress=self.progress)

        for runinfo in self.input_runinfo_dicts.values():
            sample_profile_db = dbops.ProfileDatabase(runinfo['profile_db'], quiet=True)
            sample_variable_nts_table = sample_profile_db.db.get_table_as_list_of_tuples(tables.variable_nts_table_name, tables.variable_nts_table_structure)
            sample_profile_db.disconnect()

            for tpl in sample_variable_nts_table:
                entry = tuple([variable_nts_table.next_id(tables.variable_nts_table_name)] + list(tpl[1:]))
                variable_nts_table.db_entries.append(entry)

        variable_nts_table.store()


    def merge_variable_aas_tables(self):
        self.is_all_samples_have_it('AA_frequencies_table')

        variable_aas_table = dbops.TableForAAFrequencies(self.profile_db_path, progress=self.progress)

        for runinfo in self.input_runinfo_dicts.values():
            sample_profile_db = dbops.ProfileDatabase(runinfo['profile_db'], quiet=True)
            sample_variable_aas_table = sample_profile_db.db.get_table_as_list_of_tuples(tables.variable_aas_table_name, tables.variable_aas_table_structure)
            sample_profile_db.disconnect()

            for tpl in sample_variable_aas_table:
                entry = tuple([variable_aas_table.next_id(tables.variable_aas_table_name)] + list(tpl[1:]))
                variable_aas_table.db_entries.append(entry)

        variable_aas_table.store()


    def merge_gene_coverages_tables(self):
        self.is_all_samples_have_it('gene_coverages_table')

        # create an instance from genes
        gene_coverages_table = dbops.TableForGeneCoverages(self.profile_db_path, progress=self.progress)

        # fill "genes" instance from all samples
        for runinfo in self.input_runinfo_dicts.values():
            sample_id = runinfo['sample_id']

            sample_profile_db = dbops.ProfileDatabase(runinfo['profile_db'], quiet=True)
            sample_gene_profiles = sample_profile_db.db.get_table_as_dict(tables.gene_coverages_table_name, tables.gene_coverages_table_structure)
            for g in sample_gene_profiles.values():
                gene_coverages_table.add_gene_entry(g['gene_callers_id'], g['sample_id'], g['mean_coverage'] * self.normalization_multiplier[sample_id])
            sample_profile_db.disconnect()

        gene_coverages_table.store()


    def merge_split_coverage_data(self):
        self.is_all_samples_have_it('split_coverage_values')

        output_file_path = os.path.join(self.output_directory, 'AUXILIARY-DATA.h5')
        merged_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(output_file_path, self.contigs_db_hash, create_new=True)

        # fill coverages in from all samples
        for runinfo in self.input_runinfo_dicts.values():
            input_file_path = os.path.join(os.path.dirname(runinfo['profile_db']), 'AUXILIARY-DATA.h5')
            sample_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(input_file_path, self.contigs_db_hash)

            for split_name in self.split_names:
                coverages_dict = sample_split_coverage_values.get(split_name)
                for sample_name in coverages_dict:
                    merged_split_coverage_values.append(split_name, sample_name, coverages_dict[sample_name])

            sample_split_coverage_values.close()

        merged_split_coverage_values.close()


    def set_normalization_multiplier(self):
        # WARNING: here we normalize gene coverages based on total number of reads mapped per sample.
        # this is not the best way to do it. a better way probably required all reads obtained from each
        # run, yet even that would wrongly assume equal eukaryotic contamination, etc. normalization is a bitch.
        num_reads_mapped_per_sample = {}
        for runinfo in self.input_runinfo_dicts.values():
            sample_profile_db = anvio.db.DB(runinfo['profile_db'], anvio.__profile__version__)
            num_reads_mapped_per_sample[runinfo['sample_id']] = int(sample_profile_db.get_meta_value('total_reads_mapped'))
            sample_profile_db.disconnect()

        smallest_sample_size = min(num_reads_mapped_per_sample.values())

        if smallest_sample_size == 0:
            raise ConfigError, "It seems at least one of the samples you are trying to merge has zero hits. Here is a\
                                list of all samples and number of mapped reads they have: %s." \
                                    % ', '.join(['"%s": %s' % (s, pp(num_reads_mapped_per_sample[s])) for s in num_reads_mapped_per_sample])

        for sample_id in num_reads_mapped_per_sample:
            self.normalization_multiplier[sample_id] = smallest_sample_size * 1.0 / num_reads_mapped_per_sample[sample_id]

        PRETTY = lambda x: ', '.join(['%s: %.2f' % (s, x[s]) for s in x])
        self.run.warning("anvio just set the normalization values for each sample based on\
                          how many mapped reads they contained. All normalized coverages\
                          will use this information: %s" % PRETTY(self.normalization_multiplier))


    def merge(self):
        self.sanity_check()
        self.set_sample_id()

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=self.overwrite_output_destinations)

        # init profile database
        self.profile_db_path = os.path.join(self.output_directory, 'PROFILE.db')

        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        C = lambda x: self.input_runinfo_dicts.values()[0][x]
        self.contigs_db_hash = C('contigs_db_hash')
        self.min_contig_length = C('min_contig_length')
        self.num_contigs = C('num_contigs')
        self.num_splits = C('num_splits')
        self.total_reads_mapped = C('total_reads_mapped')
        self.min_coverage_for_variability = C('min_coverage_for_variability')
        self.report_variability_full = C('report_variability_full')
        self.gene_coverages_computed = C('gene_coverages_computed')
        self.AA_frequencies_profiled = C('profile_AA_frequencies')
        self.SNVs_profiled = not C('skip_SNV_profiling')
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

        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': ','.join(self.merged_sample_ids),
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
                       'total_reads_mapped': self.total_reads_mapped,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'report_variability_full': self.report_variability_full,
                       'contigs_db_hash': self.contigs_db_hash,
                       'gene_coverages_computed': self.gene_coverages_computed}
        profile_db.create(meta_values)

        # get view data information for both contigs and splits:
        self.atomic_data_fields, self.atomic_data_for_each_run = self.read_atomic_data_tables()
        self.split_parents = self.get_split_parents()

        self.run.info('profiler_version', anvio.__profile__version__)
        self.run.info('output_dir', self.output_directory)
        self.run.info('sample_id', self.sample_id)
        self.run.info('profile_db', self.profile_db_path)
        self.run.info('merged', True)
        self.run.info('contigs_db_hash', self.contigs_db_hash)
        self.run.info('merged_sample_ids', self.merged_sample_ids)
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('num_runs_processed', len(self.merged_sample_ids))
        self.run.info('clustering_performed', not self.skip_hierarchical_clustering)

        self.set_normalization_multiplier()

        self.progress.new('Merging gene coverages tables')
        self.merge_gene_coverages_tables()
        self.progress.end()

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

        # store everything
        runinfo_serialized = os.path.join(self.output_directory, 'RUNINFO.mcp')
        self.run.info('runinfo', runinfo_serialized)
        self.run.store_info_dict(runinfo_serialized, strip_prefix=self.output_directory)

        # run CONCOCT, if otherwise is not requested:
        if not self.skip_concoct_binning and __CONCOCT_IS_AVAILABLE__:
            self.bin_contigs_concoct()

        self.run.quit()


    def get_normalized_coverage_of_split(self, target, sample_id, split_name):
        return self.atomic_data_for_each_run[target][sample_id][split_name]['mean_coverage_Q2Q3'] * self.normalization_multiplier[sample_id]


    def get_max_normalized_ratio_of_split(self, target, split_name):
        denominator = float(max(self.normalized_coverages[target][split_name].values()))

        d = {}
        for sample_id in self.merged_sample_ids:
            d[sample_id] = (self.normalized_coverages[target][split_name][sample_id] / denominator) if denominator else 0

        return d


    def get_relative_abundance_of_split(self, target, sample_id, split_name):
        denominator = float(sum(self.normalized_coverages[target][split_name].values()))
        return self.normalized_coverages[target][split_name][sample_id] / denominator if denominator else 0


    def gen_view_data_tables_from_atomic_data(self):
        essential_fields = [f for f in self.atomic_data_fields if constants.IS_ESSENTIAL_FIELD(f)]
        auxiliary_fields = [f for f in self.atomic_data_fields if constants.IS_AUXILIARY_FIELD(f)]

        # setting standard view table structure and types
        view_table_structure = ['contig'] + self.merged_sample_ids + auxiliary_fields
        view_table_types = ['text'] + ['numeric'] * len(self.merged_sample_ids) + ['text']

        # generate a dictionary for normalized coverage of each contig across samples per target
        self.normalized_coverages = {'contigs': {}, 'splits': {}}
        for target in ['contigs', 'splits']:
            for split_name in self.split_names:
                self.normalized_coverages[target][split_name] = {}
                for sample_id in self.merged_sample_ids:
                    self.normalized_coverages[target][split_name][sample_id] = self.get_normalized_coverage_of_split(target, sample_id, split_name)

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

                    for sample_id in self.merged_sample_ids:
                        if essential_field == 'normalized_coverage':
                            data_dict[split_name][sample_id] = self.normalized_coverages[target][split_name][sample_id]
                        elif essential_field == 'max_normalized_ratio':
                            data_dict[split_name][sample_id] = self.max_normalized_ratios[target][split_name][sample_id]
                        elif essential_field == 'relative_abundance':
                            data_dict[split_name][sample_id] = self.get_relative_abundance_of_split(target, sample_id, split_name)
                        else:
                            data_dict[split_name][sample_id] = self.atomic_data_for_each_run[target][sample_id][split_name][essential_field]

                # time to store the data for this view in the profile database
                table_name = '_'.join([essential_field, target])
                dbops.TablesForViews(self.profile_db_path).create_new_view(
                                                data_dict=data_dict,
                                                table_name=table_name,
                                                table_structure=view_table_structure,
                                                table_types=view_table_types,
                                                view_name=essential_field if target == 'splits' else None)

        # if SNVs were not profiled, remove all entries from variability tables:
        if not self.SNVs_profiled:
            dbops.TablesForViews(self.profile_db_path).remove(view_name='variability', table_names_to_blank=['variability_splits', 'variability_contigs'])

        self.progress.end()


    def bin_contigs_concoct(self):
        if not __CONCOCT_IS_AVAILABLE__:
            self.run.warning('The CONCOCT module is not available. Skipping the unsupervised binning\
                              with CONCOCT.')
            return

        class Args:
            pass

        args = Args()
        args.profile_db = self.profile_db_path
        args.contigs_db = self.contigs_db_path
        args.debug = self.debug

        c = concoct.CONCOCT(args)
        c.cluster()
        c.store_clusters_in_db()


    def cluster_contigs_anvio(self):
        # clustering of contigs is done for each configuration file under static/clusterconfigs/merged directory;
        # at this point we don't care what those recipes really require because we already merged and generated
        # every data file that may be required.

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

                dbops.add_hierarchical_clustering_to_db(self.profile_db_path, config_name, newick, distance=distance, linkage=linkage, make_default=config_name == constants.merged_default, run=self.run)


    def get_split_parents(self):
        parents = {}
        m = self.atomic_data_for_each_run['splits'].values()[0]
        for split in m:
            parents[split] = m[split]["__parent__"]
        return parents


    def read_atomic_data_tables(self):
        """reads atomic data for contigs and splits from the database into a dict"""
        atomic_data_table_for_each_run = {}

        for target in ['contigs', 'splits']:
            atomic_data_table_for_each_run[target] = {}

            target_table = 'atomic_data_%s' % target

            for r in self.input_runinfo_dicts.values():
                db = anvio.db.DB(r['profile_db'], dbops.get_required_version_for_db(r['profile_db']))
                atomic_data_table_for_each_run[target][r['sample_id']] = db.get_table_as_dict(target_table)

        atomic_data_table_fields = db.get_table_structure('atomic_data_splits')
        db.disconnect()

        return atomic_data_table_fields, atomic_data_table_for_each_run
