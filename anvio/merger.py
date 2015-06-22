# -*- coding: utf-8
"""The library to merge multiple profiles.

The default client of this library is under bin/anvi-merge"""


import os

import anvio
import anvio.contig
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.tables as tables
import anvio.dictio as dictio
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths

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
except ImportError, e:
    run.warning("The CONCOCT module could not be imported :( Anvi'o will still be able to perform\
                 the merging, however, the unsupervised binning results will not be available to\
                 you in the resulting merged profile database. This is what the module was upset about:\
                 '''%s''', in case you would like to fix this problem." % e)


class MultipleRuns:
    def __init__(self, args, run = run, progress = progress):
        self.progress = progress
        self.run = run 

        self.sample_id = args.sample_id
        self.merged_sample_ids = []
        self.input_runinfo_dicts = {}
        self.input_runinfo_paths = args.input
        self.split_names = None
        self.normalization_multiplier = {}
        self.profiles = []
        self.output_directory = args.output_dir
        self.skip_hierarchical_clustering = args.skip_hierarchical_clustering
        self.skip_concoct_binning = args.skip_concoct_binning
        self.skip_merging_summaries = args.skip_merging_summaries

        self.annotation_db_path = args.annotation_db_path
        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['merged']

        self.database_paths = {'ANNOTATION.db': self.annotation_db_path}

        self.overwrite_output_destinations = args.overwrite_output_destinations
        self.debug = args.debug


    def read_runinfo_dict(self, path):
        runinfo = dictio.read_serialized_object(path)
        sample_id = runinfo['sample_id']

        if not sample_id in self.merged_sample_ids:
            self.merged_sample_ids.append(sample_id)
        self.input_runinfo_dicts[sample_id] = runinfo

        input_dir = os.path.dirname(os.path.abspath(path))
        runinfo['input_dir'] = input_dir
        runinfo['profile_db'] = os.path.join(input_dir, 'PROFILE.db')
        runinfo['profile_summary_index'] = os.path.join(input_dir, 'SUMMARY.cp')

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
        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists = self.overwrite_output_destinations)

        if not len(self.input_runinfo_paths) > 1:
            raise ConfigError, "You need to provide at least 2 RUNINFO.cp files for this program\
                                           to be useful."

        if not self.annotation_db_path:
            raise ConfigError, "You must provide an annotation database for this operation."
        if not os.path.exists(self.annotation_db_path):
            raise ConfigError, "anvio couldn't find the annotation database where you said it would be :/"

        missing = [p for p in self.input_runinfo_paths if not os.path.exists(p)]
        if missing:
            raise ConfigError, "%s not found: %s." % ('Some files are' if len(missing) > 1 else "File is",
                                                                 ', '.join(missing))

        self.read_runinfo_dicts()

        if [True for v in self.input_runinfo_dicts.values() if v['merged']]:
            raise ConfigError, "This is very cute, but you can't merge already merged runs. anvio can only merge\
                                      individual profiles (which are generated through anvi-profile program). Sorry."

        for k, p in [('total_length', 'Number of nucleotides described'),
                     ('num_contigs', 'Number of contigs'),
                     ('num_splits', 'Number of splits'),
                     ('split_length', 'Split length (-L)'),
                     ('min_contig_length', 'Minimum contig length (-M)'),
                     ('min_mean_coverage', 'Minimum mean coverage (-C)'),
                     ('min_coverage_for_variability', 'Minimum coverage to report variability (-V)')]:
            v = set([r[k] for r in self.input_runinfo_dicts.values()])
            if len(v) > 1:
                raise ConfigError, "%s is not identical for all runs to be merged, which is a \
                                          deal breaker. You need to profile all runs to be merged with\
                                          identical parameters :/" % p

            # so we carry over this information into the runinfo dict for merged runs:
            self.run.info(k, v.pop())

        # get split names from one of the profile databases. split names must be identical across all 
        self.split_names = sorted(list(self.get_split_names(self.input_runinfo_dicts.values()[0]['profile_db'])))

        # make sure all runs were profiled using the same annotation database (if one used):
        sample_runinfos = self.input_runinfo_dicts.values()
        hashes_for_profile_dbs = set([r['annotation_hash'] for r in sample_runinfos])
        if len(hashes_for_profile_dbs) != 1:
            if None in hashes_for_profile_dbs:
                raise ConfigError, "It seems there is at least one run in the mix that was profiled using an\
                                          annotation database, and at least one other that was profiled without using\
                                          one. This is not good. All runs must be profiled using the same annotation\
                                          database, or all runs must be profiled without an annotation database :/"
            else:
                raise ConfigError, "It seems these runs were profiled using different annotation databases (or\
                                          different versions of the same annotation database). All runs must be\
                                          profiled using the same annotation database, or all runs must be profiled\
                                          without an annotation database :/"

        # make sure annotation hash that is common across runs is also identical to the annotation database
        annotation_db = dbops.AnnotationDatabase(self.annotation_db_path, quiet = True)
        annotation_db_hash = annotation_db.meta['annotation_hash']
        annotation_db.disconnect()

        if list(hashes_for_profile_dbs)[0] != annotation_db_hash:
            raise ConfigError, "The annotation database you provided, which is identified with hash '%s', does\
                                      not seem to match the run profiles you are trying to merge, which share the\
                                      hash identifier of '%s'. What's up with that?" % (annotation_db_hash, hashes_for_profile_dbs[0])


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
        presence = [runinfo[runinfo_variable] for runinfo in self.input_runinfo_dicts.values() if runinfo.has_key(runinfo_variable)]

        if not len(presence):
            raise ConfigError, "Something is wrong. Your profiles are not compatible with the merger. Please\
                                      re-profile everything you are trying to merge, and run merger again."

        if presence.count(True) == len(self.input_runinfo_dicts.values()):
            self.run.info(runinfo_variable, True, quiet = True)
        elif presence.count(False) == len(self.input_runinfo_dicts.values()):
            # none has it
            self.run.info(runinfo_variable, False, quiet = True)
            return
        else:
            # some weird shit must have happened.
            raise ConfigError, "anvio is confused. While merging multiple runs, it seem some of the runs have the\
                                      '%s' in their PROFILE.db's, and others do not. This should never happen. Probably\
                                      your best bet is to profile everything (with of course using the same parameters)\
                                      from scratch :/"

    def merge_variable_positions_tables(self):
        self.is_all_samples_have_it('variable_positions_table')

        variable_positions_table = dbops.TableForVariability(self.profile_db_path, anvio.__profile__version__, progress = self.progress)

        for runinfo in self.input_runinfo_dicts.values():
            sample_profile_db = dbops.ProfileDatabase(runinfo['profile_db'], quiet = True)
            sample_variable_positions_table = sample_profile_db.db.get_table_as_list_of_tuples(tables.variable_positions_table_name, tables.variable_positions_table_structure)
            sample_profile_db.disconnect()

            for tpl in sample_variable_positions_table:
                entry = tuple([variable_positions_table.next_id(tables.variable_positions_table_name)] + list(tpl[1:]))
                variable_positions_table.db_entries.append(entry)

        variable_positions_table.store()


    def merge_gene_coverages_tables(self):
        self.is_all_samples_have_it('gene_coverages_table')

        # create an instance from genes
        gene_coverages_table = dbops.TableForGeneCoverages(self.profile_db_path, anvio.__profile__version__, progress = self.progress)

        # fill "genes" instance from all samples
        for runinfo in self.input_runinfo_dicts.values():
            sample_id = runinfo['sample_id']

            sample_profile_db = dbops.ProfileDatabase(runinfo['profile_db'], quiet = True)
            sample_gene_profiles = sample_profile_db.db.get_table_as_dict(tables.gene_coverages_table_name, tables.gene_coverages_table_structure)
            for g in sample_gene_profiles.values():
                gene_coverages_table.add_gene_entry(g['prot'], g['sample_id'], g['mean_coverage'] * self.normalization_multiplier[sample_id])
            sample_profile_db.disconnect()

        gene_coverages_table.store()


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

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists = self.overwrite_output_destinations)

        # init profile database
        self.profile_db_path = os.path.join(self.output_directory, 'PROFILE.db')

        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        self.annotation_hash = self.input_runinfo_dicts.values()[0]['annotation_hash']
        self.min_contig_length = self.input_runinfo_dicts.values()[0]['min_contig_length']
        self.num_contigs = self.input_runinfo_dicts.values()[0]['num_contigs']
        self.num_splits = self.input_runinfo_dicts.values()[0]['num_splits']
        self.min_coverage_for_variability = self.input_runinfo_dicts.values()[0]['min_coverage_for_variability']
        self.total_length = self.input_runinfo_dicts.values()[0]['total_length']
        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': ','.join(self.merged_sample_ids),
                       'merged': True,
                       'contigs_clustered': not self.skip_hierarchical_clustering,
                       'default_view': 'mean_coverage',
                       'min_contig_length': self.min_contig_length,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'num_contigs': self.num_contigs,
                       'num_splits': self.num_splits,
                       'total_length': self.total_length,
                       'annotation_hash': self.annotation_hash}
        profile_db.create(meta_values)

        # get metadata information for both contigs and splits:
        self.metadata_fields, self.metadata_for_each_run = self.read_metadata_tables()
        self.split_parents = self.get_split_parents()

        self.run.info('profiler_version', anvio.__profile__version__)
        self.run.info('output_dir', self.output_directory)
        self.run.info('sample_id', self.sample_id)
        self.run.info('profile_db', self.profile_db_path)
        self.run.info('merged', True)
        self.run.info('annotation_hash', self.annotation_hash)
        self.run.info('merged_sample_ids', self.merged_sample_ids)
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('num_runs_processed', len(self.merged_sample_ids))
        #self.run.info('num_splits_found', pp(len(self.contigs.values()[0])))
        #self.run.info('contigs_total_length', pp(sum([len(s) for s in self.contigs.values()[0]])))
        self.run.info('clustering_performed', not self.skip_hierarchical_clustering)

        self.set_normalization_multiplier()

        self.progress.new('Merging gene coverages tables')
        self.merge_gene_coverages_tables()
        self.progress.end()

        self.progress.new('Merging variable positions tables')
        self.merge_variable_positions_tables()
        self.progress.end()

        if not self.skip_merging_summaries:
            self.progress.new('Generating merged summary')
            summary_dir, profile_summary_index = self.merge_split_summaries()
            self.progress.end()
            self.run.info('profile_summary_dir', summary_dir)
            self.run.info('profile_summary_index', profile_summary_index)

        # critical part:
        self.merge_metadata_files()

        # We cluster? Note: the check is being done in the function!
        self.cluster_contigs_anvio()

        self.progress.end()

        # store everything
        runinfo_serialized = os.path.join(self.output_directory, 'RUNINFO.mcp')
        self.run.info('runinfo', runinfo_serialized)
        self.run.store_info_dict(runinfo_serialized, strip_prefix = self.output_directory)

        # run CONCOCT, if otherwise is not requested:
        if not self.skip_concoct_binning and __CONCOCT_IS_AVAILABLE__:
            self.bin_contigs_concoct()

        self.run.quit()


    def get_normalized_coverage_of_split(self, target, sample_id, split_name):
        return self.metadata_for_each_run[target][sample_id][split_name]['normalized_coverage'] * self.normalization_multiplier[sample_id]


    def get_max_normalized_ratio_of_split(self, target, split_name):
        denominator = float(max(self.normalized_coverages[target][split_name].values()))

        d = {}
        for sample_id in self.merged_sample_ids:
            d[sample_id] = (self.normalized_coverages[target][split_name][sample_id] / denominator) if denominator else 0

        return d


    def get_relative_abundance_of_split(self, target, sample_id, split_name):
        denominator = float(sum(self.normalized_coverages[target][split_name].values()))
        return self.normalized_coverages[target][split_name][sample_id] / denominator if denominator else 0


    def merge_metadata_files(self):
        essential_fields = [f for f in self.metadata_fields if constants.IS_ESSENTIAL_FIELD(f)]
        auxiliary_fields = [f for f in self.metadata_fields if constants.IS_AUXILIARY_FIELD(f)]

        views_table = dbops.TableForViews(self.profile_db_path, anvio.__profile__version__, progress = self.progress)

        # setting standard metadata table structure and types
        merged_mtable_structure = ['contig'] + self.merged_sample_ids + auxiliary_fields
        merged_mtable_types = ['text'] + ['numeric'] * len(self.merged_sample_ids) + ['text']

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

        self.progress.new('Generating metadata tables')
        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet = True)
        for target in ['contigs', 'splits']:
            for essential_field in essential_fields:
                self.progress.update('Processing %s for %s ...' % (essential_field, target))

                target_table = '_'.join([essential_field, target])

                m = {}
                for split_name in self.split_names:
                    m[split_name] = {'__parent__': self.split_parents[split_name]}

                    for sample_id in self.merged_sample_ids:
                        if essential_field == 'normalized_coverage':
                            m[split_name][sample_id] = self.normalized_coverages[target][split_name][sample_id]
                        elif essential_field == 'max_normalized_ratio':
                            m[split_name][sample_id] = self.max_normalized_ratios[target][split_name][sample_id]
                        elif essential_field == 'relative_abundance':
                            m[split_name][sample_id] = self.get_relative_abundance_of_split(target, sample_id, split_name)
                        else:
                            m[split_name][sample_id] = self.metadata_for_each_run[target][sample_id][split_name][essential_field]

                # variable 'm' for the essential field is now ready to be its own table:
                profile_db.db.create_table(target_table, merged_mtable_structure, merged_mtable_types)
                db_entries = [tuple([split_name] + [m[split_name][h] for h in merged_mtable_structure[1:]]) for split_name in self.split_names]
                profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (target_table, ','.join(['?'] * len(merged_mtable_structure))), db_entries)

                if target == 'splits':
                    views_table.append(essential_field, target_table)

        profile_db.disconnect()
        self.progress.end()

        # store views in the database
        views_table.store()


    def bin_contigs_concoct(self):
        if not __CONCOCT_IS_AVAILABLE__:
            self.run.warning('The CONCOCT module is not available. Skipping the unsupervised binning\
                              with CONCOCT.')
            return

        class Args:
            pass

        args = Args()
        args.profile_db = self.profile_db_path
        args.annotation_db = self.annotation_db_path
        args.debug = self.debug

        c = concoct.CONCOCT(args)
        c.cluster()
        c.store_clusters_in_db()


    def cluster_contigs_anvio(self):
        # clustering of contigs is done for each configuration file under static/clusterconfigs/merged directory;
        # at this point we don't care what those recipes really require because we already merged and generated
        # every metadata file that may be required.
        clusterings = []

        if not self.skip_hierarchical_clustering:
            for config_name in self.clustering_configs:
                config_path = self.clustering_configs[config_name]

                config = ClusteringConfiguration(config_path, self.output_directory, db_paths = self.database_paths, row_ids_of_interest = self.split_names)

                try:
                    newick = clustering.order_contigs_simple(config, progress = self.progress)
                except Exception as e:
                    self.run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
                    self.progress.end()
                    continue

                clusterings.append(config_name)
                db_entries = tuple([config_name, newick])

                profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet = True)
                profile_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % tables.clusterings_table_name, db_entries)
                profile_db.disconnect()

        self.run.info('available_clusterings', clusterings)
        self.run.info('default_clustering', constants.merged_default)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('default_clustering', constants.merged_default)
        profile_db.db.set_meta_value('available_clusterings', ','.join(clusterings))
        profile_db.disconnect()


    def merge_split_summaries(self):
        merged_summary_index = {}
        merged_summary_index_path = os.path.join(self.output_directory, 'SUMMARY.cp')
        summary_dir = filesnpaths.gen_output_directory(os.path.join(self.output_directory, 'SUMMARY'), delete_if_exists = True)


        # read all index files per run into a dict here, so the access is easier from within
        # the for loop below
        run_sum_indices = {}
        for runinfo  in self.input_runinfo_dicts.values():
            run_sum_indices[runinfo['sample_id']] = dictio.read_serialized_object(runinfo['profile_summary_index'])

        for i in range(0, len(self.split_names)):
            self.progress.update('merging summaries for splits %s of %s' % (i + 1, len(self.split_names)))
            split_name = self.split_names[i]

            merged_summary = {}
            for runinfo in self.input_runinfo_dicts.values(): 
                run_split_summary = dictio.read_serialized_object(os.path.join(runinfo['input_dir'], run_sum_indices[runinfo['sample_id']][split_name]))
                merged_summary[runinfo['sample_id']] = run_split_summary[runinfo['sample_id']]

            merged_split_summary_path = os.path.join(summary_dir, os.path.basename(run_sum_indices[runinfo['sample_id']][split_name]))
            dictio.write_serialized_object(merged_summary, merged_split_summary_path)
            merged_summary_index[split_name] = merged_split_summary_path

        self.progress.update('Serializing merged split summary index ...')
        dictio.write_serialized_object(dictio.strip_prefix_from_dict_values(merged_summary_index, self.output_directory),\
                                           merged_summary_index_path)

        return summary_dir, merged_summary_index_path


    def get_split_parents(self):
        parents = {}
        m = self.metadata_for_each_run['splits'].values()[0]
        for split in m:
            parents[split] = m[split]["__parent__"]
        return parents


    def read_metadata_tables(self):
        """reads metadata files of contigs and splits into a dict"""
        metadata_for_each_run = {}

        for target in ['contigs', 'splits']:
            metadata_for_each_run[target] = {}

            target_table = 'metadata_%s' % target

            for r in self.input_runinfo_dicts.values():
                db = anvio.db.DB(r['profile_db'], anvio.__profile__version__)
                metadata_for_each_run[target][r['sample_id']] = db.get_table_as_dict(target_table)

        metadata_fields = db.get_table_structure('metadata_splits')
        db.disconnect()

        return metadata_fields, metadata_for_each_run


    def get_split_names(self, profile_db_path):
        profile_db = dbops.ProfileDatabase(profile_db_path)
        split_names = profile_db.db.get_single_column_from_table('metadata_splits', 'contig')
        profile_db.disconnect()

        return split_names
