# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import math
import hashlib

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.summarizer as summarizer
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.drivers.blast import BLAST
from anvio.drivers.muscle import Muscle
from anvio.drivers.diamond import Diamond
from anvio.drivers.mcl import MCL
from anvio.clusteringconfuguration import ClusteringConfiguration

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


class GenomeStorage(object):
    def __init__(self, args=None, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        self.storage_path = None
        self.genomes_storage = None
        self.genome_names_to_focus = None

        self.genomes = {}
        self.hash_to_genome_name = {}

        self.functions_are_available = False
        self.function_annotation_sources = set([])


    def load_genomes_descriptions(self, skip_functions=False):
        """Reads internal and external genome files, populates self.genomes"""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        input_file_for_internal_genomes = A('internal_genomes')
        input_file_for_external_genomes = A('external_genomes')

        fields_for_internal_genomes_input = ['name', 'bin_id', 'collection_id', 'profile_db_path', 'contigs_db_path']
        fields_for_external_genomes_input = ['name', 'contigs_db_path']

        internal_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_for_internal_genomes, expected_fields=fields_for_internal_genomes_input) if input_file_for_internal_genomes else {}
        external_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_for_external_genomes, expected_fields=fields_for_external_genomes_input) if input_file_for_external_genomes else {}

        self.internal_genome_names = internal_genomes_dict.keys()
        self.external_genome_names = external_genomes_dict.keys()

        if not self.internal_genome_names and not self.external_genome_names:
            raise ConfigError, "You in fact tried to create a genomes storage file without providing any internal or external genome\
                                descriptions! You got 5 anvi'o points for being awesome, but this is not gonna work since you really\
                                need to provide at least one of those descriptions :/"

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

                # while we are going through all genomes and reconstructing self.genomes for the first time,
                # let's add the 'name' attribute in it as well.'
                self.genomes[genome_name]['name'] = genome_name

        # add hashes for each genome in the self.genomes dict. this will allow us to see whether the HDF file already contains
        # all the information we need.
        for genome_name in self.external_genome_names:
            self.genomes[genome_name]['genome_hash'] = self.get_genome_hash_for_external_genome(self.genomes[genome_name])
        for genome_name in self.internal_genome_names:
            self.genomes[genome_name]['genome_hash'] = self.get_genome_hash_for_internal_genome(self.genomes[genome_name])

        # if the client is not interested in functions, skip the rest.
        if skip_functions:
            self.functions_are_available = False
            return

        # check whether function calls are available for all genomes involved, and whether function sources for each genome is identical
        function_annotation_sources_per_genome = {}
        all_function_annotation_sources_observed = set([])
        for genome_name in self.genomes:
            g = self.genomes[genome_name]
            contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
            sources = contigs_db.meta['gene_function_sources']
            contigs_db.disconnect()

            if sources:
                function_annotation_sources_per_genome[genome_name] = sources
                all_function_annotation_sources_observed.update(sources)

        if not len(all_function_annotation_sources_observed):
            self.run.warning("None of your genomes seem to have any functional annotation. No biggie. Things will continue to work. But\
                              then your genomes have no functional annotation. It is sad.")
        else:
            # this guy down below fills in the self.function_annotation_sources with function annotation sources
            # that are common to all genomes.
            for sources in function_annotation_sources_per_genome.values():
                if not sources:
                    continue

                if not(self.function_annotation_sources):
                    self.function_annotation_sources.update(sources)
                else:
                    self.function_annotation_sources = self.function_annotation_sources.intersection(sources)

            function_annotation_sources_some_genomes_miss = all_function_annotation_sources_observed.difference(self.function_annotation_sources)

            if not len(self.function_annotation_sources):
                # none of the functions are common
                self.run.warning("Although some of your genomes had some functional annotations, none of them were common to all genomes :/\
                                  Anvi'o will continue working with them, but you will have no functions available to you downstream. Just\
                                  so you know, these are the annotation sources observed at least once in at least one of your genomes: '%s'" % \
                                                                    (', '.join(all_function_annotation_sources_observed)))
                self.functions_are_available = False
            else:
                self.functions_are_available = True

                # good. here we know some functions are available, but let's get some further understanding, and report it to the user, you know,
                # because we're nice:
                if len(function_annotation_sources_some_genomes_miss):
                    # some functions were missing from some genomes
                    self.run.warning("Anvi'o has good news and bad news for you (very balanced, as usual). The good news is that there are some\
                                      funciton annotation sources that are common to all of your genomes, and they will be used whenever\
                                      it will be appropriate. Here they are: '%s'. The bad news is you had more functiona annotation sources,\
                                      but they were not common to all genomes. Here they are so you can say your goodbyes to them (because\
                                      they will not be used): '%s'" % \
                                            (', '.join(self.function_annotation_sources), ', '.join(function_annotation_sources_some_genomes_miss)))
                else:
                    # every function ever observed is common to all genomes.
                    self.run.warning("Good news! Anvi'o found all these functions that are common to all of your genomes and will use them for\
                                      downstream analyses and is very proud of you: '%s'." % (', '.join(self.function_annotation_sources)), lc='green')


    def init_genomes_data_storage(self):
        """Initializes an existing genomes storage by reading everything about genomes of interest"""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.storage_path = A('genomes_storage')
        self.genome_names_to_focus = A('genome_names')

        # let's take care of the genome names to focus, if there are any, first. 
        if self.genome_names_to_focus:
            if filesnpaths.is_file_exists(self.genome_names_to_focus, dont_raise=True):
                self.genome_names_to_focus = utils.get_column_data_from_TAB_delim_file(self.genome_names_to_focus, column_indices=[0])[0]
            else:
                self.genome_names_to_focus = [g.strip() for g in self.genome_names_to_focus.split(',')]

            self.run.warning("A subset of genome names is found, and anvi'o will focus only on to those.")

        filesnpaths.is_proper_genomes_storage_file(self.storage_path)

        self.genomes_storage = auxiliarydataops.GenomesDataStorage(self.storage_path, db_hash = None, genome_names_to_focus=self.genome_names_to_focus, ignore_hash = True)
        self.genomes_storage_hash = self.genomes_storage.get_storage_hash()

        self.genomes = self.genomes_storage.get_genomes_dict()

        self.external_genome_names = [g for g in self.genomes if self.genomes[g]['external_genome']]
        self.internal_genome_names = [g for g in self.genomes if not self.genomes[g]['external_genome']]

        for genome_name in self.genomes:
            self.hash_to_genome_name[self.genomes[genome_name]['genome_hash']] = genome_name


    def get_functions_dict_from_contigs_db(self, contigs_db_path):
        if not self.functions_are_available:
            return {}

        class Args: pass
        args = Args()
        args.contigs_db = contigs_db_path
        contigs_super = dbops.ContigsSuperclass(args, r=anvio.terminal.Run(verbose=False))
        contigs_super.init_functions(requested_sources=self.function_annotation_sources)

        return contigs_super.gene_function_calls_dict


    def sanity_check(self):
        """Make sure self.genomes is good to go"""

        # make sure genes are called in every contigs db:
        genomes_missing_gene_calls = [g for g in self.genomes if not self.genomes[g]['genes_are_called']]
        if len(genomes_missing_gene_calls):
            raise ConfigError, 'Genes must have been called during the generation of contigs database for this workflow to work. However,\
                                these external genomes do not have gene calls: %s' % (', '.join(genomes_missing_gene_calls))

        # if two contigs db has the same hash, we are kinda f'd:
        if len(set([self.genomes[genome_name]['genome_hash'] for genome_name in self.external_genome_names])) != len(self.external_genome_names):
            raise ConfigError, 'Not all hash values are unique across all contig databases you provided. Something\
                                very fishy is going on :/'


        if len(set([self.genomes[genome_name]['genome_hash'] for genome_name in self.internal_genome_names])) != len(self.internal_genome_names):
            raise ConfigError, "Not all hash values are unique across internal genomes. This is almost impossible to happen unless something very\
                                wrong with your workflow :/ Please let the developers know if you can't figure this one out"

        # make sure HMMs for SCGs were run for every contigs db:
        genomes_missing_hmms_for_scgs =  [g for g in self.genomes if not self.genomes[g]['hmms_for_scgs_were_run']]
        if len(genomes_missing_hmms_for_scgs):
            if len(genomes_missing_hmms_for_scgs) == len(self.genomes):
                raise ConfigError, "The contigs databases you are using for this analysis are missing HMMs for single-copy core genes. In other words,\
                                    you don't seem to have run `anvi-run-hmms` on them. Although it is perfectly legal to have anvi'o contigs databases\
                                    without HMMs run on SCGs, the current pangenomic workflow does not want to deal with this :( Sorry!"
            else:
                raise ConfigError, "Some of the genomes you have for this analysis are missing HMM hits for SCGs (%d of %d of them, to be precise). You\
                                    can run `anvi-run-hmms` on them to recover from this. Here is the list: %s" % \
                                                    (len(genomes_missing_hmms_for_scgs), len(self.genomes), ','.join(genomes_missing_hmms_for_scgs))

        # make sure genome names are not funny (since they are going to end up being db variables soon)
        [utils.is_this_name_OK_for_database('genome name "%s"' % genome_name, genome_name) for genome_name in self.genomes]


    def create_genomes_data_storage(self):
        """Creates an HDF5 file storing all genome related information for later access."""

        # some simple checks
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.storage_path = A('output_file')
        if not self.storage_path:
            self.storage_path = "GENOMES.h5"
        else:
            if not self.storage_path.endswith('-GENOMES.h5'):
                raise ConfigError, "The genomes storage file must end with '-GENOMES.h5'. Anvi'o developers do know how ridiculous\
                                    this requirement sounds like, but if you have seen the things they did, you would totally\
                                    understand why this is necessary."

        filesnpaths.is_output_file_writable(self.storage_path)

        # let's read those internal and external genome files the user sent.
        self.load_genomes_descriptions()

        # this will populate self.genomes with relevant data that can be learned about these genomes such as 'avg_gene_length',
        # 'num_splits', 'num_contigs', 'num_genes', 'percent_redundancy', 'gene_caller_ids', 'total_length', 'partial_gene_calls',
        # 'percent_complete', 'num_genes_per_kb', 'gc_content'.
        self.init_internal_genomes()
        self.init_external_genomes()

        # make sure it is OK to go with self.genomes
        self.sanity_check()

        # here we create a signature for the storage itself by concatenating all hash values from all genomes. even if one
        # split is added or removed to any of these genomes will change this signature. since we will tie this information
        # to the profile database we will generate for the pangenome analysis, even if one split is added or removed from any
        # of the genomes will make sure that the profile databases from this storage and storage itself are not compatible:
        storage_hash = hashlib.sha224('_'.join(self.genomes[genome_name]['genome_hash'] for genome_name in self.genomes)).hexdigest()[0:8]

        self.genomes_storage = auxiliarydataops.GenomesDataStorage(self.storage_path, storage_hash, create_new=True)
        self.genomes_storage.fp.attrs['functions_are_available'] = self.functions_are_available

        # some silly stuff for later fun
        num_gene_calls_added_total = 0
        num_partial_gene_calls_total = 0

        # main loop
        genome_names_to_go_through = sorted(self.genomes.keys())
        for genome_name in genome_names_to_go_through:
            self.progress.new('Initializing genomes')
            self.progress.update('%s ...' % genome_name)
            g = self.genomes[genome_name]

            self.genomes_storage.add_genome(genome_name, g)

            num_gene_calls_added = 0
            num_partial_gene_calls = 0

            contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
            protein_sequences_dict = contigs_db.db.get_table_as_dict(t.gene_protein_sequences_table_name)
            contigs_db.disconnect()

            functions_dict = self.get_functions_dict_from_contigs_db(g['contigs_db_path'])

            for gene_caller_id in g['gene_caller_ids']:
                partial_gene_call = gene_caller_id in g['partial_gene_calls']

                functions = []
                if gene_caller_id in functions_dict:
                    for annotation_source in self.function_annotation_sources:
                        if annotation_source in functions_dict[gene_caller_id]:
                            annotation_tuple = functions_dict[gene_caller_id][annotation_source]
                            if annotation_tuple:
                                functions.append((annotation_source, '%s|||%s' % (annotation_tuple[0], annotation_tuple[1])),)

                self.genomes_storage.add_gene_call_data(genome_name,
                                                        gene_caller_id,
                                                        sequence=protein_sequences_dict[gene_caller_id]['sequence'],
                                                        partial=partial_gene_call,
                                                        functions=functions)

                num_gene_calls_added += 1
                if partial_gene_call:
                    num_partial_gene_calls += 1

            self.progress.end()

            if genome_names_to_go_through.index(genome_name) == 0:
                self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)), cut_after=120, nl_before = 1)
            elif genome_names_to_go_through.index(genome_name) == len(genome_names_to_go_through) -1:
                self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)), cut_after=120, nl_after = 1)
            else:
                self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)), cut_after=120)

            num_gene_calls_added_total += num_gene_calls_added
            num_partial_gene_calls_total += num_partial_gene_calls


        self.run.info('The new genomes storage', '%s (signature: %s)' % (self.storage_path, storage_hash))
        self.run.info('Number of genomes', '%s (internal: %s, external: %s)' % (pp(len(self.genomes)), pp(len(self.internal_genome_names)), pp(len(self.external_genome_names))))
        self.run.info('Number of gene calls', '%s' % pp(num_gene_calls_added_total))
        self.run.info('Number of partial gene calls', '%s' % pp(num_partial_gene_calls_total))

        self.genomes_storage.close()


    def init_external_genomes(self):
        self.progress.new('Initializing external genomes')
        for genome_name in self.external_genome_names:
            c = self.genomes[genome_name]
            c['external_genome'] = True

            self.progress.update('working on %s' % (genome_name))

            contigs_db_summary = summarizer.get_contigs_db_info_dict(c['contigs_db_path'])

            for key in contigs_db_summary:
                c[key] = contigs_db_summary[key]

        self.progress.end()

        self.run.info('External genomes', '%d found.' % len(self.external_genome_names))


    def get_unique_profile_db_path_to_internal_genome_name_dict(self):
        """Returns a dictionary to bind all genome names that originate from the same profile db"""

        unique_profile_db_path_to_internal_genome_name = {}

        for profile_path in set([self.genomes[g]['profile_db_path'] for g in self.internal_genome_names]):
            unique_profile_db_path_to_internal_genome_name[profile_path] = [g for g in self.internal_genome_names if self.genomes[g]['profile_db_path'] == profile_path]

        return unique_profile_db_path_to_internal_genome_name


    def init_internal_genomes(self):
        self.progress.new('Initializing internal genomes')

        # to not initialize things over and over again:
        unique_profile_db_path_to_internal_genome_name = self.get_unique_profile_db_path_to_internal_genome_name_dict()

        for profile_db_path in unique_profile_db_path_to_internal_genome_name:
            self.collections = ccollections.Collections()
            self.collections.populate_collections_dict(profile_db_path)

            for genome_name in unique_profile_db_path_to_internal_genome_name[profile_db_path]:
                self.progress.update('working on %s' % (genome_name))
                c = self.genomes[genome_name]
                c['external_genome'] = False

                dbops.is_profile_db_and_contigs_db_compatible(c['profile_db_path'], c['contigs_db_path'])

                split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(c)

                # here we are using the get_contigs_db_info_dict function WITH split names we found in the collection
                # which returns a partial summary from the contigs database focusing only those splits. a small workaround
                # to be able to use the same funciton for bins in collections:
                summary_from_contigs_db_summary = summarizer.get_contigs_db_info_dict(c['contigs_db_path'], split_names=split_names_of_interest)
                for key in summary_from_contigs_db_summary:
                    c[key] = summary_from_contigs_db_summary[key]

        self.progress.end()

        self.run.info('Internal genomes', '%d have been initialized.' % len(self.internal_genome_names))


    def get_genome_hash_for_external_genome(self, entry):
        dbops.is_contigs_db(entry['contigs_db_path'])
        contigs_db = dbops.ContigsDatabase(entry['contigs_db_path'])
        genome_hash = contigs_db.meta['contigs_db_hash']
        contigs_db.disconnect()

        return genome_hash


    def get_genome_hash_for_internal_genome(self, entry):
        dbops.is_contigs_db(entry['contigs_db_path'])
        split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(entry)
        contigs_db = dbops.ContigsDatabase(entry['contigs_db_path'])
        genome_hash = hashlib.sha224('_'.join([''.join(split_names_of_interest), contigs_db.meta['contigs_db_hash']])).hexdigest()[0:12]
        contigs_db.disconnect()

        return genome_hash


    def get_split_names_of_interest_for_internal_genome(self, entry):
        dbops.is_profile_db(entry['profile_db_path'])
        # get splits of interest:
        class Args: pass
        args = Args()
        args.profile_db = entry['profile_db_path']
        args.collection_name = entry['collection_id']
        args.bin_id = entry['bin_id']

        split_names_of_interest = list(ccollections.GetSplitNamesInBins(args).get_split_names_only())

        if not len(split_names_of_interest):
            raise ConfigError, "There are 0 splits defined for bin id %s in collection %s..." % (entry['bin_id'], entry['collection_id'])

        return split_names_of_interest


class Pangenome(GenomeStorage):
    def __init__(self, args=None, run=run, progress=progress):
        GenomeStorage.__init__(self, args, run, progress)
        self.init_genomes_data_storage()

        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.project_name = A('project_name')
        self.output_dir = A('output_dir')
        self.num_threads = A('num_threads')
        self.skip_alignments = A('skip_alignments')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')
        self.min_percent_identity = A('min_percent_identity')
        self.PC_min_occurrence = A('min_occurrence')
        self.mcl_inflation = A('mcl_inflation')
        self.sensitive = A('sensitive')
        self.maxbit = A('maxbit')
        self.use_ncbi_blast = A('use_ncbi_blast')
        self.exclude_partial_gene_calls = A('exclude_partial_gene_calls')

        # when it is time to organize PCs
        self.linkage = A('linkage') or constants.linkage_method_default
        self.distance = A('distance') or constants.distance_metric_default

        self.log_file_path = None

        # to be filled during init:
        self.protein_sequences_dict = {}
        self.view_data = {}
        self.view_data_presence_absence = {}
        self.additional_view_data = {}


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
                       'maxbit': self.maxbit,
                       'exclude_partial_gene_calls': self.exclude_partial_gene_calls,
                       'gene_alignments_computed': False if self.skip_alignments else True,
                       'genomes_storage_hash': self.genomes_storage_hash,
                       'project_name': self.project_name
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
            raise ConfigError, "Please set a project name, and be prepared to see it around as (1) anvi'o will use\
                                that name to set the output directory and to name various output files such as the\
                                databases that will be generated at the end of the process. If you set your own output\
                                directory name, you can have multiple projects in it and all of those projects can use\
                                the same intermediate files whenever possible."

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

        if not isinstance(self.maxbit, float):
            raise ConfigError, "maxbit value must be of type float :("

        if self.maxbit < 0 or self.maxbit > 1:
            raise ConfigError, "Well. maxbit must be between 0 and 1. Yes. Very boring."

        if not isinstance(self.min_percent_identity, float):
            raise ConfigError, "Minimum percent identity value must be of type float :("

        if self.min_percent_identity < 0 or self.min_percent_identity > 100:
            raise ConfigError, "Minimum percent identity must be between 0%% and 100%%. Although your %.2f%% is\
                                pretty cute, too." % self.min_percent_identity


        if len([c for c in self.genomes.values() if 'genome_hash' not in c]):
            raise ConfigError, "self.genomes does not seem to be a properly formatted dictionary for\
                                the anvi'o class Pangenome."

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
            except Exception as e:
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


    def process_protein_clusters(self, protein_clusters_dict):
        self.progress.new('Generating view data')
        self.progress.update('...')

        def store_file(data, path, headers=None):
            if not headers:
                headers = ['contig'] + sorted(data.values()[0].keys())

            utils.store_dict_as_TAB_delimited_file(data, path, headers=headers)

            return path

        PCs = protein_clusters_dict.keys()

        for PC in PCs:
            self.view_data[PC] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.view_data_presence_absence[PC] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.additional_view_data[PC] = {'num_genes_in_pc': 0, 'num_genomes_pc_has_hits': 0, 'SCG': 0}

            for gene_entry in protein_clusters_dict[PC]:
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

        for PC in PCs:
            if PC not in PCs_of_interest:
                self.view_data.pop(PC)
                self.view_data_presence_absence.pop(PC)
                self.additional_view_data.pop(PC)

        if self.PC_min_occurrence > 1:
            self.run.info('PCs min occurrence', '%d (the filter removed %s PCs)' % (self.PC_min_occurrence, (len(protein_clusters_dict) - len(PCs_of_interest))))

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
        #             CHEATING THE SYSTEM FOR AN ENHANCED CLUSTERING CONFIGURATION
        ########################################################################################
        # so we want to use the clustering configurations for pan genomomic analyses to order
        # protein clusters. however, we want to add something into the clustering configuraiton
        # file, which depends on the number of genomes we have. this addition is 'num_genomes_pc_has_hits'
        # data, which pulls together protein clusters that are distributed across genomes similarly based
        # on this extra bit of inofrmation. becasue the clustering configurations framework in anvi'o
        # does not allow us to have variable information in these recipes, we are going to generate one
        # on the fly to have a more capable one.

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

            # use it to generate a clustering configuration instance:
            clustering_configuration = ClusteringConfiguration(enhanced_config_path, self.output_dir, db_paths={'PAN.db': self.pan_db_path})

            try:
                clustering_id, newick = clustering.order_contigs_simple(clustering_configuration, distance=self.distance, linkage=self.linkage, progress=self.progress)
            except Exception as e:
                self.run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
                self.progress.end()
                continue

            _, distance, linkage = clustering_id.split(':')

            dbops.add_hierarchical_clustering_to_db(self.pan_db_path,
                                                    config_name,
                                                    newick,
                                                    distance=distance,
                                                    linkage=linkage,
                                                    make_default=config_name == constants.pan_default,
                                                    run=self.run)


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

        for h in ['percent_complete', 'percent_redundancy']:
            if h in self.genomes.values()[0]:
                headers.append(h)

        headers.extend(['gc_content', 'num_genes', 'avg_gene_length', 'num_genes_per_kb'])

        for c in self.genomes.values():
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
            raise ConfigError, "Well, MCL likes its inflation parameter in 'float' form..."

        if self.mcl_inflation > 100 or self.mcl_inflation < 0.1:
            raise ConfigError, "MCL inflation parameter should have a reasonable value :/ Like between 0.1 and 100.0."

        if not isinstance(self.genomes, type({})):
            raise ConfigError, "self.genomes must be a dict. Anvi'o needs an adult :("

        if len(self.genomes) < 2:
            raise ConfigError, "There must be at least two genomes for this workflow to work. You have like '%d' of them :/" \
                    % len(self.genomes)

        if not self.skip_alignments:
            try:
                Muscle()
            except ConfigError, e:
                raise ConfigError, "It seems things are not quite in order. Anvi'o does not know what is wrong, but you will\
                                    see the actual error message in second. You can either try to address whatever causes this problem,\
                                    or you can use the `--skip-alignments` flag to avoid it completely (although the latter would clearly\
                                    be the easiest way out of this, anvi'o may stop thinking so highly of you if you choose to do that,\
                                    just FYI). Here is the actual error you got: %s." % str(e).replace('\n', '').replace('  ', ' ')

        self.check_params()

        self.run.log_file_path = self.log_file_path
        self.run.info('Args', (str(self.args)), quiet=True)


    def store_protein_clusters(self, protein_clusters_dict):
        self.progress.new('Storing protein clusters in the database')
        self.progress.update('...')

        table_for_protein_clusters = dbops.TableForProteinClusters(self.pan_db_path, run=self.run, progress=self.progress)

        num_genes_in_protein_clusters = 0
        for pc_name in protein_clusters_dict:
            for gene_entry in protein_clusters_dict[pc_name]:
                table_for_protein_clusters.add(gene_entry)
                num_genes_in_protein_clusters += 1

        self.progress.end()

        table_for_protein_clusters.store()

        pan_db = dbops.PanDatabase(self.pan_db_path, quiet=True)
        pan_db.db.set_meta_value('num_protein_clusters', len(protein_clusters_dict))
        pan_db.db.set_meta_value('num_genes_in_protein_clusters', num_genes_in_protein_clusters)
        pan_db.disconnect()

        self.run.info('protein clusters info', '%d PCs stored in the database' % len(protein_clusters_dict))


    def gen_protein_clusters_dict_from_mcl_clusters(self, mcl_clusters):
        self.progress.new('Generating the protein clusters dictionary from raw MCL clusters')
        self.progress.update('...')

        protein_clusters_dict = {}

        for PC in mcl_clusters:
            protein_clusters_dict[PC] = []

            for entry_hash, gene_caller_id in [e.split('_') for e in mcl_clusters[PC]]:
                try:
                    genome_name = self.hash_to_genome_name[entry_hash]
                except KeyError:
                    self.progress.end()
                    raise ConfigError, "Something horrible happened. This can only happen if you started a new analysis with\
                                        additional genomes without cleaning the previous work directory. Sounds familiar?"

                protein_clusters_dict[PC].append({'gene_caller_id': int(gene_caller_id), 'protein_cluster_id': PC, 'genome_name': genome_name, 'alignment_summary': ''})

        self.progress.end()

        return protein_clusters_dict


    def compute_alignments_for_PCs(self, protein_clusters_dict):
        if self.skip_alignments:
            self.run.warning('Skipping gene alignments.')
            return

        r = terminal.Run()
        r.verbose = False

        muscle = Muscle(run=r)

        self.progress.new('Aligning genes in protein sequences')
        self.progress.update('...')
        pc_names = protein_clusters_dict.keys()
        num_pcs = len(pc_names)
        for i in range(0, num_pcs):
            self.progress.update('%d of %d' % (i, num_pcs)) if i % 10 == 0 else None
            pc_name = pc_names[i]

            if len(protein_clusters_dict[pc_name]) == 1:
                # this sequence is a singleton and does not need alignment
                continue

            gene_sequences_in_pc = []
            for gene_entry in protein_clusters_dict[pc_name]:
                sequence = self.genomes_storage.get_gene_sequence(gene_entry['genome_name'], gene_entry['gene_caller_id'])
                gene_sequences_in_pc.append(('%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id']), sequence),)

            # alignment
            alignments = muscle.run_muscle_stdin(gene_sequences_in_pc)

            for gene_entry in protein_clusters_dict[pc_name]:
                gene_entry['alignment_summary'] = utils.summarize_alignment(alignments['%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id'])])

        self.progress.end()


    def process(self):
        # check sanity
        self.sanity_check()

        # gen pan_db
        self.generate_pan_db()

        # get all protein sequences:
        combined_proteins_FASTA_path = self.get_output_file_path('combined-proteins.fa')
        unique_proteins_FASTA_path, unique_proteins_names_dict = self.genomes_storage.gen_combined_protein_sequences_FASTA(combined_proteins_FASTA_path, exclude_partial_gene_calls=self.exclude_partial_gene_calls)

        # run search
        blastall_results = self.run_search(unique_proteins_FASTA_path, unique_proteins_names_dict)

        # generate MCL input from filtered blastall_results
        mcl_input_file_path = self.gen_mcl_input(blastall_results)

        # get clusters from MCL
        mcl_clusters = self.run_mcl(mcl_input_file_path)

        # we have the raw protein clusters dict, but we need to re-format it for following steps
        protein_clusters_dict = self.gen_protein_clusters_dict_from_mcl_clusters(mcl_clusters)
        del mcl_clusters

        # compute alignments for genes within each PC (or don't)
        self.compute_alignments_for_PCs(protein_clusters_dict)

        # store protein clusters dict into the db
        self.store_protein_clusters(protein_clusters_dict)

        # populate the pan db with results
        self.process_protein_clusters(protein_clusters_dict)

        # gen samples info and order files
        self.gen_samples_db()

        # done
        self.run.info('log file', self.run.log_file_path)
        self.run.quit()
