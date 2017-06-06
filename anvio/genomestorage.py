# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module to dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sene of internal or external genome descriptions is also welcome.
"""

import os
import hashlib

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError


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

        self.internal_genome_names = list(internal_genomes_dict.keys())
        self.external_genome_names = list(external_genomes_dict.keys())

        if not self.internal_genome_names and not self.external_genome_names:
            raise ConfigError("You in fact tried to create a genomes storage file without providing any internal or external genome\
                                descriptions! You got 5 anvi'o points for being awesome, but this is not gonna work since you really\
                                need to provide at least one of those descriptions :/")

        if len(self.internal_genome_names) + len(self.external_genome_names) != len(set(self.internal_genome_names + self.external_genome_names)):
            raise ConfigError("Each entry both in internal and external genome descriptions should have a unique 'name'. This does not\
                                seem to be the case with your input :/")

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
        genomes_with_no_functional_annotation = []
        function_annotation_sources_per_genome = {}
        all_function_annotation_sources_observed = set([])
        for genome_name in self.genomes:
            g = self.genomes[genome_name]
            contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
            sources = contigs_db.meta['gene_function_sources']
            contigs_db.disconnect()

            if not sources:
                genomes_with_no_functional_annotation.append(genome_name)
            else:
                function_annotation_sources_per_genome[genome_name] = sources
                all_function_annotation_sources_observed.update(sources)

        if genomes_with_no_functional_annotation:
            if len(genomes_with_no_functional_annotation) == len(self.genomes):
                self.run.warning("None of your genomes seem to have any functional annotation. No biggie. Things will continue to work. But\
                                  then your genomes have no functional annotation. SAD.")
            else:
                self.run.warning("Some of your genomes (%d of the %d, to be precise) seem to have no functional annotation. Since this workflow\
                                  can only use matching functional annotations across all genomes involved, having even one genome without\
                                  any functions means that there will be no matching function across all. Things will continue to work, but\
                                  you will have no functions at the end for your protein clusters." % \
                                                (len(genomes_with_no_functional_annotation), len(self.genomes)))

            # make sure it is clear.
            function_annotation_sources_per_genome = {}
            all_function_annotation_sources_observed = set([])
        elif not len(all_function_annotation_sources_observed):
            self.run.warning("None of your genomes seem to have any functional annotation. No biggie. Things will continue to work. But\
                              then your genomes have no functional annotation. It is sad.")
        else:
            # this guy down below fills in the self.function_annotation_sources with function annotation sources
            # that are common to all genomes.
            for sources in list(function_annotation_sources_per_genome.values()):
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

        if not self.storage_path:
            raise ConfigError("Anvi'o genomes storage is speaking. Someone called the init function,\
                               yet there is nothing to initialize since genome storage path variable\
                               (args.genomes_storage) is None. If you are an end user, please make sure\
                               you provide the genomes storage paramater to whatever program you were\
                               running. If you are a developer, you probably already figured what is\
                               wrong. If you are a cat, you need to send us an e-mail immediately.")

        # let's take care of the genome names to focus, if there are any, first. 
        if self.genome_names_to_focus:
            if filesnpaths.is_file_exists(self.genome_names_to_focus, dont_raise=True):
                self.genome_names_to_focus = utils.get_column_data_from_TAB_delim_file(self.genome_names_to_focus, column_indices=[0], expected_number_of_fields=1)[0]
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
            raise ConfigError('Genes must have been called during the generation of contigs database for this workflow to work. However,\
                                these external genomes do not have gene calls: %s' % (', '.join(genomes_missing_gene_calls)))

        # if two contigs db has the same hash, we are kinda f'd:
        if len(set([self.genomes[genome_name]['genome_hash'] for genome_name in self.external_genome_names])) != len(self.external_genome_names):
            raise ConfigError('Not all hash values are unique across all contig databases you provided. Something\
                                very fishy is going on :/')


        if len(set([self.genomes[genome_name]['genome_hash'] for genome_name in self.internal_genome_names])) != len(self.internal_genome_names):
            raise ConfigError("Not all hash values are unique across internal genomes. This is almost impossible to happen unless something very\
                                wrong with your workflow :/ Please let the developers know if you can't figure this one out")

        # make sure HMMs for SCGs were run for every contigs db:
        genomes_missing_hmms_for_scgs =  [g for g in self.genomes if not self.genomes[g]['hmms_for_scgs_were_run']]
        if len(genomes_missing_hmms_for_scgs):
            if len(genomes_missing_hmms_for_scgs) == len(self.genomes):
                raise ConfigError("The contigs databases you are using for this analysis are missing HMMs for single-copy core genes. In other words,\
                                    you don't seem to have run `anvi-run-hmms` on them. Although it is perfectly legal to have anvi'o contigs databases\
                                    without HMMs run on SCGs, the current pangenomic workflow does not want to deal with this :( Sorry!")
            else:
                raise ConfigError("Some of the genomes you have for this analysis are missing HMM hits for SCGs (%d of %d of them, to be precise). You\
                                    can run `anvi-run-hmms` on them to recover from this. Here is the list: %s" % \
                                                    (len(genomes_missing_hmms_for_scgs), len(self.genomes), ','.join(genomes_missing_hmms_for_scgs)))

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
                raise ConfigError("The genomes storage file must end with '-GENOMES.h5'. Anvi'o developers do know how ridiculous\
                                    this requirement sounds like, but if you have seen the things they did, you would totally\
                                    understand why this is necessary.")

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
        storage_hash = hashlib.sha224('_'.join(self.genomes[genome_name]['genome_hash'] for genome_name in self.genomes).encode('utf-8')).hexdigest()[0:8]

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
        genome_hash = hashlib.sha224('_'.join([''.join(split_names_of_interest), contigs_db.meta['contigs_db_hash']]).encode('utf-8')).hexdigest()[0:12]
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
            raise ConfigError("There are 0 splits defined for bin id %s in collection %s..." % (entry['bin_id'], entry['collection_id']))

        return split_names_of_interest


