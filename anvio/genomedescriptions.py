# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module for dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sense of internal or external genome descriptions is also welcome.
"""

import os
import sys
import copy
import hashlib
import argparse

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.ccollections as ccollections
import anvio.filesnpaths as filesnpaths

import pandas as pd

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class GenomeDescriptions(object):
    def __init__(self, args=None, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        self.genomes = {}
        self.internal_genomes_dict = None
        self.external_genomes_dict = None

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.just_do_it = A('just_do_it')
        self.functions_are_available = False
        self.function_annotation_sources = set([])

        self.input_file_for_internal_genomes = A('internal_genomes')
        self.input_file_for_external_genomes = A('external_genomes')
        self.skip_checking_genome_hashes = A('skip_checking_genome_hashes')
        self.list_hmm_sources = A('list_hmm_sources')          # <<< these two are look out of place, but if the args requests
        self.list_available_gene_names = A('list_available_gene_names') #     for information about HMMs, this is the bets place to set them
        self.gene_caller = A('gene_caller')

        if self.input_file_for_internal_genomes or self.input_file_for_external_genomes:
            self.read_genome_paths_from_input_files()

        # see `self.is_proper_db` function for these two variables:
        self.contigs_dbs_found = set([])
        self.profile_dbs_found = set([])

        self.external_genomes_with_identical_hashes = {}
        self.internal_genomes_with_identical_hashes = {}


    def is_proper_db(self, db_path, db_type):
        """Check if contigs db or profile db is OK.

        A given contigs database have multiple entries in an internal or external genomes file.
        The same goes for a profile database. One of the things we have to do during initialization
        is to check whether each entry in int/ext genomes files is associated with a legitimate
        contigs or profile databases. This function fills up `self.contigs_dbs_found` and
        `self.profile_dbs_found` variables every time it is called. If there is alread an entry
        it returns True without any additional function call.

        Parameters
        ==========
            db_path: str
                path to the database
            db_type: str in ['contigs', 'profile']
                whether a given database is assumed to be a contigs or profile database

        """
        db_types_factory = {'contigs': {'check': utils.is_contigs_db, 'variable': self.contigs_dbs_found},
                            'profile': {'check': utils.is_profile_db, 'variable': self.profile_dbs_found}}

        if db_type not in db_types_factory:
            raise ConfigError("is_proper_db :: wrong `db_type` :/ Pick either: %s" % ', '.join(db_types_factory))

        if db_path in db_types_factory[db_type]['variable']:
            return True
        else:
            db_types_factory[db_type]['check'](db_path)
            db_types_factory[db_type]['variable'].add(db_path)


    def names_check(self):
        i, n = list(self.internal_genomes_dict.keys() if self.internal_genomes_dict else []), \
               list(self.external_genomes_dict.keys() if self.external_genomes_dict else [])

        if not i and not n:
            raise ConfigError("You actually managed to get all the way down here in the code without actually providing any internal "
                              "or external genome files! You got 5 anvi'o points for being awesome. But this is not gonna work since "
                              "you really need to provide at least one of those files so anvi'o takes away 4 of those points :/ The "
                              "anvi'o giveth, and the anvi'o taketh away. Enjoy your point.")

        if len(i) + len(n) != len(set(i + n)):
            raise ConfigError("Each entry both in internal and external genome descriptions should have a unique 'name'. This does not "
                               "seem to be the case with your input :/")


    def read_genome_paths_from_input_files(self):
        """Reads internal and external genome files, populates self.genomes"""

        fields_for_internal_genomes_input = ['name', 'bin_id', 'collection_id', 'profile_db_path', 'contigs_db_path']
        fields_for_external_genomes_input = ['name', 'contigs_db_path']

        self.internal_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(self.input_file_for_internal_genomes, expected_fields=fields_for_internal_genomes_input) if self.input_file_for_internal_genomes else {}
        self.external_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(self.input_file_for_external_genomes, expected_fields=fields_for_external_genomes_input) if self.input_file_for_external_genomes else {}


    def list_HMM_info_and_quit(self):
        hmm_sources_in_all_genomes = self.get_HMM_sources_common_to_all_genomes()

        # since we know hmm sources in `hmm_sources_in_all_genomes` are common to all genomes,
        # we could use any of those genomes to learn about the specifics of them. here we take
        # the first one from `self.genomes`
        contigs_db = dbops.ContigsDatabase(list(self.genomes.values())[0]['contigs_db_path'])
        hmm_sources_info = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
        contigs_db.disconnect()

        if self.list_hmm_sources or self.list_available_gene_names:
            if not len(hmm_sources_in_all_genomes):
                raise ConfigError("There are no HMM sources among your external genomes that occur in every genome :/")


        if self.list_hmm_sources:
            self.run.warning(None, 'HMM SOURCES COMMON TO ALL %d GENOMES' % (len(self.genomes)), lc='yellow')
            for source in hmm_sources_in_all_genomes:
                s = hmm_sources_info[source]
                self.run.info_single('%s [type: %s] [num genes: %d]' % (source, s['search_type'], len(s['genes'])))
            sys.exit(0)

        if self.list_available_gene_names:
            self.run.warning(None, 'GENES IN HMM SOURCES COMMON TO ALL %d GENOMES' % (len(self.genomes)), lc='yellow')
            for source in hmm_sources_in_all_genomes:
                s = hmm_sources_info[source]
                gene_names = ', '.join(sorted([g.strip() for g in s['genes'].split(',')]))
                self.run.info_single('%s [type: %s]: %s' % (source, s['search_type'], gene_names), nl_after = 2)
            sys.exit(0)


    def get_HMM_sources_common_to_all_genomes(self):
        """Returns True if all HMM sources in all genomes are comparable"""

        hmm_sources_info_per_genome = {}

        # first recover hmm sources info per genome
        for genome_name in self.genomes:
            if 'hmm_sources_info' not in self.genomes[genome_name]:
                # someone did not run the expensive `init` function. but we can recover this
                # here quitte cheaply
                contigs_db = dbops.ContigsDatabase(self.genomes[genome_name]['contigs_db_path'])
                hmm_sources_info = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
            else:
                hmm_sources_info = self.genomes[genome_name]['hmm_sources_info']

            hmm_sources_info_per_genome[genome_name] = hmm_sources_info

        hmm_sources_found = set([])
        for genome_name in self.genomes:
            [hmm_sources_found.add(s) for s in hmm_sources_info.keys()]

        # find out hmm_sources that occur in all genomes
        hmm_sources_in_all_genomes = copy.deepcopy(hmm_sources_found)
        for genome_name in self.genomes:
            for hmm_source in hmm_sources_found:
                if hmm_source not in hmm_sources_info_per_genome[genome_name] and hmm_source in hmm_sources_in_all_genomes:
                    hmm_sources_in_all_genomes.remove(hmm_source)

        return hmm_sources_in_all_genomes


    def load_genomes_descriptions(self, skip_functions=False, init=True):
        """Load genome descriptions from int/ext genome dictionaries"""

        # start with a sanity check to make sure name are distinct
        self.names_check()

        self.internal_genome_names = list(self.internal_genomes_dict.keys())
        self.external_genome_names = list(self.external_genomes_dict.keys())

        # let us know if the user did not want a full init.
        self.full_init = init

        # convert relative paths to absolute paths and MERGE internal and external genomes into self.genomes:
        for source, input_file in [(self.external_genomes_dict, self.input_file_for_external_genomes),
                                   (self.internal_genomes_dict, self.input_file_for_internal_genomes)]:
            for genome_name in source:
                self.genomes[genome_name] = source[genome_name]
                for db_path_var in ['contigs_db_path', 'profile_db_path']:
                    if db_path_var not in self.genomes[genome_name]:
                        continue
                    path = self.genomes[genome_name][db_path_var]

                    if not path:
                        raise ConfigError("Bad news: anvi'o was loading genome desriptions, and it run into an empty path for "
                                          "the genome %s. How did this happen? HOW? :(" % genome_name)

                    if not path.startswith('/'):
                        self.genomes[genome_name][db_path_var] = os.path.abspath(os.path.join(os.path.dirname(input_file), path))

                # while we are going through all genomes and reconstructing self.genomes for the first time,
                # let's add the 'name' attribute in it as well.'
                self.genomes[genome_name]['name'] = genome_name

        # add hashes for each genome in the self.genomes dict.
        self.genome_hash_to_genome_name = {}

        self.progress.new('Setting up genome hash dicts', progress_total_items=len(self.genomes))
        for genome_name in self.external_genome_names:
            self.progress.update("working on %s (external)" % (genome_name), increment=True)
            g_hash = str(self.get_genome_hash_for_external_genome(self.genomes[genome_name]))
            self.genomes[genome_name]['genome_hash'] = g_hash
            self.genome_hash_to_genome_name[g_hash] = genome_name
        for genome_name in self.internal_genome_names:
            self.progress.update("working on %s (internal)" % (genome_name), increment=True)
            g_hash = str(self.get_genome_hash_for_internal_genome(self.genomes[genome_name]))
            self.genomes[genome_name]['genome_hash'] = g_hash
            self.genome_hash_to_genome_name[g_hash] = genome_name
        self.progress.end()

        # if the user wanted anvi'o to not care about checking genome hashes and we ended up
        # finding genomes with identical hashes, let them know
        if self.skip_checking_genome_hashes and (len(self.internal_genomes_with_identical_hashes) or len(self.external_genomes_with_identical_hashes)):
            self.run.warning("While processing internal and/or external genomes files you have provided, "
                             "anvi'o found genomes with identical hashes (which means they were practically "
                             "identical to each other). But since you have instructed anvi'o to ignore that "
                             "it is now continuing with the flow (even %d hashes for your internal genomes and %d) "
                             "hashes for your external gneomes appeared more than once). See below the genome names "
                             "with identical hashes:" % (len(self.internal_genomes_with_identical_hashes),
                                                         len(self.external_genomes_with_identical_hashes)),
                                                         overwrite_verbose=True)

            for _t, _d in [('Internal', self.internal_genomes_with_identical_hashes), ('External', self.external_genomes_with_identical_hashes)]:
                all_genome_hashes = list(_d.keys())
                for genome_hash in all_genome_hashes:
                    self.run.info("%s genomes with hash %s" % (_t, genome_hash), "%s" % ", ".join(_d[genome_hash]),
                                  overwrite_verbose=True,
                                  nl_after = 1 if genome_hash == all_genome_hashes[-1] else 0,
                                  lc='red')

        # if the client is not interested in functions, skip the rest.
        if skip_functions:
            self.functions_are_available = False
        else:
            self.init_functions()

        # this will populate self.genomes with relevant data that can be learned about these genomes such as 'avg_gene_length',
        # 'num_splits', 'num_contigs', 'num_genes', 'percent_redundancy', 'gene_caller_ids', 'total_length', 'partial_gene_calls',
        # 'percent_completion', 'num_genes_per_kb', 'gc_content'.
        if self.full_init:
            self.init_internal_genomes()
            self.init_external_genomes()
        else:
            # init will do everything. but it is very expensive. if the user does not want to
            # init all the bulky stuff, we still can give them the contents of the meta tables.
            for genome_name in self.genomes:
                g = self.genomes[genome_name]
                contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
                for key in contigs_db.meta:
                    g[key] = contigs_db.meta[key]

        # make sure it is OK to go with self.genomes
        self.sanity_check()


    def get_functions_and_sequences_dicts_from_contigs_db(self, genome_name, requested_source_list=None, return_only_functions=False):
        """This function fetches dictionaries of functions, AA sequences, and DNA sequences for a particular genome.

        PARAMETERS
        ==========
        genome_name, str
            the genome name you want data for
        requested_source_list, list
            the functional annotation sources you want data for. If not provided, data will be fetched for all sources in
            self.function_annotation_sources
        return_only_functions, bool
            Return only functions, and don't bother with sequences

        RETURNS
        =======
        function_calls_dict : dictionary of function annotations
        aa_sequences_dict : dictionary of corresponding amino acid sequences
        dna_sequences_dict : dictionary of corresponding nucleotide sequences
        """

        if not requested_source_list:
            requested_source_list = list(self.function_annotation_sources)

        g = self.genomes[genome_name]

        args = argparse.Namespace()
        args.contigs_db = g['contigs_db_path']

        # we are about to initialize the contigs super, but before that, we need to make sure that
        # the class will know about the splits that describe this genome in the contigs database
        # IF it is an internal genome. otherwise we will end up gathering all the functions in
        # gontigs database for it.
        if genome_name in self.internal_genome_names:
            args.split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(g)

        contigs_super = dbops.ContigsSuperclass(args, r=anvio.terminal.Run(verbose=False))

        if self.functions_are_available:
            contigs_super.init_functions(requested_sources=requested_source_list)
            function_calls_dict = contigs_super.gene_function_calls_dict
        else:
            function_calls_dict = {}

        if return_only_functions:
            return (function_calls_dict, None, None)

        # get dna sequences
        gene_caller_ids_list, dna_sequences_dict = contigs_super.get_sequences_for_gene_callers_ids(gene_caller_ids_list=list(g['gene_caller_ids']))

        # get amino acid sequences.
        # FIXME: this should be done in the contigs super.
        contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
        aa_sequences_dict = contigs_db.db.get_table_as_dict(t.gene_amino_acid_sequences_table_name)
        contigs_db.disconnect()

        return (function_calls_dict, aa_sequences_dict, dna_sequences_dict)


    def init_functions(self):
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
                self.run.warning("None of your genomes seem to have any functional annotation. No biggie. Things will continue to work. But "
                                 "then your genomes have no functional annotation. SAD.")
            else:
                self.run.warning("Some of your genomes (%d of the %d, to be precise) seem to have no functional annotation. Since this workflow "
                                 "can only use matching functional annotations across all genomes involved, having even one genome without "
                                 "any functions means that there will be no matching function across all. Things will continue to work, but "
                                 "you will have no functions at the end for your gene clusters." % \
                                                (len(genomes_with_no_functional_annotation), len(self.genomes)))

            # make sure it is clear.
            function_annotation_sources_per_genome = {}
            all_function_annotation_sources_observed = set([])
        elif not len(all_function_annotation_sources_observed):
            self.run.warning("None of your genomes seem to have any functional annotation. No biggie. Things will continue to work. But "
                             "then your genomes have no functional annotation. It is sad.")
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
                self.run.warning("Although some of your genomes had some functional annotations, none of them were common to all genomes :/ "
                                 "Anvi'o will continue working with them, but you will have no functions available to you downstream. Just "
                                 "so you know, these are the annotation sources observed at least once in at least one of your genomes: '%s'" % \
                                                                    (', '.join(all_function_annotation_sources_observed)))
                self.functions_are_available = False
            else:
                self.functions_are_available = True

                # good. here we know some functions are available, but let's get some further understanding, and report it to the user, you know,
                # because we're nice:
                if len(function_annotation_sources_some_genomes_miss):
                    # some functions were missing from some genomes
                    self.run.warning("Anvi'o has good news and bad news for you (very balanced, as usual). The good news is that there are some "
                                     "functional annotation sources that are common to all of your genomes, and they will be used whenever "
                                     "it will be appropriate. Here they are: '%s'. The bad news is you had more functiona annotation sources, "
                                     "but they were not common to all genomes. Here they are so you can say your goodbyes to them (because "
                                     "they will not be used): '%s'" % \
                                            (', '.join(self.function_annotation_sources), ', '.join(function_annotation_sources_some_genomes_miss)))
                else:
                    # every function ever observed is common to all genomes.
                    self.run.warning("Good news! Anvi'o found all these functions that are common to all of your genomes and will use them for "
                                     "downstream analyses and is very proud of you: '%s'." % (', '.join(self.function_annotation_sources)), lc='green')


    def get_genome_hash_for_external_genome(self, entry):
        self.is_proper_db(entry['contigs_db_path'], db_type='contigs')
        genome_hash = db.DB(entry['contigs_db_path'], None, ignore_version=True).get_meta_value('contigs_db_hash')

        if genome_hash in self.genome_hash_to_genome_name:
            if self.skip_checking_genome_hashes:
                if genome_hash in self.external_genomes_with_identical_hashes:
                    self.external_genomes_with_identical_hashes[genome_hash].add(entry['name'])
                    self.external_genomes_with_identical_hashes[genome_hash].add(self.genome_hash_to_genome_name[genome_hash])
                else:
                    self.external_genomes_with_identical_hashes[genome_hash] = set([self.genome_hash_to_genome_name[genome_hash], entry['name']])
            else:
                self.progress.reset()
                raise ConfigError("While working on your external genomes, anvi'o realized that genome %s and %s seem to have the same hash. "
                                  "If you are aware of this and/or if you would like anvi'o to not check genome hashes, please use the flag "
                                  "`--skip-checking-genome-hashes`." % (self.genome_hash_to_genome_name[genome_hash], entry['name']))

        return genome_hash


    def get_genome_hash_for_internal_genome(self, entry):
        self.is_proper_db(entry['contigs_db_path'], db_type='contigs')
        split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(entry)
        contigs_db_hash = db.DB(entry['contigs_db_path'], None, ignore_version=True).get_meta_value('contigs_db_hash')
        genome_hash = hashlib.sha224('_'.join([''.join(split_names_of_interest), contigs_db_hash]).encode('utf-8')).hexdigest()[0:12]

        if genome_hash in self.genome_hash_to_genome_name:
            if self.skip_checking_genome_hashes:
                if genome_hash in self.internal_genomes_with_identical_hashes:
                    self.internal_genomes_with_identical_hashes[genome_hash].add(entry['name'])
                    self.internal_genomes_with_identical_hashes[genome_hash].add(self.genome_hash_to_genome_name[genome_hash])
                else:
                    self.internal_genomes_with_identical_hashes[genome_hash] = set([self.genome_hash_to_genome_name[genome_hash], entry['name']])
            else:
                self.progress.reset()
                genome_1, genome_2 = self.genome_hash_to_genome_name[genome_hash], entry['name']
                raise ConfigError("According to hash values anvi'o has been generating for your internal genomes, not all genomes you have seem to be uniuqe. "
                                  "It is most likely you unintentionally listed the same information for different genome names. If you would like "
                                  "to double check, genome %s (in '%s') and genome %s (in '%s') seem to have the same hash (so they are basically the same genomes). "
                                  "If you are aware of this and/or if you would like anvi'o to not check genome hashes, please use the flag "
                                  "`--skip-checking-genome-hashes`." % (genome_1,
                                                               self.genomes[genome_1]['collection_id'],
                                                               genome_2,
                                                               self.genomes[genome_2]['collection_id']))

        return genome_hash


    def init_external_genomes(self):
        from anvio.summarizer import ContigSummarizer

        self.progress.new('Initializing external genomes', progress_total_items=len(self.external_genome_names))
        for genome_name in self.external_genome_names:
            c = self.genomes[genome_name]
            c['external_genome'] = True

            self.progress.update('working on %s' % (genome_name), increment=True)

            contigs_db_summary = ContigSummarizer(c['contigs_db_path']).get_contigs_db_info_dict(gene_caller_to_use=self.gene_caller)

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
        from anvio.summarizer import ContigSummarizer

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

                utils.is_profile_db_and_contigs_db_compatible(c['profile_db_path'], c['contigs_db_path'])

                split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(c)

                # here we are using the get_contigs_db_info_dict function WITH split names we found in the collection
                # which returns a partial summary from the contigs database focusing only those splits. a small workaround
                # to be able to use the same function for bins in collections:
                contigs_summary = ContigSummarizer(c['contigs_db_path'])
                summary_from_contigs_db_summary = contigs_summary.get_contigs_db_info_dict(split_names=split_names_of_interest,
                                                                                           gene_caller_to_use=self.gene_caller)

                for key in summary_from_contigs_db_summary:
                    c[key] = summary_from_contigs_db_summary[key]

        self.progress.end()

        self.run.info('Internal genomes', '%d have been initialized.' % len(self.internal_genome_names))


    def get_split_names_of_interest_for_internal_genome(self, entry):
        self.is_proper_db(entry['profile_db_path'], db_type='profile')
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


    def sanity_check(self):
        """Make sure self.genomes is good to go"""

        self.progress.new('Sanity checks')

        # depending on whether args requested such behavior.
        self.progress.update("...")
        self.list_HMM_info_and_quit()

        # make sure genes are called in every contigs db:
        self.progress.update("Checking gene calls ..")
        genomes_missing_gene_calls = [g for g in self.genomes if not self.genomes[g]['genes_are_called']]
        if len(genomes_missing_gene_calls):
            self.progress.end()
            raise ConfigError('Genes must have been called during the generation of contigs database for this workflow to work. However,\
                                these external genomes do not have gene calls: %s' % (', '.join(genomes_missing_gene_calls)))

        if not self.full_init:
            # if this is not full init, stop the sanity check here.
            self.progress.end()
            self.run.warning("You (or the programmer) requested genome descriptions for your internal and/or external "
                             "genomes to be loaded _without_ a 'full init'. There is nothing for you to be concerned. "
                             "This is just a friendly reminder to make sure you know that if something goes terribly "
                             "wrong later (like your computer sets itself on fire), this may be the reason.")

            return

        self.progress.update("Checking HMMs and SCGs ..")
        # make sure HMMs for SCGs were run for every contigs db:
        genomes_missing_hmms_for_scgs =  [g for g in self.genomes if not self.genomes[g]['hmms_for_scgs_were_run']]
        if len(genomes_missing_hmms_for_scgs):
            if len(genomes_missing_hmms_for_scgs) == len(self.genomes):
                self.progress.reset()
                self.run.warning("The contigs databases you are using for this analysis are missing HMMs for single-copy core genes. "
                                 "Maybe you haven't run `anvi-run-hmms` on your contigs database, or they didn't contain any hits. "
                                 "It is perfectly legal to have anvi'o contigs databases without HMMs or SCGs for things to work, "
                                 "but we wanted to give you heads up so you can have your 'aha' moment if you see funny things in "
                                 "the interface.")
            else:
                self.progress.end()
                raise ConfigError("Some of the genomes you have for this analysis are missing HMM hits for SCGs (%d of %d of them, to be precise). You "
                                   "can run `anvi-run-hmms` on them to recover from this. Here is the list: %s" % \
                                                    (len(genomes_missing_hmms_for_scgs), len(self.genomes), ','.join(genomes_missing_hmms_for_scgs)))

        # make sure genome names are not funny (since they are going to end up being db variables soon)
        self.progress.update("Checking genome names ..")
        [utils.is_this_name_OK_for_database('genome name "%s"' % genome_name, genome_name) for genome_name in self.genomes]

        # figure out whether there are genomes with gene calls that are NOT processed
        self.progress.update("Checking gene calls that are not processed ..")
        genomes_with_non_reported_gene_calls_from_other_gene_callers = []
        for genome_name in self.genomes:
            if self.genomes[genome_name]['gene_calls_from_other_gene_callers']:
                genomes_with_non_reported_gene_calls_from_other_gene_callers.append(genome_name)

        if len(genomes_with_non_reported_gene_calls_from_other_gene_callers):
            info = []
            for genome_name in genomes_with_non_reported_gene_calls_from_other_gene_callers:
                info.append('%s (%s)' % (genome_name,
                                         ', '.join(['%d gene calls by "%s"' % (tpl[1], tpl[0]) for \
                                                         tpl in self.genomes[genome_name]['gene_calls_from_other_gene_callers'].items()])))

            gene_caller = list(self.genomes.values())[0]['gene_caller']
            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning("Some of your genomes had gene calls identified by gene callers other than "
                                 "the gene caller anvi'o used, which was set to '%s' either by default, or because you asked for it. "
                                 "The following genomes contained genes that were not processed (this may be exactly what you expect "
                                 "to happen, but if was not, you may need to use the `--gene-caller` flag to make sure anvi'o is using "
                                 "the gene caller it should be using): %s." % \
                                                (gene_caller, ', '.join(info)), header="PLEASE READ CAREFULLY", lc='green')
            else:
                self.progress.reset()
                self.run.warning("Some of your genomes had gene calls identified by gene callers other than "
                                 "the anvi'o default, '%s', and will not be processed. Use the `--debug` flag "
                                 "if this sounds important and you would like to see more of this message." % \
                                                (gene_caller), header="JUST FYI", lc='green')

        # check whether every genome has at least one gene call.
        self.progress.update("Making sure each genome has at least one gene call ..")
        genomes_with_no_gene_calls = [g for g in self.genomes if not self.genomes[g]['num_genes']]
        if len(genomes_with_no_gene_calls):
            self.progress.reset()
            if len(genomes_with_no_gene_calls) == len(self.genomes):
                raise ConfigError("None of your genomes seem to have a gene call, which is a typical error you get if you are working "
                                  "with contigs databases with external gene calls. You can solve it by looking at the output of the "
                                  "program `anvi-db-info` for a given contigs database in your collection, and use one of the gene "
                                  "caller sources listed in the output using the `--gene-caller` parameter.")
            else:
                raise ConfigError(f"Well, {len(genomes_with_no_gene_calls)} of your {len(self.genomes)} genomes seems to have 0 gene calls. "
                                  f"We can't think of any reason to include genomes that contain no gene calls into a genomes storage, "
                                  f"hence, we are going to stop here and ask you to remove these genomes from your analysis first: "
                                  f"{', '.join(genomes_with_no_gene_calls)}. If you think this is happening because you didn't set "
                                  f"the right source for gene calls, you can always take a look at what is available in a given "
                                  f"contigs database by running the program `anvi-db-info`.")
        self.progress.end()


    def functional_enrichment_stats(self):
        """This function runs Amy Willis's enrichment script to compute functional enrichment in groups of genomes.

        It first prepares an input file describing the occurrence of functions in each group, then feeds this file to
        the script.

        Based on similar function for pangenomes in summarizer.py - Credits to Alon
        """

        # sanity check for R packages
        package_dict = utils.get_required_packages_for_enrichment_test()
        utils.check_R_packages_are_installed(package_dict)

        # output wrangling
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        output_file_path = A('output_file')
        tmp_functional_occurrence_file = filesnpaths.get_temp_file_path()

        enrichment_file_path = output_file_path
        if not enrichment_file_path:
            # if no output was requested it means a programmer is calling this function
            # in that case, we will use a tmp file for the enrichment output
            enrichment_file_path = filesnpaths.get_temp_file_path()
        if filesnpaths.is_file_exists(enrichment_file_path, dont_raise=True):
            raise ConfigError(f"The file {enrichment_file_path} already exists and anvi'o doesn't want to overwrite it "
                               "without express permission. If you think overwriting it is okay, then please `rm` the "
                               "existing file before you run this program again.")

        # this is a hack (credits to Alon) which will make the functional occurrence method print to the tmp file
        # (we need this because the functional occurrence method prints to the args.output_file)
        self.args.output_file = tmp_functional_occurrence_file
        self.run.info("Temporary functional occurrence file: ", tmp_functional_occurrence_file)

        # get functional occurrence input file
        self.functional_occurrence_stats()

        # run enrichment script
        cmd = 'anvi-script-enrichment-stats --input %s --output %s' % (tmp_functional_occurrence_file,
                                                                                      output_file_path)
        log_file = filesnpaths.get_temp_file_path()
        self.run.info("Enrichment script log file", log_file)
        self.progress.new('Functional enrichment analysis')
        self.progress.update("Running Amy's enrichment")
        utils.run_command(cmd, log_file)
        self.progress.end()

        if not filesnpaths.is_file_exists(enrichment_file_path, dont_raise=True):
            raise ConfigError("Something went wrong during the functional enrichment analysis. "
                              f"We don't know what happened, but this log file could contain some clues: {log_file}")

        if filesnpaths.is_file_empty(enrichment_file_path):
            raise ConfigError("Something went wrong during the functional enrichment analysis. "
                              "An output file was created, but it is empty... We don't know why this happened, "
                              f"but this log file could contain some clues: {log_file}")

        # if everything went okay, we remove the log file
        if anvio.DEBUG:
            self.run.warning("Just so you know, because you ran with --debug, we are keeping some temporary files around "
                             f"for you to look at - the enrichment script's log file at {log_file} and the functional "
                             f"occurrence output at {tmp_functional_occurrence_file}. You may want to remove them once you "
                             "are done with them.", lc='green', header="JUST FYI")
        else:
            self.run.warning("The temporary files generated by this program have been removed from your system. If something "
                             "seems strange and you want to take a look at those, please re-run this program with --debug flag.",
                             lc='green', header="JUST FYI")
            os.remove(log_file)
            os.remove(tmp_functional_occurrence_file)

        if not output_file_path:
            # if a programmer called this function then we return a dict
            return utils.get_TAB_delimited_file_as_dictionary(enrichment_file_path)


    def functional_occurrence_stats(self):
        """This function summarizes functional occurrence for groups of genomes.

        If an output file is provided, the functional occurrence dictionary is written to that file.
        Otherwise, the dictionary is returned.

        Based on similar function for pangenomes in summarizer.py - Credits to Alon
        """

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        output_file_path = A('output_file')
        functional_annotation_source = A('annotation_source')
        list_functional_annotation_sources = A('list_annotation_sources')
        functional_occurrence_table_output = A('functional_occurrence_table_output')
        include_ungrouped = A('include_ungrouped')

        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)
        if functional_occurrence_table_output:
            filesnpaths.is_output_file_writable(functional_occurrence_table_output)

        if not self.functions_are_available:
            raise ConfigError("We are sorry to tell you that there are no functions available for your genomes (or, if some of your genomes"
                              "are annotated, there are at least no functional sources common to all your genomes, hence functional "
                              "occurrence cannot be computed. :/")

        if list_functional_annotation_sources:
            self.run.info('Available functional annotation sources', ', '.join(self.function_annotation_sources))
            sys.exit()

        if not functional_annotation_source:
            raise ConfigError("You haven't provided a functional annotation source to make sense of functional "
                              "occurrence stats in your genomes. These are the available annotation sources "
                              f"that are common to all genomes, so pick one: {self.function_annotation_sources}.")

        if functional_annotation_source not in self.function_annotation_sources:
            sources_string = ", ".join(self.function_annotation_sources)
            raise ConfigError(f"Your favorite functional annotation source '{functional_annotation_source}' does not seem to be "
                              "among one of the sources that are available to you. Here are the ones you should choose from: "
                              f"{sources_string}.")

        # get the groups
        genomes_to_groups_dict = {}
        for g in self.genomes:
            if 'group' not in self.genomes[g].keys():
                raise ConfigError("Groups are not defined in at least one of your input files. To rectify this you must ensure that "
                "your input file has a column called 'group'. To help you figure out where this information is missing, here is a genome "
                f" that does not have a group: {g}")
            else:
                genomes_to_groups_dict[g] = self.genomes[g]['group']

        values_that_are_not_none = [s for s in genomes_to_groups_dict.values() if s is not None]
        if not values_that_are_not_none:
            raise ConfigError("Your group column(s) contains only values of type None. "
                              "This is probably a mistake, surely you didn't mean to provide no groups. "
                              "Do you think this is a mistake on our part? Let us know.")

        groups = set([str(genomes_to_groups_dict[g]) for g in genomes_to_groups_dict.keys() if \
                            (genomes_to_groups_dict[g] is not None or include_ungrouped)])

        groups_to_genomes_dict = {}
        for grp in groups:
            groups_to_genomes_dict[grp] = set([genome for genome in genomes_to_groups_dict.keys() if str(genomes_to_groups_dict[genome]) == grp])

        # give the user some info before we continue
        if not anvio.QUIET:
            self.run.info('Functional annotation source', functional_annotation_source)
            self.run.info('Groups', ', '.join(groups))
            self.run.info('Include ungrouped', include_ungrouped)

        self.progress.new('Functional occurrence analysis')
        self.progress.update('Getting functions from database')

        functions_summary_dict = {}
        for g in self.genomes:
            func, aa, dna  = self.get_functions_and_sequences_dicts_from_contigs_db(g, requested_source_list=[functional_annotation_source])
            functions_summary_dict[g] = func

        # get a dictionary of function occurrences per genome
        self.progress.update('Counting functional occurrence')
        function_occurrence_df, function_occurrence_dict = self.get_occurrence_of_functions_in_genomes(functions_summary_dict)

        if functional_occurrence_table_output:
            function_occurrence_df.astype(int).transpose().to_csv(functional_occurrence_table_output, sep='\t')
            if not anvio.QUIET:
                self.progress.reset()
                self.run.info('Occurrence frequency of functions:', functional_occurrence_table_output)

        function_occurrence_table = {}

        # populate the number of genomes per group 'N'
        groups_few_genomes = []
        for grp in groups:
            function_occurrence_table[grp] = {}
            function_occurrence_table[grp]['N'] = len(groups_to_genomes_dict[grp])
            if function_occurrence_table[grp]['N'] < 8:
                groups_few_genomes.append(grp)

        # warn user if they have a low number of genomes per group
        if groups_few_genomes:
            self.progress.reset()
            groups_string = ", ".join(groups_few_genomes)
            self.run.warning("Some of your groups have very few genomes in them, so if you are running functional enrichment, the statistical test may not be very reliable. "
                             "The minimal number of genomes in a group for the test to be reliable depends on a number of factors, "
                             "but we recommend proceeding with great caution because the following groups have fewer than 8 genomes: "
                             f"{groups_string}.")

        self.progress.update("Generating the input table for functional enrichment analysis")
        # get the # of genomes in each group in which function occurs
        function_occurrence_in_groups_df = function_occurrence_df.astype(bool).astype(int)
        # add a group column to the dataframe
        function_occurrence_in_groups_df['group'] = function_occurrence_in_groups_df.index.map(lambda x: str(genomes_to_groups_dict[x]))
        functions_in_groups = function_occurrence_in_groups_df.groupby('group').sum()

        functional_occurrence_summary_dict = {}
        for f in function_occurrence_dict.keys():
            # get proportion of function occurrence in groups 'p'
            for grp in groups:
                function_occurrence_table[grp]['p'] = functions_in_groups.loc[grp, f] / function_occurrence_table[grp]['N']
            function_occurrence_table_df = pd.DataFrame.from_dict(function_occurrence_table, orient='index')

            functional_occurrence_summary_dict[f] = {}
            present_in_genomes = [g for g in self.genomes if function_occurrence_dict[f][g]]
            functional_occurrence_summary_dict[f]["genomes"] = ", ".join(present_in_genomes)
            functional_occurrence_summary_dict[f]["accession"] = ", ".join(list(function_occurrence_dict[f]["accession"]))
            for grp in groups:
                functional_occurrence_summary_dict[f]['p_' + grp] = function_occurrence_table[grp]['p']
                functional_occurrence_summary_dict[f]['N_' + grp] = function_occurrence_table[grp]['N']
            enriched_groups_vector = utils.get_enriched_groups(function_occurrence_table_df['p'].values,
                                                               function_occurrence_table_df['N'].values)
            c_dict = dict(zip(function_occurrence_table_df['p'].index, range(len(function_occurrence_table_df['p'].index))))
            associated_groups = [grp for grp in groups if enriched_groups_vector[c_dict[grp]]]
            functional_occurrence_summary_dict[f]['associated_groups'] = ", ".join(associated_groups)

        if output_file_path:
            self.progress.update('Generating the output file')

            # Sort the columns the way we want them
            columns = [functional_annotation_source, 'accession', 'genomes', 'associated_groups']
            columns.extend([s + grp for s in ['p_', 'N_'] for grp in groups])
            utils.store_dict_as_TAB_delimited_file(functional_occurrence_summary_dict, output_file_path, headers=columns)

        self.progress.end()


    def get_occurrence_of_functions_in_genomes(self, genome_to_func_summary_dict):
        """Here we convert a dictionary of function annotations in each genome to a dictionary of counts per function.

        PARAMETERS
        ==========
        genome_to_func_summary_dict : multi-level dict
            The format of this dictionary is
            (accession, annotation, e_value) = genome_to_func_summary_dict[genome_name][gene_caller_id][annotation_source]

        RETURNS
        =======
        func_occurrence_dict :
            dictionary of function annotation counts in each genome. Its format is
            count of gene calls with function = func_occurrence_dict[function][genome]
            (set of accession numbers with this annotation) = func_occurrence_dict[function]['accession']
        func_occurrence_dataframe : dataframe
            dataframe version of the above dictionary
        """

        func_occurrence_dict = {}
        for g, genes_to_func in genome_to_func_summary_dict.items():
            for gc_id, func_annotations in genes_to_func.items():
                for source, annotation_tuple in func_annotations.items():
                    acc, function, eval = annotation_tuple

                    ## okay, so here is an annoying thing. Currently gene calls with multiple annotations from the same source
                    ## have these annotations separated by '!!!'. If we encounter this, we need to split it into multiple annotations.
                    all_accessions = acc.split('!!!')
                    all_annotations = function.split('!!!')
                    if len(all_accessions) != len(all_annotations):
                        acc_str = ", ".join(all_accessions)
                        annot_str = ", ".join(all_annotations)
                        raise ConfigError(f"Somethin' is drastically wrong here. In your genome {g} we found a gene call (id is {gc_id}) "
                                          f"with multiple annotations, but the number of accession numbers is not equal to the number of "
                                          f"functional annotations. Take a look at the list of accessions: {acc_str} and compare to the list of "
                                          f"annotations: {annot_str}. Do you know what is wrong? If not, please contact the developers for assistance.")

                    for accession, func in zip(all_accessions, all_annotations):
                        # by keying with the function name, we can automatically merge some accessions with same functional annotation. yay us.
                        if func not in func_occurrence_dict:
                            func_occurrence_dict[func] = {}
                            func_occurrence_dict[func]['accession'] = set([accession])
                            for genome_name in self.genomes:
                                func_occurrence_dict[func][genome_name] = 0
                            func_occurrence_dict[func][g] += 1
                        else:
                            # a functional annotation could have multiple accessions
                            # for example, K00844 and K12407 are both hexokinase/glucokinase
                            func_occurrence_dict[func]['accession'].add(accession)
                            func_occurrence_dict[func][g] += 1

        func_occurrence_dataframe = pd.DataFrame.from_dict(func_occurrence_dict)
        func_occurrence_dataframe.drop('accession', inplace=True)

        return func_occurrence_dataframe, func_occurrence_dict

        self.progress.end()


class MetagenomeDescriptions(object):
    def __init__(self, args=None, run=run, progress=progress, enforce_single_profiles=True):
        self.args = args
        self.run = run
        self.progress = progress
        self.enforce_single_profiles = enforce_single_profiles

        self.metagenomes = {}
        self.metagenomes_dict = None
        self.profile_dbs_available = False

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.input_file_for_metagenomes = A('metagenomes')

        if self.input_file_for_metagenomes:
            self.read_paths_from_input_file()


    def names_check(self):
        names = utils.get_column_data_from_TAB_delim_file(self.input_file_for_metagenomes, [0])[0][1:]

        if len(names) != len(set(names)):
            raise ConfigError("Each entry in your metagenomes file must be unique :/")


    def read_paths_from_input_file(self):
        """Reads metagenome files, populates self.metagenomes"""

        columns = utils.get_columns_of_TAB_delim_file(self.input_file_for_metagenomes)

        if 'profile_db_path' in columns:
            fields_for_metagenomes_input = ['name', 'contigs_db_path', 'profile_db_path']
            self.profile_dbs_available = True
        else:
            fields_for_metagenomes_input = ['name', 'contigs_db_path']
            self.profile_dbs_available = False

        self.metagenomes_dict = utils.get_TAB_delimited_file_as_dictionary(self.input_file_for_metagenomes, expected_fields=fields_for_metagenomes_input) if self.input_file_for_metagenomes else {}


    def load_metagenome_descriptions(self, skip_functions=False, init=True):
        """Load metagenome descriptions"""

        # start with a sanity check to make sure name are distinct
        self.names_check()

        self.metagenome_names = list(self.metagenomes_dict.keys())

        for metagenome_name in self.metagenomes_dict:
            self.metagenomes[metagenome_name] = self.metagenomes_dict[metagenome_name]
            for db_path_var in ['contigs_db_path', 'profile_db_path']:
                if db_path_var not in self.metagenomes[metagenome_name]:
                    continue
                path = self.metagenomes[metagenome_name][db_path_var]

                if not path:
                    raise ConfigError("Bad news: anvi'o was loading metagenome desriptions, and it run into an empty path for "
                                      "the metagenome %s. How did this happen? HOW? :(" % metagenome_name)

                if not path.startswith('/'):
                    self.metagenomes[metagenome_name][db_path_var] = os.path.abspath(os.path.join(os.path.dirname(self.input_file_for_metagenomes), path))

            # while we are going through all genomes and reconstructing self.metagenomes for the first time,
            # let's add the 'name' attribute in it as well.'
            self.metagenomes[metagenome_name]['name'] = metagenome_name

        # add hashes for each metagenome in the self.metagenomes dict.
        self.metagenome_hash_to_metagenome_name = {}
        for metagenome_name in self.metagenome_names:
            g_hash = self.get_metagenome_hash(self.metagenomes[metagenome_name])
            self.metagenomes[metagenome_name]['metagenome_hash'] = g_hash
            self.metagenome_hash_to_metagenome_name[g_hash] = metagenome_name

        for metagenome_name in self.metagenomes:
            g = self.metagenomes[metagenome_name]
            contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
            for key in contigs_db.meta:
                g[key] = contigs_db.meta[key]

        # make sure it is OK to go with self.genomes
        self.sanity_check()


    def get_metagenome_hash(self, entry):
        utils.is_contigs_db(entry['contigs_db_path'])
        contigs_db_hash = db.DB(entry['contigs_db_path'], None, ignore_version=True).get_meta_value('contigs_db_hash')

        return contigs_db_hash


    def sanity_check(self):
        """Make sure self.metagenomes is good to go"""

        if self.profile_dbs_available and self.enforce_single_profiles:
            non_single_profiles = [m for m in self.metagenomes if utils.is_profile_db_merged(self.metagenomes[m]['profile_db_path'])
                                                                  and not utils.is_blank_profile(self.metagenomes[m]['profile_db_path'])]

            if len(non_single_profiles):
                raise ConfigError("All profile databases associated with your metagenomes must be single profiles :( Here "
                                  "is a list of them that are not: '%s'." % (', '.join(non_single_profiles)))

        # make sure genes are called in every contigs db:
        metagenomes_missing_gene_calls = [g for g in self.metagenomes if not self.metagenomes[g]['genes_are_called']]
        if len(metagenomes_missing_gene_calls):
            raise ConfigError('Genes must have been called during the generation of contigs database for this workflow to work. However,\
                               these metagenomes do not have gene calls: %s' % (', '.join(metagenomes_missing_gene_calls)))

        # if two contigs db has the same hash, we are kinda f'd:
        if len(set([self.metagenomes[metagenome_name]['metagenome_hash'] for metagenome_name in self.metagenome_names])) != len(self.metagenome_names):
            raise ConfigError('Not all hash values are unique across all contig databases you provided. Something '
                               'very fishy is going on :/')

        # make sure genome names are not funny (since they are going to end up being db variables soon)
        [utils.is_this_name_OK_for_database('metagenome name "%s"' % metagenome_name, metagenome_name) for metagenome_name in self.metagenomes]


class AggregateFunctions:
    """Get functions from anywhere.

    The purpose of this class is to collect functions from many distinct databases,
    including external genomes, internal genomes, and genomes storage, and report a
    single set of disctionaries that give access to the presence/absence and frequency
    of all functions annotated by a single source.
    """

    def __init__(self, args, skip_sanity_check=False, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.genomes_storage_path = A('genomes_storage')
        self.function_annotation_source = A('annotation_source')
        self.external_genomes_path = A('external_genomes')
        self.internal_genomes_path = A('internal_genomes')
        self.min_occurrence = A('min_occurrence') or 1
        self.aggregate_based_on_accession = A('aggregate_based_on_accession')
        self.aggregate_using_best_hit = A('aggregate_using_best_hit')

        # these are some primary data structures this class reports
        self.accession_to_function_name_dict = {}
        self.accessions_to_genomes_dict_frequency = {}
        self.accessions_to_genomes_dict_presence_absence = {}
        self.genome_names_considered_for_functional_mode = set({})

        self.sanity_checked = False
        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        if not self.function_annotation_source:
            raise ConfigError("When you think about it, this mode can be useful only if someone requests a "
                              "an annotation source to be used for aggregating all the information from all "
                              "the genomes. Someoen didn't specify any function annotation source :/")

        if self.min_occurrence and not isinstance(self.min_occurrence, int):
            raise ConfigError("Obviously, --min-occurrence must be an integer.")

        if self.min_occurrence < 1:
            raise ConfigError(f"What do you have in mind when you say I want my functions to occur in at least {self.min_occurrence} genomes?")

        self.sanity_checked = True


    def update_accession_dicts(self, genome_name, accession, function):
        """Modify accession dicts with new things.

        This function is necessary to avoid redundant code to handle function dicts of different kinds
        we get from int/external genomes and genome storage. a redesign of the genome storage will likely
        fix this problem in the future by unifying how functions are collected from different anvi'o databases.
        """

        if genome_name not in self.genome_names_considered_for_functional_mode:
            self.genome_names_considered_for_functional_mode.add(genome_name)

        if not accession or not len(accession):
            raise ConfigError(f"Anvi'o is very sorry to tell you that the function annotation source you "
                              f"have chosen for this, '{self.function_annotation_source}', seem to incluce "
                              f"functions with no accession IDs. Here is one example function with no "
                              f"accession id: '{function}'. You will have to choose another function "
                              f"annotation source :(")

        if not function:
            raise ConfigError(f"It saddens anvi'o to let you know that there are some function names in "
                              f"'{self.function_annotation_source}' that clearly are blank. Here is an "
                              f"example accession ID that has a blank function name: '{accession}'. You "
                              f"will need to choose another function annotation source, or someohow fix "
                              f"this by using a combination of `anvi-export-functions` and "
                              f"`anvi-import-functions` :(")

        if self.aggregate_using_best_hit:
            accession = accession.split('!!!')[0]
            function = function.split('!!!')[0]
        else:
            pass

        if self.aggregate_based_on_accession:
            pass
        else:
            function, accession = accession, function

        if accession not in self.accessions_to_genomes_dict_frequency:
            self.accessions_to_genomes_dict_frequency[accession] = Counter({})
            self.accessions_to_genomes_dict_presence_absence[accession] = Counter({})


        self.accessions_to_genomes_dict_frequency[accession][genome_name] += 1
        self.accessions_to_genomes_dict_presence_absence[accession][genome_name] = 1

        if accession not in self.accession_to_function_name_dict:
            self.accession_to_function_name_dict[accession] = {self.function_annotation_source: function}

        return


    def _init_functions_from_int_ext_genomes(self):
        if not self.external_genomes_path and not self.internal_genomes_path:
            return

        g = GenomeDescriptions(self.args, run=terminal.Run(verbose=False))
        g.load_genomes_descriptions()
        g.init_functions()

        for genome_name in g.genomes:
            gene_functions_in_genome_dict, _, _= g.get_functions_and_sequences_dicts_from_contigs_db(genome_name, requested_source_list=[self.function_annotation_source], return_only_functions=True)
            # reminder, an entry in gene_functions_in_genome_dict looks like this:
            # 2985: {'COG20_PATHWAY': ('COG0073!!!COG0143', 'Aminoacyl-tRNA synthetases', 0)}
            for entry in gene_functions_in_genome_dict.values():
                accession, function, e_value = entry[self.function_annotation_source]

                self.update_accession_dicts(genome_name, accession, function)


    def _init_functions_from_genomes_storage(self):
        if not self.genomes_storage_path:
            return

        from anvio.genomestorage import GenomeStorage
        g = GenomeStorage(storage_path=self.genomes_storage_path, function_annotation_sources=[self.function_annotation_source], run=terminal.Run(verbose=False), progress=self.progress, skip_init=True)

        # make sure we are not overwriting existing genome names in int or ext genomes:
        genome_names_in_storage_db = g.db.get_single_column_from_table(t.genome_info_table_name, 'genome_name', unique=True)
        already_in_the_dict = [g for g in genome_names_in_storage_db if g in self.genome_names_considered_for_functional_mode]
        if len(already_in_the_dict):
            raise ConfigError(f"Anvi'o is not happy because there are some genome names that occur both in the "
                              f"genome storage and among those that are specified through internal or external "
                              f"genomes files. Here they are: {', '.join(already_in_the_dict)}.")

        gene_functions_in_genomes_dict, _ = g.get_gene_functions_in_genomes_dict()
        for entry in gene_functions_in_genomes_dict.values():
            # an entry in gene_functions_in_genomes_dict looks lke this:
            # 72645: {'genome_name': 'B_lactis_BF052', 'gene_callers_id': 443, 'source': 'COG20_PATHWAY', 'accession': 'COG0207', 'function': 'Thymidylate biosynthesis', 'e_value': 4.2e-149}
            genome_name, accession, function = entry['genome_name'], entry['accession'], entry['function']

            self.update_accession_dicts(genome_name, accession, function)


    def init(self):
        # learn functions, generate `self.views`. First start, with internal and
        # external genomes.
        # here we will update the same two dictionaries with the informaiton in the genome
        # storage
        self._init_functions_from_int_ext_genomes()
        self._init_functions_from_genomes_storage()

        if self.min_occurrence:
            num_occurrence_of_accessions = [(c, sum(self.accessions_to_genomes_dict_presence_absence[c].values())) for c in self.accessions_to_genomes_dict_presence_absence]
            accessions_to_remove = [accession for (accession, frequency) in num_occurrence_of_accessions if frequency < self.min_occurrence]

            if len(accessions_to_remove):
                for accession in accessions_to_remove:
                    self.accession_to_function_name_dict.pop(accession)
                    self.accessions_to_genomes_dict_frequency.pop(accession)
                    self.accessions_to_genomes_dict_presence_absence.pop(accession)

                self.run.warning(f"As per your request, anvi'o removed {len(accessions_to_remove)} accession IDs "
                                 f"of {self.function_annotation_source} from downstream analyses sicne they occurred "
                                 f"in less than {self.min_occurrence} genomes.")
