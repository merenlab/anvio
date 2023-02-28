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
P = terminal.pluralize

class GenomeDescriptions(object):
    def __init__(self, args=None, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        self.genomes = {}
        self.internal_genomes_dict = None
        self.external_genomes_dict = None
        self.initialized = False

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.just_do_it = A('just_do_it')
        self.functions_are_available = False
        self.function_annotation_sources = set([])
        self.function_annotation_sources_some_genomes_miss = set([])

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
        if not self.initialized:
            self.load_genomes_descriptions()

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
            for source in sorted(hmm_sources_in_all_genomes):
                s = hmm_sources_info[source]
                self.run.info_single('%s [type: %s] [num genes: %d]' % (source, s['search_type'], len(s['genes'])))
            sys.exit(0)

        if self.list_available_gene_names:
            self.run.warning(None, 'GENES IN HMM SOURCES COMMON TO ALL %d GENOMES' % (len(self.genomes)), lc='yellow')
            for source in sorted(hmm_sources_in_all_genomes):
                s = hmm_sources_info[source]
                gene_names = ', '.join(sorted([g.strip() for g in s['genes'].split(',')]))
                self.run.info_single('%s [type: %s]: %s' % (source, s['search_type'], gene_names), nl_after = 2)
            sys.exit(0)


    def get_HMM_sources_common_to_all_genomes(self, sources_that_must_be_common=None):
        """Returns True if all HMM sources in all genomes are comparable.

        If you send a list of sources via the parameter `sources_sources_that_must_be_common`, it will also
        ensure that they are common to all genomes, or will throw an exception.
        """

        self.progress.new('Identifying HMM sources common to all genomes', progress_total_items=len(self.genomes))

        hmm_sources_info_per_genome = {}

        # first recover hmm sources info per genome
        for genome_name in self.genomes:
            self.progress.update(f"Processing '{genome_name}' ...", increment=True)
            if 'hmm_sources_info' not in self.genomes[genome_name]:
                # someone did not run the expensive `init` function. but we can recover this
                # here quitte cheaply
                contigs_db = dbops.ContigsDatabase(self.genomes[genome_name]['contigs_db_path'])
                hmm_sources_info = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
            else:
                hmm_sources_info = self.genomes[genome_name]['hmm_sources_info']

            hmm_sources_info_per_genome[genome_name] = hmm_sources_info

        self.progress.update('Finalizing ...')
        hmm_sources_found = set([])
        for genome_name in self.genomes:
            [hmm_sources_found.add(s) for s in hmm_sources_info.keys()]

        # find out hmm_sources that occur in all genomes
        hmm_sources_in_all_genomes = copy.deepcopy(hmm_sources_found)
        for genome_name in self.genomes:
            for hmm_source in hmm_sources_found:
                if hmm_source not in hmm_sources_info_per_genome[genome_name] and hmm_source in hmm_sources_in_all_genomes:
                    hmm_sources_in_all_genomes.remove(hmm_source)

        self.progress.end()

        if sources_that_must_be_common:
            hmm_sources_missing = [s for s in sources_that_must_be_common if s not in hmm_sources_in_all_genomes]
            if len(hmm_sources_missing):
                genomes_missing_any = set([])

                for genome_name in self.genomes:
                    for hmm_source in sources_that_must_be_common:
                        if hmm_source not in hmm_sources_info_per_genome[genome_name].keys():
                            genomes_missing_any.add(genome_name)

                raise ConfigError(f"Bad news. {P('An HMM source', len(hmm_sources_missing), alt='Some HMM sources')} "
                                  f"you really need ({', '.join(hmm_sources_missing)}) {P('is', len(hmm_sources_missing), alt='are')} "
                                  f"missing from some of your contigs databases :/ Here is the list: [{', '.join(genomes_missing_any)}].")


        return hmm_sources_in_all_genomes


    def load_genomes_descriptions(self, skip_functions=False, init=True, skip_sanity_check=False):
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
            self.progress.new('Initializing meta information for genomes', progress_total_items=len(self.genomes))
            self.progress.update('...')
            
            for genome_name in self.genomes:
                self.progress.update(f"Working on '{genome_name}' ...", increment=True)
                g = self.genomes[genome_name]
                contigs_db = dbops.ContigsDatabase(g['contigs_db_path'])
                for key in contigs_db.meta:
                    g[key] = contigs_db.meta[key]
            self.progress.end()

        # we are done hre.
        self.initialized = True

        # make sure it is OK to go with self.genomes
        if not skip_sanity_check:
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

            self.function_annotation_sources_some_genomes_miss = all_function_annotation_sources_observed.difference(self.function_annotation_sources)

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
                if len(self.function_annotation_sources_some_genomes_miss):
                    # some functions were missing from some genomes
                    self.run.warning("Anvi'o has good news and bad news for you (very balanced, as usual). The good news is that there are some "
                                     "functional annotation sources that are common to all of your genomes, and they will be used whenever "
                                     "it will be appropriate. Here they are: '%s'. The bad news is you had more function annotation sources, "
                                     "but they were not common to all genomes. Here they are so you can say your goodbyes to them (because "
                                     "they will not be used): '%s'" % \
                                            (', '.join(self.function_annotation_sources), ', '.join(self.function_annotation_sources_some_genomes_miss)))
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

        # depending on whether args requested such behavior.
        if self.list_hmm_sources:
            self.list_HMM_info_and_quit()

        self.progress.new('Sanity checks')

        # make sure genes are called in every contigs db:
        self.progress.update("Checking gene calls ..")
        genomes_missing_gene_calls = [g for g in self.genomes if not self.genomes[g]['genes_are_called']]
        if len(genomes_missing_gene_calls):
            self.progress.end()
            raise ConfigError('Genes must have been called during the generation of contigs database for this workflow to work. However,\
                                these external genomes do not have gene calls: %s' % (', '.join(genomes_missing_gene_calls)))

        # make sure genome names are not funny (since they are going to end up being db variables soon)
        self.progress.update("Checking genome names ..")
        [utils.is_this_name_OK_for_database('genome name "%s"' % genome_name, genome_name) for genome_name in self.genomes]

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
                self.progress.reset()
                self.run.warning(f"Some of the genomes you have for this analysis are missing HMM hits for SCGs ({len(genomes_missing_hmms_for_scgs)} "
                                 f"of {len(self.genomes)} of them, to be precise). This may or may not be critical for your downstream analyses, but "
                                 f"if you would like to be certain or run into a bigger problem with the later steps of your analysis, you can always "
                                 f"run `anvi-run-hmms` on your contigs databases to call your single-copy core genes. For the sake of brevity, here "
                                 f"is the list of genomes missing SCGs: {','.join(genomes_missing_hmms_for_scgs)}.", header="SCGs ARE NOT CALLED", lc="yellow")


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


    def load_metagenome_descriptions(self, skip_functions=False, init=True, skip_sanity_check=False):
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

        if not skip_sanity_check:
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


    def get_functions_and_sequences_dicts_from_contigs_db(self, metagenome_name, requested_source_list=None, return_only_functions=False):
        """This function fetches dictionaries of functions, AA sequences, and DNA sequences for a particular metagenome.

        PARAMETERS
        ==========
        metagenome_name, str
            the metagenome name you want data for
        requested_source_list, list
            the functional annotation sources you want data for.
        return_only_functions, bool
            Return only functions, and don't bother with sequences

        RETURNS
        =======
        function_calls_dict : dictionary of function annotations
        aa_sequences_dict : dictionary of corresponding amino acid sequences
        dna_sequences_dict : dictionary of corresponding nucleotide sequences
        """

        if not requested_source_list:
            raise ConfigError("Someone is trying to call the `get_functions_and_sequences_dicts_from_contigs_db()` from the MetagenomeDescriptions "
                              "class without providing a list of annotation sources to use. This won't work unfortunately.")

        g = self.metagenomes[metagenome_name]

        if 'gene_function_sources' not in g:
            raise ConfigError(f"Uh oh. The metagenome you are trying to load functions from ({metagenome_name}) does not appear to have "
                              f"functional annotation sources initialized.")
        g_sources = g['gene_function_sources']
        srcs_missing_from_g = []
        for s in requested_source_list:
            if s not in g_sources:
                srcs_missing_from_g.append(s)
        if srcs_missing_from_g:
            req_str = ", ".join(requested_source_list)
            avail_str = ", ".join(g_sources)
            miss_str = ", ".join(srcs_missing_from_g)
            raise ConfigError(f"Whoops. Some of the functional annotation sources you requested for metagenome {metagenome_name} are "
                              f"not in its contigs database. These are the sources that you asked for: {req_str}. And these are the "
                              f"sources that are AVAILABLE for this metagenome: {avail_str}. So, the following sources are MISSING: "
                              f"{miss_str}")

        args = argparse.Namespace()
        args.contigs_db = g['contigs_db_path']

        contigs_super = dbops.ContigsSuperclass(args, r=anvio.terminal.Run(verbose=False))
        contigs_super.init_functions(requested_sources=requested_source_list)
        function_calls_dict = contigs_super.gene_function_calls_dict

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


class AggregateGenomes(object):
    """Aggregate information related to a group of genomes from anywhere.
    The purpose of this class is to collect all relevant information about
    a group of genomes comprehensively for downstream analyses. The primary
    client of this class is genome view.
    Please note that despite similar names, the use of this class and the class
    AggregateFunctions is quite different.
    """

    def __init__(self, args, run=run, progress=progress, skip_init=False):
        self.run = run
        self.progress = progress

        self.args = args
        self.initialized = False

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.external_genomes_file = A('external_genomes')
        self.internal_genomes_file = A('internal_genomes')
        self.pan_db_path = A('pan_db')

        if self.pan_db_path:
            utils.is_pan_db(self.pan_db_path)


        self.genome_descriptions = GenomeDescriptions(args, progress=terminal.Progress(verbose=False))
        self.genome_descriptions.load_genomes_descriptions()

        # critical items for genome view bottleroutes:
        self.genomes = {}
        self.gene_associations = {}

        # things to fill in optionally
        self.continuous_data_layers = {'layers': [], 'data': {}}

        # let's have this ready for convenience:
        self.genome_names = list(self.genome_descriptions.genomes.keys())

        if not skip_init:
            self.init()


    def init(self, populate_continuous_data=True):
        """Learn everything about genomes of interest.
        Calling this funciton will populate multiple critical dictionaries this class
        designed to give access to, including `self.genomes` and `self.gene_associations`.
        """

        # before going through all the gneomes, we will recover gene associations. the reason
        # we want to do it early on to make sure if there are incompatibilities between genome
        # names of interest and reosurces provided to learn gene-gene associations (such as the
        # pan database), they are discovered earlier than later:
        self.gene_associations = self.get_gene_associations()

      # now we will go through genomes of interest, and build our gigantic dictionary
        for genome_name in self.genome_names:
            self.genomes[genome_name] = {}

            # learn all about genes:
            self.genomes[genome_name]['genes'] = self.get_genes_dict(genome_name)

            # learn all about contigs:
            self.genomes[genome_name]['contigs'] = self.get_contigs_dict(genome_name)

        self.initialized = True

        if populate_continuous_data:
            self.populate_genome_continuous_data_layers()


    def get_gene_associations(self):
        """Recovers gene assoctiations through gene clusters found in a pan database.
        FIXME/TODO: also recover gene associations through user-provided input files.
        """

        d = {}

        # if we have a pan database, we will start with that, first, and learn gene clusters in it:
        if self.pan_db_path:
            pan_db = dbops.PanDatabase(self.pan_db_path)

            genome_names_missing_in_pan_db = [g for g in self.genome_names if g not in pan_db.genomes]
            if len(genome_names_missing_in_pan_db):
                raise ConfigError(f"You have provided a pan database to recover assocaitions between genes across "
                                  f"your genomes, but not all genome names in your list of genomes occur in this "
                                  f"pan database :/ Here is the list of genome names that you are missing: "
                                  f"{', '.join(genome_names_missing_in_pan_db)}")
            else:
                self.run.warning("Anvi'o found each of the genome name you are interested in the pan database you "
                                 "have provided. Which means the gene cluster information will be recovered for "
                                 "downstream magic.", header="PAN DATABASE LOOKS GOOD 🚀", lc="green")

            pan_db.disconnect()


            pan_db = dbops.PanSuperclass(self.args)
            pan_db.init_gene_clusters()

            d['anvio-pangenome'] = {}
            d['anvio-pangenome']['gene-cluster-name-to-genomes-and-genes'] = pan_db.gene_clusters
            d['anvio-pangenome']['genome-and-gene-names-to-gene-clusters'] = pan_db.gene_callers_id_to_gene_cluster

        return d


    def get_genes_dict(self, genome_name):
        """Learn everything about genes in a genome"""

        contigs_db_path = self.genome_descriptions.genomes[genome_name]['contigs_db_path']

        d = {}

        # learn gene calls, start-stop positions, and so on
        d['gene_calls'] = db.DB(contigs_db_path, None, ignore_version=True).smart_get(t.genes_in_contigs_table_name, 'gene_callers_id', self.genome_descriptions.genomes[genome_name]['gene_caller_ids'])

        # learn gene functions as well as gene amino acid and DNA sequences
        d['functions'], d['aa'], d['dna'] = self.genome_descriptions.get_functions_and_sequences_dicts_from_contigs_db(genome_name)

        return d


    def populate_continuous_GC_content_data(self):
        """Add sliding window GC-content change per contig"""

        if not self.initialized:
            raise ConfigError("You can't populate continuous data layers unless the class is properly initialized.")

        if 'GC_content' not in self.continuous_data_layers['layers']:
            self.continuous_data_layers['layers'].append('GC_content')

        self.progress.new('Populating continuous data', progress_total_items=len(self.genomes))
        for genome_name in self.genomes:
            self.progress.update(f"GC-content for {genome_name}", increment=True)

            if genome_name not in self.continuous_data_layers['data']:
                self.continuous_data_layers['data'][genome_name] = {}

            self.continuous_data_layers['data'][genome_name]['GC_content'] = {}

            for contig_name in self.genomes[genome_name]['contigs']['info']:
                contig_sequence = self.genomes[genome_name]['contigs']['dna'][contig_name]['sequence']
                self.continuous_data_layers['data'][genome_name]['GC_content'][contig_name] = utils.get_GC_content_for_sequence_as_an_array(contig_sequence)

        self.progress.end()


    def populate_genome_continuous_data_layers(self, skip_GC_content=False):
        """Function to populate continuous data layers.
        Calling this function will poupulate the `self.continuous_data_layers` data structure, which will
        look something like this:
            >>> {
            >>>   "layers": [
            >>>     "GC_content",
            >>>     "Another_layer"
            >>>   ],
            >>>   "data": {
            >>>     "GENOME_01": {
            >>>       "GC_content": {
            >>>         "CONTIG_A": [(...)],
            >>>         "CONTIG_B": [(...)],
            >>>          (...)
            >>>       },
            >>>       "Another_layer": {
            >>>         "CONTIG_A": [(...)],
            >>>         "CONTIG_B": [(...)],
            >>>          (...)
            >>>       }
            >>>     },
            >>>     "GENOME_02": {
            >>>       "GC_content": {
            >>>         "CONTIG_C": [(...)],
            >>>         "CONTIG_D": [(...)],
            >>>          (...)
            >>>       },
            >>>       "Another_layer": {
            >>>         "CONTIG_C": [(...)],
            >>>         "CONTIG_D": [(...)],
            >>>          (...)
            >>>       }
            >>>     },
            >>>     "GENOME_03": {
            >>>       "GC_content": {
            >>>         "CONTIG_E": [(...)],
            >>>         "CONTIG_F": [(...)],
            >>>          (...)
            >>>       },
            >>>       "Another_layer": {
            >>>         "CONTIG_E": [(...)],
            >>>         "CONTIG_F": [(...)],
            >>>          (...)
            >>>       }
            >>>     }
            >>>   }
            >>> }
        Please note: The number of data points in continuous data layers for a given contig will always
        be equal or less than the number of nucleotides in the same contig. When displaying htese data,
        one should always assume that the first data point matches to teh first nucleotide position, and
        cut the information in the layer short when necessary.
        """

        if not self.initialized:
            raise ConfigError("You can't populate continuous data layers unless the class is properly initialized.")

        # reset in case it was previously populated
        self.continuous_data_layers = {'layers': [], 'data': {}}

        for genome_name in self.genomes:
            self.continuous_data_layers['data'][genome_name] = {}

        # add GC-content
        if not skip_GC_content:
            self.populate_continuous_GC_content_data()


    def get_contigs_dict(self, genome_name):
        """Learn everything about contigs associated with a genome"""

        contigs_db_path = self.genome_descriptions.genomes[genome_name]['contigs_db_path']

        d = {}

        # get contig sequences
        if genome_name in self.genome_descriptions.internal_genome_names:
            raise ConfigError("This is not yet implemented :( Someone needs to find a rapid way to get contig "
                              "associated with a set of gene caller ids without having to create an instance "
                              "of ContigsSuperclass :/")
        else:
            d['info'] = db.DB(contigs_db_path, None, ignore_version=True).get_table_as_dict(t.contigs_info_table_name)
            d['dna'] = db.DB(contigs_db_path, None, ignore_version=True).get_table_as_dict(t.contig_sequences_table_name)

        return d


class AggregateFunctions:
    """Aggregate functions from anywhere.

    The purpose of this class is to collect functions from many distinct databases,
    including external genomes, internal genomes, and/or a genomes storage, and report
    a set of dictionaries that give access to the presence/absence and frequency
    of all functions annotated by a single source.

    One fancy function in AggregateFunctions is `report_functions_per_group_stats`. For
    instance, one could initiate the class in the following way to get a functions per
    group stats output file for functional enrichment analysis:

        >>> import argparse
        >>> import anvio.genomedescriptions as g
        >>> args = argparse.Namespace(external_genomes=external_genomes_path, annotation_source='KOfam')
        >>> groups = {'adoles': ['B_adolescentis', 'B_adolescentis_1_11', 'B_adolescentis_22L', 'B_adolescentis_6', 'B_adolescentis_ATCC_15703'],
                      'lactis': ['B_animalis', 'B_lactis_AD011', 'B_lactis_ATCC_27673', 'B_lactis_B420', 'B_lactis_BB_12', 'B_lactis_BF052'],
                      'longum': ['B_longum', 'B_longum_AH1206', 'B_longum_BBMN68', 'B_longum_BORI', 'B_longum_CCUG30698', 'B_longum_GT15']}
        >>> facc = g.AggregateFunctions(args, layer_groups=groups)
        >>> facc.report_functions_per_group_stats(output_file_path)

    For other uses of this class, see the `functional` mode in the interactive class
    which is accessed by anvi-display-functions, or anvi-script-gen-functions-per-group-stats-output
    that gives access to `report_functions_per_group_stats` function to generate the
    functions across groups output.

    Paremeters
    ==========
    args : argparse.Namespace object
        See the class header for options.
    layer_groups : dict
        When provided, the class will recognize that genomes belong to distinct groups
        and will prepare grouped frequency and presence-absence dicts, as well. This
        can ALTERNATIVELY be defined through a TAB-delimited input file passed through
        args.
    """

    def __init__(self, args, layer_groups=None, skip_sanity_check=False, skip_init=False, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.genomes_storage_path = A('genomes_storage')
        self.external_genomes_path = A('external_genomes')
        self.internal_genomes_path = A('internal_genomes')
        self.function_annotation_source = A('annotation_source')
        self.min_occurrence = A('min_occurrence') or 1
        self.aggregate_based_on_accession = A('aggregate_based_on_accession') or False
        self.aggregate_using_all_hits = A('aggregate_using_all_hits') or False
        self.layer_groups_input_file_path = A('groups_txt') or False
        self.print_genome_names_and_quit = A('print_genome_names_and_quit') or False
        self.functional_occurrence_table_output_path = A('functional_occurrence_table_output')
        self.functional_enrichment_output_path = A('output_file')

        # -----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----
        # these are some primary data structures this class reports

        # remember, 'key' here can be accession ids, or functio names
        # depending on `self.aggregate_based_on_accession`
        self.hash_to_key = {}

        # this variable makes sure even if functions are aggregated
        # using accession ids, there is a way to resolve function
        # names that correspond to each item. This is going to be a
        # life saver while trying to summarize things
        self.hash_to_function_dict = {}

        # distribution of 'keys' (i.e., accession ids or functions)
        # across genomes based on the frequency of observation or
        # presence absence. Having two disctionaries for this sounds
        # stupid at first (because it is), since the presence/absence
        # data can always be recovered from the frequency data. but
        # the momory fingerprint of these dicts are always going to be
        # nothing and will save additional steps in places where an
        # instance is used.
        self.functions_across_layers_frequency = {}
        self.functions_across_layers_presence_absence = {}

        # just like the previous two, but rather than genome names as
        # layers, these dicts will be tracking 'gruops' as defined by
        # the member variable, `self.layer_groups`. please note that
        # the presence-absence dictionary will not be composed of
        # binary variables, but will have the sum of presence absence
        # of a given function across all layers in a given group
        self.functions_across_groups_frequency = {}
        self.functions_across_groups_presence_absence = {}

        # two additional dicts to alwys be able to convert accession
        # ids to function names and vice versa. remember, an accession
        # id will resolve to a single function name (unless the provider
        # of function names did not screw things up), but a function
        # name can resolve to multiple accession ids, hence the latter
        # dict will contain will contain for each of its keys a list.
        self.accession_id_to_function_dict = {}
        self.function_to_accession_ids_dict = {}

        # to keep track of all layer names and where they are coming from.
        self.layer_names_considered = set({})
        self.layer_names_from_internal_genomes = []
        self.layer_names_from_external_genomes = []
        self.layer_names_from_genomes_storage = []

        # if there are 'groups' defined through the `layer_groups` variable
        # or through a `self.layer_groups_input_file_path`, this class will
        # automatically perform a functional enrichment analysis and will report
        # its output in the following dictionary.
        self.functional_enrichment_stats_dict = None
        # -----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----

        # this will summarize what happened in a text form.
        self.summary_markdown = None

        # -----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----
        # Here we will quickly deal with layer groups during the initialization of the class
        # this section of the init will establish a propeer `self.layer_groups` variable for
        # later use.
        self.layer_name_to_group_name = {}
        self.layer_groups_defined = False
        self.layer_groups = None

        if layer_groups or self.layer_groups_input_file_path:
            if layer_groups and not isinstance(layer_groups, dict):
                raise ConfigError("The variable `layer_groups` is supposed to be of type `dict`.")

            if self.layer_groups_input_file_path and layer_groups:
                raise ConfigError("You can either specify layer groups by passing a dictionary, or "
                                  "you can use the `layer_groups_input_file_path` argument, but not "
                                  "both :/")

            if self.layer_groups_input_file_path:
                self.layer_name_to_group_name, self.layer_groups = utils.get_groups_txt_file_as_dict(self.layer_groups_input_file_path)
            elif layer_groups:
                self.layer_groups = layer_groups

                # Finally, last AND least, a small helper dictionary we will use if there are groups defined
                # by the user:
                if self.layer_groups:
                    for group_name in self.layer_groups:
                        for layer_name in self.layer_groups[group_name]:
                            self.layer_name_to_group_name[layer_name] = group_name

            if self.layer_groups:
                self.layer_groups_defined = True

            group_names = sorted(list(self.layer_groups.keys()))
            self.run.info('Groups found and parsed', ', '.join(group_names))
        # -----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----

        self.key_hash_prefix = f"{'acc_' if self.aggregate_based_on_accession else 'func_'}"
        self.K = lambda: 'accession ID' if self.aggregate_based_on_accession else 'function'
        self.V = lambda: 'function' if self.aggregate_based_on_accession else 'accession ID'

        self.sanity_checked = False
        self.initialized = False

        if not skip_sanity_check:
            self.sanity_check()

        if not skip_init:
            self.init()


    def init(self):
        if self.initialized:
            raise ConfigError("Soemone already called the init function on this instance. You can't do it again. "
                              "Go get your own instance :(")

        # populate main dictionaries
        self._init_functions_from_int_ext_genomes()
        self._init_functions_from_genomes_storage()

        # show the user what genome names are being consdiered for this analysis
        if self.print_genome_names_and_quit:
            self.run.info(f"Genome names found (n={len(self.layer_names_considered)})", ' '.join(self.layer_names_considered))
            sys.exit()

        self._populate_group_dicts() # <-- this has to be called after all genomes are initialized

        if self.min_occurrence:
            num_occurrence_of_keys = [(c, sum(self.functions_across_layers_presence_absence[c].values())) for c in self.functions_across_layers_presence_absence]
            keys_to_remove = set([key for (key, frequency) in num_occurrence_of_keys if frequency < self.min_occurrence])

            # the following if block takes care of cleaning up both `self.accession_id_to_function_dict` and
            # `self.function_to_accession_ids_dict` dicts in a mindful fashion in addition to the cleanup of
            # all other dicts.
            if len(keys_to_remove):
                for key in keys_to_remove:
                    self.functions_across_layers_frequency.pop(key)
                    self.functions_across_layers_presence_absence.pop(key)

                    if self.layer_groups_defined:
                        self.functions_across_groups_frequency.pop(key)
                        self.functions_across_groups_presence_absence.pop(key)

                    self.hash_to_key.pop(key) if key in self.hash_to_key else None

                    # these are the trick ones since how this step should be handled will depend on
                    # what is key and what is value in the instance configuration:
                    if self.aggregate_based_on_accession:
                        if key in self.hash_to_key:
                            accession = self.hash_to_key[key]
                            function = self.accession_id_to_function_dict[accession]
                            self.accession_id_to_function_dict.pop(accession)
                            self.function_to_accession_ids_dict[function].pop(accession)
                            if not len(self.function_to_accession_ids_dict[function]):
                                self.function_to_accession_ids_dict.pop(function)
                    else:
                        if key in self.function_to_accession_ids_dict:
                            accessions = self.function_to_accession_ids_dict[key]
                            self.function_to_accession_ids_dict.pop(key)
                            for accession in accessions:
                                self.accession_id_to_function_dict.pop(accession)

                self.run.warning(f"As per your request, anvi'o removed {len(keys_to_remove)} {self.K()}s found in"
                                 f"{self.function_annotation_source} from downstream analyses since they occurred "
                                 f"in less than {self.min_occurrence} genomes.")

        if self.layer_groups:
            self.do_functional_enrichment_analysis()

        self.update_summary_markdown()

        self.initialized = True


    def sanity_check(self):
        if self.functional_occurrence_table_output_path:
            filesnpaths.is_output_file_writable(self.functional_occurrence_table_output_path)

        if self.functional_enrichment_output_path:
            filesnpaths.is_output_file_writable(self.functional_enrichment_output_path)

        if not self.function_annotation_source:
            raise ConfigError("Someone didn't specify any function annotation source and ended up in a bad "
                              "place. Aggregating functions require a source for functional annotations.")

        if not self.external_genomes_path and not self.internal_genomes_path and not self.genomes_storage_path:
            raise ConfigError("You must provide at least one source of genomes to this class :/")

        if self.min_occurrence and not isinstance(self.min_occurrence, int):
            raise ConfigError("Obviously, --min-occurrence must be an integer.")

        if self.min_occurrence < 1:
            raise ConfigError(f"What do you have in mind when you say I want my functions to occur in at least {self.min_occurrence} genomes?")

        if self.layer_groups_defined:
            groups_with_single_layers = set([])

            if not len(self.layer_groups) > 1:
                raise ConfigError("Layer groups must have two or more groups.")

            for layer_group in self.layer_groups:
                if not isinstance(self.layer_groups[layer_group], list):
                    raise ConfigError("Each layer group must be composed list of layer names :(")

                if not len(self.layer_groups[layer_group]) > 1:
                    groups_with_single_layers.add(layer_group)

                if len(set(self.layer_groups[layer_group])) != len(self.layer_groups[layer_group]):
                    raise ConfigError("Items in each layer group must be unique :/")

            if len(groups_with_single_layers):
                self.run.warning(f"In an ideal world, each group would describe at least two layer names. It is not "
                                 f"the case for {P('this group', len(groups_with_single_layers), alt='these groups')}: "
                                 f"{', '.join(groups_with_single_layers)}. That is OK and anvi'o will continue with this "
                                 f"analysis, but if something goes wrong with your stats or whatever, you will remember "
                                 f"this moment and go like, \"Hmm. That's why my adjusted q-values are like one point zero 🤔\".")

            # sanity check 3000 -- no joker shall pass:
            list_of_layer_names_lists = list(self.layer_groups.values())
            for i in range(0, len(list_of_layer_names_lists) - 1):
                for j in range(i + 1, len(list_of_layer_names_lists)):
                    co_occurring_names = set(list_of_layer_names_lists[i]).intersection(set(list_of_layer_names_lists[j]))
                    if len(co_occurring_names):
                        raise ConfigError(f"Layer names should occur in only one group, but AS YOU CAN GUESS BY NOW, that is not the case "
                                          f"with your groups :/ At the least, the layer name '{co_occurring_names.pop()}' occurs in more than "
                                          f"one group.")

        self.sanity_checked = True


    def update_combined_functions_dicts(self, genome_name, accession, function):
        """Modify accession dicts with new things.

        This function is necessary to avoid redundant code to handle function dicts of different kinds
        we get from int/external genomes and genome storage. a redesign of the genome storage will likely
        fix this problem in the future by unifying how functions are collected from different anvi'o databases.
        """

        if genome_name not in self.layer_names_considered:
            self.layer_names_considered.add(genome_name)

        key, value = (accession, function) if self.aggregate_based_on_accession else (function, accession)

        if not key or not len(key):
            raise ConfigError(f"Anvi'o is very sorry to tell you that the annotation source you have chosen "
                              f"here, '{self.function_annotation_source}', seem to include "
                              f"{self.V()}s with no {self.K()}s. Here is one example function with no "
                              f"{self.K()}: '{value}'. You will have to choose another annotation source :(")

        if not function:
            raise ConfigError(f"It saddens anvi'o to let you know that there are some {self.V}s in "
                              f"'{self.function_annotation_source}' that clearly are blank. Here is an "
                              f"example {self.K()} that has a blank {self.V()}: '{key}'. Either you "
                              f"need to choose another annotation source, or fix this problem with the "
                              f"existing annotations from {self.function_annotation_source} by using a "
                              f"combination of `anvi-export-functions` and `anvi-import-functions` (which "
                              f"is totally doable and you certainly can do it).")

        if self.aggregate_using_all_hits:
            pass
        else:
            key = key.split('!!!')[0]
            value = value.split('!!!')[0]

            # we wish to keep track of actual accessions and functions, too:
            accession, function = (key, value) if self.aggregate_based_on_accession else (value, key)

        # from now on we will only work with hashes of our keys, whether the keys here are function names or
        # accession ids as defined by self.aggregate_based_on_accession boolean. this is a necessary
        # complexity because function names are free text, can be very long, include weird characters. and
        # when we cluster data with those keys, some characters will be replaced with others or otherwise
        # they will break the newick file format and so on. using key hashes will make sure we don't get
        # screwed by that, and we will always use the lookup dict `self.hash_to_key` to find out what was
        # our key (and that's exactly what the next few lines of code do here):
        if key in self.hash_to_key:
            key_hash = self.hash_to_key[key]
        else:
            key_hash = self.key_hash_prefix + hashlib.sha224(key.encode('utf-8')).hexdigest()[0:12]
            self.hash_to_key[key] = key_hash
            self.hash_to_function_dict[key_hash] = {self.function_annotation_source: function}

        # --
        if key_hash not in self.functions_across_layers_frequency:
            self.functions_across_layers_frequency[key_hash] = Counter({})
            self.functions_across_layers_presence_absence[key_hash] = Counter({})

        self.functions_across_layers_frequency[key_hash][genome_name] += 1
        self.functions_across_layers_presence_absence[key_hash][genome_name] = 1

        if accession not in self.accession_id_to_function_dict:
            self.accession_id_to_function_dict[accession] = {self.function_annotation_source: function}

        if function not in self.function_to_accession_ids_dict:
            self.function_to_accession_ids_dict[function] = {self.function_annotation_source: set([accession])}
        else:
            self.function_to_accession_ids_dict[function][self.function_annotation_source].add(accession)

        return


    def check_layer_names(self, layer_names=[]):
        if not isinstance(layer_names, list):
            raise ConfigError("`layer_names` must be of type list :/")

        already_in_the_dict = [g for g in layer_names if g in self.layer_names_considered]
        if len(already_in_the_dict):
            raise ConfigError(f"Anvi'o is not happy because there are some genome or metagenome names that are not unique "
                              f"across all input databases :/ Here is an example: {already_in_the_dict[0]}.")
        else:
            # you good fam
            pass


    def _init_functions_from_int_ext_genomes(self):
        if not self.external_genomes_path and not self.internal_genomes_path:
            return

        g = GenomeDescriptions(self.args, run=terminal.Run(verbose=False))
        g.load_genomes_descriptions(skip_sanity_check=True)
        g.init_functions()

        self.layer_names_from_internal_genomes = copy.deepcopy(g.internal_genome_names)
        self.layer_names_from_external_genomes = copy.deepcopy(g.external_genome_names)

        for genome_name in g.genomes:
            self.check_layer_names([genome_name])

            gene_functions_in_genome_dict, _, _= g.get_functions_and_sequences_dicts_from_contigs_db(genome_name, requested_source_list=[self.function_annotation_source], return_only_functions=True)
            # reminder, an entry in gene_functions_in_genome_dict looks like this:
            # 2985: {'COG20_PATHWAY': ('COG0073!!!COG0143', 'Aminoacyl-tRNA synthetases', 0)}
            for entry in gene_functions_in_genome_dict.values():
                accession, function, e_value = entry[self.function_annotation_source]

                self.update_combined_functions_dicts(genome_name, accession, function)


    def _init_functions_from_genomes_storage(self):
        if not self.genomes_storage_path:
            return

        from anvio.genomestorage import GenomeStorage
        g = GenomeStorage(storage_path=self.genomes_storage_path, function_annotation_sources=[self.function_annotation_source], run=terminal.Run(verbose=False), progress=self.progress, skip_init=True)

        # make sure we are not overwriting existing genome names in int or ext genomes:
        genome_names_in_storage_db = list(g.db.get_single_column_from_table(t.genome_info_table_name, 'genome_name', unique=True))

        self.layer_names_from_genomes_storage = copy.deepcopy(genome_names_in_storage_db)

        self.check_layer_names(genome_names_in_storage_db)

        gene_functions_in_genomes_dict, _ = g.get_gene_functions_in_genomes_dict()
        for entry in gene_functions_in_genomes_dict.values():
            # an entry in gene_functions_in_genomes_dict looks lke this:
            # 72645: {'genome_name': 'B_lactis_BF052', 'gene_callers_id': 443, 'source': 'COG20_PATHWAY', 'accession': 'COG0207', 'function': 'Thymidylate biosynthesis', 'e_value': 4.2e-149}
            genome_name, accession, function = entry['genome_name'], entry['accession'], entry['function']

            self.update_combined_functions_dicts(genome_name, accession, function)


    def _populate_group_dicts(self):
        if not self.layer_groups_defined:
            return

        # let's check if layer names from the groups file has anything to do with
        # the layer names that are considered at this point.
        layer_names_with_group_association = [l for l in self.layer_names_considered if l in self.layer_name_to_group_name]
        if not len(layer_names_with_group_association):
            raise ConfigError(f"Something is wrong here :/ You provide some group associations for your genomes, however, "
                              f"the genome names anvi'o gathered from your databases does not seem to have anything to do "
                              f"with those names. You probably need to fix this. Here is an example genome name from your "
                              f"group associations: '{list(self.layer_name_to_group_name.keys())[0]}'. In comparison, here "
                              f"Here is an example genome name from your anvi'o databases: '{self.layer_names_considered.pop()}'.")

        for key_hash in self.functions_across_layers_frequency:
            if key_hash not in self.functions_across_groups_frequency:
                self.functions_across_groups_frequency[key_hash] = Counter({})
                self.functions_across_groups_presence_absence[key_hash] = Counter({})

            for layer_name in self.layer_names_considered:
                if layer_name in self.layer_name_to_group_name:
                    if layer_name in self.functions_across_layers_frequency[key_hash]:
                        group_name = self.layer_name_to_group_name[layer_name]
                        self.functions_across_groups_frequency[key_hash][group_name] += self.functions_across_layers_frequency[key_hash][layer_name]
                        self.functions_across_groups_presence_absence[key_hash][group_name] += 1


    def update_summary_markdown(self):
        G = lambda x: '\n'.join(['* %s' % l for l in x]) if len(x) else "*None :/*"

        self.summary_markdown = (f"### Quick overview\nUsing the function annotation source **'{self.function_annotation_source}'**, anvi'o "
                                 f"aggregated **{len(self.hash_to_key)} unique {self.K()}s** that occurred in **at least {self.min_occurrence}** "
                                 f"of the total {len(self.layer_names_considered)} *layers*. Here we use the term 'layer' instead of 'genomes', "
                                 f"since what anvi'o assumes to be a genome might be a metagenome depending on the contigs database you have provided. "
                                 f"If you know what you have, feel free to replace the term 'layer' with 'genome' in your mind.")

        self.summary_markdown += (f"\n\n**Internal genomes** ({P('layer', len(self.layer_names_from_internal_genomes))}):\n\n{G(self.layer_names_from_internal_genomes)}"
                                  f"\n\n**External genomes** ({P('layer', len(self.layer_names_from_external_genomes))}):\n\n{G(self.layer_names_from_external_genomes)}"
                                  f"\n\n**Genomes storage** ({P('layer', len(self.layer_names_from_genomes_storage))}):\n\n{G(self.layer_names_from_genomes_storage)}")

        if self.layer_groups:
            self.summary_markdown += (f"\n\n### Functional enrichment analysis\nAnvi'o also performed a functional enrichment analysis based on {len(self.layer_groups)} "
                                      f"groups defined by the user. The results of this analysis is shown in your view with the these additional layers: "
                                      f"'enrichment_score', 'unadjusted_p_value', 'adjusted_q_value', 'associated_groups'. You can learn more about the "
                                      f"details of the meaning of these columns [here](https://merenlab.org/software/anvio/help/main/artifacts/functional-enrichment-txt/). "
                                      f"Here is how those groups were defined:")
            for group_name in self.layer_groups:
                self.summary_markdown += (f"\n\n**Group '{group_name}'** ({P('layer', len(self.layer_groups[group_name]))}):\n\n{G(self.layer_groups[group_name])}")


    def do_functional_enrichment_analysis(self):
        """Performs functional enrichment analysis if user defined layer groups.

        This function fills in the the variable `self.functional_enrichment_stats_dict` so
        the downstream analyses can use it to do fancy things.
        """

        if not self.layer_groups:
            raise ConfigError("Functional enrichment analysis requires the `self.layer_groups` to be "
                              "initialized. But someone called this function without first intializing "
                              "groups :/ ")

        # FIXME: this is kind of a stupid design. we create this directory even if the user has declared output file names
        #        for both teh functional occurrence table output and functional enrichment output.
        output_directory = filesnpaths.get_temp_directory_path()

        if not self.functional_occurrence_table_output_path:
            self.functional_occurrence_table_output_path = os.path.join(output_directory, 'FUNC_OCCURENCE_STATS.txt')

        if not self.functional_enrichment_output_path:
            self.functional_enrichment_output_path = os.path.join(output_directory, 'FUNC_ENRICHMENT_OUTPUT.txt')

        # get functions per group stats file
        self.report_functions_per_group_stats(self.functional_occurrence_table_output_path, quiet=True)

        # run the enrichment analysis
        self.functional_enrichment_stats_dict = utils.run_functional_enrichment_stats(functional_occurrence_stats_input_file_path=self.functional_occurrence_table_output_path,
                                                                                      enrichment_output_file_path=self.functional_enrichment_output_path,
                                                                                      run=self.run,
                                                                                      progress=self.progress)


    def report_functions_across_genomes(self, output_file_prefix, quiet=False):
        """Reports text files for functions across genomes data"""

        output_file_path_for_frequency_view = f"{os.path.abspath(output_file_prefix)}-FREQUENCY.txt"
        output_file_path_for_presence_absence_view = f"{os.path.abspath(output_file_prefix)}-PRESENCE-ABSENCE.txt"

        filesnpaths.is_output_file_writable(output_file_path_for_frequency_view)
        filesnpaths.is_output_file_writable(output_file_path_for_presence_absence_view)

        with open(output_file_path_for_frequency_view, 'w') as frequency_output, open(output_file_path_for_presence_absence_view, 'w') as presence_absence_output:
            layer_names = sorted(list(self.layer_names_considered))
            frequency_output.write('\t'.join(['key'] + layer_names + [self.function_annotation_source]) + '\n')
            presence_absence_output.write('\t'.join(['key'] + layer_names + [self.function_annotation_source]) + '\n')

            for key in self.functions_across_layers_frequency:
                function = self.hash_to_function_dict[key][self.function_annotation_source]

                frequency_data = [f"{self.functions_across_layers_frequency[key][l] if l in self.functions_across_layers_frequency[key] else 0}" for l in layer_names]
                frequency_output.write('\t'.join([key] + frequency_data + [function]) + '\n')

                presence_absence_data = [f"{self.functions_across_layers_presence_absence[key][l] if l in self.functions_across_layers_presence_absence[key] else 0}" for l in layer_names]
                presence_absence_output.write('\t'.join([key] + presence_absence_data + [function]) + '\n')

        if not quiet:
            self.run.info('Functions across genomes (frequency)', output_file_path_for_frequency_view)
            self.run.info('Functions across genomes (presence/absence)', output_file_path_for_presence_absence_view)


    def report_functions_per_group_stats(self, output_file_path, quiet=False):
        """A function to summarize functional occurrence for groups of genomes.

        Please note that this function will not report functions that are associated
        with ALL groups.
        """

        filesnpaths.is_output_file_writable(output_file_path)

        if not self.layer_groups_defined:
            raise ConfigError("No groups seem to have been defined. This function is useless without :/")

        group_names = sorted(list(self.layer_groups.keys()))

        num_groups = len(group_names)

        group_counts = dict([(g, len(self.layer_groups[g])) for g in group_names])

        d = {}
        num_skipped = 0 # keep track of how many functions we skip reporting on

        key_hash_represents = "accession" if self.aggregate_based_on_accession else "function"
        self.run.info(f"Number of {self.function_annotation_source} {key_hash_represents}s found across {len(group_names)} groups", len(self.functions_across_groups_presence_absence.keys()))

        for key_hash in self.functions_across_groups_presence_absence:
            # learn which groups are associated with this function
            associated_groups = [g for g in group_names if self.functions_across_groups_presence_absence[key_hash][g]]

            # if the function is associated with all groups, simply skip that entry
            if len(associated_groups) == num_groups:
                num_skipped += 1
                continue

            function = self.hash_to_function_dict[key_hash][self.function_annotation_source]

            d[key_hash] = {}
            d[key_hash]['function'] = function
            d[key_hash]['accession'] = ','.join(self.function_to_accession_ids_dict[function][self.function_annotation_source])
            d[key_hash]['associated_groups'] = ','.join(associated_groups)

            for group_name in group_names:
                d[key_hash][f"N_{group_name}"] = group_counts[group_name]
                if group_name in self.functions_across_groups_presence_absence[key_hash]:
                    d[key_hash][f"p_{group_name}"] = self.functions_across_groups_presence_absence[key_hash][group_name] / group_counts[group_name]
                else:
                    d[key_hash][f"p_{group_name}"] = 0

        self.run.info(f"Number of {self.function_annotation_source} {key_hash_represents}s associated with all groups and SKIPPED", num_skipped)
        self.run.info(f"Number of {self.function_annotation_source} {key_hash_represents}s in final occurrence table", len(d))

        if not len(d):
            raise ConfigError("Something weird is happening here :( It seems every single function across your genomes "
                              "is associated with all groups you have defined. There is nothing much anvi'o can work with "
                              "here. If you think this is a mistake, please let us know.")

        if len(d) < 2:
            raise ConfigError("Oh, dear. It seems only one function is differentially present across the genome "
                              "groups you have defined. There is nothing much anvi'o can work with "
                              "here. If you think this is a mistake, please let us know.")

        static_column_names = ['key', 'function', 'accession', 'associated_groups']
        dynamic_column_names = []
        [dynamic_column_names.extend([f'p_{g}', f'N_{g}']) for g in group_names]

        utils.store_dict_as_TAB_delimited_file(d, output_file_path, headers=static_column_names+dynamic_column_names)

        if not quiet:
            self.run.info('Functions per group stats file', output_file_path)
