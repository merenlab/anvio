# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module to dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sene of internal or external genome descriptions is also welcome.
"""

import os
import sys
import copy
import hashlib
import argparse

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.ccollections as ccollections

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
        self.functions_are_available = False
        self.function_annotation_sources = set([])

        self.input_file_for_internal_genomes = A('internal_genomes')
        self.input_file_for_external_genomes = A('external_genomes')
        self.list_hmm_sources = A('list_hmm_sources')          # <<< these two are look out of place, but if the args requests
        self.list_available_gene_names = A('list_available_gene_names') #     for information about HMMs, this is the bets place to set them
        self.gene_caller = A('gene_caller')

        if self.input_file_for_internal_genomes or self.input_file_for_external_genomes:
            self.read_genome_paths_from_input_files()


    def names_check(self):
        i, n = list(self.internal_genomes_dict.keys()), list(self.external_genomes_dict.keys())

        if not i and not n:
            raise ConfigError("You in fact tried to create a genomes storage file without providing any internal or external genome\
                                descriptions! You got 5 anvi'o points for being awesome, but this is not gonna work since you really\
                                need to provide at least one of those descriptions :/")

        if len(i) + len(n) != len(set(i + n)):
            raise ConfigError("Each entry both in internal and external genome descriptions should have a unique 'name'. This does not\
                                seem to be the case with your input :/")


    def read_genome_paths_from_input_files(self):
        """Reads internal and external genome files, populates self.genomes"""

        fields_for_internal_genomes_input = ['name', 'bin_id', 'collection_id', 'profile_db_path', 'contigs_db_path']
        fields_for_external_genomes_input = ['name', 'contigs_db_path']

        self.internal_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(self.input_file_for_internal_genomes, expected_fields=fields_for_internal_genomes_input) if self.input_file_for_internal_genomes else {}
        self.external_genomes_dict = utils.get_TAB_delimited_file_as_dictionary(self.input_file_for_external_genomes, expected_fields=fields_for_external_genomes_input) if self.input_file_for_external_genomes else {}


    def list_HMM_info_and_quit(self):
        hmm_sources_in_all_genomes = self.get_HMM_sources_common_to_all_genomes(dont_raise=True)

        # since we know hmm sources in `hmm_sources_in_all_genomes` are common to all genomes,
        # we could use any of those genomes to learn about the specifics of them. here we take
        # the first one from `self.genomes`
        hmm_sources_info = dbops.ContigsDatabase(list(self.genomes.values())[0]['contigs_db_path']).db.get_table_as_dict(t.hmm_hits_info_table_name)

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


    def get_HMM_sources_common_to_all_genomes(self, dont_raise=False):
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

        if not len(hmm_sources_in_all_genomes):
            if dont_raise:
                return None

            raise ConfigError("There are no HMM sources among your external genomes that occur in every genome :/")

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
        for source, input_file in [(self.external_genomes_dict, self.input_file_for_external_genomes), (self.internal_genomes_dict, self.input_file_for_internal_genomes)]:
            for genome_name in source:
                self.genomes[genome_name] = source[genome_name]
                for db_path_var in ['contigs_db_path', 'profile_db_path']:
                    if db_path_var not in self.genomes[genome_name]:
                        continue
                    path = self.genomes[genome_name][db_path_var]

                    if not path:
                        raise ConfigError("Bad news: anvi'o was loading genome desriptions, and it run into an empty path for\
                                           the genome %s. How did this happen? HOW? :(" % genome_name)

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


    def get_functions_and_sequences_dicts_from_contigs_db(self, genome_name):
        g = self.genomes[genome_name]

        args = argparse.Namespace(contigs_db=g['contigs_db_path'])
        contigs_super = dbops.ContigsSuperclass(args, r=anvio.terminal.Run(verbose=False))

        if self.functions_are_available:
            contigs_super.init_functions(requested_sources=list(self.function_annotation_sources))
            function_calls_dict = contigs_super.gene_function_calls_dict
        else:
            function_calls_dict = {}

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
                self.run.warning("None of your genomes seem to have any functional annotation. No biggie. Things will continue to work. But\
                                  then your genomes have no functional annotation. SAD.")
            else:
                self.run.warning("Some of your genomes (%d of the %d, to be precise) seem to have no functional annotation. Since this workflow\
                                  can only use matching functional annotations across all genomes involved, having even one genome without\
                                  any functions means that there will be no matching function across all. Things will continue to work, but\
                                  you will have no functions at the end for your gene clusters." % \
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


    def get_genome_hash_for_external_genome(self, entry):
        utils.is_contigs_db(entry['contigs_db_path'])
        contigs_db = dbops.ContigsDatabase(entry['contigs_db_path'])
        genome_hash = contigs_db.meta['contigs_db_hash']
        contigs_db.disconnect()

        return genome_hash


    def get_genome_hash_for_internal_genome(self, entry):
        utils.is_contigs_db(entry['contigs_db_path'])
        split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(entry)
        contigs_db = dbops.ContigsDatabase(entry['contigs_db_path'])
        genome_hash = hashlib.sha224('_'.join([''.join(split_names_of_interest), contigs_db.meta['contigs_db_hash']]).encode('utf-8')).hexdigest()[0:12]
        contigs_db.disconnect()

        return genome_hash


    def init_external_genomes(self):
        self.progress.new('Initializing external genomes')
        for genome_name in self.external_genome_names:
            c = self.genomes[genome_name]
            c['external_genome'] = True

            self.progress.update('working on %s' % (genome_name))

            contigs_db_summary = summarizer.ContigSummarizer(c['contigs_db_path']).get_contigs_db_info_dict(gene_caller_to_use=self.gene_caller)

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

                utils.is_profile_db_and_contigs_db_compatible(c['profile_db_path'], c['contigs_db_path'])

                split_names_of_interest = self.get_split_names_of_interest_for_internal_genome(c)

                # here we are using the get_contigs_db_info_dict function WITH split names we found in the collection
                # which returns a partial summary from the contigs database focusing only those splits. a small workaround
                # to be able to use the same funciton for bins in collections:
                contigs_summary = summarizer.ContigSummarizer(c['contigs_db_path'])
                summary_from_contigs_db_summary = contigs_summary.get_contigs_db_info_dict(split_names=split_names_of_interest,
                                                                                           gene_caller_to_use=self.gene_caller)

                for key in summary_from_contigs_db_summary:
                    c[key] = summary_from_contigs_db_summary[key]

        self.progress.end()

        self.run.info('Internal genomes', '%d have been initialized.' % len(self.internal_genome_names))


    def get_split_names_of_interest_for_internal_genome(self, entry):
        utils.is_profile_db(entry['profile_db_path'])
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
        self.list_HMM_info_and_quit()

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

        if not self.full_init:
            # if this is not full init, stop the sanity check here.
            self.run.warning("You (or the programmer) requested genome descriptions for your internal and/or external\
                              genomes to be loaded without a 'full init'. There is nothing for you to be concerned.\
                              This is just a friendly reminder to make sure if something goes terribly wrong (like your\
                              computer sets itself on fire), this may be the reason.")
            return

        # make sure HMMs for SCGs were run for every contigs db:
        genomes_missing_hmms_for_scgs =  [g for g in self.genomes if not self.genomes[g]['hmms_for_scgs_were_run']]
        if len(genomes_missing_hmms_for_scgs):
            if len(genomes_missing_hmms_for_scgs) == len(self.genomes):
                self.run.warning("The contigs databases you are using for this analysis are missing HMMs for single-copy core genes.\
                                  Maybe you haven't run `anvi-run-hmms` on your contigs database, or they didn't contain any hits.\
                                  It is perfectly legal to have anvi'o contigs databases without HMMs or SCGs for things to work,\
                                  but we wanted to give you heads up so you can have your 'aha' moment if you see funny things in\
                                  the interface.")
            else:
                raise ConfigError("Some of the genomes you have for this analysis are missing HMM hits for SCGs (%d of %d of them, to be precise). You\
                                    can run `anvi-run-hmms` on them to recover from this. Here is the list: %s" % \
                                                    (len(genomes_missing_hmms_for_scgs), len(self.genomes), ','.join(genomes_missing_hmms_for_scgs)))

        # make sure genome names are not funny (since they are going to end up being db variables soon)
        [utils.is_this_name_OK_for_database('genome name "%s"' % genome_name, genome_name) for genome_name in self.genomes]

        # figure out whether there are genomes with gene calls that are NOT processed
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
            self.run.warning("Some of your genomes had gene calls identified by gene callers other than\
                              the gene caller anvi'o used, which was set to '%s' either by default, or because you asked for it.\
                              The following genomes contained genes that were not processed (this may be exactly what you expect\
                              to happen, but if was not, you may need to use the `--gene-caller` flag to make sure anvi'o is using\
                              the gene caller it should be using): %s." % \
                                            (gene_caller, ', '.join(info)), header="PLEASE READ CAREFULLY", lc='green')

        # check whether every genome has at least one gene call.
        genomes_with_no_gene_calls = [g for g in self.genomes if not self.genomes[g]['num_genes']]
        if len(genomes_with_no_gene_calls):
            raise ConfigError("Well, %d of your %d genomes had 0 gene calls. We can't think of any reason to include genomes that\
                               contain no gene calls into a genomes, hence, we are going to stop here and ask you to remove these\
                               genomes from your analysis first: %s. If you think this is a dumb thing to do, and they should be\
                               in the genomes storage for reasons you know and we don't, please get in touch with us, and we will\
                               be happy to reconsider." % (len(genomes_with_no_gene_calls), len(self.genomes), ', '.join(genomes_with_no_gene_calls)))
