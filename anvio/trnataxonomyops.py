#!/usr/bin/env python3
# -*- coding: utf-8
"""
Classes to use local databases for tRNA taxonomy to affiliate tRNA seqeunces in anvi'o
databases with taxon names.
"""

import os
import sys
import copy
import shutil
import hashlib
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sps

from collections import OrderedDict, Counter

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.sharedtaxonomyops as sharedtaxonomy
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.genomedescriptions import MetagenomeDescriptions
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.xxxtaxonomy import PopulateContigsDatabaseWithTaxonomy
from anvio.dbops import ContigsSuperclass, ContigsDatabase, ProfileSuperclass, ProfileDatabase


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print

HASH = lambda d: str(hashlib.sha224(''.join([str(d[level]) for level in constants.levels_of_taxonomy]).encode('utf-8')).hexdigest()[0:8])

class TRNATaxonomyContext(object):
    """The purpose of this base class is ot define file paths and constants for trna taxonomy ops."""

    def __init__(self, trna_taxonomy_data_dir=None, scgs_taxonomy_remote_database_url=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        # hard-coded GTDB variables. poor design, but I don't think we are going do need an
        # alternative to GTDB.
        self.target_database = "GTDB"
        self.target_database_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"

        # some variables from anvi'o constants
        self.hmm_source_for_trna_genes = constants.default_hmm_source_for_trna_genes
        self.default_trna_taxonomy_data_dir = constants.default_trna_taxonomy_data_dir
        self.default_anticodons_for_taxonomy = constants.default_anticodons_for_taxonomy
        self.levels_of_taxonomy = constants.levels_of_taxonomy

        # these are all the user accessible paths. defaults will serve well for all applications,
        # but these can be used for debugging.
        self.tRNA_taxonomy_data_dir = (os.path.abspath(trna_taxonomy_data_dir) if trna_taxonomy_data_dir else None) or (os.path.join(self.default_trna_taxonomy_data_dir, self.target_database))
        self.accession_to_taxonomy_file_path = os.path.join(self.tRNA_taxonomy_data_dir, 'ACCESSION_TO_TAXONOMY.txt.gz')
        self.database_version_file_path = os.path.join(self.tRNA_taxonomy_data_dir, 'VERSION')
        self.search_databases_dir_path = os.path.join(self.tRNA_taxonomy_data_dir, 'ANTICODON_SEARCH_DATABASES')
        self.target_database_URL = scgs_taxonomy_remote_database_url or self.target_database_URL

        # some dictionaries for convenience. we set them up here, but the proper place to sanity check
        # them may be somewhere else. for instance, when this class is inheritded by SetupLocalTRNATaxonomyData
        # the paths will not point to an actual file, but when it is inherited by PopulateContigsDatabaseWithTRNATaxonomy,
        # they better point to actual files.
        self.anticodons = dict([(anticodon, {'db': os.path.join(self.search_databases_dir_path, anticodon)}) for anticodon in self.default_anticodons_for_taxonomy])

        # set version for ctx, so we know what version of the databases are on disk
        if os.path.exists(self.database_version_file_path):
            self.trna_taxonomy_database_version = open(self.database_version_file_path).readline().strip()
        else:
            self.trna_taxonomy_database_version = None

        self.accession_to_taxonomy_dict = xxxtaxonomy.get_accession_to_taxonomy_dict(self.accession_to_taxonomy_file_path, self.levels_of_taxonomy, self.progress, self.run)

# here we create an instance for the module. the idea is to overwrite it if
# it is necessary to overwrite some of the defaults
ctx = TRNATaxonomyContext()


class SanityCheck(object):
    def __init__(self):
        if self.skip_sanity_check:
            self.run.warning("We are skipping all sanity checks :( Dangerous stuff is happening.")
        else:
            self.sanity_check()


    def sanity_check(self):
        if sorted(constants.default_anticodons_for_taxonomy) != sorted(self.ctx.default_anticodons_for_taxonomy):
            raise ConfigError("Oh no. The anticodons designated to be used for all tRNA taxonomy tasks in the constants.py "
                              "are not the same names described in locally known HMMs to remote FASTA files "
                              "conversion table definedd in SetupLocalTRNATaxonomyData module. If this makes zero "
                              "sense to you please ask a developer.")

        if not self.ctx.tRNA_taxonomy_data_dir:
            raise ConfigError("`SetupLocalTRNATaxonomyData` class is upset because it was inherited without "
                              "a directory for tRNA taxonomy data to be stored :( This variable can't be None.")

        if self.user_taxonomic_level and self.user_taxonomic_level not in constants.levels_of_taxonomy:
            raise ConfigError("The taxonomic level %s is not a level anvi'o knows about. Here is the list of "
                              "taxonomic levels anvi'o recognizes: %s" % (', '.join(constants.levels_of_taxonomy)))

        # sanity checks specific to classes start below
        if self.__class__.__name__ in ['SetupLocalTRNATaxonomyData']:
            pass

        if self.__class__.__name__ in ['SetupLocalTRNATaxonomyData', 'PopulateContigsDatabaseWithTRNATaxonomy']:
            if self.user_taxonomic_level:
                raise ConfigError("There is no need to set a taxonomic level while working with the class SetupLocalTRNATaxonomyData "
                                  "or PopulateContigsDatabaseWithTRNATaxonomy. Something fishy is going on :/")

        if self.__class__.__name__ in ['PopulateContigsDatabaseWithTRNATaxonomy', 'TRNATaxonomyEstimatorSingle', 'TRNATaxonomyEstimatorMulti']:
            if not os.path.exists(self.ctx.tRNA_taxonomy_data_dir):
                raise ConfigError("Anvi'o could not find the data directory for the tRNA taxonomy setup. If you have "
                                  "a non-default location for your tRNA taxonomy databases, please use the parameter "
                                  "`--trna-taxonomy-data-dir` parameter). Anvi'o tried to find your files here: '%s'" % (self.ctx.tRNA_taxonomy_data_dir))

            if not os.path.exists(self.ctx.accession_to_taxonomy_file_path) or not os.path.exists(self.ctx.database_version_file_path):
                raise ConfigError("While your tRNA taxonomy data dir seems to be in place, it is missing at least one critical "
                                  "file. This is not someting you can add or remove as this file is distributed with anvi'o "
                                  "releases :( Please get in touch with a developer, or fix it if you are one.")

            if not os.path.exists(self.ctx.accession_to_taxonomy_file_path):
                raise ConfigError("While your tRNA taxonomy data dir seems to be in place, it is missing at least one critical "
                                  "file (in this case, the file to resolve accession IDs to taxon names). This is not someting "
                                  "you can add or remove as this file is distributed with anvi'o releases :( Please get in touch "
                                  "with a developer, or fix it if you are one.")

            ###########################################################
            # PopulateContigsDatabaseWithTRNATaxonomy
            ###########################################################
            if self.__class__.__name__ in ['PopulateContigsDatabaseWithTRNATaxonomy']:
                for prefix in ['.nhr', '.nin', '.nsq']:
                    missing_anticodon_databases = [anticodon for anticodon in self.ctx.anticodons if not os.path.exists(self.ctx.anticodons[anticodon]['db'] + '.nhr')]
                    if len(missing_anticodon_databases):
                        raise ConfigError("OK. It is very likley that if you run `anvi-setup-scg-taxonomy` first you will be golden. "
                                          "Because even though anvi'o found the directory for taxonomy headquarters, "
                                          "your setup seems to be missing %d of %d databases required for everything to work "
                                          "with the current genes configuration of this class (sources say this is a record, FYI)." % \
                                                    (len(missing_anticodon_databases), len(self.ctx.anticodons)))

            ###########################################################
            # TRNATaxonomyEstimatorSingle
            #
            # Note: if something down below complains about a paramter
            #       because that actually belongs to the multi estimator
            #       class, you may need to set it to null in the class
            #       TRNATaxonomyArgs for single estimator
            #       initiation if clause
            ###########################################################
            if self.__class__.__name__ in ['TRNATaxonomyEstimatorSingle']:
                if self.metagenomes:
                    raise ConfigError("Taxonomy estimation classes have been initiated with a single contigs database, but your "
                            "arguments also include input for metagenomes. It is a no no. Please choose either. ")

                if self.output_file_prefix:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file path "
                                      "than an output file prefix.")

                if self.output_file_path:
                    filesnpaths.is_output_file_writable(self.output_file_path)

                if self.raw_output or self.matrix_format:
                    raise ConfigError("Haha in this mode you can't ask for the raw output or matrix format .. yet (we know that "
                                      "the parameter space of this program is like a mine field and we are very upset about it "
                                      "as well).")

                if not self.contigs_db_path:
                    raise ConfigError("For these things to work, you need to provide a contigs database for the anvi'o SCG "
                                      "taxonomy workflow :(")

                utils.is_contigs_db(self.contigs_db_path)

                trna_taxonomy_was_run = ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet).meta['trna_taxonomy_was_run']
                trna_taxonomy_database_version = ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet).meta['trna_taxonomy_database_version']
                if not trna_taxonomy_was_run:
                    raise ConfigError("It seems the SCG taxonomy tables were not populated in this contigs database :/ Luckily it "
                                      "is easy to fix that. Please see the program `anvi-run-trna-taxonomy`.")

                if trna_taxonomy_database_version != self.ctx.trna_taxonomy_database_version:
                    self.progress.reset()
                    self.run.warning("The SCG taxonomy database on your computer has a different version (%s) than the SCG taxonomy information "
                                     "stored in your contigs database (%s). This is not a problem and things will most likely continue to work "
                                     "fine, but we wanted to let you know. You can get rid of this warning by re-running `anvi-run-trna-taxonomy` "
                                     "on your database." % (self.ctx.trna_taxonomy_database_version, trna_taxonomy_database_version))

                if self.profile_db_path:
                    utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

                if self.collection_name and not self.profile_db_path:
                    raise ConfigError("If you are asking anvi'o to estimate taxonomy using a collection, you must also provide "
                                      "a profile database to this program.")

                if self.metagenome_mode and self.collection_name:
                    raise ConfigError("You can't ask anvi'o to treat your contigs database as a metagenome and also give it a "
                                      "collection.")

                if self.scg_name_for_metagenome_mode and not self.metagenome_mode:
                    raise ConfigError("If you are not running in `--metagenome-mode`, there is no use to define a SCG name for "
                                      "this mode :/")

                if self.scg_name_for_metagenome_mode and self.scg_name_for_metagenome_mode not in self.ctx.anticodons:
                    raise ConfigError("We understand that you wish to work with '%s' to study the taxonomic make up of your contigs "
                                      "database in metagenome mode. But then this gene is not one of those anvi'o recognizes as "
                                      "suitable SCGs to do that. Here is a list for you to choose from: '%s'." \
                                                            % (self.scg_name_for_metagenome_mode, ', '.join(self.ctx.anticodons.keys())))

                if self.compute_scg_coverages and not self.profile_db_path:
                    raise ConfigError("The flag `--compute-scg-coverages` is only good if there is a non-blank profile database around "
                                      "from which anvi'o can learn coverage statistics of genes across one or more samples :/")

                if self.profile_db_path and self.metagenome_mode and not self.compute_scg_coverages:
                    raise ConfigError("You have a profile database and you have asked anvi'o to estimate taxonomy in metagenome mode, "
                                      "but you are not asking anvi'o to compute SCG coverages which doesn't make much sense :/ Removing "
                                      "the profile database from this command or addint the flag `--compute-scg-coverages` would have "
                                      "made much more sense.")

                if self.profile_db_path and not self.metagenome_mode and not self.collection_name:
                    raise ConfigError("You have a profile database, and you are not in metagenome mode. In this case anvi'o will try to "
                                      "estimate coverages of SCGs in bins after estimating their taxonomy, but for that, you need to "
                                      "also provide a collection name. You can see what collections are available in your profile database "
                                      "you can use the program `anvi-show-collections-and-bins`, and then use the parameter "
                                      "`--collection-name` to tell anvi'o which one to use.")

                if self.update_profile_db_with_taxonomy:
                    if not self.metagenome_mode:
                        raise ConfigError("Updating the profile database with taxonomy layer data is only possible in metagenome "
                                          "mode :/ And not only that, you should also instruct anvi'o to compute single-copy core "
                                          "gene coverages.")

                    if not self.compute_scg_coverages:
                        raise ConfigError("You wish to update the profile database with taxonomy, but this will not work if anvi'o "
                                          "is computing coverages values of SCGs across samples (pro tip: you can ask anvi'o to do "
                                          "it by adding the flag `--compute-scg-coverages` to your command line).")

            ###########################################################
            # TRNATaxonomyEstimatorMulti
            ###########################################################
            if self.__class__.__name__ in ['TRNATaxonomyEstimatorMulti']:
                if self.args.contigs_db or self.args.profile_db:
                    raise ConfigError("Taxonomy estimation classes have been initiated with files for metagenomes, but your arguments "
                                      "include also a single contigs or profile database path. You make anvi'o nervous. "
                                      "Please run this program either with a metagenomes file or contigs/profile databases.")

                if self.output_file_path:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file prefix rather "
                                      "than an output file path. Anvi'o will use your prefix and will generate many files that start "
                                      "with that prefix but ends with different names for each taxonomic level.")

                if not self.output_file_prefix:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file prefix :/")

                if self.raw_output and self.matrix_format:
                    raise ConfigError("Please don't request anvi'o to report the output both in raw and matrix format. Anvi'o shall "
                                      "not be confused :(")

                if self.output_file_prefix:
                    filesnpaths.is_output_file_writable(self.output_file_prefix)


class TRNATaxonomyArgs(object):
    def __init__(self, args, format_args_for_single_estimator=False):
        """A base class to fill in common arguments for SCG Taxonomy classes.

        The purpose of this class is to reduce the complexity of setting member variables
        for various classes in this module so we can get away with a single multi-talented
        sanity check base class without any complaints regarding missing member variables.

        Parameters
        ==========
        format_args_for_single_estimator: bool
            This is a special case where an args instance is generated to be passed to the
            single estimator from within multi estimator. More specifically, the multi estimator
            class is nothing but one that iterates through all metagenomes
            given to it using the single estimator class. So it needs to create instances of
            single estimators, and collect results at upstream. The problem is, if a single
            estimtor is initiated with the args of a multi estimator, the sanity check will
            go haywire. This flag nullifies most common offenders.
        """

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A('output_file')
        self.output_file_prefix = A('output_file_prefix')
        self.just_do_it = A('just_do_it')
        self.simplify_taxonomy_information = A('simplify_taxonomy_information')
        self.metagenome_mode = True if A('metagenome_mode') else False
        self.scg_name_for_metagenome_mode = A('scg_name_for_metagenome_mode')
        self.compute_scg_coverages = A('compute_scg_coverages')
        self.report_scg_frequencies_path = A('report_scg_frequencies')
        self.metagenomes = A('metagenomes')
        self.user_taxonomic_level = A('taxonomic_level')
        self.matrix_format = A('matrix_format')
        self.raw_output = A('raw_output')

        if format_args_for_single_estimator:
            # so you're here to get an args instance to fool a single estimator class.
            # very cute. we shall make that happen.
            self.metagenomes = None
            self.output_file_path = None
            self.output_file_prefix = None
            self.matrix_format = None
            self.raw_output = None

        self.skip_sanity_check = A('skip_sanity_check')


class TRNATaxonomyEstimatorMulti(TRNATaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        """Iterate through metagenome descriptions using TRNATaxonomyEstimatorSingle"""

        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        TRNATaxonomyArgs.__init__(self, self.args)

        # set your context
        self.ctx = ctx

        # intiate sanity check
        SanityCheck.__init__(self)

        self.metagenomes = None
        self.profile_dbs_available = False


    def init_metagenomes(self):
        self.progress.new("Initializing contigs DBs")
        self.progress.update("...")
        g = MetagenomeDescriptions(self.args, run=run_quiet, progress=self.progress)
        g.load_metagenome_descriptions()

        # NOTE some enforced flags here.
        self.compute_scg_coverages = g.profile_dbs_available
        if not self.metagenome_mode and g.profile_dbs_available:
            self.metagenome_mode = True

        metagenomes_without_scg_taxonomy = [m for m in g.metagenomes if not g.metagenomes[m]['scg_taxonomy_was_run']]
        if metagenomes_without_scg_taxonomy:
            if len(metagenomes_without_scg_taxonomy) == len(g.metagenomes):
                self.progress.end()
                raise ConfigError("Surprise! None of the %d genomes had no SCG taxonomy information." % len(g.metagenomes))
            else:
                self.progress.end()
                raise ConfigError("%d of your %d genomes has no SCG taxonomy information. Here is the list: '%s'." % \
                        (len(metagenomes_without_scg_taxonomy), len(g.metagenomes), ', '.join(metagenomes_without_scg_taxonomy)))

        # if we are here, it means SCGs were run for all metagenomes. here we will quickly check if versions agree
        # with each other and with the installed version of SCG taxonomy database
        trna_taxonomy_database_versions_in_metagenomes = [g.metagenomes[m]['trna_taxonomy_database_version'] for m in g.metagenomes]
        if len(set(trna_taxonomy_database_versions_in_metagenomes)) > 1:
            self.progress.reset()
            self.run.warning("Please note that not all SCG taxonomy database versions across your metagenomes are identical. "
                             "This means the program `anvi-run-trna-taxonomy` was run on these database across different versions of "
                             "the source SCG taxonomy database. This is OK and things will continue to work, but you should consider "
                             "the fact that taxonomy estimations coming from different versions of the database may not be comparable "
                             "anymore depending on what has changed between different versions of the database. If your purpose is not "
                             "to compare different versions of the database, and if you would like to ensure consistency, you can re-run "
                             "`anvi-run-trna-taxonomy` on contigs databases that have a different version than what is installed on your "
                             "system, which is '%s' (if you run `anvi-db-info` on any contigs database you can learn the SCG database "
                             "version of it). Anvi'o found these versions across your metagenomes: '%s'." % \
                                        (self.ctx.trna_taxonomy_database_version, ', '.join(list(set(trna_taxonomy_database_versions_in_metagenomes)))))
        elif trna_taxonomy_database_versions_in_metagenomes[0] != self.ctx.trna_taxonomy_database_version:
            self.progress.reset()
            self.run.warning("While all of your metagenomes agree with each other and have the SCG taxonomy database version of %s, "
                              "this version differs from what is installed on your system, which is %s. If you don't do anything, "
                              "things will continue to work. But if you would like to get rid of this warning you will need to "
                              "re-run the program `anvi-run-trna-taxonomy` on each one of them 😬" % \
                                        (trna_taxonomy_database_versions_in_metagenomes[0], self.ctx.trna_taxonomy_database_version))

        self.metagenomes = copy.deepcopy(g.metagenomes)
        self.metagenome_names = copy.deepcopy(g.metagenome_names)
        self.profile_dbs_available = g.profile_dbs_available

        self.progress.end()


    def estimate(self):
        if not self.metagenomes:
            self.init_metagenomes()
            self.run.info("Num metagenomes", len(self.metagenome_names))

        self.run.info("Taxonomic level of interest", self.user_taxonomic_level or "(None specified by the user, so 'all levels')")
        self.run.info("Output file prefix", self.output_file_prefix)
        self.run.info("Output in matrix format", self.matrix_format)
        self.run.info("Output raw data", self.raw_output)
        self.run.info("SCG coverages will be computed?", self.compute_scg_coverages)

        if self.report_scg_frequencies_path:
            self.report_scg_frequencies_as_TAB_delimited_file()
            return

        if not self.scg_name_for_metagenome_mode:
            self.scg_name_for_metagenome_mode = self.get_best_scg_name_for_metagenome_mode()

            self.run.warning("Please not that anvi'o just set the SCG for metagenome mode as '%s' since it was the most "
                             "frequent SCG occurring across all %d contigs databases involved in this analysis. But this "
                             "is nothing more than some heuristic for your convenience, and we strongly advice you to "
                             "run this program with the parameter `--report-scg-frequencies` and examine the output "
                             "to see if there is a better choice. As you can imagine, the most frequent SCG may not be "
                             "the one that is more common across all genomes or metagenomes you are interested." % \
                                                (self.scg_name_for_metagenome_mode, len(self.metagenomes)))

            self.run.info("SCG [determined by anvi'o]", self.scg_name_for_metagenome_mode, nl_after=1, mc="green")
        else:
            self.run.info("SCG [chosen by the user]", self.scg_name_for_metagenome_mode, nl_after=1, mc="green")

        scg_taxonomy_super_dict_multi = self.get_scg_taxonomy_super_dict_for_metagenomes()

        self.store_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)


    def get_print_friendly_scg_taxonomy_super_dict_multi(self, scg_taxonomy_super_dict_multi, as_data_frame=False):
        """Extract a more print-friendly data structure from `scg_taxonomy_super_dict_multi`

        Returns
        =======
        d: dict
            This is a dictionary that summarizes taxonomy and coverages per unique SCG in
            a given metagenome and looks like this:

            ----8<----8<----8<----8<----8<----8<----8<----8<----
            {
              "USA0001": {
                "Ribosomal_L16_69898": {
                  "gene_callers_id": 69898,
                  "gene_name": "Ribosomal_L16",
                  "accession": "CONSENSUS",
                  "percent_identity": "100.0",
                  "t_domain": "Bacteria",
                  "t_phylum": "Bacteroidota",
                  "t_class": "Bacteroidia",
                  "t_order": "Bacteroidales",
                  "t_family": "Rikenellaceae",
                  "t_genus": "Alistipes",
                  "t_species": "Alistipes shahii",
                  "tax_hash": "7e3cc7fc",
                  "coverages": {
                    "USA0001_01": 59.13119533527697
                  }
                },
                "Ribosomal_L16_55413": {
                  "gene_callers_id": 55413,
                  "gene_name": "Ribosomal_L16",
                  "accession": "CONSENSUS",
                  "percent_identity": "100.0",
                  "t_domain": "Bacteria",
                  "t_phylum": "Bacteroidota",
                  "t_class": "Bacteroidia",
                  "t_order": "Bacteroidales",
                  "t_family": "Tannerellaceae",
                  "t_genus": "Parabacteroides",
                  "t_species": null,
                  "tax_hash": "d16f9017",
                  "coverages": {
                    "USA0001_01": 48.82825484764543
                  }
                }

                (...)

            }
            ----8<----8<----8<----8<----8<----8<----8<----8<----

        """

        scg_taxonomy_super_dict_multi_print_friendly = {}

        for metagenome_name in scg_taxonomy_super_dict_multi:
            args = TRNATaxonomyArgs(self.args, format_args_for_single_estimator=True)
            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            if self.metagenome_mode:
                args.metagenome_mode = True
            else:
                args.metagenome_mode = False

            if self.profile_dbs_available:
                args.profile_db = self.metagenomes[metagenome_name]['profile_db_path']
                args.compute_scg_coverages = True

            d = TRNATaxonomyEstimatorSingle(args, run=run_quiet).get_print_friendly_scg_taxonomy_super_dict(scg_taxonomy_super_dict_multi[metagenome_name])

            # NOTE: what is happening down below might look stupid, because it really is. items of `d` here
            # contain a dictionary with the key 'coverages' where the coverage value of a given gene or bin
            # is kept per sample. to simplify downstream data wrangling operations, here we replace the
            # coverages dictionary with a single 'coverage': value entry to `d`. The
            # only reason we can do this here is becasue we only allow the user to use single profiles
            # with their metagenome descriptions and we are certain that there will be only a single entry in the
            # coverages dictionary.
            if self.compute_scg_coverages:
                for item in d:
                    coverages_list = list(d[item]['coverages'].items())

                    if len(coverages_list) > 1:
                        raise ConfigError("The codebase is not ready to handle this :(")

                    d[item]['coverage'] = coverages_list[0][1]
                    d[item].pop('coverages')

            scg_taxonomy_super_dict_multi_print_friendly[metagenome_name] = d

        if as_data_frame:
            return self.print_friendly_scg_taxonomy_super_dict_multi_to_data_frame(scg_taxonomy_super_dict_multi_print_friendly)
        else:
            return scg_taxonomy_super_dict_multi_print_friendly


    def print_friendly_scg_taxonomy_super_dict_multi_to_data_frame(self, scg_taxonomy_super_dict_multi_print_friendly):
        """Take a `scg_taxonomy_super_dict_multi_print_friendly`, and turn it into a neat data frame.

        Not working with dataframes from the get go sounds like a silly thing to do, but Python
        objects have been the way we passed around all the data (including to the interactive
        interface) so it feels weird to change that completely now. But obviously there are many
        advantages of working with dataframes, so the purpose of this function is to benefit from
        those in the context of multi estimator class.

        If the global `anvio.DEBUG` is True (via `--debug`), the dataframe is stored as a pickle
        object to play.

        Parameters
        ==========
        scg_taxonomy_super_dict_multi_print_friendly: dict
             This data structure is defined in `get_print_friendly_scg_taxonomy_super_dict_multi`

        Returns
        =======
        DF: pandas dataframe
            A neatly sorted dataframe. Looks like this:

            ----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<------
                metagenome_name                         taxon     coverage taxonomic_level
            0           USA0001                      Bacteria  1077.277012        t_domain
            1           USA0001               Unknown_domains     0.000000        t_domain
            2           USA0003                      Bacteria  2668.593534        t_domain
            3           USA0003               Unknown_domains     9.177677        t_domain
            4           USA0001                  Bacteroidota   932.478308        t_phylum
            5           USA0001                    Firmicutes   126.837011        t_phylum
            6           USA0001                Proteobacteria    17.961694        t_phylum
            7           USA0001                 Unknown_phyla     0.000000        t_phylum
            8           USA0001             Verrucomicrobiota     0.000000        t_phylum
            9           USA0003                    Firmicutes  1602.189160        t_phylum
            10          USA0003                  Bacteroidota  1034.471363        t_phylum
            11          USA0003             Verrucomicrobiota    22.668447        t_phylum
            12          USA0003                Proteobacteria     9.264563        t_phylum
            13          USA0003                 Unknown_phyla     9.177677        t_phylum
            14          USA0001                   Bacteroidia   932.478308         t_class
            15          USA0001                    Clostridia   116.441241         t_class
            16          USA0001           Gammaproteobacteria    17.961694         t_class
            17          USA0001                 Negativicutes    10.395770         t_class
            18          USA0001                       Bacilli     0.000000         t_class
            19          USA0001                 Lentisphaeria     0.000000         t_class
            ..              ...                           ...          ...             ...
            296         USA0003           CAG-302 sp001916775    11.435762       t_species
            297         USA0003   Paramuribaculum sp001689565    11.149206       t_species
            298         USA0003          CAG-1435 sp000433775    11.078189       t_species
            299         USA0003           CAG-269 sp001916005     9.848404       t_species
            300         USA0003            UBA737 sp002431945     7.497512       t_species
            301         USA0003           CAG-345 sp000433315     6.733333       t_species
            302         USA0003           UBA1829 sp002405835     6.684840       t_species
            303         USA0003           TF01-11 sp001916135     6.060109       t_species
            304         USA0003       Barnesiella sp002161555     0.000000       t_species
            305         USA0003           CAG-115 sp000432175     0.000000       t_species
            306         USA0003           CAG-452 sp000434035     0.000000       t_species
            307         USA0003       Duncaniella sp002494015     0.000000       t_species
            308         USA0003   Paramuribaculum sp001689535     0.000000       t_species
            309         USA0003  Succiniclasticum sp002342505     0.000000       t_species
            ----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<------

        """

        self.progress.new("Data dict to dataframe")
        self.progress.update("Bleep very complex stuff bloop")

        taxonomic_levels = [self.user_taxonomic_level] if self.user_taxonomic_level else self.ctx.levels_of_taxonomy

        df = pd.DataFrame.from_dict({(i,j): scg_taxonomy_super_dict_multi_print_friendly[i][j]
                                for i in scg_taxonomy_super_dict_multi_print_friendly.keys()
                                for j in scg_taxonomy_super_dict_multi_print_friendly[i].keys()}, orient='index')

        df.reset_index(inplace=True)
        df.rename(columns={"level_0": "metagenome_name"}, inplace=True)
        df.drop(['level_1'], axis=1, inplace=True)

        if not self.compute_scg_coverages:
            # NOTE: this is a bit critical. if the user are working with external genomes where there is no
            # coverage information is availble for SCGs, the following code will explode for obvious
            # reasons (which include the fact that `scg_taxonomy_super_dict_multi` will not have an entry
            # for `coverage`). What we will do here is a little trick. First, we will add a coverage column
            # to `df`. And at the very end, we will replace it with `times_observed`.
            df["coverage"] = 1

        DFx = pd.DataFrame(columns=['metagenome_name', 'taxon', 'coverage', 'taxonomic_level'])
        for taxonomic_level in taxonomic_levels:
            x = df.fillna(constants.levels_of_taxonomy_unknown).groupby(['metagenome_name', taxonomic_level])['coverage'].sum().reset_index()
            x.columns = ['metagenome_name', 'taxon', 'coverage']
            x['taxonomic_level'] = taxonomic_level

            # adding for each metagenome the missing taxon names
            list_of_unique_taxon_names = x['taxon'].unique()
            for metagenome_name in x['metagenome_name'].unique():
                known_taxa = set(x[x['metagenome_name'] == metagenome_name]['taxon'])
                missing_taxa = [taxon for taxon in list_of_unique_taxon_names if taxon not in known_taxa]
                for taxon in missing_taxa:
                    row = {'metagenome_name': metagenome_name, 'taxon': taxon, 'coverage': 0.0, 'taxonomic_level': taxonomic_level}
                    x = x.append(row, ignore_index=True)

            DFx = DFx.append(x, ignore_index=True)

        DFx.sort_values(by=['metagenome_name', 'taxonomic_level', 'coverage'], ascending=[True, True, False], inplace=True)

        DF = pd.DataFrame(columns=['metagenome_name', 'taxon', 'coverage', 'taxonomic_level'])
        for taxonomic_level in constants.levels_of_taxonomy:
            DF = DF.append(DFx[DFx['taxonomic_level'] == taxonomic_level], ignore_index=True)

        del DFx

        if anvio.DEBUG:
            import pickle
            with open('DataFrame.pickle', 'wb') as output:
                pickle.dump(DF, output)

            self.progress.reset()
            self.run.info_single("The `--debug` flag made anvi'o to dump the data frame that was generated by the "
                                 "function `print_friendly_scg_taxonomy_super_dict_multi_to_data_frame` (lol) as "
                                 "a picke object stored in file name 'DataFrame.pickle'. If you would like to play "
                                 "with it, you can start an ipython session, and run the following: import pandas "
                                 "as pd; import pickle; df = pickle.load(open('DataFrame.pickle', 'rb'))",
                                 nl_before=1, nl_after=1)

        if not self.compute_scg_coverages:
            DF.rename(columns={"coverage": "times_observed"}, inplace = True)
            DF['times_observed'] = DF['times_observed'].astype(int)

        self.progress.end()

        return DF


    def store_scg_taxonomy_super_dict_multi(self, scg_taxonomy_super_dict_multi):
        if self.raw_output:
            self.store_scg_taxonomy_super_dict_raw(scg_taxonomy_super_dict_multi)
        else:
            if self.matrix_format:
                self.store_scg_taxonomy_super_dict_multi_matrix_format(scg_taxonomy_super_dict_multi)
            else:
                self.store_scg_taxonomy_super_dict_multi_long_format(scg_taxonomy_super_dict_multi)


    def store_scg_taxonomy_super_dict_multi_long_format(self, scg_taxonomy_super_dict_multi):
        df = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi, as_data_frame=True)

        output_file_path = self.output_file_prefix + '-LONG-FORMAT.txt'
        df.to_csv(output_file_path, index=True, index_label="entry_id", sep='\t')

        self.run.info("Long-format output", output_file_path)


    def store_scg_taxonomy_super_dict_raw(self, scg_taxonomy_super_dict_multi):
        d = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)

        taxonomic_levels = [self.user_taxonomic_level] if self.user_taxonomic_level else self.ctx.levels_of_taxonomy

        header = ['metagenome_name', 'gene_name', 'gene_callers_id', 'percent_identity']

        if self.compute_scg_coverages:
            header += ['coverage']

        header += taxonomic_levels

        output_file_path = self.output_file_prefix + '-RAW-LONG-FORMAT.txt'
        with open(output_file_path, 'w') as output:
            output.write('\t'.join(header) + '\n')

            for metagenome_name in d:
                for gene_name in d[metagenome_name]:
                    output.write('\t'.join([metagenome_name] + [str(d[metagenome_name][gene_name][h]) for h in header[1:]]) + '\n')

        self.run.info("Raw output", output_file_path)


    def store_scg_taxonomy_super_dict_multi_matrix_format(self, scg_taxonomy_super_dict_multi):
        df = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi, as_data_frame=True)

        taxonomic_levels = [self.user_taxonomic_level] if self.user_taxonomic_level else self.ctx.levels_of_taxonomy

        for taxonomic_level in taxonomic_levels:
            dfx = df[df['taxonomic_level'] == taxonomic_level]
            dfx.set_index(['metagenome_name', 'taxon'], inplace=True)

            if self.compute_scg_coverages:
                matrix = sps.coo_matrix((dfx.coverage, (dfx.index.codes[0], dfx.index.codes[1]))).todense().tolist()
            else:
                matrix = sps.coo_matrix((dfx.times_observed, (dfx.index.codes[0], dfx.index.codes[1]))).todense().tolist()

            rows = dfx.index.levels[0].tolist()
            cols = ['taxon'] + dfx.index.levels[1].tolist()

            output_file_path = '%s-%s-MATRIX.txt' % (self.output_file_prefix, taxonomic_level)
            temp_file_path = filesnpaths.get_temp_file_path()
            with open(temp_file_path, 'w') as output:
                output.write('\t'.join(cols) + '\n')
                for i in range(0, len(matrix)):
                    output.write('\t'.join([rows[i]] + ['%.2f' % c for c in matrix[i]]) + '\n')

            utils.transpose_tab_delimited_file(temp_file_path, output_file_path, remove_after=True)

            self.run.info('Output matrix for "%s"' % taxonomic_level, output_file_path)


    def report_scg_frequencies_as_TAB_delimited_file(self):
        scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies = self.get_scg_frequencies()

        utils.store_dict_as_TAB_delimited_file(scg_frequencies, self.report_scg_frequencies_path, headers=['genome'] + scgs_ordered_based_on_frequency, keys_order=contigs_dbs_ordered_based_on_num_scgs)

        self.run.info('SCG frequencies across contigs dbs', self.report_scg_frequencies_path)


    def get_best_scg_name_for_metagenome_mode(self):
        """Identifies the best SCG to use for taxonomy.

        Although in reality 'best' means pretty much nothing here as this
        function simply returns the most frequent SCG. In an ideal world
        there would be mulitple modes, one that takes into consideration
        the prevalence of each SCG across each `self.genome`. If we ever
        implement that, let's leave a FIXME here only to be remove by
        that hero.
        """

        scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies = self.get_scg_frequencies()

        return scgs_ordered_based_on_frequency[0]


    def get_scg_taxonomy_super_dict_for_metagenomes(self):
        """Generates the taxonomy super dict, the primary data structure for this class.

        The difference between this function and `get_scg_taxonomy_super_dict` in `TRNATaxonomyEstimatorSingle`
        is that this one aggregates multiple `scg_taxonomy_super_dict` instances into a larger dictionary.

        Returns
        =======
        scg_taxonomy_super_dict: dict
            This is a complex dictionary that looks like this:

            ----8<----8<----8<----8<----8<----8<----8<----8<----
                {
                  "USA0001": { # <- name given for this entry in internal/external genomes file
                    "taxonomy": {
                      "USA0001_01": {  # <- project name in contigs database
                        "scgs": {
                          "69898": {
                            "gene_callers_id": 69898,
                            "gene_name": "Ribosomal_L16",
                            "accession": "CONSENSUS",
                            "percent_identity": "100.0",
                            "t_domain": "Bacteria",
                            "t_phylum": "Bacteroidota",
                            "t_class": "Bacteroidia",
                            "t_order": "Bacteroidales",
                            "t_family": "Rikenellaceae",
                            "t_genus": "Alistipes",
                            "t_species": "Alistipes shahii",
                            "tax_hash": "7e3cc7fc"
                          },
                          "55413": {
                            "gene_callers_id": 55413,
                            "gene_name": "Ribosomal_L16",
                            "accession": "CONSENSUS",
                            "percent_identity": "100.0",
                            "t_domain": "Bacteria",
                            "t_phylum": "Bacteroidota",
                            "t_class": "Bacteroidia",
                            "t_order": "Bacteroidales",
                            "t_family": "Tannerellaceae",
                            "t_genus": "Parabacteroides",
                            "t_species": null,
                            "tax_hash": "d16f9017"
                          }
                        },
                        "metagenome_mode": true
                      }
                    },
                    "coverages": {
                      "69898": {
                        "USA0001_01": 18.146179401993354
                      },
                      "55413": {
                        "USA0001_01": 48.82825484764543
                      }
                    }
                  },

                  (...)

                }
            ----8<----8<----8<----8<----8<----8<----8<----8<----

            Coverages will be there only if `--compute-scg-coverages` flag was passed to the class.
        """

        scg_taxonomy_super_dict = {}

        if self.profile_dbs_available:
            self.run.info_single("Your metagenome file DOES contain profile databases, so anvi'o will turn on `--metagenome-mode`, "
                                 "set the SCG name to %s, and turn on `--compute-scg-coverages` flag. If this doesn't make sense, "
                                 "please adjust your input parameters." % (self.scg_name_for_metagenome_mode), nl_after=1)
        else:
            if self.metagenome_mode:
                self.run.info_single("Your metagenome file DOES NOT contain profile databases, but you asked anvi'o to estimate SCG "
                                     "taxonomy in metagenome mode. So be it. SCG name is set to %s." \
                                        % (self.scg_name_for_metagenome_mode), nl_after=1)
            else:
                self.run.info_single("Your metagenome file DOES NOT contain profile databases, and you haven't asked anvi'o to "
                                     "work in `--metagenome-mode`. Your contigs databases will be treated as geomes rather than "
                                     "metagenomes.", nl_after=1)

        self.progress.new("Recovering tax super dict", progress_total_items=len(self.metagenome_names))
        total_num_metagenomes = len(self.metagenome_names)
        for metagenome_name in self.metagenome_names:
            args = TRNATaxonomyArgs(self.args, format_args_for_single_estimator=True)

            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            if self.profile_dbs_available:
                args.metagenome_mode = True
                args.profile_db = self.metagenomes[metagenome_name]['profile_db_path']
                args.compute_scg_coverages = True
                args.scg_name_for_metagenome_mode = self.scg_name_for_metagenome_mode
            else:
                if self.metagenome_mode:
                    args.metagenome_mode = True
                    args.scg_name_for_metagenome_mode = self.scg_name_for_metagenome_mode
                else:
                    args.metagenome_mode = False
                    args.scg_name_for_metagenome_mode = None

            self.progress.update("[%d of %d] %s" % (self.progress.progress_current_item + 1, total_num_metagenomes, metagenome_name))
            scg_taxonomy_super_dict[metagenome_name] = TRNATaxonomyEstimatorSingle(args, progress=progress_quiet, run=run_quiet).get_scg_taxonomy_super_dict()

            self.progress.increment()

        self.progress.end()

        return scg_taxonomy_super_dict


    def get_scg_frequencies(self):
        """Calculates the SCG frequencies across `self.metagenomes`."""

        scg_frequencies = {}

        self.progress.new("Learning SCG frequencies across genomes", progress_total_items=len(self.metagenomes))
        for metagenome_name in self.metagenomes:
            self.progress.update("%s ..." % metagenome_name)
            scg_frequencies[metagenome_name] = {}

            args = TRNATaxonomyArgs(self.args, format_args_for_single_estimator=True)
            args.compute_scg_coverages = False
            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            e = TRNATaxonomyEstimatorSingle(args, progress=progress_quiet, run=run_quiet)
            for scg_name in self.ctx.default_anticodons_for_taxonomy:
                scg_frequencies[metagenome_name][scg_name] = e.frequency_of_scgs_with_taxonomy[scg_name]

            self.progress.increment()

        self.progress.update("Finalizing data")

        scg_frequencies_across_contigs_dbs = [(scg_name, sum([scg_frequencies[genome_name][scg_name] for genome_name in scg_frequencies])) for scg_name in self.ctx.default_anticodons_for_taxonomy]
        scgs_ordered_based_on_frequency = [frequency_tuple[0] for frequency_tuple in sorted(scg_frequencies_across_contigs_dbs, key = lambda x: x[1], reverse=True)]

        num_scgs_for_each_contigs_db = [(genome_name, sum(scg_frequencies[genome_name].values())) for genome_name in scg_frequencies]
        contigs_dbs_ordered_based_on_num_scgs = [frequency_tuple[0] for frequency_tuple in sorted(num_scgs_for_each_contigs_db, key = lambda x: x[1], reverse=True)]

        self.progress.end()

        return scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies


class TRNATaxonomyEstimatorSingle(TRNATaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.update_profile_db_with_taxonomy = A('update_profile_db_with_taxonomy')
        self.bin_id = A('bin_id')

        TRNATaxonomyArgs.__init__(self, self.args)

        self.ctx = ctx

        SanityCheck.__init__(self)

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Profile DB', self.profile_db_path)
        self.run.info('Metagenome mode', self.metagenome_mode)
        if self.metagenome_mode:
            self.run.info('SCG for metagenome', self.scg_name_for_metagenome_mode)

        # these dictionaries that will be initiated later
        self.contigs_db_project_name = "Unknown"
        self.scg_name_to_gene_caller_id_dict = {}
        self.frequency_of_scgs_with_taxonomy = {}
        self.gene_callers_id_to_scg_taxonomy_dict = {}
        self.split_name_to_gene_caller_ids_dict = {}
        self.gene_callers_id_to_split_name_dict = {}
        self.sample_names_in_profile_db = None

        self.initialized = False

        if not skip_init:
            self.init()


    def init(self):
        self.init_scg_data()

        if self.report_scg_frequencies_path:
            with open(self.report_scg_frequencies_path, 'w') as output:
                for scg_name, frequency in self.frequency_of_scgs_with_taxonomy.items():
                    output.write("%s\t%d\n" % (scg_name, frequency))

            self.run.info('SCG frequencies within the contigs db', self.report_scg_frequencies_path, nl_before=1)
            sys.exit()

        if self.profile_db_path:
            self.sample_names_in_profile_db = ProfileDatabase(self.profile_db_path).samples

        self.initialized = True


    def init_scg_data(self):
        """Initialize SCG taxonomy for the entire contigs database"""

        if not self.contigs_db_path:
            return None

        self.progress.new('Initializing')
        self.progress.update('SCG taxonomy dictionary')

        for scg_name in self.ctx.anticodons:
            self.scg_name_to_gene_caller_id_dict[scg_name] = set([])

        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        self.contigs_db_project_name = contigs_db.meta['project_name']
        scg_taxonomy_table = contigs_db.db.get_table_as_dict(t.scg_taxonomy_table_name)
        genes_in_splits = contigs_db.db.get_some_columns_from_table(t.genes_in_splits_table_name, "split, gene_callers_id")
        min_contig_length_in_contigs_db = contigs_db.db.get_max_value_in_column(t.contigs_info_table_name, "length", return_min_instead=True)
        contigs_db.disconnect()

        # this is important. before we begin, we need to filter out gene caller ids and splits from main dictionaries if
        # they shouldn't be there. read the warning below to see the utility of this step.
        if self.profile_db_path and self.compute_scg_coverages:
            split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.profile_db_path))
            split_names_in_contigs_db = set([tpl[0] for tpl in genes_in_splits])
            splits_missing_in_profile_db = split_names_in_contigs_db.difference(split_names_in_profile_db)

            min_contig_length_in_profile_db = ProfileDatabase(self.profile_db_path).meta['min_contig_length']

            if len(splits_missing_in_profile_db):
                self.progress.reset()
                self.run.warning("Please note that anvi'o found %s splits in your contigs database. But only %s of them "
                                 "appeared ot be in the profile database. As a result, anvi'o will now remove the %s splits "
                                 "that occur only in the contigs db from all downstream analyses here (if you didn't use the flag "
                                 "`--compute-scg-coverages` this wouldn't have been necessary, but with the current settings "
                                 "this is really the best for everyone). Where is this difference coming from though? Well. This "
                                 "is often the case because the 'minimum contig length parameter' set during the `anvi-profile` "
                                 "step can exclude many contigs from downstream analyses (often for good reasons, too). For "
                                 "instance, in your case the minimum contig length goes as low as %s nts in your contigs database. "
                                 "Yet, the minimum contig length set in the profile databaes is %s nts. Hence the difference. Anvi'o "
                                 "hopes that this explaines some things." % (pp(len(split_names_in_contigs_db)),
                                                                             pp(len(split_names_in_profile_db)),
                                                                             pp(len(splits_missing_in_profile_db)),
                                                                             pp(min_contig_length_in_contigs_db),
                                                                             pp(min_contig_length_in_profile_db)))

                self.progress.update("Removing %s splits missing form the profile db" % pp(len(splits_missing_in_profile_db)))
                genes_in_splits = [tpl for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db]

                # so now we know the final list of split names and gene caller ids as they are stored in the updated
                # `genes_in_splits` variable. time to clean up the `scg_taxonomy_table` dictionary as well.
                final_set_of_gene_caller_ids = set([tpl[1] for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db])
                entry_ids_to_remove = [entry for entry in scg_taxonomy_table if scg_taxonomy_table[entry]['gene_callers_id'] not in final_set_of_gene_caller_ids]
                [scg_taxonomy_table.pop(e) for e in entry_ids_to_remove]

        # NOTE: This will modify the taxonomy strings read from the contigs database. see the
        # function header for `trim_taxonomy_dict_entry` for more information.
        if self.simplify_taxonomy_information:
            self.progress.update('SCG taxonomy dicts ... trimming main')
            for key in scg_taxonomy_table:
                scg_taxonomy_table[key] = self.trim_taxonomy_dict_entry(scg_taxonomy_table[key])

        self.progress.update('SCG taxonomy dicts ... building g->tax')
        for entry in scg_taxonomy_table.values():
            gene_callers_id = entry['gene_callers_id']
            self.gene_callers_id_to_scg_taxonomy_dict[gene_callers_id] = entry

        self.progress.update('SCG taxonomy dicts ... building tax->g')
        for entry in self.gene_callers_id_to_scg_taxonomy_dict.values():
            scg_gene_name = entry['gene_name']
            gene_callers_id = entry['gene_callers_id']
            self.scg_name_to_gene_caller_id_dict[scg_gene_name].add(gene_callers_id)

        self.progress.update('SCG taxonomy dicts ... building s->g + g->s')
        for split_name, gene_callers_id in genes_in_splits:
            if gene_callers_id not in self.gene_callers_id_to_scg_taxonomy_dict:
                continue

            if split_name not in self.split_name_to_gene_caller_ids_dict:
                self.split_name_to_gene_caller_ids_dict[split_name] = set()

            self.split_name_to_gene_caller_ids_dict[split_name].add(gene_callers_id)
            self.gene_callers_id_to_split_name_dict[gene_callers_id] = split_name

        self.progress.end()

        self.frequency_of_scgs_with_taxonomy = OrderedDict(sorted([(g, len(self.scg_name_to_gene_caller_id_dict[g])) for g in self.scg_name_to_gene_caller_id_dict], key = lambda x: x[1], reverse=True))

        if self.metagenome_mode or anvio.DEBUG:
            self.run.info_single("A total of %s single-copy core genes with taxonomic affiliations were successfully initialized "
                                 "from the contigs database 🎉 Following shows the frequency of these SCGs: %s." % \
                                            (pp(len(self.gene_callers_id_to_scg_taxonomy_dict)),
                                             ', '.join(["%s (%d)" % (g, self.frequency_of_scgs_with_taxonomy[g]) \
                                                                for g in self.frequency_of_scgs_with_taxonomy])), nl_before=1)


    def trim_taxonomy_dict_entry(self, taxonomy_dict_entry):
        """ Remove excess information from taxonomy information.

        The purpose of this is to give an option to the user to simplify GTDB names, that
        will have a text information for every level of taxonomy depending on what branches
        genomes fit, but it is not always helpful to the user. Such as this one:

             t_domain Bacteria
             t_phylum Firmicutes
             t_class Clostridia
             t_order Monoglobales
             t_family UBA1381
             t_genus CAG-41
             t_species CAG-41 sp900066215

         in this case the user may want to get this instead:

             t_domain Bacteria
             t_phylum Firmicutes
             t_class Clostridia
             t_order Monoglobales
             t_family None
             t_genus None
             t_species None

         So this function will take a taxonomy dict entry , and will return a simplified
         version of it if trimming is applicable.

        Paremeters
        ==========
        taxonomy_dict_entry: dict
            a dictionary that contains keys for all taxon names. such as this one:
                {[...],
                 't_domain': 'Bacteria',
                 't_phylum': 'Firmicutes',
                 't_class': 'Clostridia',
                 't_order': 'Oscillospirales',
                 't_family': 'Acutalibacteraceae',
                 't_genus': 'Ruminococcus',
                 't_species': 'Ruminococcus sp002491825'
                }
         """

        # for optimization, these letters should all have three characters. if that
        # behavior needs to change, the code down below must be updates. the purpose
        # of this is not to have a comprehensive list of EVERY single GTDB-specific
        # clade designations, but to make sure teh vast majority of names are covered.
        GTDB_specific_clade_prefixes = ['CAG', 'GCA', 'UBA', 'FUL', 'PAL', '2-0', 'Fen', 'RF3', 'TAN']

        taxonomic_levels_to_nullify = []

        for taxonomic_level in self.ctx.levels_of_taxonomy[::-1]:
            if not taxonomy_dict_entry[taxonomic_level]:
                continue

            if taxonomic_level == 't_species':
                species_name = taxonomy_dict_entry[taxonomic_level].split(' ')[1]
                try:
                    int(species_name[2])
                    taxonomic_levels_to_nullify.append(taxonomic_level)
                except:
                    None
            else:
                if taxonomy_dict_entry[taxonomic_level][0:3] in GTDB_specific_clade_prefixes:
                    taxonomic_levels_to_nullify.append(taxonomic_level)

        # this is the best way to make sure we are not going to nullify order, but leave behind a family name.
        if taxonomic_levels_to_nullify:
            level_below_which_to_nullify = min([self.ctx.levels_of_taxonomy.index(l) for l in taxonomic_levels_to_nullify])
            for taxonomic_level in self.ctx.levels_of_taxonomy[level_below_which_to_nullify:]:
                taxonomy_dict_entry[taxonomic_level] = None

        return taxonomy_dict_entry


    def get_blank_hit_template_dict(self):
        hit = {}

        for level in self.ctx.levels_of_taxonomy[::-1]:
            hit[level] = None

        return hit


    def get_consensus_taxonomy(self, scg_taxonomy_dict):
        """Takes in a scg_taxonomy_dict, returns a final taxonomic string that summarize all"""

        if not len(scg_taxonomy_dict):
            return dict([(l, None) for l in self.ctx.levels_of_taxonomy])

        pd.set_option('mode.chained_assignment', None)

        scg_hits = list([v for v in scg_taxonomy_dict.values() if v['t_domain']])

        if not len(scg_hits):
            return self.get_blank_hit_template_dict()

        df = pd.DataFrame.from_records(scg_hits)

        # we have already stored a unique hash for taxonomy strings. here we will figure out most frequent
        # hash values in the df
        tax_hash_counts = df['tax_hash'].value_counts()
        tax_hash_df = tax_hash_counts.rename_axis('tax_hash').reset_index(name='frequency')
        max_frequency = tax_hash_df.frequency.max()
        tax_hash_df_most_frequent = tax_hash_df[tax_hash_df.frequency == max_frequency]

        if len(tax_hash_df_most_frequent.index) == 1:
            # if there is only a single winner, we're golden
            winner_tax_hash = tax_hash_df_most_frequent.tax_hash[0]

            # get the consensus hit based on the winner hash
            consensus_hit = df[df.tax_hash == winner_tax_hash].head(1)

            # turn it into a Python dict before returning
            return consensus_hit.to_dict('records')[0]
        else:
            # if there are competing hashes, we need to be more careful to decide
            # which taxonomic level should we use to cut things off.
            consensus_hit = self.get_blank_hit_template_dict()
            for level in self.ctx.levels_of_taxonomy[::-1]:
                if len(df[level].unique()) == 1:
                    consensus_hit[level] = df[level].unique()[0]

            return consensus_hit


    def print_scg_taxonomy_hits_in_splits(self, hits, bin_name=None):
        self.progress.reset()
        self.run.warning(None, header='Hits for %s' % (bin_name if bin_name else "a bunch of splits"), lc="green")

        if len(hits):
            header = ['SCG', 'gene', 'pct id', 'taxonomy']
            table = []

            for hit in hits:
                taxon_text = ' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy])

                # if the hit we are working on sent here as 'consensus', we will color it up a bit so it shows up
                # more clearly in the debug output.
                if hit['gene_name'] == 'CONSENSUS':
                    taxon_text = terminal.c(taxon_text, color='red')

                    for field_name in ['gene_name', 'percent_identity', 'gene_callers_id']:
                        hit[field_name] = terminal.c(hit[field_name], color='red')

                table.append([hit['gene_name'], str(hit['gene_callers_id']), str(hit['percent_identity']), taxon_text])

            anvio.TABULATE(table, header)
        else:
            self.run.info_single("No hits :/")


    def get_scg_taxonomy_dict(self, gene_caller_ids, bin_name=None):
        scg_taxonomy_dict = {}

        improper_gene_caller_ids = [g for g in gene_caller_ids if g not in self.gene_callers_id_to_scg_taxonomy_dict]
        if improper_gene_caller_ids:
            raise ConfigError("Something weird is going on. Somehow anvi'o has a bunch of gene caller ids for which it is "
                              "supposed to estimate taxonomy. However, %d of them do not occur in a key dictionary. The code "
                              "here does not know what to suggest :( Apologies." % len(improper_gene_caller_ids))

        for gene_callers_id in gene_caller_ids:
            scg_taxonomy_dict[gene_callers_id] = self.gene_callers_id_to_scg_taxonomy_dict[gene_callers_id]
            scg_taxonomy_dict[gene_callers_id]["tax_hash"] = HASH(self.gene_callers_id_to_scg_taxonomy_dict[gene_callers_id])

        return scg_taxonomy_dict


    def estimate_for_list_of_splits(self, split_names=None, bin_name=None):
        """Estimate SCG taxonomy for a bunch of splits that belong to a single population.

           The purpose of this function is to to do critical things: identify SCGs we use for taxonomy in `split_names`,
           and generate a consensus taxonomy with the assumption that these are coming from splits that represents a
           single population.

           It will return a dictionary with multiple items, including a dictionary that contains the final consensus\
           taxonomy, another one that includes every SCG and their raw associations with taxon names (from which the\
           consensus taxonomy was computed), as well as information about how many SCGs were analyzed and supported the\
           consesnus.
        """

        if self.metagenome_mode:
            raise ConfigError("Someone is attempting to estimate taxonomy for a set of splits using a class inherited in "
                              "`metagenome mode`. If you are a programmer please note that it is best to use the member "
                              "function `estimate` directly.")

        consensus_taxonomy = None

        gene_caller_ids_of_interest = self.get_gene_caller_ids_for_splits(split_names)
        scg_taxonomy_dict = self.get_scg_taxonomy_dict(gene_caller_ids_of_interest)

        try:
            consensus_taxonomy = self.get_consensus_taxonomy(scg_taxonomy_dict)
            consensus_taxonomy['gene_name'] = 'CONSENSUS'
            consensus_taxonomy['percent_identity'] = '--'
            consensus_taxonomy['gene_callers_id'] = '--'

        except Exception as e:
            self.print_scg_taxonomy_hits_in_splits(list(scg_taxonomy_dict.values()))

            raise ConfigError("While trying to sort out the consensus taxonomy for %s anvi'o failed :( The list of SCG taxon hits that "
                              "caused the failure is printed in your terminal. But the actual error message that came from the depths "
                              "of the codebase was this: '%s'." % (('the bin "%s"' % bin_name) if bin_name else 'a bunch of splits', e))

        if anvio.DEBUG:
            self.print_scg_taxonomy_hits_in_splits(list(scg_taxonomy_dict.values()) + [consensus_taxonomy], bin_name)

        # set some useful information. `total_scgs` is the number of SCGs with taxonomy found in the collection of splits. the
        # `supporting_scgs` communicate how many of them supports the consensus taxonomy fully
        total_scgs = len(scg_taxonomy_dict)
        supporting_scgs = 0

        consensus_taxonomy_levels_occupied = [level for level in self.ctx.levels_of_taxonomy if consensus_taxonomy[level]]
        consensus_taxonomy_str = ' / '.join([consensus_taxonomy[level] for level in consensus_taxonomy_levels_occupied])

        for scg_taxonomy_hit in scg_taxonomy_dict.values():
            scg_taxonomy_hit_str = ' / '.join([str(scg_taxonomy_hit[level]) for level in consensus_taxonomy_levels_occupied])

            if scg_taxonomy_hit_str == consensus_taxonomy_str:
                scg_taxonomy_hit['supporting_consensus'] = True
                supporting_scgs += 1
            else:
                scg_taxonomy_hit['supporting_consensus'] = False

        return {'consensus_taxonomy': consensus_taxonomy,
                'scgs': scg_taxonomy_dict,
                'total_scgs': total_scgs,
                'supporting_scgs': supporting_scgs,
                'metagenome_mode': False}


    def estimate_for_bins_in_collection(self):
        bins_taxonomy_dict = {}

        bin_name_to_split_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.run.info_single("%s split names associated with %s bins of in collection '%s' have been "
                             "successfully recovered 🎊" % (pp(sum([len(v) for v in bin_name_to_split_names_dict.values()])),
                                                           pp(len(bin_name_to_split_names_dict)),
                                                           self.collection_name), nl_before=1)

        for bin_name in bin_name_to_split_names_dict:
            split_names = bin_name_to_split_names_dict[bin_name]
            bins_taxonomy_dict[bin_name] = self.estimate_for_list_of_splits(split_names, bin_name)

        return bins_taxonomy_dict


    def estimate_for_contigs_db_for_genome(self):
        contigs_db_taxonomy_dict = {}

        scg_frequencies = self.frequency_of_scgs_with_taxonomy.values()
        if len([sf for sf in scg_frequencies if sf > 1]) * 100 / len(scg_frequencies) > 20:
            if self.just_do_it:
                self.run.warning("Because you asked anvi'o to just do it, it will do it, but you seem to have too much contamination "
                                 "in this contigs database for it to represent a genome. So probably taxonomy estimations are all "
                                 "garbage, but hey, at least it runs?")
            else:
                raise ConfigError("Because you haven't used the `--metagenome-mode` flag, anvi'o was trying to treat your contigs "
                                  "database as a genome. But there seems to be too much redundancy of single-copy core genes in this "
                                  "contigs database to assign taxonomy with any confidence :/ A more proper way to do this is to use the "
                                  "`--metagenome-mode` flag. Or you can also tell anvi'o to `--just-do-it`. It is your computer after "
                                  "all :( But you should still be aware that in that case you would likely get a completely irrelevant "
                                  "answer from this program.")

        splits_in_contigs_database = self.split_name_to_gene_caller_ids_dict.keys()
        contigs_db_taxonomy_dict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(split_names=splits_in_contigs_database,
                                                                                                  bin_name=self.contigs_db_project_name)
        return contigs_db_taxonomy_dict


    def estimate_for_contigs_db_for_metagenome(self):
        """Treat a given contigs database as a metagenome.

           This function deserves some attention. It relies on a single SCG to estimate the composition of a metagenome.
           For instance, its sister function, `estimate_for_contigs_db_for_genome`, works with a list of splits that are
           assumed to belong to the same genome. In which case a consensus taxonomy learned from all SCGs is most
           appropriate. In this case, however, we don't know which split will go together, hence, we can't pull together
           SCGs to learn a consensus taxonomy for independent populations in the metagenome. The best we can do is to stick
           with a single SCG with the hope that (1) it will cut through as many populations as possible and (2) will have
           reasonable power to resolve taxonomy all by itself. These independent assumptions will both work in some cases
           and both fail in others.
        """

        # we first need to decide which SCG we should use to survey taxonomy
        most_frequent_scg = next(iter(self.frequency_of_scgs_with_taxonomy))
        if self.scg_name_for_metagenome_mode:
            frequency_of_user_chosen_scg = self.frequency_of_scgs_with_taxonomy[self.scg_name_for_metagenome_mode]
            frequency_of_most_frequent_scg = self.frequency_of_scgs_with_taxonomy[most_frequent_scg]

            if frequency_of_user_chosen_scg < frequency_of_most_frequent_scg:
                additional_note = " And just so you know, there is another SCG that was observed more times (i.e., %s; %d times)\
                                   in this metagenome compared to yours (i.e., %d times). You're the boss, of course." %\
                                            (most_frequent_scg, frequency_of_most_frequent_scg, frequency_of_user_chosen_scg)
            else:
                additional_note = ""

            self.run.warning("As per your request anvi'o set '%s' to be THE single-copy core gene to survey your metagenome for its "
                             "taxonomic composition.%s" % (self.scg_name_for_metagenome_mode, additional_note))
        else:
            self.scg_name_for_metagenome_mode = most_frequent_scg

            self.run.warning("Anvi'o automatically set '%s' to be THE single-copy core gene to survey your metagenome for its "
                             "taxonomic composition. If you are not happy with that, you could change it with the parameter "
                             "`--scg-name-for-metagenome-mode`." % (self.scg_name_for_metagenome_mode))

        gene_caller_ids_of_interest = self.scg_name_to_gene_caller_id_dict[self.scg_name_for_metagenome_mode]
        scg_taxonomy_dict = self.get_scg_taxonomy_dict(gene_caller_ids=gene_caller_ids_of_interest,
                                                       bin_name=self.contigs_db_project_name)

        return {self.contigs_db_project_name: {'scgs': scg_taxonomy_dict,
                                               'metagenome_mode': True}}


    def get_scg_taxonomy_super_dict(self):
        """Function that returns the `scg_taxonomy_super_dict`.

           `scg_taxonomy_super_dict` contains a wealth of information regarding samples, SCGs,
           SCG taxonomic affiliations, consensus taxonomy, and coverages of SCGs across samples.
        """
        scg_taxonomy_super_dict = {}

        if not self.initialized:
            self.init()

        if self.profile_db_path and not self.metagenome_mode:
            scg_taxonomy_super_dict['taxonomy'] = self.estimate_for_bins_in_collection()
        elif not self.profile_db_path and not self.metagenome_mode:
            scg_taxonomy_super_dict['taxonomy'] = self.estimate_for_contigs_db_for_genome()
        elif self.metagenome_mode:
            scg_taxonomy_super_dict['taxonomy'] = self.estimate_for_contigs_db_for_metagenome()
        else:
            raise ConfigError("This class doesn't know how to deal with that yet :/")

        if self.compute_scg_coverages and self.metagenome_mode:
            scg_taxonomy_super_dict['coverages'] = self.get_scg_coverages_across_samples_dict_in_metagenome_mode(scg_taxonomy_super_dict)
        elif self.compute_scg_coverages and not self.metagenome_mode:
            scg_taxonomy_super_dict['coverages'] = self.get_scg_coverages_across_samples_dict_in_genome_mode(scg_taxonomy_super_dict)
        else:
            scg_taxonomy_super_dict['coverages'] = None

        return scg_taxonomy_super_dict


    def estimate(self):
        scg_taxonomy_super_dict = self.get_scg_taxonomy_super_dict()

        if self.update_profile_db_with_taxonomy:
            self.add_taxonomy_as_additional_layer_data(scg_taxonomy_super_dict)

        self.print_scg_taxonomy_super_dict(scg_taxonomy_super_dict)

        if self.output_file_path:
            self.store_scg_taxonomy_super_dict(scg_taxonomy_super_dict)


    def print_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        self.progress.reset()

        if self.collection_name:
            self.run.warning(None, header='Estimated taxonomy for collection "%s"' % self.collection_name, lc="green")
        elif self.metagenome_mode:
            self.run.warning(None, header='Taxa in metagenome "%s"' % self.contigs_db_project_name, lc="green")
        else:
            self.run.warning(None, header='Estimated taxonomy for "%s"' % self.contigs_db_project_name, lc="green")

        d = self.get_print_friendly_scg_taxonomy_super_dict(scg_taxonomy_super_dict)

        ordered_bin_names = sorted(list(d.keys()))

        if self.metagenome_mode:
            header = ['percent_identity', 'taxonomy']
        else:
            header = ['', 'total_scgs', 'supporting_scgs', 'taxonomy']

        # if we are in `--compute-scg-coverages` mode, and more than 5 sample names, we are in trouble since they will\
        # unlikely fit into the display while printing them. so here we will cut it to make sure things look OK.
        samples_not_shown = 0
        sample_names_to_display = None
        if self.compute_scg_coverages:
            sample_names_to_display = sorted(self.sample_names_in_profile_db)[0:5]
            samples_not_shown = sorted(self.sample_names_in_profile_db)[5:]

            header += sample_names_to_display

            if samples_not_shown:
                header += ['... %d more' % len(samples_not_shown)]

            # since we know coverages and sample names, we have a chance here to order the output
            # based on coverage. so let's do that.
            if self.metagenome_mode:
                sorted_bin_coverage_tuples = sorted([(bin_name, sum([d[bin_name]['coverages'][sample_name] for sample_name in self.sample_names_in_profile_db])) for bin_name in d], key=lambda x: x[1], reverse=True)
            else:
                sorted_bin_coverage_tuples = sorted([(bin_name, sum([(d[bin_name]['coverages'][sample_name] if d[bin_name]['supporting_scgs'] else 0) for sample_name in self.sample_names_in_profile_db])) for bin_name in d], key=lambda x: x[1], reverse=True)
            ordered_bin_names = [tpl[0] for tpl in sorted_bin_coverage_tuples]


        table = []
        for bin_name in ordered_bin_names:
            bin_data = d[bin_name]

            # set the taxonomy text depending on how much room we have. if there are sample coverages, keep it simple,
            # otherwise show the entire taxonomy text.
            if self.compute_scg_coverages:
                taxon_text_l = ['(%s) %s' % (l.split('_')[1][0], bin_data[l]) for l in self.ctx.levels_of_taxonomy[::-1] if bin_data[l]]
                taxon_text = taxon_text_l[0] if taxon_text_l else '(NA) NA'
            else:
                taxon_text = ' / '.join([bin_data[l] if bin_data[l] else '' for l in self.ctx.levels_of_taxonomy])

            # setting up the table columns here.
            if self.metagenome_mode:
                row = [bin_name, str(bin_data['percent_identity']), taxon_text]
            else:
                row = [bin_name, str(bin_data['total_scgs']), str(bin_data['supporting_scgs']), taxon_text]

            # if there are coverages, add samples to the display too
            if self.compute_scg_coverages:
                row += [d[bin_name]['coverages'][sample_name] for sample_name in sample_names_to_display]

            if samples_not_shown:
                row += ['... %d more' % len(samples_not_shown)]

            table.append(row)

        # if we are not in metagenome mode let's sort the output table based on total and
        # supporting SCGs
        if not self.metagenome_mode:
            table = sorted(table, key=lambda x: (int(x[1]), int(x[2])), reverse=True)

        anvio.TABULATE(table, header)


    def store_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        d = self.get_print_friendly_scg_taxonomy_super_dict(scg_taxonomy_super_dict)

        if self.metagenome_mode:
            headers = ['scg_name', 'percent_identity']
        else:
            headers = ['bin_name', 'total_scgs', 'supporting_scgs']

        headers += self.ctx.levels_of_taxonomy

        if self.compute_scg_coverages:
            headers_for_samples = sorted(self.sample_names_in_profile_db)
        else:
            headers_for_samples = []

        with open(self.output_file_path, 'w') as output:
            output.write('\t'.join(headers + headers_for_samples) + '\n')
            for item in d:
                line = [item] + [d[item][h] for h in headers[1:]]

                if self.compute_scg_coverages:
                    for sample_name in headers_for_samples:
                        line.append(d[item]['coverages'][sample_name])

                output.write('\t'.join([str(f) for f in line]) + '\n')

        self.run.info("Output file", self.output_file_path, nl_before=1)


    def get_print_friendly_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        d = {}

        if self.metagenome_mode:
            for scg_hit in scg_taxonomy_super_dict['taxonomy'][self.contigs_db_project_name]['scgs'].values():
                scg_hit_name = '%s_%d' % (scg_hit['gene_name'], scg_hit['gene_callers_id'])
                d[scg_hit_name] = scg_hit

                if self.compute_scg_coverages:
                    d[scg_hit_name]['coverages'] = scg_taxonomy_super_dict['coverages'][scg_hit['gene_callers_id']]
        else:
            for bin_name in scg_taxonomy_super_dict['taxonomy']:
                d[bin_name] = scg_taxonomy_super_dict['taxonomy'][bin_name]['consensus_taxonomy']
                d[bin_name]['total_scgs'] = scg_taxonomy_super_dict['taxonomy'][bin_name]['total_scgs']
                d[bin_name]['supporting_scgs'] = scg_taxonomy_super_dict['taxonomy'][bin_name]['supporting_scgs']

                if self.compute_scg_coverages:
                    d[bin_name]['coverages'] = scg_taxonomy_super_dict['coverages'][bin_name]

        return d


    def add_taxonomy_as_additional_layer_data(self, scg_taxonomy_super_dict):
        """A function that adds taxonomy to additional data tables of a given profile
           database. This will only work in metagenome mode."""

        if not self.metagenome_mode or not self.compute_scg_coverages:
            return

        self.progress.new("Adding summary taxonomy for samples")
        self.progress.update('...')

        scgs_dict = list(scg_taxonomy_super_dict['taxonomy'].values())[0]['scgs']

        # at this stage each scgs_dict entry will look like this, and most critically will
        # have the same SCG:
        #
        # "7660": {
        #       "gene_callers_id": 7660,
        #       "gene_name": "Ribosomal_S6",
        #       "accession": "CONSENSUS",
        #       "percent_identity": "98.9",
        #       "t_domain": "Bacteria",
        #       "t_phylum": "Firmicutes",
        #       "t_class": "Bacilli",
        #       "t_order": "Staphylococcales",
        #       "t_family": "Staphylococcaceae",
        #       "t_genus": "Staphylococcus",
        #       "t_species": null,
        #       "tax_hash": "b310c392"
        #     },
        #
        # this will enable us to learn the which SCG has been used to calculate coverage
        # information by only looking at a single entry:
        scg_name = list(scgs_dict.values())[0]['gene_name']

        # the might for loop to go through all taxonomic levels one by one
        for level in self.ctx.levels_of_taxonomy[::-1]:
            # setting the data group early on:
            data_group = '%s_%s' % (scg_name, level[2:])
            self.progress.update('Working on %s-level data' % level)
            data_dict = {}
            data_keys_list = set([])
            for sample_name in self.sample_names_in_profile_db:
                data_dict[sample_name] = Counter()
                for gene_callers_id in scgs_dict:
                    # starting with a tiny hack to fill in missing values. here we first find
                    # the most highly resolved level of taxonomy that is not null for this
                    # particular scg taxonomy
                    i = 0
                    for i in range(self.ctx.levels_of_taxonomy.index(level), 0, -1):
                        if scgs_dict[gene_callers_id][self.ctx.levels_of_taxonomy[i]]:
                            break

                    # just some abbreviations
                    l = self.ctx.levels_of_taxonomy[i][2:]
                    m = scgs_dict[gene_callers_id][self.ctx.levels_of_taxonomy[i]]

                    # if the best level we found in the previous step is matching to the level
                    # set by the main for loop, we're good to go with that name:
                    if level == self.ctx.levels_of_taxonomy[i]:
                        taxon_name = m
                    # otherwise we will try to replace that None name with something that is more
                    # sensible:
                    else:
                        taxon_name = "Unknown_%s_%s_%d" % (l, m, gene_callers_id)

                    # a key that will turn these data into stacked bar charts in the interface once
                    # they are added to the database:
                    key = '%s!%s' % (data_group, taxon_name)

                    # step where we add up all the values for each identical taxon names as we build
                    # the data dictionary:
                    data_dict[sample_name][key] += scg_taxonomy_super_dict['coverages'][gene_callers_id][sample_name]
                    data_keys_list.add(key)

            # next few lines demonstrate the power of anvi'o quite nicely:
            self.progress.update("Updating additional data tables...")
            args = argparse.Namespace(profile_db=self.profile_db_path, target_data_group=data_group, just_do_it=True)
            T = TableForLayerAdditionalData(args, r=run_quiet, p=progress_quiet)
            T.add(data_dict, list(data_keys_list))

            self.progress.reset()

            self.run.info_single("%s level taxonomy is added to the profile database." % (level.capitalize()))

        self.progress.end()


    def get_gene_caller_ids_for_splits(self, split_names_list):
        """Returns gene caller ids found in a list of splits"""

        gene_caller_ids_for_splits = set([])
        for split_name in split_names_list:
            if split_name in self.split_name_to_gene_caller_ids_dict:
                gene_caller_ids_for_splits.update(self.split_name_to_gene_caller_ids_dict[split_name])

        return gene_caller_ids_for_splits


    def get_split_names_for_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        """Returns a list of split names associated with SCGs found in a scg_taxonomy_super_dict."""

        if 'scgs' not in list(scg_taxonomy_super_dict['taxonomy'].values())[0]:
            raise ConfigError("Someone called this function with something that doesn't look like the kind "
                              "of input data it was expecting (sorry for the vagueness of the message, but "
                              "anvi'o hopes that will be able to find out why it is happening).")

        split_names = set([])

        for entry_name in scg_taxonomy_super_dict['taxonomy']:
            for gene_callers_id in scg_taxonomy_super_dict['taxonomy'][entry_name]['scgs']:
                split_names.add(self.gene_callers_id_to_split_name_dict[gene_callers_id])

        return split_names


    def get_scg_coverages_across_samples_dict_in_genome_mode(self, scg_taxonomy_super_dict):
        self.progress.reset()
        self.run.info_single("Anvi'o will now attempt to recover SCG coverages in GENOME MODE from the profile "
                             "database, which contains %d samples." % (len(self.sample_names_in_profile_db)), nl_before=1, nl_after=1)

        scg_coverages_across_samples_dict = self.get_scg_coverages_across_samples_dict(scg_taxonomy_super_dict)

        bin_avg_coverages_across_samples_dict = {}
        for bin_name in scg_taxonomy_super_dict['taxonomy']:
            bin_avg_coverages_across_samples_dict[bin_name] = dict([(sample_name, None) for sample_name in self.sample_names_in_profile_db])
            for sample_name in self.sample_names_in_profile_db:
                average_coverage_across_samples = [scg_coverages_across_samples_dict[gene_callers_id][sample_name] for gene_callers_id in scg_taxonomy_super_dict['taxonomy'][bin_name]['scgs']]
                if average_coverage_across_samples:
                    bin_avg_coverages_across_samples_dict[bin_name][sample_name] = np.mean(average_coverage_across_samples)

        self.run.warning("Anvi'o has just finished recovering SCG coverages from the profile database to estimate "
                         "the average coverage of your bins across your samples. Please note that anvi'o SCG taxonomy "
                         "framework is using only %d SCGs to estimate taxonomy. Which means, even a highly complete bin "
                         "may be missing all of them. In which case, the coverage of that bin will be `None` across all "
                         "your samples. The best way to prevent any misleading insights is take these results with a "
                         "huge grain of salt, and use the `anvi-summarize` output for critical applications." % len(self.ctx.anticodons),
                         header="FRIENDLY REMINDER", lc="blue")

        return bin_avg_coverages_across_samples_dict


    def get_scg_coverages_across_samples_dict_in_metagenome_mode(self, scg_taxonomy_super_dict):
        """Get SCG coverages in metagenome mode."""

        if not self.metagenome_mode:
            raise ConfigError("You're calling the wrong function. Your class is not in metagenome mode.")

        self.progress.reset()
        self.run.info_single("Anvi'o will now attempt to recover SCG coverages from the profile database, which "
                             "contains %d samples." % (len(self.sample_names_in_profile_db)), nl_before=1, nl_after=1)

        return self.get_scg_coverages_across_samples_dict(scg_taxonomy_super_dict)


    def get_scg_coverages_across_samples_dict(self, scg_taxonomy_super_dict):
        """Get SCG coverages"""
        scg_coverages_across_samples_dict = {}

        self.progress.new('Recovering coverages')
        self.progress.update('Learning all split names affiliated with SCGs ..')
        split_names_of_interest = self.get_split_names_for_scg_taxonomy_super_dict(scg_taxonomy_super_dict)
        self.progress.end()

        # initialize split coverages for splits that have anything to do with our SCGs
        args = copy.deepcopy(self.args)
        args.split_names_of_interest = split_names_of_interest
        args.collection_name = None
        profile_db = ProfileSuperclass(args, p=self.progress, r=run_quiet)
        profile_db.init_split_coverage_values_per_nt_dict()

        # recover all gene caller ids that occur in our taxonomy estimation dictionary
        # and ge their coverage stats from the profile super
        gene_caller_ids_of_interest = set([])
        for bin_name in scg_taxonomy_super_dict['taxonomy']:
            for gene_callers_id in scg_taxonomy_super_dict['taxonomy'][bin_name]['scgs']:
                gene_caller_ids_of_interest.add(gene_callers_id)

        # at this point we have everything. splits of interest are loaded in memory in `profile_db`, and we know
        # which gene caller ids we are interested in recovering coverages for. the way to access to gene coverages
        # is a bit convoluted in the dbops for historical reasons, but it is quite straightforward. the most
        # weird part is that we need a copy of a contigs super. so we will start with that:
        self.progress.new("Recovering SCG coverages")
        self.progress.update("Initiating the contigs super class")
        contigs_db = ContigsSuperclass(self.args, r=run_quiet, p=progress_quiet)

        for split_name in split_names_of_interest:
            self.progress.update("Working with %s" % split_name)
            # note for the curious: yes, here we are sending the same gene caller ids of interest over and over to
            # the `get_gene_level_coverage_stats` for each split, but that function is smart enough to not spend any
            # time on those gene caller ids that do not occur in the split name we are interested in.
            all_scg_stats_in_split = profile_db.get_gene_level_coverage_stats(split_name, contigs_db, gene_caller_ids_of_interest=gene_caller_ids_of_interest)

            for scg_stats in all_scg_stats_in_split.values():
                for entry in scg_stats.values():
                    gene_callers_id = int(entry['gene_callers_id'])
                    sample_name = entry['sample_name']
                    coverage = entry['non_outlier_mean_coverage']

                    if gene_callers_id not in scg_coverages_across_samples_dict:
                        scg_coverages_across_samples_dict[gene_callers_id] = dict([(sample_name, 0) for sample_name in self.sample_names_in_profile_db])

                    scg_coverages_across_samples_dict[gene_callers_id][sample_name] = coverage

        self.progress.end()

        return scg_coverages_across_samples_dict


class SetupLocalTRNATaxonomyData(TRNATaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        TRNATaxonomyArgs.__init__(self, self.args)

        # user accessible variables
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.redo_databases = A("redo_databases")
        self.num_threads = A('num_threads')

        self.ctx = ctx

        SanityCheck.__init__(self)


    def setup(self):
        """This function downloads all GTDB files necessary to setup the SCG databases anvi'o will rely upon.

           In addition to downloading the original files, the setup will make sure everything, including the
           DIAMOND search databases are in place.
        """

        self.run.warning("Please remember that the data anvi'o uses for SCG taxonomy is a "
                         "courtesy of The Genome Taxonomy Database (GTDB), an initiative to establish a "
                         "standardised microbial taxonomy based on genome phylogeny, primarly funded by "
                         "tax payers in Australia. Please don't forget to cite the original work, "
                         "doi:10.1038/nbt.4229 by Parks et al to explicitly mention the source of databases "
                         "anvi'o relies upon to estimate genome level taxonomy. If you are not sure how it "
                         "should look like in your methods sections, anvi'o developers will be happy to "
                         "help you if you can't find any published example to get inspiration.", lc = 'yellow')

        # let's try to figure out what do we want to do here. download files or just create/re-create
        # databases.
        an_anticodon_fasta, an_anticodon_database = list(self.ctx.anticodons.values())[0]['db'] + '.gz', list(self.ctx.anticodons.values())[0]['db'] + '.nhr'

        if os.path.exists(an_anticodon_fasta) and os.path.exists(an_anticodon_database) and self.redo_databases:
            self.run.warning("Anvi'o is removing all the previous databases so it can regenerate them from their "
                             "ashes.")

            db_paths = [v['db'] for v in self.ctx.anticodons.values()]
            for db_path in db_paths:
                os.remove(db_path) if os.path.exists(db_path) else None

        elif os.path.exists(an_anticodon_fasta) and os.path.exists(an_anticodon_database) and not self.redo_databases:
            raise ConfigError("It seems you have both the FASTA files and the search databases for anvi'o SCG taxonomy "
                              "in place. If you want to regenerate the databases for some reason, please use the flag "
                              "`--redo-databases`")

        elif os.path.exists(an_anticodon_fasta) and not os.path.exists(an_anticodon_database):
            self.run.warning("Anvi'o found your FASTA files in place, but not the databases. Now it will generate all "
                             "the search databases using the existing FASTA files.")

        self.create_search_databases()

        if not anvio.DEBUG:
            self.clean_up()


    def download_and_format_files(self):
        """There is no downloading files for this module.

        Essentially this is a manual operation. Please see the README.md file.
        """

        pass


    def create_search_databases(self):
        """Creates all the search databases"""

        self.progress.new("Creating search databases")
        self.progress.update("Removing any database that still exists in the output directory...")
        for prefix in ['.nhr', '.nin', '.nsq']:
            [os.remove(database_path) for database_path in [s['db'] + prefix for s in self.ctx.anticodons.values()] if os.path.exists(database_path)]

        # compresssing and decompressing FASTA files changes their hash and make them look like
        # modified in git. to avoid that, we will do the database generation in a temporary directory.
        temp_dir = filesnpaths.get_temp_directory_path()

        self.progress.update("Copying FASTA files to %s ..." % (temp_dir))
        # the following line basically returns a dictionary that shows the new path
        # of the FASTA file under temp_dir for a given anticodon .. apologies for the
        # incomprehensible list comprehension
        new_paths = dict([(os.path.basename(fasta_path), shutil.copy((fasta_path + '.gz'), os.path.join(temp_dir, os.path.basename(fasta_path) + '.gz'))) for fasta_path in [s['db'] for s in self.ctx.anticodons.values()]])

        missing_FASTA_files = [anticodon for anticodon in self.ctx.anticodons if not os.path.exists(new_paths[anticodon])]
        if len(missing_FASTA_files):
            raise ConfigError("Weird news :( Anvi'o is missing some FASTA files that were supposed to be somewhere. Since this "
                              "can't be your fault, it is not easy to advice what could be the solution to this. If you are not "
                              "an anvi'o programmer working on this problem this very moment, please get in touch with one.")

        self.progress.update("Decompressing FASTA files in %s" % (temp_dir))
        new_paths = dict([(anticodon, utils.gzip_decompress_file(new_paths[anticodon], keep_original=False)) for anticodon in new_paths])

        for anticodon in self.ctx.anticodons:
            self.progress.update("Working on %s in %d threads" % (anticodon, self.num_threads))

            FASTA_file_path_for_anticodon = new_paths[anticodon]

            # create a BLAST search database for `FASTA_file_path_for_anticodon`
            blast = BLAST(query_fasta=FASTA_file_path_for_anticodon, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
            blast.log_file_path = os.path.join(os.path.dirname(FASTA_file_path_for_anticodon), '%s.log' % anticodon)
            blast.makedb(dbtype='nucl')

            for prefix in ['.nhr', '.nin', '.nsq']:
                if not os.path.exists(FASTA_file_path_for_anticodon + prefix):
                    raise ConfigError("Something went wrong and BLAST did not create the database file it was supposed to "
                                      "for %s :(" % anticodon)
                else:
                    shutil.move(FASTA_file_path_for_anticodon + prefix, os.path.dirname(self.ctx.anticodons[anticodon]['db']))

        shutil.rmtree(temp_dir)

        self.progress.end()
        self.run.info_single("Every FASTA is now turned into a fancy search database. It means you are now allowed to run "
                             "`anvi-run-trna-taxonomy` on anvi'o contigs databases. This workflow is very new, and there are "
                             "caveats to it just like every other computational approach you use to make sense of complex 'omics "
                             "data. To better understand those caveats you should read our online documentation a bit. If you see "
                             "things that concerns you, please let anvi'o developers know. They love bad news. If you get good "
                             "results from this workflow, thank to those who contributed to the GTDB.", nl_after=1, mc="green")


    def clean_up(self):
        pass


class PopulateContigsDatabaseWithTRNATaxonomy(TRNATaxonomyArgs, SanityCheck, PopulateContigsDatabaseWithTaxonomy):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        TRNATaxonomyArgs.__init__(self, self.args)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.write_buffer_size = int(A('write_buffer_size') if A('write_buffer_size') is not None else 1000)
        self.contigs_db_path = A('contigs_db')
        self.num_parallel_processes = int(A('num_parallel_processes')) if A('num_parallel_processes') else 1
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        self.ctx = ctx

        self.max_target_seqs = 20
        self.evalue = float(A('e_value')) if A('e_value') else 1e-05
        self.min_pct_id = float(A('min_percent_identity')) if A('min_percent_identity') else 90

        SanityCheck.__init__(self)

        self.focus = 'trnas'
        PopulateContigsDatabaseWithTaxonomy.__init__(self, self.args)
