#!/usr/bin/env python3
# -*- coding: utf-8
"""
Classes to setup remote SCG databases in local, use local databases to affiliate SCGs in anvi'o
contigs databases with taxon names, and estimate taxonomy for genomes and metagneomes.
"""

import os
import glob
import copy
import shutil
import pandas as pd
import scipy.sparse as sps

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.dbops import ContigsDatabase
from anvio.drivers.diamond import Diamond
from anvio.genomedescriptions import MetagenomeDescriptions

from anvio.taxonomyops import AccessionIdToTaxonomy
from anvio.taxonomyops import TaxonomyEstimatorSingle
from anvio.taxonomyops import PopulateContigsDatabaseWithTaxonomy

run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


# if you need to change this, you're in trouble :) not really, but yes, you are..
locally_known_SCG_names = ['Ribosomal_S2',
                           'Ribosomal_S3_C',
                           'Ribosomal_S6',
                           'Ribosomal_S7',
                           'Ribosomal_S8',
                           'Ribosomal_S9',
                           'Ribosomal_S11',
                           'Ribosomal_S20p',
                           'Ribosomal_L1',
                           'Ribosomal_L2',
                           'Ribosomal_L3',
                           'Ribosomal_L4',
                           'Ribosomal_L6',
                           'Ribosomal_L9_C',
                           'Ribosomal_L13',
                           'Ribosomal_L16',
                           'Ribosomal_L17',
                           'Ribosomal_L20',
                           'Ribosomal_L21p',
                           'Ribosomal_L22',
                           'ribosomal_L24',
                           'Ribosomal_L27A']


class SCGTaxonomyContext(AccessionIdToTaxonomy):
    """The purpose of this base class is ot define file paths and constants for all single-copy
       core gene taxonomy operations.
    """
    def __init__(self, scgs_taxonomy_data_dir=None, database_release=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress
        self.focus = "scgs"

        # hard-coded GTDB variables. poor design, but I don't think we are going do need an
        # alternative to GTDB.
        self.target_database_name = "GTDB"
        self.target_database_release = database_release or 'v95.0'

        self.target_database_releases_yaml = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/GTDB-RELEASES.yaml')
        self.target_database_releases = utils.get_yaml_as_dict(self.target_database_releases_yaml)

        if self.target_database_release not in self.target_database_releases:
            raise ConfigError(f"Anvi'o doesn't know a release called {self.target_database_release}. It knows about "
                              f"{', '.join(self.target_database_releases)}, though.")

        self.target_database = self.target_database_releases[self.target_database_release]

        # some variables from anvi'o constants
        self.hmm_source_for_scg_taxonomy = constants.default_hmm_source_for_scg_taxonomy
        self.default_scgs_taxonomy_data_dir = constants.default_scgs_taxonomy_data_dir
        self.default_scgs_for_taxonomy = constants.default_scgs_for_taxonomy
        self.levels_of_taxonomy = constants.levels_of_taxonomy

        # these are all the user accessible paths. defaults will serve well for all applications,
        # but these can be used for debugging.
        self.SCGs_taxonomy_data_dir = (os.path.abspath(scgs_taxonomy_data_dir) if scgs_taxonomy_data_dir else None) or (os.path.join(self.default_scgs_taxonomy_data_dir, self.target_database_name))
        self.msa_individual_genes_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'MSA_OF_INDIVIDUAL_SCGs')
        self.accession_to_taxonomy_file_path = os.path.join(self.SCGs_taxonomy_data_dir, 'ACCESSION_TO_TAXONOMY.txt.gz')
        self.database_version_file_path = os.path.join(self.SCGs_taxonomy_data_dir, 'VERSION')
        self.search_databases_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'SCG_SEARCH_DATABASES')

        # some dictionaries for convenience. we set them up here, but the proper place to sanity check
        # them may be somewhere else. for instance, when this class is inheritded by SetupLocalSCGTaxonomyData
        # the paths will not point to an actual file, but when it is inherited by PopulateContigsDatabaseWithSCGTaxonomy,
        # they better point to actual files.
        self.SCGs = dict([(SCG, {'db': os.path.join(self.search_databases_dir_path, SCG + '.dmnd'), 'fasta': os.path.join(self.search_databases_dir_path, SCG)}) for SCG in self.default_scgs_for_taxonomy])

        self.accession_to_taxonomy_dict = {}

        # set version for ctx, so we know what version of the databases are on disk
        if os.path.exists(self.database_version_file_path):
            self.scg_taxonomy_database_version = open(self.database_version_file_path).readline().strip()
        else:
            self.scg_taxonomy_database_version = None

        # populate `self.accession_to_taxonomy_dict`
        AccessionIdToTaxonomy.__init__(self)


# here we create an instance for the module. the idea is to overwrite it if
# it is necessary to overwrite some of the defaults
ctx = SCGTaxonomyContext()


class SanityCheck(object):
    def __init__(self):
        if self.skip_sanity_check:
            self.run.warning("We are skipping all sanity checks :( Dangerous stuff is happening.")
        else:
            self.sanity_check()


    def sanity_check(self):
        if sorted(list(locally_known_SCG_names)) != sorted(self.ctx.default_scgs_for_taxonomy):
            raise ConfigError("Oh no. The SCGs designated to be used for all SCG taxonomy tasks in the constants.py "
                              "are not the same names described in locally known HMMs to remote FASTA files "
                              "conversion table definedd in SetupLocalSCGTaxonomyData module. If this makes zero "
                              "sense to you please ask a developer.")

        if not self.ctx.SCGs_taxonomy_data_dir:
            raise ConfigError("`SetupLocalSCGTaxonomyData` class is upset because it was inherited without "
                              "a directory for SCG taxonomy data to be stored :( This variable can't be None.")

        if self.user_taxonomic_level and self.user_taxonomic_level not in constants.levels_of_taxonomy:
            raise ConfigError("The taxonomic level %s is not a level anvi'o knows about. Here is the list of "
                              "taxonomic levels anvi'o recognizes: %s" % (', '.join(constants.levels_of_taxonomy)))

        # sanity checks specific to classes start below
        if self.__class__.__name__ in ['SetupLocalSCGTaxonomyData']:
            if self.reset and self.redo_databases:
                raise ConfigError("You can't ask anvi'o to both `--reset` and `--redo-databases` at the same time. Well. "
                                  "You can, but then this happens :/")

        if self.__class__.__name__ in ['SetupLocalSCGTaxonomyData', 'PopulateContigsDatabaseWithSCGTaxonomy']:
            if self.user_taxonomic_level:
                raise ConfigError("There is no need to set a taxonomic level while working with the class SetupLocalSCGTaxonomyData "
                                  "or PopulateContigsDatabaseWithSCGTaxonomy. Something fishy is going on :/")

        if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy', 'SCGTaxonomyEstimatorSingle', 'SCGTaxonomyEstimatorMulti']:
            if not os.path.exists(self.ctx.SCGs_taxonomy_data_dir):
                raise ConfigError("Anvi'o could not find the data directory for the single-copy core genes taxonomy "
                                  "setup. You may need to run `anvi-setup-scg-taxonomy`, or provide a directory path "
                                  "where SCG databases are set up. This is the current path anvi'o is considering (which "
                                  "can be changed via the `--scgs-taxonomy-data-dir` parameter): '%s'" % (self.ctx.SCGs_taxonomy_data_dir))

            if not os.path.exists(self.ctx.accession_to_taxonomy_file_path):
                raise ConfigError("While your SCG taxonomy data dir seems to be in place, it is missing at least one critical "
                                  "file (in this case, the file to resolve accession IDs to taxon names). You may need to run "
                                  "the program `anvi-setup-scg-taxonomy` with the `--reset` flag to set things right again.")

            filesnpaths.is_output_file_writable(self.all_hits_output_file_path, ok_if_exists=False) if self.all_hits_output_file_path else None

            filesnpaths.is_output_file_writable(self.per_scg_output_file) if self.per_scg_output_file else None


            ###########################################################
            # PopulateContigsDatabaseWithSCGTaxonomy
            ###########################################################
            if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy']:
                missing_SCG_databases = [SCG for SCG in self.ctx.SCGs if not os.path.exists(self.ctx.SCGs[SCG]['db'])]
                if len(missing_SCG_databases):
                    raise ConfigError("OK. It is very likley that if you run `anvi-setup-scg-taxonomy` first you will be golden. "
                                      "Because even though anvi'o found the directory for taxonomy headquarters, "
                                      "your setup seems to be missing %d of %d databases required for everything to work "
                                      "with the current genes configuration of this class (sources say this is a record, FYI)." % \
                                                (len(missing_SCG_databases), len(self.ctx.SCGs)))

            ###########################################################
            # SCGTaxonomyEstimatorSingle
            #
            # Note: if something down below complains about a paramter
            #       because that actually belongs to the multi estimator
            #       class, you may need to set it to null in the class
            #       SCGTaxonomyArgs for single estimator
            #       initiation if clause
            ###########################################################
            if self.__class__.__name__ in ['SCGTaxonomyEstimatorSingle']:
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

                scg_taxonomy_was_run = ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet).meta['scg_taxonomy_was_run']
                scg_taxonomy_database_version = ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet).meta['scg_taxonomy_database_version']
                if not scg_taxonomy_was_run:
                    raise ConfigError("It seems the SCG taxonomy tables were not populated in this contigs database :/ Luckily it "
                                      "is easy to fix that. Please see the program `anvi-run-scg-taxonomy`.")

                if scg_taxonomy_database_version != self.ctx.target_database_release:
                    self.progress.reset()
                    self.run.warning("The SCG taxonomy database on your computer has a different version (%s) than the SCG taxonomy information "
                                     "stored in your contigs database (%s). This is not a problem and things will most likely continue to work "
                                     "fine, but we wanted to let you know. You can get rid of this warning by re-running `anvi-run-scg-taxonomy` "
                                     "on your database." % (self.ctx.target_database_release, scg_taxonomy_database_version))

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

                if self.scg_name_for_metagenome_mode and self.scg_name_for_metagenome_mode not in self.ctx.SCGs:
                    raise ConfigError("We understand that you wish to work with '%s' to study the taxonomic make up of your contigs "
                                      "database in metagenome mode. But then this gene is not one of those anvi'o recognizes as "
                                      "suitable SCGs to do that. Here is a list for you to choose from: '%s'." \
                                                            % (self.scg_name_for_metagenome_mode, ', '.join(self.ctx.SCGs.keys())))

                if self.compute_scg_coverages and not self.profile_db_path:
                    raise ConfigError("The flag `--compute-scg-coverages` is only good if there is a non-blank profile database around "
                                      "from which anvi'o can learn coverage statistics of genes across one or more samples :/")

                if self.profile_db_path and self.metagenome_mode and not self.compute_scg_coverages:
                    raise ConfigError("You have a profile database and you have asked anvi'o to estimate taxonomy in metagenome mode, "
                                      "but you are not asking anvi'o to compute SCG coverages which doesn't make much sense :/ Removing "
                                      "the profile database from this command or adding the flag `--compute-scg-coverages` would have "
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
                                          "is NOT computing coverages values of SCGs across samples (pro tip: you can ask anvi'o to do "
                                          "it by adding the flag `--compute-scg-coverages` to your command line).")

            ###########################################################
            # SCGTaxonomyEstimatorMulti
            ###########################################################
            if self.__class__.__name__ in ['SCGTaxonomyEstimatorMulti']:
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


class SCGTaxonomyArgs(object):
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
        self.per_scg_output_file = A('per_scg_output_file')
        self.all_hits_output_file_path = A('all_hits_output_file')
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


class SCGTaxonomyEstimatorMulti(SCGTaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        """Iterate through metagenome descriptions using SCGTaxonomyEstimatorSingle"""

        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        SCGTaxonomyArgs.__init__(self, self.args)

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
        scg_taxonomy_database_versions_in_metagenomes = [g.metagenomes[m]['scg_taxonomy_database_version'] for m in g.metagenomes]
        if len(set(scg_taxonomy_database_versions_in_metagenomes)) > 1:
            self.progress.reset()
            self.run.warning("Please note that not all SCG taxonomy database versions across your metagenomes are identical. "
                             "This means the program `anvi-run-scg-taxonomy` was run on these database across different versions of "
                             "the source SCG taxonomy database. This is OK and things will continue to work, but you should consider "
                             "the fact that taxonomy estimations coming from different versions of the database may not be comparable "
                             "anymore depending on what has changed between different versions of the database. If your purpose is not "
                             "to compare different versions of the database, and if you would like to ensure consistency, you can re-run "
                             "`anvi-run-scg-taxonomy` on contigs databases that have a different version than what is installed on your "
                             "system, which is '%s' (if you run `anvi-db-info` on any contigs database you can learn the SCG database "
                             "version of it). Anvi'o found these versions across your metagenomes: '%s'." % \
                                        (self.ctx.target_database_release, ', '.join(list(set(scg_taxonomy_database_versions_in_metagenomes)))))
        elif scg_taxonomy_database_versions_in_metagenomes[0] != self.ctx.target_database_release:
            self.progress.reset()
            self.run.warning("While all of your metagenomes agree with each other and have the SCG taxonomy database version of %s, "
                              "this version differs from what is installed on your system, which is %s. If you don't do anything, "
                              "things will continue to work. But if you would like to get rid of this warning you will need to "
                              "re-run the program `anvi-run-scg-taxonomy` on each one of them 😬" % \
                                        (scg_taxonomy_database_versions_in_metagenomes[0], self.ctx.target_database_release))

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
            args = SCGTaxonomyArgs(self.args, format_args_for_single_estimator=True)
            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            if self.metagenome_mode:
                args.metagenome_mode = True
            else:
                args.metagenome_mode = False

            if self.profile_dbs_available:
                args.profile_db = self.metagenomes[metagenome_name]['profile_db_path']
                args.compute_scg_coverages = True

            d = SCGTaxonomyEstimatorSingle(args, run=run_quiet).get_print_friendly_items_taxonomy_super_dict(scg_taxonomy_super_dict_multi[metagenome_name])

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

        The difference between this function and `get_scg_taxonomy_super_dict` in `SCGTaxonomyEstimatorSingle`
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
            args = SCGTaxonomyArgs(self.args, format_args_for_single_estimator=True)

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
            scg_taxonomy_super_dict[metagenome_name] = SCGTaxonomyEstimatorSingle(args, progress=progress_quiet, run=run_quiet).get_items_taxonomy_super_dict()

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

            args = SCGTaxonomyArgs(self.args, format_args_for_single_estimator=True)
            args.compute_scg_coverages = False
            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            e = SCGTaxonomyEstimatorSingle(args, progress=progress_quiet, run=run_quiet)
            for scg_name in self.ctx.default_scgs_for_taxonomy:
                scg_frequencies[metagenome_name][scg_name] = e.frequency_of_items_with_taxonomy[scg_name]

            self.progress.increment()

        self.progress.update("Finalizing data")

        scg_frequencies_across_contigs_dbs = [(scg_name, sum([scg_frequencies[genome_name][scg_name] for genome_name in scg_frequencies])) for scg_name in self.ctx.default_scgs_for_taxonomy]
        scgs_ordered_based_on_frequency = [frequency_tuple[0] for frequency_tuple in sorted(scg_frequencies_across_contigs_dbs, key = lambda x: x[1], reverse=True)]

        num_scgs_for_each_contigs_db = [(genome_name, sum(scg_frequencies[genome_name].values())) for genome_name in scg_frequencies]
        contigs_dbs_ordered_based_on_num_scgs = [frequency_tuple[0] for frequency_tuple in sorted(num_scgs_for_each_contigs_db, key = lambda x: x[1], reverse=True)]

        self.progress.end()

        return scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies


class SCGTaxonomyEstimatorSingle(SCGTaxonomyArgs, SanityCheck, TaxonomyEstimatorSingle):
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

        SCGTaxonomyArgs.__init__(self, self.args)

        self.ctx = ctx

        SanityCheck.__init__(self)

        TaxonomyEstimatorSingle.__init__(self, skip_init=skip_init)


class SetupLocalSCGTaxonomyData(SCGTaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        SCGTaxonomyArgs.__init__(self, self.args)

        # user accessible variables
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.reset = A("reset") # complete start-over (downloading everything from GTDB)
        self.redo_databases = A("redo_databases") # just redo the databaes
        self.num_threads = A('num_threads')
        self.gtdb_release = A('gtdb_release')

        global ctx

        if self.gtdb_release:
            # re-initializing the context with the right release
            ctx = SCGTaxonomyContext(database_release=self.gtdb_release)

        self.ctx = ctx

        SanityCheck.__init__(self)


    def setup(self):
        """This function downloads all GTDB files necessary to setup the SCG databases anvi'o will rely upon.

           In addition to downloading the original files, the setup will make sure everything, including the
           DIAMOND search databases are in place.
        """

        if not anvio.DEBUG:
            self.check_initial_directory_structure()

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
        a_scg_fasta, a_scg_database = list(self.ctx.SCGs.values())[0]['fasta'] + '.gz', list(self.ctx.SCGs.values())[0]['db']

        if os.path.exists(a_scg_fasta) and os.path.exists(a_scg_database) and self.redo_databases:
            self.run.warning("Anvi'o is removing all the previous databases so it can regenerate them from their "
                             "ashes.")

            db_paths = [v['db'] for v in self.ctx.SCGs.values()]
            for db_path in db_paths:
                os.remove(db_path) if os.path.exists(db_path) else None

        elif os.path.exists(a_scg_fasta) and os.path.exists(a_scg_database) and not self.reset:
            raise ConfigError("It seems you have both the FASTA files and the search databases for anvi'o SCG taxonomy "
                              "in place. If you want to start from scratch and download everything from GTDB clean, use "
                              "the flag `--reset`.")

        elif self.reset:
            self.run.info("Local directory to setup", self.ctx.SCGs_taxonomy_data_dir)
            self.run.info("Reset the directory first", self.reset, mc="red")
            self.run.info("Remote database", self.ctx.target_database_name, nl_before=1, mc="green")
            self.run.info("Database version", self.ctx.target_database_release, mc="green")
            self.run.info("Base URL", self.ctx.target_database['base_url'])

            self.progress.new("%s setup" % self.ctx.target_database_name)
            self.progress.update("Reading the VERSION file...")
            content = utils.get_remote_file_content('/'.join([self.ctx.target_database['base_url'], self.ctx.target_database['files']['VERSION']]))
            version, release_date  = content.strip().split('\n')[0].strip(), content.strip().split('\n')[2].strip()
            self.progress.end()

            self.run.info("%s release found" % self.ctx.target_database_name, "%s (%s)" % (version, release_date), mc="green")

            self.download_and_format_files()

        elif os.path.exists(a_scg_fasta) and not os.path.exists(a_scg_database):
            self.run.warning("Anvi'o found your FASTA files in place, but not the databases. Now it will generate all "
                             "the search databases using the existing FASTA files.")

        self.create_search_databases()

        if not anvio.DEBUG:
            self.clean_up()


    def check_initial_directory_structure(self):
        if os.path.exists(self.ctx.SCGs_taxonomy_data_dir):
            if self.reset:
                shutil.rmtree(self.ctx.SCGs_taxonomy_data_dir)
                self.run.warning('The existing directory for SCG taxonomy data dir has been removed. Just so you know.')
                filesnpaths.gen_output_directory(self.ctx.SCGs_taxonomy_data_dir)
        else:
            filesnpaths.gen_output_directory(self.ctx.SCGs_taxonomy_data_dir)


    def download_and_format_files(self):
        temp_accession_to_taxonomy_file_path = '.gz'.join(self.ctx.accession_to_taxonomy_file_path.split('.gz')[:-1])

        # let's be 100% sure.
        os.remove(self.ctx.accession_to_taxonomy_file_path) if os.path.exists(self.ctx.accession_to_taxonomy_file_path) else None
        os.remove(temp_accession_to_taxonomy_file_path) if os.path.exists(temp_accession_to_taxonomy_file_path) else None

        for file_key in self.ctx.target_database['files']:
            remote_file_url = '/'.join([self.ctx.target_database['base_url'], self.ctx.target_database['files'][file_key]])
            local_file_path = os.path.join(self.ctx.SCGs_taxonomy_data_dir, file_key)

            utils.download_file(remote_file_url, local_file_path, progress=self.progress, run=self.run)

            if file_key in ['MSA_ARCHAEA.tar.gz', 'MSA_BACTERIA.tar.gz']:
                self.progress.new("Downloaded file patrol")
                self.progress.update("Unpacking file '%s'..." % os.path.basename(local_file_path))
                shutil.unpack_archive(local_file_path, extract_dir=self.ctx.msa_individual_genes_dir_path)
                os.remove(local_file_path)
                self.progress.end()

            if file_key in ['TAX_ARCHAEA.tsv', 'TAX_BACTERIA.tsv']:
                with open(temp_accession_to_taxonomy_file_path, 'a') as f:
                    f.write(open(local_file_path).read())
                    os.remove(local_file_path)

        # gzip ACCESSION_TO_TAXONOMY.txt, so it stays in its forever resting place
        # note: the following line will also remove the temporary file.
        utils.gzip_compress_file(temp_accession_to_taxonomy_file_path, self.ctx.accession_to_taxonomy_file_path, keep_original=False)

        # NEXT. learn paths for all FASTA files downloaded and unpacked
        fasta_file_paths = glob.glob(self.ctx.msa_individual_genes_dir_path + '/*.faa')

        if not fasta_file_paths:
            raise ConfigError("Something weird happened while anvi'o was trying to take care of the files "
                              "it downloaded from GTDB. Please let a developer know about this unless it is "
                              "not already reported in our issue tracker at Github :(")

        # files are done, but some of the FASTA files contain alignments solely composed of
        # gap characters :/ we will have to remove them to avoid fuck-ups in downstream
        # analyses
        self.progress.new("Clean up")
        for fasta_file_path in fasta_file_paths:
            self.progress.update("Looking for only-gap sequences from '%s'..." % os.path.basename(fasta_file_path))
            total_num_sequences, num_sequences_removed = utils.remove_sequences_with_only_gaps_from_fasta(fasta_file_path, fasta_file_path + '_CLEAN.fa', inplace=True)

            if num_sequences_removed:
                self.progress.reset()
                self.run.info_single('%d of %d seq in %s were all gaps and removed.' % (num_sequences_removed, total_num_sequences, os.path.basename(fasta_file_path)))

        # Start formatting things.
        self.progress.update("Checking output directory to store files ...")
        filesnpaths.is_output_dir_writable(os.path.dirname(self.ctx.search_databases_dir_path))
        filesnpaths.gen_output_directory(self.ctx.search_databases_dir_path, delete_if_exists=True, dont_warn=True)

        # We will be working with the files downloaded in whatever directory before. The first thing is to check
        # whether whether FASTA files in the directory are suitable for the conversion
        self.progress.update("Checking the conversion dict and FASTA files ...")
        msa_individual_gene_names_required = []
        for SCG in locally_known_SCG_names:
            msa_individual_gene_names_required.extend(self.ctx.target_database['genes'][SCG])

        fasta_file_paths = glob.glob(self.ctx.msa_individual_genes_dir_path + '/*.faa')
        msa_individual_gene_names_downloaded = [os.path.basename(f) for f in fasta_file_paths]

        missing_msa_gene_names = [n for n in msa_individual_gene_names_required if n not in msa_individual_gene_names_downloaded]
        if missing_msa_gene_names:
            self.progress.reset()
            raise ConfigError("Big trouble :( Anvi'o uses a hard-coded dictionary to convert locally known "
                              "HMM models to FASTA files reported by GTDB project. It seems something has changed "
                              "and %d of the FASTA files expected to be in the download directory are not there. "
                              "Here is that list: '%s'. Someone needs to update the codebase by reading the "
                              "appropriate documentation. If you are a user, you can't do much at this point but "
                              "contacting the developers :( Anvi'o will keep the directory that contains all the "
                              "downloaded files to update the conversion dictionary. Here is the full path to the "
                              "output: %s" % (len(missing_msa_gene_names), ', '.join(missing_msa_gene_names), self.ctx.msa_individual_genes_dir_path))
        else:
            self.progress.reset()
            self.run.info_single("Good news! The conversion dict and the FASTA files it requires seem to be in place. "
                                 "Anvi'o is now ready to to merge %d FASTA files that correspond to %d SCGs, and "
                                 "create individual search databases for them." % \
                                        (len(msa_individual_gene_names_required), len(locally_known_SCG_names)), nl_before=1, nl_after=1, mc="green")

        # Merge FASTA files that should be merged. This is defined in the conversion dictionary.
        for SCG in locally_known_SCG_names:
            self.progress.update("Working on %s ..." % (SCG))

            files_to_concatenate = [os.path.join(self.ctx.msa_individual_genes_dir_path, f) for f in self.ctx.target_database['genes'][SCG]]
            FASTA_file_for_SCG = os.path.join(self.ctx.search_databases_dir_path, SCG)

            # concatenate from the dictionary into the new destination
            utils.concatenate_files(FASTA_file_for_SCG, files_to_concatenate)

            # compress the FASTA file
            utils.gzip_compress_file(FASTA_file_for_SCG)

        self.progress.end()


    def create_search_databases(self):
        """Creates all the search databases"""

        self.progress.new("Creating search databases")
        self.progress.update("Removing any database that still exists in the output directory...")
        [os.remove(database_path) for database_path in [s['db'] for s in self.ctx.SCGs.values()] if os.path.exists(database_path)]

        # compresssing and decompressing FASTA files changes their hash and make them look like
        # modified in git. to avoid that, we will do the database generation in a temporary directory.
        temp_dir = filesnpaths.get_temp_directory_path()

        self.progress.update("Copying FASTA files to %s ..." % (temp_dir))
        # the following line basically returns a dictionary that shows the new path
        # of the FASTA file under temp_dir for a given SCG .. apologies for the
        # incomprehensible list comprehension
        new_paths = dict([(os.path.basename(fasta_path), shutil.copy((fasta_path + '.gz'), os.path.join(temp_dir, os.path.basename(fasta_path) + '.gz'))) for fasta_path in [s['fasta'] for s in self.ctx.SCGs.values()]])

        missing_FASTA_files = [SCG for SCG in self.ctx.SCGs if not os.path.exists(new_paths[SCG])]
        if len(missing_FASTA_files):
            raise ConfigError("Weird news :( Anvi'o is missing some FASTA files that were supposed to be somewhere. Since this "
                              "can't be your fault, it is not easy to advice what could be the solution to this. But you can "
                              "always try to re-run `anvi-setup-scg-taxonomy` with `--reset` flag.")

        self.progress.update("Decompressing FASTA files in %s" % (temp_dir))
        new_paths = dict([(SCG, utils.gzip_decompress_file(new_paths[SCG], keep_original=False)) for SCG in new_paths])

        # Merge FASTA files that should be merged. This is defined in the conversion dictionary.
        for SCG in self.ctx.SCGs:
            self.progress.update("Working on %s in %d threads" % (SCG, self.num_threads))

            FASTA_file_path_for_SCG = new_paths[SCG]

            # create a diamond search database for `FASTA_file_path_for_SCG`
            diamond = Diamond(query_fasta=FASTA_file_path_for_SCG, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
            diamond.makedb(output_file_path=self.ctx.SCGs[SCG]['db'])

            if not os.path.exists(self.ctx.SCGs[SCG]['db']):
                raise ConfigError("Something went wrong and DIAMOND did not create the database file it was supposed to "
                                  "for %s :(" % SCG)

        shutil.rmtree(temp_dir)

        self.progress.end()
        self.run.info_single("Every FASTA is now turned into a fancy search database. It means you are now allowed to run "
                             "`anvi-run-scg-taxonomy` on anvi'o contigs databases. This workflow is very new, and there are "
                             "caveats to it just like every other computational approach you use to make sense of complex 'omics "
                             "data. To better understand those caveats you should read our online documentation a bit. If you see "
                             "things that concerns you, please let anvi'o developers know. They love bad news. If you get good "
                             "results from this workflow, thank to those who contributed to the GTDB.", nl_after=1, mc="green")


    def clean_up(self):
        for file_path in [os.path.join(self.ctx.SCGs_taxonomy_data_dir, 'diamond-log-file.txt')]:
            if os.path.exists(file_path):
                os.remove(file_path)

        for dir_path in [self.ctx.msa_individual_genes_dir_path]:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)


class PopulateContigsDatabaseWithSCGTaxonomy(SCGTaxonomyArgs, SanityCheck, PopulateContigsDatabaseWithTaxonomy):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        SCGTaxonomyArgs.__init__(self, self.args)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.write_buffer_size = int(A('write_buffer_size') if A('write_buffer_size') is not None else 1000)
        self.contigs_db_path = A('contigs_db')
        self.num_parallel_processes = int(A('num_parallel_processes')) if A('num_parallel_processes') else 1
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        self.ctx = ctx

        self.max_target_seqs = int(A('max_num_target_sequences')) or 20
        self.evalue = float(A('e_value')) if A('e_value') else 1e-05
        self.min_pct_id = float(A('min_percent_identity')) if A('min_percent_identity') else 90

        SanityCheck.__init__(self)

        PopulateContigsDatabaseWithTaxonomy.__init__(self, self.args)
