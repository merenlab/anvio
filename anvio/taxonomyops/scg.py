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
from anvio.dbops import ContigsDatabase, ContigsSuperclass
from anvio.drivers.diamond import Diamond
from anvio.genomedescriptions import GenomeDescriptions, MetagenomeDescriptions

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
locally_known_SCG_names = ['Ribosomal_L1',
                           'Ribosomal_L13',
                           'Ribosomal_L14',
                           'Ribosomal_L16',
                           'Ribosomal_L17',
                           'Ribosomal_L19',
                           'Ribosomal_L2',
                           'Ribosomal_L20',
                           'Ribosomal_L21p',
                           'Ribosomal_L22',
                           'Ribosomal_L27A',
                           'Ribosomal_L3',
                           'Ribosomal_L4',
                           'Ribosomal_L5',
                           'Ribosomal_S11',
                           'Ribosomal_S15',
                           'Ribosomal_S16',
                           'Ribosomal_S2',
                           'Ribosomal_S6',
                           'Ribosomal_S7',
                           'Ribosomal_S8',
                           'Ribosomal_S9']


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

        # some variables from anvi'o constants
        self.hmm_source_for_scg_taxonomy = constants.default_hmm_source_for_scg_taxonomy
        self.default_scgs_taxonomy_data_dir = constants.default_scgs_taxonomy_data_dir
        self.default_scgs_for_taxonomy = constants.default_scgs_for_taxonomy
        self.levels_of_taxonomy = constants.levels_of_taxonomy

        # these are all the user accessible paths. defaults will serve well for all applications,
        # but these can be used for debugging.
        self.SCGs_taxonomy_data_dir = (os.path.abspath(scgs_taxonomy_data_dir) if scgs_taxonomy_data_dir else None) or (os.path.join(self.default_scgs_taxonomy_data_dir, self.target_database_name))
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
            raise ConfigError("Oh no. The default SCG names to be used for all SCG taxonomy tasks seem to differ "
                              "from those names for which you have names described in locally known HMMs to remote FASTA files "
                              "conversion table defined in SetupLocalSCGTaxonomyData module. If this makes zero "
                              "sense to you please ask a developer.")

        if not self.ctx.SCGs_taxonomy_data_dir:
            raise ConfigError("`SetupLocalSCGTaxonomyData` class is upset because it was inherited without "
                              "a directory for SCG taxonomy data to be stored :( This variable can't be None.")

        if self.user_taxonomic_level and self.user_taxonomic_level not in constants.levels_of_taxonomy:
            raise ConfigError("The taxonomic level %s is not a level anvi'o knows about. Here is the list of "
                              "taxonomic levels anvi'o recognizes: %s" % (', '.join(constants.levels_of_taxonomy)))

        # sanity checks specific to classes start below
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
                                  "file (in this case, the file to resolve accession IDs to taxon names). You may need to "
                                  "restore your anvi'o Git repository. Maybe reach out to us for help on Github or Discord.")

            filesnpaths.is_output_file_writable(self.all_hits_output_file_path, ok_if_exists=False) if self.all_hits_output_file_path else None

            filesnpaths.is_output_file_writable(self.per_scg_output_file) if self.per_scg_output_file else None


            ###########################################################
            # PopulateContigsDatabaseWithSCGTaxonomy
            ###########################################################
            if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy']:
                missing_SCG_databases = [SCG for SCG in self.ctx.SCGs if not os.path.exists(self.ctx.SCGs[SCG]['db'])]
                if len(missing_SCG_databases):
                    raise ConfigError("We have a problem, Houston. Even though anvi'o found the directory for taxonomy headquarters, "
                                      "your setup seems to be missing %d of %d databases required for everything to work "
                                      "properly :/ The good news? This problem will very likely go away if you run the program "
                                      "`anvi-setup-scg-taxonomy` and you will be golden." % (len(missing_SCG_databases), len(self.ctx.SCGs)))

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
                            "arguments also include an input file for multiple (meta)genomes. It is a no no. Please choose either. ")

                if self.output_file_prefix:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file path "
                                      "than an output file prefix.")

                if self.output_file_path:
                    filesnpaths.is_output_file_writable(self.output_file_path)

                if self.sequences_file_path_prefix:
                    filesnpaths.is_output_file_writable(self.sequences_file_path_prefix)

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

                if scg_taxonomy_database_version != self.ctx.scg_taxonomy_database_version:
                    self.progress.reset()
                    raise ConfigError("The SCG taxonomy database version on your computer (%s) is different than the SCG taxonomy database "
                                      "version to populate your contigs database (%s). Please re-run the program `anvi-run-scg-taxonomy` "
                                      "on your contigs-db." % (self.ctx.scg_taxonomy_database_version, scg_taxonomy_database_version))

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
                    raise ConfigError("Taxonomy estimation classes have been initiated with files for multiple (meta)genomes, but "
                                      "your arguments include also a single contigs or profile database path. You make anvi'o nervous. "
                                      "Please run this program either with a (meta)genomes file or individual contigs/profile databases.")
                
                if self.args.external_genomes and self.args.metagenomes:
                    raise ConfigError("More than one input file type (external genomes AND metagenomes) has been given to the "
                                      "taxonomy estimation classes. Please run this program with only one input type at a time.")

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
        self.sequences_file_path_prefix = A('report_scg_sequences_file_prefix')
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
            self.sequences_file_path_prefix = None
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


    def init_external_genomes(self):
        g = GenomeDescriptions(self.args, run=run_quiet, progress=self.progress)
        g.load_genomes_descriptions(skip_functions=True, init=False)

        # NOTE some enforced flags
        self.profile_dbs_available = False # we don't load profile dbs for external genomes, so this has to be off
        if self.metagenome_mode:
            self.metagenome_mode = False

        genomes_without_scg_taxonomy = [x for x in g.genomes if not g.genomes[x]['scg_taxonomy_was_run']]
        if genomes_without_scg_taxonomy:
            n_genomes = len(g.genomes)
            n_without_tax = len(genomes_without_scg_taxonomy)
            if n_without_tax == n_genomes:
                self.progress.end()
                raise ConfigError(f"Surprise! All of the {n_genomes} genomes had no SCG taxonomy information. You need "
                                  f"to run `anvi-run-scg-taxonomy` on all of them before you continue.")
            else:
                self.progress.end()
                g_str = ', '.join(genomes_without_scg_taxonomy)
                raise ConfigError(f"{n_without_tax} of your {n_genomes} genomes has no SCG taxonomy information. Here is the list of "
                                  f"genomes you need to run `anvi-run-scg-taxonomy` on: '{g_str}'.")

        # check if SCG versions agree with each other and with installed version
        scg_taxonomy_database_versions_in_genomes = [g.genomes[x]['scg_taxonomy_database_version'] for x in g.genomes]
        if len(set(scg_taxonomy_database_versions_in_genomes)) > 1:
            self.progress.reset()
            self.run.warning("Please note that not all SCG taxonomy database versions across your genomes are identical. "
                             "This means the program `anvi-run-scg-taxonomy` was run on these database across different versions of "
                             "the source SCG taxonomy database. This is OK and things will continue to work, but you should consider "
                             "the fact that taxonomy estimations coming from different versions of the database may not be comparable "
                             "anymore depending on what has changed between different versions of the database. If your purpose is not "
                             "to compare different versions of the database, and if you would like to ensure consistency, you can re-run "
                             "`anvi-run-scg-taxonomy` on contigs databases that have a different version than what is installed on your "
                             "system, which is '%s' (if you run `anvi-db-info` on any contigs database you can learn the SCG database "
                             "version of it). Anvi'o found these versions across your genomes: '%s'." % \
                                        (self.ctx.scg_taxonomy_database_version, ', '.join(list(set(scg_taxonomy_database_versions_in_genomes)))))
        elif scg_taxonomy_database_versions_in_genomes[0] != self.ctx.scg_taxonomy_database_version:
            self.progress.reset()
            self.run.warning("While all of your genomes agree with each other and have the SCG taxonomy database version of %s, "
                              "this version differs from what is installed on your system, which is %s. If you don't do anything, "
                              "things will continue to work. But if you would like to get rid of this warning you will need to "
                              "re-run the program `anvi-run-scg-taxonomy` on each one of them 😬" % \
                                        (scg_taxonomy_database_versions_in_genomes[0], self.ctx.scg_taxonomy_database_version))

        # we keep these attribute names the same as for metagenomes so that we don't have to duplicate every function later
        self.metagenomes = copy.deepcopy(g.genomes)
        self.metagenome_names = copy.deepcopy(g.external_genome_names)


    def init_metagenomes(self):
        self.progress.new("Initializing contigs DBs for metagenomes")
        self.progress.update("...")
        g = MetagenomeDescriptions(self.args, run=run_quiet, progress=self.progress)
        g.load_metagenome_descriptions(skip_sanity_check=True)

        # NOTE some enforced flags here.
        self.compute_scg_coverages = g.profile_dbs_available
        if not self.metagenome_mode and g.profile_dbs_available:
            self.metagenome_mode = True

        metagenomes_without_scg_taxonomy = [m for m in g.metagenomes if not g.metagenomes[m]['scg_taxonomy_was_run']]
        if metagenomes_without_scg_taxonomy:
            n_metagenomes = len(g.metagenomes)
            n_without_tax = len(metagenomes_without_scg_taxonomy)
            if n_without_tax == n_metagenomes:
                self.progress.end()
                raise ConfigError(f"Surprise! None of the {n_metagenomes} metagenomes had SCG taxonomy information.")
            else:
                self.progress.end()
                m_str = ', '.join(metagenomes_without_scg_taxonomy)
                raise ConfigError(f"{n_without_tax} of your {n_metagenomes} metagenomes has no SCG taxonomy information. "
                                  f"Here is the list of metagenomes you need to run `anvi-run-scg-taxonomy` on:: '{m_str}'.")

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
                                        (self.ctx.scg_taxonomy_database_version, ', '.join(list(set(scg_taxonomy_database_versions_in_metagenomes)))))
        elif scg_taxonomy_database_versions_in_metagenomes[0] != self.ctx.scg_taxonomy_database_version:
            self.progress.reset()
            self.run.warning("While all of your metagenomes agree with each other and have the SCG taxonomy database version of %s, "
                              "this version differs from what is installed on your system, which is %s. If you don't do anything, "
                              "things will continue to work. But if you would like to get rid of this warning you will need to "
                              "re-run the program `anvi-run-scg-taxonomy` on each one of them 😬" % \
                                        (scg_taxonomy_database_versions_in_metagenomes[0], self.ctx.scg_taxonomy_database_version))

        self.metagenomes = copy.deepcopy(g.metagenomes)
        self.metagenome_names = copy.deepcopy(g.metagenome_names)
        self.profile_dbs_available = g.profile_dbs_available

        self.progress.end()


    def estimate_for_genomes(self):
        if not self.metagenomes:
            self.init_external_genomes()
            self.run.info("Num genomes", len(self.metagenome_names))

        self.run.info("Taxonomic level of interest", self.user_taxonomic_level or "(None specified by the user, so 'all levels')")
        if self.output_file_path:
            self.run.info("Output file path", self.output_file_path)
        if self.output_file_prefix:
            self.run.info("Output file prefix", self.output_file_prefix)
            self.run.info("Output in matrix format", self.matrix_format)
        self.run.info("Output raw data", self.raw_output)
        self.run.info("SCG coverages will be computed?", self.compute_scg_coverages)

        if self.report_scg_frequencies_path:
            self.report_scg_frequencies_as_TAB_delimited_file()
            return

        scg_taxonomy_super_dict_multi = self.get_scg_taxonomy_super_dict_for_metagenomes()

        if self.sequences_file_path_prefix:
            self.store_sequences_for_items_multi(scg_taxonomy_super_dict_multi)

        if self.output_file_path:
            self.store_scg_taxonomy_super_dict_multi_output_file(scg_taxonomy_super_dict_multi)
        if self.output_file_prefix:
            self.store_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)


    def estimate_for_metagenomes(self):
        if not self.metagenomes:
            self.init_metagenomes()
            self.run.info("Num metagenomes", len(self.metagenome_names))

        self.run.info("Taxonomic level of interest", self.user_taxonomic_level or "(None specified by the user, so 'all levels')")
        if self.output_file_path:
            self.run.info("Output file path", self.output_file_path)
        if self.output_file_prefix:
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

        if self.sequences_file_path_prefix:
            self.store_sequences_for_items_multi(scg_taxonomy_super_dict_multi)

        if self.output_file_path:
            self.store_scg_taxonomy_super_dict_multi_output_file(scg_taxonomy_super_dict_multi)
        if self.output_file_prefix:
            self.store_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)


    def estimate(self):
        
        if self.args.metagenomes:
            self.estimate_for_metagenomes()
        elif self.args.external_genomes:
            self.estimate_for_genomes()
        else:
            raise ConfigError("Anvi'o is not sure how things got to this point, but somehow we find ourselves without input for the "
                              "SCGTaxonomyEstimatorMulti class's estimate() function.")


    def store_sequences_for_items_multi(self, scg_taxonomy_super_dict_multi):
        """Report sequences for items if possible"""

        if self.ctx.focus != 'scgs':
            raise ConfigError("This function is only tested in SCGs mode. If you need to report "
                              "sequences for taxonomy items reported in other foci, please get in "
                              "touch with anvi'o developers.")

        if not self.scg_name_for_metagenome_mode:
            raise ConfigError("You can't ask anvi'o to store sequences for SCGs unless you are "
                              "working with a specific SCG name :(")

        d = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)

        import argparse

        dna_sequences_output_file_path = self.sequences_file_path_prefix + '_DNA.fa'
        amino_acid_sequences_output_file_path = self.sequences_file_path_prefix + '_AA.fa'

        if self.just_do_it:
            pass
        elif os.path.exists(dna_sequences_output_file_path) or os.path.exists(amino_acid_sequences_output_file_path):
            raise ConfigError(f"Anvi'o has detected you already have files {dna_sequences_output_file_path} or {amino_acid_sequences_output_file_path} "
                               "and does not want to overwrite them. If you do want to overwrite them, then run the `--just-do-it` flag.")

        aa_sequences_output = open(amino_acid_sequences_output_file_path,'w')
        dna_sequences_output = open(dna_sequences_output_file_path,'w')

        for metagenome_name, v in d.items():
            contigs_db_path = self.metagenomes[metagenome_name]['contigs_db_path']

            args_for_contigsDB = argparse.Namespace()
            args_for_contigsDB.contigs_db = contigs_db_path
            c = ContigsSuperclass(args_for_contigsDB, r=run_quiet)

            gene_caller_ids = [values['gene_callers_id'] for values in v.values()]

            gene_caller_ids_list, sequences_dict = c.get_sequences_for_gene_callers_ids(gene_caller_ids, include_aa_sequences=True)

            if not len(gene_caller_ids_list):
                raise ConfigError("Something that should have never happened, happened :/ Please re-run the same command with "
                                  "`--debug` and send the Traceback to an anvi'o developer.")

            with open(amino_acid_sequences_output_file_path, 'a+') as aa_sequences_output, open(dna_sequences_output_file_path, 'a+') as dna_sequences_output:
                for header, entry in d[metagenome_name].items():
                    dna_sequence = sequences_dict[entry['gene_callers_id']]['sequence']
                    amino_acid_sequence = sequences_dict[entry['gene_callers_id']]['aa_sequence']

                    aa_sequences_output.write(f">{header}\n{amino_acid_sequence}\n")
                    dna_sequences_output.write(f">{header}\n{dna_sequence}\n")

        self.run.info("DNA sequences for SCGs", dna_sequences_output_file_path, nl_before=1)
        self.run.info("AA sequences for SCGs", amino_acid_sequences_output_file_path)


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
            args.name = metagenome_name

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
                        raise ConfigError("We found more than one coverage value per SCG, which means you gave anvi'o a merged profile "
                        "database (containing multiple samples) associated with your contigs database. The codebase is not ready to "
                        "handle this :(  You need to provide single profiles instead of merged ones.")

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

        header = ['identifier', 'metagenome_name', 'gene_name', 'gene_callers_id', 'percent_identity']

        if self.compute_scg_coverages:
            header += ['coverage']

        header += taxonomic_levels

        output_file_path = self.output_file_prefix + '-RAW-LONG-FORMAT.txt'
        with open(output_file_path, 'w') as output:
            output.write('\t'.join(header) + '\n')

            for metagenome_name in d:
                for gene_name in d[metagenome_name]:
                    output.write('\t'.join([gene_name] + [metagenome_name] + [str(d[metagenome_name][gene_name][h]) for h in header[2:]]) + '\n')

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

    
    def store_scg_taxonomy_super_dict_multi_output_file(self, scg_taxonomy_super_dict_multi):
        """Generates an output just like TaxonomyEstimatorSingle.store_items_taxonomy_super_dict(), but for multiple inputs."""
        d = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)

        headers = ['name', 'total_scgs', 'supporting_scgs']
        headers += self.ctx.levels_of_taxonomy

        with open(self.output_file_path, 'w') as output:
            output.write('\t'.join(headers) + '\n')
            for metagenome in d:
                for name in d[metagenome]:
                    # should be only one element in the innermost dictionary
                    line = [name] + [d[metagenome][name][h] for h in headers[1:]]

                    output.write('\t'.join([str(f) for f in line]) + '\n')

        self.run.info("Output file", self.output_file_path, nl_before=1)


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
                self.run.info_single("Your (meta)genome file DOES NOT contain profile databases, and you haven't asked anvi'o to "
                                     "work in `--metagenome-mode`. Your contigs databases will be treated as genomes rather than "
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
        self.name = A('name')

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
        self.num_threads = A('num_threads')
        self.SCGs_taxonomy_data_dir = A('scgs_taxonomy_data_dir')

        global ctx

        self.ctx = ctx

        SanityCheck.__init__(self)


    def setup(self):
        """The setup will make sure SCG FASTA files and DIAMOND search databases are in place."""

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


        self.run.warning("Anvi'o will now attempt to generate the SCG search databases from scratch. Fingers crossed.")

        self.create_search_databases()

        if not anvio.DEBUG:
            self.clean_up()


    def check_initial_directory_structure(self):
        if not os.path.exists(self.ctx.SCGs_taxonomy_data_dir):
            raise ConfigError(f"Bad news. Anvi'o cannot find the directory for SCG taxonomy data, which should be at "
                              f"{self.ctx.SCGs_taxonomy_data_dir}. Fixing this problem may require restoring your "
                              f"anvi'o Git repository (or perhaps cloning it from scratch). "
                              f"Since this shouldn't happen, please consider reaching out to us via anvi'o Discord so "
                              f"we can help you sort this out.")


    def create_search_databases(self):
        """Creates all the search databases"""

        # let's first check if FASTA files are in place
        missing_FASTA_files = [(s['fasta'] + '.gz') for s in self.ctx.SCGs.values() if not os.path.exists(s['fasta'] + '.gz')]
        if len(missing_FASTA_files):
            self.run.warning(None, header="MISSING FASTA FILES")
            for missing_FASTA_file in missing_FASTA_files:
                self.run.info_single(f"{missing_FASTA_file}", cut_after=None, mc='red')
            raise ConfigError("The source FASTA files (from which the setup builds databases) do not seem to be in place :/ "
                              "Anvi'o is as confused as you are regarding why this could be the case :( Please consider "
                              "reaching out to us via anvi'o Discord so we can help you sort this out.")

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
                              "can't be your fault, it is not easy to advise what could be the solution to this. You can "
                              "try to run `git status` within your copy of the anvi'o repository to see if some critical files "
                              "were deleted and can be restored. If not, please feel free to reach out to the developers for help.")

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
        self.SCGs_taxonomy_data_dir = A('scgs_taxonomy_data_dir')

        global ctx

        # if the user has specified a GTDB release or taxonomy data directory,
        # we're getting a new context.
        if self.SCGs_taxonomy_data_dir:
            # re-initializing the context with the right release
            ctx = SCGTaxonomyContext(scgs_taxonomy_data_dir=self.SCGs_taxonomy_data_dir)

        self.ctx = ctx

        self.max_target_seqs = int(A('max_num_target_sequences')) if A('max_num_target_sequences') else 20
        self.evalue = float(A('e_value')) if A('e_value') else 1e-05
        self.min_pct_id = float(A('min_percent_identity')) if A('min_percent_identity') else 90

        SanityCheck.__init__(self)

        PopulateContigsDatabaseWithTaxonomy.__init__(self, self.args)
