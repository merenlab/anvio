#!/usr/bin/env python
# -*- coding: utf-8
"""Main KEGG metabolism estimator classes."""

import os
import re
import copy
import json
import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.dbinfo import DBInfo
from anvio.errors import ConfigError
from anvio.genomedescriptions import MetagenomeDescriptions, GenomeDescriptions

from anvio.metabolism.modulesdb import ModulesDatabase
from anvio.metabolism.algorithms import KeggEstimationAlgorithms
from anvio.metabolism.dbaccess import KeggEstimatorArgs, KeggDataLoader
from anvio.metabolism.constants import OUTPUT_MODES, MODULE_METADATA_HEADERS, KO_METADATA_HEADERS, STEP_METADATA_HEADERS

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Iva Veseli"
__email__ = "iveseli@uchicago.edu"

run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print
P = terminal.pluralize


class KeggMetabolismEstimator(KeggEstimatorArgs, KeggDataLoader, KeggEstimationAlgorithms):
    """ Class for reconstructing/estimating metabolism for a SINGLE contigs DB based on hits to KEGG databases.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-estimate-metabolism
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.pan_db_path = A('pan_db')
        self.genomes_storage_path = A('genomes_storage')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file = A('bin_ids_file')
        self.contigs_db_project_name = "Unknown"
        self.database_name = A('database_name')
        self.multi_mode = True if A('multi_mode') else False

        # This can be initialized later if necessary using init_gene_coverage()
        self.profile_db = None
        # This can be initialized later if necessary by setup_ko_dict()
        self.ko_dict = {}
        self.stray_ko_dict = {}

        # INIT BASE CLASSES
        KeggEstimatorArgs.__init__(self, self.args)
        KeggDataLoader.__init__(self, self.args, self.run, self.progress)
        KeggEstimationAlgorithms.__init__(self, self.run, self.progress)

        self.name_header = None
        if self.metagenome_mode:
            self.name_header = 'contig_name'
        elif self.profile_db_path and self.collection_name and not self.metagenome_mode:
            self.name_header = 'bin_name'
        elif self.pan_db_path:
            self.name_header = 'gene_cluster_bin_name'
        else:
            self.name_header = 'genome_name'

        # update available modes and headers with appropriate genome/bin/metagenome identifier
        for m in self.available_modes:
            if m != 'modules_custom' and self.name_header not in self.available_modes[m]['headers']:
                self.available_modes[m]['headers'].insert(1, self.name_header)
            if self.metagenome_mode and self.available_modes[m]['headers'] and 'contig' in self.available_modes[m]['headers']:
                # avoid duplicate columns since contig_name is the name_header in metagenome_mode
                self.available_modes[m]['headers'].remove('contig')
        self.available_headers[self.name_header] = {
                                        'cdict_key': None,
                                        'mode_type' : 'all',
                                        'description': "Name of genome/bin/metagenome in which we find gene annotations (hits) and/or modules"
                                        }
        if self.pan_db_path:
            self.update_available_headers_for_pan()

        # INPUT OPTIONS SANITY CHECKS
        if not self.estimate_from_json and not self.contigs_db_path and self.enzymes_of_interest_df is None and not self.pan_db_path:
            raise ConfigError("NO INPUT PROVIDED. Please use the `-h` flag to see possible input options.")

        # incompatible input options
        if (self.contigs_db_path and (self.pan_db_path or self.enzymes_of_interest_df is not None)) or \
           (self.enzymes_of_interest_df is not None and self.pan_db_path):
            raise ConfigError("MULTIPLE INPUT OPTIONS DETECTED. Please check your parameters. You cannot provide more than one "
                             "of the following: a contigs database, an enzymes-txt file, or a pangenome database.")

        if self.enzymes_of_interest_df is not None:
            if None in self.enzymes_of_interest_df["enzyme_accession"].to_list():
                if self.enzymes_txt:
                    raise ConfigError("It appears that your enzymes-txt file contains one or more lines with no enzyme accession. "
                                      "Please fix this and try again.")
                else:
                    raise ConfigError("Dear programmer, the enzymes dataframe you have sent here includes enzymes with no "
                                      "accession IDs. This is no bueno.")

            if any([not enzyme_accession.startswith('K') for enzyme_accession in self.enzymes_of_interest_df["enzyme_accession"].to_list()]):
                raise ConfigError("It appears that the list of enzymes this function received includes those that do not look like "
                                  "the kind of enzyme accession IDs anvi'o is used to working with (i.e. K00001, K12345, etc). "
                                  "Please check your input.")

        if self.only_user_modules and not self.user_input_dir:
            raise ConfigError("You can only use the flag --only-user-modules if you provide a --user-modules directory.")

        self.bin_ids_to_process = None
        if self.bin_id and self.bin_ids_file:
            raise ConfigError("You have provided anvi'o with both the individual bin id %s and a file with bin ids (%s). "
                              "Please make up your mind. Which one do you want an estimate for? :)" % (self.bin_id, self.bin_ids_file))
        elif self.bin_id:
            self.bin_ids_to_process = [self.bin_id]
        elif self.bin_ids_file:
            filesnpaths.is_file_exists(self.bin_ids_file)
            self.bin_ids_to_process = [line.strip() for line in open(self.bin_ids_file).readlines()]

        # required with collection/bin input
        if (self.bin_id or self.bin_ids_file or self.collection_name) and not self.profile_db_path and not self.pan_db_path:
            raise ConfigError("You have requested metabolism estimation for a bin or set of bins, but you haven't provided "
                              "a profile database or pan database. Unfortunately, this just does not work. Please try again.")
        # required with profile db input
        if self.profile_db_path and not (self.collection_name or self.add_coverage or self.metagenome_mode):
            raise ConfigError("If you provide a profile DB, you should also provide either a collection name (to estimate metabolism "
                              "on a collection of bins) or use the --add-coverage flag (so that coverage info goes into the output "
                              "files), or both. Otherwise the profile DB is useless.")
        # required/forbidden with pangenome input
        if self.pan_db_path and not self.genomes_storage_path:
            raise ConfigError("You have provided a pan database but not its associated genomes storage database. Please give the "
                             "path to the genomes storage db using the `-g` flag.")
        if self.pan_db_path and not self.collection_name:
            raise ConfigError("You need to provide a collection name when you estimate metabolism on a pangenome. If you don't "
                              "already have a collection of gene clusters in the pan database, please make one first. Then provide "
                              "the collection to this program; you can find the collection name parameter in the INPUT #2 section "
                              "of the `-h` output.")
        if self.pan_db_path and (self.add_copy_number or self.add_coverage):
            raise ConfigError("The flags --add-copy-number or --add-coverage do not work for pangenome input.")
        # required/forbidden with JSON estimation
        if self.store_json_without_estimation and not self.json_output_file_path:
            raise ConfigError("Whoops. You seem to want to store the metabolism dictionary in a JSON file, but you haven't provided the name of that file. "
                              "Please use the --get-raw-data-as-json flag to do so.")
        if self.store_json_without_estimation and self.estimate_from_json:
            raise ConfigError("It is impossible to both estimate metabolism from JSON data and produce a JSON file without estimation at the same time... "
                              "anvi'o is judging you SO hard right now.")

        if self.profile_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)
        if self.pan_db_path:
            utils.is_pan_db_and_genomes_storage_db_compatible(self.pan_db_path, self.genomes_storage_path)

        if self.add_coverage:
            if self.enzymes_of_interest_df is None:
                if not self.profile_db_path:
                    raise ConfigError("Adding coverage values requires a profile database. Please provide one if you can. :)")
                if utils.is_blank_profile(self.profile_db_path):
                    raise ConfigError("You have provided a blank profile database, which sadly will not contain any coverage "
                                      "values, so the --add-coverage flag will not work.")

            self.add_gene_coverage_to_headers_list()

        if self.add_copy_number:
            self.available_modes["module_paths"]["headers"].extend(["num_complete_copies_of_path"])
            self.available_modes["module_steps"]["headers"].extend(["step_copy_number"])
            self.available_modes["modules"]["headers"].extend(["pathwise_copy_number", "stepwise_copy_number", "per_step_copy_numbers"])
            self.available_headers["num_complete_copies_of_path"] = {'cdict_key': None,
                                                       'mode_type': 'modules',
                                                       'description': "Number of complete copies of the path through the module"
                                                       }
            self.available_headers["step_copy_number"] = {'cdict_key': None,
                                                       'mode_type': 'modules',
                                                       'description': "Number of copies of the step"
                                                       }
            self.available_headers["pathwise_copy_number"] = {'cdict_key': None,
                                                       'mode_type': 'modules',
                                                       'description': "Pathwise module copy number, as in the maximum number of complete copies considering all the paths of highest completeness"
                                                       }
            self.available_headers["stepwise_copy_number"] = {'cdict_key': None,
                                                       'mode_type': 'modules',
                                                       'description': "Stepwise module copy number, as in the minimum copy number of all top-level steps in the module"
                                                       }
            self.available_headers["per_step_copy_numbers"] = {'cdict_key': None,
                                                       'mode_type': 'modules',
                                                       'description': "Number of copies of each top-level step in the module (the minimum of these is the stepwise module copy number)"
                                                       }

        # OUTPUT OPTIONS SANITY CHECKS
        if anvio.DEBUG:
            self.run.info("Output Modes", ", ".join(self.output_modes))
            self.run.info("Module completeness threshold", self.module_completion_threshold)
            self.run.info("Only complete modules included in output", self.only_complete)
            self.run.info("Zero-completeness modules excluded from output", self.exclude_zero_modules)
        illegal_modes = set(self.output_modes).difference(set(self.available_modes.keys()))
        if illegal_modes:
            raise ConfigError(f"You have requested some output modes that we cannot handle. The offending modes "
                              f"are: {', '.join(illegal_modes)}. Please use the flag --list-available-modes to see which ones are acceptable.")
        if self.custom_output_headers and "modules_custom" not in self.output_modes:
            raise ConfigError("You seem to have provided a list of custom headers without actually requesting a 'custom' output "
                              "mode. We think perhaps you missed something, so we are stopping you right there.")
        if "modules_custom" in self.output_modes and not self.custom_output_headers:
            raise ConfigError("You have requested a 'custom' output mode, but haven't told us what headers to include in that output. "
                              "You should be using the --custom-output-headers flag to do this.")
        if self.custom_output_headers:
            if anvio.DEBUG:
                self.run.info("Custom Output Headers", ", ".join(self.custom_output_headers))
            illegal_headers = set(self.custom_output_headers).difference(set(self.available_headers.keys()))
            if illegal_headers:
                raise ConfigError(f"You have requested some output headers that we cannot handle. The offending ones "
                                  f"are: {', '.join(illegal_headers)}. Please use the flag --list-available-output-headers to see which ones are acceptable.")

            # check if any headers requested for modules_custom mode are reserved for KOfams mode
            if "modules_custom" in self.output_modes:
                for header in self.custom_output_headers:
                    if self.available_headers[header]['mode_type'] != "modules" and self.available_headers[header]['mode_type'] != "all":
                        raise ConfigError(f"Oh dear. You requested the 'modules_custom' output mode, but gave us a header ({header}) "
                                          "that is suitable only for %s mode(s). Not good." % (self.available_headers[header]['mode_type']))

        outputs_require_ko_dict = [m for m in self.output_modes if self.available_modes[m]['data_dict'] == 'kofams']
        output_string = ", ".join(outputs_require_ko_dict)
        if self.estimate_from_json and len(outputs_require_ko_dict):
            raise ConfigError(f"You have requested to estimate metabolism from a JSON file and produce the following KOfam "
                              f"hit output mode(s): {output_string}. Unforunately, this is not possible because our JSON "
                              f"estimation function does not currently produce the required data for KOfam hit output. "
                              f"Please instead request some modules-oriented output mode(s) for your JSON input.")

        if self.matrix_format:
            raise ConfigError("You have asked for output in matrix format, but unfortunately this currently only works in "
                              "multi-mode. Please give this program an input file contining multiple bins or contigs databases instead "
                              "of the single contigs database that you have provided. We are very sorry for any inconvenience.")

        # let user know what they told anvi'o to work on
        if self.contigs_db_path:
            self.run.info("Contigs DB", self.contigs_db_path, quiet=self.quiet)
        if self.profile_db_path:
            self.run.info("Profile DB", self.profile_db_path, quiet=self.quiet)
        if self.pan_db_path:
            self.run.info("Pan DB", self.pan_db_path, quiet=self.quiet)
            self.run.info("Genomes Storage DB", self.genomes_storage_path, quiet=self.quiet)
        if self.collection_name:
            self.run.info('Collection', self.collection_name, quiet=self.quiet)
        if self.bin_id:
            self.run.info('Bin ID', self.bin_id, quiet=self.quiet)
        elif self.bin_ids_file:
            self.run.info('Bin IDs file', self.bin_ids_file, quiet=self.quiet)
        if self.enzymes_txt:
            self.run.info("Enzymes txt file", self.enzymes_txt, quiet=self.quiet)
        elif self.enzymes_of_interest_df is not None:
            self.run.info("Enzymes of interest", f"{', '.join(self.enzymes_of_interest_df['enzyme_accession'].tolist())}", quiet=self.quiet)


        estimation_mode = "Genome (or metagenome assembly)"
        if self.profile_db_path and self.collection_name:
            if not self.metagenome_mode:
                estimation_mode = "Bins in a metagenome"
            else:
                estimation_mode = "Individual contigs within a collection in a metagenome"
        elif self.metagenome_mode:
            estimation_mode = "Individual contigs in a metagenome"
        elif self.enzymes_of_interest_df is not None:
            estimation_mode = "List of enzymes"
        elif self.pan_db_path:
            estimation_mode = "Gene cluster bins in a pangenome"

        self.run.info('Mode (what we are estimating metabolism for)', estimation_mode, quiet=self.quiet)

        # a warning for high memory usage with metagenome mode in certain situations
        if self.metagenome_mode and (self.matrix_format or self.json_output_file_path):
            self.run.warning("ALERT! You are running this program in --metagenome-mode, which can have a VERY LARGE "
                             "memory footprint when used with --matrix-format or --get-raw-data-as-json, since both "
                             "of those options require storing all the per-contig data in memory. You have been warned. "
                             "The OOM-Killer may strike.")

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)
            # here we load the contigs DB just for sanity check purposes.
            # We will need to load it again later just before accessing data to avoid SQLite error that comes from different processes accessing the DB
            from anvio.dbops import ContigsDatabase # <- import here to avoid circular import
            contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
            self.contigs_db_project_name = contigs_db.meta['project_name']
        elif self.enzymes_txt:
            self.contigs_db_project_name = os.path.basename(self.enzymes_txt).replace(".", "_")
        elif self.enzymes_of_interest_df is not None:
            self.contigs_db_project_name = 'user_defined_enzymes'
        elif self.estimate_from_json:
            self.contigs_db_project_name = "json_input"
        else:
            raise ConfigError("This piece of code ended up at a place it should have never ended up at :( We need attention "
                              "from a programmer here.")

        # LOAD KEGG DATA
        if not self.only_user_modules:
            # citation output for KEGG data
            if not self.quiet:
                self.run.warning("Anvi'o will reconstruct metabolism for modules in the KEGG MODULE database, as described in "
                                 "Kanehisa and Goto et al (doi:10.1093/nar/gkr988). When you publish your findings, "
                                 "please do not forget to properly credit this work.", lc='green', header="CITATION")

            # init the enzyme accession to function definition dictionary
            # (henceforth referred to as the KO dict, even though it doesn't only contain KOs for user data)
            self.setup_ko_dict(exclude_threshold=self.exclude_kos_no_threshold)
            if self.include_stray_kos:
                self.setup_stray_ko_dict(add_entries_to_regular_ko_dict=False)
            annotation_source_set = set(['KOfam'])

            # check for kegg modules db
            if not os.path.exists(self.kegg_modules_db_path):
                raise ConfigError(f"It appears that a KEGG modules database ({self.kegg_modules_db_path}) does not exist in the provided data directory. "
                                  f"Perhaps you need to specify a different data directory using --kegg-data-dir. Or perhaps you didn't run "
                                  f"`anvi-setup-kegg-data`, though we are not sure how you got to this point in that case."
                                  f"But fine. Hopefully you now know what you need to do to make this message go away.")

            if self.contigs_db_path:
                # sanity check that contigs db was annotated with same version of MODULES.db that will be used for metabolism estimation
                if 'modules_db_hash' not in contigs_db.meta:
                    raise ConfigError("Based on the contigs DB metadata, the contigs DB that you are working with has not been annotated with hits to the "
                                      "KOfam database, so there are no KOs to estimate metabolism from. Please run `anvi-run-kegg-kofams` on this contigs DB "
                                      "before you attempt to run this script again.")
                contigs_db_mod_hash = contigs_db.meta['modules_db_hash']

                kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, args=self.args, quiet=self.quiet, run=self.run)
                mod_db_hash = kegg_modules_db.db.get_meta_value('hash')
                kegg_modules_db.disconnect()

                if contigs_db_mod_hash == "only_KOfams_were_annotated":
                    if not self.just_do_it:
                        raise ConfigError("The contigs DB that you are working with has only been annotated with KOfams, and not with a modules database. "
                                        "Since the KEGG data directory used for that annotation did not contain the modules database, we have no way of "
                                        "knowing if the set of KOfams used for annotation matches to the set of KOfams associated with your current "
                                        "modules database. Theoretically, we can still estimate metabolism even if there is a mismatch, but you risk "
                                        "getting erroneous results since 1) KOs used to define the pathways could be missing from your collection, "
                                        "and 2) KO functions could have been changed such that your KOs don't correspond to the real enzymes required "
                                        "for the pathways. If you are willing to take this risk, you can restart this program with the --just-do-it "
                                        "flag and move on with your life. But if you really want to do things properly, you should re-annotate your "
                                        "contigs database with `anvi-run-kegg-kofams`, using a KEGG data directory that includes a modules database.")
                    else:
                        self.run.warning("ALERT. ALERT. The contigs DB does not include a modules database hash, which means we can't "
                                         "tell if it was annotated with set of KOfams that match to the current modules database. Since you "
                                         "have used the --just-do-it flag, we will assume you know what you are doing. But please keep in "
                                         "mind that the metabolism estimation results could be wrong due to mismatches between the modules "
                                         "database and your set of KOfams.")

                elif contigs_db_mod_hash != mod_db_hash:
                    raise ConfigError(f"The contigs DB that you are working with has been annotated with a different version of the MODULES.db "
                                      f"than you are working with now. Basically, this means that the annotations are not compatible with the "
                                      f"metabolism data to be used for estimation. There are several ways this can happen. Please visit the "
                                      f"following URL to get help for your particular situation (copy and paste the full URL into your browser): "
                                      f"https://anvio.org/help/main/programs/anvi-estimate-metabolism/#help-im-getting-version-errors . "
                                      f"For those who need this information, the Modules DB used to annotate this contigs database has the "
                                      f"following hash: {contigs_db_mod_hash}. And the hash of the current Modules DB is: {mod_db_hash}")
        else: # USER data only
            annotation_source_set = set([])
            self.kegg_modules_db_path = None # we nullify this just in case


        # LOAD USER DATA
        if self.user_input_dir:
                # check for user modules db
                if not os.path.exists(self.user_modules_db_path):
                    raise ConfigError(f"It appears that a USER-DEFINED modules database ({self.user_modules_db_path}) does not exist in the provided data directory. "
                                      f"Perhaps you need to specify a different data directory using --user-modules. Or perhaps you didn't run "
                                      f"`anvi-setup-user-modules`. Either way, you're still awesome. Have a great day ;)")

                # sanity check that contigs db contains all necessary functional sources for user data
                user_modules_db = ModulesDatabase(self.user_modules_db_path, args=self.args, quiet=self.quiet)
                modules_db_sources = set(user_modules_db.db.get_meta_value('annotation_sources').split(','))

                if self.contigs_db_path:
                    contigs_db_sources = set(contigs_db.meta['gene_function_sources'])
                    source_in_modules_not_contigs = modules_db_sources.difference(contigs_db_sources)

                    if source_in_modules_not_contigs:
                        missing_sources = ", ".join(source_in_modules_not_contigs)
                        raise ConfigError(f"Your contigs database is missing one or more functional annotation sources that are "
                                        f"required for the modules in the database at {self.user_modules_db_path}. You will have to "
                                        f"annotate the contigs DB with these sources (or import them using `anvi-import-functions`) "
                                        f"before running this program again. Here are the missing sources: {missing_sources}")

                # expand annotation source set to include those in user db
                annotation_source_set.update(modules_db_sources)

                # we now have to add any enzymes from the user's modules db to the ko dict
                user_kos = user_modules_db.get_ko_function_dict()
                for k in user_kos:
                    if k not in self.ko_dict:
                        self.ko_dict[k] = user_kos[k]
                user_modules_db.disconnect()

        if self.contigs_db_path:
            contigs_db.disconnect()

        self.annotation_sources_to_use = list(annotation_source_set)

        # tell user what metabolism data we are using
        if self.user_input_dir:
            if self.only_user_modules:
                self.run.info('Metabolism data', "USER only")
            else:
                self.run.info('Metabolism data', "KEGG + USER-DEFINED")
        else:
            self.run.info('Metabolism data', "KEGG only")

        # sanity check for annotation sources in pangenome
        if self.genomes_storage_path:
            available_sources = DBInfo(self.genomes_storage_path).get_functional_annotation_sources()
            missing_sources = []
            for source in self.annotation_sources_to_use:
                if source not in available_sources:
                    missing_sources.append(source)
            if missing_sources:
                miss_str = ", ".join(missing_sources)
                raise ConfigError(f"The following functional annotation sources are required for metabolism "
                                  f"estimation on the chosen metabolism data, but are missing from your genome "
                                  f"storage database: {miss_str}. You'll need to figure out which genomes in the db "
                                  f"are missing those sources, and annotate them before re-making the genomes storage.")


    def update_available_headers_for_pan(self):
        """This function updates the available headers dictionary for pangenome-specific headers."""

        # in modules mode, we replace 'gene_caller_ids_in_module' with 'gene_clusters_in_module'
        self.available_headers['gene_clusters_in_module'] = {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of gene cluster IDs that contribute to a module"
                        }
        self.available_headers.pop('gene_caller_ids_in_module')

        self.available_modes['modules']["headers"] = ['gene_clusters_in_module' if x == 'gene_caller_ids_in_module' else x for x in self.available_modes['modules']["headers"]]


    def list_output_modes(self):
        """This function prints out the available output modes for the metabolism estimation script."""
        self.run.warning(None, header="AVAILABLE OUTPUT MODES", lc="green")

        for mode, mode_meta in self.available_modes.items():
            self.run.info(mode, mode_meta['description'])


    def list_output_headers(self):
        """This function prints out the available output headers for the 'custom' output mode"""
        self.run.warning(None, header="AVAILABLE OUTPUT HEADERS", lc="green")

        for header, header_meta in self.available_headers.items():
            desc_str = header_meta['description']
            type_str = header_meta['mode_type']
            mode_str = "output modes" if header_meta['mode_type'] == 'all' else "output mode"
            self.run.info(header, f"{desc_str} [{type_str} {mode_str}]")

######### ATOMIC ESTIMATION FUNCTIONS #########

    def add_gene_coverage_to_headers_list(self):
        """Updates the headers lists for relevant output modes with coverage and detection column headers.

        If a profile DB was provided, it is initialized in this function in order to get access to the sample names that will
        be part of the available coverage/detection headers.
        """

        # obtain list of sample names
        if self.enzymes_of_interest_df is not None: # in this case the input name is already determined by the initialization method
            samples_list = [self.contigs_db_project_name]

        else:
            if not self.profile_db:
                self.args.skip_consider_gene_dbs = True

                from anvio.dbops import ProfileSuperclass # <- import here to avoid circular import
                self.profile_db = ProfileSuperclass(self.args)

            samples_list = self.profile_db.p_meta['samples']

        # obtain lists of all the headers we will need to add.
        # there will be one column per sample for both coverage and detection (for individual genes and for module averages)
        kofam_hits_coverage_headers = []
        kofam_hits_detection_headers = []
        modules_coverage_headers = []
        modules_detection_headers = []

        for s in samples_list:
            # we update the available header list so that these additional headers pass the sanity checks
            kofam_hits_coverage_headers.append(s + "_coverage")
            self.available_headers[s + "_coverage"] = {'cdict_key': None,
                                                       'mode_type': 'kofams',
                                                       'description': f"Mean coverage of gene in sample {s}"
                                                       }
            kofam_hits_detection_headers.append(s + "_detection")
            self.available_headers[s + "_detection"] = {'cdict_key': None,
                                                        'mode_type': 'kofams',
                                                        'description': f"Detection of gene in sample {s}"
                                                        }
            modules_coverage_headers.extend([s + "_gene_coverages", s + "_avg_coverage"])
            self.available_headers[s + "_gene_coverages"] = {'cdict_key': None,
                                                             'mode_type': 'modules',
                                                             'description': f"Comma-separated coverage values for each gene in module in sample {s}"
                                                             }
            self.available_headers[s + "_avg_coverage"] = {'cdict_key': None,
                                                           'mode_type': 'modules',
                                                           'description': f"Average coverage of all genes in module in sample {s}"
                                                           }
            modules_detection_headers.extend([s + "_gene_detection", s + "_avg_detection"])
            self.available_headers[s + "_gene_detection"] = {'cdict_key': None,
                                                             'mode_type': 'modules',
                                                             'description': f"Comma-separated detection values for each gene in module in sample {s}"
                                                             }
            self.available_headers[s + "_avg_detection"] = {'cdict_key': None,
                                                            'mode_type': 'modules',
                                                            'description': f"Average detection of all genes in module in sample {s}"
                                                            }

        # we update the header list for the affected modes
        self.available_modes["hits"]["headers"].extend(kofam_hits_coverage_headers + kofam_hits_detection_headers)
        self.available_modes["modules"]["headers"].extend(modules_coverage_headers + modules_detection_headers)

######### ESTIMATION DRIVER FUNCTIONS #########

    def estimate_for_genome(self, kofam_gene_split_contig):
        """This is the metabolism estimation function for a contigs DB that contains a single genome.

        Assuming this contigs DB contains only one genome, it sends all of the splits and their kofam hits to the atomic
        estimation function for processing. It then returns the metabolism and ko completion dictionaries for the genome, wrapped in the superdict format.

        PARAMETERS
        ==========
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering

        RETURNS
        =======
        genome_metabolism_dict : dictionary of dictionary of dictionaries
            dictionary mapping genome name to its metabolism completeness dictionary
        genome_ko_superdict : dictionary of dictionary of dictionaries
            maps genome name to its KOfam hit dictionary
        """

        genome_metabolism_superdict = {}
        genome_ko_superdict = {}

        # since all hits belong to one genome, we can take the UNIQUE splits from all the hits
        splits_in_genome = list(set([tpl[2] for tpl in kofam_gene_split_contig]))
        metabolism_dict_for_genome, ko_dict_for_genome = self.mark_kos_present_for_list_of_splits(kofam_gene_split_contig,
                                                                                                  self.all_modules_in_db,
                                                                                                  self.all_kos_in_db,
                                                                                                  self.ko_dict,
                                                                                                  self.stray_ko_dict,
                                                                                                  self.include_stray_kos,
                                                                                                  self.exclude_kos_no_threshold,
                                                                                                  self.ignore_unknown_kos,
                                                                                                  split_list=splits_in_genome,
                                                                                                  bin_name=self.contigs_db_project_name)

        if not self.store_json_without_estimation:
            genome_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(metabolism_dict_for_genome,
                                                                                                         bin_name=self.contigs_db_project_name,
                                                                                                         all_modules_in_db=self.all_modules_in_db,
                                                                                                         module_paths_dict=self.module_paths_dict,
                                                                                                         module_completion_threshold=self.module_completion_threshold,
                                                                                                         exclude_dashed_reactions=self.exclude_dashed_reactions,
                                                                                                         add_coverage=self.add_coverage,
                                                                                                         quiet=self.quiet,
                                                                                                         genome_mode=True)
            genome_ko_superdict[self.contigs_db_project_name] = ko_dict_for_genome
        else:
            genome_metabolism_superdict[self.contigs_db_project_name] = metabolism_dict_for_genome
            genome_ko_superdict[self.contigs_db_project_name] = ko_dict_for_genome

        # append to file
        self.append_kegg_metabolism_superdicts(genome_metabolism_superdict, genome_ko_superdict)

        return genome_metabolism_superdict, genome_ko_superdict


    def estimate_for_bins_in_collection(self, kofam_gene_split_contig):
        """
        This function calls metabolism estimation for every bin the user requests.

        PARAMETERS
        ==========
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering

        RETURNS
        =======
        bins_metabolism_superdict : dictionary of dictionary of dictionaries
            dictionary mapping bin name to its metabolism completeness dictionary
        bins_ko_superdict : dictionary of dictionary of dictionaries
            dictionary mapping bin name to its KOfam hits dictionary
        """

        bins_metabolism_superdict = {}
        bins_ko_superdict = {}

        bin_name_to_split_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
        num_bins = len(bin_name_to_split_names_dict)
        self.run.info_single("%s split names associated with %s bins in collection '%s' have been "
                             "successfully recovered 🎊" % (pp(sum([len(v) for v in bin_name_to_split_names_dict.values()])),
                                                           pp(num_bins),
                                                           self.collection_name), nl_before=1, nl_after=1)

        self.progress.new("Estimating metabolism for each bin", progress_total_items=num_bins)

        for bin_name in bin_name_to_split_names_dict:
            self.progress.update("[%d of %d] %s" % (self.progress.progress_current_item + 1, num_bins, bin_name))

            splits_in_bin = bin_name_to_split_names_dict[bin_name]
            ko_in_bin = [tpl for tpl in kofam_gene_split_contig if tpl[2] in splits_in_bin]

            metabolism_dict_for_bin, ko_dict_for_bin = self.mark_kos_present_for_list_of_splits(ko_in_bin,
                                                                                                self.all_modules_in_db,
                                                                                                self.all_kos_in_db,
                                                                                                self.ko_dict,
                                                                                                self.stray_ko_dict,
                                                                                                self.include_stray_kos,
                                                                                                self.exclude_kos_no_threshold,
                                                                                                self.ignore_unknown_kos,
                                                                                                split_list=splits_in_bin,
                                                                                                bin_name=bin_name)

            if not self.store_json_without_estimation:
                bins_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(metabolism_dict_for_bin,
                                                                                       bin_name=bin_name,
                                                                                       all_modules_in_db=self.all_modules_in_db,
                                                                                       module_paths_dict=self.module_paths_dict,
                                                                                       module_completion_threshold=self.module_completion_threshold,
                                                                                       exclude_dashed_reactions=self.exclude_dashed_reactions,
                                                                                       add_coverage=self.add_coverage,
                                                                                       quiet=self.quiet,
                                                                                       genome_mode=False)
                single_bin_module_superdict = {bin_name: bins_metabolism_superdict[bin_name]}
                bins_ko_superdict[bin_name] = ko_dict_for_bin
            else:
                bins_metabolism_superdict[bin_name] = metabolism_dict_for_bin
                bins_ko_superdict[bin_name] = ko_dict_for_bin
                single_bin_module_superdict = {bin_name: metabolism_dict_for_bin}

            # append individual bin to file
            single_bin_ko_superdict = {bin_name: ko_dict_for_bin}
            self.append_kegg_metabolism_superdicts(single_bin_module_superdict, single_bin_ko_superdict)

            self.progress.increment()
            self.progress.reset()

        self.progress.end()

        return bins_metabolism_superdict, bins_ko_superdict


    def estimate_for_contigs_db_for_metagenome(self, kofam_gene_split_contig, return_superdicts=False):
        """This function handles metabolism estimation for an entire metagenome.

        We treat each contig in the metagenome to be its own 'bin' or 'genome' and estimate
        metabolism separately for each one.

        PARAMETERS
        ==========
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering
        return_superdicts : Boolean
            whether or not to return the superdicts. False by default to save on memory.

        RETURNS
        =======
        metagenome_metabolism_superdict : dictionary of dictionary of dictionaries
            dictionary mapping metagenome name to its metabolism completeness dictionary
            (will be empty dictionary if return_superdicts is False)
        metagenome_ko_superdict : dictionary of dictionary of dictionaries
            dictionary mapping metagenome name to its KOfam hits dictionary
            (will be empty dictionary if return_superdicts is False)
        """

        metagenome_metabolism_superdict = {}
        metagenome_ko_superdict = {}

        contigs_in_metagenome = list(set([tpl[3] for tpl in kofam_gene_split_contig]))
        num_contigs = len(contigs_in_metagenome)

        self.progress.new("Estimating metabolism for each contig in metagenome", progress_total_items=num_contigs)

        for contig in contigs_in_metagenome:
            self.progress.update("[%d of %d] %s" % (self.progress.progress_current_item + 1, num_contigs, contig))

            # get unique split names associated with this contig
            splits_in_contig = list(set([tpl[2] for tpl in kofam_gene_split_contig if tpl[3] == contig]))
            if anvio.DEBUG:
                self.run.info_single(f"{len(splits_in_contig)} splits recovered from contig {contig} ✌")
            ko_in_contig = [tpl for tpl in kofam_gene_split_contig if tpl[2] in splits_in_contig]
            metabolism_dict_for_contig, ko_dict_for_contig = self.mark_kos_present_for_list_of_splits(ko_in_contig,
                                                                                                      self.all_modules_in_db,
                                                                                                      self.all_kos_in_db,
                                                                                                      self.ko_dict,
                                                                                                      self.stray_ko_dict,
                                                                                                      self.include_stray_kos,
                                                                                                      self.exclude_kos_no_threshold,
                                                                                                      self.ignore_unknown_kos,
                                                                                                      split_list=splits_in_contig,
                                                                                                      bin_name=contig)


            if not self.store_json_without_estimation:
                single_contig_module_superdict = {contig: self.estimate_for_list_of_splits(metabolism_dict_for_contig,
                                                                                           bin_name=contig,
                                                                                           all_modules_in_db=self.all_modules_in_db,
                                                                                           module_paths_dict=self.module_paths_dict,
                                                                                           module_completion_threshold=self.module_completion_threshold,
                                                                                           exclude_dashed_reactions=self.exclude_dashed_reactions,
                                                                                           add_coverage=self.add_coverage,
                                                                                           quiet=self.quiet,
                                                                                           genome_mode=False)}
                if return_superdicts:
                    metagenome_metabolism_superdict[contig] = single_contig_module_superdict[contig]
                    metagenome_ko_superdict[contig] = ko_dict_for_contig
            else:
                if not return_superdicts:
                    raise ConfigError("Uh oh. Someone requested JSON-formatted estimation data from estimate_for_contigs_db_for_metagenome() "
                                      "without setting 'return_superdicts' parameter to True. ")
                metagenome_metabolism_superdict[contig] = metabolism_dict_for_contig
                metagenome_ko_superdict[contig] = ko_dict_for_contig
                single_contig_module_superdict = {contig: metabolism_dict_for_contig}


            # append individual contig to file
            single_contig_ko_superdict = {contig: ko_dict_for_contig}
            self.append_kegg_metabolism_superdicts(single_contig_module_superdict, single_contig_ko_superdict)

            self.progress.increment()
            self.progress.reset()

        self.progress.end()

        return metagenome_metabolism_superdict, metagenome_ko_superdict


    def estimate_metabolism_from_json_data(self):
        """This function runs the estimation functions on data obtained from a provided JSON file.

        Does NOT currently produce KO hits output.
        """

        self.run.info("JSON input file", self.estimate_from_json)

        filesnpaths.is_file_json_formatted(self.estimate_from_json)
        kegg_metabolism_superdict = json.load(open(self.estimate_from_json), parse_int=int)
        if ('USER' in kegg_metabolism_superdict['data_sources'] and not self.user_input_dir):
            raise ConfigError("You provided a JSON file generated from USER data, but you "
                              "did not specify which data directory to use with the `--user-modules` flag.")
        if (kegg_metabolism_superdict['data_sources'] == 'KEGG' and self.user_input_dir):
            raise ConfigError(f"You provided a JSON file generated from {kegg_metabolism_superdict['data_source']} data only, but then "
                              f"you provided us with a USER metabolism data directory. You should not use the `--user-modules` flag for this file.")

        if 'KEGG' in kegg_metabolism_superdict['data_sources']:
            kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, args=self.args, quiet=self.quiet)
            mod_db_hash = kegg_modules_db.db.get_meta_value('hash')
            kegg_modules_db.disconnect()

            if mod_db_hash != kegg_metabolism_superdict['kegg_modules_db_hash']:
                raise ConfigError(f"The modules database in the data directory you provided (or the default KEGG data directory, if you didn't "
                                  f"provide anything) has a different hash than the one used to generate this JSON input file. You probably need "
                                  f"to specify a different data directory so that we can use the modules DB with a matching hash. FYI, the hash in "
                                  f"the JSON file is {kegg_metabolism_superdict['kegg_modules_db_hash']} and the hash in the current modules DB "
                                  f"(at path `{self.kegg_modules_db_path}`) is {mod_db_hash}.")

        if self.user_input_dir:
            user_modules_db = ModulesDatabase(self.user_modules_db_path, args=self.args, quiet=self.quiet)
            mod_db_hash = user_modules_db.db.get_meta_value('hash')
            user_modules_db.disconnect()

            if mod_db_hash != kegg_metabolism_superdict['user_modules_db_hash']:
                raise ConfigError(f"The modules database in the data directory you provided with --user-modules "
                                  f"has a different hash than the one used to generate this JSON input file. You probably need "
                                  f"to specify a different data directory so that we can use the modules DB with a matching hash. FYI, the hash in "
                                  f"the JSON file is {kegg_metabolism_superdict['user_modules_db_hash']} and the hash in the current modules DB "
                                  f"(at path `{self.user_modules_db_path}`) is {mod_db_hash}.")

        new_kegg_metabolism_superdict = {}

        expected_keys_for_module = {"gene_caller_ids", "kofam_hits", "genes_to_contigs", "contigs_to_genes"}
        bins_found = []
        additional_keys = set([])

        self.init_data_from_modules_db()

        for bin_name, meta_dict_for_bin in kegg_metabolism_superdict.items():
            if bin_name in ['data_sources', 'kegg_modules_db_hash', 'user_modules_db_hash']:
                continue
            else:
                bins_found.append(bin_name)

            for mod, mod_dict in meta_dict_for_bin.items():
                # verify that dict contains the necessary keys for estimation
                if not expected_keys_for_module.issubset(set(mod_dict.keys())):
                    missing_keys = expected_keys_for_module.difference(set(mod_dict.keys()))
                    raise ConfigError("Your JSON file is incorrectly formatted for metabolism estimation. We expect the following keys: %s. "
                                      "However, we didn't find some of them for module %s in %s. Here are the missing keys: %s"
                                      % (expected_keys_for_module, mod, bin_name, missing_keys))

                additional_keys = additional_keys.union(set(mod_dict.keys()).difference(expected_keys_for_module))

                # convert certain lists to sets
                mod_dict['gene_caller_ids'] = set(mod_dict['gene_caller_ids'])
                for contig, gene_list in mod_dict['contigs_to_genes'].items():
                    mod_dict['contigs_to_genes'][contig] = set(gene_list)
                mod_dict['genes_to_contigs'] = {int(g):c for g,c in mod_dict['genes_to_contigs'].items()}
                mod_dict['warnings'] = set(mod_dict['warnings'])

            new_kegg_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(meta_dict_for_bin,
                                                                                       bin_name=bin_name,
                                                                                       all_modules_in_db=self.all_modules_in_db,
                                                                                       module_paths_dict=self.module_paths_dict,
                                                                                       module_completion_threshold=self.module_completion_threshold,
                                                                                       exclude_dashed_reactions=self.exclude_dashed_reactions,
                                                                                       add_coverage=self.add_coverage,
                                                                                       quiet=self.quiet,
                                                                                       genome_mode=False)
            single_bin_module_superdict = {bin_name: new_kegg_metabolism_superdict[bin_name]}
            self.append_kegg_metabolism_superdicts(single_bin_module_superdict, ko_superdict_for_list_of_splits={})


        if not self.quiet and additional_keys:
            self.run.warning("Just to let you know, we found the following module-level keys in your JSON file that were totally ignored during metabolism estimation "
                             "(no harm was done by including them): %s" % (additional_keys))

        self.run.info("Bins/genomes/metagenomes found", ", ".join(bins_found))
        return new_kegg_metabolism_superdict


    def estimate_metabolism_for_enzymes_of_interest(self):
        """Estimates metabolism on a set of enzymes provided by the user.

        This function assumes that all enzymes in the file are coming from a single genome, and is effectively the
        same as the estimate_for_genome() function.

        Requires the self.enzymes_of_interest_df attribute to have been established (either through self.enzymes_txt parameter,
        or by populating the self.enzymes_of_interest_df upon initialization of the args class -- see the KeggEstimatorArgs
        for details). In this mode, we make fake splits and contigs to match the expected input to the atomic functions, and
        the contigs_db_project_name attribute has been set (previously) to the name of the enzyme txt file (OR simply to
        'user_defined_enzymes', depending on the initialization).

        RETURNS
        =======
        enzyme_metabolism_superdict : dictionary of dictionary of dictionaries
            dictionary mapping the name of the enzyme txt file to its metabolism completeness dictionary
        enzyme_ko_superdict : dictionary of dictionary of dictionaries
            maps the name of the enzyme txt file to its KOfam hit dictionary
        """

        kofam_gene_split_contig = []
        # no splits or contigs here
        for gene_call_id, ko in zip(self.enzymes_of_interest_df["gene_id"], self.enzymes_of_interest_df["enzyme_accession"]):
            kofam_gene_split_contig.append((ko,gene_call_id,"NA","NA"))

        enzyme_metabolism_superdict = {}
        enzyme_ko_superdict = {}

        metabolism_dict_for_genome,ko_dict_for_genome = self.mark_kos_present_for_list_of_splits(kofam_gene_split_contig,
                                                                                                 self.all_modules_in_db,
                                                                                                 self.all_kos_in_db,
                                                                                                 self.ko_dict,
                                                                                                 self.stray_ko_dict,
                                                                                                 self.include_stray_kos,
                                                                                                 self.exclude_kos_no_threshold,
                                                                                                 self.ignore_unknown_kos,
                                                                                                 bin_name=self.contigs_db_project_name)

        if not self.store_json_without_estimation:
            enzyme_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(metabolism_dict_for_genome,
                                                                                                         bin_name=self.contigs_db_project_name,
                                                                                                         all_modules_in_db=self.all_modules_in_db,
                                                                                                         module_paths_dict=self.module_paths_dict,
                                                                                                         module_completion_threshold=self.module_completion_threshold,
                                                                                                         exclude_dashed_reactions=self.exclude_dashed_reactions,
                                                                                                         add_coverage=self.add_coverage,
                                                                                                         quiet=self.quiet,
                                                                                                         genome_mode=True)
            enzyme_ko_superdict[self.contigs_db_project_name] = ko_dict_for_genome
        else:
            enzyme_metabolism_superdict[self.contigs_db_project_name] = metabolism_dict_for_genome
            enzyme_ko_superdict[self.contigs_db_project_name] = ko_dict_for_genome

        # append to file
        self.append_kegg_metabolism_superdicts(enzyme_metabolism_superdict, enzyme_ko_superdict)

        return enzyme_metabolism_superdict, enzyme_ko_superdict


    def estimate_metabolism_for_pangenome_bins(self, enzyme_cluster_split_contig, cluster_collection):
        """Estimates metabolism individually on each bin in a pangenome.

        PARAMETERS
        ==========
        enzyme_cluster_split_contig : list
            (enzyme_accession, gene_cluster_id, split, contig) tuples in which split and contig are both NAs

        cluster_collection : dictionary
            maps bin names in the collection to the list of gene clusters in each bin
        """

        gc_bins_metabolism_superdict = {}
        gc_bins_ko_superdict = {}
        num_bins = len(cluster_collection)

        self.progress.new("Estimating metabolism for each bin of gene clusters", progress_total_items=num_bins)

        for bin_name, gc_list in cluster_collection.items():
            self.progress.update("[%d of %d] %s" % (self.progress.progress_current_item + 1, num_bins, bin_name))

            enzymes_in_bin = [tpl for tpl in enzyme_cluster_split_contig if tpl[1] in gc_list]
            metabolism_dict_for_bin, ko_dict_for_bin = self.mark_kos_present_for_list_of_splits(enzymes_in_bin,
                                                                                                self.all_modules_in_db,
                                                                                                self.all_kos_in_db,
                                                                                                self.ko_dict,
                                                                                                self.stray_ko_dict,
                                                                                                self.include_stray_kos,
                                                                                                self.exclude_kos_no_threshold,
                                                                                                self.ignore_unknown_kos,
                                                                                                bin_name=bin_name)

            if not self.store_json_without_estimation:
                gc_bins_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(metabolism_dict_for_bin,
                                                                                          bin_name=bin_name,
                                                                                          all_modules_in_db=self.all_modules_in_db,
                                                                                          module_paths_dict=self.module_paths_dict,
                                                                                          module_completion_threshold=self.module_completion_threshold,
                                                                                          exclude_dashed_reactions=self.exclude_dashed_reactions,
                                                                                          add_coverage=self.add_coverage,
                                                                                          quiet=self.quiet,
                                                                                          genome_mode=False)
                single_bin_module_superdict = {bin_name: gc_bins_metabolism_superdict[bin_name]}
                gc_bins_ko_superdict[bin_name] = ko_dict_for_bin
            else:
                gc_bins_metabolism_superdict[bin_name] = metabolism_dict_for_bin
                single_bin_module_superdict = {bin_name: metabolism_dict_for_bin}
                gc_bins_ko_superdict[bin_name] = ko_dict_for_bin

            # append individual bin to file
            single_bin_ko_superdict = {bin_name: ko_dict_for_bin}
            self.append_kegg_metabolism_superdicts(single_bin_module_superdict, single_bin_ko_superdict)

            self.progress.increment()
            self.progress.reset()

        self.progress.end()

        return gc_bins_metabolism_superdict, gc_bins_ko_superdict


    def estimate_metabolism(self, skip_storing_data=False, output_files_dictionary=None, return_superdicts=False,
                            return_subset_for_matrix_format=False, all_modules_in_db=None, all_kos_in_db=None,
                            module_paths_dict=None, prune_superdicts=False):
        """This is the driver function for estimating metabolism for a single contigs DB.

        It will decide what to do based on whether the input DB is a genome, metagenome, or pangenome.
        It usually avoids returning the metabolism data to save on memory (as this data is typically appended to
        files immediately), but this behavior can be changed by setting return_superdicts to True (for the entire
        modules/ko superdictionaries) or return_subset_for_matrix_format to True (for a subset of these dicts that
        multi-estimators need for matrix output generation).

        PARAMETERS
        ==========
        skip_storing_data : boolean
            set to True if we don't want the metabolism data dictionary to be stored as a file (useful when using this function
            for on-the-fly visualization or for generating matrix format output from a multi estimator class)
        output_files_dictionary : dictionary of mode, AppendableFile object pairs
            contains an initialized AppendableFile object to append output to for each output mode
            (used in multi-mode to direct all output from several estimators to the same files)
        return_superdicts : boolean
            set to True if you want the kegg_metabolism_superdict and kofam_hits_superdict to be returned.
            we don't return these by default to save on memory
        return_subset_for_matrix_format : boolean
            set to True if you want subsets of the superdicts to be returned: one subdict for module completeness scores, one
            subdict for module presence/absence, and one subdict for KO hits. Used for matrix format output.
        all_modules_in_db : dictionary
            if this function is called from the KeggMetabolismEstimatorMulti class, this parameter contains the module information
            loaded from the MODULES.db in init_data_from_modules_db(). Otherwise, it is None and this function will have to call
            init_data_from_modules_db()
        all_kos_in_db : dictionary
            This is the same deal as the all_modules_in_db param - it should only have a value if passed from the
            KeggMetabolismEstimatorMulti class
        module_paths_dict : dictionary
            Again, same thing as all_modules_in_db param - only provided if passed from KeggMetabolismEstimatorMulti
        prune_superdicts : bool
            when True, the function will prune superdicts to ensure they contain data only for metabolic modules
            whose `pathwise_percent_complete` score is over 0.0 and enzymes that are represented by at least one gene

        RETURNS
        =======
        kegg_metabolism_superdict : dictionary of dictionaries of dictionaries
            a complex data structure containing the metabolism estimation data for each genome/bin in the contigs DB
            (only returned if return_superdicts is True)
        kofam_hits_superdict : dictionary of dictionaries of dictionaries
            a complex data structure containing the KOfam hits information for each genome/bin in the contigs DB
            (only returned if return_superdicts is True)
        """

        kegg_metabolism_superdict = {}
        kofam_hits_superdict = {}

        if skip_storing_data or self.write_dict_to_json:
            self.output_file_dict = {} # if this object is empty, no output will be generated
        else:
            if output_files_dictionary:
                self.output_file_dict = output_files_dictionary
            else:
                self.output_file_dict = self.setup_output_for_appending()

        # recovering superdicts for metabolic modules and kofams
        if self.estimate_from_json:
            kegg_metabolism_superdict = self.estimate_metabolism_from_json_data()
        else:
            # we either get the modules DB info from the previous class, or we have to initialize it here (unless that already happened)
            if all_modules_in_db:
                self.all_modules_in_db = all_modules_in_db
                self.all_kos_in_db = all_kos_in_db
                self.module_paths_dict = module_paths_dict
            elif not self.all_modules_in_db:
                self.init_data_from_modules_db()

            if self.enzymes_of_interest_df is not None:
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_metabolism_for_enzymes_of_interest()
            elif self.pan_db_path:
                gene_cluster_collections = ccollections.Collections()
                gene_cluster_collections.populate_collections_dict(self.pan_db_path)
                if self.collection_name not in gene_cluster_collections.collections_dict:
                    c_str = ', '.join(gene_cluster_collections.collections_dict.keys())
                    raise ConfigError(f"The collection name you provided ('{self.collection_name}') is not valid for this "
                                      f"Pan DB. Here are the collections in this database: {c_str}")
                collection_dict = gene_cluster_collections.get_collection_dict(self.collection_name)

                all_gene_clusters_in_collection = []
                for bin_name, gene_cluster_list in collection_dict.items():
                    all_gene_clusters_in_collection += gene_cluster_list

                kofam_hits_info = self.init_hits_for_pangenome(self.pan_db_path, self.genomes_storage_path,
                                                             all_gene_clusters_in_collection, self.annotation_sources_to_use)
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_metabolism_for_pangenome_bins(kofam_hits_info, collection_dict)
            else:
                kofam_hits_info = self.init_hits_and_splits(self.contigs_db_path, self.profile_db_path,
                                                          self.collection_name, self.annotation_sources_to_use,
                                                          user_input_dir=self.user_input_dir,
                                                          all_kos_in_db=self.all_kos_in_db)

                if self.add_coverage:
                    self.init_gene_coverage(self.profile_db, gcids_for_kofam_hits={int(tpl[1]) for tpl in kofam_hits_info})

                if self.profile_db_path and self.collection_name and not self.metagenome_mode:
                    kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_for_bins_in_collection(kofam_hits_info)
                elif not self.collection_name and not self.metagenome_mode:
                    self.genome_mode = True
                    kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_for_genome(kofam_hits_info)
                elif self.metagenome_mode:
                    kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_for_contigs_db_for_metagenome(kofam_hits_info, return_superdicts=return_superdicts)
                else:
                    raise ConfigError("This class doesn't know how to deal with that yet :/")

        # now we have our superdicts, the first order of buisness is to check if the user is interested in pruned versions of
        # these dicts
        if prune_superdicts:
            # the next two lines simply redefine the contents of (1) kegg_metabolism_superdict based on the second level
            # entries with pathwise_percent_complete > 0.0, and (2) kofam_hits_superdict based on the second level
            # entries with more than one gene call:
            kegg_metabolism_superdict = {source: {module_name: values for module_name, values in entries.items() if values.get("pathwise_percent_complete") > 0.0} for source, entries in kegg_metabolism_superdict.items()}
            kofam_hits_superdict = {source: {kofam_name: values for kofam_name, values in entries.items() if len(values["gene_caller_ids"]) > 0} for source, entries in kofam_hits_superdict.items()}

        # we take care of the JSON output if requested
        if self.write_dict_to_json:
            self.store_metabolism_superdict_as_json(kegg_metabolism_superdict, self.json_output_file_path)

        # more housekeeping for outputs
        if not self.multi_mode:
            for mode, file_object in self.output_file_dict.items():
                file_object.close()

        # 'returning stuff' stage
        if return_superdicts:
            # if the programmer asked the superdicts to be returned, we will now extend them with metabolic
            # module NAME and CLASS information since it can't hurt to have those for downstream analyses.
            for sample_name in kegg_metabolism_superdict:
                for module_name in kegg_metabolism_superdict[sample_name]:
                    kegg_metabolism_superdict[sample_name][module_name]['NAME'] = self.all_modules_in_db[module_name]['NAME']
                    kegg_metabolism_superdict[sample_name][module_name]['CLASS'] = self.all_modules_in_db[module_name]['CLASS']

            return kegg_metabolism_superdict, kofam_hits_superdict
        elif return_subset_for_matrix_format:
            # if we are generating matrix output, we need a limited subset of this data downstream
            # so in this case, we can extract and return smaller dictionaries for module completeness,
            # module presence/absence, and KO hits.
            return self.generate_subsets_for_matrix_format(kegg_metabolism_superdict, kofam_hits_superdict, only_complete_modules=self.only_complete)
        else:
          # otherwise we return nothing at all
          return

######### OUTPUT DICTIONARY FUNCTIONS #########

    def add_common_elements_to_output_dict_for_module_in_bin(self, bin, mnum, c_dict, headers_to_include, headers_in_c_dict, d):
        """This function fills in the provided modules dictionary with data common to all module-related output modes.

        It's designed to be called from the function generate_output_dict_for_modules(), which sets up the modules dictionary
        and handles the self.modules_unique_id key. For that reason, it is best understood by taking a look at that function
        first. Most parameters are named to be consistent with the variables in that function.

        The provided modules dictionary will be modified in-place.

        PARAMETERS
        ==========
        bin : str
            current bin name
        mnum : str
            current module number
        c_dict : dictionary
            fourth level of module completion dictionary corresponding to current module
        headers_to_include : list
            which headers to include in the output dictionary
        headers_in_c_dict : list
            headers that we can include directly from the c_dict without further processing
        d : dict
            the modules output dictionary that needs to be added to
        """

        # fetch module info
        metadata_dict = self.get_module_metadata_dictionary(mnum)
        definition_list = self.all_modules_in_db[mnum]["DEFINITION"]
        if not isinstance(definition_list, list):
            definition_list = [definition_list]
        module_def = '"' + " ".join(definition_list) + '"'

        # top-level keys and keys not in superdict
        if self.name_header in headers_to_include:
            d[self.modules_unique_id][self.name_header] = bin
        if "db_name" in headers_to_include:
            d[self.modules_unique_id]["db_name"] = self.database_name
        if 'module' in headers_to_include:
            d[self.modules_unique_id]['module'] = mnum

        # module-specific info
        if "module_name" in headers_to_include:
            d[self.modules_unique_id]["module_name"] = metadata_dict["module_name"]
        if "module_class" in headers_to_include:
            d[self.modules_unique_id]["module_class"] = metadata_dict["module_class"]
        if "module_category" in headers_to_include:
            d[self.modules_unique_id]["module_category"] = metadata_dict["module_category"]
        if "module_subcategory" in headers_to_include:
            d[self.modules_unique_id]["module_subcategory"] = metadata_dict["module_subcategory"]
        if "module_definition" in headers_to_include:
            d[self.modules_unique_id]["module_definition"] = module_def
        if "module_substrates" in headers_to_include:
            if self.all_modules_in_db[mnum]['substrate_list']:
                d[self.modules_unique_id]["module_substrates"] = ",".join(self.all_modules_in_db[mnum]['substrate_list'])
            else:
                d[self.modules_unique_id]["module_substrates"] = "None"
        if "module_products" in headers_to_include:
            if self.all_modules_in_db[mnum]['product_list']:
                d[self.modules_unique_id]["module_products"] = ",".join(self.all_modules_in_db[mnum]['product_list'])
            else:
                d[self.modules_unique_id]["module_products"] = "None"
        if "module_intermediates" in headers_to_include:
            if self.all_modules_in_db[mnum]['intermediate_list']:
                d[self.modules_unique_id]["module_intermediates"] = ",".join(self.all_modules_in_db[mnum]['intermediate_list'])
            else:
                d[self.modules_unique_id]["module_intermediates"] = "None"

        # comma-separated lists of KOs and gene calls in module
        kos_in_mod = sorted(c_dict['kofam_hits'].keys())
        # gene call list should be in same order as KO list
        gcids_in_mod = []
        kos_in_mod_list = []
        if kos_in_mod:
            for ko in kos_in_mod:
                gcids_in_mod += [str(x) for x in c_dict["kofam_hits"][ko]]
                kos_in_mod_list += [ko for x in c_dict["kofam_hits"][ko]]
        if "enzyme_hits_in_module" in headers_to_include:
            d[self.modules_unique_id]["enzyme_hits_in_module"] = ",".join(kos_in_mod_list)
        if "gene_caller_ids_in_module" in headers_to_include:
            d[self.modules_unique_id]["gene_caller_ids_in_module"] = ",".join(gcids_in_mod)
        if "gene_clusters_in_module" in headers_to_include:
            d[self.modules_unique_id]["gene_clusters_in_module"] = ",".join(gcids_in_mod)

        # comma-separated list of warnings
        if "warnings" in headers_to_include:
            if not c_dict["warnings"]:
                d[self.modules_unique_id]["warnings"] = "None"
            else:
                d[self.modules_unique_id]["warnings"] = ",".join(sorted(c_dict["warnings"]))

        # list of enzymes unique to module
        unique_enzymes = sorted(list(c_dict["unique_to_this_module"]))
        if "enzymes_unique_to_module" in headers_to_include:
            if unique_enzymes:
                d[self.modules_unique_id]["enzymes_unique_to_module"] = ",".join(unique_enzymes)
            else:
                d[self.modules_unique_id]["enzymes_unique_to_module"] = "No enzymes unique to module"
        if "unique_enzymes_hit_counts" in headers_to_include:
            if unique_enzymes:
                hit_count_list = []
                for e in unique_enzymes:
                    hit_count_list.append(str(len(c_dict["kofam_hits"][e])))
                d[self.modules_unique_id]["unique_enzymes_hit_counts"] = ",".join(hit_count_list)
            else:
                d[self.modules_unique_id]["unique_enzymes_hit_counts"] = "NA"

        # everything else at c_dict level
        for h in headers_in_c_dict:
            if h not in self.available_headers.keys():
                raise ConfigError("Requested header %s not available." % (h))
            h_cdict_key = self.available_headers[h]['cdict_key']
            if not h_cdict_key:
                raise ConfigError("We don't know the corresponding key in metabolism completeness dict for header %s." % (h))

            value = c_dict[h_cdict_key]
            if isinstance(value, list):
                if not value:
                    value = "None"
                else:
                    value = ",".join(value)
            d[self.modules_unique_id][h] = value

        # add module copy number if requested
        if self.add_copy_number:
            # pathwise: we take the maximum copy number of all the paths of highest completeness
            if "pathwise_copy_number" in headers_to_include:
                d[self.modules_unique_id]["pathwise_copy_number"] = c_dict["pathwise_copy_number"]

            # stepwise copy number
            if "stepwise_copy_number" in headers_to_include:
                d[self.modules_unique_id]["stepwise_copy_number"] = c_dict["stepwise_copy_number"]
            if "per_step_copy_numbers" in headers_to_include:
                step_copy_numbers = []
                for step in c_dict["top_level_step_info"]:
                    step_copy_numbers.append(str(c_dict["top_level_step_info"][step]["copy_number"]))
                d[self.modules_unique_id]["per_step_copy_numbers"] = ",".join(step_copy_numbers)

        # add coverage if requested
        if self.add_coverage:
            for s in self.coverage_sample_list:
                sample_cov_header = s + "_gene_coverages"
                sample_det_header = s + "_gene_detection"
                sample_avg_cov_header = s + "_avg_coverage"
                sample_avg_det_header = s + "_avg_detection"

                gene_coverages_in_mod = []
                gene_detection_in_mod = []
                for gc in gcids_in_mod:
                    if self.enzymes_of_interest_df is not None:
                        gc_idx = gc
                    else:
                        gc_idx = int(gc)
                    gene_coverages_in_mod.append(c_dict["genes_to_coverage"][s][gc_idx])
                    gene_detection_in_mod.append(c_dict["genes_to_detection"][s][gc_idx])

                d[self.modules_unique_id][sample_cov_header] = ",".join([str(c) for c in gene_coverages_in_mod])
                d[self.modules_unique_id][sample_det_header] = ",".join([str(d) for d in gene_detection_in_mod])
                d[self.modules_unique_id][sample_avg_cov_header] = c_dict["average_coverage_per_sample"][s]
                d[self.modules_unique_id][sample_avg_det_header] = c_dict["average_detection_per_sample"][s]


    def generate_output_dict_for_modules(self, kegg_superdict, headers_to_include=None, only_complete_modules=False,
                                               exclude_zero_completeness=True):
        """This dictionary converts the metabolism superdict to a two-level dict containing desired headers for output.

        The metabolism superdict is a three-to-four-level dictionary. The first three levels are: genomes/metagenomes/bins, modules, and module completion information.
        The module completion dictionary also has some dictionaries in it, and those make up the fourth level.
        The structure of the module completion dictionary is like this example:
        {mnum: {"gene_caller_ids": set([132, 133, 431, 6777])
                "kofam_hits": {'K00033' : [431, 6777],
                                'K01057' : [133],
                                'K00036' : [132] },
                "genes_to_contigs": {132: 0,
                                     133: 0,
                                     431: 2,
                                    6777: 1 },
                "contigs_to_genes": { 0: set([132, 133]),
                                      1: set(6777),
                                      2: set(431) },}
                "pathway_completeness":     [0.66, 0.66, ...]
                "present_nonessential_kos":      []
                "most_complete_paths":           [['K00033','K01057','K02222'], ['K00033','K01057','K00036'], ...]
                "pathwise_percent_complete":              0.66
                "pathwise_is_complete":                      False
                (.....)
                                      }

        To distill this information into one line, we need to convert the dictionary on-the-fly to a dict of dicts,
        where each bin-module, bin-module-path, or bin-module-step (depending on output mode) is keyed by an arbitrary integer.
        There will be a lot of redundant information in the rows.

        PARAMETERS
        ==========
        kegg_superdict : dictionary of dictionaries of dictionaries
            The metabolism superdict containing KO hit and KEGG module information for each bin/genome/metagenome

        headers_to_include : list
            Which headers to include in the output dictionary

        only_complete_modules : boolean
            If True, we only put information into the output dictionary for modules whose completeness is above the threshold

        exclude_zero_completeness : boolean
            If True, we don't put modules with a 0 completeness score in the dictionary

        RETURNS
        =======
        d : dictionary of dictionaries
            The output dictionary whose format is compatible for printing to a tab-delimited file
        """

        if not headers_to_include:
            headers_to_include = set(OUTPUT_MODES['modules']['headers'])
        else:
            headers_to_include = set(headers_to_include)

        # make sure all requested headers are available
        avail_headers = set(self.available_headers.keys())
        illegal_headers = headers_to_include.difference(avail_headers)
        if illegal_headers:
            raise ConfigError("Some unavailable headers were requested. These include: %s" % (", ".join(illegal_headers)))

        module_level_headers = set(["module_name", "module_class", "module_category", "module_subcategory", "module_definition",
                                    "module_substrates", "module_products", "module_intermediates", "warnings", "enzymes_unique_to_module",
                                    "unique_enzymes_hit_counts"])
        path_level_headers = set(["path_id", "path", "path_completeness", "num_complete_copies_of_path", "annotated_enzymes_in_path"])
        step_level_headers = set(["step_id", "step", "step_completeness", "step_copy_number"])

        requested_path_info = headers_to_include.intersection(path_level_headers)
        requested_step_info = headers_to_include.intersection(step_level_headers)
        if requested_path_info and requested_step_info:
            raise ConfigError(f"Oh, bother. It seems you have requested both path-level headers and step-level headers for your modules-related "
                              f"output. Unfortunately, these two types of information are incompatible and cannot be put into the same output "
                              f"file due to the way we internally organize the data. Sorry. If you want both types of information, we recommend "
                              f"requesting both in separate files using `--output-modes module_paths,module_steps`. If you absolutely need custom "
                              f"formatting for the output files, then you can run this program twice and each time provide either only path-level "
                              f"headers or only step-level headers to the `--custom-output-headers` flag. To help you out, here are the path-level "
                              f"headers that you requested: {', '.join(requested_path_info)}. And here are the step-level headers that you requested: "
                              f"{', '.join(requested_step_info)}")

        keys_not_in_superdict = set([h for h in self.available_headers.keys() if self.available_headers[h]['cdict_key'] is None])

        remaining_headers = headers_to_include.difference(keys_not_in_superdict)
        remaining_headers = remaining_headers.difference(module_level_headers)
        remaining_headers = remaining_headers.difference(path_level_headers)
        remaining_headers = remaining_headers.difference(step_level_headers)

        # convert to two-level dict where unique id keys for a dictionary of information for each bin/module pair
        d = {}
        if not self.modules_unique_id:
             self.modules_unique_id = 0

        for bin, mod_dict in kegg_superdict.items():
            for mnum, c_dict in mod_dict.items():
                if anvio.DEBUG:
                    self.run.info("Generating output for module", mnum)

                if only_complete_modules and not (c_dict["pathwise_is_complete"] or c_dict["stepwise_is_complete"]):
                    continue
                if exclude_zero_completeness and (c_dict["pathwise_percent_complete"] == 0 and c_dict["stepwise_completeness"] == 0):
                    continue

                # handle path-level information (ie, for module_paths mode)
                if headers_to_include.intersection(path_level_headers):
                    for p_index, p in enumerate(self.module_paths_dict[mnum]):
                        d[self.modules_unique_id] = {}
                        self.add_common_elements_to_output_dict_for_module_in_bin(bin, mnum, c_dict, headers_to_include, remaining_headers, d)

                        # path-specific info
                        if "path_id" in headers_to_include:
                            d[self.modules_unique_id]["path_id"] = p_index
                        if "path" in headers_to_include:
                            d[self.modules_unique_id]["path"] = ",".join(p)
                        if "path_completeness" in headers_to_include:
                            d[self.modules_unique_id]["path_completeness"] = c_dict["pathway_completeness"][p_index]
                        if "annotated_enzymes_in_path" in headers_to_include:
                            annotated = []
                            for accession in p:
                                # handle enzyme components
                                if '+' in accession or '-' in accession:
                                    components = re.split(r'\+|\-', accession)
                                    for c in components:
                                        if c in c_dict['kofam_hits'].keys():
                                            annotated.append(c)
                                        else:
                                            annotated.append(f"[MISSING {c}]")
                                else:
                                    if (accession in self.all_modules_in_db and mod_dict[accession]["pathwise_is_complete"]) or \
                                    (accession in c_dict['kofam_hits'].keys()):
                                        annotated.append(accession)
                                    else:
                                        annotated.append(f"[MISSING {accession}]")
                            d[self.modules_unique_id]["annotated_enzymes_in_path"] = ",".join(annotated)

                        # add path-level redundancy if requested
                        if self.add_copy_number:
                            d[self.modules_unique_id]["num_complete_copies_of_path"] = c_dict["num_complete_copies_of_all_paths"][p_index]

                        self.modules_unique_id += 1

                # handle step-level information (ie, for module_steps mode)
                elif headers_to_include.intersection(step_level_headers):
                    for s_index, step_dict in c_dict["top_level_step_info"].items():
                        d[self.modules_unique_id] = {}
                        self.add_common_elements_to_output_dict_for_module_in_bin(bin, mnum, c_dict, headers_to_include, remaining_headers, d)

                        # step-specific info
                        if "step_id" in headers_to_include:
                            d[self.modules_unique_id]["step_id"] = s_index
                        if "step" in headers_to_include:
                            d[self.modules_unique_id]["step"] = step_dict["step_definition"]
                        if "step_completeness" in headers_to_include:
                            d[self.modules_unique_id]["step_completeness"] = 1 if step_dict["complete"] else 0

                        # add step-level redundancy if requested
                        if self.add_copy_number:
                            d[self.modules_unique_id]["step_copy_number"] = step_dict["copy_number"]

                        self.modules_unique_id += 1

                # handle module-level information (ie, for modules mode)
                else:
                    d[self.modules_unique_id] = {}
                    self.add_common_elements_to_output_dict_for_module_in_bin(bin, mnum, c_dict, headers_to_include, remaining_headers, d)

                    self.modules_unique_id += 1

        return d


    def generate_output_dict_for_kofams(self, ko_superdict, headers_to_include=None):
        """This dictionary converts the kofam hits superdict to a two-level dict containing desired headers for output.

        The kofam hits superdict is a three-to-four-level dictionary. The first three levels are: genomes/metagenomes/bins, KOs, and KO hit information.
        The KO hit dictionary also has some dictionaries in it, and those make up the fourth level.
        The structure of the KO hit dictionary is like this example:
        {ko: {"gene_caller_ids" : set([431, 6777]),
              "modules" : ["M00001", "M00555"],                 **Can be None if KO does not belong to any KEGG modules
              "genes_to_contigs": { 431: 2,
                                   6777: 1 },
              "contigs_to_genes": { 1: set(6777),
                                    2: set(431) }}}

        To distill this information into one line, we need to convert the dictionary on-the-fly to a dict of dicts,
        where each bin-ko-gene_caller_id is keyed by an arbitrary integer.

        PARAMETERS
        ==========
        ko_superdict : dictionary of dictionaries of dictionaries
            The metabolism superdict containing all KO hits in each bin/genome/metagenome

        headers_to_include : list
            Which headers to include in the output dictionary

        RETURNS
        =======
        d : dictionary of dictionaries
            The output dictionary whose format is compatible for printing to a tab-delimited file
        """

        # use the kofam_hits output mode header set by default
        if not headers_to_include:
            headers_to_include = set(OUTPUT_MODES["hits"]["headers"])
        else:
            headers_to_include = set(headers_to_include)

        d = {}
        if not self.ko_unique_id:
            self.ko_unique_id = 0

        for bin, ko_dict in ko_superdict.items():
            for ko, k_dict in ko_dict.items():
                if anvio.DEBUG:
                    self.run.info("Generating output for KO", ko)

                metadata_dict = self.get_ko_metadata_dictionary(ko, dont_fail_if_not_found=True)

                for gc_id in k_dict["gene_caller_ids"]:
                    d[self.ko_unique_id] = {}

                    if self.name_header in headers_to_include:
                        d[self.ko_unique_id][self.name_header] = bin
                    if "db_name" in headers_to_include:
                        d[self.ko_unique_id]["db_name"] = self.database_name
                    if "enzyme" in headers_to_include:
                        d[self.ko_unique_id]["enzyme"] = ko
                    if "gene_caller_id" in headers_to_include:
                        d[self.ko_unique_id]["gene_caller_id"] = gc_id
                    if "contig" in headers_to_include:
                        d[self.ko_unique_id]["contig"] = k_dict["genes_to_contigs"][gc_id]
                    if "modules_with_enzyme" in headers_to_include:
                        d[self.ko_unique_id]["modules_with_enzyme"] = metadata_dict["modules_with_enzyme"]
                    if "enzyme_definition" in headers_to_include:
                        d[self.ko_unique_id]["enzyme_definition"] = metadata_dict["enzyme_definition"]
                    if "warnings" in headers_to_include:
                        if not k_dict["warnings"]:
                            d[self.ko_unique_id]["warnings"] = "None"
                        else:
                            d[self.ko_unique_id]["warnings"] = ",".join(sorted(k_dict["warnings"]))

                    if self.add_coverage:
                        if self.enzymes_of_interest_df is not None:
                            for s in self.coverage_sample_list:
                                sample_cov_header = s + "_coverage"
                                d[self.ko_unique_id][sample_cov_header] = self.enzymes_of_interest_df[self.enzymes_of_interest_df["gene_id"] == gc_id]["coverage"].values[0]
                                sample_det_header = s + "_detection"
                                d[self.ko_unique_id][sample_det_header] = self.enzymes_of_interest_df[self.enzymes_of_interest_df["gene_id"] == gc_id]["detection"].values[0]
                        else:
                            if not self.profile_db:
                                raise ConfigError("We're sorry that you came all this way, but it seems the profile database has "
                                                  "not been initialized, therefore we cannot add coverage values to your output. "
                                                  "This is likely a bug or developer mistake. It is a sad day :(")

                            for s in self.coverage_sample_list:
                                sample_cov_header = s + "_coverage"
                                d[self.ko_unique_id][sample_cov_header] = self.profile_db.gene_level_coverage_stats_dict[gc_id][s]['mean_coverage']
                                sample_det_header = s + "_detection"
                                d[self.ko_unique_id][sample_det_header] = self.profile_db.gene_level_coverage_stats_dict[gc_id][s]['detection']

                    self.ko_unique_id += 1

        return d


    def generate_subsets_for_matrix_format(self, module_superdict, ko_hits_superdict, only_complete_modules=False):
        """Here we extract and return subsets of data from the superdicts, for matrix formatted output.

        The subsets of data that we need are: module completeness scores, module presence/absence, KO hit frequency,
                                              module top-level step completeness
        If --add-copy-number was provided, we also need copy numbers for modules and module steps.

        Each of these is put into a dictionary (one for modules, one for ko hits, one for module steps) and returned.

        PARAMETERS
        ==========
        module_superdict, ko_hits_superdict: multi-level dictionaries
            The superdicts containing module/enzyme information for each bin/genome/metagenome

        only_complete_modules : boolean
            If True, we only put information into the output dictionary for modules whose completeness is above the threshold
            in at least one sample in the matrix. Note that this affects all module-related info, not just completeness scores.
        """

        mod_completeness_presence_subdict = {}
        ko_hits_subdict = {}
        steps_subdict = {}

        if only_complete_modules: # here we determine which modules to include in the output
            set_of_modules_to_include = set([])

            for mod in self.all_modules_in_db:
                for bin, mod_dict in module_superdict.items():
                    if mod in mod_dict and (mod_dict[mod]["pathwise_is_complete"] or mod_dict[mod]["stepwise_is_complete"]):
                        set_of_modules_to_include.add(mod)
                        break

        for bin, mod_dict in module_superdict.items():
            mod_completeness_presence_subdict[bin] = {}
            steps_subdict[bin] = {}
            for mnum, c_dict in mod_dict.items():

                if only_complete_modules:
                    if mnum not in set_of_modules_to_include:
                        continue

                mod_completeness_presence_subdict[bin][mnum] = {}
                mod_completeness_presence_subdict[bin][mnum]["pathwise_percent_complete"] = c_dict["pathwise_percent_complete"]
                mod_completeness_presence_subdict[bin][mnum]["pathwise_is_complete"] = c_dict["pathwise_is_complete"]
                mod_completeness_presence_subdict[bin][mnum]["stepwise_completeness"] = c_dict["stepwise_completeness"]
                mod_completeness_presence_subdict[bin][mnum]["stepwise_is_complete"] = c_dict["stepwise_is_complete"]

                if self.add_copy_number:
                    if c_dict["num_complete_copies_of_most_complete_paths"]:
                        mod_completeness_presence_subdict[bin][mnum]['pathwise_copy_number'] = max(c_dict["num_complete_copies_of_most_complete_paths"])
                    else:
                        mod_completeness_presence_subdict[bin][mnum]['pathwise_copy_number'] = 'NA'
                    mod_completeness_presence_subdict[bin][mnum]['stepwise_copy_number'] = c_dict["stepwise_copy_number"]

                for step_id, step_dict in c_dict["top_level_step_info"].items():
                    step_key = mnum + "_" + f"{step_id:02d}"
                    steps_subdict[bin][step_key] = {}
                    steps_subdict[bin][step_key]["step_is_complete"] = step_dict["complete"]
                    # we include step metadata in this dictionary directly (rather than trying to access it out later using a function)
                    steps_subdict[bin][step_key]["step_definition"] = step_dict["step_definition"]
                    if self.add_copy_number:
                        steps_subdict[bin][step_key]["step_copy_number"] = step_dict["copy_number"]

        for bin, ko_dict in ko_hits_superdict.items():
            ko_hits_subdict[bin] = {}
            for knum, k_dict in ko_dict.items():
                ko_hits_subdict[bin][knum] = {}
                ko_hits_subdict[bin][knum]['num_hits'] = len(k_dict['gene_caller_ids']) # number of hits to this KO in the bin

        return mod_completeness_presence_subdict, ko_hits_subdict, steps_subdict

######### OUTPUT GENERATION FUNCTIONS #########

    def append_kegg_metabolism_superdicts(self, module_superdict_for_list_of_splits, ko_superdict_for_list_of_splits):
        """This function appends the metabolism superdicts (for a single genome, bin, or contig in metagenome) to existing files
        for each output mode.

        It appends to the initialized AppendableFile objects in self.output_file_dict.

        This is an alternative to store_kegg_metabolism_superdicts(), which prints the entire metabolism superdicts for all
        genomes/bins/contigs in metagenome at once.

        NOTE: in cases where output should not be generated (ie, self.skip_storing_data is true), the self.output_file_dict
        object will be empty, so the loop in this function will never run. :)
        """

        for mode, file_obj in self.output_file_dict.items():
            header_list = self.available_modes[mode]["headers"]
            if not header_list:
                raise ConfigError("Oh, dear. You've come all this way only to realize that we don't know which headers to use "
                                  "for the %s output mode. Something is terribly wrong, and it is probably a developer's fault. :("
                                  % (mode))
            if self.available_modes[mode]["data_dict"] == 'modules':
                output_dict = self.generate_output_dict_for_modules(module_superdict_for_list_of_splits, headers_to_include=header_list, \
                                                                    only_complete_modules=self.only_complete, \
                                                                    exclude_zero_completeness=self.exclude_zero_modules)
            elif self.available_modes[mode]["data_dict"] == 'kofams':
                output_dict = self.generate_output_dict_for_kofams(ko_superdict_for_list_of_splits, headers_to_include=header_list)
            else:
                raise ConfigError(f"Uh oh. You've requested to generate output from the {self.available_modes[mode]['data_dict']} "
                                  "data dictionary, but we don't know about that one.")

            file_obj.append(output_dict, headers=['key'] + header_list, do_not_write_key_column=True)

            if anvio.DEBUG:
                self.run.warning(f"Appended metabolism dictionary to {file_obj.path}" ,
                                header='DEBUG OUTPUT', lc='yellow')


    def store_kegg_metabolism_superdicts(self, module_superdict, ko_superdict):
        """This function writes the metabolism superdicts (in their entirety) to tab-delimited files depending
        on which output the user requested.

        The user can request a variety of output 'modes', and for each of these modes we look up the details on the output
        format which are stored in self.available_modes, use that information to generate a dictionary of dictionaries,
        and store that dictionary as a tab-delimited file.

        This is an alternative to append_kegg_metabolism_superdicts(), which adds to the output files
        one genome/bin/contig in metagenome at a time for better memory management.
        """

        for mode in self.output_modes:
            output_path = self.output_file_prefix + "_" + self.available_modes[mode]["output_suffix"]
            header_list = self.available_modes[mode]["headers"]
            if not header_list:
                raise ConfigError("Oh, dear. You've come all this way only to realize that we don't know which headers to use "
                                  "for the %s output mode. Something is terribly wrong, and it is probably a developer's fault. :("
                                  % (mode))
            if self.available_modes[mode]["data_dict"] == 'modules':
                output_dict = self.generate_output_dict_for_modules(module_superdict, headers_to_include=header_list, \
                                                                    only_complete_modules=self.only_complete, \
                                                                    exclude_zero_completeness=self.exclude_zero_modules)
            elif self.available_modes[mode]["data_dict"] == 'kofams':
                output_dict = self.generate_output_dict_for_kofams(ko_superdict, headers_to_include=header_list)
            else:
                raise ConfigError(f"Uh oh. You've requested to generate output from the {self.available_modes[mode]['data_dict']} "
                                  "data dictionary, but we don't know about that one.")
            utils.store_dict_as_TAB_delimited_file(output_dict, output_path, headers=['key'] + header_list, do_not_write_key_column=True)
            self.run.info("%s output file" % mode, output_path, nl_before=1)


    def store_metabolism_superdict_as_json(self, kegg_superdict, file_path):
        """This function writes the metabolism superdict into one json file."""

        def set_to_list(obj):
            if isinstance(obj, set):
                return list(obj)

        if self.user_input_dir:
            if self.only_user_modules:
                kegg_superdict['data_sources'] = 'USER'
            else:
                kegg_superdict['data_sources'] = 'USER,KEGG'
            user_modules_db = ModulesDatabase(self.user_modules_db_path, args=self.args, quiet=self.quiet)
            kegg_superdict['user_modules_db_hash'] = user_modules_db.db.get_meta_value('hash')
            user_modules_db.disconnect()
        else:
            kegg_superdict['data_sources'] = 'KEGG'
            kegg_superdict['user_modules_db_hash'] = ''

        if self.only_user_modules:
            kegg_superdict['kegg_modules_db_hash'] = ''
        else:
            kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, args=self.args, quiet=self.quiet)
            kegg_superdict['kegg_modules_db_hash'] = kegg_modules_db.db.get_meta_value('hash')
            kegg_modules_db.disconnect()

        filesnpaths.is_output_file_writable(file_path)
        open(file_path, 'w').write(json.dumps(kegg_superdict, indent=4, default=set_to_list))
        self.run.info("JSON Output", file_path)

######### INTERACTIVE VISUALIZATION FUNCTIONS #########

    def get_metabolism_data_for_visualization(self):
        """Returns a dictionary of metabolism data for visualization on the interactive interface.

        This function should be called from the interactive interface class to obtain metabolism data.
        It runs the metabolism estimation function on-the-fly to generate the data.
        It then pulls only certain keys from the resulting dictionary and returns them to the interface.
        """

        # add keys to this list to include the data in the visualization dictionary
        module_data_keys_for_visualization = ['pathwise_percent_complete']

        metabolism_dict, ko_hit_dict = self.estimate_metabolism(skip_storing_data=True, return_superdicts=True)
        data_for_visualization = {}

        for bin, mod_dict in metabolism_dict.items():
            data_for_visualization[bin] = {}
            for mnum, c_dict in mod_dict.items():
                data_for_visualization[bin][mnum] = {}
                for key, value in c_dict.items():
                    if key in module_data_keys_for_visualization:
                        data_for_visualization[bin][mnum][key] = value

        return data_for_visualization


class KeggMetabolismEstimatorMulti(KeggEstimatorArgs, KeggDataLoader):
    """Class for reconstructing/estimating metabolism for multiple contigs DBs at a time.

    Iterates through the provided DBs and estimates metabolism for each one using the KeggMetabolismEstimator class.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-estimate-metabolism
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        if args.contigs_db:
            raise ConfigError("You appear to have provided both an input text file and a contigs database, and "
                              "now anvi'o is not quite sure which one to work on. Please choose only one. :) ")

        # init the base classes for access to shared paths and such
        KeggEstimatorArgs.__init__(self, self.args)
        KeggDataLoader.__init__(self, self.args, run, progress)

        # This can be initialized later if necessary by setup_ko_dict()
        self.ko_dict = {}

        if anvio.DEBUG:
            self.run.info("Completeness threshold: multi estimator", self.module_completion_threshold)
            self.run.info("Output Modes", ", ".join(self.output_modes))
            self.run.info("Matrix format", self.matrix_format)
            self.run.info("Matrix will include metadata", self.matrix_include_metadata)
            if self.module_specific_matrices:
                self.run.info("Matrices for specific modules: ", ", ".join(self.module_specific_matrices))

        self.databases = None

        # INPUT SANITY CHECKS
        if (self.external_genomes_file and (self.internal_genomes_file or self.metagenomes_file)) \
            or (self.internal_genomes_file and (self.external_genomes_file or self.metagenomes_file)) \
            or (self.metagenomes_file and (self.external_genomes_file or self.internal_genomes_file)):
                raise ConfigError("Multiple file inputs were provided. Please choose only one at a time to make "
                                  "things easier on everybody.")

        if args.estimate_from_json or args.store_json_without_estimation or args.get_raw_data_as_json:
            raise ConfigError("You've provided some JSON parameters. We are sorry to say that these parameters don't "
                              "work for input files with multiple contigs DBs. :( ")

        if self.only_user_modules and not self.user_input_dir:
            raise ConfigError("You can only use the flag --only-user-modules if you provide a --user-modules directory.")

        # OUTPUT SANITY CHECKS
        if self.matrix_format and self.long_format_mode:
            raise ConfigError("Please request EITHER long-format output modes OR matrix format. When you ask for both "
                              "like this, anvi'o is confused. :) ")
        if self.matrix_include_metadata and not self.matrix_format:
            raise ConfigError("The option --include-metadata is only available when you also use the flag --matrix-format "
                              "to get matrix output. :) Plz try again.")
        if self.module_specific_matrices and not self.matrix_format:
            raise ConfigError("The option --module-specific-matrices is only available when you also use the flag --matrix-format "
                              "to get matrix output. :) Plz try again.")
        if self.matrix_format:
            for stat in ['completeness', 'presence', 'ko_hits']:
                matrix_output_file = '%s-%s-MATRIX.txt' % (self.output_file_prefix, stat)
                filesnpaths.is_output_file_writable(matrix_output_file, ok_if_exists=False)

        # set name header
        if self.metagenomes_file:
            self.name_header = 'contig_name'
        elif self.external_genomes_file:
            self.name_header = 'genome_name'
        elif self.internal_genomes_file:
            self.name_header = 'bin_name'

        # LOAD KEGG DATA
        if not self.only_user_modules:
            # citation output for KEGG data
            if not self.quiet:
                self.run.warning("Anvi'o will reconstruct metabolism for modules in the KEGG MODULE database, as described in "
                                 "Kanehisa and Goto et al (doi:10.1093/nar/gkr988). When you publish your findings, "
                                 "please do not forget to properly credit this work.", lc='green', header="CITATION")

            # init the enzyme accession to function definition dictionary
            # (henceforth referred to as the KO dict, even though it doesn't only contain KOs for user data)
            self.setup_ko_dict()
            annotation_source_set = set(['KOfam'])

            # check for kegg modules db
            if not os.path.exists(self.kegg_modules_db_path):
                raise ConfigError(f"It appears that a KEGG modules database ({self.kegg_modules_db_path}) does not exist in the provided data directory. "
                                  f"Perhaps you need to specify a different data directory using --kegg-data-dir. Or perhaps you didn't run "
                                  f"`anvi-setup-kegg-data`, though we are not sure how you got to this point in that case."
                                  f"But fine. Hopefully you now know what you need to do to make this message go away.")

        else: # USER data only
            annotation_source_set = set([])
            self.kegg_modules_db_path = None # we nullify this just in case


        # LOAD USER DATA
        if self.user_input_dir:
            # check for user modules db
            if not os.path.exists(self.user_modules_db_path):
                raise ConfigError(f"It appears that a USER-DEFINED modules database ({self.user_modules_db_path}) does not exist in the provided data directory. "
                                  f"Perhaps you need to specify a different data directory using --user-modules. Or perhaps you didn't run "
                                  f"`anvi-setup-user-modules`. Either way, you're still awesome. Have a great day ;)")

            # expand annotation source set to include those in user db
            user_modules_db = ModulesDatabase(self.user_modules_db_path, args=self.args, quiet=self.quiet)
            modules_db_sources = set(user_modules_db.db.get_meta_value('annotation_sources').split(','))
            annotation_source_set.update(modules_db_sources)

            # we now have to add any enzymes from the user's modules db to the ko dict
            user_kos = user_modules_db.get_ko_function_dict()
            for k in user_kos:
                if k not in self.ko_dict:
                    self.ko_dict[k] = user_kos[k]
            user_modules_db.disconnect()

        self.annotation_sources_to_use = list(annotation_source_set)

        # tell user what metabolism data we are using
        if self.user_input_dir:
            if self.only_user_modules:
                self.run.info('Metabolism data', "USER only")
            else:
                self.run.info('Metabolism data', "KEGG + USER-DEFINED")
        else:
            self.run.info('Metabolism data', "KEGG only")


    def list_output_modes(self):
        """This function prints out the available output modes for the metabolism estimation script."""

        self.run.warning(None, header="AVAILABLE OUTPUT MODES", lc="green")

        for mode, mode_meta in self.available_modes.items():
            self.run.info(mode, mode_meta['description'])


    def update_available_headers_for_multi(self):
        """This function updates the available headers dictionary to reflect all possibilities in the multiple DB case."""

        self.available_headers["db_name"] = {
                                        'cdict_key': None,
                                        'mode_type': 'all',
                                        'description': "Name of contigs DB. Always included in multi-mode output (so no need to specify on the command line)"
        }
        if self.name_header == 'genome_name':
            self.available_headers["genome_name"] = {
                                            'cdict_key': None,
                                            'mode_type': 'all',
                                            'description': "Name of genome in which we find gene annotations (hits) and/or modules"
                                            }
        elif self.name_header == 'bin_name':
            self.available_headers["bin_name"] = {
                                            'cdict_key': None,
                                            'mode_type': 'all',
                                            'description': "Name of bin in which we find gene annotations (hits) and/or modules"
                                            }
        elif self.name_header == 'contig_name':
            self.available_headers["contig_name"] = {
                                            'cdict_key': None,
                                            'mode_type': 'all',
                                            'description': "Name of contig (in a metagenome) in which we find gene annotations (hits) and/or modules"
                                            }

        # here we make sure db_name is always included in the multi-mode output
        for mode in self.available_modes:
            if self.available_modes[mode]["headers"] and "db_name" not in self.available_modes[mode]["headers"]:
                self.available_modes[mode]["headers"].insert(1, "db_name")


    def list_output_headers(self):
        """This function prints out the available output headers for the 'custom' output mode"""

        # include all possibilities for genome/bin/metagenome name in the output since we don't know which cases
        # we will find in the metagenomes file
        self.update_available_headers_for_multi()

        self.run.warning("Just so you know, if you used the flags --add-copy-number or --add-coverage, you "
                         "won't see the possible headers for these data in the list below. If you want to "
                         "include this information in a custom output file and need to know which headers "
                         "you can choose from, you should re-run this command with --list-available-output-headers "
                         "on a SINGLE sample from your input file.")

        self.run.warning(None, header="AVAILABLE OUTPUT HEADERS", lc="green")

        for header, header_meta in self.available_headers.items():
            desc_str = header_meta['description']
            type_str = header_meta['mode_type']
            mode_str = "output modes" if header_meta['mode_type'] == 'all' else "output mode"
            self.run.info(header, f"{desc_str} [{type_str} {mode_str}]")

######### DRIVER ESTIMATION FUNCTIONS -- MULTI #########

    def init_metagenomes(self):
        """This function parses the input metagenomes file and adjusts class attributes as needed"""

        g = MetagenomeDescriptions(self.args, run=self.run, progress=self.progress, enforce_single_profiles=False)
        g.load_metagenome_descriptions(skip_functions=(not self.matrix_format))

        # sanity check that all dbs are properly annotated with required sources
        for src in self.annotation_sources_to_use:
            bad_metagenomes = [v['name'] for v in g.metagenomes.values() if not v['gene_function_sources'] or src not in v['gene_function_sources']]
            if len(bad_metagenomes):
                bad_metagenomes_txt = [f"'{bad}'" for bad in bad_metagenomes]
                n = len(bad_metagenomes)
                it_or_them = P('it', n, alt='them')
                raise ConfigError(f"Bad news :/ It seems {n} of your {P('metagenome', len(g.metagenomes))} "
                                  f"{P('is', n, alt='are')} lacking any function annotations for "
                                  f"`{src}`. This means you either need to annotate {it_or_them} by running the appropriate "
                                  f"annotation program on {it_or_them}, import functional annotations into {it_or_them} from this source using "
                                  f"`anvi-import-functions`, or remove {it_or_them} from your internal and/or external genomes files "
                                  f"before re-running `anvi-estimate-metabolism. Here is the list of offenders: "
                                  f"{', '.join(bad_metagenomes_txt)}.")

        if self.matrix_format:
            for name in g.metagenomes:
                gene_functions_in_genome_dict, _, _= g.get_functions_and_sequences_dicts_from_contigs_db(name, requested_source_list=self.annotation_sources_to_use, return_only_functions=True)
                # reminder, an entry in gene_functions_in_genome_dict looks like this:
                # 4264: {'KOfam': None, 'COG20_FUNCTION': None, 'UpxZ': ('PF06603.14', 'UpxZ', 3.5e-53)}
                for gcid, func_dict in gene_functions_in_genome_dict.items():
                    for source, func_tuple in func_dict.items():
                        if func_tuple:
                            acc_string, func_def, evalue = func_tuple
                            for acc in acc_string.split('!!!'):
                                if acc not in self.ko_dict:
                                    self.ko_dict[acc] = {'definition': func_def}

        # enforce metagenome mode
        if not self.metagenome_mode:
            self.metagenome_mode = True

        self.databases = copy.deepcopy(g.metagenomes)
        self.database_names = copy.deepcopy(g.metagenome_names)


    def init_external_internal_genomes(self):
        """This function parses the input internal/external genomes file and adjusts class attributes as needed"""

        g = GenomeDescriptions(self.args, run=self.run, progress=progress_quiet)
        g.load_genomes_descriptions(skip_functions=(not self.matrix_format), init=False)

        # sanity check that all dbs are properly annotated with required sources
        for src in self.annotation_sources_to_use:
            bad_genomes = [v['name'] for v in g.genomes.values() if not v['gene_function_sources'] or src not in v['gene_function_sources']]
            if len(bad_genomes):
                bad_genomes_txt = [f"'{bad_genome}'" for bad_genome in bad_genomes]
                it_or_them = P('it', len(bad_genomes), alt='them')
                raise ConfigError(f"Bad news :/ It seems {len(bad_genomes)} of your {P('genome', len(g.genomes))} "
                                  f"{P('is', len(bad_genomes), alt='are')} lacking any function annotations for "
                                  f"`{src}`. This means you either need to annotate {it_or_them} by running the appropriate "
                                  f"annotation program on {it_or_them}, import functional annotations into {it_or_them} from this source using "
                                  f"`anvi-import-functions`, or remove {it_or_them} from your internal and/or external genomes files "
                                  f"before re-running `anvi-estimate-metabolism. Here is the list of offenders: "
                                  f"{', '.join(bad_genomes_txt)}.")

        if self.matrix_format:
            for genome_name in g.genomes:
                gene_functions_in_genome_dict, _, _= g.get_functions_and_sequences_dicts_from_contigs_db(genome_name, requested_source_list=self.annotation_sources_to_use, return_only_functions=True)
                # reminder, an entry in gene_functions_in_genome_dict looks like this:
                # 4264: {'KOfam': None, 'COG20_FUNCTION': None, 'UpxZ': ('PF06603.14', 'UpxZ', 3.5e-53)}
                for gcid, func_dict in gene_functions_in_genome_dict.items():
                    for source, func_tuple in func_dict.items():
                        if func_tuple:
                            acc_string, func_def, eval = func_tuple
                            for acc in acc_string.split('!!!'):
                                if acc not in self.ko_dict:
                                    self.ko_dict[acc] = {'definition': func_def}

        # metagenome mode must be off
        if self.metagenome_mode:
            self.metagenome_mode = False

        self.databases = copy.deepcopy(g.genomes)
        if self.external_genomes_file:
            self.database_names = copy.deepcopy(g.external_genome_names)
        else:
            self.database_names = copy.deepcopy(g.internal_genome_names)


    def get_args_for_single_estimator(self, db_name):
        """Returns args formatted for an instance of KeggMetabolismEstimator that will work on a contigs DB. Very tricksy.

        PARAMETERS
        ==========
        db_name : string
            the name of the contigs DB that the estimator will work on

        RETURNS
        =======
        args : Namespace object
            a set of arguments modified for use by a single estimator
        """

        args = KeggEstimatorArgs(self.args, format_args_for_single_estimator=True)

        if db_name not in self.databases:
            raise ConfigError("We cannot initialize a single estimator for the contigs DB '%s' because it is not in the metagenomes dictionary."
                              % (db_name))

        args.contigs_db = self.databases[db_name]['contigs_db_path']
        if 'profile_db_path' in self.databases[db_name]:
            args.profile_db = self.databases[db_name]['profile_db_path']
        if 'collection_id' in self.databases[db_name]:
            args.collection_name = self.databases[db_name]['collection_id']
        if 'bin_id' in self.databases[db_name]:
            args.bin_id = self.databases[db_name]['bin_id']

        args.metagenome_mode = self.metagenome_mode
        args.quiet = True
        args.database_name = db_name
        args.multi_mode = True
        args.include_metadata = self.matrix_include_metadata
        args.user_modules = self.user_input_dir or None
        args.only_user_modules = self.only_user_modules

        self.update_available_headers_for_multi()

        if anvio.DEBUG:
            self.run.info("Completeness threshold: single estimator", args.module_completion_threshold)
            self.run.info("Database name: single estimator", args.database_name)
            self.run.info("Matrix metadata: single estimator", args.matrix_include_metadata)
            self.run.info("User data directory: single estimator", args.user_modules)

        return args


    def get_metabolism_superdict_multi(self):
        """The function that calls metabolism on each individual contigs db.

         If we need matrix format output, it aggregates the results into one dictionary for modules, one for KOs,
         and one for steps. These dictionaries are returned.

         (Otherwise, empty dictionaries are returned.)
         """

        metabolism_super_dict = {}
        ko_hits_super_dict = {}
        module_steps_super_dict = {}

        total_num_metagenomes = len(self.database_names)
        self.progress.new("Estimating metabolism for contigs DBs", progress_total_items=total_num_metagenomes)

        if not self.matrix_format:
            files_dict = self.setup_output_for_appending()

        for metagenome_name in self.database_names:
            args = self.get_args_for_single_estimator(metagenome_name)
            self.progress.update("[%d of %d] %s" % (self.progress.progress_current_item + 1, total_num_metagenomes, metagenome_name))
            if not self.matrix_format:
                KeggMetabolismEstimator(args, progress=progress_quiet, run=run_quiet).estimate_metabolism(output_files_dictionary=files_dict)
            else:
                metabolism_super_dict[metagenome_name], ko_hits_super_dict[metagenome_name], module_steps_super_dict[metagenome_name] = KeggMetabolismEstimator(args, \
                                                                                                    progress=progress_quiet, \
                                                                                                    run=run_quiet).estimate_metabolism(skip_storing_data=True, \
                                                                                                    return_subset_for_matrix_format=True, \
                                                                                                    all_modules_in_db=self.all_modules_in_db, \
                                                                                                    all_kos_in_db=self.all_kos_in_db, \
                                                                                                    module_paths_dict=self.module_paths_dict)

            self.progress.increment()
            self.progress.reset()

        self.progress.end()

        if not self.matrix_format:
            for mode, file in files_dict.items():
                file.close()

        return metabolism_super_dict, ko_hits_super_dict, module_steps_super_dict


    def estimate_metabolism(self):
        """A driver function to run metabolism estimation on each provided contigs DB."""

        if not self.databases:
            self.progress.new("Initializing contigs DBs")
            self.progress.update("...")
            if self.metagenomes_file:
                self.progress.reset()
                self.run.info("Metagenomes file", self.metagenomes_file)
                self.init_metagenomes()
            elif self.external_genomes_file:
                self.progress.reset()
                self.run.info("External genomes file", self.external_genomes_file)
                self.init_external_internal_genomes()
            elif self.internal_genomes_file:
                self.progress.reset()
                self.run.info("Internal genomes file", self.internal_genomes_file)
                self.init_external_internal_genomes()
            else:
                self.progress.reset()
                raise ConfigError("Whooops. We are not sure how you got to this point without an input file, "
                                  "but you did, and now we have to crash becasue we cannot estimate metabolism "
                                  "without inputs. :/")
            self.progress.end()
            self.run.info("Num Contigs DBs in file", len(self.database_names))
            self.run.info('Metagenome Mode', self.metagenome_mode)

        self.init_data_from_modules_db()

        # these will be empty dictionaries unless matrix format
        kegg_metabolism_superdict_multi, ko_hits_superdict_multi, module_steps_superdict_multi = self.get_metabolism_superdict_multi()

        if self.matrix_format:
            self.store_metabolism_superdict_multi_matrix_format(kegg_metabolism_superdict_multi, ko_hits_superdict_multi, module_steps_superdict_multi)


######### OUTPUT GENERATION FUNCTIONS -- MULTI #########

    def write_stat_to_matrix(self, stat_name, stat_header, stat_key, stat_dict, item_list, stat_metadata_headers,
                             write_rows_with_all_zeros=False, comment_dictionary=None):
        """A generic function to write a statistic to a matrix file.

        Accesses the provided stat_dict and writes the statistic to a tab-delimited matrix file.
        Should work for module completeness, module presence/absence, and ko hits.

        PARAMETERS
        ==========
        stat_name : str
            Which statistic we are working on, ie "completeness". Used in output file name.
        stat_header : str
            The header for the items reporting this statistic in the matrix output, ie 'module' or 'ko'.
        stat_key : str
            The key used to access this statistic in the stat_dict
        stat_dict : dictionary
            A multi-level dictionary (a subset of metabolism estimation output) in which the statistic and
            relevant metadata can be found.
        item_list : list of str
            The row (item) names of the matrix. Ideally would be sorted for consistency in the output.
        stat_metadata_headers : list of str
            A list of the headers for metadata columns (which must also be keys for this metadata in the stat_dict)
            that will be included in the matrix output if self.matrix_include_metadata is True
        write_rows_with_all_zeros : boolean
            If true, rows with all zeros are included in the matrix. Otherwise we leave those out.
        comment_dictionary : dictionary
            A dictionary in which the item is a (str) comment and the key is the INDEX of the corresponding item (from item_list)
            that this comment should be printed before. When we reach this item in the list, the comment str will be
            printed (after a '#' character) before printing the item's line. Trailing "\n" should not be in the comment
            str but if this needs to be a multi-line comment internal "\n# " strings should separate each line.
        """

        output_file_path = '%s-%s-MATRIX.txt' % (self.output_file_prefix, stat_name)

        sample_list = list(stat_dict.keys())
        # here we figure out if there is more than one bin to work with in any given sample
        sample_columns = []
        sample_bin_list = {}
        for s in sample_list:
            bins = list(stat_dict[s].keys())
            bins.sort()
            sample_bin_list[s] = bins
            if len(bins) > 1:
                for b in bins:
                    sample_columns.append(s + "_" + b)
            else:
                sample_columns.append(s)

        if self.matrix_include_metadata:
            cols = [stat_header] + stat_metadata_headers + sample_columns
        else:
            cols = [stat_header] + sample_columns

        # we could be fancier with this, but we are not that cool
        with open(output_file_path, 'w') as output:
            output.write('\t'.join(cols) + '\n')
            cur_index = 0
            for m in item_list:
                # write comment, if necessary
                if comment_dictionary and cur_index in comment_dictionary:
                    comment_line = "# " + comment_dictionary[cur_index] + "\n"
                    output.write(comment_line)

                line = [m]

                if self.matrix_include_metadata:
                    if stat_header == 'module':
                        metadata_dict = self.get_module_metadata_dictionary(m)
                    elif stat_header == 'enzyme':
                        metadata_dict = self.get_ko_metadata_dictionary(m, dont_fail_if_not_found=True)
                    elif stat_header == 'module_step':
                        # get the step metadata from the first sample and first bin (b/c it will be the same in all of them)
                        first_samp = sample_list[0]
                        first_bin = list(stat_dict[first_samp].keys())[0]
                        metadata_dict = {"step_definition": stat_dict[first_samp][first_bin][m]["step_definition"]}
                    else:
                        raise ConfigError(f"write_stat_to_matrix() speaking. I need to access metadata for {stat_header} "
                                          "statistics but there is no function defined for this.")

                    for h in stat_metadata_headers:
                        if h not in metadata_dict.keys():
                            raise ConfigError(f"We couldn't find the key '{h}' in the metadata dictionary for {stat_header}s. "
                                              "Please check that your metadata accessor function obtains this data.")
                        line.append(metadata_dict[h])

                for s in sample_list:
                    bins = sample_bin_list[s]

                    for b in bins:
                        # if its not in the dict, we know it is zero
                        if m not in stat_dict[s][b].keys():
                            value = 0
                        else:
                            value = stat_dict[s][b][m][stat_key]

                        # handle presence/absence values as integers
                        if isinstance(value, bool):
                            line.append(int(value))
                        else:
                            line.append(value)

                if not write_rows_with_all_zeros:
                    only_numbers = [n for n in line if isinstance(n, (int, float))]
                    if sum(only_numbers) == 0:
                        continue

                output.write('\t'.join([str(f) for f in line]) + '\n')
                cur_index += 1

        self.run.info('Output matrix for "%s"' % stat_name, output_file_path)


    def store_metabolism_superdict_multi_matrix_format(self, module_superdict_multi, ko_superdict_multi, steps_superdict_multi):
        """Stores the multi-contigs DB metabolism data in several matrices.

        Contigs DBs are arranged in columns and KEGG modules/KOs are arranged in rows.
        Each module statistic (ie, completeness, presence/absence) will be in a different file.

        The parameters to this function are superdictionaries where each top-level entry is one
        of the multi-mode sample inputs (ie a metagenome, internal, or external genome) and its
        corresponding value comes from running estimate_metabolism() with return_subset_for_matrix_format=True.

        That is:
         module % completeness = module_superdict_multi[sample][bin][mnum]['pathwise_percent_complete']
         module is complete = module_superdict_multi[sample][bin][mnum]["pathwise_is_complete"]
         # hits for KO = ko_superdict_multi[sample][bin][knum]['num_hits']
         module step is complete = module_steps_superdict_multi[sample][bin][mnum_stepnum]['step_is_complete']

        If self.matrix_include_metadata was True, these superdicts will also include relevant metadata.
        """

        include_zeros = not self.exclude_zero_modules

        # module stats that each will be put in separate matrix file
        # key is the stat, value is the corresponding header in superdict
        module_matrix_stats = {"module_pathwise_completeness" : "pathwise_percent_complete",
                               "module_pathwise_presence" : "pathwise_is_complete",
                               "module_stepwise_completeness" : "stepwise_completeness",
                               "module_stepwise_presence" : "stepwise_is_complete",
                               }
        module_step_matrix_stats = {"step_completeness" : "step_is_complete"}
        if self.add_copy_number:
            module_matrix_stats["module_pathwise_copy_number"] = "pathwise_copy_number"
            module_matrix_stats["module_stepwise_copy_number"] = "stepwise_copy_number"
            module_step_matrix_stats["step_copy_number"] = "step_copy_number"

        # all samples/bins have the same modules in the dict so we can pull the item list from the first pair
        first_sample = list(module_superdict_multi.keys())[0]
        first_bin = list(module_superdict_multi[first_sample].keys())[0]
        module_list = list(module_superdict_multi[first_sample][first_bin].keys())
        module_list.sort()

        for stat, key in module_matrix_stats.items():
            self.write_stat_to_matrix(stat_name=stat, stat_header='module', stat_key=key, stat_dict=module_superdict_multi, \
                                      item_list=module_list, stat_metadata_headers=MODULE_METADATA_HEADERS, \
                                      write_rows_with_all_zeros=include_zeros)

        module_step_list = list(steps_superdict_multi[first_sample][first_bin].keys())
        module_step_list.sort()
        for stat, key in module_step_matrix_stats.items():
            self.write_stat_to_matrix(stat_name=stat, stat_header='module_step', stat_key=key, stat_dict=steps_superdict_multi, \
                                      item_list=module_step_list, stat_metadata_headers=STEP_METADATA_HEADERS, \
                                      write_rows_with_all_zeros=include_zeros)

        # now we make a KO hit count matrix
        ko_list = list(self.ko_dict.keys())
        ko_list.sort()
        self.write_stat_to_matrix(stat_name='enzyme_hits', stat_header='enzyme', stat_key='num_hits', stat_dict=ko_superdict_multi, \
                                  item_list=ko_list, stat_metadata_headers=KO_METADATA_HEADERS, \
                                  write_rows_with_all_zeros=include_zeros)

        # if necessary, make module specific KO matrices
        if self.module_specific_matrices:
            skipped_mods = []
            mods_defined_by_mods = []
            for mod in self.module_specific_matrices:
                if mod not in module_list:
                    skipped_mods.append(mod)
                    continue

                mod_def = self.all_modules_in_db[mod]['DEFINITION']
                if isinstance(mod_def, list):
                    mod_def = " ".join(mod_def)
                kos_in_mod = self.get_enzymes_from_module_definition_in_order(mod_def, self.all_modules_in_db)
                mod_big_steps = utils.split_by_delim_not_within_parens(mod_def, " ")

                if not kos_in_mod:
                    mods_defined_by_mods.append(mod)
                    continue

                if not self.no_comments:
                    # determine where to place comments containing module steps
                    step_comments = {}
                    lines_with_comment = []
                    for s in mod_big_steps:
                        split_s = re.sub("[\(\)\+\-,]", " ", s).strip()
                        split_s = split_s.split(" ")

                        # we skip making comments on steps without KOs like '--'
                        if not split_s:
                            continue
                        # what is the first KO in this step?
                        first_ko = split_s[0]

                        # figure out where this KO is in the list (it could be there multiple times)
                        first_ko_indices = [i for i, x in enumerate(kos_in_mod) if x == first_ko]

                        # if there are multiple instances of this KO, where should we put the step comment?
                        idx = first_ko_indices[0]
                        if len(first_ko_indices) > 1:
                            next_index = 0
                            while idx in lines_with_comment:
                                next_index += 1
                                idx = first_ko_indices[next_index]

                        step_comments[idx] = s
                        lines_with_comment.append(idx)

                else:
                    step_comments = None

                stat = f"{mod}_enzyme_hits"
                self.write_stat_to_matrix(stat_name=stat, stat_header="enzyme", stat_key='num_hits', stat_dict=ko_superdict_multi, \
                                          item_list=kos_in_mod, stat_metadata_headers=KO_METADATA_HEADERS, \
                                          write_rows_with_all_zeros=True, comment_dictionary=step_comments)

            if skipped_mods:
                skipped_list = ", ".join(skipped_mods)
                self.run.warning(f"We couldn't recognize the following module(s): {skipped_list}. So we didn't generate "
                                 "output matrices for them. If you used the `--only-complete` flag, its possible these modules "
                                 "were eliminated from the output due to having completeness scores below the threshold (in "
                                 "which case you could just remove `--only-complete` from your command and everything should "
                                 "work fine). Otherwise, maybe you made a typo? Or put an extra comma in somewhere?")

            if mods_defined_by_mods:
                skipped_list = ", ".join(mods_defined_by_mods)
                self.run.warning(f"The following modules were completely defined by other modules and therefore we didn't "
                                 f"generate enzyme hit matrices for them: {skipped_list}.")


## STATIC FUNCTIONS
def module_definition_to_enzyme_accessions(mod_definition):
    """Parses a module definition string into a list of enzyme accessions."""

    # anything that is not (),-+ should be converted to spaces, then we can split on the spaces to get the accessions
    mod_definition = re.sub('[\(\)\+\-,]', ' ', mod_definition).strip()
    acc_list = re.split(r'\s+', mod_definition)

    return acc_list
