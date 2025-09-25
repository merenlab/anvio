#!/usr/bin/env python
# -*- coding: utf-8
"""Database access and data loading operations for KEGG metabolism estimation."""

import pandas as pd

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.dbinfo import DBInfo
from anvio.errors import ConfigError

from anvio.metabolism.context import KeggContext
from anvio.metabolism.modulesdb import ModulesDatabase
from anvio.metabolism.constants import DEFAULT_OUTPUT_MODE, OUTPUT_MODES, OUTPUT_HEADERS, STRAY_KO_ANVIO_SUFFIX


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Iva Veseli"
__email__ = "iveseli@uchicago.edu"


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print
P = terminal.pluralize


class KeggEstimatorArgs():
    def __init__(self, args, format_args_for_single_estimator=False, run=terminal.Run(), progress=terminal.Progress()):
        """A base class to assign arguments for KeggMetabolism estimator classes.

        Parameters
        ==========
        format_args_for_single_estimator: bool
            This is a special case where an args instance is generated to be passed to the
            single estimator from within multi estimator. More specifically, the multi estimator
            class is nothing but one that iterates through all contigs DBs
            given to it using the single estimator class. So it needs to create instances of
            single estimators, and collect results at upstream. The problem is, if a single
            estimator is initiated with the args of a multi estimator, the sanity check will
            go haywire. This flag nullifies most common offenders.
        """

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.metagenome_mode = True if A('metagenome_mode') else False
        self.module_completion_threshold = A('module_completion_threshold') or 0.75
        self.output_file_prefix = A('output_file_prefix') or "metabolism"
        self.write_dict_to_json = True if A('get_raw_data_as_json') else False
        self.json_output_file_path = A('get_raw_data_as_json')
        self.store_json_without_estimation = True if A('store_json_without_estimation') else False
        self.estimate_from_json = A('estimate_from_json') or None
        self.enzymes_txt = A('enzymes_txt') or None
        self.output_modes = A('output_modes')
        self.custom_output_headers = A('custom_output_headers') or None
        self.matrix_format = True if A('matrix_format') else False
        self.matrix_include_metadata = True if A('include_metadata') else False
        self.exclude_zero_modules = False if A('include_zeros') else True
        self.only_complete = True if A('only_complete') else False
        self.add_coverage = True if A('add_coverage') else False
        self.add_copy_number = True if A('add_copy_number') else False
        self.exclude_kos_no_threshold = False if A('include_kos_not_in_kofam') else True
        self.include_stray_kos = True if A('include_stray_KOs') else False
        self.ignore_unknown_kos = True if A('ignore_unknown_KOs') else False
        self.exclude_dashed_reactions = True if A('exclude_dashed_reactions') else False
        self.module_specific_matrices = A('module_specific_matrices') or None
        self.no_comments = True if A('no_comments') else False
        self.external_genomes_file = A('external_genomes') or None
        self.internal_genomes_file = A('internal_genomes') or None
        self.metagenomes_file = A('metagenomes') or None
        self.kegg_data_dir = A('kegg_data_dir')

        self.modules_unique_id = None
        self.ko_unique_id = None
        self.genome_mode = False  ## controls some warnings output, will be set to True downstream if necessary

        # A bit of special attention for self.enzymes_of_interest_df. the purpose of this variable is to
        # give programmers a means to initialize estimator classes with a data frame of enzymes. so this
        # variable is not accessible to any of the command line interfaces, but if the args object includes
        # an enzymes_txt (which is accessible through command line interfaces), the contents of
        # self.enzymes_of_interest_df will be automatically filled by this class below. an important
        # point is that all classes that inherit KeggEstimatorArgs will only  use `self.enzymes_of_interest_df`
        # for all downstream tasks.
        self.enzymes_of_interest_df = A('enzymes_of_interest_df')

        # the below will be filled in by init_data_from_modules_db()
        self.all_modules_in_db = {}
        self.all_kos_in_db = {}
        self.module_paths_dict = {}

        # if necessary, assign 0 completion threshold, which evaluates to False above
        if A('module_completion_threshold') == 0:
            self.module_completion_threshold = 0.0

        # we use the below flag to find out if long format output was explicitly requested
        self.long_format_mode = True if self.output_modes else False
        # if it was not explicitly requested, we set the default output mode to "modules"
        if not self.output_modes:
            self.output_modes = DEFAULT_OUTPUT_MODE

        # output modes and headers that we can handle
        self.available_modes = OUTPUT_MODES
        self.available_headers = OUTPUT_HEADERS

        if format_args_for_single_estimator:
            # to fool a single estimator into passing sanity checks, nullify multi estimator args here
            self.databases = None
            self.matrix_format = False # we won't be storing data from the single estimator anyway
            self.module_specific_matrices = None

        # parse requested output modes if necessary
        if isinstance(self.output_modes, str):
            # parse requested output modes and make sure we can handle them all
            self.output_modes = self.output_modes.split(",")

        # parse requested output headers if necessary
        if self.custom_output_headers and isinstance(self.custom_output_headers, str):
            self.custom_output_headers = self.custom_output_headers.split(",")
            self.available_modes['modules_custom']['headers'] = self.custom_output_headers

        # parse specific matrix modules if necessary
        if self.module_specific_matrices:
            self.module_specific_matrices = [_m.strip() for _m in self.module_specific_matrices.split(",")]

        # load up the enzymes of interest if the user passed an enzyme_txt file
        if self.enzymes_txt:
            self.enzymes_of_interest_df = self.get_enzymes_of_interest_df()


    def setup_output_for_appending(self):
        """Initializes and returns a dictionary of AppendableFile objects, one for each output mode"""

        output_dict = {}
        for mode in self.output_modes:
            output_path = self.output_file_prefix + "_" + self.available_modes[mode]["output_suffix"]

            # check the status of this output file
            filesnpaths.is_output_file_writable(output_path, ok_if_exists=False)

            output_file_for_mode = filesnpaths.AppendableFile(output_path, append_type=dict, fail_if_file_exists=False)
            output_dict[mode] = output_file_for_mode

            self.run.info(f"Output file for {mode} mode", output_path)

        return output_dict


    def init_data_from_modules_db(self):
        """This function reads mucho data from the modules database(s) into dictionaries for later access.

        It generates the self.all_modules_in_db dictionary, which contains all data for all modules
        in the db, keyed by module number.
        It also generates the self.all_kos_in_db dictionary, which maps each enzyme in the db to its list of modules.
        Note that self.all_kos_in_db can contain module numbers, in cases when the module is a component of another module.
        It also generates the self.module_paths_dict dictionary by calling the appropriate function.

        We do this once at the start so as to reduce the number of on-the-fly database queries
        that have to happen during the estimation process.

        ## KEYS ADDED TO SELF.ALL_MODULES_IN_DB DICTIONARY
        * all keys from modules table in database, ie DEFINITION, CLASS, etc
        'MODULES_DB_SOURCE'         which database this module belongs to ("KEGG" for KEGG modules, "USER" for user-defined modules)
        'substrate_list'            list of substrate compounds (inputs to the module)
        'intermediate_list'         list of intermediate compounds
        'product_list'              list of product compounds (outputs of the module)
        'top_level_steps'           list of top-level steps in the module DEFINITION
        """

        # LOAD KEGG DATA (MODULES)
        if not self.only_user_modules:
            self.kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, args=self.args, run=run_quiet, quiet=self.quiet)
            self.all_modules_in_db = self.kegg_modules_db.get_modules_table_as_dict()
            # mark that these modules all came from KEGG
            for mod in self.all_modules_in_db:
                self.all_modules_in_db[mod]['MODULES_DB_SOURCE'] = "KEGG"

                # add compound lists into self.all_modules_in_db
                module_substrate_list, module_intermediate_list, module_product_list = self.kegg_modules_db.get_human_readable_compound_lists_for_module(mod)
                self.all_modules_in_db[mod]['substrate_list'] = module_substrate_list
                self.all_modules_in_db[mod]['intermediate_list'] = module_intermediate_list
                self.all_modules_in_db[mod]['product_list'] = module_product_list

                # initialize module paths into self.module_paths_dict
                self.module_paths_dict[mod] = self.init_paths_for_module(mod, mod_db=self.kegg_modules_db)

                # initialize top-level steps into self.all_modules_in_db
                self.all_modules_in_db[mod]['top_level_steps'] = self.kegg_modules_db.get_top_level_steps_in_module_definition(mod)

            self.kegg_modules_db.disconnect()

        # LOAD USER DATA (MODULES)
        if self.user_input_dir:
            self.user_modules_db = ModulesDatabase(self.user_modules_db_path, args=self.args, run=run_quiet, quiet=self.quiet)
            user_db_mods = self.user_modules_db.get_modules_table_as_dict()

            for mod in user_db_mods:
                # user modules cannot have the same name as a KEGG module
                if mod in self.all_modules_in_db:
                    self.user_modules_db.disconnect()
                    raise ConfigError(f"No. Nononono. Stop right there. You see, there is a module called {mod} in your user-defined "
                                      f"modules database (at {self.user_modules_db_path}) which has the same name as an existing KEGG "
                                      f"module. This is not allowed, for reasons. Please name that module differently. Append an "
                                      f"underscore and your best friend's name to it or something. Just make sure it's unique. OK? ok.")

                self.all_modules_in_db[mod] = user_db_mods[mod]
                # mark that this module came from the USER. it may be useful to know later.
                self.all_modules_in_db[mod]['MODULES_DB_SOURCE'] = "USER"

                # add compound lists into self.all_modules_in_db
                module_substrate_list, module_intermediate_list, module_product_list = self.user_modules_db.get_human_readable_compound_lists_for_module(mod)
                self.all_modules_in_db[mod]['substrate_list'] = module_substrate_list
                self.all_modules_in_db[mod]['intermediate_list'] = module_intermediate_list
                self.all_modules_in_db[mod]['product_list'] = module_product_list

                # initialize module paths into self.module_paths_dict
                self.module_paths_dict[mod] = self.init_paths_for_module(mod, mod_db=self.user_modules_db)

                # initialize top-level steps into self.all_modules_in_db
                self.all_modules_in_db[mod]['top_level_steps'] = self.user_modules_db.get_top_level_steps_in_module_definition(mod)

            self.user_modules_db.disconnect()

        # INIT ENZYMES
        for mod in self.all_modules_in_db:
            orthology = self.all_modules_in_db[mod]['ORTHOLOGY']
            if isinstance(orthology, str):
                ko_list = [orthology]
            else:
                ko_list = list(orthology.keys())
            for k in ko_list:
                if k not in self.all_kos_in_db:
                    src = self.all_modules_in_db[mod]['ANNOTATION_SOURCE'][k] if 'ANNOTATION_SOURCE' in self.all_modules_in_db[mod] else 'KOfam'
                    func = self.all_modules_in_db[mod]['ORTHOLOGY'][k] if 'ORTHOLOGY' in self.all_modules_in_db[mod] else 'UNKNOWN'
                    self.all_kos_in_db[k] = {'modules': [], 'annotation_source': src, 'function': func}
                self.all_kos_in_db[k]['modules'].append(mod)


    def init_paths_for_module(self, mnum, mod_db=None):
        """This function unrolls the module DEFINITION for the module provided and returns a list of all paths through it.

        It unrolls the module definition into a list of all possible paths, where each path is a list of atomic steps.
        Atomic steps include singular KOs, protein complexes, modules, non-essential steps, and steps without associated KOs.

        PARAMETERS
        ==========
        mnum : str
            The module to return paths for. Must be a key in the self.all_modules_in_db dictionary.
        mod_db : ModulesDatabase
            This must be a ModulesDatabase instance that we are connected to (ie, disconnect() has not yet been run on it)
            so that we can use its functions and access its data

        RETURNS
        ==========
        A list of all paths through the module
        """

        if not mod_db:
            raise ConfigError("Put yer hands in the air! You've tried to call init_paths_for_modules() without providing "
                              "a database to the mod_db parameter, and this is ILLEGAL.")

        if mnum not in self.all_modules_in_db:
            raise ConfigError(f"Something is wrong here. The function init_paths_for_modules() is trying to work on module "
                              f"{mnum}, but it is not a key in the self.all_modules_in_db dictionary.")
        module_definition = self.all_modules_in_db[mnum]["DEFINITION"]
        # the below function expects a list
        if not isinstance(module_definition, list):
            module_definition = [module_definition]

        return mod_db.unroll_module_definition(mnum, def_lines=module_definition)


    def get_module_metadata_dictionary(self, mnum):
        """Returns a dictionary of metadata for the given module.

        The dictionary must include all the metadata from MODULE_METADATA_HEADERS,
        using those headers as keys.

        Requires self.all_modules_in_db attribute to exist - subclasses will have to call init_data_from_modules_db()
        before this function.
        """

        if not self.all_modules_in_db:
            raise ConfigError("The function get_module_metadata_dictionary() requires the self.all_modules_in_db attribute to "
                              "be initialized. You need to make sure init_data_from_modules_db() is called before this function. ")

        class_data_val = self.all_modules_in_db[mnum]['CLASS']
        fields = class_data_val.split("; ")
        mnum_class_dict = {"class" : fields[0], "category" : fields[1], "subcategory" : fields[2] if len(fields) > 2 else None}

        metadata_dict = {}
        metadata_dict["module_name"] = self.all_modules_in_db[mnum]['NAME']
        metadata_dict["module_class"] = mnum_class_dict["class"]
        metadata_dict["module_category"] = mnum_class_dict["category"]
        metadata_dict["module_subcategory"] = mnum_class_dict["subcategory"]
        return metadata_dict


    def get_ko_metadata_dictionary(self, knum, dont_fail_if_not_found=False):
        """Returns a dictionary of metadata for the given KO.

        The dictionary must include all the metadata from KO_METADATA_HEADERS,
        using those headers as keys.

        Requires self.all_kos_in_db attribute to exist - subclasses will have to call init_data_from_modules_db()
        before this function.
        """

        if not self.all_kos_in_db:
            raise ConfigError("The function get_ko_metadata_dictionary() requires the self.all_kos_in_db attribute to "
                              "be initialized. You need to make sure init_data_from_modules_db() is called before this function. ")

        mod_list = self.all_kos_in_db[knum]['modules'] if knum in self.all_kos_in_db else None
        if mod_list:
            mod_list_str = ",".join(mod_list)
        else:
            mod_list_str = "None"

        metadata_dict = {}
        metadata_dict["modules_with_enzyme"] = mod_list_str

        if knum in self.ko_dict:
            metadata_dict["enzyme_definition"] = self.ko_dict[knum]['definition']
        elif self.include_stray_kos and self.stray_ko_dict and knum in self.stray_ko_dict or f"{knum}{STRAY_KO_ANVIO_SUFFIX}" in self.stray_ko_dict:
            # if we can't find the enzyme in the KO dictionary, try to find it in the stray KO dictionary (if it exists)
            if knum not in self.stray_ko_dict and f"{knum}{STRAY_KO_ANVIO_SUFFIX}" in self.stray_ko_dict:
                knum = f"{knum}{STRAY_KO_ANVIO_SUFFIX}"
            metadata_dict["enzyme_definition"] = self.stray_ko_dict[knum]['definition']
        else:
            # if we still can't find the enzyme, try to find it in the database
            if knum in self.all_kos_in_db and 'function' in self.all_kos_in_db[knum]:
                metadata_dict["enzyme_definition"] = self.all_kos_in_db[knum]['function']
            elif dont_fail_if_not_found:
                self.run.warning(f"The enzyme {knum} was not found in the metabolism data, so we are unable to determine "
                                 f"its functional annotation. You will see 'UNKNOWN' for this enzyme in any outputs describing "
                                 f"its function.")
                metadata_dict["enzyme_definition"] = "UNKNOWN"
            else:
                raise ConfigError("Something is mysteriously wrong. You are seeking metadata "
                                  f"for enzyme {knum} but this enzyme is not in the enzyme dictionary "
                                  "(self.ko_dict, or (self.stray_ko_dict) in some cases). This should never have happened.")

        return metadata_dict


    def get_enzymes_of_interest_df(self):
        """This function loads and sanity checks an enzymes txt file, and returns it as a pandas dataframe.

        RETURNS
        =======
        enzyme_df : Pandas DataFrame
            contains the information in the enzymes txt file
        """

        self.progress.new("Loading enzymes-txt file...")
        expected_fields = ['gene_id', 'enzyme_accession', 'source']
        enzyme_df = pd.read_csv(self.enzymes_txt, sep="\t")
        self.progress.end()

        self.run.info("Number of genes loaded from enzymes-txt file", enzyme_df.shape[0])

        # sanity check for required columns
        missing = []
        for f in expected_fields:
            if f not in enzyme_df.columns:
                missing.append(f)
        if missing:
            miss_str = ", ".join(missing)
            exp_str = ", ".join(expected_fields)
            raise ConfigError(f"Your enzymes-txt file ({self.enzymes_txt}) is missing some required columns. "
                              f"The columns it needs to have are: {exp_str}. And the missing column(s) include: {miss_str}")

        # warning about extra columns
        used_cols = expected_fields + ['coverage', 'detection']
        extra_cols = []
        for c in enzyme_df.columns:
            if c not in used_cols:
                extra_cols.append(c)
        if extra_cols:
            e_str = ", ".join(extra_cols)
            self.run.warning(f"Just so you know, your input enzymes-txt file contained some columns of data that we are not "
                             f"going to use. This isn't an issue or anything, just an FYI. We're ignoring the following field(s): {e_str}")

        # if cov/det columns are not in the file, we explicitly turn off flag to add this data to output
        if self.add_coverage and ('coverage' not in enzyme_df.columns or 'detection' not in enzyme_df.columns):
            self.run.warning("You requested coverage/detection values to be added to the output files, but your "
                             "input file does not seem to contain either a 'coverage' column or a 'detection' column, or both. "
                             "Since we don't have this data, --add-coverage will not work, so we are turning this "
                             "flag off. Sorry ¯\_(ツ)_/¯")
            self.add_coverage = False
            # remove coverage headers from the list so we don't try to access them later
            kofam_hits_coverage_headers = [self.contigs_db_project_name + "_coverage", self.contigs_db_project_name + "_detection"]
            modules_coverage_headers = [self.contigs_db_project_name + "_gene_coverages", self.contigs_db_project_name + "_avg_coverage",
                                        self.contigs_db_project_name + "_gene_detection", self.contigs_db_project_name + "_avg_detection"]
            for h in kofam_hits_coverage_headers:
                if h in self.available_modes["hits"]["headers"]:
                    self.available_modes["hits"]["headers"].remove(h)
            for h in modules_coverage_headers:
                if h in self.available_modes["modules"]["headers"]:
                    self.available_modes["modules"]["headers"].remove(h)

        return enzyme_df


class KeggDataLoader(KeggContext):
    """Handles data loading operations for KEGG metabolism estimation."""

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress
        KeggContext.__init__(self, args)

    def init_hits_and_splits(self, contigs_db_path, profile_db_path=None, collection_name=None,
                            annotation_sources=['KOfam'], splits_to_use=None,
                            user_input_dir=None, all_kos_in_db=None):
        """This function loads KOfam hits, gene calls, splits, and contigs from the contigs DB.

        We will need the hits with their KO numbers (accessions) so that we can go through the MODULES.db and determine
        which steps are present in each module. And we will need the other information so that we can determine which hits belong
        to which genomes/bins when we are handling multiple of these, and for help in computing redundancy.
        This function gets this info as a list of tuples (one tuple per kofam hit), and it makes sure that these lists don't include
        hits that we shouldn't be considering.

        PARAMETERS
        ==========
        annotation_sources : list
            which functional annotation sources to obtain gene calls from. Should at least contain 'Kofam' for
            default usage. Adding other sources may be necessary when working with user-defined metabolic modules.

        RETURNS
        =======
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering
        """

        self.progress.new("Loading split data from contigs DB")
        split_names_in_contigs_db = set(utils.get_all_item_names_from_the_database(contigs_db_path))
        if splits_to_use is None:
            splits_to_use = split_names_in_contigs_db

        # first, resolve differences in splits between profile and contigs db
        if profile_db_path:
            self.progress.update("Loading split data from profile DB")
            # if we were given a blank profile, we will assume we want all splits and pull all splits from the contigs DB
            if utils.is_blank_profile(profile_db_path):
                self.progress.reset()
                self.run.warning("You seem to have provided a blank profile. No worries, we can still estimate metabolism "
                                 "for you. But we cannot load splits from the profile DB, so instead we are assuming that "
                                 "you are interested in ALL splits and we will load those from the contigs database.")
            else:
                split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(profile_db_path))
                splits_missing_in_profile_db = split_names_in_contigs_db.difference(split_names_in_profile_db)

                if len(splits_missing_in_profile_db):
                    min_contig_length_in_profile_db = pp(int(DBInfo('PROFILE.db', expecting='profile').get_self_table()['min_contig_length']))
                    num_splits_contig = pp(len(split_names_in_contigs_db))
                    num_splits_profile = pp(len(split_names_in_profile_db))
                    num_missing = pp(len(splits_missing_in_profile_db))
                    self.progress.reset()
                    self.run.warning(f"Please note that anvi'o found {num_splits_contig} splits in your contigs database. "
                                     f"But only {num_splits_profile} of them appear in the profile database. As a result, "
                                     f"anvi'o will now remove the {num_missing} splits that occur only in the contigs db "
                                     f"from all downstream analyses. Where is this difference coming from though? Well. This "
                                     f"is often the case because the 'minimum contig length parameter' set during the `anvi-profile` "
                                     f"step can exclude many contigs from downstream analyses (often for good reasons, too). For "
                                     f"instance, in your case the minimum contig length set in the profile database is "
                                     f"{min_contig_length_in_profile_db} nts. Anvi'o hopes that this explains some things.")
                    splits_to_use = split_names_in_profile_db

            # second, if we are working with a collection, we can limit the splits to use with those from the collection
            if collection_name:
                splits_to_use = ccollections.GetSplitNamesInBins(self.args).get_split_names_only()
                self.progress.reset()
                self.run.warning(f"Since a collection name was provided, we will only work with gene calls "
                                 f"from the subset of {len(splits_to_use)} splits in the collection for the "
                                 f"purposes of estimating metabolism.")

        self.progress.update('Loading gene call data from contigs DB')
        from anvio.dbops import ContigsDatabase # <- import here to avoid circular import
        contigs_db = ContigsDatabase(contigs_db_path, run=self.run, progress=self.progress)

        split_list = ','.join(["'%s'" % split_name for split_name in splits_to_use])
        splits_where_clause = f'''split IN ({split_list})'''
        genes_in_splits = contigs_db.db.get_some_columns_from_table(t.genes_in_splits_table_name, "gene_callers_id, split",
                                                                    where_clause=splits_where_clause)

        gene_list = ','.join(["'%s'" % gcid for gcid,split in genes_in_splits])
        contigs_where_clause = f'''gene_callers_id IN ({gene_list})'''
        genes_in_contigs = contigs_db.db.get_some_columns_from_table(t.genes_in_contigs_table_name, "gene_callers_id, contig",
                                                                     where_clause=contigs_where_clause)

        source_list = ','.join(["'%s'" % src for src in annotation_sources])
        hits_where_clause = f'''source IN ({source_list}) AND gene_callers_id IN ({gene_list})'''
        kofam_hits = contigs_db.db.get_some_columns_from_table(t.gene_function_calls_table_name, "gene_callers_id, accession, function",
                                                               where_clause=hits_where_clause)

        contigs_db.disconnect()

        # combine the information for each gene call into neat tuples for returning
        # each gene call is only on one split of one contig, so we can convert these lists of tuples into dictionaries for easy access
        # but some gene calls have multiple kofam hits (and some kofams have multiple gene calls), so we must keep the tuple structure for those
        self.progress.update("Organizing gene call data")
        gene_calls_splits_dict = {tpl[0] : tpl[1] for tpl in genes_in_splits}
        gene_calls_contigs_dict = {tpl[0] : tpl[1] for tpl in genes_in_contigs}
        assert len(gene_calls_splits_dict.keys()) == len(genes_in_contigs)

        kofam_gene_split_contig = []
        for gene_call_id, ko, func in kofam_hits:
            # some genes have multiple annotations that we need to split
            for annotation in ko.split('!!!'):
                kofam_gene_split_contig.append((annotation, gene_call_id, gene_calls_splits_dict[gene_call_id], gene_calls_contigs_dict[gene_call_id]))

                # for user data, the enzymes in some loaded gene calls may not yet be in the ko dict, so we add them in here.
                if user_input_dir:
                    if not hasattr(self, 'ko_dict') or not self.ko_dict:
                        raise ConfigError("Uh oh. The code is currently trying to add gene annotations to self.ko_dict, but this attribute does not "
                                          "exist! You need to fix this.")
                    if annotation not in self.ko_dict:
                        self.ko_dict[annotation] = {'definition': func}

        self.progress.update("Done")
        self.progress.end()

        sources_str = ", ".join(annotation_sources)
        self.run.info("Annotation sources used", sources_str)
        self.run.info("Gene calls from these sources", "%d found" % len(kofam_hits), quiet=self.quiet)

        if not self.quiet and not len(kofam_hits):
            self.run.warning(f"Hmmm. No gene calls from any of the following annotation sources were found in this contigs DB: {sources_str}. "
                             f"The consequence is that all metabolism estimate outputs will be empty. This is fine, and could even be biologically "
                             f"correct. But we thought we'd mention it just in case you thought it was weird. Other, technical reasons that this could "
                             f"have happened include: 1) you didn't annotate with `anvi-run-kegg-kofams` or another annotation program, or "
                             "2) you imported functional annotations but the 'source' did not match those in the list above.")

        return kofam_gene_split_contig


    def init_gene_coverage(self, profile_db, gcids_for_kofam_hits):
        """This function initializes gene coverage/detection values from the provided profile DB.

        The profile DB should be already initialized for this to work (currently add_gene_coverage_to_headers_list()
        handles this). The reason we split the initalization of the profile db from the initialization of gene
        coverage/detection values is so that we only work on the set of gene calls with KOfam hits rather than all
        genes in the contigs DB.

        PARAMETERS
        ==========
        gcids_for_kofam_hits : set
            The gene caller ids for all genes with KOfam hits in the contigs DB
        """

        if not profile_db:
            raise ConfigError("A profile DB has not yet been initialized, so init_gene_coverage() will not work. "
                              "If you are a programmer, you should probably either 1) call this function after "
                              "add_gene_coverage_to_headers_list() or 2) extend this function so that it initializes "
                              "the profile db. If you are not a programmer, you should probably find one :) ")
        self.run.info_single("Since the --add-coverage flag was provided, we are now loading the relevant "
                             "coverage information from the provided profile database.")
        profile_db.init_gene_level_coverage_stats_dicts(gene_caller_ids_of_interest=gcids_for_kofam_hits)


    def init_hits_for_pangenome(self, pan_db_path, genomes_storage_path, gene_cluster_list, annotation_sources_to_use):
        """This function loads enzyme annotations from the pangenome for use by downstream metabolism estimation.

        For each gene cluster, it takes the most common function from each annotation source relevant to the modules.

        PARAMETERS
        ==========
        gene_cluster_list : list
            which gene cluster IDs to load from the pan DB

        RETURNS
        =======
        enzyme_cluster_split_contig : list
            (enzyme_accession, gene_cluster_id, split, contig) tuples in which split and contig are both NAs
        """

        from anvio.dbops import PanSuperclass # <- import here to avoid circular import

        pan_super = PanSuperclass(self.args)
        pan_super.init_gene_clusters(gene_cluster_ids_to_focus = gene_cluster_list)
        pan_super.init_gene_clusters_functions_summary_dict(source_list = annotation_sources_to_use, gene_clusters_of_interest = gene_cluster_list)

        enzyme_cluster_split_contig = []
        # no splits or contigs here
        for cluster_id in gene_cluster_list:
            for source in annotation_sources_to_use:
                if source in pan_super.gene_clusters_functions_summary_dict[cluster_id]:
                    acc = pan_super.gene_clusters_functions_summary_dict[cluster_id][source]['accession']
                    if acc: # avoid introducing 'None' values here
                        enzyme_cluster_split_contig.append((acc,cluster_id,"NA","NA"))

        return enzyme_cluster_split_contig
