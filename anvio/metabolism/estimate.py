#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains Kegg related classes."""

import os
import re
import copy
import json
import statistics

import pandas as pd
from scipy import stats

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.dbinfo import DBInfo
from anvio.errors import ConfigError
from anvio.genomedescriptions import MetagenomeDescriptions, GenomeDescriptions

from anvio.metabolism.context import KeggContext
from anvio.metabolism.modulesdb import ModulesDatabase
from anvio.metabolism.constants import DEFAULT_OUTPUT_MODE, OUTPUT_MODES, OUTPUT_HEADERS, STRAY_KO_ANVIO_SUFFIX, STEP_METADATA_HEADERS, KO_METADATA_HEADERS, MODULE_METADATA_HEADERS

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


    def split_module_path_into_individual_essential_components(self, path):
        """Given a list of atomic steps in a module, this function returns a list of each essential individual enzyme.

        When an atomic step is an enzyme complex (ie K01657+K01658), we need to split the complex into its individual components
        When an atomic step contains non-ssential components (ie K00765-K02502), we need to remove the nonessential components from the list
        When there are both nonessential and essential components, we need to remove the non-essential ones first and then split the essential ones

        PARAMETERS
        ==========
        path : list
            Each element in list is an atomic step, which can include enzyme complexes and non-essential components

        RETURNS
        ==========
        new_path : list
            Each element in list is a single enzyme
        """

        new_path = []
        for x in path:
            if '+' and '-' in x:
                # first remove the nonessentials, then add in the essential components
                idx_of_nonessential = x.index('-')
                new_x = x[:idx_of_nonessential]
                individual_components = new_x.split('+')
                new_path.extend(individual_components)
            elif '+' in x:
                # split essential components and add individually to list
                individual_components = x.split('+')
                new_path.extend(individual_components)
            elif '-' in x and x != '--':
                # remove nonessentials
                idx_of_nonessential = x.index('-')
                new_x = x[:idx_of_nonessential]
                new_path.append(new_x)
            else:
                new_path.append(x)

        return new_path


    def get_enzymes_from_module_definition_in_order(self, mod_definition):
        """Given a module DEFINITION string, this function parses out the enzyme accessions in order of appearance.

        PARAMETERS
        ==========
        mod_definition : a string or list of strings containing the module DEFINITION lines

        RETURNS
        ==========

        """

        if isinstance(mod_definition, list):
            mod_definition = " ".join(mod_definition)

        acc_list = module_definition_to_enzyme_accessions(mod_definition)

        # remove anything that is not an enzyme and sanity check for weird characters
        mods_to_remove = set()
        for a in acc_list:
            if a in self.all_modules_in_db:
                mods_to_remove.add(a)
            if re.match('[^a-zA-Z0-9_\.]', a):
                raise ConfigError(f"The get_enzymes_from_module_definition_in_order() function found an enzyme accession that looks a bit funny. "
                                  f"Possibly this is a failure of our parsing strategy, or maybe the enzyme accession just has unexpected characters "
                                  f"in it. We don't know what module it is, but the weird enzyme is {a}. If you think that accession looks perfectly "
                                  f"fine, you should reach out to the developers and have them fix this function to accomodate the accession. Or, you "
                                  f"could just rename the enzyme?")
        if mods_to_remove:
            for m in mods_to_remove:
                acc_list.remove(m)

        return acc_list


    def remove_nonessential_enzymes_from_module_step(self, step_string):
        """This functions removes nonessential enzymes from a step definition string and returns the resulting definition.

        It is intended to be called on top-level steps of a module definition (not on the full module definition).

        A nonessential enzyme is any accession (ie, '-K11024') or group of accessions (ie, -(K00242,K18859,K18860))
        that is marked with a leading '-' character. This function finds the '-' characters in the given string and
        removes any subsequent accessions. The resulting definition with only essential components is returned.

        If a step does not contain nonessential enzymes, the original definition is returned.

        Likewise, '--' is a special case containing the '-' character which actually indicates a step that has no enzyme profile.
        It seems to always be present as its own step (ie, '--' is the entire step definition string), so we return the original
        definition in this case. However, in case there is an internal '--' within a more complicated definition, this function
        ignores the part of the string that includes it and processes the remainder of the string before re-joining the two parts.
        It is not able to do this for steps with more than one internal '--', which would require multiple splits and joins, so
        this case results in an error. Note that when self.exclude_dashed_reactions is True, we instead remove '--' entirely.

        PARAMETERS
        ==========
        step_string: str
            A string containing the definition of one step from a module

        RETURNS
        =======
        step_string: str
            The same string, with nonessential enzyme accessions (if any) removed.
        """

        if step_string == '--' and self.exclude_dashed_reactions:
            return ""
        if step_string != '--' and '-' in step_string:
            saw_double_dash = False             # a Boolean to indicate if we found '--' within the step definition
            str_prior_to_double_dash = None     # if we find '--', this variable stores the string that occurs prior to and including this '--'
            while '-' in step_string:
                idx = step_string.index('-')
                if step_string[idx+1] == '-': # internal '--' case
                    if saw_double_dash:
                        raise ConfigError("Unfortunately, a step containing >1 internal instances of '--' was passed to the function "
                                          "remove_nonessential_enzymes_from_module_step(). This function is not currently able to handle this "
                                          "situation. Please contact a developer and ask them to turn this into a smarter function. :) ")
                    saw_double_dash = True
                    if self.exclude_dashed_reactions: # remove the internal '--'
                        str_prior_to_double_dash = step_string[:idx]
                        step_string = step_string[idx+3:] # also remove the space after it
                    else:
                        str_prior_to_double_dash = step_string[:idx+2]
                        step_string = step_string[idx+2:] # continue processing the remainder of the string
                    continue
                elif step_string[idx+1] == '(': # group to eliminate
                    parens_index = idx+1
                    while step_string[parens_index] != ')':
                        parens_index += 1
                    step_string = step_string[:idx] + step_string[parens_index+1:]
                else: # single non-essential enzyme
                    punctuation_index = idx+1
                    while punctuation_index < len(step_string) and step_string[punctuation_index] not in [')','(','+',' ',',']:
                        punctuation_index += 1
                    step_string = step_string[:idx] + step_string[punctuation_index:]

            # if we found an internal '--', we put the two pieces of the definition back together
            if saw_double_dash:
                step_string = str_prior_to_double_dash + step_string

        return step_string.strip()


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

        # check and warning for enzymes not in self.all_kos_in_db
        enzymes_not_in_modules = list(enzyme_df[~enzyme_df["enzyme_accession"].isin(self.all_kos_in_db.keys())]['enzyme_accession'].unique())
        if self.include_stray_kos:
            enzymes_not_in_modules = [e for e in enzymes_not_in_modules if e not in self.stray_ko_dict]
        if enzymes_not_in_modules:
            example = enzymes_not_in_modules[0]
            self.run.warning(f"FYI, some enzymes in the 'enzyme_accession' column of your input enzymes-txt file do not belong to any "
                             f"metabolic modules (that we know about). These enzymes will be ignored for the purposes of estimating module "
                             f"completeness, but should still appear in enzyme-related outputs (if those were requested). In case you are "
                             f"curious, here is one example: {example}")

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


class KeggMetabolismEstimator(KeggContext, KeggEstimatorArgs):
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
        KeggContext.__init__(self, self.args)

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
            raise ConfigError("You have requested some output modes that we cannot handle. The offending modes "
                              "are: %s. Please use the flag --list-available-modes to see which ones are acceptable."
                              % (", ".join(illegal_modes)))
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
                raise ConfigError("You have requested some output headers that we cannot handle. The offending ones "
                                  "are: %s. Please use the flag --list-available-output-headers to see which ones are acceptable."
                                  % (", ".join(illegal_headers)))

            # check if any headers requested for modules_custom mode are reserved for KOfams mode
            if "modules_custom" in self.output_modes:
                for header in self.custom_output_headers:
                    if self.available_headers[header]['mode_type'] != "modules" and self.available_headers[header]['mode_type'] != "all":
                        raise ConfigError(f"Oh dear. You requested the 'modules_custom' output mode, but gave us a header ({header}) "
                                          "that is suitable only for %s mode(s). Not good." % (self.available_headers[header]['mode_type']))

        outputs_require_ko_dict = [m for m in self.output_modes if self.available_modes[m]['data_dict'] == 'kofams']
        output_string = ", ".join(outputs_require_ko_dict)
        if self.estimate_from_json and len(outputs_require_ko_dict):
            raise ConfigError("You have requested to estimate metabolism from a JSON file and produce the following KOfam hit "
                              f"output mode(s): {output_string}. Unforunately, this is not possible because "
                              "our JSON estimation function does not currently produce the required data for KOfam hit output. "
                              "Please instead request some modules-oriented output mode(s) for your JSON input.")


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

    def init_hits_and_splits(self, annotation_sources=['KOfam']):
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
        split_names_in_contigs_db = set(utils.get_all_item_names_from_the_database(self.contigs_db_path))
        splits_to_use = split_names_in_contigs_db

        # first, resolve differences in splits between profile and contigs db
        if self.profile_db_path:
            self.progress.update("Loading split data from profile DB")
            # if we were given a blank profile, we will assume we want all splits and pull all splits from the contigs DB
            if utils.is_blank_profile(self.profile_db_path):
                self.progress.reset()
                self.run.warning("You seem to have provided a blank profile. No worries, we can still estimate metabolism "
                                 "for you. But we cannot load splits from the profile DB, so instead we are assuming that "
                                 "you are interested in ALL splits and we will load those from the contigs database.")
            else:
                split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.profile_db_path))
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
            if self.collection_name:
                splits_to_use = ccollections.GetSplitNamesInBins(self.args).get_split_names_only()
                self.progress.reset()
                self.run.warning(f"Since a collection name was provided, we will only work with gene calls "
                                 f"from the subset of {len(splits_to_use)} splits in the collection for the "
                                 f"purposes of estimating metabolism.")

        self.progress.update('Loading gene call data from contigs DB')
        from anvio.dbops import ContigsDatabase # <- import here to avoid circular import
        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)

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
                if self.user_input_dir:
                    if not self.ko_dict:
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


    def init_gene_coverage(self, gcids_for_kofam_hits):
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

        if not self.profile_db:
            raise ConfigError("A profile DB has not yet been initialized, so init_gene_coverage() will not work. "
                              "If you are a programmer, you should probably either 1) call this function after "
                              "add_gene_coverage_to_headers_list() or 2) extend this function so that it initializes "
                              "the profile db. If you are not a programmer, you should probably find one :) ")
        self.run.info_single("Since the --add-coverage flag was provided, we are now loading the relevant "
                             "coverage information from the provided profile database.")
        self.profile_db.init_gene_level_coverage_stats_dicts(gene_caller_ids_of_interest=gcids_for_kofam_hits)


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

    def mark_kos_present_for_list_of_splits(self, kofam_hits_in_splits, split_list=None, bin_name=None):
        """This function generates two bin-level dictionaries of dictionaries to store metabolism data.

        The first dictionary of dictionaries is the module completeness dictionary, which associates modules with the KOs
        that are present in the bin for each module.

        The structure of the dictionary is like this example:
        {mnum: {"gene_caller_ids" : set([132, 133, 431, 6777]),
                "kofam_hits" : {'K00033' : [431, 6777],
                                'K01057' : [133],
                                'K00036' : [132] },
                "genes_to_contigs": {132: 0,
                                     133: 0,
                                     431: 2,
                                    6777: 1 },
                "contigs_to_genes": { 0: set([132, 133]),
                                      1: set(6777),
                                      2: set(431) },}}
        This dictionary will be expanded later by other functions.

        The second dictionary of dictionaries is the KOfam hit dictionary, which stores all of the KOfam hits in the bin
        regardless of whether they are in a KEGG module or not.

        The structure of the dictionary is like this example:
        {ko: {"gene_caller_ids" : set([431, 6777]),
              "modules" : ["M00001", "M00555"],                 **Can be None if KO does not belong to any KEGG modules
              "genes_to_contigs": { 431: 2,
                                   6777: 1 },
              "contigs_to_genes": { 1: set(6777),
                                    2: set(431) }}}


        PARAMETERS
        ==========
        kofam_hits_in_splits : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering
        split_list : list
            splits we are considering, this is only for debugging output
        bin_name : str
            name of the bin containing these splits, this is only for debugging output

        RETURNS
        =======
        bin_level_module_dict : dictionary of dictionaries
            initialized module completeness dictionary for the list of splits (genome, metagenome, or bin) provided
        bin_level_ko_dict : dictionary of dictionaries
            dictionary of ko hits within the list of splits provided
        """

        bin_level_module_dict = {}
        bin_level_ko_dict = {}

        if anvio.DEBUG:
            self.run.info("Marking KOs present for bin", bin_name)
            if split_list:
                num_splits = len(split_list)
            else:
                num_splits = "None"
            self.run.info("Number of splits", num_splits)

        # initialize all modules with empty lists and dicts for kos, gene calls
        modules = self.all_modules_in_db.keys()
        all_kos = self.all_kos_in_db.keys()
        for mnum in modules:
            bin_level_module_dict[mnum] = {"gene_caller_ids" : set(),
                                           "kofam_hits" : {},
                                           "genes_to_contigs" : {},
                                           "contigs_to_genes" : {},
                                           "unique_to_this_module": set(),
                                           "warnings" : set()
                                          }

        for knum in all_kos:
            """
            We can only add warnings about missing KOfam profiles because for other annotation sources, we don't
            have a way to know if profiles are missing. But for KOfams with missing profiles, this step is necessary
            so that we don't add the enzyme to the bin_level_ko_dict, because later this will cause problems since
            the enzyme is not in self.ko_dict

            Furthermore, this can only be done when we are using both KEGG data and user data (ie, not --only-user-modules)
            because we need access to the self.ko_dict
            """

            if knum.startswith("M"):
                continue

            if (not self.only_user_modules
                and self.all_kos_in_db[knum]['annotation_source'] == 'KOfam'
                and knum not in self.ko_dict
                and knum not in self.stray_ko_dict
                and f"{knum}{STRAY_KO_ANVIO_SUFFIX}" not in self.stray_ko_dict
                and self.exclude_kos_no_threshold):

                mods_it_is_in = self.all_kos_in_db[knum]['modules']
                if mods_it_is_in:
                    if anvio.DEBUG:
                        mods_str = ", ".join(mods_it_is_in)
                        self.run.warning(f"Oh dear. We do not appear to have a KOfam profile for {knum}. This means "
                                        "that any modules this KO belongs to can never be fully complete (this includes "
                                        f"{mods_str}). ")
                    for m in mods_it_is_in:
                        bin_level_module_dict[m]["warnings"].add(f"No KOfam profile for {knum}")

                if anvio.DEBUG:
                    if self.exclude_kos_no_threshold:
                        self.run.warning(f"We cannot find an entry for KO {knum} in the `ko_list.txt` file downloaded "
                                         f"from KEGG. What this means is that you are somehow using KOfam annotations "
                                         f"that are different from the current version of KOfam on your computer (this can "
                                         f"happen with --enzymes-txt input). Because we are not considering these annotations, "
                                         f"you may get KeyErrors downstream. You can force the inclusion of these KOfams by "
                                         f"re-running this program with the --include-kos-not-in-kofam flag.")

                continue

            bin_level_ko_dict[knum] = {"gene_caller_ids" : set(),
                                     "modules" : None,
                                     "genes_to_contigs" : {},
                                     "contigs_to_genes" : {},
                                     "warnings" : set()
                                     }

        kos_not_in_modules = []

        for ko, gene_call_id, split, contig in kofam_hits_in_splits:
            # make sure we can count annotations to anvi'o versions of stray KO models by using KEGG's original accession
            is_anvio_version = False
            if ko.endswith(STRAY_KO_ANVIO_SUFFIX):
                is_anvio_version = True
                ko = ko.replace(STRAY_KO_ANVIO_SUFFIX, "")

            if ko not in self.all_kos_in_db:

                kos_not_in_modules.append(ko)

                # KOs that are not in modules will not be initialized above in the ko hit dictionary, so we add them here if we haven't already
                if ko not in bin_level_ko_dict:
                    bin_level_ko_dict[ko] = {"gene_caller_ids" : set(),
                                             "modules" : None,
                                             "genes_to_contigs" : {},
                                             "contigs_to_genes" : {},
                                             "warnings" : set()
                                             }
            else:

                # if we are missing the KO from the dictionary at this point, we should fail nicely instead of with a KeyError
                if ko not in bin_level_ko_dict:
                    if self.ignore_unknown_kos:
                        continue
                    raise ConfigError(f"We cannot find the KEGG enzyme {ko} in the dictionary of enzyme hits, even though this enzyme is "
                                      f"annotated in your data. There are 3 main ways this can happen: (1) you are using --enzymes-txt input "
                                      f"that includes KOs that are different from the set used for annotation with `anvi-run-kegg-kofams`, "
                                      f"(2) your contigs database was annotated with `anvi-run-kegg-kofams --include-stray-KOs` but you didn't "
                                      f"use the `--include-stray-KOs` flag for `anvi-estimate-metabolism`, or (3) you imported external annotations "
                                      f"with the source name `KOfam` that include KOs not in the KEGG data directory currently being used. "
                                      f"You have a few options to get around this error depending on which case applies to your situation. "
                                      f"If this is case (2) and you want to include these enzymes in the analysis, then re-run `anvi-estimate-metabolism` "
                                      f"with the flag `--include-stray-KOs`. If this is case (1) or (3) and you want to include these enzymes in the "
                                      f"analysis, then re-run `anvi-estimate-metabolism` with the flag `--include-kos-not-in-kofam`. And no matter "
                                      f"what the situation is, if you want to IGNORE these unknown annotations for the purposes of estimating "
                                      f"metabolism, you can re-run `anvi-estimate-metablism` with the flag `--ignore-unknown-KOs`. If this message "
                                     f"made you worry, you could also re-do your annotations or remove these unknown enzymes from your input "
                                     f"--enzymes-txt file to be on the safe side.")

                present_in_mods = self.all_kos_in_db[ko]['modules']
                bin_level_ko_dict[ko]["modules"] = present_in_mods

                # keep track of enzymes unique to this module
                is_unique = False
                if len(present_in_mods) == 1:
                    is_unique = True

                for m in present_in_mods:
                    bin_level_module_dict[m]["gene_caller_ids"].add(gene_call_id)
                    if ko in bin_level_module_dict[m]["kofam_hits"] and gene_call_id not in bin_level_module_dict[m]["kofam_hits"][ko]:
                        bin_level_module_dict[m]["kofam_hits"][ko].append(gene_call_id)
                    else:
                        bin_level_module_dict[m]["kofam_hits"][ko] = [gene_call_id]
                    bin_level_module_dict[m]["genes_to_contigs"][gene_call_id] = contig
                    if contig in bin_level_module_dict[m]["contigs_to_genes"]:
                        bin_level_module_dict[m]["contigs_to_genes"][contig].add(gene_call_id)
                    else:
                        bin_level_module_dict[m]["contigs_to_genes"][contig] = set([gene_call_id])

                    # make a special list for the enzymes that are unique
                    if is_unique:
                        bin_level_module_dict[m]["unique_to_this_module"].add(ko)
                    # warn the user if this enzyme is shared between multiple modules
                    else:
                        mod_str = "/".join(present_in_mods)
                        bin_level_module_dict[m]["warnings"].add(f"{ko} is present in multiple modules: {mod_str}")

                    # point out use of anvi'o-specific KO models
                    if is_anvio_version:
                        warning = f"used '{ko}{STRAY_KO_ANVIO_SUFFIX}' model to annotate {ko}"
                        bin_level_module_dict[m]["warnings"].add(warning)
                        bin_level_ko_dict[ko]["warnings"].add(warning)


            bin_level_ko_dict[ko]["gene_caller_ids"].add(gene_call_id)
            bin_level_ko_dict[ko]["genes_to_contigs"][gene_call_id] = contig
            if contig in bin_level_ko_dict[ko]["contigs_to_genes"]:
                bin_level_ko_dict[ko]["contigs_to_genes"][contig].add(gene_call_id)
            else:
                bin_level_ko_dict[ko]["contigs_to_genes"][contig] = set([gene_call_id])

        if anvio.DEBUG:
            self.run.info("Gene calls processed", "%d in bin" % len(kofam_hits_in_splits))
            if kos_not_in_modules:
                self.run.warning("Just so you know, the following enzymes did not belong to any modules in the MODULES.db: %s"
                % ", ".join(kos_not_in_modules))

        return bin_level_module_dict, bin_level_ko_dict


    def compute_stepwise_module_completeness_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the stepwise completeness of the specified module within the given bin dictionary.

        It uses only the "top-level" steps of the module definition, which are the steps that you get when you first
        split the module definition on a space. Each "top-level" step is comprised of one or more enzymes that either
        work together or serve as alternatives to each other. In this calculation, we ignore the possible combinations
        of enzymes and simply decide whether or not a "top-level" step is complete or not. Then the module completeness
        is computed as the number of complete "top-level" steps divided by the total number of steps.

        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "stepwise_completeness"         the stepwise completeness of the module
        "stepwise_is_complete"          whether the module completeness falls over the completeness threshold
        "top_level_step_info"           a dictionary of each top level step
                                            keyed by integer from 0 to # of top level steps.
                                            inner dict contains the following keys:
                                            'step_definition' (string)
                                            'complete' (Boolean)
                                            'includes_modules' (Boolean)
                                            'included_module_list' (list of strings)

        RETURNS
        =======
        over_complete_threshold : boolean
            whether or not the module is considered "complete" overall based on the threshold fraction of completeness
        """

        top_level_steps = self.all_modules_in_db[mnum]['top_level_steps']
        num_steps = len(top_level_steps)
        num_complete = 0
        num_nonessential_steps = 0

        present_list_for_mnum = meta_dict_for_bin[mnum]["kofam_hits"].keys()

        meta_dict_for_bin[mnum]['top_level_step_info'] = {}

        for i, step in enumerate(top_level_steps):
            step_is_present_condition_statement = ""
            cur_index = 0  # current position in the step DEFINITION
            step_includes_modules = False
            included_module_list = []

            while cur_index < len(step):
                if step[cur_index] in ['(',')']:
                    step_is_present_condition_statement += step[cur_index]
                    cur_index += 1

                elif step[cur_index] == ",":
                    step_is_present_condition_statement += " or "
                    cur_index += 1

                elif step[cur_index] == "+" or step[cur_index] == ' ':
                    step_is_present_condition_statement += " and "
                    cur_index += 1

                elif step[cur_index] == "-":
                    # '--' no associated enzyme case, by default False (assumed incomplete)
                    if step[cur_index+1] == "-":
                        if self.exclude_dashed_reactions: # skip it instead
                            cur_index += 3 # skip over both '-' AND the following space
                        else:
                            step_is_present_condition_statement += "False"
                            cur_index += 2 # skip over both '-', the next character should be a space or end of DEFINITION line

                        if anvio.DEBUG:
                            self.run.warning(f"While estimating the stepwise completeness of KEGG module {mnum}, anvi'o saw "
                                             f"'--' in the module DEFINITION. This indicates a step in the pathway that has no "
                                             f"associated enzyme. By default, anvi'o marks steps like these incomplete, *unless* "
                                             f"you are using the flag --exclude-dashed-reactions. But if you aren't using that flag, "
                                             f"it is possible that this module might be falsely considered incomplete. So it may be in your "
                                             f"interest to take a closer look at module {mnum}.")
                        if cur_index < len(step) and step[cur_index] != " ":
                            raise ConfigError(f"Serious, serious parsing sadness is happening. We just processed a '--' in "
                                              f"a DEFINITION line for module {mnum} but did not see a space afterwards. Instead, "
                                              f"we found {step[cur_index+1]}. WHAT DO WE DO NOW?")

                    # a whole set of nonessential KOs - skip all of them
                    elif step[cur_index+1] == "(":
                        while step[cur_index] != ")":
                            cur_index += 1
                        cur_index += 1 # skip over the ')'

                    # anything else that follows a '-' should be an enzyme or enzyme component, and should be skipped
                    else:
                        # find the next space or '-' or the end of the step
                        while cur_index+1 < len(step) and (step[cur_index+1] not in [' ', ',', '+', '-', '(', ')']):
                            cur_index += 1
                        cur_index += 1
                        # if we found a non-accession character, the next iteration of the loop will take care of it
                        # if we reached the end of the step and the condition statement is empty, then the entire
                        #    step is nonessential so we need to avoid counting it (taken care of later)

                else: # enzyme or module accession
                    enzyme_start_index = cur_index
                    while cur_index+1 < len(step) and step[cur_index+1] not in [' ', ',', '+', '-', '(', ')']:
                        cur_index += 1
                    enzyme_end_index = cur_index

                    accession = step[enzyme_start_index : enzyme_end_index+1]
                    # module
                    if accession in self.all_modules_in_db:
                        step_includes_modules = True
                        included_module_list.append(accession)
                        # store the module accession in the condition string to be replaced later
                        step_is_present_condition_statement += accession
                    # enzyme
                    elif accession in present_list_for_mnum:
                        step_is_present_condition_statement += "True"
                    else:
                        step_is_present_condition_statement += "False"

                    cur_index += 1

            meta_dict_for_bin[mnum]['top_level_step_info'][i] = {}
            meta_dict_for_bin[mnum]['top_level_step_info'][i]['step_definition'] = step
            meta_dict_for_bin[mnum]['top_level_step_info'][i]['includes_modules'] = step_includes_modules
            meta_dict_for_bin[mnum]['top_level_step_info'][i]['included_module_list'] = included_module_list

            # entire step was nonessential, do not count it
            if step_is_present_condition_statement == "":
                num_nonessential_steps += 1
                meta_dict_for_bin[mnum]['top_level_step_info'][i]['complete'] = "nonessential"
            elif step_includes_modules:
                # we'll eval this condition statement in a later function once all other modules have stepwise completeness
                meta_dict_for_bin[mnum]['top_level_step_info'][i]['complete'] = step_is_present_condition_statement
            else:
                step_is_present = eval(step_is_present_condition_statement)
                meta_dict_for_bin[mnum]['top_level_step_info'][i]['complete'] = step_is_present
                if step_is_present:
                    num_complete += 1

        # compute stepwise completeness as number of complete (essential) steps / number of total (essential) steps
        mod_stepwise_completeness = num_complete / (num_steps - num_nonessential_steps)
        meta_dict_for_bin[mnum]["stepwise_completeness"] = mod_stepwise_completeness

        over_complete_threshold = True if meta_dict_for_bin[mnum]["stepwise_completeness"] >= self.module_completion_threshold else False
        meta_dict_for_bin[mnum]["stepwise_is_complete"] = over_complete_threshold

        return over_complete_threshold


    def compute_pathwise_module_completeness_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the pathwise completeness of the specified module within the given bin metabolism dictionary.

        To do this, it works with the unrolled module definition: a list of all possible paths, where each path is a list of atomic steps.
        Atomic steps include singular KOs, protein complexes, modules, non-essential steps, and steps without associated KOs.
        An atomic step (or parts of a protein complex) can be considered 'present' if the corresponding KO(s) has a hit in the bin.
        For each path, the function computes the path completeness as the number of present (essential) steps divided by the number of total steps in the path.
        The module completeness is simply the highest path completeness.

        There are some special cases to consider here.
        1) Non-essential steps. These are steps that are marked with a preceding "-" to indicate that they are not required for the module to
           be considered complete. They often occur in pathways with multiple forks. What we do with these is save and count them separately as
           non-essential steps, but we do not use them in our module completeness calculations. Another thing we do is continue parsing the rest
           of the module steps as normal, even though some of them may affect steps after the non-essential one. That may eventually change.
           See comments in the code below.
        2) Steps without associated KOs. These are steps marked as "--". They may require an enzyme, but if so that enzyme is not in the KOfam
           database, so we can't know whether they are complete or not from our KOfam hits. Therefore, we assume these steps are incomplete, and
           warn the user to go back and check the module manually.
        3) Steps defined by entire modules. These steps have module numbers instead of KOs, so they require an entire module to be complete in
           order to be complete. We can't figure this out until after we've evaluated all modules, so we simply parse these steps without marking
           them complete, and later will go back to adjust the completeness score once all modules have been marked complete or not.


        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "pathway_completeness"          a list of the completeness of each pathway
        "present_nonessential_kos"      a list of non-essential KOs in the module that were found to be present
        "most_complete_paths"           a list of the paths with maximum completeness
        "pathwise_percent_complete"     the completeness of the module, which is the maximum pathway completeness
        "pathwise_is_complete"          whether the module completeness falls over the completeness threshold

        RETURNS
        =======
        over_complete_threshold : boolean
            whether or not the module is considered "complete" overall based on the threshold fraction of completeness
        has_nonessential_step : boolean
            whether or not the module contains non-essential steps. Used for warning the user about these.
        has_no_ko_step : boolean
            whether or not the module contains steps without associated KOs. Used for warning the user about these.
        defined_by_modules : boolean
            whether or not the module contains steps defined by other modules. Used for going back to adjust completeness later.
        """

        present_list_for_mnum = meta_dict_for_bin[mnum]["kofam_hits"].keys()
        if not present_list_for_mnum:
            # no KOs in this module are present
            if anvio.DEBUG:
                self.run.warning("No KOs present for module %s. Parsing for completeness is still being done to obtain module information." % mnum)

        # stuff to put in the module's dictionary
        module_nonessential_kos = [] # KOs that are present but unnecessary for module completeness

        # stuff that will be returned
        over_complete_threshold = False
        has_nonessential_step = False
        has_no_ko_step = False
        defined_by_modules = False

        meta_dict_for_bin[mnum]["pathway_completeness"] = []
        meta_dict_for_bin[mnum]["num_complete_copies_of_all_paths"] = []
        meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"] = []

        for p in self.module_paths_dict[mnum]:
            num_complete_steps_in_path = 0
            num_nonessential_steps_in_path = 0 # so that we don't count nonessential steps when computing completeness
            atomic_step_copy_number = []

            for atomic_step in p:
                # there are 5 types of atomic steps to take care of
                if any(x in atomic_step for x in ['-','+']):
                    # 1) steps without associated enzymes, ie --
                    if atomic_step == "--":
                        # when '--' in a DEFINITION line happens, it signifies a reaction step that has no associated enzyme.
                        # by default, we assume that such steps are not complete
                        has_no_ko_step = True
                        if self.exclude_dashed_reactions:
                            warning_str = "'--' step was ignored in the calculation"
                            num_nonessential_steps_in_path += 1 # this is to ensure we fix the denominator later
                        else:
                            warning_str = "'--' steps are assumed incomplete"
                            atomic_step_copy_number.append(0)

                        meta_dict_for_bin[mnum]["warnings"].add(warning_str)
                    # 2) non-essential KOs, ie -Kxxxxx
                    elif atomic_step[0] == "-" and not any(x in atomic_step[1:] for x in ['-','+']):
                        """
                        OKAY, SO HERE WE HAVE SOME POOPINESS THAT MAY NEED TO BE FIXED EVENTUALLY.
                        Basically, some DEFINITION lines have KOs that seem to be marked non-essential;
                        ie, "-K11024" in "K11023 -K11024 K11025 K11026 K11027".
                        It was difficult to decide whether we should consider only K11024, or K11024 and all following KOs, to be non-essential.
                        For instance, the module M00778 is a complex case that gave us pause - see Fiesta issue 955.
                        But for now, we have decided to just track only the one KO as a 'non-essential step', and to not include such steps in
                        the module completeness estimate.
                        """
                        ko = atomic_step[1:]
                        if ko not in module_nonessential_kos:
                            module_nonessential_kos.append(ko)
                        num_nonessential_steps_in_path += 1
                        has_nonessential_step = True

                    # 3) protein complexes, ie Kxxxxx+Kyyyyy-Kzzzzz (2 types of complex components - essential and nonessential)
                    else:
                        # split on '+' or '-'
                        pattern = re.compile('\+|\-')
                        match_idxs = []
                        for match in re.finditer(pattern, atomic_step):
                            match_idxs.append(match.start())

                        essential_components = []
                        num_matches_processed = 0
                        for i, match_idx in enumerate(match_idxs):
                            # if this is the first match, we need to handle the initial component in the complex
                            if num_matches_processed == 0:
                                essential_components.append(atomic_step[0:match_idx])

                            # handle the component after the match character
                            if i < len(match_idxs)-1:
                                next_idx = match_idxs[i+1]
                            else:
                                next_idx = len(atomic_step)
                            component_ko = atomic_step[match_idx+1:next_idx]

                            # essential component after  +
                            if atomic_step[match_idx] == '+':
                                essential_components.append(component_ko)
                            # non-essential component after '-'
                            else:
                                has_nonessential_step = True
                                if component_ko not in module_nonessential_kos:
                                    module_nonessential_kos.append(component_ko)

                            num_matches_processed += 1

                        # after processing all components of the enzyme complex, we compute the complex completeness and copy number
                        num_present_components = 0
                        component_copy_number = []
                        for c in essential_components:
                            if c in present_list_for_mnum:
                                num_present_components += 1
                                num_copies = len(meta_dict_for_bin[mnum]["kofam_hits"][c])
                            else:
                                num_copies = 0
                            component_copy_number.append(num_copies)
                        component_completeness = num_present_components / len(essential_components)
                        num_complete_steps_in_path += component_completeness

                        if component_completeness >= self.module_completion_threshold:
                            atomic_step_copy_number.append(min(component_copy_number))
                        else:
                            atomic_step_copy_number.append(0)
                else:
                    # atomic step is a single enzyme or module
                    # 4) Module numbers, ie Mxxxxx
                    if atomic_step in self.all_modules_in_db:
                        """
                        This happens when a module is defined by other modules. For example, photosynthesis module M00611 is defined as
                        (M00161,M00163) M00165 === (photosystem II or photosystem I) and calvin cycle

                        We need all the modules to have been evaluated before we can determine completeness of steps with module numbers.
                        So what we will do here is to use a flag variable to keep track of the modules that have this sort of definition
                        in a list so we can go back and evaluate completeness of steps with module numbers later.
                        """
                        defined_by_modules = True
                    # 5) regular old single enzymes, ie Kxxxxx (for KOs), COGyyyyy (for COGs), etc
                    else:
                        if atomic_step in present_list_for_mnum:
                            num_complete_steps_in_path += 1
                            num_copies = len(meta_dict_for_bin[mnum]["kofam_hits"][atomic_step])
                        else:
                            num_copies = 0

                        atomic_step_copy_number.append(num_copies)

            path_completeness = num_complete_steps_in_path / (len(p) - num_nonessential_steps_in_path)
            meta_dict_for_bin[mnum]["pathway_completeness"].append(path_completeness)

            # compute path copy number
            if defined_by_modules:
                path_copy_number = atomic_step_copy_number # save list with atomic step copy numbers to use when adjusting module copy number later
            else:
                path_copy_number = self.compute_num_complete_copies_of_path(atomic_step_copy_number)
            meta_dict_for_bin[mnum]["num_complete_copies_of_all_paths"].append(path_copy_number)

        # once all paths have been evaluated, we find the path(s) of maximum completeness and set that as the overall module completeness
        # this is not very efficient as it takes two passes over the list but okay
        meta_dict_for_bin[mnum]["pathwise_percent_complete"] = max(meta_dict_for_bin[mnum]["pathway_completeness"])
        if meta_dict_for_bin[mnum]["pathwise_percent_complete"] > 0:
            meta_dict_for_bin[mnum]["most_complete_paths"] = [self.module_paths_dict[mnum][i] for i, pc in enumerate(meta_dict_for_bin[mnum]["pathway_completeness"]) if pc == meta_dict_for_bin[mnum]["pathwise_percent_complete"]]
            if not defined_by_modules:
                meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"] = [meta_dict_for_bin[mnum]["num_complete_copies_of_all_paths"][i] for i, pc in enumerate(meta_dict_for_bin[mnum]["pathway_completeness"]) if pc == meta_dict_for_bin[mnum]["pathwise_percent_complete"]]
        else:
            meta_dict_for_bin[mnum]["most_complete_paths"] = []
            meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"] = []

        # set module copy number as the maximum copy number of the path(s) of maximum completeness
        if meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"]:
            meta_dict_for_bin[mnum]["pathwise_copy_number"] = max(meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"])
        else:
            meta_dict_for_bin[mnum]["pathwise_copy_number"] = 'NA'

        # compute proportion of unique enzymes in the module (regardless of which path(s) enzyme is in or whether enzyme is essential)
        if meta_dict_for_bin[mnum]["unique_to_this_module"]:
            num_unique_enzymes_present = 0
            num_unique_enzymes_in_mod = len(meta_dict_for_bin[mnum]["unique_to_this_module"])
            for ko in present_list_for_mnum:
                if ko in meta_dict_for_bin[mnum]["unique_to_this_module"]:
                    num_unique_enzymes_present += 1

            meta_dict_for_bin[mnum]["proportion_unique_enzymes_present"] = num_unique_enzymes_present / num_unique_enzymes_in_mod
            meta_dict_for_bin[mnum]["unique_enzymes_context_string"] = f"{num_unique_enzymes_present} of {num_unique_enzymes_in_mod} unique enzymes in module"
        else:
            meta_dict_for_bin[mnum]["proportion_unique_enzymes_present"] = "NA"
            meta_dict_for_bin[mnum]["unique_enzymes_context_string"] = "NA"


        if anvio.DEBUG and len(meta_dict_for_bin[mnum]["most_complete_paths"]) > 1:
            self.run.warning("Found %d complete paths for module %s with completeness %s. " % (len(meta_dict_for_bin[mnum]["most_complete_paths"]), mnum, meta_dict_for_bin[mnum]["pathwise_percent_complete"]),
                            header='DEBUG OUTPUT', lc='yellow')
        over_complete_threshold = True if meta_dict_for_bin[mnum]["pathwise_percent_complete"] >= self.module_completion_threshold else False
        meta_dict_for_bin[mnum]["pathwise_is_complete"] = over_complete_threshold
        meta_dict_for_bin[mnum]["present_nonessential_kos"] = module_nonessential_kos

        return over_complete_threshold, has_nonessential_step, has_no_ko_step, defined_by_modules


    def adjust_stepwise_completeness_for_bin(self, mnum, meta_dict_for_bin):
        """This function adjusts stepwise completeness of modules that are defined by other modules.

        This can only be done after all other modules have been evaluated for completeness.
        The function goes through the top-level steps established by compute_stepwise_module_completeness_for_bin()
        and re-assesses whether steps including other modules are complete. It updates the metabolism completess dictionary accordingly.

        PARAMETERS
        ==========
        mnum : string
            the module number to adjust
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dictionary for the current bin

        RETURNS
        =======
        now_complete : boolean
            whether or not the module is NOW considered "complete" overall based on the threshold fraction of completeness
        """

        num_steps = len(meta_dict_for_bin[mnum]['top_level_step_info'].keys())
        num_complete = 0
        num_nonessential_steps = 0

        for step_num, step_dict in meta_dict_for_bin[mnum]['top_level_step_info'].items():
            if step_dict['includes_modules']:
                # this condition statement has module accessions in it, we need to replace those with True/False and eval
                step_is_present_condition_statement = step_dict['complete']
                module_accessions_to_replace = step_dict['included_module_list']
                for m in module_accessions_to_replace:
                    mod_completeness = meta_dict_for_bin[m]['stepwise_completeness']
                    if mod_completeness >= self.module_completion_threshold:
                        step_is_present_condition_statement = step_is_present_condition_statement.replace(m, "True")
                    else:
                        step_is_present_condition_statement = step_is_present_condition_statement.replace(m, "False")

                # now evaluate to see if this step is complete
                step_is_present = eval(step_is_present_condition_statement)
                meta_dict_for_bin[mnum]['top_level_step_info'][step_num]['complete'] = step_is_present
                if step_is_present:
                    num_complete += 1

            else:
                if step_dict['complete'] == "nonessential":
                    num_nonessential_steps += 1
                elif step_dict['complete']:
                    num_complete += 1

        mod_stepwise_completeness = num_complete / (num_steps - num_nonessential_steps)
        meta_dict_for_bin[mnum]["stepwise_completeness"] = mod_stepwise_completeness

        now_complete = True if mod_stepwise_completeness >= self.module_completion_threshold else False
        meta_dict_for_bin[mnum]["stepwise_is_complete"] = now_complete

        return now_complete


    def adjust_pathwise_completeness_for_bin(self, mod, meta_dict_for_bin):
        """This function adjusts pathwise completeness of modules that are defined by other modules.

        This can only be done after all other modules have been evaluated for completeness.
        The function uses similar logic as compute_pathwise_module_completeness_for_bin() to re-assess whether steps defined
        by other modules are complete, and updates the metabolism completess dictionary accordingly.

        PARAMETERS
        ==========
        mod : string
            the module number to adjust
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dictionary for the current bin

        RETURNS
        =======
        now_complete : boolean
            whether or not the module is NOW considered "complete" overall based on the threshold fraction of completeness
        """

        for i in range(len(self.module_paths_dict[mod])):
            p = self.module_paths_dict[mod][i]
            num_essential_steps_in_path = 0  # note that the len(p) will include nonessential steps; we should count only essential ones
            num_complete_module_steps = 0

            # take previously computed step copy numbers. This list includes all steps except those defined by modules.
            atomic_step_copy_numbers_in_path = meta_dict_for_bin[mod]["num_complete_copies_of_all_paths"][i]
            module_copy_num_should_be_NA = False # flag to indicate whether component modules have a copy number of NA

            for atomic_step in p:
                # module step; we need to count these based on previously computed module completeness
                if atomic_step in self.all_modules_in_db:
                    num_complete_module_steps += meta_dict_for_bin[atomic_step]["pathwise_percent_complete"]
                    num_essential_steps_in_path += 1

                    if meta_dict_for_bin[atomic_step]["pathwise_copy_number"] == "NA":
                        module_copy_num_should_be_NA = True
                    else:
                        atomic_step_copy_numbers_in_path.append(meta_dict_for_bin[atomic_step]["pathwise_copy_number"])
                # non-essential KO, don't count as a step in the path
                elif atomic_step[0] == '-' and not atomic_step == "--":
                    pass
                # single enzymes, protein complexes and '--' steps; were already counted as complete by previous function
                else:
                    num_essential_steps_in_path += 1

            # now we adjust the previous pathway completeness
            old_complete_steps_in_path = meta_dict_for_bin[mod]["pathway_completeness"][i] * num_essential_steps_in_path
            adjusted_num_complete_steps_in_path = old_complete_steps_in_path + num_complete_module_steps
            meta_dict_for_bin[mod]["pathway_completeness"][i] = adjusted_num_complete_steps_in_path / num_essential_steps_in_path

            # now we adjust the path copy number
            if module_copy_num_should_be_NA:
                path_copy_number = "NA"
            else:
                path_copy_number = self.compute_num_complete_copies_of_path(atomic_step_copy_numbers_in_path)
            meta_dict_for_bin[mod]["num_complete_copies_of_all_paths"][i] = path_copy_number

        # after adjusting for all paths, adjust overall module completeness
        meta_dict_for_bin[mod]["pathwise_percent_complete"] = max(meta_dict_for_bin[mod]["pathway_completeness"])
        if meta_dict_for_bin[mod]["pathwise_percent_complete"] > 0:
            meta_dict_for_bin[mod]["most_complete_paths"] = [self.module_paths_dict[mod][i] for i, pc in enumerate(meta_dict_for_bin[mod]["pathway_completeness"]) if pc == meta_dict_for_bin[mod]["pathwise_percent_complete"]]
        else:
            meta_dict_for_bin[mod]["most_complete_paths"] = []

        now_complete = True if meta_dict_for_bin[mod]["pathwise_percent_complete"] >= self.module_completion_threshold else False
        meta_dict_for_bin[mod]["pathwise_is_complete"] = now_complete

        # and adjust overall module copy number
        if meta_dict_for_bin[mod]["num_complete_copies_of_most_complete_paths"]:
            meta_dict_for_bin[mod]["pathwise_copy_number"] = max(meta_dict_for_bin[mod]["num_complete_copies_of_most_complete_paths"])
        else:
            meta_dict_for_bin[mod]["pathwise_copy_number"] = 'NA'

        return now_complete


    def add_module_coverage(self, mod, meta_dict_for_bin):
        """This function updates the metabolism dictionary with coverage values for the given module.

        It must be called after init_gene_coverage() or add_gene_coverage_to_headers_list() so that
        the self.profile_db attribute is established.

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "genes_to_coverage"             dictionary of mean coverage in each sample for each gene
                                        coverage = meta_dict_for_bin[module]["genes_to_coverage"][sample][gcid]
        "genes_to_detection"            dictionary of detection in each sample for each gene
                                        detection = meta_dict_for_bin[module]["genes_to_detection"][sample][gcid]
        "average_coverage_per_sample"   dictionary of average mean coverage of all genes in module, per sample
                                        avg_coverage = meta_dict_for_bin[module]["average_coverage_per_sample"][sample]
        "average_detection_per_sample"  dictionary of average detection of all genes in module, per sample
                                        avg_detection = meta_dict_for_bin[module]["average_detection_per_sample"][sample]
        """

        meta_dict_for_bin[mod]["genes_to_coverage"] = {}
        meta_dict_for_bin[mod]["genes_to_detection"] = {}
        meta_dict_for_bin[mod]["average_coverage_per_sample"] = {}
        meta_dict_for_bin[mod]["average_detection_per_sample"] = {}

        if not self.enzymes_of_interest_df and not self.profile_db:
            raise ConfigError("The add_module_coverage() function cannot work without a properly initialized "
                              "profile database.")

        if self.custom_output_headers:
            # determine the specific set of samples we are interested in so we don't make the dictionary huge
            sample_set = set()
            for h in self.custom_output_headers:
                if 'coverage' in h or 'detection' in h:
                    if '_gene_coverages' in h:
                        sample = h.replace('_gene_coverages', '')
                    elif '_avg_coverage' in h:
                        sample = h.replace('_avg_coverage', '')
                    elif '_gene_detection' in h:
                        sample = h.replace('_gene_detection', '')
                    elif '_avg_detection' in h:
                        sample = h.replace('_avg_detection', '')
                    sample_set.add(sample)
            self.coverage_sample_list = list(sample_set)
        else:
            if self.enzymes_of_interest_df is not None:
                self.coverage_sample_list = [self.contigs_db_project_name]
            else:
                self.coverage_sample_list = self.profile_db.p_meta['samples']

        num_genes = len(meta_dict_for_bin[mod]["gene_caller_ids"])
        for s in self.coverage_sample_list:
            meta_dict_for_bin[mod]["genes_to_coverage"][s] = {}
            meta_dict_for_bin[mod]["genes_to_detection"][s] = {}
            coverage_sum = 0
            detection_sum = 0
            for g in meta_dict_for_bin[mod]["gene_caller_ids"]:
                if self.enzymes_of_interest_df is not None:
                    cov = self.enzymes_of_interest_df[self.enzymes_of_interest_df['gene_id'] == g]['coverage'].values[0]
                    det = self.enzymes_of_interest_df[self.enzymes_of_interest_df['gene_id'] == g]['detection'].values[0]
                else:
                    cov = self.profile_db.gene_level_coverage_stats_dict[g][s]['mean_coverage']
                    det = self.profile_db.gene_level_coverage_stats_dict[g][s]['detection']
                coverage_sum += cov
                detection_sum += det
                meta_dict_for_bin[mod]["genes_to_coverage"][s][g] = cov
                meta_dict_for_bin[mod]["genes_to_detection"][s][g] = det

            if num_genes == 0:
                meta_dict_for_bin[mod]["average_coverage_per_sample"][s] = 0
                meta_dict_for_bin[mod]["average_detection_per_sample"][s] = 0
            else:
                meta_dict_for_bin[mod]["average_coverage_per_sample"][s] = coverage_sum / num_genes
                meta_dict_for_bin[mod]["average_detection_per_sample"][s] = detection_sum / num_genes


    def estimate_for_list_of_splits(self, metabolism_dict_for_list_of_splits, bin_name=None):
        """This is the atomic metabolism estimator function, which builds up the metabolism completeness dictionary for an arbitrary list of splits.

        For example, the list of splits may represent a bin, a single isolate genome, or an entire metagenome.

        The function takes in a metabolism completeness dictionary already initialized with the relevant KOfam hits per module, and updates it
        with the individual steps and completion estimates for each module.

        PARAMETERS
        ==========
        metabolism_dict_for_list_of_splits : dictionary of dictionaries
            the metabolism completeness dictionary of dictionaries for this list of splits. It contains
            one dictionary of module steps and completion information for each module (keyed by module number).
            Calling functions should assign this dictionary to a metabolism superdict with the bin name as a key.
        bin_name : str
            the name of the bin/genome/metagenome that we are working with
        """

        pathwise_complete_mods = set([])
        stepwise_complete_mods = set([])
        mods_def_by_modules = [] # a list of modules that have module numbers in their definitions
        # modules to warn about
        mods_with_unassociated_ko = [] # a list of modules that have "--" steps without an associated KO
        mods_with_nonessential_steps = [] # a list of modules that have nonessential steps like "-K11024"

        # estimate completeness of each module
        for mod in metabolism_dict_for_list_of_splits.keys():
            # pathwise
            mod_is_complete, has_nonessential_step, has_no_ko_step, defined_by_modules \
            = self.compute_pathwise_module_completeness_for_bin(mod, metabolism_dict_for_list_of_splits)

            if mod_is_complete:
                pathwise_complete_mods.add(mod)
            if has_nonessential_step:
                mods_with_nonessential_steps.append(mod)
            if has_no_ko_step:
                mods_with_unassociated_ko.append(mod)
            if defined_by_modules:
                mods_def_by_modules.append(mod)

            # stepwise
            mod_is_complete = self.compute_stepwise_module_completeness_for_bin(mod, metabolism_dict_for_list_of_splits)
            self.compute_stepwise_module_copy_number_for_bin(mod, metabolism_dict_for_list_of_splits)

            if mod_is_complete:
                stepwise_complete_mods.add(mod)

            if self.add_coverage:
                self.add_module_coverage(mod, metabolism_dict_for_list_of_splits)

        # go back and adjust completeness/copy number of modules that are defined by other modules
        if mods_def_by_modules:
            for mod in mods_def_by_modules:
                # pathwise
                mod_is_complete = self.adjust_pathwise_completeness_for_bin(mod, metabolism_dict_for_list_of_splits)
                if mod_is_complete:
                    pathwise_complete_mods.add(mod)
                # stepwise
                mod_is_complete = self.adjust_stepwise_completeness_for_bin(mod, metabolism_dict_for_list_of_splits)
                if mod_is_complete:
                    stepwise_complete_mods.add(mod)
                self.adjust_stepwise_copy_number_for_bin(mod, metabolism_dict_for_list_of_splits)


        # estimate redundancy of each module
        for mod in metabolism_dict_for_list_of_splits.keys():
            self.compute_module_redundancy_for_bin(mod, metabolism_dict_for_list_of_splits)


        # notify user of the modules that gave some fishy results -- but only for genome mode because it's too wordy otherwise
        if not self.quiet and self.genome_mode:
            if mods_with_nonessential_steps:
                self.run.warning("Please note that anvi'o found one or more non-essential steps in the following modules: %s.   "
                                 "At this time, we are not counting these steps in our percent completion estimates."
                                 % (", ".join(mods_with_nonessential_steps)))

            if mods_with_unassociated_ko:
                self.run.warning("Just so you know, while estimating the completeness of some modules, anvi'o saw "
                                 "'--' in the module DEFINITION. This indicates a step in the pathway that has no "
                                 "associated enzyme. So we really cannot know just based on gene annotations whether or not this "
                                 "step is present. By default, anvi'o marks these steps incomplete. But they may not be, "
                                 "and as a result their modules may be falsely considered incomplete. So it may be in your "
                                 "interest to go back and take a look at these individual modules to see if you can find the "
                                 "missing enzyme in some other way. Best of luck to you. Here is the list of modules to check out: %s"
                                 % (", ".join(mods_with_unassociated_ko)))

        if anvio.DEBUG or self.genome_mode:
            self.run.info("Bin name", bin_name)
            self.run.info("Module completion threshold", self.module_completion_threshold)
            self.run.info("Number of complete modules (pathwise)", len(pathwise_complete_mods))
            self.run.info("Number of complete modules (stepwise)", len(stepwise_complete_mods))
            if pathwise_complete_mods:
                self.run.info("Pathwise complete modules", ", ".join(sorted(list(pathwise_complete_mods))))
            if stepwise_complete_mods:
                self.run.info("Stepwise complete modules", ", ".join(sorted(list(stepwise_complete_mods))))

        return metabolism_dict_for_list_of_splits

######### REDUNDANCY FUNCTIONS (UNUSED IN NON-JSON OUTPUT) #########

    def compute_naive_redundancy_for_path(self, num_ko_hits_in_path):
        """This function computes a naive redundancy measure for a module path, given the number of hits per KO in the path.

        naive redundancy = # extra hits / len(path) where a hit is "extra" if it is not the first hit to the KO.

        PARAMETERS
        ==========
        num_ko_hits_in_path : list
            stores the number of copies of each enzyme in path
        """

        extra_hits = [x - 1 if x > 1 else 0 for x in num_ko_hits_in_path]
        return sum(extra_hits)/len(num_ko_hits_in_path)


    def compute_copywise_redundancy_for_path(self, num_ko_hits_in_path, aggregation_measure="average"):
        """This function computes redundancy based on the completeness of each extra copy of a path.

        The 'base' redundancy score is determined by the number of extra copies with 100% completeness.
        The completeness measurements of all other extra copies are aggregated (using the aggregation_measure) and
        added to this 'base' redundancy to get the overall path redundancy.

        PARAMETERS
        ==========
        num_ko_hits_in_path : list
            stores the number of copies of each enzyme in path
        """

        accepted_aggregation_measures = ["average", "median", "weighted_sum", "geometric_mean"]
        extra_hits = [x - 1 if x > 1 else 0 for x in num_ko_hits_in_path]
        base_redundancy = min(extra_hits) # number of extra copies of path that are 100% complete
        extra_copy_completeness = []
        # here we get the completeness of every extra copy of the path
        for i in range((base_redundancy+1), max(extra_hits) + 1):
            num_present_kos_in_copy = len([num_hits for num_hits in extra_hits if num_hits >= i])
            extra_copy_completeness.append(num_present_kos_in_copy/len(num_ko_hits_in_path))

        aggregated_completeness = None
        if not extra_copy_completeness: # this handles the case when ALL extra copies are 100% complete
            aggregated_completeness = 0
        else:
            if aggregation_measure == "average":
                aggregated_completeness = statistics.mean(extra_copy_completeness)
            elif aggregation_measure == "median":
                aggregated_completeness = statistics.median(extra_copy_completeness)
            elif aggregation_measure == "weighted_sum":
                aggregated_completeness = 0
                for c in range(len(extra_copy_completeness)):
                    aggregated_completeness += 1/(c+1) * extra_copy_completeness[c]
            elif aggregation_measure == "geometric_mean":
                aggregated_completeness = stats.gmean(extra_copy_completeness)
            else:
                raise ConfigError("The function compute_copywise_redundancy_for_path() doesn't know how to handle the aggregation measure '%s'. "
                                  "Accepted aggregation measures include: %s " % (aggregation_measure, ", ".join(accepted_aggregation_measures)))

        return (base_redundancy + aggregated_completeness), extra_copy_completeness


    def compute_entropy_weighted_redundancy_for_bin(self, num_ko_hits_in_path):
        """This function computes naive redundancy but weights it by the entropy of the hit distribution.

        PARAMETERS
        ==========
        num_ko_hits_in_path : list
            stores the number of copies of each enzyme in path
        """

        extra_hits = [x - 1 if x > 1 else 0 for x in num_ko_hits_in_path]
        total_extra_hits = sum(extra_hits)
        num_kos = len(num_ko_hits_in_path)
        naive_redundancy = total_extra_hits/num_kos
        if all(e == 0 for e in extra_hits):
            return 0.0
        entropy = stats.entropy(extra_hits)
        max_entropy_distribution = [total_extra_hits // num_kos] * num_kos
        for i in range(total_extra_hits % num_kos):
            max_entropy_distribution[i] += 1
        max_entropy = stats.entropy(max_entropy_distribution)
        # avoid divide by 0
        max_entropy += 1e-20

        return naive_redundancy * entropy/max_entropy


    def compute_num_complete_copies_of_path(self, copy_num_of_atomic_steps):
        """This function computes the number of copies of a path that are >= x% complete,
        where x is the module completeness threshold.

        It does this based on the provided list in which each entry is the number of copies of
        each atomic step in the path.
        - first, these copy numbers are ordered (descending order)
        - then, we compute N, the number of steps needed to make the path at least X complete, where
          X is the module completeness threshold
        - finally, we loop from i=1 to the maximum number of hits. Each time, if the number x of steps
          with hit count >= i is x >= N, we add 1 to our count of path copy numbers.
        - the final count of path copy numbers is returned.

        PARAMETERS
        ==========
        copy_num_of_atomic_steps : list
            stores the number of copies of each step in path

        RETURNS
        ==========
        copy_number : int
            number of copies of path which are at least X complete, where X is module completeness threshold
        """

        import math
        path_length = len(copy_num_of_atomic_steps)
        num_enzymes_needed = math.ceil(self.module_completion_threshold * path_length)  # N
        copy_num_of_atomic_steps.sort(reverse=True)

        copy_number = 0
        for i in range(1, copy_num_of_atomic_steps[0]+1):
            x = len([h for h in copy_num_of_atomic_steps if h >= i])
            if x >= num_enzymes_needed:
                copy_number += 1

        return copy_number


    def compute_module_redundancy_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the redundancy of the specified module within the given bin metabolism dictionary.

        Each module can have multiple paths, but (in most cases) we only compute redundancy on the paths with the highest completeness
        (stored under the "most_complete_paths" key). If there are no paths in this list (which only happens when there
        are 0 KOfam hits to the module), then we do not compute redundancy.

        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        """

        meta_dict_for_bin[mnum]["naive_redundancy"] = []
        meta_dict_for_bin[mnum]["copywise_average"] = []
        meta_dict_for_bin[mnum]["copywise_completeness_distributions"] = []
        meta_dict_for_bin[mnum]["copywise_median"] = []
        meta_dict_for_bin[mnum]["copywise_weighted-sum"] = []
        meta_dict_for_bin[mnum]["copywise_geometric-mean"] = []
        meta_dict_for_bin[mnum]["entropy_weighted"] = []

        paths_of_highest_completeness = meta_dict_for_bin[mnum]["most_complete_paths"]
        if not paths_of_highest_completeness:
            return

        for p in paths_of_highest_completeness:
            p = self.split_module_path_into_individual_essential_components(p)
            num_hits_per_kofam = [len(meta_dict_for_bin[mnum]["kofam_hits"][k]) if k in meta_dict_for_bin[mnum]["kofam_hits"] else 0 for k in p]

            # for now, we will try a bunch of different redundancy calculations and put them all into the dictionary until we find the ones we like
            meta_dict_for_bin[mnum]["naive_redundancy"].append(self.compute_naive_redundancy_for_path(num_hits_per_kofam))
            cw_avg_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="average")
            meta_dict_for_bin[mnum]["copywise_average"].append(cw_avg_redundancy)
            meta_dict_for_bin[mnum]["copywise_completeness_distributions"].append(copy_completeness_distribution)
            cw_med_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="median")
            meta_dict_for_bin[mnum]["copywise_median"].append(cw_med_redundancy)
            cw_ws_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="weighted_sum")
            meta_dict_for_bin[mnum]["copywise_weighted-sum"].append(cw_ws_redundancy)
            cw_gm_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="geometric_mean")
            meta_dict_for_bin[mnum]["copywise_geometric-mean"].append(cw_gm_redundancy)
            meta_dict_for_bin[mnum]["entropy_weighted"].append(self.compute_entropy_weighted_redundancy_for_bin(num_hits_per_kofam))

        return


    def get_step_copy_number(self, step_string, enzyme_hit_counts):
        """This function recursively calculates the copy number of a step in a module.

        It parses the definition string of a step and recurses as needed to compute copy number of
        substeps. Copy numbers of any substeps are mathematically combined to obtain a copy number for
        the step as a whole.

        The key base case in the recursion is an individual enzyme accession, for which the copy number
        is simply the number of times it is annotated in the sample (which we obtain from the enzyme_hit_counts dictionary).

        Combining copy numbers works as follows: If two enzymes (or substeps) are connected by an AND, then we need both, so
        we take the minimum copy number of both of them. If they are connected by an OR, then we can use either, so we can sum
        their copy numbers. In doing this, we follow correct order of operations, as established by any parentheses in the step definition.

        In short, this function accomplishes the same thing as modifying the step definition by replacing spaces and '+' signs with min()
        operations, replacing commas with + operations, and replacing enzyme accessions with their corresponding hit counts; then returning
        the value obtained by evaluating the resulting arithmetic expression.

        Some steps are defined by other modules. When module accessions are found, we initially treat them as having a copy number of 0, but
        we re-compute the copy number of the module later once we have the overall copy number of all other modules (and then we use the
        component module's copy number in the calculation instead).

        PARAMETERS
        ==========
        step_string : str
            A string containing the definition of one step (or substep) from a module
        enzyme_hit_counts : dict
            Keys are enzyme accessions, values are the number of times the enzyme was annotated in the current sample

        RETURNS
        ==========
        The copy number (int) of the given step/substep
        """

        # first, eliminate non-essential KOs from the step definition so they won't be considered
        step_string = self.remove_nonessential_enzymes_from_module_step(step_string)

        # sometimes a step will have commas outside parentheses. In this case, we need to split those first for proper order of operations
        step_list = utils.split_by_delim_not_within_parens(step_string, ",")
        if len(step_list) > 1:
            added_step_count = 0
            # call recursively on each split
            for s in step_list:
                added_step_count += self.get_step_copy_number(s, enzyme_hit_counts)
            # combine results using addition and return
            return added_step_count

        # complex case - parentheses surround substeps, which need to be counted recursively and appropriately combined
        if '(' in step_string:
            open_parens_idx = step_string.index('(') # first (outermost) parenthesis
            close_parens_idx = None

            # find matching parenthesis
            i = open_parens_idx + 1
            parens_level = 1
            while not close_parens_idx:
                if step_string[i] == '(':
                    parens_level += 1
                if step_string[i] == ')':
                    parens_level -= 1

                    if parens_level == 0:
                        close_parens_idx = i
                i += 1

            # call recursively on string within outermost parentheses
            sub_step = step_string[open_parens_idx+1:close_parens_idx]
            sub_copy_num = self.get_step_copy_number(sub_step, enzyme_hit_counts)

            # parse the rest of the string and combine with the copy number of the stuff within parentheses
            step_copy_num = None
            # handle anything prior to parentheses
            if open_parens_idx > 0:
                previous_str = step_string[:open_parens_idx]

                previous_steps = previous_str[:-1]
                prev_copy = self.get_step_copy_number(previous_steps, enzyme_hit_counts)

                combo_element = previous_str[-1]
                if combo_element == ',': # OR
                    step_copy_num = (prev_copy + sub_copy_num)
                if combo_element == ' ' or combo_element == '+': # AND
                    step_copy_num = min(prev_copy,sub_copy_num)

            # handle anything following parentheses
            if close_parens_idx < len(step_string) - 1:
                post_steps = step_string[close_parens_idx+2:]
                post_copy = self.get_step_copy_number(post_steps, enzyme_hit_counts)

                combo_element = step_string[close_parens_idx+1]
                if step_copy_num is None:
                    # no previous clause, so we only combine the parenthetical clause and what comes after
                    if combo_element == ',': # OR
                        step_copy_num = (sub_copy_num + post_copy)
                    if combo_element == ' ' or combo_element == '+': # AND
                        step_copy_num = min(sub_copy_num,post_copy)
                else:
                    # we have to combine the post clause with the already-combined previous clause
                    # and parenthetical clause
                    if combo_element == ',': # OR
                        step_copy_num += post_copy
                    if combo_element == ' ' or combo_element == '+': # AND
                        step_copy_num = min(step_copy_num,post_copy)

            # handle edge case where parentheses circles entire step
            if (open_parens_idx == 0) and (close_parens_idx == len(step_string) - 1):
                step_copy_num = sub_copy_num

            return step_copy_num

        # simple case - no substeps within parentheses
        else:
            if ',' in step_string: # OR - combine copy numbers using addition
                or_splits = step_string.split(',')
                added_step_count = 0
                for s in or_splits:
                    added_step_count += self.get_step_copy_number(s, enzyme_hit_counts)
                return added_step_count

            elif ' ' in step_string or '+' in step_string: # AND - combine copy numbers using min()
                and_splits = step_string.replace('+', ' ').split(' ')
                min_step_count = None
                for s in and_splits:
                    s_count = self.get_step_copy_number(s, enzyme_hit_counts)
                    if min_step_count is None:
                        min_step_count = s_count # make first step the minimum
                    min_step_count = min(min_step_count, s_count)
                return min_step_count

            # base cases
            elif '-' in step_string:
                if step_string == '--': # no KO profile => no copy number (unless user wants to ignore these, in which case
                        return 0        # they were already removed by remove_nonessential_enzymes_from_module_step() above)
                else: # contains non-essential KO, should never happen because we eliminated them above
                    raise ConfigError(f"Something is very wrong, because the get_step_copy_number() function found a nonessential "
                                      f"enzyme in the step definition {step_string}")
            elif step_string == '': # entire step was nonessential KO, skip computation
                return None
            else: # accession
                if step_string in self.all_modules_in_db: # module
                    if step_string in enzyme_hit_counts: # we are currently adjusting, and know the module copy number
                        return enzyme_hit_counts[step_string]
                    else: # return 0 for now, will be adjusted later
                        return 0
                else: # enzyme
                    if step_string not in enzyme_hit_counts:
                        return 0
                    return enzyme_hit_counts[step_string]


    def are_enzymes_indirect_alternatives_within_step(self, enzyme_list: list, step: str):
        """An overly simplistic function to determine whether the relationship between the provided alternative
        enzymes in the given step is indirect.

        To do this, it simply walks through the step definition string to determine whether each pair of enzymes is separated by
        a character symbolizing a more complex relationship. That is, they are not separated only by commas and other enzymes (which
        indicates a direct relationship, as in the two enzymes are synonymous in the context of the metabolic pathway).

        For example, within the step (((K01657+K01658,K13503,K13501,K01656) K00766),K13497), the direct alternatives include
        K13503, K13501, and K01656. K01657 and K01658 are indirect alternatives to each other because they are two
        components of the same enzyme, while K01658 and K00766 are indirect because they catalyze two separate reactions in
        an alternative branch of the step.

        This algorithm is not perfect at identifying all indirect relationships - for instance, given K01658 and K13503 it will
        wrongly suggest they are direct alternatives. However, it is meant to be used only for identifying putative edge cases
        for the `get_dereplicated_enzyme_hits_for_step_in_module()` function, and it works well enough for that.

        PARAMETERS
        ==========
        enzyme_list : list of enzyme accessions
            the alternative enzymes to process
        step : string
            the definition string of the relevant step

        RETURNS
        =======
        contains_indirect : Boolean
            True if the list of provided enzymes contains those that are indirect alternatives within the given step.
        """

        enzyme_data = {e : {'index': step.index(e),
                                    'direct_alts': [],
                                    'indirect_alts': []} for e in enzyme_list}

        contains_indirect = False
        # get enzyme-specific list of alternatives
        for e in enzyme_list:
            for z in enzyme_list:
                if e != z:
                    e_index = enzyme_data[e]['index']
                    z_index = enzyme_data[z]['index']
                    indirect_alternatives = False

                    # indirect alts have a space, parentheses, or plus/minus sign between them
                    for c in step[min(e_index, z_index):max(e_index, z_index)]:
                        if c in [' ', '(', ')', '+', '-']:
                            indirect_alternatives = True

                    if indirect_alternatives:
                        enzyme_data[e]['indirect_alts'].append(z)
                        contains_indirect = True
                    else:
                        enzyme_data[e]['direct_alts'].append(z)

        return contains_indirect


    def get_dereplicated_enzyme_hits_for_step_in_module(self, meta_dict_for_mnum: dict, step_to_focus_on: str, mnum: str):
        """This function returns a dictionary of enzyme accessions matched to the number of hits, with duplicate hits to the
        same gene removed, for the provided step in a metabolic pathway.

        Depreplicating the gene calls is necessary because the same gene can be annotated with multiple alternative enzymes for the
        same reaction, and we don't want these annotations to be double-counted in the stepwise copy number calculation.

        PARAMETERS
        ==========
        meta_dict_for_mnum : dictionary of dictionaries
            metabolism completeness dict for the current bin and metabolic module
        step_to_focus_on : string
            which step in the module to resolve alternative enzymes for, passed as a definition string for the step.
        mnum : string
            module ID (used only for warning output)

        RETURNS
        =======
        derep_enzyme_hits : dictionary
            matches enzyme accession to number of hits to unique genes
        """

        derep_enzyme_hits = {k : len(meta_dict_for_mnum["kofam_hits"][k]) for k in meta_dict_for_mnum["kofam_hits"] if k in step_to_focus_on}

        # map gene caller IDs to enzyme accessions
        gene_calls_to_enzymes = {gcid : [] for gcid in meta_dict_for_mnum['gene_caller_ids']}
        for enzyme, gene_list in meta_dict_for_mnum['kofam_hits'].items():
            for g in gene_list:
                if enzyme in step_to_focus_on:
                    gene_calls_to_enzymes[g].append(enzyme)

        for gcid, enzymes in gene_calls_to_enzymes.items():
            if len(enzymes) > 1:
                # simple solution (only works well for enzymes that are direct alternatives)
                # for each duplicated gene, we arbitrarily keep only the hit to the first enzyme
                # and for all other annotations, we reduce the count of hits by one
                for acc in enzymes[1:]:
                    derep_enzyme_hits[acc] -= 1

                if self.are_enzymes_indirect_alternatives_within_step(enzymes, step_to_focus_on) and self.add_copy_number:
                    enz_str = ", ".join(enzymes)
                    self.run.warning(f"The gene call {gcid} has multiple annotations to alternative enzymes "
                                     f"within the same step of a metabolic pathway ({enz_str}), and these enzymes "
                                     f"unfortunately have a complex relationship. The affected module is {mnum}, and "
                                     f"here is the step in question: {step_to_focus_on}. We arbitrarily kept only one of "
                                     f"the annotations to this gene in order to avoid inflating the step's copy number, "
                                     f"but due to the complex relationship between these alternatives, this could mean "
                                     f"that the copy number for this step is actually too low. Please heed this warning "
                                     f"and double check the stepwise copy number results for {mnum} and other pathways "
                                     f"containing gene call {gcid}.")

        return derep_enzyme_hits


    def compute_stepwise_module_copy_number_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the copy number of the specified module within the given bin metabolism dictionary.

        It goes through the top-level steps established by compute_stepwise_module_completeness_for_bin() and determines the
        copy number of each step. Then, the overall module copy number is calculated as the minimum copy number of all steps.

        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "stepwise_copy_number"         the stepwise copy number of the module

        [keys added in "top_level_step_info" dictionary]
            "copy_number"              the copy number of an individual step
        """

        all_step_copy_nums = []
        for key in meta_dict_for_bin[mnum]["top_level_step_info"]:
            if not meta_dict_for_bin[mnum]["top_level_step_info"][key]["includes_modules"]:
                step_string = meta_dict_for_bin[mnum]["top_level_step_info"][key]["step_definition"]
                enzyme_hits_dict = self.get_dereplicated_enzyme_hits_for_step_in_module(meta_dict_for_bin[mnum], step_string, mnum)

                step_copy_num = self.get_step_copy_number(step_string, enzyme_hits_dict)
                meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"] = step_copy_num
                if step_copy_num is not None: # avoid taking minimum of None values (from non-essential steps)
                    all_step_copy_nums.append(step_copy_num)

        if all_step_copy_nums:
            module_stepwise_copy_num = min(all_step_copy_nums)
        else:
            module_stepwise_copy_num = None
        meta_dict_for_bin[mnum]["stepwise_copy_number"] = module_stepwise_copy_num


    def adjust_stepwise_copy_number_for_bin(self, mnum, meta_dict_for_bin):
        """This function adjusts stepwise copy number of modules that are defined by other modules.

        This can only be done after all other modules have had their copy numbers calculated and added to the metabolism dictionary
        by the function compute_stepwise_module_copy_number_for_bin().

        The function goes through the top-level steps in the module and re-computes copy number for steps that include other modules.
        Then it re-calculates the overall module copy number as the minimum copy number of all steps. It updates the metabolism completess
        dictionary accordingly.

        PARAMETERS
        ==========
        mnum : string
            the module number to adjust
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dictionary for the current bin
        """

        enzyme_hits_dict = {k : len(meta_dict_for_bin[mnum]["kofam_hits"][k]) for k in meta_dict_for_bin[mnum]["kofam_hits"] }

        all_step_copy_nums = []
        for key in meta_dict_for_bin[mnum]["top_level_step_info"]:
            # re-calculate ONLY for steps with modules in definition
            if meta_dict_for_bin[mnum]["top_level_step_info"][key]["includes_modules"]:
                step_string = meta_dict_for_bin[mnum]["top_level_step_info"][key]["step_definition"]

                for included_module in meta_dict_for_bin[mnum]["top_level_step_info"][key]["included_module_list"]:
                    enzyme_hits_dict[included_module] = meta_dict_for_bin[included_module]["stepwise_copy_number"]

                step_copy_num = self.get_step_copy_number(step_string, enzyme_hits_dict)
                meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"] = step_copy_num

            # take minimum over all steps, even those not defined by modules
            if meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"] is not None: # avoid taking minimum of None values (from non-essential steps)
                all_step_copy_nums.append(meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"])

        module_stepwise_copy_num = min(all_step_copy_nums)
        meta_dict_for_bin[mnum]["stepwise_copy_number"] = module_stepwise_copy_num


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
        metabolism_dict_for_genome, ko_dict_for_genome = self.mark_kos_present_for_list_of_splits(kofam_gene_split_contig, split_list=splits_in_genome,
                                                                                                    bin_name=self.contigs_db_project_name)
        if not self.store_json_without_estimation:
            genome_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(metabolism_dict_for_genome, bin_name=self.contigs_db_project_name)
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

            metabolism_dict_for_bin, ko_dict_for_bin = self.mark_kos_present_for_list_of_splits(ko_in_bin, split_list=splits_in_bin, bin_name=bin_name)

            if not self.store_json_without_estimation:
                bins_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(metabolism_dict_for_bin, bin_name=bin_name)
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
            metabolism_dict_for_contig, ko_dict_for_contig = self.mark_kos_present_for_list_of_splits(ko_in_contig, split_list=splits_in_contig, bin_name=contig)


            if not self.store_json_without_estimation:
                single_contig_module_superdict = {contig: self.estimate_for_list_of_splits(metabolism_dict_for_contig, bin_name=contig)}
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

            new_kegg_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(meta_dict_for_bin, bin_name=bin_name)
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
                                                                                                 bin_name=self.contigs_db_project_name)
        if not self.store_json_without_estimation:
            enzyme_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(metabolism_dict_for_genome,
                                                                                                         bin_name=self.contigs_db_project_name)
            enzyme_ko_superdict[self.contigs_db_project_name] = ko_dict_for_genome
        else:
            enzyme_metabolism_superdict[self.contigs_db_project_name] = metabolism_dict_for_genome
            enzyme_ko_superdict[self.contigs_db_project_name] = ko_dict_for_genome

        # append to file
        self.append_kegg_metabolism_superdicts(enzyme_metabolism_superdict, enzyme_ko_superdict)

        return enzyme_metabolism_superdict, enzyme_ko_superdict


    def init_hits_for_pangenome(self, gene_cluster_list: list):
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
        pan_super.init_gene_clusters_functions_summary_dict(source_list = self.annotation_sources_to_use, gene_clusters_of_interest = gene_cluster_list)


        enzyme_cluster_split_contig = []
        # no splits or contigs here
        for cluster_id in gene_cluster_list:
            for source in self.annotation_sources_to_use:
                if source in pan_super.gene_clusters_functions_summary_dict[cluster_id]:
                    acc = pan_super.gene_clusters_functions_summary_dict[cluster_id][source]['accession']
                    if acc: # avoid introducing 'None' values here
                        enzyme_cluster_split_contig.append((acc,cluster_id,"NA","NA"))

        return enzyme_cluster_split_contig


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
            metabolism_dict_for_bin, ko_dict_for_bin = self.mark_kos_present_for_list_of_splits(enzymes_in_bin, bin_name=bin_name)

            if not self.store_json_without_estimation:
                gc_bins_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(metabolism_dict_for_bin, bin_name=bin_name)
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

                kofam_hits_info = self.init_hits_for_pangenome(gene_cluster_list = all_gene_clusters_in_collection)
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_metabolism_for_pangenome_bins(kofam_hits_info, collection_dict)
            else:
                kofam_hits_info = self.init_hits_and_splits(annotation_sources=self.annotation_sources_to_use)

                if self.add_coverage:
                    self.init_gene_coverage(gcids_for_kofam_hits={int(tpl[1]) for tpl in kofam_hits_info})

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
            self.store_metabolism_superdict_as_json(kegg_metabolism_superdict, self.json_output_file_path + ".json")

        # more housekeeping for outputs
        if not self.multi_mode:
            for mode, file_object in self.output_file_dict.items():
                file_object.close()

        # at this point, if we are generating long-format output, the data has already been appended to files
        # so we needn't keep it in memory. We don't return it, unless the programmer wants us to.
        if return_superdicts:
            return kegg_metabolism_superdict, kofam_hits_superdict
        # on the other hand, if we are generating matrix output, we need a limited subset of this data downstream
        # so in this case, we can extract and return smaller dictionaries for module completeness, module presence/absence,
        # and KO hits.
        elif return_subset_for_matrix_format:
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


class KeggMetabolismEstimatorMulti(KeggContext, KeggEstimatorArgs):
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
        KeggContext.__init__(self, args)
        KeggEstimatorArgs.__init__(self, self.args)

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
                kos_in_mod = self.get_enzymes_from_module_definition_in_order(mod_def)
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

