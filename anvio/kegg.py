#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains Kegg related classes."""

import os
import shutil
import glob
import re
import copy
import statistics
import json
import time
import hashlib
import collections
import pandas as pd
import numpy as np
import multiprocessing as mp

from scipy import stats
from typing import Dict, List, Tuple, Union

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.tables as t
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.drivers.muscle import Muscle
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions
from anvio.dbops import ContigsSuperclass, ContigsDatabase, ProfileSuperclass, ProfileDatabase, PanSuperclass
from anvio.genomedescriptions import MetagenomeDescriptions, GenomeDescriptions
from anvio.dbinfo import DBInfo


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Iva Veseli"
__email__ = "iveseli@uchicago.edu"


run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print
P = terminal.pluralize


"""Some critical constants for metabolism estimation output formatting.

The below dictionary defines the possible output modes.
- 'output_suffix' should be unique to a mode so that multiple output modes can be used at once
- 'data_dict' indicates which data dictionary is used for generating the output (modules or kofams)
- 'headers' list describes which information to include in the output file (see OUTPUT_HEADERS dict below for more info)
- 'description' is what is printed when --list-available-modes parameter is used
"""
OUTPUT_MODES = {'modules': {
                    'output_suffix': "modules.txt",
                    'data_dict': "modules",
                    'headers': ["module", "module_name", "module_class", "module_category",
                                "module_subcategory", "module_definition",
                                "stepwise_module_completeness", "stepwise_module_is_complete",
                                "pathwise_module_completeness", "pathwise_module_is_complete",
                                "proportion_unique_enzymes_present", "enzymes_unique_to_module", "unique_enzymes_hit_counts",
                                "enzyme_hits_in_module", "gene_caller_ids_in_module", "warnings"],
                    'description': "Information on metabolic modules"
                    },
                'modules_custom': {
                    'output_suffix': "modules_custom.txt",
                    'data_dict': "modules",
                    'headers': None,
                    'description': "A custom tab-delimited output file where you choose the included modules data using --custom-output-headers"
                    },
                'module_paths': {
                                    'output_suffix': "module_paths.txt",
                                    'data_dict': "modules",
                                    'headers': ["module", "pathwise_module_completeness", "pathwise_module_is_complete",
                                                "path_id", "path", "path_completeness", "annotated_enzymes_in_path"],
                                    'description': "Information on each possible path (complete set of enzymes) in a module"
                                    },
                'module_steps': {
                                    'output_suffix': "module_steps.txt",
                                    'data_dict': "modules",
                                    'headers': ["module", "stepwise_module_completeness", "stepwise_module_is_complete",
                                                "step_id", "step", "step_completeness"],
                                    'description': "Information on each top-level step in a module"
                                    },
                'hits': {
                    'output_suffix': "hits.txt",
                    'data_dict': "kofams",
                    'headers': ["enzyme", "gene_caller_id", "contig", "modules_with_enzyme", "enzyme_definition"],
                    'description': "Information on all enzyme annotations in the contigs DB, regardless of module membership"
                    },
                }
"""
The below dictionary describes the type of information we can output
- the dictionary key corresponds to the header's key in the output dictionary (ie, as returned from generate_output_dict_for_modules() function)
- 'cdict_key' is the header's key in modules or kofams data dictionary (if any)
- 'mode_type' indicates which category of output modes (modules or kofams) this header can be used for. If both, this is 'all'
- 'description' is printed when --list-available-output-headers parameter is used
"""
OUTPUT_HEADERS = {'module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Module number"
                        },
                  'stepwise_module_is_complete' : {
                        'cdict_key': "stepwise_is_complete",
                        'mode_type': 'modules',
                        'description': "Whether a module is considered complete or not based on its STEPWISE percent completeness and the completeness threshold"
                        },
                  'stepwise_module_completeness' : {
                        'cdict_key': 'stepwise_completeness',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a module, computed as the number of complete steps divided by the number of total steps "
                                       "(where 'steps' are determined by splitting the module definition on the space character)"
                        },
                  'pathwise_module_is_complete' : {
                        'cdict_key': "pathwise_is_complete",
                        'mode_type': 'modules',
                        'description': "Whether a module is considered complete or not based on its PATHWISE percent completeness and the completeness threshold"
                        },
                  'pathwise_module_completeness' : {
                        'cdict_key': 'pathwise_percent_complete',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a module, computed as maximum completeness of all possible combinations of enzymes ('paths') in the definition"
                        },
                  'enzymes_unique_to_module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A list of enzymes that only belong to this module (ie, they are not members of multiple modules)"
                        },
                  'unique_enzymes_hit_counts' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "How many times each unique enzyme appears in the sample (order of counts corresponds to list in `enzymes_unique_to_module` field)"
                        },
                  'proportion_unique_enzymes_present' : {
                        'cdict_key': 'proportion_unique_enzymes_present',
                        'mode_type': 'modules',
                        'description': "Proportion of enzymes unique to this one module that are present in the sample"
                        },
                  'unique_enzymes_context_string' : {
                        'cdict_key': 'unique_enzymes_context_string',
                        'mode_type': 'modules',
                        'description': "Describes the unique enzymes contributing to the `proportion_unique_enzymes_present` field"
                        },
                  'module_name' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Name/description of a module"
                        },
                  'module_class' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism class of a module"
                        },
                  'module_category' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism category of a module"
                        },
                  'module_subcategory' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism subcategory of a module"
                        },
                  'module_definition' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Definition string of a module. Describes the metabolic pathway "
                                       "in terms of the enzymes (KOs, COGs, etc) that belong to the module."
                        },
                  'module_substrates' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of compounds that serve as initial input to the metabolic pathway "
                                       "(that is, substrate(s) to the initial reaction(s) in the pathway)"
                        },
                  'module_products' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of compounds that serve as final output from the metabolic pathway "
                                       "(that is, product(s) of the final reaction(s) in the pathway)"
                        },
                  'module_intermediates' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of compounds that are intermediates in the metabolic pathway "
                                       "(compounds that are both outputs and inputs of reaction(s) in the pathway)"
                        },
                  'gene_caller_ids_in_module': {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of gene caller IDs of enzymes that contribute to a module"
                        },
                  'gene_caller_id': {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': "Gene caller ID of a single enzyme in the contigs DB"
                        },
                  'enzyme_hits_in_module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of enzyme annotations that contribute to a module"
                        },
                  'enzyme_hit' : {
                        'cdict_key': 'kofam_hits',
                        'mode_type': 'kofams',
                        'description': "Enzyme identifier for a single annotation (KO, COG, etc)"
                        },
                  'contig' : {
                        'cdict_key': 'genes_to_contigs',
                        'mode_type': 'kofams',
                        'description': "Contig that an enzyme annotation is found on"
                        },
                  'path_id' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Integer ID for a path through a module. Has no real meaning and is used for data organization"
                        },
                  'path' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A path through a module: a linear sequence of enzymes that together represent each metabolic step "
                                       "in the module (most modules have several of these due to enzyme redundancy)"
                        },
                  'path_completeness' : {
                        'cdict_key': 'pathway_completeness',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a given path through a module"
                        },
                  'annotated_enzymes_in_path' : {
                        'cdict_key': 'annotated_enzymes_in_path',
                        'mode_type': 'modules',
                        'description': "Shows which enzymes in the path are annotated in your sample, and which are missing"
                        },
                  'step_id' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Integer ID for a top-level step in a module. Has no real meaning and is used for data organization"
                        },
                  'step' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A 'top-level' step in a module, represented by one or more possible enzymes that can catalyze "
                                       "a logical part of the metabolic pathway (usually one reaction)"
                        },
                  'step_completeness' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Binary completeness of a given 'top-level' step in a module"
                        },
                  'warnings' : {
                        'cdict_key': 'warnings',
                        'mode_type': 'modules',
                        'description': "This column holds a comma-separated list of notes about things that might affect completeness "
                                       "estimates for a module, such as missing enzyme profiles."
                        },
                  'enzyme' : {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': 'Identifier for an enzyme that is annotated in your database(s), ie a KO or COG number'
                        },
                  'modules_with_enzyme': {
                        'cdict_key': 'modules',
                        'mode_type': 'kofams',
                        'description': 'A comma-separated list of modules that the enzyme belongs to'
                        },
                  'enzyme_definition': {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': 'The functional annotation associated with the enzyme'
                        },
                  }

DEFAULT_OUTPUT_MODE = 'modules'
STRAY_KO_ANVIO_SUFFIX = "_anvio_version"

# global metadata header lists for matrix format
# if you want to add something here, don't forget to add it to the dictionary in the corresponding
# get_XXX_metadata_dictionary() function
MODULE_METADATA_HEADERS = ["module_name", "module_class", "module_category", "module_subcategory"]
KO_METADATA_HEADERS = ["enzyme_definition", "modules_with_enzyme"]
# Exception: if you add to this list, you must add it in the steps_subdict in generate_subsets_for_matrix_format()
# and to the relevant step metadata clause in write_stat_to_matrix()
STEP_METADATA_HEADERS = ["step_definition"]

# Global and overview map IDs have certain ranges of numbers.
GLOBAL_MAP_ID_PATTERN = re.compile(r'\d{1}11\d{2}')
OVERVIEW_MAP_ID_PATTERN = re.compile(r'\d{1}12\d{2}')


class KeggContext(object):
    """The purpose of this base class is to define shared functions and file paths for all KEGG operations."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # default data directory will be called KEGG and will store the KEGG Module data as well
        self.default_kegg_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')
        self.kegg_data_dir = os.path.abspath(A('kegg_data_dir') or self.default_kegg_dir)
        self.user_input_dir = os.path.abspath(A('user_modules')) if A('user_modules') else None
        self.only_user_modules = A('only_user_modules')
        self.orphan_data_dir = os.path.join(self.kegg_data_dir, "orphan_data")
        self.kegg_module_data_dir = os.path.join(self.kegg_data_dir, "modules")
        self.kegg_hmm_data_dir = os.path.join(self.kegg_data_dir, "HMMs")
        self.pathway_data_dir = os.path.join(self.kegg_data_dir, "pathways")
        self.brite_data_dir = os.path.join(self.kegg_data_dir, "BRITE")
        self.binary_relation_data_dir = os.path.join(self.kegg_data_dir, "binary_relations")

        # The 'KEGG/map_images' directory has a structure of nested directories. 'map_images'
        # contains 'png' for image files and 'kgml' for XML mapping files. Within both 'png' and
        # 'kgml' are directories, '1x' and '2x', for lower and higher resolution maps. 'png/1x'
        # contains 5 directories of image files highlighting different things: 'map', 'ko', 'ec',
        # 'rn', and 'org'. 'png/2x' contains 1 directory, 'map', as higher resolution images are
        # only available for manually drawn maps. 'kgml/1x' and 'kgml/2x' each contain 4 directories
        # of XML files that allow modification of different lower and higher resolution maps: 'ko',
        # 'ec', 'rn', and 'org'.
        self.map_image_data_dir = os.path.join(self.kegg_data_dir, "map_images")
        self.png_dir = os.path.join(self.map_image_data_dir, "png")
        self.kgml_dir = os.path.join(self.map_image_data_dir, "kgml")
        self.png_1x_dir = os.path.join(self.png_dir, "1x")
        self.png_2x_dir = os.path.join(self.png_dir, "2x")
        self.png_1x_map_dir = os.path.join(self.png_1x_dir, "map")
        self.png_1x_ko_dir = os.path.join(self.png_1x_dir, "ko")
        self.png_1x_ec_dir = os.path.join(self.png_1x_dir, "ec")
        self.png_1x_rn_dir = os.path.join(self.png_1x_dir, "rn")
        self.png_1x_org_dir = os.path.join(self.png_1x_dir, "org")
        self.png_2x_map_dir = os.path.join(self.png_2x_dir, "map")
        self.kgml_1x_dir = os.path.join(self.kgml_dir, "1x")
        self.kgml_2x_dir = os.path.join(self.kgml_dir, "2x")
        self.kgml_1x_ko_dir = os.path.join(self.kgml_1x_dir, "ko")
        self.kgml_1x_ec_dir = os.path.join(self.kgml_1x_dir, "ec")
        self.kgml_1x_rn_dir = os.path.join(self.kgml_1x_dir, "rn")
        self.kgml_1x_org_dir = os.path.join(self.kgml_1x_dir, "org")
        self.kgml_2x_ko_dir = os.path.join(self.kgml_2x_dir, "ko")
        self.kgml_2x_ec_dir = os.path.join(self.kgml_2x_dir, "ec")
        self.kgml_2x_rn_dir = os.path.join(self.kgml_2x_dir, "rn")
        self.kgml_2x_org_dir = os.path.join(self.kgml_2x_dir, "org")

        self.quiet = A('quiet') or False
        self.just_do_it = A('just_do_it')

        # shared variables for all KEGG subclasses
        self.kofam_hmm_file_path = os.path.join(self.kegg_hmm_data_dir, "Kofam.hmm") # file containing concatenated KOfam hmms
        self.stray_ko_hmm_file_path = os.path.join(self.orphan_data_dir, "anvio_hmm_profiles_for_stray_KOs.hmm") # anvi'o-generated concatenated hmms for stray KOs
        self.stray_ko_hmms_from_kegg = os.path.join(self.orphan_data_dir, "hmm_profiles_with_kofams_with_no_threshold.hmm") # original concatentated hmms for stray KOs
        self.ko_list_file_path = os.path.join(self.kegg_data_dir, "ko_list.txt")
        self.stray_ko_thresholds_file = os.path.join(self.orphan_data_dir, "estimated_thresholds_for_stray_kos.txt")
        self.kegg_module_file = os.path.join(self.kegg_data_dir, "modules.keg")
        self.kegg_pathway_file = os.path.join(self.kegg_data_dir, "pathways.keg")
        self.kegg_brite_hierarchies_file = os.path.join(self.kegg_data_dir, "hierarchies.json")
        self.kegg_modules_db_path = os.path.join(self.kegg_data_dir, "MODULES.db")
        self.kegg_binary_relation_files = {('KO', 'EC'): "ko2ec.xl", ('KO', 'RN'): "ko2rn.xl"}
        self.kegg_pathway_list_file = os.path.join(self.kegg_data_dir, "pathway_list.tsv")
        self.kegg_map_image_kgml_file = os.path.join(self.kegg_data_dir, "map_kgml.tsv")

        if self.user_input_dir:
            self.user_module_data_dir = os.path.join(self.user_input_dir, "modules")
            self.user_modules_db_path = os.path.join(self.user_input_dir, "USER_MODULES.db")

        # sanity check for incompatible arguments
        if A('kegg_data_dir') and A('only_user_modules'):
            raise ConfigError("The options --kegg-data-dir and --only-user-modules are incompatible. Please figure out which one you "
                              "want and try again :)")

        # sanity check to prevent automatic overwriting of non-default kegg data dir
        if self.__class__.__name__ in ['KeggSetup'] and not self.user_input_dir:
            if os.path.exists(self.kegg_data_dir) and self.kegg_data_dir != self.default_kegg_dir:
                raise ConfigError(f"You are attempting to set up KEGG in a non-default data directory ({self.kegg_data_dir}) which already exists. "
                                  f"To avoid automatically deleting a directory that may be important to you, anvi'o refuses to get rid of "
                                  f"directories that have been specified with --kegg-data-dir. If you really want to get rid of this "
                                  f"directory and replace it with the KEGG archive data, then please remove the directory yourself using "
                                  f"a command like `rm -r {self.kegg_data_dir}`. We are sorry to make you go through this extra trouble, but it really is "
                                  f"the safest way to handle things.")


    def setup_ko_dict(self, exclude_threshold=True, suppress_warnings=False):
        """The purpose of this function is to process the ko_list file into usable form by KEGG sub-classes.

        The ko_list file (which is downloaded along with the KOfam HMM profiles) contains important
        information for each KEGG Orthology number (KO, or knum), incuding pre-defined scoring thresholds
        for limiting HMM hits and annotation information.

        It looks something like this:

        knum    threshold    score_type    profile_type    F-measure    nseq    nseq_used    alen    mlen    eff_nseq    re/pos    definition
        K00001    329.57    domain    trim    0.231663    1473    1069    1798    371    17.12    0.590    alcohol dehydrogenase [EC:1.1.1.1]

        Since this information is useful for both the setup process (we need to know all the knums) and HMM process,
        all KEGG subclasses need to have access to this dictionary.

        This is a dictionary (indexed by knum) of dictionaries(indexed by column name).
        Here is an example of the dictionary structure:
        self.ko_dict["K00001"]["threshold"] = 329.57

        PARAMETERS
        ==========
        exclude_threshold : Boolean
            If this is true, we remove KOs without a bitscore threshold from the ko_dict
        suppress_warnings : Boolean
            If this is true, we don't print the warning message about stray KOs
        """

        self.ko_dict = utils.get_TAB_delimited_file_as_dictionary(self.ko_list_file_path)
        self.ko_skip_list, self.ko_no_threshold_list = self.get_ko_skip_list()

        # if we are currently setting up KOfams, we should generate a text file with the ko_list entries
        # of the KOs that have no scoring threshold
        if self.__class__.__name__ in ['KeggSetup', 'KOfamDownload']:
            stray_ko_dict = {ko:self.ko_dict[ko] for ko in self.ko_skip_list}
            stray_ko_dict.update({ko:self.ko_dict[ko] for ko in self.ko_no_threshold_list})

            if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
                raise ConfigError(f"Hmm. Something is out of order. The orphan data directory {self.orphan_data_dir} does not exist "
                                  f"yet, but it needs to in order for the setup_ko_dict() function to work.")
            stray_ko_path = os.path.join(self.orphan_data_dir, "kofams_with_no_threshold.txt")
            stray_ko_headers = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos", "definition"]
            utils.store_dict_as_TAB_delimited_file(stray_ko_dict, stray_ko_path, key_header="knum", headers=stray_ko_headers)

        [self.ko_dict.pop(ko) for ko in self.ko_skip_list]
        if exclude_threshold:
            [self.ko_dict.pop(ko) for ko in self.ko_no_threshold_list]
        else:
            if not suppress_warnings:
                self.run.warning("FYI, we are including KOfams that do not have a bitscore threshold in the analysis.")


    def setup_stray_ko_dict(self, add_entries_to_regular_ko_dict=False):
        """This class sets up a dictionary of predicted bit score thresholds for stray KOs, if possible.

        Those predicted thresholds are generated during `anvi-setup-kegg-data --include-stray-KOs`
        (see KOfamDownload.process_all_stray_kos()), and are stored in a file that looks like this:

        knum	threshold	score_type	definition
        K11700	800.4	full	poly(A) RNA polymerase Cid12 [EC:2.7.7.19]
        K14747_anvio_version	1054.2	full	benzoylacetate-CoA ligase [EC:6.2.1.-]

        The dictionary structure is identical to that of self.ko_dict. Note that the `knum` column can contain
        normal KEGG Ortholog accessions (for KOs whose HMMs we haven't updated) and accessions that end with
        STRAY_KO_ANVIO_SUFFIX (for KOs that we created new models for).

        If thresholds have not been predicted, then this function throws an error.

        Parameters
        ==========
        add_entries_to_regular_ko_dict : Boolean
            If True, we don't create a separate self.stray_ko_dict but instead add the stray KOs to the
            regular self.ko_dict attribute. Useful if you don't need to keep the two sets separate.
        """

        if os.path.exists(self.stray_ko_thresholds_file):
            if add_entries_to_regular_ko_dict:
                stray_kos = utils.get_TAB_delimited_file_as_dictionary(self.stray_ko_thresholds_file)
                self.ko_dict.update(stray_kos)
            else:
                self.stray_ko_dict = utils.get_TAB_delimited_file_as_dictionary(self.stray_ko_thresholds_file)
        else:
            raise ConfigError(f"You've requested to include stray KO models in your analysis, but we cannot find the "
                              f"estimated bit score thresholds for these models, which can be generated during "
                              f"`anvi-setup-kegg-data` and stored at the following path: {self.stray_ko_thresholds_file}. This "
                              f"means that `anvi-setup-kegg-data` was run without the `--include-stray-KOs` flag for the KEGG "
                              f"data directory that you are using. You have two options: 1) give up on including these models and "
                              f"re-run your command without the `--include-stray-KOs` flag, or 2) change the KEGG data that you "
                              f"are using, which could be as simple as specifying a new `--kegg-data-dir` or as complex as re-running "
                              f"`anvi-setup-kegg-data` to obtain a dataset with predicted thresholds for stray KO models.")


    def get_ko_skip_list(self):
        """The purpose of this function is to determine which KO numbers have no associated data or just no score threshold in the ko_list file.

        That is, their ko_list entries look like this, with hypens in all but the first and last columns:

        K14936    -    -    -    -    -    -    -    -    -    -    small nucleolar RNA snR191
        K15035    -    -    -    -    -    -    -    -    -    -    transfer-messenger RNA
        K15841    -    -    -    -    -    -    -    -    -    -    small regulatory RNA GlmY
        K15851    -    -    -    -    -    -    -    -    -    -    quorum regulatory RNA Qrr
        K16736    -    -    -    -    -    -    -    -    -    -    bantam
        K16863    -    -    -    -    -    -    -    -    -    -    microRNA 21

        These are RNAs.

        Or, their ko_list entries look like this, with no score threshold (but the rest of the data is not completely blank):

        K23749 - - - - 1 1 2266 2266 0.39 0.592 spectinabilin polyketide synthase system NorC [EC:2.3.1.290]

        Returns:
        skip_list  list of strings, each string is a KO number that has no associated data (ie, RNAs)
        no_threshold_list   list of strings, each string is a KO number that has no scoring threshold
        """

        col_names_to_check = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos"]
        skip_list = []
        no_threshold_list = []
        for k in self.ko_dict.keys():
            should_skip = True
            no_threshold = False
            for c in col_names_to_check:
                if not self.ko_dict[k][c] == "-":
                    should_skip = False
                    break # here we stop checking this KO num because we already found a value in our columns of interest

                if c == "threshold":
                    no_threshold = True # if we got to this line of code, there is a '-' in the threshold column
            if should_skip: # should be True unless we found a value above
                skip_list.append(k)
            elif no_threshold:
                no_threshold_list.append(k)
        return skip_list, no_threshold_list


    def invert_brite_json_dict(self, brite_dict):
        """Invert a BRITE hierarchy dict loaded from a json file into a dict keyed by KEGG entries.

        There are only two keys expected in a BRITE json file, 'name' and 'children'. The value for
        'name' is a string and the value for 'children' is a list of dicts.

        Here is an example of what the beginning of the json dict looks like for 'br08902 BRITE
        Hierarchy Files', the 'hierarchy of all existing hierarchies':
           {
             "name": "br08902",
             "children": [
               {
                 "name": "Pathway and Brite",
                 "children": [
                   {
                     "name": "Pathway maps",
                     "children": [
                       {
                         "name": "br08901  KEGG pathway maps"
                       }
                     ]
                   }, ...
        Observe that innermost dicts only have a single entry keyed by 'name'.

        Here is the corresponding entry in the returned dict for the item in the example:
            'br08901  KEGG pathway maps':
                [['br08902', 'Pathway and Brite', 'Pathway maps']]
        The value is a list of lists because an item can occur multiple times in the same hierarchy.

        PARAMETERS
        ==========
        brite_dict : dict
            dict loaded from BRITE hierarchy json file

        RETURNS
        =======
        categorization_dict : dict
            dict of entry categorizations in BRITE hierarchy
        """

        children_stack = collections.deque()
        children_stack.append(([brite_dict['name']], brite_dict['children']))
        categorization_dict = {}
        while children_stack:
            hierarchy, children_list = children_stack.popleft()
            for child_dict in children_list:
                child_name = child_dict['name']
                if 'children' in child_dict:
                    children_stack.append((hierarchy + [child_name], child_dict['children']))
                else:
                    try:
                        categorization_dict[child_name].append(hierarchy)
                    except KeyError:
                        categorization_dict[child_name] = [hierarchy]

        return categorization_dict


class KeggSetup(KeggContext):
    """Class for setting up KEGG Kofam HMM profiles and modules.

    It performs sanity checks and downloads, unpacks, and prepares the profiles for later use by `hmmscan`.
    It also downloads module files and creates the MODULES.db.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-setup-kegg-data. If using this class through the API, please
        provide a Namespace object with the Boolean 'reset' parameter.
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories when testing this class
    """

    def __init__(self, args, run=run, progress=progress, skip_init=False):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress
        self.num_threads = 1 if not A('num_threads') else A('num_threads')
        self.kegg_archive_path = A('kegg_archive')
        self.kegg_snapshot = A('kegg_snapshot')
        self.download_from_kegg = True if A('download_from_kegg') else False
        self.only_download = True if A('only_download') else False
        self.only_processing = True if A('only_processing') else False
        self.skip_init = skip_init
        self.skip_brite_hierarchies = True if A('skip_brite_hierarchies') else False
        self.skip_binary_relations = True if A('skip_binary_relations') else False
        self.skip_map_images = True if A('skip_map_images') else False

        if self.kegg_archive_path and self.download_from_kegg:
            raise ConfigError("You provided two incompatible input options, --kegg-archive and --download-from-kegg. "
                              "Please pick either just one or none of these. ")
        if self.kegg_snapshot and self.download_from_kegg or self.kegg_snapshot and self.kegg_archive_path:
            raise ConfigError("You cannot request setup from an anvi'o KEGG snapshot at the same time as from KEGG directly or from one of your "
                              "KEGG archives. Please pick just one setup option and try again.")

        if (not self.download_from_kegg) and (self.only_download or self.only_processing):
            raise ConfigError("Erm. The --only-download and --only-processing options are only valid if you are also using the --download-from-kegg "
                              "option. Sorry.")
        if self.only_download and self.only_processing:
            raise ConfigError("The --only-download and --only-processing options are incompatible. Please choose only one. Or, if you want both "
                              "download AND database setup to happen, then use only the -D flag without providing either of these two options.")


        # initializing these to None here so that it doesn't break things downstream
        self.pathway_dict = None
        self.brite_dict = None

        # init the base class
        KeggContext.__init__(self, self.args)

        # get KEGG snapshot info for default setup
        self.target_snapshot = self.kegg_snapshot or 'v2024-03-09'
        self.target_snapshot_yaml = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG-SNAPSHOTS.yaml')
        self.snapshot_dict = utils.get_yaml_as_dict(self.target_snapshot_yaml)

        if self.target_snapshot not in self.snapshot_dict.keys():
            self.run.warning(None, header="AVAILABLE KEGG SNAPSHOTS", lc="yellow")
            available_snapshots = sorted(list(self.snapshot_dict.keys()))
            for snapshot_name in available_snapshots:
                self.run.info_single(snapshot_name + (' (latest)' if snapshot_name == available_snapshots[-1] else ''))

            raise ConfigError("Whoops. The KEGG snapshot you requested is not one that is known to anvi'o. Please try again, and "
                                "this time pick from the list shown above.")

        # default download path for KEGG snapshot
        self.default_kegg_data_url = self.snapshot_dict[self.target_snapshot]['url']
        self.default_kegg_archive_file = self.snapshot_dict[self.target_snapshot]['archive_name']

        # the KEGG API URL, in case its needed downstream
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"

        if self.user_input_dir:
            self.run.warning(f"Just so you know, we will be setting up the metabolism data provided at the following "
                             f"location: '{self.user_input_dir}'. The success of this will be determined by how well you "
                             f"followed our formatting guidelines, so keep an eye out for errors below.")


        if not self.user_input_dir:

            # establish parent directory
            if self.download_from_kegg and not self.only_processing and not self.kegg_archive_path and not skip_init:
                filesnpaths.gen_output_directory(self.kegg_data_dir, delete_if_exists=args.reset)

        else: # user input setup
            filesnpaths.is_output_dir_writable(os.path.dirname(self.user_input_dir))

            self.check_user_input_dir_format()

            if not args.reset and not skip_init:
                self.is_user_database_exists()

            if args.reset:
                self.run.warning("Since you used the --reset flag, anvi'o will get rid of any existing user modules database. "
                                 "Now ye be warned.")
                paths_to_remove = [self.user_modules_db_path]
                for path in paths_to_remove:
                    if os.path.exists(path):
                        os.remove(path)
                        self.run.info("Successfully removed", path)


    def is_database_exists(self, files_to_check, fail_if_exists=True):
        """This function determines whether the user has already downloaded all required KEGG data.

        More specifically, it looks for the KEGG files that we use to learn what to download (as in
        the KEGG MODULE file) and for the existence of the data directories that are created by this
        program.

        PARAMETERS
        ==========
        files_to_check : list of file paths
            this list should contain the paths to all required KEGG data or directories. what those
            files are depends on the download mode.
        fail_if_exists : Boolean
            if this is True, this function will fail if the KEGG data already exists on the user's
            computer. If it is False, AND the user has already downloaded all required KEGG data,
            then this function will not fail. This is to enable the --only-processing option.
            Note that in this case we require all KEGG data to be pre-downloaded to avoid mixing
            older and newer KEGG data - so if this data is only partially downloaded, the function
            will raise an error even if this parameter is False.
        """

        if anvio.DEBUG:
            file_str = ", ".join(files_to_check)
            self.run.warning(f"We are looking for the following files to see if the KEGG data already "
                             f"exists on you computer: {file_str}")

        files_that_exist = []
        for f in files_to_check:
            if os.path.exists(f):
                if fail_if_exists:
                    raise ConfigError(f"It seems you already have data at {f}, please use the `--reset` flag "
                                      "or delete the KEGG data directory manually if you want to re-download KEGG data. "
                                      "See also the --only-processing option, which you can use if you already "
                                      "have all required KEGG data in that folder. (API users: skip this sanity "
                                      "check by initializing this class with `skip_init=True`)")
                else:
                    files_that_exist.append(f)

        if files_that_exist:
            exist_str = "\n".join(files_that_exist)
            # we require all data to be present. Otherwise we might produce chimeric KEGG data.
            if files_that_exist != files_to_check:
                raise ConfigError(f"We found some, but not all, required KEGG data on your computer in the KEGG "
                                  f"data directory. Since you don't have everything you need, we need you to re-download "
                                  f"everything from scratch. Please re-run this program using the --reset flag, and if "
                                  f"you were using the --only-processing option, remove that flag. :) HOWEVER, if you notice that "
                                  "KEGG BRITE data does not appear to be in the upcoming list, but you don't actually want "
                                  "to download BRITE data, then you can just add the --skip-brite-hierarchies to your previous "
                                  f"command and be on your way (ie, no --reset needed). Here is the KEGG data we found:\n{exist_str}")

            self.run.warning(f"We found already-downloaded KEGG data on your computer. Setup will continue using "
                             f"this data. However, if you think everything should be re-downloaded from scratch, please kill this program "
                             f"and restart it using the `--reset` flag. Here is the data you already have, in case you "
                             f"need to check it to make sure we are not using something that is too old:\n"
                             f"{exist_str}")

        if self.only_processing and not files_that_exist:
            raise ConfigError(f"We noticed that there is no KEGG data on your computer at {self.kegg_data_dir} even "
                              f"though you used the --only-processing option. If you don't actually have KEGG data already "
                              f"downloaded, you should get rid of the --only-processing flag and re-run this program. If you "
                              f"know that you DO have KEGG data, perhaps you gave us the wrong data directory?")


    def setup_from_archive(self):
        """This function sets up the KEGG data directory from an archive of a previously-setup KEGG data directory.

        To do so, it unpacks the archive and checks its structure and that all required components are there.
        """

        self.run.info("KEGG archive", self.kegg_archive_path)
        self.progress.new('Unzipping KEGG archive file...')
        if not self.kegg_archive_path.endswith("tar.gz"):
            self.progress.reset()
            raise ConfigError("The provided archive file %s does not appear to be an archive at all. Perhaps you passed "
                              "the wrong file to anvi'o?" % (self.kegg_archive_path))
        unpacked_archive_name = "KEGG_archive_unpacked"
        utils.tar_extract_file(self.kegg_archive_path, output_file_path=unpacked_archive_name, keep_original=True)

        self.progress.update('Checking KEGG archive structure and contents...')
        archive_is_ok = self.kegg_archive_is_ok(unpacked_archive_name)
        archive_contains_brite = self.check_archive_for_brite(unpacked_archive_name)
        archive_contains_binary_relations = self.check_archive_for_binary_relations(unpacked_archive_name)
        archive_contains_map_images = self.check_archive_for_map_images(unpacked_archive_name)
        self.progress.end()
        if archive_is_ok:
            if os.path.exists(self.kegg_data_dir):
                shutil.rmtree(self.kegg_data_dir)
            path_to_kegg_in_archive = os.path.join(unpacked_archive_name, "KEGG")
            shutil.move(path_to_kegg_in_archive, self.kegg_data_dir)
            shutil.rmtree(unpacked_archive_name)

            if not archive_contains_brite and not self.skip_brite_hierarchies:
                self.run.warning("The KEGG data archive does not contain the necessary files to set up BRITE hierarchy classification. "
                                 "This is not a problem, and KEGG set up proceeded without it. BRITE is guaranteed to be set up when "
                                 "downloading the latest version of KEGG with `anvi-setup-kegg-data`.")

            if not archive_contains_binary_relations and not self.skip_binary_relations:
                self.run.warning(
                    "The KEGG data archive does not contain the binary relation files needed for "
                    "`anvi-reaction-network`. This is not a problem, and KEGG setup proceeded "
                    "without it. Binary relation files are guaranteed to be set up when "
                    "downloading the latest version of KEGG with `anvi-setup-kegg-data`."
                )

            if not archive_contains_map_images and not self.skip_map_images:
                self.run.warning(
                    "The KEGG data archive does not contain the pathway map image files used for "
                    "pathway visualization. This is not a problem, and KEGG setup proceeded "
                    "without it. Map image files are guaranteed to be set up when downloading the "
                    "latest version of KEGG with `anvi-setup-kegg-data`."
                )

            # if necessary, warn user about migrating the modules db
            self.check_modules_db_version()

        else:
            debug_output = f"We kept the unpacked archive for you to take a look at it. It is at " \
                           f"{os.path.abspath(unpacked_archive_name)} and you may want " \
                           f"to delete it after you are done checking its contents."
            if not anvio.DEBUG:
                shutil.rmtree(unpacked_archive_name)
                debug_output = "The unpacked archive has been deleted, but you can re-run the script with the --debug " \
                               "flag to keep it if you want to see its contents."
            else:
                self.run.warning(f"The unpacked archive file {os.path.abspath(unpacked_archive_name)} was kept for "
                                 f"debugging purposes. You may want to clean it up after you are done looking through it.")

            raise ConfigError(f"SETUP FAILED. The provided archive file is missing some critical files, "
                              f"so anvi'o is unable to use it. {debug_output}")


    def check_modules_db_version(self):
        """This function checks if the MODULES.db is out of date and if so warns the user to migrate it"""

        # get current version of db
        db_conn = db.DB(self.kegg_modules_db_path, None, ignore_version=True)
        current_db_version = int(db_conn.get_meta_value('version'))
        db_conn.disconnect()

        # if modules.db is out of date, give warning
        target_version = int(anvio.tables.versions_for_db_types['modules'])
        if current_db_version != target_version:
            self.run.warning(f"Just so you know, the KEGG archive that was just set up contains an outdated MODULES.db (version: "
                             f"{current_db_version}). You may want to run `anvi-migrate` on this database before you do anything else. "
                             f"Here is the path to the database: {self.kegg_modules_db_path}")


    def check_archive_for_brite(self, unpacked_archive_path):
        """Check the archive for the BRITE directory and 'hierarchy of hierarchies' json file.

        It is ok for archives not to have these present, but let the user know.
        """

        is_brite_included = True

        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        brite_directories_and_files = [self.brite_data_dir,
                                       self.kegg_brite_hierarchies_file]
        for f in brite_directories_and_files:
            path_to_f_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(f))
            if not os.path.exists(path_to_f_in_archive) and not self.skip_brite_hierarchies:
                is_brite_included = False
                if anvio.DEBUG:
                    self.run.warning(f"The KEGG archive does not contain the following optional BRITE file or directory: {path_to_f_in_archive}")

        return is_brite_included


    def check_archive_for_binary_relations(self, unpacked_archive_path):
        """
        Check the archive for the binary relations directory and files.

        It is ok for archives not to have these present, but let the user know.
        """
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        binary_relation_data_dir = os.path.join(
            path_to_kegg_in_archive, os.path.basename(self.binary_relation_data_dir)
        )
        if os.path.isdir(binary_relation_data_dir):
            is_binary_relation_dir_included = True
        else:
            is_binary_relation_dir_included = False
            if anvio.DEBUG and not self.skip_binary_relations:
                self.run.warning(
                    "The KEGG archive does not contain the following optional binary relations "
                    f"directory needed for `anvi-reaction-network`: {binary_relation_data_dir}"
                )

        if is_binary_relation_dir_included:
            missing_files = []
            for file in self.kegg_binary_relation_files.values():
                path = os.path.join(binary_relation_data_dir, file)
                if not os.path.isfile(path):
                    missing_files.append(file)
            if anvio.DEBUG and missing_files:
                self.run.warning(
                    "The following binary relation files expected in an up-to-date anvi'o KEGG "
                    f"installation are missing from the directory, '{binary_relation_data_dir}', "
                    f"in the archive: {', '.join(missing_files)}"
                )

        return is_binary_relation_dir_included


    def check_archive_for_map_images(self, unpacked_archive_path):
        """
        Check the archive for the pathway map directory and image files.

        It is ok for archives not to have these present, but let the user know.
        """
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        map_image_data_dir = os.path.join(
            path_to_kegg_in_archive, os.path.basename(self.map_image_data_dir)
        )
        if os.path.isdir(map_image_data_dir):
            is_map_image_dir_included = True
        else:
            is_map_image_dir_included = False
            if anvio.DEBUG and not self.skip_map_images:
                self.run.warning(
                    f"The KEGG archive does not contain the following optional pathway map images "
                    f"directory, which is used in pathway visualization."
                )

        return is_map_image_dir_included


    def setup_kegg_snapshot(self):
        """This is the default setup strategy in which we unpack a specific KEGG archive.

        We do this so that everyone who uses the same release of anvi'o will also have the same default KEGG
        data, which facilitates sharing and also means they do not have to continuously re-annotate their datasets
        when KEGG is updated.

        It is essentially a special case of setting up from an archive.
        """

        if anvio.DEBUG:
            self.run.info("Downloading from: ", self.default_kegg_data_url)
            self.run.info("Downloading to: ", self.default_kegg_archive_file)
        utils.download_file(self.default_kegg_data_url, self.default_kegg_archive_file, progress=self.progress, run=self.run)

        # a hack so we can use the archive setup function
        self.kegg_archive_path = self.default_kegg_archive_file
        self.setup_from_archive()

        # if all went well, let's get rid of the archive we used and the log file
        if not anvio.DEBUG:
            os.remove(self.default_kegg_archive_file)
        else:
            self.run.warning(f"Because you used the --debug flag, the KEGG archive file at {self.default_kegg_archive_file} "
                             "has been kept. You may want to remove it later.")


    def kegg_archive_is_ok(self, unpacked_archive_path):
        """This function checks the structure and contents of an unpacked KEGG archive and returns True if it is as expected.

        Please note that we check for existence of the files that are necessary to run KEGG scripts, but we don't check the file
        formats. This means that people could technically trick this function into returning True by putting a bunch of crappy files
        with the right names/paths into the archive file. But what would be the point of that?

        We also don't care about the contents of certain folders (ie modules) because they are not being directly used
        when running KEGG scripts. In the case of modules, all the information should already be in the MODULES.db so we don't
        waste our time checking that all the module files are there. We only check that the directory is there. If later changes
        to the implementation require the direct use of the files in these folders, then this function should be updated
        to check for those.

        Parameters
        ==========
        unpacked_archive_path : str
            Path to the unpacked archive directory
        """

        is_ok = True

        # check top-level files and folders
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        expected_directories_and_files = [self.orphan_data_dir,
                                          self.kegg_module_data_dir,
                                          self.kegg_hmm_data_dir,
                                          self.ko_list_file_path,
                                          self.kegg_module_file,
                                          self.kegg_modules_db_path]
        for f in expected_directories_and_files:
            path_to_f_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(f))
            if not os.path.exists(path_to_f_in_archive):
                is_ok = False
                if anvio.DEBUG:
                    self.run.warning(f"The KEGG archive does not contain the following expected file or directory: "
                                     f"{path_to_f_in_archive}")

        # check hmm files
        path_to_hmms_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(self.kegg_hmm_data_dir))
        kofam_hmm_basename = os.path.basename(self.kofam_hmm_file_path)
        expected_hmm_files = [kofam_hmm_basename]
        for h in expected_hmm_files:
            path_to_h_in_archive = os.path.join(path_to_hmms_in_archive, h)
            if not os.path.exists(path_to_h_in_archive):
                is_ok = False
                if anvio.DEBUG:
                    self.run.warning(f"The KEGG archive does not contain the following expected hmm file: "
                                     f"{path_to_h_in_archive}")
            expected_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
            for ext in expected_extensions:
                path_to_expected_hmmpress_file = path_to_h_in_archive + ext
                if not os.path.exists(path_to_expected_hmmpress_file):
                    is_ok = False
                    if anvio.DEBUG:
                        self.run.warning(f"The KEGG archive does not contain the following expected `hmmpress` output: "
                                         f"{path_to_expected_hmmpress_file}")

        return is_ok


    def setup_all_data_from_archive_or_snapshot(self):
        """This driver function controls whether we download one of our KEGG snapshots and set that up, or
        set up directly from an archive file already on the user's computer.
        """

        if os.path.exists(self.kegg_data_dir) and not self.args.reset:
            raise ConfigError(f"The directory {self.kegg_data_dir} already exists. Are you sure you want to "
                              f"overwrite it? If yes, feel free to restart this program with the --reset flag.")

        if self.kegg_archive_path:
            self.setup_from_archive()
        else:
            self.setup_kegg_snapshot()


    def check_user_input_dir_format(self):
        """This function checks whether the user input directory exists and contains the required subfolders

        The required subfolders are:
            modules : directory containing the user's metabolic pathway definitions (as text files)
        """

        for path in [self.user_input_dir, self.user_module_data_dir]:
            if not os.path.exists(path):
                raise ConfigError(f"There is a problem with the input directory you provided. The following path does not "
                                  f"exist: '{path}'. Please make sure that your input folder exists and that it follows the "
                                  f"formatting requirements. We're sorry for asking this of you, but it really helps us make "
                                  f"sure everything will go smoothly.")

            file_list = [f for f in glob.glob(os.path.join(path, '*'))]
            if not file_list:
                raise ConfigError(f"The folder '{path}' appears to be empty, so we have no data to work with. Please make "
                                  f"sure that you have provided the correct input directory and formatted it correctly so "
                                  f"that anvi'o can find your data.")


    def is_user_database_exists(self):
        """This function checks whether user data has already been set up in the provided input directory."""

        if os.path.exists(self.user_modules_db_path):
            raise ConfigError(f"It seems you already have a user modules database installed at '{self.user_modules_db_path}', "
                              f"please use the --reset flag or delete this file manually if you want to re-generate it.")


    def process_pathway_file(self):
        """This function reads the kegg pathway map file into a dictionary. It should be called during setup to get the KEGG pathway ids so the pathways can be downloaded.

        The structure of this file is like this:

        +C	Map number
        #<h2><a href="/kegg/kegg2.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a>&nbsp; KEGG Pathway Maps</h2>
        !
        A<b>Metabolism</b>
        B  Global and overview maps
        C    01100  Metabolic pathways
        C    01110  Biosynthesis of secondary metabolites
        C    01120  Microbial metabolism in diverse environments
        C    01200  Carbon metabolism
        C    01210  2-Oxocarboxylic acid metabolism

        Initial lines can be ignored and thereafter the line's information can be determined by the one-letter code at the start.
        A = Category of Pathway Map
        B = Sub-category of Pathway Map
        C = Pathway Map identifier number and name

        Note that not all Pathway Maps that we download will have ORTHOLOGY fields. We don't exclude these here, but processing later
        will have to be aware of the fact that not all pathways will have associated KOs.

        We do, however, exclude Pathway Maps that don't have existing `koXXXXX` identifiers (these yield 404 errors when attempting to
        download them). For instance, we exclude those that start with the code 010 (chemical structure maps) or with 07 (drug structure maps).
        """

        self.pathway_dict = {}

        filesnpaths.is_file_exists(self.kegg_pathway_file)
        filesnpaths.is_file_plain_text(self.kegg_pathway_file)

        f = open(self.kegg_pathway_file, 'r')
        self.progress.new("Parsing KEGG Pathway file")

        current_category = None
        current_subcategory = None


        for line in f.readlines():
            line = line.strip('\n')
            first_char = line[0]

            # garbage lines
            if first_char in ["+", "#", "!"]:
                continue
            else:
                # Category
                if first_char == "A":
                    fields = re.split('<[^>]*>', line) # we split by the html tag here
                    current_category = fields[1]
                # Sub-category
                elif first_char == "B":
                    fields = re.split('\s{2,}', line) # don't want to split the subcategory name, so we have to split at least 2 spaces
                    current_subcategory = fields[1]
                elif first_char == "C":
                    fields = re.split('\s{2,}', line)
                    konum = "ko" + fields[1]
                    if konum[:5] != "ko010" and konum[:4] != "ko07":
                        self.pathway_dict[konum] = {"name" : fields[2], "category" : current_category, "subcategory" : current_subcategory}
                # unknown code
                else:
                    raise ConfigError("While parsing the KEGG file %s, we found an unknown line code %s. This has "
                                      "made the file unparseable. It is likely that an update to KEGG has broken "
                                      "things such that anvi'o doesn't know what is going on anymore. Sad, we know. :( "
                                      "Please contact the developers to see if this is a fixable issue, and in the "
                                      "meantime use an older version of the KEGG data directory (if you have one). "
                                      "If we cannot fix it, we may be able to provide you with a legacy KEGG "
                                      "data archive that you can use to setup KEGG with the --kegg-archive flag." % (self.kegg_pathway_file, first_char))
        self.progress.end()


    def get_accessions_from_htext_file(self, htext_file):
        """This function can read generic KEGG htext files to get a list of accessions.

        Here is one example of the file structure, taken from a COMPOUND htext file:
        +D	Biochemical compound
        #<h2><a href="/kegg/brite.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a>&nbsp; Compounds with Biological Roles</h2>
        !
        A<b>Organic acids</b>
        B  Carboxylic acids [Fig]
        C    Monocarboxylic acids
        D      C00058  Formate; Methanoate
        D      C00033  Acetate; Ethanoate
        D      C00163  Propionate; Propanoate
        D      C00246  Butyrate; Butanoate

        The +(letter) at the start of the first line indicates how many levels the hierarchy contains (often C or D). For the purpose of
        downloading other files, we only need the accessions from the lowest level. All other information is skipped when reading the file.

        PARAMETERS
        ==========
        htext_file : str
            The filename of the hierachical text file downloaded from KEGG

        RETURNS
        =======
        accession_list : list of str
            Contains the KEGG identifiers contained in the lowest hierarchy of this file.
        """

        accession_list = []

        filesnpaths.is_file_exists(htext_file)
        filesnpaths.is_file_plain_text(htext_file)

        f = open(htext_file, 'r')
        self.progress.new(f"Parsing KEGG htext file: {htext_file}")

        target_level = None
        for line in f.readlines():
            line = line.strip('\n')
            first_char = line[0]

            if first_char == '+':  # first line of the file; second character gives us target level
                target_level = line[1]
            elif first_char == target_level: # need to extract the accession, which is second field (split on spaces)
                fields = re.split('\s+', line)
                accession_list.append(fields[1])
            else: # skip everything else
                continue

        self.progress.end()

        num_acc = len(accession_list)
        self.run.info("Number of accessions found in htext file", num_acc)
        return accession_list


    def download_generic_htext(self, h_accession, download_dir="./"):
        """Downloads the KEGG htext file for the provided accession.

        PARAMETERS
        ==========
        h_accession : str
            The accession for a KEGG hierarchy file
        download_dir : str
            Path to directory where file will be downloaded. Current working directory by default.
        """

        htext_url_prefix = "https://www.genome.jp/kegg-bin/download_htext?htext="
        htext_url_suffix = ".keg&format=htext&filedir="
        htext_url = htext_url_prefix+h_accession+htext_url_suffix

        htext_file = h_accession + ".keg"
        path_to_download_to = os.path.join(download_dir,htext_file)

        if filesnpaths.is_file_exists(path_to_download_to, dont_raise=True):
            if not self.args.reset:
                raise ConfigError(f"The file at {path_to_download_to} already exists. If you are "
                                  f"sure that you want to download it, you can avoid this error message "
                                  f"by using the 'reset' parameter. Make sure that won't erase your other "
                                  f"KEGG data, though.")

        try:
            utils.download_file(htext_url, path_to_download_to, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError(f"Anvi'o failed to download the KEGG htext file for {h_accession} from {htext_url}.")

        return path_to_download_to


    def download_generic_flat_file(self, accession, download_dir="./"):
        """Downloads the flat file for the given accession from the KEGG API.

        PARAMETERS
        ==========
        accession : str
            A KEGG identifier
        download_dir : str
            Path to the directory in which to download the file. Current working directory by default.
        """

        file_path = os.path.join(download_dir, accession)
        if filesnpaths.is_file_exists(file_path, dont_raise=True):
            if not self.args.reset:
                raise ConfigError(f"The file at {file_path} already exists. If you are "
                                  f"sure that you want to download it, you can avoid this error message "
                                  f"by using the 'reset' parameter. Make sure that won't erase your other "
                                  f"KEGG data, though.")

        utils.download_file(self.kegg_rest_api_get + '/' + accession,
            file_path, progress=self.progress, run=self.run)
        # verify entire file has been downloaded
        f = open(file_path, 'r')
        f.seek(0, os.SEEK_END)
        f.seek(f.tell() - 4, os.SEEK_SET)
        last_line = f.readline().strip('\n')
        if not last_line == '///':
            raise ConfigError(f"The KEGG flat file {file_path} was not downloaded properly. We were expecting the last line in the file "
                              f"to be '///', but instead it was {last_line}. Formatting of these files may have changed on the KEGG website. "
                              f"Please contact the developers to see if this is a fixable issue.")


    def download_kegg_files_from_hierarchy(self, h_accession, download_dir="./"):
        """Given the accession of a KEGG hierarchy, this function downloads all of its flat files.

        PARAMETERS
        ==========
        h_accession : str
            The accession for a KEGG hierarchy file
        download_dir : str
            Path to the directory in which to download the files. Current working directory by default.
            (a folder to store the hierarchy's flat files will be generated in this folder)
        """

        filesnpaths.is_output_dir_writable(download_dir)

        htext_filename = self.download_generic_htext(h_accession, download_dir)
        acc_list = self.get_accessions_from_htext_file(htext_filename)

        download_dir_name = os.path.join(download_dir,h_accession)
        filesnpaths.gen_output_directory(download_dir_name, delete_if_exists=self.args.reset)

        self.run.info("KEGG Module Database URL", self.kegg_rest_api_get)
        self.run.info("Number of KEGG files to download", len(acc_list))
        self.run.info("Directory to store files", download_dir_name)

        # download all modules
        for acc in acc_list:
            self.download_generic_flat_file(acc, download_dir_name)


    def download_pathways(self):
        """This function downloads the KEGG Pathways.

        To do so, it first processes a KEGG file containing pathway and map identifiers into a dictionary via the process_pathway_file()
        function. To verify that each file has been downloaded properly, we check that the last line is '///'.
        """

        # note that this is the same as the REST API for modules - perhaps at some point this should be printed elsewhere so we don't repeat ourselves.
        self.run.info("KEGG Pathway Database URL", self.kegg_rest_api_get)

        # download the kegg pathway file, which lists all modules
        try:
            utils.download_file(self.kegg_pathway_download_path, self.kegg_pathway_file, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download the KEGG Pathway htext file from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")

        # get pathway dict
        self.process_pathway_file()
        self.run.info("Number of KEGG Pathways", len(self.pathway_dict.keys()))

        # download all pathways
        for konum in self.pathway_dict.keys():
            file_path = os.path.join(self.pathway_data_dir, konum)
            utils.download_file(self.kegg_rest_api_get + '/' + konum,
                file_path, progress=self.progress, run=self.run)
            # verify entire file has been downloaded
            f = open(file_path, 'r')
            f.seek(0, os.SEEK_END)
            f.seek(f.tell() - 4, os.SEEK_SET)
            last_line = f.readline().strip('\n')
            if not last_line == '///':
                raise ConfigError("The KEGG pathway file %s was not downloaded properly. We were expecting the last line in the file "
                                  "to be '///', but instead it was %s. Formatting of these files may have changed on the KEGG website. "
                                  "Please contact the developers to see if this is a fixable issue. If it isn't, we may be able to "
                                  "provide you with a legacy KEGG data archive that you can use to setup KEGG with the --kegg-archive flag."
                                  % (file_path, last_line))


    def extract_data_field_from_kegg_file(self, file_path, target_field):
        """This function parses a KEGG file and returns the data value associated with the given target field.

        It can work on flat-text files obtained via the REST API (ie, self.kegg_rest_api_get).
        """

        data_to_return = []

        f = open(file_path, 'r')
        current_data_name = None

        for line in f.readlines():
            line = line.strip('\n')

            fields = re.split('\s{2,}', line)
            data_vals = None
            data_def = None
            line_entries = []

            # when data name unknown, parse from first field
            if line[0] != ' ':
                current_data_name = fields[0]
            if line[0] == ' ' and not current_data_name:
                raise ConfigError(f"Uh oh. While trying to parse the KEGG file at {file_path}, we couldn't "
                "find the data field associated with the line '{line}'.")

            # note that if data name is known, first field still exists but is actually the empty string ''
            if len(fields) > 1:
                data_vals = fields[1]

            if (current_data_name == target_field) and (data_vals is not None):
                data_to_return.append(data_vals)

        f.close()

        return data_to_return


    def create_user_modules_dict(self):
        """This function establishes the self.module_dict parameter for user modules.

        It is essentially a replacement for the process_module_file() function.
        Since users will not have a modules file to process, we simply create the dictionary from the
        file names they provide for their module definitions. We don't add any dictionary values,
        but we won't need them (we hope).
        """

        user_module_list = [os.path.basename(k) for k in glob.glob(os.path.join(self.user_module_data_dir, '*'))]
        self.module_dict = {key: {} for key in user_module_list}

        # sanity check that they also have KEGG data since we need to compare module names
        if not os.path.exists(self.kegg_modules_db_path):
            raise ConfigError(f"Wait a second. We understand that you are setting up user-defined metabolism data, but "
                              f"unfortunately you need to FIRST have KEGG data set up on your computer. Why, you ask? "
                              f"Well, we need to make sure none of your module names overlap with those "
                              f"in the KEGG MODULES database. Long story short, we looked for KEGG data at "
                              f"{self.kegg_modules_db_path} but we couldn't find it. If this is the wrong place for us to be "
                              f"looking, please run this program again and use the --kegg-data-dir parameter to tell us where "
                              f"to find it.")

        # sanity check that user module names are distinct
        kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, args=self.args, quiet=True)
        kegg_mods = set(kegg_modules_db.get_all_modules_as_list())
        user_mods = set(user_module_list)
        bad_user_mods = kegg_mods.intersection(user_mods)
        if bad_user_mods:
            bad_mods_str = ", ".join(bad_user_mods)
            n = len(bad_user_mods)
            raise ConfigError(f"Hol'up a minute. You see, there {P('is a module', n, alt='are some modules')} "
                              f"in your user-defined modules data (at {self.user_module_data_dir}) which {P('has', n, alt='have')} "
                              f"the same name as an existing KEGG module. This is not allowed, for reasons. Please name {P('that module', n, alt='those modules')} "
                              f"differently. Append an underscore and your best friend's name to {P('it', n, alt='them')} or something. Just make sure it's "
                              f"unique. OK? ok. Here is the list of module names you should change: {bad_mods_str}")


    def setup_modules_db(self, db_path, module_data_directory, brite_data_directory=None, source='KEGG', skip_brite_hierarchies=False):
        """This function creates a Modules DB at the specified path."""

        if filesnpaths.is_file_exists(db_path, dont_raise=True):
            if self.overwrite_modules_db:
                os.remove(db_path)
            else:
                raise ConfigError(f"Woah there. There is already a modules database at {db_path}. If you really want to make a new modules database "
                                  f"in this folder, you should either delete the existing database yourself, or re-run this program with the "
                                  f"--overwrite-output-destinations flag. But the old database will go away forever in that case. Just making "
                                  f"sure you are aware of that, so that you have no regrets.")
        try:
            mod_db = ModulesDatabase(db_path, module_data_directory=module_data_directory, brite_data_directory=brite_data_directory, data_source=source, args=self.args, module_dictionary=self.module_dict, pathway_dictionary=self.pathway_dict, brite_dictionary=self.brite_dict, skip_brite_hierarchies=skip_brite_hierarchies, run=run, progress=progress)
            mod_db.create()
        except Exception as e:
            print(e)
            raise ConfigError("While attempting to build the MODULES.db, anvi'o encountered an error, which should be printed above. "
                              "If you look at that error and it seems like something you cannot handle, please contact the developers "
                              "for assistance. :) ")


    def setup_user_data(self):
        """This function sets up user metabolism data from the provided input directory.

        It processes the user's module files into the USER_MODULES.db.
        """

        self.create_user_modules_dict()
        self.setup_modules_db(db_path=self.user_modules_db_path, module_data_directory=self.user_module_data_dir, source='USER', skip_brite_hierarchies=True)


class KOfamDownload(KeggSetup):
    """Class for setting up KOfam HMM profiles.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to command-line programs relying on this
        class, such as `anvi-setup-kegg-data`. If using this class through the API, please
        provide a Namespace object with the Boolean 'reset' parameter.
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories
        when testing this class.
    """

    def __init__(self, args, run=run, progress=progress, skip_init=False):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress
        self.skip_init = skip_init
        self.include_stray_kos = True if A('include_stray_KOs') else False

        self.run.info_single("Info from KOfam Download")
        self.run.info("Stray KOs will be processed (`--include-stray-KOs` flag)", self.include_stray_kos)

        KeggSetup.__init__(self, self.args, skip_init=self.skip_init)

        filesnpaths.is_program_exists('hmmpress')
        if self.include_stray_kos:
            filesnpaths.is_program_exists('hmmbuild')
            filesnpaths.is_program_exists('muscle')

        # ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        # dictionary mapping downloaded file name to final decompressed file name or folder location
        self.kofam_files = {'ko_list.gz': self.ko_list_file_path, 'profiles.tar.gz': self.kegg_data_dir}

        expected_files_for_kofams = [self.ko_list_file_path]
        if self.only_processing:
            expected_files_for_kofams.append(os.path.join(self.kegg_data_dir, 'profiles.tar.gz'))
        else:
            expected_files_for_kofams.append(self.kofam_hmm_file_path)

        if not args.reset and not anvio.DEBUG and not self.skip_init:
            self.is_database_exists(expected_files_for_kofams, fail_if_exists=(not self.only_processing))

        if self.download_from_kegg and not self.only_processing and not self.kegg_archive_path and not self.skip_init:
            filesnpaths.gen_output_directory(self.kegg_hmm_data_dir, delete_if_exists=args.reset)
            filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)


    def download_profiles(self):
        """This function downloads the Kofam profiles."""

        self.run.info("Kofam Profile Database URL", self.database_url)

        try:
            for file_name in self.kofam_files.keys():
                utils.download_file(self.database_url + '/' + file_name,
                    os.path.join(self.kegg_data_dir, file_name), progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download KEGG KOfam profiles from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")


    def decompress_profiles(self):
        """This function decompresses the Kofam profiles."""

        self.progress.new('Decompressing files')
        for file_name in self.kofam_files.keys():
            self.progress.update('Decompressing file %s' % file_name)
            full_path = os.path.join(self.kegg_data_dir, file_name)

            if full_path.endswith("tar.gz"):
                utils.tar_extract_file(full_path, output_file_path=self.kofam_files[file_name], keep_original=False)
            else:
                utils.gzip_decompress_file(full_path, output_file_path=self.kofam_files[file_name], keep_original=False)

            self.progress.update("File decompressed. Yay.")
        self.progress.end()


    def confirm_downloaded_profiles(self):
        """This function verifies that all Kofam profiles have been properly downloaded.

        It is intended to be run after the files have been decompressed. The profiles directory should contain hmm files (ie, K00001.hmm);
        all KO numbers from the ko_list file (except those in ko_skip_list) should be included.

        This function must be called after setup_ko_dict() so that the self.ko_dict attribute is established.
        """

        ko_nums = self.ko_dict.keys()
        for k in ko_nums:
            if k not in self.ko_skip_list:
                hmm_path = os.path.join(self.kegg_data_dir, f"profiles/{k}.hmm")
                if not os.path.exists(hmm_path):
                    raise ConfigError(f"The KOfam HMM profile at {hmm_path} does not exist. This probably means that something went wrong "
                                      f"while downloading the KOfam database. Please run `anvi-setup-kegg-data` with the --reset "
                                      f"flag. If that still doesn't work, please contact the developers to see if the issue is fixable. "
                                      f"If it isn't, we may be able to provide you with a legacy KEGG data archive that you can use to "
                                      f"setup KEGG with the --kegg-archive flag.")


    def move_orphan_files(self):
        """This function moves the following to the orphan files directory:

            - profiles that do not have ko_list entries
            - profiles whose ko_list entries have no scoring threshold (in ko_no_threshold_list)

        And, the following profiles should not have been downloaded, but if they were then we move them, too:
            - profiles whose ko_list entries have no data at all (in ko_skip_list)
        """

        if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
            raise ConfigError(f"Hmm. Something is out of order. The orphan data directory {self.orphan_data_dir} does not exist "
                              "yet, but it needs to in order for the move_orphan_files() function to work.")

        no_kofam_path = os.path.join(self.orphan_data_dir, "hmm_profiles_with_no_kofams.hmm")
        no_kofam_file_list = []
        no_threshold_file_list = []
        no_data_path = os.path.join(self.orphan_data_dir, "hmm_profiles_with_kofams_with_no_data.hmm")
        no_data_file_list = []

        hmm_list = [k for k in glob.glob(os.path.join(self.kegg_data_dir, 'profiles/*.hmm'))]
        for hmm_file in hmm_list:
            ko = re.search('profiles/(K\d{5})\.hmm', hmm_file).group(1)
            if ko not in self.ko_dict.keys():
                if ko in self.ko_no_threshold_list:
                    no_threshold_file_list.append(hmm_file)
                elif ko in self.ko_skip_list: # these should not have been downloaded, but if they were we will move them
                    no_data_file_list.append(hmm_file)
                else:
                    no_kofam_file_list.append(hmm_file)

        # now we concatenate the orphan and stray KO hmms into the orphan data directory
        if no_kofam_file_list:
            utils.concatenate_files(no_kofam_path, no_kofam_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning(f"Please note that while anvi'o was building your databases, she found {len(no_kofam_file_list)} "
                             f"HMM profiles that did not have any matching KOfam entries. We have removed those HMM "
                             f"profiles from the final database. You can find them under the directory '{self.orphan_data_dir}'.")

        if no_threshold_file_list:
            utils.concatenate_files(self.stray_ko_hmms_from_kegg, no_threshold_file_list, remove_concatenated_files=False)
            filesnpaths.gen_output_directory(os.path.join(self.orphan_data_dir, "profiles"), delete_if_exists=True)
            for k_path in no_threshold_file_list:
                k = os.path.basename(k_path)
                # move individual profiles temporarily to the orphan data dir, so they don't get combined with the regular KOs
                # but we can still use them later if necessary for --include-stray-KOs
                os.rename(k_path, os.path.join(self.orphan_data_dir, f"profiles/{k}"))
            self.progress.reset()
            self.run.warning(f"Please note that while anvi'o was building your databases, she found {len(no_threshold_file_list)} "
                             f"KOfam entries that did not have any threshold to remove weak hits. We have removed those HMM "
                             f"profiles from the final database. You can find them under the directory '{self.orphan_data_dir}'. "
                             f"If you used the flag --include-stray-KOs, we will estimate their bit score thresholds using KEGG GENES "
                             f"data so that you can annotate these KOs downstream if you wish.")

        if no_data_file_list:
            utils.concatenate_files(no_data_path, no_data_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning(f"Please note that while anvi'o was building your databases, she found {len(no_data_file_list)} "
                             f"HMM profiles that did not have any associated data (besides an annotation) in their KOfam entries. "
                             f"We have removed those HMM profiles from the final database. You can find them under the directory "
                             f"'{self.orphan_data_dir}'.")


    def exec_hmmpress_command_on_ko_file(self, hmm_file_path, log_file_path):
        """Given a path to a set of KO HMMs and a log file path, this function executes the appropriate
        `hmmpress` command and deletes the log file afterwards if it was successful.
        """

        cmd_line = ['hmmpress', hmm_file_path]
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Hmm. There was an error while running `hmmpress` on the Kofam HMM profiles. "
                              "Check out the log file ('%s') to see what went wrong." % (log_file_path))
        else:
            # getting rid of the log file because hmmpress was successful
            os.remove(log_file_path)


    def run_hmmpress(self):
        """This function concatenates the Kofam profiles and runs hmmpress on them."""

        self.progress.new('Preparing Kofam HMM Profiles')

        self.progress.update('Verifying the Kofam directory %s contains all HMM profiles' % self.kegg_data_dir)
        self.confirm_downloaded_profiles()

        self.progress.update('Handling orphan files')
        self.move_orphan_files()

        self.progress.update('Concatenating HMM profiles into one file...')
        hmm_list = [k for k in glob.glob(os.path.join(self.kegg_data_dir, 'profiles/*.hmm'))]
        utils.concatenate_files(self.kofam_hmm_file_path, hmm_list, remove_concatenated_files=False)

        self.progress.update('Running hmmpress on KOs...')
        self.exec_hmmpress_command_on_ko_file(self.kofam_hmm_file_path, os.path.join(self.kegg_hmm_data_dir, '00_hmmpress_log.txt'))

        self.progress.end()


    def download_ko_files(self, kos_to_download, destination_dir, dont_raise=True):
        """Multi-threaded download of KEGG Orthology files.

        Parameters
        ==========
        kos_to_download: list of str
            List of KOs to download Orthology files for
        destination_dir: file path
            Where to download the files to
        dont_raise : Boolean
            If True (default), this function won't raise an error if some files failed to download.
        Returns
        =======
        undownloaded : list of str
            List of KOs that failed to download (will be empty if all were successful).
        """

        num_kos = len(kos_to_download)
        self.progress.new('Downloading KO files', progress_total_items=num_kos)

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for ko in kos_to_download:
            ko_file_path = os.path.join(destination_dir, ko)
            url = self.kegg_rest_api_get + '/' + ko
            input_queue.put((url, ko_file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < num_kos:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {num_kos} KO files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        self.progress.end()
        for worker in workers:
            worker.terminate()
        if undownloaded:
            if dont_raise:
                self.run.warning(f"Files for the following KOs failed to download despite multiple attempts: "
                                f"{', '.join(undownloaded)}. If this is unacceptable to you, you can try to "
                                f"re-run this program to see if things will work on the next try.")
            else:
                raise ConfigError(f"Files for the following KOs failed to download despite multiple attempts: "
                                  f"{', '.join(undownloaded)}. Since the function responsible for handling this was "
                                  f"told to quit should this happen, well, here we are. If skipping these failed KOs "
                                  f"is okay, you could always run this function with `dont_raise=True`.")

        return undownloaded


    def download_kegg_genes_files(self, genes_to_download, destination_dir, dont_raise=True):
        """Multi-threaded download of KEGG GENES files.

        Parameters
        ==========
        genes_to_download: list of str
            List of KEGG GENES accessions to download. Example format: "ctc:CTC_p60" (lowercase organism code, colon, gene accession)
        destination_dir: file path
            Where to download the files to
        dont_raise : Boolean
            If True (default), this function won't raise an error if some files failed to download.
        Returns
        =======
        undownloaded : List
            List of files that failed to download (will be empty if all were successful).
        """

        num_genes = len(genes_to_download)
        self.progress.new('Downloading KEGG GENES files', progress_total_items=num_genes)

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for g in genes_to_download:
            genes_file_path = os.path.join(destination_dir, g)
            url = self.kegg_rest_api_get + '/' + g
            input_queue.put((url, genes_file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < num_genes:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {num_genes} KEGG GENES files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        self.progress.end()
        for worker in workers:
            worker.terminate()
        if undownloaded:
            if dont_raise:
                self.run.warning(f"Files for the following KEGG GENES failed to download despite multiple attempts: "
                                f"{', '.join(undownloaded)}. If this is unacceptable to you, you can try to "
                                f"re-run this program to see if things will work on the next try.")
            else:
                raise ConfigError(f"Files for the following KEGG GENES failed to download despite multiple attempts: "
                                  f"{', '.join(undownloaded)}. Since the function responsible for handling this was "
                                  f"told to quit should this happen, well, here we are. If skipping these failed KOs "
                                  f"is okay, you could always run this function with `dont_raise=True`.")

        return undownloaded


    def get_kegg_gene_accessions_from_ko_files(self, ko_list, ko_file_dir):
        """Extracts KEGG GENES accessions from KO files and returns a dictionary mapping KO to its GENES.

        Parameters
        ==========
        ko_list: list of str
            List of KEGG accessions to process.
        ko_file_dir: file path
            Where the KO files are located
        Returns
        =======
        ko_to_genes : dict
            Dictionary with KOs as keys and list of KEGG GENES accessions as values
        """

        ko_to_genes = {k : [] for k in ko_list}

        for ko in ko_list:
            ko_file_path = os.path.join(ko_file_dir, ko)
            genes_acc_list = self.extract_data_field_from_kegg_file(ko_file_path, target_field="GENES")

            kegg_genes_code_list = []
            for i, acc in enumerate(genes_acc_list):
                acc_fields = acc.split(": ")            # example accession is "CTC: CTC_p60(tetX)"
                org_code = acc_fields[0].lower()        # the organism code (before the colon) needs to be converted to lowercase
                gene_name = acc_fields[1].split('(')[0] # the gene name (after the colon) needs to have anything in parentheses removed

                # sometimes we have multiple genes per organism, like this: "PSOM: 113322169 113340172"
                if ' ' in gene_name:
                    all_genes = gene_name.split(' ')
                    for g in all_genes:
                        kegg_genes_code = f"{org_code}:{g}"
                        kegg_genes_code_list.append(kegg_genes_code)
                else:
                    kegg_genes_code_list.append(f"{org_code}:{gene_name}")
            ko_to_genes[ko] = kegg_genes_code_list

        return ko_to_genes


    def kegg_gene_sequences_to_fasta_file(self, kegg_genes_files, target_fasta_file):
        """This function extracts the amino acid sequences for a list of KEGG GENES and prints them to a FASTA file.

        Parameters
        ==========
        kegg_genes_files : List of str
            List of paths to KEGG GENES file to extract sequences from
        target_fasta_file : list of str
            Path to FASTA file in which to store the sequences

        Returns
        =======
        seq_tuples : List of tuples
            Each sequence added to the FASTA file is also returned in this list, where each tuple contains
            (KEGG GENES name, amino acid sequence). Note that the seq name is taken from the name of the KEGG GENES file.
        """

        seq_tuples = []
        for i, gene_file_path in enumerate(kegg_genes_files):
            seq_name = os.path.basename(gene_file_path)
            # obtain the amino acid sequence and save it to the fasta file
            aa_sequence_data = self.extract_data_field_from_kegg_file(gene_file_path, target_field="AASEQ")

            aaseq = ""
            with open(target_fasta_file, 'a') as fasta:
                fasta.write(f">{i}\n") # we label the gene with its index because the HMMER parser expects an int, not a string, as the gene name
                for seq in aa_sequence_data[1:]: # we skip the first element, which is the sequence length
                    fasta.write(f"{seq}\n")
                    aaseq += seq
            seq_tuples.append((seq_name, aaseq))

        return seq_tuples


    def build_HMM_from_seqs(self, hmm_name, tuple_of_seqs, hmm_output_file, log_file_path):
        """This function aligns sequences and builds an HMM from them using `muscle` and `hmmbuild`.

        Parameters
        ==========
        hmm_name : str
            What to name the model (ie 'NAME' field in the .hmm file)
        tuple_of_seqs : List of (sequence name, sequence) tuples
            The sequences to align with 'muscle' to create the `hmmbuild` input.
            See anvio.drivers.muscle for example format
        hmm_output_file : str
            File path where to store the new HMM model
        log_file_path : str
            File path for the log file of `hmmbuild`
        """

        if len(tuple_of_seqs) < 2:
            raise ConfigError(f"The function build_HMM_from_seqs() can't build an alignment from less than "
                              f"2 sequences, but that is what it got. Here is the sequence (if any) passed to this "
                              f"function: {tuple_of_seqs}. No alignment, no HMM. Sorry!")

        m = Muscle(progress=progress_quiet, run=run_quiet)
        clw_alignment = m.run_stdin(tuple_of_seqs, debug=anvio.DEBUG, clustalw_format=True)

        hmmbuild_cmd_line = ['hmmbuild', '-n', hmm_name, '--informat', 'clustallike', hmm_output_file, '-'] # sending '-' in place of an alignment file so it reads from stdin
        utils.run_command_STDIN(hmmbuild_cmd_line, log_file_path, clw_alignment)

        if not os.path.exists(hmm_output_file):
            raise ConfigError(f"It seems that the `hmmbuild` command failed because there is no output model at {hmm_output_file}. "
                              f"Perhaps the log file {log_file_path} will hold some answers for you.")


    def estimate_bitscore_for_ko(self, ko, kegg_genes_for_ko, kegg_genes_fasta, ko_model_file):
        """This function estimates the bitscore of a single KEGG Ortholog.

        It runs `hmmscan` of the KO model against the provided list of its KEGG GENE
        sequences, and then computes the minimum bit score to use as a threshold for
        annotating this protein family.

        Parameters
        ==========
        ko : str
            KEGG identifier for the KO
        kegg_genes_for_ko : list of str
            List of KEGG GENE accessions that were used to generate the KO model (for sanity check and
            number of sequences)
        kegg_genes_fasta : str
            Path to FASTA file where the KEGG GENES sequences for this KO are stored
        ko_model_file : str
            File path of the .hmm file containg the KO model (doesn't need to contain only this model,
            but must be hmmpressed already)
        Returns
        =======
        threshold : float
            estimated bit score threshold for the KO's HMM. Will be None if kegg_genes_for_ko is empty
        """

        # sanity check for empty KEGG GENES list
        if not kegg_genes_for_ko:
            self.run.warning(f"The function estimate_bitscore_for_ko() received an empty list of KEGG GENES "
                             f"for {ko}, so it cannot estimate a bit score threshold. The function will return "
                             f"a threshold of `None` for this KO.")
            return None

        # we run hmmscan of the KO against its associated GENES sequences and process the hits
        target_file_dict = {'AA:GENE': kegg_genes_fasta}
        hmmer = HMMer(target_file_dict, num_threads_to_use=self.num_threads, progress=progress_quiet, run=run_quiet)
        hmm_hits_file = hmmer.run_hmmer('KO {ko}', 'AA', 'GENE', None, None, len(kegg_genes_for_ko), ko_model_file, None, None)

        if not hmm_hits_file:
            raise ConfigError(f"No HMM hits were found for the KO model {ko}. This is seriously concerning, because we were running it against "
                              f"gene sequences that were used to generate the model.")

        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', run=run_quiet)
        search_results_dict = parser.get_search_results()

        # take the minimum of hits from current KO model as bit score threshold
        all_relevant_bitscores = []
        for hit, hit_info_dict in search_results_dict.items():
            if hit_info_dict['gene_name'] == ko or hit_info_dict['gene_name'] == f"{ko}{STRAY_KO_ANVIO_SUFFIX}":
                all_relevant_bitscores.append(hit_info_dict['bit_score'])

        threshold = min(all_relevant_bitscores)
        return threshold


    def process_all_stray_kos(self):
        """This driver function processes each stray KO and creates a file of bit score thresholds for them.

        The following steps are run for each stray KO:
        1. download of its KO file
        2. identification and download of the KEGG GENES sequences in this family
        3. alignment with `muscle` and `hmmbuild` to create a new model (since KEGG GENES updates faster than KOfam models do)
        4. `hmmscan` of the new KO model against these sequences to get bit scores
        5. computing the minimum bit score to use as a threshold for annotating this family
        """

        num_strays = len(self.ko_no_threshold_list)
        self.run.info("Number of Stray KOs to process", num_strays)

        self.stray_ko_file_dir = os.path.join(self.orphan_data_dir, "00_STRAY_KO_FILES")
        self.stray_ko_genes_dir = os.path.join(self.orphan_data_dir, "01_STRAY_GENES_FILES")
        self.stray_ko_seqs_dir = os.path.join(self.orphan_data_dir, "02_STRAY_GENES_FASTA")
        self.stray_ko_hmms_dir = os.path.join(self.orphan_data_dir, "03_STRAY_KO_HMMS")
        filesnpaths.gen_output_directory(self.stray_ko_file_dir, delete_if_exists=True)
        filesnpaths.gen_output_directory(self.stray_ko_genes_dir, delete_if_exists=True)
        filesnpaths.gen_output_directory(self.stray_ko_seqs_dir, delete_if_exists=True)
        filesnpaths.gen_output_directory(self.stray_ko_hmms_dir, delete_if_exists=True)

        ko_files_not_downloaded = self.download_ko_files(self.ko_no_threshold_list, self.stray_ko_file_dir)

        ko_files_to_process = list(set(self.ko_no_threshold_list) - set(ko_files_not_downloaded))
        self.run.info("Number of Stray KO files successfully downloaded", len(ko_files_to_process))
        if ko_files_not_downloaded:
            self.run.warning(f"FYI, some stray KOs failed to download from KEGG. We will not estimate their bit score thresholds "
                             f"and you will not be able to annotate them later. Here they are: {', '.join(ko_files_not_downloaded)}")

        ko_to_gene_accessions = self.get_kegg_gene_accessions_from_ko_files(ko_files_to_process, self.stray_ko_file_dir)
        kegg_genes_to_download = []
        for acc_list in ko_to_gene_accessions.values():
            kegg_genes_to_download.extend(acc_list)

        kegg_genes_not_downloaded = self.download_kegg_genes_files(kegg_genes_to_download, self.stray_ko_genes_dir)
        kegg_genes_downloaded = list(set(kegg_genes_to_download) - set(kegg_genes_not_downloaded))
        self.run.info("Number of KEGG GENES files successfully downloaded", len(kegg_genes_downloaded))
        if kegg_genes_not_downloaded:
            self.run.warning(f"FYI, some KEGG GENES files failed to download from KEGG. They will not be used for "
                             f"estimating bit score thresholds. Here they are: {', '.join(kegg_genes_not_downloaded)}")

        self.progress.new("Extracting amino acid sequences for Stray KOs", progress_total_items=len(ko_files_to_process))
        ko_to_gene_seqs_list = {} # we'll store the sequences to align with muscle here. yes, we just stored them in a file.
        cur_num = 0
        for k in ko_files_to_process:
            self.progress.update(f"Working on {k} [{cur_num} of {len(ko_files_to_process)}]")
            self.progress.increment(increment_to=cur_num)
            downloaded_genes_list = [a for a in ko_to_gene_accessions[k] if a in kegg_genes_downloaded]
            gene_file_paths = [os.path.join(self.stray_ko_genes_dir, code) for code in downloaded_genes_list]
            ko_to_gene_seqs_list[k] = self.kegg_gene_sequences_to_fasta_file(gene_file_paths, os.path.join(self.stray_ko_seqs_dir, f"GENES_FOR_{k}.fa"))
            cur_num += 1
        self.progress.end()

        self.progress.new("Aligning genes and creating new HMMs for Stray KOs", progress_total_items=len(ko_files_to_process))
        list_of_new_HMMs = []
        hmmbuild_log = os.path.join(self.orphan_data_dir, "hmmbuild.log")
        cur_num = 0
        new_models = 0
        old_models = 0
        models_without_genes = []
        kos_with_one_gene = []
        models_with_anvio_version = []
        for k in ko_files_to_process:
            self.progress.update(f"Working on {k} [{cur_num} of {len(ko_files_to_process)}]")
            self.progress.increment(increment_to=cur_num)
            if not len(ko_to_gene_seqs_list[k]):
                models_without_genes.append(k)
            elif len(ko_to_gene_seqs_list[k]) == 1:
                # with only 1 sequence, we can't build a new model. We can try to use KEGG's since it is guaranteed to fit this sequence
                kegg_model_file = os.path.join(self.orphan_data_dir, f"profiles/{k}.hmm")
                if not os.path.exists(kegg_model_file):
                    kos_with_one_gene.append(k) # if we don't have the OG model, we just skip this one
                else:
                    list_of_new_HMMs.append(kegg_model_file)
                    old_models += 1
            else:
                hmm_model_file = os.path.join(self.stray_ko_hmms_dir, f"{k}_anvio.hmm")
                self.build_HMM_from_seqs(f"{k}{STRAY_KO_ANVIO_SUFFIX}", ko_to_gene_seqs_list[k], hmm_model_file, hmmbuild_log)
                list_of_new_HMMs.append(hmm_model_file)
                new_models += 1
                models_with_anvio_version.append(k)
            cur_num += 1
        self.progress.end()

        self.run.info("Number of Stray KOs with new HMMs built by anvi'o to incorporate potentially new KEGG GENES", new_models)
        self.run.info("Number of Stray KOs using KEGG's original HMM because the family includes only one gene sequence", old_models)
        if models_without_genes:
            self.run.warning(f"We weren't able to download any KEGG GENE sequences for some stray KOs, and therefore will not "
                             f"be able to estimate bit score threshold for these KOs. Here they are: {', '.join(models_without_genes)}")
            self.run.info("Number of KOs without downloaded gene sequences", len(models_without_genes))
            ko_files_to_process = list(set(ko_files_to_process) - set(models_without_genes))
        if kos_with_one_gene:
            self.run.warning(f"The following stray KOs had exactly one KEGG GENE sequence, so we couldn't build a new HMM for them, "
                             f"but we also couldn't find their models from KEGG, so we won't estimate bit score thresholds for them: "
                             f"{', '.join(kos_with_one_gene)}")
            self.run.info("Number of KOs without HMMs", len(kos_with_one_gene))
            ko_files_to_process = list(set(ko_files_to_process) - set(kos_with_one_gene))

        self.progress.new('Concatenating new Stray HMM files...')
        utils.concatenate_files(self.stray_ko_hmm_file_path, list_of_new_HMMs, remove_concatenated_files=False)
        self.progress.update('Running hmmpress on new Stray KO HMMs...')
        self.exec_hmmpress_command_on_ko_file(self.stray_ko_hmm_file_path, os.path.join(self.orphan_data_dir, '00_hmmpress_log.txt'))
        self.progress.end()
        self.run.info("File storing all new HMMs generated for Stray KOs", self.stray_ko_hmm_file_path)

        self.progress.new("Estimating bit score threshold for Stray KOs", progress_total_items=len(ko_files_to_process))
        threshold_dict = {}
        cur_num = 0
        for k in ko_files_to_process:
            self.progress.update(f"Working on {k} [{cur_num} of {len(ko_files_to_process)}]")
            self.progress.increment(increment_to=cur_num)
            downloaded_genes_list = [a for a in ko_to_gene_accessions[k] if a in kegg_genes_downloaded]
            threshold_dict[k] = self.estimate_bitscore_for_ko(k, kegg_genes_for_ko=downloaded_genes_list,
                                        kegg_genes_fasta=os.path.join(self.stray_ko_seqs_dir, f"GENES_FOR_{k}.fa"),
                                        ko_model_file=self.stray_ko_hmm_file_path)
            cur_num += 1
        self.progress.end()

        # we need to re-load the ko dictionary so that we have access to the definitions of the stray KOs
        # cannot do this before this point because the absence of an stray KO from this dict controls whether it is moved to the
        # stray data directory (and we want to keep the strays separate since we process them specially)
        self.setup_ko_dict(exclude_threshold=(not self.include_stray_kos), suppress_warnings=True)

        # write the thresholds to a file
        thresholds_not_none = 0
        with open(self.stray_ko_thresholds_file, 'w') as out:
            out.write("knum\tthreshold\tscore_type\tdefinition\n")
            for k, t in threshold_dict.items():
                if t:
                    model_name = k
                    if k in models_with_anvio_version:
                        model_name = f"{k}{STRAY_KO_ANVIO_SUFFIX}"
                    ko_definition = self.ko_dict[k]['definition']
                    out.write(f"{model_name}\t{t}\tfull\t{ko_definition}\n")
                    thresholds_not_none += 1
        self.run.info("File with estimated bit score thresholds", self.stray_ko_thresholds_file)
        self.run.info("Number of estimated thresholds", thresholds_not_none)

        # clean up downloaded files
        if not anvio.DEBUG:
            os.remove(hmmbuild_log)
            for d in [self.stray_ko_file_dir, self.stray_ko_genes_dir, self.stray_ko_seqs_dir, self.stray_ko_hmms_dir]:
                shutil.rmtree(d)
            self.run.warning("The KO and GENES files downloaded from KEGG for processing the stray KOs, as well as the "
                             "individual new HMM files that anvi'o generated for these KOs, are now deleted to save space. "
                             "If you want to keep them, next time run the program with `--debug`.")


    def setup_kofams(self):
        """This function downloads, decompresses, and runs `hmmpress` on KOfam profiles."""

        if not self.only_processing:
            self.download_profiles()

        if not self.only_download:
            self.decompress_profiles()
            self.setup_ko_dict() # get ko dict attribute
            self.run_hmmpress()

            if self.include_stray_kos:
                self.process_all_stray_kos()

            # there is no reason to keep the original HMM profiles around, unless we are debugging
            if not anvio.DEBUG:
                shutil.rmtree(os.path.join(self.kegg_data_dir, "profiles"))
                shutil.rmtree(os.path.join(self.orphan_data_dir, "profiles"))


class ModulesDownload(KeggSetup):
    """Class for setting up all KEGG data related to pathway prediction, namely KOfam profiles and KEGG MODULES;
    reaction networks, which require MODULES, BRITE, and binary relation files;
    and pathway map images and reference KO, EC, and RN KGML files.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to command-line programs relying on this
        class, such as `anvi-setup-kegg-data`. If using this class through the API, please
        provide a Namespace object with the Boolean 'reset' parameter.
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories
        when testing this class.
    """

    def __init__(self, args, run=run, progress=progress, skip_init=False):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress
        self.skip_init = skip_init
        self.skip_brite_hierarchies = A('skip_brite_hierarchies')
        self.skip_binary_relations = A('skip_binary_relations')
        self.skip_map_images = A('skip_map_images')
        self.overwrite_modules_db = A('overwrite_output_destinations')

        self.run.info_single("Info from MODULES Download")

        # we also need the init of the superclass
        KeggSetup.__init__(self, self.args, skip_init=self.skip_init)

        if (not self.download_from_kegg) and self.skip_brite_hierarchies:
            self.run.warning("Just so you know, the --skip-brite-hierarchies flag does not do anything (besides suppress some warning output) when used "
                             "without the -D option. You are setting up from an archived KEGG snapshot which may already include BRITE data, and if it "
                             "does, this data will not be removed. You can always check if the resulting modules database contains BRITE data by "
                             "running `anvi-db-info` on it and looking at the `is_brite_setup` value (which will be 1 if the database contains BRITE data).")

        if (not self.download_from_kegg) and self.skip_binary_relations:
            self.run.warning(
                "Just so you know, the --skip-binary-relations flag does not do anything (besides "
                "suppress some warning output) when used without the -D option. You are setting up "
                "from an archived KEGG snapshot which may already include binary relation files, "
                "and if it does, this data will not be removed. `anvi-reaction-network` depends on "
                "these files and will let you know if they're missing."
            )

        # download from KEGG option: module/pathway map htext files and API link
        self.kegg_module_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
        self.kegg_pathway_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir="
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"
        self.kegg_binary_relations_download_path = "https://www.genome.jp/kegg-bin/show?file="
        # download a json file containing all BRITE hierarchies, which can then be downloaded themselves
        self.kegg_brite_hierarchies_download_path = os.path.join(self.kegg_rest_api_get, "br:br08902/json")
        # download the list of pathways, used for processing map image files
        self.kegg_pathway_list_download_path = "https://rest.kegg.jp/list/pathway"

        # check if the data is already downloaded
        expected_files_for_modules = [self.kegg_module_file,
                                      self.kegg_module_data_dir]
        if not self.skip_brite_hierarchies:
            expected_files_for_modules.append(self.kegg_brite_hierarchies_file)
            expected_files_for_modules.append(self.brite_data_dir)
        if not self.skip_binary_relations:
            expected_files_for_modules.append(self.binary_relation_data_dir)
        if not self.skip_map_images:
            expected_files_for_modules.append(self.map_image_data_dir)

        if not args.reset and not anvio.DEBUG and not self.skip_init:
            self.is_database_exists(expected_files_for_modules, fail_if_exists=(not self.only_processing))

        # generate subfolders if necessary
        if self.download_from_kegg and not self.only_processing and not self.kegg_archive_path and not self.skip_init:
            filesnpaths.gen_output_directory(self.kegg_module_data_dir, delete_if_exists=args.reset)
            if not self.skip_brite_hierarchies:
                filesnpaths.gen_output_directory(self.brite_data_dir, delete_if_exists=args.reset)
            if not self.skip_binary_relations:
                filesnpaths.gen_output_directory(
                    self.binary_relation_data_dir, delete_if_exists=args.reset
                )
            if not self.skip_map_images:
                filesnpaths.gen_output_directory(
                    self.map_image_data_dir, delete_if_exists=args.reset
                )
                # Create subdirectories of the map image directory.
                for subdir in (
                    self.map_image_data_dir,
                    self.png_dir,
                    self.kgml_dir,
                    self.png_1x_dir,
                    self.png_2x_dir,
                    self.png_1x_map_dir,
                    self.png_1x_ko_dir,
                    self.png_1x_ec_dir,
                    self.png_1x_rn_dir,
                    self.png_1x_org_dir,
                    self.png_2x_map_dir,
                    self.kgml_1x_dir,
                    self.kgml_2x_dir,
                    self.kgml_1x_ko_dir,
                    self.kgml_1x_ec_dir,
                    self.kgml_1x_rn_dir,
                    self.kgml_1x_org_dir,
                    self.kgml_2x_ko_dir,
                    self.kgml_2x_ec_dir,
                    self.kgml_2x_rn_dir,
                    self.kgml_2x_org_dir
                ):
                    filesnpaths.gen_output_directory(subdir)

    def download_kegg_module_file(self):
        """This function downloads the KEGG module file, which tells us which module files to download."""

        # download the kegg module file, which lists all modules
        try:
            utils.download_file(self.kegg_module_download_path, self.kegg_module_file, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download the KEGG Module htext file from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")


    def process_module_file(self):
        """This function reads the kegg module file into a dictionary. It should be called during setup to get the KEGG module numbers so that KEGG modules can be downloaded.

        The structure of this file is like this:

        +D    Module
        #<h2><a href="/kegg/kegg2.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a>&nbsp; KEGG Modules</h2>
        !
        A<b>Pathway modules</b>
        B
        B  <b>Carbohydrate metabolism</b>
        C    Central carbohydrate metabolism
        D      M00001  Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate [PATH:map00010 map01200 map01100]
        D      M00002  Glycolysis, core module involving three-carbon compounds [PATH:map00010 map01200 map01230 map01100]
        D      M00003  Gluconeogenesis, oxaloacetate => fructose-6P [PATH:map00010 map00020 map01100]

        In other words, a bunch of initial lines to be ignored, and thereafter the line's information can be determined by the one-letter code at the start.
        A = Pathway modules (metabolic pathways) or signature modules (gene sets that indicate a phenotypic trait, ie toxins).
        B = Category of module (a type of metabolism for pathway modules. For signature modules, either Gene Set or Module Set)
        C = Sub-category of module
        D = Module

        """
        self.module_dict = {}

        filesnpaths.is_file_exists(self.kegg_module_file)
        filesnpaths.is_file_plain_text(self.kegg_module_file)

        f = open(self.kegg_module_file, 'r')
        self.progress.new("Parsing KEGG Module file")

        current_module_type = None
        current_category = None
        current_subcategory = None

        for line in f.readlines():
            line = line.strip('\n')
            first_char = line[0]

            # garbage lines
            if first_char in ["+", "#", "!"]:
                continue
            else:
                # module type
                if first_char == "A":
                    fields = re.split('<[^>]*>', line) # we split by the html tag here
                    current_module_type = fields[1]
                # Category
                elif first_char == "B":
                    fields = re.split('<[^>]*>', line) # we split by the html tag here
                    if len(fields) == 1: # sometimes this level has lines with only a B
                        continue
                    current_category = fields[1]
                # Sub-category
                elif first_char == "C":
                    fields = re.split('\s{2,}', line) # don't want to split the subcategory name, so we have to split at least 2 spaces
                    current_subcategory = fields[1]
                # module
                elif first_char == "D":
                    fields = re.split('\s{2,}', line)
                    mnum = fields[1]
                    self.module_dict[mnum] = {"name" : fields[2], "type" : current_module_type, "category" : current_category, "subcategory" : current_subcategory}
                # unknown code
                else:
                    raise ConfigError("While parsing the KEGG file %s, we found an unknown line code %s. This has "
                                      "made the file unparseable. It is likely that an update to KEGG has broken "
                                      "things such that anvi'o doesn't know what is going on anymore. Sad, we know. :( "
                                      "Please contact the developers to see if this is a fixable issue, and in the "
                                      "meantime use an older version of the KEGG data directory (if you have one). "
                                      "If we cannot fix it, we may be able to provide you with a legacy KEGG "
                                      "data archive that you can use to setup KEGG with the --kegg-archive flag." % (self.kegg_module_file, first_char))
        self.progress.end()


    def download_modules(self):
        """This function downloads the KEGG modules."""

        total = len(self.module_dict.keys())
        self.run.info("KEGG Module Database URL", self.kegg_rest_api_get)
        self.run.info("Number of KEGG Modules to download", total)
        self.run.info("Number of threads used for download", self.num_threads)

        self.progress.new("Downloading KEGG Module files", progress_total_items=total)
        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for mnum in self.module_dict.keys():
            file_path = os.path.join(self.kegg_module_data_dir, mnum)
            url = self.kegg_rest_api_get + '/' + mnum
            input_queue.put((url, file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < total:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {total} module files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        for worker in workers:
            worker.terminate()
        self.progress.end()

        if undownloaded:
            raise ConfigError(
                "Unfortunately, files for the following modules failed to download despite multiple attempts, "
                f"and so the database needs to be set up again: {', '.join(undownloaded)}"
            )


    def confirm_downloaded_modules(self):
        """This function verifies that all module files have been downloaded.

        It checks that there is a module file for every module in the self.module_dict dictionary;
        for that reason, it must be called after the function that creates that attribute,
        process_module_file(), has already been called. To verify that each file has been downloaded
        properly, we check that the last line is '///'.
        """

        for mnum in self.module_dict.keys():
            file_path = os.path.join(self.kegg_module_data_dir, mnum)
            if not os.path.exists(file_path):
                raise ConfigError(f"The module file for {mnum} does not exist at its expected location, {file_path}. "
                                  f"This probably means that something is wrong with your downloaded data, since this "
                                  f"module is present in the KEGG MODULE file that lists all modules you *should* have "
                                  f"on your computer. Very sorry to tell you this, but you need to re-download the KEGG "
                                  f"data. We recommend the --reset flag.")
            # verify entire file has been downloaded
            f = open(file_path, 'r')
            f.seek(0, os.SEEK_END)
            f.seek(f.tell() - 4, os.SEEK_SET)
            last_line = f.readline().strip('\n')
            if not last_line == '///':
                raise ConfigError("The KEGG module file %s was not downloaded properly. We were expecting the last line in the file "
                                  "to be '///', but instead it was %s. Formatting of these files may have changed on the KEGG website. "
                                  "Please contact the developers to see if this is a fixable issue. If it isn't, we may be able to "
                                  "provide you with a legacy KEGG data archive that you can use to setup KEGG with the --kegg-archive flag."
                                  % (file_path, last_line))
        self.run.info("Number of module files found", len(self.module_dict))


    def setup_modules_data(self):
        """This is a driver function which executes the setup process for pathway prediction and reaction network data from KEGG."""

        # FIXME: we will have to move user setup to a completely separate program at some point
        # PS. user setup related functions belong to the superclass for now
        if self.user_input_dir:
            self.setup_user_data()
        else:
            # download the data first
            # unless user requested only processing (mostly for developers and the adventurous)
            if not self.only_processing:
                self.download_kegg_module_file()
                self.process_module_file() # get module dict attribute
                self.download_modules()
                self.confirm_downloaded_modules()

                if not self.skip_brite_hierarchies:
                    self.download_brite_hierarchy_of_hierarchies()
                    self.process_brite_hierarchy_of_hierarchies() # get brite dict attribute
                    self.download_brite_hierarchies()
                    self.confirm_downloaded_brite_hierarchies()

                if not self.skip_binary_relations:
                    self.download_binary_relations()
                    self.confirm_downloaded_binary_relations()

                if not self.skip_map_images:
                    self.download_map_images()
            else:
                # get required attributes for database setup and make sure all expected files were downloaded
                self.process_module_file()
                self.confirm_downloaded_modules()

                if not self.skip_brite_hierarchies:
                    self.process_brite_hierarchy_of_hierarchies()
                    self.confirm_downloaded_brite_hierarchies()

                if not self.skip_binary_relations:
                    self.confirm_downloaded_binary_relations()

            # process the modules file into a database
            if not self.only_download:
                self.setup_modules_db(db_path=self.kegg_modules_db_path, module_data_directory=self.kegg_module_data_dir, brite_data_directory=self.brite_data_dir, skip_brite_hierarchies=self.skip_brite_hierarchies)


    ###### BRITE-related functions below ######
    def download_brite_hierarchy_of_hierarchies(self):
        """Download a json file of 'br08902', a "hierarchy of BRITE hierarchies."

        This hierarchy contains the names of other hierarchies which are subsequently used for
        downloading those hierarchy json files.
        """

        # note that this is the same as the REST API for modules and pathways - perhaps at some point this should be printed elsewhere so we don't repeat ourselves.
        self.run.info("KEGG BRITE Database URL", self.kegg_rest_api_get)

        try:
            utils.download_file(self.kegg_brite_hierarchies_download_path, self.kegg_brite_hierarchies_file, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download the KEGG BRITE hierarchies json file from the KEGG website. "
                              "Something likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")


    def process_brite_hierarchy_of_hierarchies(self):
        """Read the KEGG BRITE 'br08902' 'hierarchy of hierarchies' json file into a dictionary.

        This method is called during setup to find all BRITE hierarchies to be downloaded.
        Hierarchies of interest have accessions starting with 'ko' and classify genes/proteins.
        Excluded hierarchies include those for modules, pathways, and other systems for reactions,
        compounds, taxa, etc.

        The dictionary that is filled out, `self.brite_dict`, is keyed by the 'ko' hierarchy name
        exactly as given in the 'br08902' json file. The values are the categorizations of the
        hierarchy in 'br08902', going from most general to most specific category.

        Here is an example of an entry produced in self.brite_dict:
            'ko01000  Enzymes':
                ['Genes and Proteins', 'Protein families: metabolism']
        """

        filesnpaths.is_file_exists(self.kegg_brite_hierarchies_file)
        filesnpaths.is_file_json_formatted(self.kegg_brite_hierarchies_file)

        self.progress.new("Parsing KEGG BRITE Hierarchies file")

        brite_hierarchies_dict = json.load(open(self.kegg_brite_hierarchies_file))
        # store the names of all of the 'ko' hierarchies for genes/proteins
        self.brite_dict = {}
        hierarchies_appearing_multiple_times = []
        hierarchies_with_unrecognized_accession = []
        for hierarchy, categorizations in self.invert_brite_json_dict(brite_hierarchies_dict).items():
            # we have observed the hierarchy label to have an accession followed by two spaces followed by the hierarchy name,
            # but accommodate the possibility that the accession is separated from the name by a variable number of spaces
            split_hierarchy = hierarchy.split(' ')
            hierarchy_accession = split_hierarchy[0]
            hierarchy_name = ' '.join(split_hierarchy[1: ]).lstrip()
            if hierarchy_accession[: 2] == 'br':
                # hierarchy accessions beginning with 'br' are for reactions, compounds, taxa, etc., not genes/proteins
                continue
            elif hierarchy_accession == 'ko00002' and hierarchy_name == 'KEGG modules':
                # this hierarchy is for modules, not genes/proteins
                continue
            elif hierarchy_accession == 'ko00003' and hierarchy_name == 'KEGG reaction modules':
                # this hierarchy is also for modules
                continue

            if len(categorizations) > 1:
                hierarchies_appearing_multiple_times.append((hierarchy, len(categorizations)))

            if hierarchy_accession[: 2] != 'ko':
                hierarchies_with_unrecognized_accession.append(hierarchy)
                continue
            try:
                int(hierarchy_accession[2: 7])
            except ValueError:
                hierarchies_with_unrecognized_accession.append(hierarchy)
                continue
            self.brite_dict[hierarchy] = categorizations[0][1: ]

        error_first_part = ""
        if hierarchies_appearing_multiple_times:
            error_first_part = ("Each BRITE hierarchy should only appear once in the hierarchy of hierarchies, "
                                "but the following hierarchies appeared the given number of times: "
                                f"{', '.join([f'{hier}: {num_times}' for hier, num_times in hierarchies_appearing_multiple_times])}.")
        error_second_part = ""
        if hierarchies_with_unrecognized_accession:
            error_second_part = ("Each BRITE hierarchy accession is expected to have an accession formatted 'koXXXXX', where 'XXXXX' are five digits, "
                                 f"but the following hierarchies did not have this format: {', '.join(hierarchies_with_unrecognized_accession)}.")
        if hierarchies_appearing_multiple_times or hierarchies_with_unrecognized_accession:
            raise ConfigError("Please contact the developers to look into the following error. "
                              f"{error_first_part}{' ' if error_first_part and error_second_part else ''}{error_second_part}")

        self.progress.end()


    def download_brite_hierarchies(self):
        """This function downloads a json file for every BRITE hierarchy of interest.

        Hierarchies of interest classify genes/proteins and have accessions starting with 'ko'.
        """

        total = len(self.brite_dict)
        self.run.info("Number of BRITE hierarchies to download", total)
        self.progress.new("Downloading BRITE files", progress_total_items=total)
        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        unexpected_hierarchies = []
        for hierarchy in self.brite_dict:
            hierarchy_accession = hierarchy[: 7]
            brite_system = hierarchy_accession[: 2]
            if brite_system != 'ko':
                unexpected_hierarchies.append(hierarchy)
            if not unexpected_hierarchies:
                file_path = os.path.join(self.brite_data_dir, hierarchy_accession)
                url = self.kegg_rest_api_get + '/br:' + hierarchy_accession + '/json'
                input_queue.put((url, file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < total:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {total} files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        for worker in workers:
            worker.terminate()
        if undownloaded:
            raise ConfigError(
                "Unfortunately, files for the following BRITE hierarchies failed to download despite multiple attempts, "
                f"and so the database needs to be set up again: {', '.join(undownloaded)}"
            )
        self.progress.end()

        if unexpected_hierarchies:
            raise ConfigError("Accessions for BRITE hierarchies of genes/proteins should begin with 'ko'. "
                              f"Hierarchies were found that defy our assumptions; please contact a developer to investigate this: '{', '.join(unexpected_hierarchies)}'.")


    def confirm_downloaded_brite_hierarchies(self):
        """This function verifies that all BRITE hierarchy files have been downloaded.

        It checks that there is a hierarchy file for every hierarchy in the self.brite_dict dictionary;
        for that reason, it must be called after the function that creates that attribute,
        process_brite_hierarchy_of_hierarchies(), has already been called.
        """

        for hierarchy in self.brite_dict.keys():
            hierarchy_accession = hierarchy[: 7]
            file_path = os.path.join(self.brite_data_dir, hierarchy_accession)
            if not os.path.exists(file_path):
                raise ConfigError(f"The BRITE hierarchy file for {hierarchy} does not exist at its expected location, {file_path}. "
                                  f"This probably means that something is wrong with your downloaded data, since this "
                                  f"hierarchy is present in the file that lists all BRITE hierarchies you *should* have "
                                  f"on your computer. Very sorry to tell you this, but you need to re-download the KEGG "
                                  f"data. We recommend the --reset flag.")
            # verify that the whole json file was downloaded
            filesnpaths.is_file_json_formatted(file_path)
        self.run.info("Number of BRITE hierarchy files found", len(self.brite_dict))


    ###### Binary relations-related functions below ######
    def download_binary_relations(self):
        """
        Download binary relations files relating the accession of a type of KEGG data, such as KOs,
        to related accessions of another type of data, such as EC numbers.
        """
        for file in self.kegg_binary_relation_files.values():
            url = f'{self.kegg_binary_relations_download_path}{file}'
            dest = os.path.join(self.binary_relation_data_dir, file)
            try:
                utils.download_file(url, dest, progress=self.progress, run=self.run)
            except Exception as e:
                print(e)
                raise ConfigError(
                    f"Anvi'o failed to download the KEGG binary relations file, '{file}', from the "
                    "KEGG website. Something likely changed on the KEGG end. Please contact the "
                    "developers to see if this is a fixable issue. If it isn't, we may be able to "
                    "provide you with a legacy KEGG data archive that you can use to set up KEGG "
                    "with the --kegg-archive flag."
                )


    def confirm_downloaded_binary_relations(self):
        """Verify that all expected binary relations files were downloaded."""
        missing_files = []
        for file in self.kegg_binary_relation_files.values():
            path = os.path.join(self.binary_relation_data_dir, file)
            if not os.path.exists(path):
                missing_files.append(file)
        if missing_files:
            raise ConfigError(
                "The following binary relation files were not found in the expected directory, "
                f"'{self.binary_relation_data_dir}', so the KEGG data should be re-downloaded: "
                f"{', '.join(missing_files)}"
            )
        self.run.info(
            "Number of KEGG binary relations files found", len(self.kegg_binary_relation_files)
        )


    ###### Pathway map image-related functions below ######
    def download_map_images(
        self,
        add_global_reaction_line_width: Union[float, None] = 6.0,
        global_compound_circle_diameter: Union[float, None] = 17.0
    ) -> None:
        """
        Download reference pathway map image files and associated KGML files.
        
        Only download maps with at least one reference KGML file, since the purpose is to be able to
        modify maps with data, and KGML files are required to customize maps. Write a table
        indicating which KO, EC, and RN KGML files are available for every map available in KEGG,
        including those not downloaded due to an absence of KGML files.
        
        Different sets of "global" and non-global "standard" and "overview" map images are
        downloaded. The following global map images are downloaded: 1x and 2x resolution images with
        filenames starting "map", and 1x images starting "ko", "ec", and "rn". The "ko", "ec", and
        "rn" global maps color reactions with accessions in each of the KEGG KO, EC, and RN
        databases, respectively, and the "map" global maps color reactions with accessions in any of
        these databases. Non-global 1x and 2x resolution map images starting with "map" are
        downloaded. KGML files, which are tailored to the position of features in 1x maps, are
        copied to rescale features to match 2x image files.
        
        Parameters
        ==========
        add_global_reaction_line_width : Union[float, None], 6.0
            If not None, modify downloaded global map KGML files to add a width attribute to
            reaction line graphics elements. The default value of 6 (in the 1x resolution maps, 12
            in the 2x resolution maps) is just wide enough for the lines drawn from the KGML file to
            cover up the lines in the base map image.
            
        global_compound_circle_diameter : Union[float, None], 17.0
            If not None, modify downloaded global map KGML files to adjust the size of compound
            circle graphics elements. The argument value is used as the width and height attributes
            in 1x resolution maps, with twice the value used in 2x resolution maps. The default
            value of 17 is just wide enough for the circles rendered from the KGML file to cover up
            the circles in the base map image.
        """
        # Download a table from KEGG listing all available pathways.
        try:
            utils.download_file(
                self.kegg_pathway_list_download_path,
                self.kegg_pathway_list_file,
                progress=self.progress,
                run=self.run
            )
        except Exception as e:
            print(e)
            raise ConfigError(
                "Anvi'o failed to download a list of pathways from the KEGG website. Something "
                "likely changed on the KEGG end. Please contact the developers to see if this is a "
                "fixable issue."
            )
        pathway_table = pd.read_csv(
            self.kegg_pathway_list_file, sep='\t', header=None, names=['id', 'name']
        )
        
        # Determine the maximum number of map image files that may be downloaded (image files are
        # only downloaded if a corresponding KGML file is available). 5 versions of each global map
        # are downloaded: 1x and 2x "map" files and 1x "ko", "ec", and "rn" files. 2 versions of
        # each non-global map may be downloaded: 1x and 2x "map" files.
        global_map_count = sum(
            1 if re.match(GLOBAL_MAP_ID_PATTERN, pathway_id[-5:]) else 0
            for pathway_id in pathway_table['id']
        )
        nonglobal_map_count = len(pathway_table) - global_map_count
        total_dl_count = global_map_count * 5 + nonglobal_map_count * 2
        self.run.info_single(
            f"Up to {total_dl_count} map images will be downloaded. \"Up to\" because only maps "
            f"found to have associated reference KGML files are downloaded. {self.num_threads} "
            "cores (threads) will be used in downloading.",
            nl_before=1
        )
        
        # Start the worker threads for downloading map image and KGML files.
        self.progress.new("Downloading KEGG pathway map files")
        self.progress.update("0 pathway maps downloaded")
        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for pathway_id in pathway_table['id']:
            input_queue.put({
                'pathway_number': pathway_id[3:],
                'url_stem': self.kegg_rest_api_get,
                'data_dir': self.map_image_data_dir
            })
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(
                target=_download_pathway_image_files_worker, args=(input_queue, output_queue)
            )
            workers.append(worker)
            worker.start()
            
        # Process the output of download threads. The threads should return items equal to the
        # maximum number of image files that may be downloaded. Wait for threads until this number
        # of items is reached.
        successful_dls: List[str] = []
        failed_dls: List[str] = []
        # Record the paths of KGML files that need to be rescaled to fit 2x resolution images.
        kgml_paths: List[str] = []
        # Record which of types of KGML files ('KO', 'EC', 'RN') are available for each downloaded
        # pathway map image.
        kgml_availability: Dict[str, Dict[str, int]] = {}
        processed_count = 0
        while processed_count < total_dl_count:
            # For each pathway, a dictionary is returned with keys indicating each type of possible
            # map image and KGML file that can be downloaded, and length-2 list values containing
            # 1) the possible filepath and 2) an integer value indicating the success or failure
            # type of the download.
            output: Dict[str, List[str, int]] = output_queue.get()
            pathway_id = os.path.splitext(os.path.basename(output['png_1x_map'][0]))[0]
            image_keys = ['png_1x_map', 'png_2x_map']
            if re.match(GLOBAL_MAP_ID_PATTERN, pathway_id[-5:]):
                image_keys += ['png_1x_ko', 'png_1x_ec', 'png_1x_rn']
            for image_key in image_keys:
                if output[image_key][1] == 0:
                    # This occurs when there were connection errors preventing any KGML files from
                    # being downloaded, so PNG file downloads were not attempted.
                    failed_dls.append(output[image_key][0])
                elif output[image_key][1] == 1:
                    successful_dls.append(output[image_key][0])
                    self.progress.update(f"{len(successful_dls)} pathway maps downloaded")
                elif output[image_key][1] == 2:
                    # This indicates that the PNG file was unavailable for download. It should have
                    # been available given KEGG's pathway list.
                    failed_dls.append(output[image_key][0])
                elif output[image_key][1] == 3:
                    # This occurs when connection errors prevented the PNG file from being
                    # downloaded.
                    failed_dls.append(output[image_key][0])
                elif output[image_key][1] == 4:
                    # This indicates that the program did not attempt to download the PNG file
                    # because there is no KGML file available, e.g., drug maps have no KO, EC, and
                    # RN KGML files available.
                    pass
                # Record KGML files associated with 2x resolution images. These need to be rescaled.
                if image_key == 'png_2x_map':
                    for kgml_key in ('kgml_ko', 'kgml_ec', 'kgml_rn'):
                        if output[kgml_key][1] == 1:
                            kgml_paths.append(output[kgml_key][0])
                processed_count += 1
            # Record data that goes into the table of KGML availability for each pathway.
            kgml_availability[pathway_id] = pathway_kgml_availability = {}
            for pathway_org in ('ko', 'ec', 'rn'):
                if output[f'kgml_{pathway_org}'][1] == 1:
                    pathway_kgml_availability[pathway_org.upper()] = 1
                else:
                    pathway_kgml_availability[pathway_org.upper()] = 0
                    
        # Downloading is complete. Kill the worker threads.
        for worker in workers:
            worker.terminate()
        self.progress.end()
        
        # Raise an exception when expected files failed to download. Report the failed files by
        # pathway ID.
        if failed_dls:
            failed_dl_groups: Dict[str, List[str]] = {}
            for failed_dl in failed_dls:
                failed_filename = os.path.basename(failed_dl)
                pathway_number = os.path.splitext(failed_filename)[0][3:]
                try:
                    failed_dl_groups[pathway_number].append(failed_filename)
                except KeyError:
                    failed_dl_groups[pathway_number] = [failed_filename]
            failed_message = ''
            for pathway_number, failed_filenames in failed_dl_groups.items():
                failed_message += f"map{pathway_number}: {', '.join(failed_filenames)}; "
            failed_message = failed_message[:-2]
            raise ConfigError(
                "Unfortunately, files (in parentheses) for the following pathway maps failed to "
                "download despite multiple attempts, and so the database needs to be set up again: "
                f"{failed_message}"
            )
        self.run.info("Number of downloaded map images", len(successful_dls))
        
        # Add reaction line widths to global map KGML files.
        if add_global_reaction_line_width is not None:
            self._add_global_kgml_reaction_line_widths(add_global_reaction_line_width)

        # Rescale compound circles in global map KGML files.
        if global_compound_circle_diameter is not None:
            self._change_global_kgml_compound_circle_diameters(global_compound_circle_diameter)

        # Create rescaled KGML files to fit 2x resolution map images.
        self.progress.new(
            "Creating map KGML files rescaled to 2x resolution",
            progress_total_items=len(kgml_paths)
        )
        # This import can't happen at the module level due to a circular import.
        import anvio.kgml as kgml
        xml_ops = kgml.XMLOps()
        rescaled_count = 0
        for input_path in kgml_paths:
            self.progress.update(f"{rescaled_count} / {len(kgml_paths)} KGML files rescaled")
            kgml_id: str = os.path.splitext(os.path.basename(input_path))[0]
            pathway_org = kgml_id[:-5]
            pathway_number = kgml_id[-5:]
            pathway = xml_ops.load(input_path)
            pathway.scale_graphics(2)
            if pathway_org == 'ko':
                kgml_dir = self.kgml_2x_ko_dir
            elif pathway_org == 'ec':
                kgml_dir = self.kgml_2x_ec_dir
            elif pathway_org == 'rn':
                kgml_dir = self.kgml_2x_rn_dir
            else:
                raise AssertionError(
                    "Only KGML files for pathway IDs starting with 'ko', 'ec', and 'rn' should "
                    f"have been downloaded. The ID, '{kgml_id}', is not recognized."
                )
            output_path = os.path.join(kgml_dir, f'{kgml_id}.xml')
            xml_ops.write(pathway, output_path)
            rescaled_count += 1
        self.progress.end()
        
        # Write a table of the KGML files available for each map image.
        pd.DataFrame.from_dict(
            kgml_availability, orient='index', columns=['KO', 'EC', 'RN']
        ).sort_index().to_csv(self.kegg_map_image_kgml_file, sep='\t')
        
    def _add_global_kgml_reaction_line_widths(self, width: float) -> None:
        """
        Add reaction line widths to newly downloaded KGML files for global maps. Width attributes
        are not in the files.
        
        Parameters
        ==========
        width : float
            Width value to add.
        """
        assert width > 0
        
        # This import can't happen at the module level due to a circular import.
        import anvio.kgml as kgml
        xml_ops = kgml.XMLOps()
        
        for entry_type, kgml_dir in zip(
            ('ortholog', 'enzyme', 'reaction'),
            (self.kgml_1x_ko_dir, self.kgml_1x_ec_dir, self.kgml_1x_rn_dir)
        ):
            for kgml_path in glob.glob(os.path.join(kgml_dir, '*.xml')):
                if re.match(
                    GLOBAL_MAP_ID_PATTERN, os.path.splitext(os.path.basename(kgml_path))[0][-5:]
                ):
                    pathway = xml_ops.load(kgml_path)
                    for entry in pathway.get_entries(entry_type=entry_type):
                        for uuid in entry.children['graphics']:
                            graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                            graphics.width = width
                    xml_ops.write(pathway, kgml_path)
                    
    def _change_global_kgml_compound_circle_diameters(self, diameter: float) -> None:
        """
        Change the diameters of compound circles in KGML files for global maps. The purpose of this
        is to fully cover circles in base map images with circles rendered from KGML files.
        
        Parameters
        ==========
        diameter : float
            New diameter of compound cirles.
        """
        assert diameter > 0
        
        # This import can't happen at the module level due to a circular import.
        import anvio.kgml as kgml
        xml_ops = kgml.XMLOps()
        
        for kgml_dir in (self.kgml_1x_ko_dir, self.kgml_1x_ec_dir, self.kgml_1x_rn_dir):
            for kgml_path in glob.glob(os.path.join(kgml_dir, '*.xml')):
                if re.match(
                    GLOBAL_MAP_ID_PATTERN, os.path.splitext(os.path.basename(kgml_path))[0][-5:]
                ):
                    pathway = xml_ops.load(kgml_path)
                    for entry in pathway.get_entries(entry_type='compound'):
                        for uuid in entry.children['graphics']:
                            graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                            width = graphics.width
                            height = graphics.height
                            if width is not None:
                                graphics.width = diameter
                            if height is not None:
                                graphics.height = diameter
                    xml_ops.write(pathway, kgml_path)


class RunKOfams(KeggContext):
    """Class for running `hmmscan` against the KOfam database and adding the resulting hits to contigs DB for later metabolism prediction.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-run-kegg-kofams
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.num_threads = A('num_threads')
        self.hmm_program = A('hmmer_program') or 'hmmsearch'
        self.include_stray_kos = True if A('include_stray_KOs') else False
        self.keep_all_hits = True if A('keep_all_hits') else False
        self.log_bitscores = True if A('log_bitscores') else False
        self.skip_bitscore_heuristic = True if A('skip_bitscore_heuristic') else False
        self.bitscore_heuristic_e_value = A('heuristic_e_value')
        self.bitscore_heuristic_bitscore_fraction = A('heuristic_bitscore_fraction')
        self.skip_brite_hierarchies = A('skip_brite_hierarchies')
        self.ko_dict = None # should be set up by setup_ko_dict()
        self.stray_ko_dict = None # should be set up by setup_stray_ko_dict(), if possible

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists(self.hmm_program)

        # verify that Kofam HMM profiles have been set up
        if not os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError(f"Anvi'o is unable to find any KEGG files around :/ It is likely you need to first run the program "
                              f"`anvi-setup-kegg-data` to set things up. If you already have run it, but instructed anvi'o to "
                              f"store the output to a specific directory, then instead of running `anvi-setup-kegg-data` again, "
                              f"you simply need to specify the location of the KEGG data using the flag `--kegg-data-dir`. Just for "
                              f"your information, anvi'o was looking for the KEGG data here: {self.kegg_data_dir}")

        utils.is_contigs_db(self.contigs_db_path)
        filesnpaths.is_output_file_writable(self.contigs_db_path)

        # reminder to be a good citizen
        self.run.warning("Anvi'o will annotate your database with the KEGG KOfam database, as described in "
                         "Aramaki et al (doi:10.1093/bioinformatics/btz859) When you publish your findings, "
                         "please do not forget to properly credit this work.", lc='green', header="CITATION")

        self.setup_ko_dict() # read the ko_list file into self.ko_dict
        self.run.info("Stray KOs will be annotated", self.include_stray_kos)
        if self.include_stray_kos:
            self.setup_stray_ko_dict()
            self.run.warning("Please note! Because you used the flag `--include-stray-KOs`, anvi'o will annotate "
                             "your genes with KO models that do not come with a bit score threshold defined by KEGG. "
                             "We have generated new models and estimated rather conservative thresholds for them ourselves. To learn "
                             "how we did that, please read the documentation page for `anvi-setup-kegg-data`: "
                             "https://anvio.org/help/main/programs/anvi-setup-kegg-data/#what-are-stray-kos-and-what-happens-when-i-include-them")

        # load existing kegg modules db, if one exists
        if os.path.exists(self.kegg_modules_db_path):
            self.kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, module_data_directory=self.kegg_module_data_dir, brite_data_directory=self.brite_data_dir, skip_brite_hierarchies=self.skip_brite_hierarchies, args=self.args)

            if not self.skip_brite_hierarchies and not self.kegg_modules_db.db.get_meta_value('is_brite_setup'):
                self.run.warning("The KEGG Modules database does not contain BRITE hierarchy data, "
                                 "which could very well be useful to you. BRITE is guaranteed to be set up "
                                 "when downloading the latest version of KEGG with `anvi-setup-kegg-data`.")
        else:
            self.run.warning("No modules database was found in the KEGG data directory you specified. This is fine, but "
                             "you will not get functional annotations related to KEGG MODULES or BRITE hierarchies in your "
                             "contigs database. If you want to include these annotations later, you will have to rerun this "
                             "program with a data directory including a modules database (which you can obtain by running "
                             "`anvi-setup-kegg-data` again with the right mode(s).")
            self.kegg_modules_db = None


    def check_hash_in_contigs_db(self):
        """Checks the contigs DB self table to make sure it was not already annotated"""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        current_module_hash_in_contigs_db = contigs_db.db.get_meta_value('modules_db_hash', return_none_if_not_in_table=True)

        if current_module_hash_in_contigs_db and not self.just_do_it:
            contigs_db.disconnect()
            raise ConfigError("The contigs database (%s) has already been annotated with KOfam hits. If you really want to "
                              "overwrite these annotations with new ones, please re-run the command with the flag --just-do-it. "
                              "For those who need this information, the Modules DB used to annotate this contigs database previously "
                              "had the following hash: %s" % (self.contigs_db_path, current_module_hash_in_contigs_db))


    def set_hash_in_contigs_db(self):
        """Modifies the contigs DB self table to indicate which MODULES.db has been used to annotate it."""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db')

        hash_to_add = "only_KOfams_were_annotated"
        if self.kegg_modules_db:
            hash_to_add = self.kegg_modules_db.db.get_meta_value('hash')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        contigs_db.db.set_meta_value('modules_db_hash', hash_to_add)
        contigs_db.disconnect()


    def get_annotation_from_ko_dict(self, knum, ok_if_missing_from_dict=False):
        """Returns the functional annotation of the provided KO number.

        Parameters
        ==========
        knum : str
            The KO number for which to get an annotation for
        ok_if_missing_from_dict : bool
            If false, not finding the KO will raise an error. If true, the function will quietly return an "Unknown" annotation string for the missing KO

        Returns
        =======
        annotation : str
        """

        if not self.ko_dict:
            raise ConfigError("Oops! The ko_list file has not been properly loaded, so get_annotation_from_ko_dict() is "
                              "extremely displeased and unable to function properly. Please refrain from calling this "
                              "function until after setup_ko_dict() has been called.")
        if self.include_stray_kos and not self.stray_ko_dict:
            raise ConfigError("Oops! The bit score thresholds for stray KOs have not been properly loaded, so "
                              "get_annotation_from_ko_dict() is unable to work properly. If you plan to use "
                              "--include-stray-KOs, then make sure you run the setup_stray_ko_dict() function before "
                              "calling this one.")

        ret_value = None
        if knum in self.ko_dict:
            ret_value = self.ko_dict[knum]['definition']
        elif self.include_stray_kos and knum in self.stray_ko_dict:
            ret_value = self.stray_ko_dict[knum]['definition']
        else:
            if ok_if_missing_from_dict:
                return "Unknown function with KO num %s" % knum
            else:
                raise ConfigError("It seems %s found a KO number that does not exist "
                                  "in the KOfam ko_list file: %s" % (self.hmm_program, knum))

        return ret_value


    def parse_kofam_hits(self, hits_dict, hits_label = "KOfam", next_key=None):
        """This function applies bitscore thresholding (if requested) to establish the self.functions_dict
        which can then be used to store annotations in the contigs DB.

        If self.keep_all_hits is True, all hits will be added to the self.functions_dict regardless of bitscore
        threshold.

        Note that the input hits_dict contains bitscores, but the self.functions_dict does not (because the DB
        tables do not have a column for it, at least at the time of writing this).

        PARAMETERS
        ===========
        hits_dict : dictionary
            The output from the hmmsearch parser, which should contain all hits (ie, weak hits not yet removed)
        hits_label : str
            A label for the set of hits we are working on, used to keep sets separate from each other and to enable
            us to match gene caller ids to hits in different dictionaries later
        next_key : int
            The next integer key that is available for adding functions to self.functions_dict. If None is provided,
            the keys will start at 0.

        RETURNS
        ========
        counter : int
            The number of functions added to self.functions_dict. Useful for downstream functions that want to
            add to this dictionary, since it is the next available integer key.
        """

        total_num_hits = len(hits_dict.values())
        starting_annotations_in_dict = len(self.functions_dict.keys())
        self.progress.new(f"Parsing {hits_label} hits", progress_total_items=total_num_hits)
        counter = 0
        if next_key:
            counter = next_key
        num_hits_removed = 0
        cur_num_hit = 0
        for hit_key,hmm_hit in hits_dict.items():
            cur_num_hit += 1
            knum = hmm_hit['gene_name']
            gcid = hmm_hit['gene_callers_id']
            keep = False

            if cur_num_hit % 1000 == 0:
                self.progress.update("Removing weak hits [%d of %d KOs]" % (cur_num_hit, total_num_hits))
                self.progress.increment(increment_to=cur_num_hit)

            # later, we will need to quickly access the hits for each gene call. So we map gcids to the keys in the raw hits dictionary
            if gcid not in self.gcids_to_hits_dict:
                self.gcids_to_hits_dict[gcid] = {hits_label : [hit_key]}
            else:
                if hits_label not in self.gcids_to_hits_dict[gcid]:
                    self.gcids_to_hits_dict[gcid][hits_label] = [hit_key]
                else:
                    self.gcids_to_hits_dict[gcid][hits_label].append(hit_key)

            if (knum not in self.ko_dict and (self.stray_ko_dict is not None and knum not in self.stray_ko_dict)) or \
                (knum not in self.ko_dict and self.stray_ko_dict is None):
                self.progress.reset()
                raise ConfigError(f"Something went wrong while parsing the {hits_label} HMM hits. It seems that KO "
                                  f"{knum} is not in the noise cutoff dictionary for KOs. That means we do "
                                  "not know how to distinguish strong hits from weak ones for this KO. "
                                  "Anvi'o will fail now :( Please contact a developer about this error to "
                                  "get this mess fixed. ")
            # if hit is above the bitscore threshold, we will keep it
            if knum in self.ko_dict:
                if self.ko_dict[knum]['score_type'] == 'domain':
                    if hmm_hit['domain_bit_score'] >= float(self.ko_dict[knum]['threshold']):
                        keep = True
                elif self.ko_dict[knum]['score_type'] == 'full':
                    if hmm_hit['bit_score'] >= float(self.ko_dict[knum]['threshold']):
                        keep = True
                else:
                    self.progress.reset()
                    raise ConfigError(f"The KO noise cutoff dictionary for {knum} has a strange score type which "
                                    f"is unknown to anvi'o: {self.ko_dict[knum]['score_type']}")
            elif knum in self.stray_ko_dict:
                if self.stray_ko_dict[knum]['score_type'] == 'domain':
                    if hmm_hit['domain_bit_score'] >= float(self.stray_ko_dict[knum]['threshold']):
                        keep = True
                elif self.stray_ko_dict[knum]['score_type'] == 'full':
                    if hmm_hit['bit_score'] >= float(self.stray_ko_dict[knum]['threshold']):
                        keep = True
                else:
                    self.progress.reset()
                    raise ConfigError(f"The KO noise cutoff dictionary for the stray KO {knum} has a strange score type which "
                                    f"is unknown to anvi'o: {self.stray_ko_dict[knum]['score_type']}")
            else:
                raise ConfigError(f"We cannot find KO {knum} in either self.ko_dict or in self.stray_ko_dict. This is likely a "
                                  f"problem for the developers.")

            if keep or self.keep_all_hits:
                self.functions_dict[counter] = {
                    'gene_callers_id': gcid,
                    'source': 'KOfam',
                    'accession': knum,
                    'function': self.get_annotation_from_ko_dict(knum, ok_if_missing_from_dict=True),
                    'e_value': hmm_hit['e_value'],
                }
                # later, we will need to know if a particular gene call has hits or not. So here we are just saving for each
                # gene caller id the keys for its corresponding hits in the function dictionary.
                if gcid not in self.gcids_to_functions_dict:
                    self.gcids_to_functions_dict[gcid] = [counter]
                else:
                    self.gcids_to_functions_dict[gcid].append(counter)

                # add associated KEGG module information to database
                mods = None
                if self.kegg_modules_db:
                    mods = self.kegg_modules_db.get_modules_for_knum(knum)
                    names = self.kegg_modules_db.get_module_names_for_knum(knum)
                    classes = self.kegg_modules_db.get_module_classes_for_knum_as_list(knum)

                if mods:
                    mod_annotation = "!!!".join(mods)
                    mod_class_annotation = "!!!".join(classes) # why do we split by '!!!'? Because that is how it is done in COGs. So so sorry. :'(
                    mod_name_annotation = ""

                    for mod in mods:
                        if mod_name_annotation:
                            mod_name_annotation += "!!!" + names[mod]
                        else:
                            mod_name_annotation = names[mod]

                    self.kegg_module_names_dict[counter] = {
                        'gene_callers_id': gcid,
                        'source': 'KEGG_Module',
                        'accession': mod_annotation,
                        'function': mod_name_annotation,
                        'e_value': None,
                    }
                    self.kegg_module_classes_dict[counter] = {
                        'gene_callers_id': gcid,
                        'source': 'KEGG_Class',
                        'accession': mod_annotation,
                        'function': mod_class_annotation,
                        'e_value': None,
                    }

                if self.kegg_modules_db and not self.skip_brite_hierarchies:
                    # get BRITE categorization information in the form to be added to the contigs database
                    ortholog_categorizations_dict = self.get_ortholog_categorizations_dict(knum, gcid)
                    if ortholog_categorizations_dict:
                        self.kegg_brite_categorizations_dict[counter] = ortholog_categorizations_dict

                counter += 1
            else:
                num_hits_removed += 1

        self.progress.end()
        ending_annotations_in_dict = len(self.functions_dict.keys())
        self.run.info(f"Number of weak hits removed by {hits_label} parser", num_hits_removed)
        self.run.info(f"Number of annotations added for {hits_label}", ending_annotations_in_dict - starting_annotations_in_dict)

        return counter


    def update_dict_for_genes_with_missing_annotations(self, gcids_list, super_hits_dict, next_key):
        """This function adds functional annotations for genes with missing hits to the dictionary.

        The reason this is necessary is that the bitscore thresholds can be too stringent, causing
        us to miss legitimate annotations. To find these annotations, we adopt the following heuristic:
            For every gene without a KOfam annotation, we examine all the hits with an e-value below X
            and a bitscore above Y% of the threshold. If those hits are all to a unique KO profile,
            then we annotate the gene call with that KO.

            X is self.bitscore_heuristic_e_value, Y is self.bitscore_heuristic_bitscore_fraction

        For reasons that are hopefully obvious, this function must be called after parse_kofam_hits(),
        which establishes the self.functions_dict attribute.

        PARAMETERS
        ===========
        gcids_list : list
            The list of gene caller ids in the contigs database. We will use this to figure out which
            genes have no annotations
        super_hits_dict : dictionary
            A two-level dictionary in which keys are the labels for each set of hits and values are the dictionary output
            from the hmmsearch parser, which should contain all hits from the set (ie, weak hits not yet removed)
        next_key : int
            The next integer key that is available for adding functions to self.functions_dict
        """

        self.run.warning("Anvi'o will now re-visit genes without KOfam annotations to see if potentially valid "
                         "functional annotations were missed. These genes will be annotated with a KO only if "
                         f"all KOfam hits to this gene with e-value <= {self.bitscore_heuristic_e_value} and bitscore > "
                         f"({self.bitscore_heuristic_bitscore_fraction} * KEGG threshold) are hits to the same KO. Just "
                         "so you know what is going on here. If this sounds like A Very Bad Idea to you, then please "
                         "feel free to turn off this behavior with the flag --skip-bitscore-heuristic or to change "
                         "the e-value/bitscore parameters (see the help page for more info).")

        num_annotations_added = 0
        num_stray_KOs_added = 0
        total_num_genes = len(gcids_list)
        self.progress.new("Relaxing bitscore threshold", progress_total_items=total_num_genes)

        # for each gene call, check for annotation in self.functions_dict
        current_gene_num = 0
        for gcid in gcids_list:
            current_gene_num += 1
            if current_gene_num % 1000 == 0:
                self.progress.update("Adding back decent hits [%d of %d gene calls]" % (current_gene_num, total_num_genes))
                self.progress.increment(increment_to=current_gene_num)

            if gcid not in self.gcids_to_functions_dict:
                decent_hit_kos = set()
                best_e_value = 100 # just an arbitrary positive value that will be larger than any evalue
                best_hit_key = None
                best_hit_label = None

                # if no annotation, get all hits for gene caller id from the hits dictionaries
                if gcid in self.gcids_to_hits_dict:
                    for hit_label in self.gcids_to_hits_dict[gcid]:
                        for hit_key in self.gcids_to_hits_dict[gcid][hit_label]:
                            knum = super_hits_dict[hit_label][hit_key]['gene_name']
                            ko_threshold = None
                            ko_score_type = None
                            if knum in self.ko_dict:
                                ko_threshold = float(self.ko_dict[knum]['threshold'])
                                ko_score_type = self.ko_dict[knum]['score_type']
                            elif self.include_stray_kos and knum in self.stray_ko_dict:
                                ko_threshold = float(self.stray_ko_dict[knum]['threshold'])
                                ko_score_type = self.stray_ko_dict[knum]['score_type']
                            else:
                                raise ConfigError(f"panik. the function update_dict_for_genes_with_missing_annotations() "
                                                  f"cannot find the bit score threshold for {knum}.")

                            # get set of hits that fit specified heuristic parameters
                            if ko_score_type == 'domain':
                                hit_bitscore = super_hits_dict[hit_label][hit_key]['domain_bit_score']
                            elif ko_score_type == 'full':
                                hit_bitscore = super_hits_dict[hit_label][hit_key]['bit_score']
                            if super_hits_dict[hit_label][hit_key]['e_value'] <= self.bitscore_heuristic_e_value and hit_bitscore > (self.bitscore_heuristic_bitscore_fraction * ko_threshold):
                                decent_hit_kos.add(knum)
                                # keep track of hit with lowest e-value we've seen so far
                                if super_hits_dict[hit_label][hit_key]['e_value'] <= best_e_value:
                                    best_e_value = super_hits_dict[hit_label][hit_key]['e_value']
                                    best_hit_key = hit_key
                                    best_hit_label = hit_label

                    # if unique KO, add annotation with best e-value to self.functions_dict
                    if len(decent_hit_kos) == 1:
                        best_knum = super_hits_dict[best_hit_label][best_hit_key]['gene_name']
                        ## TODO: WE NEED A GENERIC FUNCTION FOR THIS SINCE IT IS SAME AS ABOVE
                        self.functions_dict[next_key] = {
                            'gene_callers_id': gcid,
                            'source': 'KOfam',
                            'accession': best_knum,
                            'function': self.get_annotation_from_ko_dict(best_knum, ok_if_missing_from_dict=True),
                            'e_value': super_hits_dict[best_hit_label][best_hit_key]['e_value'],
                        }
                        # we may never access this downstream but let's add to it to be consistent
                        self.gcids_to_functions_dict[gcid] = [next_key]

                        # track how many of stray KOs are added back
                        if best_hit_label == "Stray KO":
                            num_stray_KOs_added += 1

                        # add associated KEGG module information to database
                        mods = None
                        if self.kegg_modules_db:
                            mods = self.kegg_modules_db.get_modules_for_knum(best_knum)
                            names = self.kegg_modules_db.get_module_names_for_knum(best_knum)
                            classes = self.kegg_modules_db.get_module_classes_for_knum_as_list(best_knum)

                        if mods:
                            mod_annotation = "!!!".join(mods)
                            mod_class_annotation = "!!!".join(classes) # why do we split by '!!!'? Because that is how it is done in COGs. So so sorry. :'(
                            mod_name_annotation = ""

                            for mod in mods:
                                if mod_name_annotation:
                                    mod_name_annotation += "!!!" + names[mod]
                                else:
                                    mod_name_annotation = names[mod]

                            self.kegg_module_names_dict[next_key] = {
                                'gene_callers_id': gcid,
                                'source': 'KEGG_Module',
                                'accession': mod_annotation,
                                'function': mod_name_annotation,
                                'e_value': None,
                            }
                            self.kegg_module_classes_dict[next_key] = {
                                'gene_callers_id': gcid,
                                'source': 'KEGG_Class',
                                'accession': mod_annotation,
                                'function': mod_class_annotation,
                                'e_value': None,
                            }

                        if self.kegg_modules_db and not self.skip_brite_hierarchies:
                            # get BRITE categorization information in the form to be added to the contigs database
                            ortholog_categorizations_dict = self.get_ortholog_categorizations_dict(best_knum, gcid)
                            if ortholog_categorizations_dict:
                                self.kegg_brite_categorizations_dict[next_key] = ortholog_categorizations_dict

                        next_key += 1
                        num_annotations_added += 1

        self.progress.end()
        self.run.info("Number of decent hits added back after relaxing bitscore threshold", num_annotations_added)
        if self.include_stray_kos:
            self.run.info("... of these, number of regular KOs is", num_annotations_added - num_stray_KOs_added)
            self.run.info("... of these, number of stray KOs is", num_stray_KOs_added)
        self.run.info("Total number of hits in annotation dictionary after adding these back", len(self.functions_dict.keys()))


    def get_ortholog_categorizations_dict(self, ortholog_accession, gene_callers_id=None):
        """Return a dictionary of ortholog BRITE categorizations.

        The dictionary is formatted to represent a row of the `gene_functions` table in the contigs database.
        """

        ortholog_brite_dict = self.kegg_modules_db.get_ortholog_brite_categorizations(ortholog_accession)
        if not ortholog_brite_dict:
            return None

        # the following explains the format of BRITE "accession" and "function" entries in the
        # table. Orthologs can be in multiple hierarchies, and can be categorized in a hierarchy
        # multiple times. Each categorization of the ortholog in a hierarchy is separated by "!!!",
        # and each category in the categorization is separated by ">>>". The hierarchy name is
        # placed at the beginning of each categorization. Perhaps the name of hierarchy "ko00001",
        # which is "KEGG Orthology (KO)", should not be placed at the beginning of categorizations
        # due to its uninformativeness, but for the sake of consistency, the format is maintained
        # for this hierarchy. For example, "K01647 citrate synthase" produces the following
        # "accession" and "function" strings:
        # "ko00001!!!ko00001!!!ko01000"
        # "KEGG Orthology (KO)>>>09100 Metabolism>>>09101 Carbohydrate metabolism>>>00020 Citrate cycle (TCA cycle)!!!
        #  KEGG Orthology (KO)>>>09100 Metabolism>>>09101 Carbohydrate metabolism>>>00630 Glyoxylate and dicarboxylate metabolism!!!
        #  Enzymes>>>2. Transferases>>>2.3 Acyltransferases>>>2.3.3 Acyl groups converted into alkyl groups on transfer>>>2.3.3.1 citrate (Si)-synthase"
        hierarchy_accession = ""
        categorizations = ""
        categorization_dicts = list(ortholog_brite_dict.values())
        for categorization_dict in categorization_dicts[: -1]:
            hierarchy_accession += f"{categorization_dict['hierarchy_accession']}!!!"
            categorizations += f"{categorization_dict['hierarchy_name']}>>>{categorization_dict['categorization']}!!!"
        last_categorization_dict = categorization_dicts[-1]
        hierarchy_accession += last_categorization_dict['hierarchy_accession']
        categorizations += f"{last_categorization_dict['hierarchy_name']}>>>{last_categorization_dict['categorization']}"

        ortholog_categorizations_dict = {'gene_callers_id': gene_callers_id,
                                         'source': 'KEGG_BRITE',
                                         'accession': hierarchy_accession,
                                         'function': categorizations,
                                         'e_value': None}

        return ortholog_categorizations_dict


    def store_annotations_in_db(self):
        """Takes the dictionary of function annotations (already parsed, if necessary) and puts them in the DB.

        Should be called after the function that parses the HMM hits and creates self.functions_dict :) which is
        parse_kofam_hits()
        """

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        if self.functions_dict:
            gene_function_calls_table.create(self.functions_dict)
            if self.kegg_module_names_dict:
                gene_function_calls_table.create(self.kegg_module_names_dict)
            if self.kegg_module_classes_dict:
                gene_function_calls_table.create(self.kegg_module_classes_dict)
            if self.kegg_brite_categorizations_dict:
                gene_function_calls_table.create(self.kegg_brite_categorizations_dict)
        else:
            self.run.warning("There are no KOfam hits to add to the database. Returning empty handed, "
                             "but still adding KOfam as a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})


    def process_kofam_hmms(self):
        """This is a driver function for running HMMs against the KOfam database and processing the hits into the provided contigs DB."""

        tmp_directory_path = filesnpaths.get_temp_directory_path()
        contigs_db = ContigsSuperclass(self.args) # initialize contigs db
        # we will need the gene caller ids later
        all_gcids_in_contigs_db = contigs_db.genes_in_contigs_dict.keys()

        # safety check for previous annotations so that people don't overwrite those if they don't want to
        self.check_hash_in_contigs_db()

        # get AA sequences as FASTA
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_files_dict['AA:GENE'],
                                                      simple_headers=True,
                                                      report_aa_sequences=True)

        # run hmmscan
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer('KOfam', 'AA', 'GENE', None, None, len(self.ko_dict), self.kofam_hmm_file_path, None, None)

        has_stray_hits = False
        stray_hits_file = None
        if self.include_stray_kos:
            ohmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
            stray_hits_file = ohmmer.run_hmmer('Stray KOs', 'AA', 'GENE', None, None, len(self.stray_ko_dict), self.stray_ko_hmm_file_path, None, None)
            has_stray_hits = True if stray_hits_file else False

        if not hmm_hits_file and not has_stray_hits:
            run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                             "now anvi'o will add KOfam as a functional source with no hits, clean the temporary directories "
                             "and gracefully quit.", nl_before=1, nl_after=1)
            if not anvio.DEBUG:
                shutil.rmtree(tmp_directory_path)
                hmmer.clean_tmp_dirs()
            else:
                self.run.warning("Because you ran this script with the --debug flag, anvi'o will not clean up the temporary "
                                 "directories located at %s and %s. Please be responsible for cleaning up this directory yourself "
                                 "after you are finished debugging :)" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
            gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})
            return

        # set up some attributes that we'll need later
        self.functions_dict = {}
        self.kegg_module_names_dict = {}
        self.kegg_module_classes_dict = {}
        self.kegg_brite_categorizations_dict = {}
        self.gcids_to_hits_dict = {}
        self.gcids_to_functions_dict = {}
        super_hits_dict = {} # will store the hits from each set of HMMs

        # parse hmmscan output
        self.run.warning('', header='HMM hit parsing for KOfams', lc='green')
        self.run.info("HMM output table", hmm_hits_file)
        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()
        next_key_in_functions_dict = self.parse_kofam_hits(search_results_dict)
        super_hits_dict["KOfam"] = search_results_dict
        self.run.info("Current number of annotations in functions dictionary", len(self.functions_dict))

        if has_stray_hits:
            self.run.warning('', header='HMM hit parsing for Stray KOs', lc='green')
            self.run.info("HMM output table", stray_hits_file)
            oparser = parser_modules['search']['hmmer_table_output'](stray_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
            stray_search_results = oparser.get_search_results()
            next_key_in_functions_dict = self.parse_kofam_hits(stray_search_results, hits_label = "Stray KO", next_key=next_key_in_functions_dict)
            super_hits_dict["Stray KO"] = stray_search_results
            self.run.info("Current number of annotations in functions dictionary", len(self.functions_dict))

        if not self.skip_bitscore_heuristic:
            self.update_dict_for_genes_with_missing_annotations(all_gcids_in_contigs_db, super_hits_dict, next_key=next_key_in_functions_dict)

        # add functions and KEGG modules info to database
        self.store_annotations_in_db()

        # If requested, store bit scores of each hit in file
        if self.log_bitscores:
            self.bitscore_log_file = os.path.splitext(os.path.basename(self.contigs_db_path))[0] + "_bitscores.txt"
            anvio.utils.store_dict_as_TAB_delimited_file(search_results_dict, self.bitscore_log_file, key_header='entry_id')
            self.run.info("Bit score information file: ", self.bitscore_log_file)

        # mark contigs db with hash of modules.db content for version tracking
        self.set_hash_in_contigs_db()

        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                        "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single("Cleaning up the temp directory (you can use `--debug` if you would "
                            "like to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()


class KeggEstimatorArgs():
    def __init__(self, args, format_args_for_single_estimator=False, run=run, progress=progress):
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
        self.module_specific_matrices = A('module_specific_matrices') or None
        self.no_comments = True if A('no_comments') else False
        self.external_genomes_file = A('external_genomes') or None
        self.internal_genomes_file = A('internal_genomes') or None
        self.metagenomes_file = A('metagenomes') or None
        self.kegg_data_dir = A('kegg_data_dir')
        self.modules_unique_id = None
        self.ko_unique_id = None
        self.genome_mode = False  ## controls some warnings output, will be set to True downstream if necessary

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


    def setup_output_for_appending(self):
        """Initializes and returns a dictionary of AppendableFile objects, one for each output mode"""

        output_dict = {}
        for mode in self.output_modes:
            output_path = self.output_file_prefix + "_" + self.available_modes[mode]["output_suffix"]
            if filesnpaths.is_file_exists(output_path, dont_raise=True):
                raise ConfigError("It seems like output files with your requested prefix already exist, for "
                                  f"example: {output_path}. Please delete the existing files or provide a "
                                  "different output prefix.")
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

        self.all_modules_in_db = {}
        self.all_kos_in_db = {}
        self.module_paths_dict = {}

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

                # if we have our own versions of any stray KOs, then we include them here to enable lookups downstream
                k_anvio = f"{k}{STRAY_KO_ANVIO_SUFFIX}"
                if self.include_stray_kos and k_anvio in self.ko_dict:
                    if k_anvio not in self.all_kos_in_db:
                        src = 'KOfam'
                        func = self.all_modules_in_db[mod]['ORTHOLOGY'][k] if 'ORTHOLOGY' in self.all_modules_in_db[mod] else self.ko_dict[k_anvio]['definition']
                        self.all_kos_in_db[k_anvio] = {'modules': [], 'annotation_source': src, 'function': func}
                    self.all_kos_in_db[k_anvio]['modules'].append(mod)


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
        this case results in an error.

        PARAMETERS
        ==========
        step_string: str
            A string containing the definition of one step from a module

        RETURNS
        =======
        step_string: str
            The same string, with nonessential enzyme accessions (if any) removed.
        """

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
        elif self.include_stray_kos and self.stray_ko_dict and knum in self.stray_ko_dict:
            # if we can't find the enzyme in the KO dictionary, try to find it in the stray KO dictionary (if it exists)
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


class KeggMetabolismEstimator(KeggContext, KeggEstimatorArgs):
    """ Class for reconstructing/estimating metabolism for a SINGLE contigs DB based on hits to KEGG databases.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-estimate-metabolism
    """

    def __init__(self, args, run=run, progress=progress):
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

        if self.enzymes_txt:
            self.contigs_db_project_name = os.path.basename(self.enzymes_txt).replace(".", "_")

        # INPUT OPTIONS SANITY CHECKS
        if not self.estimate_from_json and not self.contigs_db_path and not self.enzymes_txt and not self.pan_db_path:
            raise ConfigError("NO INPUT PROVIDED. Please use the `-h` flag to see possible input options.")
        # incompatible input options
        if (self.contigs_db_path and (self.pan_db_path or self.enzymes_txt)) or \
           (self.enzymes_txt and self.pan_db_path):
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
            if not self.enzymes_txt:
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
            run.info("Output Modes", ", ".join(self.output_modes))
            run.info("Module completeness threshold", self.module_completion_threshold)
            run.info("Only complete modules included in output", self.only_complete)
            run.info("Zero-completeness modules excluded from output", self.exclude_zero_modules)
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

        estimation_mode = "Genome (or metagenome assembly)"
        if self.profile_db_path and self.collection_name:
            if not self.metagenome_mode:
                estimation_mode = "Bins in a metagenome"
            else:
                estimation_mode = "Individual contigs within a collection in a metagenome"
        elif self.metagenome_mode:
            estimation_mode = "Individual contigs in a metagenome"
        elif self.enzymes_txt:
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
            contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
            self.contigs_db_project_name = contigs_db.meta['project_name']


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
                self.setup_stray_ko_dict(add_entries_to_regular_ko_dict=True)
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
        run.warning(None, header="AVAILABLE OUTPUT MODES", lc="green")

        for mode, mode_meta in self.available_modes.items():
            self.run.info(mode, mode_meta['description'])


    def list_output_headers(self):
        """This function prints out the available output headers for the 'custom' output mode"""
        run.warning(None, header="AVAILABLE OUTPUT HEADERS", lc="green")

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
                    min_contig_length_in_profile_db = pp(ProfileDatabase(self.profile_db_path).meta['min_contig_length'])
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
        if self.enzymes_txt: # for this input the sample name is just the name of the input file (dots converted to underscores)
            samples_list = [self.contigs_db_project_name]

        else:
            if not self.profile_db:
                self.args.skip_consider_gene_dbs = True
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
            if not self.only_user_modules and self.all_kos_in_db[knum]['annotation_source'] == 'KOfam' and knum not in self.ko_dict \
                        and (f"{knum}{STRAY_KO_ANVIO_SUFFIX}" not in self.ko_dict) and self.exclude_kos_no_threshold:
                mods_it_is_in = self.all_kos_in_db[knum]['modules']
                if mods_it_is_in:
                    if anvio.DEBUG:
                        mods_str = ", ".join(mods_it_is_in)
                        self.run.warning(f"Oh dear. We do not appear to have a KOfam profile for {knum}. This means "
                                        "that any modules this KO belongs to can never be fully complete (this includes "
                                        f"{mods_str}). ")
                    for m in mods_it_is_in:
                        if knum[0] != 'M':
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
                                     "contigs_to_genes" : {}
                                     }

        kos_not_in_modules = []
        for ko, gene_call_id, split, contig in kofam_hits_in_splits:
            if ko not in self.all_kos_in_db:
                kos_not_in_modules.append(ko)
                # KOs that are not in modules will not be initialized above in the ko hit dictionary, so we add them here if we haven't already
                if ko not in bin_level_ko_dict:
                    bin_level_ko_dict[ko] = {"gene_caller_ids" : set(),
                                             "modules" : None,
                                             "genes_to_contigs" : {},
                                             "contigs_to_genes" : {}
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

                # make sure we can count annotations to anvi'o versions of stray KO models by using KEGG's original accession
                is_anvio_version = False
                if ko.endswith(STRAY_KO_ANVIO_SUFFIX):
                    is_anvio_version = True
                    ko = ko.replace(STRAY_KO_ANVIO_SUFFIX, "")

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
                        bin_level_module_dict[m]["warnings"].add(f"used '{ko}{STRAY_KO_ANVIO_SUFFIX}' model to annotate {ko}")

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
                    # '--' no associated enzyme case, always False (assumed incomplete)
                    if step[cur_index+1] == "-":
                        step_is_present_condition_statement += "False"
                        cur_index += 2 # skip over both '-', the next character should be a space or end of DEFINITION line

                        if anvio.DEBUG:
                            self.run.warning(f"While estimating the stepwise completeness of KEGG module {mnum}, anvi'o saw "
                                             f"'--' in the module DEFINITION. This indicates a step in the pathway that has no "
                                             f"associated enzyme. By default, anvi'o is marking this step incomplete. But it may not be, "
                                             f"and as a result this module might be falsely considered incomplete. So it may be in your "
                                             f"interest to take a closer look at this individual module.")
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
                        # we assume that such steps are not complete
                        has_no_ko_step = True
                        warning_str = "'--' steps are assumed incomplete"
                        meta_dict_for_bin[mnum]["warnings"].add(warning_str)
                        atomic_step_copy_number.append(0)
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

        was_already_complete = meta_dict_for_bin[mnum]["stepwise_is_complete"]
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

        was_already_complete = meta_dict_for_bin[mod]["pathwise_is_complete"]
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

        if not self.enzymes_txt and not self.profile_db:
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
            if self.enzymes_txt:
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
                if self.enzymes_txt:
                    cov = self.enzymes_txt_data[self.enzymes_txt_data['gene_id'] == g]['coverage'].values[0]
                    det = self.enzymes_txt_data[self.enzymes_txt_data['gene_id'] == g]['detection'].values[0]
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
                post_str = step_string[close_parens_idx+1:]

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
                if step_string == '--': # no KO profile => no copy number
                    return 0
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
            raise ConfigError(f"You provided a JSON file generated from USER data, but you "
                              f"did not specify which data directory to use with the `--user-modules` flag.")
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


    def load_data_from_enzymes_txt(self):
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
            self.run.warning("Just so you know, your input enzymes-txt file contained some columns of data that we are not "
                             "going to use. This isn't an issue or anything, just an FYI. We're ignoring the following field(s): {e_str}")

        # check and warning for enzymes not in self.all_kos_in_db
        enzymes_not_in_modules = list(enzyme_df[~enzyme_df["enzyme_accession"].isin(self.all_kos_in_db.keys())]['enzyme_accession'].unique())
        if enzymes_not_in_modules:
            example = enzymes_not_in_modules[0]
            self.run.warning(f"FYI, some enzymes in the 'enzyme_accession' column of your input enzymes-txt file do not belong to any "
                             f"metabolic modules (that we know about). These enzymes will be ignored for the purposes of estimating module "
                             f"completeness, but should still appear in enzyme-related outputs (if those were requested). In case you are "
                             f"curious, here is one example (run this program with --debug to get a full list): {example}")

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


    def estimate_metabolism_from_enzymes_txt(self):
        """Estimates metabolism on a set of enzymes provided in a text file.

        This function assumes that all enzymes in the file are coming from a single genome, and is effectively the
        same as the estimate_for_genome() function.

        Requires the self.enzymes_txt_data attribute to have been established (ie, by loading the self.enzymes_txt file).
        We make fake splits and contigs to match the expected input to the atomic functions, and the contigs_db_project_name
        attribute has been set (previously) to the name of the enzyme txt file

        RETURNS
        =======
        enzyme_metabolism_superdict : dictionary of dictionary of dictionaries
            dictionary mapping the name of the enzyme txt file to its metabolism completeness dictionary
        enzyme_ko_superdict : dictionary of dictionary of dictionaries
            maps the name of the enzyme txt file to its KOfam hit dictionary
        """

        kofam_gene_split_contig = []
        # no splits or contigs here
        for gene_call_id, ko in zip(self.enzymes_txt_data["gene_id"], self.enzymes_txt_data["enzyme_accession"]):
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
                            return_subset_for_matrix_format=False, all_modules_in_db=None, all_kos_in_db=None, module_paths_dict=None):
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

        if self.estimate_from_json:
            kegg_metabolism_superdict = self.estimate_metabolism_from_json_data()
        else:
            # we either get the modules DB info from the previous class, or we have to initialize it here
            if all_modules_in_db:
                self.all_modules_in_db = all_modules_in_db
                self.all_kos_in_db = all_kos_in_db
                self.module_paths_dict = module_paths_dict
            else:
                self.init_data_from_modules_db()

            if self.enzymes_txt:
                self.enzymes_txt_data = self.load_data_from_enzymes_txt()
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_metabolism_from_enzymes_txt()
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

        if self.write_dict_to_json:
            self.store_metabolism_superdict_as_json(kegg_metabolism_superdict, self.json_output_file_path + ".json")

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
                    if self.enzymes_txt:
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

                    if self.add_coverage:
                        if self.enzymes_txt:
                            for s in self.coverage_sample_list:
                                sample_cov_header = s + "_coverage"
                                d[self.ko_unique_id][sample_cov_header] = self.enzymes_txt_data[self.enzymes_txt_data["gene_id"] == gc_id]["coverage"].values[0]
                                sample_det_header = s + "_detection"
                                d[self.ko_unique_id][sample_det_header] = self.enzymes_txt_data[self.enzymes_txt_data["gene_id"] == gc_id]["detection"].values[0]
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

    def __init__(self, args, run=run, progress=progress):
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
                if filesnpaths.is_file_exists(matrix_output_file, dont_raise=True):
                    raise ConfigError(f"Uh oh... there is already matrix output (such as {matrix_output_file}) "
                                      "using this file prefix in the current directory. Please either remove these files "
                                      "or give us a different prefix (with -O).")

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

        run.warning(None, header="AVAILABLE OUTPUT MODES", lc="green")

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
                                 "output matrices for them. Maybe you made a typo? Or put an extra comma in somewhere?")

            if mods_defined_by_mods:
                skipped_list = ", ".join(mods_defined_by_mods)
                self.run.warning(f"The following modules were completely defined by other modules and therefore we didn't "
                                 f"generate enzyme hit matrices for them: {skipped_list}.")


class ModulesDatabase(KeggContext):
    """To create or access a Modules DB.

    This DB should be created in the Kegg Data folder during KEGG setup, and will be populated with information from the
    Kegg Module files.

    If you want to load an existing database from the python terminal, all you need is the path to the database and an
    empty args object:
    ```
    >>> import argparse
    >>> from anvio import kegg
    >>> path_to_db = "YOUR/PATH/HERE/MODULES.db"
    >>> args = argparse.Namespace()
    >>> kegg.ModulesDatabase(path_to_db, args)
    ```
    """

    def __init__(self, db_path, args, module_data_directory=None, brite_data_directory=None, module_dictionary=None, pathway_dictionary=None, brite_dictionary=None, data_source='KEGG', skip_brite_hierarchies=False, run=run, progress=progress, quiet=False):
        self.db = None
        self.db_path = db_path
        self.module_data_directory = module_data_directory # only required for create()
        self.brite_data_directory = brite_data_directory # only required for create()
        self.module_dict = module_dictionary
        self.pathway_dict = pathway_dictionary
        self.brite_dict = brite_dictionary
        self.data_source = data_source
        # BRITE setup can be skipped to allow newer versions of the database to be consistent with older versions that lacked BRITE
        self.skip_brite_hierarchies = skip_brite_hierarchies
        self.run = run
        self.progress = progress
        self.quiet = quiet

        # keep track of functional annotation sources needed for the modules in this db
        self.annotation_sources = set()
        if self.data_source == 'KEGG':
            self.annotation_sources.add('KOfam')

        if anvio.DEBUG:
            self.run.info("Modules DB quiet param", self.quiet)

        # init the base class for access to shared paths and such
        KeggContext.__init__(self, args)

        # modules table info
        self.module_table_name = t.module_table_name
        self.module_table_structure = t.module_table_structure
        self.module_table_types = t.module_table_types

        # pathway maps table info
        self.pathway_table_name = t.pathway_table_name
        self.pathway_table_structure = t.pathway_table_structure
        self.pathway_table_types = t.pathway_table_types

        # BRITE hierarchies table info
        self.brite_table_name = t.brite_table_name
        self.brite_table_structure = t.brite_table_structure
        self.brite_table_types = t.brite_table_types

        if os.path.exists(self.db_path):
            utils.is_kegg_modules_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__kegg_modules_version__, new_database=False)

            if not self.quiet:
                self.run.info("Modules database", f"An existing database, {self.db_path}, has been loaded.", quiet=self.quiet)
                self.run.info("Modules", f"{self.db.get_meta_value('num_modules')} found", quiet=self.quiet)
                self.run.info("BRITE KO hierarchies", f"{self.db.get_meta_value('num_brite_hierarchies')} found", quiet=self.quiet)

        else:
            # if self.module_dict is None, then we tried to initialize the DB outside of setup
            if not self.module_dict:
                raise ConfigError("ERROR - a new ModulesDatabase() cannot be initialized without providing a modules dictionary. This "
                                  "usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-data` may fix this.")

            if not self.skip_brite_hierarchies and not self.brite_dict:
                raise ConfigError("ERROR - a new ModulesDatabase() cannot be initialized without providing a BRITE dictionary. This "
                                  "usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-data` may fix this.")

######### DB GENERATION FUNCTIONS #########

    def touch(self):
        """Creates an empty Modules database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db.
        """

        # sanity check to avoid overriding previous Modules DB
        # this will probably never happen as long as this function is called through the setup script, but we check just in case
        if os.path.exists(self.db_path):
            raise ConfigError("A modules database at %s already exists. Please use the --reset flag when you restart the setup "
                              "if you really want to get rid of this one and make a new one." % (self.db_path))


        self.db = db.DB(self.db_path, anvio.__kegg_modules_version__, new_database=True)

        self.db.create_table(self.module_table_name, self.module_table_structure, self.module_table_types)
        self.db.create_table(self.brite_table_name, self.brite_table_structure, self.brite_table_types)


    def data_vals_sanity_check(self, data_vals, current_data_name, current_module_num):
        """This function checks if the data values were correctly parsed from a line in a KEGG module file.

        This is a sadly necessary step because some KEGG module file lines are problematic and don't follow the right format (ie, 2+ spaces
        between different fields). So here we check if the values that we parsed look like they are the right format, without any extra bits.
        Each data name (ORTHOLOGY, DEFINITION, etc) has a different format to check for.

        Note that we don't check the following data name types: NAME, CLASS, REFERENCE

        WARNING: The error checking and correction is by no means perfect and may well fail when KEGG is next updated. :(

        PARAMETERS
        ==========
        data_vals : str
            the data values field (split from the kegg module line)
        current_data_name : str
            which data name we are working on. It should never be None because we should have already figured this out by parsing the line.
        current_module_num : str
            which module we are working on. We need this to keep track of which modules throw parsing errors.

        RETURNS
        =======
        is_ok : bool
            whether the values look correctly formatted or not
        """

        is_ok = True
        is_corrected = False
        corrected_vals = None
        corrected_def = None

        if not current_data_name:
            raise ConfigError("data_vals_sanity_check() cannot be performed when the current data name is None. Something was not right "
                              "when parsing the KEGG module line.")
        elif current_data_name == "ENTRY":
            # example format: M00175
            if data_vals[0] != 'M' or len(data_vals) != 6:
                is_ok = False
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)
        elif current_data_name == "DEFINITION":
            # example format: (K01647,K05942) (K01681,K01682) (K00031,K00030) (K00164+K00658+K00382,K00174+K00175-K00177-K00176)
            # another example: (M00161,M00163) M00165
            knums = [x for x in re.split('\(|\)|,| |\+|-',data_vals) if x]
            for k in knums:
                if k[0] not in ['K','M'] or len(k) != 6:
                    is_ok = False
            if not is_ok: # this goes here to avoid counting multiple errors for the same line
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)
        elif current_data_name == "ORTHOLOGY":
            # example format: K00234,K00235,K00236,K00237
            # more complex example: (K00163,K00161+K00162)+K00627+K00382-K13997
            # another example:  (M00161         [ie, from (M00161  Photosystem II)]
            knums = [x for x in re.split('\(|\)|,|\+|-', data_vals) if x]
            for k in knums:
                if k[0] not in ['K','M'] or len(k) != 6:
                    is_ok = False
            # try to fix it by splitting on first space
            if not is_ok:
                self.parsing_error_dict['bad_line_splitting'].append(current_module_num)
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                # double check that we don't have a knum in the new definition
                if re.match("K\d{5}",corrected_def):
                    corrected_vals = "".join([corrected_vals,corrected_def])
                    corrected_def = None
                is_corrected = True
        elif current_data_name == "PATHWAY":
            # example format: map00020
            if data_vals[0:3] != "map" or len(data_vals) != 8:
                is_ok = False
                self.parsing_error_dict['bad_line_splitting'].append(current_module_num)
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                is_corrected = True
        elif current_data_name == "REACTION":
            # example format: R01899+R00268,R00267,R00709
            rnums = [x for x in re.split(',|\+', data_vals) if x]
            for r in rnums:
                if r[0] != 'R' or len(r) != 6:
                    is_ok = False
            if not is_ok:
                self.parsing_error_dict['bad_line_splitting'].append(current_module_num)
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                is_corrected = True
        elif current_data_name == "COMPOUND":
            # example format: C00024
            if data_vals[0] not in ['C','G'] or len(data_vals) != 6:
                is_ok = False
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)
        elif current_data_name == "RMODULE":
            # example format: RM003
            if data_vals[0:2] != "RM" or len(data_vals) != 5:
                is_ok = False
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)


        if not is_ok and not is_corrected:
            self.num_uncorrected_errors += 1
            if self.just_do_it:
                self.progress.reset()
                self.run.warning("While parsing, anvi'o found an uncorrectable issue with a KEGG Module line in module %s, but "
                                 "since you used the --just-do-it flag, anvi'o will quietly ignore this issue and add the line "
                                 "to the MODULES.db anyway. Please be warned that this may break things downstream. In case you "
                                 "are interested, the line causing this issue has data name %s and data value %s."
                                 % (current_module_num, current_data_name, data_vals))
                is_ok = True # let's pretend that everything is alright so that the next function will take the original parsed values

            else:
                raise ConfigError("While parsing, anvi'o found an uncorrectable issue with a KEGG Module line in module %s. The "
                                  "current data name is %s, here is the incorrectly-formatted data value field: %s. If you think "
                                  "this is totally fine and want to ignore errors like this, please re-run the setup with the "
                                  "--just-do-it flag. But if you choose to do that of course we are obliged to inform you that things "
                                  "may eventually break as a result." % (current_module_num, current_data_name, data_vals))

        if is_corrected:
            self.num_corrected_errors += 1
            if anvio.DEBUG and not self.quiet:
                self.progress.reset()
                self.run.warning("While parsing a KEGG Module line, we found an issue with the formatting. We did our very best to parse "
                                 "the line correctly, but please check that it looks right to you by examining the following values.")
                self.run.info("Incorrectly parsed data value field", data_vals)
                self.run.info("Corrected data values", corrected_vals)
                self.run.info("Corrected data definition", corrected_def)

        return is_ok, corrected_vals, corrected_def


    def parse_kegg_modules_line(self, line, current_module, line_num=None, current_data_name=None, error_dictionary=None):
        """This function parses information from one line of a KEGG module file.

        These files have fields separated by 2 or more spaces. Fields can include data name (not always), data value (always), and data definition (not always).
        Lines for pathway module files can have between 1 and 4 fields, but in fact the only situation where there should be 4 lines is the ENTRY data,
        which for some inexplicable reason has multiple spaces between "Pathway" and "Module" in the data definition field. We can safely ignore this last "Module", I think.

        Some lines will have multiple entities in the data_value field (ie, multiple KOs or reaction numbers) and will be split into multiple db entries.

        PARAMETERS
        ==========
        line : str
            the line to parse
        current_module : str
            which module we are working on. We need this to keep track of which modules throw parsing errors
        line_num : int
            which line number we are working on. We need this to keep track of which entities come from the same line of the file.
        current_data_name : str
            which data name we are working on. If this is None, we need to parse this info from the first field in the line.

        RETURNS
        =======
        line_entries : list
            tuples, each containing information for one db entry, namely data name, data value, data definition, and line number.
            Not all parts of the db entry will be included (module num, for instance), so this information must be parsed and combined with
            the missing information before being added to the database.
        """

        if anvio.DEBUG:
            self.progress.reset()
            self.run.info("[DEBUG] Parsing line", line, mc='red', lc='yellow')
        fields = re.split('\s{2,}', line)
        data_vals = None
        data_def = None
        line_entries = []

        # when data name unknown, parse from first field
        if not current_data_name:
            # sanity check: if line starts with space then there is no data name field and we should have passed a current_data_name
            if line[0] == ' ':
                raise ConfigError("Oh, please. Some silly developer (you know who you are) has tried to call parse_kegg_modules_line() on "
                                  "a line without a data name field, and forgot to give it the current data name. Shame on you, go fix "
                                  "this. (For reference here is the line: %s)" % (line))

            current_data_name = fields[0]
        # note that if data name is known, first field still exists but is actually the empty string ''
        # so no matter the situation, data value is field 1 (index 0) and data definition (if any) is field 2 (index 1)
        # the only exception is that sometimes there is nothing in the data definition field (REFERENCE lines sometimes do this)
        if len(fields) > 1:
            data_vals = fields[1]

            if self.data_source == 'KEGG':
                # need to sanity check data value field because SOME modules don't follow the 2-space separation formatting
                vals_are_okay, corrected_vals, corrected_def = self.data_vals_sanity_check(data_vals, current_data_name, current_module)
            else:
                # TODO: USER data sanity check
                vals_are_okay = True # for now we always assume USER data is formatted properly

            if vals_are_okay and len(fields) > 2: # not all lines have a definition field
                data_def = fields[2]
            elif not vals_are_okay:
                data_vals = corrected_vals
                data_def = corrected_def
        else: # only the data name was in the line
            # these are the data types that we don't care if they have an empty line
            data_types_can_be_empty = ['REFERENCE', 'AUTHORS', 'TITLE', 'JOURNAL']
            if current_data_name in data_types_can_be_empty or self.just_do_it:
                if anvio.DEBUG:
                    self.run.warning(f"While parsing module {current_module} we found an empty {current_data_name} line. "
                                     "We think it is okay and it probably won't cause issues downstream.",
                                     header="DEBUG OUTPUT", lc='yellow')
            else:
                raise ConfigError(f"While parsing module {current_module} we found an empty {current_data_name} line. "
                                  "We are quitting here so you can check it, because this data type might be important. "
                                  "However, if you disagree, you can re-run the setup with --just-do-it and we will quietly "
                                  "incorporate this empty line into the MODULES.db (you may also need the --reset flag when you re-run). ")

        # some types of information may need to be split into multiple db entries
        data_types_to_split = ["ORTHOLOGY","REACTION"] # lines that fall under these categories need to have data_vals split on comma
        if current_data_name in data_types_to_split:
            # here we should NOT split on any commas within parentheses
            vals = [x for x in re.split('\(|\)|,|\+|-', data_vals) if x]
            for val in vals:
                line_entries.append((current_data_name, val, data_def, line_num))
        else:
            line_entries.append((current_data_name, data_vals, data_def, line_num))

        return line_entries


    def create(self):
        """Creates the Modules DB"""

        if not self.module_data_directory:
            raise ConfigError("Some dumb programmer forgot to provide a module_data_directory parameter value to the ModulesDatabase "
                              "class. The DB can't be created unless it knows where the modules are... Get yourself together.")
        if not self.skip_brite_hierarchies and not self.brite_data_directory:
            raise ConfigError("Some dumb programmer forgot to provide a brite_data_directory parameter value to the ModulesDatabase "
                              "class. The DB can't be created unless it knows where the BRITE hierarchies are... Get yourself together.")

        self.touch()

        self.progress.new("Loading %s modules into Modules DB..." % len(self.module_dict.keys()))

        # sanity check that we setup the modules previously.
        # It shouldn't be a problem since this function should only be called during the setup process after modules download, but just in case.
        if not os.path.exists(self.module_data_directory) or len(self.module_dict.keys()) == 0:
            raise ConfigError("Appparently, the module data files were not correctly setup and now all sorts of things are broken. The "
                              "Modules DB cannot be created from broken things. BTW, this error is not supposed to happen to anyone "
                              "except maybe developers, so if you do not fall into that category you are likely in deep doo-doo. "
                              "Maybe re-running setup with --reset will work? (if not, you probably should email/Discord/telepathically "
                              "cry out for help to the developers). Anyway, if this helps make things any clearer, the number of modules "
                              "in the module dictionary is currently %s" % len(self.module_dict.keys()))

        # init the Modules table
        mod_table = ModulesTable(self.module_table_name)

        # keep track of errors encountered while parsing
        self.parsing_error_dict = {"bad_line_splitting" : [], "bad_kegg_code_format" : []}
        self.num_corrected_errors = 0
        self.num_uncorrected_errors = 0
        # keep track of BRITE parsing information
        self.num_hierarchies_parsed = 0
        self.num_brite_categorizations = 0

        num_modules_parsed = 0
        line_number = 0
        for mnum in self.module_dict.keys():
            self.progress.update("Parsing Module %s" % mnum)
            mod_file_path = os.path.join(self.module_data_directory, mnum)
            f = open(mod_file_path, 'r')

            prev_data_name_field = None
            module_has_annotation_source = False
            orthology_to_annotation_source = {}  # for sanity check that each enzyme has an annotation source
            # for sanity check that each enzyme in definition has an orthology line
            mod_definition = []
            orth_list = set()
            for line in f.readlines():
                line = line.strip('\n')
                line_number += 1

                # check for last line ///. We don't want to send the last line to the parsing function because it will break.
                # we also check here that the line is not entirely blank (this happens sometimes in KEGG modules, inexplicably)
                if not line == '///' and re.search(r"\S+", line):
                    # parse the line into a tuple
                    entries_tuple_list = None
                    # here is the tricky bit about parsing these files. Not all lines start with the data_name field; those that don't start with a space.
                    # if this is the case, we need to tell the parsing function what the previous data_name field has been.
                    if line[0] == ' ':
                        entries_tuple_list = self.parse_kegg_modules_line(line, mnum, line_number, prev_data_name_field)
                    else:
                        entries_tuple_list = self.parse_kegg_modules_line(line, mnum, line_number)

                    prev_data_name_field = entries_tuple_list[0][0]

                    for name, val, definition, line in entries_tuple_list:
                        # there is one situation in which we want to ignore the entry, and that is Modules appearing in the ORTHOLOGY category, like so:
                        # (M00531  Assimilatory nitrate reduction, nitrate => ammonia)
                        if not (name == "ORTHOLOGY" and val[0] == '('):
                            # append_and_store will collect db entries and store every 10000 at a time
                            mod_table.append_and_store(self.db, mnum, name, val, definition, line)
                        else:
                            line -= 1

                        # save definition for later sanity check
                        if name == 'DEFINITION':
                            mod_definition.append(val)
                        if name == 'ORTHOLOGY':
                            orth_list.add(val)
                        # keep track of distinct annotation sources for user modules
                        if self.data_source != 'KEGG' and name == "ORTHOLOGY":
                            orthology_to_annotation_source[val] = None
                        if self.data_source != 'KEGG' and name == "ANNOTATION_SOURCE":
                            self.annotation_sources.add(definition)
                            module_has_annotation_source = True
                            if val not in orthology_to_annotation_source:
                                raise ConfigError(f"Woah. While parsing module {mnum} we found an ANNOTATION_SOURCE for "
                                                  f"an enzyme, {val}, that does not have an ORTHOLOGY line. Please check the "
                                                  f"module file and make sure that 1) ORTHOLOGY comes before ANNOTATION_SOURCE, "
                                                  f"and 2) each enzyme with an ORTHOLOGY also has a corresponding ANNOTATION_SOURCE.")
                            orthology_to_annotation_source[val] = definition

                        # sanity check for user-defined CLASS value
                        if self.data_source != 'KEGG' and name == "CLASS":
                            if len(val.split(";")) != 3:
                                raise ConfigError(f"The module {mnum} appears to have an invalid CLASS value. That value should be "
                                                  f"a string with a class, category, and sub-category separated by semi-colons (for a "
                                                  f"total of two semi-colons in the string). Instead, it is this: {val}")

                f.close()

            # every enzyme in the module definition needs an orthology line
            mod_definition = " ".join(mod_definition)
            # anything that is not (),-+ should be converted to spaces, then we can split on the spaces to get the accessions
            mod_definition = re.sub('[\(\)\+\-,]', ' ', mod_definition).strip()
            acc_list = re.split(r'\s+', mod_definition)
            accessions_in_def = set(acc_list)
            # remove any accession that is for a module (ie, not an enzyme)
            enzymes_in_def = set([acc for acc in accessions_in_def if acc not in self.module_dict.keys()])
            enzymes_without_orth = enzymes_in_def.difference(orth_list)
            if self.data_source != 'KEGG' and enzymes_without_orth:
                bad_list = ", ".join(enzymes_without_orth)
                n = len(enzymes_without_orth)
                raise ConfigError(f"So, there is a thing. And that thing is that there {P('is an enzyme', n, alt='are some enzymes')} "
                                  f"in the DEFINITION string for module {mnum} that {P('does', n, alt='do')} not have an ORTHOLOGY line, "
                                  f"and {P('it', n, alt='they')} really should have one. Here {P('it is', n, alt='they are')}: {bad_list}")

            # every user module needs at least one annotation source
            if self.data_source != 'KEGG' and not module_has_annotation_source:
                os.remove(self.db_path)
                raise ConfigError(f"While parsing user module {mnum}, we noticed that it does not have a single "
                                  f"'ANNOTATION_SOURCE' field. We are sorry to tell you that this is not okay, "
                                  f"because every user-defined module requires at least one of those fields to tell "
                                  f"anvi'o where to find its gene annotations. So you should go and take a look at "
                                  f"the module file at {mod_file_path} and add one 'ANNOTATION_SOURCE' line for each "
                                  f"gene in the module definition, before re-trying this setup program. Thank you!")
            # every enzyme needs an annotation source
            if self.data_source != 'KEGG' and not all(orthology_to_annotation_source.values()):
                nones = [k for k,v in orthology_to_annotation_source.items() if not v]
                nones_str = ", ".join(nones)
                n = len(nones)
                raise ConfigError(f"*Dalek noises* EXTERMINATE! EXTERMINATE! It seems that your module {mnum} contains "
                                  f"{P('an enzyme that is', n, alt='some enzymes that are')} missing an ANNOTATION_SOURCE line. "
                                  f"Please go back in time and fix this. Here {P('is the enzyme', n, alt='are the enzymes')} in "
                                  f"question here: {nones_str}")

            num_modules_parsed += 1
        # once we are done parsing all modules, we store whatever db entries remain in the db_entries list
        # this is necessary because append_and_store() above only stores every 10000 entries
        self.progress.update("Storing final batch of module entries into DB")
        mod_table.store(self.db)

        self.progress.end()

        # warn user about parsing errors
        if self.num_corrected_errors > 0 or self.num_uncorrected_errors > 0:
            if anvio.DEBUG:
                self.run.warning("Several parsing errors were encountered while building the Modules DB. "
                                 "Below you will see which modules threw each type of parsing error. Note that modules which "
                                 "threw multiple errors will occur in the list as many times as it threw each error.")
                self.run.info("Bad line splitting (usually due to rogue or missing spaces)", self.parsing_error_dict["bad_line_splitting"])
                self.run.info("Bad KEGG code format (not corrected; possibly problematic)", self.parsing_error_dict["bad_kegg_code_format"])
            else: # less verbose
                self.run.warning("First things first - don't panic. Several parsing errors were encountered while building the Modules DB. "
                                 "But that is probably okay, because if you got to this point it is likely that we already fixed all of them "
                                 "ourselves. So don't worry too much. Below you will see how many of each type of error was encountered. If "
                                 "you would like to see which modules threw these errors, please re-run the setup using the --debug flag (you "
                                 "will also probably need the --reset or --overwrite-output-destinations flag). When doing so, you will also "
                                 "see which lines caused issues; this can be a lot of output, so you can suppress the line-specific output with "
                                 "the `--quiet` flag if that makes things easier to read. So, in summary: You can probably ignore this warning. "
                                 "But if you want more info: run setup again with `--reset --debug --quiet` to see exactly which modules had "
                                 "issues, or run with `--reset --debug` to see exactly which lines in which modules had issues. Anvi'o developers "
                                 "thank you for your attention and patience 😇")
                self.run.info("Bad line splitting (usually due to rogue or missing spaces)", len(self.parsing_error_dict["bad_line_splitting"]))
                self.run.info("Bad KEGG code format (usually not correctable)", len(self.parsing_error_dict["bad_kegg_code_format"]))

        if not self.annotation_sources:
            os.remove(self.db_path)
            raise ConfigError("We're not sure how you made it this far without having any annotation sources defined in your module files, "
                              "because we should have noticed this while parsing them. But it happened, and here we are. You need to go add "
                              "'ANNOTATION_SOURCE' fields to those module files, and then re-do this setup.")
        annotation_source_list = ",".join(list(self.annotation_sources))

        self.populate_brite_table()

        # give some run info
        self.run.info('Modules database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of modules', num_modules_parsed, quiet=self.quiet)
        self.run.info('Number of module entries', mod_table.get_total_entries(), quiet=self.quiet)
        self.run.info('Number of module parsing errors (corrected)', self.num_corrected_errors, quiet=self.quiet)
        self.run.info('Number of module parsing errors (uncorrected)', self.num_uncorrected_errors, quiet=self.quiet)
        self.run.info('Annotation sources required for estimation', ", ".join(self.annotation_sources))
        if not self.skip_brite_hierarchies and self.brite_dict:
            self.run.info('Number of BRITE hierarchies', self.num_hierarchies_parsed, quiet=self.quiet)
            self.run.info('Number of ortholog BRITE categorizations', self.num_brite_categorizations, quiet=self.quiet)

        # record some useful metadata
        self.db.set_meta_value('db_type', 'modules')
        self.db.set_meta_value('data_source', self.data_source)
        self.db.set_meta_value('annotation_sources', annotation_source_list)
        self.db.set_meta_value('num_modules', num_modules_parsed)
        self.db.set_meta_value('total_module_entries', mod_table.get_total_entries())
        if not self.skip_brite_hierarchies and self.brite_dict:
            self.db.set_meta_value('is_brite_setup', True)
            self.db.set_meta_value('num_brite_hierarchies', self.num_hierarchies_parsed)
            self.db.set_meta_value('total_brite_entries', self.num_brite_categorizations)
        else:
            self.db.set_meta_value('is_brite_setup', False)
            self.db.set_meta_value('num_brite_hierarchies', None)
            self.db.set_meta_value('total_brite_entries', None)
        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value('hash', self.get_db_content_hash())
        self.db.set_meta_value('version', t.metabolic_modules_db_version)

        self.db.disconnect()


    def populate_brite_table(self):
        if self.skip_brite_hierarchies or not self.brite_dict:
            return

        self.progress.new("Loading BRITE hierarchies into Modules DB...")

        # init the BRITE table
        brite_table = BriteTable(self.brite_table_name)

        num_hierarchies_parsed = 0
        unrecognized_items = []
        for hierarchy in self.brite_dict:
            self.progress.update(f"Parsing BRITE hierarchy '{hierarchy}'")
            hierarchy_accession = hierarchy[: 7] # the validity of the hierarchy accession was checked in `KeggSetup.process_brite_hierarchy_of_hierarchies`
            hierarchy_name = hierarchy[7: ].lstrip()

            brite_file_path = os.path.join(self.brite_data_directory, hierarchy_accession)
            for ortholog, categorizations in self.invert_brite_json_dict(json.load(open(brite_file_path))).items():
                split_ortholog = ortholog.split(' ')
                ortholog_accession = split_ortholog[0]

                # record items in the hierarchy that do not have expected ortholog accessions formatted KXXXXX
                if len(ortholog_accession) != 6:
                    unrecognized_items.append(f'{hierarchy}: {ortholog}')
                    continue
                if ortholog_accession[0] != 'K':
                    unrecognized_items.append(f'{hierarchy}: {ortholog}')
                    continue
                try:
                    int(ortholog_accession[1: ])
                except ValueError:
                    unrecognized_items.append(f'{hierarchy}: {ortholog}')
                    continue

                ortholog_name = ' '.join(split_ortholog[1: ]).lstrip()

                # process each of the ortholog's categorizations in the hierarchy
                for categorization in categorizations:
                    if hierarchy_accession == 'ko00001':
                        # the expected top-level classes of the "ko00001  KEGG Orthology (KO)"
                        # hierarchy are "09100 Metabolism", "09120 Genetic Information Processing",
                        # "09130 Environmental Information Processing", "09140 Cellular Processes",
                        # "09150 Organismal Systems", "09160 Human Diseases", "09180 Brite
                        # Hierarchies", and "09190 Not Included in Pathway or Brite". "09180 Brite
                        # Hierarchies" is a representation of other hierarchies, such as "01001
                        # Protein kinases [BR:ko01001]", as totally flat categories, with all
                        # subcategories flattened out. We download and process json files for these
                        # other hierarchies separately. Therefore ignore entries in "09180 Brite
                        # Hierarchies". "09150 Organismal Systems" and "09160 Human Diseases" are
                        # also ignored due to their focus on human genes. The value of "09190 Not
                        # Included in Pathway or Brite" is debatable, but certain proteins are only
                        # found in this category, such as bacterial circadian clock proteins
                        # (classified under "09193 Unclassified: signaling and cellular processes"
                        # >>> "99995 Signaling proteins"), so it is retained.
                        category_accession = categorization[1].split(' ')[0]
                        if category_accession == '09180' or category_accession == '09150' or category_accession == '09160':
                            continue
                    brite_table.append_and_store(self.db, hierarchy_accession, hierarchy_name, ortholog_accession, ortholog_name, '>>>'.join(categorization[1: ])) # ignore the first category, the accession of the hierarchy itself
            num_hierarchies_parsed += 1
        self.num_hierarchies_parsed = num_hierarchies_parsed
        self.num_brite_categorizations = brite_table.get_total_entries()

        if unrecognized_items and anvio.DEBUG:
            self.run.warning("We attempted to parse some names of items in hierarchies as orthologs, "
                             "but ignored them since they did not start with an accession formatted 'KXXXXX', where 'XXXXX' are five digits. "
                             f"The following entries are formatted as '<hierarchy>: <ignored item>': {', '.join(set(unrecognized_items))}")

        # once we are done parsing all hierarchies, we store whatever db entries remain in the db_entries list
        # this is necessary because append_and_store() above only stores every 10000 entries
        self.progress.update("Storing final batch of BRITE entries into DB")
        brite_table.store(self.db)

        self.progress.end()


    def disconnect(self):
        self.db.disconnect()

######### SELF TABLE ACCESS FUNCTIONS #########

    def get_days_since_creation(self):
        """Returns the time (in days) since MODULES.db was created.

        The time units are seconds, and there are 60*60*24 = 86400 seconds per day,
        so we do the appropriate division to get the time in days.
        """

        return round((time.time() - float(self.db.get_meta_value('creation_date'))) / 86400)


    def get_db_content_hash(self):
        """Compute hash of all KOs and module numbers present in the db (used for tracking major changes to db content with future KEGG updates)"""
        mods = self.get_all_modules_as_list()
        mods.sort()
        orths = self.get_all_knums_as_list()
        orths.sort()
        mods_and_orths = mods + orths
        mods_and_orths = "".join(mods_and_orths)
        return str(hashlib.sha224(mods_and_orths.encode('utf-8')).hexdigest())[0:12]

######### MODULES TABLE ACCESS FUNCTIONS #########

    def get_modules_table_as_dict(self, data_names_of_interest=[]):
        """This function loads the modules table and returns it as a dictionary keyed by module.

        Every data_name for the module (NAME, DEFINITION, ORTHOLOGY, etc) will become a key in the
        inner dictionary, and the corresponding data_value will become its value. For data_names that
        have multiple values (ORTHOLOGY, COMPOUND, REACTION, etc), the value becomes yet another dictionary
        of data_value -> data_definition key-value pairs.

        The one exception is DEFINITION lines - there are occasionally multiple of these for one module but
        these must be returned as a list, not as a dictionary.

        PARAMETERS
        ==========
        data_names_of_interest : list of str
            the returned dictionary will contain only data names from this list. If the list is empty,
            all data names are returned.

        RETURNS
        =======
        module_dictionary : dict of dicts
            data for each module in the modules table. Outer dictionary is keyed by module number and
            inner dictionary is keyed by data name
        """

        if data_names_of_interest:
            data_names_list = [f"'{n}'" for n in data_names_of_interest]
            where_clause_string = f"data_name in ({','.join(data_names_list)})"
            # this WILL fail if you ask for a data name that doesn't exist, so know your data before you query
            dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True)
        else:
            dict_from_mod_table = self.db.get_table_as_dict(self.module_table_name, row_num_as_key=True)
        # the returned dictionary is keyed by an arbitrary integer, and each value is a dict containing one row from the modules table
        # ex of one row in this dict: 0: {'module': 'M00001', 'data_name': 'ENTRY', 'data_value': 'M00001', 'data_definition': 'Pathway', 'line': 1}

        # now we convert this to a per-module dictionary
        module_dictionary = {}
        for entry in dict_from_mod_table:
            mod = dict_from_mod_table[entry]['module']
            data_name = dict_from_mod_table[entry]['data_name']
            data_value = dict_from_mod_table[entry]['data_value']
            data_definition = dict_from_mod_table[entry]['data_definition']

            if mod not in module_dictionary:
                module_dictionary[mod] = {}

            if data_name not in module_dictionary[mod]:
                if not data_definition:
                    module_dictionary[mod][data_name] = data_value
                else:
                    module_dictionary[mod][data_name] = {data_value : data_definition}
            else:
                # place multiple definition lines into list
                if data_name == "DEFINITION":
                    if isinstance(module_dictionary[mod][data_name], list):
                        module_dictionary[mod][data_name].append(data_value)
                    else:
                        existing_val = module_dictionary[mod][data_name]
                        module_dictionary[mod][data_name] = [existing_val, data_value]
                else:
                    # data_value -> data_definition dictionary
                    if isinstance(module_dictionary[mod][data_name], dict):
                        if data_value in module_dictionary[mod][data_name]:
                            module_dictionary[mod][data_name][data_value] += " / " + data_definition
                        else:
                            module_dictionary[mod][data_name][data_value] = data_definition
                    else:
                        existing_val = module_dictionary[mod][data_name]
                        existing_val_data_def = dict_from_mod_table[entry-1]['data_definition']
                        module_dictionary[mod][data_name] = {existing_val: existing_val_data_def, data_value: data_definition}

        return module_dictionary


    def get_data_value_entries_for_module_by_data_name(self, module_num, data_name):
        """This function returns data_value elements from the modules table for the specified module and data_name pair.

        All elements corresponding to the pair (ie, M00001 and ORTHOLOGY) will be returned.
        The function relies on the db.get_some_rows_from_table_as_dict() function to first fetch all rows corresponding \
        to a particular model, and then parses the resulting dictionary to find all the elements with the given data_name field.

        PARAMETERS
        ==========
        module_num : str
            the module to fetch data for
        data_name : str
            which data_name field we want

        RETURNS
        =======
        data_values_to_ret : list of str
            the data_values corresponding to the module/data_name pair
        """

        where_clause_string = "module = '%s'" % (module_num)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True)
        # the returned dictionary is keyed by an arbitrary integer, and each value is a dict containing one row from the modules table
        # ex of one row in this dict: 0: {'module': 'M00001', 'data_name': 'ENTRY', 'data_value': 'M00001', 'data_definition': 'Pathway', 'line': 1}
        data_values_to_ret = []
        for key in dict_from_mod_table.keys():
            if dict_from_mod_table[key]['data_name'] == data_name:
                data_values_to_ret.append(dict_from_mod_table[key]['data_value'])

        if not data_values_to_ret:
            self.run.warning("Just so you know, we tried to fetch data from the Modules database for the data_name field %s "
                             "and module %s, but didn't come up with anything, so an empty list is being returned. This may "
                             "cause errors down the line, and if so we're very sorry for that.")

        return data_values_to_ret


    def get_data_definition_entries_for_module_by_data_name(self, module_num, data_name):
        """This function returns data_definition elements from the modules table for the specified module and data_name pair.

        All elements corresponding to the pair (ie, M00001 and ORTHOLOGY) will be returned.
        The function relies on the db.get_some_rows_from_table_as_dict() function to first fetch all rows corresponding \
        to a particular model, and then parses the resulting dictionary to find all the elements with the given data_name field.

        PARAMETERS
        ==========
        module_num : str
            the module to fetch data for
        data_name : str
            which data_name field we want

        RETURNS
        =======
        data_defs_to_ret : list of str
            the data_definitions corresponding to the module/data_name pair
        """

        where_clause_string = "module = '%s'" % (module_num)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True)

        data_defs_to_ret = []
        for key in dict_from_mod_table.keys():
            if dict_from_mod_table[key]['data_name'] == data_name:
                data_defs_to_ret.append(dict_from_mod_table[key]['data_definition'])

        if not data_defs_to_ret and anvio.DEBUG:
            self.run.warning("Just so you know, we tried to fetch data definitions from the Modules database for the data_name field %s "
                             "and module %s, but didn't come up with anything, so an empty list is being returned. This may "
                             "cause errors down the line, and if so we're very sorry for that.")

        return data_defs_to_ret


    def get_all_modules_as_list(self):
        """This function returns a list of all modules in the DB."""
        return self.db.get_single_column_from_table(self.module_table_name, 'module', unique=True)


    def get_all_knums_as_list(self):
        """This function returns a list of all KO numbers in the DB."""
        where_clause_string = "data_name = 'ORTHOLOGY'"
        return self.db.get_single_column_from_table(self.module_table_name, 'data_value', unique=True, where_clause=where_clause_string)


    def get_ko_function_dict(self):
        """This function returns a 2-level dictionary keyed by KO number. The inner dict contains a 'definition' field
           that stores the KO's functional annotation.

           This method effectively builds a partial KO dict similar to the one produced by the setup_ko_dict() function,
           but containing only the KOs/gene annotations from the modules DB. When being used for user-defined metabolism,
           it should be expanded later with gene annotations from the contigs database(s) being worked on.
        """
        where_clause_string = "data_name = 'ORTHOLOGY'"
        kos_and_functions = self.db.get_some_columns_from_table(self.module_table_name, "data_value,data_definition", unique=True, where_clause=where_clause_string)
        ko_func_dict = {}
        for k,f in kos_and_functions:
            if k not in ko_func_dict:
                ko_func_dict[k] = {'definition': f }
        return ko_func_dict


    def get_modules_for_knum(self, knum):
        """This function returns a list of modules that the given KO belongs to."""
        where_clause_string = "data_value = '%s'" % (knum)
        return self.db.get_single_column_from_table(self.module_table_name, 'module', unique=True, where_clause=where_clause_string)


    def get_module_classes_for_knum_as_dict(self, knum):
        """This function returns the classes for the modules that a given KO belongs to in a dictionary of dictionaries keyed by module number."""
        mods = self.get_modules_for_knum(knum)
        all_mods_classes_dict = {}
        for mnum in mods:
            all_mods_classes_dict[mnum] = self.get_kegg_module_class_dict(mnum)
        return all_mods_classes_dict


    def get_module_classes_for_knum_as_list(self, knum):
        """This function returns the classes for the modules that a given KO belongs to as a list of strings."""
        mods = self.get_modules_for_knum(knum)
        all_mods_classes_list = []
        for mnum in mods:
            mod_class = self.get_data_value_entries_for_module_by_data_name(mnum, "CLASS")[0]
            all_mods_classes_list.append(mod_class)
        return all_mods_classes_list


    def get_ortholog_brite_categorizations(self, ortholog_accession):
        """Return a list of the BRITE categorizations of the ortholog."""
        where_clause_string = f"ortholog_accession = '{ortholog_accession}'"
        return self.db.get_some_rows_from_table_as_dict(self.brite_table_name, where_clause=where_clause_string, error_if_no_data=False, row_num_as_key=True)


    def get_module_name(self, mnum):
        """This function returns the name of the specified KEGG module."""

        # there should only be one NAME per module, so we return the first list element
        return self.get_data_value_entries_for_module_by_data_name(mnum, "NAME")[0]


    def get_module_names_for_knum(self, knum):
        """This function returns all names of each KEGG module that the given KO belongs to in a dictionary keyed by module number."""
        mods = self.get_modules_for_knum(knum)
        module_names = {}
        for mnum in mods:
            module_names[mnum] = self.get_module_name(mnum)
        return module_names


    def parse_kegg_class_value(self, class_data_val):
        """This function takes a data_value string for the CLASS field in the modules table and parses it into a dictionary.

        The data_value string of CLASS fields should look something like this: Pathway modules; Amino acid metabolism; Lysine metabolism
        so they can be parsed into 3 parts: class, category, and subcategory.
        """

        fields = class_data_val.split("; ")
        class_dict = {"class" : fields[0], "category" : fields[1], "subcategory" : fields[2] if len(fields) > 2 else None}
        return class_dict


    def get_kegg_module_class_dict(self, mnum, class_value=None):
        """This function returns a dictionary of values in the CLASS field for a specific module

        It really exists only for convenience to put together the data fetch and parsing functions.

        PARAMETERS
        ==========
        mnum : str
            the module number
        class_value : str
            The 'CLASS' string for the module. This parameter is optional, and if it is not provided,
            the 'CLASS' value will be queried from the modules DB.
        """

        if not class_value:
            # there should only be one CLASS line per module, so we extract the first list element
            class_value = self.get_data_value_entries_for_module_by_data_name(mnum, "CLASS")[0]
        return self.parse_kegg_class_value(class_value)


    def get_kegg_module_definition(self, mnum):
        """This function returns module DEFINITION fields as one string"""

        def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        def_lines = [l.strip() for l in def_lines] # sometimes there are stray spaces at the end of the string that will mess us up later
        return " ".join(def_lines)


    def get_ko_definition_from_modules_table(self, ko_num):
        """This function returns the definition for the given KO from the modules data table.

        Note that the modules table will only contain information for KOs that belong to modules, so this
        function returns None for those KOs that are not in modules. If your use case depends on accessing
        definitions for all KOs, you are better off calling KeggContext.setup_ko_dict() and taking the
        definition from that dictionary.
        """

        where_clause_string = "data_name = 'ORTHOLOGY' AND data_value = '%s'" % (ko_num)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True, error_if_no_data=False)
        if not dict_from_mod_table:
            self.run.warning("get_ko_definition() speaking: No ORTHOLOGY entry found for KO %s - returning None."
                            % (ko_num))
            return None
        else:
            # there could be several rows for the same KO in different modules, but each definition should be
            # the same or similar, so we arbitrarily return the first one
            return dict_from_mod_table[0]['data_definition']

    
    def get_ko_reactions_from_modules_table(self, ko_num):
        """This function returns the KEGG reaction ID for the given KO from its data definition entry in the modules table.

        Reactions are indicated within brackets of the data definition entry, like these: [RN:R05339] or [RN:R01538 R03033]. 
        This function parses all reactions out of the entry and returns a list in which each reaction ID number is prefixed by 
        the standard KEGG reaction indicator 'RN:', as in ["RN:R01538", "RN:R03033"].
        
        Note that the modules table will only contain information for KOs that belong to modules, so this function 
        returns None for those KOs that are not in modules.
        """

        definition_line = self.get_ko_definition_from_modules_table(ko_num)
        if not definition_line: # this KO was not in the modules db
            return None
        
        def_fields = definition_line.split('[') # the last split should start with RN: and end with ]
        for f in def_fields:
            if f.startswith("RN:"):
                react_ids = f[3:-1].split(' ') # extract the ID numbers (without the initial RN: part or the final ])
                return ["RN:" + id for id in react_ids]



    def get_kos_in_module(self, mnum):
        """This function returns a list of KOs in the given module.

        It does this by parsing the ORTHOLOGY lines in the modules database. However,
        please note that these KOs are not always in the same order as the module
        definition, and may even contain duplicate entries for a KO. A good example
        of this is http://rest.kegg.jp/get/M00091 (K00551 is in two ORTHOLOGY lines)
        and http://rest.kegg.jp/get/M00176 (see KOs in the first top-level step). If
        this will be a problem, you should use the function get_kos_from_module_definition()
        instead.
        """

        return self.get_data_value_entries_for_module_by_data_name(mnum, "ORTHOLOGY")


    def get_kos_from_module_definition(self, mnum):
        """This function returns a list of KOs in the given module, in order of the DEFINITION.

        An alternative to get_kos_in_module().
        """

        mod_def = self.get_kegg_module_definition(mnum)
        ko_list = []
        k_indices = [x for x, v in enumerate(mod_def) if v == 'K']
        for idx in k_indices:
            ko_list.append(mod_def[idx:idx+6])

        return ko_list


    def get_kegg_module_compound_lists(self, mnum):
        """This function returns a list of substrates, a list of intermediates, and a list of products for the given module.

        We define 'substrate' to be any compound that is an input to but not an output from reactions in the module pathway.
        Likewise, a 'product' is any compound that is an output from but not an input to reactions in the module pathway.
        'Intermediate' is a compound that is both an input to and and output from reactions in the pathway.

        Note that this function refers to compounds by their KEGG identifier (format is 'C#####' where # is a digit).
        A separate function is used to convert these lists to human-readable compound names.

        RETURNS
        =======
        substrates : list
            Compounds that are only inputs to the module's metabolic pathway
        intermediates : list
            Compounds that are both outputs and inputs in the module's metabolic reactions
        products : list
            Compunds that are only outputs from the module's metabolic pathway
        """

        reactions_list = self.get_data_definition_entries_for_module_by_data_name(mnum, "REACTION")
        if not reactions_list:
            if anvio.DEBUG:
                self.run.warning(f"No REACTION entries found for module {mnum}, so no compounds will be returned by "
                                 "get_kegg_module_compound_lists()")

        inputs = set([])
        outputs = set([])

        for rxn_string in reactions_list:
            if '<->' in rxn_string:
                split_rxn = rxn_string.split('<->')
            else:
                split_rxn = rxn_string.split('->')
            if len(split_rxn) != 2:
                raise ConfigError(f"get_kegg_module_compound_lists('{mnum}') ran into an issue splitting the reaction {rxn_string}"
                                  "into 2 parts. Here is what the split looks like: {split_rxn}")
            rxn_inputs = [x.strip() for x in split_rxn[0].split('+')]
            rxn_outputs = [x.strip() for x in split_rxn[1].split('+')]
            # as of December 2023, some REACTION lines now include stoichiometry, like this:
            # C00390 + 2 C00125 -> C00399 + 2 C00126 + 2 C00080
            # to make sure we don't keep the stoichiometric numbers, we need to split each compound substring one more time
            # on a space and keep only the last element
            rxn_inputs = [x.split()[-1] for x in rxn_inputs]
            rxn_outputs = [x.split()[-1] for x in rxn_outputs]
            inputs = inputs.union(set(rxn_inputs))
            outputs = outputs.union(set(rxn_outputs))

        substrates = inputs.difference(outputs)
        products = outputs.difference(inputs)
        intermediates = inputs.intersection(outputs)

        return list(substrates), list(intermediates), list(products)


    def get_compound_dict_for_module(self, mnum, raise_error_if_no_data=False):
        """This function returns a dictionary mapping compound identifiers to their human-readable name for the given module

        If the module has no compounds, this function will either raise an error or return an empty dictionary depending on raise_error_if_no_data.
        If a compound doesn't have a human-readable name, then the compound identifier is used as the 'name'

        PARAMETERS
        ==========
        mnum : str
            module number to get compounds for
        raise_error_if_no_data : bool
            whether to quit all things if we don't get what we want
        """

        where_clause_string = "data_name = 'COMPOUND' AND module = '%s'" % (mnum)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True, error_if_no_data=raise_error_if_no_data)
        compound_dict = {}
        for key,row in dict_from_mod_table.items():
            compound = row['data_value']
            compound_name = row['data_definition']
            # if compound has no human-readable name in the database, we use the compound ID after all
            if not compound_name:
                compound_name = compound
            compound_dict[compound] = compound_name

        return compound_dict


    def get_human_readable_compound_lists_for_module(self, mnum):
        """This function returns a human-readable list of substrates, a list of intermediates, and a list of products for the given module.

        We define 'substrate' to be any compound that is an input to but not an output from reactions in the module pathway.
        Likewise, a 'product' is any compound that is an output from but not an input to reactions in the module pathway.
        'Intermediate' is a compound that is both an input to and and output from reactions in the pathway.

        RETURNS
        =======
        substrate_name_list : list of str
            List of substrate compounds
        intermediate_name_list : list of str
            List of intermediate compounds
        product_name_list : list of str
            List of product compounds
        """
        compound_to_name_dict = self.get_compound_dict_for_module(mnum)
        substrate_compounds, intermediate_compounds, product_compounds = self.get_kegg_module_compound_lists(mnum)

        substrate_name_list = [compound_to_name_dict[c] for c in substrate_compounds]
        intermediate_name_list = [compound_to_name_dict[c] for c in intermediate_compounds]
        product_name_list = [compound_to_name_dict[c] for c in product_compounds]

        return substrate_name_list, intermediate_name_list, product_name_list

######### BRITE TABLE ACCESS FUNCTIONS #########

    def get_brite_table_as_ortholog_dict(self, ortholog_accessions_of_interest=None, hierarchy_accessions_of_interest=None, category_substrings_of_interest=None, case_insensitive_substrings=False, use_ortholog_accessions_as_keys=False):
        """Load the BRITE hierarchies table as a dictionary keyed by ortholog.

        The returned dictionary contains each hierarchy and each categorization within the hierarchy
        in which the ortholog is found.

        The returned dictionary is structured as follows:
            {
                (<ortholog 1 accession>, <ortholog 1 name>):
                    {
                        (<hierarchy A accession>, <hierarchy A name>):
                            [
                                [(<category i accession>, <category i name>), (<category j accession>, <category j name>), ...],
                                [(<category k accession>, <category k name>), (<category l accession>, <category l name>), ...],
                                ...
                            ],
                        (<hierarchy B accession>, <hierarchy B name>):
                            [
                                [(<category x accession>, <category x name>), (<category y accession>, <category y name>), ...],
                                ...
                            ],
                        ...
                    },
                (<ortholog 2 accession>, <ortholog 2 name>):
                    {...},
                ...
            }

        Here is an example of the entry for arginyl-tRNA synthetase:
        ('K01887', 'RARS, argS; arginyl-tRNA synthetase [EC:6.1.1.19]'):
            {
                ('ko00001', 'KEGG Orthology (KO)'):
                    [
                        [('09120', 'Genetic Information Processing'), ('09122', 'Translation'), ('00970', 'Aminoacyl-tRNA biosynthesis')]
                    ],
                ('ko01000', 'Enzymes'):
                    [
                        [('6.', 'Ligases'), ('6.1', 'Forming carbon-oxygen bonds'), ('6.1.1', 'Ligases forming aminoacyl-tRNA and related compounds'), ('6.1.1.19', 'arginine---tRNA ligase')]
                    ],
                ('ko01007', 'Amino acid related enzymes'):
                    [
                        [('', 'Aminoacyl-tRNA synthetase'), ('', 'Class I (G)')]
                    ],
                ('ko03016', 'Transfer RNA biogenesis'):
                    [
                        [('', 'Eukaryotic type'), ('', 'Aminoacyl-tRNA synthetases (AARSs)'), ('', 'Multi-aminoacyl-tRNA synthetase complex (MSC)')],
                        [('', 'Prokaryotic type'), ('', 'Aminoacyl-tRNA synthetases (AARSs)'), ('', 'Other AARSs')]
                    ],
                ('ko03029', 'Mitochondrial biogenesis'):
                    [
                        [('', 'Mitochondrial DNA transcription, translation, and replication factors'), ('', 'Mitochondrial transcription and translation factors'), ('', 'Other mitochondrial DNA transcription and translation factors')]
                    ]
            }

        Keys and list items are split by accession and description, even in the absence of an
        accession for a category in the hierarchy. Given the hierarchies that are used in
        construction of the Modules database, only two are known to contain category "accessions."
        "ko01000 Enzyme" hierarchy categories yield EC number accessions, and "k00001 KEGG Orthology
        (KO)" hierarchy categories yield five digit accessions.

        Categorization lists proceed from most general to most specific level.

        Filtration with `hierarchy_accessions_of_interest` and `category_substrings_of_interest`
        returns the orthologs in the hierarchies and matched categories of interest, and also
        reduces the returned dictionary to the selected hierarchies and categorizations with matched
        categories. In the example of arginyl-tRNA synthetase, if `hierarchy_accessions_of_interest`
        is ['ko03016'] and `category_substrings_of_interest` is ['aminoacyl-tRNA synthetase'], then
        the returned dictionary becomes:
        ('K01887', 'RARS, argS; arginyl-tRNA synthetase [EC:6.1.1.19]'):
            {
                ('ko03016', 'Transfer RNA biogenesis'):
                    [
                        [('', 'Eukaryotic type'), ('', 'Aminoacyl-tRNA synthetases (AARSs)'), ('', 'Multi-aminoacyl-tRNA synthetase complex (MSC)')],
                        [('', 'Prokaryotic type'), ('', 'Aminoacyl-tRNA synthetases (AARSs)'), ('', 'Other AARSs')]
                    ]
            }

        PARAMETERS
        ==========
        ortholog_accessions_of_interest : list, None
            filters results to orthologs of interest

        hierarchy_accessions_of_interest : list, None
            filters results to hierarchies of interest

        category_substrings_of_interest : list, None
            filters results to categories containing substrings of interest

        case_insensitive_substrings : bool, False
            changes category substring search to be case insensitive

        use_ortholog_accessions_as_keys : bool, False
            ortholog keys of returned dictionary are accession strings rather than tuples

        RETURNS
        =======
        ortholog_dict : dict
            dictionary of ortholog BRITE categorizations
        """

        if ortholog_accessions_of_interest or hierarchy_accessions_of_interest:
            # filter table by orthologs, hierarchies, and category substrings of interest
            where_clause_string = ""

            if ortholog_accessions_of_interest:
                ortholog_list = [f"'{knum}'" for knum in ortholog_accessions_of_interest]
                where_clause_string += f"ortholog_accession IN ({','.join(ortholog_list)})"

            if hierarchy_accessions_of_interest:
                hierarchy_list = [f"'{konum}'" for konum in hierarchy_accessions_of_interest]
                if where_clause_string:
                    where_clause_string += " AND "
                where_clause_string += f"hierarchy_accession IN ({','.join(hierarchy_list)})"

            if category_substrings_of_interest:
                if where_clause_string:
                    where_clause_string += " AND ("
                for substring in category_substrings_of_interest:
                    if case_insensitive_substrings:
                        where_clause_string += f"UPPER(categorization) LIKE UPPER('%{substring}%') OR "
                    else:
                        where_clause_string += f"categorization LIKE '%{substring}%' OR "
                where_clause_string = where_clause_string[: -4] + ")"

            # this WILL fail if you ask for a data name that doesn't exist, so know your data before you query
            dict_from_brite_table = self.db.get_some_rows_from_table_as_dict(self.brite_table_name, where_clause_string, row_num_as_key=True)
        else:
            dict_from_brite_table = self.db.get_table_as_dict(self.brite_table_name, row_num_as_key=True)

        # the returned dict is keyed by an arbitrary integer, and each value is a dict containing one row from the BRITE table, e.g.,
        # 0: {'hierarchy_accession': 'ko00001',
        #     'hierarchy_name': 'KEGG Orthology (KO)',
        #     'ortholog_accession': 'K00001',
        #     'ortholog_name': 'E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]',
        #     'categorization': '09100 Metabolism>>>09101 Carbohydrate metabolism>>>00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'}

        # now we convert this to a per-ortholog dict
        ortholog_dict = {}
        for entry_dict in dict_from_brite_table.values():
            ortholog_accession = entry_dict['ortholog_accession']
            ortholog_name = entry_dict['ortholog_name']
            hierarchy_accession = entry_dict['hierarchy_accession']
            hierarchy_name = entry_dict['hierarchy_name']
            categorization = entry_dict['categorization']

            if use_ortholog_accessions_as_keys:
                ortholog_key = ortholog_accession
            else:
                ortholog_key = (ortholog_accession, ortholog_name)
            try:
                hierarchy_dict = ortholog_dict[ortholog_key]
            except KeyError:
                ortholog_dict[ortholog_key] = hierarchy_dict = {}

            hierarchy_key = (hierarchy_accession, hierarchy_name)
            try:
                category_list = hierarchy_dict[hierarchy_key]
            except KeyError:
                hierarchy_dict[hierarchy_key] = category_list = []

            categories = categorization.split('>>>')
            if hierarchy_accession == 'ko00001' or hierarchy_accession == 'ko01000':
                # the hierarchies, "ko00001 KEGG Orthology (KO)" and "ko01000 Enzymes", should have "accessions" for each category
                parsed_categories = []
                for category in categories:
                    split_category = category.split(' ')
                    parsed_categories.append((split_category[0], ' '.join(split_category[1: ])))
            else:
                parsed_categories = [('', category) for category in categories]
            category_list.append(parsed_categories)

        return ortholog_dict


    def get_brite_table_as_hierarchy_dict(self, hierarchy_accessions_of_interest=None, level_cutoff=None, collapse_keys=False, collapse_mixed_branches=True):
        """Load the BRITE hierarchies table as a dictionary keyed by hierarchy.

        The returned dictionary contains the category structure of the hierarchy and a set of
        orthologs in each categorization.

        With `collapse_keys` set to the default of False and `collapse_mixed_branches` set to the
        default of True, the returned dictionary is structured as follows:
            {
                (<hierarchy A accession>, <hierarchy A name>):
                    {
                        (<level 1 category i accession>, <level 1 category i name>):
                            {
                                (<level 2 category j accession>, <level 2 category j name>):
                                    set([(<ortholog 1 accession>, <ortholog 1 name>), (<ortholog 2 accession>, <ortholog 2 name>), ...]),
                                (<level 2 category k accession>, <level 2 category k name>):
                                    set([(<ortholog 3 accession>, <ortholog 3 name>), (<ortholog 4 accession>, <ortholog 4 name>), ...]),
                                (<level 2 category l accession>, <level 2 category l name>):
                                    {...},
                                ...
                            },
                        (<level 1 category m accession>, <level 1 category m name>):
                            {...},
                        ...
                    },
                    {
                        (<level 1 category n accession>, <level 1 category n name>):
                            {...},
                        ...
                    },
                (<hierarchy B accession>, <hierarchy B name>):
                    {...},
                ...
            }

        Here is an example of the entry for the "Ribosome" hierarchy:
            ('ko03011', 'Ribosome'):
                {
                    ('', 'Ribosomal proteins'):
                        {
                            ('', 'Eukaryotes'):
                                {
                                    ('', 'Small subunit'):
                                        set([('K02981', 'RP-S2e, RPS2; small subunit ribosomal protein S2e'), ('K02985', 'RP-S3e, RPS3; small subunit ribosomal protein S3e'), ...]),
                                    ('', 'Large subunit'):
                                        set([('K02925', 'RP-L3e, RPL3; large subunit ribosomal protein L3e'), ('K02930', 'RP-L4e, RPL4; large subunit ribosomal protein L4e'), ...])
                                },
                            ('', 'Mitochondria/ Chloroplast'):
                                {
                                    ('', 'Small subunit'):
                                        {...}
                                    ('', 'Large subunit'):
                                        {...}
                                },
                            ('', 'Bacteria'):
                                {
                                    ...
                                },
                            ('', 'Archaea'):
                                {
                                    ...
                                }
                        },
                    ('', 'Ribosomal RNAs'):
                        {
                            ('Eukaryotes'):
                                set([('K01979', 'SSUrRNA; small subunit ribosomal RNA'), ('K01982', 'LSUrRNA; large subunit ribosomal RNA'), ...]),
                            ('Prokaryotes'):
                                set([('K01985', '5SrRNA, rrf; 5S ribosomal RNA'), ('K01977', '16SrRNA, rrs; 16S ribosomal RNA'), ('K01980', '23SrRNA, rrl; 23S ribosomal RNA')])
                        }
                }

        Keys and set items are split by accession and description, even in the absence of an
        accession for a category in the hierarchy, as seen in the example. The only hierarchies
        expected to contain category "accessions" are "ko00001 Gene Ontology (KO)" and "ko01000
        Enzymes".

        Levels of the hierarchy can be collapsed with the `level_cutoff` argument. If `level_cutoff`
        is a positive number, it is measured down from the top level of the hierarchy (level 1). If
        `level_cutoff` is a negative number, it is measured up from the bottom-most level of the
        hierarchy. Since different branches of the hierarchy can have different depths, the negative
        number is converted to a positive number given the deepest branch: in the "Ribosome"
        example, the "Ribosomal proteins" branch has 3 levels, and the "Ribosomal RNAs" branch has 2
        levels, level -2 would be measured against the "Ribosomal proteins" branch and be converted
        to level 1. If the negative parameterization would eliminate all levels of the hierarchy,
        the cutoff is set to level 1.

        Example: `level_cutoff` is set to 1 or -2, so levels below "Ribosomal proteins" and
        "Ribosomal RNAs" are removed:
            ('ko03011', 'Ribosome'):
                {
                    ('', 'Ribosomal proteins'):
                        set([('K02981', 'RP-S2e, RPS2; small subunit ribosomal protein S2e'), ('K02985', 'RP-S3e, RPS3; small subunit ribosomal protein S3e'), ...]),
                    ('', 'Ribosomal RNAs'):
                        set([('K01979', 'SSUrRNA; small subunit ribosomal RNA'), ('K01982', 'LSUrRNA; large subunit ribosomal RNA'), ...])
                }

        Example: `level_cutoff` is set to -1, only removing "Small subunit" and "Large subunit"
        levels under "Ribosomal proteins" but not levels under "Ribosomal RNAs".
            ('ko03011', 'Ribosome'):
                {
                    ('', 'Ribosomal proteins'):
                        {
                            ('', 'Eukaryotes'):
                                set([('K02981', 'RP-S2e, RPS2; small subunit ribosomal protein S2e'), ('K02985', 'RP-S3e, RPS3; small subunit ribosomal protein S3e'), ...]),
                            ('', 'Mitochondria/ Chloroplasts'):
                                set([...]),
                            ('', 'Bacteria'):
                                set([...]),
                            ('', 'Archaea'):
                                set([...])
                        },
                    ('', 'Ribosomal RNAs'):
                        {
                            ('', 'Eukaryotes'):
                                set([...]),
                            ('', 'Prokaryotes'):
                                set([...])
                        }

        Dictionary nesting can be simplified by `collapse_keys`. Each ortholog set is keyed by a
        single categorization tuple. Applied to the "Ribosome" hierarchy:
            ('ko03011', 'Ribosome'):
                {
                    (('', 'Ribosomal proteins'), ('', 'Eukaryotes'), ('', 'Small subunit')):
                        set([('K02981', 'RP-S2e, RPS2; small subunit ribosomal protein S2e'), ('K02985', 'RP-S3e, RPS3; small subunit ribosomal protein S3e'), ...]),
                    (('', 'Ribosomal proteins'), ('', 'Eukaryotes'), ('', 'Large subunit')):
                        set([('K02925', 'RP-L3e, RPL3; large subunit ribosomal protein L3e'), ('K02930', 'RP-L4e, RPL4; large subunit ribosomal protein L4e'), ...]),
                    (('', 'Ribosomal proteins'), ('', 'Mitochondria/ Chloroplasts'), ('', 'Small subunit')):
                        set([...]),
                    ...
                }

        To this point, we have ignored the possibility that a category containing ortholog "leaves"
        can also contain additional category "branches". An example of this is in the "RNases"
        category of the "Ribosome biogenesis" hierarchy. The category contains RNases such as
        "K14812  NGL2; RNA exonuclease NGL2 [EC:3.1.-.-]", but also includes a category of "RNase
        MRP" subunits, including "K01164  POP1; ribonuclease P\/MRP protein subunit POP1
        [EC:3.1.26.5]". With the parameter, `collapse_mixed_branches`, set to the default of True,
        categories in such "mixed" branches are collapsed out of existence: subunit orthologs are
        placed in "RNases" rather than "RNase MRP", which is removed.

        Mixed branches can be preserved by setting `collapse_mixed_branches` to False. This also
        changes the structure of the returned dictionary. Now, rather than a category key mapping to
        EITHER a dict or a set, each category key maps to a tuple of (1) an ortholog set and (2) a
        category dict. Without `collapse keys`, the returned dictionary is structured as follows:
            {
                (<hierarchy A accession>, <hierarchy A name>):
                    {
                        (<level 1 category i accession>, <level 1 category i name>):
                            (
                                set([...]),
                                {
                                    (<level 2 category j accession>, <level 2 category j name>):
                                        (
                                            set([(<ortholog 1 accession>, <ortholog 1 name>), (<ortholog 2 accession>, <ortholog 2 name>), ...]),
                                            {...}
                                        ),
                                    (<level 2 category k accession>, <level 2 category k name>):
                                        (
                                            set([(<ortholog 3 accession>, <ortholog 3 name>), (<ortholog 4 accession>, <ortholog 4 name>), ...]),
                                            {}
                                        ),
                                    (<level 2 category l accession>, <level 2 category l name>):
                                        (
                                            set([]),
                                            {...},
                                        )
                                    ...
                                }
                            ),
                        (<level 1 category m accession>, <level 1 category m name>):
                            (
                                set([...]),
                                {...}
                            ),
                        ...
                (<hierarchy B accession>, <hierarchy B name>):
                    {...},
                ...
            }

        With `collapse_keys` and without `collapse_mixed_branches`, the format of the returned
        dictionary is the same as with `collapse_keys` and `collapse_mixed_branches`. The only
        difference is that there can be entries for orthologs in a category and orthologs in a
        subcategory: the categorization key tuples are the same for such entries up to the
        subcategory elements.

        PARAMETERS
        ==========
        hierarchy_accessions_of_interest : list, None
            filters results to hierarchies of interest

        level_cutoff : int, None
            collapse branches below the level cutoff, with a positive level measured top-down and a
            negative level measured bottom-up from the deepest branch of the tree

        collapse_keys: bool, False
            eliminate category nesting, keying each ortholog set by a single categorization tuple

        collapse_mixed_branches : bool, True
            collapse category branches off a category node that also ends in ortholog leaves

        RETURNS
        =======
        hierarchy_dict : dict
            dictionary of BRITE hierarchies
        """

        if hierarchy_accessions_of_interest:
            hierarchy_list = [f"'{konum}'" for konum in hierarchy_accessions_of_interest]
            where_clause_string = f"hierarchy_accession IN ({','.join(hierarchy_list)})"
            # this WILL fail if you ask for a data name that doesn't exist, so know your data before you query
            dict_from_brite_table = self.db.get_some_rows_from_table_as_dict(self.brite_table_name, where_clause_string, row_num_as_key=True)
        else:
            dict_from_brite_table = self.db.get_table_as_dict(self.brite_table_name, row_num_as_key=True)

        # the returned dict is keyed by an arbitrary integer, and each value is a dict containing one row from the BRITE table, e.g.,
        # 0: {'hierarchy_accession': 'ko01000',
        #     'hierarchy_name': 'Enzymes',
        #     'ortholog_accession': 'K00001',
        #     'ortholog_name': 'E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]',
        #     'categorization': '1. Oxidoreductases>>>1.1  Acting on the CH-OH group of donors>>>1.1.1  With NAD+ or NADP+ as acceptor>>>1.1.1.1  alcohol dehydrogenase'}

        # find the maximum depth of each hierarchy
        max_depth_dict = self.get_brite_max_depth_dict(dict_from_brite_table)

        if level_cutoff == 0 or type(level_cutoff) != int:
            raise ConfigError("`level_cutoff` must be a nonzero integer.")

        # set the level cutoff for each hierarchy
        if level_cutoff > 0:
            topdown_level_cutoff_dict = {hierarchy_accession: min(level_cutoff, max_depth) for hierarchy_accession, max_depth in max_depth_dict.items()}
        elif level_cutoff < 0:
            # find the positive level corresponding to the negative level cutoff for each
            # hierarchy, ensuring that at least one category remains per hierarchy
            topdown_level_cutoff_dict = {hierarchy_accession: max(max_depth + level_cutoff, 1) for hierarchy_accession, max_depth in max_depth_dict.items()}
        else:
            topdown_level_cutoff_dict = max_depth_dict

        # hierarchy level cutoffs can be affected by collapsing subcategories of mixed categories
        if collapse_mixed_branches:
            topdown_level_cutoff_dict = self.get_brite_topdown_level_cutoff_dict_ignoring_subcategories_of_mixed_categories(topdown_level_cutoff_dict, dict_from_brite_table)

        # create the per-hierarchy dict
        hierarchy_dict = {}
        for row_id, entry_dict in dict_from_brite_table.items():
            ortholog_accession = entry_dict['ortholog_accession']
            ortholog_name = entry_dict['ortholog_name']
            hierarchy_accession = entry_dict['hierarchy_accession']
            hierarchy_name = entry_dict['hierarchy_name']
            categorization = entry_dict['categorization']
            categories = categorization.split('>>>')

            # the hierarchies, "ko00001 KEGG Orthology (KO)" and "ko01000 Enzymes", should have "accessions" for each category
            if hierarchy_accession == 'ko00001' or hierarchy_accession == 'ko01000':
                parsed_categories = []
                for category in categories:
                    split_category = category.split(' ')
                    parsed_categories.append((split_category[0], ' '.join(split_category[1: ])))
            else:
                parsed_categories = [('', category) for category in categories]

            # make the top-level category dict for the hierarchy
            try:
                category_dict = hierarchy_dict[(hierarchy_accession, hierarchy_name)]
            except KeyError:
                hierarchy_dict[(hierarchy_accession, hierarchy_name)] = category_dict = {}

            topdown_level_cutoff = topdown_level_cutoff_dict[hierarchy_accession]

            if collapse_mixed_branches:
                # each value of a category dict is either a set or a dict
                if collapse_keys:
                    key = tuple(parsed_categories[: topdown_level_cutoff])
                    try:
                        ortholog_set = category_dict[key]
                    except KeyError:
                        category_dict[key] = ortholog_set = set()
                    ortholog_set.add((ortholog_accession, ortholog_name))
                else:
                    num_categories = len(parsed_categories)
                    for level, category in enumerate(parsed_categories, 1):
                        if level == topdown_level_cutoff:
                            try:
                                ortholog_set = category_dict[category]
                            except KeyError:
                                category_dict[category] = ortholog_set = set()
                            ortholog_set.add((ortholog_accession, ortholog_name))
                            break
                        try:
                            category_dict = category_dict[category]
                        except KeyError:
                            category_dict[category] = category_dict = {}
            else:
                # each value of a category dict is a tuple containing a set and a dict
                if collapse_keys:
                    key = tuple(parsed_categories[: topdown_level_cutoff])
                    try:
                        ortholog_set = category_dict[key]
                    except KeyError:
                        category_dict[key] = ortholog_set = set()
                    ortholog_set.add((ortholog_accession, ortholog_name))
                else:
                    num_categories = len(parsed_categories)
                    for level, category in enumerate(parsed_categories, 1):
                        if level == topdown_level_cutoff:
                            try:
                                category_tuple = category_dict[category]
                            except KeyError:
                                category_dict[category] = category_tuple = (set(), {})
                            category_tuple[0].add((ortholog_accession, ortholog_name))
                            break
                        try:
                            category_dict = category_dict[category][1]
                        except KeyError:
                            category_dict[category] = category_tuple = (set(), {})
                            category_dict = category_tuple[1]

        return hierarchy_dict


    def get_brite_max_depth_dict(self, dict_from_brite_table=None):
        """Return a dictionary of the maximum depths of BRITE hierarchies.

        "Maximum depth" is greatest number of levels of categorization of an ortholog in the
        hierarchy.

        By default, without `dict_from_brite_table`, every hierarchy in the database is analyzed.
        With that argument, only the hierarchies represented in the dict are analyzed.

        PARAMETERS
        ==========
        dict_from_brite_table : dict
            contains BRITE table rows of interest, as returned, for example, by `get_some_rows_from_table_as_dict`

        RETURNS
        =======
        max_depth_dict : dict
            relates hierarchy accession keys to maximum depths
        """
        if not dict_from_brite_table:
            dict_from_brite_table = self.db.get_table_as_dict(self.brite_table_name, row_num_as_key=True)

        max_depth_dict = {}
        for entry_dict in dict_from_brite_table.values():
            hierarchy_accession = entry_dict['hierarchy_accession']
            categories = entry_dict['categorization'].split('>>>')
            try:
                current_max_depth = max_depth_dict[hierarchy_accession]
            except KeyError:
                current_max_depth = 0
            max_depth_dict[hierarchy_accession] = max(current_max_depth, len(categories))

        return max_depth_dict


    def get_brite_depth_dict_ignoring_subcategories_of_mixed_categories(self, dict_from_brite_table=None, input_depth_dict=None):
        """Return hierarchy depths disregarding bottom levels of hierarchy that only exist due to subcategories of mixed categories.

        "Mixed" categories contain both subcategories and orthologs.

        Example: The depth of a hierarchy is 4, but this is only due to a subcategory of a depth 3
        mixed category. After collapsing this subcategory, hierarchy depth is reduced to 3.

        With default arguments, all BRITE table entries and maximum hierarchy depths are considered.

        PARAMETERS
        ==========
        dict_from_brite_table : dict, None
            contains BRITE table rows of interest, as returned, for example, by
            `self.db.get_some_rows_from_table_as_dict`

        input_depth_dict : dict, None
            contains hierarchy depths of interest, as returned, for example, by
            `self.get_brite_max_depth_dict`

        RETURNS
        =======
        depth_dict : dict, None
            contains adjusted hierarchy depths
        """
        if not dict_from_brite_table:
            dict_from_brite_table = self.db.get_table_as_dict(self.brite_table_name, row_num_as_key=True)

        if input_depth_dict:
            depth_dict = copy.deepcopy(input_depth_dict)
        else:
            depth_dict = self.get_brite_max_depth_dict(dict_from_brite_table)

        deep_categorizations = {} # depth of categorization greater than or equal to level cutoff
        shallow_categorizations = {} # depth of categorization less than level cutoff
        hierarchy_culprits = set([]) # hierarchies containing subcategories of mixed categories responsible for depth exceeding level cutoff
        is_checked = False
        # keep trimming hierarchy "culprits" until the deepest level is not a subcategory of a mixed category
        while hierarchy_culprits or not is_checked:
            for entry_dict in dict_from_brite_table.values():
                hierarchy_accession = entry_dict['hierarchy_accession']
                if is_checked and hierarchy_accession not in hierarchy_culprits:
                    continue
                categorization = entry_dict['categorization']
                categories = categorization.split('>>>')
                if len(categories) >= depth_dict[hierarchy_accession]:
                    try:
                        deep_categorizations[hierarchy_accession].add(categorization)
                    except KeyError:
                        deep_categorizations[hierarchy_accession] = set([categorization])
                else:
                    try:
                        shallow_categorizations[hierarchy_accession].add(categorization)
                    except KeyError:
                        shallow_categorizations[hierarchy_accession] = set([categorization])

            # interrogate the "deep" categorizations in each hierarchy
            for hierarchy_accession, hierarchy_deep_categorizations in deep_categorizations.items():
                # if all categorizations are deep, then the hierarchy cannot be a "culprit"
                try:
                    hierarchy_shallow_categorizations = shallow_categorizations[hierarchy_accession]
                except KeyError:
                    continue

                # compare deep categorization strings with all shallow categorizations
                # strings, checking if the shallow categorization is a substring of the deep
                # categorization, and thus the deep categorization is a subcategory of a
                # mixed category
                subcategorization_culprits = []
                for deep_categorization in hierarchy_deep_categorizations:
                    for shallow_categorization in hierarchy_shallow_categorizations:
                        if deep_categorization[: len(shallow_categorization)] == shallow_categorization:
                            if deep_categorization[len(shallow_categorization): len(shallow_categorization) + 3] == '>>>':
                                subcategorization_culprits.append(deep_categorization)
                                break

                # remove culprit subcategories
                for deep_categorization in subcategorization_culprits:
                    hierarchy_deep_categorizations.remove(deep_categorization)

                hierarchy_culprits.add(hierarchy_accession)
                if hierarchy_deep_categorizations:
                    # there are no more subcategory culprits, but the hierarchy is still
                    # deeper than the level cutoff: the deep categories will be collapsed as
                    # per normal
                    hierarchy_culprits.remove(hierarchy_accession)
                else:
                    # there are no more subcategory culprits, so the level cutoff can simply be reduced by 1
                    depth_dict[hierarchy_accession] = depth_dict[hierarchy_accession] - 1
            is_checked = True

        return depth_dict


    def list_brite_hierarchies(self, as_accessions=False, as_tuples=False):
        """List all BRITE hierarchies in the database.

        PARAMETERS
        ==========
        as_accessions : bool, False
            return list of hierarchy accessions

        as_tuples : bool, False
            return list of tuples of hierarchy accessions and names

        RETURNS
        =======
        hierarchy_entries : list
            database BRITE hierarchies, formatted "<accession> <name>" with default parameterization
        """

        if as_accessions:
            hierarchy_entries = self.db.get_single_column_from_table(self.brite_table_name, 'hierarchy_accession', unique=True)
        else:
            hierarchy_entries = self.db.get_some_columns_from_table(self.brite_table_name, 'hierarchy_accession, hierarchy_name', unique=True)
            if not as_tuples:
                hierarchy_entries = [f"{accession} {name}" for accession, name in hierarchy_entries]

        return hierarchy_entries


######### MODULE DEFINITION UNROLLING FUNCTIONS #########

    def get_top_level_steps_in_module_definition(self, mnum):
        """This function access the DEFINITION line of a KEGG Module and returns the top-level steps as a list

        A 'top-level' step is one that you get by splitting on spaces (but not spaces in parentheses) just once -
        ie, the 'first layer' when unrolling the module.
        """

        def_string = self.get_kegg_module_definition(mnum)
        return utils.split_by_delim_not_within_parens(def_string, " ")


    def unroll_module_definition(self, mnum, def_lines = None):
        """This function accesses the DEFINITION line of a KEGG Module, unrolls it into all possible paths through the module, and
        returns the list of all paths.

        This is a driver for the recursive functions that do the actual unrolling of each definition line.

        PARAMETERS
        ==========
        mnum : str
            module number
        def_lines : list of str
            The DEFINITION lines for the module. This parameter is optional, and if it is not passed, the module
            definition will be looked up from the modules DB.
        """

        if not def_lines:
            def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        combined_def_line = ""
        for d in def_lines:
            d = d.strip()
            combined_def_line += d + " "
        combined_def_line = combined_def_line.strip()
        def_line_paths = self.recursive_definition_unroller(combined_def_line)

        return def_line_paths


    def recursive_definition_unroller(self, step):
        """This function recursively splits a module definition into its components.

        First, the definition is split into its component steps (separated by spaces).
        Each step is either an atomic step (a single KO, module number, '--', or nonessential KO starting with '-'),
        a protein complex, or a compound step.

        Atomic steps are used to extend each path that has been found so far. Protein complexes are split into
        their respective components, which may be split further by the split_paths() function to find all possible
        alternative complexes, before being used to extend each path. Compound steps are split and recursively processed
        by the split_paths() function before the resulting downstream paths are used to extend each path.

        PARAMETERS
        ==========
        step : str
            step definition to split into component steps as necessary

        RETURNS
        =======
        paths_list : list
            all paths that the input step has been unrolled into
        """

        split_steps = utils.split_by_delim_not_within_parens(step, " ")
        paths_list = [[]]  # list to save all paths, with initial empty path list to extend from
        for s in split_steps:
            # base case: step is a ko, mnum, non-essential step, or '--'
            if (len(s) == 6 and s[0] == "K") or (len(s) == 6 and s[0] == "M") or (s == "--") or (len(s) == 7 and s[0] == "-"):
                for p in paths_list:
                    p.extend([s])
            else:
                if s[0] == "(" and s[-1] == ")":
                    # here we try splitting to see if removing the outer parentheses will make the definition become unbalanced
                    # (the only way to figure this out is to try it because regex cannot handle nested parentheses)
                    comma_substeps = utils.split_by_delim_not_within_parens(s[1:-1], ",")
                    if not comma_substeps: # if it doesn't work, try without removing surrounding parentheses
                        comma_substeps = utils.split_by_delim_not_within_parens(s, ",")
                    space_substeps = utils.split_by_delim_not_within_parens(s[1:-1], " ")
                    if not space_substeps:
                        space_substeps = utils.split_by_delim_not_within_parens(s, " ")
                else:
                    comma_substeps = utils.split_by_delim_not_within_parens(s, ",")
                    space_substeps = utils.split_by_delim_not_within_parens(s, " ")

                # complex case: no commas OR spaces outside parentheses so this is a protein complex rather than a compound step
                if len(comma_substeps) == 1 and len(space_substeps) == 1:
                    complex_components, delimiters = utils.split_by_delim_not_within_parens(s, ["+","-"], return_delims=True)
                    complex_strs = [""]

                    # reconstruct the complex (and any alternate possible complexes) while keeping the +/- structure the same
                    for i in range(len(complex_components)):
                        c = complex_components[i]
                        if c[0] == '(':
                            alts = self.split_path(c)
                            new_complex_strs = []
                            for a in alts:
                                if len(a) > 1:
                                    raise ConfigError("Uh oh. recursive_definition_unroller() speaking. We found a protein complex with more "
                                                      "than one KO per alternative option here: %s" % s)
                                for cs in complex_strs:
                                    extended_complex = cs + a[0]
                                    new_complex_strs.append(extended_complex)
                            complex_strs = new_complex_strs
                        else:
                            for j in range(len(complex_strs)):
                                complex_strs[j] += c

                        if i < len(delimiters):
                            for j in range(len(complex_strs)):
                                complex_strs[j] += delimiters[i]

                    new_paths_list = []
                    for cs in complex_strs:
                        for p in paths_list:
                            p_copy = copy.copy(p)
                            p_copy.extend([cs])
                            new_paths_list.append(p_copy)
                    paths_list = new_paths_list

                # compound step case:
                else:
                    alts = self.split_path(s)
                    new_paths_list = []
                    for a in alts:
                        for p in paths_list:
                            p_copy = copy.copy(p)
                            p_copy.extend(a)
                            new_paths_list.append(p_copy)
                    paths_list = new_paths_list

        return paths_list


    def split_path(self, step):
        """This function handles compound steps that should be split into multiple alternative paths.

        It first splits the input step into substeps, and then since each substep could be its own mini-definition,
        it recursively calls the definition unrolling function to parse it. The list of all alternative paths
        that can be made from this step is returned.
        """

        if step[0] == "(" and step[-1] == ")":
            substeps = utils.split_by_delim_not_within_parens(step[1:-1], ",")
            if not substeps: # if it doesn't work, try without removing surrounding parentheses
                substeps = utils.split_by_delim_not_within_parens(step, ",")
        else:
            substeps = utils.split_by_delim_not_within_parens(step, ",")

        alt_path_list = []
        for s in substeps:
            alt_paths_from_substep = self.recursive_definition_unroller(s)
            for a in alt_paths_from_substep:
                alt_path_list.append(a)

        return alt_path_list


class ModulesTable:
    """This class defines operations for creating the KEGG Modules table in Modules.db"""

    def __init__(self, mod_table_name = None):
        """"""
        self.db_entries = []
        self.total_entries = 0

        if mod_table_name:
            self.module_table_name = mod_table_name
        else:
            raise ConfigError("Beep Beep. Warning. ModulesTable was initialized without knowing its own name.")


    def append_and_store(self, db, module_num, data_name, data_value, data_definition=None, line_num=None):
        """This function handles collects db entries (as tuples) into a list, and once we have 10,000 of them it stores that set into the Modules table.

        The db_entries list is cleared after each store so that future stores don't add duplicate entries to the table.
        """

        db_entry = tuple([module_num, data_name, data_value, data_definition, line_num])
        self.db_entries.append(db_entry)
        self.total_entries += 1

        # we can store chunks of 5000 at a time, so we don't want over 10,000 entries.
        if len(self.db_entries) >= 10000:
            self.store(db)
            self.db_entries = []


    def store(self, db):
        if len(self.db_entries):
            db._exec_many('''INSERT INTO %s VALUES (%s)''' % (self.module_table_name, (','.join(['?'] * len(self.db_entries[0])))), self.db_entries)


    def get_total_entries(self):
        return self.total_entries


class BriteTable:
    """This class defines operations for creating the KEGG BRITE hierarchies table in Modules.db"""

    def __init__(self, brite_table_name = None):
        """"""
        self.db_entries = []
        self.total_entries = 0

        if brite_table_name:
            self.brite_table_name = brite_table_name
        else:
            raise ConfigError("Beep Beep. Warning. BriteTable was initialized without knowing its own name.")


    def append_and_store(self, db, hierarchy_accession, hierarchy_name, ortholog_accession, ortholog_name, categorization):
        """This function handles collects db entries (as tuples) into a list, and once we have 10,000 of them it stores that set into the Modules table.

        The db_entries list is cleared after each store so that future stores don't add duplicate entries to the table.
        """

        db_entry = tuple([hierarchy_accession, hierarchy_name, ortholog_accession, ortholog_name, categorization])
        self.db_entries.append(db_entry)
        self.total_entries += 1

        # we can store chunks of 5000 at a time, so we don't want over 10,000 entries.
        if len(self.db_entries) >= 10000:
            self.store(db)
            self.db_entries = []


    def store(self, db):
        if len(self.db_entries):
            db._exec_many('''INSERT INTO %s VALUES (%s)''' % (self.brite_table_name, (','.join(['?'] * len(self.db_entries[0])))), self.db_entries)


    def get_total_entries(self):
        return self.total_entries


class KeggModuleEnrichment(KeggContext):
    """This class is a driver for anvi-script-enrichment-stats for modules input.

    It takes in the modules mode output from anvi-estimate-metabolism, formats it for the enrichment script,
    and runs the script.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-compute-functional-enrichment
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.modules_txt = A('modules_txt')
        self.groups_txt = A('groups_txt')
        self.sample_header_in_modules_txt = A('sample_header') or 'db_name'
        self.module_completion_threshold = A('module_completion_threshold') or 0.75
        self.output_file_path = A('output_file')
        self.include_missing = True if A('include_samples_missing_from_groups_txt') else False
        self.use_stepwise_completeness = A('use_stepwise_completeness')

        # init the base class
        KeggContext.__init__(self, self.args)

        # if necessary, assign 0 completion threshold, which evaluates to False above
        if A('module_completion_threshold') == 0:
            self.module_completion_threshold = 0.0
            if not self.quiet:
                self.run.warning("Your completion threshold is set to 0, which will make the enrichment results MEANINGLESS. Why? Because "
                                 "with a threshold this low, every module will be considered present in every single sample and therefore will "
                                 "be equally present in every single group. So you should stop what you are doing RIGHT NOW.")
            if not self.just_do_it:
                raise ConfigError("We are stopping you right there because your completion threshold is 0 and that will make the enrichment results "
                                  "meaningless (see the warnings above, if you haven't suppressed them with --quiet). But if you really really really "
                                  "want to do it, you can run again with --just-do-it and then we won't stop you. You have the right to be without meaning.")

        # sanity checkses my precious
        if not self.modules_txt:
            raise ConfigError("To compute module enrichment, you must provide a modules-txt file (aka modules mode output from "
                              "`anvi-estimate-metabolism`).")
        if not self.groups_txt:
            raise ConfigError("To compute module enrichment, you must provide a groups-txt file mapping each sample to a group.")

        filesnpaths.is_file_exists(self.modules_txt)
        filesnpaths.is_file_plain_text(self.modules_txt)
        filesnpaths.is_file_exists(self.groups_txt)
        filesnpaths.is_file_plain_text(self.groups_txt)

        if filesnpaths.is_file_exists(self.output_file_path, dont_raise=True):
            raise ConfigError(f"Whoops... we almost overwrote the existing output file {self.output_file_path}. But we stopped just in time. "
                               "If you really want us to replace the contents of that file with new enrichment results, then remove this "
                               "file before you run this program again.")
        filesnpaths.is_output_file_writable(self.output_file_path)

        if not self.quiet:
            self.run.info("modules-txt input file", self.modules_txt)
            self.run.info("groups-txt input file", self.groups_txt)
            self.run.info("sample column in modules-txt", self.sample_header_in_modules_txt)
            self.run.info("module completion threshold", self.module_completion_threshold)


    def get_enrichment_input(self, output_file_path):
        """This function converts modules mode output into input for anvi-script-enrichment-stats

        The input format for anvi-script-enrichment-stats is described in a comment at the top of that script, and here is
        how we get the values for each column:
        The first column, 'KEGG_MODULE', and second column 'accession', are already in the modules mode output as 'module_name'
        and 'module', respectively.
        The 'N_*' columns are the total number of samples in each group.
        For each module, this function determines which samples the module is 'present' in according to the specified completion threshold.
        This determines the list of samples for the 'sample_ids' column as well as the 'p_*' proportions for each group of samples.
        Finally, the fourth column, 'associated_groups', is computed from the 'p_*' proportions and 'N_*' totals.

        PARAMETERS
        ==========
        output_file_path : str
            a file path where we will store the (temporary) input file for the enrichment script
        """

        filesnpaths.is_output_file_writable(output_file_path)

        # read the files into dataframes
        modules_df = pd.read_csv(self.modules_txt, sep='\t')

        completeness_header = 'pathwise_module_completeness'
        if self.use_stepwise_completeness:
            completeness_header = 'stepwise_module_completeness'
        self.progress.reset()
        self.run.info("Completeness score being used for determining sample presence", completeness_header)

        # make sure we have all the columns we need in modules mode output, since this output can be customized
        required_modules_txt_headers = ['module', 'module_name', completeness_header]
        missing_headers = []
        for h in required_modules_txt_headers:
            if h not in modules_df.columns:
                missing_headers.append(h)
        if missing_headers:
            missing_string = ", ".join(missing_headers)
            self.progress.reset()
            raise ConfigError("We cannot go on! *dramatic sweep*   We trust that you have provided us with "
                              "modules mode output, but unfortunately the modules-txt input does not contain "
                              f"the following required headers: {missing_string}   Please re-generate your "
                              "modules-txt to include these before trying again.")

        if 'unique_id' in modules_df.columns:
            modules_df = modules_df.drop(columns=['unique_id'])

        # samples column sanity check - this column will become the index
        if self.sample_header_in_modules_txt not in modules_df.columns:
            col_list = ", ".join(modules_df.columns)
            self.progress.reset()
            raise ConfigError(f"You have specified that your sample names are in the column with header '{self.sample_header_in_modules_txt}' "
                               "in the modules-txt file, but that column does not exist. :( Please figure out which column is right and submit "
                               "it using the --sample-header parameter. Just so you know, the columns in modules-txt that you can choose from "
                               f"are: {col_list}")

        samples_to_groups_dict, groups_to_samples_dict = utils.get_groups_txt_file_as_dict(self.groups_txt, include_missing_samples_is_true=self.include_missing)

        # make sure the samples all have a group
        samples_with_none_group = []
        for s,g in samples_to_groups_dict.items():
            if not g:
                samples_with_none_group.append(s)

        for s in samples_with_none_group:
            samples_to_groups_dict.pop(s)

        if samples_with_none_group:
            self.progress.reset()
            none_group_str = ", ".join(samples_with_none_group)
            self.run.warning("Some samples in your groups-txt did not have a group, and we will ignore those samples. If you "
                                 "want them to be included in the analysis, you need to fix the groups-txt to have a group for "
                                 "these samples. Anyway. Here are the samples we will be ignoring: "
                                 f"{none_group_str}")

        # sanity check for mismatch between modules-txt and groups-txt
        sample_names_in_modules_txt = set(modules_df[self.sample_header_in_modules_txt].unique())
        sample_names_in_groups_txt = set(samples_to_groups_dict.keys())
        samples_missing_in_groups_txt = sample_names_in_modules_txt.difference(sample_names_in_groups_txt)
        samples_missing_in_modules_txt = sample_names_in_groups_txt.difference(sample_names_in_modules_txt)
        if anvio.DEBUG:
            self.run.info("Samples in modules-txt", ", ".join(list(sample_names_in_modules_txt)))
            self.run.info("Samples in groups-txt", ", ".join(list(sample_names_in_groups_txt)))
            self.run.info("Missing samples from groups-txt", ", ".join(list(samples_missing_in_groups_txt)))
            self.run.info("Missing samples from modules-txt", ", ".join(list(samples_missing_in_modules_txt)))

        if samples_missing_in_groups_txt:
            missing_samples_str = ", ".join(samples_missing_in_groups_txt)
            if not self.include_missing:
                self.progress.reset()
                self.run.warning(f"Your groups-txt file does not contain some samples present in your modules-txt ({self.sample_header_in_modules_txt} "
                                "column). Since you have not elected to --include-samples-missing-from-groups-txt, we are not going to take these samples into consideration at all. "
                                "Here are the samples that we will be ignoring: "
                                f"{missing_samples_str}")
                # drop the samples that are not in groups-txt
                modules_df = modules_df[~modules_df[self.sample_header_in_modules_txt].isin(list(samples_missing_in_groups_txt))]
                if anvio.DEBUG:
                    self.run.info("Samples remaining in modules-txt dataframe after removing ungrouped", ", ".join(modules_df[self.sample_header_in_modules_txt].unique()))

            else:
                self.progress.reset()
                self.run.warning(f"Your groups-txt file does not contain some samples present in your modules-txt ({self.sample_header_in_modules_txt} "
                                "column). Since you have chosen to --include-samples-missing-from-groups-txt, for the purposes of this analysis we will now consider all of "
                                "these samples to belong to one group called 'UNGROUPED'."
                                "Here are the {len(samples_missing_in_groups_txt)} UNGROUPED samples that we will consider as one big happy family: "
                                f"{missing_samples_str}")
                # add those samples to the UNGROUPED group
                ungrouped_samples = list(samples_missing_in_groups_txt)
                for s in ungrouped_samples:
                    samples_to_groups_dict[s] = 'UNGROUPED'

        if samples_missing_in_modules_txt:
            missing_samples_str = ", ".join(samples_missing_in_modules_txt)
            if not self.just_do_it:
                self.progress.reset()
                raise ConfigError(f"Your modules-txt file ({self.sample_header_in_modules_txt} column) does not contain some samples that "
                                 "are present in your groups-txt. This is not necessarily a huge deal, it's just that those samples will "
                                 "not be included in the enrichment analysis because, well, you don't have any module information for them. "
                                 "If all of the missing samples belong to groups you don't care about at all, then feel free to ignore this "
                                 "message and re-run using --just-do-it. But if you do care about those groups, you'd better fix this because "
                                 "the enrichment results for those groups will be wrong. Here are the samples in question: "
                                  f"{missing_samples_str}")
            else:
                self.progress.reset()
                self.run.warning(f"Your modules-txt file ({self.sample_header_in_modules_txt} column) does not contain some samples that "
                                 "are present in your groups-txt. This is not necessarily a huge deal, it's just that those samples will "
                                 "not be included in the enrichment analysis because, well, you don't have any module information for them. "
                                 "Since you have used the --just-do-it parameter, we assume you don't care about this and are going to keep "
                                 "going anyway. We hope you know what you are doing :) Here are the samples in question: "
                                  f"{missing_samples_str}")
                # drop the samples that are not in modules-txt
                for s in list(samples_missing_in_modules_txt):
                    samples_to_groups_dict.pop(s)
                if anvio.DEBUG:
                    self.run.info("Samples remaining in groups-txt dataframe after removing ungrouped", ", ".join(samples_to_groups_dict.keys()))


        modules_df.set_index(self.sample_header_in_modules_txt, inplace=True)
        sample_groups_df = pd.DataFrame.from_dict(samples_to_groups_dict, orient="index", columns=['group'])

        # convert modules mode output to enrichment input
        N_values = sample_groups_df['group'].value_counts()
        group_list = N_values.keys()
        module_list = modules_df['module'].unique()

        output_dict = {}
        header_list = ['MODULE', 'accession', 'sample_ids', 'associated_groups']
        for c in group_list:
            header_list.append(f"p_{c}")
            header_list.append(f"N_{c}")

        for mod_num in module_list:
            query_string = f"module == '{mod_num}' and {completeness_header} >= {self.module_completion_threshold}"
            samples_with_mod_df = modules_df.query(query_string)
            if samples_with_mod_df.shape[0] == 0:
                continue
            # if we are working with module data from metagenomes, we may have multiple complete copies of the module in
            # the same sample. We drop these duplicates before proceeding.
            duplicates = samples_with_mod_df.index.duplicated()
            samples_with_mod_df = samples_with_mod_df[~duplicates]

            # we need to explicitly ignore samples without a group here, because they were taken out of sample_groups_df
            # and if only ungrouped samples end up having this module, we will get an index error
            samples_with_mod_list = list(samples_with_mod_df.index)
            for s in samples_with_none_group:
                if s in samples_with_mod_list:
                    samples_with_mod_list.remove(s)
            if len(samples_with_mod_list) == 0:
                continue

            mod_name = samples_with_mod_df['module_name'][0]
            output_dict[mod_name] = {}
            output_dict[mod_name]['MODULE'] = mod_name
            output_dict[mod_name]['accession'] = mod_num
            output_dict[mod_name]['sample_ids'] = ','.join(samples_with_mod_list)
            sample_group_subset = sample_groups_df.loc[samples_with_mod_list]
            p_values = sample_group_subset['group'].value_counts()

            # we need the categories p and N values to be in the same order for finding associated groups
            p_vector = np.array([])
            N_vector = np.array([])
            for c in group_list:
                if c not in p_values.index:
                    p_values[c] = 0
                p_vector = np.append(p_vector, p_values[c]/N_values[c])
                N_vector = np.append(N_vector, N_values[c])

            # compute associated groups for functional enrichment
            enriched_groups_vector = utils.get_enriched_groups(p_vector, N_vector)

            associated_groups = [c for i,c in enumerate(group_list) if enriched_groups_vector[i]]
            output_dict[mod_name]['associated_groups'] = ','.join(associated_groups)

            for c in group_list:
                output_dict[mod_name]["p_%s" % c] = p_values[c]/N_values[c]
                output_dict[mod_name]["N_%s" % c] = N_values[c]

        utils.store_dict_as_TAB_delimited_file(output_dict, output_file_path, key_header='accession', headers=header_list)


    def run_enrichment_stats(self):
        """This function is the driver for running the enrichment script on the modules data."""

        self.progress.new('Enrichment analysis')

        self.progress.update('Converting modules mode output into input for enrichment script')
        enrichment_input_path = filesnpaths.get_temp_file_path()

        if anvio.DEBUG:
            self.progress.reset()
            self.run.info("Temporary input file for enrichment script", enrichment_input_path)

        self.get_enrichment_input(enrichment_input_path)

        self.progress.end()

        # run the enrichment analysis
        enrichment_stats = utils.run_functional_enrichment_stats(enrichment_input_path,
                                                                 self.output_file_path,
                                                                 run=self.run,
                                                                 progress=self.progress)

        return enrichment_stats


## STATIC FUNCTIONS
def module_definition_to_enzyme_accessions(mod_definition):
    """Parses a module definition string into a list of enzyme accessions."""

    # anything that is not (),-+ should be converted to spaces, then we can split on the spaces to get the accessions
    mod_definition = re.sub('[\(\)\+\-,]', ' ', mod_definition).strip()
    acc_list = re.split(r'\s+', mod_definition)

    return acc_list


def _download_worker(
    input_queue: mp.Queue,
    output_queue: mp.Queue,
    max_num_tries: int = 100,
    wait_secs: float = 10.0
) -> None:
    """
    Multiprocessing worker to download files from a queue.

    Parameters
    ==========
    input_queue : multiprocessing.Queue
        Queue of length-two iterables of the URL and local path for each file to download.

    output_queue : multiprocessing.Queue
        Queue in which the success of each download operation is recorded, with True put in the
        output queue if the download succeeded and the local path from the input queue put in the
        output queue if the download failed (after exceeding the maximum number of tries).

    max_num_tries : int, 100
        The maximum number of times to try downloading a file (in case of a connection reset).

    wait_secs : float, 10.0
        The number of seconds to wait between each file download attempt.

    Returns
    =======
    None
    """
    while True:
        url, path = input_queue.get()
        num_tries = 0
        while True:
            try:
                utils.download_file(url, path)
                output = True
                break
            except (ConfigError, ConnectionResetError) as e:
                num_tries += 1
                if num_tries > max_num_tries:
                    output = path
                    break
                time.sleep(wait_secs)
        output_queue.put(output)


def download_org_pathway_image_files(
    pathway_name: str,
    data_dir: str,
    kegg_rest_api_get: str = 'http://rest.kegg.jp/get'
) -> Tuple[str, str]:
    """
    Download an organism-specific pathway map and associated KGML file.
    
    Parameters
    ==========
    pathway_name : str
        This ID has 2 parts: the first 3 org characters are specific to the organism, such as 'eco'
        for E. coli, and the last 5 digits identify the pathway, such as '00010'.
    
    data_dir : str
        Path to KEGG data directory set up by anvi'o with the necessary subdirectory structure.
    
    kegg_rest_api_get : str, 'http://rest.kegg.jp/get'
        KEGG API URL for downloading files.
    
    Returns
    =======
    Tuple[str, str]
        Pathway PNG image and KGML XML filepaths of downloaded files.
    """
    png_url = f'{kegg_rest_api_get}/{pathway_name}/image'
    kgml_url = f'{kegg_rest_api_get}/{pathway_name}/kgml'
    
    png_path = os.path.join(data_dir, 'png', '1x', 'org', f'{pathway_name}.png')
    kgml_path = os.path.join(data_dir, 'kgml', '1x', 'org', f'{pathway_name}.xml')
    
    utils.download_file(png_url, png_path)
    utils.download_file(kgml_url, kgml_path)
    
    return (png_path, kgml_path)

def _download_pathway_image_files_worker(
    input_queue: mp.Queue,
    output_queue: mp.Queue,
    max_num_tries: int = 100,
    wait_secs: float = 10.0
) -> None:
    """
    Multiprocessing worker to download pathway maps and associated KGML files given a pathway ID.

    Parameters
    ==========
    input_queue : multiprocessing.Queue
        Queue of input data stored in dictionaries formatted as follows, with values being strings.
        {
            'pathway_number': <last 5-digit part of the pathway ID>,
            'url_stem': <URL stem for KEGG downloads>,
            'data_dir': <KEGG map data directory with proper subdirectory structure>
        }
        Here is a description of the required subdirectory structure of the data directory. It must
        contain subdirectories 'png' and 'kgml', within each of which are subdirectories '1x' and
        '2x'. Within 'png/1x' are 5 directories, 'map', 'ko', 'ec', 'rn', and 'org'. Within 'png/2x'
        is one directory, 'map'. Within 'kgml/1x' and 'kgml/2x' are 4 directories, 'ko', 'ec', 'rn',
        and 'org'.

    output_queue : multiprocessing.Queue
        Queue of output data stored in dictionaries formatted as follows, with values being length-2
        lists of 1) the target download filepath and 2) an integer indicating what happened with the
        download. A value of 0 indicates that there was no attempt at downloading the file because
        the program did not need to try, e.g., for non-global maps, 'ko', 'ec', and 'rn' map images
        are not downloaded; also, if there was a connection error in trying to download a KGML file,
        then the associated map image files did not need to be downloaded. A value of 1 indicates
        that the file downloaded successfully. A value of 2 indicates that the file was unavailable
        for download, e.g., there is no KGML RN file available for the pathway. A value of 3
        indicates that there was a connection error preventing download. A value of 4 indicates that
        there was no attempt to download because the program found that other requisite files were
        unavailable, e.g., a map image is not downloaded if it has no reference KGML files
        associated with it.
        {
            'png_1x_map': [<filepath>, <integer>],
            'png_2x_map': [<filepath>, <integer>],
            'png_1x_ko': [<filepath>, <integer>],
            'png_1x_ec': [<filepath>, <integer>],
            'png_1x_rn': [<filepath>, <integer>],
            'kgml_ko': [<filepath>, <integer>],
            'kgml_ec': [<filepath>, <integer>],
            'kgml_rn': [<filepath>, <integer>]
        }

    max_num_tries : int, 10
        The maximum number of times to try downloading a file (in case of a connection reset).

    wait_secs : float, 10.0
        The number of seconds to wait between each file download attempt.

    Returns
    =======
    None
    """
    while True:
        input = input_queue.get()
        pathway_number: str = input['pathway_number']
        url: str = input['url_stem']
        data_dir: str = input['data_dir']
        
        png_1x_map_url = f'{url}/map{pathway_number}/image'
        png_2x_map_url = f'{url}/map{pathway_number}/image2x'
        png_1x_ko_url = f'{url}/ko{pathway_number}/image'
        png_1x_ec_url = f'{url}/ec{pathway_number}/image'
        png_1x_rn_url = f'{url}/rn{pathway_number}/image'
        kgml_ko_url = f'{url}/ko{pathway_number}/kgml'
        kgml_ec_url = f'{url}/ec{pathway_number}/kgml'
        kgml_rn_url = f'{url}/rn{pathway_number}/kgml'
        
        png_1x_map_path = os.path.join(data_dir, 'png', '1x', 'map', f'map{pathway_number}.png')
        png_2x_map_path = os.path.join(data_dir, 'png', '2x', 'map', f'map{pathway_number}.png')
        png_1x_ko_path = os.path.join(data_dir, 'png', '1x', 'ko', f'ko{pathway_number}.png')
        png_1x_ec_path = os.path.join(data_dir, 'png', '1x', 'ec', f'ec{pathway_number}.png')
        png_1x_rn_path = os.path.join(data_dir, 'png', '1x', 'rn', f'rn{pathway_number}.png')
        kgml_ko_path = os.path.join(data_dir, 'kgml', '1x', 'ko', f'ko{pathway_number}.xml')
        kgml_ec_path = os.path.join(data_dir, 'kgml', '1x', 'ec', f'ec{pathway_number}.xml')
        kgml_rn_path = os.path.join(data_dir, 'kgml', '1x', 'rn', f'rn{pathway_number}.xml')
        
        output: Dict[str, List[str, int]] = {
            'png_1x_map': [png_1x_map_path, 0],
            'png_2x_map': [png_2x_map_path, 0],
            'png_1x_ko': [png_1x_ko_path, 0],
            'png_1x_ec': [png_1x_ec_path, 0],
            'png_1x_rn': [png_1x_rn_path, 0],
            'kgml_ko': [kgml_ko_path, 0],
            'kgml_ec': [kgml_ec_path, 0],
            'kgml_rn': [kgml_rn_path, 0]
        }
        
        if re.match(GLOBAL_MAP_ID_PATTERN, pathway_number):
            is_global_map = True
        else:
            is_global_map = False
        
        # First try to download KGML files for the pathway. Map images are only downloaded if there
        # is at least 1 KGML file associated with it.
        max_tries_exceeded = False
        for key, kgml_url, kgml_path in (
            ('kgml_ko', kgml_ko_url, kgml_ko_path),
            ('kgml_ec', kgml_ec_url, kgml_ec_path),
            ('kgml_rn', kgml_rn_url, kgml_rn_path)
        ):
            num_tries = 0
            while True:
                try:
                    utils.download_file(kgml_url, kgml_path)
                    output[key][1] = 1
                    break
                except ConnectionResetError:
                    num_tries += 1
                    if num_tries > max_num_tries:
                        max_tries_exceeded = True
                        output[key][1] = 3
                        break
                    time.sleep(wait_secs)
                except ConfigError as e:
                    if 'HTTP Error 404' in str(e):
                        output[key][1] = 2
                        break
                    else:
                        num_tries += 1
                        if num_tries > max_num_tries:
                            max_tries_exceeded = True
                            output[key][1] = 3
                            break
                        time.sleep(wait_secs)
                        
        if max_tries_exceeded:
            # Connection errors prevented at least 1 of the KO, EC, or RN KGML files from being
            # downloaded, so it remains unknown if these files are actually available for the
            # pathway map.
            output_queue.put(output)
            continue
        elif output['kgml_ko'][1] == 2 and output['kgml_ec'][1] == 2 and output['kgml_rn'][1] == 2:
            # No KO, EC, and RN KGML files are available for the pathway map. For instance, this is
            # the case for drug maps with KEGG IDs starting with 'map07', such as 'map07011',
            # 'Penicillins'.
            output['png_1x_map'][1] = 4
            output['png_2x_map'][1] = 4
            if is_global_map:
                output['png_1x_ko'][1] = 4
                output['png_1x_ec'][1] = 4
                output['png_1x_rn'][1] = 4
            output_queue.put(output)
            continue
        
        dl_items = [
            ('png_1x_map', png_1x_map_url, png_1x_map_path),
            ('png_2x_map', png_2x_map_url, png_2x_map_path)
        ]
        if is_global_map:
            if output['kgml_ko'][1] == 1:
                dl_items.append(('png_1x_ko', png_1x_ko_url, png_1x_ko_path))
            elif output['kgml_ko'][1] == 2:
                output['png_1x_ko'][1] = 4
                
            if output['kgml_ec'][1] == 1:
                dl_items.append(('png_1x_ec', png_1x_ec_url, png_1x_ec_path))
            elif output['kgml_ec'][1] == 2:
                output['png_1x_ec'][1] = 4
                
            if output['kgml_rn'][1] == 1:
                dl_items.append(('png_1x_rn', png_1x_rn_url, png_1x_rn_path))
            elif output['kgml_rn'][1] == 2:
                output['png_1x_rn'][1] = 4
        for key, image_url, image_path in dl_items:
            num_tries = 0
            while True:
                try:
                    utils.download_file(image_url, image_path)
                    output[key][1] = 1
                    break
                except ConnectionResetError:
                    num_tries += 1
                    if num_tries > max_num_tries:
                        output[key][1] = 3
                        break
                    time.sleep(wait_secs)
                except ConfigError as e:
                    if 'HTTP Error 404' in str(e):
                        output[key][1] = 2
                        break
                    else:
                        num_tries += 1
                        if num_tries > max_num_tries:
                            output[key][1] = 3
                            break
                        time.sleep(wait_secs)
        output_queue.put(output)
    