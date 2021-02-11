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
import pandas as pd
import numpy as np
import scipy.sparse as sps
from scipy import stats

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.tables as t
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions
from anvio.dbops import ContigsSuperclass, ContigsDatabase, ProfileDatabase
from anvio.genomedescriptions import MetagenomeDescriptions, GenomeDescriptions


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

"""Some critical constants for metabolism estimation output formatting."""
# dict containing possible output modes
# output_suffix should be unique to a mode so that multiple output modes can be used at once
# data_dict indicates which data dictionary is used for generating the output (modules or kofams)
# headers list describes which information to include in the output file; see OUTPUT_HEADERS dict below for more info
# description is what is printed when --list-available-modes parameter is used
OUTPUT_MODES = {'kofam_hits_in_modules': {
                    'output_suffix': "kofam_hits_in_modules.txt",
                    'data_dict': "modules",
                    'headers': ["unique_id", "kegg_module", "module_is_complete",
                                "module_completeness", "path_id", "path", "path_completeness",
                                "kofam_hit", "gene_caller_id", "contig"],
                    'description': "Information on each KOfam hit that belongs to a KEGG module"
                    },
                'modules': {
                    'output_suffix': "modules.txt",
                    'data_dict': "modules",
                    'headers': ["unique_id", "kegg_module", "module_name", "module_class", "module_category",
                                "module_subcategory", "module_definition", "module_completeness", "module_is_complete",
                                "kofam_hits_in_module", "gene_caller_ids_in_module", "warnings"],
                    'description': "Information on KEGG modules"
                    },
                'modules_custom': {
                    'output_suffix': "modules_custom.txt",
                    'data_dict': "modules",
                    'headers': None,
                    'description': "A custom tab-delimited output file where you choose the included KEGG modules data using --custom-output-headers"
                    },
                'kofam_hits': {
                    'output_suffix': "kofam_hits.txt",
                    'data_dict': "kofams",
                    'headers': ["unique_id", "ko", "gene_caller_id", "contig", "modules_with_ko", "ko_definition"],
                    'description': "Information on all KOfam hits in the contigs DB, regardless of KEGG module membership"
                    },
                }
# dict containing matrix headers of information that we can output in custom mode
# key corresponds to the header's key in output dictionary (returned from generate_output_dict_for_modules() function)
# cdict_key is the header's key in modules or kofams data dictionary (if any)
# mode_type indicates which category of output modes (modules or kofams) this header can be used for. If both, this is 'all'
# description is printed when --list-available-output-headers parameter is used
OUTPUT_HEADERS = {'unique_id' : {
                        'cdict_key': None,
                        'mode_type': 'all',
                        'description': "Just an integer that keeps our data organized. No real meaning here. Always included in output, so no need to specify it on the command line"
                        },
                  'kegg_module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "KEGG module number"
                        },
                  'module_is_complete' : {
                        'cdict_key': 'complete',
                        'mode_type': 'modules',
                        'description': "Whether a KEGG module is considered complete or not based on its percent completeness and the completeness threshold"
                        },
                  'module_completeness' : {
                        'cdict_key': 'percent_complete',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a KEGG module"
                        },
                  'module_name' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "English name/description of a KEGG module"
                        },
                  'module_class' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism class of a KEGG module"
                        },
                  'module_category' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism category of a KEGG module"
                        },
                  'module_subcategory' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism subcategory of a KEGG module"
                        },
                  'module_definition' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "KEGG-formatted definition of a KEGG module. Describes the metabolic pathway "
                                       "in terms of the KOS that belong to the module"
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
                        'description': "Comma-separated list of compounds that are intermediates the metabolic pathway "
                                       "(compounds that are both outputs and inputs of reaction(s) in the pathway)"
                        },
                  'gene_caller_ids_in_module': {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of gene caller IDs of KOfam hits in a module"
                        },
                  'gene_caller_id': {
                        'cdict_key': None,
                        'mode_type': 'all',
                        'description': "Gene caller ID of a single KOfam hit in the contigs DB. If you choose this header, each "
                                       "line in the output file will be a KOfam hit"
                        },
                  'kofam_hits_in_module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of KOfam hits in a module"
                        },
                  'kofam_hit' : {
                        'cdict_key': 'kofam_hits',
                        'mode_type': 'modules',
                        'description': "KO number of a single KOfam hit. If you choose this header, each line in the output file "
                                       "will be a KOfam hit"
                        },
                  'contig' : {
                        'cdict_key': 'genes_to_contigs',
                        'mode_type': 'all',
                        'description': "Contig that a KOfam hit is found on. If you choose this header, each line in the output "
                                       "file will be a KOfam hit"
                        },
                  'path_id' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Integer ID for a path through a KEGG module. No real meaning and just for data organization. "
                                       "If you choose this header, each line in the output file will be a KOfam hit"
                        },
                  'path' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A path through a KEGG module (a linear sequence of KOs that together represent each metabolic step "
                                       "in the module. Most modules have several of these due to KO redundancy). If you choose this header, "
                                       "each line in the output file will be a KOfam hit"
                        },
                  'path_completeness' : {
                        'cdict_key': 'pathway_completeness',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a particular path through a KEGG module. If you choose this header, each line "
                                       "in the output file will be a KOfam hit"
                        },
                  'warnings' : {
                        'cdict_key': 'warnings',
                        'mode_type': 'modules',
                        'description': "If we are missing a KOfam profile for one of the KOs in a module, there will be a note in this column. "
                        },
                  'ko' : {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': 'KEGG Orthology (KO) number of a KOfam hit'
                        },
                  'modules_with_ko': {
                        'cdict_key': 'modules',
                        'mode_type': 'kofams',
                        'description': 'A comma-separated list of modules that the KO belongs to'
                        },
                  'ko_definition': {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': 'The functional annotation associated with the KO number'
                        },
                  }

# global metadata header lists for matrix format
# if you want to add something here, don't forget to add it to the dictionary in the corresponding
# get_XXX_metadata_dictionary() function
MODULE_METADATA_HEADERS = ["module_name", "module_class", "module_category", "module_subcategory"]
KO_METADATA_HEADERS = ["ko_definition", "modules_with_ko"]



class KeggContext(object):
    """The purpose of this base class is to define shared functions and file paths for all KEGG operations."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # default data directory will be called KEGG and will store the KEGG Module data as well
        self.default_kegg_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')
        self.kegg_data_dir = A('kegg_data_dir') or self.default_kegg_dir
        self.orphan_data_dir = os.path.join(self.kegg_data_dir, "orphan_data")
        self.module_data_dir = os.path.join(self.kegg_data_dir, "modules")
        self.hmm_data_dir = os.path.join(self.kegg_data_dir, "HMMs")
        self.pathway_data_dir = os.path.join(self.kegg_data_dir, "pathways")
        self.quiet = A('quiet') or False
        self.just_do_it = A('just_do_it')

        # shared variables for all KEGG subclasses
        self.kofam_hmm_file_path = os.path.join(self.hmm_data_dir, "Kofam.hmm") # file containing concatenated KOfam hmms
        self.ko_list_file_path = os.path.join(self.kegg_data_dir, "ko_list.txt")
        self.kegg_module_file = os.path.join(self.kegg_data_dir, "modules.keg")
        self.kegg_pathway_file = os.path.join(self.kegg_data_dir, "pathways.keg")
        self.kegg_modules_db_path = os.path.join(self.kegg_data_dir, "MODULES.db")

        # sanity check to prevent automatic overwriting of non-default kegg data dir
        if A('reset') and A('kegg_data_dir'):
            raise ConfigError("You are attempting to run KEGG setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --kegg-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with KEGG data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.kegg_data_dir, self.kegg_data_dir))


    def setup_ko_dict(self):
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
        """

        self.ko_dict = utils.get_TAB_delimited_file_as_dictionary(self.ko_list_file_path)
        self.ko_skip_list, self.ko_no_threshold_list = self.get_ko_skip_list()

        # if we are currently setting up KEGG, we should generate a text file with the ko_list entries
        # of the KOs that have no scoring threshold
        if self.__class__.__name__ in ['KeggSetup']:
            orphan_ko_dict = {ko:self.ko_dict[ko] for ko in self.ko_skip_list}
            orphan_ko_dict.update({ko:self.ko_dict[ko] for ko in self.ko_no_threshold_list})

            if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
                raise ConfigError("Hmm. Something is out of order. The orphan data directory %s does not exist "
                                  "yet, but it needs to in order for the setup_ko_dict() function to work." % self.orphan_data_dir)
            orphan_ko_path = os.path.join(self.orphan_data_dir, "01_ko_fams_with_no_threshold.txt")
            orphan_ko_headers = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos", "definition"]
            utils.store_dict_as_TAB_delimited_file(orphan_ko_dict, orphan_ko_path, key_header="knum", headers=orphan_ko_headers)

        [self.ko_dict.pop(ko) for ko in self.ko_skip_list]
        [self.ko_dict.pop(ko) for ko in self.ko_no_threshold_list]


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


    def get_module_metadata_dictionary(self, mnum):
        """Returns a dictionary of metadata for the given module.

        The dictionary must include all the metadata from MODULE_METADATA_HEADERS,
        using those headers as keys.
        """

        mnum_class_dict = self.kegg_modules_db.get_kegg_module_class_dict(mnum)

        metadata_dict = {}
        metadata_dict["module_name"] = self.kegg_modules_db.get_module_name(mnum)
        metadata_dict["module_class"] = mnum_class_dict["class"]
        metadata_dict["module_category"] = mnum_class_dict["category"]
        metadata_dict["module_subcategory"] = mnum_class_dict["subcategory"]
        return metadata_dict


    def get_ko_metadata_dictionary(self, knum):
        """Returns a dictionary of metadata for the given KO.

        The dictionary must include all the metadata from KO_METADATA_HEADERS,
        using those headers as keys.
        """

        mod_list = self.kegg_modules_db.get_modules_for_knum(knum)
        if mod_list:
            mod_list_str = ",".join(mod_list)
        else:
            mod_list_str = "None"

        if knum not in self.ko_dict:
            raise ConfigError("Something is mysteriously wrong. Your contigs database "
                              f"has an annotation for KO {knum} but this KO is not in "
                              "the KO dictionary. This should never have happened.")

        metadata_dict = {}
        metadata_dict["ko_definition"] = self.ko_dict[knum]['definition']
        metadata_dict["modules_with_ko"] = mod_list_str
        return metadata_dict


class KeggSetup(KeggContext):
    """Class for setting up KEGG Kofam HMM profiles and modules.

    It performs sanity checks and downloads, unpacks, and prepares the profiles for later use by `hmmscan`.
    It also downloads module files and creates the MODULES.db.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-setup-kegg-kofams
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories when testing this class
    """

    def __init__(self, args, run=run, progress=progress, skip_init=False):
        self.args = args
        self.run = run
        self.progress = progress
        self.kegg_archive_path = args.kegg_archive
        self.download_from_kegg = True if args.download_from_kegg else False
        self.kegg_snapshot = args.kegg_snapshot

        if self.kegg_archive_path and self.download_from_kegg:
            raise ConfigError("You provided two incompatible input options, --kegg-archive and --download-from-kegg. "
                              "Please pick either just one or none of these. ")
        if self.kegg_snapshot and self.download_from_kegg or self.kegg_snapshot and self.kegg_archive_path:
            raise ConfigError("You cannot request setup from an anvi'o KEGG snapshot at the same time as from KEGG directly or from one of your "
                              "KEGG archives. Please pick just one setup option and try again.")

        # initializing this to None here so that it doesn't break things downstream
        self.pathway_dict = None

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists('hmmpress')

        # this is to avoid a strange os.path.dirname() bug that returns nothing if the input doesn't look like a path
        if '/' not in self.kegg_data_dir:
            self.kegg_data_dir += '/'
        filesnpaths.is_output_dir_writable(os.path.dirname(self.kegg_data_dir))

        if not args.reset and not anvio.DEBUG and not skip_init:
            self.is_database_exists()

        if self.download_from_kegg and not self.kegg_archive_path and not skip_init:
            filesnpaths.gen_output_directory(self.kegg_data_dir, delete_if_exists=args.reset)
            filesnpaths.gen_output_directory(self.hmm_data_dir, delete_if_exists=args.reset)
            filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)
            filesnpaths.gen_output_directory(self.module_data_dir, delete_if_exists=args.reset)
            filesnpaths.gen_output_directory(self.pathway_data_dir, delete_if_exists=args.reset)

        # get KEGG snapshot info for default setup
        self.target_snapshot = self.kegg_snapshot or 'v2020-12-23'
        self.target_snapshot_yaml = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG-SNAPSHOTS.yaml')
        self.snapshot_dict = utils.get_yaml_as_dict(self.target_snapshot_yaml)

        if self.target_snapshot not in self.snapshot_dict.keys():
            self.run.warning(None, header="AVAILABLE KEGG SNAPSHOTS", lc="yellow")
            available_snapshots = sorted(list(self.snapshot_dict.keys()))
            for snapshot_name in available_snapshots:
                self.run.info_single(snapshot_name + (' (latest)' if snapshot_name == available_snapshots[-1] else ''))

            raise ConfigError("Whoops. The KEGG snapshot you requested is not one that is known to anvi'o. Please try again, and "
                              f"this time pick from the list shown above.")

        # default download path for KEGG snapshot
        self.default_kegg_data_url = self.snapshot_dict[self.target_snapshot]['url']
        self.default_kegg_archive_file = self.snapshot_dict[self.target_snapshot]['archive_name']

        # download from KEGG option: ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        # dictionary mapping downloaded file name to final decompressed file name or folder location
        self.files = {'ko_list.gz': self.ko_list_file_path, 'profiles.tar.gz': self.kegg_data_dir}

        # download from KEGG option: module/pathway map htext files and API link
        self.kegg_module_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
        self.kegg_pathway_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir="
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"


    def is_database_exists(self):
        """This function determines whether the user has already downloaded the Kofam HMM profiles and KEGG modules."""

        if os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("It seems you already have KOfam HMM profiles installed in '%s', please use the --reset flag "
                              "or delete this directory manually if you want to re-download it." % self.kegg_data_dir)

        if os.path.exists(self.kegg_module_file):
            raise ConfigError("Interestingly, though KOfam HMM profiles are not installed on your system, KEGG module "
                              "information seems to have been already downloaded in %s. Please use the --reset flag or "
                              "delete this directory manually to let this script re-download everything from scratch."
                              % self.kegg_data_dir)

        if os.path.exists(self.kegg_pathway_file):
            raise ConfigError("Interestingly, though KOfam HMM profiles are not installed on your system, KEGG pathway "
                              "information seems to have been already downloaded in %s. Please use the --reset flag or "
                              "delete this directory manually to let this script re-download everything from scratch."
                              % self.kegg_data_dir)

        if os.path.exists(self.module_data_dir):
            raise ConfigError("It seems the KEGG module directory %s already exists on your system. This is even more "
                              "strange because Kofam HMM profiles have not been downloaded. We suggest you to use the "
                              "--reset flag or delete the KEGG directory (%s) manually to download everything from scratch."
                              % (self.module_data_dir, self.kegg_data_dir))

        if os.path.exists(self.pathway_data_dir):
            raise ConfigError("It seems the KEGG pathway directory %s already exists on your system. This is even more "
                              "strange because Kofam HMM profiles have not been downloaded. We suggest you to use the "
                              "--reset flag or delete the KEGG directory (%s) manually to download everything from scratch."
                              % (self.pathway_data_dir, self.kegg_data_dir))


    def download_profiles(self):
        """This function downloads the Kofam profiles."""

        self.run.info("Kofam Profile Database URL", self.database_url)

        try:
            for file_name in self.files.keys():
                utils.download_file(self.database_url + '/' + file_name,
                    os.path.join(self.kegg_data_dir, file_name), progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download KEGG KOfam profiles from the KEGG website. Something "
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

        f = open(self.kegg_module_file, 'rU')
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

        f = open(self.kegg_pathway_file, 'rU')
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


    def download_modules(self):
        """This function downloads the KEGG modules.

        To do so, it also processes the KEGG module file into a dictionary via the process_module_file() function.
        To verify that each file has been downloaded properly, we check that the last line is '///'.
        """

        self.run.info("KEGG Module Database URL", self.kegg_rest_api_get)

        # download the kegg module file, which lists all modules
        try:
            utils.download_file(self.kegg_module_download_path, self.kegg_module_file, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download the KEGG Module htext file from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")

        # get module dict
        self.process_module_file()
        self.run.info("Number of KEGG Modules", len(self.module_dict.keys()))

        # download all modules
        for mnum in self.module_dict.keys():
            file_path = os.path.join(self.module_data_dir, mnum)
            utils.download_file(self.kegg_rest_api_get + '/' + mnum,
                file_path, progress=self.progress, run=self.run)
            # verify entire file has been downloaded
            f = open(file_path, 'rU')
            f.seek(0, os.SEEK_END)
            f.seek(f.tell() - 4, os.SEEK_SET)
            last_line = f.readline().strip('\n')
            if not last_line == '///':
                raise ConfigError("The KEGG module file %s was not downloaded properly. We were expecting the last line in the file "
                                  "to be '///', but instead it was %s. Formatting of these files may have changed on the KEGG website. "
                                  "Please contact the developers to see if this is a fixable issue. If it isn't, we may be able to "
                                  "provide you with a legacy KEGG data archive that you can use to setup KEGG with the --kegg-archive flag."
                                  % (file_path, last_line))


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
            f = open(file_path, 'rU')
            f.seek(0, os.SEEK_END)
            f.seek(f.tell() - 4, os.SEEK_SET)
            last_line = f.readline().strip('\n')
            if not last_line == '///':
                raise ConfigError("The KEGG pathway file %s was not downloaded properly. We were expecting the last line in the file "
                                  "to be '///', but instead it was %s. Formatting of these files may have changed on the KEGG website. "
                                  "Please contact the developers to see if this is a fixable issue. If it isn't, we may be able to "
                                  "provide you with a legacy KEGG data archive that you can use to setup KEGG with the --kegg-archive flag."
                                  % (file_path, last_line))


    def decompress_files(self):
        """This function decompresses the Kofam profiles."""

        self.progress.new('Decompressing files')
        for file_name in self.files.keys():
            self.progress.update('Decompressing file %s' % file_name)
            full_path = os.path.join(self.kegg_data_dir, file_name)

            if full_path.endswith("tar.gz"):
                utils.tar_extract_file(full_path, output_file_path=self.files[file_name], keep_original=False)
            else:
                utils.gzip_decompress_file(full_path, output_file_path=self.files[file_name], keep_original=False)

            self.progress.update("File decompressed. Yay.")
        self.progress.end()


    def confirm_downloaded_profiles(self):
        """This function verifies that all Kofam profiles have been properly downloaded.

        It is intended to be run after the files have been decompressed. The profiles directory should contain hmm files from K00001.hmm to
        K23763.hmm with some exceptions; all KO numbers from ko_list file (except those in ko_skip_list) should be included.
        """

        ko_nums = self.ko_dict.keys()
        for k in ko_nums:
            if k not in self.ko_skip_list:
                hmm_path = os.path.join(self.kegg_data_dir, "profiles/%s.hmm" % k)
                if not os.path.exists(hmm_path):
                    raise ConfigError("The KOfam HMM profile at %s does not exist. This probably means that something went wrong "
                                      "while downloading the KOfam database. Please run `anvi-setup-kegg-kofams` with the --reset "
                                      "flag. If that still doesn't work, please contact the developers to see if the issue is fixable. "
                                      "If it isn't, we may be able to provide you with a legacy KEGG data archive that you can use to "
                                      "setup KEGG with the --kegg-archive flag." % (hmm_path))


    def move_orphan_files(self):
        """This function moves the following to the orphan files directory:

            - profiles that do not have ko_list entries
            - profiles whose ko_list entries have no scoring threshold (in ko_no_threshold_list)

        And, the following profiles should not have been downloaded, but if they were then we move them, too:
            - profiles whose ko_list entries have no data at all (in ko_skip_list)
        """

        if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
            raise ConfigError("Hmm. Something is out of order. The orphan data directory %s does not exist "
                              "yet, but it needs to in order for the move_orphan_files() function to work." % self.orphan_data_dir)

        no_kofam_path = os.path.join(self.orphan_data_dir, "00_hmm_profiles_with_no_ko_fams.hmm")
        no_kofam_file_list = []
        no_threshold_path = os.path.join(self.orphan_data_dir, "02_hmm_profiles_with_ko_fams_with_no_threshold.hmm")
        no_threshold_file_list = []
        no_data_path = os.path.join(self.orphan_data_dir, "03_hmm_profiles_with_ko_fams_with_no_data.hmm")
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

        # now we concatenate the orphan KO hmms into the orphan data directory
        if no_kofam_file_list:
            utils.concatenate_files(no_kofam_path, no_kofam_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning("Please note that while anvi'o was building your databases, she found %d "
                             "HMM profiles that did not have any matching KOfam entries. We have removed those HMM "
                             "profiles from the final database. You can find them under the directory '%s'."
                             % (len(no_kofam_file_list), self.orphan_data_dir))

        if no_threshold_file_list:
            utils.concatenate_files(no_threshold_path, no_threshold_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning("Please note that while anvi'o was building your databases, she found %d "
                             "KOfam entries that did not have any threshold to remove weak hits. We have removed those HMM "
                             "profiles from the final database. You can find them under the directory '%s'."
                             % (len(no_threshold_file_list), self.orphan_data_dir))

        if no_data_file_list:
            utils.concatenate_files(no_data_path, no_data_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning("Please note that while anvi'o was building your databases, she found %d "
                             "HMM profiles that did not have any associated data (besides an annotation) in their KOfam entries. "
                             "We have removed those HMM profiles from the final database. You can find them under the directory '%s'."
                             % (len(no_data_file_list), self.orphan_data_dir))


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

        # there is no reason to keep the original HMM profiles around, unless we are debugging
        if not anvio.DEBUG:
            shutil.rmtree((os.path.join(self.kegg_data_dir, "profiles")))

        self.progress.update('Running hmmpress...')
        cmd_line = ['hmmpress', self.kofam_hmm_file_path]
        log_file_path = os.path.join(self.hmm_data_dir, '00_hmmpress_log.txt')
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Hmm. There was an error while running `hmmpress` on the Kofam HMM profiles. "
                              "Check out the log file ('%s') to see what went wrong." % (log_file_path))
        else:
            # getting rid of the log file because hmmpress was successful
            os.remove(log_file_path)

        self.progress.end()


    def setup_modules_db(self):
        """This function creates the Modules DB from the KEGG Module files."""

        try:
            mod_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args, module_dictionary=self.module_dict, pathway_dictionary=self.pathway_dict, run=run, progress=progress)
            mod_db.create()
        except Exception as e:
            print(e)
            raise ConfigError("While attempting to build the MODULES.db, anvi'o encountered an error, which should be printed above. "
                              "If you look at that error and it seems like something you cannot handle, please contact the developers "
                              "for assistance, as it may be possible that a recent update to KEGG has broken our setup process. If that "
                              "is the case, we may be able to provide you with a legacy KEGG data archive that you can use to set up "
                              "KEGG with the --kegg-archive flag, while we try to fix things. :) ")


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
        """

        is_ok = True

        # check top-level files and folders
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        expected_directories_and_files = [self.orphan_data_dir,
                                          self.module_data_dir,
                                          self.hmm_data_dir,
                                          #self.pathway_data_dir,   #TODO: uncomment me when we start incorporating pathways
                                          self.ko_list_file_path,
                                          self.kegg_module_file,
                                          #self.kegg_pathway_file,  #TODO: uncomment me when we start incorporating pathways
                                          self.kegg_modules_db_path]
        for f in expected_directories_and_files:
            path_to_f_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(f))
            if not os.path.exists(path_to_f_in_archive):
                is_ok = False
                if anvio.DEBUG:
                    self.run.warning("The KEGG archive does not contain the folllowing expected file or directory: %s"
                                     % (path_to_f_in_archive))

        # check hmm files
        path_to_hmms_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(self.hmm_data_dir))
        kofam_hmm_basename = os.path.basename(self.kofam_hmm_file_path)
        expected_hmm_files = [kofam_hmm_basename]
        for h in expected_hmm_files:
            path_to_h_in_archive = os.path.join(path_to_hmms_in_archive, h)
            if not os.path.exists(path_to_h_in_archive):
                is_ok = False
                if anvio.DEBUG:
                    self.run.warning("The KEGG archive does not contain the folllowing expected hmm file: %s"
                                     % (path_to_h_in_archive))
            expected_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
            for ext in expected_extensions:
                path_to_expected_hmmpress_file = path_to_h_in_archive + ext
                if not os.path.exists(path_to_expected_hmmpress_file):
                    is_ok = False
                    if anvio.DEBUG:
                        self.run.warning("The KEGG archive does not contain the folllowing expected `hmmpress` output: %s"
                                         % (path_to_expected_hmmpress_file))

        return is_ok


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
                             "{current_db_version}). You may want to run `anvi-migrate` on this database before you do anything else. "
                             f"Here is the path to the database: {self.kegg_modules_db_path}")


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
        self.progress.end()
        if archive_is_ok:
            if os.path.exists(self.kegg_data_dir) and self.kegg_data_dir != self.default_kegg_dir:
                raise ConfigError("You are attempting to set up KEGG from a KEGG data archive in a non-default data directory (%s) which already exists. "
                                  "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to get rid of "
                                  "directories that have been specified with --kegg-data-dir. If you really want to get rid of this "
                                  "directory and replace it with the KEGG archive data, then please remove the directory yourself using "
                                  "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                                  "the safest way to handle things." % (self.kegg_data_dir, self.kegg_data_dir))
            elif os.path.exists(self.kegg_data_dir):
                shutil.rmtree(self.kegg_data_dir)
            path_to_kegg_in_archive = os.path.join(unpacked_archive_name, "KEGG")
            shutil.move(path_to_kegg_in_archive, self.kegg_data_dir)
            shutil.rmtree(unpacked_archive_name)

            # if necessary, warn user about migrating the modules db
            self.check_modules_db_version()

        else:
            debug_output = "We kept the unpacked archive for you to take a look at it. It is at %s and you may want " \
                           "to delete it after you are done checking its contents." % os.path.abspath(unpacked_archive_name)
            if not anvio.DEBUG:
                shutil.rmtree(unpacked_archive_name)
                debug_output = "The unpacked archive has been deleted, but you can re-run the script with the --debug " \
                               "flag to keep it if you want to see its contents."
            else:
                self.run.warning("The unpacked archive file %s was kept for debugging purposes. You may want to "
                                 "clean it up after you are done looking through it." % (os.path.abspath(unpacked_archive_name)))
            raise ConfigError("The provided archive file %s does not appear to be a KEGG data directory, so anvi'o is unable "
                              "to use it. %s" % (self.kegg_archive_path, debug_output))


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


    def setup_data(self):
        """This is a driver function which executes the KEGG setup process.

        It downloads, decompresses, and hmmpresses the KOfam profiles.
        It also downloads and processes the KEGG Module files into the MODULES.db.
        """

        if self.kegg_archive_path:
            self.setup_from_archive()
        elif self.download_from_kegg:
            # mostly for developers and the adventurous
            self.download_profiles()
            self.decompress_files()
            self.download_modules()
            #self.download_pathways()   # This is commented out because we do not do anything with pathways downstream, but we will in the future.
            self.setup_ko_dict()
            self.run_hmmpress()
            self.setup_modules_db()
        else:
            # the default, set up from frozen KEGG release
            self.setup_kegg_snapshot()


class KeggRunHMMs(KeggContext):
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
        self.keep_all_hits = True if A('keep_all_hits') else False
        self.log_bitscores = True if A('log_bitscores') else False
        self.ko_dict = None # should be set up by setup_ko_dict()

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists(self.hmm_program)

        # verify that Kofam HMM profiles have been set up
        if not os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("Anvi'o is unable to find the Kofam.hmm file at %s. This can happen one of two ways. Either you "
                              "didn't specify the correct KEGG data directory using the flag --kegg-data-dir, or you haven't "
                              "yet set up the Kofam data by running `anvi-setup-kegg-kofams`. Hopefully you now know what to do "
                              "to fix this problem. :) " % self.hmm_data_dir)

        utils.is_contigs_db(self.contigs_db_path)

        self.setup_ko_dict() # read the ko_list file into self.ko_dict

        # load existing kegg modules db
        self.kegg_modules_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args)

        # reminder to be a good citizen
        self.run.warning("Anvi'o will annotate your database with the KEGG KOfam database, as described in "
                         "Aramaki et al (doi:10.1093/bioinformatics/btz859) When you publish your findings, "
                         "please do not forget to properly credit this work.", lc='green', header="CITATION")


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

        contigs_db = ContigsDatabase(self.contigs_db_path)
        contigs_db.db.set_meta_value('modules_db_hash', self.kegg_modules_db.db.get_meta_value('hash'))
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

        if not knum in self.ko_dict:
            if ok_if_missing_from_dict:
                return "Unknown function with KO num %s" % knum
            else:
                raise ConfigError("It seems %s found a KO number that does not exist "
                                  "in the KOfam ko_list file: %s" % (self.hmm_program, knum))

        return self.ko_dict[knum]['definition']


    def parse_kofam_hits(self, hits_dict):
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
        """

        total_num_hits = len(hits_dict.values())
        self.progress.new("Parsing KOfam hits", progress_total_items=total_num_hits)
        self.functions_dict = {}
        self.kegg_module_names_dict = {}
        self.kegg_module_classes_dict = {}
        self.gcids_to_hits_dict = {}
        self.gcids_to_functions_dict = {}
        counter = 0
        num_hits_removed = 0
        for hit_key,hmm_hit in hits_dict.items():
            knum = hmm_hit['gene_name']
            gcid = hmm_hit['gene_callers_id']
            keep = False

            self.progress.update("Removing weak hits for %s [%d of %d]" % (knum, self.progress.progress_current_item + 1, total_num_hits))

            # later, we will need to quickly access the hits for each gene call. So we map gcids to the keys in the raw hits dictionary
            if gcid not in self.gcids_to_hits_dict:
                self.gcids_to_hits_dict[gcid] = [hit_key]
            else:
                self.gcids_to_hits_dict[gcid].append(hit_key)

            if knum not in self.ko_dict:
                self.progress.reset()
                raise ConfigError("Something went wrong while parsing the KOfam HMM hits. It seems that KO "
                                  f"{knum} is not in the noise cutoff dictionary for KOs. That means we do "
                                  "not know how to distinguish strong hits from weak ones for this KO. "
                                  "Anvi'o will fail now :( Please contact a developer about this error to "
                                  "get this mess fixed. ")
            # if hit is above the bitscore threshold, we will keep it
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

                counter += 1
            else:
                num_hits_removed += 1

            self.progress.increment()
            self.progress.reset()

        self.progress.end()
        self.run.info("Number of weak hits removed by KOfam parser", num_hits_removed)
        self.run.info("Number of hits remaining in annotation dict ", len(self.functions_dict.keys()))


    def update_dict_for_genes_with_missing_annotations(self, gcids_list, hits_dict):
        """This function adds functional annotations for genes with missing hits to the dictionary.

        The reason this is necessary is that the bitscore thresholds can be too stringent, causing
        us to miss legitimate annotations. To find these annotations, we adopt the following heuristic:
            For every gene without a KOfam annotation, we examine all the hits with an e-value below X
            and a bitscore above Y% of the threshold. If those hits are all to a unique KO profile,
            then we annotate the gene call with that KO.

        For reasons that are hopefully obvious, this function must be called after parse_kofam_hits(),
        which establishes the self.functions_dict attribute.

        PARAMETERS
        ===========
        gcids_list : list
            The list of gene caller ids in the contigs database. We will use this to figure out which
            genes have no annotations
        hits_dict : dictionary
            The output from the hmmsearch parser, which should contain all hits (ie, weak hits not yet removed)
        """

        # for each gene call, check for annotation in self.functions_dict
        for gcid in gcids_list:
            if gcid not in self.gcids_to_functions_dict:
                # if no annotation, get all hits for gene caller id from hits_dict

        # get set of hits that fit parameters
        # if unique KO, add annotation with best e-value to self.functions_dict


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

        if not hmm_hits_file:
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

        # parse hmmscan output
        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        # add functions and KEGG modules info to database
        self.parse_kofam_hits(search_results_dict)
        self.update_dict_for_genes_with_missing_annotations(all_gcids_in_contigs_db, search_results_dict)
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
        self.output_file_prefix = A('output_file_prefix') or "kegg-metabolism"
        self.write_dict_to_json = True if A('get_raw_data_as_json') else False
        self.json_output_file_path = A('get_raw_data_as_json')
        self.store_json_without_estimation = True if A('store_json_without_estimation') else False
        self.estimate_from_json = A('estimate_from_json') or None
        self.output_modes = A('kegg_output_modes') or A('output_modes') or "modules"
        self.custom_output_headers = A('custom_output_headers') or None
        self.matrix_format = True if A('matrix_format') else False
        self.matrix_include_metadata = True if A('include_metadata') else False
        self.exclude_zero_modules = False if A('include_zeros') else True
        self.only_complete = True if A('only_complete') else False
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

            if "unique_id" not in self.custom_output_headers:
                self.custom_output_headers = ["unique_id"] + self.custom_output_headers
            elif self.custom_output_headers.index("unique_id") != 0:
                self.custom_output_headers.remove("unique_id")
                self.custom_output_headers = ["unique_id"] + self.custom_output_headers
            self.available_modes['modules_custom']['headers'] = self.custom_output_headers

        # parse specific matrix modules if necessary
        if self.module_specific_matrices:
            self.module_specific_matrices = self.module_specific_matrices.split(",")


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
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file = A('bin_ids_file')
        self.contigs_db_project_name = "Unknown"
        self.database_name = A('database_name')
        self.multi_mode = True if A('multi_mode') else False

        KeggEstimatorArgs.__init__(self, self.args)

        self.name_header = None
        if self.metagenome_mode:
            self.name_header = 'contig_name'
        elif self.profile_db_path and not self.metagenome_mode:
            self.name_header = 'bin_name'
        elif not self.profile_db_path and not self.metagenome_mode:
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
                                        'description': "Name of genome/bin/metagenome in which we find KOfam hits and/or KEGG modules"
                                        }

        # input options sanity checks
        if not self.estimate_from_json and not self.contigs_db_path:
            raise ConfigError("NO INPUT PROVIDED. You must provide (at least) a contigs database or genomes file to this program, unless you are using the --estimate-from-json "
                              "flag, in which case you must provide a JSON-formatted file.")

        self.bin_ids_to_process = None
        if self.bin_id and self.bin_ids_file:
            raise ConfigError("You have provided anvi'o with both the individual bin id %s and a file with bin ids (%s). "
                              "Please make up your mind. Which one do you want an estimate for? :)" % (self.bin_id, self.bin_ids_file))
        elif self.bin_id:
            self.bin_ids_to_process = [self.bin_id]
        elif self.bin_ids_file:
            filesnpaths.is_file_exists(self.bin_ids_file)
            self.bin_ids_to_process = [line.strip() for line in open(self.bin_ids_file).readlines()]


        if (self.bin_id or self.bin_ids_file or self.collection_name) and not self.profile_db_path:
            raise ConfigError("You have requested metabolism estimation for a bin or set of bins, but you haven't provided "
                              "a profiles database. Unfortunately, this just does not work. Please try again.")

        if not self.metagenome_mode and self.profile_db_path and not self.collection_name:
            raise ConfigError("If you provide a profiles DB, you should also provide a collection name.")

        if self.store_json_without_estimation and not self.json_output_file_path:
            raise ConfigError("Whoops. You seem to want to store the metabolism dictionary in a JSON file, but you haven't provided the name of that file. "
                              "Please use the --get-raw-data-as-json flag to do so.")
        if self.store_json_without_estimation and self.estimate_from_json:
            raise ConfigError("It is impossible to both estimate metabolism from JSON data and produce a JSON file without estimation at the same time... "
                              "anvi'o is judging you SO hard right now.")

        if self.profile_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        # output options sanity checks
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


        # init the base class
        KeggContext.__init__(self, self.args)

        # init the KO dictionary
        self.setup_ko_dict()

        if not self.estimate_from_json:
            utils.is_contigs_db(self.contigs_db_path)

        # load existing kegg modules db
        if not os.path.exists(self.kegg_modules_db_path):
            raise ConfigError("It appears that a modules database (%s) does not exist in the KEGG data directory %s. "
                              "Perhaps you need to specify a different KEGG directory using --kegg-data-dir. Or perhaps you didn't run "
                              "`anvi-setup-kegg-kofams`, though we are not sure how you got to this point in that case "
                              "since you also cannot run `anvi-run-kegg-kofams` without first having run KEGG setup. But fine. Hopefully "
                              "you now know what you need to do to make this message go away." % ("MODULES.db", self.kegg_data_dir))
        kegg_modules_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args, quiet=self.quiet)

        if not self.estimate_from_json:
            # here we load the contigs DB just for sanity check purposes.
            # We will need to load it again later just before accessing data to avoid SQLite error that comes from different processes accessing the DB
            contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
            self.contigs_db_project_name = contigs_db.meta['project_name']

            # sanity check that contigs db was annotated with same version of MODULES.db that will be used for metabolism estimation
            if 'modules_db_hash' not in contigs_db.meta:
                raise ConfigError("Based on the contigs DB metadata, the contigs DB that you are working with has not been annotated with hits to the "
                                  "KOfam database, so there are no KOs to estimate metabolism from. Please run `anvi-run-kegg-kofams` on this contigs DB "
                                  "before you attempt to run this script again.")
            contigs_db_mod_hash = contigs_db.meta['modules_db_hash']
            mod_db_hash = kegg_modules_db.db.get_meta_value('hash')
            if contigs_db_mod_hash != mod_db_hash:
                raise ConfigError("The contigs DB that you are working with has been annotated with a different version of the MODULES.db than you are working with now. "
                                  "Perhaps you updated your KEGG setup after running `anvi-run-kegg-kofams` on this contigs DB? Or maybe you have multiple KEGG data "
                                  "directories set up on your computer, and the one you are using now is different from the one that you used for `anvi-run-kegg-kofams`? "
                                  "Well. The solution to the first problem is to re-run `anvi-run-kegg-kofams` on the contigs DB (%s) using the updated MODULES.db "
                                  "(located in the KEGG data directory %s). The solution to the second problem is to specify the appropriate KEGG data directory using "
                                  "the --kegg-data-dir flag. If neither of those things make this work, then you should contact the developers to see if they can help you "
                                  "figure this out. For those who need this information, the Modules DB used to annotate this contigs database previously had the "
                                  "following hash: %s. And the hash of the current Modules DB is: %s" % (self.contigs_db_path, self.kegg_data_dir, contigs_db_mod_hash, mod_db_hash))
            contigs_db.disconnect()
        kegg_modules_db.disconnect()

        if not self.quiet:
            self.run.warning("Anvi'o will reconstruct metabolism for modules in the KEGG MODULE database, as described in "
                             "Muto et al (doi:10.1021/ci3005379). When you publish your findings, "
                             "please do not forget to properly credit this work.", lc='green', header="CITATION")


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

    def init_hits_and_splits(self):
        """This function loads KOfam hits, gene calls, splits, and contigs from the contigs DB.

        We will need the hits with their KO numbers (accessions) so that we can go through the MODULES.db and determine
        which steps are present in each module. And we will need the other information so that we can determine which hits belong
        to which genomes/bins when we are handling multiple of these, and for help in computing redundancy.
        This function gets this info as a list of tuples (one tuple per kofam hit), and it makes sure that these lists don't include
        hits that we shouldn't be considering.

        RETURNS
        =======
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering
        """

        self.progress.new('Loading data from Contigs DB')
        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        genes_in_splits = contigs_db.db.get_some_columns_from_table(t.genes_in_splits_table_name, "gene_callers_id, split")
        genes_in_contigs = contigs_db.db.get_some_columns_from_table(t.genes_in_contigs_table_name, "gene_callers_id, contig")
        kofam_hits = contigs_db.db.get_some_columns_from_table(t.gene_function_calls_table_name, "gene_callers_id, accession",
                                                               where_clause="source = 'KOfam'")
        min_contig_length_in_contigs_db = contigs_db.db.get_max_value_in_column(t.contigs_info_table_name, "length", return_min_instead=True)
        contigs_db.disconnect()

        # get rid of gene calls that are not associated with KOfam hits.
        all_gene_calls_in_splits = set([tpl[0] for tpl in genes_in_splits])
        gene_calls_with_kofam_hits = set([tpl[0] for tpl in kofam_hits])
        gene_calls_without_kofam_hits = all_gene_calls_in_splits.difference(gene_calls_with_kofam_hits)

        if gene_calls_without_kofam_hits:
            self.progress.update("Removing %s gene calls without KOfam hits" % len(gene_calls_without_kofam_hits))
            genes_in_splits = [tpl for tpl in genes_in_splits if tpl[0] not in gene_calls_without_kofam_hits]
            genes_in_contigs = [tpl for tpl in genes_in_contigs if tpl[0] not in gene_calls_without_kofam_hits]
            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning("The following gene calls in your contigs DB were removed from consideration as they "
                                 "do not have any hits to the KOfam database: %s" % (gene_calls_without_kofam_hits))


        # get rid of splits and contigs (and their associated gene calls) that are not in the profile DB
        if self.profile_db_path:
            # if we were given a blank profile, we will assume we want all splits and pull all splits from the contigs DB
            if utils.is_blank_profile(self.profile_db_path):
                self.run.warning("You seem to have provided a blank profile. No worries, we can still estimate metabolism for you. "
                                 "But we cannot load splits from the profile DB, so instead we are assuming that you are interested in "
                                 "ALL splits and we will load those from the contigs database.")
                split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.contigs_db_path))
            else:
                split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.profile_db_path))
            split_names_in_contigs_db = set([tpl[1] for tpl in genes_in_splits])
            splits_missing_in_profile_db = split_names_in_contigs_db.difference(split_names_in_profile_db)

            min_contig_length_in_profile_db = ProfileDatabase(self.profile_db_path).meta['min_contig_length']

            if len(splits_missing_in_profile_db):
                self.progress.reset()
                self.run.warning("Please note that anvi'o found %s splits in your contigs database with KOfam hits. But only %s of them "
                                 "appear in the profile database. As a result, anvi'o will now remove the %s splits with KOfam hits "
                                 "that occur only in the contigs db from all downstream analyses. Where is this difference coming from though? "
                                 "Well. This is often the case because the 'minimum contig length parameter' set during the `anvi-profile` "
                                 "step can exclude many contigs from downstream analyses (often for good reasons, too). For "
                                 "instance, in your case the minimum contig length goes as low as %s nts in your contigs database. "
                                 "Yet, the minimum contig length set in the profile databaes is %s nts. Hence the difference. Anvi'o "
                                 "hopes that this explaines some things." % (pp(len(split_names_in_contigs_db)),
                                                                             pp(len(split_names_in_profile_db)),
                                                                             pp(len(splits_missing_in_profile_db)),
                                                                             pp(min_contig_length_in_contigs_db),
                                                                             pp(min_contig_length_in_profile_db)))

                self.progress.update("Removing %s splits (and associated gene calls) that were missing from the profile db" % pp(len(splits_missing_in_profile_db)))
                genes_in_splits = [tpl for tpl in genes_in_splits if tpl[1] not in splits_missing_in_profile_db]
                remaining_gene_calls = [tpl[0] for tpl in genes_in_splits]
                genes_in_contigs = [tpl for tpl in genes_in_contigs if tpl[0] in remaining_gene_calls]
                kofam_hits = [tpl for tpl in kofam_hits if tpl[0] in remaining_gene_calls]

        # combine the information for each gene call into neat tuples for returning
        # each gene call is only on one split of one contig, so we can convert these lists of tuples into dictionaries for easy access
        # but some gene calls have multiple kofam hits (and some kofams have multiple gene calls), so we must keep the tuple structure for those
        self.progress.update("Organizing KOfam hit data")
        gene_calls_splits_dict = {tpl[0] : tpl[1] for tpl in genes_in_splits}
        gene_calls_contigs_dict = {tpl[0] : tpl[1] for tpl in genes_in_contigs}
        assert len(gene_calls_splits_dict.keys()) == len(genes_in_contigs)

        kofam_gene_split_contig = []
        for gene_call_id, ko in kofam_hits:
            kofam_gene_split_contig.append((ko, gene_call_id, gene_calls_splits_dict[gene_call_id], gene_calls_contigs_dict[gene_call_id]))

        self.progress.update("Done")
        self.progress.end()

        self.run.info("Contigs DB", self.contigs_db_path, quiet=self.quiet)
        self.run.info("KOfam hits", "%d found" % len(kofam_hits), quiet=self.quiet)
        self.run.info("Profile DB", self.profile_db_path, quiet=self.quiet)
        self.run.info('Metagenome mode', self.metagenome_mode)
        if self.collection_name:
            self.run.info('Collection', self.collection_name)
        if self.bin_id:
            self.run.info('Bin ID', self.bin_id)
        elif self.bin_ids_file:
            self.run.info('Bin IDs file', self.bin_ids_file)

        if not self.quiet and not len(kofam_hits):
            self.run.warning("Hmmm. No KOfam hits were found in this contigs DB, so all metabolism estimate outputs will be empty. This is fine, and "
                             "could even be biologically correct. But we thought we'd mention it just in case you thought it was weird. "
                             "Other, technical reasons that this could have happened include: 1) you didn't annotate with `anvi-run-kegg-kofams` "
                             "and 2) you imported KEGG functional annotations but the 'source' was not 'KOfam'.")

        return kofam_gene_split_contig


    def init_paths_for_modules(self):
        """This function unrolls the module DEFINITION for each module and places it in an attribute variable for
        all downstream functions to access.

        It unrolls the module definition into a list of all possible paths, where each path is a list of atomic steps.
        Atomic steps include singular KOs, protein complexes, modules, non-essential steps, and steps without associated KOs.
        """

        self.module_paths_dict = {}
        modules = self.kegg_modules_db.get_all_modules_as_list()
        for m in modules:
            self.module_paths_dict[m] = self.kegg_modules_db.unroll_module_definition(m)


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
            self.run.info("Number of splits", len(split_list))

        # initialize all modules with empty lists and dicts for kos, gene calls
        modules = self.kegg_modules_db.get_all_modules_as_list()
        all_kos = self.kegg_modules_db.get_all_knums_as_list()
        for mnum in modules:
            bin_level_module_dict[mnum] = {"gene_caller_ids" : set(),
                                           "kofam_hits" : {},
                                           "genes_to_contigs" : {},
                                           "contigs_to_genes" : {},
                                           "warnings" : []
                                          }
        for knum in all_kos:
            if knum not in self.ko_dict:
                mods_it_is_in = self.kegg_modules_db.get_modules_for_knum(knum)
                if mods_it_is_in:
                    if anvio.DEBUG:
                        mods_str = ", ".join(mods_it_is_in)
                        self.run.warning(f"Oh dear. We do not appear to have a KOfam profile for {knum}. This means "
                                        "that any modules this KO belongs to can never be fully complete (this includes "
                                        f"{mods_str}). ")
                    for m in mods_it_is_in:
                        bin_level_module_dict[m]["warnings"].append(f"No KOfam profile for {knum}")
                continue

            bin_level_ko_dict[knum] = {"gene_caller_ids" : set(),
                                     "modules" : None,
                                     "genes_to_contigs" : {},
                                     "contigs_to_genes" : {}
                                     }

        kos_not_in_modules = []
        for ko, gene_call_id, split, contig in kofam_hits_in_splits:
            present_in_mods = self.kegg_modules_db.get_modules_for_knum(ko)
            if not present_in_mods:
                kos_not_in_modules.append(ko)
                # KOs that are not in modules will not be initialized above in the ko hit dictionary, so we add them here if we haven't already
                if ko not in bin_level_ko_dict:
                    bin_level_ko_dict[ko] = {"gene_caller_ids" : set(),
                                             "modules" : None,
                                             "genes_to_contigs" : {},
                                             "contigs_to_genes" : {}
                                             }
            else:
                bin_level_ko_dict[ko]["modules"] = present_in_mods
            for m in present_in_mods:
                bin_level_module_dict[m]["gene_caller_ids"].add(gene_call_id)
                if ko in bin_level_module_dict[m]["kofam_hits"]:
                    bin_level_module_dict[m]["kofam_hits"][ko].append(gene_call_id)
                else:
                    bin_level_module_dict[m]["kofam_hits"][ko] = [gene_call_id]
                bin_level_module_dict[m]["genes_to_contigs"][gene_call_id] = contig
                if contig in bin_level_module_dict[m]["contigs_to_genes"]:
                    bin_level_module_dict[m]["contigs_to_genes"][contig].add(gene_call_id)
                else:
                    bin_level_module_dict[m]["contigs_to_genes"][contig] = set([gene_call_id])

            bin_level_ko_dict[ko]["gene_caller_ids"].add(gene_call_id)
            bin_level_ko_dict[ko]["genes_to_contigs"][gene_call_id] = contig
            if contig in bin_level_ko_dict[ko]["contigs_to_genes"]:
                bin_level_ko_dict[ko]["contigs_to_genes"][contig].add(gene_call_id)
            else:
                bin_level_ko_dict[ko]["contigs_to_genes"][contig] = set([gene_call_id])

        if anvio.DEBUG:
            self.run.info("KOs processed", "%d in bin" % len(kofam_hits_in_splits))
            if kos_not_in_modules:
                self.run.warning("Just so you know, the following KOfam hits did not belong to any KEGG modules in the MODULES.db: %s"
                % ", ".join(kos_not_in_modules))

        return bin_level_module_dict, bin_level_ko_dict


    def compute_module_completeness_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the completeness of the specified module within the given bin metabolism dictionary.

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
        "percent_complete"              the completeness of the module, which is the maximum pathway completeness
        "complete"                      whether the module completeness falls over the completeness threshold

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

        for p in self.module_paths_dict[mnum]:
            num_complete_steps_in_path = 0
            num_nonessential_steps_in_path = 0 # so that we don't count nonessential steps when computing completeness
            for atomic_step in p:
                # there are 5 types of atomic steps to take care of
                # 1) regular old single KOs, ie Kxxxxx
                if atomic_step[0] == "K" and len(atomic_step) == 6:
                    if atomic_step in present_list_for_mnum:
                        num_complete_steps_in_path += 1
                # 2) protein complexes, ie Kxxxxx+Kyyyyy-Kzzzzz (2 types of complex components - essential and nonessential)
                elif atomic_step[0] == "K" and (atomic_step[6] == "+" or atomic_step[6] == "-"):
                    idx = 6
                    essential_components = [atomic_step[0:idx]]
                    while idx < len(atomic_step):
                        component_ko = atomic_step[idx+1:idx+7]
                        if atomic_step[idx] == "+":
                            essential_components.append(component_ko)
                        else:
                            has_nonessential_step = True
                            if component_ko not in module_nonessential_kos:
                                module_nonessential_kos.append(component_ko)
                        idx += 7

                    num_present_components = 0
                    for c in essential_components:
                        if c in present_list_for_mnum:
                            num_present_components += 1
                    component_completeness = num_present_components / len(essential_components)
                    num_complete_steps_in_path += component_completeness
                # 3) non-essential KOs, ie -Kxxxxx
                elif atomic_step[0:2] == "-K" and len(atomic_step) == 7:
                    """
                    OKAY, SO HERE WE HAVE SOME POOPINESS THAT MAY NEED TO BE FIXED EVENTUALLY.
                    Basically, some DEFINITION lines have KOs that seem to be marked non-essential;
                    ie, "-K11024" in "K11023 -K11024 K11025 K11026 K11027".
                    It was difficult to decide whether we should consider only K11024, or K11024 and all following KOs, to be non-essential.
                    For instance, the module M00778 is a complex case that gave us pause - see Fiesta issue 955.
                    But for now, we have decided to just track only the one KO as a 'non-essential step', and to not include such steps in
                    the module completeness estimate.
                    """
                    if atomic_step[1:] not in module_nonessential_kos:
                        module_nonessential_kos.append(atomic_step[1:])
                    num_nonessential_steps_in_path += 1
                    has_nonessential_step = True
                # 4) steps without associated KOs, ie --
                elif atomic_step == "--":
                    # when '--' in a DEFINITION line happens, it signifies a reaction step that has no associated KO.
                    # we assume that such steps are not complete,  because we really can't know if it is from the KOfam hits alone
                    has_no_ko_step = True
                    warning_str = "'--' steps are assumed incomplete"
                    if warning_str not in meta_dict_for_bin[mnum]["warnings"]:
                        meta_dict_for_bin[mnum]["warnings"].append(warning_str)
                # 5) Module numbers, ie Mxxxxx
                elif atomic_step[0] == "M" and len(atomic_step) == 6:
                    """
                    This happens when a module is defined by other modules. For example, photosynthesis module M00611 is defined as
                    (M00161,M00163) M00165 === (photosystem II or photosystem I) and calvin cycle

                    We need all the modules to have been evaluated before we can determine completeness of steps with module numbers.
                    So what we will do here is to use a flag variable to keep track of the modules that have this sort of definition
                    in a list so we can go back and evaluate completeness of steps with module numbers later.
                    """
                    defined_by_modules = True
                else:
                    raise ConfigError("Well. While estimating completeness for module %s, we found an atomic step in the pathway that we "
                                        "are not quite sure what to do with. Here it is: %s" % (mnum, atomic_step))


            path_completeness = num_complete_steps_in_path / (len(p) - num_nonessential_steps_in_path)
            meta_dict_for_bin[mnum]["pathway_completeness"].append(path_completeness)

        # once all paths have been evaluated, we find the path(s) of maximum completeness and set that as the overall module completeness
        # this is not very efficient as it takes two passes over the list but okay
        meta_dict_for_bin[mnum]["percent_complete"] = max(meta_dict_for_bin[mnum]["pathway_completeness"])
        if meta_dict_for_bin[mnum]["percent_complete"] > 0:
            meta_dict_for_bin[mnum]["most_complete_paths"] = [self.module_paths_dict[mnum][i] for i, pc in enumerate(meta_dict_for_bin[mnum]["pathway_completeness"]) if pc == meta_dict_for_bin[mnum]["percent_complete"]]
        else:
            meta_dict_for_bin[mnum]["most_complete_paths"] = []


        if anvio.DEBUG and len(meta_dict_for_bin[mnum]["most_complete_paths"]) > 1:
            self.run.warning("Found %d complete paths for module %s with completeness %s. " % (len(meta_dict_for_bin[mnum]["most_complete_paths"]), mnum, meta_dict_for_bin[mnum]["percent_complete"]),
                            header='DEBUG OUTPUT', lc='yellow')
        over_complete_threshold = True if meta_dict_for_bin[mnum]["percent_complete"] >= self.module_completion_threshold else False
        meta_dict_for_bin[mnum]["complete"] = over_complete_threshold
        meta_dict_for_bin[mnum]["present_nonessential_kos"] = module_nonessential_kos
        if over_complete_threshold:
            meta_dict_for_bin["num_complete_modules"] += 1

        return over_complete_threshold, has_nonessential_step, has_no_ko_step, defined_by_modules


    def adjust_module_completeness_for_bin(self, mod, meta_dict_for_bin):
        """This function adjusts completeness of modules that are defined by other modules.

        This can only be done after all other modules have been evaluated for completeness.
        The function uses similar logic as compute_module_completeness_for_bin() to re-assess whether steps defined
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

            for atomic_step in p:
                # single KOs and protein complexes and '--' steps; were already counted as complete by previous function
                if atomic_step[0] == "K" or atomic_step == "--":
                    num_essential_steps_in_path += 1
                # non-essential KO, don't count as a step in the path
                elif atomic_step[0:2] == "-K" and len(atomic_step) == 7:
                    pass
                # module step; we need to count these based on previously computed module completeness
                elif atomic_step[0] == "M" and len(atomic_step) == 6:
                    num_complete_module_steps += meta_dict_for_bin[atomic_step]["percent_complete"]
                    num_essential_steps_in_path += 1
                else:
                    raise ConfigError("Well. While adjusting completeness estimates for module %s, we found an atomic step in the pathway that we "
                                      "are not quite sure what to do with. Here it is: %s" % (mod, atomic_step))

            # now we adjust the previous pathway completeness
            old_complete_steps_in_path = meta_dict_for_bin[mod]["pathway_completeness"][i] * num_essential_steps_in_path
            adjusted_num_complete_steps_in_path = old_complete_steps_in_path + num_complete_module_steps
            meta_dict_for_bin[mod]["pathway_completeness"][i] = adjusted_num_complete_steps_in_path / num_essential_steps_in_path

        # after adjusting for all paths, adjust overall module completeness
        meta_dict_for_bin[mod]["percent_complete"] = max(meta_dict_for_bin[mod]["pathway_completeness"])
        if meta_dict_for_bin[mod]["percent_complete"] > 0:
            meta_dict_for_bin[mod]["most_complete_paths"] = [self.module_paths_dict[mod][i] for i, pc in enumerate(meta_dict_for_bin[mod]["pathway_completeness"]) if pc == meta_dict_for_bin[mod]["percent_complete"]]
        else:
            meta_dict_for_bin[mod]["most_complete_paths"] = []

        was_already_complete = meta_dict_for_bin[mod]["complete"]
        now_complete = True if meta_dict_for_bin[mod]["percent_complete"] >= self.module_completion_threshold else False
        meta_dict_for_bin[mod]["complete"] = now_complete
        if now_complete and not was_already_complete:
            meta_dict_for_bin["num_complete_modules"] += 1

        return now_complete


    def estimate_for_list_of_splits(self, metabolism_dict_for_list_of_splits, bin_name=None):
        """This is the atomic metabolism estimator function, which builds up the metabolism completeness dictionary for an arbitrary list of splits.

        For example, the list of splits may represent a bin, a single isolate genome, or an entire metagenome.

        The function takes in a metabolism completeness dictionary already initialized with the relevant KOfam hits per module, and updates it
        with the individual steps and completion estimates for each module.

        PARAMETERS
        ==========
        metabolism_dict_for_list_of_splits : dictionary of dictionaries
            the metabolism completeness dictionary of dictionaries for this list of splits. It contains
            one dictionary of module steps and completion information for each module (keyed by module number),
            as well as one key num_complete_modules that tracks the number of complete modules found in these splits.
            Calling functions should assign this dictionary to a metabolism superdict with the bin name as a key.
        bin_name : str
            the name of the bin/genome/metagenome that we are working with
        """

        metabolism_dict_for_list_of_splits["num_complete_modules"] = 0

        complete_mods = []
        mods_def_by_modules = [] # a list of modules that have module numbers in their definitions
        # modules to warn about
        mods_with_unassociated_ko = [] # a list of modules that have "--" steps without an associated KO
        mods_with_nonessential_steps = [] # a list of modules that have nonessential steps like "-K11024"

        # estimate completeness of each module
        for mod in metabolism_dict_for_list_of_splits.keys():
            if mod == "num_complete_modules":
                continue
            mod_is_complete, has_nonessential_step, has_no_ko_step, defined_by_modules \
            = self.compute_module_completeness_for_bin(mod, metabolism_dict_for_list_of_splits)

            if mod_is_complete:
                complete_mods.append(mod)
            if has_nonessential_step:
                mods_with_nonessential_steps.append(mod)
            if has_no_ko_step:
                mods_with_unassociated_ko.append(mod)
            if defined_by_modules:
                mods_def_by_modules.append(mod)

        # go back and adjust completeness of modules that are defined by other modules
        if mods_def_by_modules:
            for mod in mods_def_by_modules:
                mod_is_complete = self.adjust_module_completeness_for_bin(mod, metabolism_dict_for_list_of_splits)

                if mod_is_complete:
                    complete_mods.append(mod)


        # estimate redundancy of each module
        for mod in metabolism_dict_for_list_of_splits.keys():
            if mod == "num_complete_modules":
                continue

            self.compute_module_redundancy_for_bin(mod, metabolism_dict_for_list_of_splits)


        # notify user of the modules that gave some fishy results -- but only for genome mode because it's too wordy otherwise
        if not self.quiet and self.genome_mode:
            if mods_with_nonessential_steps:
                self.run.warning("Please note that anvi'o found one or more non-essential steps in the following KEGG modules: %s.   "
                                 "At this time, we are not counting these steps in our percent completion estimates."
                                 % (", ".join(mods_with_nonessential_steps)))

            if mods_with_unassociated_ko:
                self.run.warning("Just so you know, while estimating the completeness of some KEGG modules, anvi'o saw "
                                 "'--' in the module DEFINITION. This indicates a step in the pathway that has no "
                                 "associated KO. So we really cannot know just based on KOfam hits whether or not this "
                                 "step is present. By default, anvi'o marks these steps incomplete. But they may not be, "
                                 "and as a result their modules may be falsely considered incomplete. So it may be in your "
                                 "interest to go back and take a look at these individual modules to see if you can find the "
                                 "missing enzyme in some other way. Best of luck to you. Here is the list of modules to check out: %s"
                                 % (", ".join(mods_with_unassociated_ko)))

        if anvio.DEBUG or self.genome_mode:
            self.run.info("Bin name", bin_name)
            self.run.info("Module completion threshold", self.module_completion_threshold)
            self.run.info("Number of complete modules", metabolism_dict_for_list_of_splits["num_complete_modules"])
            if complete_mods:
                self.run.info("Complete modules", ", ".join(complete_mods))

        return metabolism_dict_for_list_of_splits

######### REDUNDANCY FUNCTIONS (UNUSED) #########

    def compute_naive_redundancy_for_path(self, num_ko_hits_in_path_dict):
        """This function computes a naive redundancy measure for a module path, given the number of hits per KO in the path.

        naive redundancy = # extra hits / len(path) where a hit is "extra" if it is not the first hit to the KO.
        """

        extra_hits = [num_ko_hits_in_path_dict[ko] - 1 if num_ko_hits_in_path_dict[ko] > 1 else 0 for ko in num_ko_hits_in_path_dict]
        return sum(extra_hits)/len(num_ko_hits_in_path_dict.keys())


    def compute_copywise_redundancy_for_path(self, num_ko_hits_in_path_dict, aggregation_measure="average"):
        """This function computes redundancy based on the completeness of each extra copy of a path.

        The 'base' redundancy score is determined by the number of extra copies with 100% completeness.
        The completeness measurements of all other extra copies are aggregated (using the aggregation_measure) and
        added to this 'base' redundancy to get the overall path redundancy.
        """

        accepted_aggregation_measures = ["average", "median", "weighted_sum", "geometric_mean"]
        extra_hits = [num_ko_hits_in_path_dict[ko] - 1 if num_ko_hits_in_path_dict[ko] > 1 else 0 for ko in num_ko_hits_in_path_dict]
        base_redundancy = min(extra_hits) # number of extra copies of path that are 100% complete
        extra_copy_completeness = []
        # here we get the completeness of every extra copy of the path
        for i in range((base_redundancy+1), max(extra_hits) + 1):
            num_present_kos_in_copy = len([num_hits for num_hits in extra_hits if num_hits >= i])
            extra_copy_completeness.append(num_present_kos_in_copy/len(num_ko_hits_in_path_dict.keys()))

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


    def compute_entropy_weighted_redundancy_for_bin(self, num_ko_hits_in_path_dict):
        """This function computes naive redundancy but weights it by the entropy of the hit distribution."""

        extra_hits = [num_ko_hits_in_path_dict[ko] - 1 if num_ko_hits_in_path_dict[ko] > 1 else 0 for ko in num_ko_hits_in_path_dict]
        total_extra_hits = sum(extra_hits)
        num_kos = len(num_ko_hits_in_path_dict.keys())
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


    def compute_module_redundancy_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the redundancy of the specified module within the given bin metabolism dictionary.

        Each module can have multiple paths, but we only compute redundancy on the paths with the highest completeness
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
        meta_dict_for_bin[mnum]["entropy_weighted"] = []

        paths_of_highest_completeness = meta_dict_for_bin[mnum]["most_complete_paths"]
        if not paths_of_highest_completeness:
            # put zero values in dict wherever necessary
            return

        for p in paths_of_highest_completeness:
            kofam_hits_in_path = { ko : meta_dict_for_bin[mnum]["kofam_hits"][ko] for ko in meta_dict_for_bin[mnum]["kofam_hits"].keys() if ko in p }
            num_hits_per_kofam = { ko : len(kofam_hits_in_path[ko]) for ko in kofam_hits_in_path.keys() }
            for ko in p:
                if ko not in num_hits_per_kofam:
                    num_hits_per_kofam[ko] = 0

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
            meta_dict_for_bin[mnum]["copywise_weighted-sum"].append(cw_gm_redundancy)
            meta_dict_for_bin[mnum]["entropy_weighted"].append(self.compute_entropy_weighted_redundancy_for_bin(num_hits_per_kofam))

        return

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
                                                           self.collection_name), nl_before=1)

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


    def estimate_for_contigs_db_for_metagenome(self, kofam_gene_split_contig):
        """This function handles metabolism estimation for an entire metagenome.

        We treat each contig in the metagenome to be its own 'bin' or 'genome' and estimate
        metabolism separately for each one.

        PARAMETERS
        ==========
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering

        RETURNS
        =======
        metagenome_metabolism_superdict : dictionary of dictionary of dictionaries
            dictionary mapping metagenome name to its metabolism completeness dictionary
        metagenome_ko_superdict : dictionary of dictionary of dictionaries
            dictionary mapping metagenome name to its KOfam hits dictionary
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
                metagenome_metabolism_superdict[contig] = self.estimate_for_list_of_splits(metabolism_dict_for_contig, bin_name=contig)
                single_contig_module_superdict = {contig: metagenome_metabolism_superdict[contig]}
                metagenome_ko_superdict[contig] = ko_dict_for_contig
            else:
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
        new_kegg_metabolism_superdict = {}

        expected_keys_for_module = {"gene_caller_ids", "kofam_hits", "genes_to_contigs", "contigs_to_genes"}
        bins_found = []
        additional_keys = set([])

        self.init_paths_for_modules()

        for bin_name, meta_dict_for_bin in kegg_metabolism_superdict.items():
            bins_found.append(bin_name)
            for mod, mod_dict in meta_dict_for_bin.items():
                if mod == "num_complete_modules":
                    self.run.warning("Your JSON file appears to have been generated from data that already contains metabolic module completeness information. "
                                     "We say this because the key 'num_complete_modules' was found. This isn't a problem; however you should know that anvi'o "
                                     "won't take any of the existing estimation information into account. The only module-level keys that will be used from this file "
                                     "are: %s" % (expected_keys_for_module))
                    continue
                # verify that dict contains the necessary keys for estimation
                if not expected_keys_for_module.issubset(set(mod_dict.keys())):
                    missing_keys = expected_keys_for_module.difference(set(mod_dict.keys()))
                    raise ConfigError("Your JSON file is incorrectly formatted for metabolism estimation. We expect the following keys: %s. "
                                      "However, we didn't find some of them for module %s in %s. Here are the missing keys: %s"
                                      % (expected_keys_for_module, mod, bin_name, missing_keys))

                additional_keys = additional_keys.union(set(mod_dict.keys()).difference(expected_keys_for_module))

                # convert gene_caller_ids and contigs_to_genes lists to sets
                mod_dict['gene_caller_ids'] = set(mod_dict['gene_caller_ids'])
                for contig, gene_list in mod_dict['contigs_to_genes'].items():
                    mod_dict['contigs_to_genes'][contig] = set(gene_list)
                mod_dict['genes_to_contigs'] = {int(g):c for g,c in mod_dict['genes_to_contigs'].items()}

            new_kegg_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(meta_dict_for_bin, bin_name=bin_name)
            single_bin_module_superdict = {bin_name: new_kegg_metabolism_superdict[bin_name]}
            self.append_kegg_metabolism_superdicts(single_bin_module_superdict, ko_superdict_for_list_of_splits={})


        if not self.quiet and additional_keys:
            self.run.warning("Just to let you know, we found the following module-level keys in your JSON file that were totally ignored during metabolism estimation "
                             "(no harm was done by including them): %s" % (additional_keys))

        self.run.info("Bins/genomes/metagenomes found", ", ".join(bins_found))
        return new_kegg_metabolism_superdict


    def estimate_metabolism(self, skip_storing_data=False, output_files_dictionary=None, return_superdicts=False,
                            return_subset_for_matrix_format=False):
        """This is the driver function for estimating metabolism for a single contigs DB.

        It will decide what to do based on whether the input contigs DB is a genome or metagenome.
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

        self.kegg_modules_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args, run=run_quiet, quiet=self.quiet)

        if skip_storing_data or self.write_dict_to_json:
            self.output_file_dict = {}
        else:
            if output_files_dictionary:
                self.output_file_dict = output_files_dictionary
            else:
                self.output_file_dict = self.setup_output_for_appending()

        if self.estimate_from_json:
            kegg_metabolism_superdict = self.estimate_metabolism_from_json_data()
        else:

            kofam_hits_info = self.init_hits_and_splits()
            self.init_paths_for_modules()

            if self.profile_db_path and not self.metagenome_mode:
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_for_bins_in_collection(kofam_hits_info)
            elif not self.profile_db_path and not self.metagenome_mode:
                self.genome_mode = True
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_for_genome(kofam_hits_info)
            elif self.metagenome_mode:
                kegg_metabolism_superdict, kofam_hits_superdict = self.estimate_for_contigs_db_for_metagenome(kofam_hits_info)
            else:
                raise ConfigError("This class doesn't know how to deal with that yet :/")

        if self.write_dict_to_json:
            self.store_metabolism_superdict_as_json(kegg_metabolism_superdict, self.json_output_file_path + ".json")

        self.kegg_modules_db.disconnect()
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
            return self.generate_subsets_for_matrix_format(kegg_metabolism_superdict, kofam_hits_superdict)
        # otherwise we return nothing at all
        return

######### OUTPUT DICTIONARY FUNCTIONS #########

    def generate_output_dict_for_modules(self, kegg_superdict, headers_to_include=None, only_complete_modules=False, exclude_zero_completeness=True):
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
                "percent_complete":              0.66
                "complete":                      False
                                      }

        To distill this information into one line, we need to convert the dictionary on-the-fly to a dict of dicts,
        where each bin-module-path-kofam_hit-gene_caller_id is keyed by an arbitrary integer. There will be a lot of redundant information
        in the rows.

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

        keys_not_in_superdict = set(["unique_id", "genome_name", "bin_name", "contig_name", "kegg_module", "db_name",
                                     "kofam_hits_in_module", "gene_caller_ids_in_module"])
        module_level_headers = set(["module_name", "module_class", "module_category", "module_subcategory", "module_definition",
                                    "module_substrates", "module_products", "module_intermediates"])
        path_and_ko_level_headers = set(["path_id", "path", "path_completeness", "kofam_hit", "gene_caller_id", "contig"])
        remaining_headers = headers_to_include.difference(keys_not_in_superdict)
        remaining_headers = remaining_headers.difference(module_level_headers)
        remaining_headers = remaining_headers.difference(path_and_ko_level_headers)

        # convert to two-level dict where unique id keys for a dictionary of information for each bin/module pair
        d = {}
        if not self.modules_unique_id:
             self.modules_unique_id = 0
        """
        ### FIXME ###
        The unique_id problematic to deal with when we are generating multiple output files because this self.modules_unique_id
        variable is shared between all output files. Since we append to each output file in turn, this means that within an
        output file the unique_id column will have 'jumps' in value (such that a row's unique_id is not always 1 greater than the
        previous row's unique_id). This does not really cause problems and is not an urgent fix, per se, but it looks weird.
        ### FIXME ###
        """

        for bin, mod_dict in kegg_superdict.items():
            for mnum, c_dict in mod_dict.items():
                if mnum == "num_complete_modules":
                    continue

                if anvio.DEBUG:
                    self.run.info("Generating output for module", mnum)

                if only_complete_modules and not c_dict["complete"]:
                    continue
                if exclude_zero_completeness and c_dict["percent_complete"] == 0:
                    continue

                # fetch module info from db
                metadata_dict = self.get_module_metadata_dictionary(mnum)
                module_def = '"' + self.kegg_modules_db.get_kegg_module_definition(mnum) + '"'
                module_substrate_list, module_intermediate_list, module_product_list = self.kegg_modules_db.get_human_readable_compound_lists_for_module(mnum)

                # handle path- and ko-level information
                if headers_to_include.intersection(path_and_ko_level_headers):
                    for p_index in range(len(self.module_paths_dict[mnum])):
                        p = self.module_paths_dict[mnum][p_index]

                        # handle ko-level information
                        for ko in c_dict['kofam_hits']:

                            # some paths include protein complexes, so we must look for KO within these protein complexes as well
                            kos_in_path = set([])
                            for ko_or_complex in p:
                                split_kos = re.split('\+|\-', ko_or_complex)
                                kos_in_path.update(split_kos)
                            if ko not in kos_in_path:
                                continue

                            for gc_id in c_dict["kofam_hits"][ko]:
                                d[self.modules_unique_id] = {}

                                # kofam hit specific info
                                if "kofam_hit" in headers_to_include:
                                    d[self.modules_unique_id]["kofam_hit"] = ko
                                if "gene_caller_id" in headers_to_include:
                                    d[self.modules_unique_id]["gene_caller_id"] = gc_id
                                if "contig" in headers_to_include:
                                    d[self.modules_unique_id]["contig"] = c_dict["genes_to_contigs"][gc_id]

                                # repeated information for each hit
                                # path specific info
                                if "path_id" in headers_to_include:
                                    d[self.modules_unique_id]["path_id"] = p_index
                                if "path" in headers_to_include:
                                    d[self.modules_unique_id]["path"] = ",".join(p)
                                if "path_completeness" in headers_to_include:
                                    d[self.modules_unique_id]["path_completeness"] = c_dict["pathway_completeness"][p_index]

                                # top-level keys and keys not in superdict
                                if self.name_header in headers_to_include:
                                    d[self.modules_unique_id][self.name_header] = bin
                                if "db_name" in headers_to_include:
                                    d[self.modules_unique_id]["db_name"] = self.database_name
                                if "kegg_module" in headers_to_include:
                                    d[self.modules_unique_id]["kegg_module"] = mnum

                                # module specific info
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
                                    if module_substrate_list:
                                        d[self.modules_unique_id]["module_substrates"] = ",".join(module_substrate_list)
                                    else:
                                        d[self.modules_unique_id]["module_substrates"] = "None"
                                if "module_products" in headers_to_include:
                                    if module_product_list:
                                        d[self.modules_unique_id]["module_products"] = ",".join(module_product_list)
                                    else:
                                        d[self.modules_unique_id]["module_products"] = "None"
                                if "module_intermediates" in headers_to_include:
                                    if module_intermediate_list:
                                        d[self.modules_unique_id]["module_intermediates"] = ",".join(module_intermediate_list)
                                    else:
                                        d[self.modules_unique_id]["module_intermediates"] = "None"

                                # comma-separated lists of KOs and gene calls in module
                                if "kofam_hits_in_module" in headers_to_include:
                                    kos_in_mod = c_dict['kofam_hits'].keys()
                                    d[self.modules_unique_id]["kofam_hits_in_module"] = ",".join(kos_in_mod)
                                if "gene_caller_ids_in_module" in headers_to_include:
                                    gcids_in_mod = c_dict['genes_to_contigs'].keys()
                                    gcids_in_mod = [str(x) for x in gcids_in_mod]
                                    d[self.modules_unique_id]["gene_caller_ids_in_module"] = ",".join(gcids_in_mod)

                                # everything else at c_dict level
                                for h in remaining_headers:
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

                                self.modules_unique_id += 1
                else:
                    d[self.modules_unique_id] = {}

                    # top-level keys and keys not in superdict
                    if self.name_header in headers_to_include:
                        d[self.modules_unique_id][self.name_header] = bin
                    if "db_name" in headers_to_include:
                        d[self.modules_unique_id]["db_name"] = self.database_name
                    if "kegg_module" in headers_to_include:
                        d[self.modules_unique_id]["kegg_module"] = mnum

                    # module specific info
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
                        if module_substrate_list:
                            d[self.modules_unique_id]["module_substrates"] = ",".join(module_substrate_list)
                        else:
                            d[self.modules_unique_id]["module_substrates"] = "None"
                    if "module_products" in headers_to_include:
                        if module_product_list:
                            d[self.modules_unique_id]["module_products"] = ",".join(module_product_list)
                        else:
                            d[self.modules_unique_id]["module_products"] = "None"
                    if "module_intermediates" in headers_to_include:
                        if module_intermediate_list:
                            d[self.modules_unique_id]["module_intermediates"] = ",".join(module_intermediate_list)
                        else:
                            d[self.modules_unique_id]["module_intermediates"] = "None"

                    # comma-separated lists of KOs and gene calls in module
                    if "kofam_hits_in_module" in headers_to_include:
                        kos_in_mod = c_dict['kofam_hits'].keys()
                        d[self.modules_unique_id]["kofam_hits_in_module"] = ",".join(kos_in_mod)
                    if "gene_caller_ids_in_module" in headers_to_include:
                        gcids_in_mod = c_dict['genes_to_contigs'].keys()
                        gcids_in_mod = [str(x) for x in gcids_in_mod]
                        d[self.modules_unique_id]["gene_caller_ids_in_module"] = ",".join(gcids_in_mod)

                    # everything else at c_dict level
                    for h in remaining_headers:
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
            headers_to_include = set(OUTPUT_MODES["kofam_hits"]["headers"])
        else:
            headers_to_include = set(headers_to_include)

        d = {}
        if not self.ko_unique_id:
            self.ko_unique_id = 0
        """
        ### FIXME ###
        See note in previous function about the weirdness of the unique_id with multiple output files.
        ### FIXME ###
        """

        for bin, ko_dict in ko_superdict.items():
            for ko, k_dict in ko_dict.items():
                if anvio.DEBUG:
                    self.run.info("Generating output for KO", ko)

                metadata_dict = self.get_ko_metadata_dictionary(ko)

                for gc_id in k_dict["gene_caller_ids"]:
                    d[self.ko_unique_id] = {}

                    if self.name_header in headers_to_include:
                        d[self.ko_unique_id][self.name_header] = bin
                    if "db_name" in headers_to_include:
                        d[self.ko_unique_id]["db_name"] = self.database_name
                    if "ko" in headers_to_include:
                        d[self.ko_unique_id]["ko"] = ko
                    if "gene_caller_id" in headers_to_include:
                        d[self.ko_unique_id]["gene_caller_id"] = gc_id
                    if "contig" in headers_to_include:
                        d[self.ko_unique_id]["contig"] = k_dict["genes_to_contigs"][gc_id]
                    if "modules_with_ko" in headers_to_include:
                        d[self.ko_unique_id]["modules_with_ko"] = metadata_dict["modules_with_ko"]
                    if "ko_definition" in headers_to_include:
                        d[self.ko_unique_id]["ko_definition"] = metadata_dict["ko_definition"]

                    self.ko_unique_id += 1

        return d


    def generate_subsets_for_matrix_format(self, module_superdict, ko_hits_superdict):
        """Here we extract and return three subsets of data from the superdicts, for matrix formatted output.

        The subsets of data that we need are: module completeness scores, module presence/absence, and KO hit frequency.
        Each of these is put into a dictionary (one for modules, one for ko hits) and returned.
        """

        mod_completeness_presence_subdict = {}
        ko_hits_subdict = {}

        for bin, mod_dict in module_superdict.items():
            mod_completeness_presence_subdict[bin] = {}
            for mnum, c_dict in mod_dict.items():
                if mnum == "num_complete_modules":
                    continue
                mod_completeness_presence_subdict[bin][mnum] = {}
                mod_completeness_presence_subdict[bin][mnum]['percent_complete'] = c_dict['percent_complete']
                mod_completeness_presence_subdict[bin][mnum]['complete'] = c_dict['complete']

        for bin, ko_dict in ko_hits_superdict.items():
            ko_hits_subdict[bin] = {}
            for knum, k_dict in ko_dict.items():
                ko_hits_subdict[bin][knum] = {}
                ko_hits_subdict[bin][knum]['num_hits'] = len(k_dict['gene_caller_ids']) # number of hits to this KO in the bin

        return mod_completeness_presence_subdict, ko_hits_subdict

######### OUTPUT GENERATION FUNCTIONS #########

    def append_kegg_metabolism_superdicts(self, module_superdict_for_list_of_splits, ko_superdict_for_list_of_splits):
        """This function appends the metabolism superdicts (for a single genome, bin, or contig in metagenome) to existing files
        for each output mode.

        It appends to the initialized AppendableFile objects in self.output_file_dict.

        This is an alternative to store_kegg_metabolism_superdicts(), which prints the entire metabolism superdicts for all
        genomes/bins/contigs in metagenome at once.
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

            file_obj.append(output_dict, key_header="unique_id", headers=header_list)

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
            utils.store_dict_as_TAB_delimited_file(output_dict, output_path, key_header="unique_id", headers=header_list)
            self.run.info("%s output file" % mode, output_path, nl_before=1)


    def store_metabolism_superdict_as_json(self, kegg_superdict, file_path):
        """This function writes the metabolism superdict into one json file."""

        def set_to_list(obj):
            if isinstance(obj, set):
                return list(obj)

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
        module_data_keys_for_visualization = ['percent_complete']

        metabolism_dict, ko_hit_dict = self.estimate_metabolism(skip_storing_data=True, return_superdicts=True)
        data_for_visualization = {}

        for bin, mod_dict in metabolism_dict.items():
            data_for_visualization[bin] = {}
            for mnum, c_dict in mod_dict.items():
                if mnum == "num_complete_modules":
                    continue

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

        # init the base class for access to shared paths and such
        KeggContext.__init__(self, args)

        KeggEstimatorArgs.__init__(self, self.args)

        if anvio.DEBUG:
            self.run.info("Completeness threshold: multi estimator", self.module_completion_threshold)
            self.run.info("Output Modes", ", ".join(self.output_modes))
            self.run.info("Matrix format", self.matrix_format)
            self.run.info("Matrix will include metadata", self.matrix_include_metadata)
            if self.module_specific_matrices:
                self.run.info("Matrices for specific modules: ", ", ".join(self.module_specific_matrices))

        self.databases = None

        # input sanity checks
        if (self.external_genomes_file and (self.internal_genomes_file or self.metagenomes_file)) \
            or (self.internal_genomes_file and (self.external_genomes_file or self.metagenomes_file)) \
            or (self.metagenomes_file and (self.external_genomes_file or self.internal_genomes_file)):
                raise ConfigError("Multiple file inputs were provided. Please choose only one at a time to make "
                                  "things easier on everybody.")

        if args.estimate_from_json or args.store_json_without_estimation or args.get_raw_data_as_json:
            raise ConfigError("You've provided some JSON parameters. We are sorry to say that these parameters don't "
                              "work for input files with multiple contigs DBs. :( ")

        # output sanity checks
        if self.matrix_format and args.kegg_output_modes:
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

        if not self.quiet:
            self.run.warning("Anvi'o will reconstruct metabolism for modules in the KEGG MODULE database, as described in "
                             "Muto et al (doi:10.1021/ci3005379). When you publish your findings, "
                             "please do not forget to properly credit this work.", lc='green', header="CITATION")


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
                                            'description': "Name of genome in which we find KOfam hits and/or KEGG modules"
                                            }
        elif self.name_header == 'bin_name':
            self.available_headers["bin_name"] = {
                                            'cdict_key': None,
                                            'mode_type': 'all',
                                            'description': "Name of bin in which we find KOfam hits and/or KEGG modules"
                                            }
        elif self.name_header == 'contig_name':
            self.available_headers["contig_name"] = {
                                            'cdict_key': None,
                                            'mode_type': 'all',
                                            'description': "Name of contig (in a metagenome) in which we find KOfam hits and/or KEGG modules"
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

        run.warning(None, header="AVAILABLE OUTPUT HEADERS", lc="green")

        for header, header_meta in self.available_headers.items():
            desc_str = header_meta['description']
            type_str = header_meta['mode_type']
            mode_str = "output modes" if header_meta['mode_type'] == 'all' else "output mode"
            self.run.info(header, f"{desc_str} [{type_str} {mode_str}]")

######### DRIVER ESTIMATION FUNCTIONS -- MULTI #########

    def init_metagenomes(self):
        """This function parses the input metagenomes file and adjusts class attributes as needed"""

        g = MetagenomeDescriptions(self.args, run=self.run, progress=self.progress, enforce_single_profiles=False)
        g.load_metagenome_descriptions()

        # enforce metagenome mode
        if not self.metagenome_mode:
            self.metagenome_mode = True

        self.databases = copy.deepcopy(g.metagenomes)
        self.database_names = copy.deepcopy(g.metagenome_names)


    def init_external_internal_genomes(self):
        """This function parses the input internal/external genomes file and adjusts class attributes as needed"""

        g = GenomeDescriptions(self.args, run=self.run, progress=progress_quiet)
        g.load_genomes_descriptions(skip_functions=True, init=False)

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

        self.update_available_headers_for_multi()

        if anvio.DEBUG:
            self.run.info("Completeness threshold: single estimator", args.module_completion_threshold)
            self.run.info("Database name: single estimator", args.database_name)
            self.run.info("Matrix metadata: single estimator", args.matrix_include_metadata)

        return args


    def get_metabolism_superdict_multi(self):
        """The function that calls metabolism on each individual contigs db.

         If we need matrix format output, it aggregates the results into one dictionary for modules
         and one for KOs, and returns these. (Otherwise, empty dictionaries are returned.)
         """

        metabolism_super_dict = {}
        ko_hits_super_dict = {}

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
                metabolism_super_dict[metagenome_name], ko_hits_super_dict[metagenome_name] = KeggMetabolismEstimator(args, progress=progress_quiet, run=run_quiet).estimate_metabolism(skip_storing_data=True, return_subset_for_matrix_format=True)

            self.progress.increment()
            self.progress.reset()

        self.progress.end()

        if not self.matrix_format:
            for mode, file in files_dict.items():
                file.close()

        return metabolism_super_dict, ko_hits_super_dict


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

        # these will be empty dictionaries unless matrix format
        kegg_metabolism_superdict_multi, ko_hits_superdict_multi = self.get_metabolism_superdict_multi()

        if self.matrix_format:
            self.store_metabolism_superdict_multi_matrix_format(kegg_metabolism_superdict_multi, ko_hits_superdict_multi)

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
                    elif stat_header == 'KO':
                        metadata_dict = self.get_ko_metadata_dictionary(m)
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


    def store_metabolism_superdict_multi_matrix_format(self, module_superdict_multi, ko_superdict_multi):
        """Stores the multi-contigs DB metabolism data in several matrices.

        Contigs DBs are arranged in columns and KEGG modules/KOs are arranged in rows.
        Each module statistic (ie, completeness, presence/absence) will be in a different file.

        The parameters to this function are superdictionaries where each top-level entry is one
        of the multi-mode sample inputs (ie a metagenome, internal, or external genome) and its
        corresponding value comes from running estimate_metabolism() with return_subset_for_matrix_format=True.

        That is:
         module % completeness = module_superdict_multi[sample][bin][mnum]['percent_complete']
         module is complete = module_superdict_multi[sample][bin][mnum]['complete']
         # hits for KO = ko_superdict_multi[sample][bin][knum]['num_hits']

        If self.matrix_include_metadata was True, these superdicts will also include relevant metadata.
        """

        # we need this for metadata
        self.kegg_modules_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args)
        include_zeros = not self.exclude_zero_modules

        # module stats that each will be put in separate matrix file
        # key is the stat, value is the corresponding header in superdict
        module_matrix_stats = {"completeness" : "percent_complete", "presence" : "complete"}
        # all samples/bins have the same modules in the dict so we can pull the item list from the first pair
        first_sample = list(module_superdict_multi.keys())[0]
        first_bin = list(module_superdict_multi[first_sample].keys())[0]
        module_list = list(module_superdict_multi[first_sample][first_bin].keys())
        module_list.sort()

        for stat, key in module_matrix_stats.items():
            self.write_stat_to_matrix(stat_name=stat, stat_header='module', stat_key=key, stat_dict=module_superdict_multi, \
                                      item_list=module_list, stat_metadata_headers=MODULE_METADATA_HEADERS, \
                                      write_rows_with_all_zeros=include_zeros)

        # now we make a KO hit count matrix
        self.setup_ko_dict()
        ko_list = list(self.ko_dict.keys())
        ko_list.sort()
        self.write_stat_to_matrix(stat_name='ko_hits', stat_header='KO', stat_key='num_hits', stat_dict=ko_superdict_multi, \
                                  item_list=ko_list, stat_metadata_headers=KO_METADATA_HEADERS, \
                                  write_rows_with_all_zeros=include_zeros)

        # if necessary, make module specific KO matrices
        if self.module_specific_matrices:
            skipped_mods = []
            for mod in self.module_specific_matrices:
                if mod not in module_list:
                    skipped_mods.append(mod)
                    continue

                kos_in_mod = self.kegg_modules_db.get_kos_from_module_definition(mod)
                mod_big_steps = self.kegg_modules_db.get_top_level_steps_in_module_definition(mod)

                if not self.no_comments:
                    # determine where to place comments containing module steps
                    step_comments = {}
                    lines_with_comment = []
                    for s in mod_big_steps:
                        # what is the first KO in this step?
                        first_k = s.find("K")

                        # we skip making comments on steps without KOs like '--'
                        if first_k < 0:
                            continue

                        # figure out where this KO is in the list
                        first_ko = s[first_k:first_k+6]
                        first_ko_indices = [i for i, x in enumerate(kos_in_mod) if x == first_ko]
                        if not first_ko_indices:
                            raise ConfigError(f"Something went wrong while writing a comment for step '{s}' in the "
                                              f"matrix for {mod}. We couldn't find the first KO, {first_ko}, in the "
                                              "KO list for this module.")

                        # where should we put the step comment?
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

                stat = f"{mod}_ko_hits"
                self.write_stat_to_matrix(stat_name=stat, stat_header="KO", stat_key='num_hits', stat_dict=ko_superdict_multi, \
                                          item_list=kos_in_mod, stat_metadata_headers=KO_METADATA_HEADERS, \
                                          write_rows_with_all_zeros=True, comment_dictionary=step_comments)

            if skipped_mods:
                skipped_list = ", ".join(skipped_mods)
                self.run.warning(f"We couldn't recognize the following module(s): {skipped_list}. So we didn't generate "
                                 "output matrices for them. Maybe you made a typo? Or put an extra comma in somewhere?")

        self.kegg_modules_db.disconnect()


class KeggModulesDatabase(KeggContext):
    """To create or access a Modules DB.

    This DB should be created in the Kegg Data folder during KEGG setup, and will be populated with information from the
    Kegg Module files.
    """

    def __init__(self, db_path, args, module_dictionary=None, pathway_dictionary=None, run=run, progress=progress, quiet=False):
        self.db = None
        self.db_path = db_path
        self.module_dict = module_dictionary
        self.pathway_dict = pathway_dictionary
        self.run = run
        self.progress = progress
        self.quiet = quiet

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

        if os.path.exists(self.db_path):
            utils.is_kegg_modules_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__kegg_modules_version__, new_database=False)

            if not self.quiet:
                self.run.info('Modules database', 'An existing database, %s, has been loaded.' % self.db_path, quiet=self.quiet)
                self.run.info('Kegg Modules', '%d found' % self.db.get_meta_value('num_modules'), quiet=self.quiet)

        else:
            # if self.module_dict is None, then we tried to initialize the DB outside of setup
            if not self.module_dict:
                raise ConfigError("ERROR - a new KeggModulesDatabase() cannot be initialized without providing a modules dictionary. This "
                                  "usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-kofams` may fix this.")
            # This is commented out because we are not yet using pathways. But it should be uncommented when we get to the point of using them :)
            # if not self.pathway_dict:
            #     raise ConfigError("ERROR - a new KeggModulesDatabase() cannot be initialized without providing a pathway dictionary. This "
            #                       "usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-kofams` may fix this.")

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
            # need to sanity check data value field because SOME modules don't follow the 2-space separation formatting
            vals_are_okay, corrected_vals, corrected_def = self.data_vals_sanity_check(data_vals, current_data_name, current_module)

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

        self.touch()

        self.progress.new("Loading %s KEGG modules into Modules DB..." % len(self.module_dict.keys()))

        # sanity check that we setup the modules previously.
        # It shouldn't be a problem since this function should only be called during the setup process after modules download, but just in case.
        if not os.path.exists(self.module_data_dir) or len(self.module_dict.keys()) == 0:
            raise ConfigError("Appparently, the Kegg Modules were not correctly setup and now all sorts of things are broken. The "
                              "Modules DB cannot be created from broken things. BTW, this error is not supposed to happen to anyone "
                              "except maybe developers, so if you do not fall into that category you are likely in deep doo-doo. "
                              "Maybe re-running setup with --reset will work? (if not, you probably should email/Slack/telepathically "
                              "cry out for help to the developers). Anyway, if this helps make things any clearer, the number of modules "
                              "in the module dictionary is currently %s" % len(self.module_dict.keys()))

        # init the Modules table
        mod_table = KeggModulesTable(self.module_table_name)

        # keep track of errors encountered while parsing
        self.parsing_error_dict = {"bad_line_splitting" : [], "bad_kegg_code_format" : []}
        self.num_corrected_errors = 0
        self.num_uncorrected_errors = 0

        num_modules_parsed = 0
        line_number = 0
        for mnum in self.module_dict.keys():
            self.progress.update("Parsing KEGG Module %s" % mnum)
            mod_file_path = os.path.join(self.module_data_dir, mnum)
            f = open(mod_file_path, 'rU')

            prev_data_name_field = None
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

                f.close()

            num_modules_parsed += 1
        # once we are done parsing all modules, we store whatever db entries remain in the db_entries list
        # this is necessary because append_and_store() above only stores every 10000 entries
        self.progress.update("Storing final batch of module entries into DB")
        mod_table.store(self.db)

        self.progress.end()

        # warn user about parsing errors
        if anvio.DEBUG:
            self.run.warning("Several parsing errors were encountered while building the KEGG Modules DB. "
                             "Below you will see which modules threw each type of parsing error. Note that modules which "
                             "threw multiple errors will occur in the list as many times as it threw each error.")
            self.run.info("Bad line splitting (usually due to rogue or missing spaces)", self.parsing_error_dict["bad_line_splitting"])
            self.run.info("Bad KEGG code format (not corrected; possibly problematic)", self.parsing_error_dict["bad_kegg_code_format"])
        else: # less verbose
            self.run.warning("First things first - don't panic. Several parsing errors were encountered while building the KEGG Modules DB. "
                             "But that is probably okay, because if you got to this point it is likely that we already fixed all of them "
                             "ourselves. So don't worry too much. Below you will see how many of each type of error was encountered. If "
                             "you would like to see which modules threw these errors, please re-run the setup using the `--debug` flag (you "
                             "will also probably need the `--reset` flag). When doing so, you will also see which lines caused issues; this "
                             "can be a lot of output, so you can suppress the line-specific output with the `--quiet` flag if that makes things "
                             "easier to read. So, in summary: You can probably ignore this warning. But if you want more info: run setup again "
                             "with `--reset --debug --quiet` to see exactly which modules had issues, or run with `--reset --debug` to see exactly "
                             "which lines in which modules had issues. Anvi'o developers thank you for your attention and patience 😇")
            self.run.info("Bad line splitting (usually due to rogue or missing spaces)", len(self.parsing_error_dict["bad_line_splitting"]))
            self.run.info("Bad KEGG code format (usually not correctable)", len(self.parsing_error_dict["bad_kegg_code_format"]))

        # give some run info
        self.run.info('Modules database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of KEGG modules', num_modules_parsed, quiet=self.quiet)
        self.run.info('Number of entries', mod_table.get_total_entries(), quiet=self.quiet)
        self.run.info('Number of parsing errors (corrected)', self.num_corrected_errors, quiet=self.quiet)
        self.run.info('Number of parsing errors (uncorrected)', self.num_uncorrected_errors, quiet=self.quiet)

        # record some useful metadata
        self.db.set_meta_value('db_type', 'modules')
        self.db.set_meta_value('num_modules', num_modules_parsed)
        self.db.set_meta_value('total_entries', mod_table.get_total_entries())
        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value('hash', self.get_db_content_hash())
        self.db.set_meta_value('version', t.metabolic_modules_db_version)

        self.db.disconnect()


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
            self.run.warning("Just so you know, we tried to fetch data from the KEGG Modules database for the data_name field %s "
                             "and KEGG module %s, but didn't come up with anything, so an empty list is being returned. This may "
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
            self.run.warning("Just so you know, we tried to fetch data definitions from the KEGG Modules database for the data_name field %s "
                             "and KEGG module %s, but didn't come up with anything, so an empty list is being returned. This may "
                             "cause errors down the line, and if so we're very sorry for that.")

        return data_defs_to_ret


    def get_all_modules_as_list(self):
        """This function returns a list of all modules in the DB."""
        return self.db.get_single_column_from_table(self.module_table_name, 'module', unique=True)


    def get_all_knums_as_list(self):
        """This function returns a list of all KO numbers in the DB."""
        where_clause_string = "data_name = 'ORTHOLOGY'"
        return self.db.get_single_column_from_table(self.module_table_name, 'data_value', unique=True, where_clause=where_clause_string)


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


    def get_kegg_module_class_dict(self, mnum):
        """This function returns a dictionary of values in the CLASS field for a specific module

        It really exists only for convenience to put together the data fetch and parsing functions.
        """

        # there should only be one CLASS line per module, so we extract the first list element
        class_value = self.get_data_value_entries_for_module_by_data_name(mnum, "CLASS")[0]
        return self.parse_kegg_class_value(class_value)


    def get_kegg_module_definition(self, mnum):
        """This function returns module DEFINITION fields as one string"""

        def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
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

        """
        compound_to_name_dict = self.get_compound_dict_for_module(mnum)
        substrate_compounds, intermediate_compounds, product_compounds = self.get_kegg_module_compound_lists(mnum)

        substrate_name_list = [compound_to_name_dict[c] for c in substrate_compounds]
        intermediate_name_list = [compound_to_name_dict[c] for c in intermediate_compounds]
        product_name_list = [compound_to_name_dict[c] for c in product_compounds]

        return substrate_name_list, intermediate_name_list, product_name_list

######### MODULE DEFINITION UNROLLING FUNCTIONS #########
    def get_top_level_steps_in_module_definition(self, mnum):
        """This function access the DEFINITION line of a KEGG Module and returns the top-level steps as a list

        A 'top-level' step is one that you get by splitting on spaces (but not spaces in parentheses) just once -
        ie, the 'first layer' when unrolling the module.
        """

        def_string = self.get_kegg_module_definition(mnum)
        return self.split_by_delim_not_within_parens(def_string, " ")


    def unroll_module_definition(self, mnum):
        """This function accesses the DEFINITION line of a KEGG Module, unrolls it into all possible paths through the module, and
        returns the list of all paths.

        This is a driver for the recursive functions that do the actual unrolling of each definition line.
        """

        def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        combined_def_line = ""
        for d in def_lines:
            d = d.strip()
            combined_def_line += d + " "
        combined_def_line = combined_def_line.strip()
        def_line_paths = self.recursive_definition_unroller(combined_def_line)

        return def_line_paths


    def split_by_delim_not_within_parens(self, d, delims, return_delims=False):
        """Takes a string, and splits it on the given delimiter(s) as long as the delimeter is not within parentheses.

        This function exists because regular expressions don't handle nested parentheses very well. It is used in the
        recursive module definition unrolling functions to split module steps, but it is generically written in case
        it could have other uses in the future.

        The function can also be used to determine if the parentheses in the string are unbalanced (it will return False
        instead of the list of splits in this situation)

        PARAMETERS
        ==========
        d : str
            string to split
        delims : str or list of str
            a single delimiter, or a list of delimiters, to split on
        return_delims : boolean
            if this is true then the list of delimiters found between each split is also returned

        RETURNS
        =======
        If parentheses are unbalanced in the string, this function returns False. Otherwise:
        splits : list
            strings that were split from d
        delim_list : list
            delimiters that were found between each split (only returned if return_delims is True)
        """

        parens_level = 0
        last_split_index = 0
        splits = []
        delim_list = []
        for i in range(len(d)):
            # only split if not within parentheses
            if d[i] in delims and parens_level == 0:
                splits.append(d[last_split_index:i])
                delim_list.append(d[i])
                last_split_index = i + 1 # we add 1 here to skip the space
            elif d[i] == "(":
                parens_level += 1
            elif d[i] == ")":
                parens_level -= 1

            # if parentheses become unbalanced, return False to indicate this
            if parens_level < 0:
                return False
        splits.append(d[last_split_index:len(d)])

        if return_delims:
            return splits, delim_list
        return splits


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

        split_steps = self.split_by_delim_not_within_parens(step, " ")
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
                    comma_substeps = self.split_by_delim_not_within_parens(s[1:-1], ",")
                    if not comma_substeps: # if it doesn't work, try without removing surrounding parentheses
                        comma_substeps = self.split_by_delim_not_within_parens(s, ",")
                    space_substeps = self.split_by_delim_not_within_parens(s[1:-1], " ")
                    if not space_substeps:
                        space_substeps = self.split_by_delim_not_within_parens(s, " ")
                else:
                    comma_substeps = self.split_by_delim_not_within_parens(s, ",")
                    space_substeps = self.split_by_delim_not_within_parens(s, " ")

                # complex case: no commas OR spaces outside parentheses so this is a protein complex rather than a compound step
                if len(comma_substeps) == 1 and len(space_substeps) == 1:
                    complex_components, delimiters = self.split_by_delim_not_within_parens(s, ["+","-"], return_delims=True)
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
            substeps = self.split_by_delim_not_within_parens(step[1:-1], ",")
            if not substeps: # if it doesn't work, try without removing surrounding parentheses
                substeps = self.split_by_delim_not_within_parens(step, ",")
        else:
            substeps = self.split_by_delim_not_within_parens(step, ",")

        alt_path_list = []
        for s in substeps:
            alt_paths_from_substep = self.recursive_definition_unroller(s)
            for a in alt_paths_from_substep:
                alt_path_list.append(a)

        return alt_path_list


class KeggModulesTable:
    """This class defines operations for creating the KEGG Modules table in Modules.db"""

    def __init__(self, mod_table_name = None):
        """"""
        self.db_entries = []
        self.total_entries = 0

        if mod_table_name:
            self.module_table_name = mod_table_name
        else:
            raise ConfigError("Beep Beep. Warning. KeggModulesTable was initialized without knowing its own name.")


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
        self.include_ungrouped = True if A('include_ungrouped') else False
        self.include_missing = True if A('include_samples_missing_from_groups_txt') else False

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
        and 'kegg_module', respectively.
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

        # make sure we have all the columns we need in modules mode output, since this output can be customized
        required_modules_txt_headers = ['kegg_module', 'module_completeness', 'module_name']
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

        required_groups_txt_headers = ['sample', 'group']
        sample_groups_dict = utils.get_TAB_delimited_file_as_dictionary(self.groups_txt, expected_fields=required_groups_txt_headers)
        samples_to_groups_dict = {samp : sample_groups_dict[samp]['group'] for samp in sample_groups_dict.keys()}

        # make sure the samples all have a group
        samples_with_none_group = []
        for s,g in samples_to_groups_dict.items():
            if not g:
                samples_with_none_group.append(s)
                if self.include_ungrouped:
                    samples_to_groups_dict[s] = 'UNGROUPED'

        if not self.include_ungrouped:
            for s in samples_with_none_group:
                samples_to_groups_dict.pop(s)

        if samples_with_none_group:
            self.progress.reset()
            none_group_str = ", ".join(samples_with_none_group)
            if self.include_ungrouped:
                self.run.warning("Some samples in your groups-txt did not have a group, but since you elected to --include-ungrouped, "
                                 "we will consider all of those samples to belong to one group called 'UNGROUPED'. Here are those "
                                 f"UNGROUPED samples: {none_group_str}")
            else:
                self.run.warning("Some samples in your groups-txt did not have a group, and we will ignore those samples. If you "
                                 "want them to be included in the analysis (but without assigning a group), you can simply re-run "
                                 "this program with the --include-ungrouped flag. Now. Here are the samples we will be ignoring: "
                                 f"{none_group_str}")

        # sanity check for mismatch between modules-txt and groups-txt
        sample_names_in_modules_txt = set(modules_df[self.sample_header_in_modules_txt].unique())
        sample_names_in_groups_txt = set(sample_groups_dict.keys())
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
                                "these samples to belong to one group called 'UNGROUPED'. If you wish to ignore these samples instead, please run again "
                                "without the --include-ungrouped parameter. "
                                "Here are the UNGROUPED samples that we will consider as one big happy family: "
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
        module_list = modules_df['kegg_module'].unique()

        output_dict = {}
        header_list = ['KEGG_MODULE', 'accession', 'sample_ids', 'associated_groups']
        for c in group_list:
            header_list.append(f"p_{c}")
            header_list.append(f"N_{c}")

        for mod_num in module_list:
            query_string = f"kegg_module == '{mod_num}' and module_completeness >= {self.module_completion_threshold}"
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
            if not self.include_ungrouped:
                for s in samples_with_none_group:
                    if s in samples_with_mod_list:
                        samples_with_mod_list.remove(s)
            if len(samples_with_mod_list) == 0:
                continue

            mod_name = samples_with_mod_df['module_name'][0]
            output_dict[mod_name] = {}
            output_dict[mod_name]['KEGG_MODULE'] = mod_name
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

        cmd = f"anvi-script-enrichment-stats --input {enrichment_input_path} --output {self.output_file_path}"
        log_file_path = filesnpaths.get_temp_file_path()

        self.progress.update("Running Amy's enrichment")
        utils.run_command(cmd, log_file_path)
        self.progress.end()
        if not filesnpaths.is_file_exists(self.output_file_path, dont_raise=True):
            raise ConfigError("It looks like something went wrong during the enrichment analysis. "
                              f"We don't know what happened, but this log file could contain some clues: {log_file_path}")

        if filesnpaths.is_file_empty(self.output_file_path):
            raise ConfigError("It looks like something went wrong during the functional enrichment analysis. "
                              "An output file was created, but it is empty "
                              f"We don't know why this happened, but this log file could contain some clues: {log_file_path}")


        self.run.info('Enrichment results', self.output_file_path)

        if not anvio.DEBUG:
            # get rid of the temporary files now that we are done
            os.remove(enrichment_input_path)
            os.remove(log_file_path)
        else:
            self.run.info('Enrichment log file:', log_file_path)
            self.run.warning("Because you ran this script with the --debug flag, the temporary files are kept. Please "
                             "consider cleaning them up when you are done taking a look at them. Here they are: "
                             f"{enrichment_input_path}, {log_file_path}")
