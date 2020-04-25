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
from anvio.constants import KEGG_SETUP_INTERVAL


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Iva Veseli"
__email__ = "iveseli@uchicago.edu"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class KeggContext(object):
    """The purpose of this base class is to define shared functions and file paths for all KEGG operations."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # default data directory will be called KEGG and will store the KEGG Module data as well
        self.kegg_data_dir = A('kegg_data_dir') or os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')
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


class KeggSetup(KeggContext):
    """Class for setting up KEGG Kofam HMM profiles and modules.

    It performs sanity checks and downloads, unpacks, and prepares the profiles for later use by `hmmscan`.
    It also downloads module files and creates the MODULES.db.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-setup-kegg-kofams
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists('hmmpress')

        filesnpaths.is_output_dir_writable(os.path.dirname(self.kegg_data_dir))

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.kegg_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.hmm_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.module_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.pathway_data_dir, delete_if_exists=args.reset)

        # ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        # dictionary mapping downloaded file name to final decompressed file name or folder location
        self.files = {'ko_list.gz': self.ko_list_file_path, 'profiles.tar.gz': self.kegg_data_dir}

        # Kegg module text files
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

        for file_name in self.files.keys():
            utils.download_file(self.database_url + '/' + file_name,
                os.path.join(self.kegg_data_dir, file_name), progress=self.progress, run=self.run)


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
                                      "meantime use an older version of the KEGG data directory (if you have one)." % (self.kegg_module_file, first_char))
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
                                      "meantime use an older version of the KEGG data directory (if you have one)." % (self.kegg_pathway_file, first_char))
        self.progress.end()


    def download_modules(self):
        """This function downloads the KEGG modules.

        To do so, it also processes the KEGG module file into a dictionary via the process_module_file() function.
        To verify that each file has been downloaded properly, we check that the last line is '///'.
        """

        self.run.info("KEGG Module Database URL", self.kegg_rest_api_get)

        # download the kegg module file, which lists all modules
        utils.download_file(self.kegg_module_download_path, self.kegg_module_file, progress=self.progress, run=self.run)

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
                                  "to be '///', but instead it was %s." % (file_path, last_line))


    def download_pathways(self):
        """This function downloads the KEGG Pathways.

        To do so, it first processes a KEGG file containing pathway and map identifiers into a dictionary via the process_pathway_file()
        function. To verify that each file has been downloaded properly, we check that the last line is '///'.
        """

        # note that this is the same as the REST API for modules - perhaps at some point this should be printed elsewhere so we don't repeat ourselves.
        self.run.info("KEGG Pathway Database URL", self.kegg_rest_api_get)

        # download the kegg pathway file, which lists all modules
        utils.download_file(self.kegg_pathway_download_path, self.kegg_pathway_file, progress=self.progress, run=self.run)

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
                                  "to be '///', but instead it was %s." % (file_path, last_line))


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
                                      "flag." % (hmm_path))


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
        remove_old_files = not anvio.DEBUG # if we are running in debug mode, we will not remove the individual hmm files after concatenation
        if no_kofam_file_list:
            utils.concatenate_files(no_kofam_path, no_kofam_file_list, remove_concatenated_files=remove_old_files)
            self.progress.reset()
            self.run.warning("Please note that while anvi'o was building your databases, she found %d "
                             "HMM profiles that did not have any matching KOfam entries. We have removed those HMM "
                             "profiles from the final database. You can find them under the directory '%s'."
                             % (len(no_kofam_file_list), self.orphan_data_dir))

        if no_threshold_file_list:
            utils.concatenate_files(no_threshold_path, no_threshold_file_list, remove_concatenated_files=remove_old_files)
            self.progress.reset()
            self.run.warning("Please note that while anvi'o was building your databases, she found %d "
                             "KOfam entries that did not have any threshold to remove weak hits. We have removed those HMM "
                             "profiles from the final database. You can find them under the directory '%s'."
                             % (len(no_threshold_file_list), self.orphan_data_dir))

        if no_data_file_list:
            utils.concatenate_files(no_data_path, no_data_file_list, remove_concatenated_files=remove_old_files)
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
        """This function creates the Modules DB from the Kegg Module files."""

        mod_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args, module_dictionary=self.module_dict, run=run, progress=progress)
        mod_db.create()


    def setup_profiles(self):
        """This is a driver function which executes the KEGG setup process.

        It downloads, decompresses, and hmmpresses the KOfam profiles.
        It also downloads and processes the KEGG Module files into the MODULES.db.
        """

        self.download_profiles()
        self.decompress_files()
        self.download_modules()
        self.setup_ko_dict()
        self.run_hmmpress()
        self.setup_modules_db()


class KeggRunHMMs(KeggContext):
    """ Class for running `hmmscan` against the KOfam database and adding the resulting hits to contigs DB for later metabolism prediction.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-run-kegg-kofams
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.contigs_db_path = args.contigs_db
        self.num_threads = args.num_threads
        self.hmm_program = args.hmmer_program or 'hmmsearch'
        self.ko_dict = None # should be set up by setup_ko_dict()

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists(self.hmm_program)

        # verify that Kofam HMM profiles have been set up
        if not os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("Anvi'o is unable to find the Kofam.hmm file at %s. This can happen one of two ways. Either you "
                              "didn't specify the correct KEGG data directory using the flag --kegg-data-dir, or you haven't "
                              "yet set up the Kofam data by running `anvi-setup-kegg-kofams`. Hopefully you now know what to do "
                              "to fix this problem. :) " % self.kegg_data_dir)

        utils.is_contigs_db(self.contigs_db_path)

        self.setup_ko_dict() # read the ko_list file into self.ko_dict

        # load existing kegg modules db
        self.kegg_modules_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args)


    def set_hash_in_contigs_db(self):
        """Modifies the contigs DB self table to indicate which MODULES.db has been used to annotate it."""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        current_module_hash_in_contigs_db = contigs_db.db.get_meta_value('modules_db_hash', return_none_if_not_in_table=True)

        if current_module_hash_in_contigs_db and not self.just_do_it:
            raise ConfigError("The contigs database (%s) has already been annotated with KOfam hits. If you really want to "
                              "overwrite these annotations with new ones, please re-run the command with the flag --just-do-it. "
                              "For those who need this information, the Modules DB used to annotate this contigs database previously "
                              "had the following hash: %s" % (self.contigs_db_path, current_module_hash_in_contigs_db))

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


    def process_kofam_hmms(self):
        """This is a driver function for running HMMs against the KOfam database and processing the hits into the provided contigs DB."""

        tmp_directory_path = filesnpaths.get_temp_directory_path()
        contigs_db = ContigsSuperclass(self.args) # initialize contigs db

        # mark contigs db with hash of modules.db content for version tracking
        # this function also includes a safety check for previous annotations so that people don't overwrite those if they don't want to
        self.set_hash_in_contigs_db()

        # get AA sequences as FASTA
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['AA:GENE'],
                                                                   simple_headers=True,
                                                                   rna_alphabet=False,
                                                                   report_aa_sequences=True)

        # run hmmscan
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmscan('KOfam', 'AA', 'GENE', None, None, len(self.ko_dict), self.kofam_hmm_file_path, None, None)

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

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
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})
            return

        # parse hmmscan output
        parser = parser_modules['search']['hmmscan'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results(noise_cutoff_dict=self.ko_dict)

        # add functions and KEGG modules info to database
        functions_dict = {}
        kegg_module_names_dict = {}
        kegg_module_classes_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            knum = hmm_hit['gene_name']
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'KOfam',
                'accession': knum,
                'function': self.get_annotation_from_ko_dict(hmm_hit['gene_name'], ok_if_missing_from_dict=True),
                'e_value': hmm_hit['e_value'],
            }

            # add associated KEGG module information to database
            mods = self.kegg_modules_db.get_modules_for_knum(knum)
            names = self.kegg_modules_db.get_module_names_for_knum(knum)
            classes = self.kegg_modules_db.get_module_classes_for_knum_as_list(knum)

            # FIXME? some KOs are not associated with modules. Should we report this?
            if mods:
                mod_annotation = "!!!".join(mods)
                mod_class_annotation = "!!!".join(classes) # why do we split by '!!!'? Because that is how it is done in COGs. So so sorry. :'(
                mod_name_annotation = ""

                for mod in mods:
                    if mod_name_annotation:
                        mod_name_annotation += "!!!" + names[mod]
                    else:
                        mod_name_annotation = names[mod]

                kegg_module_names_dict[counter] = {
                    'gene_callers_id': hmm_hit['gene_callers_id'],
                    'source': 'KEGG_Module',
                    'accession': mod_annotation,
                    'function': mod_name_annotation,
                    'e_value': None,
                }
                kegg_module_classes_dict[counter] = {
                    'gene_callers_id': hmm_hit['gene_callers_id'],
                    'source': 'KEGG_Class',
                    'accession': mod_annotation,
                    'function': mod_class_annotation,
                    'e_value': None,
                }

            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
            gene_function_calls_table.create(kegg_module_names_dict)
            gene_function_calls_table.create(kegg_module_classes_dict)
        else:
            self.run.warning("KOfam class has no hits to process. Returning empty handed, but still adding KOfam as "
                             "a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})


        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                        "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single("Cleaning up the temp directory (you can use `--debug` if you would "
                            "like to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()


class KeggMetabolismEstimator(KeggContext):
    """ Class for reconstructing/estimating metabolism based on hits to KEGG databases.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-estimate-kegg-metabolism
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
        self.metagenome_mode = True if A('metagenome_mode') else False
        self.completeness_threshold = A('module_completion_threshold') or 0.75
        self.output_file_prefix = A('output_file_prefix') or "kegg-metabolism"
        self.contigs_db_project_name = "Unknown"
        self.write_dict_to_json = True if A('get_raw_data_as_json') else False
        self.json_output_file_path = A('get_raw_data_as_json')
        self.store_json_without_estimation = True if A('store_json_without_estimation') else False
        self.estimate_from_json = A('estimate_from_json') or None

        if not self.estimate_from_json and not self.contigs_db_path:
            raise ConfigError("NO INPUT PROVIDED. You must provide (at least) a contigs database to this program, unless you are using the --estimate-from-json "
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

        if self.bin_id or self.bin_ids_file or self.collection_name and not self.profile_db_path:
            raise ConfigError("You have requested metabolism estimation for a bin or set of bins, but you haven't provided "
                              "a profiles database. Unfortunately, this just does not work. Please try again.")

        if self.profile_db_path and not self.collection_name:
            raise ConfigError("If you provide a profiles DB, you should also provide a collection name.")

        if self.store_json_without_estimation and not self.json_output_file_path:
            raise ConfigError("Whoops. You seem to want to store the metabolism dictionary in a JSON file, but you haven't provided the name of that file. "
                              "Please use the --get-raw-data-as-json flag to do so.")
        if self.store_json_without_estimation and self.estimate_from_json:
            raise ConfigError("It is impossible to both estimate metabolism from JSON data and produce a JSON file without estimation at the same time... "
                              "anvi'o is judging you SO hard right now.")


        # init the base class
        KeggContext.__init__(self, self.args)

        if not self.estimate_from_json:
            utils.is_contigs_db(self.contigs_db_path)

        # load existing kegg modules db
        if not os.path.exists(self.kegg_modules_db_path):
            raise ConfigError("It appears that a modules database (%s) does not exist in the KEGG data directory %s. "
                              "Perhaps you need to specify a different KEGG directory using --kegg-data-dir. Or perhaps you didn't run "
                              "`anvi-setup-kegg-kofams`, though we are not sure how you got to this point in that case "
                              "since you also cannot run `anvi-run-kegg-kofams` without first having run KEGG setup. But fine. Hopefully "
                              "you now know what you need to do to make this message go away." % ("MODULES.db", self.kegg_data_dir))
        self.kegg_modules_db = KeggModulesDatabase(self.kegg_modules_db_path, args=self.args)

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
        self.contigs_db_project_name = contigs_db.meta['project_name']

        # sanity check that contigs db was annotated with same version of MODULES.db that will be used for metabolism estimation
        contigs_db_mod_hash = contigs_db.meta['modules_db_hash']
        mod_db_hash = self.kegg_modules_db.db.get_meta_value('hash')
        if contigs_db_mod_hash != mod_db_hash:
            raise ConfigError("The contigs DB that you are working with has been annotated with a different version of the MODULES.db than you are working with now. "
                                "Perhaps you updated your KEGG setup after running `anvi-run-kegg-kofams` on this contigs DB? Or maybe you have multiple KEGG data "
                                "directories set up on your computer, and the one you are using now is different from the one that you used for `anvi-run-kegg-kofams`? "
                                "Well. The solution to the first problem is to re-run `anvi-run-kegg-kofams` on the contigs DB (%s) using the updated MODULES.db "
                                "(located in the KEGG data directory %s). The solution to the second problem is to specify the appropriate KEGG data directory using "
                                "the --kegg-data-dir flag. If neither of those things make this work, then you should contact the developers to see if they can help you "
                                "figure this out." % (self.contigs_db_path, self.kegg_data_dir))

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
        assert len(gene_calls_splits_dict.keys()) == len(genes_in_splits)
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

        return kofam_gene_split_contig


    def mark_kos_present_for_list_of_splits(self, kofam_hits_in_splits, split_list=None, bin_name=None):
        """This function generates a bin-level dictionary of dictionaries, which associates modules with the KOs
        that are present in the bin for each module.

        The structure of the dictionary is like this example:
        {mnum: {"gene_caller_ids" : set([132, 133, 431, 6777])
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
            initialized metabolism completeness dictionary for the list of splits (genome, metagenome, or bin) provided
        """

        bin_level_module_dict = {}

        if anvio.DEBUG:
            self.run.info("Marking KOs present for bin", bin_name)
            self.run.info("Number of splits", len(split_list))

        # initialize all modules with empty lists and dicts for kos, gene calls
        modules = self.kegg_modules_db.get_all_modules_as_list()
        for mnum in modules:
            bin_level_module_dict[mnum] = {"gene_caller_ids" : set(),
                                           "kofam_hits" : {},
                                           "genes_to_contigs" : {},
                                           "contigs_to_genes" : {}
                                          }

        kos_not_in_modules = []
        for ko, gene_call_id, split, contig in kofam_hits_in_splits:
            present_in_mods = self.kegg_modules_db.get_modules_for_knum(ko)
            if not present_in_mods:
                kos_not_in_modules.append(ko)
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

        # TODO: at some point I think we should save these KOs somewhere so that the user can look at them manually
        if anvio.DEBUG:
            self.run.info("KOs processed", "%d in bin" % len(kofam_hits_in_splits))
            if kos_not_in_modules:
                self.run.warning("Just so you know, the following KOfam hits did not belong to any KEGG modules in the MODULES.db: %s"
                % ", ".join(kos_not_in_modules))

        return bin_level_module_dict


    def compute_module_completeness_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the completeness of the specified module within the given bin metabolism dictionary.

        To do this, it unrolls the module definition into a list of all possible paths, where each path is a list of atomic steps.
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
        "paths"                         a list of all possible paths (each is a list of atomic) through the module DEFINITION
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

        # unroll the module definition to get all possible paths
        meta_dict_for_bin[mnum]["paths"] = self.kegg_modules_db.unroll_module_definition(mnum)
        meta_dict_for_bin[mnum]["pathway_completeness"] = []

        for p in meta_dict_for_bin[mnum]["paths"]:
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
            meta_dict_for_bin[mnum]["most_complete_paths"] = [meta_dict_for_bin[mnum]["paths"][i] for i, pc in enumerate(meta_dict_for_bin[mnum]["pathway_completeness"]) if pc == meta_dict_for_bin[mnum]["percent_complete"]]
        else:
            meta_dict_for_bin[mnum]["most_complete_paths"] = []


        if anvio.DEBUG and len(meta_dict_for_bin[mnum]["most_complete_paths"]) > 1:
            self.run.warning("Found %d complete paths for module %s with completeness %s. " % (len(meta_dict_for_bin[mnum]["most_complete_paths"]), mnum, meta_dict_for_bin[mnum]["percent_complete"]),
                            header='DEBUG OUTPUT', lc='yellow')
        over_complete_threshold = True if meta_dict_for_bin[mnum]["percent_complete"] >= self.completeness_threshold else False
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

        for i in range(len(meta_dict_for_bin[mod]["paths"])):
            p = meta_dict_for_bin[mod]["paths"][i]
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
            meta_dict_for_bin[mod]["most_complete_paths"] = [meta_dict_for_bin[mod]["paths"][i] for i, pc in enumerate(meta_dict_for_bin[mod]["pathway_completeness"]) if pc == meta_dict_for_bin[mod]["percent_complete"]]
        else:
            meta_dict_for_bin[mod]["most_complete_paths"] = []

        now_complete = True if meta_dict_for_bin[mod]["percent_complete"] >= self.completeness_threshold else False
        meta_dict_for_bin[mod]["complete"] = now_complete
        if now_complete:
            meta_dict_for_bin["num_complete_modules"] += 1

        return now_complete


    def compute_naive_redundancy_for_path(self, num_ko_hits_in_path_dict):
        """This function computes a naive redundancy measure for a module path, given the number of hits per KO in the path.

        naive redundancy = # extra hits / len(path) where a hit is "extra" if it is not the first hit to the KO.
        """

        extra_hits = [num_ko_hits_in_path_dict[ko] - 1 for ko in num_ko_hits_in_path_dict if num_ko_hits_in_path_dict[ko] > 1 ]
        return sum(extra_hits)/len(num_ko_hits_in_path_dict.keys())


    def compute_copywise_redundancy_for_path(self, num_ko_hits_in_path_dict, aggregation_measure="average"):
        """This function computes redundancy based on the completeness of each extra copy of a path.

        The 'base' redundancy score is determined by the number of extra copies with 100% completeness.
        The completeness measurements of all other extra copies are aggregated (using the aggregation_measure) and
        added to this 'base' redundancy to get the overall path redundancy.
        """

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
            elif aggregation_measure == "knee":
                raise ConfigError("aggregation measure 'knee' not implemented yet")
            else:
                raise ConfigError("The function compute_copywise_redundancy_for_path() doesn't know how to handle the aggregation measure '%s'", aggregation_measure)

        return (base_redundancy + aggregated_completeness), extra_copy_completeness


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

        return


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


        # notify user of the modules that gave some fishy results
        if not self.quiet:
            if mods_with_nonessential_steps:
                self.run.warning("Please note that anvi'o found one or more non-essential steps in the following KEGG modules: %s.   "
                                 "At this time, we are not counting these steps in our percent completion estimates. But we still kept track of which "
                                 "of these non-essential steps were found to be complete. You can see this information in the output file."
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

        self.run.info("Bin name", bin_name)
        self.run.info("Module completion threshold", self.completeness_threshold)
        self.run.info("Number of complete modules", metabolism_dict_for_list_of_splits["num_complete_modules"])
        if complete_mods:
            self.run.info("Complete modules", ", ".join(complete_mods))

        return metabolism_dict_for_list_of_splits


    def estimate_for_genome(self, kofam_gene_split_contig):
        """This is the metabolism estimation function for a contigs DB that contains a single genome.

        Assuming this contigs DB contains only one genome, it sends all of the splits and their kofam hits to the atomic
        estimation function for processing. It then returns the metabolism completion dictionary for the genome, wrapped in the superdict format.

        PARAMETERS
        ==========
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering

        RETURNS
        =======
        genome_metabolism_dict : dictionary of dictionary of dictionaries
            dictionary mapping genome name to its metabolism completeness dictionary
        """

        genome_metabolism_superdict = {}
        # since all hits belong to one genome, we can take the UNIQUE splits from all the hits
        splits_in_genome = list(set([tpl[2] for tpl in kofam_gene_split_contig]))
        metabolism_dict_for_genome = self.mark_kos_present_for_list_of_splits(kofam_gene_split_contig, split_list=splits_in_genome,
                                                                                                    bin_name=self.contigs_db_project_name)
        if not self.store_json_without_estimation:
            genome_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(metabolism_dict_for_genome, bin_name=self.contigs_db_project_name)
        else:
            genome_metabolism_superdict[self.contigs_db_project_name] = metabolism_dict_for_genome

        return genome_metabolism_superdict


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
        """

        bins_metabolism_superdict = {}

        bin_name_to_split_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.run.info_single("%s split names associated with %s bins of in collection '%s' have been "
                             "successfully recovered 🎊" % (pp(sum([len(v) for v in bin_name_to_split_names_dict.values()])),
                                                           pp(len(bin_name_to_split_names_dict)),
                                                           self.collection_name), nl_before=1)

        for bin_name in bin_name_to_split_names_dict:
            splits_in_bin = bin_name_to_split_names_dict[bin_name]
            ko_in_bin = [tpl for tpl in kofam_gene_split_contig if tpl[2] in splits_in_bin]
            metabolism_dict_for_bin = self.mark_kos_present_for_list_of_splits(ko_in_bin, split_list=splits_in_bin, bin_name=bin_name)

            if not self.store_json_without_estimation:
                bins_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(metabolism_dict_for_bin, bin_name=bin_name)
            else:
                bins_metabolism_superdict[bin_name] = metabolism_dict_for_bin

        return bins_metabolism_superdict


    def estimate_for_contigs_db_for_metagenome(self, kofam_gene_split_contig):
        """This function handles metabolism estimation for an entire metagenome.

        Similar to isolate genomes, we treat the entire metagenome as one big 'bin'. This means that there
        will be a large amount of redundancy (repeated pathways) due to the presence of multiple populations
        in the metagenome.

        In fact, because we essentially consider the metagenome to be one big genome, this function is exactly the same
        as estimate_for_genome(). Why is it a separate function? Well, because we may eventually want to do something
        differently here.

        PARAMETERS
        ==========
        kofam_gene_split_contig : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering

        RETURNS
        =======
        metagenome_metabolism_superdict : dictionary of dictionary of dictionaries
            dictionary mapping metagenome name to its metabolism completeness dictionary
        """

        metagenome_metabolism_superdict = {}
        # since we consider all the hits in the metagenome collectively, we can take the UNIQUE splits from all the hits
        splits_in_metagenome = list(set([tpl[2] for tpl in kofam_gene_split_contig]))
        metabolism_dict_for_metagenome = self.mark_kos_present_for_list_of_splits(kofam_gene_split_contig, split_list=splits_in_metagenome,
                                                                                                    bin_name=self.contigs_db_project_name)
        if not self.store_json_without_estimation:
            metagenome_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(metabolism_dict_for_metagenome, bin_name=self.contigs_db_project_name)
        else:
            metagenome_metabolism_superdict[self.contigs_db_project_name] = metabolism_dict_for_metagenome

        return metagenome_metabolism_superdict


    def estimate_metabolism_from_json_data(self):
        """This function runs the estimation functions on data obtained from a provided JSON file"""

        self.run.info("JSON input file", self.estimate_from_json)

        filesnpaths.is_file_json_formatted(self.estimate_from_json)
        kegg_metabolism_superdict = json.load(open(self.estimate_from_json), parse_int=int)
        new_kegg_metabolism_superdict = {}

        expected_keys_for_module = {"gene_caller_ids", "kofam_hits", "genes_to_contigs", "contigs_to_genes"}
        bins_found = []
        additional_keys = set([])

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


        if not self.quiet and additional_keys:
            self.run.warning("Just to let you know, we found the following module-level keys in your JSON file that were totally ignored during metabolism estimation "
                             "(no harm was done by including them): %s" % (additional_keys))

        self.run.info("Bins/genomes/metagenomes found", ", ".join(bins_found))
        return new_kegg_metabolism_superdict


    def estimate_metabolism(self):
        """This is the driver function for estimating metabolism.

        It will decide what to do based on whether the input contigs DB is a genome or metagenome.
        It returns the metabolism superdict which contains a metabolism completion dictionary for each genome/bin in the contigs db.
        The metabolism completion dictionary is keyed by KEGG module number, with a few exceptions for summary data (ie, 'num_complete_modules').
        """

        kegg_metabolism_superdict = {}

        if self.estimate_from_json:
            kegg_metabolism_superdict = self.estimate_metabolism_from_json_data()
        else:

            kofam_hits_info = self.init_hits_and_splits()

            if self.profile_db_path and not self.metagenome_mode:
                kegg_metabolism_superdict = self.estimate_for_bins_in_collection(kofam_hits_info)
            elif not self.profile_db_path and not self.metagenome_mode:
                kegg_metabolism_superdict = self.estimate_for_genome(kofam_hits_info)
            elif self.metagenome_mode:
                kegg_metabolism_superdict = self.estimate_for_contigs_db_for_metagenome(kofam_hits_info)
            else:
                raise ConfigError("This class doesn't know how to deal with that yet :/")

        if not self.store_json_without_estimation:
            self.store_kegg_metabolism_superdict(kegg_metabolism_superdict)
        if self.write_dict_to_json:
            self.store_metabolism_superdict_as_json(kegg_metabolism_superdict, self.json_output_file_path + ".json")


    def store_kegg_metabolism_superdict(self, kegg_superdict):
        """This function writes the metabolism superdict to a tab-delimited file, and also generates a file summarizing the complete modules.

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
                "paths":             [['K00033','K01057','K02222'], ['K00033','K01057','K00036'], ...]
                "pathway_completeness":     [0.66, 0.66, ...]
                "present_nonessential_kos":      []
                "most_complete_paths":           [['K00033','K01057','K02222'], ['K00033','K01057','K00036'], ...]
                "percent_complete":              0.66
                "complete":                      False
                                      }

        To distill this information into one line, we need to convert the dictionary on-the-fly to a dict of dicts,
        where each bin-module-path-kofam_hit-gene_caller_id is keyed by an arbitrary integer. There will be a lot of redundant information
        in the rows.

        The complete modules summary file includes only a portion of the information in the metabolism dictionary. Its purpose is to give the user
        quick access to the complete modules in each bin. Every bin-module pair in this file is keyed by an arbitrary integer (with no relation to the
        id in the other file).
        """

        hits_output_path = self.output_file_prefix + "-all_kofam_hits.txt"
        complete_module_summary_path = self.output_file_prefix + "-complete_modules.txt"

        name_header = None
        if self.profile_db_path and not self.metagenome_mode:
            name_header = "bin_name"
        elif not self.profile_db_path and not self.metagenome_mode:
            name_header = "genome_name"
        elif self.metagenome_mode:
            name_header = "metagenome_name"

        header_list = ["unique_id", name_header, "kegg_module", "module_is_complete", "module_completeness",
        "path_id", "path", "path_completeness", "kofam_hit_in_path", "gene_caller_id", "contig"]
        summary_header_list = ["unique_id", name_header, "kegg_module","module_completeness", "module_name", "module_class",
        "module_category", "module_subcategory"]

        d = {}
        cm_summary = {}
        unique_id = 0
        summary_unique_id = 0
        for bin, mod_dict in kegg_superdict.items():
            for mnum, c_dict in mod_dict.items():
                if mnum == "num_complete_modules":
                    continue


                if c_dict["complete"]:
                    cm_summary[summary_unique_id] = {}
                    cm_summary[summary_unique_id][name_header] = bin
                    cm_summary[summary_unique_id]["kegg_module"] = mnum
                    cm_summary[summary_unique_id]["module_completeness"] = c_dict["percent_complete"]
                    cm_summary[summary_unique_id]["module_name"] = self.kegg_modules_db.get_module_name(mnum)
                    mnum_class_dict = self.kegg_modules_db.get_kegg_module_class_dict(mnum)
                    cm_summary[summary_unique_id]["module_class"] = mnum_class_dict["class"]
                    cm_summary[summary_unique_id]["module_category"] = mnum_class_dict["category"]
                    cm_summary[summary_unique_id]["module_subcategory"] = mnum_class_dict["subcategory"]

                    summary_unique_id += 1

                for p_index in range(len(c_dict['paths'])):
                    p = c_dict['paths'][p_index]

                    for ko in c_dict['kofam_hits']:
                        if ko not in p:
                            continue

                        for gc_id in c_dict["kofam_hits"][ko]:
                            d[unique_id] = {}
                            d[unique_id][name_header] = bin
                            d[unique_id]["kegg_module"] = mnum
                            d[unique_id]["module_is_complete"] = c_dict["complete"]
                            d[unique_id]["module_completeness"] = c_dict["percent_complete"]
                            d[unique_id]["path_id"] = p_index
                            d[unique_id]["path"] = ",".join(p)
                            d[unique_id]["path_completeness"] = c_dict["pathway_completeness"][p_index]
                            d[unique_id]["kofam_hit_in_path"] = ko
                            d[unique_id]["gene_caller_id"] = gc_id
                            d[unique_id]["contig"] = c_dict["genes_to_contigs"][gc_id]

                            unique_id += 1

        utils.store_dict_as_TAB_delimited_file(d, hits_output_path, key_header="unique_id", headers=header_list)
        self.run.info("Kofam hits output file", hits_output_path, nl_before=1)
        utils.store_dict_as_TAB_delimited_file(cm_summary, complete_module_summary_path, key_header="unique_id", headers=summary_header_list)
        self.run.info("Complete modules summary file", complete_module_summary_path)


    def store_metabolism_superdict_as_json(self, kegg_superdict, file_path):
        """This function writes the metabolism superdict into one json file."""

        def set_to_list(obj):
            if isinstance(obj, set):
                return list(obj)

        filesnpaths.is_output_file_writable(file_path)
        open(file_path, 'w').write(json.dumps(kegg_superdict, indent=4, default=set_to_list))
        self.run.info("JSON Output", file_path)


class KeggModulesDatabase(KeggContext):
    """To create or access a Modules DB.

    This DB should be created in the Kegg Data folder during KEGG setup, and will be populated with information from the
    Kegg Module files.
    """

    def __init__(self, db_path, args, module_dictionary=None, run=run, progress=progress, quiet=False):
        self.db = None
        self.db_path = db_path
        self.module_dict = module_dictionary
        self.run = run
        self.progress = progress
        self.quiet = quiet

        # init the base class for access to shared paths and such
        KeggContext.__init__(self, args)

        # modules table info
        # I wonder if these should be moved to the tables __init__.py at some point?
        self.module_table_name = "kegg_modules"
        self.module_table_structure = ['module', 'data_name', 'data_value', 'data_definition', 'line']
        self.module_table_types     = [ 'str'  ,   'str'    ,     'str'   ,       'str'      ,'numeric' ]

        if os.path.exists(self.db_path):
            utils.is_kegg_modules_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__kegg_modules_version__, new_database=False)

            self.run.info('Modules database', 'An existing database, %s, has been loaded.' % self.db_path, quiet=self.quiet)
            self.run.info('Kegg Modules', '%d found' % self.db.get_meta_value('num_modules'), quiet=self.quiet)

            days_since_created = self.get_days_since_creation()
            if not self.quiet and days_since_created >= KEGG_SETUP_INTERVAL:
                self.run.warning("Just a friendly PSA here: it has been at least %s days since the MODULES.db was created (%s days to be exact). "
                                 "It is entirely possible that KEGG has been updated since then, so perhaps it is a good idea to re-run "
                                 "anvi-setup-kegg-kofams to be sure that you are working with the latest KEGG data. No pressure, though. If you do "
                                 "want to reset your KEGG setup, we STRONGLY encourage saving a copy of your current KEGG data directory, just "
                                 "in case there was an update that breaks everything and you need to go back to your previous KEGG setup. Don't say we "
                                 "didn't warn you. And we will even be so nice as to tell you that your current KEGG data directory is %s"
                                 % (KEGG_SETUP_INTERVAL, days_since_created, self.kegg_data_dir))
        else:
            # if self.module_dict is None, then we tried to initialize the DB outside of setup
            if not self.module_dict:
                raise ConfigError("ERROR - a new KeggModulesDatabase() cannot be initialized without providing a modules dictionary. This "
                                  "usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-kofams` may fix this.")

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
        Lines for pathway module files can have between 2 and 4 fields, but in fact the only situation where there should be 4 lines is the ENTRY data,
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
        # so no matter the situation, data value is field 1 and data definition (if any) is field 2
        data_vals = fields[1]
        # need to sanity check data value field because SOME modules don't follow the 2-space separation formatting
        vals_are_okay, corrected_vals, corrected_def = self.data_vals_sanity_check(data_vals, current_data_name, current_module)

        if vals_are_okay and len(fields) > 2: # not all lines have a definition field
            data_def = fields[2]
        elif not vals_are_okay:
            data_vals = corrected_vals
            data_def = corrected_def

        # some types of information may need to be split into multiple db entries
        data_types_to_split = ["ORTHOLOGY","REACTION"] # lines that fall under these categories need to have data_vals split on comma
        if current_data_name in data_types_to_split:
            # here we should NOT split on any commas within parentheses
            for val in re.split(',(?!.*\))', data_vals):
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
                             "you would like to see which modules threw these errors, please re-run the setup using the --debug flag (you "
                             "will also probably need the --reset flag). When doing so, you will also see which lines caused issues; this "
                             "can be a lot of output, so you can suppress the line-specific output with the --quiet flag if that makes things "
                             "easier to read. So, in summary: You can probably ignore this warning. But if you want more info: run setup again "
                             "with --reset --debug --quiet to see exactly which modules had issues, or run with --reset --debug to see exactly "
                             "which lines in which modules had issues. Now, here is a kiss for you because you have been so patient and good with anvi'o 😚")
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

        self.db.disconnect()


    def get_days_since_creation(self):
        """Returns the time (in days) since MODULES.db was created"""
        return (time.time() - float(self.db.get_meta_value('creation_date'))) / 3600


    def get_db_content_hash(self):
        """Compute hash of all KOs and module numbers present in the db (used for tracking major changes to db content with future KEGG updates)"""
        mods = self.get_all_modules_as_list()
        mods.sort()
        orths = self.get_all_knums_as_list()
        orths.sort()
        mods_and_orths = mods + orths
        mods_and_orths = "".join(mods_and_orths)
        return str(hashlib.sha224(mods_and_orths.encode('utf-8')).hexdigest())[0:12]


    # KEGG Modules Table functions for data access and parsing start below
    # ====================================================================
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


    def unroll_module_definition(self, mnum):
        """This function accesses the DEFINITION line of a KEGG Module, unrolls it into all possible paths through the module, and
        returns the list of all paths.

        This is a driver for the recursive functions that do the actual unrolling of each definition line.
        """

        all_paths = [[]]
        def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        for d in def_lines:
            d = d.strip()
            def_line_paths = self.recursive_definition_unroller(d)
            new_paths_list = []
            for a in def_line_paths:
                for p in all_paths:
                    p_copy = copy.copy(p)
                    p_copy.extend(a)
                    new_paths_list.append(p_copy)
            all_paths = new_paths_list

        return all_paths


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
        Each step is either an atomic step (a single KO, module number '--', or nonessential KO starting with '-'),
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
