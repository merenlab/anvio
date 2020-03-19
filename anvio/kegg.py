#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains Kegg related classes."""

import os
import gzip
import shutil
import requests
import glob
import re

import anvio
import anvio.db as db
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.tables as t
import anvio.ccollections as ccollections

from anvio.errors import ConfigError, FilesNPathsError
from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions
from anvio.dbops import ContigsSuperclass, ContigsDatabase, ProfileSuperclass, ProfileDatabase


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
        self.quiet = A('quiet') or False
        self.just_do_it = A('just_do_it')

        # shared variables for all KEGG subclasses
        self.kofam_hmm_file_path = os.path.join(self.kegg_data_dir, "Kofam.hmm") # file containing concatenated KOfam hmms
        self.ko_list_file_path = os.path.join(self.kegg_data_dir, "ko_list")
        self.kegg_module_file = os.path.join(self.kegg_data_dir, "ko00002.keg")


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
                raise ConfigError("Hmm. Something is out of order. The orphan data directory %s does not exist \
                yet, but it needs to in order for the setup_ko_dict() function to work." % self.orphan_data_dir)
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
        skip_list  list of strings, each string is a KO number
        no_threshold_list   list of strings, each string is a KO number
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

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.kegg_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.module_data_dir, delete_if_exists=args.reset)

        # ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        self.files = ['ko_list.gz', 'profiles.tar.gz']

        # Kegg module text files
        self.kegg_module_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"


    def is_database_exists(self):
        """This function determines whether the user has already downloaded the Kofam HMM profiles and KEGG modules."""

        if os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("It seems you already have KOfam HMM profiles installed in '%s', please use --reset flag if you want to re-download it." % self.kegg_data_dir)

        if os.path.exists(self.kegg_module_file):
            raise ConfigError("Interestingly, though KOfam HMM profiles are not installed on your system, KEGG module information seems to have been \
            already downloaded in %s. Please use the --reset flag to re-download everything from scratch." % self.kegg_data_dir)

        if os.path.exists(self.module_data_dir):
            raise ConfigError("It seems the KEGG module directory %s already exists on your system. This is even more strange because Kofam HMM \
            profiles have not been downloaded. We suggest you to use the --reset flag to download everything from scratch." % self.module_data_dir)

    def download_profiles(self):
        """This function downloads the Kofam profiles."""

        self.run.info("Kofam Profile Database URL", self.database_url)

        for file_name in self.files:
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
                    raise ConfigError("While parsing the KEGG file %s, we found an unknown line code %s. This has \
                    made the file unparseable. Sad. :(" % (self.kegg_module_file, first_char))
        self.progress.end()

    def download_modules(self):
        """This function downloads the KEGG modules.

        To do so, it also processes the KEGG module file into a dictionary via the
        process_module_file() function. To verify that each file has been downloaded properly, we check that the last line is '///'.
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
                raise ConfigError("The KEGG module file %s was not downloaded properly. We were expecting the last line in the file \
                to be '///', but instead it was %s." % (file_path, last_line))


    def decompress_files(self):
        """This function decompresses the Kofam profiles."""

        for file_name in self.files:
            self.progress.new('Decompressing file %s' % file_name)
            full_path = os.path.join(self.kegg_data_dir, file_name)

            if full_path.endswith("tar.gz"):
                utils.tar_extract_file(full_path, output_file_path = self.kegg_data_dir, keep_original=False)
            else:
                utils.gzip_decompress_file(full_path, keep_original=False)

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
                    raise ConfigError("The KOfam HMM profile at %s does not exist. This probably means that something went wrong \
                                    while downloading the KOfam database. Please run `anvi-setup-kegg-kofams` with the --reset \
                                    flag." % (hmm_path))

    def move_orphan_files(self):
        """This function moves the following to the orphan files directory:

            - profiles that do not have ko_list entries
            - profiles whose ko_list entries have no scoring threshold (in ko_no_threshold_list)

        And, the following profiles should not have been downloaded, but we check if they exist and move any that do:
            - profiles whose ko_list entries have no data at all (in ko_skip_list)
        """

        if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
            raise ConfigError("Hmm. Something is out of order. The orphan data directory %s does not exist \
            yet, but it needs to in order for the move_orphan_files() function to work." % self.orphan_data_dir)

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
            self.run.warning("Please note that while anvi'o was building your databases, she found %d \
                            HMM profiles that did not have any matching KOfam entries. We have removed those HMM \
                            profiles from the final database. You can find them under the directory '%s'."
                            % (len(no_kofam_file_list), self.orphan_data_dir))

        if no_threshold_file_list:
            utils.concatenate_files(no_threshold_path, no_threshold_file_list, remove_concatenated_files=remove_old_files)
            self.run.warning("Please note that while anvi'o was building your databases, she found %d \
                            KOfam entries that did not have any threshold to remove weak hits. We have removed those HMM \
                            profiles from the final database. You can find them under the directory '%s'."
                            % (len(no_threshold_file_list), self.orphan_data_dir))
        if no_data_file_list:
            utils.concatenate_files(no_data_path, no_data_file_list, remove_concatenated_files=remove_old_files)
            self.run.warning("Please note that while anvi'o was building your databases, she found %d \
                            HMM profiles that did not have any associated data (besides an annotation) in their KOfam entries. \
                            We have removed those HMM profiles from the final database. You can find them under the directory '%s'."
                            % (len(no_data_file_list), self.orphan_data_dir))


    def run_hmmpress(self):
        """This function concatenates the Kofam profiles and runs hmmpress on them."""

        self.progress.new('Preparing Kofam HMM Profiles')
        log_file_path = os.path.join(self.kegg_data_dir, '00_hmmpress_log.txt')

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
        log_file_path = os.path.join(self.kegg_data_dir, '00_hmmpress_log.txt')
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Hmm. There was an error while running `hmmpress` on the Kofam HMM profiles. \
                                Check out the log file ('%s') to see what went wrong." % (log_file_path))
        else:
            # getting rid of the log file because hmmpress was successful
            os.remove(log_file_path)

        self.progress.end()

    def setup_modules_db(self):
        """This function creates the Modules DB from the Kegg Module files. """

        mod_db = KeggModulesDatabase(os.path.join(self.kegg_data_dir, "MODULES.db"), args=self.args, module_dictionary=self.module_dict, run=run, progress=progress)
        mod_db.create()


    def setup_profiles(self):
        """This is a driver function which executes the KEGG setup process by downloading, decompressing, and hmmpressing the profiles."""

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
        self.ko_dict = None # should be set up by setup_ko_dict()

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists('hmmscan')

        # verify that Kofam HMM profiles have been set up
        if not os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("Anvi'o is unable to find the Kofam.hmm file at %s. This can happen one of two ways. Either you \
                                didn't specify the correct KEGG data directory using the flag --kegg-data-dir, or you haven't \
                                yet set up the Kofam data by running `anvi-setup-kegg-kofams`. Hopefully you now know what to do \
                                to fix this problem. :) " % self.kegg_data_dir)

        utils.is_contigs_db(self.contigs_db_path)

        self.setup_ko_dict() # read the ko_list file into self.ko_dict

        # load existing kegg modules db
        self.kegg_modules_db = KeggModulesDatabase(os.path.join(self.kegg_data_dir, "MODULES.db"), args=self.args)

    def get_annotation_from_ko_dict(self, knum, ok_if_missing_from_dict=False):
        if not self.ko_dict:
            raise ConfigError("Oops! The ko_list file has not been properly loaded, so get_annotation_from_ko_dict() is \
                                extremely displeased and unable to function properly. Please refrain from calling this \
                                function until after setup_ko_dict() has been called.")

        if not knum in self.ko_dict:
            if ok_if_missing_from_dict:
                return "Unknown function with KO num %s" % knum
            else:
                raise ConfigError("It seems hmmscan found a KO number that does not exist\
                                   in the KOfam ko_list file: %s" % knum)

        return self.ko_dict[knum]['definition']

    def process_kofam_hmms(self):
        """This is a driver function for running HMMs against the KOfam database and processing the hits into the provided contigs DB"""

        tmp_directory_path = filesnpaths.get_temp_directory_path()
        contigs_db = ContigsSuperclass(self.args) # initialize contigs db

        # get AA sequences as FASTA
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['AA:GENE'],
                                                                   simple_headers=True,
                                                                   rna_alphabet=False,
                                                                   report_aa_sequences=True)

        # run hmmscan
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads)
        hmm_hits_file = hmmer.run_hmmscan('KOfam', 'AA', 'GENE', None, None, len(self.ko_dict), self.kofam_hmm_file_path, None, None)

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        if not hmm_hits_file:
            run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But\
                             now anvi'o will add KOfam as a functional source with no hits, clean the temporary directories\
                             and gracefully quit.", nl_before=1, nl_after=1)
            if not anvio.DEBUG:
                shutil.rmtree(tmp_directory_path)
                hmmer.clean_tmp_dirs()
            else:
                self.run.warning("Because you ran this script with the --debug flag, anvi'o will not clean up the temporary\
                directories located at %s and %s. Please be responsible for cleaning up this directory yourself \
                after you are finished debugging :)" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})
            return

        # parse hmmscan output
        parser = parser_modules['search']['hmmscan'](hmm_hits_file, alphabet='AA', context='GENE')
        search_results_dict = parser.get_search_results(ko_list_dict=self.ko_dict)

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
                mod_annotation = "\n".join(mods)
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
            self.run.warning("KOfam class has no hits to process. Returning empty handed, but still adding KOfam as \
                              a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})

        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up\
                         later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single('Cleaning up the temp directory (you can use `--debug` if you would\
                             like to keep it for testing purposes)', nl_before=1, nl_after=1)
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
        self.completeness_threshold = A('module-completion-threshold') or 0.75
        self.output_file_path = A('output_file') or "kegg-metabolism.txt"
        self.contigs_db_project_name = "Unknown"

        self.bin_ids_to_process = None
        if self.bin_id and self.bin_ids_file:
            raise ConfigError("You have provided anvi'o with both the individual bin id %s and a file with bin ids (%s). \
            Please make up your mind. Which one do you want an estimate for? :)" % (self.bin_id, self.bin_ids_file))
        elif self.bin_id:
            self.bin_ids_to_process = [self.bin_id]
        elif self.bin_ids_file:
            filesnpaths.is_file_exists(self.bin_ids_file)
            self.bin_ids_to_process = [line.strip() for line in open(self.bin_ids_file).readlines()]

        if self.profile_db_path and not self.collection_name:
            raise ConfigError("If you provide a profiles DB, you should also provide a collection name.")


        # init the base class
        KeggContext.__init__(self, self.args)

        utils.is_contigs_db(self.contigs_db_path)

        # load existing kegg modules db
        if not os.path.exists(os.path.join(self.kegg_data_dir, "MODULES.db")):
            raise ConfigError("It appears that a modules database (%s) does not exist in the KEGG data directory %s. \
            Perhaps you need to specify a different KEGG directory using --kegg-data-dir. Or perhaps you didn't run \
            `anvi-setup-kegg-kofams`, though we are not sure how you got to this point in that case \
            since you also cannot run `anvi-run-kegg-kofams` without first having run KEGG setup. But fine. Hopefully \
            you now know what you need to do to make this message go away." % ("MODULES.db", self.kegg_data_dir))
        self.kegg_modules_db = KeggModulesDatabase(os.path.join(self.kegg_data_dir, "MODULES.db"), args=self.args)

    def init_hits_and_splits(self):
        """This function loads splits and KOfam hits from the contigs DB.

        We will need the hits with their KO numbers (accessions) so that we can go through the MODULES.db and determine
        which steps are present in each module. And we will need the splits so that we can determine which hits belong
        to which genomes/bins when we are handling multiple of these. This function gets these hits and splits (as lists
        of tuples), and it makes sure that these lists don't include hits/splits we shouldn't be considering.

        RETURNS
        =======
        kofam_hits          list of (gene_call_id, ko_num) tuples
        genes_in_splits     list of (split, gene_call_id) tuples
        """

        self.progress.new('Loading')
        self.progress.update('Contigs DB')
        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        self.contigs_db_project_name = contigs_db.meta['project_name']
        self.progress.update('Splits')
        genes_in_splits = contigs_db.db.get_some_columns_from_table(t.genes_in_splits_table_name, "split, gene_callers_id")
        genes_in_contigs = contigs_db.db.get_some_columns_from_table(t.genes_in_contigs_table_name, "contig, gene_callers_id")
        self.progress.update('KOfam hits')
        kofam_hits = contigs_db.db.get_some_columns_from_table(t.gene_function_calls_table_name, "gene_callers_id, accession",
                                                where_clause="source = 'KOfam'")
        min_contig_length_in_contigs_db = contigs_db.db.get_max_value_in_column(t.contigs_info_table_name, "length", return_min_instead=True)
        contigs_db.disconnect()

        # get rid of gene calls in genes_in_splits that are not associated with KOfam hits.
        # Perhaps this is not a necessary step. But it makes me feel clean.
        all_gene_calls_in_splits = set([tpl[1] for tpl in genes_in_splits])
        gene_calls_with_kofam_hits = set([tpl[0] for tpl in kofam_hits])
        gene_calls_without_kofam_hits = all_gene_calls_in_splits.difference(gene_calls_with_kofam_hits)

        if gene_calls_without_kofam_hits:
            self.progress.update("Removing %s gene calls without KOfam hits" % len(gene_calls_without_kofam_hits))
            genes_in_splits = [tpl for tpl in genes_in_splits if tpl[1] not in gene_calls_without_kofam_hits]
            if anvio.DEBUG:
                self.run.warning("The following gene calls in your contigs DB were removed from consideration as they \
                do not have any hits to the KOfam database: %s" % (gene_calls_without_kofam_hits))

        # get rid of splits (and their associated gene calls) that are not in the profile DB
        if self.profile_db_path:
            split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.profile_db_path))
            split_names_in_contigs_db = set([tpl[0] for tpl in genes_in_splits])
            splits_missing_in_profile_db = split_names_in_contigs_db.difference(split_names_in_profile_db)

            min_contig_length_in_profile_db = ProfileDatabase(self.profile_db_path).meta['min_contig_length']

            if len(splits_missing_in_profile_db):
                self.progress.reset()
                self.run.warning("Please note that anvi'o found %s splits in your contigs database with KOfam hits. But only %s of them "
                                 "appear in the profile database. As a result, anvi'o will now remove the %s splits with KOfam hits"
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
                genes_in_splits = [tpl for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db]
                remaining_gene_calls = [tpl[1] for tpl in genes_in_splits]
                kofam_hits = [tpl for tpl in kofam_hits if tpl[0] in remaining_gene_calls]

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

        return kofam_hits, genes_in_splits

    def mark_kos_present_for_list_of_splits(self, kofam_hits_in_splits, split_list=None, bin_name=None):
        """This function generates a bin-level dictionary of dictionaries, which associates modules with the list of KOs
        that are present in the bin for each module.

        The structure of the dictionary is like this:
        {mnum: {present_kos: [knum1, knum2, ....]}}
        Why do we need an inner dictionary with just one list? Well. This dictionary will be expanded later by other functions, not to worry.

        PARAMETERS
        ==========
        kofam_hits_in_splits        list of KO numbers that are hits in the current list of splits
        split_list                  list of splits we are considering, this is only for debugging output
        bin_name                    name of the bin containing these splits, this is only for debugging output

        RETURNS
        =======
        bin_level_module_dict       dict of dicts that maps module number to dictionary of KOs present in the splits for that module
        """

        bin_level_module_dict = {}

        if anvio.DEBUG:
            self.run.info("Marking KOs present for bin", bin_name)
            self.run.info("Number of splits", len(split_list))

        # initialize all modules with empty presence list
        modules = self.kegg_modules_db.get_all_modules_as_list()
        for mnum in modules:
            bin_level_module_dict[mnum] = {"present_kos" : []}

        kos_not_in_modules = []
        for ko in kofam_hits_in_splits:
            present_in_mods = self.kegg_modules_db.get_modules_for_knum(ko)
            if not present_in_mods:
                kos_not_in_modules.append(ko)
            for m in present_in_mods:
                bin_level_module_dict[m]["present_kos"].append(ko)

        if anvio.DEBUG:
            self.run.info("KOs processed", "%d in bin" % len(kofam_hits_in_splits))
            if kos_not_in_modules:
                self.run.warning("Just so you know, the following KOfam hits did not belong to any KEGG modules in the MODULES.db: %s"
                % ", ".join(kos_not_in_modules))

        return bin_level_module_dict

    def compute_module_completeness_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the completeness of the specified module.

        This requires some parsing of the module DEFINITION fields. In these fields, we have the following:
        "Kxxxxx"    (KO numbers) indicating which enzyme contributes to a step in the module
        "Mxxxxx"    (module numbers) indicating that the module encompasses another module. This is rare. See note below.
        " "         (spaces) separating module steps; indicating an AND operation
        ","         (commas) separating alternatives (which can be singular KOs or entire pathways); indicating an OR operation
        "()"        (parentheses) enclosing comma-separated alternatives
        "+"         (plus sign) indicating the following KO is a necessary component of a complex; indicating an AND operation
        "-"         (minus sign) indicating the following KO is non-essential in a complex; so in other words we don't care if it is there

        What we will do is build a condition statement out of each step which will evaulate to True if the step can be considered present based
        on the available KOs in the current genome/bin.
        For example, suppose we have a step like: (K13937,((K00036,K19243) (K01057,K07404)))
        This will be parsed into the condition statement: (K13937 OR ((K00036 OR K19243) AND (K01057 OR K07404)))
        where the KOs will be replaced by True if they are present and False otherwise.

        While we are parsing, we save the individual module steps in lists (ie, one for all steps, one for complete steps) for easy access later.
        Afterwards we compute the completeness of the module based on the specified completion threshold.
        Then, we return a bunch of information about the completeness of the module, which can then be placed into the module completeness dictionary.

        There are 3 special cases to consider here.
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
        mnum                    string, module number to work on
        meta_dict_for_bin       metabolism completeness dict for the current bin, to be modified in-place

        VARIABLES FOR UPDATING METABOLISM COMPLETENESS DICT
        =======
        module_step_list                list of strings, each string is an individual step in the module (may have sub-steps if there are alternate pathways)
        module_complete_steps           list of strings, each string is a step in the module that is considered complete based on KO availability
        module_nonessential_steps       list of strings, each string is a step in the module that doesn't count for completeness estimates
        module_complete_nonessential_steps          list of strings, each string is a non-essential step that is considered complete based on KO availability
        module_total_steps                  int, the total number of steps in the module
        module_num_complete_steps           int, the number of complete steps in the module
        module_num_nonessential_steps       int, the total number of nonessential steps in the module
        module_num_complete_nonessential_steps      int, the number of nonessential steps in the module that were found to be complete
        module_completeness             float, a decimal indicating the fraction of complete steps in the module

        RETURNS
        =======
        over_complete_threshold         boolean, whether or not the module is considered "complete" overall based on the threshold fraction of completeness
        has_nonessential_step           boolean, whether or not the module contains non-essential steps. Used for warning the user about these.
        has_no_ko_step                  boolean, whether or not the module contains steps without associated KOs. Used for warning the user about these.
        defined_by_modules              boolean, whether or not the module contains steps defined by other modules. Used for going back to adjust completeness later.
        """

        present_list_for_mnum = meta_dict_for_bin[mnum]["present_kos"]
        if not present_list_for_mnum:
            # no KOs in this module are present
            if anvio.DEBUG:
                self.run.warning("No KOs present for module %s. Parsing for completeness is still being done to obtain module steps." % mnum)

        # module information to return
        module_step_list = [] # while we are at it, we'll remember what the (essential) steps are
        module_complete_steps = [] # and what the complete steps are
        module_nonessential_steps = [] # steps that aren't necessary for module completeness
        module_complete_nonessential_steps = [] # and those nonessential steps which we find are complete
        module_total_steps = 0
        module_num_complete_steps = 0
        module_num_nonessential_steps = 0
        module_num_complete_nonessential_steps = 0
        has_nonessential_step = False
        has_no_ko_step = False
        defined_by_modules = False

        def_lines = self.kegg_modules_db.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        for d in def_lines:
            d = d.strip()
            cur_index = 0  # current position in the DEFINITION line
            parens_level = 0 # how deep we are in nested parentheses
            step_is_present_condition_statement = ""
            last_step_end_index = 0

            while cur_index < len(d):
                if d[cur_index] == "K": # we have found a KO
                    ko = d[cur_index:cur_index+6]
                    defined_by_modules = False  # reset this flag just in case KO-defined step comes after a module-defined step
                    if ko in present_list_for_mnum:
                        step_is_present_condition_statement += "True"
                    else:
                        step_is_present_condition_statement += "False"
                    cur_index += 6

                elif d[cur_index] == "(":
                    parens_level += 1
                    step_is_present_condition_statement += "("
                    cur_index += 1

                elif d[cur_index] == ")":
                    parens_level -= 1
                    step_is_present_condition_statement += ")"
                    cur_index += 1

                elif d[cur_index] == ",":
                    step_is_present_condition_statement += " or "
                    cur_index += 1

                elif d[cur_index] == "+":
                    step_is_present_condition_statement += " and "
                    cur_index += 1

                elif d[cur_index] == "-":
                    # either a singular KO or a set of KOs in parentheses can follow this character
                    # since the following KO(s) are non-essential in the complex, we skip over them to ignore them
                    # unless this is its own step, in which case we consider the whole step non-essential

                    # singular nonessential KO
                    if d[cur_index+1] == "K":
                        nonessential_ko = d[cur_index+1:cur_index+7]
                        cur_index += 7
                        """
                        OKAY, SO HERE WE HAVE SOME POOPINESS THAT MAY NEED TO BE FIXED EVENTUALLY.
                        Basically, some DEFINITION lines have KOs that seem to be marked non-essential;
                        ie, "-K11024" in "K11023 -K11024 K11025 K11026 K11027".
                        It was difficult to decide whether we should consider only K11024, or K11024 and all following KOs, to be non-essential.
                        For instance, the module M00778 is a complex case that gave us pause - see Fiesta issue 955.
                        But for now, we have decided to just track only the one KO as a 'non-essential step', and to not include such steps in
                        the module completeness estimate.
                        """
                        # if this is the first KO in the step and we find a space after this KO, then we have found a non-essential step
                        if step_is_present_condition_statement == "" and (cur_index == len(d) or d[cur_index] == " "):
                            has_nonessential_step = True
                            module_nonessential_steps.append(d[last_step_end_index:cur_index])
                            module_num_nonessential_steps += 1

                            if nonessential_ko in present_list_for_mnum:
                                module_complete_nonessential_steps.append(d[last_step_end_index:cur_index])
                                module_num_complete_nonessential_steps += 1

                            # reset for next step
                            last_step_end_index = cur_index + 1
                            cur_index += 1

                    # a whole set of nonessential KOs
                    elif d[cur_index+1] == "(":
                        while d[cur_index] != ")":
                            cur_index += 1
                        cur_index += 1 # skip over the ')'

                    # the '--' (no KO) situation
                    elif d[cur_index+1] == "-":
                        # when '--' in a DEFINITION line happens, it signifies a reaction step that has no associated KO.
                        # we assume that such steps are not complete,  because we really can't know if it is from the KOfam hits alone
                        has_no_ko_step = True
                        step_is_present_condition_statement += "False"
                        cur_index += 2 # skip over both '-', the next character should be a space or end of DEFINITION line

                        if cur_index < len(d) and d[cur_index] != " ":
                            raise ConfigError("Serious, serious parsing sadness is happening. We just processed a '--' in "
                                              "a DEFINITION line for module %s, but did not see a space afterwards. Instead, we found %s. "
                                              "WHAT DO WE DO NOW?" % (mnum, d[cur_index+1]))
                    # anything else that follows a '-'
                    else:
                        raise ConfigError("The following character follows a '-' in the DEFINITION line for module %s "
                        "and we just don't know what to do: %s" % (mnum, d[cur_index+1]))

                elif d[cur_index] == " ":
                    # if we are outside of parentheses, we are done processing the current step
                    if parens_level == 0:
                        module_step_list.append(d[last_step_end_index:cur_index])
                        module_total_steps += 1
                        # we do not evaluate completeness of this step yet if it is defined by other modules
                        if not defined_by_modules:
                            step_is_present = eval(step_is_present_condition_statement)
                            if step_is_present:
                                module_complete_steps.append(d[last_step_end_index:cur_index])
                                module_num_complete_steps += 1
                        # reset for next step
                        step_is_present_condition_statement = ""
                        last_step_end_index = cur_index + 1
                        cur_index += 1
                    # otherwise, we are processing an alternative path so AND is required
                    else:
                        step_is_present_condition_statement += " and "
                        cur_index += 1

                elif d[cur_index] == "M":
                    """
                    This happens when a module is defined by other modules. For example, photosynthesis module M00611 is defined as
                    (M00161,M00163) M00165 === (photosystem II or photosystem I) and calvin cycle

                    We need all the modules to have been evaluated before we can determine completeness of steps with module numbers.
                    So what we will do here is just add the step to the appropriate lists without evaluating completeness, and use a
                    flag variable to keep track of the modules that have this sort of definition in a list so we can go back and
                    evaluate completeness of steps with module numbers later.
                    """
                    defined_by_modules = True
                    cur_index += 6

                else:
                    raise ConfigError("While parsing the DEFINITION field for module %s, (which is %s), anvi'o found the following character "
                                        "that she didn't understand: %s. Unfortunately, this means we cannot determine the module "
                                        "completeness. For context, here is the current index in the DEFINITION line: %s and the "
                                        "surrounding characters: %s" % (mnum, d, d[cur_index], cur_index, d[cur_index-5:cur_index+6]))

            # once we have processed the whole line, we still need to eval the last step.
            # Unless we already did (this can happen with non-essential steps), which we check by seeing if the condition statement is empty
            # However, if this step is defined by modules, the condition statement will be empty, but we still need to save the step
            if step_is_present_condition_statement != "" or defined_by_modules:
                module_step_list.append(d[last_step_end_index:cur_index])
                module_total_steps += 1
                if not defined_by_modules:
                    step_is_present = eval(step_is_present_condition_statement)
                    if step_is_present:
                        module_complete_steps.append(d[last_step_end_index:cur_index])
                        module_num_complete_steps += 1

        # once we have processed all DEFINITION lines, we can compute the overall completeness
        module_completeness = module_num_complete_steps / module_total_steps
        over_complete_threshold = True if module_completeness >= self.completeness_threshold else False

        # instead of returning everything, we update the metabolism completeness dictionary in place
        meta_dict_for_bin[mnum]["step_list"] = module_step_list
        meta_dict_for_bin[mnum]["complete_step_list"] = module_complete_steps
        meta_dict_for_bin[mnum]["nonessential_step_list"] = module_nonessential_steps
        meta_dict_for_bin[mnum]["complete_nonessential_step_list"]= module_complete_nonessential_steps
        meta_dict_for_bin[mnum]["num_steps"] = module_total_steps
        meta_dict_for_bin[mnum]["num_complete_steps"] = module_num_complete_steps
        meta_dict_for_bin[mnum]["num_nonessential_steps"] = module_num_nonessential_steps
        meta_dict_for_bin[mnum]["num_complete_nonessential_steps"] = module_num_complete_nonessential_steps
        meta_dict_for_bin[mnum]["percent_complete"] = module_completeness
        meta_dict_for_bin[mnum]["complete"] = over_complete_threshold
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
        mod                 string, the module number to adjust
        meta_dict_for_bin   metabolism completeness dictionary for the current bin

        RETURNS
        =======
        now_complete        boolean, whether or not the module is NOW considered "complete" overall based on the threshold fraction of completeness
        """

        for step in meta_dict_for_bin[mod]["step_list"]:
            cur_index = 0  # current position in the step definition
            parens_level = 0 # how deep we are in nested parentheses
            step_is_present_condition_statement = ""
            is_ko_step = False
            while cur_index < len(step):
                # we have found a KO so we can ignore this step; it has already been counted as complete or not
                if step[cur_index] == "K":
                    is_ko_step = True
                    break

                # we have found a module so we must evaluate this steps's completeness by checking if the module is complete
                elif step[cur_index] == "M":
                    mnum = step[cur_index:cur_index+6]
                    if meta_dict_for_bin[mnum]["complete"]:
                        step_is_present_condition_statement += "True"
                    else:
                        step_is_present_condition_statement += "False"
                    cur_index += 6

                elif step[cur_index] == "(":
                    parens_level += 1
                    step_is_present_condition_statement += "("
                    cur_index += 1

                elif step[cur_index] == ")":
                    parens_level -= 1
                    step_is_present_condition_statement += ")"
                    cur_index += 1

                elif step[cur_index] == ",":
                    step_is_present_condition_statement += " or "
                    cur_index += 1

                elif step[cur_index] == " ":
                    # if we are outside of parentheses, something is wrong because this should all be just one step
                    if parens_level == 0:
                        raise ConfigError("Much parsing sadness. We thought we were re-evaluating the completeness of just one step in "
                                          "module %s (step: %s), but we found a space that seems to indicate another step. HALP." % (mod, step))
                    # otherwise, we are processing an alternative path so AND is required
                    else:
                        step_is_present_condition_statement += " and "
                        cur_index += 1

                else:
                    raise ConfigError("While correcting completeness for module %s, (step %s), anvi'o found the following character "
                                        "that she didn't understand: %s. Unfortunately, this means we cannot determine the module "
                                        "completeness. For context, here is the current index in the DEFINITION line: %s and the "
                                        "surrounding characters: %s" % (mod, step, step[cur_index], cur_index, step[cur_index-5:cur_index+6]))
            # once we have processed everything, we need to re-evaluate the step (provided its not a KO step that has already been evaluated)
            if not is_ko_step:
                step_is_present = eval(step_is_present_condition_statement)
                if step_is_present:
                    meta_dict_for_bin[mod]["complete_step_list"].append(step)
                    meta_dict_for_bin[mod]["num_complete_steps"] += 1

        # now, we recalculate module completeness
        meta_dict_for_bin[mod]["percent_complete"] = meta_dict_for_bin[mod]["num_complete_steps"] / meta_dict_for_bin[mod]["num_steps"]
        now_complete = True if meta_dict_for_bin[mod]["percent_complete"] >= self.completeness_threshold else False
        meta_dict_for_bin[mod]["complete"] = now_complete
        if now_complete:
            meta_dict_for_bin["num_complete_modules"] += 1

        return now_complete

    def estimate_for_list_of_splits(self, ko_hits_in_splits, splits=None, bin_name=None):
        """This is the atomic metabolism estimator function, which builds a metabolism completeness dictionary for an arbitrary list of splits.

        For example, the list of splits may represent a bin or a single isolate genome.
        The metabolism completeness dictionary is first initialized to contain the KOs that are present in the genome for each KEGG module.
        It is later updated with the individual steps and completion estimates for each module.

        PARAMETERS
        ==========
        ko_hits_in_splits       a list of KO numbers indicating the KOfam hits that have occurred in this list of splits
        splits                  a list of splits identifiers
        bin_name                the name of the bin that we are working with

        RETURNS
        =======
        metabolism_dict_for_list_of_splits      the metabolism completeness dictionary of dictionaries for this list of splits. It contains
                                                one dictionary of module steps and completion information for each module (keyed by module number),
                                                as well as one key num_complete_modules that tracks the number of complete modules found in these splits.
                                                Calling functions should assign this dictionary to a metabolism superdict with the bin name as a key.
        """

        metabolism_dict_for_list_of_splits = self.mark_kos_present_for_list_of_splits(ko_hits_in_splits, split_list=splits,
                                                                                                    bin_name=bin_name)
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


    def estimate_for_genome(self, kofam_hits, genes_in_splits):
        """This is the metabolism estimation function for a contigs DB that contains a single genome.

        Assuming this contigs DB contains only one genome, it sends all of the splits and their kofam hits to the atomic
        estimation function for processing. It then returns the metabolism completion dictionary for the genome, wrapped in the superdict format.

        PARAMETERS
        ==========
        kofam_hits          list of (gene_call_id, ko_num) tuples, all belong to this single genome
        genes_in_splits     list of (split, gene_call_id) tuples, all belong to this single genome

        RETURNS
        =======
        genome_metabolism_dict      dictionary mapping genome name to its metabolism completeness dictionary
        """

        genome_metabolism_superdict = {}
        # get list of KOs only - since all splits belong to one genome, we can take all the hits
        ko_in_genome = [tpl[1] for tpl in kofam_hits]
        splits_in_genome = [tpl[0] for tpl in genes_in_splits]

        genome_metabolism_superdict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(ko_in_genome, splits=splits_in_genome, bin_name=self.contigs_db_project_name)

        return genome_metabolism_superdict


    def estimate_for_bins_in_collection(self, kofam_hits, genes_in_splits):
        bins_metabolism_superdict = {}

        bin_name_to_split_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.run.info_single("%s split names associated with %s bins of in collection '%s' have been "
                             "successfully recovered 🎊" % (pp(sum([len(v) for v in bin_name_to_split_names_dict.values()])),
                                                           pp(len(bin_name_to_split_names_dict)),
                                                           self.collection_name), nl_before=1)

        for bin_name in bin_name_to_split_names_dict:
            splits_in_bin = bin_name_to_split_names_dict[bin_name]
            genes_in_bin = [tpl[1] for tpl in genes_in_splits if tpl[0] in splits_in_bin]
            ko_in_bin = [tpl[1] for tpl in kofam_hits if tpl[0] in genes_in_bin]
            bins_metabolism_superdict[bin_name] = self.estimate_for_list_of_splits(ko_in_bin, splits=splits_in_bin, bin_name=bin_name)

        return bins_metabolism_superdict


    def estimate_metabolism(self):
        """This is the driver function for estimating metabolism.

        It will decide what to do based on whether the input contigs DB is a genome or metagenome.
        It returns the metabolism superdict which contains a metabolism completion dictionary for each genome/bin in the contigs db.
        The metabolism completion dictionary is keyed by KEGG module number.
        """

        hits_to_consider, splits_to_consider = self.init_hits_and_splits()

        kegg_metabolism_superdict = {}

        if self.profile_db_path and not self.metagenome_mode:
            kegg_metabolism_superdict = self.estimate_for_bins_in_collection(hits_to_consider, splits_to_consider)
        elif not self.profile_db_path and not self.metagenome_mode:
            kegg_metabolism_superdict = self.estimate_for_genome(hits_to_consider, splits_to_consider)
        elif self.profile_db_path and self.metagenome_mode:
            raise ConfigError("This class doesn't know how to deal with that yet :/")
            # metagenome, with profiling
            #self.estimate_for_contigs_db_for_metagenome()
        elif not self.profile_db_path and self.metagenome_mode:
            raise ConfigError("This class doesn't know how to deal with that yet :/")
            # metagenome without profiling
            #self.estimate_for_contigs_db_for_metagenome()
        else:
            raise ConfigError("This class doesn't know how to deal with that yet :/")

        self.store_kegg_metabolism_superdict(kegg_metabolism_superdict)


    def store_kegg_metabolism_superdict(self, kegg_superdict):
        """This function writes the metabolism superdict to a tab-delimited file.

        The metabolism superdict is a three-level dictionary (genomes/bins, modules, and module completion information).
        To distill this information into one line, we need to convert the dictionary on-the-fly to a dict of dicts,
        where each genome/bin-module pair is keyed by an arbitrary integer.
        """

        d = {}
        i = 0
        for bin, mod_dict in kegg_superdict.items():
            for mnum, c_dict in mod_dict.items():
                if mnum == "num_complete_modules":
                    continue
                d[i] = c_dict
                d[i]["bin_name"] = bin
                d[i]["kegg_module"] = mnum
                i += 1

        utils.store_dict_as_TAB_delimited_file(d, self.output_file_path, key_header="unique_id")
        self.run.info("Output file", self.output_file_path, nl_before=1)


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
        else:
            # if self.module_dict is None, then we tried to initialize the DB outside of setup
            if not self.module_dict:
                raise ConfigError("ERROR - a new KeggModulesDatabase() cannot be initialized without providing a modules dictionary. This \
                usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-kofams` may fix this.")

    def touch(self):
        """Creates an empty Modules database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db.
        """

        # sanity check to avoid overriding previous Modules DB
        # this will probably never happen as long as this function is called through the setup script, but we check just in case
        if os.path.exists(self.db_path):
            raise ConfigError("A modules database at %s already exists. Please use the --reset flag when you restart the setup \
            if you really want to get rid of this one and make a new one." % (self.db_path))


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
        data_vals           str, the data values field (split from the kegg module line)
        current_data_name   str, which data name we are working on. It should never be None because we should have already figured this out by parsing the line.
        current_module_num  str, which module we are working on. We need this to keep track of which modules throw parsing errors.

        RETURNS
        =======
        is_ok               bool, whether the values look correctly formatted or not
        """

        is_ok = True
        is_corrected = False
        corrected_vals = None
        corrected_def = None

        if not current_data_name:
            raise ConfigError("data_vals_sanity_check() cannot be performed when the current data name is None. Something was not right when parsing the KEGG \
            module line.")
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
                self.run.warning("While parsing, anvi'o found an uncorrectable issue with a KEGG Module line in module %s, but since you used the --just-do-it flag, \
                anvi'o will quietly ignore this issue and add the line to the MODULES.db anyway. Please be warned that this may break things downstream. \
                In case you are interested, the line causing this issue has data name %s and data value %s" % (current_module_num, current_data_name, data_vals))
                is_ok = True # let's pretend that everything is alright so that the next function will take the original parsed values
            else:
                raise ConfigError("While parsing, anvi'o found an uncorrectable issue with a KEGG Module line in module %s. The current data name is %s, \
                here is the incorrectly-formatted data value field: %s. If you think this is totally fine and want to ignore errors like this, please \
                re-run the setup with the --just-do-it flag. But if you choose to do that of course we are obliged to inform you that things may eventually \
                break as a result." % (current_module_num, current_data_name, data_vals))

        if is_corrected:
            self.num_corrected_errors += 1
            if anvio.DEBUG and not self.quiet:
                self.run.warning("While parsing a KEGG Module line, we found an issue with the formatting. We did our very best to parse the line \
                correctly, but please check that it looks right to you by examining the following values.")
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
        line                 str, the line to parse
        current_module       str, which module we are working on. We need this to keep track of which modules throw parsing errors
        line_num             int, which line number we are working on. We need this to keep track of which entities come from the same line of the file.
        current_data_name    str, which data name we are working on. If this is None, we need to parse this info from the first field in the line.

        RETURNS
        =======
        line_entries        a list of tuples, each containing information for one db entry, namely data name, data value, data definition, and line number.
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
                raise ConfigError("Oh, please. Some silly developer (you know who you are) has tried to call parse_kegg_modules_line() on \
                a line without a data name field, and forgot to give it the current data name. Shame on you, go fix this. (For reference here \
                is the line: %s)" % (line))

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
            raise ConfigError("Appparently, the Kegg Modules were not correctly setup and now all sorts of things are broken. The \
             Modules DB cannot be created from broken things. BTW, this error is not supposed to happen to anyone except maybe developers, so \
             if you do not fall into that category you are likely in deep doo-doo. Maybe re-running setup with --reset will work? (if not, you \
             probably should email/Slack/telepathically cry out for help to the developers). Anyway, if this helps make things any clearer, \
             the number of modules in the module dictionary is currently %s" % len(self.module_dict.keys()))

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
                        # append_and_store will collect db entries and store every 10000 at a time
                        mod_table.append_and_store(self.db, mnum, name, val, definition, line)

                f.close()

            num_modules_parsed += 1
        # once we are done parsing all modules, we store whatever db entries remain in the db_entries list
        # this is necessary because append_and_store() above only stores every 10000 entries
        self.progress.update("Storing final batch of module entries into DB")
        mod_table.store(self.db)

        self.progress.end()

        # warn user about parsing errors
        if anvio.DEBUG:
            self.run.warning("Several parsing errors were encountered while building the KEGG Modules DB. \
            Below you will see which modules threw each type of parsing error. Note that modules which threw multiple \
            errors will occur in the list as many times as it threw each error.")
            self.run.info("Bad line splitting (usually due to rogue or missing spaces)", self.parsing_error_dict["bad_line_splitting"])
            self.run.info("Bad KEGG code format (not corrected; possibly problematic)", self.parsing_error_dict["bad_kegg_code_format"])
        else: # less verbose
            self.run.warning("First things first - don't panic. Several parsing errors were encountered while building the KEGG Modules DB. But that \
            is probably okay, because if you got to this point it is likely that we already fixed all of them ourselves. So don't worry too much. \
            Below you will see how many of each type of error was encountered. If you would like to see which modules threw these errors, please \
            re-run the setup using the --debug flag (you will also probably need the --reset flag). When doing so, you will also see which lines \
            caused issues; this can be a lot of output, so you can suppress the line-specific output with the --quiet flag if that makes things easier to read. \
            So, in summary: You can probably ignore this warning. But if you want more info: \
            run setup again with --reset --debug --quiet to see exactly which modules had issues, or \
            run --reset --debug to see exactly which lines in which modules had issues. \
            Now, here is a kiss for you because you have been so patient and good with anvi'o 😚")
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

        self.db.disconnect()


    # KEGG Modules Table functions for data access and parsing start below
    # ====================================================================
    def get_data_value_entries_for_module_by_data_name(self, module_num, data_name):
        """This function returns data_value elements from the modules table for the specified module and data_name pair.

        All elements corresponding to the pair (ie, M00001 and ORTHOLOGY) will be returned.
        The function relies on the db.get_some_rows_from_table_as_dict() function to first fetch all rows corresponding \
        to a particular model, and then parses the resulting dictionary to find all the elements with the given data_name field.

        PARAMETERS
        ==========
        module_num           str, the module to fetch data for
        data_name            str, which data_name field we want

        RETURNS
        =======
        data_values_to_ret   list of str, the data_values corresponding to the module/data_name pair
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
            self.run.warning("Just so you know, we tried to fetch data from the KEGG Modules database for the data_name field %s and KEGG module %s, \
            but didn't come up with anything, so an empty list is being returned. This may cause errors down the line, and if so we're very sorry for that.")

        return data_values_to_ret

    def get_all_modules_as_list(self):
        """This function returns a list of all modules in the DB."""
        return self.db.get_single_column_from_table(self.module_table_name, 'module', unique=True)

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
        where_clause_string = "module = '%s'" % (mnum)
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
