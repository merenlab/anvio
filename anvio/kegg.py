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

from anvio.errors import ConfigError, FilesNPathsError
from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions

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
    """The purpose of this base class is to define shared functions and file paths for all KOfam operations."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # default directory will be called KEGG and will store the KEGG Module data as well
        self.kofam_data_dir = A('kofam_data_dir') or os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')
        self.orphan_data_dir = os.path.join(self.kofam_data_dir, "orphan_data")
        self.module_data_dir = os.path.join(self.kofam_data_dir, "modules")

        # shared variables for all KOfam subclasses
        self.kofam_hmm_file_path = os.path.join(self.kofam_data_dir, "Kofam.hmm") # file containing concatenated KOfam hmms
        self.ko_list_file_path = os.path.join(self.kofam_data_dir, "ko_list")
        self.kegg_module_file = os.path.join(self.kofam_data_dir, "ko00002.keg")


    def setup_ko_dict(self):
        """The purpose of this function is to process the ko_list file into usable form by Kofam sub-classes.

        The ko_list file (which is downloaded along with the KOfam HMM profiles) contains important
        information for each KEGG Orthology number (KO, or knum), incuding pre-defined scoring thresholds
        for limiting HMM hits and annotation information.

        It looks something like this:

        knum    threshold    score_type    profile_type    F-measure    nseq    nseq_used    alen    mlen    eff_nseq    re/pos    definition
        K00001    329.57    domain    trim    0.231663    1473    1069    1798    371    17.12    0.590    alcohol dehydrogenase [EC:1.1.1.1]

        Since this information is useful for both the setup process (we need to know all the knums) and HMM process,
        all Kofam subclasses need to have access to this dictionary.

        This is a dictionary (indexed by knum) of dictionaries(indexed by column name).
        Here is an example of the dictionary structure:
        self.ko_dict["K00001"]["threshold"] = 329.57
        """

        self.ko_dict = utils.get_TAB_delimited_file_as_dictionary(self.ko_list_file_path)
        self.ko_skip_list, self.ko_no_threshold_list = self.get_ko_skip_list()

        # if we are currently setting up KOfams, we should generate a text file with the ko_list entries
        # of the KOs that have no scoring threshold
        if self.__class__.__name__ in ['KofamSetup']:
            orphan_ko_dict = {ko:self.ko_dict[ko] for ko in self.ko_skip_list}
            orphan_ko_dict.update({ko:self.ko_dict[ko] for ko in self.ko_no_threshold_list})

            if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
                raise ConfigError("Hmm. Something is out of order. The orphan data directory %s does not exist \
                yet, but it needs to in order for the setup_ko_dict() function to work." % self.orphan_data_dir)
            orphan_ko_path = os.path.join(self.orphan_data_dir, "01_ko_fams_with_no_threshold.txt")
            orphan_ko_headers = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos", "definition"]
            utils.store_dict_as_TAB_delimited_file(orphan_ko_dict, orphan_ko_path, key_header="knum", headers=orphan_ko_headers)

        # here we remove KOs from the dictionary if they are in the skip list or no threshold list
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
    """Class for setting up KEGG Kofam HMM profiles. It performs sanity checks and downloads, unpacks, and prepares the profiles for later use by `hmmscan`.

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

        filesnpaths.gen_output_directory(self.kofam_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.module_data_dir, delete_if_exists=args.reset)

        # ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        self.files = ['ko_list.gz', 'profiles.tar.gz']

        # Kegg module text file
        self.kegg_module_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"


    def is_database_exists(self):
        """This function determines whether the user has already downloaded the Kofam HMM profiles and KEGG modules."""

        if os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("It seems you already have KOfam HMM profiles installed in '%s', please use --reset flag if you want to re-download it." % self.kofam_data_dir)

        if os.path.exists(self.kegg_module_file):
            raise ConfigError("Interestingly, though KOfam HMM profiles are not installed on your system, KEGG module information seems to have been \
            already downloaded in %s. Please use the --reset flag to re-download everything from scratch." % self.kofam_data_dir)

        if os.path.exists(self.module_data_dir):
            raise ConfigError("It seems the KEGG module directory %s already exists on your system. This is even more strange because Kofam HMM \
            profiles have not been downloaded. We suggest you to use the --reset flag to download everything from scratch." % self.module_data_dir)

    def download_profiles(self):
        """This function downloads the Kofam profiles."""

        self.run.info("Kofam Profile Database URL", self.database_url)

        for file_name in self.files:
            utils.download_file(self.database_url + '/' + file_name,
                os.path.join(self.kofam_data_dir, file_name), progress=self.progress, run=self.run)

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
            full_path = os.path.join(self.kofam_data_dir, file_name)

            if full_path.endswith("tar.gz"): # extract tar file instead of doing gzip
                utils.tar_extract_file(full_path, output_file_path = self.kofam_data_dir, keep_original=False)
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
                hmm_path = os.path.join(self.kofam_data_dir, "profiles/%s.hmm" % k)
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

        hmm_list = [k for k in glob.glob(os.path.join(self.kofam_data_dir, 'profiles/*.hmm'))]
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
        log_file_path = os.path.join(self.kofam_data_dir, '00_hmmpress_log.txt')

        self.progress.update('Verifying the Kofam directory %s contains all HMM profiles' % self.kofam_data_dir)
        self.confirm_downloaded_profiles()

        self.progress.update('Handling orphan files')
        self.move_orphan_files()

        self.progress.update('Concatenating HMM profiles into one file...')
        hmm_list = [k for k in glob.glob(os.path.join(self.kofam_data_dir, 'profiles/*.hmm'))]
        utils.concatenate_files(self.kofam_hmm_file_path, hmm_list, remove_concatenated_files=False)

        # there is no reason to keep the original HMM profiles around, unless we are debugging
        if not anvio.DEBUG:
            shutil.rmtree((os.path.join(self.kofam_data_dir, "profiles")))

        self.progress.update('Running hmmpress...')
        cmd_line = ['hmmpress', self.kofam_hmm_file_path]
        log_file_path = os.path.join(self.kofam_data_dir, '00_hmmpress_log.txt')
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

        mod_db = KeggModulesDatabase(os.path.join(self.kofam_data_dir, "MODULES.db"), module_dictionary=self.module_dict, args=self.args, run=run, progress=progress)
        mod_db.create()


    def setup_profiles(self):
        """This is a driver function which executes the Kofam setup process by downloading, decompressing, and hmmpressing the profiles."""

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
                                didn't specify the correct Kofam data directory using the flag --kofam-data-dir, or you haven't \
                                yet set up the Kofam data by running `anvi-setup-kegg-kofams`. Hopefully you now know what to do \
                                to fix this problem. :) " % self.kofam_data_dir)

        utils.is_contigs_db(self.contigs_db_path)

        self.setup_ko_dict() # read the ko_list file into self.ko_dict

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
        contigs_db = dbops.ContigsSuperclass(self.args) # initialize contigs db

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

        # add functions to database
        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'KOfam',
                'accession': hmm_hit['gene_name'],
                'function': self.get_annotation_from_ko_dict(hmm_hit['gene_name'], ok_if_missing_from_dict=True),
                'e_value': hmm_hit['e_value'],
            }

            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
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

class KeggModulesDatabase(KeggContext):
    """To create or access a Modules DB.

    This DB should be created in the Kegg Data folder during Kofam setup, and will be populated with information from the
    Kegg Module files.
    """

    def __init__(self, db_path, module_dictionary, args, run=run, progress=progress, quiet=False):
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

        ## here we should call init function if the db exists

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

    def data_vals_sanity_check(self, data_vals, current_data_name):
        """This function checks if the data values were correctly parsed from a line in a KEGG module file.

        This is a sadly necessary step because some KEGG module file lines are problematic and don't follow the right format (ie, 2+ spaces
        between different fields). So here we check if the values that we parsed look like they are the right format, without any extra bits.
        Each data name (ORTHOLOGY, DEFINITION, etc) has a different format to check for.

        Note that we don't check the following data name types: NAME, CLASS, REFERENCE

        PARAMETERS
        ==========
        data_vals           str, the data values field (split from the kegg module line)
        current_data_name   str, which data name we are working on. It should never be None because we should have already figured this out by parsing the line.

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
        elif current_data_name == "DEFINITION":
            # example format: (K01647,K05942) (K01681,K01682) (K00031,K00030) (K00164+K00658+K00382,K00174+K00175-K00177-K00176)
            # another example: (M00161,M00163) M00165
            knums = [x for x in re.split('\(|\)|,| |\+|-',data_vals) if x]
            for k in knums:
                if k[0] not in ['K','M'] or len(k) != 6:
                    is_ok = False
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
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                is_corrected = True
        elif current_data_name == "COMPOUND":
            # example format: C00024
            if data_vals[0] not in ['C','G'] or len(data_vals) != 6:
                is_ok = False
        elif current_data_name == "RMODULE":
            # example format: RM003
            if data_vals[0:2] != "RM" or len(data_vals) != 5:
                is_ok = False


        if not is_ok and not is_corrected:
            # in production, this should not end with an error. This raises an error for now just so I can easily find errors that I haven't implemented
            # correction for yet
            raise ConfigError("Found an issue with a KEGG Module line. Data values incorrectly parsed. Current data name is %s, here is the \
            incorrectly-formatted data value field: %s" % (current_data_name, data_vals))

        if is_corrected:
            self.run.warning("While parsing a KEGG Module line, we found an issue with the formatting. We did our very best to parse the line \
            correctly, but please check that it looks right to you by examining the following values.")
            self.run.info("Incorrectly parsed data value field", data_vals)
            self.run.info("Corrected data values", corrected_vals)
            self.run.info("Corrected data definition", corrected_def)

        return is_ok, corrected_vals, corrected_def


    def parse_kegg_modules_line(self, line, line_num = None, current_data_name=None):
        """This function parses information from one line of a KEGG module file.

        These files have fields separated by 2 or more spaces. Fields can include data name (not always), data value (always), and data definition (not always).
        Lines for pathway module files can have between 2 and 4 fields, but in fact the only situation where there should be 4 lines is the ENTRY data,
        which for some inexplicable reason has multiple spaces between "Pathway" and "Module" in the data definition field. We can safely ignore this last "Module", I think.

        Some lines will have multiple entities in the data_value field (ie, multiple KOs or reaction numbers) and will be split into multiple db entries.

        PARAMETERS
        ==========
        line                 str, the line to parse
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

        # data name unknown, parse from first field
        if not current_data_name:
            # sanity check: if line starts with space then there is no data name field and we should have passed a current_data_name
            if line[0] == ' ':
                raise ConfigError("Oh, please. Some silly developer (you know who you are) has tried to call parse_kegg_modules_line() on \
                a line without a data name field, and forgot to give it the current data name. Shame on you, go fix this. (For reference here \
                is the line: %s)" % (line))

            current_data_name = fields[0]
        # note that if data name is known, first field still exists but is actually the empty string ''
        # so no matter which situation, data value is field 1 and data definition (if any) is field 2
        data_vals = fields[1]
        # need to sanity check data value field because SOME modules don't follow the 2-space separation formatting
        vals_are_okay, corrected_vals, corrected_def = self.data_vals_sanity_check(data_vals, current_data_name)

        if vals_are_okay and len(fields) > 2: # not all lines have a definition field
            data_def = fields[2]

        # some types of information may need to be split into multiple db entries
        data_types_to_split = ["ORTHOLOGY","REACTION"] # lines that fall under these categories need to have data_vals split on comma
        if current_data_name in data_types_to_split:
            for val in data_vals.split(','):
                line_entries.append((current_data_name, val, data_def, line_num))
        else: # just send what we found without splitting the line
            line_entries.append((current_data_name, data_vals, data_def, line_num))

        # still need to figure out what to do about REFERENCE info type (includes AUTHORS, TITLE, JOURNAL) - do we want this?
        return line_entries


    def create(self):
        """Creates the Modules DB"""

        self.touch()

        self.progress.new("Loading KEGG modules into Modules DB...")

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
                if not line == '///':
                    # parse the line into a tuple
                    entries_tuple_list = None
                    # here is the tricky bit about parsing these files. Not all lines start with the data_name field; those that don't start with a space.
                    # if this is the case, we need to tell the parsing function what the previous data_name field has been.
                    if line[0] == ' ':
                        entries_tuple_list = self.parse_kegg_modules_line(line, line_number, prev_data_name_field)
                    else:
                        entries_tuple_list = self.parse_kegg_modules_line(line, line_number)

                    # update prev_data_name_field; use the first (and perhaps only) entry by default
                    prev_data_name_field = entries_tuple_list[0][0]

                    # unpack that tuple info
                    for name, val, definition, line in entries_tuple_list:
                        # call append_and_store which will collect db entries and store every 10000 at a time
                        mod_table.append_and_store(self.db, mnum, name, val, definition, line)



            num_modules_parsed += 1
        self.progress.end()

        # give some run info
        self.run.info('Modules database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of KEGG modules', num_modules_parsed, quiet=self.quiet)
        self.run.info('Number of entries', mod_table.get_total_entries(), quiet=self.quiet)

        # record some useful metadata
        self.db.set_meta_value('db_type', 'modules')
        self.db.set_meta_value('num_modules', num_modules_parsed)
        self.db.set_meta_value('total_entries', mod_table.get_total_entries())

        self.db.disconnect()

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
