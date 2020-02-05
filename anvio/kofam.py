#!/usr/bin/env python
# -*- coding: utf-8
"""
    This file contains KofamSetup and Kofam classes.

"""

import os
import gzip
import shutil
import requests
import glob
import re

import anvio
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


class KofamContext(object):
    """
    The purpose of this base class is to define shared functions and file paths for all KOfam operations.
    """
    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # default directory will be called KEGG and will store the KEGG Module data as well
        self.kofam_data_dir = A('kofam_data_dir') or os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')
        self.orphan_data_dir = os.path.join(self.kofam_data_dir, "orphan_data")

        # shared variables for all KOfam subclasses
        self.kofam_hmm_file_path = os.path.join(self.kofam_data_dir, "Kofam.hmm") # file containing concatenated KOfam hmms
        self.ko_list_file_path = os.path.join(self.kofam_data_dir, "ko_list")

    def setup_ko_dict(self):
        """
        The purpose of this function is to process the ko_list file into usable form by Kofam sub-classes.

        The ko_list file (which is downloaded along with the KOfam HMM profiles) contains important
        information for each KEGG Orthology number (KO, or knum), incuding pre-defined scoring thresholds
        for limiting HMM hits and annotation information.

        It looks something like this:

        knum	threshold	score_type	profile_type	F-measure	nseq	nseq_used	alen	mlen	eff_nseq	re/pos	definition
        K00001	329.57	domain	trim	0.231663	1473	1069	1798	371	17.12	0.590	alcohol dehydrogenase [EC:1.1.1.1]

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
        """
        The purpose of this function is to determine which KO numbers have no associated data or just no score threshold in the ko_list file.
        That is, their ko_list entries look like this, with hypens in all but the first and last columns:

        K14936	-	-	-	-	-	-	-	-	-	-	small nucleolar RNA snR191
        K15035	-	-	-	-	-	-	-	-	-	-	transfer-messenger RNA
        K15841	-	-	-	-	-	-	-	-	-	-	small regulatory RNA GlmY
        K15851	-	-	-	-	-	-	-	-	-	-	quorum regulatory RNA Qrr
        K16736	-	-	-	-	-	-	-	-	-	-	bantam
        K16863	-	-	-	-	-	-	-	-	-	-	microRNA 21

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

class KofamSetup(KofamContext):
    """ Class for setting up KEGG Kofam HMM profiles. It performs sanity checks and downloads, unpacks, and prepares
    the profiles for later use by `hmmscan`.

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
        KofamContext.__init__(self, self.args)

        filesnpaths.is_program_exists('hmmpress')

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.kofam_data_dir, delete_if_exists=args.reset)
        filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)

        # ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        self.files = ['ko_list.gz', 'profiles.tar.gz']


    def is_database_exists(self):
        """This function determines whether the user has already downloaded the Kofam HMM profiles."""
        if os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError("It seems you already have KOfam HMM profiles installed in '%s', please use --reset flag if you want to re-download it." % self.kofam_data_dir)

    def download(self):
        """This function downloads the Kofam profiles."""
        self.run.info("Database URL", self.database_url)

        for file_name in self.files:
            utils.download_file(self.database_url + '/' + file_name,
                os.path.join(self.kofam_data_dir, file_name), progress=self.progress, run=self.run)


    def decompress_files(self):
        """This function decompresses the Kofam profiles."""
        for file_name in self.files:
            full_path = os.path.join(self.kofam_data_dir, file_name)

            if full_path.endswith("tar.gz"): # extract tar file instead of doing gzip
                utils.tar_extract_file(full_path, output_file_path = self.kofam_data_dir, keep_original=False)
            else:
                utils.gzip_decompress_file(full_path, keep_original=False)

    def confirm_downloaded_files(self):
        """This function verifies that all Kofam profiles have been properly downloaded. It is intended to be run
        after the files have been decompressed. The profiles directory should contain hmm files from K00001.hmm to
        K23763.hmm with some exceptions; all KO numbers from ko_list file (except those in ko_skip_list) should be
        included."""
        ko_nums = self.ko_dict.keys()
        for k in ko_nums:
            if k not in self.ko_skip_list:
                hmm_path = os.path.join(self.kofam_data_dir, "profiles/%s.hmm" % k)
                if not os.path.exists(hmm_path):
                    raise ConfigError("The KOfam HMM profile at %s does not exist. This probably means that something went wrong \
                                    while downloading the KOfam database. Please run `anvi-setup-kegg-kofams` with the --reset \
                                    flag." % (hmm_path))

    def move_orphan_files(self):
        """
        This function moves the following to the orphan files directory:
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

        self.progress.update('Verifying that the Kofam directory at %s contains all HMM profiles' % self.kofam_data_dir)
        self.confirm_downloaded_files()

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

    def setup_profiles(self):
        """This is a driver function which executes the Kofam setup process by downloading, decompressing, and hmmpressing the profiles."""
        self.download()
        self.decompress_files()
        self.setup_ko_dict()
        self.run_hmmpress()

class KofamRunHMMs(KofamContext):
    """ Class for running `hmmscan` against the KOfam database and adding the resulting hits to contigs DBs
    for later metabolism prediction.

    Parameters
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
        KofamContext.__init__(self, self.args)

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
        """This is a driver function for running HMMs against the KOfam database and processing the hits into the
        provided contigs DB"""

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
