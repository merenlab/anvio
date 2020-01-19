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
        self.ko_skip_list = self.get_ko_skip_list()

    def get_ko_skip_list(self):
        """
        The purpose of this function is to determine which KO numbers have no associated data in the ko_list file.
        That is, their ko_list entries look like this, with hypens in all but the first and last columns:

        K14936	-	-	-	-	-	-	-	-	-	-	small nucleolar RNA snR191
        K15035	-	-	-	-	-	-	-	-	-	-	transfer-messenger RNA
        K15841	-	-	-	-	-	-	-	-	-	-	small regulatory RNA GlmY
        K15851	-	-	-	-	-	-	-	-	-	-	quorum regulatory RNA Qrr
        K16736	-	-	-	-	-	-	-	-	-	-	bantam
        K16863	-	-	-	-	-	-	-	-	-	-	microRNA 21

        These are RNAs.

        Returns: skip_list  list of strings, each string is a KO number
        """
        col_names_to_check = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos"]
        skip_list = []
        for k in self.ko_dict.keys():
            should_skip = True
            for c in col_names_to_check:
                if not self.ko_dict[k][c] == "-":
                    should_skip = False
                    break # here we stop checking this KO num because we already found a value in our columns of interest
            if should_skip: # should be True unless we found a value above
                skip_list.append(k)
        return skip_list

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
        K23763.hmm with some exceptions; all KO numbers from ko_list file should be included."""
        ko_nums = self.ko_dict.keys()
        for k in ko_nums:
            if k not in self.ko_skip_list:
                hmm_path = os.path.join(self.kofam_data_dir, "profiles/%s.hmm" % k)
                if not os.path.exists(hmm_path):
                    raise ConfigError("The KOfam HMM profile at %s does not exist. This probably means that something went wrong \
                                    while downloading the KOfam database. Please run `anvi-setup-kegg-kofams` with the --reset \
                                    flag." % (hmm_path))


    def run_hmmpress(self):
        """This function concatenates the Kofam profiles and runs hmmpress on them."""
        self.progress.new('Preparing Kofam HMM Profiles')
        log_file_path = os.path.join(self.kofam_data_dir, '00_hmmpress_log.txt')

        self.progress.update('Verifying that the Kofam directory at %s contains all HMM profiles' % self.kofam_data_dir)
        self.confirm_downloaded_files()

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
        hmm_hits_file = hmmer.run_hmmscan('KOfam', 'AA', 'GENE', None, None, len(self.ko_dict), self.kofam_hmm_file_path, None, '--cut_ga')
        ## TODO: here we have an issue. the number of genes in our HMM model (len(self.ko_dict) is wrong, because ko_dict (derived from ko_list)
        ## does not contain the same number of entries as there were .hmm files originally in the downloaded KOfam profiles...
        ## As far as I can tell this argument is not used for much except for printing out the value, but it would be nice if it was right

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        if not hmm_hits_file:
            run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But\
                             now anvi'o will add KOfam as a functional source with no hits, clean the temporary directories\
                             and gracefully quit.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})
            return

        # parse hmmscan output
        parser = parser_modules['search']['hmmscan'](hmm_hits_file, alphabet='AA', context='GENE')
        search_results_dict = parser.get_search_results()
