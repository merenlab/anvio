#!/usr/bin/env python
# -*- coding: utf-8
"""
    This file contains KofamSetup and Kofam classes.

"""

import os
import gzip
import shutil
import requests

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class KofamSetup(object):
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
        self.kofam_data_dir = args.kofam_data_dir

        filesnpaths.is_program_exists('hmmpress')

        # default directory will be called KEGG and will store the KEGG Module data as well
        if not self.kofam_data_dir:
            self.kofam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')

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
        if os.path.exists(os.path.join(self.kofam_data_dir, 'K00001.hmm')): # TODO: update this after determining final structure
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
        K23763.hmm with some exceptions."""
        skip_list = [17, 47, 56, 80, 92, 110] # the KO profiles that don't exist, based on ko_list
        for k in range(1,23764): # there is likely a better way to do this. Perhaps we should process the ko_list file into a dict first
            if k not in skip_list:
                hmm_path = os.path.join(self.kofam_data_dir, "profiles/K%05d.hmm" % k)
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

        self.progress.end()

    def setup_profiles(self):
        """This is a driver function which executes the Kofam setup process by downloading, decompressing, and hmmpressing the profiles."""
        self.download()
        self.decompress_files()
        # TODO: add concatenation and hmmpress
