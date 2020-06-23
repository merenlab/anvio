#!/usr/bin/env python
# -*- coding: utf-8
"""Setup and utilize Interacdome.

The Interacdome is from Mona Singh's group at Princeton. In brief, they attribute empirical
residue-by-residue binding frequencies for Pfam families. Only families with members that have
crystallized structures with ligand co-complexes are used, since otherwise there are no binding
frequencies to attribute. The assumption is that if a residue is close to a bound ligand, it gets a
non-zero weight for being important for binding.

References
==========
- https://interacdome.princeton.edu/
- Systematic domain-based aggregation of protein structures highlights DNA-, RNA- and other
  ligand-binding positions.  Shilpa Nadimapalli Kobren and Mona Singh. Nucleic Acids Research (2019)
  47: 582-593.
"""

import anvio
import anvio.pfam as pfam
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

import os
import argparse


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class InteracdomeSetup(object):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        """Setup a Pfam database for anvi'o

        Parameters
        ==========
        args : argparse.Namespace
            See `bin/anvi-setup-interacdome` for available arguments
        """

        self.run = run
        self.progress = progress
        self.interacdome_data_dir = args.interacdome_data_dir

        self.interacdome_files = {
            'representable_interactions.txt': 'https://interacdome.princeton.edu/session/344807c477180230032a2a3807ba47c3/download/downloadBP?w=',
            'confident_interactions.txt': 'https://interacdome.princeton.edu/session/344807c477180230032a2a3807ba47c3/download/downloadConfidentBP?w=',
        }

        if self.interacdome_data_dir and args.reset:
            raise ConfigError("You are attempting to run Interacdome setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --interacdome-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with Interacdome data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.interacdome_data_dir, self.interacdome_data_dir))

        if not self.interacdome_data_dir:
            self.interacdome_data_dir = constants.default_interacdome_data_path

        self.run.info('Data directory', self.interacdome_data_dir)
        self.run.info('Reset contents', args.reset)

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.interacdome_data_dir)))

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.interacdome_data_dir, delete_if_exists=args.reset)


    def is_database_exists(self):
        """Raise ConfigError if database exists

        Currently, this primitively decides if the database exists by looking at whether Pfam-A.hmm
        or Pfam-A.hmm.gz exists.
        """

        if (os.path.exists(os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm') or 
            os.path.exists(os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm.gz')))):
            raise ConfigError("It seems you already have the Interacdome data downloaded in '%s', please "
                              "use --reset flag if you want to re-download it." % self.interacdome_data_dir)


    def setup(self):
        self.run.info_single('Downloading Interacdome tables', nl_after=1, nl_before=1)
        self.download_interacdome_files()

        self.run.info_single('Downloading associated Pfam HMM profiles', nl_after=1, nl_before=1)
        self.download_pfam_subset()


    def download_pfam_subset(self):
        """Setup the pfam data subset used by interacdome

        Currently, interacdome only works for pfam version 31.0, so that is the version downloaded here.
        After downloading, the pfam hmm is filtered to only include those in the interacdome.
        """

        pfam_args = argparse.Namespace(
            pfam_data_dir=self.interacdome_data_dir,
            pfam_version='31.0',
            reset=False,
        )

        pfam_setup = pfam.PfamSetup(pfam_args)
        pfam_setup.get_remote_version()
        pfam_setup.download()

        # filter pfam FIXME
        pass


    def download_interacdome_files(self):
        """Download the confident and representable non-redundant Interacdome datasets

        These datasets can be found at the interacdome webpage: https://interacdome.princeton.edu/
        """

        for path, url in self.interacdome_files.items():
            utils.download_file(
                url,
                os.path.join(self.interacdome_data_dir, path),
                check_certificate=False,
                progress=self.progress,
                run=self.run
            )
