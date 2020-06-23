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
import pandas as pd


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class InteracdomeTableData(object):
    """Manages query and access of the interacdome tabular data sets"""

    def __init__(self, kind='representable', interacdome_data_dir=None):
        self.kind = kind
        self.interacdome_data_dir = interacdome_data_dir

        if not self.interacdome_data_dir:
            self.interacdome_data_dir = constants.default_interacdome_data_path

        self.files = {
            'representable': os.path.join(self.interacdome_data_dir, 'representable_interactions.txt'),
            'confident': os.path.join(self.interacdome_data_dir, 'confident_interactions.txt'),
        }

        if kind not in self.files:
            raise ConfigError("Unknown kind '%s' of Interacdome data. Known kinds: %s" \
                              % (kind, ','.join(list(self.files.keys()))))

        self.filepath = self.files[self.kind]
        filesnpaths.is_file_exists(self.filepath)


    def get_as_dataframe(self):
        """Return the dataset as a dataframe"""

        df = pd.read_csv(os.path.join(self.interacdome_data_dir, self.filepath), sep='\t', comment='#')
        df.index = df['pfam_id'].str.split('_', n=1, expand=True)[0]

        return df


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

        self.pfam_setup = None

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

        self.run.warning('', header='Setting up Interacdome', lc='yellow')
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
        """The main method of this class. Sets up the Interacdome data directory for usage"""

        self.run.warning('', header='Downloading Interacdome tables', lc='yellow')
        self.download_interacdome_files()

        self.run.warning('', header='Downloading associated Pfam HMM profiles', lc='yellow')
        self.download_pfam()

        self.run.warning('', header='Filtering Pfam HMM profiles', lc='yellow')
        self.filter_pfam()


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


    def download_pfam(self):
        """Setup the pfam data subset used by interacdome

        Currently, interacdome only works for pfam version 31.0, so that is the version downloaded here.
        """

        pfam_args = argparse.Namespace(
            pfam_data_dir=self.interacdome_data_dir,
            pfam_version='31.0',
            reset=False,
        )

        self.pfam_setup = pfam.PfamSetup(pfam_args)
        self.pfam_setup.get_remote_version()
        self.pfam_setup.download(hmmpress_files=False)


    def load_interacdome(self, kind='representable'):
        """Loads the representable interacdome dataset as pandas df"""

        data = InteracdomeTableData(kind=kind, interacdome_data_dir=self.interacdome_data_dir)
        return data.get_as_dataframe()


    def get_interacdome_pfam_accessions(self):
        """Get the representable interacdome Pfam accessions"""

        return set(self.load_interacdome(kind='representable').index.tolist())


    def filter_pfam(self):
        """Filter Pfam data according to whether the ACC is in the Interacdome dataset"""

        interacdome_pfam_accessions = self.get_interacdome_pfam_accessions()

        hmm_profiles = pfam.HMMProfile(os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm'))
        hmm_profiles.filter(by='ACC', subset=interacdome_pfam_accessions, filepath=None)

        # hmmpresses the new .hmm
        self.pfam_setup.hmmpress_files()

        # We also filter out the Pfam-A.clans.tsv file, since it is used as a catalog
        clans_file = os.path.join(self.interacdome_data_dir, 'Pfam-A.clans.tsv')
        clans = pd.read_csv(clans_file, sep='\t', header=None)
        clans[clans[0].isin(interacdome_pfam_accessions)].to_csv(clans_file, sep='\t', index=False, header=False)


