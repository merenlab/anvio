#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains classes related to metabolism estimation, especially for user-defined metabolic pathways."""

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.utils as utils

from anvio.errors import ConfigError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Iva Veseli"
__email__ = "iva.veseli@hifmb.de"

run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)

class PathwayYAML:
    """A class to store definitions of metabolic pathways from YAML-formatted data. 

    The data must be provided either in the form of a YAML file, or as a dictionary.
    
    PARAMETERS
    ==========
    file_path : string
        path to input YAML file
    path_dict : dictionary
        dictionary containing the pathway info, ie from an already-parsed YAML file
    """

    def __init__(self, file_path=None, path_dict=None):

        if not file_path and not path_dict:
            raise ConfigError("The PathwayYAML class requires an input. You should provide either a path to a YAML file or a "
                              "dictionary of pathway information.")

        if file_path:
            self.file_path = file_path
            filesnpaths.is_file_exists(self.file_path)
            pathway_dict = utils.get_yaml_as_dict(self.file_path)
            if anvio.DEBUG:
                run.info_single(f"The pathway file {self.file_path} has been successfully loaded.")
        elif path_dict:
            self.file_path = None
            pathway_dict = path_dict
        self.id = list(pathway_dict.keys())[0]
        self.dict = pathway_dict[self.id]

    #### ACCESSOR FUNCTIONS ####
    def get_pathway_dict(self):
        """Returns the pathway dictionary."""

        return self.dict

    def print_pathway(self):
        """Prints the pathway to the terminal. Mainly for debugging output."""

        print(self.dict)

    #### SANITY CHECKING FUNCTIONS ####
