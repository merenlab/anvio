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
    """A YAML-formatted file to store definitions of metabolic pathways.
    
    PARAMETERS
    ==========
    file_path : string
        path to input YAML file (required)
    """

    def __init__(self, file_path):
        filesnpaths.is_file_exists(file_path)
        pathway_dict = utils.get_yaml_as_dict(file_path)

        run.info_single(f"The pathway file {file_path} has been successfully loaded.")

