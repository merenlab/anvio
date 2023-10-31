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

## REQUIRED FIELDS FOR A PATHWAY FILE
## AND THEIR ACCEPTABLE VALUES (IF APPLICABLE)
REQ_FIELDS = {'type': {'accepted_vals': ['pathway', 'enzyme set'], 'data_type': str},
              'source': {'accepted_vals': ['KEGG', 'user'], 'data_type': str},
              'functional_definition': {'data_type': list},
              'functions': {'data_type': dict},
            }

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
        
        num_pathways_in_file = len(pathway_dict.keys())
        if anvio.DEBUG:
            run.info("Number of pathways found", num_pathways_in_file)

        if num_pathways_in_file > 1:
            # check if user forgot to nest the pathway info within a dictionary keyed by pathway ID
            for key in pathway_dict.keys():
                if key in list(REQ_FIELDS.keys()):
                    raise ConfigError(f"We noticed that the top-level of the data dictionary contains pathway information, like '{key}'. "
                                      f"The first level of this dictionary should be the pathway ID (such as 'M00001'). Please nest the "
                                      f"pathway info one level lower in your input data.")
            raise ConfigError("We found more than one pathway in the input data, but we are not yet equipped to handle that :/ ")
        
        self.id = list(pathway_dict.keys())[0]
        self.dict = pathway_dict[self.id]

        self.check_required_fields()
    #### ACCESSOR FUNCTIONS ####
    def get_pathway_dict(self):
        """Returns the pathway dictionary."""

        return self.dict

    def print_pathway(self):
        """Prints the pathway to the terminal. Mainly for debugging output."""

        print(self.dict)

    #### SANITY CHECKING FUNCTIONS ####
    def check_required_fields(self):
        """Ensure that all required fields are defined in the pathway file."""

        for field, requirements in REQ_FIELDS.items():
            if field not in self.dict:
                raise ConfigError(f"The required field '{field}' was not found in the input data.")
            if 'accepted_vals' in requirements and (self.dict[field] not in requirements['accepted_vals']):
                raise ConfigError(f"There is an issue with the definition of {field} in the input pathway data. "
                                  f"This field has the value '{self.dict[field]}', but only the following values are accepted: "
                                  f"{', '.join(requirements['accepted_vals'])}")

            if not isinstance(self.dict[field], requirements['data_type']):
                raise ConfigError(f"There is an issue with the definition of {field} in the input pathway data. "
                                  f"This field is supposed to be of type '{requirements['data_type']}' but instead it is a {type(self.dict[field])}.")
                

