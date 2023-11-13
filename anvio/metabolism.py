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
P = terminal.pluralize

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
    file_path : str
        path to input YAML file
    path_dict : dict
        dictionary containing the pathway info, ie from an already-parsed YAML file
    """

    def __init__(self, file_path: str = None, path_dict: dict = None):

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

        self.parse_pathway_data()

        self.print_pathway()


    #### ACCESSOR FUNCTIONS ####
    def get_pathway_dict(self):
        """Returns the pathway dictionary."""

        return self.dict

    def print_pathway(self):
        """Prints the pathway to the terminal. Mainly for debugging output."""

        print(self.dict)


    #### PARSING AND SANITY CHECKING FUNCTIONS ####
    def check_required_fields(self):
        """Ensures that all required fields are defined in the pathway file."""

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

    def parse_accessions_from_definition(self, passed_definition: str = None):
        """Returns a list of functional accessions from the definition string. 
        
        If no definition is passed using the `passed_definition` variable, this function works on the definition 
        stored in `self.definition_string`.

        PARAMETERS
        ==========
        passed_definition : str
            The definition string to parse (optional)
        """
        
        def_to_parse = self.definition_string
        if passed_definition:
            def_to_parse = passed_definition

        substrs_to_remove = ['(',')', 'and', 'or']
        for r in substrs_to_remove:
            def_to_parse = def_to_parse.replace(r, '')
        acc_list = def_to_parse.split()
        return acc_list
                
    def crosscheck_definition_and_functions(self):
        """Verifies that every accession in the functional definition has a corresponding entry in the functions dict."""

        acc_missing_from_functions_field = []
        for a in self.parse_accessions_from_definition():
            if a not in self.functions:
                acc_missing_from_functions_field.append(a)
        
        if acc_missing_from_functions_field:
            acc_str = ", ".join(acc_missing_from_functions_field)
            n = len(acc_missing_from_functions_field)
            raise ConfigError(f"There {P('is a function', n, alt='are some functions')} in the functional definition for pathway "
                              f"{self.id} that {P('is', n, alt='are')} missing from the 'functions' field in the input data. Please "
                              f"add an entry for each accession into the 'functions' field. Here {P('is the accession', n, alt='are the accessions')} "
                              f"that you should add before trying to parse the pathway data again: {acc_str}")

    def function_sanity_checks(self):
        """Checks the formatting of each individual function within the self.functions dictionary.
        
        Each function requires at least an annotation 'source'.
        """

        for accession, function_data in self.functions.items():
            if 'source' not in function_data or function_data['source'] is None:
                raise ConfigError(f"The function {accession} does not have an associated annotation source. Please add "
                                  f"one using the key 'source' to ensure that we will be able to find the annotations "
                                  f"corresponding to this function later. If you don't know the annotation source for "
                                  f"this function, this webpage might help: https://anvio.org/help/main/artifacts/user-modules-data/#1-find-the-enzymes")

    def parse_pathway_data(self):
        """Parses the pathway dictionary to obtain the attributes of the pathway. 
        
        This function also performs sanity checking on the data inside critical fields.
        """

        self.check_required_fields()

        # attributes from required fields
        self.source = self.dict['source']
        self.type = self.dict['type']
        self.functional_definition = self.dict['functional_definition']
        self.functions = self.dict['functions']

        # attributes from optional fields that may be important downstream
        self.name = self.dict['name'] if 'name' in self.dict else 'UNDEFINED'

        # create other useful forms of the input data
        self.definition_string = " and ".join(['('+step+')' for step in self.functional_definition])

        # sanity checks
        self.crosscheck_definition_and_functions()
        self.function_sanity_checks()
