"""Some critical constants for metabolism estimation output formatting.

The below dictionary defines the possible output modes.
- 'output_suffix' should be unique to a mode so that multiple output modes can be used at once
- 'data_dict' indicates which data dictionary is used for generating the output (modules or kofams)
- 'headers' list describes which information to include in the output file (see OUTPUT_HEADERS dict below for more info)
- 'description' is what is printed when --list-available-modes parameter is used
"""

import re

OUTPUT_MODES = {'modules': {
                    'output_suffix': "modules.txt",
                    'data_dict': "modules",
                    'headers': ["module", "module_name", "module_class", "module_category",
                                "module_subcategory", "module_definition",
                                "stepwise_module_completeness", "stepwise_module_is_complete",
                                "pathwise_module_completeness", "pathwise_module_is_complete",
                                "proportion_unique_enzymes_present", "enzymes_unique_to_module", "unique_enzymes_hit_counts",
                                "enzyme_hits_in_module", "gene_caller_ids_in_module", "warnings"],
                    'description': "Information on metabolic modules"
                    },
                'modules_custom': {
                    'output_suffix': "modules_custom.txt",
                    'data_dict': "modules",
                    'headers': None,
                    'description': "A custom tab-delimited output file where you choose the included modules data using --custom-output-headers"
                    },
                'module_paths': {
                                    'output_suffix': "module_paths.txt",
                                    'data_dict': "modules",
                                    'headers': ["module", "pathwise_module_completeness", "pathwise_module_is_complete",
                                                "path_id", "path", "path_completeness", "annotated_enzymes_in_path"],
                                    'description': "Information on each possible path (complete set of enzymes) in a module"
                                    },
                'module_steps': {
                                    'output_suffix': "module_steps.txt",
                                    'data_dict': "modules",
                                    'headers': ["module", "stepwise_module_completeness", "stepwise_module_is_complete",
                                                "step_id", "step", "step_completeness"],
                                    'description': "Information on each top-level step in a module"
                                    },
                'hits': {
                    'output_suffix': "hits.txt",
                    'data_dict': "kofams",
                    'headers': ["enzyme", "gene_caller_id", "contig", "modules_with_enzyme", "enzyme_definition", "warnings"],
                    'description': "Information on all enzyme annotations in the contigs DB, regardless of module membership"
                    },
                }
"""
The below dictionary describes the type of information we can output
- the dictionary key corresponds to the header's key in the output dictionary (ie, as returned from generate_output_dict_for_modules() function)
- 'cdict_key' is the header's key in modules or kofams data dictionary (if any)
- 'mode_type' indicates which category of output modes (modules or kofams) this header can be used for. If both, this is 'all'
- 'description' is printed when --list-available-output-headers parameter is used
"""
OUTPUT_HEADERS = {'module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Module number"
                        },
                  'stepwise_module_is_complete' : {
                        'cdict_key': "stepwise_is_complete",
                        'mode_type': 'modules',
                        'description': "Whether a module is considered complete or not based on its STEPWISE percent completeness and the completeness threshold"
                        },
                  'stepwise_module_completeness' : {
                        'cdict_key': 'stepwise_completeness',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a module, computed as the number of complete steps divided by the number of total steps "
                                       "(where 'steps' are determined by splitting the module definition on the space character)"
                        },
                  'pathwise_module_is_complete' : {
                        'cdict_key': "pathwise_is_complete",
                        'mode_type': 'modules',
                        'description': "Whether a module is considered complete or not based on its PATHWISE percent completeness and the completeness threshold"
                        },
                  'pathwise_module_completeness' : {
                        'cdict_key': 'pathwise_percent_complete',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a module, computed as maximum completeness of all possible combinations of enzymes ('paths') in the definition"
                        },
                  'enzymes_unique_to_module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A list of enzymes that only belong to this module (ie, they are not members of multiple modules)"
                        },
                  'unique_enzymes_hit_counts' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "How many times each unique enzyme appears in the sample (order of counts corresponds to list in `enzymes_unique_to_module` field)"
                        },
                  'proportion_unique_enzymes_present' : {
                        'cdict_key': 'proportion_unique_enzymes_present',
                        'mode_type': 'modules',
                        'description': "Proportion of enzymes unique to this one module that are present in the sample"
                        },
                  'unique_enzymes_context_string' : {
                        'cdict_key': 'unique_enzymes_context_string',
                        'mode_type': 'modules',
                        'description': "Describes the unique enzymes contributing to the `proportion_unique_enzymes_present` field"
                        },
                  'module_name' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Name/description of a module"
                        },
                  'module_class' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism class of a module"
                        },
                  'module_category' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism category of a module"
                        },
                  'module_subcategory' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Metabolism subcategory of a module"
                        },
                  'module_definition' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Definition string of a module. Describes the metabolic pathway "
                                       "in terms of the enzymes (KOs, COGs, etc) that belong to the module."
                        },
                  'module_substrates' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of compounds that serve as initial input to the metabolic pathway "
                                       "(that is, substrate(s) to the initial reaction(s) in the pathway)"
                        },
                  'module_products' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of compounds that serve as final output from the metabolic pathway "
                                       "(that is, product(s) of the final reaction(s) in the pathway)"
                        },
                  'module_intermediates' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of compounds that are intermediates in the metabolic pathway "
                                       "(compounds that are both outputs and inputs of reaction(s) in the pathway)"
                        },
                  'gene_caller_ids_in_module': {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of gene caller IDs of enzymes that contribute to a module"
                        },
                  'gene_caller_id': {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': "Gene caller ID of a single enzyme in the contigs DB"
                        },
                  'enzyme_hits_in_module' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Comma-separated list of enzyme annotations that contribute to a module"
                        },
                  'enzyme_hit' : {
                        'cdict_key': 'kofam_hits',
                        'mode_type': 'kofams',
                        'description': "Enzyme identifier for a single annotation (KO, COG, etc)"
                        },
                  'contig' : {
                        'cdict_key': 'genes_to_contigs',
                        'mode_type': 'kofams',
                        'description': "Contig that an enzyme annotation is found on"
                        },
                  'path_id' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Integer ID for a path through a module. Has no real meaning and is used for data organization"
                        },
                  'path' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A path through a module: a linear sequence of enzymes that together represent each metabolic step "
                                       "in the module (most modules have several of these due to enzyme redundancy)"
                        },
                  'path_completeness' : {
                        'cdict_key': 'pathway_completeness',
                        'mode_type': 'modules',
                        'description': "Percent completeness of a given path through a module"
                        },
                  'annotated_enzymes_in_path' : {
                        'cdict_key': 'annotated_enzymes_in_path',
                        'mode_type': 'modules',
                        'description': "Shows which enzymes in the path are annotated in your sample, and which are missing"
                        },
                  'step_id' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Integer ID for a top-level step in a module. Has no real meaning and is used for data organization"
                        },
                  'step' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "A 'top-level' step in a module, represented by one or more possible enzymes that can catalyze "
                                       "a logical part of the metabolic pathway (usually one reaction)"
                        },
                  'step_completeness' : {
                        'cdict_key': None,
                        'mode_type': 'modules',
                        'description': "Binary completeness of a given 'top-level' step in a module"
                        },
                  'warnings' : {
                        'cdict_key': 'warnings',
                        'mode_type': 'all',
                        'description': "This column holds a comma-separated list of notes about things that might affect completeness "
                                       "estimates for a module, such as missing enzyme profiles."
                        },
                  'enzyme' : {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': 'Identifier for an enzyme that is annotated in your database(s), ie a KO or COG number'
                        },
                  'modules_with_enzyme': {
                        'cdict_key': 'modules',
                        'mode_type': 'kofams',
                        'description': 'A comma-separated list of modules that the enzyme belongs to'
                        },
                  'enzyme_definition': {
                        'cdict_key': None,
                        'mode_type': 'kofams',
                        'description': 'The functional annotation associated with the enzyme'
                        },
                  }

DEFAULT_OUTPUT_MODE = 'modules'
STRAY_KO_ANVIO_SUFFIX = "_anvio_version"

# global metadata header lists for matrix format
# if you want to add something here, don't forget to add it to the dictionary in the corresponding
# get_XXX_metadata_dictionary() function
MODULE_METADATA_HEADERS = ["module_name", "module_class", "module_category", "module_subcategory"]
KO_METADATA_HEADERS = ["enzyme_definition", "modules_with_enzyme"]
# Exception: if you add to this list, you must add it in the steps_subdict in generate_subsets_for_matrix_format()
# and to the relevant step metadata clause in write_stat_to_matrix()
STEP_METADATA_HEADERS = ["step_definition"]




