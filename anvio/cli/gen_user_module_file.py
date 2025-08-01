#!/usr/bin/env python
# -*- coding: utf-8
"""
This script converts a file with a list of KOs (one KO per line) into a modules file. The file should be tab-delimited and
contain 3 columns: the KO identifier, the orthology, and the annotation source. If annotation source is KOfam, the orthology
field can be left blank because the script can look it up.

Users can pass the ENTRY, NAME, DEFINITION, and CLASS lines as parameters.
"""

import sys
import os

import anvio
import anvio.kegg as kegg
import anvio.utils as utils

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva']
__provides__ = ["user-modules-data"]
__requires__ = ["enzymes-list-for-module"]
__description__ = ("This script generates a user-defined module file from a tab-delimited file of enzymes and other input parameters.")


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    enzyme_txt_path = A('enzymes_list_for_module')
    entry = A('module_entry')
    name = A('module_name')
    class_str = A('module_class')
    def_str =  A('module_definition')
    kegg_dir_path = A('kegg_data_dir') or os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')

    # LOADING DATA AND SANITY CHECKS
    # sanity check entry is just one word
    if len(entry.split(' ')) > 1:
        raise ConfigError("You did not provide a valid ENTRY ID. It should be just one word.")
    # sanity check module file does not yet exist
    if os.path.exists("./" + entry):
        raise ConfigError(f"A module file for {entry} already exists in this directory. Please delete that file first if you want to re-generate it.")
    # sanity check class format
    if len(class_str.split(";")) != 3:
        raise ConfigError("You did not provide a valid CLASS string. It should be three sections separated by semi-colons.")

    # load enzymes txt
    enzyme_dict = utils.get_TAB_delimited_file_as_dictionary(enzyme_txt_path, expected_fields=['enzyme','source','orthology'])

    if def_str:
        #sanity check all enzymes from definition are in file
        enzyme_set = set(kegg.module_definition_to_enzyme_accessions(def_str))
        dict_set = set(enzyme_dict.keys())
        missing_enzymes = enzyme_set.difference(dict_set)
        if missing_enzymes:
            missing_str = ", ".join(missing_enzymes)
            raise ConfigError(f"Oh no... anvi'o found a discrepancy between the enzymes-txt and the module definition you provided. "
                              f"All enzymes in the definition should be in the text file. Here are the one(s) that are not: {missing_str}")
    else:   # we create definition from enzymes in file
        def_str = " ".join(enzyme_dict.keys())

    # if any orthology line is blank, access db
    enzyme_has_blank_orthology = [x for x in enzyme_dict if not enzyme_dict[x]['orthology']]
    if enzyme_has_blank_orthology:
        annotations_sources_for_blank_orthology = set([enzyme_dict[x]['source'] for x in enzyme_has_blank_orthology])
        acceptable_sources = set(['KOfam'])
        unacceptable_sources = annotations_sources_for_blank_orthology.difference(acceptable_sources)
        if unacceptable_sources:
            ok_str = ", ".join(acceptable_sources)
            bad_str = ", ".join(unacceptable_sources)
            raise ConfigError(f"Some enzymes in your input file have an empty orthology field, but are not from an annotation source "
                              f"that we can automatically fill in. Currently, we are only able to look up orthology for enzymes from "
                              f"source(s): {ok_str}. But the enzymes with blank orthology are coming from source(s): {bad_str}")

        ko_file_path = os.path.join(kegg_dir_path, "ko_list.txt")
        if not os.path.exists(ko_file_path):
            raise ConfigError(f"There does not appear to be a KO file in the KEGG data directory that you provided. This is "
                              f"unfortunately an essential file for filling in orthology lines for KOfam enzymes, so anvi'o has "
                              f"to fail now. You can fix this by either passing in a proper KEGG data directory, or filling in the "
                              f"empty 'orthology' fields in the enzymes-txt file yourself. (Just a reminder, we were looking for the "
                              f"KO file at this path: {ko_file_path})")
        ko_dict = utils.get_TAB_delimited_file_as_dictionary(ko_file_path, expected_fields = ['definition'])

    # WRITE MODULE FILE
    orthology_lines = []
    source_lines = []
    for acc,d in enzyme_dict.items():
        # first line gets the data name
        if len(orthology_lines) == 0:
            o_line_prefix = "ORTHOLOGY   "
            a_line_prefix = "ANNOTATION_SOURCE  "
        else:
            o_line_prefix = a_line_prefix = "            "
            a_line_prefix = "                   "

        if acc in enzyme_has_blank_orthology:
            ortho = ko_dict[acc]['definition']
        else:
            ortho = d['orthology']
        src = d['source']

        orthology_lines.append(o_line_prefix + acc + "  " + ortho)
        source_lines.append(a_line_prefix + acc + "  " + src)


    all_lines = [f"ENTRY       {entry}",
              f"NAME        {name}",
              f"DEFINITION  {def_str}"] + \
              orthology_lines + \
             [f"CLASS       {class_str}"] + \
              source_lines + \
             ['///' + '\n']

    with open(entry, 'w') as f:
        f.write("\n".join(all_lines))


def get_args():
    parser = ArgumentParser(description=__description__)

    groupR = parser.add_argument_group('REQUIRED INPUT', "Everything you NEED to make a module definition.")
    groupR.add_argument(*anvio.A('enzymes-list-for-module'), **anvio.K('enzymes-list-for-module', {'required': True,}))
    groupR.add_argument(*anvio.A('module-entry'), **anvio.K('module-entry'))
    groupR.add_argument(*anvio.A('module-name'), **anvio.K('module-name'))
    groupR.add_argument(*anvio.A('module-class'), **anvio.K('module-class'))

    groupD = parser.add_argument_group('OPTIONAL - MODULE DEFINITION', "If you don't pass a module definition, anvi'o will "
                                                                       "make each enzyme in the enzymes-list-for-module file "
                                                                       "a distinct step in the module. Useful if you are not "
                                                                       "defining a metabolic pathway, but rather a "
                                                                       "characteristic set of genes (ie, transporters, "
                                                                       "tRNA modification enzymes, etc).")
    groupD.add_argument(*anvio.A('module-definition'), **anvio.K('module-definition'))

    groupK = parser.add_argument_group('OPTIONAL - KEGG DATA', "If you want anvi'o to fill in the ORTHOLOGY of "
                                                               "enzymes from source 'KOfam', you might need to "
                                                               "tell her where to find any non-default KEGG data. ")
    groupK.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
