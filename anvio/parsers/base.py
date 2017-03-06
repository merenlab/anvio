# -*- coding: utf-8
"""
    Base class for parsers to take care of the boring stuff.
"""

import os
import hashlib

import anvio.tables as t
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.utils import get_TAB_delimited_file_as_dictionary as get_dict
from anvio.utils import get_FASTA_file_as_dictionary as get_dict_f


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class Parser(object):
    def __init__(self, annotation_source, input_file_paths, files_expected={}, files_structure={}):
        self.annotation_source = annotation_source
        self.input_file_paths = input_file_paths
        self.files_expected = files_expected
        self.files_structure = files_structure
        self.input_file_names = [os.path.basename(p) for p in input_file_paths]
        self.paths = {}
        self.dicts = {}

        if len(input_file_paths) != len(files_expected):
            raise ConfigError("This parser (%s) requires %d file(s), but %d of them were sent. This class is now\
                                confused :/" % (self.annotation_source, len(files_expected), len(input_file_paths)))

        if sorted(files_expected.keys()) != sorted(files_structure.keys()):
            raise ConfigError("Items in files_expected and files_structure must match.")

        missing_files = []
        for f in list(self.files_expected.values()):
            if os.path.basename(f) not in self.input_file_names:
                missing_files.append(f)
        if missing_files:
            if sorted(missing_files) == sorted(self.files_expected.values()):
                raise ConfigError("%s parser requires these file(s): %s. Please refer to the documentation if you\
                                    don't know how to generate them" % (self.annotation_source,
                                                                        ', '.join(list(self.files_expected.values()))))

            raise ConfigError("%s parser requires %d files (%s). %s missing from your input: %s"\
                                     % (self.annotation_source,
                                        len(self.files_expected),
                                        ', '.join(list(self.files_expected.values())),
                                        "These files were" if len(missing_files) > 1 else "This file was",
                                        ", ".join(missing_files)))

        for alias in self.files_expected:
            for i in range(0, len(self.input_file_names)):
                file_name = self.input_file_names[i]
                if os.path.basename(self.files_expected[alias]) == file_name:
                    self.paths[alias] = self.input_file_paths[i]

        for alias in self.files_expected:
            f = self.files_structure[alias]
            if 'type' in f:
                if f['type'] == 'fasta':
                    self.dicts[alias] = get_dict_f(self.paths[alias])
                else:
                    raise ConfigError("Parser class does not know about file type '%s' :/" % f['type'])
            else:
                # then it is tab-delimited
                no_header = f['no_header'] if 'no_header' in f else False
                separator = f['separator'] if 'separator' in f else '\t'
                indexing_field = f['indexing_field'] if 'indexing_field' in f else 0
                self.dicts[alias] = get_dict(self.paths[alias], no_header=no_header,
                                             column_names=self.files_structure[alias]['col_names'],
                                             column_mapping=self.files_structure[alias]['col_mapping'],
                                             indexing_field=indexing_field, separator=separator,
                                             ascii_only=True)

class TaxonomyHelper(object):
    """This is the class that takes an annotations dictionary, and returns
       genes_taxonomy, and taxon names dicts.

       annotations dictionary must contain t_species, t_genus,
       t_family, t_class, t_order, and t_phylum as keys for each
       gene call:

            annotations_dict[gene_call_id]['t_species'] = 'Bifidobacterium longum'

       the purpose of this class is to return to dictionaries that removes the
       redundancy of taxon names, by creating a dict that contains each unique
       taaxonomy found, and a secondary dict that associates each gene call with
       the names dictionary.
       """
    def __init__(self, annotations_dict, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress
        self.annotations_dict = annotations_dict


    def get_genes_taxonomy_and_taxon_names_dicts(self):
        """Return genes_taxonomy_dict, and taxon_names_dict"""
        genes_taxonomy_dict = {}
        taxon_names_dict = {}
        hash_to_taxon_name_id = {}
        taxon_name_id_counter = 1

        taxon_names = t.taxon_names_table_structure[1:]

        for gene_callers_id in self.annotations_dict:
            t_hash = hashlib.sha224(''.join(self.annotations_dict[gene_callers_id][taxon] or '' for taxon in taxon_names).encode('utf-8')).hexdigest()

            if t_hash in hash_to_taxon_name_id:
                taxon_name_id = hash_to_taxon_name_id[t_hash]
            else:
                hash_to_taxon_name_id[t_hash] = taxon_name_id_counter
                taxon_names_dict[taxon_name_id_counter] = self.annotations_dict[gene_callers_id]
                taxon_name_id = taxon_name_id_counter

                taxon_name_id_counter += 1

            genes_taxonomy_dict[gene_callers_id] = taxon_name_id

        return genes_taxonomy_dict, taxon_names_dict


if __name__ == "__main__":
    pass
