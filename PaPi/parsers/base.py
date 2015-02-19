#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2015, PaPi Developers
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

"""
    Base class for parsers to take care of the boring stuff.
"""

import os

from PaPi.utils import ConfigError
from PaPi.utils import get_TAB_delimited_file_as_dictionary as get_dict
from PaPi.utils import get_FASTA_file_as_dictionary as get_dict_f


class Parser(object):
    def __init__(self, annotation_source, input_file_paths, files_expected = {}, files_structure = {}):
        self.annotation_source = annotation_source
        self.input_file_paths = input_file_paths
        self.files_expected = files_expected
        self.files_structure = files_structure
        self.input_file_names = [os.path.basename(p) for p in input_file_paths]
        self.paths = {}
        self.dicts = {}

        if sorted(files_expected.keys()) != sorted(files_structure.keys()):
            raise ConfigError, "Items in files_expected and files_structure must match."

        missing_files = []
        for f in self.files_expected.values():
            if os.path.basename(f) not in self.input_file_names:
                missing_files.append(f)
        if missing_files:
            if sorted(missing_files) == sorted(self.files_expected.values()):
                raise ConfigError, "%s parser requires these file(s): %s. Please refer to the documentation if you\
                                    don't know how to generate them" % (self.annotation_source,
                                                                        ', '.join(self.files_expected.values()))

            raise ConfigError, "%s parser requires %d files (%s). %s missing from your input: %s"\
                                     % (self.annotation_source,
                                        len(self.files_expected),
                                        ', '.join(self.files_expected.values()),
                                        "These files were" if len(missing_files) > 1 else "This file was",
                                        ", ".join(missing_files))

        for alias in self.files_expected:
            for i in range(0, len(self.input_file_names)):
                file_name = self.input_file_names[i]
                if os.path.basename(self.files_expected[alias]) == file_name:
                    self.paths[alias] = self.input_file_paths[i]

        for alias in self.files_expected:
            f = self.files_structure[alias]
            if f.has_key('type'):
                if f['type'] == 'fasta':
                    self.dicts[alias] = get_dict_f(self.paths[alias])
                else:
                    raise ConfigError, "Parser class does not know about file type '%s' :/" % f['type']
            else:
                # then it is tab-delimited
                no_header = f['no_header'] if f.has_key('no_header') else False
                separator = f['separator'] if f.has_key('separator') else '\t'
                indexing_field = f['indexing_field'] if f.has_key('indexing_field') else 0
                self.dicts[alias] = get_dict(self.paths[alias], no_header = no_header,
                                             column_names = self.files_structure[alias]['col_names'],
                                             column_mapping = self.files_structure[alias]['col_mapping'],
                                             indexing_field = indexing_field, separator = separator)


if __name__ == "__main__":
    pass
