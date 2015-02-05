#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.


import os
import PaPi.filesnpaths as filesnpaths

from PaPi.utils import ConfigError
from PaPi.utils import get_TAB_delimited_file_as_dictionary as get_dict
from PaPi.utils import store_dict_as_TAB_delimited_file as store_dict
from PaPi.parsers.base import Parser


class MyRastCMDLine_DO_NOT_USE(Parser):
    """
    OK. This class was parsing the output of this command:
    
        svr_assign_to_dna_using_figfams < ../contigs.fa > svr_assign_to_dna_using_figfams.txt

    svr_assign_to_dna_using_figfams.txt looked like this:

        204_10M_MERGED.PERFECT.gz.keep_contig_878    530    204_10M_MERGED.PERFECT.gz.keep_contig_878_3443_6817    Carbamoyl-phosphate synthase large chain (EC 6.3.5.5)
        204_10M_MERGED.PERFECT.gz.keep_contig_878    428    204_10M_MERGED.PERFECT.gz.keep_contig_878_12284_14530    Helicase PriA essential for oriC/DnaA-independent DNA replication    Bifidobacterium adolescentis ATCC 15703
        204_10M_MERGED.PERFECT.gz.keep_contig_878    271    204_10M_MERGED.PERFECT.gz.keep_contig_878_18914_20476    Pup ligase PafA' paralog, possible component of postulated heterodimer PafA-PafA'
        204_10M_MERGED.PERFECT.gz.keep_contig_878    316    204_10M_MERGED.PERFECT.gz.keep_contig_878_21745_23202    Pup ligase PafA, possible component of postulated heterodimer PafA-PafA'    Bifidobacterium adolescentis ATCC 15703
        (...)

    which was great, because we had almost everything we needed to know about our contigs. however,
    this output did not contain information about open reading frames without a known function. so we
    had to go back to the shitty implementation. If we can find a way to 
    """

    def __init__(self, input_file_paths, annotation_table_structure):
        files_expected = {'svr_output': 'svr_assign_to_dna_using_figfams.txt'}

        files_structure = {'svr_output': 
                                {'col_names': ['contig', 'field1', 'prot', 'function', 't_species'],
                                 'col_mapping': [str, int, str, str, str],
                                 'indexing_field': 2}}

        self.annotation_table_structure = annotation_table_structure
        Parser.__init__(self, 'MyRastCMDLine', input_file_paths, files_expected, files_structure)


    def get_annotations_dict(self):
        annotations_dict = {}

        counter = 1
        for key in self.dicts['svr_output']:
            entry = {}
            prot = 'prot_%.12d' % counter
            counter += 1

            for field in self.annotation_table_structure:
                entry[field] = None

            start, stop = [t for t in key.split('_')[-2:]]
            start, stop = int(start), int(stop)
            entry['start'], entry['stop'], entry['direction'] = (start, stop, 'f') if start < stop else (stop, start, 'r')
            entry['contig'] = self.dicts['svr_output'][key]['contig']
            entry['function'] = self.dicts['svr_output'][key]['function']

            t_species_str = self.dicts['svr_output'][key]['t_species']
            if len(t_species_str.split()) > 2:
                t_species_str = ' '.join(t_species_str.split()[0:2])
            entry['t_species'] = t_species_str

            annotations_dict[prot] = entry

        return annotations_dict

