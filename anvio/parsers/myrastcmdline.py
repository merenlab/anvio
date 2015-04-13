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


from anvio.parsers.base import Parser


class MyRastCMDLine(Parser):
    def __init__(self, input_file_paths, genes_table_structure):
        files_expected = {'functions': 'svr_assign_using_figfams.txt', 'genes': 'svr_call_pegs.txt'}

        files_structure = {'functions': 
                                {'col_names': ['t_species', 'field2', 'prot', 'function'],
                                 'col_mapping': [str, int, str, str],
                                 'indexing_field': 2},
                           'genes': 
                                {'type': 'fasta'},}

        self.genes_table_structure = genes_table_structure
        Parser.__init__(self, 'MyRastCMDLine', input_file_paths, files_expected, files_structure)


    def get_annotations_dict(self):
        # start with the fasta dict to identify start and stop. fasta file looked like this, as a reminder:
        #
        #    >prot_00001 D23-1_contig_1_127_1215
        #    MTDCSXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #    >prot_00002 D23-1_contig_1_1281_1499
        #    MKDNAERKAKRRIFLXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #    >prot_00003 D23-1_contig_1_2630_1626
        #    MAKQKIRIRLKAYDHRVIDQSAEKIVETAKRSGADVSGPIPLPTE
        #    >prot_00004 D23-1_contig_1_2987_3376
        #    MKNGPKLSLALIGIFLILCEFFYGIPFLGATFILSFGWQPLIFNA
        #    >prot_00005 D23-1_contig_1_4901_3567
        #    MKNYFQFDKYGTNFKREILGGITTFLSMAYILAVNPQVLSLAGVK
        #    >prot_00006 D23-1_contig_1_7149_5014
        #    MKSLILAEKPSVARDIADALQINQKRNGYFENNQYIVTWALGHLV
        #    >prot_00007 D23-1_contig_1_7266_8231
        #    MLISLLTFISVEILYNKSNKKYGGNDMSIVQLYDITQIKSFIEHS

        annotations_dict = {}

        for key in self.dicts['genes']:
            entry = {}
            for field in self.genes_table_structure:
                entry[field] = None

            prot, remainder = key.split()
            contig = '_'.join(remainder.split('_')[:-2])
            start, stop = [t for t in remainder.split('_')[-2:]]
            start, stop = int(start), int(stop)
            entry['start'], entry['stop'], entry['direction'] = (start, stop, 'f') if start < stop else (stop, start, 'r')
            entry['contig'] = contig

            annotations_dict[prot] = entry

        for prot in self.dicts['functions']:
            d = self.dicts['functions'][prot]

            t_species_str = d['t_species']
            if len(t_species_str.split()) > 2:
                t_species_str = ' '.join(t_species_str.split()[0:2])
            annotations_dict[prot]['t_species'] = t_species_str
            annotations_dict[prot]['function'] = d['function']

        return annotations_dict

