# -*- coding: utf-8

"""
Parser for RAST output (the suboptimal version that does now know about hypothetical proteins)

This class is parsing the output of this command:

    svr_assign_to_dna_using_figfams < ../contigs.fa > svr_assign_to_dna_using_figfams.txt

svr_assign_to_dna_using_figfams.txt looks like this:

    204_10M_MERGED.PERFECT.gz.keep_contig_878    530    204_10M_MERGED.PERFECT.gz.keep_contig_878_3443_6817    Carbamoyl-phosphate synthase large chain (EC 6.3.5.5)
    204_10M_MERGED.PERFECT.gz.keep_contig_878    428    204_10M_MERGED.PERFECT.gz.keep_contig_878_12284_14530    Helicase PriA essential for oriC/DnaA-independent DNA replication    Bifidobacterium adolescentis ATCC 15703
    204_10M_MERGED.PERFECT.gz.keep_contig_878    271    204_10M_MERGED.PERFECT.gz.keep_contig_878_18914_20476    Pup ligase PafA' paralog, possible component of postulated heterodimer PafA-PafA'
    204_10M_MERGED.PERFECT.gz.keep_contig_878    316    204_10M_MERGED.PERFECT.gz.keep_contig_878_21745_23202    Pup ligase PafA, possible component of postulated heterodimer PafA-PafA'    Bifidobacterium adolescentis ATCC 15703
    (...)

we almost have everything we need to know about our contigs in this output. however,
this output does not contain information about open reading frames without a known function
(hence the "suboptimal" version)."""


from anvio.parsers.base import Parser


class MyRastCMDLine_DO_NOT_USE(Parser):

    def __init__(self, input_file_paths, genes_table_structure):
        files_expected = {'svr_output': 'svr_assign_to_dna_using_figfams.txt'}

        files_structure = {'svr_output': 
                                {'col_names': ['contig', 'field1', 'prot', 'function', 't_species'],
                                 'col_mapping': [str, int, str, str, str],
                                 'indexing_field': 2}}

        self.genes_table_structure = genes_table_structure
        Parser.__init__(self, 'MyRastCMDLine', input_file_paths, files_expected, files_structure)


    def get_annotations_dict(self):
        annotations_dict = {}

        counter = 1
        for key in self.dicts['svr_output']:
            entry = {}
            prot = 'prot_%.12d' % counter
            counter += 1

            for field in self.genes_table_structure:
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

