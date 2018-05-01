# -*- coding: utf-8

"""Parser for HMMer's hmmscan output"""

import anvio

from anvio.errors import ConfigError
from anvio.parsers.base import Parser


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class HMMScan(Parser):
    def __init__(self, hmm_scan_hits_txt, alphabet='AA', context='GENE'):
        self.alphabet = alphabet
        self.context = context

        files_expected = {'hits': hmm_scan_hits_txt}

        if self.context == "GENE":
            # see the HMMER user guide for details of the fields for AA sequence search, and DNA sequence search.
            col_names = ['gene_name', 'gene_hmm_id', 'gene_callers_id', 'f', 'e_value', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
            col_mapping = [str, str, int, str, float, str, str, str, str, str, str, str, str, str, str, str, str, str]
        elif self.context == "CONTIG" and (self.alphabet == "DNA" or self.alphabet == "RNA"):
            # 'hmm_target', 'hmm_acc', 'query_id', 'query_acc', 'hmm_from', 'hmm_to', 'alignment_from', 'alignment_to', 'envelope_from', 'envelope_to', 'seq_len', 'strand', 'e_value', 'score', 'bias', 'desc']
            col_names = ['gene_name', 'gene_hmm_id', 'contig_name', 'f', 'hmm_from', 'hmm_to', 'alignment_from', 'alignment_to', 'envelope_from', 'envelope_to', 'f', 'f', 'e_value', 'f', 'f', 'f']
            col_mapping = [str, str, str, str, str, str, int, int, int, int, str, str, float, str, str, str]
        else:
            raise ConfigError("HMMScan driver is confused. Yor context and alphaet pair ('%s' and '%s')\
                               does not seem to be implemented in the parser module. If you think this is\
                               not a mistake on your part, please get in touch with the anvi'o developers\
                               and watch them fix it like actual pros." % (self.context, self.alphabet))

        files_structure = {'hits':
                                {'col_names': col_names,
                                 'col_mapping': col_mapping,
                                 'indexing_field': -1,
                                 'no_header': True
                                 },
                        }

        Parser.__init__(self, 'HMMScan', [hmm_scan_hits_txt], files_expected, files_structure)


    def get_search_results(self):
        annotations_dict = {}

        # this is the stuff we are going to try to fill with this:
        # search_table_structure = ['entry_id', 'source', 'alphabet', 'contig', 'gene_callers_id' 'gene_name', 'gene_hmm_id', 'e_value']

        entry_id = 0
        for hit in list(self.dicts['hits'].values()):
            if self.context == 'GENE':
                entry = {'entry_id': entry_id,
                         'gene_name': hit['gene_name'],
                         'gene_hmm_id': hit['gene_hmm_id'],
                         'gene_callers_id': hit['gene_callers_id'],
                         'e_value': hit['e_value']}
            elif self.context == 'CONTIG' and (self.alphabet == 'DNA' or self.alphabet == 'RNA'):
                entry = {'entry_id': entry_id,
                         'gene_name': hit['gene_name'],
                         'gene_hmm_id': hit['gene_hmm_id'],
                         'contig_name': hit['contig_name'],
                         'start': hit['alignment_from'],
                         'stop': hit['alignment_to'],
                         'e_value': hit['e_value']}
            else:
                raise ConfigError("Anvi'o does not know how to parse %s:%s" % (self.alphabet, self.context))

            entry_id += 1
            annotations_dict[entry_id] = entry

        return annotations_dict
