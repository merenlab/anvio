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

run = anvio.terminal.Run()


class HMMScan(Parser):
    def __init__(self, hmm_scan_hits_txt, alphabet='AA', context='GENE'):
        self.alphabet = alphabet
        self.context = context

        self.run = run

        files_expected = {'hits': hmm_scan_hits_txt}

        if self.context == "GENE":
            # see the HMMER user guide for details of the fields for AA sequence search, and DNA sequence search.
            #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description
            col_names = ['gene_name', 'gene_hmm_id', 'gene_callers_id', 'f', 'e_value', 'bit_score', 'f', 'f', 'dom_bit_score', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
            col_mapping = [str, str, int, str, float, float, str, str, float, str, str, str, str, str, str, str, str, str]
        elif self.context == "CONTIG" and (self.alphabet == "DNA" or self.alphabet == "RNA"):
            # 'hmm_target', 'hmm_acc', 'query_id', 'query_acc', 'hmm_from', 'hmm_to', 'alignment_from', 'alignment_to', 'envelope_from', 'envelope_to', 'seq_len', 'strand', 'e_value', 'score', 'bias', 'desc']
            col_names = ['gene_name', 'gene_hmm_id', 'contig_name', 'f', 'hmm_from', 'hmm_to', 'alignment_from', 'alignment_to', 'envelope_from', 'envelope_to', 'f', 'f', 'e_value', 'f', 'f', 'f']
            col_mapping = [str, str, str, str, str, str, int, int, int, int, str, str, float, str, str, str]
        else:
            raise ConfigError("HMMScan driver is confused. Yor context and alphaet pair ('%s' and '%s') "
                              "does not seem to be implemented in the parser module. If you think this is "
                              "not a mistake on your part, please get in touch with the anvi'o developers "
                              "and watch them fix it like actual pros." % (self.context, self.alphabet))

        files_structure = {'hits':
                                {'col_names': col_names,
                                 'col_mapping': col_mapping,
                                 'indexing_field': -1,
                                 'no_header': True
                                 },
                        }

        Parser.__init__(self, 'HMMScan', [hmm_scan_hits_txt], files_expected, files_structure)


    def get_search_results(self, ko_list_dict = None):
        """This function goes through the hits provided by `hmmscan` and generates an annotation dictionary with the relevant information about each hit.

        If we are parsing Kofam hits, then this function makes sure only hits with a high enough bit score make it into the annotation dictionary.

        Parameters
        ==========
        ko_list_dict    dictionary of the ko_list file; see setup_ko_dict in kofam.py for more details

        Returns
        =======
        annotations_dict    dictionary of annotations
        """

        annotations_dict = {}

        # this is the stuff we are going to try to fill with this:
        # search_table_structure = ['entry_id', 'source', 'alphabet', 'contig', 'gene_callers_id' 'gene_name', 'gene_hmm_id', 'e_value']

        entry_id = 0
        num_hits_removed = 0 # a counter for the number of hits we don't add to the annotation dictionary
        for hit in list(self.dicts['hits'].values()):
            entry = None
            if self.context == 'GENE':
                # This is for KEGG Kofams. Here we only add the hit to the annotations_dict if the appropriate bit score is above the
                # threshold set in ko_list_dict (which is indexed by ko num, aka gene_name in the hits dict)
                if ko_list_dict and hit['gene_name'] in ko_list_dict.keys():
                    knum =  hit['gene_name']
                    score_type = ko_list_dict[knum]['score_type']
                    threshold = ko_list_dict[knum]['threshold']
                    keep = True
                    if score_type == 'full':
                        if hit['bit_score'] < float(threshold):
                            keep = False
                    elif score_type == 'domain':
                        if hit['dom_bit_score'] < float(threshold):
                            keep = False
                    else:
                        self.run.warning("Oh dear. The Kofam profile %s has a strange score_type value: %s. The only accepted values \
                        for this type are 'full' or 'domain', so anvi'o cannot parse the hits to this profile. All hits will be kept \
                        regardless of bit score. You have been warned." % (hit['gene_name'], score_type))

                    if keep:
                        entry = {'entry_id': entry_id,
                                 'gene_name': hit['gene_name'],
                                 'gene_hmm_id': hit['gene_hmm_id'],
                                 'gene_callers_id': hit['gene_callers_id'],
                                 'e_value': hit['e_value']}
                    else:
                        num_hits_removed += 1

                elif ko_list_dict and hit['gene_name'] not in ko_list_dict.keys():
                    # this should never happen, in an ideal world where everything is filled with butterflies and happiness
                    self.run.warning("Hmm. While parsing your Kofam hits, it seems the Kofam profile %s was not found in the ko_list dictionary. \
                    This should probably not ever happen, and you should contact a developer as soon as possible to figure out what \
                    is going on. But for now, anvi'o is going to keep all hits to this profile. Consider those hits with a grain of salt, \
                    as not all of them may be good." % hit['gene_name'])
                    entry = {'entry_id': entry_id,
                             'gene_name': hit['gene_name'],
                             'gene_hmm_id': hit['gene_hmm_id'],
                             'gene_callers_id': hit['gene_callers_id'],
                             'e_value': hit['e_value']}

                else:
                    # but in Pfams, we don't care, we just keep all hits
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

            if entry:
                entry_id += 1
                annotations_dict[entry_id] = entry

        self.run.info("Number of weak hits removed", num_hits_removed)
        self.run.info("Number of hits in annotation dict ", len(annotations_dict.keys()))

        return annotations_dict
