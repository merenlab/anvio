# -*- coding: utf-8

"""
    Classes to compute completeness estimates based on the information stored in search tables in the
    contigs database.
"""

from collections import Counter

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class Completeness:
    def __init__(self, contigs_db_path, source = None, run = run, progress = progress):
        # hi db
        contigs_db = dbops.ContigsDatabase(contigs_db_path)

        # read info table to get what is available in the db
        info_table = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)

        # identify and remove non-single-copy sources of hmm search results:
        non_singlecopy_sources = set([k for k in info_table.keys() if info_table[k]['search_type'] != 'singlecopy'])
        singlecopy_sources = set([k for k in info_table.keys() if info_table[k]['search_type'] == 'singlecopy'])
        for non_singlecopy_source in non_singlecopy_sources:
            info_table.pop(non_singlecopy_source)

        # read search table (which holds hmmscan hits for splits).
        self.search_table = utils.get_filtered_dict(contigs_db.db.get_table_as_dict(t.hmm_hits_splits_table_name), 'source', singlecopy_sources)

        # an example entry in self.search_table looks loke this:
        #
        # {
        #    'percentage_in_split'   : 100,
        #    'source'                : u'Campbell_et_al',
        #    'gene_unique_identifier': u'c70c1cc3025b636100fd8a910b5b7f0dd09752fc78e2a1f10ee60954',
        #    'e_value'               : 0.0013,
        #    'gene_name'             : u'UvrC_HhH_N',
        #    'split'                 : u'ANTARCTICAAQUATIC_SMPL_SITE231_3.0UMcontig18439_split_00001'
        # }
        #


        # a little convenience for potential clients:
        self.http_refs = {}
        for source_in_db in info_table:
            self.http_refs[source_in_db] = [h for h in info_table[source_in_db]['ref'].split() if h.startswith('http')][0]

        self.genes_in_db = dict([(s, info_table[s]['genes'].split(', ')) for s in info_table])


        # we're done with the db
        contigs_db.disconnect()

        self.sources = info_table.keys()

        if source:
            if source not in self.sources:
                raise ConfigError, 'Source "%s" is not one of the single-copy gene sources found in the database.' % source

            # filter out sources that are not requested
            self.sources = [source]
            self.genes_in_db = {source: self.genes_in_db[source]}
            self.search_table = utils.get_filtered_dict(self.search_table, 'source', set([source]))

        self.unique_gene_id_to_gene_name = {}
        self.splits_unique_gene_id_occurs = {}
        # these will be very useful later. trust me.
        for entry in self.search_table.values():
            if entry['gene_unique_identifier'] not in self.unique_gene_id_to_gene_name:
                self.unique_gene_id_to_gene_name[entry['gene_unique_identifier']] = entry['gene_name']

            if entry['gene_unique_identifier'] not in self.splits_unique_gene_id_occurs:
                self.splits_unique_gene_id_occurs[entry['gene_unique_identifier']] = [entry['split']]
            else:
                self.splits_unique_gene_id_occurs[entry['gene_unique_identifier']].append(entry['split'])


    def get_info_for_splits(self, split_names, min_e_value = 1e-5):
        hits = utils.get_filtered_dict(self.search_table, 'split', split_names)

        # we need to restructure 'hits' into a dictionary that gives access to sources and genes in a more direct manner
        info_dict, gene_name_to_unique_id = {}, {}
        for source in self.sources:
            info_dict[source], gene_name_to_unique_id[source] = {}, {}

        # here we go through every hit and populate 'info_dict' and 'gene_name_to_unique_id':
        for entry in hits.values():
            if entry['e_value'] > min_e_value:
                continue

            source = entry['source']
            e_value = entry['e_value']
            gene_name = entry['gene_name']
            percentage = entry['percentage_in_split']
            gene_unique_id = entry['gene_unique_identifier']

            if info_dict[source].has_key(gene_unique_id):
                info_dict[source][gene_unique_id]['percentage'] += percentage
            else:
                info_dict[source][gene_unique_id] = {}
                info_dict[source][gene_unique_id] = {'gene_name': gene_name, 'percentage': percentage, 'e_value': e_value}

            if gene_name_to_unique_id[source].has_key(gene_name):
                gene_name_to_unique_id[source][gene_name].add(gene_unique_id)
            else:
                gene_name_to_unique_id[source][gene_name] = set([gene_unique_id])

        # here we generate the results information
        results_dict = {}
        for source in self.sources:
            results_dict[source] = {}

        for source in self.sources:
            genes_count = Counter([v['gene_name'] for v in info_dict[source].values()])

            # report results
            results_dict[source]['percent_complete'] = len(genes_count) * 100.0 / len(self.genes_in_db[source])

            # report redundancy:
            genes_that_occur_multiple_times = [g for g in genes_count if genes_count[g] > 1]
            results_dict[source]['percent_redundancy'] = sum([genes_count[g] - 1 for g in genes_that_occur_multiple_times]) * 100.0 / len(self.genes_in_db[source])

            # identify splits that contribute the same single_copy_gene
            redundants = {}
            for gene_name in genes_that_occur_multiple_times:
                redundants[gene_name] = [self.splits_unique_gene_id_occurs[unique_gene_id] for unique_gene_id in gene_name_to_unique_id[source][gene_name]]
            results_dict[source]['redundants'] = redundants

        return results_dict
