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

from collections import Counter

import anvio.terminal as terminal

from anvio.parsers.base import Parser
from anvio.parsers.base import TaxonomyHelper


class Centrifuge(Parser):
    def __init__(self, input_file_paths, taxonomy_table_structure, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        self.min_hit_score = 250

        files_expected = {'report': 'centrifuge_report.tsv', 'hits': 'centrifuge_hits.tsv'}

        files_structure = {'report':
                                {'col_names': ['t_species', 'taxon_id', 'f1', 'f2', 'f3', 'f4', 'f5'],
                                 'col_mapping': [str, int, str, str, str, str, str],
                                 'indexing_field': 1},
                           'hits':
                                {'col_names': ['gene_callers_id', 'f1', 'taxon_id', 'score', 'f2', 'f3', 'f4', 'f5'],
                                 'col_mapping': [lambda x: int(x.split('|')[0]), str, int, int, str, str, str, str],
                                 'indexing_field': -1},
                          }

        self.taxonomy_table_structure = taxonomy_table_structure
        Parser.__init__(self, 'centrifuge', input_file_paths, files_expected, files_structure)


    def process(self):
        """Parse two files, returns two dicts: genes_taxonomy, taxon_names.
        
        the file self.dicts['report'] points should look like this:

            name                                     taxID    taxRank  genomeSize  numReads  numUniqueReads  abundance
            Enterobacter_cloacae                     550      species  5123200     1         1               0.0
            Bacillus_cereus                          1396     species  6217228     1         0               0.0
            Bacillus_thuringiensis                   1428     species  6903565     1         0               0.0
            Clostridium_butyricum                    1492     species  9246221     1         0               0.0

        the file self.dicts['hits'] points should look like this:
        
            readID                                                                                                               seqID          taxID    score    2ndBestScore  hitLength  queryLength  numMatches
            0|contig:204_10M_MERGED.PERFECT.gz.keep_contig_878|start:0|stop:933|direction:f|rev_compd:False|length:933           NZ_CP007443.1  1680     842724   0             933        933          1
            1|contig:204_10M_MERGED.PERFECT.gz.keep_contig_878|start:1113|stop:1677|direction:f|rev_compd:False|length:564       NZ_CP010437.1  1680     301401   301401        564        564          3
            1|contig:204_10M_MERGED.PERFECT.gz.keep_contig_878|start:1113|stop:1677|direction:f|rev_compd:False|length:564       NC_008618.1    367928   301401   301401        564        564          3
            1|contig:204_10M_MERGED.PERFECT.gz.keep_contig_878|start:1113|stop:1677|direction:f|rev_compd:False|length:564       NZ_CP007443.1  1680     301401   301401        564        564          3

        """
        annotations_dict = {}

        report = {}
        hits = {}

        for taxon_id in self.dicts['report']:
            taxon = self.dicts['report'][taxon_id]['t_species']
            report[taxon_id] = {'t_species': ' '.join(taxon.split()[0:2])}

        # we are done with this one:
        del self.dicts['report']

        self.run.info('Total num hits found', len(self.dicts['hits']))

        num_hits_below_hit_score = 0
        for hit in self.dicts['hits'].values():
            if hit['score'] < self.min_hit_score:
                num_hits_below_hit_score += 1
                continue

            gene_callers_id = hit['gene_callers_id']
            taxon_id = hit['taxon_id']

            if not taxon_id in report:
                continue

            taxon = report[taxon_id]['t_species']

            if gene_callers_id not in hits:
                hits[gene_callers_id] = Counter()

            hits[gene_callers_id][taxon] += 1

        self.run.info('Removed due to low hit score of %d' % self.min_hit_score, num_hits_below_hit_score)

        # we are done with this one too:
        del self.dicts['hits']

        removed_due_to_too_many_hits = 0
        for gene_callers_id in hits:
            counter_obj = hits[gene_callers_id]

            num_hits = len(counter_obj)
            if num_hits == 1:
                annotations_dict[gene_callers_id] = {'t_species': counter_obj.most_common()[0][0]}
            else:
                counter_obj_ordered = counter_obj.most_common()
                if counter_obj_ordered[0][1] > counter_obj_ordered[1][1]:
                    annotations_dict[gene_callers_id] = {'t_species': counter_obj.most_common()[0][0]}
                else:
                    removed_due_to_too_many_hits += 1
                    continue

            # add the other levels of taxonomy. this has to be done in a much better
            # way at some point, but here is a workaround to deal with this for now:
            annotations_dict[gene_callers_id]['t_genus'] = annotations_dict[gene_callers_id]['t_species'].split()[0]
            annotations_dict[gene_callers_id]['t_family'] = None
            annotations_dict[gene_callers_id]['t_order'] = None
            annotations_dict[gene_callers_id]['t_class'] = None
            annotations_dict[gene_callers_id]['t_phylum'] = None

        self.run.info('Removed due to too many competing hits', removed_due_to_too_many_hits)
        self.run.info('Final num hits', len(annotations_dict))

        genes_taxonomy, taxon_names = TaxonomyHelper(annotations_dict).get_genes_taxonomy_and_taxon_names_dicts()

        return (genes_taxonomy, taxon_names)
