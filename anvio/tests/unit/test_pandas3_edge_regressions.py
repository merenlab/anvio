"""Focused regressions for pandas compatibility edge cases."""

import unittest

import numpy as np
import pandas as pd

from anvio import constants
from anvio.genomictrnaseq import Affinitizer
from anvio.reactionnetwork import _select_network_rows
from anvio.taxonomyops import TaxonomyTree
from anvio.trnaseq import _anticodon_sequence_or_empty


class TestPandas3EdgeRegressions(unittest.TestCase):
    def test_contribution_stack_preserves_missing_cells_and_blank_anticodon(self):
        wide = pd.DataFrame(
            [[1.5, np.nan], [np.nan, 2.5]],
            index=pd.Index(['g1', 'g2'], name='gene_caller_id'),
            columns=['AAA', ''],
        )

        long = Affinitizer._stack_to_long(wide, 'sample-1', 'contribution')

        self.assertEqual(len(long), 4)
        self.assertIn('', set(long['anticodon']))
        self.assertEqual(long['trnaseq_sample_name'].unique().tolist(), ['sample-1'])
        self.assertEqual(long['contribution'].isna().sum(), 2)
        self.assertEqual(set(long['gene_caller_id']), {'g1', 'g2'})

    def test_missing_anticodon_sequence_becomes_blank(self):
        self.assertEqual(_anticodon_sequence_or_empty(np.nan), '')
        self.assertEqual(_anticodon_sequence_or_empty(pd.NA), '')
        self.assertEqual(_anticodon_sequence_or_empty('GTA'), 'GTA')

    def test_sparse_taxonomy_stops_at_first_missing_rank(self):
        entry = {level: None for level in constants.levels_of_taxonomy}
        entry['t_domain'] = 'Bacteria'
        entry['t_phylum'] = np.nan
        entry['t_class'] = 'must not appear below a missing rank'

        tree = TaxonomyTree([entry], max_taxonomic_level='t_class', use_unicode=False)

        domain = tree.root['children']['Bacteria']
        unknown_phylum = constants.levels_of_taxonomy_unknown['t_phylum']
        self.assertEqual(domain['children'][unknown_phylum]['count'], 1)
        self.assertNotIn('must not appear below a missing rank', str(domain))

    def test_reaction_ko_selection_keeps_membership_and_stable_order(self):
        annotations = pd.DataFrame({
            'accession': ['K00002', 'K99999', 'K00001', 'K00002'],
            'gene_callers_id': [20, 99, 10, 21],
        })

        selected = _select_network_rows(
            annotations, 'accession', {'K00002', 'K00001', 'K00003'})

        self.assertEqual(selected.index.tolist(), ['K00001', 'K00002', 'K00002'])
        self.assertEqual(selected['gene_callers_id'].tolist(), [10, 20, 21])
