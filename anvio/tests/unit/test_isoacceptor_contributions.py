# pylint: disable=line-too-long
"""
Tests for the per-isoacceptor leave-one-out contribution analysis in
`anvio.genomictrnaseq.Affinitizer.get_isoacceptor_contributions`.

The headline invariant is that the closed-form Δ(g, s, i) = β(g, s) - β₋ᵢ(g, s) the analysis
emits matches a brute-force refit of `scipy.stats.linregress` with isoacceptor `i` removed.
Aggregation correctness, NaN propagation through SE-normalized contributions, the leverage
floor on (1 - h_{ii}), and the function-mode vs. gene-mode row-index conventions are also
checked here.

These tests synthesize the inputs `Affinitizer.get_affinities` consumes and bypass `__init__`
(via `Affinitizer.__new__`) so no tRNA-seq or genomic database is required to run them.
"""

import unittest

import numpy as np
import pandas as pd
from scipy.stats import linregress
from types import SimpleNamespace

import anvio

from anvio.errors import ConfigError
from anvio.genomictrnaseq import Affinitizer


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


SAMPLES = ['s1', 's2', 's3']
ANTICODONS = ['ACA', 'ACG', 'ACC', 'ACT', 'AAA', 'AAG', 'AAT', 'AAC']
GENOMES = ['gA', 'gB']
SOURCES = ['KEGG', 'Pfam']


def _make_synthetic_inputs(rng, gene_mode=False, n_funcs=10):
    """
    Build (`isoacceptor_abund_ratios_df`, `isoacceptor_codon_weights_df`) for two genomes.

    Function mode uses two function sources per genome to exercise the function_source
    pivot-level handling; gene mode uses `gene_callers_id`.
    """
    rows = []
    for gen in GENOMES:
        for s in SAMPLES:
            for ac in ANTICODONS:
                rows.append({
                    'genome_name': gen,
                    'decoded_amino_acid': 'X',
                    'anticodon': ac,
                    'trnaseq_sample_name': s,
                    'abundance_ratio': float(np.exp(rng.normal())),
                })
    abund = pd.DataFrame(rows)

    if gene_mode:
        idx = pd.MultiIndex.from_tuples(
            [(gen, i) for gen in GENOMES for i in range(n_funcs)],
            names=['genome_name', 'gene_callers_id'])
    else:
        idx = pd.MultiIndex.from_tuples([
            (gen, src, f'K{i:05d}', f'f{i}')
            for gen in GENOMES for src in SOURCES for i in range(n_funcs)],
            names=['genome_name', 'function_source', 'function_accession', 'function_name'])
    weights = pd.DataFrame(
        rng.gamma(2.0, 1.0, size=(len(idx), len(ANTICODONS))),
        index=idx, columns=ANTICODONS)
    return abund, weights


def _make_affinitizer(gene_affinity, genome_info_dict=None):
    """Build a minimal `Affinitizer` instance that bypasses `__init__`'s DB requirements."""
    aff = Affinitizer.__new__(Affinitizer)
    aff.genome_info_dict = genome_info_dict or {gen: {} for gen in GENOMES}
    aff.run = SimpleNamespace(warning=lambda *a, **k: None, info_single=lambda *a, **k: None)
    aff.gene_affinity = gene_affinity
    return aff


class TestContributionMath(unittest.TestCase):
    """Verify the closed-form contribution matches brute-force LOO refits and that
    aggregations behave as plain pandas operations on the underlying long table."""

    def setUp(self):
        self.rng = np.random.default_rng(42)
        self.abund, self.weights = _make_synthetic_inputs(self.rng)
        self.aff = _make_affinitizer(gene_affinity=False)
        _, _, self.intermediates = self.aff.get_affinities(
            self.abund, self.weights, return_intermediates=True)
        self.contributions, self.sanity = self.aff.get_isoacceptor_contributions(
            self.intermediates)
        # Per-gene relative codon weights used to reconstruct y for brute-force refits.
        self.rel_weights = self.weights.div(self.weights.sum(axis=1), axis=0)


    def test_delta_raw_matches_brute_force_loo_refit(self):
        """For every (g, s, i) in the long-format Δ_raw table, removing isoacceptor `i` and
        refitting via `scipy.stats.linregress` should produce β₋ᵢ such that β - β₋ᵢ equals
        the reported Δ_raw within floating-point precision."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                for _, row in long_raw.iterrows():
                    inter = self.intermediates[(gen, row['trnaseq_sample_name'])]
                    gene_row = (gen, src, row['function_accession'], row['function_name'])
                    y_full = self.rel_weights.loc[gene_row, inter['anticodons']].values
                    i = inter['anticodons'].index(row['anticodon'])
                    refit = linregress(np.delete(inter['x'], i), np.delete(y_full, i))
                    brute_delta = inter['beta'].loc[gene_row] - refit.slope
                    self.assertTrue(
                        np.isclose(row['delta_raw'], brute_delta, atol=1e-10),
                        msg=f"({gene_row}, {row['trnaseq_sample_name']}, {row['anticodon']}): "
                            f"closed-form Δ={row['delta_raw']} vs. brute={brute_delta}")


    def test_delta_norm_equals_delta_raw_over_se(self):
        """Δ_norm should equal Δ_raw / SE(β) cell-for-cell, with both NaN where SE is NaN
        (n ≤ 2 or degenerate supply) -- and *not* infinite where SE is 0 (in this synthetic
        dataset SE > 0 everywhere, so we only verify the finite case)."""
        for gen in GENOMES:
            for src in SOURCES:
                lr = self.contributions[gen][src]['LONG-RAW']
                ln = self.contributions[gen][src]['LONG-NORM']
                # Both long tables share the same key columns; align and divide.
                key_cols = [c for c in lr.columns if c not in {'delta_raw'}]
                merged = lr.merge(ln, on=key_cols)
                for _, row in merged.iterrows():
                    inter = self.intermediates[(gen, row['trnaseq_sample_name'])]
                    gene_row = (gen, src, row['function_accession'], row['function_name'])
                    se = inter['se'].loc[gene_row]
                    expected = row['delta_raw'] / se if se != 0 else np.nan
                    self.assertTrue(
                        np.isclose(row['delta_norm'], expected, atol=1e-10, equal_nan=True),
                        msg=f"({gene_row}, {row['trnaseq_sample_name']}, {row['anticodon']}): "
                            f"Δ_norm={row['delta_norm']} vs. Δ_raw/SE={expected}")


    def test_per_sample_mean_matches_groupby(self):
        """`PER_SAMPLE-RAW-MEAN` should reproduce a direct groupby-mean over the long table."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                per_sample = self.contributions[gen][src]['PER_SAMPLE-RAW-MEAN']
                expected = (long_raw.groupby(
                    ['genome_name', 'function_source', 'trnaseq_sample_name', 'anticodon'],
                    dropna=False)['delta_raw'].mean().unstack('anticodon'))
                # `per_sample` is column-reindexed to the bucket's anticodon universe; align.
                expected = expected.reindex(columns=per_sample.columns)
                pd.testing.assert_frame_equal(per_sample, expected, check_names=True)


    def test_global_mean_collapses_over_both_dimensions(self):
        """`GLOBAL-RAW-MEAN` should reproduce a direct groupby-mean over (genome,
        function_source, anticodon)."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                glob = self.contributions[gen][src]['GLOBAL-RAW-MEAN']
                expected = (long_raw.groupby(
                    ['genome_name', 'function_source', 'anticodon'],
                    dropna=False)['delta_raw'].mean().unstack('anticodon'))
                expected = expected.reindex(columns=glob.columns)
                pd.testing.assert_frame_equal(glob, expected, check_names=True)


    def test_per_gene_abs_mean_matches_groupby(self):
        """`PER_GENE-RAW-ABS_MEAN` should reproduce |Δ|.mean() per (gene, anticodon)."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                per_gene = self.contributions[gen][src]['PER_GENE-RAW-ABS_MEAN']
                expected = (long_raw.groupby([
                    'genome_name', 'function_source', 'function_accession',
                    'function_name', 'anticodon'],
                    dropna=False)['delta_raw'].apply(
                        lambda v: v.abs().mean()).unstack('anticodon'))
                expected = expected.reindex(columns=per_gene.columns)
                pd.testing.assert_frame_equal(per_gene, expected, check_names=True)


    def test_function_mode_aggregations_include_function_source_level(self):
        """`function_source` must appear in the row index for every aggregated table in
        function mode -- otherwise cross-source concatenation by the writer produces
        duplicate (genome, sample) or (genome,) rows."""
        sample = self.contributions[GENOMES[0]][SOURCES[0]]
        for key in ['PER_SAMPLE-RAW-MEAN', 'PER_GENE-RAW-MEAN', 'GLOBAL-RAW-MEAN']:
            self.assertIn(
                'function_source', sample[key].index.names,
                msg=f"{key}: expected function_source in index, got {sample[key].index.names}")


class TestGeneMode(unittest.TestCase):
    """Gene-affinity mode has a different gene/function index structure (no function_source
    level). Verify the contribution outputs reflect that."""

    def setUp(self):
        rng = np.random.default_rng(7)
        self.abund, self.weights = _make_synthetic_inputs(rng, gene_mode=True)
        self.aff = _make_affinitizer(gene_affinity=True)
        _, _, self.intermediates = self.aff.get_affinities(
            self.abund, self.weights, return_intermediates=True)
        self.contributions, _ = self.aff.get_isoacceptor_contributions(self.intermediates)


    def test_single_source_bucket_named_genes(self):
        """In gene mode the source bucket is named 'genes' and is the only one."""
        for gen in GENOMES:
            self.assertEqual(list(self.contributions[gen].keys()), ['genes'])


    def test_aggregations_omit_function_source_level(self):
        """Gene mode's aggregations should NOT include function_source as a row-index level."""
        sample = self.contributions[GENOMES[0]]['genes']
        self.assertEqual(
            list(sample['PER_SAMPLE-RAW-MEAN'].index.names),
            ['genome_name', 'trnaseq_sample_name'])
        self.assertEqual(
            list(sample['PER_GENE-RAW-MEAN'].index.names),
            ['genome_name', 'gene_callers_id'])
        self.assertEqual(
            list(sample['GLOBAL-RAW-MEAN'].index.names),
            ['genome_name'])


class TestNumericalEdgeCases(unittest.TestCase):
    """Cases where the closed-form Δ would be ill-defined: high-leverage isoacceptors, SE
    exactly zero, and degenerate supply vectors."""

    def setUp(self):
        self.rng = np.random.default_rng(0)


    def test_high_leverage_isoacceptor_is_clipped_to_nan(self):
        """When one isoacceptor's log-ratio is far outside the spread of the others, its
        leverage h_{ii} approaches 1; (1 − h_{ii}) drops below `contribution_leverage_floor`
        and the corresponding Δ entries must be NaN."""
        # Stack 5 near-zero log-ratios and one extreme one so the extreme isoacceptor
        # carries ~all of Σ(x − x̄)² (h_ii ≈ 1).
        rows = []
        for s in SAMPLES:
            base = {'genome_name': 'gA', 'decoded_amino_acid': 'X', 'trnaseq_sample_name': s}
            for ac in ANTICODONS[:5]:
                rows.append({**base, 'anticodon': ac, 'abundance_ratio': 1.0 + 1e-6})
            # One isoacceptor with extreme abundance ratio.
            rows.append({**base, 'anticodon': ANTICODONS[5], 'abundance_ratio': 1e6})
        abund = pd.DataFrame(rows)
        idx = pd.MultiIndex.from_tuples(
            [('gA', 'KEGG', f'K{i:05d}', f'f{i}') for i in range(5)],
            names=['genome_name', 'function_source', 'function_accession', 'function_name'])
        weights = pd.DataFrame(
            self.rng.gamma(2.0, 1.0, size=(5, 6)), index=idx, columns=ANTICODONS[:6])
        aff = _make_affinitizer(gene_affinity=False, genome_info_dict={'gA': {}})

        _, _, inters = aff.get_affinities(abund, weights, return_intermediates=True)
        contribs, _ = aff.get_isoacceptor_contributions(inters)

        # Δ for the high-leverage isoacceptor should be NaN; Δ for the others should be
        # finite (the (1 − h_{ii}) factor for them is well above the floor).
        long_raw = contribs['gA']['KEGG']['LONG-RAW']
        ext_rows = long_raw[long_raw['anticodon'] == ANTICODONS[5]]
        self.assertGreater(len(ext_rows), 0)
        self.assertTrue(ext_rows['delta_raw'].isna().all(),
                        msg="Δ_raw for the high-leverage isoacceptor must be NaN")


    def test_perfect_fit_makes_delta_norm_nan_not_infinite(self):
        """When SE(β) is exactly zero (the residuals are zero -- a perfect line through the
        n points), Δ_norm = Δ / 0 must be normalized to NaN, not ±inf."""
        # Construct y = α + β·x exactly so residuals are zero and SE(β) = 0.
        x_supply = np.linspace(0.5, 4.0, len(ANTICODONS))  # arbitrary spread
        abund_ratios = np.power(2.0, x_supply)  # log2(abund_ratio) = x_supply
        rows = []
        for s in SAMPLES:
            for ac, ar in zip(ANTICODONS, abund_ratios):
                rows.append({
                    'genome_name': 'gA', 'decoded_amino_acid': 'X', 'anticodon': ac,
                    'trnaseq_sample_name': s, 'abundance_ratio': float(ar)})
        abund = pd.DataFrame(rows)
        # Make codon weights a perfect linear function of x_supply so the OLS fit has zero
        # residuals. Use the *relative* weights y = α + β·x; back-solve raw weights so the
        # row-sum normalization recovers those relative weights.
        intercept, slope = 0.05, 0.02
        relative_y = intercept + slope * x_supply
        relative_y = relative_y / relative_y.sum()  # normalize so the row sums to 1
        # After row-normalization in `get_affinities`, weights[ac] for this single-function
        # genome become relative_y[ac], producing a perfect linear y = a + b·x_supply fit.
        idx = pd.MultiIndex.from_tuples(
            [('gA', 'KEGG', 'K00000', 'f0')],
            names=['genome_name', 'function_source', 'function_accession', 'function_name'])
        weights = pd.DataFrame(relative_y[None, :], index=idx, columns=ANTICODONS)
        aff = _make_affinitizer(gene_affinity=False, genome_info_dict={'gA': {}})

        _, _, inters = aff.get_affinities(abund, weights, return_intermediates=True)
        # Confirm we got a zero-SE fit (within fp tolerance).
        inter = inters[('gA', 's1')]
        self.assertLess(float(inter['se'].iloc[0]), 1e-12,
                        msg="setup assumption: synthetic fit should have ~zero SE")

        contribs, _ = aff.get_isoacceptor_contributions(inters)
        long_norm = contribs['gA']['KEGG']['LONG-NORM']
        # No ±inf entries: the writer's `replace([inf, -inf], nan)` should have kicked in.
        self.assertFalse(
            np.isinf(long_norm['delta_norm']).any(),
            msg="Δ_norm must not contain ±inf when SE(β) = 0")


class TestValidationAndEmptyInput(unittest.TestCase):
    """`get_isoacceptor_contributions` validates its inputs up front so programmatic API
    users (or sanity-check-bypassed CLI invocations) get a clear error instead of a late
    one downstream."""

    def setUp(self):
        rng = np.random.default_rng(1)
        abund, weights = _make_synthetic_inputs(rng)
        self.aff = _make_affinitizer(gene_affinity=False)
        _, _, self.intermediates = self.aff.get_affinities(
            abund, weights, return_intermediates=True)


    def test_unknown_variant_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, 'unknown contribution variant'):
            self.aff.get_isoacceptor_contributions(self.intermediates, variants=['bogus'])


    def test_empty_variants_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, 'at least one contribution variant'):
            self.aff.get_isoacceptor_contributions(self.intermediates, variants=[])


    def test_aggregated_levels_without_statistics_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, 'at least one statistic'):
            self.aff.get_isoacceptor_contributions(
                self.intermediates, aggregations=['per_sample'], statistics=[])


    def test_empty_intermediates_returns_empty_results(self):
        contribs, sanity = self.aff.get_isoacceptor_contributions({})
        self.assertEqual(contribs, {})
        self.assertEqual(sanity, [])


    def test_long_only_request_skips_statistic_validation(self):
        """The aggregated-levels-without-statistics check shouldn't fire when only 'long' is
        requested -- 'long' has no statistic dimension."""
        contribs, _ = self.aff.get_isoacceptor_contributions(
            self.intermediates, aggregations=['long'], statistics=[])
        # Should produce only LONG-* keys, no PER_SAMPLE / PER_GENE / GLOBAL.
        for gen_outputs in contribs.values():
            for src_outputs in gen_outputs.values():
                for key in src_outputs:
                    self.assertTrue(
                        key.startswith('LONG-'),
                        msg=f"unexpected output key for long-only request: {key}")


if __name__ == '__main__':
    unittest.main()
