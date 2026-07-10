# pylint: disable=line-too-long
"""
Tests for the per-isoacceptor contribution decomposition in
`anvio.genomictrnaseq.Affinitizer.get_isoacceptor_contributions`.

The affinity is the uncentered, demand-weighted log-supply index
    affinity(g, s) = Σ_i demand(g, i) · log2(supply(i, s) / reference(i)),
so each isoacceptor's contribution is `demand(g, i) · log_supply(i)` and the contributions sum to
the affinity by construction. These tests check that headline invariant, the cell-level
contribution definition, SE-normalization, aggregation correctness, and the function-mode vs.
gene-mode row-index conventions.

These tests synthesize the inputs `Affinitizer.get_affinities` consumes and bypass `__init__`
(via `Affinitizer.__new__`) so no tRNA-seq or genomic database is required to run them.
"""

import unittest

import numpy as np
import pandas as pd
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
    pivot-level handling; gene mode uses `gene_caller_id`.
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
                    'log_supply_variance': float(rng.uniform(0.01, 0.5)),
                })
    abund = pd.DataFrame(rows)

    if gene_mode:
        idx = pd.MultiIndex.from_tuples(
            [(gen, i) for gen in GENOMES for i in range(n_funcs)],
            names=['genome_name', 'gene_caller_id'])
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


class TestAffinityIndex(unittest.TestCase):
    """The affinity itself is the demand-weighted sum of log-supply."""

    def setUp(self):
        self.rng = np.random.default_rng(42)
        self.abund, self.weights = _make_synthetic_inputs(self.rng)
        self.aff = _make_affinitizer(gene_affinity=False)
        self.affinities, self.stderrs, self.intermediates = self.aff.get_affinities(
            self.abund, self.weights, return_intermediates=True)

    def test_affinity_is_demand_weighted_log_supply(self):
        """affinity(g, s) == Σ_i relative_weight(g, i) · log2(ratio(i, s))."""
        relative_weights = self.weights.div(self.weights.sum(axis=1), axis=0)
        for (gen, sample), inter in self.intermediates.items():
            log_supply = inter['log_supply']  # Series indexed by anticodon
            for gene_row in inter['beta'].index:
                expected = float(
                    (relative_weights.loc[gene_row, log_supply.index] * log_supply).sum())
                self.assertTrue(
                    np.isclose(inter['beta'].loc[gene_row], expected, atol=1e-12),
                    msg=f"{gene_row} @ {sample}: β={inter['beta'].loc[gene_row]} vs {expected}")

    def test_stderr_is_demand_sq_weighted_variance(self):
        """SE(g, s) == sqrt(Σ_i relative_weight(g, i)² · Var(x_i, s)), with the per-isoacceptor
        log-supply variances supplied alongside the abundance ratios."""
        relative_weights = self.weights.div(self.weights.sum(axis=1), axis=0)
        for gen in GENOMES:
            for sample in SAMPLES:
                sub = self.abund[(self.abund['genome_name'] == gen)
                                 & (self.abund['trnaseq_sample_name'] == sample)]
                variance = sub.set_index('anticodon')['log_supply_variance']
                for gene_row in self.intermediates[(gen, sample)]['beta'].index:
                    w = relative_weights.loc[gene_row, variance.index].values
                    expected = float(np.sqrt(((w ** 2) * variance.values).sum()))
                    actual = self.stderrs.loc[gene_row, sample]
                    self.assertTrue(
                        np.isclose(actual, expected, atol=1e-12),
                        msg=f"{gene_row} @ {sample}: SE={actual} vs {expected}")


class TestContributionMath(unittest.TestCase):
    """Per-isoacceptor contributions sum to the affinity and match their definition; aggregations
    behave as plain pandas operations on the underlying long table."""

    def setUp(self):
        self.rng = np.random.default_rng(42)
        self.abund, self.weights = _make_synthetic_inputs(self.rng)
        self.relative_weights = self.weights.div(self.weights.sum(axis=1), axis=0)
        self.aff = _make_affinitizer(gene_affinity=False)
        _, _, self.intermediates = self.aff.get_affinities(
            self.abund, self.weights, return_intermediates=True)
        self.contributions = self.aff.get_isoacceptor_contributions(self.intermediates)


    def test_contribution_sums_to_affinity(self):
        """For every (g, s), Σ_i contribution(g, s, i) must equal the affinity β(g, s). This is the
        defining property of the decomposition (it holds exactly, by construction)."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                summed = long_raw.groupby(
                    ['function_accession', 'function_name', 'trnaseq_sample_name'],
                    dropna=False)['contribution'].sum()
                for (acc, name, sample), contribution_sum in summed.items():
                    beta = self.intermediates[(gen, sample)]['beta'].loc[(gen, src, acc, name)]
                    self.assertTrue(
                        np.isclose(contribution_sum, beta, atol=1e-12, equal_nan=True),
                        msg=f"({gen}, {src}, {acc}, {sample}): Σ contribution = {contribution_sum} "
                            f"vs. β = {beta}")


    def test_contribution_equals_demand_times_log_supply(self):
        """Cell-level: contribution(g, s, i) == demand(g, i) · log_supply(i, s)."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                for row in long_raw.itertuples(index=False):
                    gene_row = (gen, src, row.function_accession, row.function_name)
                    log_supply = self.intermediates[(gen, row.trnaseq_sample_name)]['log_supply']
                    expected = self.relative_weights.loc[gene_row, row.anticodon] * \
                        log_supply[row.anticodon]
                    self.assertTrue(
                        np.isclose(row.contribution, expected, atol=1e-12),
                        msg=f"({gene_row}, {row.trnaseq_sample_name}, {row.anticodon}): "
                            f"{row.contribution} vs. {expected}")


    def test_contribution_norm_equals_contribution_over_se(self):
        """contribution_norm == contribution / SE(affinity) cell-for-cell, using the count-based
        standard errors computed by `get_affinities`."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                long_norm = self.contributions[gen][src]['LONG-NORM']
                key_cols = [c for c in long_raw.columns if c != 'contribution']
                merged = long_raw.merge(long_norm, on=key_cols)
                for row in merged.itertuples(index=False):
                    se = self.intermediates[(gen, row.trnaseq_sample_name)]['se'].loc[
                        (gen, src, row.function_accession, row.function_name)]
                    expected = row.contribution / se
                    self.assertTrue(
                        np.isclose(row.contribution_norm, expected, atol=1e-9, equal_nan=True),
                        msg=f"norm={row.contribution_norm} vs contribution/SE={expected}")


    def test_per_sample_mean_matches_groupby(self):
        """`PER_SAMPLE-RAW-MEAN` should reproduce a direct groupby-mean over the long table."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                per_sample = self.contributions[gen][src]['PER_SAMPLE-RAW-MEAN']
                expected = (long_raw.groupby(
                    ['genome_name', 'function_source', 'trnaseq_sample_name', 'anticodon'],
                    dropna=False)['contribution'].mean().unstack('anticodon'))
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
                    dropna=False)['contribution'].mean().unstack('anticodon'))
                expected = expected.reindex(columns=glob.columns)
                pd.testing.assert_frame_equal(glob, expected, check_names=True)


    def test_per_gene_abs_mean_matches_groupby(self):
        """`PER_GENE-RAW-ABS_MEAN` should reproduce |contribution|.mean() per (gene, anticodon)."""
        for gen in GENOMES:
            for src in SOURCES:
                long_raw = self.contributions[gen][src]['LONG-RAW']
                per_gene = self.contributions[gen][src]['PER_GENE-RAW-ABS_MEAN']
                expected = (long_raw.groupby([
                    'genome_name', 'function_source', 'function_accession',
                    'function_name', 'anticodon'],
                    dropna=False)['contribution'].apply(
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
        self.contributions = self.aff.get_isoacceptor_contributions(self.intermediates)


    def test_single_source_bucket_named_genes(self):
        """In gene mode the source bucket is named 'genes' and is the only one."""
        for gen in GENOMES:
            self.assertEqual(list(self.contributions[gen].keys()), ['genes'])


    def test_contribution_sums_to_affinity(self):
        """Σ_i contribution == β must hold in gene mode too."""
        for gen in GENOMES:
            long_raw = self.contributions[gen]['genes']['LONG-RAW']
            summed = long_raw.groupby(
                ['gene_caller_id', 'trnaseq_sample_name'], dropna=False)['contribution'].sum()
            for (gcid, sample), contribution_sum in summed.items():
                beta = self.intermediates[(gen, sample)]['beta'].loc[(gen, gcid)]
                self.assertTrue(np.isclose(contribution_sum, beta, atol=1e-12, equal_nan=True))


    def test_aggregations_omit_function_source_level(self):
        """Gene mode's aggregations should NOT include function_source as a row-index level."""
        sample = self.contributions[GENOMES[0]]['genes']
        self.assertEqual(
            list(sample['PER_SAMPLE-RAW-MEAN'].index.names),
            ['genome_name', 'trnaseq_sample_name'])
        self.assertEqual(
            list(sample['PER_GENE-RAW-MEAN'].index.names),
            ['genome_name', 'gene_caller_id'])
        self.assertEqual(
            list(sample['GLOBAL-RAW-MEAN'].index.names),
            ['genome_name'])


class TestMissingAndEdgeCases(unittest.TestCase):
    """Index-specific edge cases: isoacceptors with no codon demand are dropped, and SE == 0
    yields NaN (not infinite) normalized contributions."""

    def setUp(self):
        self.rng = np.random.default_rng(0)


    def test_isoacceptor_without_codon_weights_is_dropped(self):
        """An isoacceptor present in the abundance table but absent from the codon-weight columns
        carries no demand and must be excluded from contributions and from the affinity."""
        abund, weights = _make_synthetic_inputs(self.rng)
        # Add an extra isoacceptor to the abundance table that has no codon-weight column.
        extra_rows = []
        for gen in GENOMES:
            for s in SAMPLES:
                extra_rows.append({
                    'genome_name': gen, 'decoded_amino_acid': 'X', 'anticodon': 'GGG',
                    'trnaseq_sample_name': s, 'abundance_ratio': 2.0,
                    'log_supply_variance': 0.1})
        abund = pd.concat([abund, pd.DataFrame(extra_rows)], ignore_index=True)
        aff = _make_affinitizer(gene_affinity=False)
        _, _, inters = aff.get_affinities(abund, weights, return_intermediates=True)
        for inter in inters.values():
            self.assertNotIn('GGG', inter['anticodons'])
        contribs = aff.get_isoacceptor_contributions(inters)
        long_raw = contribs['gA']['KEGG']['LONG-RAW']
        self.assertNotIn('GGG', set(long_raw['anticodon']))


    def test_zero_se_makes_norm_nan_not_infinite(self):
        """When an injected SE is exactly zero, contribution_norm = contribution / 0 must be
        normalized to NaN, not +/-inf."""
        abund, weights = _make_synthetic_inputs(self.rng)
        aff = _make_affinitizer(gene_affinity=False)
        _, _, inters = aff.get_affinities(abund, weights, return_intermediates=True)
        # Inject SE = 0 for every gene/function in every (genome, sample).
        zeroed = {}
        for key, inter in inters.items():
            inter = dict(inter)
            inter['se'] = pd.Series(0.0, index=inter['beta'].index)
            zeroed[key] = inter
        contribs = aff.get_isoacceptor_contributions(zeroed)
        long_norm = contribs['gA']['KEGG']['LONG-NORM']
        self.assertFalse(np.isinf(long_norm['contribution_norm']).any())


class TestValidationAndEmptyInput(unittest.TestCase):
    """`get_isoacceptor_contributions` validates its inputs up front so programmatic API
    users (or sanity-check-bypassed CLI invocations) get a clear error instead of a late one."""

    def setUp(self):
        rng = np.random.default_rng(1)
        abund, weights = _make_synthetic_inputs(rng)
        self.aff = _make_affinitizer(gene_affinity=False)
        _, _, self.intermediates = self.aff.get_affinities(
            abund, weights, return_intermediates=True)


    def test_unknown_variant_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, 'unknown contribution variant'):
            self.aff.get_isoacceptor_contributions(self.intermediates, variants=['bogus'])


    def test_removed_residual_variant_is_rejected(self):
        """'residual' is no longer a valid variant for the index decomposition."""
        with self.assertRaisesRegex(ConfigError, 'unknown contribution variant'):
            self.aff.get_isoacceptor_contributions(self.intermediates, variants=['residual'])


    def test_empty_variants_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, 'at least one contribution variant'):
            self.aff.get_isoacceptor_contributions(self.intermediates, variants=[])


    def test_aggregated_levels_without_statistics_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, 'at least one statistic'):
            self.aff.get_isoacceptor_contributions(
                self.intermediates, aggregations=['per_sample'], statistics=[])


    def test_empty_intermediates_returns_empty_results(self):
        contribs = self.aff.get_isoacceptor_contributions({})
        self.assertEqual(contribs, {})


    def test_long_only_request_skips_statistic_validation(self):
        """The aggregated-levels-without-statistics check shouldn't fire when only 'long' is
        requested -- 'long' has no statistic dimension."""
        contribs = self.aff.get_isoacceptor_contributions(
            self.intermediates, aggregations=['long'], statistics=[])
        for gen_outputs in contribs.values():
            for src_outputs in gen_outputs.values():
                for key in src_outputs:
                    self.assertTrue(
                        key.startswith('LONG-'),
                        msg=f"unexpected output key for long-only request: {key}")


class TestChooseAmbiguousSeedGenome(unittest.TestCase):
    """
    Unit tests for `Affinitizer._choose_ambiguous_seed_genome`, the genome-choosing logic of
    `--seed-assignment ambiguous_choose`. Each input is {sample: {genome: unambiguous coverage}}.
    """

    @staticmethod
    def choose(sample_genome_coverages, min_coverage_ratio=5):
        return Affinitizer._choose_ambiguous_seed_genome(
            sample_genome_coverages, min_coverage_ratio)

    def test_single_sample_clear_winner(self):
        self.assertEqual(self.choose({'s1': {'gA': 100, 'gB': 10}}), 'gA')

    def test_ratio_not_met_drops(self):
        # 100 / 30 = 3.33 < 5
        self.assertIsNone(self.choose({'s1': {'gA': 100, 'gB': 30}}))

    def test_runner_up_zero_chooses_leader(self):
        self.assertEqual(self.choose({'s1': {'gA': 50, 'gB': 0}}), 'gA')

    def test_single_candidate_chosen(self):
        self.assertEqual(self.choose({'s1': {'gA': 50}}), 'gA')

    def test_agreement_across_samples(self):
        self.assertEqual(
            self.choose({'s1': {'gA': 100, 'gB': 10}, 's2': {'gA': 80, 'gB': 5}}), 'gA')

    def test_disagreement_across_samples_drops(self):
        self.assertIsNone(
            self.choose({'s1': {'gA': 100, 'gB': 10}, 's2': {'gA': 10, 'gB': 100}}))

    def test_uninformative_sample_is_ignored(self):
        # No candidate has unambiguous coverage in s1; the decision comes from s2.
        self.assertEqual(
            self.choose({'s1': {'gA': 0, 'gB': 0}, 's2': {'gA': 100, 'gB': 10}}), 'gA')

    def test_no_informative_sample_drops(self):
        self.assertIsNone(self.choose({'s1': {'gA': 0, 'gB': 0}}))
        self.assertIsNone(self.choose({}))

    def test_one_failing_sample_drops_whole_seed(self):
        # s1 is decisive for gA, but s2 fails the ratio (100 / 30 < 5), so the seed is unassignable.
        self.assertIsNone(
            self.choose({'s1': {'gA': 100, 'gB': 10}, 's2': {'gA': 100, 'gB': 30}}))

    def test_tie_drops(self):
        # ratio = 1 < 5
        self.assertIsNone(self.choose({'s1': {'gA': 50, 'gB': 50}}))


if __name__ == '__main__':
    unittest.main()
