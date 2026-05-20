#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functional consensus analysis for gene clusters.

Given a gene cluster table and functional annotations from one or more sources,
this module evaluates within-cluster annotation consistency, cross-source
coherence, and optionally propagates the dominant annotation to unannotated
genes.

Works with any gene cluster table (pre-pan, post-pan, or entirely custom) and
any contigs databases, genomes storage, or pre-loaded annotation DataFrames.

Usage (Python API)
------------------
    from anvio.clusterannotation import GeneClusterFunctionalConsensus

    c = GeneClusterFunctionalConsensus(
        gene_clusters   = 'gene_clusters.txt',          # or a DataFrame
        annotations     = 'CONTIGS.db',                 # or a storage path / DataFrame
        annotation_sources = ['Pfam', 'COG20_FUNCTION'],
    )
    consensus_df, propagated_df = c.compute()
"""

import os
import re
import numpy as np
import pandas as pd
from itertools import combinations
from collections import Counter

import anvio
import anvio.db as anvio_db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['cvanni']


run = terminal.Run()
progress = terminal.Progress()


class GeneClusterFunctionalConsensus:
    """Evaluate functional consensus across genes within clusters.

    Parameters
    ----------
    gene_clusters : str or pd.DataFrame
        Path to a TSV or a DataFrame with columns:
        ``gene_callers_id``, ``cluster_id`` and optionally ``genome_name``.
        Accepts the canonical anvi'o gene-clusters-txt format
        (``genome_name`` / ``gene_caller_id`` / ``gene_cluster_name``) and
        normalises column names automatically.
    annotations : str, list, dict, or pd.DataFrame
        Annotation data in any of these forms:
        * ``pd.DataFrame`` with columns ``gene_callers_id``, ``source``,
          ``accession``, ``function``, ``e_value`` (+ ``genome_name`` for
          multi-genome data).
        * ``str`` — path to a single contigs database.
        * ``list`` of ``(genome_name, contigs_db_path)`` tuples.
        * ``dict`` mapping ``genome_name`` → ``contigs_db_path``.
        * ``str`` ending in ``.db`` that is a genomes storage database.
    annotation_sources : list of str, optional
        Sources to evaluate.  ``None`` uses every source present in the
        annotation data.
    min_coverage : float
        Minimum fraction of annotated genes to trust a cluster for a given
        source.  Default 0.3.
    min_jacc_sc : float
        Minimum coverage-scaled Jaccard score to classify a cluster as
        HIGH_CONSENSUS.  Default 0.6.
    min_cross_source_coherence : float
        Minimum cross-source keyword-Jaccard score before a PURE/HIGH_CONSENSUS
        cluster is downgraded to MIXED.  Default 0.3.
    propagate : bool
        If True, emit a propagated-annotations table for PURE and
        HIGH_CONSENSUS clusters.  Default True.
    pfam_data_dir : str, optional
        Path to the Pfam data directory (containing ``Pfam-A.clans.tsv``).
        Auto-detected from the anvi'o data directory when not provided.
    cogs_data_dir : str, optional
        Path to the COG data directory.  Auto-detected when not provided.
    """

    # Words that carry no functional meaning in annotation descriptions
    BIOLOGICAL_STOPWORDS = {
        'protein', 'domain', 'family', 'putative', 'hypothetical',
        'uncharacterized', 'predicted', 'conserved', 'possible', 'probable',
        'related', 'like', 'type', 'subunit', 'component', 'homolog',
        'homologue', 'superfamily', 'containing', 'binding', 'associated',
        'involved', 'function', 'unknown', 'activity', 'process', 'response',
        'factor', 'repeat', 'motif', 'region', 'site', 'group', 'class',
        'member', 'with', 'and', 'for', 'the', 'in', 'of', 'to', 'a', 'an',
        'by', 'gene', 'genes', 'bacterial', 'archaeal', 'eukaryotic',
    }

    # Strips standard N/C-terminal split-domain suffixes from Pfam function names
    _TERMINAL_RE = re.compile(
        r'[\s_\-]*(N|C|N[\-_]term(?:inal)?|C[\-_]term(?:inal)?|NT|CT|'
        r'Middle|Central|N_terminal|C_terminal)[\s_\-]*$',
        re.IGNORECASE,
    )

    # Matches EC numbers embedded in KEGG/KOfam descriptions: [EC:1.2.3.4]
    _EC_RE = re.compile(r'\[EC:([\d\.\-][\d\.\-\s,]+)\]')

    # Tokenise description text into words
    _TOKEN_RE = re.compile(r'[A-Za-z]{3,}')

    def __init__(
        self,
        gene_clusters,
        annotations,
        annotation_sources=None,
        min_coverage=0.3,
        min_jacc_sc=0.6,
        min_cross_source_coherence=0.3,
        propagate=True,
        pfam_data_dir=None,
        cogs_data_dir=None,
        run=run,
        progress=progress,
    ):
        self.run = run
        self.progress = progress

        self.min_coverage = min_coverage
        self.min_jacc_sc = min_jacc_sc
        self.min_cross_source_coherence = min_cross_source_coherence
        self.propagate = propagate

        self.progress.new('Initialising GeneClusterFunctionalConsensus')

        self.progress.update('Loading gene clusters …')
        self.gene_clusters_df = self._load_gene_clusters(gene_clusters)

        self.multi_genome = 'genome_name' in self.gene_clusters_df.columns
        self.gene_key_cols = (
            ['genome_name', 'gene_callers_id'] if self.multi_genome
            else ['gene_callers_id']
        )

        self.progress.update('Loading annotations …')
        self.annotations_df = self._load_annotations(annotations)

        # Determine which sources to evaluate
        available_sources = self.annotations_df['source'].unique().tolist()
        if annotation_sources:
            missing = [s for s in annotation_sources if s not in available_sources]
            if missing:
                raise ConfigError(
                    f"The following requested annotation sources were not found "
                    f"in the annotation data: {', '.join(missing)}. "
                    f"Available sources: {', '.join(available_sources)}"
                )
            self.annotation_sources = annotation_sources
        else:
            self.annotation_sources = available_sources

        self.run.info('Gene clusters', f"{self.gene_clusters_df['cluster_id'].nunique():,} clusters, "
                                       f"{len(self.gene_clusters_df):,} genes")
        self.run.info('Annotation sources', ', '.join(self.annotation_sources))
        self.run.info('Multi-genome mode', self.multi_genome)

        # Lazy-loaded hierarchy lookups
        self._pfam_clans = None           # accession -> clan_id
        self._cog_categories = None       # cog_id    -> [category_letter, ...]
        self._pfam_data_dir = pfam_data_dir
        self._cogs_data_dir = cogs_data_dir

        self.progress.end()

    # ------------------------------------------------------------------
    # Data loaders
    # ------------------------------------------------------------------

    def _load_gene_clusters(self, gene_clusters) -> pd.DataFrame:
        """Return a normalised DataFrame: gene_callers_id, cluster_id [, genome_name]."""
        if isinstance(gene_clusters, pd.DataFrame):
            df = gene_clusters.copy()
        elif isinstance(gene_clusters, str):
            filesnpaths.is_file_exists(gene_clusters)
            df = pd.read_csv(gene_clusters, sep='\t')
        else:
            raise ConfigError(
                "gene_clusters must be a file path (str) or a pandas DataFrame."
            )

        # Normalise column names from anvi'o gene-clusters-txt format
        rename_map = {
            'gene_caller_id': 'gene_callers_id',   # pan DB uses singular
            'gene_cluster_id': 'cluster_id',
            'gene_cluster_name': 'cluster_id',
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        required = {'gene_callers_id', 'cluster_id'}
        missing = required - set(df.columns)
        if missing:
            raise ConfigError(
                f"Gene clusters table is missing required columns: {', '.join(missing)}. "
                f"Expected at minimum: gene_callers_id (or gene_caller_id), cluster_id "
                f"(or gene_cluster_id / gene_cluster_name)."
            )

        df['gene_callers_id'] = df['gene_callers_id'].astype(int)
        return df

    def _load_annotations(self, annotations) -> pd.DataFrame:
        """Normalise any supported annotation input to a DataFrame."""
        if isinstance(annotations, pd.DataFrame):
            df = annotations.copy()
            self._validate_annotations_df(df)
            return df

        if isinstance(annotations, dict):
            annotations = list(annotations.items())

        if isinstance(annotations, list):
            frames = [
                self._load_from_contigs_db(path, genome_name=name)
                for name, path in annotations
            ]
            return pd.concat(frames, ignore_index=True)

        if isinstance(annotations, str):
            filesnpaths.is_file_exists(annotations)
            # Distinguish genomes storage from contigs database by db_type meta value
            db_type = self._get_db_type(annotations)
            if db_type == 'genomestorage':
                return self._load_from_genomes_storage(annotations)
            elif db_type == 'contigs':
                return self._load_from_contigs_db(annotations)
            else:
                raise ConfigError(
                    f"The database at '{annotations}' has type '{db_type}', "
                    f"but anvi-compute-functional-consensus expects either a "
                    f"contigs database or a genomes storage database."
                )

        raise ConfigError(
            "annotations must be a DataFrame, a contigs DB path, a genomes "
            "storage path, a list of (genome_name, contigs_db_path) tuples, "
            "or a dict mapping genome_name to contigs_db_path."
        )

    def _get_db_type(self, db_path: str) -> str:
        _db = anvio_db.DB(db_path, None, ignore_version=True)
        db_type = _db.get_meta_value('db_type', return_none_if_not_in_table=True)
        _db.disconnect()
        return db_type or 'unknown'

    def _load_from_contigs_db(self, db_path: str, genome_name=None) -> pd.DataFrame:
        filesnpaths.is_file_exists(db_path)
        _db = anvio_db.DB(db_path, None, ignore_version=True)
        df = _db.get_table_as_dataframe(
            t.gene_function_calls_table_name,
            columns_of_interest=['gene_callers_id', 'source', 'accession', 'function', 'e_value'],
            error_if_no_data=False,
        )
        _db.disconnect()

        if df is None or df.empty:
            self.run.warning(f"No functional annotations found in '{db_path}'.")
            df = pd.DataFrame(columns=['gene_callers_id', 'source', 'accession', 'function', 'e_value'])

        df['gene_callers_id'] = df['gene_callers_id'].astype(int)
        if genome_name:
            df.insert(0, 'genome_name', genome_name)
        return df

    def _load_from_genomes_storage(self, storage_path: str) -> pd.DataFrame:
        filesnpaths.is_file_exists(storage_path)
        _db = anvio_db.DB(storage_path, None, ignore_version=True)
        df = _db.get_table_as_dataframe(
            t.genome_gene_function_calls_table_name,
            columns_of_interest=['genome_name', 'gene_callers_id', 'source', 'accession', 'function', 'e_value'],
            error_if_no_data=False,
        )
        _db.disconnect()

        if df is None or df.empty:
            self.run.warning(f"No functional annotations found in genomes storage '{storage_path}'.")
            return pd.DataFrame(
                columns=['genome_name', 'gene_callers_id', 'source', 'accession', 'function', 'e_value']
            )

        df['gene_callers_id'] = df['gene_callers_id'].astype(int)
        return df

    def _validate_annotations_df(self, df: pd.DataFrame):
        required = {'gene_callers_id', 'source', 'accession', 'function'}
        missing = required - set(df.columns)
        if missing:
            raise ConfigError(
                f"Annotations DataFrame is missing required columns: {', '.join(missing)}."
            )

    # ------------------------------------------------------------------
    # Hierarchy loaders (lazy, graceful degradation)
    # ------------------------------------------------------------------

    def _load_pfam_clans(self):
        if self._pfam_clans is not None:
            return

        data_dir = self._pfam_data_dir or os.path.join(
            os.path.dirname(anvio.__file__), 'data/misc/Pfam'
        )
        clan_file = os.path.join(data_dir, 'Pfam-A.clans.tsv')

        if not filesnpaths.is_file_exists(clan_file, dont_raise=True):
            self.run.warning(
                "Pfam clan file not found; functionally_coherent sub-flag "
                "will not be set for Pfam annotations."
            )
            self._pfam_clans = {}
            return

        catalog = utils.get_TAB_delimited_file_as_dictionary(
            clan_file,
            column_names=['accession', 'clan', 'unk1', 'unk2', 'function'],
            no_header=True,
        )
        # accession may carry version suffix (PF00069.26) — strip it
        self._pfam_clans = {
            k.split('.')[0]: v['clan']
            for k, v in catalog.items()
            if v['clan']
        }
        self.run.info_single(
            f"Loaded Pfam clans for {len(self._pfam_clans):,} accessions.",
            nl_before=0,
        )

    def _load_cog_categories(self):
        if self._cog_categories is not None:
            return

        # Try to find a COG data directory via the same priority chain as cogs.py
        base_dir = (
            self._cogs_data_dir
            or os.environ.get('ANVIO_COG_DATA_DIR')
            or os.path.join(os.path.dirname(anvio.__file__), 'data/misc/COG')
        )

        # Walk versioned sub-directories (COG20, COG24, COG14 …)
        import glob
        cog_txt = None
        for version_dir in sorted(glob.glob(os.path.join(base_dir, '*COG*')), reverse=True):
            candidate = os.path.join(version_dir, 'COG.txt')
            if os.path.exists(candidate):
                cog_txt = candidate
                break

        if not cog_txt:
            self.run.warning(
                "COG data not found; functionally_coherent sub-flag will not "
                "be set for COG annotations."
            )
            self._cog_categories = {}
            return

        raw = utils.get_TAB_delimited_file_as_dictionary(
            cog_txt,
            column_names=['COG', 'categories', 'annotation'],
            no_header=True,
        )
        self._cog_categories = {
            cog: [c.strip() for c in v['categories'].split(',')]
            for cog, v in raw.items()
        }
        self.run.info_single(
            f"Loaded COG categories for {len(self._cog_categories):,} COGs.",
            nl_before=0,
        )

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _jaccard(a: frozenset, b: frozenset) -> float:
        """Jaccard similarity between two annotation sets."""
        if not a and not b:
            return 1.0
        if not a or not b:
            return 0.0
        return len(a & b) / len(a | b)

    @staticmethod
    def _normalize_terminal(name: str) -> str:
        """Strip N/C-terminal suffixes from a Pfam function name."""
        return GeneClusterFunctionalConsensus._TERMINAL_RE.sub('', name).strip()

    @staticmethod
    def _parse_multi_hit(value: str) -> list:
        """Split a !!!-joined COG multi-hit string into individual tokens."""
        if not value or value in ('None', 'nan'):
            return []
        return [v.strip() for v in str(value).split('!!!') if v.strip() and v.strip() not in ('None', 'nan')]

    @staticmethod
    def _median_pairwise_jaccard(annotation_sets: list) -> float:
        """Median pairwise Jaccard similarity over a list of frozensets."""
        n = len(annotation_sets)
        if n == 0:
            return 0.0
        if n == 1:
            return 1.0
        scores = [
            GeneClusterFunctionalConsensus._jaccard(annotation_sets[i], annotation_sets[j])
            for i, j in combinations(range(n), 2)
        ]
        return float(np.median(scores))

    def _extract_ec_numbers(self, text: str) -> frozenset:
        """Return a frozenset of EC numbers found in a description string."""
        matches = self._EC_RE.findall(text)
        ecs = set()
        for m in matches:
            for ec in re.split(r'[,\s]+', m):
                ec = ec.strip()
                if ec:
                    ecs.add(ec)
        return frozenset(ecs)

    def _extract_keywords(self, text: str) -> frozenset:
        """Return a frozenset of meaningful lowercase words from a description."""
        if not text or text in ('None', 'nan'):
            return frozenset()
        tokens = self._TOKEN_RE.findall(text.lower())
        return frozenset(t for t in tokens if t not in self.BIOLOGICAL_STOPWORDS)

    # ------------------------------------------------------------------
    # Functional coherence checks (Pfam clan / COG category)
    # ------------------------------------------------------------------

    def _check_functional_coherence(self, accessions: set, source: str) -> bool:
        """Return True if all accessions in a cluster map to the same functional class."""
        if not accessions:
            return False

        if source == 'Pfam':
            self._load_pfam_clans()
            if not self._pfam_clans:
                return False
            clans = {self._pfam_clans.get(a.split('.')[0]) for a in accessions}
            clans.discard(None)
            return len(clans) == 1

        if source.startswith('COG'):
            self._load_cog_categories()
            if not self._cog_categories:
                return False
            # Collect all category letters across all accessions
            all_cats = set()
            for acc in accessions:
                cats = self._cog_categories.get(acc, [])
                all_cats.update(cats)
            return len(all_cats) == 1

        return False

    # ------------------------------------------------------------------
    # Per-cluster, per-source statistics
    # ------------------------------------------------------------------

    def _build_annotation_sets(
        self, source_rows: pd.DataFrame, normalize_terminals: bool = False
    ) -> dict:
        """Return {gene_key: frozenset(accessions)} for genes with annotations.

        Works uniformly for both Pfam (multiple rows per gene) and COG
        (single !!!-joined row per gene).
        """
        annotation_sets = {}

        group_cols = [c for c in self.gene_key_cols if c in source_rows.columns]
        if not group_cols:
            return annotation_sets

        for key, group in source_rows.groupby(group_cols):
            accs = set()
            for _, row in group.iterrows():
                for acc in self._parse_multi_hit(str(row['accession'])):
                    if normalize_terminals and 'function' in row:
                        # normalize is applied to the description, we compare
                        # the normalised function text as the "key" rather than
                        # the accession — so we substitute the acc with the
                        # normalised description to unify N/C variants
                        func = self._normalize_terminal(
                            self._parse_multi_hit(str(row['function']))[0]
                            if self._parse_multi_hit(str(row['function'])) else acc
                        )
                        accs.add(func)
                    else:
                        accs.add(acc)
            if accs:
                annotation_sets[key] = frozenset(accs)

        return annotation_sets

    def _compute_source_stats(
        self,
        cluster_id: str,
        source: str,
        n_genes: int,
        cluster_annotations: pd.DataFrame,
    ) -> dict:
        """Compute all statistics for one cluster × source combination."""

        base = {
            'cluster_id': cluster_id,
            'source': source,
            'n_genes': n_genes,
            'n_annotated': 0,
            'coverage': 0.0,
            'dominant_accession': None,
            'dominant_function': None,
            'jacc_median_raw': None,
            'jacc_median_sc': None,
            'is_multi_domain': False,
            'architecturally_consistent': False,
            'functionally_coherent': False,
            'n_sources_evaluated': None,
            'cross_source_coherence': None,
            'cross_source_method': None,
            'cross_source_conflict': False,
            'classification': 'LOW_EVIDENCE',
        }

        # Rows for this source that have a non-null accession
        src_rows = cluster_annotations[
            (cluster_annotations['source'] == source)
            & cluster_annotations['accession'].notna()
        ]

        # Count annotated genes (unique gene keys)
        group_cols = [c for c in self.gene_key_cols if c in src_rows.columns]
        if src_rows.empty or not group_cols:
            return base

        n_annotated = src_rows[group_cols].drop_duplicates().shape[0]
        coverage = n_annotated / n_genes if n_genes else 0.0

        base.update({'n_annotated': n_annotated, 'coverage': round(coverage, 4)})

        if coverage < self.min_coverage:
            return base

        # Raw accession sets — used only for multi-domain and dominant-accession logic
        raw_annotation_sets = self._build_annotation_sets(src_rows, normalize_terminals=False)
        raw_sets_list = list(raw_annotation_sets.values())

        # Terminal-normalized sets — used for Jaccard and architectural consistency
        # (N/C-terminal Pfam split variants are unified before comparison)
        norm_annotation_sets = self._build_annotation_sets(src_rows, normalize_terminals=True)
        sets_list = list(norm_annotation_sets.values())

        # Dominant accession (most frequent individual accession)
        acc_counter: Counter = Counter()
        func_map: dict = {}
        for _, row in src_rows.iterrows():
            accs = self._parse_multi_hit(str(row['accession']))
            funcs = self._parse_multi_hit(str(row.get('function', '')))
            for i, acc in enumerate(accs):
                acc_counter[acc] += 1
                if i < len(funcs):
                    func_map[acc] = funcs[i]

        dominant_accession = acc_counter.most_common(1)[0][0] if acc_counter else None
        dominant_function = func_map.get(dominant_accession) if dominant_accession else None

        # Pairwise Jaccard (on normalized sets)
        jacc_raw = self._median_pairwise_jaccard(sets_list)
        jacc_sc = round(jacc_raw * coverage, 4)
        jacc_raw = round(jacc_raw, 4)

        # Sub-flags
        is_multi_domain = any(len(s) > 1 for s in raw_sets_list)
        architecturally_consistent = len(set(sets_list)) == 1

        # Functional coherence via clan / category hierarchy (raw accessions)
        all_accessions: set = set()
        for s in raw_sets_list:
            all_accessions.update(s)
        functionally_coherent = self._check_functional_coherence(all_accessions, source)

        # Classification
        if jacc_raw == 1.0:
            classification = 'PURE'
        elif jacc_sc >= self.min_jacc_sc:
            classification = 'HIGH_CONSENSUS'
        else:
            classification = 'MIXED'

        base.update({
            'dominant_accession': dominant_accession,
            'dominant_function': dominant_function,
            'jacc_median_raw': jacc_raw,
            'jacc_median_sc': jacc_sc,
            'is_multi_domain': is_multi_domain,
            'architecturally_consistent': architecturally_consistent,
            'functionally_coherent': functionally_coherent,
            'classification': classification,
        })
        return base

    # ------------------------------------------------------------------
    # Cross-source coherence
    # ------------------------------------------------------------------

    def _cross_source_coherence(self, evaluable_results: list) -> tuple:
        """Compare dominant annotation descriptions across sources.

        Returns (score, method, conflict) where:
        - score  : float 0-1 or None
        - method : 'ec' | 'keyword' | None
        - conflict : bool
        """
        descriptions = []
        for r in evaluable_results:
            text = r.get('dominant_function') or ''
            if text and text not in ('None', 'nan', 'unknown'):
                descriptions.append(text)

        if len(descriptions) < 2:
            return None, None, False

        # Tier 1: EC number overlap
        ec_sets = [self._extract_ec_numbers(d) for d in descriptions]
        if any(ec_sets):
            # At least one source has EC numbers
            ec_pairs = [
                (i, j) for i, j in combinations(range(len(ec_sets)), 2)
                if ec_sets[i] and ec_sets[j]
            ]
            if ec_pairs:
                scores = [
                    self._jaccard(ec_sets[i], ec_sets[j])
                    for i, j in ec_pairs
                ]
                score = float(np.median(scores))
                conflict = score < self.min_cross_source_coherence
                return round(score, 4), 'ec', conflict

        # Tier 2: keyword Jaccard
        kw_sets = [self._extract_keywords(d) for d in descriptions]
        usable = [s for s in kw_sets if s]
        if len(usable) < 2:
            return None, None, False

        scores = [
            self._jaccard(usable[i], usable[j])
            for i, j in combinations(range(len(usable)), 2)
        ]
        score = float(np.median(scores))
        conflict = score < self.min_cross_source_coherence
        return round(score, 4), 'keyword', conflict

    # ------------------------------------------------------------------
    # Propagation
    # ------------------------------------------------------------------

    def _build_propagated_table(
        self, consensus_df: pd.DataFrame, merged: pd.DataFrame
    ) -> pd.DataFrame:
        """Build annotation propagation table for PURE and HIGH_CONSENSUS clusters."""
        rows = []

        propagate_clusters = consensus_df[
            consensus_df['classification'].isin(['PURE', 'HIGH_CONSENSUS'])
            & ~consensus_df['cross_source_conflict']
        ][['cluster_id', 'source', 'dominant_accession', 'dominant_function']].copy()

        for _, rec in propagate_clusters.iterrows():
            cluster_id = rec['cluster_id']
            source = rec['source']

            cluster_genes = merged[
                merged['cluster_id'] == cluster_id
            ][self.gene_key_cols].drop_duplicates()

            annotated_genes = merged[
                (merged['cluster_id'] == cluster_id)
                & (merged['source'] == source)
                & merged['accession'].notna()
            ][self.gene_key_cols].drop_duplicates()

            annotated_keys = set(
                annotated_genes.itertuples(index=False, name=None)
            )

            for gene_row in cluster_genes.itertuples(index=False, name=None):
                gene_dict = dict(zip(self.gene_key_cols, gene_row))
                is_propagated = gene_row not in annotated_keys

                out = {
                    **gene_dict,
                    'cluster_id': cluster_id,
                    'source': source,
                    'accession': rec['dominant_accession'],
                    'function': rec['dominant_function'],
                    'is_propagated': is_propagated,
                }
                rows.append(out)

        return pd.DataFrame(rows) if rows else pd.DataFrame(
            columns=self.gene_key_cols + [
                'cluster_id', 'source', 'accession', 'function', 'is_propagated'
            ]
        )

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def compute(self) -> tuple:
        """Run the full consensus pipeline.

        Returns
        -------
        consensus_df : pd.DataFrame
            One row per cluster × source combination with all statistics and
            the final classification (LOW_EVIDENCE / PURE / HIGH_CONSENSUS /
            MIXED).
        propagated_df : pd.DataFrame
            One row per gene × source for clusters classified as PURE or
            HIGH_CONSENSUS.  ``is_propagated=True`` marks genes that received
            an annotation by propagation (were unannotated before).
        """
        self.progress.new('Computing functional consensus')

        # ---- Merge gene clusters with annotations (left join) ----------
        self.progress.update('Merging clusters with annotations …')

        src_df = self.annotations_df[
            self.annotations_df['source'].isin(self.annotation_sources)
        ]
        merged = self.gene_clusters_df.merge(src_df, on=self.gene_key_cols, how='left')

        # Pre-compute cluster sizes from the original clusters table
        cluster_sizes = self.gene_clusters_df['cluster_id'].value_counts().to_dict()

        n_clusters = len(cluster_sizes)
        results = []

        # ---- Per-cluster loop ------------------------------------------
        self.progress.update(f"Processing {n_clusters:,} clusters …")

        for cluster_id, cluster_data in merged.groupby('cluster_id'):
            n_genes = cluster_sizes.get(cluster_id, 0)
            cluster_source_results = []

            for source in self.annotation_sources:
                r = self._compute_source_stats(
                    cluster_id, source, n_genes, cluster_data
                )
                cluster_source_results.append(r)

            # ---- Cross-source coherence --------------------------------
            evaluable = [
                r for r in cluster_source_results
                if r['classification'] != 'LOW_EVIDENCE'
            ]
            n_sources_evaluated = len(evaluable)

            if n_sources_evaluated >= 2:
                coherence, method, conflict = self._cross_source_coherence(evaluable)
            else:
                coherence, method, conflict = None, None, False

            for r in cluster_source_results:
                r['n_sources_evaluated'] = n_sources_evaluated
                r['cross_source_coherence'] = coherence
                r['cross_source_method'] = method
                r['cross_source_conflict'] = conflict
                if conflict and r['classification'] in ('PURE', 'HIGH_CONSENSUS'):
                    r['classification'] = 'MIXED'

            results.extend(cluster_source_results)

        self.progress.end()

        consensus_df = pd.DataFrame(results)

        # ---- Propagation -----------------------------------------------
        if self.propagate and not consensus_df.empty:
            self.progress.new('Building propagated annotations table')
            self.progress.update('…')
            propagated_df = self._build_propagated_table(consensus_df, merged)
            self.progress.end()
        else:
            propagated_df = pd.DataFrame(
                columns=self.gene_key_cols + [
                    'cluster_id', 'source', 'accession', 'function', 'is_propagated'
                ]
            )

        # ---- Summary ---------------------------------------------------
        if not consensus_df.empty:
            counts = consensus_df['classification'].value_counts()
            self.run.info_single(
                "Classification summary (across all cluster × source pairs):",
                nl_before=1,
            )
            for label, n in counts.items():
                self.run.info_single(f"  {label}: {n:,}")

            if self.propagate:
                n_prop = int(propagated_df['is_propagated'].sum()) if not propagated_df.empty else 0
                self.run.info('Propagated annotations', f"{n_prop:,} genes received a propagated annotation")

        return consensus_df, propagated_df
