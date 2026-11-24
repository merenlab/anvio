"""
Lazy SQL-backed access to the `variable_nucleotides` table of a profile.db.

Historically, code that touched per-sample SNVs across a profile database (most
prominently `anvi-report-dgrs` on the DGR_Finder branch) loaded the entire
table into a single pandas DataFrame up-front. On typical anvi'o data this is
fine. On the read-recruitment-scale profile databases that pair this module's
arrival (hundreds of millions of SNV rows in a single merged profile.db) it is
not: the materialized frame is tens to hundreds of gigabytes and even before
the analysis starts, the parent process has saturated memory.

The observation that motivates `SNVAccessor` is that every real read site is
actually localized. Window-scan algorithms iterate per (split, sample); BLAST
hit filtering needs a single (contig, position-range); variability profiling
needs (sample, contig, position-range). None of these need the full table at
once.

`SNVAccessor` exposes exactly those three patterns. They are fast only when the
`variable_nucleotides(split_name)` index exists; that index is created on demand
via the `anvi-index-table` program rather than by default (building it on every
profile db would force an expensive migration nobody asked for). Without the
index the same queries still return correct results, but degrade to full table
scans -- so `warn_or_raise_on_missing_index()` lets a consumer decide, based on
db size, whether to proceed or to stop and point the user at `anvi-index-table`.
Returned DataFrames use the same downcast dtypes that the prior eager-load code
produced, so downstream code can adopt the accessor without changes.
"""

import os

import numpy as np
import pandas as pd

import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.constants import nucleotides


# Default columns returned by SNV queries. Mirrors what the historical
# `init_snv_table` in dgrs.py selected, plus the derived `contig_name`
# (computed in `_categorize` from `split_name`).
DEFAULT_COLUMNS_OF_INTEREST = [
    'sample_id', 'split_name', 'pos_in_contig', 'base_pos_in_codon',
    'departure_from_reference', 'reference',
] + list(nucleotides)


# Cap on how many split names go into a single `split_name IN (...)` clause.
# Far below any SQLite parameter or statement-length limit, but small enough
# that even a pathologically large bin (tens of thousands of splits) stays
# responsive.
SPLIT_IN_LIST_BATCH = 5_000


class SNVAccessor:
    """Per-split / per-region SQL-backed accessor for `variable_nucleotides`.

    Constructed against a profile.db path and a `departure_from_reference`
    threshold. All queries apply the threshold so callers always see the same
    filter the historical eager-load applied.

    Designed to be picklable across `multiprocessing.spawn` workers: only the
    path and filter scalars survive pickling; the SQLite connection and any
    cached metadata are rebuilt lazily in the consumer process.
    """

    def __init__(self, profile_db_path, departure_from_reference_threshold,
                 columns_of_interest=None):
        filesnpaths.is_file_exists(profile_db_path)
        utils.is_profile_db(profile_db_path)

        self.profile_db_path = profile_db_path
        self.threshold = float(departure_from_reference_threshold)
        self.columns_of_interest = (list(columns_of_interest)
                                    if columns_of_interest is not None
                                    else list(DEFAULT_COLUMNS_OF_INTEREST))

        self._reset_caches()


    def _reset_caches(self):
        self._db = None
        self._sample_ids = None
        self._splits_with_snvs = None
        self._contig_to_splits = None


    def __getstate__(self):
        # Don't carry the connection or derived caches across a pickle boundary.
        # The child process rebuilds them on first use.
        return {
            'profile_db_path': self.profile_db_path,
            'threshold': self.threshold,
            'columns_of_interest': self.columns_of_interest,
        }


    def __setstate__(self, state):
        self.__dict__.update(state)
        self._reset_caches()


    def __enter__(self):
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disconnect()
        return False


    def _get_db(self):
        if self._db is None:
            self._db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        return self._db


    def disconnect(self):
        if self._db is not None:
            self._db.disconnect()
            self._db = None


    @property
    def _threshold_clause(self):
        return f"departure_from_reference >= {self.threshold}"


    def is_indexed(self):
        """Whether the `variable_nucleotides(split_name)` index exists in this profile db.

        The index is created on demand via the `anvi-index-table` program. Without it, the
        per-split and per-region queries below still work but fall back to full table scans.
        """
        index_name = t.index_name_for(t.variable_nts_table_name, ['split_name'])
        db_ = self._get_db()
        response = db_.db._exec("SELECT name FROM sqlite_master WHERE type='index' AND name=?", (index_name,))
        return len(response.fetchall()) > 0


    def warn_or_raise_on_missing_index(self, run=None, size_threshold_gb=10):
        """Decide, based on profile-db size, whether a missing index is fatal or merely slow.

        If the `variable_nucleotides(split_name)` index is present, this is a no-op. If it is
        missing, the cost of the full table scans the queries fall back to scales with the
        table size, so:

          - On a profile db larger than `size_threshold_gb`, raise a `ConfigError` telling the
            user exactly how to build the index with `anvi-index-table` and re-run. Scanning a
            table this large on every query would take tens of minutes.
          - On a smaller profile db, emit a warning (if `run` is given) and return, letting the
            caller proceed: at this scale the full scans are cheap enough.

        The decision stays with the caller: it chooses whether to call this and with what
        threshold. `anvi-report-dgrs` is the intended first consumer.
        """
        if self.is_indexed():
            return

        db_size = os.path.getsize(self.profile_db_path)
        human_size = utils.human_readable_file_size(db_size)
        command = (f"anvi-index-table {self.profile_db_path} "
                   f"--table {t.variable_nts_table_name} --column split_name")

        if db_size > size_threshold_gb * 1e9:
            raise ConfigError(
                f"The profile database at '{self.profile_db_path}' is large ({human_size}) and the "
                f"'{t.variable_nts_table_name}(split_name)' index it needs for fast per-split SNV access does "
                f"not exist. Without it, each query would do a full table scan that can take tens of minutes. "
                f"Please build the index once, and then re-run this program:\n\n"
                f"    {command}\n\n"
                f"You can reclaim the disk space later by running the same command with the `--drop-index` flag.")

        if run:
            run.warning(
                f"The '{t.variable_nts_table_name}(split_name)' index does not exist in this profile database, so "
                f"SNV queries will fall back to full table scans. This profile db is small enough ({human_size}) "
                f"that anvi'o will just go ahead. If you find this step slow, you can build the index with: "
                f"{command}",
                header="NO SNV INDEX -- PROCEEDING ANYWAY", lc='yellow')


    def get_sample_ids(self):
        """Return the sorted list of sample IDs that have at least one SNV passing the threshold."""
        if self._sample_ids is None:
            db_ = self._get_db()
            rows = db_.db.get_single_column_from_table(
                t.variable_nts_table_name, 'sample_id',
                unique=True, where_clause=self._threshold_clause)
            self._sample_ids = sorted(rows)
        return list(self._sample_ids)


    def get_splits_with_snvs(self):
        """Return the set of split names that have at least one SNV passing the threshold."""
        self._ensure_split_index()
        return set(self._splits_with_snvs)


    def get_contig_to_splits_map(self):
        """Return `{contig_name: set(split_names)}` restricted to splits that carry SNVs.

        Built from a single `SELECT DISTINCT split_name` so the cost is one
        query regardless of table size. Cached after first call.
        """
        self._ensure_split_index()
        return self._contig_to_splits


    def _ensure_split_index(self):
        if self._splits_with_snvs is not None:
            return

        db_ = self._get_db()
        splits = db_.db.get_single_column_from_table(
            t.variable_nts_table_name, 'split_name',
            unique=True, where_clause=self._threshold_clause)

        self._splits_with_snvs = set(splits)
        contig_to_splits = {}
        for split_name in self._splits_with_snvs:
            contig_name = split_name.split('_split_', 1)[0]
            contig_to_splits.setdefault(contig_name, set()).add(split_name)
        self._contig_to_splits = contig_to_splits


    def get_snvs_for_splits(self, split_names):
        """Return a DataFrame of SNVs for any of the given split names.

        Used per-bin (collections mode) and per-split (standard mode). Batches
        the `IN (...)` clause to keep statements bounded; concats the per-batch
        frames at the end.
        """
        split_names = list(split_names)
        if not split_names:
            return self._empty_frame()

        chunks = []
        for i in range(0, len(split_names), SPLIT_IN_LIST_BATCH):
            chunks.append(self._fetch_for_split_in(split_names[i:i + SPLIT_IN_LIST_BATCH]))

        df = chunks[0] if len(chunks) == 1 else pd.concat(chunks, ignore_index=True)
        return self._categorize(df)


    def get_snvs_for_region(self, contig_name, start, end, sample_id=None):
        """Return SNVs in [start, end] on `contig_name`, optionally restricted to one sample.

        Looks up the splits for `contig_name` in the cached contig->splits
        map (so the WHERE clause uses the indexed `split_name` column) and
        adds a `pos_in_contig BETWEEN ? AND ?` range. Returns an empty frame
        if the contig has no SNV-carrying splits.
        """
        self._ensure_split_index()
        splits = self._contig_to_splits.get(contig_name, set())
        if not splits:
            return self._empty_frame()

        # Build the IN-list manually; the codebase elsewhere (smart_get,
        # get_table_as_dataframe) uses the same pattern of inlining quoted
        # string literals rather than parameterized placeholders.
        in_list = ','.join(f'"{s}"' for s in splits)
        where = (f"{self._threshold_clause} "
                 f"AND split_name IN ({in_list}) "
                 f"AND pos_in_contig >= {int(start)} "
                 f"AND pos_in_contig <= {int(end)}")
        if sample_id is not None:
            # Escape any embedded single-quote -- belt-and-suspenders; anvi'o
            # sample IDs are constrained elsewhere to alnum+underscore.
            sample_lit = str(sample_id).replace("'", "''")
            where += f" AND sample_id = '{sample_lit}'"

        df = self._fetch_with_where(where)
        return self._categorize(df)


    def _fetch_for_split_in(self, split_names):
        in_list = ','.join(f'"{s}"' for s in split_names)
        where = f"{self._threshold_clause} AND split_name IN ({in_list})"
        return self._fetch_with_where(where)


    def _fetch_with_where(self, where_clause):
        db_ = self._get_db()
        col_list = ', '.join(f'"{c}"' for c in self.columns_of_interest)
        # Use pd.read_sql_query rather than db.get_table_as_dataframe because
        # the latter raises on empty results by default; for this accessor an
        # empty result is a normal outcome (e.g. a region with no SNVs).
        query = (f'SELECT {col_list} FROM "{t.variable_nts_table_name}" '
                 f'WHERE {where_clause}')
        df = pd.read_sql_query(query, db_.db.conn)
        return self._downcast(df)


    def _empty_frame(self):
        return pd.DataFrame(columns=list(self.columns_of_interest) + ['contig_name'])


    def _downcast(self, df):
        """Match the dtype contract produced by the historical `init_snv_table` in dgrs.py."""
        if df.empty:
            return df
        if 'pos_in_contig' in df.columns:
            df['pos_in_contig'] = df['pos_in_contig'].astype(np.int32)
        if 'base_pos_in_codon' in df.columns:
            df['base_pos_in_codon'] = df['base_pos_in_codon'].astype(np.int8)
        if 'departure_from_reference' in df.columns:
            df['departure_from_reference'] = df['departure_from_reference'].astype(np.float32)
        for col in nucleotides:
            if col in df.columns:
                df[col] = df[col].astype(np.int32)
        return df


    def _categorize(self, df):
        """Add `contig_name` and convert string-valued columns to `category` dtype.

        Done after concatenation rather than per-batch so the resulting
        category dictionary covers the full result set instead of being merged
        from per-batch dictionaries (which is what pandas does, but at a cost).
        """
        if df.empty:
            return df
        df['contig_name'] = df['split_name'].astype(str).str.split('_split_', n=1).str[0]
        for col in ('sample_id', 'split_name', 'contig_name', 'reference'):
            if col in df.columns:
                df[col] = df[col].astype('category')
        return df
