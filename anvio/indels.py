"""Cross-sample analysis of large indels (often mobile genetic elements) from
profile-db indels tables — produces a loci × samples allele frequency matrix
for visualization with `anvi-interactive --manual`.

The pipeline:

  1. Load indel rows (INS + DEL) from one or more profile-dbs.
  2. Reconstruct DEL sequences from the contigs-db (only INS sequences are stored
     in the indels table; DEL sequences are contig-specific and not in reads).
  3. Cluster indel sequences by sequence similarity. Identity threshold is
     length-tiered (longer indels tolerate more divergence) — each tier runs
     mmseqs2 easy-cluster independently.
  4. Within each cluster, collapse events that occur at nearby contig positions
     into a single 'locus' (a unique insertion/deletion site for this family).
  5. Compute the allele frequency (count / coverage) of each locus in each sample.
  6. Emit a samples × loci matrix + row decorations + items-order Newick suitable
     for `anvi-interactive --manual`.
"""

import os
import mmap
import array
import bisect
import shutil
import statistics
import numpy as np
from collections import Counter, defaultdict

import anvio
import anvio.tables as t
import anvio.dbinfo as dbi
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.drivers.mmseqs2 import MMseqs2


# For a sample with no recorded DEL event at a locus, we use the auxiliary per-nt coverage
# to ask "is the reference allele (= the element) actually present here?". A sample is
# considered to carry the element if at least this fraction of the DELd region has ≥1x coverage.
DETECTION_THRESHOLD_FOR_PRESENCE = 0.5


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Florian Trigodet"
__email__ = "trigodet.florian@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
PL = terminal.pluralize


# Length-tiered identity defaults for clustering. Shorter indels need stricter
# identity (random matches over 500 bp are common at low identity); longer
# elements can tolerate more divergence because their length itself anchors
# the homology call.
DEFAULT_LENGTH_TIERS = [
    # (label,        max_length_exclusive, min_seq_id)
    ('short',        1000,                 0.95),
    ('mid',          5000,                 0.90),
    ('long',         float('inf'),         0.85),
]


# Preferred functional annotation sources, in priority order. Whichever one is
# present first wins for the per-locus 'gene_function' decoration.
PREFERRED_FUNCTION_SOURCES = ['KOfam', 'COG20_FUNCTION', 'Pfam', 'COG14_FUNCTION', 'KEGG_Module']


class _IndelRow:
    """Read-only, dict-like view of one row of an `_IndelStore`, so that `e['count']` and `e.get(...)` work."""

    __slots__ = ('store', 'i')
    KEYS = ('sample_id', 'profile_db', 'split_name', 'contig_name', 'pos_in_contig', 'corresponding_gene_call',
            'reference', 'type', 'sequence', 'length', 'count', 'coverage')

    def __init__(self, store, i):
        self.store, self.i = store, i


    def __getitem__(self, key):
        return self.store.value(self.i, key)


    def get(self, key, default=None):
        return self.store.value(self.i, key) if key in self.KEYS else default


    def keys(self):
        return self.KEYS


class _IndelStore:
    """Columnar storage for indel rows.

    Integer columns are numpy arrays, strings are codes into small tables, and the indel sequences live in a
    file on disk that is memory-mapped for reading. Hundreds of millions of rows then fit in tens of GB, where
    one Python dict per row (with its sequence) does not fit at all. Rows are appended while the profile-dbs
    are read, `finalize` turns the columns into numpy arrays, and `take` keeps or reorders rows.
    """

    TYPES = ('DEL', 'INS')   # code order is string order, so type codes sort like the type strings
    INT_COLS = {'sample': 'i', 'profile': 'i', 'split': 'i', 'contig': 'i', 'pos': 'q', 'gene': 'q', 'ref': 'i',
                'type': 'b', 'length': 'i', 'count': 'q', 'coverage': 'q', 'seq_off': 'q', 'seq_len': 'i'}

    def __init__(self, work_dir):
        self.seq_path = os.path.join(work_dir, 'indel_sequences.bin')
        self._fh = open(self.seq_path, 'wb')
        self._end = 0
        self._mm = None
        self._mm_size = 0
        self._cols = {k: array.array(code) for k, code in self.INT_COLS.items()}
        self.tables = {k: [] for k in ('sample', 'profile', 'split', 'contig', 'ref')}
        self._index = {k: {} for k in self.tables}
        self.n = 0


    def code(self, table, value):
        """The integer code of `value` in string table `table` (added if new)."""

        index = self._index[table]
        c = index.get(value)
        if c is None:
            c = index[value] = len(self.tables[table])
            self.tables[table].append(value)

        return c


    def _write_seq(self, seq):
        b = seq.encode() if isinstance(seq, str) else seq
        offset = self._end
        self._fh.write(b)
        self._end += len(b)

        return offset, len(b)


    def append(self, sample, profile, split, contig, pos, gene, ref, typ, seq, length, count, coverage):
        c = self._cols
        offset, seq_length = self._write_seq(seq)
        c['sample'].append(self.code('sample', sample))
        c['profile'].append(self.code('profile', profile))
        c['split'].append(self.code('split', split))
        c['contig'].append(self.code('contig', contig))
        c['pos'].append(pos)
        c['gene'].append(gene)
        c['ref'].append(self.code('ref', ref))
        c['type'].append(self.TYPES.index(typ))
        c['length'].append(length)
        c['count'].append(count)
        c['coverage'].append(coverage)
        c['seq_off'].append(offset)
        c['seq_len'].append(seq_length)
        self.n += 1


    def finalize(self):
        """Turn the append buffers into numpy arrays (call once, after the last `append`)."""

        dtypes = {'i': np.int32, 'q': np.int64, 'b': np.int8}
        for k, code in self.INT_COLS.items():
            setattr(self, k, np.frombuffer(self._cols[k], dtype=dtypes[code]).copy())
        self._cols = None
        self._fh.flush()


    def _map(self):
        if self._mm is None or self._mm_size != self._end:
            self._fh.flush()
            if self._mm is not None:
                self._mm.close()
                self._mm = None
            if self._end == 0:
                return b''
            with open(self.seq_path, 'rb') as f:
                self._mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            self._mm_size = self._end

        return self._mm


    def seq_bytes(self, i):
        seq_length = int(self.seq_len[i])
        if seq_length == 0:
            return b''
        offset = int(self.seq_off[i])

        return self._map()[offset:offset + seq_length]


    def set_seq(self, i, seq):
        self.seq_off[i], self.seq_len[i] = self._write_seq(seq)


    def take(self, order):
        """Keep and reorder rows: every column becomes column[order]."""

        for k in self.INT_COLS:
            setattr(self, k, getattr(self, k)[order])
        self.n = len(order)


    def value(self, i, key):
        if key == 'sample_id':
            return self.tables['sample'][self.sample[i]]
        if key == 'profile_db':
            return self.tables['profile'][self.profile[i]]
        if key == 'split_name':
            return self.tables['split'][self.split[i]]
        if key == 'contig_name':
            return self.tables['contig'][self.contig[i]]
        if key == 'pos_in_contig':
            return int(self.pos[i])
        if key == 'corresponding_gene_call':
            return int(self.gene[i])
        if key == 'reference':
            return self.tables['ref'][self.ref[i]]
        if key == 'type':
            return self.TYPES[self.type[i]]
        if key == 'sequence':
            return self.seq_bytes(i).decode()
        if key == 'length':
            return int(self.length[i])
        if key == 'count':
            return int(self.count[i])
        if key == 'coverage':
            return int(self.coverage[i])

        raise KeyError(key)


    def ranks(self, table):
        """rank[code] is the position of that code's string in sorted (Python string) order."""

        values = self.tables[table]
        order = sorted(range(len(values)), key=lambda j: values[j])
        r = np.empty(len(values), dtype=np.int64)
        r[order] = np.arange(len(values))

        return r


    def close(self):
        if self._mm is not None:
            self._mm.close()
            self._mm = None
        if not self._fh.closed:
            self._fh.close()


    def __len__(self):
        return self.n


    def __bool__(self):
        return self.n > 0


    def __getitem__(self, i):
        return _IndelRow(self, i)


    def __iter__(self):
        for i in range(self.n):
            yield _IndelRow(self, i)


class Indels:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x, default=None: args.__dict__[x] if x in args.__dict__ and args.__dict__[x] is not None else default

        # input
        self.contigs_db_path = A('contigs_db')
        self.profile_db_paths = A('profile_db', []) or []
        self.profile_dbs_file = A('profile_dbs_file')

        # event filters
        self.min_length = A('min_length', 500)
        self.max_length = A('max_length')
        self.event_type = A('event_type', 'both')

        # clustering
        self.min_identity_override = A('min_identity')
        self.min_coverage = A('min_coverage', 0.8)
        self.position_tolerance = A('position_tolerance', 50)
        self.min_samples = A('min_samples', 2)
        self.min_reads_per_element = A('min_reads_per_element', 5)

        # scale
        self.exclude_contig_edges = A('exclude_contig_edges', 0) or 0
        self.max_loci_for_tree = A('max_loci_for_tree', 0) or 0

        # output
        self.output_dir = A('output_dir')
        self.num_threads = A('num_threads', 1)
        self.just_do_it = A('just_do_it', False)

        # state populated by process()
        self.contigs_super = None
        self.samples = []                 # ordered list of sample names
        self.sample_to_profile_path = {}  # sample_id → profile-db path it lives in
        self.aux_data_caches = {}         # profile-db path → AuxiliaryDataForSplitCoverages (or None)
        self._split_index = {}            # contig → (split starts, ends, names), for DEL-locus aux lookups
        self.indels = []                  # list of dicts (one per accepted indel row)
        self.clusters_to_members = {}     # {cluster_id: [indel idx, ...]}
        self.cluster_representatives = {} # {cluster_id: fasta_id_of_representative}
        self.loci = {}                    # {locus_id: {...}}
        self.locus_event_types = {}       # locus_id → 'INS' or 'DEL' (dominant type)
        self.matrix = {}                  # {locus_id: {sample: presence}}
        self.items_additional = {}        # {locus_id: {decoration_col: value, ...}}
        self.items_tree_newick = None
        self.n_edge_rows = 0              # INS rows set aside by --exclude-contig-edges
        self.work_dir = None              # holds the indel sequences (and other large temporary files)


    def sanity_check(self):
        if not self.contigs_db_path:
            raise ConfigError("You need to provide a contigs database (-c / --contigs-db).")
        filesnpaths.is_file_exists(self.contigs_db_path)

        # gather profile-db paths from either -p (repeated) or --profile-dbs-file
        if self.profile_dbs_file and self.profile_db_paths:
            raise ConfigError("Pick one: either pass -p/--profile-db one or more times, OR pass "
                              "--profile-dbs-file with a file listing profile-db paths. Not both.")

        if self.profile_dbs_file:
            filesnpaths.is_file_exists(self.profile_dbs_file)
            with open(self.profile_dbs_file) as f:
                self.profile_db_paths = [line.strip() for line in f if line.strip() and not line.startswith('#')]

        if not self.profile_db_paths:
            raise ConfigError("You need at least one profile database. Use -p / --profile-db (you can pass "
                              "this flag multiple times) or --profile-dbs-file.")

        for p in self.profile_db_paths:
            filesnpaths.is_file_exists(p)

        if self.event_type not in ('INS', 'DEL', 'both'):
            raise ConfigError(f"--event-type must be one of INS, DEL, both (got '{self.event_type}').")

        if self.min_length < 1:
            raise ConfigError("--min-length must be >= 1.")

        if self.max_length is not None and self.max_length < self.min_length:
            raise ConfigError("--max-length must be >= --min-length.")

        if self.min_identity_override is not None:
            if not (0 < self.min_identity_override <= 1.0):
                raise ConfigError("--min-identity must be in (0, 1].")

        if not (0 < self.min_coverage <= 1.0):
            raise ConfigError("--min-coverage must be in (0, 1].")

        if self.position_tolerance < 0:
            raise ConfigError("--position-tolerance must be >= 0.")

        if self.min_samples < 1:
            raise ConfigError("--min-samples must be >= 1.")

        if self.min_reads_per_element < 1:
            raise ConfigError("--min-reads-per-element must be >= 1.")

        if not self.output_dir:
            raise ConfigError("You need to provide an output directory with -o / --output-dir.")

        # verify all profile-dbs match the contigs-db hash
        contigs_hash = dbi.DBInfo(self.contigs_db_path).hash
        for p in self.profile_db_paths:
            profile_hash = dbi.DBInfo(p).hash
            if profile_hash != contigs_hash:
                msg = (f"Profile-db '{p}' has contigs_db_hash '{profile_hash}', but the supplied "
                       f"contigs-db has hash '{contigs_hash}'. They were not produced together.")
                if self.just_do_it:
                    self.run.warning(msg + " You said --just-do-it, so we will press on, but expect chaos.")
                else:
                    raise ConfigError(msg + " If you really know what you are doing, you can re-run with --just-do-it.")

        # external program
        utils.is_program_exists('mmseqs')

        # ensure output dir is writable
        filesnpaths.gen_output_directory(self.output_dir, delete_if_exists=False)


    def process(self):
        """Top-level orchestration."""
        self.sanity_check()

        self.run.warning(None, header="anvi-report-indels", lc="green")
        self.run.info('Contigs database', self.contigs_db_path)
        self.run.info('Profile databases', f"{len(self.profile_db_paths)} ({', '.join(self.profile_db_paths)})")
        self.run.info('Min indel length', self.min_length)
        self.run.info('Max indel length', self.max_length if self.max_length else 'no upper limit')
        self.run.info('Event type', self.event_type)
        self.run.info('Position tolerance', f"{self.position_tolerance} bp")
        self.run.info('Min samples per locus', self.min_samples)
        self.run.info('Min reads per element', self.min_reads_per_element)

        # the indel rows keep their sequences on disk, in a work directory inside the output directory
        self.work_dir = os.path.join(self.output_dir, '.anvi-report-indels-tmp')
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)
        os.makedirs(self.work_dir)

        self._init_contigs_super()
        self.collect_indels()
        self._set_aside_contig_edge_ins()
        if not self.indels:
            raise ConfigError(f"No indels passed --min-length {self.min_length} (and other filters) "
                              f"across the supplied profile-db(s). Nothing to report.")

        self.reconstruct_del_sequences()
        self.cluster_sequences()
        self._filter_clusters_by_occurrence()
        self.define_loci()
        if not self.loci:
            raise ConfigError(f"No locus remained after collapsing events within "
                              f"--position-tolerance {self.position_tolerance} bp. Try raising "
                              f"the tolerance (or check that any indel events survived the "
                              f"upstream filters).")

        self.build_matrix()
        self._apply_min_samples_filter()
        if not self.loci:
            raise ConfigError(f"No locus has ≥ {self.min_samples} sample(s) with element "
                              f"presence > 0. Lower --min-samples or relax upstream filters "
                              f"(--min-reads-per-element, --min-length).")

        self.build_items_additional_data()
        self.build_items_tree()
        self.build_long_format()
        self.report()

        self.indels.close()
        shutil.rmtree(self.work_dir, ignore_errors=True)


    def _filter_clusters_by_occurrence(self):
        """Drop clusters whose total supporting-read count (summed `count` across all member
        indel rows, including INS and DEL events from any sample) is below
        --min-reads-per-element. A 'cluster' here represents one MGE/element family, so this
        filter removes elements that the data simply doesn't see often enough to trust."""

        if self.min_reads_per_element <= 1:
            return

        kept_clusters = {}
        kept_reps = {}
        n_clusters_dropped = 0
        n_events_dropped = 0
        n_reads_dropped = 0

        for cluster_id, member_idxs in self.clusters_to_members.items():
            total_reads = sum(self.indels[idx]['count'] for idx in member_idxs)
            if total_reads >= self.min_reads_per_element:
                kept_clusters[cluster_id] = member_idxs
                if cluster_id in self.cluster_representatives:
                    kept_reps[cluster_id] = self.cluster_representatives[cluster_id]
            else:
                n_clusters_dropped += 1
                n_events_dropped += len(member_idxs)
                n_reads_dropped += total_reads

        if n_clusters_dropped:
            self.run.info(
                f"Clusters dropped by --min-reads-per-element ({self.min_reads_per_element})",
                f"{pp(n_clusters_dropped)} cluster(s), {pp(n_events_dropped)} event(s), "
                f"{pp(n_reads_dropped)} supporting read(s)"
            )

        self.clusters_to_members = kept_clusters
        self.cluster_representatives = kept_reps

        if not self.clusters_to_members:
            raise ConfigError(
                f"All clusters were filtered out by --min-reads-per-element "
                f"{self.min_reads_per_element}. Try lowering that threshold.")


    def _init_contigs_super(self):
        """Instantiate ContigsSuperclass once and load gene calls and functions (contig sequences are
        loaded later, and only for the contigs with DEL events; see `reconstruct_del_sequences`)."""
        import argparse
        self.progress.new('Initializing contigs database')
        self.progress.update('...')
        args = argparse.Namespace(contigs_db=self.contigs_db_path)
        self.contigs_super = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        try:
            self.contigs_super.init_functions()
        except ConfigError:
            # functional annotations may not be present; that's OK.
            pass
        self._build_gene_position_index()
        self.progress.end()


    def _build_gene_position_index(self):
        """Build a per-contig sorted gene index so we can answer 'which gene overlaps position P
        on contig C?' in O(log N) instead of O(total_genes_in_db).

        ContigsSuperclass.nt_position_to_gene_caller_id() walks the entire genes_in_contigs_dict
        on every call, which is fine for a handful of lookups but blows up on a metagenome with
        hundreds of thousands of gene calls × thousands of indel loci (=> hundreds of millions of
        comparisons). This index turns each lookup into a binary search over the genes on a
        single contig, then a tiny walk over any (rare) overlapping genes.
        """
        self._gene_index = None
        if not self.contigs_super.a_meta.get('genes_are_called', False):
            return

        genes_by_contig = defaultdict(list)
        for gene_id, gene_call in self.contigs_super.genes_in_contigs_dict.items():
            genes_by_contig[gene_call['contig']].append((gene_call['start'], gene_call['stop'], gene_id))

        # Per-contig: sort by start and pre-extract the starts array for bisect lookups.
        self._gene_index = {}
        for contig, genes in genes_by_contig.items():
            genes.sort()
            self._gene_index[contig] = {
                'starts': [g[0] for g in genes],
                'genes': genes,
            }


    def _gene_caller_id_at(self, contig_name, position_in_contig):
        """Return the gene_caller_id of the gene overlapping `position_in_contig` on `contig_name`,
        or -1 if none. Uses the pre-built index from `_build_gene_position_index`."""
        if self._gene_index is None:
            return -1
        entry = self._gene_index.get(contig_name)
        if not entry:
            return -1
        # rightmost gene-start that is <= position
        idx = bisect.bisect_right(entry['starts'], position_in_contig) - 1
        # walk back over any (rare) overlapping genes whose start <= position. Stop as soon
        # as a gene starts after `position` or we've left the array.
        while idx >= 0:
            start, stop, gene_id = entry['genes'][idx]
            if start > position_in_contig:
                break
            if position_in_contig < stop:
                return gene_id
            idx -= 1
        return -1


    def collect_indels(self):
        """Load indels rows from each profile-db, apply filters, and keep them in a columnar store (sequences
        on disk). Rows are then put in a canonical, content-only order."""

        self.progress.new('Collecting indels')

        # explicit SELECT — avoids the leading entry_id that the generic
        # get_some_rows_from_table() helper prepends for this table.
        select_cols = ['sample_id', 'split_name', 'pos_in_contig', 'corresponding_gene_call',
                       'reference', 'type', 'sequence', 'length', 'count', 'coverage']
        select_clause = ', '.join(select_cols)

        where = f"length >= {int(self.min_length)}"
        if self.max_length is not None:
            where += f" AND length <= {int(self.max_length)}"
        if self.event_type in ('INS', 'DEL'):
            where += f" AND type = '{self.event_type}'"

        splits_to_contigs = {s: info['parent'] for s, info in self.contigs_super.splits_basic_info.items()}

        store = _IndelStore(self.work_dir)
        sample_order = []
        sample_seen = set()

        for path in self.profile_db_paths:
            self.progress.update(f"Reading indels from {os.path.basename(path)} ...")

            pdb = dbops.ProfileDatabase(path, quiet=True)
            if not pdb.meta.get('INDELs_profiled', 0):
                self.progress.reset()
                self.run.warning(f"Profile-db '{path}' was not profiled with INDELs enabled. Skipping.")
                pdb.db.disconnect()
                continue

            # Record every sample known to this profile-db (even those with no indel rows
            # passing our filter) so they appear as columns in the matrix.
            for s in sorted(pdb.samples):
                if s not in sample_seen:
                    sample_seen.add(s)
                    sample_order.append(s)
                self.sample_to_profile_path[s] = path

            db_basename = os.path.basename(path).replace('.db', '')

            # rows are read in batches, so that a profile-db's indels table is never held in memory at once
            try:
                cursor = pdb.db._exec(f"SELECT {select_clause} FROM {t.indels_table_name} WHERE {where}")
                while True:
                    batch = cursor.fetchmany(100000)
                    if not batch:
                        break

                    for row in batch:
                        sample_id, split_name, pos_in_contig, gene_id, reference, ev_type, sequence, length, count, coverage = row
                        contig = splits_to_contigs.get(split_name)
                        if contig is None:
                            # split is unknown to this contigs-db; skip
                            continue

                        store.append(sample_id, db_basename, split_name, contig, int(pos_in_contig),
                                     int(gene_id) if gene_id is not None else -1, reference or '', ev_type,
                                     sequence or '', int(length), int(count), int(coverage))

                        if sample_id not in sample_seen:
                            sample_seen.add(sample_id)
                            sample_order.append(sample_id)
            finally:
                pdb.db.disconnect()

        store.finalize()

        # Put rows in a canonical, content-only order. Row order in a profile-db's indels table is
        # whatever order the profiler's worker processes happened to finish contigs in, so it differs
        # between two runs on the same BAM. Every later step (row indices, mmseqs2 input order and
        # representatives, cluster and locus numbering, the running-mean locus merge) depends on this
        # order, so without sorting the report is not reproducible.
        # The key is (contig_name, pos_in_contig, type, length, sequence, sample_id, split_name, count, coverage,
        # corresponding_gene_call, reference, profile_db), with strings compared as strings: the string columns are
        # sorted through the ranks of their codes, and sequences are ranked only among the rows that tie on
        # (contig_name, pos_in_contig, type, length), where they are needed.
        self.progress.update(f"Putting {pp(len(store))} rows in canonical order ...")
        if len(store):
            contig_rank = store.ranks('contig')[store.contig]
            seq_rank = np.zeros(len(store), dtype=np.int64)
            pre = np.lexsort((store.length, store.type, store.pos, contig_rank))
            keys = np.stack([contig_rank[pre], store.pos[pre], store.type[pre], store.length[pre]])
            new_group = np.ones(len(pre), dtype=bool)
            new_group[1:] = np.any(keys[:, 1:] != keys[:, :-1], axis=0)
            starts = np.nonzero(new_group)[0]
            ends = np.append(starts[1:], len(pre))
            for start, end in zip(starts[(ends - starts) > 1], ends[(ends - starts) > 1]):
                members = pre[start:end]
                seqs = [store.seq_bytes(i) for i in members]
                rank_of = {b: r for r, b in enumerate(sorted(set(seqs)))}
                for i, b in zip(members, seqs):
                    seq_rank[i] = rank_of[b]

            order = np.lexsort((store.ranks('profile')[store.profile], store.ranks('ref')[store.ref], store.gene,
                                store.coverage, store.count, store.ranks('split')[store.split],
                                store.ranks('sample')[store.sample], seq_rank, store.length, store.type,
                                store.pos, contig_rank))
            store.take(order)

        self.indels = store
        self.samples = sample_order
        self.progress.end()

        self.run.info('Indel rows passing filters', pp(len(self.indels)))
        self.run.info('Distinct samples', pp(len(self.samples)))


    def _set_aside_contig_edge_ins(self):
        """Move INS rows within --exclude-contig-edges bp of a contig end out of the analysis.

    Notes
    =====
    - The rows set aside are written to INDELS-CONTIG-EDGE-INS.txt in the work directory, and `report` moves
      that file into the output directory.
    """

        if self.exclude_contig_edges <= 0:
            return

        n = self.exclude_contig_edges
        store = self.indels
        lengths_by_name = {contig: info['length'] for contig, info in self.contigs_super.contigs_basic_info.items()}
        contig_lengths = np.array([lengths_by_name.get(c, -1) for c in store.tables['contig']], dtype=np.int64)[store.contig]

        is_edge = ((store.type == store.TYPES.index('INS')) & (contig_lengths >= 0) &
                   ((store.pos < n) | (store.pos >= contig_lengths - n)))

        self.n_edge_rows = int(is_edge.sum())
        if self.n_edge_rows:
            edge_cols = ['sample_id', 'profile_db', 'contig_name', 'split_name', 'pos_in_contig', 'type', 'length', 'count', 'coverage']
            with open(os.path.join(self.work_dir, 'INDELS-CONTIG-EDGE-INS.txt'), 'w') as f:
                f.write('\t'.join(edge_cols) + '\n')
                for i in np.nonzero(is_edge)[0]:
                    f.write('\t'.join(str(store.value(i, c)) for c in edge_cols) + '\n')

        store.take(np.nonzero(~is_edge)[0])
        self.run.info(f'INS rows within {n} bp of a contig end (set aside)', pp(self.n_edge_rows))


    def reconstruct_del_sequences(self):
        """Fill in `sequence` for DEL rows by slicing from the contig sequence."""

        # contig sequences are needed only here, and only for the contigs that carry DEL rows: loading every
        # contig instead takes about as much memory as the contigs-db is large (~150 GB for 141 Gbp)
        store = self.indels
        del_idxs = np.nonzero(store.type == store.TYPES.index('DEL'))[0]
        del_contigs = set(store.tables['contig'][c] for c in np.unique(store.contig[del_idxs]))
        if del_contigs:
            self.progress.new('Loading contig sequences')
            self.progress.update(f"for {pp(len(del_contigs))} contigs with DEL events ...")
            self.contigs_super.init_contig_sequences(contig_names_of_interest=del_contigs)
            self.progress.end()

        n_reconstructed = 0
        n_skipped_mismatch = 0
        n_skipped_oob = 0
        keep = np.ones(len(store), dtype=bool)

        for i in del_idxs:
            contig = store.value(i, 'contig_name')
            pos = store.value(i, 'pos_in_contig')
            length = store.value(i, 'length')
            reference = store.value(i, 'reference')

            contig_seq = self.contigs_super.contig_sequences.get(contig, {}).get('sequence')
            if not contig_seq:
                n_skipped_oob += 1
                keep[i] = False
                continue

            if pos < 0 or pos + length > len(contig_seq):
                n_skipped_oob += 1
                keep[i] = False
                continue

            del_seq = contig_seq[pos:pos + length].upper()

            if reference and del_seq[:1].upper() != reference.upper():
                # reference mismatch — likely a position-semantics edge case; warn and skip
                n_skipped_mismatch += 1
                keep[i] = False
                continue

            store.set_seq(i, del_seq)
            n_reconstructed += 1

        store.take(np.nonzero(keep)[0])

        if n_reconstructed:
            self.run.info('Reconstructed DEL sequences', pp(n_reconstructed))
        if n_skipped_mismatch:
            self.run.warning(f"{pp(n_skipped_mismatch)} DEL rows were dropped because the reconstructed first "
                             f"base did not match the stored `reference` value.")
        if n_skipped_oob:
            self.run.warning(f"{pp(n_skipped_oob)} DEL rows were dropped because their position fell outside "
                             f"the contig's known bounds (this should be rare; check your contigs-db).")


    def cluster_sequences(self):
        """Cluster indel sequences with mmseqs2. Identity threshold is length-tiered unless --min-identity overrides.
        Tiers whose shortest sequence falls below the mmseqs2 nucleotide k-mer size (15) fall back to exact-match
        clustering, which is the only sensible thing to do at those lengths anyway."""

        # split indels into tiers by length (unless --min-identity is set)
        if self.min_identity_override is not None:
            tiers = [('all', float('inf'), self.min_identity_override)]
        else:
            tiers = DEFAULT_LENGTH_TIERS

        # bin indels by tier: the first tier whose max_length_exclusive is above the length, else the last one
        lengths = self.indels.length
        tier_of = np.full(len(self.indels), len(tiers) - 1, dtype=np.int64)
        assigned = np.zeros(len(self.indels), dtype=bool)
        for k, (label, max_excl, min_id) in enumerate(tiers):
            hit = (~assigned) & (lengths < max_excl)
            tier_of[hit] = k
            assigned |= hit
        bins = {(label, min_id): [int(i) for i in np.nonzero(tier_of == k)[0]] for k, (label, _, min_id) in enumerate(tiers)}

        self.run.warning(None, header="CLUSTERING (length-tiered)", lc="green")

        all_clusters = {}
        all_reps = {}
        cluster_counter = 0
        work_dir = filesnpaths.get_temp_directory_path()

        MMSEQS_MIN_LEN = 15  # default nucleotide k-mer size in mmseqs2

        # iterate tiers in their defined order (not dict insertion order) so cluster numbering is canonical
        for label, _, min_id in tiers:
            idxs = bins.get((label, min_id), [])
            if not idxs:
                continue

            tier_lengths = [self.indels[i]['length'] for i in idxs]
            use_mmseqs = (min(tier_lengths) >= MMSEQS_MIN_LEN) and (len(idxs) > 1)

            method = 'mmseqs2' if use_mmseqs else 'exact-match'
            self.run.info(f"Tier '{label}'", f"{pp(len(idxs))} indels @ min_seq_id={min_id} ({method})")

            if not use_mmseqs:
                # exact-match clustering: group by (type, sequence) — useful when sequences are too short for k-mer search
                exact_groups = defaultdict(list)
                for idx in idxs:
                    e = self.indels[idx]
                    key = (e['type'], e['sequence'].upper())
                    exact_groups[key].append(idx)

                for key, members in exact_groups.items():
                    cluster_counter += 1
                    cluster_id = f"C_{cluster_counter:08d}"
                    all_clusters[cluster_id] = members
                    rep_idx = members[0]
                    e = self.indels[rep_idx]
                    all_reps[cluster_id] = f"{rep_idx}|{e['sample_id']}|{e['contig_name']}|{e['pos_in_contig']}|{e['type']}|{e['length']}"
                continue

            fasta_path = os.path.join(work_dir, f"indels_{label}.fa")
            with open(fasta_path, 'w') as fa:
                for idx in idxs:
                    e = self.indels[idx]
                    seq_id = f"{idx}|{e['sample_id']}|{e['contig_name']}|{e['pos_in_contig']}|{e['type']}|{e['length']}"
                    fa.write(f">{seq_id}\n{e['sequence']}\n")

            m = MMseqs2(fasta_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False),
                        num_threads=self.num_threads)
            m.min_seq_id = min_id
            m.coverage = self.min_coverage
            m.cov_mode = 0   # bidirectional coverage

            sub_clusters = m.get_clusters_dict(name_prefix=f"TIER{label.upper()}")

            # member ids look like '<idx>|sample|contig|pos|type|length' — pull the leading idx back out.
            # Number clusters by their smallest member index (canonical once rows are sorted), and keep
            # members sorted, instead of relying on the order mmseqs2 output was parsed in.
            parsed = []
            for sub_cid, members in sub_clusters.items():
                member_idxs = sorted(int(name.split('|', 1)[0]) for name in members)
                parsed.append((member_idxs, m.representatives.get(sub_cid)))
            for member_idxs, rep in sorted(parsed, key=lambda x: x[0][0]):
                cluster_counter += 1
                cluster_id = f"C_{cluster_counter:08d}"
                all_clusters[cluster_id] = member_idxs
                all_reps[cluster_id] = rep

            m.cleanup()

        self.clusters_to_members = all_clusters
        self.cluster_representatives = all_reps

        self.run.info('Total clusters', pp(len(self.clusters_to_members)))


    def define_loci(self):
        """Within each cluster, collapse events by contig + position window into loci."""

        loci = {}
        locus_counter = 0

        for cluster_id in sorted(self.clusters_to_members.keys()):
            member_idxs = self.clusters_to_members[cluster_id]

            # group by contig
            by_contig = defaultdict(list)
            for idx in member_idxs:
                by_contig[self.indels[idx]['contig_name']].append(idx)

            for contig in sorted(by_contig):
                idxs = by_contig[contig]
                # sort by pos_in_contig (ties by canonical row index), sweep, merge while within
                # tolerance of running mean
                idxs.sort(key=lambda i: (self.indels[i]['pos_in_contig'], i))

                current_members = []
                current_positions = []
                current_sum = 0  # running sum of current_positions (ints, so the mean is exact)

                def flush():
                    nonlocal locus_counter, current_members, current_positions, current_sum
                    if not current_members:
                        return
                    locus_counter += 1
                    locus_id = f"locus_{locus_counter:08d}"
                    loci[locus_id] = {
                        'cluster_id': cluster_id,
                        'contig_name': contig,
                        'positions': list(current_positions),
                        'mean_pos': int(round(current_sum / len(current_positions))),
                        'min_pos': min(current_positions),
                        'max_pos': max(current_positions),
                        'members': list(current_members),
                    }
                    current_members = []
                    current_positions = []
                    current_sum = 0

                for idx in idxs:
                    pos = self.indels[idx]['pos_in_contig']
                    if current_positions:
                        running_mean = current_sum / len(current_positions)
                        if pos - running_mean > self.position_tolerance:
                            flush()
                    current_members.append(idx)
                    current_positions.append(pos)
                    current_sum += pos

                flush()

        self.loci = loci
        self.run.info('Candidate loci defined', pp(len(self.loci)))


    def _open_aux_data_caches(self):
        """Open one AuxiliaryDataForSplitCoverages per profile-db so we can query per-nt coverage
        for samples that have no indel event at a DEL locus (where we still need to know whether
        the reference allele -- i.e. the element -- is present in that sample)."""

        contigs_hash = dbi.DBInfo(self.contigs_db_path).hash
        for path in self.profile_db_paths:
            if path in self.aux_data_caches:
                continue
            aux_path = os.path.join(os.path.dirname(path), 'AUXILIARY-DATA.db')
            if not os.path.exists(aux_path):
                self.run.warning(f"No AUXILIARY-DATA.db next to '{path}'. Samples in that profile-db "
                                 f"without DEL events will be reported as element-absent rather than "
                                 f"checked for actual coverage in the DEL region.")
                self.aux_data_caches[path] = None
                continue
            try:
                self.aux_data_caches[path] = auxiliarydataops.AuxiliaryDataForSplitCoverages(
                    aux_path, contigs_hash, ignore_hash=self.just_do_it)
            except Exception as e:
                self.run.warning(f"Could not open auxiliary data at '{aux_path}': {e}. "
                                 f"Falling back to element-absent default for missing DEL cells.")
                self.aux_data_caches[path] = None


    def _close_aux_data_caches(self):
        for aux in self.aux_data_caches.values():
            if aux is not None:
                try:
                    aux.close()
                except Exception:
                    pass
        self.aux_data_caches = {}


    def _build_split_index(self):
        """Index the splits of every contig by their start positions.

    Returns
    =======
    index : dict
        {contig_name: (starts, ends, split_names)}, where each list is sorted by split start, so
        the split that holds a given position can be found by bisection (see `_split_at`).
    """

        by_contig = defaultdict(list)
        for split_name, info in self.contigs_super.splits_basic_info.items():
            by_contig[info['parent']].append((info['start'], info['end'], split_name))

        index = {}
        for contig, splits in by_contig.items():
            splits.sort()
            index[contig] = ([s[0] for s in splits], [s[1] for s in splits], [s[2] for s in splits])

        return index


    def _split_at(self, contig, pos):
        """Find the split of a contig that contains a given position.

    Parameters
    ==========
    contig : str
        Contig name.
    pos : int
        Position in the contig.

    Returns
    =======
    split : tuple
        (split_name, split_start), or (None, None) if no split of `contig` contains `pos`.
    """

        entry = self._split_index.get(contig)
        if entry is None:
            return None, None

        starts, ends, names = entry
        i = bisect.bisect_right(starts, pos) - 1
        if i >= 0 and pos < ends[i]:
            return names[i], starts[i]

        return None, None


    def _detection_in_slice(self, coverage_array, rel_start, length):
        """Compute the fraction of a region with at least 1x coverage.

    Parameters
    ==========
    coverage_array : numpy.ndarray
        Per-nucleotide coverage of a split in one sample.
    rel_start : int
        Start of the region relative to the split start.
    length : int
        Length of the region.

    Returns
    =======
    detection : float
        Fraction of the `length` positions with ≥1x coverage.

    Notes
    =====
    - Detection is computed against the requested length rather than the slice length, so
      positions past the split boundary count as zero-coverage (for 20 kb splits and sub-kb
      indels this is essentially never tripped).
    """

        if length <= 0:
            return 0.0

        cov_slice = coverage_array[rel_start:rel_start + length]
        return float(np.count_nonzero(cov_slice)) / float(length)


    def _resolve_del_detection_queries(self, queries):
        """Answer the auxiliary-data coverage lookups `build_matrix` queued for DEL loci.

    A queued cell becomes 1.0 if its sample covers ≥ DETECTION_THRESHOLD_FOR_PRESENCE of the
    DELd region, and otherwise keeps the 0.0 it already holds.

    Parameters
    ==========
    queries : dict
        {(profile_db_path, split_name): [(locus_id, sample, rel_start, length), ...]}

    Notes
    =====
    - Lookups are answered grouped by profile-db and split, so each coverage blob is read and
      decompressed once, and no blob is kept once its split is done.
    """

        n_blobs = len(queries)
        self.progress.new('Checking DEL-region coverage', progress_total_items=n_blobs)
        self.progress.update('...')
        update_every = max(1, n_blobs // 200)

        splits_by_profile = defaultdict(list)
        for profile_path, split_name in queries:
            splits_by_profile[profile_path].append(split_name)

        done = 0
        for profile_path in sorted(splits_by_profile):
            aux = self.aux_data_caches.get(profile_path)
            for split_name in sorted(splits_by_profile[profile_path]):
                if done % update_every == 0:
                    self.progress.increment(increment_to=done)
                    self.progress.update(f"split {pp(done)} of {pp(n_blobs)} ...")
                done += 1

                try:
                    split_cov = aux.get(split_name)
                except Exception:
                    split_cov = {}

                for locus_id, sample, rel_start, length in queries[(profile_path, split_name)]:
                    # lookup unavailable -> the cell keeps the conservative default of 0.0
                    if not split_cov or sample not in split_cov:
                        continue

                    if self._detection_in_slice(split_cov[sample], rel_start, length) >= DETECTION_THRESHOLD_FOR_PRESENCE:
                        self.matrix[locus_id][sample] = 1.0

        self.progress.end()


    def build_matrix(self):
        """For each (sample, locus) cell, report the fraction of reads carrying the MGE *element*
        at this locus -- i.e., element presence -- on a consistent 0..1 scale, regardless of
        whether the event is recorded as an INS (element in read, not in reference) or as a DEL
        (element in reference, not in read).

        Per cell:
          - Sample has events at locus:
              * INS-dominant locus: presence = max(count/cov) across the locus's positions.
              * DEL-dominant locus: presence = 1 - max(count/cov)  (i.e., fraction of reads still
                carrying the reference allele = the element).
          - Sample has no events at locus:
              * INS-dominant locus: presence = 0.0 (no signal -> element absent).
              * DEL-dominant locus: look up the sample's per-nt coverage over the DELd region in
                the auxiliary-data db. If ≥ DETECTION_THRESHOLD_FOR_PRESENCE of positions have
                ≥1x coverage, presence = 1.0 (element is there; reads simply don't deviate);
                otherwise presence = 0.0 (no coverage = no signal we can interpret).

        Note: anvi-profile records one indels row per *unique* inserted/deleted sequence at a
        position; sequencing-error variants of the same biological event end up as many rows
        with count=1 that share the same per-position `coverage` value. We therefore sum count
        across rows at the same (sample, position) but use the position's coverage once."""

        # Determine the dominant event type per locus
        for locus_id, info in self.loci.items():
            type_counts = Counter(self.indels[idx]['type'] for idx in info['members'])
            ins_n = type_counts.get('INS', 0)
            del_n = type_counts.get('DEL', 0)
            # tie-breaker: INS wins, matching the simpler interpretation
            self.locus_event_types[locus_id] = 'DEL' if del_n > ins_n else 'INS'

        # If any locus is DEL-dominant, we'll need the aux-data lookups
        any_del_locus = any(t == 'DEL' for t in self.locus_event_types.values())
        if any_del_locus:
            self._open_aux_data_caches()
            self._split_index = self._build_split_index()

        # aux-coverage lookups for DEL-locus cells without events, answered after the loop grouped by
        # (profile-db, split): {(profile_path, split_name): [(locus_id, sample, rel_start, length), ...]}
        del_detection_queries = defaultdict(list)

        n_loci = len(self.loci)
        self.progress.new('Building presence matrix', progress_total_items=n_loci)
        self.progress.update('...')
        update_every = max(1, n_loci // 200)

        for i, (locus_id, info) in enumerate(self.loci.items()):
            if i % update_every == 0:
                self.progress.increment(increment_to=i)
                self.progress.update(f"locus {pp(i)} of {pp(n_loci)} ...")

            locus_type = self.locus_event_types[locus_id]

            # sample -> position -> [supporting_reads_total, coverage_at_position]
            per_sample_per_pos = defaultdict(lambda: defaultdict(lambda: [0, 0]))
            for idx in info['members']:
                e = self.indels[idx]
                slot = per_sample_per_pos[e['sample_id']][e['pos_in_contig']]
                slot[0] += e['count']
                slot[1] = e['coverage']

            # representative (pos, length) used for aux-coverage detection lookups; the DEL
            # region varies slightly with sequencing noise, so we take median across member events
            if locus_type == 'DEL':
                del_lengths = [self.indels[idx]['length'] for idx in info['members'] if self.indels[idx]['type'] == 'DEL']
                rep_length = int(statistics.median(del_lengths)) if del_lengths else 0
                rep_pos = info['mean_pos']
                split_name, split_start = self._split_at(info['contig_name'], rep_pos)
            else:
                rep_length = 0
                rep_pos = info['mean_pos']

            row = {}
            for sample in self.samples:
                positions = per_sample_per_pos.get(sample)

                if positions:
                    freqs = [(count / cov) if cov > 0 else 0.0 for count, cov in positions.values()]
                    freq = max(freqs) if freqs else 0.0
                    if locus_type == 'DEL':
                        row[sample] = max(0.0, 1.0 - freq)
                    else:
                        row[sample] = freq
                else:
                    # no event in this sample at this locus
                    if locus_type == 'INS':
                        row[sample] = 0.0
                    else:
                        # DEL locus: did this sample actually have coverage in the DELd region? 0.0
                        # until the aux lookup (queued here, answered below) says so; it stays 0.0 if
                        # the lookup is unavailable (conservative default)
                        row[sample] = 0.0
                        profile_path = self.sample_to_profile_path.get(sample)
                        if split_name is not None and profile_path and self.aux_data_caches.get(profile_path) is not None:
                            del_detection_queries[(profile_path, split_name)].append(
                                (locus_id, sample, rep_pos - split_start, rep_length))

            self.matrix[locus_id] = row

        self.progress.end()

        if del_detection_queries:
            self._resolve_del_detection_queries(del_detection_queries)

        if any_del_locus:
            self._close_aux_data_caches()
            self._split_index = {}


    def _apply_min_samples_filter(self):
        """Drop loci where fewer than --min-samples samples have element presence > 0 in the
        matrix. This counts samples with non-zero presence in the FINAL matrix regardless of
        whether they contributed an indel event or had their presence inferred from auxiliary
        coverage detection, giving INS-dominant and DEL-dominant loci a symmetric treatment
        (`the element appears in at least N samples`)."""

        if self.min_samples <= 1:
            return

        n_dropped = 0
        for locus_id in list(self.loci.keys()):
            row = self.matrix.get(locus_id, {})
            n_with_element = sum(1 for v in row.values() if v > 0)
            if n_with_element < self.min_samples:
                del self.loci[locus_id]
                del self.matrix[locus_id]
                if locus_id in self.locus_event_types:
                    del self.locus_event_types[locus_id]
                n_dropped += 1

        if n_dropped:
            self.run.info(
                f"Loci dropped by --min-samples ({self.min_samples})",
                f"{pp(n_dropped)} (locus needs ≥{self.min_samples} sample(s) with element presence > 0)"
            )
        self.run.info('Final loci', pp(len(self.loci)))


    def _gene_function_for(self, gene_caller_id):
        """Return (source, accession, function_text) using the preferred source order."""
        if gene_caller_id is None or gene_caller_id < 0:
            return ('', '', '')

        calls = getattr(self.contigs_super, 'gene_function_calls_dict', {})
        gene_calls = calls.get(gene_caller_id)
        if not gene_calls:
            return ('', '', '')

        order = PREFERRED_FUNCTION_SOURCES + [s for s in gene_calls.keys() if s not in PREFERRED_FUNCTION_SOURCES]
        for source in order:
            entry = gene_calls.get(source)
            if entry:
                accession, function_text = entry[0], entry[1]
                # function_text may carry multiple !!!-separated hits; keep the first
                return (source, accession or '', (function_text or '').split('!!!')[0])
        return ('', '', '')


    def build_items_additional_data(self):
        """Build per-locus decoration columns for `--items-additional-data`."""

        has_genes_called = self.contigs_super.a_meta.get('genes_are_called', False)

        n_loci = len(self.loci)
        self.progress.new('Building per-locus decoration data', progress_total_items=n_loci)
        self.progress.update('...')

        # per-cluster values are computed once, over every member of the cluster, so that all loci of a
        # cluster report the same value (a per-locus median would give one cluster several different ones)
        loci_per_cluster = Counter(info['cluster_id'] for info in self.loci.values())
        cluster_median_length = {cid: int(statistics.median(self.indels[i]['length'] for i in self.clusters_to_members[cid]))
                                 for cid in loci_per_cluster}

        update_every = max(1, n_loci // 200)
        for i, (locus_id, info) in enumerate(self.loci.items()):
            if i % update_every == 0:
                self.progress.increment(increment_to=i)
                self.progress.update(f"locus {pp(i)} of {pp(n_loci)} ...")

            members = [self.indels[idx] for idx in info['members']]
            sample_counts = Counter(m['sample_id'] for m in members)
            lengths = [m['length'] for m in members]
            n_supporting_reads = sum(m['count'] for m in members)

            mean_pos = info['mean_pos']
            gene_caller_id = self._gene_caller_id_at(info['contig_name'], mean_pos) if has_genes_called else -1

            source, accession, function = self._gene_function_for(gene_caller_id)

            self.items_additional[locus_id] = {
                'cluster_id': info['cluster_id'],
                'locus_type': self.locus_event_types.get(locus_id, ''),
                'contig_name': info['contig_name'],
                'mean_pos_in_contig': mean_pos,
                'pos_range': f"{info['min_pos']}-{info['max_pos']}",
                'cluster_size_median_nt': cluster_median_length[info['cluster_id']],
                'n_loci_in_cluster': loci_per_cluster[info['cluster_id']],
                'locus_size_median_nt': int(statistics.median(lengths)),
                'n_events': len(members),
                'n_supporting_reads': n_supporting_reads,
                'n_samples_supporting': len(sample_counts),
                'in_orf': 'yes' if gene_caller_id >= 0 else 'no',
                'gene_caller_id': gene_caller_id if gene_caller_id >= 0 else '',
                'gene_function_source': source,
                'gene_function_accession': accession,
                'gene_function': function,
                'representative_id': self.cluster_representatives.get(info['cluster_id'], ''),
            }

        self.progress.end()


    def build_items_tree(self):
        """Hierarchical clustering of the matrix → Newick string."""
        if len(self.matrix) < 3:
            self.items_tree_newick = None
            return

        if self.max_loci_for_tree and len(self.matrix) > self.max_loci_for_tree:
            self.run.warning(f"Skipping the items-order tree: {pp(len(self.matrix))} loci is more than "
                             f"--max-loci-for-tree {pp(self.max_loci_for_tree)} (its hierarchical clustering "
                             f"needs memory that grows with the square of the number of loci).")
            self.items_tree_newick = None
            return

        self.progress.new('Computing items-order tree')
        self.progress.update(f"hierarchical clustering of {pp(len(self.matrix))} loci × {pp(len(self.samples))} samples ...")
        self.items_tree_newick = clustering.get_newick_tree_data_for_dict(
            self.matrix, distance='euclidean', linkage='ward', zero_fill_missing=True)
        self.progress.end()


    def build_long_format(self):
        """Map every indel row to its locus for the per-event audit trail, which `report` writes row by row
        straight from the indel store (one dict per row does not fit in memory for very large inputs)."""

        self._idx_to_locus = {}
        for locus_id, info in self.loci.items():
            for idx in info['members']:
                self._idx_to_locus[idx] = locus_id


    def _cluster_id_for_indel(self, idx):
        # build a reverse lookup once on first call
        if not hasattr(self, '_idx_to_cluster'):
            self._idx_to_cluster = {}
            for cid, members in self.clusters_to_members.items():
                for i in members:
                    self._idx_to_cluster[i] = cid
        return self._idx_to_cluster.get(idx, '')


    def report(self):
        """Write all output files into `self.output_dir` and print the next-step invocation."""

        view_path = os.path.join(self.output_dir, 'INDELS-VIEW-DATA.txt')
        items_add_path = os.path.join(self.output_dir, 'INDELS-ITEMS-ADDITIONAL.txt')
        tree_path = os.path.join(self.output_dir, 'INDELS-ITEMS-ORDER.newick')
        long_path = os.path.join(self.output_dir, 'INDELS-LOCI-LONG.txt')
        clusters_fa_path = os.path.join(self.output_dir, 'INDELS-CLUSTERS.fa')
        clusters_txt_path = os.path.join(self.output_dir, 'INDELS-CLUSTERS.txt')

        self.progress.new('Writing outputs')

        # view-data matrix: rows=loci, columns=samples
        self.progress.update('view-data matrix ...')
        with open(view_path, 'w') as f:
            f.write('item\t' + '\t'.join(self.samples) + '\n')
            for locus_id in sorted(self.matrix.keys()):
                row = self.matrix[locus_id]
                f.write(locus_id + '\t' + '\t'.join(f"{row[s]:.6f}" for s in self.samples) + '\n')

        # items-additional-data
        decoration_cols = [
            'cluster_id', 'locus_type', 'contig_name', 'mean_pos_in_contig', 'pos_range',
            'cluster_size_median_nt', 'n_loci_in_cluster', 'locus_size_median_nt', 'n_events', 'n_supporting_reads', 'n_samples_supporting',
            'in_orf', 'gene_caller_id',
            'gene_function_source', 'gene_function_accession', 'gene_function',
            'representative_id',
        ]
        self.progress.update('items-additional data ...')
        utils.store_dict_as_TAB_delimited_file(self.items_additional, items_add_path,
                                               headers=['item'] + decoration_cols)

        # items-order newick
        self.progress.update('items-order tree ...')
        if self.items_tree_newick:
            with open(tree_path, 'w') as f:
                f.write(self.items_tree_newick)
        else:
            tree_path = None

        # long-format audit
        self.progress.update('long-format audit ...')
        long_cols = ['indel_index', 'cluster_id', 'locus_id', 'sample_id', 'profile_db',
                     'contig_name', 'split_name', 'pos_in_contig', 'type', 'length',
                     'count', 'coverage', 'allele_frequency',
                     'corresponding_gene_call', 'reference']
        with open(long_path, 'w') as f:
            f.write('\t'.join(long_cols) + '\n')
            n_events = len(self.indels)
            self.progress.update('long-format audit ...')
            for idx in range(n_events):
                e = self.indels[idx]
                row = {c: e[c] for c in long_cols if c in _IndelRow.KEYS}
                row['indel_index'] = idx
                row['cluster_id'] = self._cluster_id_for_indel(idx)
                row['locus_id'] = self._idx_to_locus.get(idx, '')
                row['allele_frequency'] = (row['count'] / row['coverage']) if row['coverage'] > 0 else 0.0
                f.write('\t'.join(str(row.get(c, '')) for c in long_cols) + '\n')

        # INS rows set aside by --exclude-contig-edges
        edge_path = None
        if self.n_edge_rows:
            edge_path = os.path.join(self.output_dir, 'INDELS-CONTIG-EDGE-INS.txt')
            shutil.move(os.path.join(self.work_dir, 'INDELS-CONTIG-EDGE-INS.txt'), edge_path)

        # cluster-representative FASTA — one entry per final cluster (i.e. clusters that have at
        # least one locus surviving every filter). Header carries metadata for downstream BLAST /
        # HMM lookups; sequence is the mmseqs2 (or exact-match) cluster representative's actual
        # indel sequence, pulled from self.indels via the leading idx encoded in the rep id.
        self.progress.update('cluster-representative FASTA ...')
        final_cluster_ids = sorted({info['cluster_id'] for info in self.loci.values()})

        # one-pass count of loci per cluster so we don't do O(N_loci × N_clusters) scans below
        loci_per_cluster = Counter(info['cluster_id'] for info in self.loci.values())

        # both cluster files are written cluster by cluster: holding the metadata and representative
        # sequences of every final cluster at once does not fit in memory for tens of millions of clusters
        cluster_cols = ['length', 'cluster_type', 'n_supporting_reads', 'n_loci', 'n_samples_supporting', 'representative_id']
        with open(clusters_txt_path, 'w') as txt, open(clusters_fa_path, 'w') as fa:
            txt.write('\t'.join(['cluster_id'] + cluster_cols) + '\n')

            for cid in final_cluster_ids:
                member_idxs = self.clusters_to_members.get(cid, [])
                if not member_idxs:
                    continue

                rep_id = self.cluster_representatives.get(cid, '')
                rep_seq = ''
                if rep_id:
                    try:
                        rep_idx = int(rep_id.split('|', 1)[0])
                        rep_seq = self.indels[rep_idx].get('sequence', '') or ''
                    except (ValueError, IndexError):
                        pass

                # clusters without a representative sequence are left out of both files
                if not rep_seq:
                    continue

                type_counts = Counter(self.indels[i]['type'] for i in member_idxs)
                ins_n = type_counts.get('INS', 0)
                del_n = type_counts.get('DEL', 0)
                if ins_n and del_n:
                    cluster_type = 'mixed'
                elif ins_n:
                    cluster_type = 'INS'
                else:
                    cluster_type = 'DEL'

                meta = {'length': len(rep_seq),
                        'cluster_type': cluster_type,
                        'n_supporting_reads': sum(self.indels[i]['count'] for i in member_idxs),
                        'n_loci': loci_per_cluster.get(cid, 0),
                        'n_samples_supporting': len({self.indels[i]['sample_id'] for i in member_idxs}),
                        'representative_id': rep_id}
                txt.write('\t'.join([str(cid)] + [str(meta[c]) for c in cluster_cols]) + '\n')

                # the defline is just the cluster id so that this file can go straight into
                # `anvi-gen-contigs-database`; the rest of the metadata is in INDELS-CLUSTERS.txt
                fa.write(f">{cid}\n")
                for i in range(0, len(rep_seq), 70):
                    fa.write(rep_seq[i:i + 70] + "\n")

        self.progress.end()

        self.run.warning(None, header="OUTPUT", lc="green")
        self.run.info('View data (matrix)', view_path, mc='green')
        self.run.info('Items additional data', items_add_path, mc='green')
        if tree_path:
            self.run.info('Items-order tree (newick)', tree_path, mc='green')
        else:
            self.run.info('Items-order tree (newick)', 'skipped (< 3 loci, or more than --max-loci-for-tree)', mc='yellow')
        self.run.info('Per-event audit (long format)', long_path, mc='green')
        if edge_path:
            self.run.info('Contig-edge INS rows (set aside)', edge_path, mc='green')
        self.run.info('Cluster representatives (FASTA)', clusters_fa_path, mc='green')
        self.run.info('Cluster representatives (metadata)', clusters_txt_path, mc='green')

        manual_profile = os.path.join(self.output_dir, 'manual-profile.db')
        invocation = [
            'anvi-interactive --manual-mode',
            f'  -p {manual_profile}',
            f'  --view-data {view_path}',
            f'  --additional-layers {items_add_path}',
        ]
        if tree_path:
            invocation.append(f'  --tree {tree_path}')

        self.run.warning("\n".join(invocation),
                         header="Next step: visualize with anvi-interactive --manual",
                         lc="cyan", raw=True)
