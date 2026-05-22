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

        # output
        self.output_dir = A('output_dir')
        self.num_threads = A('num_threads', 1)
        self.just_do_it = A('just_do_it', False)

        # state populated by process()
        self.contigs_super = None
        self.samples = []                 # ordered list of sample names
        self.sample_to_profile_path = {}  # sample_id → profile-db path it lives in
        self.aux_data_caches = {}         # profile-db path → AuxiliaryDataForSplitCoverages (or None)
        self.indels = []                  # list of dicts (one per accepted indel row)
        self.clusters_to_members = {}     # {cluster_id: [indel idx, ...]}
        self.cluster_representatives = {} # {cluster_id: fasta_id_of_representative}
        self.loci = {}                    # {locus_id: {...}}
        self.locus_event_types = {}       # locus_id → 'INS' or 'DEL' (dominant type)
        self.matrix = {}                  # {locus_id: {sample: presence}}
        self.items_additional = {}        # {locus_id: {decoration_col: value, ...}}
        self.items_tree_newick = None
        self.long_format_rows = []


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

        self._init_contigs_super()
        self.collect_indels()
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
        """Instantiate ContigsSuperclass once and pre-load contig sequences + gene calls/functions."""
        import argparse
        self.progress.new('Initializing contigs database')
        self.progress.update('...')
        args = argparse.Namespace(contigs_db=self.contigs_db_path)
        self.contigs_super = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        self.contigs_super.init_contig_sequences()
        try:
            self.contigs_super.init_functions()
        except ConfigError:
            # functional annotations may not be present; that's OK.
            pass
        self.progress.end()


    def collect_indels(self):
        """Load indels rows from each profile-db, apply filters, store as list of dicts."""

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

        all_rows = []
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

            try:
                cursor = pdb.db._exec(f"SELECT {select_clause} FROM {t.indels_table_name} WHERE {where}")
                rows = cursor.fetchall()
            except Exception:
                rows = []
            finally:
                pdb.db.disconnect()

            db_basename = os.path.basename(path).replace('.db', '')

            for row in rows:
                sample_id, split_name, pos_in_contig, gene_id, reference, ev_type, sequence, length, count, coverage = row
                contig = splits_to_contigs.get(split_name)
                if contig is None:
                    # split is unknown to this contigs-db; skip
                    continue

                entry = {
                    'sample_id': sample_id,
                    'profile_db': db_basename,
                    'split_name': split_name,
                    'contig_name': contig,
                    'pos_in_contig': int(pos_in_contig),
                    'corresponding_gene_call': int(gene_id) if gene_id is not None else -1,
                    'reference': reference or '',
                    'type': ev_type,
                    'sequence': sequence or '',
                    'length': int(length),
                    'count': int(count),
                    'coverage': int(coverage),
                }
                all_rows.append(entry)

                if sample_id not in sample_seen:
                    sample_seen.add(sample_id)
                    sample_order.append(sample_id)

        self.indels = all_rows
        self.samples = sample_order
        self.progress.end()

        self.run.info('Indel rows passing filters', pp(len(self.indels)))
        self.run.info('Distinct samples', pp(len(self.samples)))


    def reconstruct_del_sequences(self):
        """Fill in `sequence` for DEL rows by slicing from the contig sequence."""

        n_reconstructed = 0
        n_skipped_mismatch = 0
        n_skipped_oob = 0
        kept = []

        for entry in self.indels:
            if entry['type'] != 'DEL':
                kept.append(entry)
                continue

            contig = entry['contig_name']
            pos = entry['pos_in_contig']
            length = entry['length']

            contig_seq = self.contigs_super.contig_sequences.get(contig, {}).get('sequence')
            if not contig_seq:
                n_skipped_oob += 1
                continue

            if pos < 0 or pos + length > len(contig_seq):
                n_skipped_oob += 1
                continue

            del_seq = contig_seq[pos:pos + length].upper()

            if entry['reference'] and del_seq[:1].upper() != entry['reference'].upper():
                # reference mismatch — likely a position-semantics edge case; warn and skip
                n_skipped_mismatch += 1
                continue

            entry['sequence'] = del_seq
            kept.append(entry)
            n_reconstructed += 1

        self.indels = kept

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

        # bin indels by tier
        bins = defaultdict(list)
        for idx, entry in enumerate(self.indels):
            for label, max_excl, min_id in tiers:
                if entry['length'] < max_excl:
                    bins[(label, min_id)].append(idx)
                    break
            else:
                bins[(tiers[-1][0], tiers[-1][2])].append(idx)

        self.run.warning(None, header="CLUSTERING (length-tiered)", lc="green")

        all_clusters = {}
        all_reps = {}
        cluster_counter = 0
        work_dir = filesnpaths.get_temp_directory_path()

        MMSEQS_MIN_LEN = 15  # default nucleotide k-mer size in mmseqs2

        for (label, min_id), idxs in bins.items():
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

            for sub_cid, members in sub_clusters.items():
                cluster_counter += 1
                cluster_id = f"C_{cluster_counter:08d}"
                # member ids look like '<idx>|sample|contig|pos|type|length' — pull the leading idx back out
                member_idxs = [int(name.split('|', 1)[0]) for name in members]
                all_clusters[cluster_id] = member_idxs
                all_reps[cluster_id] = m.representatives.get(sub_cid)

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

            for contig, idxs in by_contig.items():
                # sort by pos_in_contig, sweep, merge while within tolerance of running mean
                idxs.sort(key=lambda i: self.indels[i]['pos_in_contig'])

                current_members = []
                current_positions = []

                def flush():
                    nonlocal locus_counter, current_members, current_positions
                    if not current_members:
                        return
                    locus_counter += 1
                    locus_id = f"locus_{locus_counter:08d}"
                    loci[locus_id] = {
                        'cluster_id': cluster_id,
                        'contig_name': contig,
                        'positions': list(current_positions),
                        'mean_pos': int(round(sum(current_positions) / len(current_positions))),
                        'min_pos': min(current_positions),
                        'max_pos': max(current_positions),
                        'members': list(current_members),
                    }
                    current_members = []
                    current_positions = []

                for idx in idxs:
                    pos = self.indels[idx]['pos_in_contig']
                    if current_positions:
                        running_mean = sum(current_positions) / len(current_positions)
                        if pos - running_mean > self.position_tolerance:
                            flush()
                    current_members.append(idx)
                    current_positions.append(pos)

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


    def _detection_in_region(self, sample, contig, pos_start, length):
        """Return the fraction of positions in [pos_start, pos_start+length) on `contig`
        that have ≥1x coverage in `sample`, looked up via the auxiliary-data db.

        Returns None if the lookup can't be performed (no aux data, unknown sample, etc.)."""

        profile_path = self.sample_to_profile_path.get(sample)
        if not profile_path:
            return None
        aux = self.aux_data_caches.get(profile_path)
        if aux is None:
            return None

        # Find the split containing pos_start
        split_name = None
        split_start = None
        for sn, sinfo in self.contigs_super.splits_basic_info.items():
            if sinfo['parent'] == contig and sinfo['start'] <= pos_start < sinfo['end']:
                split_name = sn
                split_start = sinfo['start']
                break
        if split_name is None:
            return None

        try:
            split_cov = aux.get(split_name)
        except Exception:
            return None
        if not split_cov or sample not in split_cov:
            return None

        rel_start = pos_start - split_start
        cov_slice = split_cov[sample][rel_start:rel_start + length]
        # detection is computed against the requested length, not the slice length, so positions
        # past the split boundary count as zero-coverage. For 20 kb splits and sub-kb indels
        # this is essentially never tripped.
        if length <= 0:
            return 0.0
        return float(np.count_nonzero(cov_slice)) / float(length)


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

        for locus_id, info in self.loci.items():
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
                        # DEL locus: did this sample actually have coverage in the DELd region?
                        detection = self._detection_in_region(sample, info['contig_name'], rep_pos, rep_length)
                        if detection is None:
                            # aux lookup unavailable -> conservative default of 0.0
                            row[sample] = 0.0
                        elif detection >= DETECTION_THRESHOLD_FOR_PRESENCE:
                            row[sample] = 1.0
                        else:
                            row[sample] = 0.0

            self.matrix[locus_id] = row

        if any_del_locus:
            self._close_aux_data_caches()


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

        for locus_id, info in self.loci.items():
            members = [self.indels[i] for i in info['members']]
            sample_counts = Counter(m['sample_id'] for m in members)
            lengths = [m['length'] for m in members]
            n_supporting_reads = sum(m['count'] for m in members)

            mean_pos = info['mean_pos']
            gene_caller_id = -1
            if has_genes_called:
                try:
                    gene_caller_id = self.contigs_super.nt_position_to_gene_caller_id(info['contig_name'], mean_pos)
                except ConfigError:
                    gene_caller_id = -1

            source, accession, function = self._gene_function_for(gene_caller_id)

            self.items_additional[locus_id] = {
                'cluster_id': info['cluster_id'],
                'locus_type': self.locus_event_types.get(locus_id, ''),
                'contig_name': info['contig_name'],
                'mean_pos_in_contig': mean_pos,
                'pos_range': f"{info['min_pos']}-{info['max_pos']}",
                'cluster_size_median_nt': int(statistics.median(lengths)),
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


    def build_items_tree(self):
        """Hierarchical clustering of the matrix → Newick string."""
        if len(self.matrix) < 3:
            self.items_tree_newick = None
            return
        self.items_tree_newick = clustering.get_newick_tree_data_for_dict(
            self.matrix, distance='euclidean', linkage='ward', zero_fill_missing=True)


    def build_long_format(self):
        """Per-event audit trail."""
        # locus lookup: idx -> locus_id
        idx_to_locus = {}
        for locus_id, info in self.loci.items():
            for idx in info['members']:
                idx_to_locus[idx] = locus_id

        for idx, e in enumerate(self.indels):
            row = dict(e)
            row.pop('sequence', None)  # sequence is too big for the audit file
            row['indel_index'] = idx
            row['cluster_id'] = self._cluster_id_for_indel(idx)
            row['locus_id'] = idx_to_locus.get(idx, '')
            row['allele_frequency'] = (e['count'] / e['coverage']) if e['coverage'] > 0 else 0.0
            self.long_format_rows.append(row)


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

        # view-data matrix: rows=loci, columns=samples
        with open(view_path, 'w') as f:
            f.write('item\t' + '\t'.join(self.samples) + '\n')
            for locus_id in sorted(self.matrix.keys()):
                row = self.matrix[locus_id]
                f.write(locus_id + '\t' + '\t'.join(f"{row[s]:.6f}" for s in self.samples) + '\n')

        # items-additional-data
        decoration_cols = [
            'cluster_id', 'locus_type', 'contig_name', 'mean_pos_in_contig', 'pos_range',
            'cluster_size_median_nt', 'n_events', 'n_supporting_reads', 'n_samples_supporting',
            'in_orf', 'gene_caller_id',
            'gene_function_source', 'gene_function_accession', 'gene_function',
            'representative_id',
        ]
        utils.store_dict_as_TAB_delimited_file(self.items_additional, items_add_path,
                                               headers=['item'] + decoration_cols)

        # items-order newick
        if self.items_tree_newick:
            with open(tree_path, 'w') as f:
                f.write(self.items_tree_newick)
        else:
            tree_path = None

        # long-format audit
        long_cols = ['indel_index', 'cluster_id', 'locus_id', 'sample_id', 'profile_db',
                     'contig_name', 'split_name', 'pos_in_contig', 'type', 'length',
                     'count', 'coverage', 'allele_frequency',
                     'corresponding_gene_call', 'reference']
        with open(long_path, 'w') as f:
            f.write('\t'.join(long_cols) + '\n')
            for r in self.long_format_rows:
                f.write('\t'.join(str(r.get(c, '')) for c in long_cols) + '\n')

        # cluster-representative FASTA — one entry per final cluster (i.e. clusters that have at
        # least one locus surviving every filter). Header carries metadata for downstream BLAST /
        # HMM lookups; sequence is the mmseqs2 (or exact-match) cluster representative's actual
        # indel sequence, pulled from self.indels via the leading idx encoded in the rep id.
        final_cluster_ids = sorted({info['cluster_id'] for info in self.loci.values()})
        cluster_meta = {}
        for cid in final_cluster_ids:
            member_idxs = self.clusters_to_members.get(cid, [])
            if not member_idxs:
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
            samples = {self.indels[i]['sample_id'] for i in member_idxs}
            n_loci = sum(1 for info in self.loci.values() if info['cluster_id'] == cid)
            n_supporting_reads = sum(self.indels[i]['count'] for i in member_idxs)
            rep_id = self.cluster_representatives.get(cid, '')
            rep_seq = ''
            rep_length = 0
            if rep_id:
                try:
                    rep_idx = int(rep_id.split('|', 1)[0])
                    rep_seq = self.indels[rep_idx].get('sequence', '') or ''
                    rep_length = len(rep_seq)
                except (ValueError, IndexError):
                    pass
            cluster_meta[cid] = {
                'representative_id': rep_id,
                'sequence': rep_seq,
                'length': rep_length,
                'cluster_type': cluster_type,
                'n_supporting_reads': n_supporting_reads,
                'n_loci': n_loci,
                'n_samples_supporting': len(samples),
            }

        with open(clusters_fa_path, 'w') as f:
            for cid in final_cluster_ids:
                m = cluster_meta.get(cid)
                if not m or not m['sequence']:
                    continue
                header = (f">{cid} length={m['length']} cluster_type={m['cluster_type']} "
                          f"n_supporting_reads={m['n_supporting_reads']} "
                          f"n_loci={m['n_loci']} n_samples={m['n_samples_supporting']} "
                          f"representative={m['representative_id']}")
                f.write(header + "\n")
                seq = m['sequence']
                for i in range(0, len(seq), 70):
                    f.write(seq[i:i + 70] + "\n")

        self.run.warning(None, header="OUTPUT", lc="green")
        self.run.info('View data (matrix)', view_path, mc='green')
        self.run.info('Items additional data', items_add_path, mc='green')
        if tree_path:
            self.run.info('Items-order tree (newick)', tree_path, mc='green')
        else:
            self.run.info('Items-order tree (newick)', 'skipped (< 3 loci)', mc='yellow')
        self.run.info('Per-event audit (long format)', long_path, mc='green')
        self.run.info('Cluster representatives (FASTA)', clusters_fa_path, mc='green')

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
