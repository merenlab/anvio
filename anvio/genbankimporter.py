# -*- coding: utf-8
"""GenBank → contigs-database sequence-features importer.

`anvi-import-genbank-features` is the user-facing CLI; this module hosts the
implementation. The two-pass design keeps relationship resolution simple by
guaranteeing every candidate parent has been parsed before any child is asked
to find one. All database writes happen in a single transaction at the end of
pass 2 via `TablesForSequenceFeatures.populate_features`.
"""

import os
import re
import gzip
import io
import hashlib

from Bio import SeqIO
from Bio.SeqFeature import (CompoundLocation, BeforePosition, AfterPosition, BetweenPosition)

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.tables as t
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tables.sequencefeatures import TablesForSequenceFeatures, BUILTIN_FEATURE_TYPES


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10']
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


VALID_SOURCE_NAME_RE = re.compile(r'^[A-Za-z0-9_-]+$')
BUILTIN_FEATURE_TYPE_NAMES = {row[0] for row in BUILTIN_FEATURE_TYPES}

# Qualifiers that get promoted to dedicated columns and therefore should not also
# be written to `feature_qualifiers`. Only applies to CDS features.
CDS_DEDICATED_QUALIFIERS = {'translation', 'codon_start', 'transl_table'}

# Feature-type pairs (child_type → ordered list of parent-type preferences plus the
# relationship label that pairing produces). Only types we resolve relationships for.
PARENT_RULES = {
    'CDS':    [('mRNA', 'part_of'), ('gene', 'derives_from')],
    'mRNA':   [('gene', 'part_of')],
    'exon':   [('mRNA', 'part_of')],
    'intron': [('mRNA', 'part_of')],
}


def compute_feature_id(contig, feature_type, source, start, stop, direction):
    """Compute the 16-character hex `feature_id` for one row.

    The hash input is the TAB-joined tuple. The TAB is mandatory — concatenating
    fields without a separator would let adjacent coordinates blur into hash
    collisions (start=12,stop=345 vs start=123,stop=45 both produce '12345').
    Before hashing we assert no field contains TAB itself; silent hash divergence
    would be a far worse failure mode than an explicit error.
    """

    direction_str = direction if direction is not None else ''
    fields = [contig, feature_type, source, str(start), str(stop), direction_str]
    for field in fields:
        if '\t' in field:
            raise ConfigError(f"compute_feature_id refuses to hash an input field that contains a TAB character: '{field!r}'. "
                              f"This would create silent hash divergence. Please report this as an anvi'o bug if a TAB ever appears in a contig name, "
                              f"feature type, or source name.")
    return hashlib.sha224('\t'.join(fields).encode('utf-8')).hexdigest()[:16]


class GenbankFeatureImporter:
    """Two-pass GenBank feature importer.

    Pass 1 parses the entire GenBank file into an in-memory representation
    (features, qualifiers, CDS-specific data) and dedups hash collisions
    that arise from duplicate rows in the input file. Pass 2 reconciles each
    single-segment gene against the existing `genes_in_contigs` gene calls,
    populates `external_id` from `locus_tag` where available, resolves parent
    relationships per contig, and calls `TablesForSequenceFeatures.populate_features`
    once to commit everything atomically.

    Validation and warnings live in this class (not in the SQL layer):
    direction in {'f','r',None}, partial flags in {0,1,None}, `codon_start` in
    {1,2,3,None}, LOCUS-to-contig matching, sequence-length parity, and the
    structural invariants for multi-segment groups.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path   = A('contigs_db')
        self.input_genbank_path = A('input_genbank')
        self.source_name        = A('source_name') or 'genbank_import'
        # --force is captured by anvi'o's global flag handling at module-import time
        self.force = anvio.FORCE

        if not self.contigs_db_path:
            raise ConfigError("Anvi'o needs a contigs database to import features into. Please provide one with `-c / --contigs-db`.")
        if not self.input_genbank_path:
            raise ConfigError("Anvi'o needs a GenBank file as input. Please provide one with `-i / --input-genbank`.")
        if not VALID_SOURCE_NAME_RE.match(self.source_name):
            raise ConfigError(f"The --source-name value '{self.source_name}' is not acceptable. Source names must match the regex "
                              f"`^[A-Za-z0-9_-]+$` because the value is hashed into every `feature_id` with TAB as a separator, "
                              f"and allowing arbitrary characters could let a user accidentally create separator-collisions. "
                              f"Pick something like 'genbank_import' or 'ncbi-refseq'.")

        filesnpaths.is_file_exists(self.input_genbank_path)
        utils.is_contigs_db(self.contigs_db_path)

        # state populated by setup / pass 1 / pass 2
        self.contig_lengths = {}
        self.contig_is_circular = {}
        self.gene_calls_lookup = {}

        self.features = []                  # list of dicts ready for populate_features
        self.relationships = []             # list of dicts ready for populate_features
        self.qualifiers = []                # list of dicts ready for populate_features
        self.cds_specific = []              # list of dicts ready for populate_features

        # per-feature metadata kept around for pass 2 relationship resolution
        # (not written to the database directly — derived from source SeqFeatures
        # in pass 1 and consumed in pass 2)
        self.feature_meta = {}              # feature_id → {locus_tag, gene_q, contig, ftype, canonical_fid, group_extent}
        self.groups = {}                    # canonical_fid → {'segments': [fids in transcription order],
                                            #                  'contig', 'feature_type', 'min_start', 'max_stop',
                                            #                  'locus_tag', 'gene_q'}

        # summary counters
        self.counts = {
            'features_inserted': 0,
            'relationships_inserted': 0,
            'qualifiers_inserted': 0,
            'new_non_builtin_types': 0,
            'duplicates_skipped': 0,
            'skipped_no_extent': 0,
            'skipped_malformed_origin_linear': 0,
            'skipped_multiseg_undefined_strand': 0,
            'gcid_matched': 0,
            'gcid_unmatched': 0,
            'cds_without_parents': 0,
            'invalid_codon_start_nulled': 0,
        }
        self.types_seen = {}                # feature_type → count for summary breakdown
        self.new_non_builtin_types = set()


    def process(self):
        """Top-level entry point. Runs setup, pass 1, pass 2, then the summary."""

        self.progress.new("Importing GenBank features")
        try:
            self.progress.update("Loading contigs database metadata...")
            self._setup()
            self.progress.update("Parsing GenBank file (pass 1)...")
            self._parse_genbank_pass1()
            self.progress.update("Resolving relationships (pass 2)...")
            self._resolve_pass2()
            self.progress.update("Writing features to the database...")
            self._write()
        finally:
            self.progress.end()

        self._emit_summary()


    # -----------------------------------------------------------------
    # Setup
    # -----------------------------------------------------------------

    def _setup(self):
        """Open the contigs database, read contig metadata + existing gene calls, and
        decide whether `--force` is needed up front (we abort here rather than after
        pass 1 if rows for the chosen source already exist)."""

        # Use the high-level access class for version handling and meta. We will not
        # hold the connection open across pass 1; ContigsSuperclass / TablesForSequenceFeatures
        # open their own connections as needed.
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))

        try:
            contig_rows = contigs_db.db.get_table_as_dict(t.contigs_info_table_name)
            for contig_name, row in contig_rows.items():
                self.contig_lengths[contig_name] = row['length']

            gene_rows = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name, error_if_no_data=False)
            for _gcid_unused, gene in gene_rows.items():
                # the first column of genes_in_contigs is gene_callers_id; get_table_as_dict
                # keys by that, but we also need the gcid as a value field. The dict value
                # already contains every other column (contig, start, stop, direction, ...).
                # We need the gcid; reach into the row's natural key.
                pass
            # Rebuild gene_calls_lookup the right way — gcid is the dict key, the rest is value.
            for gcid, gene in gene_rows.items():
                key = (gene['contig'], int(gene['start']), int(gene['stop']), gene['direction'])
                self.gene_calls_lookup[key] = int(gcid)

            existing_in_source = contigs_db.db.get_single_column_from_table(
                t.contigs_sequence_features_table_name,
                'feature_id',
                where_clause=f"source = '{self.source_name}'",
            )
        finally:
            contigs_db.disconnect()

        if existing_in_source and not self.force:
            raise ConfigError(f"There are already {len(existing_in_source)} sequence-feature rows in this contigs database with "
                              f"`source = '{self.source_name}'`. Anvi'o refuses to silently merge or overwrite them. Either pass "
                              f"`--force` to delete and re-import this source (rows belonging to OTHER source names are unaffected "
                              f"by --force), or pick a different `--source-name`.")


    # -----------------------------------------------------------------
    # Pass 1: parse the GenBank file
    # -----------------------------------------------------------------

    def _parse_genbank_pass1(self):
        """Buffer the GenBank records, validate LOCUS↔contig parity, then iterate features."""

        records = list(self._iter_genbank_records())
        if not records:
            raise ConfigError(f"The GenBank file at '{self.input_genbank_path}' did not parse into any records. It may be empty, "
                              f"malformed, or in a format Biopython cannot read.")

        unmatched_loci = [r.name for r in records if r.name not in self.contig_lengths]
        if unmatched_loci:
            raise ConfigError(f"The following GenBank LOCUS name(s) do not correspond to any contig in your contigs database "
                              f"(matched against `contigs_basic_info`): {', '.join(unmatched_loci)}. Importing zero features from "
                              f"those records would silently swallow data, so anvi'o is bailing out rather than continuing. "
                              f"Either fix the GenBank file to match your contigs DB exactly, or generate a new contigs DB from "
                              f"the same FASTA underlying this GenBank file.")

        for record in records:
            expected_len = self.contig_lengths[record.name]
            actual_len = len(record.seq)
            if actual_len > 0 and actual_len != expected_len:
                self.run.warning(f"GenBank LOCUS '{record.name}' reports a sequence length of {actual_len:,} bp but the contigs "
                                 f"database has the same contig at {expected_len:,} bp. Coordinates may still be meaningful in the "
                                 f"matching prefix, but this is suspicious — the FASTA used to build the contigs DB may differ from "
                                 f"the GenBank's sequence.")

        for record in records:
            self.contig_is_circular[record.name] = (record.annotations.get('topology', '').lower() == 'circular')
            self._iter_features_in_record(record)


    def _iter_genbank_records(self):
        """Yield SeqRecord objects from the input GenBank file, transparently handling gzip."""

        if self.input_genbank_path.endswith('.gz'):
            handle = io.TextIOWrapper(gzip.open(self.input_genbank_path, 'rb'), encoding='utf-8')
        else:
            handle = open(self.input_genbank_path, 'r')
        try:
            for record in SeqIO.parse(handle, 'genbank'):
                yield record
        finally:
            handle.close()


    def _iter_features_in_record(self, record):
        contig = record.name
        contig_length = self.contig_lengths[contig]
        is_circular = self.contig_is_circular[contig]

        for seqfeature in record.features:
            ftype = seqfeature.type
            if ftype == 'source':
                continue

            loc = seqfeature.location
            if loc is None:
                continue

            # reject BetweenPosition endpoints (`100^101` style) and zero-extent features
            if isinstance(loc.start, BetweenPosition) or isinstance(loc.end, BetweenPosition):
                self.run.warning(f"Skipping a feature ({ftype}) on contig '{contig}' whose location uses BetweenPosition "
                                 f"(`100^101` notation). Such features have no biological extent we can represent.")
                self.counts['skipped_no_extent'] += 1
                continue

            parts = list(loc.parts) if hasattr(loc, 'parts') else [loc]
            if any(int(p.start) == int(p.end) for p in parts):
                self.run.warning(f"Skipping a feature ({ftype}) on contig '{contig}' that resolves to zero genomic extent.")
                self.counts['skipped_no_extent'] += 1
                continue

            direction = self._biopython_strand_to_direction(loc.strand)
            multi_segment = len(parts) > 1

            if multi_segment and direction is None:
                self.run.warning(f"Skipping a multi-segment feature ({ftype}) on contig '{contig}' that has an undefined strand. "
                                 f"Transcription order is undefined for such features.")
                self.counts['skipped_multiseg_undefined_strand'] += 1
                continue

            origin_crossing = self._is_origin_crossing(parts, contig_length, is_circular)
            if multi_segment and not is_circular and self._is_origin_crossing(parts, contig_length, True):
                # malformed: join() spans the origin but the contig is linear
                self.run.warning(f"Skipping a feature ({ftype}) on contig '{contig}': its `join(...)` location looks like an "
                                 f"origin-crossing feature but the contig is linear. This is almost certainly malformed.")
                self.counts['skipped_malformed_origin_linear'] += 1
                continue

            ordered_parts = self._transcription_order(parts, direction, origin_crossing)

            # partial flags: derived from the whole feature's 5' and 3' endpoints, then
            # replicated on every segment row.
            partial_5p, partial_3p = self._partial_flags(ordered_parts, direction)

            # extract source qualifiers once for the whole feature
            quals_dict = self._extract_qualifiers(seqfeature)
            locus_tag = quals_dict.get('locus_tag', [None])[0]
            gene_q   = quals_dict.get('gene',      [None])[0]

            # generate one row per segment
            new_segment_ids = []
            for idx, segment in enumerate(ordered_parts):
                seg_start = int(segment.start)
                seg_stop  = int(segment.end)
                if not (0 <= seg_start < seg_stop <= contig_length):
                    self.run.warning(f"Skipping segment {idx} of a {ftype} on contig '{contig}': coordinates [{seg_start}, {seg_stop}) "
                                     f"are out of bounds for a contig of length {contig_length}.")
                    self.counts['skipped_no_extent'] += 1
                    new_segment_ids = []
                    break
                fid = compute_feature_id(contig, ftype, self.source_name, seg_start, seg_stop, direction)
                row = {
                    'feature_id':         fid,
                    'contig':             contig,
                    'feature_type':       ftype,
                    'source':             self.source_name,
                    'start':              seg_start,
                    'stop':               seg_stop,
                    'direction':          direction,
                    'partial_fiveprime':  partial_5p,
                    'partial_threeprime': partial_3p,
                    'feature_group_id':   None,
                    'segment_order':      None,
                    'external_id':        locus_tag,
                    'gene_callers_id':    None,
                }
                if multi_segment:
                    row['segment_order'] = idx
                # we'll fill feature_group_id below after we know the canonical id
                if not self._accept_or_dedup_collision(row):
                    continue
                new_segment_ids.append(fid)

            if not new_segment_ids:
                continue

            canonical_fid = new_segment_ids[0]
            if multi_segment:
                for fid in new_segment_ids:
                    # find the just-added row in self.features and set feature_group_id
                    for row in reversed(self.features):
                        if row['feature_id'] == fid:
                            row['feature_group_id'] = canonical_fid
                            break

            # group bookkeeping
            self.groups[canonical_fid] = {
                'segments':    list(new_segment_ids),
                'contig':      contig,
                'feature_type': ftype,
                'min_start':   min(int(p.start) for p in ordered_parts),
                'max_stop':    max(int(p.end)   for p in ordered_parts),
                'locus_tag':   locus_tag,
                'gene_q':      gene_q,
            }

            for fid in new_segment_ids:
                self.feature_meta[fid] = {
                    'canonical_fid': canonical_fid,
                    'contig':        contig,
                    'feature_type':  ftype,
                }

            # per-feature qualifier rows (replicated onto every segment row)
            for fid in new_segment_ids:
                self._record_qualifiers(fid, ftype, quals_dict)

            # CDS-specific data — one row per segment; non-canonical rows get NULLs
            if ftype == 'CDS':
                self._record_cds(canonical_fid, new_segment_ids, quals_dict)

            self.types_seen[ftype] = self.types_seen.get(ftype, 0) + 1
            if ftype not in BUILTIN_FEATURE_TYPE_NAMES:
                self.new_non_builtin_types.add(ftype)


    def _biopython_strand_to_direction(self, strand):
        if strand == 1:
            return 'f'
        if strand == -1:
            return 'r'
        return None


    def _is_origin_crossing(self, parts, contig_length, is_circular):
        if not is_circular or len(parts) < 2:
            return False
        starts_at_zero = any(int(p.start) == 0 for p in parts)
        ends_at_length = any(int(p.end) == contig_length for p in parts)
        return starts_at_zero and ends_at_length


    def _transcription_order(self, parts, direction, origin_crossing):
        """Return parts in transcription order. See module docstring & prompt for rules:
        - origin-crossing on circular contig: trust the GenBank `join` order (parts as-is).
        - forward strand otherwise: sort by genomic start ascending.
        - reverse strand otherwise: sort by genomic start ascending then reverse.
        """

        if origin_crossing:
            return list(parts)
        sorted_parts = sorted(parts, key=lambda p: int(p.start))
        if direction == 'r':
            return list(reversed(sorted_parts))
        return sorted_parts


    def _partial_flags(self, ordered_parts, direction):
        """Compute (partial_fiveprime, partial_threeprime) from the whole feature's endpoints.

        Forward: 5' = leftmost.start, 3' = rightmost.end.
        Reverse: 5' = rightmost.end (becomes leftmost in transcription), 3' = leftmost.start.
        For unstranded features we only flag a partial if a BeforePosition/AfterPosition is
        present anywhere; we cannot distinguish which side is 5' vs 3'.
        """

        if not ordered_parts:
            return 0, 0

        if direction == 'f':
            five_end = ordered_parts[0].start
            three_end = ordered_parts[-1].end
        elif direction == 'r':
            five_end = ordered_parts[0].end   # rightmost.end in transcription order
            three_end = ordered_parts[-1].start
        else:
            # unstranded: flag a partial if any endpoint is fuzzy
            any_before = any(isinstance(p.start, BeforePosition) for p in ordered_parts)
            any_after  = any(isinstance(p.end,   AfterPosition)  for p in ordered_parts)
            return (1 if any_before else 0), (1 if any_after else 0)

        partial_5p = 1 if isinstance(five_end,  BeforePosition) else 0
        partial_3p = 1 if isinstance(three_end, AfterPosition)  else 0
        return partial_5p, partial_3p


    def _extract_qualifiers(self, seqfeature):
        """Strip whitespace from `translation`; preserve whitespace elsewhere; coerce to strings."""

        out = {}
        for key, values in seqfeature.qualifiers.items():
            cleaned = []
            for v in values:
                s = str(v)
                if key == 'translation':
                    s = re.sub(r'\s+', '', s)
                cleaned.append(s)
            out[key] = cleaned
        return out


    def _accept_or_dedup_collision(self, row):
        """Append `row` to self.features unless it duplicates an existing row.

        Exact-match duplicates emit a warning and are silently dropped. Hash collisions
        on different (contig, ftype, source, start, stop, direction) tuples are
        astronomically unlikely with 16-char SHA-224 truncation at first-PR scales —
        if one occurs we raise rather than silently merge data."""

        fid = row['feature_id']
        for existing in self.features:
            if existing['feature_id'] != fid:
                continue
            same = all(existing[k] == row[k] for k in ('contig', 'feature_type', 'source', 'start', 'stop', 'direction'))
            if same:
                self.run.warning(f"Duplicate feature in GenBank input: contig={row['contig']}, type={row['feature_type']}, "
                                 f"[{row['start']}, {row['stop']}), direction={row['direction']}. Dropping the duplicate.")
                self.counts['duplicates_skipped'] += 1
                return False
            raise ConfigError(f"Hash collision between two distinct features mapping to feature_id={fid}. Expected: at first-PR "
                              f"scales this is astronomically unlikely (16 hex chars of SHA-224). Either the GenBank file has "
                              f"genuinely pathological data or anvi'o has a bug. Refusing to silently merge.")
        self.features.append(row)
        return True


    def _record_qualifiers(self, fid, ftype, quals_dict):
        """Emit one `feature_qualifiers` row per (key, position). For CDSs the three
        dedicated keys (translation, codon_start, transl_table) live in CDS_features
        instead and are not duplicated here."""

        for key, values in quals_dict.items():
            if ftype == 'CDS' and key in CDS_DEDICATED_QUALIFIERS:
                continue
            for pos, value in enumerate(values):
                self.qualifiers.append({
                    'feature_id': fid,
                    'key':        key,
                    'value':      value,
                    'position':   pos,
                })


    def _record_cds(self, canonical_fid, segment_ids, quals_dict):
        """Emit one `CDS_features` row per segment. The canonical row carries the actual
        values; all other segment rows store NULL across the three dedicated columns.

        `codon_start` is validated against {1,2,3}; invalid values are nulled and warned."""

        codon_start = self._parse_codon_start(quals_dict.get('codon_start'))
        translation = quals_dict.get('translation', [None])[0]
        transl_table = self._parse_transl_table(quals_dict.get('transl_table'))

        for fid in segment_ids:
            if fid == canonical_fid:
                self.cds_specific.append({
                    'feature_id':   fid,
                    'codon_start':  codon_start if codon_start is not None else 1,
                    'translation':  translation,
                    'transl_table': transl_table if transl_table is not None else 1,
                })
            else:
                self.cds_specific.append({
                    'feature_id':   fid,
                    'codon_start':  None,
                    'translation':  None,
                    'transl_table': None,
                })


    def _parse_codon_start(self, raw):
        if raw is None:
            return None
        try:
            v = int(raw[0])
        except (ValueError, TypeError, IndexError):
            self.run.warning(f"GenBank codon_start qualifier '{raw}' is not a valid integer; storing NULL.")
            self.counts['invalid_codon_start_nulled'] += 1
            return None
        if v not in (1, 2, 3):
            self.run.warning(f"GenBank codon_start qualifier {v} is not in {{1, 2, 3}}; storing NULL.")
            self.counts['invalid_codon_start_nulled'] += 1
            return None
        return v


    def _parse_transl_table(self, raw):
        if raw is None:
            return None
        try:
            return int(raw[0])
        except (ValueError, TypeError, IndexError):
            return None


    # -----------------------------------------------------------------
    # Pass 2: relationships, GCID reconciliation, external_id
    # -----------------------------------------------------------------

    def _resolve_pass2(self):
        # GCID reconciliation only applies to single-segment gene rows.
        for row in self.features:
            if row['feature_type'] != 'gene':
                continue
            if row['feature_group_id'] is not None:
                # multi-segment gene → leave NULL; this is expected, not a data issue.
                continue
            if row['direction'] is None:
                # genes_in_contigs gene calls always have a direction; no match possible.
                continue
            key = (row['contig'], int(row['start']), int(row['stop']), row['direction'])
            gcid = self.gene_calls_lookup.get(key)
            if gcid is not None:
                row['gene_callers_id'] = gcid
                self.counts['gcid_matched'] += 1
            else:
                self.run.warning(f"No gene call matches GenBank gene row at contig={row['contig']}, [{row['start']}, {row['stop']}), "
                                 f"direction={row['direction']}. Importing with gene_callers_id = NULL.")
                self.counts['gcid_unmatched'] += 1

        # Relationship resolution — per contig, then per child group.
        groups_by_contig = {}
        for canonical_fid, info in self.groups.items():
            groups_by_contig.setdefault(info['contig'], []).append(canonical_fid)

        for contig, canonical_fids in groups_by_contig.items():
            self._resolve_relationships_on_contig(canonical_fids)


    def _resolve_relationships_on_contig(self, canonical_fids):
        """Apply the prompt's matching rules to assign parent relationships."""

        # Index by type for the contig
        by_type = {}
        for cfid in canonical_fids:
            by_type.setdefault(self.groups[cfid]['feature_type'], []).append(cfid)

        for cfid in canonical_fids:
            info = self.groups[cfid]
            ftype = info['feature_type']
            if ftype not in PARENT_RULES:
                continue

            chosen_parents = []   # list of (parent_canonical_fid, relationship)

            for parent_type, label in PARENT_RULES[ftype]:
                candidates = by_type.get(parent_type, [])
                matches = self._match_parents(info, candidates)
                if matches:
                    chosen_parents = [(p, label) for p in matches]
                    break  # only the first preference produces parents (CDS-mRNA preferred over CDS-gene)

            if not chosen_parents and ftype == 'CDS':
                self.run.warning(f"CDS '{cfid}' on contig '{info['contig']}' has no mRNA or gene parent that satisfied both the "
                                 f"qualifier-matching rule and coordinate containment. Importing with no parent links.")
                self.counts['cds_without_parents'] += 1
                continue
            if not chosen_parents:
                self.run.warning(f"{ftype} '{cfid}' on contig '{info['contig']}' has no {' or '.join(p for p, _ in PARENT_RULES[ftype])} parent.")
                continue

            for parent_cfid, label in chosen_parents:
                if parent_cfid == cfid:
                    continue  # never self-link
                # every segment of the child → one row referencing the canonical parent
                for seg_fid in info['segments']:
                    self.relationships.append({
                        'child_feature_id':  seg_fid,
                        'parent_feature_id': parent_cfid,
                        'relationship':      label,
                    })


    def _match_parents(self, child_info, candidate_cfids):
        """Apply the three-rule precedence (locus_tag → gene qualifier → coord containment),
        then filter by spanning-extent containment. Returns the list of parent canonical
        feature_ids that match. The matching rule selects candidates; containment finalizes."""

        ctag = child_info['locus_tag']
        cgene = child_info['gene_q']
        matched = []

        for pfid in candidate_cfids:
            pinfo = self.groups[pfid]
            pt    = pinfo['locus_tag']
            pg    = pinfo['gene_q']

            # the three rules, in order — first applicable rule is final for this pair
            if ctag is not None and pt is not None:
                if ctag != pt:
                    continue
                # locus_tag matches: still requires containment below
            elif cgene is not None and pg is not None:
                if cgene != pg:
                    continue
            else:
                # neither qualifier is shared: containment alone decides
                pass

            if pinfo['min_start'] <= child_info['min_start'] and pinfo['max_stop'] >= child_info['max_stop']:
                matched.append(pfid)

        return matched


    # -----------------------------------------------------------------
    # Write
    # -----------------------------------------------------------------

    def _write(self):
        tables_for_features = TablesForSequenceFeatures(self.contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        tables_for_features.populate_features(
            self.features,
            self.relationships,
            self.qualifiers,
            self.cds_specific,
            source_name=self.source_name,
            force=self.force,
        )

        self.counts['features_inserted'] = len(self.features)
        self.counts['relationships_inserted'] = len(self.relationships)
        self.counts['qualifiers_inserted'] = len(self.qualifiers)
        self.counts['new_non_builtin_types'] = len(self.new_non_builtin_types)


    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------

    def _emit_summary(self):
        self.run.info_single("Sequence-features import complete.", nl_before=1, nl_after=1, mc='green')
        self.run.info("Source name",              self.source_name)
        self.run.info("Features inserted",        self.counts['features_inserted'])
        self.run.info("Relationships inserted",   self.counts['relationships_inserted'])
        self.run.info("Qualifiers inserted",      self.counts['qualifiers_inserted'])
        self.run.info("CDS rows inserted",        len(self.cds_specific))
        self.run.info("Gene calls matched",       self.counts['gcid_matched'])
        self.run.info("Gene calls unmatched",     self.counts['gcid_unmatched'])
        self.run.info("CDS without parent links", self.counts['cds_without_parents'])
        self.run.info("Duplicates skipped",       self.counts['duplicates_skipped'])
        self.run.info("Malformed origin-crossings skipped (linear)", self.counts['skipped_malformed_origin_linear'])
        self.run.info("Multi-segment features with undefined strand skipped", self.counts['skipped_multiseg_undefined_strand'])
        self.run.info("Zero-extent / between-position features skipped",     self.counts['skipped_no_extent'])
        self.run.info("Invalid codon_start values nulled",                   self.counts['invalid_codon_start_nulled'])

        if self.new_non_builtin_types:
            self.run.info("New non-builtin feature types registered", ", ".join(sorted(self.new_non_builtin_types)))

        if self.types_seen:
            self.run.warning(None, header="By feature type", lc='green')
            for ftype in sorted(self.types_seen):
                self.run.info(ftype, self.types_seen[ftype])
