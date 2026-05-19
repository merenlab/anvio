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

# Feature types that represent transcript-source RNAs — used by step 13 to link
# them all to gene parents (the v25 implementation only linked mRNA), and by
# step 13.5 (synthesis) to find a gene's transcript-source children. The set
# excludes literal `transcript` features in step 13.5a's case 1 lookups but
# IS included here so that literal transcript features, if any, still receive
# the same treatment as the other RNA types. New transcript-source feature
# types should be added here.
TRANSCRIPT_SOURCE_TYPES = {
    'transcript', 'mRNA', 'primary_transcript', 'ncRNA', 'tRNA', 'rRNA',
    'precursor_RNA', 'misc_RNA', 'tmRNA', 'snRNA', 'snoRNA', 'scRNA', 'antisense_RNA',
}

# Transcript-source types that step 13 sub-rule 4.1 links to gene parents. We
# exclude literal `transcript` because synthesis treats it as a special case
# (shadow-transcript pattern in step 13.5a) — see Section 4 of the prompt.
TRANSCRIPT_SOURCE_TYPES_FOR_GENE_PARENTING = TRANSCRIPT_SOURCE_TYPES - {'transcript'}

# Feature-type pairs (child_type → ordered list of (candidate parent types, relationship)
# tuples). Each list entry is one precedence level; candidates from all listed parent
# types at the same level are pooled and considered together. The first precedence
# level that produces matches wins. Sub-rule 4.1 (every transcript-source type except
# `transcript` parents to `gene`) is appended after the static rules.
PARENT_RULES = {
    'CDS':    [(['mRNA'], 'part_of'), (['gene'], 'derives_from')],
    # Sub-rule 4.2: exon and intron features parent to any transcript-source feature
    # (mRNA, ncRNA, tRNA, ..., and the rare literal `transcript`). The v25 rule only
    # considered `mRNA`; this generalization lets an exon under an ncRNA, tRNA, etc.
    # be linked correctly to that RNA.
    'exon':   [(sorted(TRANSCRIPT_SOURCE_TYPES), 'part_of')],
    'intron': [(sorted(TRANSCRIPT_SOURCE_TYPES), 'part_of')],
}
# Sub-rule 4.1: every transcript-source type except literal `transcript` parents to
# `gene` with `part_of`. The literal `transcript` case is handled by the synthesis
# layer (step 13.5) via the shadow-transcript pattern.
for _ttype in TRANSCRIPT_SOURCE_TYPES_FOR_GENE_PARENTING:
    PARENT_RULES[_ttype] = [(['gene'], 'part_of')]
del _ttype


HASH_INPUT_FIELDS = ('contig', 'feature_type', 'source', 'start', 'stop', 'direction',
                     'external_id', 'derivation', 'derived_from_feature_id')


def compute_feature_id(contig, feature_type, source, start, stop, direction,
                       external_id=None, derivation=None, derived_from_feature_id=None):
    """Compute the 16-character hex `feature_id` for one row.

    The hash input is the TAB-joined tuple of nine fields. The TAB is mandatory —
    concatenating fields without a separator would let adjacent coordinates blur
    into hash collisions (start=12,stop=345 vs start=123,stop=45 both produce
    '12345'). Before hashing we assert no field contains TAB itself; silent hash
    divergence would be a far worse failure mode than an explicit error.

    The nine fields are: contig, feature_type, source, start, stop, direction,
    external_id, derivation, derived_from_feature_id. `external_id` disambiguates
    alternative-splicing isoforms whose shared exon coordinates would otherwise
    collide under the original 6-field convention. `derivation` and
    `derived_from_feature_id` are NULL for literal features (read from the
    GenBank file) and non-NULL for synthesized rows produced by step 13.5. NULL
    values are encoded as the empty string in the hash input, and the trailing
    empty TABs are appended unconditionally so the hash input string is
    well-formed for every combination of NULL/non-NULL field values.
    """

    direction_str               = direction               if direction               is not None else ''
    external_id_str             = external_id             if external_id             is not None else ''
    derivation_str              = derivation              if derivation              is not None else ''
    derived_from_feature_id_str = derived_from_feature_id if derived_from_feature_id is not None else ''
    fields = [contig, feature_type, source, str(start), str(stop), direction_str,
              external_id_str, derivation_str, derived_from_feature_id_str]
    for field in fields:
        if '\t' in field:
            raise ConfigError(f"compute_feature_id refuses to hash an input field that contains a TAB character: '{field!r}'. "
                              f"This would create silent hash divergence. Please report this as an anvi'o bug if a TAB ever appears in a contig name, "
                              f"feature type, source name, external_id (locus_tag), derivation, or derived_from_feature_id.")
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

        # step 13.5 synthesis side-tables. Populated by `_synthesize_pass2` and merged
        # into `self.features` / `self.relationships` once the per-gene loop and the
        # uniqueness check are done. Kept separate during synthesis so the duplicate
        # check operates only on the synthesized set.
        self.synthesized_features = []
        self.synthesized_relationships = []
        self.synth_counts = {
            'transcripts_by_derivation': {},
            'exons_by_derivation': {},
            'genes_without_synthesized_transcript': 0,
        }

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
            self.progress.update("Synthesizing transcript and exon hierarchy (step 13.5)...")
            self._synthesize_pass2()
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
                # literal features always have derivation=None / derived_from=None;
                # external_id comes from locus_tag (None when the source feature lacks one)
                fid = compute_feature_id(contig, ftype, self.source_name, seg_start, seg_stop, direction,
                                         external_id=locus_tag)
                row = {
                    'feature_id':              fid,
                    'contig':                  contig,
                    'feature_type':            ftype,
                    'source':                  self.source_name,
                    'start':                   seg_start,
                    'stop':                    seg_stop,
                    'direction':               direction,
                    'partial_fiveprime':       partial_5p,
                    'partial_threeprime':      partial_3p,
                    'feature_group_id':        None,
                    'segment_order':           None,
                    'external_id':             locus_tag,
                    'gene_callers_id':         None,
                    'derivation':              None,
                    'derived_from_feature_id': None,
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
        on different hash-input tuples are astronomically unlikely with 16-char SHA-224
        truncation at first-PR scales — if one occurs we raise rather than silently merge data.

        The comparison covers every field that feeds the 9-field hash convention so
        that a literal feature and another feature with the same coordinates but a
        distinct `external_id` (or `derivation`, or `derived_from_feature_id`) is
        correctly recognized as a distinct feature rather than misclassified as a
        genuine hash collision.
        """

        fid = row['feature_id']
        for existing in self.features:
            if existing['feature_id'] != fid:
                continue
            same = all(existing[k] == row[k] for k in ('contig', 'feature_type', 'source', 'start', 'stop', 'direction',
                                                       'external_id', 'derivation', 'derived_from_feature_id'))
            if same:
                self.run.warning(f"Duplicate feature in GenBank input: contig={row['contig']}, type={row['feature_type']}, "
                                 f"[{row['start']}, {row['stop']}), direction={row['direction']}, external_id={row['external_id']}. "
                                 f"Dropping the duplicate.")
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

            for parent_types, label in PARENT_RULES[ftype]:
                # candidates from every listed parent type at this precedence level are
                # pooled — important for exon/intron which considers all transcript-source
                # types together rather than cascading through them in order.
                candidates = []
                for pt in parent_types:
                    candidates.extend(by_type.get(pt, []))
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
                # flatten all candidate parent types into a unique, sorted list for the warning
                parent_type_names = sorted({pt for parent_types, _ in PARENT_RULES[ftype] for pt in parent_types})
                self.run.warning(f"{ftype} '{cfid}' on contig '{info['contig']}' has no {' or '.join(parent_type_names)} parent.")
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
    # Pass 2 step 13.5: synthesis of transcript / exon hierarchy
    # -----------------------------------------------------------------

    def _synthesize_pass2(self):
        """Synthesize `transcript` and `exon` rows per gene so the database exposes a
        uniform gene → transcript → exon hierarchy regardless of how the source GenBank
        file structured its annotations.

        Runs once per gene canonical feature. The three cases (transcript-source children
        present, only CDS children, no children) are dispatched by `_synthesize_for_gene`.
        Synthesis writes into `self.synthesized_features` and `self.synthesized_relationships`
        first; after a defensive uniqueness check (step 13.5e) the synthesized rows are
        merged into the main lists for `populate_features` to consume in step 14.
        """

        # Index relationships by parent_feature_id so per-gene synthesis can find children
        # in O(1) per lookup rather than scanning self.relationships for each gene.
        children_by_parent = {}
        for rel in self.relationships:
            children_by_parent.setdefault(rel['parent_feature_id'], []).append(rel)

        # fid → row lookup for the literal features. Used throughout synthesis for
        # coordinate / direction / external_id reads.
        features_by_fid = {row['feature_id']: row for row in self.features}

        # Process each gene canonical feature. Insertion order in self.groups is
        # deterministic (Python 3.7+ dict ordering) so synthesized rows are emitted in
        # a stable order across runs of the same input file.
        gene_canonical_fids = [fid for fid, info in self.groups.items() if info['feature_type'] == 'gene']
        for gene_canonical_fid in gene_canonical_fids:
            produced = self._synthesize_for_gene(gene_canonical_fid, children_by_parent, features_by_fid)
            if not produced:
                self.synth_counts['genes_without_synthesized_transcript'] += 1
                gene_info = self.groups[gene_canonical_fid]
                self.run.warning(f"Gene '{gene_canonical_fid}' on contig '{gene_info['contig']}' produced no synthesized "
                                 f"transcript. Every gene should receive at least one — this almost certainly indicates "
                                 f"an anvi'o bug; please report.")

        # Step 13.5e: defensive uniqueness check on synthesized feature_ids. Synthesized
        # rows are designed to hash uniquely by construction; this guard catches bugs
        # and unexpected edge cases rather than relying on the underlying UNIQUE INDEX
        # to surface a duplicate at insert time.
        seen = set()
        for row in self.synthesized_features:
            if row['feature_id'] in seen:
                self.progress.end()
                raise ConfigError(f"Step 13.5e detected two synthesized features with the same feature_id "
                                  f"'{row['feature_id']}'. Synthesized features are designed to hash uniquely "
                                  f"by construction; encountering a collision indicates a bug in the synthesis "
                                  f"layer or extremely pathological input. Refusing to silently merge.")
            seen.add(row['feature_id'])

        # Merge synthesized rows into the main lists. From here on populate_features
        # treats them as ordinary feature / relationship rows tagged by `derivation`.
        self.features.extend(self.synthesized_features)
        self.relationships.extend(self.synthesized_relationships)


    def _synthesize_for_gene(self, gene_canonical_fid, children_by_parent, features_by_fid):
        """Synthesize transcript(s) and exon(s) for one gene canonical row.

        Dispatches to one of three cases per step 13.5a:
            1. Transcript-source literal children present → one transcript per source.
            2. CDS literal children only → one transcript with the gene's coordinates.
            3. No children at all → one transcript covering the gene.

        Returns True if at least one transcript was synthesized (which should always
        be the case; the explicit return lets the caller surface a warning if not).
        """

        # Collect this gene's canonical children, grouped by feature_type. The gene's
        # canonical fid is always the parent_feature_id in feature_relationships even
        # when the gene is multi-segment (canonical-parent convention).
        canonical_children_by_type = self._gather_canonical_children(gene_canonical_fid, children_by_parent)

        # Case 1: at least one transcript-source RNA among the gene's children.
        # Order: by feature_type alphabetically, then by source min_start within type.
        # This determinism keeps synthesized-row emission order stable across runs.
        source_canonicals_in_order = []
        for ftype in sorted(TRANSCRIPT_SOURCE_TYPES):
            for fid in sorted(canonical_children_by_type.get(ftype, []),
                              key=lambda c: self.groups[c]['min_start']):
                source_canonicals_in_order.append(fid)

        if source_canonicals_in_order:
            for src_canonical in source_canonicals_in_order:
                self._synthesize_case1(gene_canonical_fid, src_canonical, children_by_parent, features_by_fid)
            return True

        # Case 2: no transcript-source children, but CDS children present.
        cds_canonicals = canonical_children_by_type.get('CDS', [])
        if cds_canonicals:
            # Tie-breaker for the unusual multi-CDS-group-under-one-gene case: pick the
            # CDS group whose spanning extent has the smallest min_start.
            best_cds = min(cds_canonicals, key=lambda c: self.groups[c]['min_start'])
            self._synthesize_case2(gene_canonical_fid, best_cds, features_by_fid)
            return True

        # Case 3: lone gene — no transcript-source, no CDS children.
        self._synthesize_case3(gene_canonical_fid, features_by_fid)
        return True


    def _gather_canonical_children(self, parent_canonical_fid, children_by_parent):
        """Return `{feature_type: [canonical_fid, ...]}` for the children of `parent_canonical_fid`.

        Multi-segment children appear once (by their canonical fid). Children whose
        canonical fid is not registered in `self.feature_meta` are skipped (defensive
        — shouldn't happen in practice since step 13 only adds rows for known features).
        """

        result = {}
        seen_canonical = set()
        for rel in children_by_parent.get(parent_canonical_fid, []):
            child_fid = rel['child_feature_id']
            meta = self.feature_meta.get(child_fid)
            if not meta:
                continue
            child_canonical = meta['canonical_fid']
            if child_canonical in seen_canonical:
                continue
            seen_canonical.add(child_canonical)
            result.setdefault(meta['feature_type'], []).append(child_canonical)
        return result


    def _synthesize_case1(self, gene_canonical_fid, src_canonical_fid, children_by_parent, features_by_fid):
        """Case 1 of step 13.5a: synthesize one transcript from a transcript-source literal."""

        src_info = self.groups[src_canonical_fid]
        src_canonical_row = features_by_fid[src_canonical_fid]
        derivation = src_info['feature_type']
        external_id = src_info['locus_tag']
        contig = src_info['contig']
        direction = src_canonical_row['direction']
        src_segments = src_info['segments']

        # Origin-crossing sources synthesize a multi-segment transcript with 1:1 segment
        # correspondence. Non-origin-crossing sources synthesize a single-segment transcript
        # spanning the whole source extent; the per-segment structure goes to exons in 13.5b(ii).
        if self._is_origin_crossing_for_canonical(src_canonical_fid, features_by_fid):
            transcript_coords = [(features_by_fid[s]['start'], features_by_fid[s]['stop']) for s in src_segments]
        else:
            transcript_coords = [(src_info['min_start'], src_info['max_stop'])]

        # Inherit the source's 5' and 3' partial flags. The 5' flag lives on the canonical
        # (segment_order=0) row; the 3' flag lives on the last segment in transcription order.
        partial_5p = src_canonical_row['partial_fiveprime']
        partial_3p = features_by_fid[src_segments[-1]]['partial_threeprime']

        transcript_canonical_fid = self._add_transcript_rows(
            gene_canonical_fid, contig, direction, transcript_coords,
            external_id, derivation, src_canonical_fid, partial_5p, partial_3p,
        )

        # Step 13.5b(i): re-parent every literal child of the source to the synthesized
        # transcript with `part_of`. The original relationship rows are not modified —
        # both coexist.
        for rel in children_by_parent.get(src_canonical_fid, []):
            self.synthesized_relationships.append({
                'child_feature_id':  rel['child_feature_id'],
                'parent_feature_id': transcript_canonical_fid,
                'relationship':      'part_of',
            })

        # Step 13.5b(ii): synthesize exons unless any literal exon is reachable from the
        # synthesized transcript (which is equivalent to: any literal exon child of the source).
        if not self._has_literal_exon_among_children(src_canonical_fid, children_by_parent):
            for src_seg_fid in src_segments:
                src_seg_row = features_by_fid[src_seg_fid]
                self._add_exon_row(transcript_canonical_fid, contig, direction,
                                   src_seg_row['start'], src_seg_row['stop'],
                                   src_seg_row['external_id'], derivation, src_seg_fid)


    def _synthesize_case2(self, gene_canonical_fid, cds_canonical_fid, features_by_fid):
        """Case 2 of step 13.5a: synthesize ONE transcript from a CDS (gene has CDS children
        but no transcript-source children). Coordinates are the GENE's, not the CDS's, since
        gene coordinates capture the annotator's intent including UTRs (in eukaryotes).
        """

        gene_info = self.groups[gene_canonical_fid]
        gene_canonical_row = features_by_fid[gene_canonical_fid]
        cds_info = self.groups[cds_canonical_fid]
        contig = gene_info['contig']
        direction = gene_canonical_row['direction']
        external_id = gene_canonical_row['external_id']

        if self._is_origin_crossing_for_canonical(gene_canonical_fid, features_by_fid):
            gene_segments = gene_info['segments']
            transcript_coords = [(features_by_fid[s]['start'], features_by_fid[s]['stop']) for s in gene_segments]
        else:
            transcript_coords = [(gene_info['min_start'], gene_info['max_stop'])]

        partial_5p = gene_canonical_row['partial_fiveprime']
        partial_3p = features_by_fid[gene_info['segments'][-1]]['partial_threeprime']

        transcript_canonical_fid = self._add_transcript_rows(
            gene_canonical_fid, contig, direction, transcript_coords,
            external_id, 'CDS', cds_canonical_fid, partial_5p, partial_3p,
        )

        # Step 13.5b(i) special case: the CDS is the *source* (not a child of one), so
        # link every CDS segment to the synthesized transcript with `part_of`. This is
        # in addition to the CDS's step-13 `derives_from gene` relationship.
        for cds_seg_fid in cds_info['segments']:
            self.synthesized_relationships.append({
                'child_feature_id':  cds_seg_fid,
                'parent_feature_id': transcript_canonical_fid,
                'relationship':      'part_of',
            })

        # Step 13.5b(ii): one exon per CDS segment. UTRs absent — `derivation='CDS'` tells
        # consumers these exons reflect only the coding portion.
        for cds_seg_fid in cds_info['segments']:
            cds_seg_row = features_by_fid[cds_seg_fid]
            self._add_exon_row(transcript_canonical_fid, contig, direction,
                               cds_seg_row['start'], cds_seg_row['stop'],
                               cds_seg_row['external_id'], 'CDS', cds_seg_fid)


    def _synthesize_case3(self, gene_canonical_fid, features_by_fid):
        """Case 3 of step 13.5a: lone gene — no transcript-source, no CDS children.

        Synthesize one transcript covering the gene and one exon per gene segment
        (almost always one, since multi-segment genes are exceedingly rare).
        """

        gene_info = self.groups[gene_canonical_fid]
        gene_canonical_row = features_by_fid[gene_canonical_fid]
        contig = gene_info['contig']
        direction = gene_canonical_row['direction']
        external_id = gene_canonical_row['external_id']
        gene_segments = gene_info['segments']

        if self._is_origin_crossing_for_canonical(gene_canonical_fid, features_by_fid):
            transcript_coords = [(features_by_fid[s]['start'], features_by_fid[s]['stop']) for s in gene_segments]
        else:
            transcript_coords = [(gene_info['min_start'], gene_info['max_stop'])]

        partial_5p = gene_canonical_row['partial_fiveprime']
        partial_3p = features_by_fid[gene_segments[-1]]['partial_threeprime']

        transcript_canonical_fid = self._add_transcript_rows(
            gene_canonical_fid, contig, direction, transcript_coords,
            external_id, 'gene', gene_canonical_fid, partial_5p, partial_3p,
        )

        # Step 13.5b(ii): one exon per gene segment. For the typical single-segment gene
        # this collapses to one exon spanning the whole gene.
        for seg_fid in gene_segments:
            seg_row = features_by_fid[seg_fid]
            self._add_exon_row(transcript_canonical_fid, contig, direction,
                               seg_row['start'], seg_row['stop'],
                               seg_row['external_id'], 'gene', seg_fid)


    def _add_transcript_rows(self, gene_canonical_fid, contig, direction, segment_coords,
                             external_id, derivation, derived_from_canonical_fid,
                             partial_5p, partial_3p):
        """Append synthesized transcript rows and their `part_of gene` relationships.

        `segment_coords` is a list of `(start, stop)` tuples in transcription order. A
        single-element list produces a single-segment transcript (NULL `feature_group_id`
        and `segment_order` per the v25 convention); a multi-element list produces a
        multi-segment transcript sharing the canonical row's fid as `feature_group_id`.
        Returns the canonical (segment_order=0) row's feature_id.
        """

        multi_segment = len(segment_coords) > 1
        seg_fids = []
        for start, stop in segment_coords:
            fid = compute_feature_id(contig, 'transcript', self.source_name, start, stop, direction,
                                     external_id=external_id, derivation=derivation,
                                     derived_from_feature_id=derived_from_canonical_fid)
            seg_fids.append(fid)
        canonical_fid = seg_fids[0]
        last_idx = len(segment_coords) - 1

        for idx, (start, stop) in enumerate(segment_coords):
            row = {
                'feature_id':              seg_fids[idx],
                'contig':                  contig,
                'feature_type':            'transcript',
                'source':                  self.source_name,
                'start':                   start,
                'stop':                    stop,
                'direction':               direction,
                'partial_fiveprime':       partial_5p if idx == 0        else 0,
                'partial_threeprime':      partial_3p if idx == last_idx else 0,
                'feature_group_id':        canonical_fid if multi_segment else None,
                'segment_order':           idx if multi_segment else None,
                'external_id':             external_id,
                'gene_callers_id':         None,
                'derivation':              derivation,
                'derived_from_feature_id': derived_from_canonical_fid,
            }
            self.synthesized_features.append(row)
            # Canonical-parent convention: every transcript segment becomes a child of the
            # gene's canonical fid (not of the gene's segments).
            self.synthesized_relationships.append({
                'child_feature_id':  seg_fids[idx],
                'parent_feature_id': gene_canonical_fid,
                'relationship':      'part_of',
            })

        self.synth_counts['transcripts_by_derivation'][derivation] = \
            self.synth_counts['transcripts_by_derivation'].get(derivation, 0) + 1

        return canonical_fid


    def _add_exon_row(self, transcript_canonical_fid, contig, direction,
                      start, stop, external_id, derivation, derived_from_fid):
        """Append one synthesized exon row, parented to the synthesized transcript.

        Synthesized exons are always single-segment in `contigs_sequence_features`. A
        multi-segment source produces multiple independent single-segment exon rows,
        each pointing back to its source segment via `derived_from_feature_id`.
        """

        fid = compute_feature_id(contig, 'exon', self.source_name, start, stop, direction,
                                 external_id=external_id, derivation=derivation,
                                 derived_from_feature_id=derived_from_fid)
        row = {
            'feature_id':              fid,
            'contig':                  contig,
            'feature_type':            'exon',
            'source':                  self.source_name,
            'start':                   start,
            'stop':                    stop,
            'direction':               direction,
            'partial_fiveprime':       0,
            'partial_threeprime':      0,
            'feature_group_id':        None,
            'segment_order':           None,
            'external_id':             external_id,
            'gene_callers_id':         None,
            'derivation':              derivation,
            'derived_from_feature_id': derived_from_fid,
        }
        self.synthesized_features.append(row)
        self.synthesized_relationships.append({
            'child_feature_id':  fid,
            'parent_feature_id': transcript_canonical_fid,
            'relationship':      'part_of',
        })
        self.synth_counts['exons_by_derivation'][derivation] = \
            self.synth_counts['exons_by_derivation'].get(derivation, 0) + 1


    def _has_literal_exon_among_children(self, parent_canonical_fid, children_by_parent):
        """Return True iff any literal `exon` feature is a child of `parent_canonical_fid`."""

        for rel in children_by_parent.get(parent_canonical_fid, []):
            meta = self.feature_meta.get(rel['child_feature_id'])
            if meta and meta['feature_type'] == 'exon':
                return True
        return False


    def _is_origin_crossing_for_canonical(self, canonical_fid, features_by_fid):
        """Detect origin-crossing for an existing canonical group by inspecting segment coords.

        Mirrors the predicate used by pass 1 (`_is_origin_crossing`) but operates on the
        already-parsed `self.features` rows rather than Biopython part objects, so it can
        be called during synthesis when only canonical fids are at hand.
        """

        info = self.groups[canonical_fid]
        contig = info['contig']
        if not self.contig_is_circular.get(contig, False):
            return False
        if len(info['segments']) < 2:
            return False
        contig_length = self.contig_lengths[contig]
        starts_at_zero = False
        ends_at_length = False
        for seg_fid in info['segments']:
            seg = features_by_fid[seg_fid]
            if seg['start'] == 0:
                starts_at_zero = True
            if seg['stop'] == contig_length:
                ends_at_length = True
        return starts_at_zero and ends_at_length


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

        # Step 13.5 synthesis stats. Reported by derivation so users can see how their
        # genes' transcript hierarchies were derived (literal mRNAs vs. CDS fallback vs.
        # lone-gene fallback). The genes-without-transcript count should always be zero —
        # we emit it unconditionally so any nonzero value is immediately visible.
        synth_t = self.synth_counts['transcripts_by_derivation']
        synth_e = self.synth_counts['exons_by_derivation']
        if synth_t or synth_e or self.synth_counts['genes_without_synthesized_transcript']:
            self.run.warning(None, header="Synthesis (step 13.5)", lc='green')
            self.run.info("Synthesized transcripts (total)", sum(synth_t.values()))
            if synth_t:
                for derivation in sorted(synth_t):
                    self.run.info(f"  from {derivation}", synth_t[derivation])
            self.run.info("Synthesized exons (total)",       sum(synth_e.values()))
            if synth_e:
                for derivation in sorted(synth_e):
                    self.run.info(f"  from {derivation}", synth_e[derivation])
            # This should be zero unless there is a bug; per-gene warnings during synthesis
            # already named the offending genes, so here we just surface the aggregate count.
            self.run.info("Genes with no synthesized transcript", self.synth_counts['genes_without_synthesized_transcript'])
