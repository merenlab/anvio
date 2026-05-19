This program imports sequence features from a %(genbank-file)s into the `contigs_sequence_features` table (and its companions: `feature_types`, `feature_relationships`, `feature_qualifiers`, and `CDS_features`) of a %(contigs-db)s. Available starting from contigs-database version 25.

Unlike %(anvi-script-process-genbank)s, which converts a GenBank file into a FASTA plus external gene calls that you then feed into %(anvi-gen-contigs-database)s, this program is **additive** — it leaves the existing `genes_in_contigs` table (and every gene-callers-id-dependent code path) untouched. Use it when you already have a contigs database and want to attach the rich annotation hierarchy from a GenBank file — including non-`gene` feature types like `mRNA`, `exon`, `intron`, and arbitrary user-defined types like `regulatory` — without rebuilding the database from scratch.

## What gets imported

Every `SeqFeature` in the GenBank file becomes one row in `contigs_sequence_features` (or N rows if the location is multi-segment, e.g. `join(...)` or `complement(join(...))`). Coordinates are stored as zero-based half-open `[start, stop)`. For multi-segment features, the row with `segment_order = 0` is the canonical row (the 5' end in transcription direction — leftmost on forward strand, rightmost on reverse strand), and every segment of the group shares a `feature_group_id` equal to the canonical row's `feature_id`. Origin-crossing features on circular contigs are stored as two segments with the join-order honored.

Each feature's GenBank qualifiers land in `feature_qualifiers` (one row per `(feature_id, key, position)`), with the exceptions of `translation`, `codon_start`, and `transl_table` on CDS features — those are promoted to dedicated columns in `CDS_features`, stored only on the canonical row of a multi-segment CDS.

Single-segment `gene` features are reconciled against the existing `genes_in_contigs` table by exact coordinate+direction match; on success, the matched `gene_callers_id` is stored on the row. Unmatched gene rows are imported with `gene_callers_id = NULL` and a warning. Multi-segment gene features are never reconciled.

Parent relationships are resolved per contig using a three-rule precedence — `locus_tag` exact match → `gene` qualifier exact match → coordinate containment alone — and link `CDS → mRNA` or `CDS → gene` (`part_of` / `derives_from` respectively), every transcript-source RNA type (`mRNA`, `ncRNA`, `tRNA`, `rRNA`, `primary_transcript`, `precursor_RNA`, `misc_RNA`, `tmRNA`, `snRNA`, `snoRNA`, `scRNA`, `antisense_RNA`) `→ gene` with `part_of`, and `exon` / `intron` `→` any of those transcript-source types with `part_of`. For multi-segment children every segment has its own row in `feature_relationships`, but the parent side always references the canonical row of the parent group.

## Synthesized `transcript` and `exon` hierarchy

GenBank files are wildly inconsistent in how they represent transcript structure: bacterial files typically have only `gene` + `CDS`; eukaryotic NCBI RefSeq files have `gene` + `mRNA(join)` + `CDS(join)`, sometimes with explicit `exon` features; non-coding genes use `ncRNA`, `tRNA`, `rRNA`, or `precursor_RNA`; pseudogenes may have `gene` + `mRNA(join)` with no `CDS`. To make canonical "find all transcripts and their exons" queries possible without file-format-specific logic, this program produces a uniform `gene → transcript → exon` hierarchy for every gene — synthesizing missing rows where needed and using literal rows in place wherever they already exist.

The literal rows and the synthesized rows share the same two tables and feature types. Two columns on `contigs_sequence_features` distinguish them:

- `derivation` is NULL for literal features (those that appeared in the GenBank file with that feature type). For synthesized features, it holds the source feature's type — e.g. `'mRNA'`, `'CDS'`, `'gene'`, `'ncRNA'`, `'transcript'`. The column answers both "is this synthesized?" (NULL vs. not) and "what was it synthesized from?".
- `derived_from_feature_id` is NULL for literal features and references the source feature's canonical `feature_id` for synthesized rows.

### Priority cascade

For each gene, transcripts are populated using this cascade:

1. **At least one transcript-source literal child of any type other than `transcript`** (`mRNA`, `ncRNA`, `tRNA`, `rRNA`, `primary_transcript`, `precursor_RNA`, `misc_RNA`, `tmRNA`, `snRNA`, `snoRNA`, `scRNA`, or `antisense_RNA`): synthesize **one transcript per source child**, with `derivation` set to that source's feature type. Multi-isoform genes get one synthesized transcript per isoform.
2. **At least one literal `transcript` child** (rare in NCBI files; some non-NCBI annotations include them): the literal transcript is already canonical — it is **promoted in place** rather than shadowed by a synthesized row. The literal keeps `derivation IS NULL` and retains the `part_of gene` relationship the parent-resolution pass gave it.
3. **Only CDS children** (no transcript children of any kind): synthesize **one transcript** with the gene's coordinates (not the CDS's — gene coords capture annotated UTRs), `derivation='CDS'`.
4. **No children at all**: synthesize **one transcript** covering the gene, `derivation='gene'`.

Steps 1 and 2 can fire together for the same gene if both kinds of children are present.

For each transcript (whether synthesized or a promoted literal), exons are populated:

- If literal `exon` features are linked to that transcript, no exon synthesis happens — literal exons fulfill the canonical hierarchy.
- Otherwise, synthesized exons mirror the segments of the transcript's source feature. For a synthesized transcript that's the literal mRNA / ncRNA / etc.; for a promoted literal transcript that's the literal transcript itself; for a CDS-derived transcript that's the CDS; for a gene-derived transcript that's the gene.

So the overall priority is: **literal exons > exons synthesized from the transcript's source segments**.

### Querying the canonical hierarchy

All transcripts and exons are queryable through `feature_type` alone — literal and synthesized rows live in the same table and are equivalent for canonical-hierarchy queries:

```sql
-- Every gene's transcript count and exon count, regardless of input style
SELECT g.feature_id            AS gene_id,
       g.external_id           AS gene_locus_tag,
       COUNT(DISTINCT t.feature_id) AS n_transcripts,
       COUNT(DISTINCT e.feature_id) AS n_exons
FROM contigs_sequence_features g
LEFT JOIN feature_relationships rt ON rt.parent_feature_id = g.feature_id AND rt.relationship = 'part_of'
LEFT JOIN contigs_sequence_features t ON t.feature_id = rt.child_feature_id AND t.feature_type = 'transcript'
LEFT JOIN feature_relationships re ON re.parent_feature_id = t.feature_id AND re.relationship = 'part_of'
LEFT JOIN contigs_sequence_features e ON e.feature_id = re.child_feature_id AND e.feature_type = 'exon'
WHERE g.feature_type = 'gene'
GROUP BY g.feature_id;
```

Filter on `derivation` only when you specifically want to separate literal-from-file rows from importer-synthesized rows:

- `WHERE feature_type = 'transcript'` — all transcripts (literal-promoted + synthesized).
- `WHERE feature_type = 'transcript' AND derivation IS NULL` — only literal transcripts that came from the input file (rare).
- `WHERE feature_type = 'transcript' AND derivation IS NOT NULL` — only synthesized transcripts.
- `WHERE feature_type = 'exon' AND derivation IS NULL` — only literal exon features.
- `WHERE feature_type = 'exon' AND derivation = 'mRNA'` — only synthesized exons derived from mRNAs.

### Biological caveats

- **Transcripts with `derivation = 'CDS'`** reflect the gene boundary, not the true transcript boundary. The UTRs are not annotated in the source CDS, so the synthesized transcript inherits the gene's `[start, stop)`. In bacteria gene and CDS coords are typically identical, so this is fine; in eukaryotic CDS-only annotations it means UTRs are unknown.
- **Exons with `derivation = 'CDS'`** cover only the coding portion of the true exons. UTRs are absent. The `derivation` column tells consumers exactly this.
- **Multi-isoform genes** produce per-isoform synthesized exon sets without deduplication. Two isoforms sharing an exon coordinate produce two distinct exon rows (different `derived_from_feature_id`, different `external_id`). For "unique exonic regions of gene X" the user must `SELECT DISTINCT start, stop, direction`; deduplication is not done at synthesis time because it would erase isoform identity.
- **Mixed transcript types under one gene** (e.g. an `ncRNA` and a `CDS` under the same gene, as in some viral genomes): the synthesized transcript is derived from the `ncRNA`; the CDS is NOT linked to it, because the parent-relationship rules link CDSes to mRNAs or genes, not to ncRNAs. Walking "CDSes under synthesized transcripts" misses such CDSes — for comprehensive coverage walk via the literal `gene` row's children.

### When synthesis cannot help

Synthesis depends on the parent-relationship resolution pass to find a gene's children. If that pass fails to link a CDS to its mRNA because the input GenBank file lacks consistent `locus_tag` qualifiers, the CDS will not be re-parented to the synthesized transcript and will remain only `derives_from` the gene. The program emits warnings about unlinked CDSes; if you see them, examine the input file for malformed or inconsistent qualifiers.

## Example: importing features from a GenBank file

After running %(anvi-script-process-genbank)s + %(anvi-gen-contigs-database)s (or starting from an existing %(contigs-db)s that was built some other way), you can add GenBank features like this:

{{ codestart }}
anvi-import-genbank-features -c %(contigs-db)s \
                             -i %(genbank-file)s
{{ codestop }}

Behind the scenes the program parses the GenBank file twice: once to buffer all features and validate LOCUS names against the contigs DB, and a second time to resolve parent relationships and write everything atomically. A failed import (e.g. a duplicate feature_id) rolls back the entire transaction; the database is either fully updated or not touched at all.

## --source-name and multi-source coexistence

Every row imported by this program carries a `source` value (default: `genbank_import`). You can override it with `--source-name`:

{{ codestart }}
anvi-import-genbank-features -c %(contigs-db)s \
                             -i %(genbank-file)s \
                             --source-name 'refseq_v100'
{{ codestop }}

`--source-name` must match the regex `^[A-Za-z0-9_-]+$` (this is enforced because the value is hashed into every `feature_id` with TAB as a separator; allowing arbitrary characters would risk separator collisions).

Multiple `--source-name` values can coexist in the same database — running the program a second time with a different source name is additive and does not require `--force`. This lets you keep, e.g., a stable RefSeq import side-by-side with an experimental annotation pipeline's output.

## --force semantics

Running the program with the same `--source-name` against a database that already has rows from that source will fail with a clear error message — anvi'o refuses to silently merge or overwrite. To replace the existing rows of one source while leaving rows from other sources untouched, add `--force`:

{{ codestart }}
anvi-import-genbank-features -c %(contigs-db)s \
                             -i %(genbank-file)s \
                             --source-name 'refseq_v100' \
                             --force
{{ codestop }}

Because `feature_id` generation is deterministic (a 16-character SHA-224 truncation of the TAB-joined input fields), replaying the same GenBank file with `--force` produces the exact same set of `feature_id`s — so external systems that have cached `feature_id`s for cross-referencing remain consistent across runs.

## Related programs and artifacts

- %(anvi-script-process-genbank)s — the older path that converts a GenBank into FASTA + external gene calls for %(anvi-gen-contigs-database)s. The two programs are complementary: `anvi-script-process-genbank` populates `genes_in_contigs` from the GenBank's gene calls; `anvi-import-genbank-features` populates the new sequence-features tables (and reconciles back to the gene calls when possible).
- %(anvi-gen-contigs-database)s — generates the contigs database the new tables live in.
- %(anvi-import-functions)s — imports gene-level functional annotations into the existing `gene_functions` table; orthogonal to and unaffected by `anvi-import-genbank-features`.
- %(contigs-db)s — the artifact this program writes into.
