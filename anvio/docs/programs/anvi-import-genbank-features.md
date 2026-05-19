This program imports sequence features from a %(genbank-file)s into the `contigs_sequence_features` table (and its companions: `feature_types`, `feature_relationships`, `feature_qualifiers`, and `CDS_features`) of a %(contigs-db)s. Available starting from contigs-database version 25.

Unlike %(anvi-script-process-genbank)s, which converts a GenBank file into a FASTA plus external gene calls that you then feed into %(anvi-gen-contigs-database)s, this program is **additive** — it leaves the existing `genes_in_contigs` table (and every gene-callers-id-dependent code path) untouched. Use it when you already have a contigs database and want to attach the rich annotation hierarchy from a GenBank file — including non-`gene` feature types like `mRNA`, `exon`, `intron`, and arbitrary user-defined types like `regulatory` — without rebuilding the database from scratch.

## What gets imported

Every `SeqFeature` in the GenBank file becomes one row in `contigs_sequence_features` (or N rows if the location is multi-segment, e.g. `join(...)` or `complement(join(...))`). Coordinates are stored as zero-based half-open `[start, stop)`. For multi-segment features, the row with `segment_order = 0` is the canonical row (the 5' end in transcription direction — leftmost on forward strand, rightmost on reverse strand), and every segment of the group shares a `feature_group_id` equal to the canonical row's `feature_id`. Origin-crossing features on circular contigs are stored as two segments with the join-order honored.

Each feature's GenBank qualifiers land in `feature_qualifiers` (one row per `(feature_id, key, position)`), with the exceptions of `translation`, `codon_start`, and `transl_table` on CDS features — those are promoted to dedicated columns in `CDS_features`, stored only on the canonical row of a multi-segment CDS.

Single-segment `gene` features are reconciled against the existing `genes_in_contigs` table by exact coordinate+direction match; on success, the matched `gene_callers_id` is stored on the row. Unmatched gene rows are imported with `gene_callers_id = NULL` and a warning. Multi-segment gene features are never reconciled.

Parent relationships are resolved per contig using a three-rule precedence — `locus_tag` exact match → `gene` qualifier exact match → coordinate containment alone — and link `CDS → mRNA` or `CDS → gene` (`part_of` / `derives_from` respectively), `mRNA → gene`, `exon → mRNA`, and `intron → mRNA`. For multi-segment children every segment has its own row in `feature_relationships`, but the parent side always references the canonical row of the parent group.

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
