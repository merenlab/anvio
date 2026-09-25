A tab-delimited text file describing known tRNA modification enzymes, the genomic functions that encode them, and the tRNA positions they target.

This file is provided by the user as input to %(anvi-get-trna-isoacceptor-frequencies)s. It defines which modification enzymes to search for in a genome and which tRNA positions and isoacceptors to associate with each enzyme.

Only enzymes that modify the anticodon wobble position (position 34) into inosine or lysidine currently affect the output of %(anvi-get-trna-isoacceptor-frequencies)s, since these are the only two wobble modifications its decoding-efficiency model understands. Rows for other enzymes and positions (e.g., a wobble uridine or a position-37 modification) are valid entries in this file, but are ignored by that program with a warning.

## Required columns

| Column | Description |
| --- | --- |
| `modifying_enzyme_name` | Short name of the modification enzyme (e.g., `tadA`, `tilS`) |
| `modification` | Modification produced. Use `inosine` (or `I`) for an A34-to-inosine wobble modification, and `lysidine` (or `L`, `k2C`) for a C34-to-lysidine wobble modification at tRNA-Ile2; any other value is accepted but not actionable by %(anvi-get-trna-isoacceptor-frequencies)s |
| `canonical_position` | Canonical tRNA position targeted by the enzyme (integer, e.g., `34`, `37`) |
| `expected_reference` | Nucleotide in the unmodified tRNA at this position (e.g., `A`, `C`) |
| `isoacceptor_specificity` | `Specific` if the enzyme targets a single isoacceptor; `Non-specific` otherwise |
| `aa` | Three-letter amino acid abbreviation (or isotype label, e.g. `Ile2`) for the targeted isoacceptor; `N/A` for non-specific enzymes |
| `anticodon` | Anticodon sequence of the targeted isoacceptor (e.g., `ACG`); `N/A` for non-specific enzymes |
| `function_accession` | Accession of the genomic function that encodes the enzyme (e.g., a COG accession) |
| `function_source` | Source of the function annotation (e.g., `COG20_FUNCTION`) |

## Example

| modifying_enzyme_name | modification | canonical_position | expected_reference | isoacceptor_specificity | aa | anticodon | function_accession | function_source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tadA | inosine | 34 | A | Specific | Arg | ACG | COG0590 | COG20_FUNCTION |
| tilS | lysidine | 34 | C | Specific | Ile2 | CAT | COG0037 | COG20_FUNCTION |
| miaA | i6A | 37 | A | Specific | Phe | GAA | COG0324 | COG20_FUNCTION |
