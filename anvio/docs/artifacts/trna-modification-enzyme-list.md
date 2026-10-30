A tab-delimited text file describing known tRNA modification enzymes, the genomic functions that encode them, and the tRNA positions they target.

This file is provided by the user as input to %(anvi-export-genome-specific-trna-modifications)s. It defines which modification enzymes to search for in the genomes and which tRNA positions and isoacceptors to associate with each enzyme.

## Required columns

| Column | Description |
| --- | --- |
| `modifying_enzyme_name` | Short name of the modification enzyme (e.g., `cmoA`, `mnmA`) |
| `modification` | Modification produced (e.g., `cmo5U`, `s2U`, `m1G`) |
| `canonical_position` | Canonical tRNA position targeted by the enzyme (integer, e.g., `34`, `37`) |
| `expected_reference` | Nucleotide in the unmodified tRNA at this position (e.g., `T`, `A`, `G`) |
| `isoacceptor_specificity` | `Specific` if the enzyme targets particular isoacceptors; `Non-specific` otherwise |
| `aa` | Three-letter amino acid abbreviation for the targeted isoacceptor (e.g., `Leu`); `N/A` for non-specific enzymes |
| `anticodon` | Anticodon sequence of the targeted isoacceptor (e.g., `TAG`); `N/A` for non-specific enzymes |
| `function_accession` | Accession of the genomic function (e.g., a COG accession) |
| `function_source` | Source of the function annotation (e.g., `COG24_FUNCTION`) |

## Example

| modifying_enzyme_name | modification | canonical_position | expected_reference | isoacceptor_specificity | aa | anticodon | function_accession | function_source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cmoA | cmo5U | 34 | T | Specific | Leu | TAG | COG2226 | COG24_FUNCTION |
| mnmA | s2U | 34 | T | Specific | Gln | TTG | COG0482 | COG24_FUNCTION |
| miaA | i6A | 37 | A | Specific | Phe | GAA | COG0324 | COG24_FUNCTION |
| trm61 | m1A | 57 | A | Non-specific | N/A | N/A | COG2519 | COG24_FUNCTION |
