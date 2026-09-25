Get tRNA isoacceptor frequency statistics from genomes, genes, and functions.

This program converts codon usage in genes or functions into tRNA isoacceptor (anticodon-level) decoding demand, accounting for wobble pairing and, optionally, wobble modification.


## Overview

For each isoacceptor decoded by a tRNA gene found in a genome via %(anvi-scan-trnas)s, this program sums the codon frequencies of the codons that isoacceptor can decode, weighted by decoding efficiency (the same bacterial wobble weights of Sabi and Tuller, 2014, used by %(anvi-compute-trnaseq-functional-affinity)s). If a %(trna-modification-enzyme-list)s is provided and the corresponding tRNA modifying enzyme's gene is found among the genome's annotated functions, the isoacceptor's anticodon is treated as wobble-modified (expanding, in most cases, the set of codons it can decode) before this weighting is applied. For example, tadA and TilS are major tRNA modifying enzymes to consider.

Unlike %(anvi-compute-trnaseq-functional-affinity)s, which relates tRNA-seq seed abundances to codon usage, this program works from a %(contigs-db)s alone. It infers wobble modification state (A34 to inosine, C34 to lysidine at tRNA-Ile2) from whether the gene encoding the modifying enzyme is present in the genome, rather than from sequencing data.

## Basic commands

### Gene frequencies, no tRNA modifying enzyme inference

{{ codestart }}
anvi-get-trna-isoacceptor-frequencies -c %(contigs-db)s \
                                      -o path/to/output.txt
{{ codestop }}

Without `--trna-modification-enzyme-list`, every isoacceptor is treated as canonically unmodified.

### Gene frequencies with tRNA modifying enzyme inference

{{ codestart }}
anvi-get-trna-isoacceptor-frequencies -c %(contigs-db)s \
                                      --trna-modification-enzyme-list path/to/enzyme_list.txt \
                                      -o path/to/output.txt
{{ codestop }}

### Function frequencies

{{ codestart }}
anvi-get-trna-isoacceptor-frequencies -c %(contigs-db)s \
                                      --function-sources \
                                      --function-table-output path/to/function_output.txt
{{ codestop }}

### Multiple genomes

{{ codestart }}
anvi-get-trna-isoacceptor-frequencies -e %(external-genomes)s \
                                      --trna-modification-enzyme-list path/to/enzyme_list.txt \
                                      -o path/to/output.txt
{{ codestop }}

## Frequency statistics

| Get | Options |
| --- | ------- |
| Absolute weighted decoding demand (default) | |
| Relative composition (proportions of each row's total) | `--relative` |
| Summed demand across genes, one row per genome | `--sum` |
| Averaged demand across genes, one row per genome | `--average` |
| Each function's/gene's demand as a fold-change against its genome's average | `--relative --deviation-from-genome-average` |

`--deviation-from-genome-average` requires `--relative`: absolute weighted decoding demand scales with how many codons a gene or function has, so without `--relative`, a longer gene would appear to deviate from the genome average on every isoacceptor purely because of its length, not because its isoacceptor usage is actually biased. It also cannot be combined with `--sum`/`--average`, since collapsing to one row per genome and then comparing that row to the genome average is a trivial, always-1.0 comparison.

No pseudocount is added before this division. An isoacceptor with zero demand in a gene but a nonzero genome average simply divides to `0.0`. If an isoacceptor's value and its genome's average are *both* exactly zero (0/0), the result is `NaN` unless `--infinity-to-zero` is also given, in which case it becomes `0.0` -- but an isoacceptor that is genuinely absent from a genome (no tRNA gene for it at all) is always left as `NaN`, regardless of `--infinity-to-zero`.

## Inputs

The %(contigs-db)s must have been processed by %(anvi-scan-trnas)s (or %(anvi-run-hmms)s with `--also-scan-trnas`) so that it has tRNA gene calls to identify which isoacceptors exist in the genome. To use `--trna-modification-enzyme-list`, the genome must also have been annotated with gene functions (e.g., by `anvi-run-ncbi-cogs`) using the same source specified in the enzyme list.

## Decoding weights

By default, decoding efficiencies are the bacterial mean wobble values of Sabi and Tuller (2014). Use `--get-default-decoding-weights` to write this table as a template, edit it, and pass it back in with `--decoding-weights-txt`.
