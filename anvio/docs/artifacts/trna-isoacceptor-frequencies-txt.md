This file contains **the tRNA isoacceptor (anticodon-level) decoding demand implied by codon usage in your %(contigs-db)s**, as computed by %(anvi-get-trna-isoacceptor-frequencies)s.

This is a tab-delimited table where each column represents an isoacceptor (formatted as `{amino_acid}_{anticodon}`, e.g. `Arg_ICG`) and each row represents a gene or a function, depending on how the program was run. A value is the sum of that gene's or function's codon frequencies, weighted by how efficiently the isoacceptor's anticodon decodes each codon (accounting for wobble pairing and, when a %(trna-modification-enzyme-list)s is provided, wobble modification state inferred from the presence of the relevant modifying enzyme).

An isoacceptor column is only present for a genome if that genome actually has a tRNA gene for it (found via %(anvi-scan-trnas)s). When multiple genomes are analyzed together, a genome lacking a given isoacceptor's tRNA gene entirely will have `NaN` (not `0`) in that column, since `0` would incorrectly imply the isoacceptor exists but is never used.

Values are absolute weighted decoding demand by default. `--relative` reports proportions of each row's total; `--sum`/`--average` collapse to one row per genome; `--relative --deviation-from-genome-average` reports each row as a fold-change against that genome's own average isoacceptor profile (requiring `--relative` so the comparison isn't just driven by gene length). No pseudocount is added before this division: a 0/0 quotient (an isoacceptor with zero demand everywhere in the genome) produces `NaN` unless `--infinity-to-zero` is also given, in which case it becomes `0.0` -- an isoacceptor genuinely absent from a genome is always left as `NaN` either way.

### Example

    genome_name  gene_caller_id  Arg_ACG  Arg_ICG  Ile2_CAT ...
        HIMB83             1            0        2         1
        HIMB83             2            1        0         0
        .
        .
        .
