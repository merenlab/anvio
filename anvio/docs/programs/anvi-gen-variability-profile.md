
This program takes the variability data stored within a %(profile-db)s and compiles it from across samples into a single matrix that comprehensively describes your SNVs, SCVs or SAAVs (a %(variability-profile-txt)s).

This program is described in detail in [this blog post](http://merenlab.org/2015/07/20/analyzing-variability/#the-anvio-way), which also covers the biological motivation, output column definitions, and example use cases.

## A note on default filtering of variability signal during profiling

Before this program even runs, anvi'o applies a dynamic, coverage-dependent filter during the profiling step (i.e., when you run %(anvi-profile)s on your BAM files) to reduce the reporting of variation due to sequencing errors and other factors that are difficult to distinguish frmo noise. Specifically, the base frequencies observed at a position in the read recruitment data are only stored in the %(profile-db)s *if* `departure_from_reference` value for that position, which quantifies the fraction of reads that differ from the reference nucleotide at a given position, meets or exceeds a coverage-dependent minimum threshold defined by:

$$y = \left(\frac{1}{b}\right)^{x^{\frac{1}{b}} - m} + c$$

where $x$ is the coverage depth and the model parameters are $b = 2$, $m = 1.45$, $c = 0.05$. This gives a threshold of approximately 0.17 at 20X coverage, 0.07 at 50X, and asymptotically approaches 0.05 at very high coverage. Positions that do not meet this threshold are silently excluded from the database and will never appear in the output of this program.

The user can bypass this filter entirely and store all observed variation regardless of frequency by including `--report-variability-full` when running %(anvi-profile)s. But the use of this parameter will dramatically increase database size, and should only be used with extreme care and attention (this is a very politically correct way to say "please do not use this parameter unless you are certain that this is what you need").

## Basic usage

Here is a basic run with no bells or whistles:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C DEFAULT \
                             -b EVERYTHING
{{ codestop }}

Note that this program requires you to specify a subset of the data to focus on, so to work with everything in the databases, first run %(anvi-script-add-default-collection)s and use the resulting %(collection)s and %(bin)s as shown above.

You can also specify an output file path, which is useful when running multiple times with different `--engine` settings:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --output-file /path/to/your/variability.txt
{{ codestop }}

## Choosing what to analyze

### Focusing on a subset of the input

Instead of a collection and bin, there are three alternatives:

1. Provide gene caller IDs directly or as a file:

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 --gene-caller-ids 1,2,3
    {{ codestop }}

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 --genes-of-interest %(genes-of-interest-txt)s
    {{ codestop }}

2. Provide a %(splits-txt)s to focus on a specific set of splits:

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 --splits-of-interest %(splits-txt)s
    {{ codestop }}

3. Provide any %(collection)s and %(bin)s:

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 -C %(collection)s \
                                 -b %(bin)s
    {{ codestop }}

### Restricting to specific samples

You can limit the analysis to a subset of samples by providing a file with one sample name per line:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --samples-of-interest my_samples.txt
{{ codestop }}

where `my_samples.txt` looks like:

{{ codestart }}
DAY_17A
DAY_18A
DAY_22A
...
{{ codestop }}

### Excluding intergenic positions

For nucleotide-level analyses, you can restrict output to positions that fall within gene calls:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --exclude-intergenic
{{ codestop }}

## Choosing the variability type (engine)

Which type of variants you analyze depends on the `--engine` parameter:

| Engine | Variant type | Requires |
|--------|-------------|---------|
| `NT` (default) | Single nucleotide variants (SNVs) | standard profiling |
| `CDN` | Single codon variants (SCVs) | `--profile-SCVs` at profiling time |
| `AA` | Single amino acid variants (SAAVs) | `--profile-SCVs` at profiling time |

To analyze SAAVs:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine AA
{{ codestop }}

To analyze SCVs:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine CDN
{{ codestop }}

## Adding structural annotations

You can add structural annotations by providing a %(structure-db)s:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C DEFAULT \
                             -b EVERYTHING \
                             -s %(structure-db)s
{{ codestop }}

When a %(structure-db)s is provided, you can also limit your analysis to only genes that have structures in the database:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -s %(structure-db)s \
                             --only-if-structure
{{ codestop }}

## Filtering the output

### By departure from reference or consensus

`departure_from_reference` is the fraction of reads at a position that differ from the reference nucleotide. `departure_from_consensus` is similar but measured against the most frequent allele in that sample. You can set minimum and maximum bounds on either:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --min-departure-from-reference 0.05 \
                             --max-departure-from-reference 0.90
{{ codestop }}

### By minimum occurrence across samples

To keep only positions that are variable in at least N samples (useful for reducing stochastic noise):

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --min-occurrence 3
{{ codestop }}

### By coverage across all samples

To remove any position that has insufficient coverage in even one sample (requires `--quince-mode` to have coverage data for all samples at all positions):

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --quince-mode \
                             --min-coverage-in-each-sample 10
{{ codestop }}

### By number of positions per split

To randomly subsample variable positions and keep at most N per split:

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --num-positions-from-each-split 100
{{ codestop }}

## Special reporting modes

### --quince-mode

By default, if a position is variable in only some samples, the other samples will have no entry for that position in the output. With `--quince-mode`, the program goes back to the raw data and fills in allele frequencies for every sample at every reported position, even those where no variation was detected. This is essential for statistical approaches that require a complete matrix.

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --quince-mode
{{ codestop }}

Note that `--quince-mode` substantially increases runtime and output file size.

### --kiefl-mode

When using `--engine AA` or `--engine CDN`, the default behavior reports only positions that had detectable variation during profiling. With `--kiefl-mode`, all positions in the analyzed genes are reported, with invariant positions given a reference allele frequency of 1. This is useful for analyses that need a complete picture of every codon or amino acid position, not just the variable ones. Incompatible with `--quince-mode`.

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine CDN \
                             --kiefl-mode
{{ codestop }}

This flag was added in this [pull request](https://github.com/merenlab/anvio/pull/1794) where you can read about the tests performed to validate its behavior.

## Additional output columns

### Contig and split names

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --include-contig-names \
                             --include-split-names
{{ codestop }}

Contig names are excluded by default since they can nearly double file size.

### Gene-level coverage statistics

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --compute-gene-coverage-stats
{{ codestop }}

This appends per-gene coverage statistics to each row. It is computationally expensive and off by default.

### Per-site pN/pS values (CDN engine only)

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine CDN \
                             --include-site-pnps
{{ codestop }}

This adds 12 columns of per-site synonymous and nonsynonymous substitution information, computed relative to the reference, the consensus, and the most common consensus across samples.

### Additional data from the database

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine AA \
                             --include-additional-data
{{ codestop }}

This appends any data stored in the `amino_acid_additional_data` table as extra columns. Currently only supported for the `AA` engine.
