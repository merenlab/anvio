A tab-delimited text file linking tRNA-seq modification profiles to the genomes and modification enzymes responsible for them.

This file is produced by %(anvi-export-genome-specific-trna-modifications)s. Each row corresponds to a tRNA seed detected in a specific sample, at a specific canonical position, in the context of a specific genomic tRNA gene and the modification enzyme encoded by that genome.

## Columns

| Column | Description |
| --- | --- |
| `modifying_enzyme_name` | Short name of the modification enzyme |
| `modification` | Modification produced |
| `function_name` | Full name of the genomic function encoding the enzyme |
| `function_accession` | Accession of the genomic function |
| `gene_callers_id` | Gene callers ID of the tRNA seed in the tRNA-seq contigs database |
| `contig_name` | Contig name of the tRNA seed in the tRNA-seq contigs database |
| `anticodon` | Anticodon of the tRNA seed |
| `aa` | Amino acid of the tRNA seed |
| `domain` | Taxonomic domain of the seed |
| `phylum` | Taxonomic phylum |
| `class` | Taxonomic class |
| `order` | Taxonomic order |
| `family` | Taxonomic family |
| `genus` | Taxonomic genus |
| `species` | Taxonomic species |
| `taxon_percent_id` | Percent identity of the taxonomic hit |
| `genome_name` | Name of the external genome containing the modification enzyme and tRNA gene |
| `mean_coverage` | Mean specific coverage of the tRNA seed across all positions in this sample |
| `relative_mean_coverage` | Mean coverage of the seed relative to all seeds in this sample |
| `relative_discriminator_coverage` | Discriminator position coverage relative to all seeds in this sample |
| `coverage_at_position` | Coverage at the canonical modification position (sum of A+C+G+T for SNV rows; position-specific coverage from %(seeds-specific-txt)s for no-SNV rows); `NA` if below the minimum coverage threshold |
| `seed_position` | 0-based index of the modification position in the seed sequence; `NA` if no SNV was detected at this position in any sample |
| `ordinal_name` | Structural name of the position (e.g., `anticodon_loop_3`); `NA` if no SNV detected |
| `ordinal_position` | Ordinal index of the position across all possible tRNA positions; `NA` if no SNV detected |
| `canonical_position` | Canonical tRNA position number (e.g., `34`, `37`) |
| `reference_abundant_nucleotide` | Most abundant nucleotide observed at this position in this sample; set to `reference_genome` for no-SNV rows |
| `reference_genome` | Unmodified nucleotide at this position according to the tRNA gene sequence in the genome |
| `sample_name` | tRNA-seq sample name |
| `A` | Count of A reads at this position; `NA` for no-SNV rows |
| `C` | Count of C reads; `NA` for no-SNV rows |
| `G` | Count of G reads; `NA` for no-SNV rows |
| `T` | Count of T reads; `NA` for no-SNV rows |
| `modified.fraction.reference_genome` | Fraction of reads differing from the genomic reference nucleotide (modification signal); `NA` if coverage is below threshold; `0.0` for no-SNV rows with sufficient coverage |
| `modified.fraction.reference_abundant_nucleotide` | Fraction of reads differing from the most abundant observed nucleotide; `NA` if below threshold; `0.0` for no-SNV rows with sufficient coverage |

## Notes

Rows in this table fall into two categories:

- **SNV rows**: a substitution was detected at this position in this sample by `anvi-merge-trnaseq`, indicating a modification signal. Raw nucleotide counts (A, C, G, T) and modification fractions are reported directly from %(modifications-txt)s.
- **No-SNV rows**: the seed maps to a genome encoding the enzyme, but no substitution was detected at this position in this sample. The tRNA gene is present and the modification site has the expected reference nucleotide, but the modification was not detected (fraction = 0.0). `coverage_at_position` is drawn from the position-specific coverage column in %(seeds-specific-txt)s. Coverage below the minimum threshold is reported as `NA`.

The companion output from `--enzyme-distribution-output` reports how many genes encoding each enzyme are present in each genome.
