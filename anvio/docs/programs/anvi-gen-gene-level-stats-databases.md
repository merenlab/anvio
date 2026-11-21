This program **generates a %(genes-db)s, which stores the coverage and detection values for all of the genes in your %(contigs-db)s.**

This information is usually calculated when it's needed (for example when running %(anvi-interactive)s in genes mode), but this program lets you break this process into two steps. This way, you can easily change the parameters of %(anvi-interactive)s without having to recalculate the gene-level statistics.

Given a %(contigs-db)s and %(profile-db)s pair, as well as a %(collection)s, this program will calculate the stats for the genes in each of your %(bin)ss and give each bin its own %(genes-db)s that includes this information.

For example, if a %(collection)s called `GENE_COLLECTION` contained the bins `bin_0001`, `bin_0002`, and `bin_0003` and you ran:

{{ codestart }}
anvi-gen-gene-level-stats-databases -c %(contigs-db)s \
                                    -p %(profile-db)s \
                                    -C GENE_COLLECTION
{{ codestop }}

Then it will create a directory called `GENES` next to your %(profile-db)s that contains three %(genes-db)ss called `GENE_COLLECTION-bin_0001.db`, `GENE_COLLECTION-bin_0002.db`, and `GENE_COLLECTION-bin_0003.db`. Once these databases are there, running %(anvi-interactive)s in gene mode for any of these bins will use them directly, instead of computing gene-level stats from scratch.

You also have the option to provide a list of %(bin)s (either as a file or as a string) to analyze instead of every bin in a %(collection)s.

If genes databases already exist for some of your bins, anvi'o will ask you to use the flag `--just-do-it` before it does anything. With that flag, it will keep each existing genes database that was generated with the same parameters and the same splits, and regenerate the ones that weren't.

### Normalized gene coverages

By default, this program also computes five normalizations of the `non_outlier_mean_coverage` values of each gene in each sample, and stores them in the genes database. They become additional views when you run %(anvi-interactive)s in gene mode:

- **log1p**: log(x + 1) of coverage values, which reduces their dynamic range.
- **cov_per_million**: mean coverage depth per million mapped reads in each sample, which accounts for differences in sequencing depth across samples. Please note that this is not CPM (counts per million), since it is based on coverage depth rather than read counts.
- **zscore_raw**, **zscore_log1p**, and **zscore_cov_per_million**: per-gene z-scores of raw, log1p, and cov_per_million values across samples, which show how many standard deviations a gene's coverage in a given sample deviates from that gene's mean coverage across all samples.

Z-scores are computed for each gene separately, so they are not comparable between genes within the same sample. With transcriptomic data, they reflect each gene's own expression pattern. With metagenomic data, genes of a single population tend to follow the abundance of that population, so z-scores mostly help spot genes whose coverage patterns deviate from the rest of their genome (such as those that are found in multiple copies, or are absent from some samples). In either case, setting a __center__ value of 0 for z-score views in the %(interactive)s interface will display them as diverging values.

If an existing genes database does not have normalized coverages yet (e.g., because it was generated with an earlier version of anvi'o), this program will compute them using the gene-level stats already stored in it.

You can skip this step with the flag `--skip-computing-gene-level-normalized-coverages`, in which case these data will not be generated, and these views will not be available:

{{ codestart }}
anvi-gen-gene-level-stats-databases -c %(contigs-db)s \
                                    -p %(profile-db)s \
                                    -C %(collection)s \
                                    --skip-computing-gene-level-normalized-coverages
{{ codestop }}

{:.notice}
Z-scores are not meaningful when there is only a single sample. So, if your %(profile-db)s has only one sample, anvi'o will refuse to create a new genes database until you add `--skip-computing-gene-level-normalized-coverages` to your command, so it is clear that these values will not be available.

{:.notice}
Computing coverage per million mapped reads requires the number of reads mapped to each sample, which some very old %(profile-db)ss may not have. In that case anvi'o will skip computing normalized coverages with a warning.

### Other Parameters

You can also change the definition of an outlier nucleotide position (with `--outliers-threshold`), or whether nucleotide positions with zero coverage should be treated as outliers (with `--zeros-are-outliers`).
