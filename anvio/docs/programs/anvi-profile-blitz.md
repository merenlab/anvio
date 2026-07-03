This program **produces a %(bam-stats-txt)s from one or more %(bam-file)s given a %(contigs-db)s**. It is designed to serve people who only need to process read recruitment data stored in a %(bam-file)s to recover coverage and detection statistics as well as the number of mapped reads (along with other statistics) for their genes, contigs, and/or genomes. It will report what's going on nicely with memory usage information and estimated time of completion:

[![anvi-profile-blitz](../../images/anvi-profile-blitz.png){:.center-img}](../../images/anvi-profile-blitz.png)

There are other programs in anvi'o software ecosystem that are similar to this one:

* %(anvi-profile)s also takes a %(bam-file)s and profiles it. **They both require a %(contigs-db)s**. But while %(anvi-profile)s produces a %(single-profile-db)s for downstream analyses in anvi'o, %(anvi-profile-blitz)s produces text files for downstream analyses by the user (via R, Python, or other solutions). In contrast to %(anvi-profile)s, %(anvi-profile-blitz)s is orders of magnitude faster with similar memory usage.

* %(anvi-script-get-coverage-from-bam)s also takes a %(bam-file)s and profiles it. **They both produce text output files.** But while %(anvi-script-get-coverage-from-bam)s does not require a %(contigs-db)s, %(anvi-profile-blitz)s requires one to work. They will both run very rapidly, %(anvi-script-get-coverage-from-bam)s will work with a much smaller amount of memory.

## Output files

For output file formats, please see %(bam-stats-txt)s.

## Running

You can use this program with one or more BAM files to recover minimal or extended statistics for contigs or genes in a %(contigs-db)s.

{:.warning}
Since the program will not be able to ensure the %(contigs-db)s was generated from the same %(contigs-fasta)s that was used for read recruitment that resulted in %(bam-file)ss for analysis, you can make serious mistakes if you mix up your workflow and start profiling BAM files that have nothing to do with a %(contigs-db)s. If you make a mistake like that, in the best case scenario you will get an empty output file because the program will skip all contigs with non-matching name. In the worst case scenario you will get a file if some names in %(contigs-db)s incorrectly matches to some names in the %(bam-file)s. While this warning may be confusing, you can avoid all these if you use the SAME FASTA FILE both as reference for read recruitment and as input for %(anvi-gen-contigs-database)s.

### Contigs mode, default output

Profile contigs, produce a default output:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   -o OUTPUT.txt
{{ codestop }}

This example is with a single BAM file, but you can also have multiple BAM files as a parameter by using wildcards,

{{ codestart }}
anvi-profile-blitz *.bam \
                   -c %(contigs-db)s \
                   -o OUTPUT.txt
{{ codestop }}

or by providing multiple paths:

{{ codestart }}
anvi-profile-blitz /path/to/SAMPLE-01.bam \
                   /path/to/SAMPLE-02.bam \
                   /another/path/to/SAMPLE-03.bam
                   -c %(contigs-db)s \
                   -o OUTPUT.txt
{{ codestop }}

### Contigs mode, minimal output

Profile contigs, produce a minimal output. This is the fastest option:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --report-minimal \
                   -o OUTPUT.txt
{{ codestop }}

### Genome/bin mode via collections

If you have a %(collection-txt)s, you can summarize coverage/detection statistics per bin/genome instead of per contig.

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   -C %(collection-txt)s \
                   -o OUTPUT.txt
{{ codestop }}

{:.note}
Bin-level detection and coverage are computed by treating all contigs in a bin as a single genome (nucleotide arrays are concatenated under the hood). With `--report-minimal`, these statistics are streamed to keep memory usage close to contig mode even for large bins.

### Genes mode, default output

Instead of contigs, profile genes, produce a default output:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --gene-mode \
                   -o OUTPUT.txt
{{ codestop }}

`--report-minimal` will behave the same, and produce minimal output.

### Genes mode, on a subset of genes

You don't want to profile all genes in your %(contigs-db)s? Then you should tell anvi'o which specific genes you are interested in by either providing a comma-separated list of gene-caller IDs:

{{ codestart }}
# the following will profile 4 specific gene calls
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --gene-mode \
                   -o OUTPUT.txt \
                   --gene-caller-ids 5,13,74,203
{{ codestop }}

_Or_ by providing a %(genes-of-interest-txt)s file:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --gene-mode \
                   -o OUTPUT.txt \
                   --genes-of-interest %(genes-of-interest-txt)s
{{ codestop }}

## Modifying DisCov parameters

For genomes or contigs, %(anvi-profile-blitz)s can compute the Distribution of Coverage (DisCov) score alongside other coverage statistics. DisCov combines a spread score _S_ (proportion of windows with coverage) and an evenness score _E_ (proportion of covered bases within a fold-range of the median nonzero coverage). See %(discov-stats)s for a full description of the metric and parameters.

The default parameters are designed to work well across typical metagenomics data, but you can adjust them:

**Window sizing for _S_**

The spread score (_S_) is computed by dividing the input sequence into non-overlapping windows. By default, contig-level stats use a window length equal to 1% of the contig length (minimum 300 bp). Genome/bin-level stats (when using `--collection-txt`) use a fixed 1,000 bp window. You can override these defaults with the following parameters:

* `--window-length INT` — use a fixed window size in bp for all sequences
* `--window-length-as-percentage FLOAT` — set window length as a percentage of each sequence's length
* `--min-window-length INT` — set a minimum window length floor when using percentage mode

**Fold-range for _E_**

The evenness score (_E_) counts bases with coverage between a lower and upper fold-multiple of the median nonzero coverage. The default fold-range endpoints are 0.5x and 2.0x. To adjust:

* `--foldrange-lower FLOAT` — lower bound (default: 0.5)
* `--foldrange-upper FLOAT` — upper bound (default: 2.0)

**Combining _S_ and _E_**

You can adjust how _S_ and _E_ are combined by changing either their weights or the overall DisCov formula:

* `--alpha FLOAT` — weight of _S_ relative to _E_, in [0, 1] (default: 0.5)
* `--discov-formula STRING` — `linear` (DisCov = _αS_ + (1-_α_)_E_) or `geometric` (DisCov = _S_^_α_ × _E_^(1-_α_)) (default: `linear`)

**Window-level output**

To inspect the per-window values used to compute _S_ (useful for debugging or visualization), add the `--gen-window-level-output` flag:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   -o OUTPUT.txt \
                   --gen-window-level-output
{{ codestop }}

This produces an additional file named `OUTPUT-WINDOWS.txt` with per-window start/stop positions, coverage presence, and base counts within the fold-range.

Note that computing DisCov is not compatible with `--gene-mode` or `--report-minimal`.

## Performance

The memory use will be correlated linaerly with the size of the %(contigs-db)s, but once everything is loaded, the memory usage will not increase substantially over time.

With the flag `--report-minimal`, %(anvi-profile-blitz)s profiled on a laptop computer 100,000 contigs that contained 1 billion nts in 6 minutes and used  ~300 Mb memory. This contigs database had 1.5 million genes, and memory usage increased to 1.7 Gb when %(anvi-profile-blitz)s run in `--gene-mode`. The flag `--gene-mode` does not change time complexity dramatically.

Anvi'o has this program because [Emile Faure](https://twitter.com/faureemile) presented us with a challenge: Emile had a ~140 Gb anvi'o %(contigs-db)s that contained nearly 70 million contig sequences from over 200 single-assembled metagenomes, and wanted to learn the coverages of each gene in the contigs database in 200 metagenomes individually. Yet the combination of %(anvi-profile)s and %(anvi-summarize)s jobs would take **more than 40 days** to complete. Since all Emile needed was to learn the coverages from BAM files, we implemented %(anvi-profile-blitz)s to skip the profiling step. The run took **8 hours to compute and report coverage values for 175 million genes in 70 million contigs**, and the memory use remained below 200 Gb.
