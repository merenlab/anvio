This program takes advantage of paired-end read technology and uses reverse-forward (`RF`) read pairs to determine whether a given contig in a given %(bam-file)s is circular, linear, or indeterminate by making use of the insert size distributions in the entirety of the %(bam-file)s.

### Example run

{{ codestart }}
anvi-report-circularity ANY-SORTED-BAM-FILE.bam \
                        -o CIRCULARITY-REPORT.txt
{{ codestop }}

### How it works

The following describes what %(anvi-report-circularity)s does step by step while shedding light on the purpose and relevance of the command-line parameters.

* **Estimate insert-size statistics**. The run starts by trying to establish an understanding of the insert size statistics by sampling up to `--max-num-pairs-for-is-est` forward-reverse (`FR`) pairs from across all contigs (this is to make sure the program does not process ALL reads in a %(bam-file)s). Then, it computes the median insert size for the library, `M`, and its [median absolute deviation](https://en.wikipedia.org/wiki/Median_absolute_deviation), `MAD`. For this to work, a minimum of `--min-pairs-for-stats` `FR` pairs is required, otherwise we can't have an `M` and so the program stops with an error already. `MAD` should NEVER be zero (since it would indicate that you have a hyper-uniform library, but it is not possible), but in case you are running the program with simulated data, the program simply assumes `MAD` equals 1.0 and moves on.

{:.notice}
If you prepared your genomic or shotgun metagenomic libraries without a proper size selection step, this program will unlikely work well since `M` will be irrelevant for most pairs, and `MAD` will be relatively useless.

* **Scan each contig**. Then the program proceeds to read in the %(bam-file)s for each contig to count `FR` pairs, and keep track of reads near contig edges. If the length of a contig, `L` is shorter than `2 × M`, they will be ignore ad this step as their miserable length would make it impossible for the downstream code to confidently determine if they are circular or not. These contigs are marked as `indeterminate` in the output file with a warning flag.

* **Compute the expected number of RFs**. The algorithm then calculates the expected number of `RF` pairs for a given contig of length `L` to be circular using the following logic,
  
   `E[RF] = N_FR × M / (L - M)`  

   where `N_FR` is the number of `FR` pairs on the contig. This reflects the chance that an `RF` pair spans a circular junction when linearized.

{:.notice}
The efficacy of the step above heavily depends on your assembler's restraint to not add stupid extras such as non-existent tandem repeats to contig ends (as they generally do with circular entities) :/

* **Score RFs that supports circularity**. For each `RF` pair, the algorithm then computes the 'circular insert size', let's call it `circular_insert`. Circular insert size is a new insert size value calculated by the formula `left_end + (L - right_start)` and simply answers the question, "*what the insert size of this `RF` paired end would have been __if__ the contig was indeed circular in the environment?". A given `RF` is one that 'supports' circularity, if `|circular_insert - M| ≤ T`, where `T` is tolerance, and it is calculated using `MAD` and the user-defined tolerance factor: `T = --insert-tolerance-factor × MAD`. Once all `RF` paris are concidered, the overall circularity support for a given contig is the fraction of observed `RF` pairs that are supporting.

{:.notice}
This step will most likely miss pro-viruses that occur in high-coverage chromosomes *and* found as circular genomes in capsids.

* **Assess edge coherence**. This is not immediately relevant to the circularity calculation, but it is something good to have in the report: the proportion of paired-end reads at the extremeties of contigs where both mates map *into* the contig. If both mates in every single paired-end read that occur at the contig edge, which is defined by thos that are within `M` nucleotides form either of the contig edges, the edge coherence is maximum (which may indicate that a given contig is linear, and fragment of no one). If the is a large number of paired-end reads in these edges with their mate mapping elsewhere or nowhere, then the edge coherence is minimum (which may indicate that the contig is a fragment of something).

* **Report**. The program concludes by generating an output file that lists for each contig their status, support scores, coverage counts, and warning flags.

### Details of the decision making

The program makes use of a combination of user-defined thresholds and the `min_required` variable which is equivalent to `--min-supporting-pairs` or `--expected-fraction-threshold × E[RF]`, whichever is larger.

Given all these, the final decision for contig status is made based on the rules below, which are evaluated in this order:

* **Circular** if supporting `RF` pairs ≥ `min_required` **and** circularity support ≥ `--circularity-confidence-threshold`.
* **Circular** if supporting `RF` pairs ≥ `--min-supporting-pairs` **and** ≥80%% of observed `RF` pairs are supporting.
* **Linear** if no `RF` pairs are observed **and** there are at least 100 `FR` pairs on the contig.
* **Linear** if supporting `RF` pairs < `min_required` **and** circularity support < `--circularity-confidence-threshold`.
* Otherwise **indeterminate** (flagged as ambiguous evidence). Contigs with no `FR` pairs, very short lengths (`L < 2M`), or extremely low edge coverage are also reported as indeterminate with somewhat informative flags.
