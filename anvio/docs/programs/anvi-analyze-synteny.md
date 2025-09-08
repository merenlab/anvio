Briefly, %(anvi-analyze-synteny)s counts %(ngrams)s by converting contigs into strings of annotations for a given user-defined source of gene annotation. A source annotation for %(functions)s **must** be provided to create %(ngrams)s, upon which anvi'o will use a sliding window of size `N` to deconstruct the loci of interest into %(ngrams)s and count their frequencies.

### Run for a given function annotation source

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     --annotation-source %(functions)s \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}

For instance, if you have run %(anvi-run-ncbi-cogs)s on each %(contigs-db)s you have used to generate your %(genomes-storage-db)s, your `--annotation-source` can be `NCBI_COGS`:

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     --annotation-source NCBI_COGS \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}


### Handling genes with unknown functions 

By default, %(anvi-analyze-synteny)s will ignore genes with unknown functions based on the annotation source of interest. However, this can be circumvented either by providing a %(pan-db)s, so the program would use gene cluster identities as function names:

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     -p %(pan-db)s \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}

or by explicitly asking the program to consider unknown functions, in which case the program would not discard ngrams that include genes without functions:

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     --annotation-source %(functions)s \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s \
                     --analyze-unknown-functions
{{ codestop }}

The disadvantage of the latter strategy is that since all genes with unknown functions will be considered the same, the frequency of ngrams that contain genes with unknown functions may be inflated in your final results.

### Run with multiple annotations

If multiple gene annotation sources are provided (i.e., a pangenome for gene clusters identities as well as a functional annotation source), the user must define which annotation source will be used to create the %(ngrams)s using the parameter `--ngram-source`. The resulting %(ngrams)s will then be re-annotated with the second annotation source and also reported. 

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     -p %(pan-db)s \
                     --annotation-source %(functions)s \
                     --ngram-source gene_clusters \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}

### Test cases for developers

If you are following the anvi'o master branch on your computer, you can create a test case for this program.

First, go to any work directory, and run the following commands:

``` bash
anvi-self-test --suite metagenomics-full \
               --output-dir TEST_OUTPUT
```

Run one or more alternative scenarios and check output files:

```
anvi-analyze-synteny -g TEST_OUTPUT/TEST-GENOMES.db \
                     --annotation-source COG20_FUNCTION \
                     --ngram-window-range 2:3 \
                     -o TEST_OUTPUT/synteny_output_no_unknowns.tsv

anvi-analyze-synteny -g TEST_OUTPUT/TEST-GENOMES.db \
                     --annotation-source COG20_FUNCTION \
                     --ngram-window-range 2:3 \
                     -o TEST_OUTPUT/synteny_output_with_unknowns.tsv \
                     --analyze-unknown-functions

anvi-analyze-synteny -g TEST_OUTPUT/TEST-GENOMES.db \
                     --annotation-source COG20_FUNCTION \
                     --ngram-window-range 2:3 \
                     -o TEST_OUTPUT/tsv.txt \
                     --analyze-unknown-functions
```
