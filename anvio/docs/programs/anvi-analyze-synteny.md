The %(anvi-analyze-synteny)s program quantifies %(ngrams)s by transforming contigs into strings of annotations based on a user-specified gene annotation source. A functional annotation source for %(functions)s **must** be provided to generate %(ngrams)s. The program employs a sliding window of size `N` to deconstruct genomic loci of interest into %(ngrams)s and calculate their occurrence frequencies.

### Run for a given function annotation source

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     --annotation-source %(functions)s \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}

For instance, if you have executed %(anvi-run-ncbi-cogs)s on each %(contigs-db)s used to generate your %(genomes-storage-db)s, your `--annotation-source` parameter can be specified as `NCBI_COGS`:

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     --annotation-source NCBI_COGS \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}


### Handling genes with unknown functions 

By default, %(anvi-analyze-synteny)s excludes genes with unknown functions based on the specified annotation source. However, this behavior can be modified through two alternative approaches. First, by providing a %(pan-db)s, which enables the program to utilize gene cluster identities as functional annotations:

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     -p %(pan-db)s \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s
{{ codestop }}

Alternatively, you can explicitly instruct the program to consider genes with unknown functions, which will include ngrams containing functionally unannotated genes in the analysis:

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db)s \
                     --annotation-source %(functions)s \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)s \
                     --analyze-unknown-functions
{{ codestop }}

The primary limitation of this latter approach is that all genes lacking functional annotations are treated as identical entities, which may artificially inflate the frequency of ngrams containing unannotated genes in your final results.

### Run with multiple annotations

When multiple gene annotation sources are provided (such as a pangenome database for gene cluster identities in addition to a functional annotation source), you must specify which annotation source will be used to construct the %(ngrams)s using the `--ngram-source` parameter. The resulting %(ngrams)s will subsequently be re-annotated with the secondary annotation source and reported accordingly. 

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

First, navigate to any working directory and execute the following commands:

``` bash
anvi-self-test --suite metagenomics-full \
               --output-dir TEST_OUTPUT
```

Execute one or more alternative scenarios and examine the output files:

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
