The input for %(anvi-analyze-synteny) is a %(genomes-storage-db) which should contain one or more similar loci or genomes, and an annotation source. 

Briefly, %(anvi-analyze-synteny) counts %(ngrams) by converting contigs into strings of annotations with the user-defined annotation source. An target annotation source **must** be designated to create %(ngrams) (e.g. gene_clusters and COG_FUNCTION). Next, it uses a sliding window of size N to deconstruct the loci of interest into %(ngrams). By default, %(anvi-analyze-synteny) will ignore genes that are not annotated. 



### Run with gene-cluster annotations from pangenome

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db) \
                     -p %(pan-db) \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)
{{ codestop }}

### Run with NCBI COGs annotations

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db) \
                     --annotation-source COG_FUNCTION \
                     --ngram-window-range 2:3 \
                     -o %(ngrams) 
{{ codestop }}

### Run with multiple annotations

If multiple annotation sources are provided, the user must define which annotation source will be used to create the %(ngrams) using the parameter `--ngram-source`. The resulting %(ngrams) will then be re-annotated with the second annotation source and also reported. 

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db) \
                     -p %(pan-db) \
                     --annotation-source COG_FUNCTION \
                     --ngram-source gene_clusters \
                     --ngram-window-range 2:3 \
                     -o %(ngrams)
{{ codestop }}

### Record genes without annotations

%(ngrams) that contain genes that are not annotated can be recorded by using the parameter `--analyze-unknown-functions`

{{ codestart }}
anvi-analyze-synteny -g %(genomes-storage-db) \
                     --annotation-source COG_FUNCTION \
                     --ngram-window-range 2:3 \
                     -o %(ngrams) \
                     --analyze-unknown-functions
{{ codestop }}
