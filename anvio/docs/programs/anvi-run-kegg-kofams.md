Essentially, this program uses the KEGG database to annotate functions and metabolic pathways in a %(contigs-db)s. More specifically, %(anvi-run-kegg-kofams)s annotates a %(contigs-db)s with HMM hits from KOfam, a database of KEGG Orthologs (KOs). You must set up these HMMs on your computer using %(anvi-setup-kegg-kofams)s before you can use this program.

Briefly, what this program does is extract all the gene calls from the %(contigs-db)s and checks each one for hits to the KOfam HMM profiles in your %(kegg-db)s. This can be time-consuming given that the number of HMM profiles is quite large, and especially time-consuming if the number of genes in the %(contigs-db)s is also large. Multi-threading is a good idea if you have the computational capability to do so.

Many HMM hits will be found, most of them weak. The weak hits will by default be eliminated according to the score thresholds provided by KEGG; that is, only hits with scores above the threshold for a given KO profile will be annotated in the %(contigs-db)s. It is perfectly normal to notice that the number of raw hits found is many, many times larger than the number of annotated KO hits in your database.

In the %(contigs-db)s functions table, annotated KO hits (%(kegg-functions)s) will have the source `KOfam`.

Running this program is a pre-requisite for metabolism estimation with %(anvi-estimate-metabolism)s. Note that if you are planning to run metabolism estimation, it must be run with the same %(kegg-db)s that is used in this program to annotate KOfam hits.

### Standard usage

{{ codestart }}
anvi-run-kegg-kofams -c CONTIGS.db
{{ codestop }}

### Use a specific non-default KEGG data directory

{{ codestart }}
anvi-run-kegg-kofams -c CONTIGS.db --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

### Run with multiple threads

{{ codestart }}
anvi-run-kegg-kofams -c CONTIGS.db -T 4
{{ codestop }}

### Use a different HMMER program
By default, %(anvi-run-kegg-kofams)s uses `hmmsearch` to find KO hits. If for some reason you would rather use a different program (`hmmscan` is also currently supported), you can do so.

{{ codestart }}
anvi-run-kegg-kofams -c CONTIGS.db --hmmer-program hmmscan
{{ codestop }}

### Keep all HMM hits
Usually, this program parses out weak HMM hits and keeps only those that are above the score threshold for a given KO. If you would like to turn off this behavior and keep all hits (there will be _a lot_ of weak ones), you can follow the example below:

{{ codestart }}
anvi-run-kegg-kofams -c CONTIGS.db --keep-all-hits
{{ codestop }}
