Essentially, this program uses the KEGG database to annotate functions and metabolic pathways in a %(contigs-db)s. More specifically, %(anvi-run-kegg-kofams)s annotates a %(contigs-db)s with HMM hits from KOfam, a database of KEGG Orthologs (KOs). You must set up these HMMs on your computer using %(anvi-setup-kegg-kofams)s before you can use this program. Membership of KOfam functions in KEGG metabolic MODULES and BRITE hierarchies is also stored in the %(contigs-db)s.

Running this program is a pre-requisite for metabolism estimation with %(anvi-estimate-metabolism)s. Note that if you are planning to run metabolism estimation, it must be run with the same %(kegg-data)s that is used in this program to annotate KOfam hits.

## How does it work?
**1) Run an HMM search against KOfam**
Briefly, what this program does is extract all the gene calls from the %(contigs-db)s and checks each one for hits to the KOfam HMM profiles in your %(kegg-data)s. This can be time-consuming given that the number of HMM profiles is quite large, even more so if the number of genes in the %(contigs-db)s is also large. Multi-threading is a good idea if you have the computational capability to do so.

**2) Eliminate weak hits based on bitscore**
Many HMM hits will be found, most of them weak. The weak hits will by default be eliminated according to the bitscore thresholds provided by KEGG; that is, hits with bitscores below the threshold for a given KO profile will be discarded, and those with bitscores above the threshold will be annotated in the %(contigs-db)s. It is perfectly normal to notice that the number of raw hits found is many, many times larger than the number of annotated KO hits in your database.

**3) Add back valid hits that were missed**
There is one issue with this practice of removing _all_ KOfam hits below the KEGG bitscore threshold for a given profile. We (and others) have noticed that the KEGG thresholds can sometimes be too stringent, eliminating hits that are actually valid annotations. To solve this problem, we
have implemented the following heuristic for relaxing the bitscore thresholds and annotating genes that would otherwise go without a valid KO annotation:

For every gene without a KOfam annotation, we examine all the hits with an e-value below `X` and a bitscore above `Y` percent of the threshold. If those hits are all to a unique KOfam profile, then we annotate the gene call with that KO.

`X` and `Y` are parameters that can be modified (see below), but by default the e-value threshold (`X`) is 1e-05 and the bitscore fraction (`Y`) is 0.5.

Please note that this strategy is just a heuristic. We have tried to pick default parameters that seemed reasonable but by no means have we comprehensively tested and optimized them. This is why X and Y are mutable so that you can explore different values and see how they work for your data. It is always a good idea to double-check your annotations to make sure they are reasonable and as stringent as you'd like them to be. In addition, if you do not feel comfortable using this heuristic at all, you can always turn this behavior off and rely solely on KEGG's bitscore thresholds. :)

**3) Put annotations in the database**
In the %(contigs-db)s functions table, annotated KO hits (%(kegg-functions)s) will have the source `KOfam`. Metabolic Modules and BRITE functional classifications containing these functions also have entries in the table, with sources labeled `KEGG_Module` and `KEGG_BRITE`. BRITE classification will not occur if %(anvi-setup-kegg-kofams)s was not set up with BRITE data (see the artifact for that program to see how to include BRITE).

## Standard usage

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s
{{ codestop }}

## Use a specific non-default KEGG data directory
If you have previously setup your KEGG data directory using `--kegg-data-dir` (see %(anvi-setup-kegg-kofams)s), or have moved the KEGG data directory that you wish to use to a non-default location (maybe you like keeping the older versions around when you update, we don't know how you roll), then you may need to specify where to find the KEGG data so that this program can use the right one. In that case, this is how you do it:

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s \
                     --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

## Run with multiple threads

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s -T 4
{{ codestop }}

## Use a different HMMER program
By default, %(anvi-run-kegg-kofams)s uses `hmmsearch` to find KO hits. If for some reason you would rather use a different program (`hmmscan` is also currently supported), you can do so.

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s \
                     --hmmer-program hmmscan
{{ codestop }}

## Keep all HMM hits
Usually, this program parses out weak HMM hits and keeps only those that are above the score threshold for a given KO. If you would like to turn off this behavior and keep all hits (there will be _a lot_ of weak ones), you can follow the example below:

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s \
                     --keep-all-hits
{{ codestop }}

## Modify the bitscore relaxation heuristic
As described above, this program does its best to avoid missing valid annotations by relaxing the bitscore threshold for genes without any annotations. For such a gene, hits with e-value <= X and bitscore > (Y * KEGG threshold) that are all hits to the same KOfam profile are used to annotate the gene with that KO.

### Skip this heuristic entirely
If you don't want any previously-eliminated hits to be used for annotation, you can skip this heuristic by using the flag `--skip-bitscore-heuristic`. Then, _only_ hits with bitscores above the KEGG-provided threshold for a given KO will be used for annotation.

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s \
                     --skip-bitscore-heuristic
{{ codestop }}

### Modify the heuristic parameters
If our default values are too stringent or not stringent enough for your tastes, you can change them! The e-value threshold (X, default: 1e-05) can be set using `-E` or `--heuristic-e-value` and the bitscore fraction (Y, default: 0.50) can be set using `-H` or `--heuristic-bitscore-fraction`. Like so:

{{ codestart }}
anvi-run-kegg-kofams -c %(contigs-db)s \
                     -E 1e-15 \
                     -H 0.90
{{ codestop }}
