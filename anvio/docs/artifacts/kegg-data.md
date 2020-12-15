A **directory of data** downloaded from the [KEGG database resource](https://www.kegg.jp/) for use in function annotation and metabolism estimation.

It is created by running the program %(anvi-setup-kegg-kofams)s. Not everything from KEGG is included in this directory, only the information relevant to downstream programs. The most critical components of this directory are KOfam HMM profiles and the %(modules-db)s which contains information on metabolic pathways as described in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html).

Programs that rely on this data directory include %(anvi-run-kegg-kofams)s and %(anvi-estimate-metabolism)s.
