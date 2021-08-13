[Kegg Orthology](https://www.genome.jp/kegg/ko.html) (KO) functional annotations, produced by finding HMM hits to the KEGG KOfam database.

You can annotate a %(contigs-db)s with these KEGG functions by running %(anvi-run-kegg-kofams)s. They will be added to the gene functions table under the source 'KOfam'.

Another program that relies on these annotations is %(anvi-estimate-metabolism)s, which uses them to determine presence and completeness of metabolic pathways that are defined by KOs.
