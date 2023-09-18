This program allows you to remove low quality HMM alignments from a %(hmm-source)s in a %(contigs-db)s with HMM alignment parameters such as model-coverage (query-coverage) and gene-coverage (target-coverage), or by removing partial genes (i.e., genes that are not partial and that start with a start codon and end with a stop codon). Briefly, the program will remove all records from an %(hmm-source)s in the %(hmm-hits)s, then import a new %(hmm-hits)s table into the %(contigs-db)s that was filtered to your specifications.

## Filter with HMM alignment parameters

Similar to query coverage in BLAST, we can also use HMM alignment coverage to help determine if an hmm-hit is homologous. A small coverage value means only a small proportion of the query/target is aligning. Before anvi'o can filter out %(hmm-hits)s with alignment coverage, you must run %(anvi-run-hmms)s and report a domain hits table by including `--domain-hits-table` flag in your command:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71 \
              --hmmer-output-dir path/to/dir
              --domain-hits-table
{{ codestop }}

After the command above, your HMM hits will be stored in your %(contigs-db)s as usual. However, with the domain hits table, you can filter out hits from your %(contigs-db)s using thresholds for model or gene coverage of each hit i.e. you can filter out %(hmm-hits)s where the profile HMM and gene align well to each other.

For example, following the command above, the command below will remove %(hmm-hits)s from your %(contigs-db)s for profile HMMs that had less than 90%% coverage of the target genes:

{{ codestart }}
anvi-script-filter-hmm-hits-table -c %(contigs-db)s \
                                  --hmm-source Bacteria_71 \
                                  --domain-hits-table path/to/dir/hmm.domtable \
                                  --model-coverage 0.9
{{ codestop }}

### HMMs with multiple hits to one gene

Some HMM profiles align multiple times to the same gene at different coordinates. The program `anvi-script-filter-hmm-hits-table` by default will use only one of those domain hits table records which could represent very little alignment coverage. To combine the domain hits table records into one hit and thus increasing alignment coverage, use the parameter `--merge-partial-hits-within-X-nts`. Briefly, if you give the parameter `--merge-partial-hits-within-X-nts` 300, `anvi-script-filter-hmm-hits-table` will merge all hits to the same gene in the domain hits table that have coordinates within 300 nucleotides of each other.  

{{ codestart }}
anvi-script-filter-hmm-hits-table -c %(contigs-db)s \
                                  --hmm-source Bacteria_71 \
                                  --domain-hits-table path/to/dir/hmm.domtable \
                                  --model-coverage 0.9 \
                                  --merge-partial-hits-within-X-nts
{{ codestop }}

{:.notice}
The input domtblout file for %(anvi-script-filter-hmm-hits-table)s will be saved as `hmm.domtable.orig` and the output, filtered version will be saved as `hmm.domtable`. If you decide to change the coverage filtering threshold or `--merge-partial-hits-within-X-nts`, be sure to change the path for `--domain-hits-table`  to `hmm.domtable.orig`.

## Filter out hmm-hits from partial genes

HMMs are able to detect partial genes (i.e., genes that are not partial and that start with a start codon and end with a stop codon) with good alignment coverage and homology statistics. However, partial genes can lead to spurious phylogenetic branches and/or inflate the number of observed populations or functions in a given set of genomes/metagenomes. Using `--filter-out-partial-gene-calls`, you can remove partial gene hmm-hits.

{{ codestart }}
anvi-script-filter-hmm-hits-table -c %(contigs-db)s \
                                  --hmm-source Bacteria_71 \
                                  --domain-hits-table path/to/dir/hmm.domtable \
                                  --filter-out-partial-gene-calls
{{ codestop }}