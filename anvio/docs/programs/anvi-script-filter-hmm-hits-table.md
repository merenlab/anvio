This program allows you to remove low quality HMM alignments from a %(hmm-source)s in a %(contigs-db)s by leveraging HMM alignment parameters such as model-coverage (query-coverage) and gene-coverage (target-coverage) calculated from a %(hmm-hits)s. Briefly, the program will remove all records from an %(hmm-source)s in the %(hmm-hits)s, then import a new %(hmm-hits)s table into the %(contigs-db)s that was filtered to your specifications.

For this, you first need to have %(anvi-run-hmms)s to ask HMMER to report a domain hits table by including `--domain-hits-table` flag in your command:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71 \
              --hmmer-output-dir path/to/dir
              --domain-hits-table
{{ codestop }}

After the command above, your HMM hits will be stored in your %(contigs-db)s as usual. However, with the availability of the domain hits table, you can filter out hits from your contigs database using thresholds for model or gene coverage of each hit i.e. you can filter out %(hmm-hits)s where the profile HMM and gene align well to each other.

For example, following the command above, the command below will remove %(hmm-hits)s from your %(contigs-db)s for profile HMMs that had less than 90%% coverage of the target genes:

{{ codestart }}
anvi-script-filter-hmm-hits-table -c %(contigs-db)s \
                                  --hmm-source Bacteria_71 \
                                  --domain-hits-table path/to/dir/hmm.domtable \
                                  --model-coverage 0.9
{{ codestop }}

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