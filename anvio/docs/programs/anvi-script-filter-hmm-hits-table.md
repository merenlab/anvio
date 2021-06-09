This program is for filtering a %(hmm-source)s from a %(hmm-hits)s in a %(contigs-db)s using HMM alignment parameters such as query-coverage and target-coverage. Briefly, the program will remove all records from an %(hmm-source)s in the %(hmm-hits)s then import a new %(hmm-hits)s table into the %(contigs-db)s that was filtered to your specifications. At the moment, this tool is only designed to work with `hmmsearch` with protein sequences. The `--domtblout` can be produced from running %(anvi-run-hmms)s with your %(hmm-source)s of interest and using the `--domtblout` parameter AND  `hmmsearch` as the program.

For this, you first need to ask %(anvi-run-hmms)s to ask HMMER to report a domain hits table by including `--hmm-domain-tblout-path` flag in your command:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71 \
              --hmm-domain-tblout-path DOMTABLE.txt
{{ codestop }}

At the end of this run, your HMM hits will be stored in your contigs database as usual. But with the availability of the domain hits table from this run, you can filter out hits from your contigs database using thresholds for query or target coverage of each hit.

For instance following the command above, the command below will remove HMM hits from your contigs database for genes that had less than 90%% coverage of the target:

{{ codestart }}
anvi-script-filter-hmm-hits-table -c %(contigs-db)s \
                                  --hmm-source Bacteria_71 \
                                  --domtblout DOMTABLE.txt \
                                  --target-coverage 0.9
{{ codestop }}
