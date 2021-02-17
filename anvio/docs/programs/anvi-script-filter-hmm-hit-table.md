This program is for filtering a %(hmm-source)s from a %(hmm-hits)s in a %(contigs-db)s using HMM alignment parameters such as query-coverage and target-coverage. Briefly, the program will remove all records from an %(hmm-source)s in the %(hmm-hits)s then import a new %(hmm-hits)s table into the %(contigs-db)s that was filtered to your specifications. At the moment, this tool is only designed to work with `hmmsearch` with protein sequences. The `--domtblout` can be produced from running `anvi-run-hmms` with your %(hmm-source)s of interest and using the `--domtblout` parameter AND  `hmmsearch` as the program.

To list available %(hmm-source)s in a database:

{{ codestart }}
anvi-script-filter-hmm-hit-table -c %(contigs-db)s \
                                 --list-hmm-sources
{{ codestop }}

Make `--domtblout` with `anvi-run-hmms`

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71 \
              --just-do-it \
              --hmmer-program hmmsearch \
              --hmm-domain-tblout-path hmmsearch
{{ codestop }}

Then you can filter out hits using query or target coverage! Here's an example where we can filter out hmm_hits with a target-coverage less than 90%

{{ codestart }}
anvi-script-filter-hmm-hit-table -c %(contigs-db)s \
                                 --hmm-source %(hmm-source)s \
                                 --domtblout hmmsearch_domtable \
                                 --target-coverage 0.9
{{ codestop }}
