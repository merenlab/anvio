This program is for filtering %(hmm-hits)s from a %(contigs-db) using HMM alignment parameters such as query-coverage and target-coverage. Briefly, the program will remove the %(hmm-hits)s from a %(contigs-db) then import a new %(hmm-hits)s table into the %(contigs-db) that was filtered to your specifications. At the moment, this tool is only designed to work with `hmmsearch` with protein sequences which outputs the domain-table-output file from `hmmsearch` (--domtblout).

To list available %(hmm-source)s in a database:

{{ codestart }}
anvi-script-filter-hmm-hit-table -c %(contigs-db)s \
                                 --list-hmm-sources
{{ codestop }}

Then you can filter out hits using query or target coverage! Here's an example where we can filter out hmm_hits with a target-coverage less than 90%

{{ codestart }}
anvi-script-filter-hmm-hit-table -c %(contigs-db)s \
                                 --hmm-source %(hmm-source)s \
                                 --domtblout hmmsearch_domtable \
                                 --target-coverage 0.9
{{ codestop }}
