This program removes %(hmm-hits)s from a %(contigs-db)s. This allows you to repopulate functional annotations with a different source or program, or simply remove data that may be cluttering the interface.

It is generally advisable to export your information before deletion as a precautionary measure. The HMM hits will appear in most displays, so if you have already run %(anvi-summarize)s, you should have this information preserved. 

To list available %(hmm-source)ss in a database, execute:

{{ codestart }}
anvi-delete-hmms -c %(contigs-db)s \
                 --list-hmm-sources
{{ codestop }}

You can then remove %(hmm-hits)s from a specific source using the command:

{{ codestart }}
anvi-delete-hmms -c %(contigs-db)s \
                 --hmm-source %(hmm-source)s 
{{ codestop }}
