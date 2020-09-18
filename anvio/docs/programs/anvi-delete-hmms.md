This program, as implied by the name, is used to delete a %(hmm-hits)s from a %(contigs-db)s. This way, you can repopulate the function annotations with a different source or program or just delete data that's clogging up the interface.

It is generally a good idea to export your information before deleting it, just in case. The HMM hits will show up in most displays, so if you've already run %(anvi-summarize)s, you should be good. 

To list available %(hmm-source)ss in a database, call 

{{ codestart }}
anvi-delete-hmms -c %(contigs-db)s \
                 --list-hmm-sources
{{ codestop }}

Then, you can easily delete %(hmm-hits)s from a specific source with the command

{{ codestart }}
anvi-delete-hmms -c %(contigs-db)s \
                 --hmm-source %(hmm-source)s 
{{ codestop }}
