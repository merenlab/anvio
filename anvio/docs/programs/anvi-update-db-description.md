This program allows you to update the description of any anvi'o database with the push of a button (and the writing of an updated description). 

This descirption helps make UIs a little prettier by showing up when you run programs like %(anvi-interactive)s and %(anvi-summarize)s. 

Simply write out the description that you would prefer in a plain text file (with markdown syntax) and use this program to update the description of any %(pan-db)s, %(profile-db)s, %(contigs-db)s, or %(genomes-storage-db)s: 

{{ codestart }}
anvi-update-db-description --description my_description.txt \
                           %(contigs-db)s
{{ codestop }}

