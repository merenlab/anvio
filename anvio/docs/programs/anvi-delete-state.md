This program, as implied by the name, is used to delete a %(state)s from a %(pan-db)s or %(profile-db)s. This way, you can remove states that are clogging up the state list in the interface. 

It is generally a good idea to export your state before deleting it, just in case ((anvi-export-state)s).

To list available %(state)ss in a database, call 

{{ codestart }}
anvi-delete-state -p %(pan-db)s \
                 --list-states
{{ codestop }}

Then, you can easily delete a %(state)s with the command

{{ codestart }}
anvi-delete-hmms -p %(profile-db)s \
                 -s %(state)s 
{{ codestop }}
