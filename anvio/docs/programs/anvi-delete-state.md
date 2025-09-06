This program removes a %(state)s from a %(pan-db)s or %(profile-db)s. This allows you to remove states that may be cluttering the state list in the interface. 

It is generally advisable to export your state before deletion as a precautionary measure (%(anvi-export-state)s).

To list available %(state)ss in a database, execute:

{{ codestart }}
anvi-delete-state -p %(pan-db)s \
                 --list-states
{{ codestop }}

You can then remove a %(state)s using the command:

{{ codestart }}
anvi-delete-state -p %(profile-db)s \
                  -s %(state)s 
{{ codestop }}
