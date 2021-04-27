This program allows you to export a %(state)s from a %(pan-db)s or %(profile-db)s. The output of this is a %(state-json)s, which you can import into another anvi'o project with %(anvi-import-state)s. 

You can run this program on a %(profile-db)s or %(pan-db)s as follows: 

{{ codestart }}
anvi-export-state -s %(state)s \
                  -p %(profile-db)s  \
                  -o path/to/output
{{ codestop }}

To list the collections available in this database, you can run 

{{ codestart }}
anvi-export-state -p %(pan-db)s \
                  --list-states
{{ codestop }}
