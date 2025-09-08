This program allows you to import a %(state)s from a %(state-json)s.

You can run this program on a %(profile-db)s or %(pan-db)s like so: 

{{ codestart }}
anvi-import-state -p %(profile-db)s \
                  -s %(state-json)s \
                  -n MY_STATE
{{ codestop }}

This will import the state described in your %(state-json)s into your %(profile-db)s with the name `MY_STATE`. 
