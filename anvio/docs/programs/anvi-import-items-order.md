This program, as one might think, allows you to import a %(misc-data-items-order-txt)s to describe a specific order of items stored in a %(profile-db)s, %(pan-db)s, or %(genes-db)s.

{{ codestart }}
anvi-import-items-order -p %(profile-db)s \
                        -i %(misc-data-items-order-txt)s
{{ codestop }}

It may also be nice to give it a good name, so that it's easy to find in the interface.

{{ codestart }}
anvi-import-items-order -p %(profile-db)s \
                        -i %(misc-data-items-order-txt)s \
                        --name ORDER_NAME
{{ codestop }}
