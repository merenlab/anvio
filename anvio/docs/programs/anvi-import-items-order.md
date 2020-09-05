This program, as one might think, allows you to import a %(misc-data-items-order-txt)s, so that you can view it as a %(misc-data-items-order)s. 

You can import an item order into a %(profile-db)s, %(pan-db)s, or %(genes-db)s as follows: 

{{ codestart }}
anvi-import-items-order -p %(profile-db)s \
                        --i %(misc-data-items-order-txt)s
{{ codestop }}

It may also be nice to give it a good name, so that it's easy to find in the interface. 

{{ codestart }}
anvi-import-items-order -p %(profile-db)s \
                        --i %(misc-data-items-order-txt)s \
                        --name my_favorite_tree
{{ codestop }}
