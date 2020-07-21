This program, as one might think, allows you to export a %(collection)s. This allows you to take your binning results elsewhere (including into another Anvi'o project with the command %(anvi-import-collection)s). 

You can run this program on a %(profile-db)s or %(pan-db)s as follows: 

{{ codestart }}
anvi-export-collection -C my_favorite_collection \
                        -p %(profile-db)s 
{{ codestop }}

This will give you a %(collection-txt)s file that describes the collection `my_favorite_collection`. 

To list the collections available in this database, you can run 

{{ codestart }}
anvi-export-collection -p %(pan-db)s \
                        --list-colllections
{{ codestop }}

You can also add the flag `--include-unbinned` to have all unbinned contigs in the database show up at the end of your %(collection-txt)s file in a bin titled `UNBINNED`. 
