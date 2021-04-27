This program, as implied by the name, is used to delete a %(collection)s from a %(profile-db)s. 

When you do this, you'll lose the collection forever, as well as the %(bin)ss within it. It is generally a good idea to export your binning effort into a %(collection-txt)s using %(anvi-export-collection)s before deleting it, just to be safe. 

To list available collections in a database, call 

{{ codestart }}
anvi-delete-collection -p %(profile-db)s \
                       --list-collections
{{ codestop }}

Then, you can easily delete a collection with the command

{{ codestart }}
anvi-delete-collection -p %(profile-db)s \
                       -C %(collection)s
{{ codestop }}
