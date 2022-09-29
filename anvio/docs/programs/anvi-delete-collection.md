This program, as implied by the name, will delete a %(collection)s from a %(profile-db)s or a %(pan-db)s. 

Using this program will delete the %(collection)s and every %(bin)s it describes from the database forever and without any additional warning. If you are not sure whether you may need a given collection later, it may be a good idea to export your binning effort into a %(collection-txt)s using %(anvi-export-collection)s before deleting it, just to be safe. 

To list available collections in a database you can use the program %(anvi-show-collections-and-bins)s. When you know which collection you wish to remove, you can run the program on the target collection name:

{{ codestart }}
anvi-delete-collection -p %(profile-db)s \
                       -C %(collection)s
{{ codestop }}
