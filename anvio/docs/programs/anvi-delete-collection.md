This program removes a %(collection)s from a %(profile-db)s or a %(pan-db)s. 

Executing this program will permanently delete the %(collection)s and all %(bin)ss it contains from the database without any additional confirmation prompts. If you are uncertain whether you may need a given collection for future analysis, it is recommended to export your binning effort into a %(collection-txt)s using %(anvi-export-collection)s before deletion, as a precautionary measure. 

To list available collections in a database, you can use the program %(anvi-show-collections-and-bins)s. Once you have identified the collection you wish to remove, you can execute the program on the target collection:

{{ codestart }}
anvi-delete-collection -p %(profile-db)s \
                       -C %(collection)s
{{ codestop }}
