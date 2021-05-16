This program tells you about the %(collection)ss within a %(profile-db)s or %(pan-db)s. 

Just run it like so 

{{ codestart }}
anvi-show-collections-and-bins -p %(profile-db)s 
{{ codestop }}

and Anvi'o will output to your console the following information for each of the %(collection)ss in the database: 

* The name and ID of the collection
* The number of %(bin)ss within the collection, and each of their names
* The number of splits contained within those bins 
