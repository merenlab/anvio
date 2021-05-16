This program adds a new %(collection)s and %(bin)s to your %(pan-db)s or %(profile-db)s and %(contigs-db)s pair. This collection and bin will both contain all of your contigs. 

This way, you can perform collection and bin specfic operations without having to bin anything yourself. For example, running %(anvi-interactive)s in gene-mode requires you to specify a collection and bin (as is done [in the Infant Gut Tutorial](http://merenlab.org/tutorials/infant-gut/#the-gene-mode-studying-distribution-patterns-at-the-gene-level)). 

By default, the collection is named `DEFAULT` and the bin is named `EWVERYTHING`, but you can change these names with the `-C` and `-b` parameters respectively. 

Here is an example run on a %(profile-db)s: 

{{ codestart }}
anvi-script-add-default-collection -c %(contigs-db)s \ 
                                   -p %(profile-db)s \ 
                                   -C MY_COLLECTION \
                                   -b MY_BIN 
{{ codestop }}

Once this is run, your profile database will contain a collection called `MY_COLLECTION` with a single bin (called `MY_BIN`) which contains all of your contigs. 
