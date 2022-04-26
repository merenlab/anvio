This program adds a 'default' %(collection)s and %(bin)s to your %(pan-db)s or %(profile-db)s and %(contigs-db)s that describes every item in your database.

This way, you can perform anvi'o tasks that require a collection or a bin even if you do not have a particular collection for your data, or all items in your database represent a meaningful bin (such as every contig in a %(contigs-db)s that represents a single genome).

As an example, see this program in action in the [Infant Gut Tutorial](http://merenlab.org/tutorials/infant-gut/#the-gene-mode-studying-distribution-patterns-at-the-gene-level) where it is used to run %(anvi-interactive)s on a genome in 'gene mode'.

Run in its simples form,

{{ codestart }}
anvi-script-add-default-collection -c %(contigs-db)s \
                                   -p %(profile-db)s
{{ codestop }}

the program will add a new collection into the profile database named `DEFAULT`, which will contain a single bin that describes all items in the database named `EVERYTHING`. You can set these default names to your liking using additional parameters:

{{ codestart }}
anvi-script-add-default-collection -c %(contigs-db)s \
                                   -p %(profile-db)s \
                                   -C MY_COLLECTION \
                                   -b MY_BIN
{{ codestop }}

Also see related programs, %(anvi-show-collections-and-bins)s and %(anvi-delete-collection)s.
