Essentially, a collection **is a group of %(bin)ss**.

You can use collections to represent all kinds of things. The default collection (that you'll get by running %(anvi-script-add-default-collection)s) is called DEFAULT and contains all of your contigs. However, you can put any group of bins into their own collection, and use that to limit what you're analyzing downstream. A ton of anvi'o programs are able to take in a bin or a collection so that you don't have to analyze your entire contigs-db when you just want to look at one section. 

To look at the collections contained within an Anvi'o database, just run %(anvi-show-collections-and-bins)s. Or, to look at the completion estimates for all bins within a collection, use %(anvi-estimate-genome-completeness)s.

To view the content within a collection, use %(anvi-summarize)s.

### Examples of collections in action

If your bins represent MAGs, you could use a collection to group related MAGs together. For example, you could group together genomes from the same population, or from the same taxonomic order. Or your could go even wider and have a collection for all of the Archaea genomes in your sample. You're not limited by taxonomy either. You could go wild and have a collection for all of your bins that you suspect are prophages. 

