This artifact describes the **order of the branches for a visualization**.

This is the tree in the central section of the default anvi'o interactive interface, or the order that your contigs (or genes or bins, depending on the mode) display in when in the circles. In order words, it is to %(misc-data-items)s as %(misc-data-layer-orders)s is to %(misc-data-layers)s: a description not of the items themselves, but of what order they go in on the interface. 

Most often, this is a phylogenetic tree, but it can also just describe what order to display the branches or items in. 

As of now, this is an provided by programs that generate a tree of this kind, including %(anvi-pan-genome)s, %(anvi-merge)s, and %(anvi-profile)s. When you run these programs, they will put this information into your resulting %(pan-db)s or %(profile-db)s. 

However, you can also export this information to give to a fellow Anvi'o user or import this information if you have your own phylogenetic tree or desired order for your contigs. This is especially useful if you want to perform manual binning. 

To import your own order for your items, use %(anvi-import-items-order)s. To export this information, use %(anvi-export-items-order)s. 


