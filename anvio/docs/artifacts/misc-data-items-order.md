This artifact describes the **order of items in visualization tasks**.

In anvi'o, main display items (such as 'gene clusters' in a pan database, 'contigs' in a profile database, etc) can be ordered either by a NEWICK formatted tree (such as a phylogenetic tree or a hierarchical clustering dendrogram), or by an array (such as a flat list of item names).

When a NEWICK tree is used to order items, it will appear as the tree in the central section of the default anvi'o interactive interface. When a flat list of items are provided to order items, the central display where a tree appears will be blank and the displayed items will still be ordered according to the list. In order words, items order is to %(misc-data-items)s as %(misc-data-layer-orders)s is to %(misc-data-layers)s: a description not of the items themselves, but of what order they go in on the interface. 

Anvi'o programs such as %(anvi-pan-genome)s, %(anvi-merge)s, and %(anvi-profile)s automatically generate NEWICK-formatted items order if possible (i.e., if you have less than 20,000 items). When you run these programs, they will put this information into your resulting %(pan-db)s or %(profile-db)s. 

You can also export this information to give to a fellow anvi'o user or import this information if you have your own phylogenetic tree or desired order for your contigs.

You can use %(anvi-import-items-order)s to import specific orders for your items, or %(anvi-export-items-order)s to export this information.