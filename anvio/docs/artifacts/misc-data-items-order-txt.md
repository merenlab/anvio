This is a text file that **contains the information for a %(misc-data-items-order)s**, used for importing into and exporting this information from your anvi'o project.

## NEWICK order

If you intend to import a tree order, the contents of your file should look something like this (but probably much more complicated depending on the number of items in your anvi'o database): 

```
(contig_4, ((contig_1, contig_2), contig_3))
```

When a NEWICK order is imported into an anvi'o project, the contigs will be displayed in the order `contig_4, contig_1, contig_2, contig_3`, and the following tree will be generated in the interface:

```
    contig_4    contig_1    contig_2    contig_3
        |           |           |           |
        |           -------------           |
        |                 |                 |
        |                 -------------------
        |                           |
        -----------------------------
                    |
                    |
```

## LIST order

Alternative to the NEWICK order, you can provide a list of items in flat form. For instance, if you want to order your items this way, your text file should look like the following, where each line contains a single item name in your database:

```
contig_4
contig_1
contig_2
contig_3
```

{:.warning}
After importing an order into a database, you may need to specifically select that order in the interactive interface through the "Item orders" dropbox and re-draw your display to change the default order.
