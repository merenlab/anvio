This is the section of your %(profile-db)s/%(pan-db)s that contains custom additional information about each of the items in the central section of the interactive interface. When you run %(anvi-interactive)s, this data will appear as additional concentric circles. 

As also defined in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology), this type of data will include information about each item (whether that's a contig, gene, or bin). This data is either numerical or categorical and can be imported into another database from a %(misc-data-items-txt)s using %(anvi-import-misc-data)s. It is also displayed when you run %(anvi-show-misc-data)s and can be exported or deleted with %(anvi-export-misc-data)s and %(anvi-delete-misc-data)s respectively. 

To change the order that the items are displayed in, take a look at %(anvi-import-items-order)s.

For example, this information could describe whether or not each bin reached a certain completion threshold, the e-score of the function annotation on each gene, or different categories that the total length of a contig could fall into (1-1.5 kb, 1.5-2 kb, 2-2.5 kb, and so on). 
