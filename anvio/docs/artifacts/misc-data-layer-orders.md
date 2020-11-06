This is the section of your %(profile-db)s/%(pan-db)s that contains custom additional information about the order that your layers are displayed in and the tree that relates them to each other . When you run %(anvi-interactive)s, this data will determine what order the concentric circles are displayed in, as well as the tree that appears above the %(misc-data-layers)s graphs.

As also defined in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology), this type of data will include information about how your layers are related to each other and determines their order. This data is stored as a tree that is displayed in the top-right. 

This data can be imported from a %(misc-data-layer-orders-txt)s artifact with %(anvi-import-misc-data)s.  It is also displayed when you run %(anvi-show-misc-data)s and can be exported or deleted with %(anvi-export-misc-data)s and %(anvi-delete-misc-data)s respectively. 

For example, you could use this tree to indicate and group together samples that came from the same geographic location, samples that came from the same donor, samples of the same type,  samples collected with the same collection method, and so on. 

This is also used to import the taxonomy information at the end of [the pangenomics and phylogenomics workflow](http://merenlab.org/2017/06/07/phylogenomics/#pangenomic--phylogenomics). 
