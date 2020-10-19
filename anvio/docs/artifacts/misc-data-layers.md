This is the section of your %(profile-db)s/%(pan-db)s that contains custom additional information about each of the layers of the interactive interface (usually displayed as the concentric circles). When you run %(anvi-interactive)s, this data will appear as additional graphs in line with your layers, similar to how the sample names are displayed at the top. 

As also defined in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology), this type of data will include information about each layer of the interface (usually representing your samples). This data is either numerical or categorical and can be imported into another database from a %(misc-data-layers-txt)s using %(anvi-import-misc-data)s. It is also displayed when you run %(anvi-show-misc-data)s and can be exported or deleted with %(anvi-export-misc-data)s and %(anvi-delete-misc-data)s respectively. 

If you would like to change the order that your layers are displayed, take a look at %(misc-data-layer-orders)s. Or, if you want to specifically import taxnomic information for your layers (if applicable), check out %(anvi-import-taxonomy-for-layers)s.

For example, this information could describe the salinity of a series of ocean samples, the continent your samples were taken in, or which of several collection methods was used. 
