This program lets you associate your layers with taxonomic information through a %(single-profile-db)s. 

This information is displayed in the interactive interface at the same place as %(misc-data-layers)s, which is point (4) on [this page](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology). 

If instead you want the layers to *represent* taxonomic ranks, then you'll want to take a look at [this tutorial on phylogenomics](http://merenlab.org/2017/06/07/phylogenomics/).

Usually, the layers describe separate samples. However, when working with only one sample, you may break up different aspects of that sample to be represented in each layer, hence why you might want to associate them with taxonomic information. 

To run this program, simply provide a %(layer-taxonomy-txt)s

{{ codestart }}
anvi-import-taxonomy-for-layers -p %(single-profile-db)s \
                                -i %(layer-taxonomy-txt)s 
{{ codestop }}

You also have the option to change the minimum abundance cut off using `--min-abundance`. The default value is 0.1 percent. 
