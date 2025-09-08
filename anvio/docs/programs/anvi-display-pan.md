This program **displays the contents of a %(pan-db)s in the [anvi'o interactive interface](http://merenlab.org/2016/02/27/the-anvio-interactive-interface//#using-the-anvio-interactive-interface), much like %(anvi-interactive)s.**

Like you can see in the [pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#displaying-the-pan-genome), this opens a window of the interactive interface where each item is a gene cluster and each layer represents one of your genomes. 

### A general run 

You can run it with only two parameters: 

{{ codestart }}
anvi-display-pan -p %(pan-db)s \
                 -g %(genomes-storage-db)s 
{{ codestop }}

There are several default layer orders to choose from, including organizing based on gene cluster presence/absense or gene cluster frequency. These will both group your core gene clusters and singletons separately. 

Beyond that, there are many different settings you can change in the side panel of the interface and you can import various additional data (primarily with the program %(anvi-import-misc-data)s). Once you're happy with the data displayed in the interface (and the prettiness of that data), you can  save those preferences in a %(state)s. 

### I want MORE data displayed 

There are several other data types you can additionally choose to display in this program. Namely, you can add:

- a title (very fancy I know) using `--title` 
- a NEWICK formatted tree (or import it as a %(misc-data-items-order-txt)s with %(anvi-import-items-order)s or as a %(misc-data-layer-orders)s with %(anvi-import-misc-data)s). 
- view data in a tab-delimited file
- an additional view (provide this in a tab-delimited matrix where each column corresponds to a sample and each row corresponds to a gene cluster)
- an additional layer in the form of a %(misc-data-layers-txt)s (or import it into your %(pan-db)s with %(anvi-import-misc-data)s

### How to minimize mouse clicks 

Wondering how to autoload specific aspects of the interface? You're in the right place. 

You have the option to specify quite a few aspects of the interface through the parameters to save you those sweet mouse clicks. 

- You can specify which view to start the interface with. Check which views are available with `--list-views`. 
- You can load a specific %(state)s (either a previous state or a state that you've imported with %(anvi-import-state)s). Check which states are available with the flag `--list-states`. 
- You can also load a specific %(collection)s with `--collection-autoload`. To check which collections are availible, use `--list-collections`. 

### Other parameters 

You can also skip processes like intializing functions or automatically ordering your items to save time, as well as configure the server to your heart's content. 
