Anvi-interactive opens the Anvi'o interactive interface, which is one of the most sophisticated parts of Anvi'o. 

The most widely-known view gives you beautiful concentric circles of data, but the interface has many forms and vast functionality, from manual metagenomic binning (check out anvi-interactive in a metagenomic workflow [here](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-interactive)) to examining various stats about your data. However, in case you don't like circles, you can also display your data in a rectangle (as seen [here](http://merenlab.org/tutorials/interactive-interface/#lets-go-all-corners)). 

In fact, the interface has many of its own blog posts, including a pretty comprehensive introductory tutorial [here](http://merenlab.org/tutorials/interactive-interface/) and a breakdown of its data types  [here](http://merenlab.org/2016/02/27/the-anvio-interactive-interface/). 

Here, we'll go through *some* things that the Anvi'o interactive interface is capable of through this program. 

## Running anvi-interactive on a profile database

One of the simplest ways to run the interactive interface (espeically useful for manual binning) is just providing a profile database and contigs database:

{{ codestart }}
anvi-interactive -p %(profile-db)s \ 
            -c %(contigs-db)s
{{ codestop }}

For the central tree to display correctly, you'll need to have run hierarchical clustering at some point while making your profile database (either during %(anvi-merge)s, or, if this is a %(single-profile-db)s, while running %(anvi-profile)s). 

You'll get lovely concentric circles, each filled with data that was contained in your databases and that you are now free to interact with. See the page for the %(interface)s for more information. 

### Running on a specific collection 

You can also run anvi-interactive on a specific collection. When doing this, Anvi'o will calculate various information about each of your bins, and display the interface. Each leaf of your plot will not represent a contig, but a bin within your collection. 

{{ codestart }}
anvi-interactive -p %(profile-db)s \ 
            -c %(contigs-db)s \
            -C %(collection)s
{{ codestop }}

Since clustering is done here, you can also customize the linkage method and distance metric if desired.

See the note on this mode in the [metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-interactive) for more information. 

If instead you want to run the standard anvi'o interface, but only looking at contigs within a specific collection, use the tag `--colection-autoload`. 

### Looking at *genes* instead of bins or contigs

You can also start the interface in "gene mode," in which each leaf of the central tree is not a split or contig (or bin like in "collection mode"). A lot of the same functionality is availble, including looking at detection, coverage, inspection, and sequence functionality. However, you cannot store states or collections in this mode. 

## Manual Inputs: I want to provide my own (non-Anvi'o) data!

You can do that with the flag `--manual-mode` and then by providing the following types of files: 

- a %(fasta)s file
- a tab-delimated view data file
- a NEWICK formatted tree
- a flat file containing the order of the items you want to display

When doing this kind of run, if you don't provide a profile database, Anvi'o will create an empty one for you. 

## Visualization Settings

In Anvi'o, the visualization settings at a given time are called a %(state)s. 

To open the interface in a specific state, you can use the `--state-autoload` flag or by importing a state using %(anvi-import-state)s. 

You can also customize various aspects of the interactive interface. For example, you can change the preselected view, title, and taxnomic level displayed (for example, showing the class name instead of the genus name). You can also hide outlier single nucleotide variations or open only a specific collection. 

### Adding additional information 

You can add any additional layers of your choice using the parameter `--additional-layers` and providing a file containing the information you want to appear as another layer. You could also choose to split non-singlecopy gene HMM hits into their own layer with the `--split-hmm-layers` parameter. 

If you want to add an entirely new view to the interface, you can do that too, as long as you provide a file containing all split names and their associated values. For more information, see the parameter `--additional-view`. 

You can also provide the manual inputs even if you're using an Anvi'o database. For example, if you provide your own NEWICK formatted tree and it will be displayed instead of the one in your database. 

## Other things 

### Viewing your data

You can use this program to look at the available information in your databases, which is very convenient. For example, you can view all of the available

- views (using `--show-views`)
- states (using `--show-states`)
-  collections (using `--list-collections`)

### Note for power users 

You can also configure the server to your heart's content, skip function call initizations, and change any of the output paths. 
