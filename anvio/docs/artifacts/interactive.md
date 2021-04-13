This page describes general properties of anvi'o interactive displays and programs that offer anvi'o interactive artifacts.

## Terminology

Anvi'o uses a simple terminology to address various aspects of interactive displays it produces, such as items, layers, views, orders, and so on. The purpose of this section is to provide some insights into these terminology using the figure below:

![an anvi'o display](../../images/interactive_interface/anvio_display_template.png){:.center-img}

Even though the figure is a product of %(anvi-display-pan)s, the general terminology does not change across different interfaces, including the default visualizations of %(anvi-interactive)s. Here are the descriptions of numbered areas in the figure:

* The tree denoted by **(1)** shows the organization of each `item`. Items could be contigs, gene clusters, bins, genes, or anything else depending on which mode the anvi'o interactive interface was initiated. The structure that orders items and denoted by **(1)** in the figure can be a phylogenetic or phylogenomic tree, or a dendrogram produced by a hierarchical clustering algorithm. In addition, there may be nothing there, if the user has requested or set a linear items order through %(misc-data-items-order)s.
* Each concentric circle underneath the number **(2)** is called a `layer` and the data shown for items and layers as a whole is called a `view`. A **layer** can be a genome, a metagenome, or anything else depending on which mode the anvi'o interactive was initiated. The **view** is like a data table where a datum is set for each **item** in each **layer**. The view data is typically computed by anvi’o and stored in pan databases by %(anvi-pan-genome)s or profile databases by %(anvi-profile)s. The user add another view to the relevant combo box in the interface by providing a TAB-delimited file to %(anvi-interactive)s through the command line argument `--additional-view`, or add new layers to extend these vies with additional data through %(misc-data-items)s.
* The tree denoted by **(3)** shows a specific ordering of layers. Anvi'o will compute various layer orders automatically based on available **view** depending on the analysis or visualization mode, and users can extend available **layer orders** through %(misc-data-layer-orders)s.
* What is shown by **(4)** is the additional data for layers. the user can extend this section with additional information on layers using the %(misc-data-layers)s.

The orchestrated use of %(anvi-import-misc-data)s, %(anvi-export-misc-data)s, and %(anvi-delete-misc-data)s provides a powerful framework to decorate items or layers in a display and enhance visualization of complex data. Please take a look at the following article on how to extend anvi'o displays:

* [https://merenlab.org/2017/12/11/additional-data-tables/](https://merenlab.org/2017/12/11/additional-data-tables/)

## Programs that give interactive access

If you're new to the anvi'o interactive interface, you'll probably want to check out [this tutorial for beginners](http://merenlab.org/tutorials/interactive-interface/) or the other resources on the  %(anvi-interactive)s page.

However, there are more interfaces available in anvi'o than just that one, so let's list them out:

- %(anvi-display-structure)s lets you examine specific protein structures, along with SCV and SAAVs within it. (It even has [its own software page.](http://merenlab.org/software/anvio-structure/). It's kind of a big deal.)

- %(anvi-display-contigs-stats)s shows you various stats about the contigs within a %(contigs-db)s, such as their hmm-hits, lengths, N and L statistics, and so on.

- %(anvi-display-functions)s lets you quickly browse the functional pool for a given set of genomes or metagenomes.

- %(anvi-display-metabolism)s is still under development but will allow you to interactively view metabolism estimation data using %(anvi-estimate-metabolism)s under the hood.

- %(anvi-display-pan)s displays information about the gene clusters that are stored in a %(pan-db)s. It lets you easily view your core and accessory genes, and can even be turned into a metapangenome through importing additional data tables.

- %(anvi-inspect)s lets you look at a single split across your samples, as well as the genes identified within it. This interface can also be opened from the %(anvi-interactive)s interface by asking for details about a specific split.

- %(anvi-interactive)s displays the information in a %(profile-db)s. It lets you view the distribution of your contigs across your samples, manually bin metagenomic data into MAGSs (and refine those bins with %(anvi-refine)s), and much more. You can also use this to look at your genes instead of your contigs or [examine the genomes after a phylogenomic anlysis](http://merenlab.org/2017/06/07/phylogenomics/). Just look at that program page for a glimpse of this program's amazingness.

- %(anvi-script-snvs-to-interactive)s lets you view a comprehensive summary of the SNVs, SCVs, and SAAVs within your contigs.

<!-- wireframe of documentation layout 
     each subtopic should link to relevant tutorials/workflows that highlight focused portion of the interactive interface

     pull in content from 
     https://merenlab.org/2016/02/27/the-anvio-interactive-interface/
     https://merenlab.org/tutorials/interactive-interface/

     https://merenlab.org/software/anvio/help/main/programs/anvi-inspect/ 
-->
- Settings Panel 
    - Main
    - Layers
    - Bins
    - Legends
    - Search 
    - hamburger dropdown menu 
    - Bottom Left
        - State stuff 
        - zoom, download, upload 

- Mouse Context Menu
    - layers
    - tree 

- Mouse Panel 

- Inspect Page 
    - accessed from context menu and terminal 

- Glossary of keyboard shortcuts 
- Glossry of helpful hints

# An overview of the display

The interactive interface has two major areas of interaction: the space for visualization on the right, and the settings panel. The settings panel gives access to various controls to work with the data visualized, and improve the presentataion of it.

# Settings panel

At the bottom of the settings panel there is a section with tiny controls that are available in all tabs. Through these controls you can,
- __Create or refresh__ the display when necessary using the draw button (some changes require you to do that),
- __Zoom in, zoom out, and center__ the display.
- __Download your display as an SVG file.__

## Main tab 
- Display
    - content 
- Layers
    - content
## Layers tab 
Through the layers tab you can,

- __Change general settings for the tree__ (i.e., switching between circle or rectengular displays, changing tree radius or width), __and layers__ (i.e., editing layer margins, or activating custom layer margins).
- __Load or save states__ to store all visual settings, or load a previously saved state.
- __Customize individual__ layers by switching between different __display modes__ depending on the layer type (i.e., ‘text’ or ‘color’ mode for categorical layers, or ‘bar’ or ‘intensity’ mode for numerical layers), __set normalization__ (i.e., ‘square-root’, or ‘log’ normalization), __minimum, and maximum cutoff__ values for numerical layers, or set __layer height__, and __layer margin__ (i.e., its distance from the previous layer).
- Use the __multi-selector__ at the bottom to change settings for multiple layers at once.

<!-- blog post links to 'working with SVGs blog post, good info, should incorporate? -->
## Samples tab 
Samples tab is for the additional data you provide the interface through a samples database (see samples order and samples infomration sections above). Through this layer you can,
- __Change the order__ of layers using automatically-generated or user-provided orders of layers using the Sample order combo box,
- __Customize individual samples information entries.__
Changes in this tab can be reflected to the current display without re-drawing the entire tree unless the sample order is changed.
## Bins tab 
Anvi’o allows you to create selections of items shown in the display (whether they are metagenomic contigs, 16S rRNA tags, or any other type of information). Bins tab allow you to maintain these selections. Any selection on the tree will be added to active bin in this tab (the state radio button next to a bin defines its activity). Through this tab you can,
- __Create or delete bins, set bin names, change the color of a given bin__, or sort bins based on their name, the number of units they carry, or completion and contamination estimates (completion / contamination estimates are only computed for genomic or metagenomic analyses).
- View the __number of selected units__ in a given bin, and see the __list of names in the selection__ by clicking the button that shows the number of units described in the bin.
- __Store a collection of bins__, or __load a previously stored collection.__
## Legends tab 
- content 
## Search tab 
It does what the name suggests. Using this tab you can,
- __Build expressions to search items__ visualized in the main display.
- __Highlight matches__, and __append__ them to, or __remove__ them from the __selected bin__ in the Bins tab.




# The right panels
## Mouse panel
The mouse panel displays the value of items underneath the mouse pointer while the user browse the tree.

Displaying the numerical or categorical value of an item shown on the tree is not an easy task. We originally thought that displaying pop-up windows would solve it, but besides the great overhead, it often became a nuisance while browsing parts of the tree. We could show those pop-up displays only when use clicks on the tree, however click-behavior is much more appropriate to add or remove individual items from a bin, hence, it wasn’t the best solution either. So we came up with the ‘mouse panel’. You have a better idea? I am not surprised! We would love to try improve your experience: please enter an issue, and let’s discuss.

## News panel 
- some info here
## Description panel 
- some info here




