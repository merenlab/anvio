This page describes general properties of anvi'o interactive displays and programs that offer anvi'o interactive artifacts.

## Terminology

Anvi'o uses a simple terminology to address various aspects of interactive displays it produces. The same terminology applies regardless of which program initiated the interface -- whether it is %(anvi-interactive)s, %(anvi-display-pan)s, or any other display program. Here is a quick-reference table:

| Term | Definition |
|------|------------|
| **Items** | The individual elements displayed in the central tree/dendrogram. Depending on the context, items can be contigs, gene clusters, genes, samples, or any other unit of data. Items are the *rows* of your data matrix. |
| **Layers** | The concentric rings of data that surround the central tree. Each layer represents a different data dimension -- for instance, a sample’s coverage, a genome in a pangenome, or a categorical annotation. Layers are the *columns* of your data matrix. |
| **View** | A specific representation of your data. Different views show different aspects of the same dataset (e.g., mean coverage vs. detection in a metagenomic context). Views are stored in the database and can be extended with `--additional-view`. |
| **Items order** | The tree or dendrogram at the center of the display that organizes items. It can be a phylogenetic/phylogenomic tree, a hierarchical clustering dendrogram, or a linear order. Users can provide custom orders through %(misc-data-items-order)s. |
| **Layer order** | The tree or dendrogram that organizes layers around the display. Users can extend available layer orders through %(misc-data-layer-orders)s. |
| **Bins** | User-defined selections (groups) of items. You can create bins by clicking on branches of the tree. Bins are central to tasks like genome binning in metagenomics or selecting gene clusters in pangenomics. |
| **Collection** | A named set of bins that can be stored in the database and recalled later. |
| **State** | A saved configuration of all visual settings -- layer order, colors, normalization, drawing type, etc. States are stored in the database and can be exported/imported with %(anvi-export-state)s and %(anvi-import-state)s. |
| **Additional data** | Extra information associated with items or layers that is not part of the core data matrix. Managed with %(anvi-import-misc-data)s, %(anvi-export-misc-data)s, and %(anvi-delete-misc-data)s. |

Now let’s see where each of these concepts lives in an actual display.

### Items and layers

![Items and layers in an anvi'o display](../../images/interactive_interface/terms_items_and_layers.gif){:.center-img .width-70}

**Items** are the units arranged around the center of the display. What they represent depends entirely on the context: contigs in a metagenomic analysis, gene clusters in a pangenome, genomes in a phylogenomic tree, or anything else. **Layers** are the concentric rings of data surrounding the items tree. Each layer shows a data value for every item, like coverage of a contig in a given sample, or presence/absence of a gene cluster in a given genome.

Together, items and layers form a data matrix: items are the rows, layers are the columns.

### Items order and layer order (dendrograms)

![The two dendrograms in an anvi'o display](../../images/interactive_interface/terms_dendrograms.gif){:.center-img .width-70}

The **items order** is the tree or dendrogram at the center of the display. It determines how items are organized: it can be a phylogenetic tree, a hierarchical clustering dendrogram, or a simple linear order. The **layer order** is the smaller tree on the outside of the display that organizes layers. Anvi'o computes layer orders automatically, but you can provide your own through %(misc-data-layer-orders)s.

### Views

![Two views of the same data in an anvi'o display](../../images/interactive_interface/terms_views.gif){:.center-img .width-70}

A **view** defines what numerical values are displayed in the layers. For instance, in a metagenomics dataset, the data that can be viewed are mean coverage, detection, or variability; and with pangenome, the options are presence/absence or frequencies of gene clusters. Views are stored in the database and can be switched using the **Data** dropdown in the settings panel. You can also add custom views with the `--additional-view` flag.

### Bins

![Bins in an anvi'o display](../../images/interactive_interface/terms_bins.gif){:.center-img .width-70}

**Bins** are user-defined groups of items. You create them by clicking on branches of the dendrogram. In a metagenomic context, bins typically represent metagenome-assembled genomes (MAGs). In a pangenomic context, they can represent groups of gene clusters (core genome, accessory genome, etc.). Bins are organized into **collections** that can be stored in the database and used by downstream programs like %(anvi-summarize)s.

### Additional data for items and layers

![Additional data in an anvi'o display](../../images/interactive_interface/terms_additional_data.gif){:.center-img .width-70}

Beyond the core view data, you can decorate your display with **additional data** for both items and layers. Additional data for items appears as extra layers on the outside of the display (e.g., taxonomy assignments, GC content, bin colors). Additional data for layers appears next to the layer order dendrogram (e.g., sample metadata like body site, or genome metadata like species name). These are managed through %(anvi-import-misc-data)s.

See [this article](https://merenlab.org/2017/12/11/additional-data-tables/) for a detailed guide on extending anvi'o displays with additional data tables.

## Programs that give interactive access

Several anvi'o programs produce an interactive display. The most commonly used ones are:

- %(anvi-interactive)s is the main interactive interface. It displays information from a %(profile-db)s and %(contigs-db)s, and supports manual binning, gene mode, collection mode, and a fully flexible manual mode for visualizing any data.

- %(anvi-display-pan)s displays pangenomic data from a %(pan-db)s. It shows gene cluster distributions across genomes and supports gene cluster inspection and binning.

- %(anvi-inspect)s provides a detailed, nucleotide-level view of a single contig across samples, including coverage, SNVs, and gene calls.

- %(anvi-display-contigs-stats)s shows summary statistics for contigs in a %(contigs-db)s.

- %(anvi-display-functions)s lets you browse the functional pool for a set of genomes or metagenomes.

- %(anvi-display-metabolism)s provides an interactive view of metabolism estimation data.

- %(anvi-display-structure)s lets you examine protein structures along with SCVs and SAAVs.

And a few others, including %(anvi-script-snvs-to-interactive)s.

## Artifacts that give interactive access

- %(contig-inspection)s shows detailed contig information (coverage, SNVs, gene calls).

- %(gene-cluster-inspection)s lets you examine individual gene clusters (alignments, genomic context, functional annotations).

## An overview of the display

The interactive interface has two major areas of interaction:

* The space for visualization in the middle area,
* The Settings panel on the left of the screen,

We will spend most of our time in the `Settings` panel.

### Settings panel

<img align="left" height="350px" src="../../images/interactive_interface/interactive_interface_tabs.png" style="margin-right: 20px;">

If closed, the settings panel can be opened by clicking on the little button on the left-middle part of your browser. When opened, you will see multiple tabs:

But before we start talking about these tabs, it is worthwhile to mention that at the bottom of the settings panel you will find a section with tiny controls that are available in all tabs:

Through these controls you can,

* __Create or refresh__ the display when necessary using the Draw button (or just press ***D***),

* __Load or Save__ the current state of the display using the buttons below the State section,

* __Download your display as an SVG file__ (or just press ***E***),

* __Center__ the display.

💡  Before talking about Tabs, here's a  *useful tip*: you can conveniently switch between tabs using your keypad. For instance, pressing `1` will navigate you to the Main tab, `2` to the Options tab, and so forth.

For more tips you can check out the [Tips + Tricks](#interactive-interface-tips--tricks) section at the end of this page.

OK. Let's talk about each tab you will find in the settings panel.


### Main Tab

This is one of the most frequently used tabs in the interface, and there are multiple sections in it (keeps growing over time, so things may be missing here).

![an anvi'o main tab](../../images/interactive_interface/interactive-settings-panel-main.png){:.center-img}

- #### Items subsection
  - Provides high level options for adjusting _drawing type_, _items order_ and _data type_.
  - Display part of the interface is where you can adjust individual layer attributes like _color_, display _type_, _height_  and _min/max_ values. Click + drag each layer to rearrange how layers are ordered. Or _edit attributes for multiple lpngayers_ as well.

- #### Layers subsection
  - __Change general settings for the tree__ (i.e., switching between circle or rectangular displays, changing tree radius or width), __and layers__ (i.e., editing layer margins, or activating custom layer margins).
  - __Load or save states__ to store all visual settings, or load a previously saved state.
  - Display usage is same as Items subsection.
  - Change the order of layers using automatically-generated or user-provided orders of layers using the Sample order combo box.
  - Customize individual samples information entries. Changes in this tab can be reflected to the current display without re-drawing the entire tree unless the sample order is changed.

- #### Legends subsection
  - This subsection enables users to easily change individual or batch legend colors for any of their additional data items.d

### Options tab

In the Options tab, just like in the Main tab, we have two sections: Items and Layers, providing a separation for these components.,

![an anvi'o options items](../../images/interactive_interface/interactive-settings-options-items.png){:.center-img}

* **Bins Selection subsection**. Allows us to customize the bins.
* **Cosmetics subsection**. _Margin and Background Opacity_ adjustments on the chart.
* **Dendrogram subsection**. _Radius_ and _Angle_ , and _Edge length normalization_ adjustments for the dendrogram.
* **Branch support subsection**. Settings for displaying _bootstrap values_ on the dendrogram.
* **Selections subsection**. To adjust _height_, _grid_ and/or _shade_ display, as well as selection _name_ settings.
* **Performance subsection**. Whether the SVG output is optimized for performance or granularity (very advanced stuff).

![an anvi'o options items](../../images/interactive_interface/interactive-settings-options-layers.png){:.center-img}

* **Tree/Dendogram subsection**. _Height and Edge length normalization_ adjustments on the Tree.
* **Label subsection**. Settings for maximum font size on the Layer.

Mastering these in the Options Tab will minimize the post-processing of your anvi'o figures for high-quality and good-looking publication ready images.

### Bins tab

Anvi'o allows you to create selections of items shown in the display (whether they are contigs, gene clusters, or any other type of data shown in the display). Bins tab allow you to maintain these selections. Any selection on the tree will be added to active bin in this tab (the state radio button next to a bin defines its activity). Through this tab you can,

![an anvi'o options items](../../images/interactive_interface/interactive-settings-bins-tabs.png){:.center-img}

- __Create or delete bins, set bin names, change the color of a given bin__, or sort bins based on their name, the number of units they carry, or completion and contamination estimates (completion / contamination estimates are only computed for genomic or metagenomic analyses).

- View the __number of selected units__ in a given bin, and see the __list of names in the selection__ by clicking the button that shows the number of units described in the bin.

- __Store a collection of bins__, or __load a previously stored collection.__

### Data tab

The data tab displays the value of items underneath the mouse pointer while the user browse the tree.

Displaying the numerical or categorical value of an item shown on the tree is not an easy task. We originally thought that displaying pop-up windows would solve it, but besides the great overhead, it often became a nuisance while browsing parts of the tree. We could show those pop-up displays only when use clicks on the tree, however click-behavior is much more appropriate to add or remove individual items from a bin, hence, it wasn’t the best solution either. So we came up with the ‘data tab’. You have a better idea? I am not surprised! We would love to try improve your experience: please enter an issue, and let’s discuss.

### Search tab

It does what the name suggests. Using this tab you can,

- __Build expressions to search items__ visualized in the main display.

- __Highlight matches__, and __append__ them to, or __remove__ them from the __selected bin__ in the Bins tab.

### News tab

The news panel provides information and external links tracking major Anvi'o releases and development updates.

### Notes tab

- The notes tab is a flexible, multipurpose space where users can,
- Store notes, comments, and any other stray items related to their project, in a feature-rich markdown environment.
- Display context, references, reproducibility instructions, and any other salient details for published figures.
- This section moved under Settings panel. `Settings > Notes`

![The Description panel in action](../../images/interactive_interface/interactive-settings-description-panel.png){:.center-img}

## Interactive interface tips + tricks

Here are some small conveniences that may help the interface serve you better (we are happy to expand these little tricks with your suggestions).

* You can zoom to a section of the display by making a rectangular selection of the area __while the pressing the shift button.__

* You can click an entire branch to add items into the selected bin, and remove them by __right-clicking__ a branch.

* If you click a branch __while pressing the `Command` or `CTRL` button__, it will create a new bin, and add the content of the selection into that bin.

* Tired of selecting items for binning one by one? __right-click__ on an item and select __Mark item as 'range start'__ to set an 'in point', then __right-click__ on another item and select __Add items in range to active bin__ or __Remove items in range from any bin__ to manipulate many items with few clicks. Nice!

* By pressing `1`,`2`,`3`,`4`,`5`,`6`,`7` and `8` you can go between Main, Options, Bins, Data, Notes, Search, News and Anvi'o tabs!

## Keyboard shortcuts

The interactive interface recognizes a handful of keyboard shortcuts to help speed up your workflow

- The `S` key toggles the Settings panel
- The `D` key triggers a redraw of your visualization
- The `T` key toggles showing the Title panel
- Keys `1` through `8` will toggle between tabs within the Settings panel, granted the Settings panel is currently shown.
- `CTRL`+`Z` and `CTRL`+`SHIFT`+`Z` will undo or redo bin actions, respectively.





