%(anvi-draw-kegg-pathways)s draws %(kegg-pathway-map)s files incorporating data from anvi'o databases. The visualization of user data in the context of KEGG's curated biochemical pathways can reveal patterns in metabolism.

## Setup

There are hundreds of pathway maps, listed and categorized [here](https://www.genome.jp/kegg/pathway.html). %(anvi-setup-kegg-data)s downloads the maps that have corresponding [files](https://www.kegg.jp/kegg/xml/) that allow elements of the map to be modified. Make sure to use the `--kegg-snapshot` option with the newest snapshot of %(kegg-data)s that includes pathway map files, `v2024-08-30`.

{{ codestart }}
anvi-setup-kegg-data --kegg-snapshot v2024-08-30
{{ codestop }}

### Download newest available files

Alternatively, KEGG data including maps can be set up not from a snapshot but by downloading the newest files available from KEGG with the `-D` flag. (`-T` downloads faster by using multiple threads, here 5.)

{{ codestart }}
anvi-setup-kegg-data -D -T 5
{{ codestop }}

### Install in non-default location

At the moment, KEGG data that includes maps does _not_ include stray KOs (see %(anvi-setup-kegg-data)s) due to changes in the available files. To preserve KEGG data that is already set up, for this reason or another, the new snapshot or download can be placed in a non-default location with `--kegg-data-dir`.

{{ codestart }}
anvi-setup-kegg-data --kegg-snapshot v2024-08-30 --kegg-data-dir path/to/other/directory
{{ codestop }}

`anvi-draw-kegg-pathways` requires the `--kegg-dir` argument to seek KEGG data in a non-default location.

## Pathway selection

By default, this program draws the maps that contain data of interest, e.g., KO gene annotations in a %(contigs-db)s. For _all_ maps in %(kegg-data)s to be drawn, including those that don't contain data of interest, use the flag, `--draw-bare-maps`.

The optional argument, `--pathway-numbers`, limits the output to maps of interest. A single ID number can be provided, e.g., '00010' for 'Glycolysis / Gluconeogenesis', or multiple numbers can be listed, e.g., '00010 00020'. Regular expressions can also be provided, e.g., '011.. 012..', where '.' represents any character: here the set of numbers given by '011..' corresponds to 'global' maps and '012..' to 'overview' maps.

## KO occurrence

Gene sequences in anvi'o databases can be annotated with KEGG Orthologs (KOs): see %(anvi-run-kegg-kofams)s. A KO indicates functional capabilities of the gene product. KO data from one or more contigs databases or a pan database can be mapped using the `--ko` flag, enabling investigation of the metabolic capabilities of individual organisms or multiple organisms, including community samples. Reactions associated with KOs are colored on the pathway maps.

### Single contigs database

Here is the basic command to draw KO data from a single %(contigs-db)s.

{{ codestart }}
anvi-draw-kegg-pathways --contigs-dbs %(contigs-db)s -o output_dir --ko
{{ codestop }}

Here are three maps drawn with this command from a bacterial genomic contigs database. The map in the upper left, "00010 Glycolysis / Gluconeogenesis," is a "standard" map, in which boxes are associated with reactions and sets of KOs. The map in the upper right, "01200 Carbon metabolism," is a metabolic "overview" map; these have numerical IDs in the range '012XX'. Reaction arrows are associated with KOs and are colored and widened given KO data. The bottom map, "01100 Metabolic pathways," is a "global" metabolic map; these have numerical IDs in the range '011XX'. Reaction lines are associated with KOs are colored given KO data. In all maps, circles are colored if the compound they represent is involved in reactions that are also colored.

![Three maps showing KOs from a single contigs database](../../images/anvi-draw-kegg-pathways/kos_single_contigs_db.png)

#### Set color

The default color can be changed with the option, `--set-color`, which takes either a color hex code (e.g., "#FF0000" for red) or "original" for the original color scheme of the reference map. Global maps are especially colorful, with reactions varying in color across the map. It is necessary to enclose a color hex code argument value in quotation marks, as "#" otherwise causes the rest of the command to be ignored as a comment.

{{ codestart }}
anvi-draw-kegg-pathways --contigs-dbs %(contigs-db)s -o output_dir --pathway-numbers 00010 --ko --set-color "#2986cc"
{{ codestop }}

![Change color to blue](../../images/anvi-draw-kegg-pathways/kos_color_blue.png)

{{ codestart }}
anvi-draw-kegg-pathways --contigs-dbs %(contigs-db)s -o output_dir --pathway-numbers 00010 01100 01200 --ko --set-color original
{{ codestop }}

![Use original color scheme](../../images/anvi-draw-kegg-pathways/kos_color_original.png)

### Multiple contigs databases

The KO content of multiple contigs databases can be compared. Database file paths can be provided directly on the command line or in an %(external-genomes)s text file.

{{ codestart }}
anvi-draw-kegg-pathways --contigs-dbs %(contigs-db)s_1 %(contigs-db)s_2 ... %(contigs-db)s_N -o output_dir --ko
{{ codestop }}

{{ codestart }}
anvi-draw-kegg-pathways --external-genomes %(external-genomes)s -o output_dir --ko
{{ codestop }}

The images in this section show data from contigs databases of bacterial genomes from different strains of the same species.

#### Color by database

When comparing a small number of databases (realistically, two or three), reactions can be colored by their occurrence across databases, with each color representing a different database or combination of databases. A colorbar key is drawn in a separate file in the output directory, 'colorbar.pdf'. Compound circles are imparted the color of the associated reaction found in the greatest number of databases.

![Three maps showing KOs from three contigs databases](../../images/anvi-draw-kegg-pathways/kos_three_contigs_dbs.png)

#### Color by count

When comparing a larger number of databases, it makes more sense to color reactions by the number of databases in which they occur using a sequential colormap than by database or combination of databases using a qualitative colormap. By default, coloring explicitly by database automatically applies to three or fewer databases, whereas coloring by database count applies to four or more databases. The user can override this default with the argument, `--colormap-scheme`, which accepts the values "by_database" and "by_count" (e.g., the user may have three databases but wish to color reactions by database count, and so would specify `--colormap-scheme by_count`).

![Three maps showing KOs from six contigs databases](../../images/anvi-draw-kegg-pathways/kos_six_contigs_dbs.png)

#### Reverse colormap

Changing the colormap can draw attention to different information on maps. When coloring by count, the default sequential colormap, "plasma_r", goes from dark to light colors: reactions shared among all of the contigs databases are assigned the darkest color, and reactions unique to a single database are assigned the lightest color. The colormap can be reversed to accentuate unshared reactions in the darkest colors and shared reactions in the lightest colors. Reversing the colormap is accomplished with the option, `--colormap plasma 0.1 0.9`. Note that Matplotlib colormap names differing by "_r" are reversed versions. The second and third numerical `--colormap` values are not mandatory, but can be provided to trim a fraction of the colormap from each end to eliminate the lightest and darkest colors. By default with "plasma_r", the limits are "0.1 0.9". Just changing the colormap (e.g., `--colormap plasma`) removes the limits (i.e., "0.0 1.0"), so exactly reversing the colormap requires that the same limits be specified. The flag, `--reverse-overlay`, should also be used to reverse the default drawing order such that unshared reactions appear above shared reactions, which is especially important in cluttered global maps.

{{ codestart }}
anvi-draw-kegg-pathways --external-genomes %(external-genomes)s -o output_dir --ko --colormap plasma 0.1 0.9 --reverse-overlay
{{ codestop }}

![Emphasize unshared reactions with reversed coloring](../../images/anvi-draw-kegg-pathways/kos_reverse_colormap.png)

#### Showing individual database maps

Coloring by count obviously obscures the individual databases containing different reactions. However, options are provided to enable investigation of the distribution of reactions across databases. Standalone map files showing the presence/absence of reactions in individual contigs databases can be drawn with the flag, `--draw-individual-files`. Maps for individual databases can also be drawn alongside the map of all databases, facilitating comparisons, using the flag, `--draw-grid`. The following command would draw individual map files and grid files with a reverse colormap to emphasize unshared reactions.

{{ codestart }}
anvi-draw-kegg-pathways --external-genomes %(external-genomes)s -o output_dir --draw-grid --draw-individual-files --ko --colormap plasma 0.1 0.9 --reverse-overlay
{{ codestop }}

The following map grid reveals unique aspects of galactose metabolism among six related genomes.

![Map grid](../../images/anvi-draw-kegg-pathways/kos_database_grid.png)

### Pangenomic database

Pangenomes are treated similarly to multiple contigs databases. Rather than comparing the occurrence of KOs across contigs databases, consensus KO annotations of gene clusters are compared across genomes in a pangenomic database. Here is the basic structure of the command.

{{ codestart }}
anvi-draw-kegg-pathways -p %(pan-db)s -g %(genomes-storage-db)s -o output_dir --ko
{{ codestop }}

Pangenomes were constructed from 12 strains of two related species. The following maps were produced with a basic command.

![Three maps showing KOs from a pangenome](../../images/anvi-draw-kegg-pathways/kos_pan.png)

As with the comparison of contigs databases, it can be useful to reverse the colormap and create map grids to compare the KO content of genomes in the pangenome.

{{ codestart }}
anvi-draw-kegg-pathways -p %(pan-db)s -g %(genomes-storage-db)s -o output_dir --draw-grid --ko --colormap plasma 0.1 0.9 --reverse-overlay
{{ codestop }}

The following map grid reveals certain differences between the strains, and particularly the two species, in carbohydrate metabolism, with *faecalis* enriched in enzymes for xylose metabolism (towards the bottom of the map), and *faecium* enriched in enzymes for uronate metabolism (towards the top of the map).

![Pangenomic map grid](../../images/anvi-draw-kegg-pathways/kos_pan_grid.png)

