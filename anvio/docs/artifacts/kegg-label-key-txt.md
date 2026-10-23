This artifact is a TAB-delimited text file written by %(anvi-draw-kegg-pathways)s. It lists the reaction boxes drawn without their labels on %(kegg-pathway-map)s files, because of the bars inside them.

With `--reaction-bar-column`, a map whose bar scale has values on both sides of 0 draws bars that grow both ways from a gray zero line. The label of a box fills its upper half, where it would hide the zero line and the ends of short bars, so the boxes of the default `fill` style are drawn without their labels. A box whose bar runs up its height, with `--reaction-bar-orientation vertical`, is drawn without its label too, since the top of its bar would run behind the label. This table says what each of these boxes stands for. The program writes it as `label_key_reactions.txt` at the top of the output directory, beside the bar keys. One file covers every map of the run. It is written where the bars are vertical, or where the bar scale of the `unified` map, or the scale shared by the maps of the samples or groups, has values on both sides of 0. Boxes of the `strip` style keep their labels, so a run of strips writes no table. See %(anvi-draw-kegg-pathways)s for how bars are drawn.

## File format

The file has a header line and these columns:

|Column|Description|
|:--|:--|
|`map`|The number of the KEGG map, such as `00010`.|
|`drawn_map`|Which map of that KEGG map was drawn. This is `unified` for the `unified` map, or the name of a sample or group for the map of that sample or group. A %(kegg-reaction-txt)s without a `sample` column draws only the `unified` map, unless a %(kegg-compound-txt)s with a `sample` column draws maps of samples too. Its boxes are then the same on every map, and each map lists them.|
|`entry_id`|The ID of the box's entry in the map's KGML file.|
|`x`|The x coordinate of the center of the box in the map's KGML file.|
|`y`|The y coordinate of the center of the box.|
|`label`|The label of the box in the KGML file, such as `K01623...`. This is the label that was not drawn.|
|`accession`|The KOs or KEGG reactions of the box that the file has on the drawn map, in KGML order, joined by commas without spaces, such as `K01622,K16306`. On the map of a sample, these are the accessions with rows in that sample. On the map of a group or the `unified` map, they are those of any of its samples. A file without a `sample` column gives every map all of its accessions.|
|`bar`|The bar value of the box on the drawn map. On the map of a group or the `unified` map of a file with a `sample` column, it is the summary of the bar values in the samples or groups.|

## Rows

Each row is one box on one drawn map. The rows are sorted by `map`, then by `drawn_map`, then by entry ID. The `unified` map comes first, and the maps of the samples or groups follow in the order in which they are drawn. Rows name maps, not files. A map drawn only as a panel of a grid is listed under its sample or group, like any other.

A box is listed on a drawn map only where its bar is vertical or its bar scale has values on both sides of 0. The `unified` map and the maps of the samples or groups have bar scales of their own, so a box can be listed on one and labeled on the other. A box left uncolored for want of a bar value keeps its label and is not listed. Boxes are listed whether or not `--reaction-labels` labels standard maps.

## Examples

Here are rows for three boxes of map `00010`, from a file of log ratios in three samples, with `--reaction-bar-sample-summary mean`. The bar scale of the `unified` map runs from -1 to 3.5, and the scale shared by the maps of the samples from -1 to 5. On the `unified` map, the box of entry 18 shows the mean of its bar values across the samples. On the map of `SAMPLE_1`, only two of its KOs have rows.

|map|drawn_map|entry_id|x|y|label|accession|bar|
|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|unified|18|966|814|K01623...|K01624,K01622,K16306|3.4166666666666665|
|00010|unified|61|916|968|K00134...|K00134|-0.5|
|00010|unified|159|1148|1592|K01007|K01007|-1.0|
|00010|SAMPLE_1|18|966|814|K01623...|K01622,K16306|2.75|
