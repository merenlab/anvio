This artifact is a TAB-delimited text file written by %(anvi-draw-kegg-pathways)s. It lists the slices of the reaction elements drawn on %(kegg-pathway-map)s files, or the slices of the rings that color their borders.

With `--reaction-slice-by`, a box or line that stands for several KOs, KEGG reactions, or genes is drawn in slices, one for each accession or gene with a value on the map. Slices are too narrow to label, and a box label names only the first KO of the element. This table says what each slice stands for. The program writes it as `slice_key_reactions.txt` at the top of the output directory, beside the colorbars. One file covers every map of the run. Options that arrange the map files, such as `--categorize-files` and `--collate-files-by-map`, do not move it. See %(anvi-draw-kegg-pathways)s for how slices are drawn and how their shares are found.

With `--reaction-border-slice-by`, the ring that colors the border of a box is divided among the accessions or genes with border values. The program then writes a second table of the same form, `slice_key_reactions_border.txt`, beside the first. It lists the slices of the rings, as described in [Slices of rings](#slices-of-rings).

## File format

The file has a header line and these columns:

|Column|Description|
|:--|:--|
|`map`|The number of the KEGG map, such as `00010`.|
|`drawn_map`|Which map of that KEGG map was drawn. This is `unified` for the `unified` map, or the name of a sample or group for the map of that sample or group. The one map drawn from a %(kegg-reaction-txt)s without a `sample` column is `unified`.|
|`entry_id`|The ID of the element's entry in the map's KGML file.|
|`x`|The x coordinate of the element in the map's KGML file. For a box, this is its center. For a line, this is the point near the middle of its first piece where Biopython puts the label of a line.|
|`y`|The y coordinate of the same point.|
|`label`|The label of the element in the KGML file, such as `K01623...`.|
|`slice`|The position of the slice, counting from 1. The slices of a box count from the left. The slices of a line count from its first point in the KGML file, within each cycle of dashes. This is blank for an element drawn whole for a reason.|
|`accession`|The KO or KEGG reaction that the slice stands for. For a slice of a gene (`--reaction-slice-by gene`), this lists the gene's KOs or KEGG reactions in the element, in KGML order, joined by commas without spaces, such as `K01623,K01624`. These are the accessions of the element for which the gene has rows anywhere in the file, so a gene has the same label on every map. On a map where the gene has rows for only some of them, all are still listed.|
|`gene_id`|The gene that the slice stands for, with `--reaction-slice-by gene`. This is blank for a slice of an accession.|
|`value`|The value of the accession or gene on the drawn map. This is a number, or a categorical value with `--reaction-value-type categorical`. A gene with rows for two or more accessions of the element takes the mean of these rows in a sample by default. On the map of a group or the `unified` map, the value is the summary of the values in the samples or groups.|
|`share`|The share of the element's width or length that the slice takes. The shares of an element on a drawn map sum to 1.|
|`color`|The color of the slice, as a hex code. For an element drawn whole, it is the color of the whole element.|
|`drawn_whole`|Blank for a slice that is drawn. For an element drawn whole, it gives the reason, `over the cap` for a box or `line too short` for a line. A ring can also be drawn whole for the reason `fill not sliced`.|
|`bar`|Only in the table of a layer with a bar column (`--reaction-bar-column`). The bar value of the accession or gene on the drawn map, which its column of the box shows as a bar. A box drawn whole lists the bar values of the slices it would have had. The column is blank where the slice has no bar value, for the slices of lines, which take no bars, and on a map that shows no bars, such as the `unified` map without `--reaction-bar-sample-summary`. A column without a bar value is white in the default style, and keeps its color without a strip in the `strip` style.|

## Rows

Each row is one slice on one drawn map. An element on the `unified` map and on the maps of three samples has rows for each of the four maps. Its slices can differ from map to map, since an accession or gene gets a slice only on the maps where it has a value. The same position can then hold different accessions or genes on different maps.

An element with one slice is drawn whole in the color of that slice. It has one row, with `slice` 1 and `share` 1.0.

An element with more slices than fit is drawn whole in the color of a stand-in value. A box is drawn whole where a slice of equal width would be narrower than 8 points. In a layer with a border column, the slices of a box fit inside its ring, which leaves 85 of its 92 points. A line is drawn whole where any of its pieces would give a slice of equal length shorter than the line is wide. Only the part of the piece that compound circles leave visible counts. Such an element has a row for each slice it would have had. These rows have a blank `slice`, the share and value that each slice would have had, and the color of the whole element. `drawn_whole` gives the reason. The stand-in value itself is not listed.

Elements left uncolored have no rows. These include an element whose slices with values all have a weight of 0, and an element drawn whole whose stand-in value is undefined.

The rows are sorted by `map`, then by `drawn_map`, then by entry ID, and then by `slice`. The `unified` map comes first, and the maps of the samples or groups follow in the order in which they are drawn. Rows name maps, not files. A map drawn only as a panel of a grid is listed under its sample or group, like any other. The file is written whenever slicing is on. It holds only the header if no slice was drawn.

## Examples

Here are rows for one box of map `00010`, drawn in equal slices on the maps of two samples. The box has two slices on the map of `SAMPLE_1`. On the map of `SAMPLE_2`, only K16306 has a value, so it fills the box.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|SAMPLE_1|18|966|814|K01623...|1|K01622||36.5|0.5|#a01a9c||
|00010|SAMPLE_1|18|966|814|K01623...|2|K16306||15.5|0.5|#f79143||
|00010|SAMPLE_2|18|966|814|K01623...|1|K16306||30.5|1.0|#c13b82||

Here are rows for the same box from a run sliced by gene, on the maps of two other samples, `control` and `treated`. The rows come from the example file for slices of genes in %(anvi-draw-kegg-pathways)s. `gene_1` has rows for K01623 and K01624, so it is one slice. On the map of `treated`, its rows are 1.2 for K01623 and 0.8 for K01624, so its value is their mean, 1.0. `gene_2` has a row for K01623 only.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|control|18|966|814|K01623...|1|K01623,K01624|gene_1|0.2|0.5|#db5c68||
|00010|control|18|966|814|K01623...|2|K01623|gene_2|0.5|0.5|#a31e9a||
|00010|treated|18|966|814|K01623...|1|K01623,K01624|gene_1|1.0|0.5|#0d0887||
|00010|treated|18|966|814|K01623...|2|K01623|gene_2|-0.4|0.5|#f0f921||

Here are the first three of the 14 rows of a box drawn whole, from a run without samples. The box would have had 14 slices, which is more than the 11 that fit. The whole box takes one color.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00190|unified|255|960|772|K00234...||K00234||16.5|0.09090909090909091|#c33d80|over the cap|
|00190|unified|255|960|772|K00234...||K00235||7.5|0.12121212121212122|#c33d80|over the cap|
|00190|unified|255|960|772|K00234...||K00236||21.5|0.030303030303030304|#c33d80|over the cap|

## Columns of bars

With `--reaction-bar-column`, the table of the layer has the column `bar` at the end. Each column of a box shows the bar of its own accession or gene, so the bar value belongs to the slice of its row. A box drawn whole shows one bar, the mean of the bar values of its slices, weighted by their shares, and that mean is not listed. Here are rows for the box of the first example from a file with `expression` as its value column and `coverage` as its bar column, on the `unified` map, where `--reaction-sample-summary mean` and `--reaction-bar-sample-summary mean` pool the samples.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|bar|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|unified|18|966|814|K01623...|1|K01623||6.0|0.25|#f79143||13.333333333333334|
|00010|unified|18|966|814|K01623...|2|K01624||7.666666666666667|0.25|#e56a5d||21.0|
|00010|unified|18|966|814|K01623...|3|K01622||16.0|0.25|#0d0887||17.25|
|00010|unified|18|966|814|K01623...|4|K16306||2.5|0.25|#f0f921||15.0|

A file whose only column of values is its bar column colors its fills by presence. Its bars are sliced on their own, among the accessions or genes with bar values, and its table lists these columns. Their `value` is blank, since the fill has none, and their `color` is the one color of the fill. Here are rows of the same box from such a file, with `coverage` as its bar column and no samples.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|bar|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|unified|18|966|814|K01623...|1|K01624|||0.3333333333333333|#2ca02c||33.5|
|00010|unified|18|966|814|K01623...|2|K01622|||0.3333333333333333|#2ca02c||36.5|
|00010|unified|18|966|814|K01623...|3|K16306|||0.3333333333333333|#2ca02c||23.0|

The tables of layers without a bar column, and the tables of rings, have no `bar` column.

## Slices of rings

The table of the slices of rings, `slice_key_reactions_border.txt`, has the same columns. `value` and `color` are those of the border column, and `drawn_whole` gives why a ring is drawn whole. The rings of a run are either aligned with the slices of the fill or divided on their own, and the positions and shares depend on which.

- **Aligned with the fill.** Where `--reaction-slice-by` and `--reaction-border-slice-by` name the same unit, and the file has a value column for the fill, each column of the fill takes the border value of its own accession or gene, as bands above and below it. A row lists such a column. Its `slice` and `share` are those of the column in `slice_key_reactions.txt`. A column without a border value has no row.
- **Divided on its own.** Otherwise, each ring is divided among its accessions or genes with border values. The slices count clockwise from the top left corner of the box, and `share` is the share of the ring that a slice takes.

A ring of one slice is drawn whole in the color of that slice, with `slice` 1 and `share` 1.0. Rings can also be drawn whole for these reasons, in the color of a stand-in value:

- `over the cap`. An aligned ring is drawn whole where its fill has more slices than fit. Its rows list the units with both a slice of the fill and a border value, with the shares of the fill. A ring divided on its own is drawn whole where it would have more than 11 slices.
- `fill not sliced`. An aligned ring is drawn whole where the fill of its box is not sliced, since the map shows presence or the box has no fill value there. Its rows list every unit with a border value, with equal shares.

Rings left out have no rows. These include a ring whose slices with values all have a weight of 0, and a ring drawn whole whose stand-in value is undefined.

Here are rows of the aligned rings of the box of the first example, on the map of `SAMPLE_1`. The fill was sliced by KO, with shares from a column of weights, and the border column holds the categorical values `up` and `down`.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|SAMPLE_1|18|966|814|K01623...|1|K01624||up|0.47619047619047616|#d62728||
|00010|SAMPLE_1|18|966|814|K01623...|2|K01622||down|0.047619047619047616|#1f77b4||
|00010|SAMPLE_1|18|966|814|K01623...|3|K16306||up|0.47619047619047616|#d62728||

Here are rows of the same box on the `unified` map of that run, where the fills show presence. The fill is not sliced, so the ring is whole. Its color is that of `up`, which `plurality_order` chooses among the three KOs.

|map|drawn_map|entry_id|x|y|label|slice|accession|gene_id|value|share|color|drawn_whole|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|00010|unified|18|966|814|K01623...||K01624||up|0.3333333333333333|#d62728|fill not sliced|
|00010|unified|18|966|814|K01623...||K01622||down|0.3333333333333333|#d62728|fill not sliced|
|00010|unified|18|966|814|K01623...||K16306||up|0.3333333333333333|#d62728|fill not sliced|
