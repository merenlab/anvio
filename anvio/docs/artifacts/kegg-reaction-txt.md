This artifact is a TAB-delimited text file that can be used to color the **reaction layer** of %(kegg-pathway-map)s files drawn by %(anvi-draw-kegg-pathways)s.

The "reaction layer" contains reaction lines in global and overview maps and boxes in standard maps. This file, passed to %(anvi-draw-kegg-pathways)s with the `--reaction-txt` option, can be paired with a %(kegg-compound-txt)s file (`--compound-txt`) to color the reaction and compound layers of a map at once. One or more %(contigs-db)ss, a %(pan-db)s, or a metabolic model file (see %(reaction-network-json)s) can alternatively be used by the program to color the reaction layer.

## File format

The file must have a header line with a required `accession` column; the other columns are optional.

Each column name may appear only once, and names that start with `__` are reserved for anvi'o's own columns. Double quotes around a whole field or name, as R writes them, are dropped. Any other quote is kept as text, so a number that holds one is refused.

|Column|Required|Description|
|:--|:--|:--|
|`accession`|YES|A KEGG ID. Accessions in the file must **all** be **KEGG Ortholog (KO) IDs** (`K#####`) or **KEGG reaction IDs** (`R#####`) — a single file cannot mix the two. Both KOs and Reactions address the reaction layer, so combining them can clash with two scales on the same map elements.|
|`gene_id`|NO|An identifier for the gene a value came from. The "gene" can be anything that carries the annotation, such as a gene, a transcript, a protein, or a gene cluster, and any text works as its ID. Several genes may carry the same accession, and their values are combined into that accession's value (by default their `sum` — see `--reaction-gene-aggregation`). With `--reaction-slice-by gene`, each gene is drawn as a slice of its own instead, and `--reaction-border-slice-by gene` does the same for the rings of a border column. Every row then needs a gene ID, and one gene ID must name the same gene in every sample (see [Drawing elements in slices](#drawing-elements-in-slices)).|
|`sample`|NO|A sample of origin, enabling comparison across samples (and, with a %(groups-txt)s, across groups of samples). If this column is present, every row must have a sample.|
|*a value column*|NO|A SINGLE column with any name of your choosing, for coloring by a value instead of by presence/absence. Its values are numbers by default, for coloring by a continuous quantitative value. Its meaning is up to you (coverage, expression, flux, ...). With `--reaction-value-type categorical`, its values are categorical values instead, such as `up` and `down` (see below). If this column is present, every row must have a value. This column is **auto-detected** as not being named `accession`, `gene_id`, or `sample`, and not being named by `--reaction-slice-weights` or `--reaction-border-column`.|
|*a weight column*|NO|A column of numbers that sizes the slices of reaction elements drawn with `--reaction-slice-by` (see below), and the slices of rings that `--reaction-border-slice-by` divides on their own. It can have any name. It is not auto-detected, so name it with `--reaction-slice-weights`. Without that option, the file would have two candidate value columns and be refused. Every row needs a weight, which is a finite number of 0 or more. The file still needs a separate value column, which colors the slices, unless the weights size only the slices of rings.|
|*a border column*|NO|A column whose values color the borders of reaction boxes on standard maps. Each border is a ring just inside the black outline of a box, and the value column still colors the fill. The column can have any name. It is not auto-detected, so name it with `--reaction-border-column`. Without that option, a file with a value column would have two candidate value columns and be refused. It cannot be `accession`, `gene_id`, `sample`, or the weight column. Its values are numbers by default, and each must be a finite number. With `--reaction-border-value-type categorical`, its values are categorical values instead. A row can leave its border cell empty, and then it has no border value. A border column with no values at all is refused. A file with a border column is read as text, so `NA` is a name in every column, including `sample` and `gene_id`. Reaction lines show no borders. See the example below.|

To color by **presence** of accessions, do not include a value column. A border column is not a value column, so a file with a border column and no value column colors its fills by presence.

With a value column or a border column, no two rows may describe the same thing — the same accession, gene, and sample. Repeated rows are ambiguous, as they can be replicate measurements to average or separate contributions to add, so combine them yourself in whichever way suits your data. Note that several *different* genes carrying one accession are not repeats, and are the ordinary case of the same type of protein being encoded by different genes.

A reaction element on a map is defined by one or more KEGG reactions (~23%% of elements have multiple), which may be carried out by one or more KOs. A file of **KO** accessions addresses the reaction layer by KO ID, while a file of **KEGG reaction** accessions addresses it by reaction ID — useful when your values are per-reaction, as with metabolic model fluxes. Since a single line or box can stand for several KOs or reactions, its color is the aggregate of whichever of its accessions are in the file (`--reaction-accession-aggregation`; see below). With `--reaction-slice-by`, the element is drawn in slices instead, one for each of its accessions or genes with a value (see [Drawing elements in slices](#drawing-elements-in-slices)).

## Coloring by presence or quantitative values

The inclusion of different columns results in different coloring modes.

**No sample column, no value column**

Reactions are colored a single presence color, green by default — set it with `--reaction-color`, or use `--original-color` for the reference map's colors.

**Sample column, no value column**

Reactions are colored by sample (or, with a %(groups-txt)s, group) count or membership, exactly as when comparing multiple %(contigs-db)ss or genomes of a %(pan-db)s. Counts are drawn either in discrete bands, one per count (`--reaction-sample-summary count`), or as a gradient from the lowest count to the highest (`--reaction-sample-summary count_continuous`), the latter being the only one that works when there are more samples than the colormap has distinguishable colors — without `--reaction-sample-summary`, discrete automatically switches to continuous given enough samples. By default, the scale stops at the highest count actually observed rather than at the number of samples, so that the colors spread over the counts that occur; use `--count-scale-max` to change the maximum value setting.

With `--reaction-color` or `--original-color`, the single color or original colors of the reference map override dynamic coloring, simply showing whether a reaction is present in any of the samples.

**Value column**

Reactions are colored by the continuous value through a sequential colormap (`--reaction-colormap`, default `plasma_r`). The colorbar is labeled by the value column's header. A cyclic colormap called `clocktime` was defined in anvi'o to show a value column of clock times in hours. Give clock times `--reaction-value-period 24`: each color scale of values then runs from 0 to 24 h and is colored by `clocktime`. Each reduction of the values takes their circular mean, `circular_mean`, so that 23.5 h and 0.5 h average to 0 h, not 12 h. See %(anvi-draw-kegg-pathways)s.

Per-gene values are aggregated to a per-accession value, and a map element's constituent KO and reaction accessions to a per-element value, by `--reaction-gene-aggregation` and `--reaction-accession-aggregation` (both `sum` by default) — these reductions happen within each sample.

With a `sample` column, `--draw-individual-files` and/or `--draw-grid` color maps of value column data that share a single `colorbar_reactions_samples.pdf` so that samples are comparable on the same scale. How the `unified` map, and with %(groups-txt)s the per-group maps, summarize samples is a separate choice (see below).

Each of these scales spans exactly the values it is given, which a few extreme elements can stretch until the rest cannot be told apart. `--reaction-value-limits` bounds the `unified` map's scale and `--reaction-individual-value-limits` the scale its per-sample or per-group maps share. These parameters take an argument of two values, a minimum and a maximum, either of which can be `none` to leave that end open while setting the other. A limit sets its end of the scale whether or not any value reaches it. The colorbar labels each limit. A limit that values lie past is marked `≤` or `≥`, since its color then stands for that value or anything past it.

A scale also sits wherever its own values leave it, which for a signed quantity means the neutral middle color of a diverging colormap lands somewhere other than zero. `--reaction-value-center` and `--reaction-individual-value-center` put a value at the middle of those two scales: used as a bare flag each centers its scale on 0, and each takes a number as an argument to center it on that instead. The scale is only ever widened to make room, running the same distance either side of the center and staying linear, so nothing is clipped that was not clipped already.

The two scales are colored by one colormap, `--reaction-colormap`, unless `--reaction-individual-colormap` gives the per-sample or per-group scale one of its own. That is worth doing where the two show quantities of different kinds: with `--reaction-sample-summary std` the `unified` map shows how much samples disagree rather than how large their values are, so drawing it in the colors that show magnitude on the per-sample maps invites mistaking one scale for the other.

## Coloring by categorical values

A value column can hold **categorical values** instead of numbers. These are names, such as `up`, `down` and `same`, or the names of clusters. Give `--reaction-value-type categorical` to read the column this way. Values written as numbers, such as `1` and `2`, are names too. A value is read as text, so `NA` is a value like any other. Every row needs a value. Each accession can have one value per gene and sample, so repeated rows are refused, as they are for numbers.

Each value takes a color of its own. By default, the values take colors in sorted order. Where every value is a number, such as `-1`, `0` and `1`, the values are sorted as numbers. Otherwise the digits in the values are compared as numbers, so `cluster 2` comes before `cluster 10`. Up to 9 values take the colors of Matplotlib's `tab10` without its gray. The 10th to 18th values take the light colors of `tab20` without its light gray. Gray is left out because global maps draw the reactions they do not highlight in gray. To choose the colors yourself, give a %(kegg-value-colors-txt)s file with `--reaction-value-colors`. The order of its rows is then the order of the values. More than 18 values need such a file, except in a grouped run without a sample rule, where no map shows the values. The legend lists the values in order, with the first value at the top.

**Choosing one value where several meet.** Each row has one value. Several values can still meet in two places within a sample. The genes of one accession can have different values. The KOs or reactions of one map element can have different values. A rule chooses one value in each place. `--reaction-gene-aggregation` sets the rule for the genes of an accession, and `--reaction-accession-aggregation` the rule for the accessions of an element. Each gene or accession gives one vote.

|Rule|What it does|
|:--|:--|
|`plurality_order`|The default. Takes the most common value, even with less than half of the votes. A tie goes to the value that comes first in the order of the values.|
|`plurality_random`|Takes the most common value. A tie is broken at random.|
|`plurality`|Takes the most common value. A tie gives no value.|
|`majority`|Takes the value with more than half of the votes. Without one, there is no value.|
|`order`|Takes the value that comes first in the order of the values, whatever the votes.|
|`random`|Takes one of the values at random, whatever the votes.|

The order of the values is the row order of the colors file, or else the sorted order described above. If a colors file lists `up` before `same`, `order` chooses `up` for a map element whose KOs have both. A random choice comes from `--seed` (42 by default) and from what the value is chosen for, such as the KOs of one element in one sample. The same element then gets the same value on every map, and in every run with the same seed. The run reports how many choices each rule made, and how: by plurality, by majority, by value order, at random, or without a value.

`plurality` and `majority` can give no value. An accession given no value has no value in that sample, so it does not vote in its map elements. A map element given no value is left uncolored. That looks like an element the file lacks, so a warning says how many there are. These two rules suit a map that should show only clear choices. With two voters, `majority` means that both agree.

On global and overview maps drawn from reactions alone, each compound takes a color from the reactions that touch it. The accession rule chooses among their values. `plurality` and `majority` can leave a compound without a value. It then keeps the color of a compound that no colored reaction touches. A warning says how many such compounds there are. `--reaction-reverse-overlay` puts the last value on top where reactions overlap. It does not change which value a rule chooses.

**Samples.** With a `sample` column, the map of each sample shows that sample's values. Draw these maps with `--draw-individual-files` and/or `--draw-grid`. They share one legend, `colorbar_reactions_samples.pdf`. A value in one sample need not mean the same thing in another sample. For example, cluster 1 of one sample need not be cluster 1 of another. The `unified` map therefore shows presence across the samples by default. It counts the samples that contain each element, however few samples there are. A run that draws only the `unified` map then colors nothing by the values. anvi'o warns about it. Coloring by which samples contain an element would use the `tab10` colors of the values. A color would then mean a sample on one map and a value on another. `--reaction-sample-summary membership` still asks for it. Where a value does mean the same thing in every sample, give a rule to `--reaction-sample-summary`. The rule then chooses one value per element across the samples, and the `unified` map shows it. Each sample with a value for the element gives one vote. A sample without one does not vote.

**Groups of samples.** With a %(groups-txt)s, maps are drawn for groups rather than samples. By default, each group's map counts the group's samples, and the `unified` map counts the groups. Such a run colors nothing by the values, and anvi'o warns about it. A rule given to `--reaction-sample-summary` chooses one value per element among the samples of each group. Each group's map then shows that value. A rule given to `--reaction-group-summary` then chooses one value across the groups for the `unified` map. It needs a rule for the samples, since the groups otherwise have no values. Each sample with a value gives one vote in its group. Each group with a value gives one vote across the groups, whatever its size. With few voters, ties are common. With two groups, every disagreement between them is a tie, so `plurality_order` settles it by the order of the values. `majority` among the samples of a group shows a value only where more than half of the samples with a value agree. A sample without a value does not count, so one sample with a value is enough. Where the `unified` map counts the groups, `--group-threshold` decides when a group contains an element, as for any layer.

The options that shape a scale of numbers are refused for categorical values. These are the value limits and centers, `--reaction-value-period`, `--reaction-element-normalization` and `--reaction-individual-colormap`. `--reaction-colormap` colors only the presence on the `unified` map.

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --reaction-value-type categorical \
                        --draw-grid \
                        -o output_dir
{{ codestop }}

## Summarizing across samples and groups

There are four distinct reductions of values, each with its own option, plus one rescaling that is not a reduction. The first two reductions apply whenever the file has a value column. The next two, and the rescaling, need a `sample` column. Slices drawn with `--reaction-slice-by` or `--reaction-border-slice-by` add two reductions of weights, in the last two rows of the table.

|Level|Option|What it does|
|:--|:--|:--|
|genes of an accession|`--reaction-gene-aggregation` (an aggregation; `sum` by default; or a rule for categorical values)|the genes carrying an accession → that accession's value. With `--reaction-slice-by gene`, each gene keeps its own value, and this option is refused.|
|accessions of a map element|`--reaction-accession-aggregation` (an aggregation; `sum` by default; or a rule for categorical values)|the constituent accessions of a map element → that element's value. With `--reaction-slice-by accession`, each slice keeps the value of its own accession instead. With `--reaction-slice-by gene`, the rows of one gene for the accessions of an element → that gene's value, with `mean` as the default. An aggregation that changes a single value, such as `std`, is then refused.|
|across samples|`--reaction-sample-summary` (`count`, `count_continuous`, `membership`, an aggregation, or a rule for categorical values)|a set of samples → one value or one presence value per map element, or per slice with `--reaction-slice-by`|
|across groups|`--reaction-group-summary` (`count`, `count_continuous`, `membership`, an aggregation, or a rule for categorical values)|the groups of a %(groups-txt)s → one value or one presence value per map element, or per slice with `--reaction-slice-by`|
|each sample's or group's value for an element|`--reaction-element-normalization` (a normalization; none by default)|an element's value in one sample or group → that value rescaled against the element's values across all samples or groups|
|weights of the genes of an accession|`--reaction-slice-weight-aggregation` (`sum` by default, or `mean`)|the weights of the genes carrying an accession → that accession's weight, which sizes its slice. This applies only where `--reaction-slice-by accession` slices the elements or `--reaction-border-slice-by accession` slices the rings.|
|shares across samples and groups|`--reaction-slice-share-summary` (`share_mean` by default, or `pooled_sum`)|the shares of a slice in a set of samples, or in the groups → its share on the map of a group or the `unified` map|

An **aggregation** is `sum` (default), `mean`, `max`, `min`, `median`, `std` — or any other unsuggested pandas aggregation that reduces a series of numbers to one number, such as `var` or `sem`. Names that transform rather than reduce (`cumsum`) or that only a grouping offers (`first`) are rejected. Where an aggregation is undefined for the values available, as `std` is for a single value, the affected elements are left uncolored. `circular_mean` averages values that repeat after a period, such as clock times, and requires `--reaction-value-period`. `angular_deviation` measures how far such values are spread apart, on the `unified` map only.

A summary pools the values of a **map element**, not of an accession. An element's value in a sample is the aggregate of its accessions in that sample. Without `--reaction-element-normalization`, this is the value that the sample's own map draws, so the sample summary pools what the maps of the individual samples show. With groups, the group summary then pools the values of the groups. With `--reaction-slice-by`, a summary pools the values of each **slice** instead, which are the values of its one accession or gene.

The sample summary drives the `unified` map when there are no groups, and each per-group map (produced when using `--draw-individual-files` or `--draw-grid`) when there are. The group summary drives the `unified` map when there are groups.

Note that with groups, the per-group map either shows the count of that group's sample — the default case — or, with a `--reaction-sample-summary` aggregation argument, the samples' pooled value; `--reaction-sample-summary` is therefore not permitted with a presence name (`count`, `count_continuous` or `membership`) when using groups. For categorical values, a rule given to `--reaction-sample-summary` colors each group's map by the value that it chooses among the group's samples.

The sample and group summaries in the `unified` maps both default to **presence**, even with a quantitative value column. The default is `membership` for 3 or fewer samples or groups. Above that, it is `count`. Where there are more samples or groups than the colormap has distinguishable colors, it is `count_continuous`. Presence is the default because it is meaningful across any set of samples. Pooling values is only meaningful across commensurable samples. For example, it makes sense to average transcript abundances across samples from replicate conditions, but not from unrelated conditions. Pooling values across samples is triggered by providing an aggregation argument to `--reaction-sample-summary` or `--reaction-group-summary`. `std` at this level maps how much samples disagree. A layer of categorical values counts the samples by default, however few there are.

A layer can show **presence in the overview and values per sample or group**. With the default presence summary, the `unified` map gets a colorbar of sample or group counts or memberships, while per-sample maps keep a continuous colorbar. The sample and group levels can be treated independently — replicate samples can be averaged within each condition (group) while the `unified` map summarizing all conditions shows in how many conditions each reaction occurs, with a group containing a reaction when at least half of its samples do:

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --groups-txt %(groups-txt)s \
                        --group-threshold 0.5 \
                        --reaction-sample-summary mean \
                        --reaction-group-summary count \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

## Normalizing an element across samples

With a `sample` column and a value column, each reaction element in a sample map is colored by value aggregated across accessions in the element. `--reaction-element-normalization` instead colors each element relative to all of the samples. Each element is rescaled against its own values across the samples, so under `relative_to_mean` a value of +0.25 means the sample has 25%% more of that element than it averages across samples, and -0.31 means 31%% less.

This cannot be worked out in the file itself before analyzing the maps, since each element in each map is related to one or more accessions — the same accession factors into different rescaled element values on each map.

anvi'o encodes a number of normalizations with names recognized as an argument: `relative_to_mean`, `relative_to_median`, `difference_from_mean`, `difference_from_median`, `log2_ratio_to_mean`, `log2_ratio_to_median`, `z_score`, `rank`, `fraction_of_max`, and `fraction_of_total`. With `--reaction-value-period`, the only normalization is `difference_from_circular_mean`, the signed offset from the element's circular mean. A different name is taken to be a pandas Series method that transforms each value into a new value, such as `abs`.

Only the maps of individual samples or groups are rescaled, so element normalization needs `--draw-individual-files` and/or `--draw-grid`. The `unified` map summarizes the unnormalized values, so is unaffected by normalization.

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --reaction-element-normalization relative_to_mean \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

## Drawing elements in slices

With `--reaction-slice-by accession`, %(anvi-draw-kegg-pathways)s draws a reaction element in slices, one for each of its accessions with a value on the map. A file of KOs is sliced by KO, and a file of KEGG reactions by KEGG reaction. Each slice takes the color of its accession's value. Slicing needs a value column, of numbers or of categorical values. Maps that show presence are drawn whole.

Without weights, the slices of an element are equal. A weight column, named with `--reaction-slice-weights`, sizes each slice in proportion to its weight instead. When slicing by accession with a `gene_id` column, the weights of the genes of an accession are added up, or averaged with `--reaction-slice-weight-aggregation mean`. The maps of groups and the `unified` map pool the shares of the samples, as `--reaction-slice-share-summary` sets.

|accession|sample|log2_ratio|coverage|
|:--|:--|:--|:--|
|K01623|treated|1.2|30.5|
|K01624|treated|-0.4|10.2|
|K01623|control|0.3|28.0|
|K01624|control|0.1|12.9|

The `log2_ratio` column colors the slices, and `coverage` sizes them. K01623 and K01624 share a box on map `00010`, so that box has two slices on the map of each sample. On the map of `treated`, K01623 takes 0.75 of the box.

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

With `--reaction-slice-by gene`, each gene of the `gene_id` column gets a slice instead, and every row needs a gene ID. A gene with rows for two or more accessions of an element is one slice. In each sample, its value is the mean of these rows, or their circular mean with `--reaction-value-period`. A gene with one row takes the value of that row. Without a period, `--reaction-accession-aggregation` can name another aggregation that returns a single value unchanged, such as `sum` or `max`. An aggregation that changes a single value, such as `std` or `count`, is refused. For categorical values, the accession rule chooses among the rows. The weight of the gene is the mean of the weights of these rows. Each gene is a slice of its own, so `--reaction-gene-aggregation` is refused. `--reaction-slice-weight-aggregation` is refused too, unless the rings of a border column are sliced by accession.

The maps of groups and the `unified` map pool the values of each gene across samples, so one gene ID must name the same gene in every sample. Gene caller IDs start at 0 in every %(contigs-db)s, so IDs from separate assemblies must be made unique, for example with a sample prefix. Otherwise, unrelated genes are pooled as one. On the drawn maps that pool samples, a warning names the gene IDs that have no accession in common between their rows in two samples. A gene of one contigs database has the same accessions in every sample, so such an ID likely names two genes. A gene that only lacks rows for some of its accessions in a sample is not named. Another warning names the maps where no gene ID occurs in two or more samples. Nothing is pooled there.

With a border column, the columns of a sliced box fit inside its ring. Inside the ring, a box is 85 points wide rather than 92. A box in a layer with a border column therefore holds up to 10 slices rather than 11. The ring stays whole unless `--reaction-border-slice-by` divides it among the accessions or genes with border values. Where `--reaction-slice-by` names the same unit, each column of the fill takes the border value of its own accession or gene, as bands above and below it. Otherwise, each ring is divided on its own, clockwise from the top left corner, and the weight column sizes its slices. The slices of rings are listed in `slice_key_reactions_border.txt`.

The program writes a %(kegg-slice-key-txt)s that lists each slice drawn on each map. See %(anvi-draw-kegg-pathways)s for which maps are sliced, how the shares of samples are pooled, how slices are drawn and ordered, and how to read sliced maps.

## Examples

### KO coverage (quantitative)

Color reactions by the coverage of their KOs, aggregated from genes. The `coverage` column is auto-detected as the value column.

|accession|gene_id|coverage|
|:--|:--|:--|
|K00844|gene_1|12.4|
|K12407|gene_2|3.1|
|K00845|gene_3|20.0|

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        -o output_dir
{{ codestop }}

### Reaction flux (quantitative, reaction IDs)

Color reactions directly by a per-reaction value, such as a flux from a metabolic model. Use reaction (`R#####`) accessions rather than KO accessions.

|accession|flux|
|:--|:--|
|R00200|8.7|
|R00299|1.2|
|R01068|4.5|

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        -o output_dir
{{ codestop }}

### Expression status per sample (categorical)

Color each sample's reactions by a status of its KOs. The `status` column holds names, which `--reaction-value-type categorical` reads as categorical values. A colors file gives each status a color of its own.

|accession|sample|status|
|:--|:--|:--|
|K00844|treated|up|
|K00845|treated|same|
|K00844|control|down|

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --reaction-value-type categorical \
                        --reaction-value-colors %(kegg-value-colors-txt)s \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

### KO coverage with a status border (quantitative fill, categorical border)

Color the fill of each box by the coverage of its KOs, and its border by a status. `coverage` is auto-detected as the value column. `status` is named as the border column, and `--reaction-border-value-type categorical` reads its values as categorical values. The last row leaves its status empty, so K01810 has no border value.

|accession|gene_id|sample|coverage|status|
|:--|:--|:--|:--|:--|
|K00844|gene_1|treated|12.5|up|
|K00844|gene_1|control|3.1|down|
|K00845|gene_2|treated|8.0|same|
|K00845|gene_2|control|7.7|same|
|K01810|gene_3|treated|4.2||

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --reaction-border-sample-summary plurality_order \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

The map of each sample shows that sample's coverage in the fills and its status in the borders. The border summaries have no default, so the `unified` map shows borders here only because `--reaction-border-sample-summary` gives a rule. Its fills show presence across the samples, as they do by default. Each border option, such as `--reaction-border-colormap`, works like the option of the fill whose name lacks `border-`, but some defaults differ. The border aggregations default to `mean` for numbers, rather than `sum`. Numeric border values take `viridis` from 0.2 to 1, and categorical border values take the 7 colors of Matplotlib's `Dark2` without its gray. The borders get keys of their own, such as `colorbar_reactions_border.pdf` and `colorbar_reactions_border_samples.pdf`. `--reaction-border-slice-by` divides each ring among the accessions or genes with border values (see [Drawing elements in slices](#drawing-elements-in-slices)). See %(anvi-draw-kegg-pathways)s for the border options, their defaults, and the maps that show borders.

### KO presence (single color)

With no value column and no `sample` column, present reactions are drawn one color.

|accession|
|:--|
|K00844|
|K00845|
|K12407|

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        -o output_dir
{{ codestop }}

### Comparing samples

Add a `sample` column to compare origins (genomes, metagenomes, transcriptomes, ...) on the same maps, as when comparing multiple contigs databases. Without a value column, an accession is present in a sample only when a row places it there (below, K00844 is in genome_A only, K00845 in both); add a value column to color each sample's map by its own magnitude. Combine with a %(groups-txt)s to compare groups of samples, and set `--reaction-sample-summary`/`--reaction-group-summary` to choose whether the summary views show presence or pooled values. See %(anvi-draw-kegg-pathways)s for the full set of sample, group, and colormap options.

|accession|sample|
|:--|:--|
|K00844|genome_A|
|K00845|genome_A|
|K00845|genome_B|

### With a compound layer

Pair this file with a %(kegg-compound-txt)s file to color both layers on the same maps, each by its own mode. For example, the presence of KOs in a genome can be shown alongside quantitative metabolomics of compounds.

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --compound-txt %(kegg-compound-txt)s \
                        -o output_dir
{{ codestop }}
