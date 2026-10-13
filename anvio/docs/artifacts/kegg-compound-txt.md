This artifact is a TAB-delimited text file that can be used to color the **compound layer** of %(kegg-pathway-map)s files drawn by %(anvi-draw-kegg-pathways)s.

The "compound layer" contains compound circles in maps. This file, passed to %(anvi-draw-kegg-pathways)s with the `--compound-txt` option, can be paired with a %(kegg-reaction-txt)s file (`--reaction-txt`) to color the compound and reaction layers of a map at once.

## File format

The file must have a header line with a required `accession` column; the other columns are optional.

|Column|Required|Description|
|:--|:--|:--|
|`accession`|YES|A KEGG compound ID (`C#####`).|
|`sample`|NO|A sample of origin, enabling comparison across samples (and, with a %(groups-txt)s, across groups of samples). If this column is present, every row must have a sample.|
|*a value column*|NO|A SINGLE column with any name of your choosing, for coloring by a value instead of by presence/absence. Its values are numbers by default, for coloring by a continuous quantitative value. Its meaning is up to you (concentration, intensity, ...). With `--compound-value-type categorical`, its values are categorical values instead, such as `up` and `down` (see below). If this column is present, every row must have a value. This column is **auto-detected** as not being named `accession` or `sample`.|

To color by **presence** of accessions, do not include a value column.

With a value column, no two rows may name the same accession in the same sample. Combine repeated measurements, or the separate ions of the same molecule, in whichever way suits your data.

A compound element on a map is defined by one or more KEGG compounds (~2%% of elements have multiple). Since a single compound circle can stand for several compounds, its color is the aggregate of whichever of its accessions are in the file (`--compound-accession-aggregation`; see below).

## Coloring by presence or quantitative values

The inclusion of different columns results in different coloring modes.

**No sample column, no value column**

Compounds are colored a single presence color, pink by default — set it with `--compound-color`.

On global and overview maps, if reaction layer data supplied by %(kegg-reaction-txt)s is colored with `--original-color`, the compounds of those reactions are automatically colored to match, just as they appear in the reference maps.

**Sample column, no value column**

Compounds are colored by sample (or, using a %(groups-txt)s, group) count or membership. Counts are drawn either in discrete bands, one per count (`--compound-sample-summary count`), or as a gradient from the lowest count to the highest (`--compound-sample-summary count_continuous`), the latter being the only one that works when there are more samples than the colormap has distinguishable colors — without `--compound-sample-summary`, discrete automatically switches to continuous given enough samples. By default, the scale stops at the highest count actually observed rather than at the number of samples, so that the colors spread over the counts that occur; use `--count-scale-max` to change the maximum value setting.

**Value column**

Compounds are colored by the continuous value through a sequential colormap (`--compound-colormap`, default `plasma_r`). The colorbar is labeled by the value column's header. A cyclic colormap called `clocktime` was defined in anvi'o to show a value column of clock times in hours. Give clock times `--compound-value-period 24`: each color scale of values then runs from 0 to 24 h and is colored by `clocktime`. Each reduction of the values takes their circular mean, `circular_mean`, so that 23.5 h and 0.5 h average to 0 h, not 12 h. See %(anvi-draw-kegg-pathways)s.

A map element's constituent compound accessions are aggregated to a per-element value by `--compound-accession-aggregation` (`sum` by default) — this reduction happens within each sample.

With a `sample` column, `--draw-individual-files` and/or `--draw-grid` color maps of value column data that share a single `colorbar_compounds_samples.pdf` so that samples are comparable on the same scale. How the `unified` map, and with %(groups-txt)s the per-group maps, summarize samples is a separate choice (see below).

Each of these scales spans exactly the values it is given, which a few extreme elements can stretch until the rest cannot be told apart. `--compound-value-limits` bounds the `unified` map's scale and `--compound-category-value-limits` the scale its per-sample or per-group maps share. These parameters take an argument of two values, a minimum and a maximum, either of which can be `none` to leave that end open while setting the other. A limit sets its end of the scale whether or not any value reaches it. The colorbar labels each limit. A limit that values lie past is marked `≤` or `≥`, since its color then stands for that value or anything past it.

A scale also sits wherever its own values leave it, which for a signed quantity means the neutral middle color of a diverging colormap lands somewhere other than zero. `--compound-value-center` and `--compound-category-value-center` put a value at the middle of those two scales: used as a bare flag each centers its scale on 0, and each takes a number as an argument to center it on that instead. The scale is only ever widened to make room, running the same distance either side of the center and staying linear, so nothing is clipped that was not clipped already.

The two scales are colored by one colormap, `--compound-colormap`, unless `--compound-category-colormap` gives the per-sample or per-group scale one of its own. That is worth doing where the two show quantities of different kinds: with `--compound-sample-summary std` the `unified` map shows how much samples disagree rather than how large their values are, so drawing it in the colors that show magnitude on the per-sample maps invites mistaking one scale for the other.

## Coloring by categorical values

A value column can hold **categorical values** instead of numbers. These are names, such as `up`, `down` and `same`, or the names of clusters. Give `--compound-value-type categorical` to read the column this way. Values written as numbers, such as `1` and `2`, are names too. A value is read as text, so `NA` is a value like any other. Every row needs a value. Each compound can have one value per sample, so repeated rows are refused, as they are for numbers.

Each value takes a color of its own. By default, the values take colors in sorted order. Where every value is a number, such as `-1`, `0` and `1`, the values are sorted as numbers. Otherwise the digits in the values are compared as numbers, so `cluster 2` comes before `cluster 10`. Up to 9 values take the colors of Matplotlib's `tab10` without its gray. The 10th to 18th values take the light colors of `tab20` without its light gray. Gray is left out because global maps draw the compounds they do not highlight in gray. To choose the colors yourself, give a %(kegg-value-colors-txt)s file with `--compound-value-colors`. The order of its rows is then the order of the values. More than 18 values need such a file. The legend lists the values in order, with the first value at the top.

**Choosing one value where several meet.** A circle on a map occasionally stands for more than one compound, and these can have different values. `--compound-accession-aggregation` then takes a rule that chooses one of them. Each compound gives one vote. The rules are `plurality_order` (the default), `plurality_random`, `plurality`, `majority`, `order` and `random`, as described in %(kegg-reaction-txt)s. `plurality` and `majority` can give no value. The circle is then left uncolored. The order of the values is the row order of the colors file, or else sorted order. A random choice comes from `--seed` (42 by default). `--compound-reverse-overlay` puts the last value on top where circles overlap.

**Samples.** With a `sample` column, the map of each sample shows that sample's values. Draw these maps with `--draw-individual-files` and/or `--draw-grid`. They share one legend, `colorbar_compounds_samples.pdf`. A value in one sample need not mean the same thing in another sample. So by default the `unified` map shows presence across the samples. It counts the samples that contain each element, however few samples there are. Coloring by which samples contain an element would use the `tab10` colors of the values. A color would then mean a sample on one map and a value on another. `--compound-sample-summary membership` still asks for it. Where a value does mean the same thing in every sample, give a rule to `--compound-sample-summary`. The rule then chooses one value per compound across the samples, and the `unified` map shows it.

**Groups of samples.** With a %(groups-txt)s, maps are drawn for groups rather than samples. By default, each group's map counts the group's samples, and the `unified` map counts the groups. Such a run colors nothing by the values, and anvi'o warns about it. A rule given to `--compound-sample-summary` chooses one value per element among the samples of each group. Each group's map then shows that value. A rule given to `--compound-group-summary` then chooses one value across the groups for the `unified` map. It needs a rule for the samples, since the groups otherwise have no values. Each sample with a value gives one vote in its group. Each group with a value gives one vote across the groups, whatever its size. With few voters, ties are common. With two groups, every disagreement between them is a tie, so `plurality_order` settles it by the order of the values. `majority` among the samples of a group shows a value only where more than half of the samples with a value agree. A sample without a value does not count. So one sample with a value is enough. Where the `unified` map counts the groups, `--group-threshold` decides when a group contains an element, as for any layer.

The options that shape a scale of numbers are refused for categorical values. These are the value limits and centers, `--compound-value-period`, `--compound-element-normalization` and `--compound-category-colormap`. `--compound-colormap` colors only the presence on the `unified` map.

## Summarizing across samples and groups

There are three distinct reductions, each with its own option, plus one rescaling that is not a reduction. The first reduction applies whenever the file has a value column; everything below it needs a `sample` column:

|Level|Option|What it does|
|:--|:--|:--|
|accessions of a map element|`--compound-accession-aggregation` (an aggregation; `sum` by default; or a rule for categorical values)|the constituent accessions of a map element → that element's value|
|across samples|`--compound-sample-summary` (`count`, `count_continuous`, `membership`, an aggregation, or a rule for categorical values)|a set of samples → one value or one presence value per map element|
|across groups|`--compound-group-summary` (`count`, `count_continuous`, `membership`, an aggregation, or a rule for categorical values)|the groups of a %(groups-txt)s → one value or one presence value per map element|
|each sample's or group's value for an element|`--compound-element-normalization` (a normalization; none by default)|an element's value in one sample or group → that value rescaled against the element's values across all samples or groups|

An **aggregation** is `sum` (default), `mean`, `max`, `min`, `median`, `std` — or any other unsuggested pandas aggregation that reduces a series of numbers to one number, such as `var` or `sem`. Names that transform rather than reduce (`cumsum`) or that only a grouping offers (`first`) are rejected. Where an aggregation is undefined for the values available, as `std` is for a single value, the affected elements are left uncolored. `circular_mean` averages values that repeat after a period, such as clock times, and requires `--compound-value-period`. `angular_deviation` measures how far such values are spread apart, on the `unified` map only.

A summary pools the values of a **map element**, not of an accession. An element's value in a sample is the aggregate of its accessions in that sample. Without `--compound-element-normalization`, this is the value that the sample's own map draws. So the sample summary pools what the maps of the individual samples show. With groups, the group summary then pools the values of the groups.

The sample summary drives the `unified` map when there are no groups, and each per-group map (produced when using `--draw-individual-files` or `--draw-grid`) when there are. The group summary drives the `unified` map when there are groups.

Note that with groups, the per-group map either shows the count of that group's sample — the default case — or, with a `--compound-sample-summary` aggregation argument, the samples' pooled value; `--compound-sample-summary` is therefore not permitted with a presence name (`count`, `count_continuous` or `membership`) when using groups. For categorical values, a rule given to `--compound-sample-summary` colors each group's map by the value that it chooses among the group's samples.

The sample and group summaries in the `unified` maps both default to **presence** (`membership` for 3 or fewer categories, `count` above that, and `count_continuous` where there are more categories than the colormap has distinguishable colors), even with a quantitative value column, because presence is meaningful across any set of samples while pooling values is only meaningful for commensurable ones. For example, it makes sense to average metabolite concentrations across samples from replicate conditions, but not from unrelated conditions. Pooling values across samples is triggered by providing an aggregation argument to `--compound-sample-summary` or `--compound-group-summary`. `std` at this level maps how much samples disagree. A layer of categorical values counts the samples by default, however few there are.

A layer can show **presence in the overview and values per sample or group**. With the default presence summary, the `unified` map gets a colorbar of sample or group counts or memberships, while per-sample maps keep a continuous colorbar. The sample and group levels can be treated independently — replicate samples can be averaged within each condition (group) while the `unified` map summarizing all conditions shows in how many conditions each compound occurs, with a group containing a compound when at least half of its samples do:

{{ codestart }}
anvi-draw-kegg-pathways --compound-txt kegg-compound.txt \
                        --groups-txt %(groups-txt)s \
                        --group-threshold 0.5 \
                        --compound-sample-summary mean \
                        --compound-group-summary count \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

## Normalizing an element across samples

With a `sample` column and a value column, each compound element in a sample map is colored by value aggregated across accessions in the element. `--compound-element-normalization` instead colors each element relative to all of the samples. Each element is rescaled against its own values across the samples, so under `relative_to_mean` a value of +0.25 means the sample has 25%% more of that element than it averages across samples, and -0.31 means 31%% less.

This cannot be worked out in the file itself before analyzing the maps, since each element in each map is related to one or more accessions — the same accession factors into different rescaled element values on each map.

anvi'o encodes a number of normalizations with names recognized as an argument: `relative_to_mean`, `relative_to_median`, `difference_from_mean`, `difference_from_median`, `log2_ratio_to_mean`, `log2_ratio_to_median`, `z_score`, `rank`, `fraction_of_max`, and `fraction_of_total`. With `--compound-value-period`, the only normalization is `difference_from_circular_mean`, the signed offset from the element's circular mean. A different name is taken to be a pandas Series method that transforms each value into a new value, such as `abs`.

Only the maps of individual samples or groups are rescaled, so element normalization needs `--draw-individual-files` and/or `--draw-grid`. The `unified` map summarizes the unnormalized values, so is unaffected by normalization.

{{ codestart }}
anvi-draw-kegg-pathways --compound-txt kegg-compound.txt \
                        --compound-element-normalization relative_to_mean \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

## Maps with compounds that aren't colored

A small number of maps (`00121`, `00621`, `01052`, `01054`) contain compounds drawn as rectangles rather than circles. The rectangles overlay chemical structures on the base maps, so anvi'o does not color them in order to leave the structures unobscured.

## Examples

### Metabolomics (quantitative)

Color compound circles by measured concentrations. The `concentration` column is auto-detected as the value column.

|accession|concentration|
|:--|:--|
|C00031|55.0|
|C00092|12.0|
|C00074|30.5|

{{ codestart }}
anvi-draw-kegg-pathways --compound-txt kegg-compound.txt \
                        -o output_dir
{{ codestop }}

### Compound presence (single color)

With no value column and no `sample` column, present compounds are drawn one color.

|accession|
|:--|
|C00031|
|C00092|
|C00074|

{{ codestart }}
anvi-draw-kegg-pathways --compound-txt kegg-compound.txt \
                        -o output_dir
{{ codestop }}

### Comparing samples

Add a `sample` column to compare metabolomes (or other origins) on the same maps. Each sample's map is colored by its own values, while `--compound-sample-summary`/`--compound-group-summary` choose whether the summary views show presence or pooled values. See %(anvi-draw-kegg-pathways)s for the full set of sample, group, and colormap options.

|accession|sample|concentration|
|:--|:--|:--|
|C00031|metabolome_A|55.0|
|C00031|metabolome_B|31.0|
|C00074|metabolome_A|30.5|

### With a reaction layer

Pair this file with a %(kegg-reaction-txt)s file to color both layers on the same maps, each by its own mode. For example, quantitative metabolomics of compounds can be shown alongside the presence of KOs in a genome.

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt %(kegg-reaction-txt)s \
                        --compound-txt kegg-compound.txt \
                        -o output_dir
{{ codestop }}
