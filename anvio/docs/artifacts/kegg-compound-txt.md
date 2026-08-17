This artifact is a TAB-delimited text file that can be used to color the **compound layer** of %(kegg-pathway-map)s files drawn by %(anvi-draw-kegg-pathways)s.

The "compound layer" contains compound circles in maps. This file, passed to %(anvi-draw-kegg-pathways)s with the `--compound-txt` option, can be paired with a %(kegg-reaction-txt)s file (`--reaction-txt`) to color the compound and reaction layers of a map at once.

## File format

The file must have a header line with a required `accession` column; the other columns are optional.

|Column|Required|Description|
|:--|:--|:--|
|`accession`|YES|A KEGG compound ID (`C#####`).|
|`sample`|NO|A sample of origin, enabling comparison across samples (and, with a %(groups-txt)s, across groups of samples). If this column is present, every row must have a sample.|
|*a value column*|NO|A SINGLE numeric column with any name of your choosing, for coloring by a continuous quantitative value instead of by presence/absence. Its meaning is up to you (concentration, intensity, ...). If this column is present, every row must have a value. This column is **auto-detected** as not being named `accession` or `sample`.|

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

Compounds are colored by the continuous value through a sequential colormap (`--compound-colormap`, default `plasma_r`). The colorbar is labeled by the value column's header.

A map element's constituent compound accessions are aggregated to a per-element value by `--reaction-accession-aggregation` (`sum` by default) — this reduction happens within each sample.

With a `sample` column, `--draw-individual-files` and/or `--draw-grid` color maps of value column data that share a single `colorbar_compounds_samples.pdf` so that samples are comparable on the same scale. How the `unified` map, and with %(groups-txt)s the per-group maps, summarize samples is a separate choice (see below).

Each of these scales spans exactly the values it is given, which a few extreme elements can stretch until the rest cannot be told apart. `--compound-value-limits` bounds the `unified` map's scale and `--compound-category-value-limits` the scale its per-sample or per-group maps share. These parameters take an argument of two values, a minimum and a maximum, either of which can be `none` to leave that end open while truncating the other. A limit only truncates where the values actually cross it, and a truncated end of the colorbar is marked `≤` or `≥`, since its color then stands for that value or anything past it.

A scale also sits wherever its own values leave it, which for a signed quantity means the neutral middle color of a diverging colormap lands somewhere other than zero. `--compound-value-center` and `--compound-category-value-center` put a value at the middle of those two scales: used as a bare flag each centers its scale on 0, and each takes a number as an argument to center it on that instead. The scale is only ever widened to make room, running the same distance either side of the center and staying linear, so nothing is clipped that was not clipped already.

The two scales are colored by one colormap, `--compound-colormap`, unless `--compound-category-colormap` gives the per-sample or per-group scale one of its own. That is worth doing where the two show quantities of different kinds: with `--compound-sample-summary std` the `unified` map shows how much samples disagree rather than how large their values are, so drawing it in the colors that show magnitude on the per-sample maps invites mistaking one scale for the other.

## Summarizing across samples and groups

There are three distinct reductions, each with its own option. The first applies whenever the file has a value column; the last two need a `sample` column:

|Level|Option|What it reduces|
|:--|:--|:--|
|accessions of a map element|`--compound-accession-aggregation` (an aggregation; `sum` by default)|the constituent accessions of a map element → that element's value|
|across samples|`--compound-sample-summary` (`count`, `count_continuous`, `membership`, or an aggregation)|a set of samples → one continuous value or one presence value per accession|
|across groups|`--compound-group-summary` (`count`, `count_continuous`, `membership`, or an aggregation)|the groups of a %(groups-txt)s → one continuous value or one presence value per accession|

An **aggregation** is `sum` (default), `mean`, `max`, `min`, `median`, `std` — or any other unsuggested pandas aggregation that reduces a series of numbers to one number, such as `var` or `sem`. Names that transform rather than reduce (`cumsum`) or that only a grouping offers (`first`) are rejected. Where an aggregation is undefined for the values available, as `std` is for a single value, the affected elements are left uncolored and a warning says how many accessions are affected.

A summary reduces only the samples (or groups) that actually contain an accession: a sample with no row for an accession is treated as not having observed it, so `mean` over three samples where only one lists the accession is that one sample's value.

The sample summary drives the `unified` map when there are no groups, and each per-group map (produced when using `--draw-individual-files` or `--draw-grid`) when there are. The group summary drives the `unified` map when there are groups.

Note that with groups, the per-group map either shows the count of that group's sample — the default case — or, with a `--compound-sample-summary` aggregation argument, the samples' pooled value; `--compound-sample-summary` is therefore not permitted with a presence name (`count`, `count_continuous` or `membership`) when using groups.

The sample and group summaries in the `unified` maps both default to **presence** (`membership` for 3 or fewer categories, `count` above that, and `count_continuous` where there are more categories than the colormap has distinguishable colors), even with a quantitative value column, because presence is meaningful across any set of samples while pooling values is only meaningful for commensurable ones. For example, it makes sense to average metabolite concentrations across samples from replicate conditions, but not from unrelated conditions. Pooling values across samples is triggered by providing an aggregation argument to `--compound-sample-summary` or `--compound-group-summary`. `std` at this level maps how much samples disagree.

A layer can be **categorical in the overview and continuous per sample or group**. With the default presence summary, the `unified` map gets a colorbar of sample or group counts or memberships, while per-sample maps keep a continuous colorbar. The sample and group levels can be treated independently — replicate samples can be averaged within each condition (group) while the `unified` map summarizing all conditions shows in how many conditions each compound occurs, with a group containing a compound when at least half of its samples do:

{{ codestart }}
anvi-draw-kegg-pathways --compound-txt kegg-compound.txt \
                        --groups-txt %(groups-txt)s \
                        --group-threshold 0.5 \
                        --compound-sample-summary mean \
                        --compound-group-summary count \
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
