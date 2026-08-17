This artifact is a TAB-delimited text file that can be used to color the **reaction layer** of %(kegg-pathway-map)s files drawn by %(anvi-draw-kegg-pathways)s.

The "reaction layer" contains reaction lines in global and overview maps and boxes in standard maps. This file, passed to %(anvi-draw-kegg-pathways)s with the `--reaction-txt` option, can be paired with a %(kegg-compound-txt)s file (`--compound-txt`) to color the reaction and compound layers of a map at once. One or more %(contigs-db)ss, a %(pan-db)s, or a metabolic model file (see %(reaction-network-json)s) can alternatively be used by the program to color the reaction layer.

## File format

The file must have a header line with a required `accession` column; the other columns are optional.

|Column|Required|Description|
|:--|:--|:--|
|`accession`|YES|A KEGG ID. Accessions in the file must **all** be **KEGG Ortholog (KO) IDs** (`K#####`) or **KEGG reaction IDs** (`R#####`) — a single file cannot mix the two. Both KOs and Reactions address the reaction layer, so combining them can clash with two scales on the same map elements.|
|`gene_id`|NO|An identifier for the gene a value came from. Several genes may carry the same accession, and their values are combined into that accession's value (by default their `sum` — see `--reaction-gene-aggregation`).|
|`sample`|NO|A sample of origin, enabling comparison across samples (and, with a %(groups-txt)s, across groups of samples). If this column is present, every row must have a sample.|
|*a value column*|NO|A SINGLE numeric column with any name of your choosing, for coloring by a continuous quantitative value instead of by presence/absence. Its meaning is up to you (coverage, expression, flux, ...). If this column is present, every row must have a value. This column is **auto-detected** as not being named `accession`, `gene_id`, or `sample`.|

To color by **presence** of accessions, do not include a value column.

With a value column, no two rows may describe the same thing — the same accession, gene, and sample. Repeated rows are ambiguous, as they can be replicate measurements to average or separate contributions to add, so combine them yourself in whichever way suits your data. Note that several *different* genes carrying one accession are not repeats, and are the ordinary case of the same type of protein being encoded by different genes.

A reaction element on a map is defined by one or more KEGG reactions (~23%% of elements have multiple), which may be carried out by one or more KOs. A file of **KO** accessions addresses the reaction layer by KO ID, while a file of **KEGG reaction** accessions addresses it by reaction ID — useful when your values are per-reaction, as with metabolic model fluxes. Since a single line or box can stand for several KOs or reactions, its color is the aggregate of whichever of its accessions are in the file (`--reaction-accession-aggregation`; see below).

## Coloring by presence or quantitative values

The inclusion of different columns results in different coloring modes.

**No sample column, no value column**

Reactions are colored a single presence color, green by default — set it with `--reaction-color`, or use `--original-color` for the reference map's colors.

**Sample column, no value column**

Reactions are colored by sample (or, with a %(groups-txt)s, group) count or membership, exactly as when comparing multiple %(contigs-db)ss or genomes of a %(pan-db)s. Counts are drawn either in discrete bands, one per count (`--reaction-sample-summary count`), or as a gradient from the lowest count to the highest (`--reaction-sample-summary count_continuous`), the latter being the only one that works when there are more samples than the colormap has distinguishable colors — without `--reaction-sample-summary`, discrete automatically switches to continuous given enough samples. By default, the scale stops at the highest count actually observed rather than at the number of samples, so that the colors spread over the counts that occur; use `--count-scale-max` to change the maximum value setting.

With `--reaction-color` or `--original-color`, the single color or original colors of the reference map override dynamic coloring, simply showing whether a reaction is present in any of the samples.

**Value column**

Reactions are colored by the continuous value through a sequential colormap (`--reaction-colormap`, default `plasma_r`). The colorbar is labeled by the value column's header.

Per-gene values are aggregated to a per-accession value, and a map element's constituent KO and reaction accessions to a per-element value, by `--reaction-gene-aggregation` and `--reaction-accession-aggregation` (both `sum` by default) — these reductions happen within each sample.

With a `sample` column, `--draw-individual-files` and/or `--draw-grid` color maps of value column data that share a single `colorbar_reactions_samples.pdf` so that samples are comparable on the same scale. How the `unified` map, and with %(groups-txt)s the per-group maps, summarize samples is a separate choice (see below).

Each of these scales spans exactly the values it is given, which a few extreme elements can stretch until the rest cannot be told apart. `--reaction-value-limits` bounds the `unified` map's scale and `--reaction-category-value-limits` the scale its per-sample or per-group maps share. These parameters take an argument of two values, a minimum and a maximum, either of which can be `none` to leave that end open while truncating the other. A limit only truncates where the values actually cross it, and a truncated end of the colorbar is marked `≤` or `≥`, since its color then stands for that value or anything past it.

A scale also sits wherever its own values leave it, which for a signed quantity means the neutral middle color of a diverging colormap lands somewhere other than zero. `--reaction-value-center` and `--reaction-category-value-center` put a value at the middle of those two scales: used as a bare flag each centers its scale on 0, and each takes a number as an argument to center it on that instead. The scale is only ever widened to make room, running the same distance either side of the center and staying linear, so nothing is clipped that was not clipped already.

The two scales are colored by one colormap, `--reaction-colormap`, unless `--reaction-category-colormap` gives the per-sample or per-group scale one of its own. That is worth doing where the two show quantities of different kinds: with `--reaction-sample-summary std` the `unified` map shows how much samples disagree rather than how large their values are, so drawing it in the colors that show magnitude on the per-sample maps invites mistaking one scale for the other.

## Summarizing across samples and groups

There are four distinct reductions, each with its own option. The first two apply whenever the file has a value column; the last two need a `sample` column:

|Level|Option|What it reduces|
|:--|:--|:--|
|genes of an accession|`--reaction-gene-aggregation` (an aggregation; `sum` by default)|the genes carrying an accession → that accession's value|
|accessions of a map element|`--reaction-accession-aggregation` (an aggregation; `sum` by default)|the constituent accessions of a map element → that element's value|
|across samples|`--reaction-sample-summary` (`count`, `count_continuous`, `membership`, or an aggregation)|a set of samples → one continuous value or one presence value per accession|
|across groups|`--reaction-group-summary` (`count`, `count_continuous`, `membership`, or an aggregation)|the groups of a %(groups-txt)s → one continuous value or one presence value per accession|

An **aggregation** is `sum` (default), `mean`, `max`, `min`, `median`, `std` — or any other unsuggested pandas aggregation that reduces a series of numbers to one number, such as `var` or `sem`. Names that transform rather than reduce (`cumsum`) or that only a grouping offers (`first`) are rejected. Where an aggregation is undefined for the values available, as `std` is for a single value, the affected elements are left uncolored and a warning says how many accessions are affected.

A summary reduces only the samples (or groups) that actually contain an accession: a sample with no row for an accession is treated as not having observed it, so `mean` over three samples where only one lists the accession is that one sample's value.

The sample summary drives the `unified` map when there are no groups, and each per-group map (produced when using `--draw-individual-files` or `--draw-grid`) when there are. The group summary drives the `unified` map when there are groups.

Note that with groups, the per-group map either shows the count of that group's sample — the default case — or, with a `--reaction-sample-summary` aggregation argument, the samples' pooled value; `--reaction-sample-summary` is therefore not permitted with a presence name (`count`, `count_continuous` or `membership`) when using groups.

The sample and group summaries in the `unified` maps both default to **presence** (`membership` for 3 or fewer categories, `count` above that, and `count_continuous` where there are more categories than the colormap has distinguishable colors), even with a quantitative value column, because presence is meaningful across any set of samples while pooling values is only meaningful for commensurable ones. For example, it makes sense to average transcript abundances across samples from replicate conditions, but not from unrelated conditions. Pooling values across samples is triggered by providing an aggregation argument to `--reaction-sample-summary` or `--reaction-group-summary`. `std` at this level maps how much samples disagree.

A layer can be **categorical in the overview and continuous per sample or group**. With the default presence summary, the `unified` map gets a colorbar of sample or group counts or memberships, while per-sample maps keep a continuous colorbar. The sample and group levels can be treated independently — replicate samples can be averaged within each condition (group) while the `unified` map summarizing all conditions shows in how many conditions each reaction occurs, with a group containing a reaction when at least half of its samples do:

{{ codestart }}
anvi-draw-kegg-pathways --reaction-txt kegg-reaction.txt \
                        --groups-txt %(groups-txt)s \
                        --group-threshold 0.5 \
                        --reaction-sample-summary mean \
                        --reaction-group-summary count \
                        --draw-individual-files \
                        -o output_dir
{{ codestop }}

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
