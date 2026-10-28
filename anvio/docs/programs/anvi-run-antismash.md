This program runs [antiSMASH](https://antismash.secondarymetabolites.org/) on a %(contigs-db)s to identify **biosynthetic gene clusters** (BGCs) that, for example, produce secondary metabolites such as antibiotics, siderophores, pigments, and quorum-sensing signals. The program stores antiSMASH's annotations in the database as %(functions)s.

{:.notice}
Just to avoid any misunderstandings and/or confusions down the line: what a gene cluster describes in this context is very different than what it describes in the context of pangenomics.

For downstream analysis, you can export the annotations with %(anvi-export-functions)s, search them with %(anvi-search-functions)s, view them in %(anvi-interactive)s, include them in %(anvi-summarize)s, and compare BGCs across different contigs databases with %(anvi-script-gen-antismash-matrix)s.

antiSMASH describes several different facets of each gene in a cluster (what it does, its role, the biosynthetic domains it carries, and which cluster it belongs to). The anvi'o functions table stores a single accession and function per row, so each of these facets becomes its own annotation source:

|**source**|**what it holds**|
|:--|:--|
|`antiSMASH`|the clearest functional name for each gene: antiSMASH's smCOG description where available (e.g. `crotonyl-CoA reductase / alcohol dehydrogenase`, with the smCOG id as the accession and its E-value), otherwise the biosynthetic domain name. A gene with several names gets several rows.|
|`antiSMASH_ROLE`|the gene's role in the cluster: `biosynthetic`, `biosynthetic-additional`, `regulatory`, `transport`, or `other`.|
|`antiSMASH_DOMAIN`|all the biosynthetic domains an individual gene in the BGC carries, each with its E-value.|
|`antiSMASH_REGION`|which BGC region the gene belongs to (the accession is a compact, assembly-unique region id; the function is the region's product/type, e.g. `terpene` or `NRPS`). Genes that sit on a region boundary count as part of that region.|

All four sources are always registered in the database after a run, even if antiSMASH found nothing for some (or all) of them, so you can always tell a database that was searched from one that was never processed.

## How does it work?

**Export the sequences and gene calls**
antiSMASH uses both sequence information and gene locations, so this program exports a FASTA file and a GFF3 file of the gene calls already stored in your %(contigs-db)s. At this point, very short contigs are filtered out before the BGC search (by default, anything shorter than 3,000 bp), because antiSMASH cannot meaningfully call clusters on short fragments.

**Gene and BGC annotation**
antiSMASH scans the sequences for the signature enzymes and domains of known biosynthetic pathways and groups them into **regions** (the clusters). The primary annotation source is **smCOGs** (secondary-metabolism Clusters of Orthologous Groups), a curated library of the protein families that recur in biosynthetic gene clusters. In contrast to COGs, smCOGs cover only the enzymes, transporters, and regulators found in and around BGCs. Each one (for example `SMCOG1006`, acyl-CoA dehydrogenase) carries a readable functional name and a hint about the gene's likely role in the cluster, and `anvi-run-antismash` uses the smCOG description as a gene's clearest name whenever one exists. Additionally, antiSMASH draws on the **TIGRFAM** and **Pfam** protein-family databases, uses **RRE-Finder** to identify RiPP-recognition elements (the signatures of RiPPs, ribosomally synthesized and post-translationally modified peptides), and runs an **Active Site Finder**. The per-gene results are imported into your %(contigs-db)s as the four sources above: the TIGRFAM, Pfam, and RRE domains are stored in the `antiSMASH_DOMAIN` source, while the active-site predictions are only visible in the browsable HTML report.

## Setting up antiSMASH

antiSMASH has to be installed in a **separate** conda environment, because its strict dependency versions conflict with anvi'o's. This is a one-time setup:

{{ codestart }}
conda create -y -n antismash python=3.10
conda activate antismash
conda install -y -c conda-forge -c bioconda antismash=7.1.0 "bcbio-gff<0.7"
download-antismash-databases
conda deactivate
{{ codestop }}

{:.notice}
The `download-antismash-databases` step fetches several gigabytes of reference data, so it can take a while and needs disk space. You only ever have to do it once.

{:.notice}
The `bcbio-gff<0.7` pin is a workaround for a current incompatibility: newer `bcbio-gff` (0.7+) needs a `biopython` version that antiSMASH 7.1.0 does not allow, and without the pin antiSMASH fails to parse the gene calls (with a `could not parse records from GFF3 file` error). You can drop it once a future antiSMASH release resolves this.

## Standard usage

To annotate a single %(contigs-db)s with BGCs, run:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s
{{ codestop }}

To annotate **many databases in one command**, list them in an %(external-genomes)s file and pass it with `-e` instead of `-c` (just like %(anvi-estimate-metabolism)s):

{{ codestart }}
anvi-run-antismash -e %(external-genomes)s \
                   --num-threads 8
{{ codestop }}

Each database in the file is annotated in place and gets its own report directory (`<db>-ANTISMASH`) next to it. Databases that already carry antiSMASH annotations are skipped, so you can safely re-run to pick up where a batch left off (add `--just-do-it` to re-annotate them anyway). The single-database output options (`--output-dir`, `--regions-output`, `--work-dir`) apply to `-c` runs only.

You normally do **not** need to tell anvi'o where antiSMASH is. Because antiSMASH lives in its own conda environment, anvi'o first looks for an `antismash` executable on your `PATH`, and failing that automatically searches your other conda environments for it.

## Telling anvi'o where antiSMASH is

If anvi'o cannot find antiSMASH on its own, or if you have several antiSMASH installations and want to use a specific one, point it at the antiSMASH installation directory (i.e. its conda environment) with `--antismash-dir`:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s \
                   --antismash-dir ~/miniforge3/envs/antismash
{{ codestop }}

## Extended output option

`anvi-run-antismash` annotates the %(contigs-db)s as described above and saves a browsable antiSMASH HTML report together with the two summary files. By default it keeps the view-only HTML report and the summary files, but does not generate the per-region GenBank files, the raw JSON, or the ClusterBlast/smCOG comparison data.

Adding `--include-detailed-output` runs antiSMASH's cluster-comparison analyses and keeps their output, together with the per-region GenBank (`.gbk`) files (a useful input to downstream tools such as [BiG-SCAPE](https://bigscape-corason.secondarymetabolites.org/)) and the raw JSON:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s \
                   --include-detailed-output
{{ codestop }}

The comparison analyses this enables are:

* **`--cb-general`** (ClusterBlast): compares each cluster against a large database of antiSMASH-predicted clusters from other genomes.
* **`--cb-knownclusters`** (KnownClusterBlast): compares against experimentally characterized clusters in [MIBiG](https://mibig.secondarymetabolites.org/) (i.e. clusters with a known product).
* **`--cb-subclusters`** (SubClusterBlast): compares against known sub-clusters that make specific precursor building blocks.
* **`--cc-mibig`**: a cluster-level comparison against the MIBiG reference set.
* **`--smcog-trees`**: builds phylogenetic trees for the smCOG gene families in each cluster.

The database annotations are unaffected by these additional outputs.

{:.notice}
The detailed output is much larger and slower to produce: on one test genome the default run produced about 2.4 MB in roughly 50 seconds, versus about 91 MB in roughly 4 minutes with `--include-detailed-output`.

## Optional flags

**1. Change the output directory.** By default the report is saved next to the %(contigs-db)s as `<contigs-db basename>-ANTISMASH`. Use `-o` / `--output-dir` to put it somewhere else:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s \
                   --output-dir MY_BGC_REPORT/
{{ codestop }}

**2. Use multiple threads.** antiSMASH can take a while, so if you have the cores to spare, give it more threads:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s -T 4
{{ codestop }}

**3. Switch the taxon to fungi.** By default anvi'o tells antiSMASH to expect a bacterial genome. If your %(contigs-db)s describes a fungal genome, switch the taxon:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s \
                   --taxon fungi
{{ codestop }}

**4. Change the minimum contig length.** Contigs shorter than 3,000 bp are dropped before the run. To change that cutoff, for instance to be stricter, set `--min-contig-length`:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s \
                   --min-contig-length 5000
{{ codestop }}

## Additional options

antiSMASH provides a series of additional annotation tools that you can also include in an `anvi-run-antismash` run:

* **`--fullhmmer`**: whole-genome Pfam annotation (essentially what %(anvi-run-pfams)s does).
* **`--clusterhmmer`**: Pfam annotation limited to genes inside clusters.
* **`--tigrfam`**: TIGRFAM family annotation of cluster genes (antiSMASH's core detection already contributes the relevant ones).
* **`--pfam2go`**: maps Pfam domains to Gene Ontology terms (report only).
* **`--tfbs`**: TFBS finder that looks for transcription-factor binding sites in clusters (regulatory context).
* **`--cassis`**: an alternative, motif-based method for predicting cluster boundaries.

To include any of them, pass them through the `--antismash-add-on` flag. Please use the `=` form shown below, otherwise the command line parser may mistake the antiSMASH flags for flags of `anvi-run-antismash` itself:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s \
                   --antismash-add-on="--fullhmmer --tfbs"
{{ codestop }}

Flags anvi'o already sets for antiSMASH on its own (such as `--cpus`, `--taxon`, `--databases`, `--output-dir`, or the ones `--include-detailed-output` turns on) are not accepted here; please use the corresponding `anvi-run-antismash` flags for those instead.

## Re-annotating a database

If the %(contigs-db)s already contains antiSMASH annotations from a previous run, the program refuses to overwrite them, so you don't clobber earlier results by accident. To re-annotate, either remove the old annotations first:

{{ codestart }}
anvi-delete-functions -c %(contigs-db)s \
                      --annotation-sources antiSMASH,antiSMASH_ROLE,antiSMASH_DOMAIN,antiSMASH_REGION
anvi-run-antismash -c %(contigs-db)s
{{ codestop }}

or pass `--just-do-it` to overwrite them in a single step:

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s --just-do-it
{{ codestop }}

With `--just-do-it`, anvi'o first removes all four antiSMASH sources from the previous run, so nothing stale survives even if the new run finds fewer (or no) clusters.

{:.notice}
Like other anvi'o programs, `anvi-run-antismash` will not write into an output directory that already exists. So if the report directory from the previous run (the default `<contigs-db basename>-ANTISMASH`, or the one you gave with `--output-dir`) is still around, either use `--output-dir` to put the new report elsewhere, or add `--force-overwrite` to have anvi'o remove the existing directory first (please use that flag with care, since it removes the directory with everything in it). The same applies to a directory you give with `--work-dir`.

{{ codestart }}
anvi-run-antismash -c %(contigs-db)s --just-do-it --force-overwrite
{{ codestop }}

## What you get

**Four annotation sources in the database.** These are the four sources described at the top of this page. You can confirm they were added with %(anvi-db-info)s.

**A browsable HTML report,** saved by default as `<contigs-db basename>-ANTISMASH`. Open its `index.html` in a web browser to explore each region interactively, exactly as you would on the antiSMASH website. Please note that contigs in the HTML report are named `anvio_c000000`, `anvio_c000001`, and so on, rather than with your contig names (see below for why, and for how to translate them back).

**A contig name mapping file,** `contig-name-mapping.txt`, inside the report directory. It translates the contig names you see in antiSMASH's own output back to the contig names in your %(contigs-db)s:

|**antismash_record_id**|**contig**|
|:--|:--|
|anvio_c000000|c_000000000001|
|anvio_c000001|c_000000000002|

Why is this necessary? antiSMASH renames any sequence whose name is longer than 16 characters (e.g., `NODE_1_length_10000_cov_5.2` becomes something like `c00001_NODE_1_..`), and once that happens anvi'o could no longer tell which of your contigs a given BGC came from. To make sure every result can be linked back to the right contig regardless of how long your contig names are (or how future versions of antiSMASH decide to handle them), anvi'o gives every contig a short, predictable name before handing it to antiSMASH, and translates these names back when it imports the results. So, the annotations in your database and the two summary files below always use your own contig names, and the only places the `anvio_c` names show up are the files antiSMASH itself generates: the HTML report and, with `--include-detailed-output`, the GenBank (`.gbk`) and JSON files (keep this in mind if you use the `.gbk` files with other tools such as BiG-SCAPE). Whenever you need to know which contig an `anvio_c` name refers to, `contig-name-mapping.txt` is where to look.

**Two self-contained TAB-delimited (`.txt`) summary files inside the report directory.** These lay everything out for direct inspection, so you never have to go back to the database to read the results.

`antismash_genes.txt` has one row per antiSMASH-annotated gene, with the function name, biosynthetic domains, smCOG, role, region, and E-value as separate columns:

|**gene_callers_id**|**contig**|**region_id**|**role**|**function**|**smcog**|**domains**|**e_value**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|698|c_000000000001|e5c0ecfbd09e_region_001|biosynthetic-additional|crotonyl-CoA reductase / alcohol dehydrogenase|SMCOG1028|ADH_N;ADH_zinc_N|7.4e-85|
|700|c_000000000001|e5c0ecfbd09e_region_001|regulatory|autoinducer-binding transcriptional regulator|SMCOG1197||2.3e-30|
|701|c_000000000001|e5c0ecfbd09e_region_001|biosynthetic|Autoind_synth||Autoind_synth|1.8e-14|

`antismash_regions.txt` has one row per BGC region, fully describing it (its compact id, contig, product/type, coordinates, and the genes it contains):

|**region_id**|**contig**|**product**|**start**|**end**|**length**|**num_genes**|**gene_callers_ids**|**contig_edge**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|e5c0ecfbd09e_region_001|c_000000000001|hserlactone|690123|711456|21333|22|689;690;691;...;710|False|

## Inspecting and using the results

Confirm the sources were added, export them, or search them:

{{ codestart }}
anvi-db-info %(contigs-db)s | grep antiSMASH
anvi-export-functions -c %(contigs-db)s --annotation-sources antiSMASH -o bgc.txt
anvi-search-functions -c %(contigs-db)s --search-terms NRPS,PKS,siderophore
{{ codestop }}

## Comparing BGCs across genomes

To compare biosynthetic potential across many genomes, first annotate them all (a single `anvi-run-antismash -e %(external-genomes)s` does the whole set, as shown above), then summarize them with %(anvi-script-gen-antismash-matrix)s. The matrix tool reads the four antiSMASH sources straight from each database, so there is nothing extra to run per database beyond the annotation itself.

Point it at the same %(external-genomes)s file and pick one of two views. A per-cluster matrix (`--regions`), where the rows are BGC product types such as `terpene` or `NRPS` and each cell counts how many clusters of that type a genome has:

{{ codestart }}
anvi-script-gen-antismash-matrix -e %(external-genomes)s \
                                 --regions \
                                 -O BGC_TYPES
{{ codestop }}

or a per-gene matrix (`--genes`), where the rows are gene signatures described by their function, role, smCOG, and domains:

{{ codestart }}
anvi-script-gen-antismash-matrix -e %(external-genomes)s \
                                 --genes \
                                 -O BGC_GENES
{{ codestop }}

Either way you get a frequency matrix and a presence/absence matrix across all of your genomes. See %(anvi-script-gen-antismash-matrix)s for the full description, including the companion per-gene file that links each matrix cell back to the individual genes behind it.
