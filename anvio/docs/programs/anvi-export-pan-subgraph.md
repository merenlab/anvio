This program exports the genomic loci that fall between two nodes of an anvi'o %(pan-graph-db)s as individual %(contigs-db)s files (one per genome).

You will typically identify the two nodes of interest by visualizing the pangenome graph with %(anvi-display-pan-graph)s, and then ask this program to extract, for each genome in the graph, the stretch of DNA that lies between the genes that correspond to those two nodes.

### Required Parameters

You will need to provide a %(pan-graph-db)s, an %(external-genomes)s file that describes the genomes that went into the graph (so anvi'o knows where to find each source %(contigs-db)s), an output directory, and a way to describe the loci to export. The loci can be described either by two graph nodes (`--graph-nodes`) or by a single region ID (`--region-id`). You must provide one of these two, but not both.

#### Describing loci with two graph nodes

Provide the two graph nodes that flank the loci you wish to export with `--graph-nodes`:

{{ codestart }}
anvi-export-pan-subgraph -p PAN-GRAPH.db \
                         -e %(external-genomes)s \
                         --graph-nodes NODE_1,NODE_2 \
                         -o OUTPUT_DIRECTORY
{{ codestop }}

#### Describing loci with a region ID

Alternatively, you can provide a single region ID with `--region-id`, and anvi'o will resolve the two boundary nodes of that region automatically:

{{ codestart }}
anvi-export-pan-subgraph -p PAN-GRAPH.db \
                         -e %(external-genomes)s \
                         --region-id 234 \
                         -o OUTPUT_DIRECTORY
{{ codestop }}

Please note that the way by which boundary nodes are resolved depends on the type of the region:

* For a **backbone region**, anvi'o uses the leftmost and rightmost nodes of the region (by position) as the boundaries.
* For a **variable region**, the first and the last nodes in the region will not necessarily present in all genomes, so anvi'o instead uses the closest eligible (non-RNA, coding) backbone nodes flanking the region on each side. As a result, the exported loci will include the variable region content plus at least two extra genes from the flanking backbone regions (one on each side), and anvi'o will print a warning describing exactly what was used. Keep this in mind if you intend to analyze the conservancy of genes in variable regions, to avoid misinterpretations.

You will typically learn which nodes or region IDs you are interested in by visualizing the pangenome graph in anvi'o using the program %(anvi-display-pan-graph)s. Region IDs will appear as you zoom in to certain parts of the pangenome graph like shown below:

[![Region IDs in pangenome graphs](../../images/anvi-export-pan-subgraph-region.png){:.center-img .width-70}](../../images/anvi-export-pan-subgraph-region.png)

#### Output

Regardless of how you describe the loci, the output directory will contain a separate %(contigs-db)s (along with a FASTA file and a blank %(profile-db)s) for the locus extracted from each genome.

### Gene caller IDs

By default, the gene caller IDs in each exported %(contigs-db)s **match those in the source %(contigs-db)s** that the locus was extracted from. This way you can trace every gene in an exported locus back to the gene it originated from in the original database.

If you would rather have the gene caller IDs reset so that they start from `0` in each output database (which was the historical behavior of this program), use the `--reset-gene-caller-ids` flag:

{{ codestart }}
anvi-export-pan-subgraph -p PAN-GRAPH.db \
                         -e %(external-genomes)s \
                         --graph-nodes NODE_1,NODE_2 \
                         -o OUTPUT_DIRECTORY \
                         --reset-gene-caller-ids
{{ codestop }}
