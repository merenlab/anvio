%(anvi-search-sequence-motifs)s will search one or more sequence motifs in applicable anvi'o databases and will report their frequency. If you have more than one motif to search, you can list them as comma-separated sequences

In this context we assume a motif is a 4 to 10 nucleotide-long string, although, anvi'o will not impose any limit to length, and will search any motif it is given along with its reverse-complement across all sequences and report frequencies.

The most primitive output is a TAB-delimited text file, but anvi'o will store frequency information also into your databases like a pro if you use the `--store-in-db` flag.

The following subsections include some examples.

## A contigs database

The minimum amount of stuff you need to run this program is a motif sequence and a %(contigs-db)s:

{{ codestart }}
anvi-search-sequence-motifs -c %(contigs-db)s \
                            --motifs ATCG,TAAAT \
                            --output-file motifs.txt
{{ codestop }}

Running this will yield an output file with as many columns as the number of sequence motifs that show their frequencies across each contig found in the %(contigs-db)s. Here is an example:

|contig_name|ATCG|TAAAT|
|:--|:--:|:--:|
|204_10M_contig_1720|101|159|
|204_10M_contig_6515|64|31|
|204_10M_contig_878|435|3|

## Contigs database + profile database

If you provide this program with a %(profile-db)s, this time it will count your motif sequences in split sequences rather than contigs,

{{ codestart }}
anvi-search-sequence-motifs -c %(contigs-db)s \
                            -p %(profile-db)s
                            --motifs ATCG,TAAAT \
                            --output-file motifs.txt
{{ codestop }}

And the output will look like this:

|split_name|ATCG|TAAAT|
|:--|:--:|:--:|
|204_10M_contig_1720_split_00001|14|22|
|204_10M_contig_1720_split_00002|2|6|
|204_10M_contig_1720_split_00003|14|23|
|204_10M_contig_1720_split_00004|8|18|
|204_10M_contig_1720_split_00005|9|17|
|204_10M_contig_1720_split_00006|19|28|
|204_10M_contig_1720_split_00007|4|8|
|204_10M_contig_1720_split_00008|31|32|
|204_10M_contig_1720_split_00009|0|5|
|204_10M_contig_6515_split_00001|7|5|
|204_10M_contig_6515_split_00002|5|2|
|204_10M_contig_6515_split_00003|5|4|
|204_10M_contig_6515_split_00004|25|8|
|204_10M_contig_6515_split_00005|6|2|
|204_10M_contig_6515_split_00006|8|3|
|204_10M_contig_6515_split_00007|3|3|
|204_10M_contig_6515_split_00008|5|3|
|204_10M_contig_878_split_00001|17|0|
|204_10M_contig_878_split_00002|14|0|
|204_10M_contig_878_split_00003|108|1|
|204_10M_contig_878_split_00004|35|0|
|204_10M_contig_878_split_00005|7|0|
|204_10M_contig_878_split_00006|18|0|
|204_10M_contig_878_split_00007|42|0|
|204_10M_contig_878_split_00008|12|1|
|204_10M_contig_878_split_00009|13|0|
|204_10M_contig_878_split_00010|18|0|
|204_10M_contig_878_split_00011|28|0|
|204_10M_contig_878_split_00012|0|1|
|204_10M_contig_878_split_00013|24|0|
|204_10M_contig_878_split_00014|11|0|
|204_10M_contig_878_split_00015|33|0|
|204_10M_contig_878_split_00016|13|0|
|204_10M_contig_878_split_00017|2|0|
|204_10M_contig_878_split_00018|40|0|

{:.notice}
This output format may enable you to bin your splits based on their motif composition and use %(anvi-import-collection)s to import them as a new collection into your profile database, or use %(anvi-matrix-to-newick)s to cluster them based on this information to organize splits in the interface based on their motif composition.

You can also store this information into your profile database using the flag `--store-in-db`. When you do that, running %(anvi-interactive)s on this profile database will include additional layers where these frequencies are displayed. Here is an example:

{{ codestart }}
%(anvi-search-sequence-motifs)s -c %(contigs-db)s \
                             -p %(profile-db)s
                             --motifs ATCG,TAAAT \
                             --store-in-db
{{ codestop }}

And this is how things will look like in the interface:

{{ codestart }}
%(anvi-interactive)s -c %(contigs-db)s \
                  -p %(profile-db)s
{{ codestop }}

[![motifs](../../images/layers_for_sequence_motifs.png){:.center-img .width-50}](../../images/layers_for_sequence_motifs.png)

Layers for sequence motif frequencies will be automatically colored to a shade of blue (although the user can change this through the %(interactive)s interface and/or through %(state)s files).

## Contigs database + genes database

Instead of a profile database, this program can also run on an anvi'o %(genes-db)s and search sequence motifs for each gene rather than split or contig sequences.