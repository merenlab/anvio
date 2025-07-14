The aim of this script is to find potential assembly errors from long read assemblers. For a comprehensive assembly error analysis, which include the use of this script, you can check [the reproducible workflow](https://merenlab.org/data/benchmarking-long-read-assemblers/) for [Assemblies of long-read metagenomes suffer from diverse forms of errors”](https://doi.org/10.1101/2025.04.22.649783) by Trigodet et al.

### Principle

This script searches for potential assembly errors in a %(bam-file)s generated from the mapping of long reads onto an assembly made using the same reads. The basic principle is that (1) every single nucleotide, and (2) consecutive pairs/triplets of nucleotides in an assembly should be covered by at least one of the long reads that was used to generate the assembly.

To find out if every nucleotide is covered by at least one read, one simply needs to find regions with 0x coverage, and that is one output of this script.

To find where two/three consecutive nucleotides are not covered by at least one read, we can leverage the clipping information reported by long read mapping software like [minimap2](https://github.com/lh3/minimap2). Clipping happens when the left- or right-most part of a read does not align to the reference. If 100%% of reads clip at the same nucleotide position, then it means that not a single read is covering at least two consecutive nucleotides. AND THAT'S NOT GOOD.

Here is an example visualized in [IGV](https://igv.org/). All reads are clipped at the same position. There is no support in the long reads that the left and right sides of the contig should be joined, suggesting that there is a misassembly here.

![clipping_example](../../images/anvi-script-find-misassemblies.png)


### Basic usage

The only input file to this script is a simple %(bam-file)s. But not any BAM file. It has to made by mapping long reads onto an assembly generated using the same reads. You also need to provide an output prefix.

{{ codestart }}
anvi-script-find-misassemblies -b sample01.bam -o result
{{ codestop }}


### Outputs

The first output is a table reporting regions in your assemblies with zero coverage. This output includes the contig in which the region occurs, the contig's length, the region's start and stop positions, and the length of the region with no coverage.

|**contig**|**length**|**range**|**range_size**|
|:--|:--|:--|:--|
|contig_001|1665603|0-498|498|
|contig_001|1665603|100500-101000|500|
|contig_001|1665603|1665106-1665603|497|

The second output is a table reporting the positions with high relative abundance of clipped reads. It includes the contig in which the position occurs, the contig's length, the position in the contig with high levels of clipping, its relative position in the contig, the total coverage at that position, the coverage of reads clipping at that position, and the ratio between these two coverage values. The 'relative position' is the position divided by the contig's length (a value from 0 to 1), which tells you whether the position is close to the contig ends or somewhere in the middle.

|**contig**|**length**|**pos**|**relative_pos**|**cov**|**clipping**|**clipping_ratio**|
|:--|:--|:--|:--|:--|:--|:--|
|contig_001|1665603|498|0.0002989908159387321|48|48|1.0|
|contig_001|1665603|500999|0.30079136504917436|120|120|1.0|
|contig_001|1665603|501000|0.30079196543233894|79|79|1.0|
|contig_001|1665603|1665105|0.9997010091840612|45|45|1.0|


### Additional parameters

By default, the script will report clipping positions if the ratio of clipping reads to the total number of reads is over 80%%. You can change that threshold with the flag `--clipping-ratio`:

{{ codestart }}
anvi-script-find-misassemblies -b sample01.bam -o result --clipping-ratio 0.6
{{ codestop }}


Another default behavior of the script is to skip the first and last 100bp of a contig (only valid for the table reporting positions of high clipping). This is because contig ends often have high proportions of clipping which have nothing to do with misassemblies. You can change that parameter with the flag `--min-dist-to-end`:

{{ codestart }}
anvi-script-find-misassemblies -b sample01.bam -o result --min-dist-to-end 0
{{ codestop }}
