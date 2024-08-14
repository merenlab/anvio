The aim of this script is to find potential assembly errors from long read assembler such as Flye, hifiasm, metaMDBG.

### Principle

This script search for potential assembly errors in a BAM file generated from the mapping of long read onto an assembly made using the same reads. The basic principle used by anvi'o is that (1) every single nucleotide, and (2) consecutive pair of nucleotide in an assembly should be covered by at least one long read that was used to generate the assembly in the first place.

To find out if every nucleotide are covered by at least one reads, one simply need to find region with 0x coverage, and that is one output of this script.

![zero_cov_example](../../images/anvi-script-find-misassemblies-01.png)

And to find where two consecutive nucleotides are not covered by at least one reads, we can leverage the clipping reported by long read mapping software like minimap2. Clipping happens when the left or right-most part of a read does not align to the reference. If a single nucleotide in the contigs is covered by 100%% of reads clipping at that position, it means that not a single read is covering at least two consecutive nucleotides. AND THAT'S NOT GOOD.

![clipping_example](../../images/anvi-script-find-misassemblies-02.png)


### Basic usage

The only input is a simple BAM file. But not any BAM file. It has to made by mapping long reads onto an assembly generated using the same reads. You also need to provide an output prefix.

{{ codestart }}
anvi-script-find-misassemblies -b sample01.bam -o result
{{ codestop }}


### Outputs

The first output is a table reporting regions in your assemblies with zero coverage. It includes the contig's length, start and stop position and length of the region with no coverage.

|**contig**|**length**|**range**|**range_size**|
|:--|:--|:--|:--|
|contig_001|1665603|0-498|498|
|contig_001|1665603|100500-101000|500|
|contig_001|1665603|1665106-1665603|497|

The second output is a table reporting the position with high relative abundance of clipped reads. It includes the contig's length, contig's position, the relative position in the contig (0 to 1), the total coverage at that position, the coverage of clipping at that position, and the ratio between these two coverage.

|**contig**|**length**|**pos**|**relative_pos**|**cov**|**clipping**|**clipping_ratio**|
|:--|:--|:--|:--|:--|:--|:--|
|contig_001|1665603|498|0.0002989908159387321|48|48|1.0|
|contig_001|1665603|500999|0.30079136504917436|120|120|1.0|
|contig_001|1665603|501000|0.30079196543233894|79|79|1.0|
|contig_001|1665603|1665105|0.9997010091840612|45|45|1.0|


### Additional parameters

By default, the script will report clipping position if the ratio of clipping reads to total number of read is over 80%%. You can change that threshold with the flag `--clipping-ratio`:

{{ codestart }}
anvi-script-find-misassemblies -b sample01.bam -o result --clipping-ratio 0.6
{{ codestop }}


Another default behaviour of the script is to skip the first and last 100bp of a contig (only valid for the clipping output). You can change that parameter with the flag `--min-dist-to-end`:

{{ codestart }}
anvi-script-find-misassemblies -b sample01.bam -o result --min-dist-to-end 500
{{ codestop }}
