This script takes a FASTQ file and a short input sequence and finds all of the short reads in your fastq file that align to your short sequence. 

The purpose of this is to get back short reads that may be extending into hypervariable regions of genomes, resulting a decreased mappability of short reads in the metagenome given a reference. You often see those areas of genomes as significant dips in coverage, and in most cases with a large number of SNVs. When you provide the downstream conserved sequence, this program allows you to take a better look at those regions at the short read level without any mapping.

To instead get short reads mapping to a gene, use %(anvi-get-short-reads-mapping-to-a-gene)s.

Here is an example run of this program with the default parameters, where the user is searching for alignments to `AAAAAAAAAAAA` in the sample named `example_sample` stored in the two fastq files `fastaq_one.fastq` and `fastq_two.fastq`: 

{{ codestart }}
anvi-script-get-short-reads-matching-something --match-sequence AAAAAAAAAAAA \
                                               -s example_sample \ 
                                               -O example_sample_AAAAAAA_results
                                               fastaq_one.fastq fastq_two.fastq
{{ codestop }}

This will output all of the matching sequences into a %(fasta)s file in the directory `example_sample_AAAAAAA_results`. 

Note that this will only report sequences where the length of the short read after the matching sequence is above a certain threshold. The default is 60. For example, if this dataset has the sequence `TTAAAAAAAAAAAAGGGGGGGGG`, this would not be included in the results, but if the sequence was followed by 60 `G` nucleotides, it would be because the length of the sequence after the match is longer than the threshold. You can change this threshold with the parameter `--min-remainder-length`.

You can also choose to stop the program after it finds a certain number of matches or report the raw sequences instead of trimming them to the relevant sections. 
