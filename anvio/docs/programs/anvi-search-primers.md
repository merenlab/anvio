This program finds all reads in a given set of FASTQ files provided as %(samples-txt)s based on user-provided primer sequences as %(primers-txt)s.

One of many potential uses of this program is to get back short reads that may be extending into hypervariable regions of genomes that often suffer from significant drops in coverage in conventional read-recruitment analyses, thus preventing any meaningful insights into coverage or variability patterns. In such situations, one can identify downstream conserved sequences (typically 15 to 25 nucleotides long) using the anvi'o interactive interface or through other means, and then provide those sequences to this program so it can find all matching sequences in a set of FASTQ files without any mapping.

{:.notice}
To instead get short reads mapping to a gene, use %(anvi-get-short-reads-mapping-to-a-gene)s.

Here is a typical command line to run it:

{{ codestart }}
anvi-search-primers --samples-txt %(samples-txt)s \
                    --primers-txt %(primers-txt)s \
                    --output-dir OUTPUT
{{ codestop }}

The %(samples-txt)s file is to list all the samples one is interested in, and the %(primers-txt)s file lists each primer sequence of interest, and their user-defined names. Each of these files can contain a single entry, or multiple ones.

This will output all of the matching sequences into three %(fasta)s files in the directory `OUTPUT`. These %(fasta)s files differ in their format and will include those that describe,

* Remainders are the downstream sequences after primer match, excluding the primer sequence.
* Primer matches are the primer-matching part of the match sequences (useful if one is working with degenerate primers and wishes to see the diversity of matching seqeunces).
* Trimmed sequences are trimmed to the shortest length (and include primer match). All matching sequences will start at the same position.
* Gapped sequences are not trimmed, but shorter ones are padded with gaps to eliminate length variation artificially.

The last two formats provide downstream possibilities to generate %(oligotypes)s and cluster short reads from an hypervariable region to estimate their diversity and oligotype proportion.

There will only be a single FASTA file in the output directory for raw sequences if the user asked only the primer matches to be reported with the flag `--only-report-primer-matches` or `--only-report-remainders`.

### For programmers

You can access to the functionality this program provides also programmatically. Here is an example:

``` python

import argparse

from anvio.sequencefeatures import PrimerSearch

# define a samples dictionary, there may be as many samples as you want
samples = {'sample_01': {'r1': 'sample_01_R1.fastq', 'r2': 'sample_01_R2.fastq'},
           'sample_02': {'r1': 'sample_02_R1.fastq', 'r2': 'sample_02_R2.fastq'}}

# define a primers dictionary, again, you may have as many primers as you
# wish
primers = {'primer_01': {'primer_sequence': 'GAGCAAAGATCATGTTTCAAAA.ACGTTC'},
           'primer_02': {'primer_sequence': 'AAGT.CTATCAGAACTTAGAGTAGAGCAC'},
           'primer_03': {'primer_sequence': 'GGCAGAAATGCCAAGT.CTATCAGAACTT'}}

# get an instance of the class, see the class header for all
# parameters.
s = PrimerSearch(argparse.Namespace(samples_dict=samples, primers_dict=primers, min_remainder_length=6))

# you can go through a for loop for each sample, or simply call
# s.process() to process all samples with all primers automatcially.
# here, though, this example will simply focus on a single sample
# to recover all primer hits, and then get sequences for a single primer
sample_dict, primers_dict = s.process_sample('sample_01')

# once primer hits are recovered, one can get any set of sequences
# of interest
sequences = s.get_sequences('primer_01', primers_dict, target='gapped')
print(sequences)
>>> ['GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGA-----------------------------------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGATGCAAC-----------------------------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGATGCAACA----------------------------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGATGCAACAA---------------------------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGATGCAACAAAGATAAGC-------------------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGATGCAACAAAGATAAGCCGCTTTTTT----------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCGTCTGATGCAACAAAGATAAGCCGCTTTTTT----------------------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCCTTTTTTGAAACACTGTTTTGGCTCTGCTCACTGAAGGCCAAAGG--------------------------------------------------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCCTTTTTTGAAACACTGTTTTGGCTCTGCTCACTGAAGGCCAAAGGAAGAGATAAATGGCTGATAATTAAAACAATGTAGAAATATTTGC------------------------',
     'GAGCAAAGATCATGTTTCAAAAGACGTTCCTTTTTTGAAACACTGTTTTGGCTCTGCTCACTGAAGGCCAAAGGAAGAGATAAATGGCTGATAATTAAAACAATGTAGAAATATTTGCACAGATGAAAAAAGCGGCTTATCT']

sequences = s.get_sequences('primer_01', primers_dict, target='trimmed')
print(sequences)
>>> ['GAAGATAGCCGTAGAAAGTGTAGAGTTTTAGGAGT',
     'AGCCGTAGAAAGTGTAGAGTTTCAGGAGTTTGGAG',
     'GCCGTAGAAAGTGTAGAGTTTTAGGAGTTTGGAGG',
     'CGTAGAAAGTGTAGAGTTTTAGGAGTTTGGAGGGG',
     'AGTGTAGAGTTTTAGGAGTTTGGAGGGGAGAATTA',
     'TTTAGGAGTTTGGAGGGGAGAATTAAGAAACGGTA',
     'TTTAGGAGTTTGGAGGGGAGAATTAAGAAACGGTA',
     'AGGGTAGAATTAAGAAACGGTAACGGTTGGTCTTG',
     'AAGAATAGTTGAAGAAGAATTATTGTATGGGAGAG',
     'TGTATGGGAGAGCAAAGATCATGTTTCAAAAGACG']

sequences = s.get_sequences('primer_01', primers_dict, target='primer_matches')
print(sequences)
>>> ['GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC',
     'GAGCAAAGATCATGTTTCAAAAGACGTTC']

s = PrimerSearch(argparse.Namespace(samples_dict=samples, primers_dict=primers, stop_after=10, min_remainder_length=6, only_keep_remainder=True))
sample_dict, primers_dict = s.process_sample('sample_01')
sequences = s.get_sequences('primer_01', primers_dict, target='remainder')
print(sequences)
>>> ['GTCTGA',
     'GTCTGA',
     'GTCTGA',
     'GTCTGA',
     'GTCTGA',
     'GTCTGA',
     'GTCTGA',
     'CTTTTT',
     'CTTTTT',
     'CTTTTT']
```
