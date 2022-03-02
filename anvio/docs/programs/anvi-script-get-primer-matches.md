This program finds all reads in a given set of FASTQ files based on user-provided primer sequences.

The primary utility of this program is to get back short reads that may be extending into hypervariable regions of genomes that often suffer from significant drops in coverage in conventional read-recruitment analyses, thus preventing any meaningful insights into coverage or variability patterns.

In these situations, one can identify downstream conserved sequences (typically 15 to 25 nucleotides long) using the anvi'o interactive interface or through other means, and then provide those sequences to this program so it can find all matching sequences in a set of FASTQ files without any mapping.

{:.notice}
To instead get short reads mapping to a gene, use %(anvi-get-short-reads-mapping-to-a-gene)s.

Here is a typical command line to run it:

{{ codestart }}
anvi-script-get-primer-matches --%(samples-txt)s samples.txt  \
                               --primer-sequences sequences.txt \
                               --output-dir OUTPUT
{{ codestop }}

The %(samples-txt)s file is to list all the samples one is interested in, and the primer sequences file lists each primer sequence of interest. Each of these files can contain a single entry, or multiple ones.

This will output all of the matching sequences into three %(fasta)s files in the directory `OUTPUT`. These %(fasta)s files differ in their format and will include those that describe,

* Raw sequences: sequences from the FASTQ files that matched to a primer where each sequence reported as is with no processing.
* Trimmed sequences: Raw sequences where the upstream of the primer sequence trimmed, as a result all matching sequences will start at the same position, and
* Gapped sequences: Trimmed sequences padded with gap characters to eliminate length variation artificially.

The last two formats provide downstream possibilities to generate %(oligotypes)s and cluster short reads from an hypervariable region to estimate their diversity and oligotype proportion. 
