This program finds all reads in a given set of FASTQ files provided as %(samples-txt)s based on user-provided primer sequences as %(primers-txt)s.

One of many potential uses of this program is to get back short reads that may be extending into hypervariable regions of genomes that often suffer from significant drops in coverage in conventional read-recruitment analyses, thus preventing any meaningful insights into coverage or variability patterns. In such situations, one can identify downstream conserved sequences (typically 15 to 25 nucleotides long) using the anvi'o interactive interface or through other means, and then provide those sequences to this program so it can find all matching sequences in a set of FASTQ files without any mapping.

{:.notice}
To instead get short reads mapping to a gene, use %(anvi-get-short-reads-mapping-to-a-gene)s.

Here is a typical command line to run it:

{{ codestart }}
anvi-script-get-primer-matches --samples-txt %(samples-txt)s \
                               --primers-txt %(primers-txt)s \
                               --output-dir OUTPUT
{{ codestop }}

The %(samples-txt)s file is to list all the samples one is interested in, and the %(primers-txt)s file lists each primer sequence of interest, and their user-defined names. Each of these files can contain a single entry, or multiple ones.

This will output all of the matching sequences into three %(fasta)s files in the directory `OUTPUT`. These %(fasta)s files differ in their format and will include those that describe,

* **Raw sequences**: Sequences from the FASTQ files that matched to a primer where each sequence reported as is with no processing.
* **Remainders**: Only downstream of the sequences after the primer match.
* **Trimmed sequences**: Raw sequences where the upstream of the primer sequence trimmed, as a result all matching sequences will start at the same position, and
* **Gapped sequences**: Trimmed sequences padded with gap characters to eliminate length variation artificially.

The last two formats provide downstream possibilities to generate %(oligotypes)s and cluster short reads from an hypervariable region to estimate their diversity and oligotype proportion.

There will only be a single FASTA file in the output directory for raw sequences if the user asked only the primer matches to be reported with the flag `--only-report-primer-matches` or `--only-report-remainders`.
