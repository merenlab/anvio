Get short reads from a %(bam-file)s in the form of %(short-reads-fasta)s).

{:.warning}
The purpose of this program is not to replace more efficient tools to recover short reads from BAM files such as `samtool`. Since it was designed to address much more subtle needs, this program may have a huge memory fingerprint for very large and numerous BAM files.

Using this program you can,

* Get all reads from one or more BAM files
* Get reads matching to contig names found in any %(bin)s in a %(collection)s
* Get reads matching to contig names found in one or more specific %(bin)ss in a %(collection)s
* Get all reads matching to a specific contig name
* Get reads matching to a specific region of a specific contig name

In addition, you can use the previously-defined fetch filters via the `--fetch-filter` parmeter to get only short reads satisfy a particular set of criteria (i.e., those that are in forward-forward or reverse-reverse orientation, those that have a template length longer than 1,000 nucleotides, and so on). For a complete set of fetch filters you can use, please see the help menu of the program.

The program can report all reads in a single file, or you can ask reads to be split into R1 and R2 files for mapping results of paired-end sequences using the flag `--split-R1-and-R2`. In this case, reads that are not paired will be reported in a file with the prefix `_UNPAIRED.fa`.

Reads reported as a FASTA will contain necessary information in their deflines to recover which BAM file, contig, sample they are from with explicit start/stop positions on the contig to which they matched.

### Getting all reads

A basic run of this program is as follows:

{{ codestart }}
anvi-get-short-reads-from-bam BAM_FILE_1.bam BAM_FILE_2.bam (...) \
                              --output-file OUTPUT.fa
{{ codestop }}

This will report all short reads found in BAM files `BAM_FILE_1.bam` and `BAM_FILE_2.bam` and store them into a single file. You can use as many BAM files as you wish.

### Narrowing the input with anvi'o files:

You can choose to only return the short reads that are contained within a %(collection)s:

{{ codestart }}
anvi-get-short-reads-from-bam BAM_FILE_1.bam BAM_FILE_2.bam \
                              -c %(contigs-db)s \
                              -p %(profile-db)s \
                              -C %(collection)s \
                              --output-file OUTPUT.fa
{{ codestop }}

Or in a bin that is described in a collection:

{{ codestart }}
anvi-get-short-reads-from-bam BAM_FILE_1.bam BAM_FILE_2.bam \
                              -c %(contigs-db)s \
                              -p %(profile-db)s \
                              -C %(collection)s \
                              -b %(bin)s \
                              --output-file OUTPUT.fa
{{ codestop }}

### Focusing on individual contigs

You can get all reads mapped to a contig:

{{ codestart }}
anvi-get-short-reads-from-bam BAM_FILE_1.bam BAM_FILE_2.bam \
                              --target-contig CONTIG_NAME \
                              --output-file OUTPUT.fa
{{ codestop }}

Or define explicit start/stop positions on it:

{{ codestart }}
anvi-get-short-reads-from-bam BAM_FILE_1.bam BAM_FILE_2.bam \
                              --target-contig CONTIG_NAME \
                              --target-region-start 100 \
                              --target-region-end 1000 \
                              --output-file OUTPUT.fa
{{ codestop }}

In this mode, the program will fetch any read that includes a nucleotide that matches to anywhere in the region defined by the user. Which means, if the user sets `--target-region-start` to `100` and `--target-region-end` to `101`, all reads that have a nuclotide mapping to the `100th` position will be returned.

### Changing the output format

You can split the output based on the directionality of paired-end reads. Adding the tag `--split-R1-and-R2` causes the program to create three separate output files: one for R1 (sequences in the forward direction), one for R2 (sequences in the reverse direction; i.e. reverse complement of R1 sequences), and one for unparied reads. When doing this, you can name these three files with a prefix by using the flag `-O`.

{{ codestart }}
anvi-get-short-reads-from-bam -o path/to/output \
                              --split-R1-and-R2 \
                              -O BAM_1_and_BAM_2 \
                              BAM_FILE_1.bam BAM_FILE_2.bam
{{ codestop }}

You can also compress the output by adding the flag `--gzip-output`.
