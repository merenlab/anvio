This program finds all short reads from (%(bam-file)s) that align to a specific gene and returns them as a %(short-reads-fasta)s.

If instead you want to extract these short reads from a FASTQ file, get your gene sequence with %(anvi-export-gene-calls)s and take a look at %(anvi-search-primers)s.

To run this program, just specify the bam files you're looking at and the gene of interest. To do this, name the %(contigs-db)s containing your gene and the gene caller ID (either directly through the parameter `--gene-caller-id` or through a file). Here is an example:

{{ codestart }}
anvi-get-short-reads-mapping-to-a-gene -c %(contigs-db)s \
                                       --gene-caller-id 2 \
                                       -i BAM_FILE_ONE.bam \
                                       -O GENE_2_MATCHES
{{ codestop }}

The output of this will be a file named `GENE_2_MATCHES_BAM_FILE_ONE.fasta` (prefix + bam file name), which will contain all short reads that aligned to gene 2 with more than 100 nucleotides.

You also have the option to provide multiple bam files; in this case, there will be an output files for each bam file inputted.

Additionally, you can change the number of nucleotides required to map to a short read for it to be reported. For example, to expand your search, you could decrease the required mapping length to 50 nucleotides, as so:

{{ codestart }}
anvi-get-short-reads-mapping-to-a-gene -c %(contigs-db)s \
                                       --gene-caller-id 2 \
                                       -i Bam_file_one.bam Bam_file_two.bam \
                                       -O GENE_2_MATCHES \
                                       --leeway 50
{{ codestop }}
