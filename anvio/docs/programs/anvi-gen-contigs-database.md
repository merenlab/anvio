The input for this program is a %(contigs-fasta)s, which should contain one or more sequences. These sequences may belong to a single genome or could be contigs obtained from an assembly.

Make sure the input file matches the requirements of a %(contigs-fasta)s. If you are planning to use the resulting contigs-db with %(anvi-profile)s, it is essential that you convert your %(fasta)s file to a properly formatted %(contigs-fasta)s *before* you perform the read recruitment.

An anvi'o contigs database will keep all the information related to your sequences: positions of open reading frames, k-mer frequencies for each contig, functional and taxonomic annotation of genes, etc. The contigs database is one of the most essential components of anvi'o.

When run on a %(contigs-fasta)s this program will,

* **Compute k-mer frequencies** for each contig (the default is `4`, but you can change it using `--kmer-size` parameter if you feel adventurous).

* **Soft-split contigs** longer than 20,000 bp into smaller ones (you can change the split size using the `--split-length` flag). When the gene calling step is not skipped, the process of splitting contigs will consider where genes are and avoid cutting genes in the middle. For very, very large assemblies this process can take a while, and you can skip it with `--skip-mindful-splitting` flag.

* **Identify open reading frames** using [Prodigal](http://prodigal.ornl.gov/), UNLESS, (1) you have used the flag `--skip-gene-calling` (no gene calls will be made) or (2) you have provided %(external-gene-calls)s.

{:.notice}
This program can work with compressed input FASTA files (i.e., the file name ends with a `.gz` extention).

### Create a contigs database from a FASTA file

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s
{{ codestop }}

### Create a contigs database with external gene calls

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s \
                          --external-gene-calls %(external-gene-calls)s
{{ codestop }}

See %(external-gene-calls)s for the description and formatting requirements of this file.

If user-provided or anvi'o-calculated amino acid sequences contain internal stop codons, anvi'o will yield an error. The following command will persist through this error:

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s \
                          --external-gene-calls %(external-gene-calls)s \
                          --ignore-internal-stop-codons
{{ codestop }}
