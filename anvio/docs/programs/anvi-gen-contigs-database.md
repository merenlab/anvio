The input for this program is a %(contigs-fasta)s, which should contain one or more sequences. These sequences may belong to a single genome or could be many contigs obtained from an assembly or a single sequence of any kind.

Make sure the input file matches the requirements of a %(contigs-fasta)s. If you are planning to use the resulting contigs-db with %(anvi-profile)s, it is essential that you convert your %(fasta)s file to a properly formatted %(contigs-fasta)s *before* you perform the read recruitment.

The contigs database is one of the most essential components of anvi'o, and a contigs database will keep all the information related to your sequences: positions of open reading frames, k-mer frequencies for each contig, functional and taxonomic annotation of genes, etc. 

When run on a %(contigs-fasta)s this program will,

* **Compute k-mer frequencies** for each contig (the default is `4`, but you can change it using `--kmer-size` parameter if you feel adventurous).

* **Soft-split contigs** longer than 20,000 bp into smaller ones (you can change the split size using the `--split-length` flag). When the gene calling step is not skipped, the process of splitting contigs will consider where genes are and avoid cutting genes in the middle. For very, very large assemblies this process can take a while, and you can skip it with `--skip-mindful-splitting` flag.

* **Identify open reading frames** using [Prodigal](http://prodigal.ornl.gov/), UNLESS, (1) you have used the flag `--skip-gene-calling` (no gene calls will be made) or (2) you have provided %(external-gene-calls)s.

{:.notice}
This program can work with compressed input FASTA files (i.e., the file name ends with a `.gz` extention).

### Create a contigs database from a FASTA file

The simplest form of this command that will give you a %(contigs-db)s is the following:

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s
{{ codestop }}

But we suggest you to explicitly assign a unique 'project name' for each contigs-db you are generating through the `--project-name` parameter. Project name becomes an idenfitier of a contigs-db for most downstream analyses, and when you are working with many contigs-db files, non-unique names for each one of them may lead to various issues. Here is an example for a single genome:

{{ codestart }}
anvi-gen-contigs-database -f Patient_6557_E_faecalis_cultivar.fa \
                          --project-name E_faecalis_P6557 \
                          -o E_faecalis_P6557.db
{{ codestop }}

and a metagenome:

{{ codestart }}
anvi-gen-contigs-database -f scaffolds.fa \
                          --project-name North_Atlantic_MGX_004 \
                          -o North_Atlantic_MGX_004.db
{{ codestop }}

There are a myriad of programs you can run on a %(contigs-db)s once it is created to add more and more layers of information on it. Please see the artifact %(contigs-db)s to see a list of steps you can follow.

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

### Changing k-mer size

You can change the k-mer size by modifying the `--kmer-size` parameter:

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s \
                          --kmer-size 3
{{ codestop }}

A word of caution: you can increase the k-mer size up to a maximum of k=5 for standard installations of anvi'o. This is because the contigs database stores the k-mer frequencies in a big table with one column per k-mer, and SQLite has an upper limit on the number of columns per table. [The default limit is 2,000 columns](https://www.sqlite.org/limits.html), which translates into an upper limit of k=5 (with 4^5 = 1,096 possible k-mers). Trying to increase `k` beyond this point will result in the following error: `sqlite3.OperationalError: too many columns on kmer_contigs`.

If you want to increase `k` even further, you can re-compile the `sqlite3` library to increase the column limit (the constant `SQLITE_MAX_COLUMN`). Note that you can only increase this limit up to a maximum of 32,000 columns, which makes k=7 (with 4^7 = 16,384 possible k-mers) the new upper limit for k-mer size.