Suppose you have downloaded some genomes from NCBI (using [this](https://github.com/kblin/ncbi-genome-download) incredibly useful program) and you have a metadata table describing those genomes. This program will convert that metadata table into some useful files, namely: a FASTA file of contig sequences, an external gene calls file, and an external functions file for each genome you have downloaded; as well as a single tab-delimited fasta-txt file (like the one shown [here](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)) describing the path to each of these files for all downloaded genomes (that you can pass directly to a snakemake workflow if you need to). Yay.

### The metadata file

The prerequisite for running this program is to have a tab-delimited metadata file containing information about each of the genomes you downloaded from NCBI. Let's say your download command started like this: `ncbi-genome-download --metadata-table ncbi_metadata.txt -t ....` So for the purposes of this usage tutorial, your metadata file is called `ncbi_metadata.txt`.

In case you are wondering, that file should have a header that looks something like this:
```
assembly_accession	bioproject	biosample	wgs_master	excluded_from_refseq	refseq_category	relation_to_type_material	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_dateasm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	local_filename
```

### Basic usage

If you run this, all the output files will show up in your current working directory.

{{ codestart }}
anvi-script-process-genbank-metadata -m ncbi_metadata.txt
{{ codestop }}

### Choosing an output directory

Alternatively, you can specify a directory in which to generate the output:

{{ codestart }}
anvi-script-process-genbank-metadata -m ncbi_metadata.txt -o DOWNLOADED_GENOMES
{{ codestop }}

### Picking a name for the fasta-txt file

The default name for the fasta-txt file is `fasta-input.txt`, but you can change that with the `--output-fasta-txt` parameter.

{{ codestart }}
anvi-script-process-genbank-metadata -m ncbi_metadata.txt --output-fasta-txt ncbi_fasta.txt
{{ codestop }}

### Make a fasta-txt without the gene calls and functions columns

The default columns in the fasta-txt file are:
```
name	path	external_gene_calls	gene_functional_annotation
```

But sometimes, you don't want your downstream snakemake workflow to use those external gene calls or functional annotations files. So to skip adding those columns into the fasta-txt file, you can use the `-E` flag:
{{ codestart }}
anvi-script-process-genbank-metadata -m ncbi_metadata.txt --output-fasta-txt ncbi_fasta.txt -E
{{ codestop }}

Then the fasta-txt will only contain a `name` column and a `path` column.
