The tRNA-seq workflow is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow run by %(anvi-run-workflow)s.

The workflow can run the following programs in order:

- [Illumina-utils](https://github.com/merenlab/illumina-utils), for merging paired-end reads and quality control
- %(anvi-script-reformat-fasta)s, for making FASTA deflines anvio-compliant
- %(anvi-trnaseq)s, for predicting tRNA sequences, structures, and modification sites in each sample
- %(anvi-merge-trnaseq)s, for predicting tRNA seed sequences and their modification sites from the set of samples
- %(anvi-run-trna-taxonomy)s, for assigning taxonomy to tRNA seeds
- %(anvi-tabulate-trnaseq)s, for generating tables of seed and modification information that are easily manipulated

## Input

The tRNA-seq workflow requires two files to run: a %(workflow-config)s config file and a %(samples-txt)s. You can obtain a 'default' config file for this workflow to further edit using the following command.

{{ codestart }}
anvi-run-workflow -w trnaseq \
                  --get-default-config config.json
{{ codestop }}

Different "rules," or steps, of the workflow can be turned on and off as needed in the config file. The workflow can be restarted at intermediate rules without rerunning prior rules that have already completed.

%(samples-txt)s will contain a list of FASTQ or FASTA files and associated information on each library. FASTQ files contain unmerged paired-end tRNA-seq reads. Reads are merged in the workflow by [Illumina-utils](https://github.com/merenlab/illumina-utils). FASTA files contain merged reads, and the initial read-merging steps in the workflow are skipped.

Here is an example tRNA-seq samples file with FASTQ inputs.

| sample | treatment | r1 | r2 | r1_prefix | r2_prefix |
| --- | --- | --- | --- | --- | --- |
| ecoli_A1_noDM | untreated | FASTQ/ecoli_A1_noDM.r1.fq.gz | FASTQ/ecoli_A1_noDM.r2.fq.gz | NNNNNN | TTCCAGT |
| ecoli_A1_DM | demethylase | FASTQ/ecoli_A1_DM.r1.fq.gz | FASTQ/ecoli_A1_DM.r2.fq.gz | NNNNNN | TCTGAGT |
| ecoli_B1_noDM | untreated | FASTQ/ecoli_B1_noDM.r1.fq.gz | FASTQ/ecoli_B1_noDM.r2.fq.gz | NNNNNN | TGGTAGT |
| ecoli_B1_DM | demethylase | FASTQ/ecoli_B1_DM.r1.fq.gz | FASTQ/ecoli_B1_DM.r2.fq.gz | NNNNNN | CTGAAGT |

The treatment column is optional. The treatment indicates a chemical application, such as demethylase, and can be used to have a bearing on seed sequence determination in %(anvi-merge-trnaseq)s. In the absence of a treatment column, all samples are assigned the same treatment, which can be specified in the `anvi_trnaseq` section of the workflow config file and defaults to `untreated`.

Read 1 and 2 prefix columns are also optional. These represent sequences that Illumina-utils should identify and trim from the start of the read. In the example, the read 1 prefix is a unique molecular identifier (UMI) of 6 random nucleotides, and the read 2 prefix is a sample barcode. Illumina-utils will discard the paired-end read if the prefix is not found. In the example, the read 1 UMI will always be found, but the read 2 barcode must match exactly.

Here is an equivalent tRNA-seq samples file with FASTA inputs.

| sample | treatment | fasta |
| --- | --- | --- |
| ecoli_A1_noDM | untreated | FASTA/ecoli_A1_noDM.fa.gz |
| ecoli_A1_DM | demethylase | FASTA/ecoli_A1_DM.fa.gz |
| ecoli_B1_noDM | untreated | FASTA/ecoli_B1_noDM.fa.gz |
| ecoli_B1_DM | demethylase | FASTA/ecoli_B1_DM.fa.gz |

Note that barcodes and other sequence prefixes should already be trimmed from FASTA sequences.
