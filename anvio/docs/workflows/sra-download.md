The `sra_download` workflow is a Snakemake workflow that downloads FASTQ files from SRA-accessions from [NCBI](https://www.ncbi.nlm.nih.gov/sra) e.g. SRR000001 and ERR000001. using [NCBI sra-tools wiki](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump), verifies the downlaod the MD5 checksum, gzips them using [pigz](https://zlib.net/pigz/), and provides a %(samples-txt)s. You will need to have these tools installed before you start.

Let's get started.

## Required input

### Configuration file

The first step is to make a %(workflow-config)s.

{{ codestart }}
anvi-run-workflow -w sra_download --get-default-config sra_download_config.json
{{ codestop }}

Here's what the %(workflow-config)s file looks like:

{{ codestart }}
$ cat sra_download_config.json
{
    "SRA_accession_list": "SRA_accession_list.txt",
    "Remove_unzipped_SRA_files": true,
    "prefetch": {
        "--max-size": "40g",
        "threads": 2
    },
    "fasterq_dump": {
        "threads": 6
    },
    "pigz": {
        "threads": 8
    },
    "output_dirs": {
        "SRA_prefetch": "01_NCBI_SRA",
        "FASTAS": "02_FASTA",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": "",
    "config_version": "3",
    "workflow_name": "sra_download"
}
{{ codestop }}

#### Modify any of the bells and whistles in the config file

{:.notice}
If this is the first time using an anvi'o Snakemake workflow, check out [Alon's blog post first](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#configjson).

Feel free to adjust anything in the config file! Here are some to consider:
- `threads`: this can be optimized for any of the steps depending on the size and number of SRA accessions you are downloaded.
- `prefetch` `--max-size`: The default is 40g but maybe you need more! For reference, this `--max-size` can download TARA Ocean metagenomes. You can use `vdb-dump --info` to learn how much the `prefetch` step will download e.g. `vdb-dump SRR000001 --info`. Read more about that [here](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump#check-the-maximum-size-limit-of-the-prefetch-tool).
- `Remove_unzipped_SRA_files`: Toggle this parameter to false if you want to keep the raw `.sra` files.

### List of SRA accessions

The input for the `sra_download` workflow is `SRA_accession_list.txt`. This contains a list of your SRA accessions you would like to download and it looks like this:

{:.warning}
All SRA accessions begin with the prefix `SRR` or `ERR` to denote their uploads to [NCBI](https://www.ncbi.nlm.nih.gov/sra) or [EBI](https://www.ebi.ac.uk/ena/browser/home) respectively.

{{ codestart }}
$ cat SRA_accession_list.txt
ERR6450080
ERR6450081
SRR5965623
{{ codestop }}

{:.warning}
The .sra files are stored in `01_NCBI_SRA/`. This directory will be deleted upon successful completion of the workflow because I don't know any use for .sra files. If you need these feel free to update the workflow.

## Start the workflow!

Here's a basic command to start the workflow:

### Run on your local computer

{{ codestart }}
anvi-run-workflow -w sra_download -c sra_download_config.json
{{ codestop }}

### Go big and use an HPC!

The power of Snakemake shines when you can leverage a High Performance Computing system to parallelize jobs. Check out the [Snakemake cluster documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#) on how to launch this workflow on your own HPC.

### MD5 Checksum verification

The `sra_download` workflow verifies that `prefetch` properly download all the `.sra` files you requested by comparing MD5 checksum values. 

If the workflow successfully runs then everything checked out! You can find the checksum messages by doing this:

{{ codestart }}
$ grep "Checksums match" 00_LOGS/*_check_md5sum.log
00_LOGS/ERR6450080_check_md5sum.log:Checksums match: a09e0eafd4c01f9685fb72efce447a2e
00_LOGS/ERR6450081_check_md5sum.log:Checksums match: 5c56ba814f649e77dac2091f1956e7aa
00_LOGS/ERR770960_check_md5sum.log:Checksums match: d1043d78938fb132f5d4b585f65ad8d7
00_LOGS/SRR5965623_check_md5sum.log:Checksums match: 95b69957f6381047763fb50782959cd8
{{ codestop }}

If a SRA accession fails to download properly, the workflow will stop and show you an error like this: 
{{ codestart }}
RuleException:
ValueError in file /Users/mschechter/github/anvio/anvio/workflows/sra_download/Snakefile, line 122:
[SRR5965623] MD5 checksum failed: expected 95b69957f6381047763fb50782959cd8, got 76057d4503f322d0cbd97897fce68872 for file 01_NCBI_SRA/SRR5965623/SRR5965623.sra
{{ codestop }}

## Common use cases

### Download sequencing files associated with an NCBI BioSample

Here is how to use the `sra_download` workflow to download all of the sequencing files from an NCBI BioSample:

1. Search for the [NCBI BioSample](https://www.ncbi.nlm.nih.gov/biosample/) under `All Databases` on the [NCBI website](https://www.ncbi.nlm.nih.gov/).
2. Under `Genomes` click `SRA`
3. Send results to Run selector by clicking `Send to:` and then `Run Selector`
4. Here you can filter for specific sequencing in the project OR you can download the `Metadata` or `Accession list` to download a text file with ALL of the SRA accesssions associated with the BioSample. Put the SRA accessions into the `SRA_accession_list.txt` and start the workflow!