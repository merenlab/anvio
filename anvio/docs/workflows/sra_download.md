The sra_download workflow is a Snakemake workflow that downloads FASTQ files from SRA-accessions using [NCBI sra-tools wiki](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) then gzips them using [pigz](https://zlib.net/pigz/). You will need to have these tools installed before you start.

{:.warning}
The workflow currently ONLY works for paired-end reads and WILL crash if an SRA accession points to any other kind of FASTQ file. Feel free to reach out if becomes important for this workflow to handle different kinds of FASTQ files is necessary. 

Let's get started.

## Required input

### Configuration file

The first step is to make a %(workflow-config)s.

```bash
anvi-run-workflow -w sra_download --get-default-config sra_download_config.json
```

Here's what the %(workflow-config)s file looks like:

```bash
$ cat sra_download_config.json
{
    "SRA_accession_list": "SRA_accession_list.txt",
    "prefetch": {
        "--max-size": "40g",
        "threads": 2
    },
    "fasterq_dump": {
        "threads": 6
    },
    "pigz": {
        "threads": 8,
        "--processes": ""
    },
    "output_dirs": {
        "SRA_prefetch": "01_NCBI_SRA",
        "FASTAS": "02_FASTA",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": "",
    "config_version": "3",
    "workflow_name": "sra_download"
```

#### Modify any of the bells and whistles in the config file

{:.notice}
If this is the first time using an anvi'o Snakemake workflow, I would check out [Alon's blog post first](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#configjson).

Feel free to adjust anything in the config file! Here are some to consider:
- `threads`: this can be optimized for any of the steps depending on the size and number of SRA accessions you are downloaded.
- `prefetch` `--max-size`: I already upped the amount from the default 40g but maybe you need more! For reference, I can download TARA Ocean metagenomes with the current parameter. You can use `vdb-dump --info` to learn how much the the `prefetch` step will download e.g. `vdb-dump SRR000001 --info`. Read more about that [here](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump#check-the-maximum-size-limit-of-the-prefetch-tool). 

### List of SRA accessions

The input for the `sra_download` workflow is `SRA_accession_list.txt`. This contains a list of your SRA accession you would like to download and it looks like this:

```bash
$ cat SRA_accession_list.txt
ERR6450080
ERR6450081
ERR6450082
```

{:.warning}
The .sra files are stored in `01_NCBI_SRA/`. This directory will be deleted upon successful completion of the workflow because I don't know any use for .sra files. If you need these feel free to update the workflow.

## Start the workflow!

Here's a basic command to start the workflow:

### Run on your local computer

```bash
anvi-run-workflow -w sra_download -c sra_download_config.json
```

### Go big and use an HPC!

The power of Snakemake shines when you can leverage a High Performance Computing system to parallize jobs. Check out the [Snakemake cluster documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#) on how to launch this workflow on your own HPC.
