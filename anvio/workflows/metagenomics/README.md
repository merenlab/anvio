# Snakemake workflow for assembly based metagenomics

# contents

- [Introduction](#introduction)
- [Standard Usage](#standard-usage)
- [Reference Mode](#reference-mode)
- [Running the workflow on a cluster](#running-the-workflow-on-a-cluster)
	- [Defining the number of threads per rule in a cluster](#defining-the-number-of-threads-per-rule-in-a-cluster)
- [The config file](#the-config-file)
	- [General configurations](#general-configurations)
		- [Output Directories](#output-directories)
		- [The "all against all" option](#the-all-against-all-option)
		- [Optional Steps](#optional-steps)
	- [Step-Specific Configurations](#step-specific-configurations)
		- [QC](#qc)
		- [reformat_fasta](#reformat_fasta)
		- [megahit](#megahit)
		- [idba_ud](#idba_ud)
		- [run_centrifuge](#run_centrifuge)
		- [anvi_run_hmms](#anvi_run_hmms)
		- [anvi_run_ncbi_cogs](#anvi_run_ncbi_cogs)
		- [samtools_view](#samtools_view)
		- [bowtie](#bowtie)
		- [anvi_profile](#anvi_profile)
		- [anvi_merge](#anvi_merge)
	- [Example config.json file](#example-configjson-file)

# Introduction

**Important note**: this pipeline was evaluated using snakemake version 3.12.0. If you are using an older version, then we suggest upgrading to the newest version.

The majority of the steps used in this pipeline are extensively described in the [anvi'o user tutorial for metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). However, in contrast to that tuturial which starts with a FASTA files of contigs and BAM files, this pipeline includes steps to get there, including quality filtering, assembly, and mapping steps. detailed below.

The default entering point to this pipeline are the unprocessed raw reads from one or more metagenomes. The default output of the pipline is an anvi'o merged profile database ready for refinement of bins (or whatever it is that you want to do with it).

The pipline includes the following steps:

1. QC of the metagenomes using [illumina-utils](https://github.com/merenlab/illumina-utils/). Including generating a report of the results of the QC.
2. (Co-)Assembly of metagenomes using [megahit](https://github.com/voutcn/megahit).
3. Generating an anvi'o CONTIGS database using [anvi-gen-contigs-database](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-gen-contigs-database).
4. Running default anvi'o HMM profiles on the CONTIGS database using [anvi-run-hmms](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-hmms) (optional step).
5. Assigning taxonomy to genes with [centrifuge](https://ccb.jhu.edu/software/centrifuge/) and importing results into the CONTIGS database using [anvi-import-taxonomy](http://merenlab.org/2016/06/18/importing-taxonomy/) (optional step).
6. Mapping short reads from metagenomes to the contigs using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
7. Profiling individual BAM files using [anvi-profile](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile).
8. Merging resulting anvi'o profile databases using [anvi-merge](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-merge).

A directed acyclic graph (DAG) describing the workflow for a mock dataset could be seen below:

![Alt text](mock_files_for_merenlab_metagenomics_pipeline/mock-dag.png?raw=true "mock-dag")


If you want to create a DAG for your dataset, simply run:

```
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile --dag | dot -Tsvg > dag.svg
```

### Using `dot` on MAC OSX

If you are using MAC OSX, you can use `dot` by installing `graphviz`, simply run `brew install graphviz`.

# Standard usage

All you need is a bunch of FASTQ files, and a `samples.txt` file. A properly formatted `samples.txt` is available [here](mock_files_for_merenlab_metagenomics_pipeline/samples.txt).

The `samples.txt` file specifies the names of your samples and which group they belong to (if you optionally would like to do multiple co-assemblies as we did when we binned the TARA Oceans project metagenomes). It also describes where to find the pair-end FASTQ files (for now, we do not support single-end FASTQ runs).

The defalt name for your samples file is `samples.txt`, but you can use a different name by specifying it in the config file (see below).

[Back to Table of Contents](#contents)

# Reference Mode
## Estimating occurence of population genomes in metagenomes

Along with assembly-based metagenomics, we often use anvi'o to explore the occurence of population genomes accross metagenomes. You can see a nice example of that here: [Please insert a nice example here. Probably the blog about DWH thingy](link-to-nice-example).
In that case, what you have is a bunch of fastq files (metagenomes) and fasta files (reference genomes), and all you need to do is to let the workflow know where to find these files, using two `.txt` files: `samples.txt`, and `references.txt`. 

`references.txt` should be a 2 column tab-separated file, where the first column specifies a reference name and the second column specifies the filepath of the fasta file for that reference. An example `references.txt` can be found [here](mock_files_for_merenlab_metagenomics_pipeline/references.txt).


The `samples.txt` stays as before, but this time the `group` column will specify for each sample, which reference should be used (aka the name of the reference as defined in `references.txt`). If the `samples.txt` files doesn't have a `group` column, then an ["all against all"](#the-all-against-all-option) mode would be provoked. Below you can see how the DAG looks like for this mode:

![alt text](mock_files_for_merenlab_metagenomics_pipeline/mock-dag-references-mode.png?raw=true "mock-dag-references-mode")

After properly formatting your `samples.txt` and `references.txt`, reference mode is initiated by adding this to your `config.json`:

```
"references_txt": "references.txt"
```

[Back to Table of Contents](#contents)

# I only want to create a banch of contigs databases

Regardless if you are running in [reference mode](#reference-mode) or not, you can decide you want to only create contigs databases and annotate them with functions, taxonomy, etc. If you want to do that then simply run the following:

```bash
snakemake -s merenlab-metagenomics-pipeline.snakefile --until annotate_contigs_database
```

This would create the contigs databases (and also run assembly if that's what is needed), and would run the annotations that were specified according to your config file.

# Running the workflow on a cluster

When submitting to a cluster, you can utilize the [snakemake cluster execution](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution). Notice that the number of threads per rule could be changed using the `config.json` file (and not by using the [cluster configuration](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) file). For more details, refer to the documentation of the configuration file below.

When submitting a workflow to a cluster, snakemake requires you to limit the number of jobs using `--jobs`. If you prefer to limit the number of threads that would be used by your workflow (for example, if you share your cluster with others and you don't want to consume all resources), then you can make use of the snakemake built-in `resources` directive. You can set the number of jobs to your limit (or to a very big number if you dont care), and use `--resources nodes=30`, if you wish to only use 30 threads. We used the word `nodes` so that to not confuse with the reserved word `threads` in snakemake.

Notice: if you don't include `--jobs` (identical to `--cores`) in your command line, then snakemake will only use one node, and will not utilize multiple nodes even if the `threads` parameter for a rule is higher than 1. This is simply the behaviour of snakemake (described [here](http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads)).

## Defining the number of threads per rule in a cluster

In order to change the number of threads per rule when running on a cluster, the following structure should be used: 

```
	"rule_name":
		"threads": number_of_threads
```

The following defaults have been set:

**rule**|**threads**
:-----:|:-----:
qc|2
megahit|11
gen\_contigs\_db|5
run\_centrifuge|5
anvi\_run\_hmms|20
anvi\_run_\ncbi\_cogs|20
bowtie\_build|4
bowtie|10
samtools\_view|4
anvi\_init\_bam|4
anvi\_profile|5

All other rules use 1 thread by default.

Notice: if you want to use multiple threads, don't forget to include `--cores` in your snakemake command line. For more information, refer to [this section](http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads) of the snakemake documentation.

## A note on cluster-config

This note is here mainly for documentation of the code, and for those of you who are interested in snakemake. The reason we decided not to use the [cluster configuration](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) file to control the number of threads per rule, is becuase certain software require the number of threads as an input (for example `megahit` and `anvi-profile`), but the cluster config file is not available for shell commands within snakemake rules. To bypass this issue we simply put the threads configuration in the `config.json`, thus available for the user to modify.

[Back to Table of Contents](#contents)

# The config file

To make changes easy and accessible for the user, we tried our best to make all relevant configurations available to the user through a `JSON` formatted config file, and thus avoiding the need to change the Snakefile. An example config file is [here](mock_files_for_merenlab_metagenomics_pipeline/config.json). There are some general configurations, and there are step specific configurations.

## General configurations

### Output directories

The default directory structure that will appear in the output directory include these directories: "00\_LOGS", "01\_QC", "02\_ASSEMBLY", "03\_CONTIGS", "04\_MAPPING", "05\_ANVIO_PROFILE", "06\_MERGED"

Don't like these names? You can specify what is the name of the folder, by providing the following information in the config file:

```
    "output_dirs":{
        "LOGS_DIR" : "00_MY_beAuTiFul_LOGS",
        "QC_DIR" : "BEST_QC_DIR_EVER",
        "ASSEMBLY_DIR" : "assemblies",
        "CONTIGS_DIR" : "/absolute/path/to/my/contigs/dir",
        "MAPPING_DIR" : "relative/path/to/my/mapping/dir",
        "PROFILE_DIR": "/I/already/did/my/profiling/and/this/is/where/you/can/find/it/",
        "MERGE_DIR": "06_Keep_Calm_and_Merge_On"
    }
```

When using the "reference mode" (see below) the default name for the `ASSEMBLY_DIR` is `02_REFERENCE_FASTA`. You can change the default name the following way:

```
    "output_dirs":{
    	"REFERENCES_DIR" : "02_REF"
    }
```

You can change all, or just some of the names of these output folders. And you can provide an absolute or a relative path for them.

[Back to Table of Contents](#contents)

### The "all against all" option

The default behaviour for this workflow is to create a contigs database for each _group_ and map (and profile, and merge) the samples that belong to that _group_. If you wish to map all samples to all contigs, use the `all_against_all` option in the config file:

```
    "all_against_all: "True"
```

If you are new to `snakemake`, you might be surprised to see how easy it is to switch between modes. All we need to do is tell the `anvi_merge` rule that we want all samples merged for each _group_, and snakemake immediatly infers that it needs to also run the extra mapping, and profiling steps. *Thank you snakemake!* (says everyone).

An updated DAG for the workflow for our mock data is available below:

![alt text](mock_files_for_merenlab_metagenomics_pipeline/mock-dag-all-against-all.png?raw=true "mock-dag-all-against-all")

A little more of a mess! But also has a beauty to it :-).

[Back to Table of Contents](#contents)

### Optional steps

The following steps are only optional:

1. Assigning taxonomy with centrifuge (default is **not** running).
2. Running hmm profiles on the contigs database (default is **running**).
3. QC for the input metagenomes (default is **running**, but if you already performed QC, you can skip this step).
4. Reformating the labels of the fasta files with `anvi-script-reformat-fasta` (default is **running**).

For more details refer to the specific documentation for these steps below.

[Back to Table of Contents](#contents)

## Step-specific configurations 

Some of the steps in the workflow have parameters with defaults that could be changed. We tried to keep things flexible and accessible for the user, but we know we didn't do everything possible. If there is something that you want to have access to and is not possible, please create an issue on [github](https://github.com/merenlab/MerenLab-workflows/issues). Or, better yet, make those changes and send us a pull request. We plan to do a better job to let you access in a flexible form all the parameters of each step, and if this is of special interest to you, you can refer to the note below regarding wrappers.

The step-specific configurations in the `config.json` file always have the following structure:
```
	"step_name":{
		"configurable_parameter": "value"
	}
```

Notice that everything has to have quotation marks (to be compatible with the JSON format).

### qc

A report with the full results of the QC for each sample is generated. Below you can see an example:

```
sample	number of pairs analyzed	total pairs passed	total pairs passed (percent of all pairs)	total pair_1 trimmed	total pair_1 trimmed (percent of all passed pairs)	total pair_2 trimmed	total pair_2 trimmed (percent of all passed pairs)	total pairs failed	total pairs failed (percent of all pairs)	pairs failed due to pair_1	pairs failed due to pair_1 (percent of all failed pairs)	pairs failed due to pair_2	pairs failed due to pair_2 (percent of all failed pairs)	pairs failed due to both	pairs failed due to both (percent of all failed pairs)	FAILED_REASON_P	FAILED_REASON_P (percent of all failed pairs)	FAILED_REASON_N	FAILED_REASON_N (percent of all failed pairs)	FAILED_REASON_C33	FAILED_REASON_C33 (percent of all failed pairs)
01_QC/BM_HC_HMP_S001_01	1787927	1655580	92.60	580936	35.09	453098	27.37	132347	7.40	51935	39.24	49249	37.21	31163	23.55	0	0	30573	23.10	101774	76.90
01_QC/BM_HC_HMP_S002_01	4718043	4224421	89.54	828338	19.61	687063	16.26	493622	10.46	207050	41.95	197028	39.91	89544	18.14	0	0	277734	56.26	215888	43.74
01_QC/BM_HC_HMP_S003_01	1483467	1378881	92.95	688143	49.91	1378654	99.98	104586	7.05	22155	21.18	64076	61.27	18355	17.55	0	0	2436	2.33	102150	97.67
```

If you already performed QC, and so wish to skip qc, then simply add this to your config file:

```
	"qc": {
		"run": false
	}
```

A nice trick worth knowing: if you only want to qc your files and then compress them (and not do anything else), simply invoke the workflow with the following command:

```
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile --until gzip_fastqs
```

To understand this better, refer to the snakemake documentation.

[Back to Table of Contents](#contents)

### reformat_fasta

In "reference mode", you may choose to skip this step, and keep your contigs names. In order to do so, add this to your config file:

```json
	"reformat_fasta": {
		"run": false
	}
```

In assembly mode, this rule is always excecuted.

[Back to Table of Contents](#contents)

### megahit

`run` - You must specify this (`run: true`), otherwise you would probably get an error message. Notice that your config file should only include one assembly software, if it includes two, you would get an error message.

`memory` (see `-m/--memory` in the megahit documentation) - The default is 0.4.

`min_contig_len` (`--min-contig-len`) - default is 1,000.

[Back to Table of Contents](#contents)

### idba_ud

`run` - You must specify this (`run: true`), otherwise you would probably get an error message. Notice that your config file should only include one assembly software, if it includes two, you would get an error message.

`min_contig` - default is 1,000.

**Important**: if you are using `idba_ud` together, and you are skipping qc, then your fastq files must be uncompressed. This is because you will simply be providing a path to the pair of fastq files, and hence there would be no easy way for us to know that they are compresed or not, and in which format they are.

Another thing to note regarding `idba_ud` is that it requires a single fasta as an input. Because of that, what we do is use `fq2fa` to merge the pair of reads of each sample to one fasta, and then we use `cat` to concatenate multiple samples for a co-assembly. The `fasta` file that is created is create as a temporary file, and is deleted once `idba_ud` finishes running. If this is annoying to you, then feel free to contact us or just hack it yourself.

### run_centrifuge

`run` - could get values of `true` or `false` (all lower case!) - to configure whether to run centrifuge or not. The default is `false`.

`db` - if you choose run centrifuge, you **must** provide the path to the database (for example `$CENTRIFUGE_BASE/p+h+v/p+h+v`).

[Back to Table of Contents](#contents)

### anvi_run_hmms

`run` - could get values of `true` or `false` (all lower case!) - to configure whether to run hmms or not. The default is `true`.

[Back to Table of Contents](#contents)

### anvi_run_ncbi_cogs

`run` - could get values of `true` or `false` (all lower case!) - to configure whether to run hmms or not. The default is `true`.

Additionaly, you can set all the parameters that are available for `anvi-run-ncbi-cogs` (the default setting for all the following parameters is to take the default setting of `anvi-run-ncbi-cogs`):

`cog_data_dir` - path to the cog data directory.

`sensitive` - flag for DIAMOND sensitivity (should be either `true` or `false`. The default is `false`).

`temporary_dir_path` - see `anvi-run-ncbi-cogs` documentation.

`search_with` - see `anvi-run-ncbi-cogs` documentation.

Example:

```JSON
	"anvi_run_ncbi_cogs":
		"run": true,
		"cog_data_dir": "/USER/COG_DIR/",
		"sensitive": true,
		"temporary_dir_path": "/USER/MY_TEMP_DIR/",
		"search_with": "blastp",
		"threads": 1
```

[Back to Table of Contents](#contents)

### samtools_view

`s` - the samtools command executed is `samtools view {additional_params} -bS {stuff} -o {stuff}`, where `additional_params` specifies what goes in place of `{additional_params}` and `{stuff}` refers to stuff handled internally by our workflow (and therefore shouldn't be messed with). You can therefore specify all options that aren't `-bS` or `-o` with `additional_params`. For example, you could set `view_flag` to be `-f 2`, or `-f 2 -q 1` (for a full list see the samtools [documentation](http://www.htslib.org/doc/samtools.html)). The default is `-F 4`.

[Back to Table of Contents](#contents)

### bowtie

`additional_params` - the bowtie2 command executed is `bowtie2 --threads {stuff} -x {stuff} -1 {stuff} -2 {stuff} {additional_params} -S {stuff}`, where `additional_params` specifies what goes in place of `{additional_params}` and `{stuff}` refers to stuff handled internally by our workflow (and therefore shouldn't be messed with). You can therefore specify all parameters that aren't `--threads`, `-x`, `-1`, `-2`, or `-S` with `additional_params`. For example, if you don't want gapped alignment (aka the reference does not recruit any reads that contain indels with respect to it), and you don't want to store unmapped reads in the SAM output file, set `additional_params` to be `--rfg 10000,10000 --no-unal` (for a full list of options see the bowtie2 [documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#options)). The default is `--no-unal`.

[Back to Table of Contents](#contents)

### anvi_profile

`min_contig_length` - see anvi-profile docummentation for `--min-contig-length`. The default is going with the default of `anvi-profile` (which is 2,500).

[Back to Table of Contents](#contents)

### anvi_merge

`skip_concoct_binning` - see the `anvi-merge` docummentation for `--skip-concoct-binning`.

[Back to Table of Contents](#contents)

## Example config.json file

So let's say I want to run centrifuge, I don't want to run hmms, and I want my minimum contig length for megahit and anvi-profile to be 500 and and 3,000 respectively. Then my config file would like like this:

```
{
	"run_centrifuge":{
		"run": true,
		"db": "$CENTRIFUGE_BASE/p+h+v/p+h+v"
	},
	"anvi_run_hmms":{
		"run": false
	},
	"anvi_profile:{
		"min_contig_length": 3000,
		"threads": 10
	},
	"megahit":{
		"min_contig_len": 500
	}
}
```

[Back to Table of Contents](#contents)

## Wrappers

*(soon)*

[Back to Table of Contents](#contents)
