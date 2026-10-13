This program allows you to run a %(workflow)s implemented by anvi'o developers for various commonly used set of steps to typically process your raw data (i.e., short reads or contigs from genomes, transcriptomes, metagenomes, metatranscriptomes, etc). Some aspects of this program is described in [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/).

For a list of currently available anvi'o workflows, please see the %(workflow)s artifact.

### Before running the workflow

Each workflow requires a %(workflow-config)s: the file that details all of the parameters for the workflow. To get the %(workflow-config)s with the default parameters, just run

{{ codestart }}
anvi-run-workflow -w %(workflow)s \
                  --get-default-config CONFIG.json
{{ codestop }}

Before running a workflow, it is also a good idea to check the required dependencies by running

{{ codestart }}
anvi-run-workflow -w %(workflow)s \
                  --list-dependencies
{{ codestop }}

### The main run

The main run of the workflow should look like this:

{{ codestart }}
anvi-run-workflow -w %(workflow)s \
                  -c CONFIG.json
                  --save-workflow-graph
{{ codestop }}

The flag `--save-workflow-graph` creates a visual representation of the anvio programs that the workflow you're running used.

You can also use the `-A` flag at the end of the parameter list to change other [Snakemake](https://snakemake.readthedocs.io/en/stable/) parameters.

### Logs and the workflow manifest

Every workflow run writes rule logs under `00_LOGS`, within a subdirectory for the workflow or named run. Logs are then organized by rule name:

{{ codestart }}
00_LOGS/<workflow-or-run-name>/<rule-name>/<job-specific-name>.log
{{ codestop }}

For instance, a metagenomics profiling job may write to `00_LOGS/metagenomics/anvi_profile/G01-S01.log`, while an SRA download checksum job may write to `00_LOGS/sra_download/check_md5sum/SRR5965623.log`.

Each workflow run also creates a tab-delimited manifest in the same workflow-specific log directory:

{{ codestart }}
00_LOGS/<workflow-or-run-name>/<workflow-name>-workflow-manifest.tsv
{{ codestop }}

This file lists the status of each Snakemake job, the rule name, the `group` and `read` wildcards when they exist, the rule log path, and the Snakemake log path when Snakemake reports one. If a workflow stops because a rule failed, this manifest is often the quickest way to find the relevant rule log.
