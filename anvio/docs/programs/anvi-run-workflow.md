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
