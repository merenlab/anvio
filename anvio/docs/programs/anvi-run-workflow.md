This program allows you to run [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows for common anvi'o processes. Some aspects of this program is described in [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/), and workflow-specific details can be found in worklfow pages.

Essentially, an anvi'o workflow will run several anvi'o programs for you in quick succession (based on a standard set of initial steps that will allow you to quickly get to a point where you can ask your data novel questions).

As of now, the available workflows are the %(contigs-workflow)s, the %(metagenomics-workflow)s, the %(pangenomics-workflow)s, the %(phylogenomics-workflow)s, the %(trnaseq-workflow)s, the %(ecophylo-workflow)s, and the %(sra-download-workflow)s.

### Before running the workflow

Each workflow requires a %(workflow-config)s: the file that details all of the parameters for the workflow. To get the %(workflow-config)s with the default parameters, just run

{{ codestart }}
anvi-run-workflow -w WORKFLOW-NAME \
                  --get-default-config CONFIG.json
{{ codestop }}

Before running a workflow, it is also a good idea to check the required dependencies by running

{{ codestart }}
anvi-run-workflow -w WORKFLOW-NAME \
                  --list-dependencies
{{ codestop }}

### The main run

The main run of the workflow should look like this:

{{ codestart }}
anvi-run-workflow -w WORKFLOW-NAME \
                  -c CONFIG.json
                  --save-workflow-graph
{{ codestop }}

The flag `--save-workflow-graph` creates a visual representation of the anvio programs that the workflow you're running used.

You can also use the `-A` flag at the end of the parameter list to change other [Snakemake](https://snakemake.readthedocs.io/en/stable/) parameters.
