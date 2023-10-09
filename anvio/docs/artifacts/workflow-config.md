A `JSON`-formated configuration file that describes steps and parameters to be considered by an anvio %(workflow)s.

You can create a default config file for a given workflow using the following command:

```
anvi-run-workflow --workflow ANVIO-WORKFLOW \
                  --get-default-config CONFIG.json
```

Following this, the file `CONFIG.json` will contain all configurable flags and parameters set to their default value for that workflow. From there, you can edit this file to your hearts content.

### What's in this file?

The config file contains three types of information:

1. **General parameters**, including the name of the workflow, the version of this config file, and links to the %(fasta-txt)s or %(samples-txt)s file)
2. **Rule specific parameters** which allow you to set the parameters on individual anvi'o programs that are run in the workflow.
3. **Output directory names** which just tell anvi'o what to name all of the intermediate and final outputs (to help keep things organized).

For example, the default config file for the [contigs workflow](../../workflows/contigs) has no rule specific parameters and looks like this:

    {
        "workflow_name": "contigs",
        "config_version": 1,
        "fasta_txt": "fasta.txt",
        "output_dirs": {
            "FASTA_DIR":   "01_FASTA_contigs_workflow",
            "CONTIGS_DIR": "02_CONTIGS_contigs_workflow",
            "LOGS_DIR":    "00_LOGS_contigs_workflow"
        }
    }

On the other hand, the default config file for the [contigs workflow](../../workflows/metagenomics) is much longer, because it has sections for each rule specific parameter. For example, its section on parameters for the program %(anvi-gen-contigs-database)s looks like this:

    "anvi_gen_contigs_database": {
       "--project-name": "{group}",
       "threads": 5,
       "--description": "",
       "--skip-gene-calling": "",
       "--ignore-internal-stop-codons": "",
       "--skip-mindful-splitting": "",
       "--contigs-fasta": "",
       "--split-length": "",
       "--kmer-size": ""
    },

Note that the empty string `""` here means that the default parameter for the program %(anvi-gen-contigs-database)s will be used.

For more details on the anvi'o snakemake workflows, please refer to [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/).

