A `JSON`-formated configuration file that describes steps and parameters to be considered by an anvio workflow, which includes %(contigs-workflow)s, %(metagenomics-workflow)s, %(pangenomics-workflow)s, %(phylogenomics-workflow)s, and %(trnaseq-workflow)s.

You can create a default config file for a gein workflow using the following command:

```
anvi-run-workflow --workflow ANVIO-WORKFLOW \
                  --get-default-config CONFIG.json
```

For details of anvi'o snakemake workflows, please refer to [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/).
