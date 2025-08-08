
This workflow is extremely useful if you have one or more %(fasta)s files that describe one or more contig sequences for your genomes or assembled metagenomes, and all you want to turn them into %(contigs-db)s files.

{:.warning}
If you have not yet run anvi'o programs %(anvi-setup-ncbi-cogs)s and %(anvi-setup-scg-taxonomy)s on your system yet, you will get a cryptic error from this workflow if you run it with the default %(workflow-config)s. You can avoid this by first running these two anvi'o programs to setup the necessary databases (which is done only once for every anvi'o installation), **or** set the rules for COG functions and/or SCG taxonomy to `run=false` explicitly.

To start things going with this workflow, first ask anvi'o to give you a default %(workflow-config)s file for the contigs workflow:

```bash
anvi-run-workflow -w contigs \
                  --get-default-config config-contigs-default.json
```

This will generate a file in your work directory called `config-contigs-default.json`. You should investigate its contents, and familiarize youself with it. It should look something like this, but much longer:
and you could examine its content to find out all possible options to tweak. We included a much simpler config file, `config-contigs.json`, in the mock data package for the sake of demonstrating how the contigs workflow works:

```json
{
    "workflow_name": "contigs",
    "config_version": "2",
    "fasta_txt": "fasta.txt",
    "output_dirs": {
        "FASTA_DIR": "01_FASTA",
        "CONTIGS_DIR": "02_CONTIGS",
        "LOGS_DIR": "00_LOGS"
    }
}
```

The only mandatory thing you need to do is to (1) manually create a %(fasta-txt)s file to describe the name and location of each FASTA file you wish to work with, and (2) make sure the `fasta_txt` variable in your %(workflow-config)s point to the location of your %(fasta-txt)s.

To see if everything looks alright, you can simply run the following command, which should generate a 'workflow graph' for you, given your config file parameters and input files:

```bash
anvi-run-workflow -w contigs \
                  -c config-contigs.json \
                  --save-workflow-graph
```

For the example config file shown above, this command will generate something similar to this:

[![DAG-contigs](../../images/workflows/contigs/DAG-contigs.png)]( ../../images/workflows/contigs/DAG-contigs.png){:.center-img .width-50}

{:.notice}
Please note that the generation of this workflow graph requires the usage of a program called [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)). If you are using MAC OSX, you can use [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) by installing [graphviz](http://www.graphviz.org/) through `brew` or `conda`.

If everything looks alright, you can run this workflow the following way:

```bash
anvi-run-workflow -w contigs \
                  -c config-contigs.json
```

If everything goes smoothly, you should see happy messages flowing on your screen, and at the end of it all you should see your contigs databases are generated and annotated properly. At the end of this process, you will have all your %(contigs-db)s files in the `02_CONTIGS` directory (as per the instructions in the config file, which you can change). You can use the program %(anvi-display-contigs-stats)s on one of them to see if everything makes sense.
