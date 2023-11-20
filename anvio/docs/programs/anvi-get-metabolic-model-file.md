This program **exports a metabolic %(reaction-network)s from a %(contigs-db)s OR a %(pan-db)s and %(genomes-storage-db)s to a %(reaction-network-json)s file** formatted for flux balance analysis.

The required input to this program is a %(contigs-db)s OR a %(pan-db)s in which a %(reaction-network)s has been stored by %(anvi-reaction-network)s. The %(pan-db)s must be accompanied by a %(genomes-storage-db)s input.

The %(reaction-network-json)s file output contains sections on the metabolites, reactions, and genes (or gene clusters) constituting the %(reaction-network)s that had been predicted from the genome (or pangenome). An "objective function" representing the biomass composition of metabolites in the ["core metabolism" of *E. coli*](http://bigg.ucsd.edu/models/e_coli_core) is automatically added as the first entry in the "reactions" section of the file and can be deleted as needed. An objective function is needed for flux balance analysis.

## Usage

%(anvi-get-metabolic-model-file)s requires a %(contigs-db)s OR a %(pan-db)s and %(genomes-storage-db)s as input, plus the path to an output %(reaction-network-json)s file.

{{ codestart }}
anvi-get-metabolic-model-file -c /path/to/contigs-db \
                              -o /path/to/ouput.json
{{ codestop }}

{{ codestart }}
anvi-get-metabolic-model-file -p /path/to/pan-db \
                              -g /path/to/genomes-storage-db \
                              -o /path/to/output.json
{{ codestop }}

An existing file at the target output location must be explicitly overwritten with the flag, `--overwrite-output-destinations`.

{{ codestart }}
anvi-get-metabolic-model-file -c /path/to/contigs-db \
                              -o /path/to/output.json \
                              --overwrite-output-destinations
{{ codestop }}

The flag, `--remove-missing-objective-metabolites` must be used to remove metabolites in the *E. coli* core biomass objective function from the %(reaction-network-json)s file if the metabolites are not produced or consumed by the predicted %(reaction-network)s. [COBRApy](https://opencobra.github.io/cobrapy/), for instance, cannot load the JSON file if metabolites in the objective function are missing from the model.

{{ codestart }}
anvi-get-metabolic-model-file -c /path/to/contigs-db \
                              -o /path/to/output.json \
                              --remove-missing-objective-metabolites
{{ codestop }}

It is possible that the gene KO annotations used to construct the stored reaction network have since been changed in the %(contigs-db)s or the %(genomes-storage-db)s. By default, without using the flag, `--ignore-changed-gene-annotations`, this program checks that the set of gene KO annotations that is currently stored was also that used in construction of the %(reaction-network)s, and raises an error if this is not the case. Use of this flag ignores that check, permitting the set of gene annotations to have changed since creation of the network.

{{ codestart }}
anvi-get-metabolic-model-file -p /path/to/contigs-db \
                              -o /path/to/output.json \
                              --ignore-changed-gene-annotations
{{ codestop }}

For a pangenomic network, the option `--record-genomes` determines which additional information is added to the output %(reaction-network-json)s file regarding genome membership. By default, genome names are recorded for gene clusters and reactions, which is equivalent to `--record-genomes cluster reaction`. 'cluster' records in the 'notes' section of each 'gene' (cluster) entry in the JSON file which genomes are part of the cluster. 'reaction' and 'metabolite', respectively, record the genomes predicted to encode enzymes associated with reaction and metabolite entries. The arguments, 'cluster', 'reaction', and 'metabolite', are valid, and are all used in the following example.

{{ codestart }}
anvi-get-metabolic-model-file -p /path/to/pan-db \
                              -g /path/to/genomes-storage-db \
                              --record-genomes cluster reaction metabolite
{{ codestop }}

The use of `--record-genomes` as a flag without any arguments prevents genome membership from being recorded at all in the %(reaction-network-json)s file.

{{ codestart }}
anvi-get-metabolic-model-file -p /path/to/pan-db \
                              -g /path/to/genomes-storage-db \
                              --record-genomes
{{ codestop }}
