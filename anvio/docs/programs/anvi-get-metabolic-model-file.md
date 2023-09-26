This program **exports a metabolic %(reaction-network)s from a %(contigs-db)s to a %(reaction-network-json)s file** suitable for inspection and flux balance analysis.

The required input to this program is a %(contigs-db)s in which a %(reaction-network)s has been stored by %(anvi-reaction-network)s.

The %(reaction-network-json)s file output contains sections on the metabolites, reactions, and genes constituting the %(reaction-network)s that had been predicted from the genome. An "objective function" representing the biomass composition of metabolites in the ["core metabolism" of *E. coli*](http://bigg.ucsd.edu/models/e_coli_core) is automatically added as the first entry in the "reactions" section of the file and can be deleted as needed. An objective function is needed for flux balance analysis.

## Usage

%(anvi-get-metabolic-model-file)s requires a %(contigs-db)s as input and the path to an output %(reaction-network-json)s file.

{{ codestart }}
anvi-get-metabolic-model-file -c /path/to/contigs-db -o /path/to/ouput-json-file
{{ codestop }}

An existing file at the target output location must be explicitly overwritten with the `-W` flag.

{{ codestart }}
anvi-get-metabolic-model-file -c /path/to/contigs-db -o /path/to/existing-file -W
{{ codestop }}

The flag, `--remove-missing-objective-metabolites` must be used to remove metabolites in the *E. coli* core biomass objective function from the output file if the metabolites are not produced or consumed by the predicted %(reaction-network)s. [COBRApy](https://opencobra.github.io/cobrapy/), for instance, cannot load the JSON file if metabolites in the objective function are missing from the genomic model.

{{ codestart }}
anvi-get-metabolic-model-file -c /path/to/contigs-db -o /path/to/output-json-file --remove-missing-objective-metabolites
{{ codestop }}
