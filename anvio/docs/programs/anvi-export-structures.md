This program exports structures from a %(structure-db)s into the globally recognized PDB format (%(protein-structure-txt)s), enabling their use in downstream analyses conducted outside of anvi'o.


To execute this program, provide a %(structure-db)s and an output path: 

{{ codestart }}
anvi-export-structures -s %(structure-db)s \
                       -o path/to/output
{{ codestop }}

You can also specify a list of gene caller IDs, either directly through the parameter `--gene-caller-ids` or through a file containing one gene caller ID per line.


