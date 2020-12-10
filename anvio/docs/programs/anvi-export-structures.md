
This program exports the structures from a %(structure-db)s into the globally understood pdb format (%(protein-structure-txt)s), so they may be used for any follow-up analyses taking place outside of anvi'o.


To run, just provide a %(structure-db)s and an output path: 

{{ codestart }}
anvi-export-structures -s %(structure-db)s \
                       -o path/to/output
{{ codestop }}

You can also provide a list of gene caller IDs, either directly through the parameter `--gene-caller-ids` or through a file with one gene caller ID per line.


