This program exports the structures in a %(structure-db)s as %(protein-structure-txt)s files. 

This allows you to view your 3D structures through other softwares that aren't focused on variable positions (like %(anvi-display-structure)s is). 

To run, just provide a %(structure-db)s and an output path:

{{ codestart }}
anvi-export-structures -s %(structure-db)s \
                       -o path/to/output
{{ codestop }}

You can also provide a list of gene caller IDs, either directly through the parameter `--gene-caller-ids` or through a file with one gene caller ID per line. 
