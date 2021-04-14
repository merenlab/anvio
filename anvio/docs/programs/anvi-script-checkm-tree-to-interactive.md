A helper script to process CheckM tree output to generate files compatible with %(anvi-interactive)s.

An example use:

{{ codestart }}
anvi-script-checkm-tree-to-interactive -t CheckM_concatenated.tree \
                                       -o OUTPUT_PATH
cd OUTPUT_PATH/
anvi-interactive -p PROFILE.db \
                 -t newick.tree \
                 -d view_data.txt \
                 --manual
{{ codestop }}
