This program is used to add additional genes to or re-run the analysis of genes already within a %(structure-db)s.

For that reason, it is very similar to %(anvi-gen-structure-database)s and the parameters used to run that program (when you first generated your %(structure-db)s) will be automatically applied when you run this program. To know what MODELLER parameters are being used, you run this program on a %(structure-db)s with the flag `--list-modeller-params`. 

To run this program, just provide a %(contigs-db)s and %(structure-db)s, and name your genes of interest (either in a file or directly). If the named genes are not already in your %(structure-db)s, they will be added to the database. 

For example, if your %(structure-db)s already contains the genes with caller-IDs 1, 2 and 3, and you run

{{ codestart }}
anvi-update-structure-database -c %(contigs-db)s \
                               -s %(structure-db)s \
                               --gene-caller-ids 1,4,5
{{ codestop }}

Then the structural analysis for genes 4 and 5 will be added to your %(structure-db)s (assuming templates are found). Gene 1 will be ignored, since it is already present.

If instead you want to re-run the structural analysis on genes that are already in your %(structure-db)s, you'll need to specify that by adding the flag `--rerun-genes`

{{ codestart }}
anvi-update-structure-database -c %(contigs-db)s \
                               -s %(structure-db)s \
                               --gene-caller-ids 1,4,5 \
                               --rerun-genes
{{ codestop }}

Now, the program will rerun the analysis for gene 1 and will still add genes 4 and 5 to the %(structure-db)s. 

Both of these runs will have the same MODELLER parameters as your run of %(anvi-gen-structure-database)s. However, to get the raw outputs, you will need to use the parameter `--dump-dir`. You can also set a specific MODELLER program with `--modeller-executable`. Parameters for multi-threading would also have to be given again.

