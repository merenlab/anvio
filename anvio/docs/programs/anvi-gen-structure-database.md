This program attempts to solve for the 3D strucutres of proteins encoded by genes in your %(contigs-db)s using DIAMOND and MODELLER. 

MODELLER first searches your sequence(s) against a database of proteins with a known structure (in Anvi'o, this is either your %(pdb-db)s or the online copy of [the RCSB database](https://www.rcsb.org/) using [DIAMOND](http://www.diamondsearch.org/index.php). After sequence alignments, the program will select a base template based on the best hits. Then, the program creates a 3D alignment for your sequence and makes final adjustments to it based off of intermolecular interactions. For more information, see [here](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#how-modeller-works). 

The output of this is a %(structure-db)s, which can be used to run %(anvi-3dev)s to visualize all of this information. You can also export your strucutres into external .pdb files (%(anvi-export-structures)s), generate the fixation index matrix (%(anvi-gen-fixation-index-matrix)s), or the variability profile (%(anvi-gen-variability-profile)s). 

### Basic run 

Here is a simple run:

{{ codestart }}
anvi-gen-structure-database -c %(contigs-db)s \
                            --gene-caller-ids 1,2,3 \
                            -o STRUCTURE.db 
{{ codestop }}

Following this, you will have the strucutures for genes 1, 2, and 3 stored in `STRUCTURE.db`. Alternatively, you can provide a file name with the gene caller-ids (one ID per line) with the flag `--genes-of-interest`. 

To indicate that you have already run %(anvi-setup-pdb-database)s to set up a local copy of representative PDB structures, provide the path to your %(pdb-db)s:

{{ codestart }}
anvi-gen-structure-database -c %(contigs-db)s \
                            --gene-caller-ids 1,2,3 \
                            --pdb-database %(pdb-db)s \
                            -o STRUCTURE.db 
{{ codestop }}

To quickly get a very rough estimate for your structures, you can run with the flag `--very-fast`. 

### Advanced Parameters

Here, we will go through a brief overview of the MODELLER parameters that you are able to change. See [this page](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#description-of-all-modeller-parameters) for more information. 

- The number of models to be simultated. The default is 1. 
- The standard deviation of atomic perturbation of the initial strucutre (i.e. how much you change the position of the atoms before fine tuning with other analysis). The default is 4.
- The MODELLER database used. The default is `pdb_95`, which can be found [here](https://salilab.org/modeller/supplemental.html). This is the same database that is downloaded by %(anvi-setup-pdb-database)s. 
- The socring function used to compare potential models. The default is `DOPE_score`.
- The normalized percent identity cutoff for a template from the database to be further considered. 
- The maximum number of templates that the program will consider. The default is 5. 
- The MODELLER program to use. The default is `mod9.19`. 

For a case study on how some of these parameters matter, see [here](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#a-quick-case-study-on-the-importance-of-key-parameters). 

You also have the option to 

- Skip the use of DSSP, which predicts the locations of beta sheets, alpha helices, etc. 
- Turn on `offline mode`, which will prevent the system from looking up sequences it did not find a match to if it cannot find one in the downloaded database. 

Additionally, to output the raw data, just provide a path to the desired directory with the flag `--dump-dir`. 


