
This program creates a %(structure-db)s either by (a) attempting to solve for the 3D structures of proteins encoded by genes in your %(contigs-db)s using DIAMOND and MODELLER, or (b) importing pre-existing structures provided by the user using an %(external-structures)s file.

### The basics of the pipeline

This section covers option (a), where the user is interested in having structures predicted for them.

DIAMOND first searches your sequence(s) against a database of proteins with a known structure.  This database is downloaded from the [Sali lab](https://salilab.org/modeller/supplemental.html), who created and maintain MODELLER, and contains all of the PDB sequences clustered at 95%% identity.

If any good hits are found, they are selected as templates, and their structures are nabbed either from [the RCSB directly](https://www.rcsb.org/), or from a local %(pdb-db)s database which you can create yourself with %(anvi-setup-pdb-database)s. Then, anvi'o passes control over to MODELLER, which creates a 3D alignment for your sequence to the template structures, and makes final adjustments to it based off of empirical distributions of bond angles. For more information, check [this blogpost](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#how-modeller-works).

The output of this program is a %(structure-db)s, which contains all of the modelled structures. Currently, the primary use of the %(structure-db)s is for interactive exploration with %(anvi-display-structure)s. You can also export your structures into external .pdb files with %(anvi-export-structures)s, or incorporate structural information in the %(variability-profile-txt)s with %(anvi-gen-variability-profile)s.

### Basic standard run

Here is a simple run: 

{{ codestart }}
anvi-gen-structure-database -c %(contigs-db)s \
                            --gene-caller-ids 1,2,3 \
                            -o STRUCTURE.db 
{{ codestop }}

Following this, you will have the structures for genes 1, 2, and 3 stored in `STRUCTURE.db`, assuming reasonable templates were found. Alternatively, you can provide a file name with the gene caller IDs (one ID per line) with the flag `--genes-of-interest`.  

If you have already run %(anvi-setup-pdb-database)s and therefore have a local copy of representative PDB structures, make sure you use it by providing the `--offline` flag. If you put it in a non-default location, provide the path to your %(pdb-db)s: 

{{ codestart }}
anvi-gen-structure-database -c %(contigs-db)s \
                            --gene-caller-ids 1,2,3 \
                            --pdb-database %(pdb-db)s \
                            -o STRUCTURE.db 
{{ codestop }}

To quickly get a very rough estimate for your structures, you can run with the flag `--very-fast`. 

### Basic import run

If you already possess structures and would like to create a %(structure-db)s for downstream anvi'o uses such as %(anvi-display-structure)s, you should create a %(external-structures)s file. Then, create the database as follows:

{{ codestart }}
anvi-gen-structure-database -c %(contigs-db)s \
                            --external-structures %(external-structures)s \
                            -o STRUCTURE.db 
{{ codestop }}

{:.notice}
Please avoid using any MODELLER-specific parameters when using this mode, as they will be silently ignored.


### Advanced Parameters

Here, we will go through a brief overview of the MODELLER parameters that you are able to change. See [this page](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#description-of-all-modeller-parameters) for more information. 

- The number of models to be simulated. The default is 1. 
- The standard deviation of atomic perturbation of the initial structure (i.e. how much you change the position of the atoms before fine tuning with other analysis). The default is 4 angstroms.
- The MODELLER database used. The default is `pdb_95`, which can be found [here](https://salilab.org/modeller/supplemental.html). This is the same database that is downloaded by %(anvi-setup-pdb-database)s.
- The scoring function used to compare potential models. The default is `DOPE_score`.
- The minimum percent identity cutoff for a template to be further considered.
- The minimum alignment fraction that the sequence is covered by the template in order to be further considered.
- The maximum number of templates that the program will consider. The default is 5. 
- The MODELLER program to use. The default is `mod9.19`, but anvi'o is somewhat intelligent and will
  look for the most recent version it can find.

For a case study on how some of these parameters matter, see [here](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#a-quick-case-study-on-the-importance-of-key-parameters). 

You also have the option to

- Skip the use of DSSP, which predicts beta sheets, alpha helices, certain bond angles, and relative
  solvent acessibility of residues.
- Output **all** the raw data, just provide a path to the desired directory with the flag `--dump-dir`.


