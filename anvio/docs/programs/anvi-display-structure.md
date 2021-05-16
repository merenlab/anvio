
This program opens an interactive interface to explore single amino acid variants (SAAVs) and single codon variants (SCVs) in the context of predicted tertiary protein structures and binding sites.  There are many example uses [here](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#display-metagenomic-sequence-variants-directly-on-predicted-structures) and you can work through an example as part of [the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures) as well.  This is an integral program of anvi'o structure, which you can learn more about [here](https://merenlab.org/software/anvio-structure/).


In short, this program enables users to explore sequence variation in the context of 3D protein structure, which reveals insight that cannot be learned from purely sequence-based approaches.


### Before running 

To run this program, you'll need to have created a %(structure-db)s which can be easily done with a %(contigs-db)s and the program %(anvi-gen-structure-database)s.


You'll also need a %(profile-db)s that was created using %(anvi-profile)s's flag `--profile-SCVs`, which means that single codon variants (SCVs) have been profiled. Very sorry if this forces you to re-profile, but as of v6.2, this is now a very expedient process.


### Basic Run 

There are two ways to provide the variability information to this program.  

The first is to provide a %(contigs-db)s and %(profile-db)s pair, and let this program calculate SAAVs and SCVs as they are requested by the interface.


{{ codestart }}
anvi-display-structure -s %(structure-db)s \
                       -p %(profile-db)s \
                       -c %(contigs-db)s 
{{ codestop }}

The second is to use %(anvi-gen-variability-profile)s to create a %(variability-profile-txt)s. This way, you pre-load all of the variability data and don't have to wait for %(anvi-display-structure)s to calculate variability on-the-fly. This option is probably most convenient in instances where you have already generated a %(variability-profile-txt)s for other reasons. If you fall into this camp, you can run %(anvi-display-structure)s as so:


{{ codestart }}
anvi-display-structure -s %(structure-db)s \
                       -c %(contigs-db)s \
                       -v %(variability-profile-txt)s
{{ codestop }}

{:.notice}
You still must provide the %(contigs-db)s used to generate the %(variability-profile-txt)s, since it contains other necessary information such as functional annotations and ligand binding predictions.  You may optionally provide a %(profile-db)s if custom sample grouping is important to you.

{:.notice}
During %(anvi-gen-variability-profile)s, if you are _only_ interested in genes that have predicted structures, you may want to run %(anvi-gen-variability-profile)s with the flag `--only-if-structure`.

### Refining your search

You have several options to refine what proteins and variants you're looking at: 

- Provide a list of gene caller IDs to only display specific genes (this can be provided either directly as a parameter or as a file with one gene caller ID per line)
- Specify the minimum departure from the consensus sequence. This is a number from 0-1 that describes the threshold for a variability position to be displayed. For example, if this is set to 0.2, then all SAAVs and SCVs where less than 20 percent of the reads vary from the consensus sequence will not be displayed.
- Specify samples of interest. Those in your %(profile-db)s or %(variability-profile-txt)s that are not in the samples of interest will be filtered out.

If you're choosing to have %(anvi-display-structure)s calculate variability on-the-fly, you can speed things up by choosing to _only_ calculate SAAVs or _only_ calculate SCVs.


### Other parameters 

Power users can also change the server configuration (i.e. set the IP address, port number, browser path, server password, etc.)


