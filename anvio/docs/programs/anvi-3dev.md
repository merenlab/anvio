This program opens the interactive interface to visualize variable positions directly on the 3D structure of a protein. There are many example uses [here](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#display-metagenomic-sequence-variants-directly-on-predicted-structures) and you can work through an example as part of [the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures). 

In short, though, this interface lets you view the predicted 3D structures of your protein (as stored in a %(structure-db)s) with the variability positions (SCVs and SAAVs, usually determined from your metagenomic data) directly mapped on. This can give you new insights, for example the solvent accessibility of individual SAAVs and the strucutral distribution of variability positions. This is espeically useful after using %(anvi-import-misc-data)s to annotate additional data into your %(contigs-db)s, such as binding sites for substrates or other enzymes. 

### Before running

To run this program, you'll need to have used your %(contigs-db)s created a %(structure-db)s with %(anvi-gen-structure-database)s. 

You'll also need a %(profile-db)s that has had its SCVs profiled. In other words, when you ran either %(anvi-profile)s or %(anvi-merge)s, you need to have added the flag `--profile-SCVs` for this to work. This may be computationally intensive, but it is also necessary to run anvi-3dev. 

### Basic Run

There are two ways to provide the variability information to this program. 

The first is to provide a %(contigs-db)s and %(profile-db)s pair, and let this program calculate the variability positions for you in the moment. 

{{ codestart }}
anvi-3dev -s %(structure-db)s \
          -p %(profile-db)s \
          -c %(contigs-db)s 
{{ codestop }}

The second is to use %(anvi-gen-variability-profile)s to create a %(variability-profile-txt)s. This way, you don't have to wait for %(anvi-3dev)s to run all of this analysis and it is easier to share your variaiblity information. 

For this purpose, you'll probably want to run %(anvi-gen-variability-profile)s with the flag `--only-if-structure` so that it only calculates varaibility proteins that can be visualized. Then you can run %(anvi-3dev)s as so:

{{ codestart }}
anvi-3dev -s %(structure-db)s \
          -v %(variability-profile-txt)s
{{ codestop }}

### Refining your search 

You have several options to refine what proteins and variable positions you're looking at: 

- Provide a list of gene caller IDs to only display specific genes (this can be provided either directly as a parameter or as a file with one gene caller ID per line)
- Provide a %(splits-txt)s to only look at specfic splits (though usually you'll want to refine your search by genes)
- Specify the minimum departure from the consensus sequence. This is a number from 0-1 that describes the threshold for a variability position to be displayed. For example, if this is set to 0.2, then all SAAVs and SCVs where less than 20 percent of the reads vary from the consensus sequence will not be displayed. 

If you're choosing to not use a %(variability-profile-txt)s and have anvi-3dev calculate your variability, you can also change a few other parameters. If you are using a %(variability-profile-txt)s, you can instead set these parameters when running %(anvi-gen-variability-profile)s
- Provide a list of samples to use to calculate variable positions 
- Specify whether to look speicfically at SCVs or SAAVs. 

### Other parameters

Power users can also change the server configuration (i.e. set the IP address, port number, browser path, server password, etc.)
