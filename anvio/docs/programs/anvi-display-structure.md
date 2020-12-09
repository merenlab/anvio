This program opens an interactive interface to explore sequence variants (SAAVs and SCVs) in the
context of tertiary protein structure and predicted binding sites. There are many example uses
[here](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#display-metagenomic-sequence-variants-directly-on-predicted-structures)
and you can work through an example as part of [the infant gut
tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures).
This is an integral program of anvi'o structure, which you can learn more about
[here](https://merenlab.org/software/anvio-structure/).

In short, this program enables users to explore sequence variation in the context of 3D protein
structure, which reveals insight that cannot be learned from purely sequence-based approaches.

### Before running

To run this program, you'll need to have created a %(structure-db)s which can be easily done with a
%(contigs-db)s and the program %(anvi-gen-structure-database)s. 

You'll also need a %(profile-db)s that was run with `--profile-SCVs`, which means that single codon
variants (SCVs) have been profiled. Very sorry for this inconvenience.


### Basic Run

There are two ways to provide the variability information to this program. 

The first is to provide a %(contigs-db)s and %(profile-db)s pair, and let this program calculate the variability positions for you in the moment. 

{{ codestart }}
anvi-display-structure -s %(structure-db)s \
          -p %(profile-db)s \
          -c %(contigs-db)s 
{{ codestop }}

The second is to use %(anvi-gen-variability-profile)s to create a %(variability-profile-txt)s. This way, you don't have to wait for %(anvi-display-structure)s to run all of this analysis and it is easier to share your variaiblity information. 

For this purpose, you'll probably want to run %(anvi-gen-variability-profile)s with the flag `--only-if-structure` so that it only calculates varaibility proteins that can be visualized. Then you can run %(anvi-display-structure)s as so:

{{ codestart }}
anvi-display-structure -s %(structure-db)s \
          -v %(variability-profile-txt)s
{{ codestop }}

### Refining your search 

You have several options to refine what proteins and variable positions you're looking at: 

- Provide a list of gene caller IDs to only display specific genes (this can be provided either directly as a parameter or as a file with one gene caller ID per line)
- Provide a %(splits-txt)s to only look at specfic splits (though usually you'll want to refine your search by genes)
- Specify the minimum departure from the consensus sequence. This is a number from 0-1 that describes the threshold for a variability position to be displayed. For example, if this is set to 0.2, then all SAAVs and SCVs where less than 20 percent of the reads vary from the consensus sequence will not be displayed. 

If you're choosing to not use a %(variability-profile-txt)s and have anvi-display-structure calculate your variability, you can also change a few other parameters. If you are using a %(variability-profile-txt)s, you can instead set these parameters when running %(anvi-gen-variability-profile)s
- Provide a list of samples to use to calculate variable positions 
- Specify whether to look speicfically at SCVs or SAAVs. 

### Other parameters

Power users can also change the server configuration (i.e. set the IP address, port number, browser path, server password, etc.)
