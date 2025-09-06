This program opens an interactive interface to explore single amino acid variants (SAAVs) and single codon variants (SCVs) in the context of predicted tertiary protein structures and binding sites. There are many example applications [here](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#display-metagenomic-sequence-variants-directly-on-predicted-structures) and you can work through an example as part of [the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures) as well. This is an integral component of anvi'o structure, which you can learn more about [here](https://merenlab.org/software/anvio-structure/).


In short, this program enables users to explore sequence variation in the context of 3D protein structure, which reveals insights that cannot be obtained from purely sequence-based approaches.


### Before running 

To execute this program, you'll need to have created a %(structure-db)s which can be easily accomplished with a %(contigs-db)s and the program %(anvi-gen-structure-database)s.


You'll also need a %(profile-db)s that was created using %(anvi-profile)s's flag `--profile-SCVs`, which means that single codon variants (SCVs) have been profiled. While this may require re-profiling your data, as of v6.2, this is now a very efficient process.


### Basic Run 

There are two approaches for providing the variability information to this program.  

The first is to provide a %(contigs-db)s and %(profile-db)s pair, and let this program calculate SAAVs and SCVs as they are requested by the interface.


{{ codestart }}
anvi-display-structure -s %(structure-db)s \
                       -p %(profile-db)s \
                       -c %(contigs-db)s 
{{ codestop }}

The second approach is to use %(anvi-gen-variability-profile)s to create a %(variability-profile-txt)s. This method pre-loads all of the variability data and eliminates the need to wait for %(anvi-display-structure)s to calculate variability on-the-fly. This option is probably most convenient when you have already generated a %(variability-profile-txt)s for other purposes. If this applies to your situation, you can execute %(anvi-display-structure)s as follows:


{{ codestart }}
anvi-display-structure -s %(structure-db)s \
                       -c %(contigs-db)s \
                       -v %(variability-profile-txt)s
{{ codestop }}

{:.notice}
You must still provide the %(contigs-db)s used to generate the %(variability-profile-txt)s, since it contains other necessary information such as functional annotations and ligand binding predictions. You may optionally provide a %(profile-db)s if custom sample grouping is important to you.

{:.notice}
During %(anvi-gen-variability-profile)s execution, if you are _only_ interested in genes that have predicted structures, you may want to run %(anvi-gen-variability-profile)s with the flag `--only-if-structure`.

### Refining your search

Several options are available to refine which proteins and variants you're examining: 

- Provide a list of gene caller IDs to display only specific genes (this can be provided either directly as a parameter or as a file with one gene caller ID per line)
- Specify the minimum departure from the consensus sequence. This is a value from 0-1 that describes the threshold for a variability position to be displayed. For example, if this is set to 0.2, then all SAAVs and SCVs where less than 20 percent of the reads vary from the consensus sequence will not be displayed.
- Specify samples of interest. Those in your %(profile-db)s or %(variability-profile-txt)s that are not in the samples of interest will be filtered out.

If you're choosing to have %(anvi-display-structure)s calculate variability on-the-fly, you can accelerate the process by choosing to calculate _only_ SAAVs or _only_ SCVs.


### Other parameters 

Power users can also modify the server configuration (i.e., set the IP address, port number, browser path, server password, etc.)


