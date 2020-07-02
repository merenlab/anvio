Stores %(hmm-hits)s for a given %(hmm-source)s in a %(contigs-db)s. In short, this is the program that will do a search for HMMs against a %(contigs-db)s and store that information into the contigs-db's %(hmm-hits)s. 

This is one of the programs that users commonly run on newly generated %(contigs-db)s, along with %(anvi-scan-trnas)s, %(anvi-run-ncbi-cogs)s, %(anvi-run-scg-taxonomy)s, and so on.

In short, this step will help annotate the genes in your %(contigs-db)s.

### Default Usage

To run this program with all default settings (against all default anvio %(hmm-source)s), you only need to provide a %(contigs-db)s. 

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB 
{{ codestop }}

### Running against a custom set of %(hmm-source)s

In order to run against your own hmm-source or a custom subset of anvi'o's hmm-sources, you have two choices. 

#### Choice 1: I have my own hmm-sources on my computer

This way the source can be completely outside of anvi'o. 

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB -H path_to_your_hmm_profile
{{ codestop }}

#### Choice 2: I prefer anvi'o's hmm-sources, but I don't need all of them.

By default, anvi'o will look through all of its hmm-sources when doing a search. If you only want to run against a specific one, you're in the right place. These are the currently available ones: "Bacteria_71" (type: singlecopy), "Archaea_76" (type: singlecopy), "Protista_83" (type: singlecopy), and "Ribosomal_RNAs" (type: Ribosomal_RNAs). See the page for %(hmm-source)s for more information. 

For example, 

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB -I Bacteria_71 
{{ codestop }}

### Other things anvi-run-hmms can do

- Add the tag `--also-scan-trnas` to basically run %(anvi-scan-trnas)s for you at the same time. It's very convientient. 
- Add the tag `--just-do-it` to hide all warnings and questions in case you don't want to deal with those.
-  There are also parameters that can help speed up the runtime of this program. However, be aware of the limits of your system, espeically if running on a SGE.  For example, you can increase the number of threads or switch to hmmsearch if you are scanning  a large umber of HMMs. For more information on that, check out [here](http://merenlab.org/software/anvio/vignette/#anvi-run-hmms). 

