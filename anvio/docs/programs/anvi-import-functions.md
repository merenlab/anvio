This program **takes in a %(functions-txt)s to annotate your %(contigs-db)s with %(functions)s.** Basically, if you have already have the gene functions for the contigs in your %(contigs-db)s available in a file, you can import them into anvi'o using this command. 

You can find a really comprehesive walkthrough of this program on [this blog post about importing functions](http://merenlab.org/2016/06/18/importing-functions/), including information about built-in anvi'o parsers for InterProScan and the EggNOG database.

If you want to overwrite any function annotations you already have, just add the tag `--drop-previous-annotations`. 

