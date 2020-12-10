
This program (much like all of the other programs that begin with `anvi-setup`) sets up a local copy of the InteracDome database for %(anvi-run-interacdome)s as well as a local copy of Pfam v31.0, which is what InteracDome is defined for. Note that anvi'o only needs this program to be run once.


Specifically, this downloads [InteracDome](https://interacdome.princeton.edu/)’s [tab-separated files](https://interacdome.princeton.edu/#tab-6136-4) and the Pfam v31.0 HMM profiles for the Pfams in your InteracDome data. This data is stored in the %(interacdome-data)s artifact. 


It's easy as 1-2-3:

{{ codestart }}
anvi-setup-interacdome
{{ codestop }}

When running this program, you can provide a path to store your InteracDome data in. The default path is `anvio/data/misc/InteracDome`; if you use a custom path, you will have to provide it to %(anvi-run-interacdome)s with the same parameter. Here is an example run: 


{{ codestart }}
anvi-setup-interacdome --interacdome-data-dir path/to/directory 
{{ codestop }}

If you want to overwrite any data that you have already downloaded (for example if you suspect something went wrong in the download), add the `--reset` flag: 

{{ codestart }}
anvi-setup-interacdome  --interacdome-data-dir path/to/directory \ 
                        --reset
{{ codestop }}

