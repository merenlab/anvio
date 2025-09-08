This program clusters the contigs stored in a %(profile-db)s using your binning algorithm of choice and stores the results in several %(bin)ss. 

This is a quick alternative to manually binning your contigs, but it might miss some details that a human doing manual binning would find. After running this, you might want to run %(anvi-summarize)s on the resulting %(collection)s to look through your bins, and, if necessary, use %(anvi-refine)s to change the contents of them. 

You have to option to use several different clustering algorithms, which you'll specify with the `driver` parameter: [concoct](https://github.com/BinPro/CONCOCT/blob/develop/doc/source/index.rst), [metabat2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/), [maxbin2](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), [dastool](https://github.com/cmks/DAS_Tool), and [binsanity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5345454/). 

So, a run of this program will look like the following:

{{ codestart }}
anvi-cluster-contigs -p %(profile-db)s \
                     -c %(contigs-db)s \ 
                     -C %(collection)s \ 
                     --driver concoct
{{ codestop }}

Once you specify an algorithm, there are many algorithm specific parameters that you can change to your liking. When this program is set up, these parameters will appear in the help menu for the algorithms that anvi'o can find. 
