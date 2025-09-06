This program performs automated clustering of contigs stored in a %(profile-db)s using the specified binning algorithm and stores the results in multiple %(bin)ss. 

This approach provides a rapid alternative to manual binning procedures, though it may overlook subtle patterns that would be detected through careful manual curation. Following execution of this program, it is recommended to run %(anvi-summarize)s on the resulting %(collection)s to review the generated bins, and if necessary, use %(anvi-refine)s to modify their composition. 

The program supports several clustering algorithms, which are specified using the `driver` parameter: [concoct](https://github.com/BinPro/CONCOCT/blob/develop/doc/source/index.rst), [metabat2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/), [maxbin2](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), [dastool](https://github.com/cmks/DAS_Tool), and [binsanity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5345454/). 

A typical execution of this program follows this format:

{{ codestart }}
anvi-cluster-contigs -p %(profile-db)s \
                     -c %(contigs-db)s \ 
                     -C %(collection)s \ 
                     --driver concoct
{{ codestop }}

Each algorithm provides specific parameters that can be customized according to your requirements. When this program is properly configured, these algorithm-specific parameters will appear in the help menu for the algorithms that anvi'o can locate on your system.
