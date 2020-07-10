This artifact **contains information about the functions of its contigs in a tab delimated text file.**

This file is formatted with a single gene per row with the following columns: 
1. The gene caller ID
2. The source (the database that you got this function data from)
3. The assession number of the gene (optional)
4. The function information
1. The e-value (optional, but helpful if you plan to use the interactive interface and have more than one function for a single gene)

For an example, check out [this lovely page](http://merenlab.org/2016/06/18/importing-functions/#simple-matrix). 

This is primarily used to import and export function information from an anvi'o project. See [this page](http://merenlab.org/2016/06/18/importing-functions/) for more information on importing function information. 

It is also the output of %(anvi-search-functions)s which searches for specific terms in your functional annotations. For example, you could search for all genes related to "ribosomes" and get a functions-txt of all of those genes. 
