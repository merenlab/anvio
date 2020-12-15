This describes the output of %(anvi-script-calculate-pn-ps-ratio)s, which calculates the pN/pS ratio for each gene in a %(contigs-db)s. 

{:.notice}
See the page for %(anvi-script-calculate-pn-ps-ratio)s for an explanation of the pN/pS ratio 

This describes a directory that contains the following three files: 

`pN_pS_ratio.txt`: a matrix of the pN/pS values where each row is a gene and each column is a metagenome. 

    gene_caller_id  ANE_004_05M             ANE_004_05M             ANE_150_05M           ...
    1248            0.024404785965479868    0.024404785965479868    0.019650249960943812
    1249            0.014399205561072536    0.014399205561072536    0.0
    ...

`SAAV_counts.txt`: a matrix that describes the number of SAAVs for each gene in each metagenome.  
        
    gene_caller_id  ANE_004_05M     ANE_004_05M     ANE_150_05M      ...
    1248            12              11              10    
    1249            1               1               0
    ...
        
`SCV_counts.txt`: a matrix that describes the number of SCVs for each gene in each metagenome.  

    gene_caller_id  ANE_004_05M     ANE_004_05M     ANE_150_05M      ...
    1248            143             154             148    
    1249            19              19              14
    ...
