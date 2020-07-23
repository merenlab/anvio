This is a text file containing **the average coverage for each contig in each sample** that was in the %(profile-db)s and %(contigs-db)s that you used when you ran %(anvi-export-splits-and-coverages)s or %(anvi-export-gene-coverage-and-detection)s. 

This is a tab-delimited file where each row describes a specific split/gene and each column describes one of your samples. Each cell contains the average coverage of that contig in that sample. 

This artifact is really only used when taking information out of anvi'o, so enjoy your coverage information :) 

### Example for splits

(the type of output you would get from %(anvi-export-splits-and-coverages)s)

    contig                  sample_1    sample_2    sample_3 ...
    Day1_contig1_split1     5.072727    4.523432    1.2343243         
    Day1_contig1_split2     6.895844    5.284812    9.3721947
    Day1_contig2_split1     2.357049    3.519150    8.2385691
    ...


### Example for genes

(the type of output you would get from %(anvi-export-gene-coverage-and-detection)s)

    key       sample_1    sample_2    sample_3 ...
    13947     10.29109    1.984394    6.8289432         
    13948     34.89584    6.284812    3.3721947
    23026     23.94938    9.239235    13.238569
    ...




