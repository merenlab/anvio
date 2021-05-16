This program converts a %(linkmers-txt)s artifact into %(oligotypes)s data.

A %(linkmers-txt)s artifact describes each of your short reads that mapped to specific target nucleotide positions in a reference contig. This program counts the total occurance of each combination in those target positions within each of your samples. 

For example, if your %(linkmers-txt)s focused on two target positions, and you ran the following:

{{ codestart }}
anvi-oligotype-linkmers -i %(linkmers-txt)s 
{{ codestop }}

The output (which by default is called `oligotype-counts-001.txt`) might look like the following:

    key         AG   CA    CG    GA    GG    TA    TG   
    sample_001  0    320   12    2     0     3     579    
    sample_002  0    142   2     0     2     10    353  
    sample_003  3    404   1     1     0     2     610   
    sample_004  0    209   6     0     1     0     240

Note that combinations with zero reads in every sample are not included. 
