This is a text file that **contains the information in a %(misc-data-items-order)s**, used for importing into and exporting this information from your anvi'o project. 

If you open this file, it will look something like this: 

    (contig_4, ((contig_1, contig_2), contig_3))
    
    
(but probably much more complicated.) If this is imported into an Anvi'o project, the contigs will be displayed in the order `contig_4, contig_1, contig_2, contig_3`, and the following tree will be generated in the interface: 
    
    contig_4    contig_1    contig_2    contig_3
        |           |           |           |
        |           -------------           |
        |                 |                 |
        |                 -------------------
        |                           |
        -----------------------------
                    |
                    |
