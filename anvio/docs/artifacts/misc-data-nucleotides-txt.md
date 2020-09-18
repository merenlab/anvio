This a tab-delimited text file that describes information contained in a %(misc-data-nucleotides)s. 

To import this information into a database, use %(anvi-import-misc-data)s. 

In this table, the first column should provide two pieces of information, both identifying a specific nucleotide position: the name of the contig the nucleotide is on, and its position on that contig. These should be separated by a colon. The following columns can contain any categorical or numerical data of your choosing.

Here is an example with very abstract data:

    item_name   categorical_data  numerical_data    data_group
    contig_1:4       group_1          4.3245         cool_data 
    contig_4:72      group_2          1.3542         cool_data
    contig_7:24      group_1          3.2526         cool_data
    ...

For a more concrete example, check out the example table for %(misc-data-amino-acids)s (which is formatted very similarly)  [here](http://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database). The second table on this page is what you would provide to %(anvi-import-misc-data)s. 
