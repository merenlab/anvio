This a tab-delimited text file that describes information contained in a %(misc-data-amino-acids)s. 

To import this information into a database, use %(anvi-import-misc-data)s. 

In this table, the first column should provide two pieces of information, both identifying a specific amino acid residue: the gene caller id (for the gene this residue is on) and its codon order in that gene. These should be separated by a colon. The following columns can contain any categorical or numerical data of your choosing.

Here is an example with very abstract data:

    item_name   categorical_data   numerical_data     data_group
      1:42            group_1          4.3245         cool_data 
      6:3             group_2          1.3542         cool_data
      9:96            group_1          3.2526         cool_data
      ...

There is another example of this table [here](http://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database). The second table on this page is what you would provide to %(anvi-import-misc-data)s, whereas the first lays out the same data more conceptually. 
