
This is a section of your %(contigs-db)s that contains custom additional information about specific amino acid residues.  

Take a look at [this blogpost](http://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database) for potential uses in the InteracDome (which will likely be added to anvi'o in v7) and the motivation behind this program.  

Similarly to other types of miscellaneous data (like %(misc-data-items)s), this information is either numerical or categorical and can be populated into a %(contigs-db)s (from a %(misc-data-amino-acids-txt)s) with %(anvi-import-misc-data)s. It is also displayed when you run %(anvi-show-misc-data)s and can be exported or deleted with %(anvi-export-misc-data)s and %(anvi-delete-misc-data)s respectively.  

For example, this could describe various key residues for binding to ligands, or residues otherwise determined to be important to the user for whatever reason.  

