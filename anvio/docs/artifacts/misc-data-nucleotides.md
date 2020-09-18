This is a section of your %(contigs-db)s that contains custom additional information about specific nucleotides. 

Take a look at [this blogpost](http://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database) for potential uses in the interactome (which will likely be added to anvi'o in v7) and the motivation behind this program. 

Similarly to other types of miscellaneous data (like %(misc-data-items)s), this information can be populated into a %(contigs-db)s with %(anvi-import-misc-data)s and can be either numerical or categorical. 

For example, this information could describe specific nucleotides that are known to be SNVs from another experiment, various key nucleotides for binding to ligands, or positions known to have other modifications (such as m1A or s4U).
