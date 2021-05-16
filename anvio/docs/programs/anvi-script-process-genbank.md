This program extracts the data from a %(genbank-file)s and converts it into anvi'o friendly artifacts: namely, a %(contigs-fasta)s, %(external-gene-calls)s and a %(functions-txt)s.

The %(contigs-fasta)s and %(external-gene-calls)s can be given to %(anvi-gen-contigs-database)s to create a %(contigs-db)s, and then you can use %(anvi-import-functions)s to bring the function data (in the %(functions-txt)s) into the database. Then you'll have all of the data in your %(genbank-file)s converted into a single %(contigs-db)s, which you can use for a variety of anvi'o analyses.

The parameters of this program entirely deal with the outputs. Besides telling the program where to put them, you can also give the function annotation source (in the %(functions-txt)s) a custom name. 
