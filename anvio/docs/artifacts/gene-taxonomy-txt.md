This is a file containing **the taxonomy information for the genes in your %(contigs-db)s**. 

You can use %(anvi-import-taxonomy-for-genes)s to integrate this information into your contigs database. See [this blog post](http://merenlab.org/2016/06/18/importing-taxonomy/) for a comprehensive tutorial. 

In its simplest form, this file is a tab-delimited text file that lists gene caller IDs and their associated taxonomy information. However, anvi'o can also parse outputs from taxonomy-based software like [Kaiju](https://github.com/bioinformatics-centre/kaiju) or [Centrifuge](https://github.com/infphilo/centrifuge). 

For example:

    gene_caller_id  t_domain     t_phylum       t_class      ...
          1         Eukarya      Chordata       Mammalia
          2         Prokarya     Bacteroidetes  Bacteroidia
          ...


