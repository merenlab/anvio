This program **generates a %(genes-db)s, which stores the coverage and detection values for all of the genes in your %(contigs-db)s.** 

This information is usually calculated when it's needed (for example when running %(anvi-interactive)s in genes mode), but this program lets you break this process into two steps. This way, you can easily change the parameters of %(anvi-interactive)s without having to recalculate the gene-level statistics. 

Given a %(contigs-db)s and %(profile-db)s pair, as well as a %(collection)s, this program will calculate the stats for the genes in each of your %(bin)ss and give each bin its own %(profile-db)s that includes this information. 

For example, if a %(collection)s called `GENE_COLLECTION` contained the bins `bin_0001`, `bin_0002`, and `bin_0003` and you ran:

{{ codestart }}
anvi-gen-gene-level-stats-databases -c %(contigs-db)s \
                                    -p %(profile-db)s \
                                    -C %(collection)s 
{{ codestop }}

Then it will create a directory called `GENES` that contains three %(profile-db)s called `GENE_COLLECTION-bin_0001.db`, `GENE_COLLECTION-bin_0002.db`, and `GENE_COLLECTION-bin_0003.db`. In terms of output, this program is similar to %(anvi-split)s: each of these databases can now be treated as self-contained anvi'o projects but they also contain the gene-level information. Thus, you then could run %(anvi-interactive)s in genes mode on one of these profile databases. 

You also have the option to provide a list of %(bin)s (either as a file or as a string) to anlyze instead of a single %(collection)s. 

### Other Parameters

You can also change the definition of an outlier nucleotide position or switch calculations to use the [INSeq/Tn-Seq](https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/in-seq-tn-seq.html) statistical methods. 
