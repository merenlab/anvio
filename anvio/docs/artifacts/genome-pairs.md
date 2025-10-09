This input file type describes pairs of genomes that go together in an analysis, and is usually paired with an %{external-genomes}s file that describes where the %{contigs-db} associated with each genome name is located.

In the context of %(anvi-predict-metabolic-exchanges)s, this file describes which pairs of genomes to predict exchanges between.

The file should be tab-delimited and contain at least two columns named `genome_1` and `genome_2`. Here is an example:

|**`genome_1`**|**`genome_2`**|
|:--|:--|
|name_of_one_genome|name_of_another_genome|
|E_coli|Pelagibacter_sp|
|another_genome|awesomegenome|

As long as the genome names match to those described in the accompanying %{external-genomes}s file, you should be good to go.