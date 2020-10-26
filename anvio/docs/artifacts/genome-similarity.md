This is the output of %(anvi-compute-genome-similarity)s and describes the level of similarity between all of the input genomes. 

{:.notice}
The output of %(anvi-compute-genome-similarity)s will only be in this structure if you did not input a %(pan-db)s. Otherwise, the data will be put directly into the additional data tables of the %(pan-db)s. 

This is a directory (named by the user) that contains both a %(dendrogram)s (NEWICK-tree) and a matrix of the similarity scores between each pair for a variety of metrics dependent on the program that you used to run %(anvi-compute-genome-similarity)s.

For example, if you used `pyANI`'s `ANIb` (the default program), the output directory will contain the following twelve files: 

-`ANIb_alignment_coverage.newick` and `ANIb_alignment_coverage.txt`

-`ANIb_full_percentage_identity.newick` and `ANIb_full_percentage_identity.txt` 

-`ANIb_hadamard.newick` and `ANIb_hadamard.txt`

-`ANIb_percentage_identity.newick` and `ANIb_percentage_identity.txt`

-`ANIb_similarity_errors.newick` and `ANIb_similarity_errors.txt`

