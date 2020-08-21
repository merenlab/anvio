The main function of `anvi-merge` is to convert multiple %(single-profile-db)ss into a single %(profile-db)s (also called a merged profile database). Basically, this takes the alignment data from each sample (each contained in its own %(single-profile-db)s) and combines them into a single database that anvi'o can look through more easily. 

### Overview: How to run anvi-merge

1. Set up your %(contigs-db)s. See that page for more information

1. Use %(anvi-profile)s to create a %(single-profile-db)s for each of your samples (formatted into a %(bam-file)s) *(Note: for each of these runs, you'll need to use the same %(contigs-db)s and parameters)*

1. Use `anvi-merge` to merge those %(single-profile-db)ss into a single database, called a %(profile-db)s. This will look something like the following:

{{ codestart }}
anvi-merge -c cool_contigs.db \
            Single_profile_db_1 Single_profile_db_2 \
            -o cool_contigs_merge
{{ codestop }}
                    
This will put all of the output files (the final %(profile-db)s as well as a %(misc-data-items-order)s which is the result of your hierarchical clustering and describes the order to display your contigs in) into the folder `cool_contigs_merge `.
    

## Other Parameters

You must give `anvi-merge` your contigs database and single profile databases. However, you can also provide more information or give addtional instructions. Use the flag `-h` at any time to display the help menu.

### Hierarchical Clustering 

#### To run or not to run? 
* Use the flag `--skip-hierarchical-clustering` to turn hierarchical clustering off. This will save on computation time, but will skip out on creating the tree of contigs at the center of the interactive interface. If you have more than 25,000 splits in the final profile, this will be set automatically. 
* Use the flag `--enforce-hierarchical-clustering` to turn hierarchical clustering back on. This will take a long time, but will produce a lovely contigs tree for the interactive interface. 

#### Additional parameters
* Provide a custom distance metric for clustering using the flag `--distance.` (The default is euclidean)
* Provide a custom linkage method for clustering using the flag `--linkage.` (The default is ward)

### Providing additional information
* Provide the sample name using the flag `-S`. If you don't, anvi'o will come up with one, but it probably won't be any good. 
* Provide a text file in markdown to describe the project using the flag `--description`. This will show up when you later use the interactive interface to analyze your profile-db. 

### Output Information
* Provide an output destination with the flag `-o`.
* Add the flag `-W` to overwrite existing files in that directory. 
