This ddatabase contains the **predicted 3D structure data** of the sequences encoded by contigs in a %(contigs-db)s. It is the result of running %(anvi-gen-structure-database)s (which uses MODELLER to predict your protein strucutres initially based on alignment to sequences with known strucutres). 

You can use this database to generate a variability profile (with %(anvi-gen-variability-profile)s) and then use that data (along with this database) to visualize your 3D structures and their relation to SCVs and SAAVs with %(anvi-display-structure)s. 

Besides this, you can export the structure data into an external `.pdb` file (using %(anvi-export-structures)s) or generation the fixation index matrix (with %(anvi-gen-fixation-index-matrix)s). 

For more information on the structure database, see [this blog post](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#the-structure-database). 

