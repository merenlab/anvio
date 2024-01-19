This program imports a protein abundance profile, such as from proteomic experiments, into a %(profile-db)s.

This program takes as input a tab-delimited file of protein abundance data and a %(profile-db)s. The tabular file must have four columns with the following names: "source", "accession", "sample", and "abundance". Each row of the table corresponds to a distinct protein abundance measurement.

- "source" is the source of the protein accessions. It must be a gene function annotation source stored in the anvi'o %(profile-db)s (available sources can be found with the program, %(anvi-db-info)s).
- "accession" is the protein ID in the annotation source. A contigs database built from a GenBank file, for example, could contain the source, "NCBI_PGAP", and the accession, "WP_011862028.1".
- "sample" is the name of the sample in which the measurement was made. It need not be the same as any nucleotide sequence samples stored in the profile database.
- "abundance" is the protein abundance value, however defined.

Once protein abundances are stored in a profile database, they can be loaded into a metabolic %(reaction-network)s for analysis in the context of biochemical pathways.
