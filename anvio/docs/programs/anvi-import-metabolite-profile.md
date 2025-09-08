This program imports a metabolite abundance profile, such as from metabolomic experiments, into a %(profile-db)s.

This program takes as input a tab-delimited file of metabolite abundance data and a %(profile-db)s. The tabular file must have three columns with the following names: "accession", "sample", and "abundance". Each row of the table corresponds to a distinct metabolite abundance measurement.

- "accession" is the ModelSEED Compound ID, e.g., "cpd00027" for D-glucose.
- "sample" is the name of the sample in which the measurement was made. It need not be the same as any nucleotide sequence samples stored in the profile database.
- "abundance" is the metabolite abundance value, however defined.

Once metabolite abundances are stored in a profile database, they can be loaded into a metabolic %(reaction-network)s for analysis in the context of biochemical pathways. Metabolites in the network are defined in terms of ModelSEED Compounds.
