An anvi’o database that contains **key information associated with your gene clusters**.

This is vital for its pangenomic analsis, hence the name. 

This is the output of the program %(anvi-pan-genome)s, which can be run after you've created a %(genomes-storage-db)s with the genomes you want to analyze. That script will identify your gene clusters and put them in a pan-db. 

Then, you can use the pan database to run a variety of other functions. You can think of it as a %(contigs-db)s focused on the gene clusters of multiple genomes. Some of these functions include %(anvi-compute-genome-similarity)s, %(anvi-analyze-synteny)s, and %(anvi-get-enriched-functions-per-pan-group)s. You can also view the data in a pan-db using %(anvi-display-pan)s. 

## Advanced information for programmers

While it is possible to read and write a given anvi'o pan database through SQLite functions directly, one can also use anvi'o libraries to initiate a pan database to read from. 

### Initiate a pan database instance

``` python
import argparse

from anvio.dbops import PanSuperclass

args = argparse.Namespace(pan_db="PAN.db", genomes_storage="GENOMES.db")

pan_db = PanSuperclass(args)

```

### Gene clusters dictionary

Once an instance from `PanSuperclass` is initiated, the following member function will give access to gene clusters:

``` pyton
pan_db.init_gene_clusters()
print(pan_db.gene_clusters)
```

```
{
  "GC_00000001": {
    "Genome_A": [19, 21],
    "Genome_B": [30, 32],
    "Genome_C": [122, 125],
    "Genome_D": [44, 42]
  },
  "GC_00000002": {
    "Genome_A": [123],
    "Genome_B": [176],
    "Genome_C": [175],
    "Genome_D": []
  },
  (...)
  "GC_00000036": {
    "Genome_A": [],
    "Genome_B": [24],
    "Genome_C": [],
    "Genome_D": []
  }
  (...)
```

Each item in this dictionary is a gene cluster describes anvi'o gene caller ids of each gene from each genome that contributes to this cluster.

### Sequences in gene clusters

```
gene_clusters_of_interest = set(["GC_00000006", "GC_00000036"])
gene_cluster_sequences = pan_db.get_sequences_for_gene_clusters(gene_cluster_names= gene_clusters_of_interest)

print(gene_cluster_sequences)
```

```
{
  "GC_00000006": {
    "Genome_A": {
      23: "MDVKKGWSGNNLND--NNNGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    },
    "Genome_B": {
      34: "MDVKKGWSGNNLND--NNNGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    },
    "Genome_C": {
      23: "MDVKKGWSGNNLNDWVNNNGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    },
    "Genome_D": {
      23: "MDVKKGWSGNNLNDWVNNAGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    }
  },
  "GC_00000036": {
    "Genome_A": {},
    "Genome_B": {
      24: "MSKRHKFKQFMKKKNLNPMNNRKKVGIILFATSIGLFFLFAFRTTYIVATGKVAGVSLKEKTA"
    },
    "Genome_C": {},
    "Genome_D": {}
  }
}
```
