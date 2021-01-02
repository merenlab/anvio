As an artifact, this describes the variability information about a single sample calculated when you ran %(anvi-profile)s. To examine variability across samples, you'll want to use this information (which is stored within your %(profile-db)s) to run %(anvi-gen-variability-profile)s. 

## Details about Variability

In the context of anvi'o, variability means divergence of environmental populations from the reference used to perform metagenomic read recruitment.

Here, the term "population" describes an assemblage of co-existing microbial genomes in an environment that are similar enough to map to the context of the same reference genome.

The variability profile of a metagenome enables studies of [microbial population genetics with anvi'o](http://merenlab.org/2015/07/20/analyzing-variability/).

There are two types of variability the program %(anvi-profile)s can characterize and store: substitutions and indels.

### Substitutions: SNVs, SCVs, SAAVs

Anvi'o can make sense of single-nucleotide variants (SNVs), single-codon variants (SCVs), and single-amino acid variants (SAAVs). See [this article](http://merenlab.org/2015/07/20/analyzing-variability) for more information.

You can learn the name of the table in which anvi'o stores this in a given %(profile-db)s by running this command in your anvi'o environment:

``` bash
python -c 'import anvio.tables as t; print(t.variable_nts_table_name)'
```

This will tell you about its structure:

``` bash
python -c 'import anvio.tables as t; print(t.variable_nts_table_structure)'
```

### Indels: insertions and deletions

Anvi'o can also characterize insertions and deletions found within an environment based on short-read recruitment results and will store in the following table:

``` bash
python -c 'import anvio.tables as t; print(t.indels_table_name)'
```

**Notes for programmers**: The convention for the start position of an insertion is defined like so:

```
    pos_in_contig ...0123456 7890123456
    reference     ...CTACTAC TACTTCATGA...
    read              TACTAC TAC
    insertion               └──ACTG
```

In this case, the start position of the insertion in the contig is 6. The insertion _follows_ the position it is defined by. This is opposite to IGV, in which the insertion _precedes_ the position it is defined by.

For deletions, there is no such ambiguity in the start position, since the deletion starts on a reference position, not in between two reference positions.
