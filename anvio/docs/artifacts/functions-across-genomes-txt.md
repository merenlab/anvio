A TAB-delimited output file that describes the distribution of functions across a group of genomes.

## General format

With this file type anvi'o may describe the presence-absence of individual functions across genomes, or the frequency of functions. The first column contains a unique `key` for each function, and the second column the annotation generated from the annotation source the user chose. The third column is optional and includes the accession ids. The remaining columns hold one value per (meta)genome.

### Frequency of function hits across genomes

Here is an example output file that shows the frequency of functions across 9 _Bifidobacterium_ genomes:

|**key**|**COG20_FUNCTION**|**B_adolescentis_1_11**|**B_adolescentis_22L**|**B_adolescentis_6**|**B_lactis_Bl_04**|**B_lactis_CNCM_I_2494**|**B_lactis_DSM_10140**|**B_longum_JDM301**|**B_longum_KACC_91563**|**B_longum_NCIMB8809**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|`func_3874a4219ffa`|Chromosomal replication initiation ATPase DnaA (DnaA) (PDB:1L8Q)|1|1|1|1|1|1|1|1|1|
|`func_1530552f6b61`|DNA polymerase III sliding clamp (beta) subunit, PCNA homolog (DnaN) (PDB:1JQJ)|1|1|1|1|1|1|1|1|1|
|`func_060988dd6d3a`|Recombinational DNA repair ATPase RecF (RecF) (PDB:5Z67)|1|1|1|1|1|1|1|1|1|
|`func_46477a763d38`|Predicted  nucleic acid-binding protein, contains Zn-ribbon domain (includes truncated derivatives)|1|1|1|1|1|1|1|1|1|
|`func_184ca83fe067`|DNA gyrase/topoisomerase IV, subunit B (GyrB) (PDB:1EI1)|2|2|2|2|2|2|2|2|2|
|`func_2ee67347f04b`|DNA gyrase/topoisomerase IV, subunit A (GyrA) (PDB:1SUU)|2|2|2|2|2|2|2|2|2|
|`func_ff160914972f`|Predicted ATPase, archaeal AAA+ ATPase superfamily|1|0|0|0|0|0|3|1|1|
|`func_598325c18ceb`|Molybdopterin or thiamine biosynthesis adenylyltransferase (ThiF) (PDB:1ZUD) (PUBMED:32239579)|1|0|0|0|0|0|1|1|1|
|`func_8ce5ff84aa42`|Predicted arabinose efflux permease AraJ, MFS family (AraJ) (PDB:4LDS)|14|13|13|10|10|10|22|17|16|
|`func_300e2c6e37e4`|Glutamate dehydrogenase/leucine dehydrogenase (GdhA) (PDB:1B3B) (PUBMED:24391520)|1|1|1|1|1|1|1|1|1|
|`func_a666bc87a03f`|Large-conductance mechanosensitive channel (MscL) (PDB:2OAR)|1|1|1|1|1|1|1|1|1|
|`func_71379d89f0c6`|UTP pyrophosphatase, metal-dependent hydrolase family (YgjP) (PDB:4JIU) (PUBMED:27941785)|1|1|1|1|1|1|1|1|1|
|`func_772efbb7cb2e`|Hemolysin-related protein, contains CBS domains, UPF0053 family (TlyC) (PDB:2NQW)|2|2|2|2|2|2|2|2|2|
|`func_5a5bf9735c8a`|Carbonic anhydrase (CynT) (PDB:1EKJ) (PUBMED:22081392)|1|1|1|1|1|1|1|1|1|
|`func_1ffc26d0034f`|Transposase (or an inactivated derivative) (IS285)|4|4|1|3|3|3|12|13|2|
|`func_0631d9fecae6`|Phosphoenolpyruvate carboxylase (Ppc) (PDB:1FIY)|2|1|1|1|1|1|2|1|1|
|`func_62068e09f9f3`|Chromosome segregation ATPase Smc (Smc) (PDB:5XG3)|2|1|1|0|0|0|1|1|1|
|`func_f4b8779f6ad8`|Peptidoglycan-binding (PGRP) domain of peptidoglycan hydrolases (PGRP) (PDB:4FET)|2|0|0|0|0|0|0|0|0|
|`func_700652d0a2dd`|Uncharacterized membrane protein YjjP, DUF1212 family (YjjP)|1|1|1|1|1|1|1|1|1|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Presence/absence of function hits across genomes

In contrast, here is the presence/absence report for the same data:

|**key**|**COG20_FUNCTION**|**B_adolescentis_1_11**|**B_adolescentis_22L**|**B_adolescentis_6**|**B_lactis_Bl_04**|**B_lactis_CNCM_I_2494**|**B_lactis_DSM_10140**|**B_longum_JDM301**|**B_longum_KACC_91563**|**B_longum_NCIMB8809**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|`func_3874a4219ffa`|Chromosomal replication initiation ATPase DnaA (DnaA) (PDB:1L8Q)|1|1|1|1|1|1|1|1|1|
|`func_1530552f6b61`|DNA polymerase III sliding clamp (beta) subunit, PCNA homolog (DnaN) (PDB:1JQJ)|1|1|1|1|1|1|1|1|1|
|`func_060988dd6d3a`|Recombinational DNA repair ATPase RecF (RecF) (PDB:5Z67)|1|1|1|1|1|1|1|1|1|
|`func_46477a763d38`|Predicted  nucleic acid-binding protein, contains Zn-ribbon domain (includes truncated derivatives)|1|1|1|1|1|1|1|1|1|
|`func_184ca83fe067`|DNA gyrase/topoisomerase IV, subunit B (GyrB) (PDB:1EI1)|1|1|1|1|1|1|1|1|1|
|`func_2ee67347f04b`|DNA gyrase/topoisomerase IV, subunit A (GyrA) (PDB:1SUU)|1|1|1|1|1|1|1|1|1|
|`func_ff160914972f`|Predicted ATPase, archaeal AAA+ ATPase superfamily|1|0|0|0|0|0|1|1|1|
|`func_598325c18ceb`|Molybdopterin or thiamine biosynthesis adenylyltransferase (ThiF) (PDB:1ZUD) (PUBMED:32239579)|1|0|0|0|0|0|1|1|1|
|`func_8ce5ff84aa42`|Predicted arabinose efflux permease AraJ, MFS family (AraJ) (PDB:4LDS)|1|1|1|1|1|1|1|1|1|
|`func_300e2c6e37e4`|Glutamate dehydrogenase/leucine dehydrogenase (GdhA) (PDB:1B3B) (PUBMED:24391520)|1|1|1|1|1|1|1|1|1|
|`func_a666bc87a03f`|Large-conductance mechanosensitive channel (MscL) (PDB:2OAR)|1|1|1|1|1|1|1|1|1|
|`func_71379d89f0c6`|UTP pyrophosphatase, metal-dependent hydrolase family (YgjP) (PDB:4JIU) (PUBMED:27941785)|1|1|1|1|1|1|1|1|1|
|`func_772efbb7cb2e`|Hemolysin-related protein, contains CBS domains, UPF0053 family (TlyC) (PDB:2NQW)|1|1|1|1|1|1|1|1|1|
|`func_5a5bf9735c8a`|Carbonic anhydrase (CynT) (PDB:1EKJ) (PUBMED:22081392)|1|1|1|1|1|1|1|1|1|
|`func_1ffc26d0034f`|Transposase (or an inactivated derivative) (IS285)|1|1|1|1|1|1|1|1|1|
|`func_0631d9fecae6`|Phosphoenolpyruvate carboxylase (Ppc) (PDB:1FIY)|1|1|1|1|1|1|1|1|1|
|`func_62068e09f9f3`|Chromosome segregation ATPase Smc (Smc) (PDB:5XG3)|1|1|1|0|0|0|1|1|1|
|`func_f4b8779f6ad8`|Peptidoglycan-binding (PGRP) domain of peptidoglycan hydrolases (PGRP) (PDB:4FET)|1|0|0|0|0|0|0|0|0|
|`func_700652d0a2dd`|Uncharacterized membrane protein YjjP, DUF1212 family (YjjP)|1|1|1|1|1|1|1|1|1|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Per-population copy number of functions across metagenomes

When the program that generates this file is run with the `--add-per-population-copy-number` flag, anvi'o additionally reports the per-population copy number (PPCN) of each function. In this file, the frequency of a function in a given metagenome is divided by the number of populations estimated to be present in that metagenome (which anvi'o estimates from single-copy core genes). This normalization makes the values comparable across metagenomes that describe communities of different sizes.

The layout is identical to the frequency file, except that the values are per-population copy numbers rather than raw counts. If anvi'o was unable to estimate the number of populations for a given metagenome, its values in this file will be reported as `NA`.
