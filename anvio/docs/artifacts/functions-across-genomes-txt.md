A TAB-delimited output file that describes the distribution of functions across a group of genomes.

## General format

With this file type anvi'o may describe the presence-absence of individual functions across genomes, or the frequency of functions. The header of the last column in a %(functions-across-genomes-txt)s file will be the name of the function annotation source from which the gene functions were recovered.

### Frequency of function hits across genomes

Here is an example output file that shows the frequency of functions across 9 _Bifidobacterium_ genomes:

|**key**|**B_adolescentis_1_11**|**B_adolescentis_22L**|**B_adolescentis_6**|**B_lactis_Bl_04**|**B_lactis_CNCM_I_2494**|**B_lactis_DSM_10140**|**B_longum_JDM301**|**B_longum_KACC_91563**|**B_longum_NCIMB8809**|**COG20_FUNCTION**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|`func_3874a4219ffa`|1|1|1|1|1|1|1|1|1|Chromosomal replication initiation ATPase DnaA (DnaA) (PDB:1L8Q)|
|`func_1530552f6b61`|1|1|1|1|1|1|1|1|1|DNA polymerase III sliding clamp (beta) subunit, PCNA homolog (DnaN) (PDB:1JQJ)|
|`func_060988dd6d3a`|1|1|1|1|1|1|1|1|1|Recombinational DNA repair ATPase RecF (RecF) (PDB:5Z67)|
|`func_46477a763d38`|1|1|1|1|1|1|1|1|1|Predicted  nucleic acid-binding protein, contains Zn-ribbon domain (includes truncated derivatives)|
|`func_184ca83fe067`|2|2|2|2|2|2|2|2|2|DNA gyrase/topoisomerase IV, subunit B (GyrB) (PDB:1EI1)|
|`func_2ee67347f04b`|2|2|2|2|2|2|2|2|2|DNA gyrase/topoisomerase IV, subunit A (GyrA) (PDB:1SUU)|
|`func_ff160914972f`|1|0|0|0|0|0|3|1|1|Predicted ATPase, archaeal AAA+ ATPase superfamily|
|`func_598325c18ceb`|1|0|0|0|0|0|1|1|1|Molybdopterin or thiamine biosynthesis adenylyltransferase (ThiF) (PDB:1ZUD) (PUBMED:32239579)|
|`func_8ce5ff84aa42`|14|13|13|10|10|10|22|17|16|Predicted arabinose efflux permease AraJ, MFS family (AraJ) (PDB:4LDS)|
|`func_300e2c6e37e4`|1|1|1|1|1|1|1|1|1|Glutamate dehydrogenase/leucine dehydrogenase (GdhA) (PDB:1B3B) (PUBMED:24391520)|
|`func_a666bc87a03f`|1|1|1|1|1|1|1|1|1|Large-conductance mechanosensitive channel (MscL) (PDB:2OAR)|
|`func_71379d89f0c6`|1|1|1|1|1|1|1|1|1|UTP pyrophosphatase, metal-dependent hydrolase family (YgjP) (PDB:4JIU) (PUBMED:27941785)|
|`func_772efbb7cb2e`|2|2|2|2|2|2|2|2|2|Hemolysin-related protein, contains CBS domains, UPF0053 family (TlyC) (PDB:2NQW)|
|`func_5a5bf9735c8a`|1|1|1|1|1|1|1|1|1|Carbonic anhydrase (CynT) (PDB:1EKJ) (PUBMED:22081392)|
|`func_1ffc26d0034f`|4|4|1|3|3|3|12|13|2|Transposase (or an inactivated derivative) (IS285)|
|`func_0631d9fecae6`|2|1|1|1|1|1|2|1|1|Phosphoenolpyruvate carboxylase (Ppc) (PDB:1FIY)|
|`func_62068e09f9f3`|2|1|1|0|0|0|1|1|1|Chromosome segregation ATPase Smc (Smc) (PDB:5XG3)|
|`func_f4b8779f6ad8`|2|0|0|0|0|0|0|0|0|Peptidoglycan-binding (PGRP) domain of peptidoglycan hydrolases (PGRP) (PDB:4FET)|
|`func_700652d0a2dd`|1|1|1|1|1|1|1|1|1|Uncharacterized membrane protein YjjP, DUF1212 family (YjjP)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Presence/absence of function hits across genomes

In contrast, here is the presence/absence report for the same data:

|**key**|**B_adolescentis_1_11**|**B_adolescentis_22L**|**B_adolescentis_6**|**B_lactis_Bl_04**|**B_lactis_CNCM_I_2494**|**B_lactis_DSM_10140**|**B_longum_JDM301**|**B_longum_KACC_91563**|**B_longum_NCIMB8809**|**COG20_FUNCTION**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|`func_3874a4219ffa`|1|1|1|1|1|1|1|1|1|Chromosomal replication initiation ATPase DnaA (DnaA) (PDB:1L8Q)|
|`func_1530552f6b61`|1|1|1|1|1|1|1|1|1|DNA polymerase III sliding clamp (beta) subunit, PCNA homolog (DnaN) (PDB:1JQJ)|
|`func_060988dd6d3a`|1|1|1|1|1|1|1|1|1|Recombinational DNA repair ATPase RecF (RecF) (PDB:5Z67)|
|`func_46477a763d38`|1|1|1|1|1|1|1|1|1|Predicted  nucleic acid-binding protein, contains Zn-ribbon domain (includes truncated derivatives)|
|`func_184ca83fe067`|1|1|1|1|1|1|1|1|1|DNA gyrase/topoisomerase IV, subunit B (GyrB) (PDB:1EI1)|
|`func_2ee67347f04b`|1|1|1|1|1|1|1|1|1|DNA gyrase/topoisomerase IV, subunit A (GyrA) (PDB:1SUU)|
|`func_ff160914972f`|1|0|0|0|0|0|1|1|1|Predicted ATPase, archaeal AAA+ ATPase superfamily|
|`func_598325c18ceb`|1|0|0|0|0|0|1|1|1|Molybdopterin or thiamine biosynthesis adenylyltransferase (ThiF) (PDB:1ZUD) (PUBMED:32239579)|
|`func_8ce5ff84aa42`|1|1|1|1|1|1|1|1|1|Predicted arabinose efflux permease AraJ, MFS family (AraJ) (PDB:4LDS)|
|`func_300e2c6e37e4`|1|1|1|1|1|1|1|1|1|Glutamate dehydrogenase/leucine dehydrogenase (GdhA) (PDB:1B3B) (PUBMED:24391520)|
|`func_a666bc87a03f`|1|1|1|1|1|1|1|1|1|Large-conductance mechanosensitive channel (MscL) (PDB:2OAR)|
|`func_71379d89f0c6`|1|1|1|1|1|1|1|1|1|UTP pyrophosphatase, metal-dependent hydrolase family (YgjP) (PDB:4JIU) (PUBMED:27941785)|
|`func_772efbb7cb2e`|1|1|1|1|1|1|1|1|1|Hemolysin-related protein, contains CBS domains, UPF0053 family (TlyC) (PDB:2NQW)|
|`func_5a5bf9735c8a`|1|1|1|1|1|1|1|1|1|Carbonic anhydrase (CynT) (PDB:1EKJ) (PUBMED:22081392)|
|`func_1ffc26d0034f`|1|1|1|1|1|1|1|1|1|Transposase (or an inactivated derivative) (IS285)|
|`func_0631d9fecae6`|1|1|1|1|1|1|1|1|1|Phosphoenolpyruvate carboxylase (Ppc) (PDB:1FIY)|
|`func_62068e09f9f3`|1|1|1|0|0|0|1|1|1|Chromosome segregation ATPase Smc (Smc) (PDB:5XG3)|
|`func_f4b8779f6ad8`|1|0|0|0|0|0|0|0|0|Peptidoglycan-binding (PGRP) domain of peptidoglycan hydrolases (PGRP) (PDB:4FET)|
|`func_700652d0a2dd`|1|1|1|1|1|1|1|1|1|Uncharacterized membrane protein YjjP, DUF1212 family (YjjP)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
