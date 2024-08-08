You can use this program to search for genes that encode specific functions in a %(contigs-db)s, %(genomes-storage-db)s or %(pan-db)s.

## Search

For example, you could search for all genes in a given %(contigs-db)s that encode a function related to the keyword 'kinase' the following way:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms kinase
{{ codestop }}

If you include `--verbose` flag in your command, %(anvi-search-functions)s will offer you a quick look at the matching function names:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms Phage \
                      --verbose
{{ codestop }}

```
Searching in Contigs Database ................: CONTIGS.db
Contigs DB ...................................: Initialized: CONTIGS.db (v. 23)
Search terms .................................: 1 found: 'Phage'
Case sensitive search? .......................: False
Exact match? .................................: False

Matches ......................................: 380 unique genes contained the search term "Phage"

Sneak peak into matching functions (up to 25)
===============================================
* Mobilome: prophages, transposons; Mobilome: prophages, transposons :: X; X within 'COG20_CATEGORY'
* Mobilome: prophages, transposons :: X within 'COG20_CATEGORY'
* Phage terminase-like protein, large subunit, contains N-terminal HTH domain (YmfN) :: COG4626 within 'COG20_FUNCTION'
* Phage-related protein (PblB) :: COG4926 within 'COG20_FUNCTION'
* Putative phage head morphogenesis protein, F of phage Mu or gp7 of SPP1 :: COG5585 within 'COG20_FUNCTION'
* Phage portal protein BeeE (BeeE) :: COG4695 within 'COG20_FUNCTION'
* Replication, recombination and repair; Mobilome: prophages, transposons :: L; X within 'COG20_CATEGORY'
* Abortive infection bacteriophage resistance protein (AbiF) :: COG4823 within 'COG20_FUNCTION'
* Phage-related holin (Lysis protein) :: COG4824 within 'COG20_FUNCTION'
* Phage terminase large subunit (XtmB) (PDB:4IDH) :: COG1783 within 'COG20_FUNCTION'
* Uncharacterized membrane protein YhgE, phage infection protein (PIP) family (YhgE) :: COG1511 within 'COG20_FUNCTION'
* Cell wall/membrane/envelope biogenesis; Mobilome: prophages, transposons :: M; X within 'COG20_CATEGORY'
* Prophage pi2 protein 07 :: COG4707 within 'COG20_FUNCTION'
* Prophage antirepressor :: COG3617 within 'COG20_FUNCTION'
* Phage antirepressor protein YoqD, KilAC domain (KilAC) :: COG3645 within 'COG20_FUNCTION'
* phage shock protein E :: K03972 within 'KOfam'
* Escherichia/Staphylococcus phage prohead protease :: K06904 within 'KOfam'
* putative DNA-invertase from lambdoid prophage Rac :: K14060 within 'KOfam'
* DNA primase, phage- or plasmid-associated :: COG3378 within 'COG20_FUNCTION'
* phage terminase large subunit :: K06909 within 'KOfam'
* Phage-related protein :: COG5412 within 'COG20_FUNCTION'
* phage terminase small subunit :: K07474 within 'KOfam'
* Phage head maturation protease :: COG3740 within 'COG20_FUNCTION'
* Phage-related replication protein YjqB, UPF0714/DUF867 family (YjqB) (PDB:3A9L) :: COG4195 within 'COG20_FUNCTION'
* phage shock protein C :: K03973 within 'KOfam'

Items additional data compatible output ......: search_results.txt
```

By default, the search is case insensitive. But you can make it case sensitive with the flag `--case-sensitive`:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms Phage \
                      --case-sensitive \
                      --verbose
{{ codestop }}

```
Searching in Contigs Database ................: CONTIGS.db
Contigs DB ...................................: Initialized: CONTIGS.db (v. 23)
Search terms .................................: 1 found: 'Phage'
Case sensitive search? .......................: True
Exact match? .................................: False

Matches ......................................: 109 unique genes contained the search term "Phage"

Sneak peak into matching functions (up to 25)
===============================================
* Phage terminase-like protein, large subunit, contains N-terminal HTH domain (YmfN) :: COG4626 within 'COG20_FUNCTION'
* Phage-related replication protein YjqB, UPF0714/DUF867 family (YjqB) (PDB:3A9L) :: COG4195 within 'COG20_FUNCTION'
* Phage antirepressor protein YoqD, KilAC domain (KilAC) :: COG3645 within 'COG20_FUNCTION'
* Phage anti-repressor protein Ant (PUBMED:27099293); Phage antirepressor protein YoqD, KilAC domain (KilAC) :: COG3561; COG3645 within 'COG20_FUNCTION'
* Phage terminase large subunit (XtmB) (PDB:4IDH) :: COG1783 within 'COG20_FUNCTION'
* Phage-related protein YomH (YomH) (PDB:2X8K) :: COG4722 within 'COG20_FUNCTION'
* Phage portal protein BeeE (BeeE) :: COG4695 within 'COG20_FUNCTION'
* Phage-related protein :: COG5412 within 'COG20_FUNCTION'
* Phage repressor protein C, contains Cro/C1-type HTH and peptisase s24 domains :: COG2932 within 'COG20_FUNCTION'
* Phage-related protein (PblB) :: COG4926 within 'COG20_FUNCTION'
* Phage shock protein PspC (stress-responsive transcriptional regulator) (PspC) :: COG1983 within 'COG20_FUNCTION'
* Murein DD-endopeptidase MepM and murein hydrolase activator NlpD, contains LysM domain (NlpD) (PDB:2GU1); Phage-related protein :: COG0739; COG5412 within 'COG20_FUNCTION'
* Phage-related holin (Lysis protein) :: COG4824 within 'COG20_FUNCTION'
* Phage head-tail adaptor (PDB:2Y3D) :: COG5614 within 'COG20_FUNCTION'
* Phage-related minor tail protein (YqbO); Phage-related protein :: COG5280; COG5412 within 'COG20_FUNCTION'
* Phage-related tail protein :: COG5283 within 'COG20_FUNCTION'
* Phage regulatory protein Rha (pRha) :: COG3646 within 'COG20_FUNCTION'
* Phage head maturation protease :: COG3740 within 'COG20_FUNCTION'
* Phage terminase, small subunit (XtmA) (PDB:4ZC3) :: COG3728 within 'COG20_FUNCTION'
* Prophage antirepressor; Phage antirepressor protein YoqD, KilAC domain (KilAC) :: COG3617; COG3645 within 'COG20_FUNCTION'
* Phage terminase, small subunit :: COG3747 within 'COG20_FUNCTION'
* Murein DD-endopeptidase MepM and murein hydrolase activator NlpD, contains LysM domain (NlpD) (PDB:2GU1); SLT domain protein (SLT) (PUBMED:15608122); Phage-related minor tail protein (YqbO); Phage-related protein :: COG0739; COG3953; COG5280; COG5412 within 'COG20_FUNCTION'
* Phage shock protein A (PspA) (PDB:4WHE) :: COG1842 within 'COG20_FUNCTION'
* Transcriptional regulator, contains XRE-family HTH domain (HipB) (PDB:1ADR); Phage repressor protein C, contains Cro/C1-type HTH and peptisase s24 domains :: COG1396; COG2932 within 'COG20_FUNCTION'

Items additional data compatible output ......: search_results.txt
```

By defult the search term is searched anywhere in the function. The `--exact-match` flag change that behavior, and will only return matches if the entire field of text matches to the search term:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms "Phage head maturation protease" \
                      --exact-match \
                      --verbose
{{ codestop }}

```
Searching in Contigs Database ................: CONTIGS.db
Contigs DB ...................................: Initialized: CONTIGS.db (v. 23)
Search terms .................................: 1 found: 'Phage head maturation protease'
Case sensitive search? .......................: False
Exact match? .................................: True

Matches ......................................: 3 unique genes contained the search term "Phage head maturation protease"

Sneak peak into matching functions (up to 25)
===============================================
* Phage head maturation protease :: COG3740 within 'COG20_FUNCTION'

Items additional data compatible output ......: search_results.txt
```

By default, all function annotation sources are included in the search, however, the parameter `--annotation-sources` can be used to limit the search to one or more specific annotation sources:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms Phage \
                      --annotation-sources KOfam \
                      --verbose
{{ codestop }}

```
Searching in Contigs Database ................: CONTIGS.db
Contigs DB ...................................: Initialized: CONTIGS.db (v. 23)
Search terms .................................: 1 found: 'Phage'
Case sensitive search? .......................: False
Exact match? .................................: False

Matches ......................................: 51 unique genes contained the search term "Phage"

Sneak peak into matching functions (up to 25)
===============================================
* putative DNA-invertase from lambdoid prophage Rac :: K14060 within 'KOfam'
* phage terminase small subunit :: K07474 within 'KOfam'
* phage shock protein E :: K03972 within 'KOfam'
* DNA polymerase bacteriophage-type [EC:2.7.7.7] :: K02334 within 'KOfam'
* phage terminase large subunit :: K06909 within 'KOfam'
* macrophage erythroblast attacher :: K18624 within 'KOfam'
* Escherichia/Staphylococcus phage prohead protease :: K06904 within 'KOfam'
* phage shock protein A :: K03969 within 'KOfam'
* phage shock protein C :: K03973 within 'KOfam'

Items additional data compatible output ......: search_results.txt
```

Similar to functions, accession IDs can be searched, as well:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms COG4824,COG5614 \
                      --verbose
{{ codestop }}

```
Searching in Contigs Database ................: CONTIGS.db
Contigs DB ...................................: Initialized: CONTIGS.db (v. 23)
Search terms .................................: 2 found: 'COG4824|COG5614'
Case sensitive search? .......................: False
Exact match? .................................: False

Matches ......................................: 8 unique genes contained the search term "COG4824"

Sneak peak into matching functions (up to 25)
===============================================
* Phage-related holin (Lysis protein) :: COG4824 within 'COG20_FUNCTION'

Matches ......................................: 1 unique genes contained the search term "COG5614"

Sneak peak into matching functions (up to 25)
===============================================
* Phage head-tail adaptor (PDB:2Y3D) :: COG5614 within 'COG20_FUNCTION'

Items additional data compatible output ......: search_results.txt
```

## Output

By default, the output will be a fairly barren, and will only show which contigs contain genes that matched to a given. This output file will be most helpful as an additional layer in the anvi'o interactive interface to quickly see the items that include genes with functions that match to the search terms.

However, generating a more comprehensive report is also an option through the parameter `--full-report`:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms Phage \
                      --full-report phage-results.txt \
                      --verbose
{{ codestop }}

```
anvi-search-functions -c CONTIGS.db --search-terms Phage --verbose --full-report phage-reuslts.txt
Searching in Contigs Database ................: CONTIGS.db
Contigs DB ...................................: Initialized: CONTIGS.db (v. 23)
Search terms .................................: 1 found: 'Phage'
Case sensitive search? .......................: False
Exact match? .................................: False

Matches ......................................: 380 unique genes contained the search term "Phage"

Sneak peak into matching functions (up to 25)
===============================================
* Uncharacterized membrane protein YhgE, phage infection protein (PIP) family (YhgE) :: COG1511 within 'COG20_FUNCTION'
* Phage-related holin (Lysis protein) :: COG4824 within 'COG20_FUNCTION'
* Mobilome: prophages, transposons :: X within 'COG20_CATEGORY'
* DNA primase, phage- or plasmid-associated :: COG3378 within 'COG20_FUNCTION'
* Phage-related protein (PblB) :: COG4926 within 'COG20_FUNCTION'
* Phage-related replication protein YjqB, UPF0714/DUF867 family (YjqB) (PDB:3A9L) :: COG4195 within 'COG20_FUNCTION'
* Uncharacterized phage-associated protein, contains DUF4065 domain (GepA) :: COG3600 within 'COG20_FUNCTION'
* phage terminase large subunit :: K06909 within 'KOfam'
* Predicted phage phi-C31 gp36 major capsid-like protein (gp36) :: COG4653 within 'COG20_FUNCTION'
* phage terminase small subunit :: K07474 within 'KOfam'
* Phage anti-repressor protein Ant (PUBMED:27099293); Phage antirepressor protein YoqD, KilAC domain (KilAC) :: COG3561; COG3645 within 'COG20_FUNCTION'
* Phage-related minor tail protein (YqbO); Phage-related protein :: COG5280; COG5412 within 'COG20_FUNCTION'
* Abortive infection bacteriophage resistance protein (AbiF) :: COG4823 within 'COG20_FUNCTION'
* Phage-related protein :: COG5412 within 'COG20_FUNCTION'
* Phage shock protein PspC (stress-responsive transcriptional regulator) (PspC) :: COG1983 within 'COG20_FUNCTION'
* Phage-related protein YomH (YomH) (PDB:2X8K) :: COG4722 within 'COG20_FUNCTION'
* Holliday junction resolvase RusA (prophage-encoded endonuclease) (RusA) (PDB:1Q8R) :: COG4570 within 'COG20_FUNCTION'
* Prophage antirepressor; Phage antirepressor protein YoqD, KilAC domain (KilAC) :: COG3617; COG3645 within 'COG20_FUNCTION'
* phage shock protein E :: K03972 within 'KOfam'
* Cell wall/membrane/envelope biogenesis; Mobilome: prophages, transposons; Mobilome: prophages, transposons; Mobilome: prophages, transposons :: M; X; X; X within 'COG20_CATEGORY'
* Phage terminase large subunit (XtmB) (PDB:4IDH) :: COG1783 within 'COG20_FUNCTION'
* DNA polymerase bacteriophage-type [EC:2.7.7.7] :: K02334 within 'KOfam'
* Mobilome: prophages, transposons; Mobilome: prophages, transposons :: X; X within 'COG20_CATEGORY'
* putative DNA-invertase from lambdoid prophage Rac :: K14060 within 'KOfam'
* phage shock protein C :: K03973 within 'KOfam'

Items additional data compatible output ......: search_results.txt
Full report ..................................: phage-reuslts.txt
```

This will result in an output file that will look like this:

|**gene_callers_id**|**source**|**accession**|**function**|**search_term**|**contigs**|
|:--|:--|:--|:--|:--|:--|
|560|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig1_split_00030|
|563|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig1_split_00030|
|756|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig1_split_00041|
|757|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig1_split_00041|
|779|COG20_FUNCTION|COG1983|Phage shock protein PspC (stress-responsive transcriptional regulator) (PspC)|Phage|Day17a_QCcontig1_split_00042|
|1550|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig4_split_00008|
|1555|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig4_split_00009|
|1556|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig4_split_00009|
|1557|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig4_split_00009|
|1604|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig4_split_00011|
|1605|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig4_split_00011|
|2580|COG20_FUNCTION|COG4653|Predicted phage phi-C31 gp36 major capsid-like protein (gp36)|Phage|Day17a_QCcontig9_split_00011|
|2580|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig9_split_00011|
|2733|COG20_FUNCTION|COG4824|Phage-related holin (Lysis protein)|Phage|Day17a_QCcontig13_split_00002|
|2733|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig13_split_00002|
|2767|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig14_split_00001|
|2825|COG20_FUNCTION|COG1511|Uncharacterized membrane protein YhgE, phage infection protein (PIP) family (YhgE)|Phage|Day17a_QCcontig15_split_00003|
|2965|COG20_FUNCTION|COG4823|Abortive infection bacteriophage resistance protein (AbiF)|Phage|Day17a_QCcontig16_split_00006|
|3029|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig16_split_00009|
|3170|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig17_split_00002|
|13008|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig276_split_00001|
|13009|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig277_split_00001|
|13022|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig279_split_00001|
|13029|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig280_split_00001|
|13106|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig283_split_00001|
|13107|COG20_FUNCTION|COG3378|DNA primase, phage- or plasmid-associated|Phage|Day17a_QCcontig283_split_00001|
|13107|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig283_split_00001|
|13116|COG20_FUNCTION|COG3728|Phage terminase, small subunit (XtmA) (PDB:4ZC3)|Phage|Day17a_QCcontig283_split_00001|
|13116|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig283_split_00001|
|13135|COG20_FUNCTION|COG2369|Uncharacterized protein, contains phage Mu head morphogenesis gpF-like domain|Phage|Day17a_QCcontig285_split_00001|
|13135|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig285_split_00001|
|13141|COG20_FUNCTION|COG2369|Uncharacterized protein, contains phage Mu head morphogenesis gpF-like domain|Phage|Day17a_QCcontig286_split_00001|
|13141|COG20_CATEGORY|X|Mobilome: prophages, transposons|Phage|Day17a_QCcontig286_split_00001|
|15448|KOfam|K06909|phage terminase large subunit|Phage|Day17a_QCcontig476_split_00003|
|13143|KOfam|K06909|phage terminase large subunit|Phage|Day17a_QCcontig286_split_00001|
|15293|KOfam|K06909|phage terminase large subunit|Phage|Day17a_QCcontig470_split_00001|
|15450|KOfam|K06909|phage terminase large subunit|Phage|Day17a_QCcontig476_split_00003|
|8686|KOfam|K03973|phage shock protein C|Phage|Day17a_QCcontig100_split_00001|
|9711|KOfam|K03973|phage shock protein C|Phage|Day17a_QCcontig125_split_00007|
|14379|KOfam|K03972|phage shock protein E|Phage|Day17a_QCcontig379_split_00001|
|27227|KOfam|K18624|macrophage erythroblast attacher|Phage|Day17a_QCcontig2462_split_00001|
|31918|KOfam|K06904|Escherichia/Staphylococcus phage prohead protease|Phage|Day17a_QCcontig4102_split_00001|
|31919|KOfam|K06904|Escherichia/Staphylococcus phage prohead protease|Phage|Day17a_QCcontig4102_split_00001|
|2561|KOfam|K02334|DNA polymerase bacteriophage-type [EC:2.7.7.7]|Phage|Day17a_QCcontig9_split_00011|
|1803|KOfam|K02334|DNA polymerase bacteriophage-type [EC:2.7.7.7]|Phage|Day17a_QCcontig4_split_00021|
|3955|KOfam|K02334|DNA polymerase bacteriophage-type [EC:2.7.7.7]|Phage|Day17a_QCcontig26_split_00010|
|7105|KOfam|K03969|phage shock protein A|Phage|Day17a_QCcontig74_split_00005|
|4518|KOfam|K07474|phage terminase small subunit|Phage|Day17a_QCcontig33_split_00009|
|4315|KOfam|K14060|putative DNA-invertase from lambdoid prophage Rac|Phage|Day17a_QCcontig30_split_00002|
|4299|KOfam|K14060|putative DNA-invertase from lambdoid prophage Rac|Phage|Day17a_QCcontig30_split_00001|
|3770|KOfam|K14060|putative DNA-invertase from lambdoid prophage Rac|Phage|Day17a_QCcontig25_split_00001|
|4519|KOfam|K06909|phage terminase large subunit|Phage|Day17a_QCcontig33_split_00009|

It is also possible to include sequences of genes into this output:

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms Phage \
                      --full-report phage-results.txt \
                      --include-sequences \
                      --verbose
{{ codestop }}
