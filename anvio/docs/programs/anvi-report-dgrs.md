## Usage
This program locates diversity-generating retroelements (DGRs for short - to save words and typing) within read mapped metagenomic or genomic samples.

In brief, the tool first searches for the variable regions of the DGR, by locating areas of high SNV density which are deemed as possible variable regions. These SNVs should by default occur in the 1st and 2nd codon base positions so that they are within an ORF. Next we search for the template regions associated with those possible variable regions, by BLASTn'ing the genome against itself to locate regions of similarity other than the disparity in the Adenine bases. Finally, the tool then locates a reverse transcriptase, this works by locating the nearest reverse transcriptase to each template region via a homological HMM search.

### First thing's first what are diversity-generating retroelements?
Diversity-generating retroelements (DGRs) are molecular evolutionary mechanisms which facilitate rapid microbial adaptation to environmental changes. DGRs are capable of changing the protein expression of open reading frames in which they act in, by generating site-specific nucleotide variability in the hypermutated target genes via non-synonymous changes.

DGRs are widespread across bacteria, archaea, and many viruses, and have been identified in more than 90 environmental contexts. This widespread mechanisms are composed of three essential elements: the **template region** (TR), **variable region** (VR), and an error-prone **reverse transcriptase** (RT). The TR and VR are nearly identical sequences of the same length, differing only at positions subject to mutagenesis. This controlled divergence generates a repertoire of protein variants, enhancing an organism’s adaptability and functional potential.

 The template region, which are often encoded in non-coding DNA regions, are used as the template for error-prone reverse transcription. Reverse transcription produces cDNA, which is then inserted into an open reading frame, with the adenine site mutagenesis, known now as the variable region. A single template region can be associated with multiple variable regions, allowing  multiple proteins of similar functions to diversify as a result of one system. Template and variable regions can be distally located to one another, and their the specifics of the mechanisms structure can also vary between organisms.

### Our philosophy for finding DGRs
Locating the template and variable region first is an integral philosophy for locating DGRs without relying on reverse transcriptase homology. In fact, the reverse transcriptase is not a requirement for the reporting of a DGR.

Another integral part of this tool, is that it relies on SNVs in the variable region and therefore, only locates DGRs that are active within the sample, much like a [weeping angel](https://en.wikipedia.org/wiki/Weeping_Angel) cannot attack what doesn't move, anvi'o cannot locate inactive DGRs.

### Pre-requisites for running this programme

Before you run this program you need to run locate the reverse transcriptases for the DGRs in your sample/s. You can do this via the anvi'o own reverse transcriptase collection which is compsed of the [Roux et al., 2021](https://doi.org/10.1038/s41467-021-23402-7), alongside 4 Pfam general RT HMMs: PF00078/RVT_1 Reverse transcriptase (RNA-dependent DNA polymerase), PF07727/RVT_2 reversevtranscriptase (RNA-dependent DNA polymerase), PF13456/RVT_3 reverse transcriptase-like – found in plants, PF13655/RVT_N N-terminal domain of reverse transcriptase (do i cite simon and zimmerly 2008 here??).

{{ codestart }}
anvi-run-hmms -c contigs.db  -I Reverse_Transcriptase
{{ codestop }}

You can of course create your own HMM and use that instead, who is anvi'o to tell you what to do?

### Essential inputs for searching for DGRs

If you want to see more diversity then your %(single-profile-db) needs to be a %(merged-profile-db) so you can look at activity across samples

* %(contigs-db)
* ideally a merged %(profile.db) - we don't discriminate against singles here either
* %(hmm-usage)

Here is an example of the standard command you could use with default parameters:

{{ codestart }}

anvi-report-dgrs -c contigs.db -p profile.db -I Reverse_Transcriptase -o dgrs-output --skip-compute-DGR-variability-profiling

{{ codestop }}

You can also run the program and compute the activity of the DGR in your samples (explained under the header: **Computing DGR variability profiling**) for this you will need the samples that your want to look for the activity across and for this anvi'o needs to know the location of these samples via a %{samples-txt}.

{{ codestart }}

anvi-report-dgrs -c contigs.db -p profile.db -I Reverse_Transcriptase -o dgrs-output --samples-txt samples.txt

{{ codestop }}

### Identifying putative variable regions

The tool bases the location of the variable regions on the principle that at the mutagenesis positions of the variable region there will be SNVs, these clusters of dense SNVs are located. The identification of these SNVs is tunable through multiple parameters. `--departure-from-reference-percentage`is the minimum departure from reference for the tool to count a SNV. The tool locates a SNV and then walks along to the next SNV the '--distance-between-snv' controls this length. These are then passed through a `--minimum-snv-density` threshold and a `--minimum-range-size`. Then any of these overlapping SNV clusters are merged together and a `--variable-buffer-length` is added to each side of a high-SNV region.

### Locating template regions

Once we have located the putative variable region sequences, we need to confirm that there is a template region in the sample, or else it is just an area with a lot of SNVs. For this we use BLASTn to BLAST the putative variable regions back across the metagenomic sample to find areas which match the characteristics of a template and variable region pair.

Then we check that these are truly a template and variable region pair. For that they have to have a `--number-of-mismatches` and those mismatches have to mainly be from one type of base (previously DGRs have only been reports to be adenine based specific for this use the flag `--only-a-bases`). The flag `--percentage-mismatch` looks at the percentage of mismatches that comes from one base type in the template region.

### Confirming DGRs

Due to the very complex nature of metagenomic samples anvi'o has a number of different ways to ensure that we are truly finding a DGR pair and not just two random sequences in the sample that match at almost every position. Having said that it is always a good idea to double check the results yourself.

The sequences that are found can be rather short in length and heavily repeated. In order to avoid heavily repeated sequences being found as false positives we use the python compiled version of [tantan](https://doi.org/10.1093/nar/gkq1212) to remove sequences that are over a `--repeat-threshold` of the fraction of that sequence that is repeated. We also check that there are multiple base types on both of the regions to ensure that they are diverse this can be changed with `--min-base-types-vr` and `--min-base-types-tr`.

At this point, the algorithm has no concept of the strand in which the regions were found on, therefore, if the mutagenesis base is thymine, this is reverse complemented to be adenine and saves in the output information table.

`anvi-report-dgrs` only searches for active DGRs, the search is based on SNVs, as we expect these to be located in the mismatching positions of the TR and VR which are presumed to correspond to the mutagenesis base. Due to this logic, we check for two main SNV related things. The first being that the SNVs have to come from the majority of the first and second codon positions as such random logic that less than a third of the SNVs would be from the 3rd codon position, variable with the flag `----snv-codon-position`. This is therefore reported in the output as `snv_at_3_codon_over_a_third`.

The second SNV related check performed is a result of expecting the places of the mismatches in the VR to be the sites of diversification, we therefore want to check that the majority of the SNVs are also at those positions, meaning that it is unexpected for SNVs to be over matching positions other than the site of mutagenesis. Of course other reasons for DGR spurious SNVs to occur over a VR due to populations variances and such. Therefore, the tool checks for SNVs in matches that are not over the mutagenesis base type matches, so in the case of adenine it would search for SNVs in matches of the VR that occur over the T, G, C bases. This is also then reported, so that the user can investigate more closely and decide on thresholds themselves using the flag `----snv-matching-proportion`.

## Locating associated Reverse Transcriptases

Next, the tool then locates a reverse transcriptase, this works by locating the nearest reverse transcriptase to each template region. these reverse transcriptases are annotated via a HMM annotation when you run `anvi-run-hmms` on your %{contigs-db}.

{{ codestart }}
anvi-run-hmms -c contigs.db  -i Reverse_Transcriptase
{{ codestop }}

Because the tool relies on the **reverse transcriptase** being on the same contig as a **template region**, these are not always found or associated with a found TR/VR pair. Therefore, a reverse transciptase is not a critical factor for the determinance of a DGR.

The user can of course run another HMM that has reverse transcriptases by creating a new HMM collection via a [user-defined HMM source](https://anvio.org/help/main/artifacts/hmm-source/), and then running this with your %(contigs-db) and then `anvi-report-dgrs`.

### Reporting genomic context around DGRs

This section of the code provides genomic context surrounding each variable region, this is to help you understand the context of the context of the variable region. The number of genes up and downstream of the variable region can be altered with the parameter `--num-genes-to-consider-in-context`. This information is then printed in the output directory as tsv file called DGR-OUTPUT_DGR_genes_found.tsv.

By default because anvi'o uses [Prodigal](https://github.com/hyattpd/Prodigal) for gene prediction curing the generation of your %{contigs-db}, unless you have stated otherwise. Therefore, for this to work, you need to have the same gene prediction tool as your contigs, you can change this using `--gene_caller_to_consider_in_context`. NB: "!!!" is a splitter for there being more than one function or accession associated with a gene

| DGR_ID  | VR_ID  | Contig                                         | Start | Stop  | Direction | Partial | Call_Type | Gene_Caller_Source | Version | Gene_Caller_ID | DNA_Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | AA_Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Length | Gene_Functions                                                                                                                                                                                                                                                 | Gene_Function_Source         | Gene_Function_Accession |
| ------- | ------ | ---------------------------------------------- | ----- | ----- | --------- | ------- | --------- | ------------------ | ------- | -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------- | ----------------------- |
| DGR_001 | VR_001 | B_fragilis_ARW016_000000000001_VR_001_rev_comp | 1718  | 2303  | r         | 0       | 1         | prodigal           | v2.6.3  | 4              | ATGTCAGAACTTAAAATTACACAAGAAAAGGTTACAGCCGCTTTTAGTGAAGCAAACGACTGTCCTAAAGCAATTAGTATTCTAACAGCTTTATTCGGAAAGCAAAAGCCGGATTATACAGATTATCACAATATCAAAACCTACGAAGATGCTTGTGAAGCAATAGGTGTAAAACCTATTGTTCGCCTACTTGTTGAAGATGAAGACGGACACAAAGAAGAAGTGGCTGATATTGCACACCTCGCCTACATCAAACTATGCACAATTGCTCGTGCATTGAACAACGATCCTAATTTTCCACGATTTACTAAAGATGAATACCGTTATACGCCGTGGTTTTATCTTTATAATCAGAAAGAAATTGATGAAATGGACGAAGAGGATCGTAATCGGCTGGTTCTTTGGGGCGGCAATGCGAGTGCCGGTGCGTATTGCGGCCTCGCTTATGCGAACTCGACTCACGCTTGGTCGTACTCGACTGCGAAGCTCGGCTCTCGCCTTGCTGTAAAATCAAGTGAAATCGCAATTTACTTTGGAGAACAATTCAAAGAATTGTGGAAAGACTTTCTGATTGGAAAAAAGTAA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | MSELKITQEKVTAAFSEANDCPKAISILTALFGKQKPDYTDYHNIKTYEDACEAIGVKPIVRLLVEDEDGHKEEVADIAHLAYIKLCTIARALNNDPNFPRFTKDEYRYTPWFYLYNQKEIDEMDEEDRNRLVLWGGNASAGAYCGLAYANSTHAWSYSTAKLGSRLAVKSSEIAIYFGEQFKELWKDFLIGKK                                                                                                                                                                                                                                                                                                                                                                                                                | 585    |                                                                                                                                                                                                                                                                |                              |                         |        |
| DGR_002 | VR_001 | SAMEA2619974_000000008595_DGR_no_RT            | 15732 | 16398 | r         | 0       | 1         | prodigal           | v2.6.3  | 70             | ATGGGCAAAAGCCTAAGACGGTTGGAACAGAAGAAGTTGCAGGTAGGCATCGGCAAACTAGCCGATGCGGGCGTGGATGTGGAACGGTGGCTCGCCATGTGTGACGCGGACGATGCAGCGATGCAGCGCCTCGTCGAGGCATGGCCAGTCAAACTCCCATCCTTCGTCTACGATGCCAAGACCGTTTGCAGCATCTTGGGTCTCCCGAACAACTGCAAGGAGGAGACACCGAAGGCTGCTTACGGAGAGGTCGTGGTTTGGTACGGTGGTTGGACACTCGGGGAGTTAGTGGCGACCGACAAGGTCGTCAACTATCTCTCCAAAGAGCGAGAGTCGTGGAAAGCCCCTCCGGGTTATTACCACACTCGCATTCCTGTCCCAGACAGCAACCGGATGACGCACGACGAGCAGGTCGCTGCTCACTTGGCGCAGCTCTACGCTGCATTCAAGGAGCTGCCCACACCCGTTGGCGCGACTGTGCTTGCTGCGCATCTCGACGCTACAGGCGAAGATCTCCTCAACGGGAACTTCTGCCGTTGCGCCGAACCGCTCCCCGACGACCGCCATGCCGGGCTCTTTGTCAGCGAGGGTCGCGTGCGCGTCAACTACTACTGGGATGGCGATCGCCGCGTCGGCGTCTTCCTTGGTGCCGCTCGGAAGGCCTGA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | MGKSLRRLEQKKLQVGIGKLADAGVDVERWLAMCDADDAAMQRLVEAWPVKLPSFVYDAKTVCSILGLPNNCKEETPKAAYGEVVVWYGGWTLGELVATDKVVNYLSKERESWKAPPGYYHTRIPVPDSNRMTHDEQVAAHLAQLYAAFKELPTPVGATVLAAHLDATGEDLLNGNFCRCAEPLPDDRHAGLFVSEGRVRVNYYWDGDRRVGVFLGAARKA                                                                                                                                                                                                                                                                                                                                                                                     | 666    |                                                                                                                                                                                                                                                                |                              |                         |        |
| DGR_003 | VR_001 | T_erythraeum_IMS101_000000000001_VR_001        | 9496  | 11278 | f         | 0       | 1         | prodigal           | v2.6.3  | 27             | ATGCAGTCTGAACAAATTCTGAGAAACCGTTATAAAGTCATTAAATCTCTAGGAACTGGAGGGTTTGGCTATACCTACTTAGCTGAAGATCTAGACTTACCCGGATATCCCAAATGCGTGGTCAAACATCTGAAGCCAAAAAGCCCGGACTCTACAGTTCTAAATATTGCTAGAAAGTTATTTCTAAGAGAAGCAGATATCTTATACAAATTAGGTAATGACTCAGACCAAATACCAAGATTATTTGCCTATTTCCAAGAACAAAGAGAATTTTATTTAGTACAAGAATATATAGAAGGGCAAGATATAAGTCGGGAGCTGACTCCGGGAAAAAAATTAAGTGAGTCTGACACTATTGCTTTACTCAAAGGAATATTAGAAGCATTAACAGTTGCTCATCAAAATAATGTAATTCACCGAGATATTAAACCCCAAAATTTGATGCGCCGTAGGTCAGACAATAAAATAGTATTAATAGACTTTGGGGCGGTAAAAGAAATAGATGTATTGACTATAAACCAACAAGGTGCAACAACTTTAACTGTTGCAGTGGGAACTCCTGGCTATATGCCCAGCGAACAAAGCAATGGTAAACCGAAGCTCAGTAGCGACATATATGCAGTAGGAATGGTAGGGATAAGAGCATTGACAGGTAAGGAACCCCAGAGTTTACTCACAGACCCCAAAACAGGAAATGTGATCTGGCGAAATGAAGCTCAAGTAAGCAACCGTCTGGCAGATATACTAGATAAGATGGTGCATGAATATTTTCCCCAGCGTTATGAAAATGCCATGGAAGTCTTAGATGTTTTACAAGAAGCGAGGAAAAAACCTAACCCACCTAAGCCAACACAGAGACCAGCGGCAGAAAGCACACTCAACGTCAACCCATCCCAACCCCGGGGTCAAGGAGGAAGTTGGTTTGCGGGGAGATTTTTCAATCAAAAACCAGCAAAAATTTCTATTCAGACATTTACTACTGTCACAGTTAACCGAAGAGGAGAAATAATAAGTCGCGCTCAAGGTCAAGCAGAAGTAATAACAGAAAATATTGGTAATGGAGTTTCTCTAGAAATGGTAAAAATTCCCGGAGGTAGGTTTTTGATGGGGTCTCCAAAGAAAGAAGCAGAAAGACTTGACAAGGAAGGTCCGCAACATTATGTAGATGTGCCAGAATTTTTGATGGGAAAATATGCAGTTACTCAAGCACAGTGGGAAGCAGTTATGGGAAATAACCCTGCTAGGTTTAAAGGTGCAAACCGTCCTGTGGAAAAAGTAAGTTGGAATGATGCGACAGAATTTTGTCGGAAACTCTCACGAATAACAGGAAAACAATATAGTCTACCCAGTGAGAGCCAATGGGAATATGCTTGTCGAGCCCGAACAACAACACCATTTTATTTTGGAGAAACTATAACGCCTGAGTTAGTTAACTACGATGGCAACTATACTTATGCTGATGCTCCGAAAGGAAAATATCGAAAAGAAACAACAGATGTGGGAATTTTTCCACCGAATGCATTTGGTTTGTACGATATGCATGGGAATGTGTACGAGTTTTGTCAGGATGTTTGGCATGAAAACTATAATGGAGCACCTACTGATGGCAGTGCCTGGGAAACTGGAGGTGATAGTAGTAGAAGAGTCTGTCGTGGCGGCTCCTGGGTCAACTATCCTGGGAGGTGTCGCTCTGCGGACCGCATCAACTATGACTCGGTCGGGGCGGACTACATCGGTATTGGTTTTCGTCTTGTGAGTTTCCCCCCCAGGACTTTTGAATAG | MQSEQILRNRYKVIKSLGTGGFGYTYLAEDLDLPGYPKCVVKHLKPKSPDSTVLNIARKLFLREADILYKLGNDSDQIPRLFAYFQEQREFYLVQEYIEGQDISRELTPGKKLSESDTIALLKGILEALTVAHQNNVIHRDIKPQNLMRRRSDNKIVLIDFGAVKEIDVLTINQQGATTLTVAVGTPGYMPSEQSNGKPKLSSDIYAVGMVGIRALTGKEPQSLLTDPKTGNVIWRNEAQVSNRLADILDKMVHEYFPQRYENAMEVLDVLQEARKKPNPPKPTQRPAAESTLNVNPSQPRGQGGSWFAGRFFNQKPAKISIQTFTTVTVNRRGEIISRAQGQAEVITENIGNGVSLEMVKIPGGRFLMGSPKKEAERLDKEGPQHYVDVPEFLMGKYAVTQAQWEAVMGNNPARFKGANRPVEKVSWNDATEFCRKLSRITGKQYSLPSESQWEYACRARTTTPFYFGETITPELVNYDGNYTYADAPKGKYRKETTDVGIFPPNAFGLYDMHGNVYEFCQDVWHENYNGAPTDGSAWETGGDSSRRVCRGGSWVNYPGRCRSADRINYDSVGADYIGIGFRLVSFPPRTFE | 1782   | Serine/threonine protein kinase (SPS1) (PDB:6G4J)!!!Formylglycine-generating enzyme, required for sulfatase activity, contains SUMF1/FGE domain (YfmG) (PDB:2AFT)Signal transduction mechanisms!!!Posttranslational modification, protein turnover, chaperones | COG20_FUNCTIONCOG20_CATEGORY | COG0515!!!COG1262T!!!O  |
| DGR_003 | VR_002 | T_erythraeum_IMS101_000000000001_VR_002        | 6606  | 8190  | f         | 0       | 1         | prodigal           | v2.6.3  | 42             | ATGACAGAAATAACAGATAAAATTATCGCTCAGAATAAAGTAGAACAATTTGTAAGTCGTTTTGAAGATTCCTATTTCTTATTAGCTTCCCATGTCGCTTTACCTCTAGTATTAACTCCAGAATTAGTCAATTATCTCAGAATTGAATTTCTCAAAAGTGAAAAAGTTCCTTGGATTTCTGAAGCTGATTTATTGTTATCTAATTTGTGTCGTCCAGTTGGTTATGAACTTTATGTAATGGATGAAGGAGTAAGAATTTATTTACTCCAAAAATTAGCAGAAAACCCAAAGTTTGGTCAGAAAAGAATTAAAGAAATTGCCCGTTTATTACTGAGCTATATTAAATATTTAGCTAAAACTAATCATTTTATTGGAGACAAAGAACTACAACTTCAGACATGGGGTGCAATGTTGTATTTAGATACAGAGAAAACTGTCAAGGAAATAGCTTTTGCTATCTCTCAATGTATAGATAAAGCAGAGTTAGCTCGTTTACTAAAAATTACCCAAGATTTTAAACAACAAATAGCTTCATCAGGACAATTTCAAGACTTTTTAGACTATGCTCAACTTTGTAATGAATTACTTAGAGAACCTGAAAAAGTAACCCTAGAAAAAATAACACATTCATACCAGGTCGCAAAGCTTGAATTAACTATACCAAAAATAGTTTTATCAAAAGAATTAGATCTAGGGGAATTATTTTATTTTGAAGTTGTAAAAGTTGATAATTTTGGCAATATTATTGAGCGTAGTCAGGAAAGCGCTAGGCAGAAAATAGAATATTTGGGTAATGAAATTAAATTGGAAATGGTTTATATTCCTGGGGGTACTTTTACCATGGGTTCTCCTAAAAGTGAAGTGGATAGCAGAGATAATGAACGGCCCCAACATCAGGTAACTGTTCCTTCCTTTTTTATGGGCAAATATCCAGTAACTCAAGGACAGTGGAAAGCGATCGCCTCTCGAACAGAATTAAAAGTAGAACGAGACCTAGACCCAGAGCCATCATATTTTAAAGAACCATACGAAAGTATAGATAGATGGCAGAGACCAGTGGAGCTGGTTAGTTGGAACCAAGCAGTAGAATTTTGTCAAAGACTATCAAAACTAACAGGAAGAAAGTACAGGTTACCCAGTGAAGCCGAATGGGAATATGCTTGTCGTGCAGGAACGAGTACACCATTCTACTTTGGAGAAACCATAACATCTGAGTTAGCTAACTATGATGGCAACCATACTTACGGCAAAGGGCCAAAAGGAGAAAATAGAGAGCAAACAACTCCTGTAGGTAAATTTCCACCAAATGCTTTTGGTTTATATGATATGCACGGAAACGTTGATGAATGGTGTGCTGACGACTGGTATGATAATTATGTAGGAGCTCCTACAGATGGAAGTGCTAGGATTAATAGTAATAAAAATAGTAATACTAATAGTACAAAAGTATTTCGTGGCGGCTCCTGGATCTACTATCCTTGGTGGTGTCGCTCTGCGCTCCGCGGCAACTTTAACTCGGACGAGGCGGACAGCATCTCTTTTGGTTTTCGTCTTGTGAGTTTCCCCCCCAGGACTTTTGAATAG                                                                                                                                                                                                       | MTEITDKIIAQNKVEQFVSRFEDSYFLLASHVALPLVLTPELVNYLRIEFLKSEKVPWISEADLLLSNLCRPVGYELYVMDEGVRIYLLQKLAENPKFGQKRIKEIARLLLSYIKYLAKTNHFIGDKELQLQTWGAMLYLDTEKTVKEIAFAISQCIDKAELARLLKITQDFKQQIASSGQFQDFLDYAQLCNELLREPEKVTLEKITHSYQVAKLELTIPKIVLSKELDLGELFYFEVVKVDNFGNIIERSQESARQKIEYLGNEIKLEMVYIPGGTFTMGSPKSEVDSRDNERPQHQVTVPSFFMGKYPVTQGQWKAIASRTELKVERDLDPEPSYFKEPYESIDRWQRPVELVSWNQAVEFCQRLSKLTGRKYRLPSEAEWEYACRAGTSTPFYFGETITSELANYDGNHTYGKGPKGENREQTTPVGKFPPNAFGLYDMHGNVDEWCADDWYDNYVGAPTDGSARINSNKNSNTNSTKVFRGGSWIYYPWWCRSALRGNFNSDEADSISFGFRLVSFPPRTFE                                                                   | 1584   | Formylglycine-generating enzyme, required for sulfatase activity, contains SUMF1/FGE domain (YfmG) (PDB:2AFT)Posttranslational modification, protein turnover, chaperones                                                                                      | COG20_FUNCTIONCOG20_CATEGORY | COG1262O                |        |

Of course it is possible to also skip this step entirely in case you somehow have no desire to understand context information via the flag `--skip-recovering-genomic-context`.

### TSV output

The main output fo th etool is a table formatted in a TSV, called DGR-OUTPUT_DGRs_found.tsv. 

Need a dgrs.txt output directory help

### HTML output

### DGRs_found
table example photo



### Computing DGR variability profiling


NB: These next steps are not within anvi'o directly. You will need to [install](https://merenlab.org/2014/08/16/installing-the-oligotyping-pipeline/) [oligotyping](https://merenlab.org/software/oligotyping/) to do this process.


Visualising the variability of the variable regions produced in the context of your samples
The user give the tools samples that were used to create the Profile.db, where they also have a samples.txt with paths to the short reads of these samples, obviously, also need the fastq files of the short reads.

All VRs are forced to be on the forward strand, this includes reverse complementing them to fit this narrative. There is then an initial primer region that is by default 12 bp in front of the start of the VR.

#### Creating different primers for every sample

## Power Users:
#### Discovery mode
Default use of the program is that it will locate SNVs based on their codon position so they need to be in ORFs and the SNVs are limited to the first and second codon position.

Discovery mode however, places no limits on where the SNVs are located. This allows for a much broader search area of investigation and goes against some of the known features of DGRs.

#### Collections mode
Searching for DGRs in metagenomes. Uses collections, more specifically the bins in the specified collection by the user and checks that the TR and VRs are found in the same bin. Because in theory the bins should each be one sample ('species' [notice the quotation marks]), therefore you could give the program a metagenome as long as you have binned groups.

Can specify a specific collection with the flag `-C` and then the program will only look through that specific collection for DGRs.

## Changeable parameters:

There are many ways to alter the behavioral specifics of this programme. You can find some commonly adjusted parameters below. For a full list of parameters, check the program help (-h) output.


#### BLASTn Arguments
- Step: The length of base pairs you would like to cut your genome/sequence into. Default = 100 (This is only for fasta file mode)
- Word Size: BLASTn word size parameter. Default = 8

#### Locating VR Arguments:
 - skip-Ns: Skip 'N' bases when searching for mismatches
 - skip-dashes: Skip '-' bases when searching for mismatches
 - gene caller: The gene caller to show gene calls if you are using a contigs.db. This is used to tell the program that you want to find the genes that your Variable Regions occurred in.
INPUT DATA

contigs-db (required): Input Contigs.db.

profile-db (required): Input Profile.db, preferably a merged profile database.

CONTIGS AND PROFILE DB INPUT ARGUMENTS

hmm-usage (required):
The name of the HMM profile run on your Contigs.db (e.g., Reverse_Transcriptase or a custom set).
Accepts comma-separated list without spaces.
Available profiles: {available_hmm_sources_pretty}.

gene-caller:
Gene caller to use for locating genes that contain Variable Regions.

COLLECTIONS MODE

collections-mode:
Enable searching only within a specified collection.

collection-name:
Name of the collection (single) to restrict the search to.

SNV CLUSTER LOCATOR PARAMETERS

departure-from-reference-percentage:
Minimum departure from reference to count an SNV. Default = 0.1.

minimum-snv-density:
Minimum proportion of SNVs over the VR window. Default = 0.2.

distance-between-snv:
Max bp distance between SNVs for joining them into a window. Default = 8.

minimum-range-size:
Minimum size of the SNV window. Default = 5.

variable-buffer-length:
bp added to each side of a high-SNV region. Default = 35.

BLASTn ARGUMENTS

word-size:
BLASTn word size. Default = 8.

LOCATING VR OPTIONS

skip-Ns: Skip 'N' bases when searching for mismatches. (default: True)

skip-dashes: Skip '-' bases when searching for mismatches. (default: True)

discovery-mode:
Use all SNVs (not only codon positions 1 & 2) to find candidate VRs.
More sensitive but less conservative. (default: False)

FILTER BLAST RESULTS FOR TR/VRs

repeat-threshold:
Max allowed proportion of masked sequence (pytantan). Default = 0.5.

num-imperfect-tandem-repeats:
Max number of imperfect repeats in TR/VR (pytrf). Default = 10.

repeat-motif-coverage:
Coverage threshold of imperfect repeat motifs. Default = 0.8.

percentage-mismatch:
Proportion of mismatching bases in TR (0.5–1.0). Default = 0.8.

number-of-mismatches:
Number of same-type mismatching bases in TR. Default = 7.

only-a-bases:
Only allow mismatches arising from A bases in the TR.

min-mismatching-base-types-vr:
Minimum different mismatch base types in the VR. Default = 2.

min-base-types-tr:
Minimum base types in the TR. Default = 2.

min-base-types-vr:
Minimum base types in the VR. Default = 2.

snv-matching-proportion:
Maximum allowed proportion of SNVs that match TR base positions.

snv-codon-position:
Max fraction of SNVs allowed in 3rd codon position. Default = 0.33.

OUTPUT DIRECTORY

output-dir: Directory to place all outputs.

parameter-output: Output a CSV file of all parameters supplied.

overwrite-output-destinations: Allow overwriting existing outputs.

REPORTING GENOMIC CONTEXT AROUND DGRs

skip-recovering-genomic-context:
Skip collecting surrounding gene context.

num-genes-to-consider-in-context:
Number of upstream/downstream genes to include. Default = 3.

COMPUTING VARIABLE REGION VARIABILITY PROFILING

samples-txt:
Input file containing r1 and r2 short-read paths used for the merged profile.

initial-variable-primer-length:
Length of the primer segment immediately before the VR. Default = 12 bp.

whole-primer-length:
Total length of the designed in-silico primer for VR mapping. Default = 65 bp.

skip-compute-DGR-variability-profiling:
Skip the computationally heavy DGR activity analysis.

CALCULATE ACTIVITY FROM KNOWN DGRs

pre-computed-dgrs:
Use existing DGRs_found.tsv to recalculate DGR activity.
(Most other parameters become irrelevant except the VR-profiling ones.)