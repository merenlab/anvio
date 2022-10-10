This is the output of %(anvi-report-inversions)s. 

### Per sample inversion report table

%(anvi-report-inversions)s searches for inversions in every single sample at a time and thus genereates a TAB-delimited table for every sample: `INVERSIONS-IN-SAMPLE_01.txt`, `INVERSIONS-IN-SAMPLE_02`, ...

Here is an example output:

|**entry_id**|**sample_name**|**contig_name**|**first_seq**|**midline**|**second_seq**|**first_start**|**first_end**|**first_oligo_primer**|**first_oligo_reference**|**second_start**|**second_end**|**second_oligo_primer**|**second_oligo_reference**|**num_mismatches**|**num_gaps**|**length**|**distance**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|c_000000000001_10214_10541|S01|c_000000000001|TGTTTCAAAAAACGTTCGT|\|\|\|\|\|\|\|\|\|\|x\|\|\|\|\|\|\|\||ACGAACGTCTTTTGAAACA|10306|10325|TCGATCAATTGATGTTTCAAAA.ACGTTCGT|CTTTTG|10493|10512|TCAGTAGTGAATGTGTTTCAAAA.ACGTTCGT|TTAATA|1|0|19|168|
|c_000000000002_10148_11135|S01|c_000000000002|TGTTTCAAAAAACGTTCGT|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\||ACGAACGTTTTTTGAAACA|10514|10533|AGATACGGTTTATGTTTCAAAAAACGTTCGT|CTTTTG|10724|10743|AGAAAAGAAGGCGTGTTTCAAAAAACGTTCGT|TCAATG|0|0|19|191|

These tables contains the following columns:

* Entry ID made with the contig's name and the start and stop position of the stretch
* The contig's name
* The first palindrome sequence
* The aligment midline
* The second palindrome sequence
* The start and stop position of the first and second palindrome sequence
* The number of mismatches
* The number of gaps
* The length of the palindrome sequence
* The distance between the first and second palindrome seqeuences, i.e. the size of the inversion
* The number of samples in which it was detected and confirmed
* The in silico primers used to compute the inversion's activity, for the first and second palindrome
* The oligo corresponding to the reference sequence

### Inversions consensus table

Anvi'o eventually create a consensus table with all the unique inversions found accross all your samples in a file called `INVERSIONS-CONSENSUS.txt`. This table has the same format as the individual sample outputs, with the 'entry ID' replaced by a unique inversion ID. It also has column reporting the samples where the inversion was detected.

The table should look like this:

|**inversion_id**|**contig_name**|**first_seq**|**midline**|**second_seq**|**first_start**|**first_end**|**second_start**|**second_end**|**num_mismatches**|**num_gaps**|**length**|**distance**|**num_samples**|**sample_names**|**first_oligo_primer**|**first_oligo_reference**|**second_oligo_primer**|**second_oligo_reference**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|INV_0001|c_000000000001|TGTTTCAAAAAACGTTCGT|\|\|\|\|\|\|\|\|\|\|x\|\|\|\|\|\|\|\||ACGAACGTCTTTTGAAACA|10306|10325|10493|10512|1|0|19|168|3|S01,S02,S03|TCGATCAATTGATGTTTCAAAA.ACGTTCGT|CTTTTG|TCAGTAGTGAATGTGTTTCAAAA.ACGTTCGT|TTAATA|
|INV_0002|c_000000000002|TGTTTCAAAAAACGTTCGT|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\||ACGAACGTTTTTTGAAACA|10514|10533|10724|10743|0|0|19|191|3|S01,S02,S03|AGATACGGTTTATGTTTCAAAAAACGTTCGT|CTTTTG|AGAAAAGAAGGCGTGTTTCAAAAAACGTTCGT|TCAATG|

### All stretches considered

Another default output table is named `ALL-STRETCHES-CONSIDERED.txt` and it reports every stretch that passed the 'Identifying regions of interest' parameters. It reports the maximum coverage of FWD/FWD and REV/REV in that stretch, per sample. It also reports the number of palindromes found and if a true inversion was confirmed.

|**entry_id**|**sequence_name**|**sample_name**|**contig_name**|**start_stop**|**max_coverage**|**num_palindromes_found**|**true_inversions_found**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|S01_c_000000000001_10214_10541|c_000000000001_10214_10541|S01|c_000000000001|10214-10541|16|4|False|
|S01_c_000000000002_10148_11135|c_000000000002_10148_11135|S01|c_000000000002|10148-11135|69|20|False|
|S02_c_000000000001_10283_10542|c_000000000001_10283_10542|S02|c_000000000001|10283-10542|11|3|False|
|S02_c_000000000002_10200_11052|c_000000000002_10200_11052|S02|c_000000000002|10200-11052|96|17|False|
|S03_c_000000000001_10033_10801|c_000000000001_10033_10801|S03|c_000000000001|10033-10801|30|12|False|
|S03_c_000000000002_10498_10764|c_000000000002_10498_10764|S03|c_000000000002|10498-10764|13|2|False|

### Surrounding genes and functions

If the user enable the reporting of the genomic context, two addition TAB-delimited tables are generated: `INVERSIONS-CONSENSUS-SURROUNDING-GENES.txt` and `INVERSIONS-CONSENSUS-SURROUNDING-FUNCTIONS.txt`.

The first table report the gene calls surrounging every inversion when possible:

|**inversion_id**|**entry_type**|**gene_callers_id**|**start**|**stop**|**direction**|**partial**|**call_type**|**source**|**version**|**contig**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|INV_0001|FIRST_PALINDROME||10306|10325||||||c_000000000001|
|INV_0001|SECOND_PALINDROME||10493|10512||||||c_000000000001|
|INV_0001|GENE|6|6818|7595|r|0|1|prodigal|v2.6.3|c_000000000001|
|INV_0001|GENE|7|7632|8976|r|0|1|prodigal|v2.6.3|c_000000000001|
|INV_0001|GENE|8|9145|9616|r|0|1|prodigal|v2.6.3|c_000000000001|
|INV_0001|GENE|9|9651|10170|r|0|1|prodigal|v2.6.3|c_000000000001|
|INV_0001|GENE|10|11311|14161|r|0|1|prodigal|v2.6.3|c_000000000001|
|INV_0001|GENE|11|14165|14495|r|0|1|prodigal|v2.6.3|c_000000000001|
|INV_0001|GENE|12|14524|16072|r|0|1|prodigal|v2.6.3|c_000000000001|

The second table report the function associated to every gene call reported in the first file:

|**inversion_id**|**gene_callers_id**|**source**|**accession**|**function**|
|:--|:--|:--|:--|:--|
|INV_0001|6|COG20_FUNCTION|COG1208|NDP-sugar pyrophosphorylase, includes eIF-2Bgamma, eIF-2Bepsilon, and LPS biosynthesis protein s (GCD1) (PDB:6JQ8)|
|INV_0001|6|COG20_CATEGORY|J|Translation, ribosomal structure and biogenesis|
|INV_0001|6|KOfam|K00978|glucose-1-phosphate cytidylyltransferase [EC:2.7.7.33]|
|INV_0001|7|COG20_FUNCTION|COG0399|dTDP-4-amino-4,6-dideoxygalactose transaminase (WecE) (PDB:4PIW) (PUBMED:15271350)|
|INV_0001|7|COG20_CATEGORY|M|Cell wall/membrane/envelope biogenesis|
|INV_0001|7|KOfam|K12452|CDP-4-dehydro-6-deoxyglucose reductase, E1 [EC:1.17.1.1]|
|INV_0001|9|COG20_FUNCTION|COG0250|Transcription termination/antitermination protein NusG (NusG) (PDB:3EWG) (PUBMED:19500594)|
|INV_0001|9|COG20_CATEGORY|K|Transcription|
|INV_0001|10|COG20_FUNCTION|COG2605|Predicted kinase related to galactokinase and mevalonate kinase (PDB:4USK)|
|INV_0001|10|COG20_CATEGORY|R|General function prediction only|
|INV_0001|10|KOfam|K05305|fucokinase [EC:2.7.1.52]|
|INV_0001|11|COG20_FUNCTION|COG3254|L-rhamnose mutarotase (RhaM) (PDB:1X8D)|
|INV_0001|11|COG20_CATEGORY|M|Cell wall/membrane/envelope biogenesis|
|INV_0001|11|KOfam|K03534|L-rhamnose mutarotase [EC:5.1.3.32]|
|INV_0001|12|COG20_FUNCTION|COG0305|Replicative DNA helicase (DnaB) (PDB:1B79)|
|INV_0001|12|COG20_CATEGORY|L|Replication, recombination and repair|
|INV_0001|12|KOfam|K02314|replicative DNA helicase [EC:3.6.4.12]|


### Inversion's activity

Finally, if the user provide R1 and R2 fastq files and enable the reporting of inversion's activity, %(anvi-report-inversions)s will generate a long-format file named `INVERSION-ACTIVITY.txt`. This file reports, for every inversion and sample, the relative proportion and read abundance of unique oligos, which either correspond to the reference contig (no inversion), or to an inversion sequence. The inversion's activity is computed and reported for both side of each inversion.

|**sample**|**inversion_id**|**oligo_primer**|**oligo**|**reference**|**frequency_count**|**relative_abundance**|
|:--|:--|:--|:--|:--|:--|:--|
|S01|INV_0001|first_oligo_primer|TTAATA|False|6|0.097|
|S01|INV_0001|first_oligo_primer|CTTTTG|True|55|0.887|
|S01|INV_0001|first_oligo_primer|CTTTTT|False|1|0.016|
|S01|INV_0001|second_oligo_primer|CTTTTG|False|11|0.169|
|S01|INV_0001|second_oligo_primer|TTAATA|True|54|0.831|
|S01|INV_0002|first_oligo_primer|TCAATG|False|37|0.587|
|S01|INV_0002|first_oligo_primer|CTTTTG|True|25|0.397|
|S01|INV_0002|first_oligo_primer|TCAATT|False|1|0.016|
|S01|INV_0002|second_oligo_primer|CTTTTG|False|53|0.609|
|S01|INV_0002|second_oligo_primer|TCAATG|True|33|0.379|
|S01|INV_0002|second_oligo_primer|TCAATC|False|1|0.011|
