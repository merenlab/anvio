By default, anvi'o uses Prodigal for gene calling when the user is generating a %(contigs-db)s. Yet, if the user provides an external gene calls file, then anvi'o does not perform gene calling, and uses this file to store the gene information into the new %(contigs-db)s.

External gene calls is a user-provided TAB-delimited file that should follow this format:

|gene_callers_id|contig|start|stop|direction|partial|call_type|source|version|
|:---:|:---|:---:|:---:|:---:|:---:|:--:|:---:|
|1|contig_01|1113|1677|f|0|1|program|v1.0|
|2|contig_01|1698|2142|f|0|1|program|v1.0|
|3|contig_01|2229|3447|f|0|1|program|v1.0|
|4|contig_01|3439|6820|r|0|1|program|v1.0|
|7|contig_01|8496|10350|r|1|1|program|v1.0|
|8|contig_02|306|1650|f|0|1|program|v1.0|
|9|contig_02|1971|3132|f|0|1|program|v1.0|
|10|contig_02|3230|4007|f|0|1|program|v1.0|
|11|contig_02|4080|5202|f|0|1|program|v1.0|
|12|contig_02|5194|5926|f|0|1|program|v1.0|
|13|contig_03|606|2514|f|0|1|program|v1.0|
|14|contig_03|2751|3207|f|0|1|program|v1.0|
|15|contig_03|3219|5616|f|0|1|program|v1.0|
|16|contig_03|5720|6233|f|0|1|program|v1.0|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Please note that while anvi'o will not perform any gene prediction, it will still try to translate DNA sequence found in start-stop positions of each gene using the standard genetic code. You can prevent that by providing your own amino acid sequences by adding an optional column to your external gene calls file, `aa_sequence`:

|gene_callers_id|contig|start|stop|direction|partial|call_type|source|version|aa_sequence|
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:--:|:--|
|1|contig_01|1113|1677|f|0|1|program|v1.0|MAQTTNDIKNGSVLNLDGQLWTVI(...)|
|2|contig_01|1698|2142|f|0|1|program|v1.0|MARSTARKRALNTLYEADEKGQDI(...)|
|3|contig_01|2229|3447|f|0|1|program|v1.0|MNQYDSEAVMFDPQDAVLVLEDGQ(...)|
|4|contig_01|3439|6820|r|0|1|program|v1.0|MPKRTDIKSVMVIGSGPIVIGQAA(...)|
|7|contig_01|8496|10350|r|1|1|program|v1.0|MMSSPSSEEVNAQRSDFGLRLSNS(...)|
|8|contig_02|306|1650|f|0|1|program|v1.0|MADSQHGRLIVLCGPAGVGKGTVL(...)|
|9|contig_02|1971|3132|f|0|1|program|v1.0|MRSAKLMNGRVFAGARALYRAAGV(...)|
|10|contig_02|3230|4007|f|0|1|program|v1.0|MAFGTEPTPTGLADPPIDDLMEHA(...)|
|11|contig_02|4080|5202|f|0|1|program|v1.0|MAELKLISAESVTEGHPDKVCDQI(...)|
|12|contig_02|5194|5926|f|0|1|program|v1.0|MRYPCIMTNEDAEQLALDGLAPRK(...)|
|13|contig_03|606|2514|f|0|1|program|v1.0|MTLTLRMEKRMKGWPGEPQMEYDV(...)|
|14|contig_03|2751|3207|f|0|1|program|v1.0|MLKVLFAGTPDVAVPSLKLLAQDT(...)|
|15|contig_03|3219|5616|f|0|1|program|v1.0|MLEQETPNIASMASLPTLSAPGLL(...)|
|16|contig_03|5720|6233|f|0|1|program|v1.0|MLESEVDMNDHDEETLASLQQAND(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Explicitly defining amino acid sequences could be particularly useful when working with eukaryotic genomes, and/or genomes that use non-standard genetic code. Sections below discuss specific information about columns of this file.

### Gene start/stop positions

Anvi'o follows the convention of string indexing and splicing that is identical to the way one does it in Python or C. This means that the index of the first nucleotide in any contig should be `0`. In other words, for a gene call that starts at the position `x`th position and ends at position `y`th position, we start counting from `x-1`, and not from `x` (but we still end at `y`). The `start` and `stop` positions in the input file should comply with this criterion. Here is an example gene in a contig:

``` bash
                 1         2         3
nt pos: 12345678901234567890123456789012 (...)
contig: NNNATGNNNNNNNNNNNNNNNNNTAGAAAAAA (...)
           |______ gene X _______|
```

The `start` and `stop` positions in the input file for this gene should be `3` and `26`, respectively. Which means, if you are trying to generate an external gene calls file from gene calls produced by a gene caller that reports start/stop positions starting with the index of `1` rather than `0`, you basically need to subtract one from the start position of every gene call for a matching anvi'o external gene calls file.

Gene `start` and `stop` positions do not care about the direction of the gene as they simply address how the gene sequence should be sliced out from a longer sequence. Whether a gene is forward or reverse is defined in the column `direction`.

{:.notice}
You can read the previous discussions regarding this behavior in [this issue](https://github.com/meren/anvio/issues/374)). Thanks for your patience!

### Call type

{:.notice}
This is a feature added after anvi'o `v6.2`. If you are using anvi'o `v6.2` or earlier, please remove `call_type` column from your external gene calls file.

The column `call_type` declares the nature of the call. It can take one of the following three integer values:

* `1`, indicates that the gene call is for a CODING gene. For gene calls marked as CODING genes, anvi'o will try to predict the proper coding frame when %(anvi-gen-contigs-database)s is run using Markov models trained on a large number of protein sequences and was first described in this [pull request](https://github.com/merenlab/anvio/pull/1428). This is the default behavior for CODING sequences regardless of whether the gene call is partial or not. However, there are two ways the user can change this: (1) by providing an amino acid sequence for the call in the `aa_sequence` column or (2) by asking %(anvi-gen-contigs-database)s to `--skip-predict-frame`.

* `2`, indicates that the gene call is for a NONCODING gene. This is used for non-coding RNAs (transfer RNAs or ribosomal RNAs). For gene calls marked as NONCODING, anvi'o will not attempt to predict an amino acid sequence (nor it will tolerate entries in the `aa_sequence` column).

* `3`, indicates that the gene call is for an UNKNOWN genomic region. This is currently reserved for experimental purposes.
