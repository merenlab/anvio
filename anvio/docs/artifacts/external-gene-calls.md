
By default, anvi'o uses Prodigal for gene calling when the user is generating a %(contigs-db)s. Yet, if the user provides an external gene calls file, then anvi'o does not perform gene calling, and uses this file to store the gene information into the new %(contigs-db)s.

External gene calls is a user-provided TAB-delimited file that should follow this format:

|gene_callers_id|contig|start|stop|direction|partial|source|version|
|:---:|:---|:---:|:---:|:---:|:---:|:---:|
|1|contig_01|1113|1677|f|0|program|v1.0|
|2|contig_01|1698|2142|f|0|program|v1.0|
|3|contig_01|2229|3447|f|0|program|v1.0|
|4|contig_01|3439|6820|r|0|program|v1.0|
|7|contig_01|8496|10350|r|1|program|v1.0|
|8|contig_02|306|1650|f|0|program|v1.0|
|9|contig_02|1971|3132|f|0|program|v1.0|
|10|contig_02|3230|4007|f|0|program|v1.0|
|11|contig_02|4080|5202|f|0|program|v1.0|
|12|contig_02|5194|5926|f|0|program|v1.0|
|13|contig_03|606|2514|f|0|program|v1.0|
|14|contig_03|2751|3207|f|0|program|v1.0|
|15|contig_03|3219|5616|f|0|program|v1.0|
|16|contig_03|5720|6233|f|0|program|v1.0|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Please note that while anvi'o will not perform any gene prediction, it will still try to translate DNA sequence found in start-stop positions of each gene using the standard genetic code. You can prevent that by providing your own amino acid sequences by adding an optional column to your external gene calls file, `aa_sequence`:

|gene_callers_id|contig|start|stop|direction|partial|source|version|aa_sequence|
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:--|
|1|contig_01|1113|1677|f|0|program|v1.0|MAQTTNDIKNGSVLNLDGQLWTVI(...)|
|2|contig_01|1698|2142|f|0|program|v1.0|MARSTARKRALNTLYEADEKGQDI(...)|
|3|contig_01|2229|3447|f|0|program|v1.0|MNQYDSEAVMFDPQDAVLVLEDGQ(...)|
|4|contig_01|3439|6820|r|0|program|v1.0|MPKRTDIKSVMVIGSGPIVIGQAA(...)|
|7|contig_01|8496|10350|r|1|program|v1.0|MMSSPSSEEVNAQRSDFGLRLSNS(...)|
|8|contig_02|306|1650|f|0|program|v1.0|MADSQHGRLIVLCGPAGVGKGTVL(...)|
|9|contig_02|1971|3132|f|0|program|v1.0|MRSAKLMNGRVFAGARALYRAAGV(...)|
|10|contig_02|3230|4007|f|0|program|v1.0|MAFGTEPTPTGLADPPIDDLMEHA(...)|
|11|contig_02|4080|5202|f|0|program|v1.0|MAELKLISAESVTEGHPDKVCDQI(...)|
|12|contig_02|5194|5926|f|0|program|v1.0|MRYPCIMTNEDAEQLALDGLAPRK(...)|
|13|contig_03|606|2514|f|0|program|v1.0|MTLTLRMEKRMKGWPGEPQMEYDV(...)|
|14|contig_03|2751|3207|f|0|program|v1.0|MLKVLFAGTPDVAVPSLKLLAQDT(...)|
|15|contig_03|3219|5616|f|0|program|v1.0|MLEQETPNIASMASLPTLSAPGLL(...)|
|16|contig_03|5720|6233|f|0|program|v1.0|MLESEVDMNDHDEETLASLQQAND(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Explicitly defining amino acid sequences could be particularly useful when working with eukaryotic genomes, or genomes that use non-standard genetic code.

{:.notice}
**Please note**: If gene calls are partial (i.e. they take the value `1` under the `partial` column), or if complete genes (i.e. they take the value `0` under the `partial` column) have lengths (`stop - start`) that are indivisible by 3, anvi'o by default does not attempt to translate their sequence, since there is ambiguity in which codon frame should be used for translation. However, during `anvi-gen-contigs-database`, one can provide `--predict-frame`, which causes anvi'o to find the best codon frame for the given gene call. This is not the default behavior, as it gives anvi'o the responsibility and power to slightly trim the `start`/`stop` that you explicitly told it to use. Read more in this pull request where the feature was [first introduced](https://github.com/merenlab/anvio/pull/1428).

{:.notice}
**Please note**: anvi'o follows the convention of string indexing and splicing that is identical to the way one does it in Python or C.

The statement above means that the index of the first nucleotide in any contig should be `0`. In other words, for a gene call that starts at the position `x`th position and ends at position `y`th position, we start counting from `x-1`, and not from `x` (but we still end at `y`. The `start` and `stop` positions in the input file should comply with this. Here is an example gene in a contig:

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


