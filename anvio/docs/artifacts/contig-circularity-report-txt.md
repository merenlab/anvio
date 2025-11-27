A report produced by %(anvi-report-circularity)s will have 12 TAB-delimited columns. Here is a snippet from an example output file:

|**`contig`**|**`contig_length`**|**`status`**|**`circularity_support`**|**`edge_coherence`**|**`num_fr_pairs`**|**`expected_rf_pairs`**|**`observed_rf_pairs`**|**`supporting_rf_pairs`**|**`num_edge_reads`**|**`num_edge_reads_incoherent`**|**`flags`**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|c_001045|6184|linear|0.0000|0.6613|3264|204.74|5|0|685|232||
|c_001380|4450|linear|0.0000|0.6452|1497|133.76|1|0|451|160||
|c_001206|4129|circular|0.7459|0.9247|2662|258.14|122|91|930|70||
|c_001137|5761|linear|0.0000|0.8947|97|6.56|0|0|19|2|low_fr_coverage (n=97);low_edge_coverage (n=19)|
|c_001925|2771|circular|0.9275|0.8053|1322|200.55|69|64|673|131||
|c_001947|2784|indeterminate|1.0000|1.0000|47|7.09|3|3|20|0|low_fr_coverage (n=47);ambiguous_evidence|
|c_002285|42768|linear|0.0000|0.9936|26431|227.51|107|0|2195|14||
|c_002729|2044|circular|0.7380|0.9527|2782|604.78|271|200|1967|93||
|c_003377|1704|linear|0.0000|0.6584|2087|568.90|1|0|2181|745||
|c_005100|16122|linear|0.0000|0.9892|9051|209.66|5|0|277|3||
|c_002482|2237|linear|0.0000|0.6667|23|4.48|0|0|9|3|low_fr_coverage (n=23);low_expected_rf_pairs (expected=4.5);low_edge_coverage (n=9)|
|c_009498|7406|linear|0.0000|0.9603|6083|315.34|61|0|378|15||
|c_009666|7209|linear|0.0000|0.9735|7456|397.64|41|0|339|9||
|c_009940|6893|linear|0.0000|0.8895|3494|195.36|16|0|371|41||
|c_007466|10285|circular|0.8000|0.9302|551|20.27|10|8|86|6||

where,

* `contig`: name of the reference sequence in the %(bam-file)s.
* `contig_length`: the actual length of the reference sequence.
* `status`: the classification result, and whether the contig determined to be 'circular', 'linear', or none ('indeterminate').
* `circularity_support`: a score that ranges from 0 to 1 and indicates the strenght of evidence that support circularity.
* `edge_coherence`: a score that ranges from 0 to 1 and indicates the fraction of edge-mapping reads whose mates also map to this contig.
* `num_fr_pairs`: number of `FR` (forward-reverse) pairs (these are the normal pairs that map like normal people do).
* `expected_rf_pairs`: the expected number of `RF` pairs *if* the contig were circular, calculated from `num_fr_pairs` and the geometry of insert size vs contig length.
* `observed_rf_pairs`: the total number of `RF`-oriented read pairs that are actually found on this contig.
* `supporting_rf_pairs`: the umber of `RF` pairs whose circular insert size fell within the expected insert size distribution (i.e., weird pairs that do look like plausible junction-spanning pairs).
* `num_edge_reads`: the number of reads mapping within 'average insert size' length of contig edges.
* `num_edge_reads_incoherent`: the number of edge reads whose mate maps to a different contig.
* `flags`: warning flags for edge cases.
