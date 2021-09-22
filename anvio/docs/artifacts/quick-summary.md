The output of %(anvi-quick-summary)s.

%(anvi-quick-summary)s summarizes read-recruitment statistics for a collection of bins across multiple samples. It produces long format output in which each row contains the (weighted) average statistics of a bin in a sample. Each statistic is summarized in a different column of the file.

Here is an example output file from this program, summarizing detection and mean_coverage_Q2Q3 data for 3 bins across multiple samples:

unique_id | bin_name | sample | detection | mean_coverage_Q2Q3
|:---|:---|:---|:---|:---|
0 | bin_1 | sample_1 | 0.015553023620503776 | 1.0272713907674214
1 | bin_2 | sample_1 | 0.0004871607502275562 | 0.0
2 | bin_3 | sample_1 | 0.0023636043452898497 | 0.0
3 | bin_1 | sample_2 | 0.015767421346662747 | 1.1101759286484367
4 | bin_2 | sample_2 | 0.0004871607502275562 | 0.0
5 | bin_3 | sample_2 | 0.001595914458984989 | 0.0
[...] | [...] |[...] |[...] |[...]
