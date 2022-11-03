A collection of TAB-delimited text files generated from the profiling of BAM files.

## Example outputs

The number of columns and their content for files that are considered artifact %(bam-stats-txt)s will be variable and depend on the user parameters set for %(anvi-profile-blitz)s.

The column names may be one of these:

* `gene_callers_id`: Unique number assigned by the gene caller during the creation of the %(contigs-db)s.
* `contig`: Contig name as appears in the %(bam-file)s and %(contigs-db)s
* `sample`: The name of the %(bam-file)s without its prefix. I.e., the value *SAMPLE-01* will appear in the *sample* column if the BAM file path was */path/to/SAMPLE-01.bam*.
* `length`: Depending on the context, the length of the gene or contig.
* `num_mapped_reads`: The actual number of short reads mapping to a gene or contig in a given sample. Useful for those who wish to do TPM/RPKM normalizations.
* `detection`: Proportion of nucleotides that have at least 1X coverage.
* `mean_cov`: Mean covearge.
* `q2q3_cov`: Mean of the coverage (inner quartiles).
* `median_cov`: Median coverage.
* `min_cov`: Minimum coverage value observed for the gene or the contig.
* `max_cov`: Minimum coverage value observed for the gene or the contig.
* `std_cov`: Standard deviation of coverage.

### Contig mode, default output

12-column TAB delimited file, where each row represents a single contig x sample pair (so the values in the first column are not unique):

|**contig**|**sample**|**length**|**gc_content**|**num_mapped_reads**|**detection**|**mean_cov**|**q2q3_cov**|**median_cov**|**min_cov**|**max_cov**|**std_cov**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|contig_878|SAMPLE-01|27538|0.608|11877|0.9995|63.61|65.2|65.0|0|107|15.21|
|contig_6515|SAMPLE-01|12315|0.446|7669|0.9985|91.51|92.5|92.0|0|195|23.98|
|contig_1720|SAMPLE-01|16856|0.312|4237|0.9993|37.85|38.25|38.0|0|56|7.961|
|contig_878|SAMPLE-02|27538|0.608|1594|0.9999|9.262|9.161|9.0|0|21|3.42|
|contig_6515|SAMPLE-02|12315|0.446|2562|0.9918|33.05|33.47|33.0|0|56|8.503|
|contig_1720|SAMPLE-02|16856|0.312|926|0.9986|8.93|8.751|9.0|0|19|3.306|
|contig_878|SAMPLE-03|27538|0.608|6395|1.0|37.32|37.21|37.0|0|75|11.46|
|contig_6515|SAMPLE-03|12315|0.446|300|0.9276|3.953|3.682|4.0|0|15|2.644|
|contig_1720|SAMPLE-03|16856|0.312|18175|1.0|178.1|178.1|178.0|1|269|29.13|

### Contig mode, minimal output:

7-column TAB delimited file, where each row represents a single contig x sample pair:


|**contig**|**sample**|**length**|**gc_content**|**num_mapped_reads**|**detection**|**mean_cov**|
|:--|:--|:--|:--|:--|:--|:--|
|contig_878|SAMPLE-01|27538|0.608|11877|0.9995|63.61|
|contig_6515|SAMPLE-01|12315|0.446|7669|0.9985|91.51|
|contig_1720|SAMPLE-01|16856|0.312|4237|0.9993|37.85|
|contig_878|SAMPLE-02|27538|0.608|1594|0.9999|9.262|
|contig_6515|SAMPLE-02|12315|0.446|2562|0.9918|33.05|
|contig_1720|SAMPLE-02|16856|0.312|926|0.9986|8.93|
|contig_878|SAMPLE-03|27538|0.608|6395|1.0|37.32|
|contig_6515|SAMPLE-03|12315|0.446|300|0.9276|3.953|
|contig_1720|SAMPLE-03|16856|0.312|18175|1.0|178.1|

### Gene mode, default output

12-column TAB delimited file, where each row represents a single gene x sample pair:

|**gene_callers_id**|**contig**|**sample**|**length**|**num_mapped_reads**|**detection**|**mean_cov**|**q2q3_cov**|**median_cov**|**min_cov**|**max_cov**|**std_cov**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|0|contig_878|SAMPLE-01|933|385|0.9871|53.97|58.9|59.0|0|85|21.18|
|1|contig_878|SAMPLE-01|564|326|1.0|66.97|66.88|66.0|35|103|20.6|
|2|contig_878|SAMPLE-01|444|318|1.0|81.72|81.59|82.0|70|95|4.88|
|3|contig_878|SAMPLE-01|1218|522|1.0|54.45|57.1|58.0|15|87|17.8|
|4|contig_878|SAMPLE-01|3381|1476|1.0|60.89|60.96|61.0|19|95|13.66|
|5|contig_878|SAMPLE-01|942|472|1.0|64.34|63.98|63.0|38|92|10.49|
|6|contig_878|SAMPLE-01|588|320|1.0|67.51|66.18|66.0|51|92|9.591|
|7|contig_878|SAMPLE-01|1854|852|1.0|62.63|63.03|63.0|31|85|10.14|
|8|contig_878|SAMPLE-01|285|195|1.0|67.43|68.41|70.0|51|80|8.741|
|9|contig_878|SAMPLE-01|1215|567|1.0|60.96|63.68|64.0|16|83|13.94|
|10|contig_878|SAMPLE-01|2250|1018|1.0|62.36|62.91|62.0|9|107|18.65|
|11|contig_878|SAMPLE-01|741|433|1.0|70.23|70.42|71.0|44|94|11.49|
|12|contig_878|SAMPLE-01|963|470|1.0|63.49|65.79|65.0|24|88|13.8|
|13|contig_878|SAMPLE-01|684|310|1.0|56.33|57.69|58.0|26|85|13.38|
|14|contig_878|SAMPLE-01|1569|724|1.0|61.79|63.8|64.0|10|95|18.3|
|15|contig_878|SAMPLE-01|1584|775|1.0|65.69|66.14|67.0|44|88|8.792|
|16|contig_878|SAMPLE-01|831|456|1.0|67.95|67.82|67.0|48|91|8.154|
|17|contig_878|SAMPLE-01|192|179|1.0|81.12|81.14|82.0|69|91|5.041|
|18|contig_878|SAMPLE-01|1467|675|1.0|60.06|60.98|62.0|25|91|14.24|
|19|contig_878|SAMPLE-01|801|456|1.0|68.31|67.72|68.0|58|86|5.945|
|20|contig_878|SAMPLE-01|360|252|1.0|71.42|72.38|74.0|53|87|8.963|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|0|contig_878|SAMPLE-02|933|58|1.0|8.17|7.548|8.0|1|21|4.7|
|1|contig_878|SAMPLE-02|564|41|1.0|8.151|8.28|8.0|2|14|3.13|
|2|contig_878|SAMPLE-02|444|42|1.0|10.84|10.31|10.0|7|15|2.128|
|3|contig_878|SAMPLE-02|1218|83|1.0|9.86|10.55|11.0|1|15|3.041|
|4|contig_878|SAMPLE-02|3381|224|1.0|10.17|9.925|10.0|3|18|2.796|
|5|contig_878|SAMPLE-02|942|78|1.0|11.07|10.98|11.0|6|17|2.34|
|6|contig_878|SAMPLE-02|588|47|1.0|10.09|9.296|9.0|5|18|3.169|
|7|contig_878|SAMPLE-02|1854|121|1.0|9.417|9.186|9.0|2|16|2.75|
|8|contig_878|SAMPLE-02|285|30|1.0|10.33|10.0|10.0|7|15|2.217|
|9|contig_878|SAMPLE-02|1215|79|1.0|9.386|8.685|9.0|3|20|3.965|
|10|contig_878|SAMPLE-02|2250|115|0.9991|7.619|7.98|8.0|0|14|2.97|
|11|contig_878|SAMPLE-02|741|58|1.0|10.82|11.18|11.0|3|16|3.067|
|12|contig_878|SAMPLE-02|963|46|1.0|6.849|7.004|7.0|2|12|2.912|
|13|contig_878|SAMPLE-02|684|36|1.0|7.281|7.029|8.0|2|14|3.168|
|14|contig_878|SAMPLE-02|1569|74|1.0|6.505|6.345|6.0|1|13|2.363|
|15|contig_878|SAMPLE-02|1584|102|1.0|9.199|9.064|9.0|4|15|2.398|
|16|contig_878|SAMPLE-02|831|60|1.0|10.59|10.73|11.0|6|15|2.31|
|17|contig_878|SAMPLE-02|192|16|1.0|8.208|8.854|9.0|3|11|2.857|
|18|contig_878|SAMPLE-02|1467|108|1.0|10.67|10.36|11.0|4|20|3.894|
|19|contig_878|SAMPLE-02|801|68|1.0|10.85|10.74|11.0|5|19|2.706|
|20|contig_878|SAMPLE-02|360|34|1.0|10.32|10.24|10.0|6|15|1.886|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Gene mode, minimal output:

6-column TAB delimited file, where each row represents a single gene x sample pair:

|gene_callers_id|contig|sample|length|detection|mean_cov|
|:--|:--:|:--:|:--:|:--:|:--:|
|0|contig_878|SAMPLE-01|933|0.9871|53.97|
|1|contig_878|SAMPLE-01|564|1.0|66.97|
|2|contig_878|SAMPLE-01|444|1.0|81.72|
|3|contig_878|SAMPLE-01|1218|1.0|54.45|
|4|contig_878|SAMPLE-01|3381|1.0|60.89|
|5|contig_878|SAMPLE-01|942|1.0|64.34|
|6|contig_878|SAMPLE-01|588|1.0|67.51|
|7|contig_878|SAMPLE-01|1854|1.0|62.63|
|8|contig_878|SAMPLE-01|285|1.0|67.43|
|9|contig_878|SAMPLE-01|1215|1.0|60.96|
|10|contig_878|SAMPLE-01|2250|1.0|62.36|
|11|contig_878|SAMPLE-01|741|1.0|70.23|
|12|contig_878|SAMPLE-01|963|1.0|63.49|
|13|contig_878|SAMPLE-01|684|1.0|56.33|
|14|contig_878|SAMPLE-01|1569|1.0|61.79|
|15|contig_878|SAMPLE-01|1584|1.0|65.69|
|16|contig_878|SAMPLE-01|831|1.0|67.95|
|17|contig_878|SAMPLE-01|192|1.0|81.12|
|18|contig_878|SAMPLE-01|1467|1.0|60.06|
|19|contig_878|SAMPLE-01|801|1.0|68.31|
|20|contig_878|SAMPLE-01|360|1.0|71.42|
|(...)|(...)|(...)|(...)|(...)|(...)|
|0|contig_878|SAMPLE-02|933|1.0|8.17|
|1|contig_878|SAMPLE-02|564|1.0|8.151|
|2|contig_878|SAMPLE-02|444|1.0|10.84|
|3|contig_878|SAMPLE-02|1218|1.0|9.86|
|4|contig_878|SAMPLE-02|3381|1.0|10.17|
|5|contig_878|SAMPLE-02|942|1.0|11.07|
|6|contig_878|SAMPLE-02|588|1.0|10.09|
|7|contig_878|SAMPLE-02|1854|1.0|9.417|
|8|contig_878|SAMPLE-02|285|1.0|10.33|
|9|contig_878|SAMPLE-02|1215|1.0|9.386|
|10|contig_878|SAMPLE-02|2250|0.9991|7.619|
|11|contig_878|SAMPLE-02|741|1.0|10.82|
|12|contig_878|SAMPLE-02|963|1.0|6.849|
|13|contig_878|SAMPLE-02|684|1.0|7.281|
|14|contig_878|SAMPLE-02|1569|1.0|6.505|
|15|contig_878|SAMPLE-02|1584|1.0|9.199|
|16|contig_878|SAMPLE-02|831|1.0|10.59|
|17|contig_878|SAMPLE-02|192|1.0|8.208|
|18|contig_878|SAMPLE-02|1467|1.0|10.67|
|19|contig_878|SAMPLE-02|801|1.0|10.85|
|20|contig_878|SAMPLE-02|360|1.0|10.32|
|(...)|(...)|(...)|(...)|(...)|(...)|

## Reproducing these output files

Examples above generated by running %(anvi-profile-blitz)s in the mini-test output directory. To reproduce them, you can run this command to generate the necessary files,

```
anvi-self-test --suite mini -o TEST
```

then go into the directory,

```
cd TEST
```

and run %(anvi-profile-blitz)s in coresponding modes.
