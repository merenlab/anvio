A **TAB-delimited** file to describe samples and paired-end FASTQ files associated with them. By doing so, this file type links sample names to raw sequencing reads.

This file type includes required and optional columns.

{:.notice}
While these required and optional columns are what anvi'o is going to look for anytime you expect to process a TAB-delimited file as %(samples-txt)s, you can have as many columns as you like in a given TAB-delimited to be used as %(samples-txt)s as long as it includes these required and optional columns.

The following three columns are **required** for this file type:

* `sample`: a single-word sample name,
* `r1`: path to the FASTQ file for pair one, and
* `r2`: path to the FASTQ file for pair two.

{:.notice}
You can also use `name` as your first column instead of `sample`.

While you can use relative paths for `r1` and `r2`, it is always better to have absolute paths to improve reproducibility.

The following is an **optional** column:

* `group`: A single-word categorical variable that assigns two or more samples into two or more groups. This is useful to co-assemble multiple samples so that you can bin them later. 

For more information, see the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt)

### Examples samples.txt file

Here is an example file:

|sample|group|r1|r2|
|:--|:--|:--|:--|
|Sample_01|WARM|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz|
|Sample_02|COLD|/path/to/YYY-02-R1.fastq.gz|/path/to/YYY-02-R2.fastq.gz|
|Sample_03|COLD|/path/to/ZZZ-03-R1.fastq.gz|/path/to/ZZZ-03-R2.fastq.gz|
