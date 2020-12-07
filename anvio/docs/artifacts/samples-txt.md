A **TAB-delimited** file to describe samples and paired-end FASTQ files associated with them. By doing so, this file type, which we commonly refer to as %(samples-txt)s, is a way to links sample names to raw sequencing reads.

This file type includes required and optional columns.

The following three columns are **required** for this file type:

* `sample`: a single-word sample name,
* `r1`: path to the FASTQ file for pair one, and
* `r2`: path to the FASTQ file for pair two.

While you can use relative paths for `r1` and `r2`, it is always better to have absolute paths to improve reproducibility.

The following is an **optional** column:

* `group`: A single-word categorical variable that assigns two or more samples into two or more groups.

### Examples samples.txt file

Here is an example file:

|sample|group|r1|r2|
|:--|:--|:--|:--|
|Sample_01|WARM|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz|
|Sample_02|COLD|/path/to/YYY-02-R1.fastq.gz|/path/to/YYY-02-R2.fastq.gz|
|Sample_03|COLD|/path/to/ZZZ-03-R1.fastq.gz|/path/to/ZZZ-03-R2.fastq.gz|