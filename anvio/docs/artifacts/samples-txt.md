A **TAB-delimited** file to describe samples and FASTQ files associated with them. By doing so, this file type links sample names to raw sequencing reads.

This file type includes required and optional columns.

{:.notice}
While these required and optional columns are what anvi'o is going to look for anytime you expect to process a TAB-delimited file as %(samples-txt)s, you can have as many columns as you like in a given TAB-delimited to be used as %(samples-txt)s as long as it includes these required and optional columns.

The only **required** column is `sample` or `name`, which should be a single-word sample name.

The other required columns will depend on the type of sequencing data you have. For example, if you have paired-end Illumina reads, the
required columns will be `r1` and `r2`, which should point to the FASTQ files for pair one and pair two, respectively. If you have single-end Illumina reads, the only required column will be `r1`.
If you have long reads, the only required column will be `lr`, which should point to the FASTQ file for long reads. Long reads must be FASTQ (`.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz`): anvi'o's long-read steps rely on per-base quality scores, so FASTA long reads are not supported.

{:.notice}
The paths to FASTQ files can be absolute or relative to the location of the samples-txt file.

### Paired-end reads example
For paired-end reads, the samples-txt file should have at least the following three columns:

|sample|r1|r2|
|:--|:--|:--|
|Sample_01|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz|
|Sample_02|/path/to/YYY-02-R1.fastq.gz|/path/to/YYY-02-R2.fastq.gz|
|Sample_03|/path/to/ZZZ-03-R1.fastq.gz|/path/to/ZZZ-03-R2.fastq.gz|

### Single-end reads example
For single-end reads, the samples-txt file should have at least the following two columns:

|sample|r1|
|:--|:--|
|Sample_01|/path/to/XXX-01-R1.fastq.gz|
|Sample_02|/path/to/YYY-02-R1.fastq.gz|
|Sample_03|/path/to/ZZZ-03-R1.fastq.gz|

### Long reads example
For long reads, the samples-txt file should have at least the following two columns:

|sample|lr|
|:--|:--|
|Sample_01|/path/to/XXX-01-lr.fastq.gz|
|Sample_02|/path/to/YYY-02-lr.fastq.gz|
|Sample_03|/path/to/ZZZ-03-lr.fastq.gz|

### Mixed reads example
If you have a mix of paired-end, single-end, and long reads, your samples-txt file should have at least the following columns:

|sample|r1|r2|lr|
|:--|:--|:--|:--|
|Sample_01|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz||
|Sample_02|/path/to/YYY-02-R1.fastq.gz||/path/to/YYY-02-lr.fastq.gz|
|Sample_03|||/path/to/ZZZ-03-lr.fastq.gz|

### Additional optional column: group

The following is an **optional** column:

* `group`: A single-word categorical variable that assigns two or more samples into two or more groups. This is useful to co-assemble multiple samples so that you can bin them later.

For more information, see the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt)

Here is an example samples.txt file with the optional `group` column in addition to the required columns for paired-end reads:

|sample|group|r1|r2|
|:--|:--|:--|:--|
|Sample_01|WARM|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz|
|Sample_02|COLD|/path/to/YYY-02-R1.fastq.gz|/path/to/YYY-02-R2.fastq.gz|
|Sample_03|COLD|/path/to/ZZZ-03-R1.fastq.gz|/path/to/ZZZ-03-R2.fastq.gz|

### Additional optional column: lr_technology

The following is an **optional** column relevant only to long reads (i.e., samples with an `lr` path):

* `lr_technology`: The long-read sequencing technology used for a sample. When you provide this column, anvi'o automatically selects the appropriate presets for each long-read tool it runs — the minimap2 mapping preset and the Flye read-type flag — so you do not have to set them by hand in your %(workflow-config)s.

The accepted values are:

|lr_technology|Description|
|:--|:--|
|`ont`|Oxford Nanopore (any kit)|
|`pb-clr`|PacBio CLR (RS II / Sequel / Sequel II)|
|`pb-hifi`|PacBio HiFi (CCS)|

A few rules govern this column:

* It is **optional**. If you omit it entirely, anvi'o falls back to the presets you set explicitly in your %(workflow-config)s (e.g. `minimap2: {"preset": ...}` and the Flye read-type flag). Anvi'o will not run these tools with their built-in defaults, so if the column is absent you must set those presets yourself or the workflow will stop with an error telling you exactly what is missing.
* It is **all-or-nothing**. If you include the column, then *every* sample that has long reads must have a value in it. Short-read-only samples should leave it blank.

Here is an example samples.txt file with the optional `lr_technology` column for a set of long-read samples:

|sample|lr|lr_technology|
|:--|:--|:--|
|Sample_01|/path/to/XXX-01-lr.fastq.gz|ont|
|Sample_02|/path/to/YYY-02-lr.fastq.gz|pb-hifi|

And here is a mixed short/long-read example, where the short-read-only sample leaves `lr_technology` blank:

|sample|r1|r2|lr|lr_technology|
|:--|:--|:--|:--|:--|
|Sample_01|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz|/path/to/XXX-01-lr.fastq.gz|ont|
|Sample_02|/path/to/YYY-02-R1.fastq.gz|/path/to/YYY-02-R2.fastq.gz|||
|Sample_03|||/path/to/ZZZ-03-lr.fastq.gz|ont|
