A **TAB-delimited** file to describe anvi'o %(single-profile-db)s and %(bam-file)s pairs along with the %(contigs-db)s used the profile the BAM file.

This file type includes required and optional columns. The following four columns are **required** for this file type:

* `name`: a single-word name for the entry.
* `contigs_db_path`: path to a %(contigs-db)s.
* `bam_file_path`: path to a %(bam-file)s.
* `profile_db_path`: path to a %(single-profile-db)s generated from the BAM file and the contigs database mentioned.

### Example

Here is an example file:

|name|contigs_db_path|profile_db_path|bam_file_path|
|:--|:--:|:--:|:--:|
|D01|CONTIGS.db|D01/PROFILE.db|D01.bam|
|R01|CONTIGS.db|R01/PROFILE.db|R01.bam|
|R02|CONTIGS.db|R02/PROFILE.db|R02.bam|

### Optional columns

In addition to the required columns shown above, you can add as many columns as you like in your file. But two of these columns will be further processed during sanity check: `r1` and `r2`, whith the expectation that these columns will include information regarding the location of the raw FASTQ files. For each row the %(bams-and-profiles-txt)s file describes, the FASTQ files must be those that were used to generate the BAM files. Here is an example file with these two additional columns: 

|name|contigs_db_path|profile_db_path|bam_file_path|r1|r2|
|:--|:--:|:--:|:--:|:--:|:--:|
|D01|CONTIGS.db|D01/PROFILE.db|D01.bam|D01-R1.fastq|D01-R2.fastq|
|R01|CONTIGS.db|R01/PROFILE.db|R01.bam|R01-R1.fastq|R01-R2.fastq|
|R02|CONTIGS.db|R02/PROFILE.db|R02.bam|R02-R1.fastq|R02-R2.fastq|

Some programs, such as %(anvi-report-inversions)s, can process the `r1` and `r2` files.
