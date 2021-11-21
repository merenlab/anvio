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
