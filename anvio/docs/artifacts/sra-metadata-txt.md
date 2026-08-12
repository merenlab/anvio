A **TAB-delimited** file describing a set of NCBI SRA sequencing runs: what kind of reads each one holds, and how much room it is going to take up.

Anvi'o builds this file for you the first time it needs it, and you are welcome to correct it.

### Why it exists

Before a workflow can download reads from the SRA, it has to know some things about them that an accession alone does not say. Whether a run is paired-end decides how many FASTQ files will come out of it. Whether it came off a short-read or a long-read instrument decides which assembler, mapper, and quality-control steps it goes through. And how big it is decides how many runs can be on disk at the same time when you have set a disk budget. All of this has to be settled before any work begins, since it determines the shape of the workflow itself.

So anvi'o asks NCBI once, writes the answers here, and reads this file every time after that.

### What is in it

|accession|read_type|lr_technology|platform|model|library_layout|spots|bases|size_mb|source|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|ERR6450080|PE-SR||ILLUMINA|Illumina MiSeq|PAIRED|132807|71988518|38|ncbi|
|SRR39242466|LR|ont|OXFORD_NANOPORE|MinION|SINGLE|781|466614|0|ncbi|
|SRR11951439|LR|pb-hifi|PACBIO_SMRT|Sequel II|SINGLE|225118|1367510760|800|user|

* `read_type` is anvi'o's own reading of the run: `PE-SR` for paired-end short reads, `LR` for long reads, and `SE-SR` for single-end short reads (which the metagenomics workflow cannot process).
* `lr_technology` is the long-read technology, when anvi'o could determine it. Nanopore runs are unambiguous. PacBio runs are not always: an instrument like the Sequel II was used for both CLR and HiFi sequencing, so anvi'o leaves this blank rather than guessing, and tells you which accessions it could not settle.
* `platform`, `model`, `library_layout`, `spots`, `bases` and `size_mb` are what NCBI reports, kept here so you can see what anvi'o based its reading on.
* `source` is `ncbi` for a row anvi'o filled in and `user` for one you have edited.

### Correcting it

Much of what NCBI knows about a run is whatever the person who submitted it typed in, so it is occasionally wrong. When it is, edit this file: change the `read_type` of a run NCBI mislabelled, or fill in the `lr_technology` of a PacBio run you happen to know the chemistry of. Set `source` to `user` on any row you touch, and anvi'o will neither overwrite it nor pester you about it again.

Anvi'o only ever looks up accessions that are *missing* from this file, so nothing you have written is at risk.

### Making one ahead of time

Anvi'o will build this file on its own when a workflow first needs it. If the machine that runs your workflows has no internet access, or you would simply like to see what anvi'o thinks before it acts on it, you can build it in advance with %(anvi-script-get-sra-metadata)s.
