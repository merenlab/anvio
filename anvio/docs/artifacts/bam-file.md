A BAM file contains **already aligned sequence data.** However, it is written in binary to save space (so it will look like jibberish if you open it). 

BAM files (and their text file cousin SAM files) are often used in 'omics analysis and are described in more detail in [this file](https://samtools.github.io/hts-specs/SAMv1.pdf), written by the developers of samtools. 

If your BAM file is not indexed, it is actually a %(raw-bam-file)s and you can run %(anvi-init-bam)s to turn it into a BAM file. You can tell if your BAM file is indexed if in the same folder as your `XXXX.bam` file, there is another file with the same name called `XXXX.bam.bai`.

As of now, no anvi'o programs will output results in BAM format, so you'll primary use BAM files to import sequence data into anvi'o. For example, in %(anvi-profile)s (which generates a %(profile-db)s), your BAM file is expected to contain the aligned short reads from your samples. 

