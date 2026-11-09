This program asks NCBI what it knows about a set of SRA run accessions and writes the answers into a %(sra-metadata-txt)s.

You do not have to run it. Any anvi'o workflow that downloads reads from the SRA builds this file itself the first time it needs one. This program exists for the two cases where doing it by hand is better: when the computer that runs your workflows cannot reach the internet, and when you would like to look over what anvi'o concluded before it acts on it.

### Starting from a list of accessions

{{ codestart }}
anvi-script-get-sra-metadata --accession-list SRA_accession_list.txt \
                             -o SRA-METADATA.txt
{{ codestop }}

where `SRA_accession_list.txt` is a file with one SRA run accession per line.

### Starting from a samples-txt

If you already have a %(samples-txt)s with an `sra_accession` column in it, point this program at that instead and it will pick the accessions out for you:

{{ codestart }}
anvi-script-get-sra-metadata --samples-txt %(samples-txt)s \
                             -o SRA-METADATA.txt
{{ codestop }}

### What it tells you

Along with the file itself you get a summary of what is in it: how many runs are paired-end short reads and how many are long reads, and how much disk space they would take up if you downloaded every one of them at once. That last number is a useful sanity check before setting `max_disk_gb` in your %(workflow-config)s.

You will also hear about anything that needs your attention — long-read runs whose sequencing chemistry could not be determined from NCBI's description, and single-end short-read runs, which the metagenomics workflow cannot process.

### Running it again

Running this program a second time on the same output file does not start over. Anvi'o reads what is already there, looks up only the accessions that are missing, and appends them. Rows you have corrected by hand are left exactly as you wrote them, which is what makes this file a reasonable place to fix NCBI metadata that is wrong.
