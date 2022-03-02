This program **converts a %(fasta)s file to a %(contigs-fasta)s.** In other words, it reformats your FASTA formatted file to meet the conditions required of a %(contigs-fasta)s, which is able to be used by other anvi'o programs.

{{ codestart }}
anvi-script-reformat-fasta %(fasta)s \
                           -o %(contigs-fasta)s \
                           --simplify-names
{{ codestop }}

{:.notice}
If you use the flag *--report-file*, it will also create a TAB-delimited file for you to keep track of which defline in the new file corresponds to which defline in the original file.

### Removing the short reads

Removing short contigs from a FASTA file will improve the performance of the %(contigs-db)s later. The example below runs the same command while also removing sequences that are shorter than 1,000 nts:

{{ codestart }}
anvi-script-reformat-fasta %(fasta)s \
                           -o %(contigs-fasta)s \
                           -l 1000 \
                           --simplify-names
{{ codestop }}

