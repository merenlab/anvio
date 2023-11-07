This program **reformats a %(bam-file)s file using a %(contig-rename-report-txt)s to a %(bam-file)s compatible with a re-formatted %(contigs-fasta)s.** If you had to run %(anvi-script-reformat-fasta)s to convert the contigs name of your FASTA file to make it compatible with anvi'o, and you already had some BAM files mapped to the original FASTA file, you will need to reformat those BAM files to make them compatible with the new FASTA file. This program will use the report file from %(anvi-script-reformat-fasta)s to do that for you.

{{ codestart }}
anvi-script-reformat-bam %(bam-file)s \
                         -l %(contig-rename-report-txt)s \
                         -o %(bam-file)s
{{ codestop }}

{:.notice}
This program is required only if you ran `anvi-script-reformat-fasta` with the `--simplify-names` flag, and will require the output file generated with the flag `--report-file`.

### Example output

```
anvi-script-reformat-bam SAMPLE.bam \
                         --list REPORT.txt \
                         --output SAMPLE-REFORMATTED.bam
```

```text
Input BAM file ...............................: SAMPLE.bam
Rename list ..................................: REPORT.txt
Overwrite output? ............................: True

WHAT WAS THERE
===============================================
Loaded REPORT.txt ............................: 3

WHAT WAS DONE
===============================================
Sequences in BAM file ........................: 3
Output BAM file created ......................: SAMPLE-REFORMATTED.bam
```
