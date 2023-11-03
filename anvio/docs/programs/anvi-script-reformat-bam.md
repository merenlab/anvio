This program **reformats a %(bam-file)s file to a %(bam-file)s compatible with a re-formatted %(contigs-fasta).** If you had to run `anvi-script-reformat-fasta` to convert the contigs name of your FASTA file to make it compatible with Anvi'o, and you already had some BAM files mapped to the original FASTA file, you will need to reformat those BAM files to make them compatible with the new FASTA file. This program will use the report file from `anvi-script-reformat-fasta` to do that for you.


{{ codestart }}
anvi-script-reformat-fasta %(bam-file)s \
                           -l %(report-file)s \
                           -o %(bam-file)s
{{ codestop }}

{:.notice}
This program is required only if you ran `anvi-script-reformat-fasta` with the `--simplify-names` switch, and will require the output file generadet with the flag `--report-file`.



### Example output

```
anvi-script-reformat-bam  sample1_raw.bam \
                           --list REPORT.txt \
                           --output sample1.bam \
                           --overwrite-output
```

```text
Input BAM file ...............................: test/input/raw.bam
Rename list ..................................: test/REPORT_FILE.txt
Overwrite output? ............................: False

WHAT WAS THERE
===============================================
Loaded test/REPORT_FILE.txt ..................: 3

WHAT WAS DONE
===============================================
Sequences in BAM file ........................: 3
Output BAM file created ......................: test/output/final.bam
```

{:.warning}
Please use the flag `--overwrite-output` with extreme caution, as it will overwrite a pre-existing BAM file with the same name as the output file.
