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

### Example output

```
anvi-script-reformat-fasta contigs.fa \
                           --simplify-names \
                           --prefix YYY \
                           --min-len 1000 \
                           --seq-type NT \
                           --overwrite-input
```

```
Input ........................................: contigs.fa
Output .......................................: (anvi'o will overwrite your input file)

WHAT WAS THERE
===============================================
Total num contigs ............................: 4,189
Total num nucleotides ........................: 35,766,167

WHAT WAS ASKED
===============================================
Simplify deflines? ...........................: Yes
Add prefix to sequence names? ................: Yes, add 'YYY'
Minimum length of contigs to keep ............: 1,000
Max % gaps allowed ...........................: 100.00%
Max num gaps allowed .........................: 1,000,000
Exclude specific sequences? ..................: No
Keep specific sequences? .....................: No
Enforce sequence type? .......................: Yes, enforce 'NT'

WHAT HAPPENED
===============================================
Contigs removed ..............................: 3,156 (75.34% of all)
Nucleotides removed ..........................: 6,121,239 (17.11% of all)
Nucleotides modified .........................: 161 (0.00045% of all)
Deflines simplified ..........................: True


* The contents of your input file have changed because you used the flag
`--overwrite-input`.

```

{:.warning}
Please use the flag `--overwrite-input` with extreme caution.
