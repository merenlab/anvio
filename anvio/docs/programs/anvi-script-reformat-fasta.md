This program **converts a %(fasta)s file to a %(contigs-fasta)s.** In other words, it reformats your FASTA formatted file to meet the conditions required of a %(contigs-fasta)s, which is able to be used by other anvi'o programs.

{{ codestart }}
anvi-script-reformat-fasta %(fasta)s \
                           -o %(contigs-fasta)s \
                           --simplify-names
{{ codestop }}

{:.notice}
If you use the flag `--report-file`, it will also create a TAB-delimited file for you to keep track of which defline in the new file corresponds to which defline in the original file.

{:.notice}
This program can work with compressed input FASTA files (i.e., the file name ends with a `.gz` extention) and will report a compressed output FASTA file (i.e., if the output file name ends with a `.gz` extension).

In addition to simplifying names, this program will allow you to do a combination of the operations that include,

* Add a prefix to sequnce names in a FASTA file,
* Remove sequences that are shorter than a specific length or only keep sequences that match to a specific length,
* Remove sequences if they contain more than a number of gap characters or exceed the precentage of gap characters you permit,
* Exclude sequences that match to a list of sequence IDs, or only keep those that match to a list of sequence IDs,
* Enforce a sequence type to replace any character with `N` for nucleotide sequences that are not A, C, T, or G, or replace any character with `X` for amino acid sequences if the character does not match any of the single-letter amino acid characters.

### Removing the short reads is important

If your FASTA file includes a lot of very short contigs, removing them may dramatically improve the performance of the generation and processing of your %(contigs-db)s. The example below runs the same command while also removing sequences that are shorter than 1,000 nts:

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
Max %% gaps allowed ...........................: 100.00%%
Max num gaps allowed .........................: 1,000,000
Exclude specific sequences? ..................: No
Keep specific sequences? .....................: No
Enforce sequence type? .......................: Yes, enforce 'NT'

WHAT HAPPENED
===============================================
Contigs removed ..............................: 3,156 (75.34%% of all)
Nucleotides removed ..........................: 6,121,239 (17.11%% of all)
Nucleotides modified .........................: 161 (0.00045%% of all)
Deflines simplified ..........................: True


* The contents of your input file have changed because you used the flag
`--overwrite-input`.

```

{:.warning}
Please use the flag `--overwrite-input` with extreme caution.
