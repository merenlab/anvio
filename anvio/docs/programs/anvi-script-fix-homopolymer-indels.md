This program takes an input %(fasta)s file with one or more sequences, then **corrects INDELs associated with homopolymer regions given a reference %(fasta)s file**, and reports edited sequences as a new %(fasta)s file.

{:.warning}
You must be extremely careful with this program since it reports edited sequences.

We developed this tool to ameliorate the large number of INDEL errors Oxford Nanopore Technology yields. When there is a high-quality reference genome, this program can align a set of input sequences to the reference, and when it sees something like this in the alignment:

```
Input sequence ......: ... CGAAAAACG ...
Reference sequence ..: ... CGAAA--CG ...
```

It can correct the input sequence to look like this, since this would indicate that the additional `A` nucleotides after `AAA` were likely due to errors from the long-read sequencing:

```
Input sequence ......: ... CGAAACG ...
```

Similarly, if the program sees a case like this:

```
Input sequence ......: ... CGAAA--CG ...
Reference sequence ..: ... CGAAAAACG ...
```

The program would correct the input sequence this way, assuming that the lacking `A` nucleotides were likely due to errors from long-read sequencing:


```
Input sequence ......: ... CGAAAAACG ...
```


## Homopolymer length

The parameter `--min-homopolymer-length` helps the program to determine what to call a homopolymer region. Please read the help menu for this parameter carefully since it is not exactly intuitive.

The value `3` would be a stringent setting. But you may want to lower it to `2` if you promise to evaluate your output carefully.

The script includes some pre-aligned test sequences for you to see how `--min-homopolymer-length` influences things. For instance, here is the output for `--min-homopolymer-length 2`:

```
anvi-script-fix-homopolymer-indels --test-run \
                                   --min-homopolymer-length 2

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCG-AAATCGATCGATCG

GAPS BEFORE REPEATS OF "A"
===============================================
Query ........................................: AAA
Reference ....................................: -AA
Resolution ...................................: {'action': 'DEL', 'positions': [12], 'nt': 'A'}

Edited query sequence ........................: ATCGATCGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAA-TCGATCGATCG

GAPS AFTER REPEATS OF "A"
===============================================
Query ........................................: AAA
Reference ....................................: AA-
Resolution ...................................: {'action': 'DEL', 'positions': [14], 'nt': 'A'}

Edited query sequence ........................: ATCGATCGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAA-TCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAAAATCGATCGATCG

GAPS AFTER REPEATS OF "A"
===============================================
Query ........................................: AA-
Reference ....................................: AAA
Resolution ...................................: {'action': 'INS', 'positions': [15], 'nt': 'A'}

Edited query sequence ........................: ATCGATCGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: A-C-AT-GATCG-AAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAAAATCGATCGATCG

GAPS BEFORE REPEATS OF "A"
===============================================
Query ........................................: -AA
Reference ....................................: AAA
Resolution ...................................: {'action': 'INS', 'positions': [9], 'nt': 'A'}

Edited query sequence ........................: ACATGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAA-TCGATCGATCG

GAPS AFTER REPEATS OF "A"
===============================================
Query ........................................: AAA
Reference ....................................: AA-
Resolution ...................................: {'action': 'DEL', 'positions': [14], 'nt': 'A'}

Edited query sequence ........................: ATCGATCGATCGAATCGATCGATCG
Query sequence edited correctly? .............: Yes
```

The output for the same input sequences with `--min-homopolymer-length 3`:


```
anvi-script-fix-homopolymer-indels --test-run \
                                   --min-homopolymer-length 3


* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCG-AAATCGATCGATCG

GAPS BEFORE REPEATS OF "A"
===============================================
Query ........................................: AAAA
Reference ....................................: -AAA
Resolution ...................................: {'action': 'DEL', 'positions': [12], 'nt': 'A'}

Edited query sequence ........................: ATCGATCGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAA-TCGATCGATCG
Edited query sequence ........................: ATCGATCGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: No

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAA-TCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAAAATCGATCGATCG

GAPS AFTER REPEATS OF "A"
===============================================
Query ........................................: AAA-
Reference ....................................: AAAA
Resolution ...................................: {'action': 'INS', 'positions': [15], 'nt': 'A'}

Edited query sequence ........................: ATCGATCGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: A-C-AT-GATCG-AAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAAAATCGATCGATCG

GAPS BEFORE REPEATS OF "A"
===============================================
Query ........................................: -AAA
Reference ....................................: AAAA
Resolution ...................................: {'action': 'INS', 'positions': [9], 'nt': 'A'}

Edited query sequence ........................: ACATGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: Yes

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAA-TCGATCGATCG
Edited query sequence ........................: ATCGATCGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: No
```

The output for the same input sequences with `--min-homopolymer-length 4`:

```
anvi-script-fix-homopolymer-indels --test-run \
                                   --min-homopolymer-length 4

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCG-AAATCGATCGATCG
Edited query sequence ........................: ATCGATCGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: No

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAA-TCGATCGATCG
Edited query sequence ........................: ATCGATCGATCGAAAATCGATCGATCG
Query sequence edited correctly? .............: No

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAA-TCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAAAATCGATCGATCG
Edited query sequence ........................: ATCGATCGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: No

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: A-C-AT-GATCG-AAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAAAATCGATCGATCG
Edited query sequence ........................: ACATGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: No

* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Query sequence ...............................: ATCGATCGATCGAAATCGATCGATCG
Reference sequence ...........................: ATCGATCGATCGAA-TCGATCGATCG
Edited query sequence ........................: ATCGATCGATCGAAATCGATCGATCG
Query sequence edited correctly? .............: No
```

## Tips and Warnings

As you correct your input sequences one round, the BLAST may produce new homopolymers. So you may want to re-run the tool by turning the output sequence into an input sequence. For instance, We had an *Akkermensia* genome reconstructed using long-read sequencing that matched to a gold-standard genome on NCBI.

Running the script the first time this way,

``` bash
anvi-script-fix-homopolymer-indels --input Akkermansia_minION.fasta \
                                   --reference Akkermansia_REFERENCE.fasta \
                                   --output Akkermansia_minION_CORRECTED.fasta
```

Produced the following output:

```
OVERALL & PER-SEQUENCE STATS
===============================================
Num input sequences ..........................: 1
Num homopolymers associated with INDELs ......: 529
Num actions ..................................: 597
Num insertions ...............................: 292
Num deletions ................................: 305

* contig_332_pilon
    - Homopolymers: 529
    - Insertions: 292
    - Deletions: 305

Corrected output FASTA .......................: Akkermansia_minION_CORRECTED.fasta
```

Then copying the output file as the input file,

```
cp Akkermansia_minION_CORRECTED.fasta Akkermansia_minION.fasta
```

And re-running the script the same way multiple times gave the following outputs:

``` bash

# (first round)

* contig_332_pilon
    - Homopolymers associated with INDELs: 29
    - Insertions: 21
    - Deletions: 10

# (second round)

* contig_332_pilon
    - Homopolymers associated with INDELs: 17
    - Insertions: 9
    - Deletions: 13

# (third round)

* contig_332_pilon
    - Homopolymers associated with INDELs: 4
    - Insertions: 4
    - Deletions: 0

# (fourth round)

* contig_332_pilon
    - Homopolymers associated with INDELs: 0
    - Insertions: 0
    - Deletions: 0
```

At the end, there were no more homopolymers associated with INDELs.

**Please consider the following points**:

* If the input and reference genomes are not closely related enough (i.e., expected ANI if there were no sequencing errors > 98%-99%), this process may yield very incorrect outcomes. But it should work great for genomes reconstructed from the same culture.

* The iterative improvement of a given input genome may reach to a 'back-and-forth' situation where there is no overall improvement, but the homopolymers associated with INDELs do not reach to 0. This happens when there are repeats in the reference genome that are identical to each other expect the number of nucleotides in homopolymers.

* You can always add `--verbose` to your command to see every single case that is considered, and resolution anvi'o reached.

* Under all circumstances, it is important to double check your results, and make sure you keep in mind that anything you see outstanding in your downstream analyses may be due to this step.