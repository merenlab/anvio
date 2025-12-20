%(anvi-reorient-genomes)s aligns circular genomes in a %(fasta-txt)s to a chosen reference genome, rotates and/or reverse-complements them to synchronize their arbitrary circulizarization beginnings with the reference, and generates new 'reoriented %(fasta)s files' for downstream use. By doing so, it maximizes the gene order across genomes, so that downstream analyses that rely on synteny can have easier time.

It reports (1) alignment quality between each pair of genomes (coverage between the two genomes and very crappy approximate ANI that no one should trust), (2) what actions were taken to reach a consensus (reverse-complement, rotations, etc), and (3) a final determination of the level of trust for outcomes. The program also prints out simple dotplots to visualize the final alignment of each genome to the reference.

The default usage is simple:

{{ codestart }}
%(anvi-reorient-genomes)s --fasta-txt %(fasta-txt)s \
                           --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

In this use case, **the program will auto-select the longest genome as reference**, and it will then *intelligently* find and rotate to an optimal starting position that is conserved across most of the input genomes mentioned in the %(fasta-txt)s file. It is intelligent becasue this step will ensure that all genomes will start at a biologically meaningful position rather than an arbitrary coordinate. %(anvi-reorient-genomes)s does it by first aligning each genome to the reference, and then using a 10,000 nt sliding window to search for a region that is shared between all genomes in a greedy fashion (once it finds a region that works for all genomes, it stops).

{:.notice}
Please note that when a reference genome is chosen automatically, the program may end up rotating the automatically chosen reference genome to a more meaningful and conserved start position to maximize agreements across all genomes.

Another option is to **declare a reference genome** to be used to orient all others:

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --reference REF-GENOME-NAME \
                      --output-dir REORIENTED-FASTA-FILES/
                      --threads 8
{{ codestop }}

If you have a reference genome that you trust (downlodaded from a reliable source or manually set the starting point using the DnaA gene yourself or something to that effect), you can ask the program to use that as a reference. In which case the program will not be tinkering with the reference, and do its best to match every other genome to it.

For the most up-to-date list of parameters and their default values, please see the help menu on your terminal.

### Example terminal output

Outputs in the terminal will help you figure out what is going on, so please keep an eye on them. 

Here is an example output for a genome that nicely aligned to the reference:

```
REORIENTING HIMB1559 (10 OF 26 TOTAL)
===============================================
Orientation outcome ..........................: TRUSTWORTHY
Applied actions ..............................: reverse-complemented query (pass1 anchor '-'), rotated 793356 bp (pass1 snap), rotated 1388641 bp (pass2 snap), rotated 1 bp (snap-to-ref0 anchor), rotated 3246 bp (final snap to ref0)
Output FASTA .................................: /path/to/HIMB1559.fa
Final alignment strand .......................: +
Start in query ...............................: 61
Start in reference ...........................: 61
Query length .................................: 1,489,298
Reference length .............................: 1,538,266
Aligned length ...............................: 1,521,299
Query coverage by alignment ..................: 100.0%%
Reference coverage by alignment ..............: 98.9%%
Approx ANI to reference ......................: 100.0%%

                       After reorientation: HIMB1559 vs HIMB1506
      ┌────────────────────────────────────────────────────────────────────────┐
  1.5M┤ ▞▞ + strand                                                       ▄▄▘  │
      │                                                              ▄▄▞▀▀     │
      │                                                         ▗▄▞▀▀          │
  1.1M┤                                                    ▗▄▄▀▀▘              │
      │                                                ▄▄▀▀▘                   │
      │                                           ▄▄▞▀▀                        │
      │                                      ▗▄▞▀▀                             │
744.6K┤                                 ▗▄▄▀▀▘                                 │
      │                             ▄▄▀▀▘                                      │
      │                        ▄▄▞▀▀                                           │
372.3K┤                   ▗▄▞▀▀                                                │
      │              ▗▄▄▀▀▘                                                    │
      │          ▄▄▀▀▘                                                         │
      │     ▄▄▞▀▀                                                              │
     0┤▄▄▞▀▀                                                                   │
      └┬─────────────────┬─────────────────┬────────────────┬─────────────────┬┘
       0              384.6K            769.1K            1.2M             1.5M
Query start                         Reference start
```

The fact that 'Start in query' and 'Start in reference' is identical is great news, and shows that the program was able to set a start position that works perfectly. The high coverage alignment between the two is also great as it shows that you brough together genomes that are quite closely related. The approximate ANI says 100%%, and it is not to be trusted. For instance, the true ANI between these two genomes in this example is about 93%%. But since we are cutting corners using the `minimap2` aligned bits between the two, we really really overestimate the ANI value. It is useful, but please don't rely on the apprximate ANI reported here for anything serious.

In contrast to that example, here is a terminal output for an alignment attempt that did not work so well:

```
REORIENTING M (5 OF 5 TOTAL)
===============================================
Orientation outcome ..........................: NOT TRUSTWORTHY
Applied actions ..............................: reverse-complemented query (pass1 anchor '-'), rotated 2481312 bp (pass1 snap), rotated 1107466 bp (pass2 snap), rotated 1 bp (snap-to-ref0 anchor), rotated 1913130 bp (final snap to ref0)
Output FASTA .................................: /path/to/M.fa
Final alignment strand .......................: -
Start in query ...............................: 2,322,697
Start in reference ...........................: 518,631
Query length .................................: 3,035,044
Reference length .............................: 3,469,552
Aligned length ...............................: 231,049
Query coverage by alignment ..................: 7.6%%
Reference coverage by alignment ..............: 6.7%%
Approx ANI to reference ......................: 94.6%%

                              After reorientation: M vs D
      ┌────────────────────────────────────────────────────────────────────────┐
    3M┤ ▞▞ + strand ▝       ▖▖   ▖ ▘        ▘     ▗ ▗ ▝▖      ▗▗▗  ▗▄▞▘▗▗     ▝│
      │ ▞▞ - strand   ▝   ▝      ▗  ▘▝▘      ▞▀  ▗ ▘▝  ▖  ▖  ▘ ▗▄▞▀▘        ▘  │
      │ ▝ ▝   ▐ ▗  ▗   ▖  ▝  ▖ ▖     ▖▝    ▝   ▘ ▘  ▘▗ ▘ ▝ ▄▄▀▀▘▗ ▚▄           │
  2.3M┤▗                    ▖       ▗                 ▗▄▄▀▀ ▗  ▗        ▘ ▖    │
      │   ▖  ▖▖           ▖▗▗  ▝ ▝▘ ▗            ▗▗▄▞▀▘▞       ▝ ▘▘     ▝    ▖ │
      │ ▝  ▟  ▘    ▝▘▗    ▝   ▗ ▘▘           ▝▄▄▀▀▘▝    ▘         ▘   ▝    ▗   │
      │   ▖▝▝▌ ▝ ▖  ▛▘▄▝▖▘▝ ▘     ▘ ▝▖ ▘ ▗▄▄▀▀▘▗ ▗ ▄▗    ▌    ▀▀ ▞▘     ▚  ▝▗  │
  1.5M┤      ▝                 ▝   ▘ ▗▄▞▀▘               ▘                     │
      │▗▗  ▖    ▘  ▝ ▗▗      ▗ ▗ ▄▄▞▀▘        ▘   ▜ ▘▝          ▖ ▝  ▖ ▗    ▗  │
      │ ▖▖ ▀  ▝▐▖▐    ▝   ▗  ▄▄▀▀ ▗▝▖ ▘    ▝     ▗▗    ▗ ▗  ▌     ▝        ▗ ▗▘│
758.8K┤    ▘▀▝▘ ▘     ▝▘▗▄▞▀▀  ▖  ▗▌▘ ▝   ▗ ▝▝      ▌▗     ▄ ▖    ▝          ▝ │
      │            ▗▄▄▞▀▘ ▝▀ ▀   ▝           ▝        ▝                      ▝ │
      │▗ ▗▖    ▗▄▄▀▀    ▗ ▐   ▝  ▘▝ ▗▖    ▗▝   ▘ ▖▝ ▄▖ ▗  ▝     ▞▘▗▗  ▘  ▗▖   ▘│
      │ ▗  ▗▄▞▀▀ ▘ ▝ ▝ ▝    ▖▛▚ ▖          ▖ ▞▘▖ ▖ ▖  ▀   ▌ ▝ ▀▗▘▖▖         ▘ ▖│
     0┤▄▄▞▀▘▘ ▝ ▗    ▘▗      ▛ ▘ ▖  ▀▝▘   ▝ ▖ ▝▝  ▗              ▖     ▝▖ ▖▘▖  │
      └┬─────────────────┬─────────────────┬────────────────┬─────────────────┬┘
       0              867.4K             1.7M             2.6M             3.5M
Query start                         Reference start
```

The program will give you something to read and think about at the very end:

```
(...)

FINAL REPORT
===============================================
Your genome reorientation task considered 5 genomes and 1 reference. Some
outcomes are not trustworthy (low alignment coverage can lead to unreliable
orientation). Please review the dotplots above and the summary below to decide
which FASTA files to use downstream.

Trustworthy ..................................: 1
    - D -> /path/to/D.fa

Somewhat OK ..................................: 1
    - G -> /path/to/G.fa

Not trustworthy ..............................: 4
    - C -> /path/to/C.fa
    - B -> /path/to/B.fa
    - F -> /path/to/F.fa
    - M -> /path/to/M.fa

Failed .......................................: 0

✓ reorient_genomes.py took 0:00:18.481794
```


### What the program does

Here is a more detailed description of what is going on behind the scenes when you press ENTER:

1. **Parse inputs and pick a reference**. Reads %(fasta-txt)s, and if `--reference` is not set, the program picks the longest genome. All FASTAs are sanity-checked (existence, FASTA format, single contig).

2. **Find optimal starting position (when reference is auto-selected)**. If you don't specify `--reference`, all genomes are considered to find the most conserved position in the reference genome. It aligns each genome to the reference, builds a coverage map showing which positions are covered by how many genomes, and identifies the position with maximum coverage (stopping early if it finds 100%% coverage). The reference is then rotated to start at this optimal position, ensuring all genomes will start at a biologically meaningful, conserved region.

3. **Initial alignment (reference vs query)**. Runs `minimap2` (preset `asm5`, with as many threads as you asked for using `--threads`) to align each query to the reference and identifies a primary anchor near reference position 0. If the anchor is on `-` strand, the query is reverse-complemented. The query is rotated so that reference position 0 maps onto the query (first snap).

4. **Second alignment and snap**. Re-aligns the rotated query, finds the primary anchor with the smallest reference start, and rotates again to bring reference 0 onto the query (second snap).

5. **Snap-to-zero with a ref0-focused anchor**. Aligns once again, picks the primary anchor closest to reference position 0, rotates, aligns once more, and applies a final snap so that reference position 0 maps to query position 0.

6. **Iterative correction for perfect alignment**. After the final snap, the program checks if the alignment truly starts at position 0 in both the query and reference. If not (e.g., due to `minimap2` soft-clipping divergent regions), it calculates the necessary rotation, applies it, and re-aligns. This iterates up to 5 times or until the genomes are perfectly aligned at position 0.

7. **Write outputs**. Copies the reference %(fasta)s to the output directory (potentially rotated if auto-selected). Writes each reoriented query %(fasta)s under the same name (and using the original extension).

8. **Report per-genome stats and dotplots**. For each genome, the program reports the orientation outcome (whether it was `TRUSTWORTHY`, `SOMEWHAT OK`, or `NOT TRUSTWORTHY` based on average alignment coverage of query and reference where less than 50%% is assumed trustworthy, betwen 50%% and 90%% is assumed somewhat OK, and over 90%% is assumed trustworthy), many other statistics, and a dotplot of the final alignment to visualize orientation quality.

9. **Final report**. Summarizes the number of genomes in different trust categories along their output FASTA paths for the user to decide which outputs are safe for downstream analyses.


### Tips, caveats, and runtime

* **Interpreting "Start in query" and "Start in reference"**: These numbers show where `minimap2` primary alignments begin. **When both values are equal** (e.g., both 61), your genomes are correctly aligned at the same biological position. The small offset is just natural sequence variation in the first few dozen base pairs since `minimap2` soft-clips divergent regions at the beginning. In fact, Meren was so concerned about it, he triple-checked this, and discovered that despite differences in offsets across genomes, the amino acid sequences of the first genes in every aligned genome was identical even though their 'start' positions in individual genomes differed. So that's that. But **when the values differ** (e.g., query=0, reference=106), the alignment may have issues and the genome might not be well-oriented.

* **Interpreting trust labels**: Low coverage (either genome covered <50%%) yields `NOT TRUSTWORTHY`. This often means the genome is too divergent or structurally different to reorient confidently. Check the dotplot and stats before using such outputs. Genomes with `TRUSTWORTHY` labels and matching start positions are safe for downstream comparative analyses.

* **Auto-selected reference and optimal start**: When you don't specify `--reference`, the program finds and rotates to a conserved position across your genomes. This is especially useful for sets of closely related genomes where you want them all to start at a biologically meaningful position (like the beginning of a conserved gene). If you want a specific genome or starting position, use `--reference` to override this behavior.

* **Approximate ANI**: This is quite a garbage number at this point, and it is no one's fault. When `minimap2` emits `dv`, ANI is computed as `(1 - dv) x 100`. Otherwise it falls back to `nmatch/alen × 100`. So treat it as a quick proxy, not an accurate ANI calculation. It is not the purpose here.

* **Circular ambiguity**: Circular genomes can align equally well at different offsets. The program applies multiple snaps and iterative corrections to align reference position 0 to query position 0, but in highly repetitive cases the true biological origin may still be ambiguous.

* **Runtime**: Each query triggers several `minimap2` runs and `seqkit` rotations, so it can take some time to converge. Then the optimal start finding step (when auto-selecting reference) adds an initial survey phase. But it takes no more than 30 seconds on a laptop computer to align 30 SAR11 genomes, so it is not that bad all things considered.

{:.notice}
YES. We will make it also work for MAGs and fragmented genomes. One step at a time.
