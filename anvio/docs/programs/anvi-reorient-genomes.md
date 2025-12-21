%(anvi-reorient-genomes)s aligns genomes listed in a %(fasta-txt)s to a chosen reference genome and reorients them to match the reference coordinate system. For **circular genomes** (complete bacterial genomes, plasmids, viral genomes), it rotates and/or reverse-complements them to synchronize their arbitrary circularization beginnings with the reference. For **fragmented genomes** (MAGs, draft assemblies), it orders and orients contigs to match the reference genome's coordinate system. By doing so, it maximizes gene order conservation across genomes, so that downstream analyses that rely on synteny can have an easier time.

It reports (1) alignment quality between each pair of genomes (coverage between the two genomes and very crappy approximate ANI that no one should trust), (2) what actions were taken to reach a consensus (reverse-complement, rotations, contig ordering, etc), and (3) a final determination of the level of trust one should have for outcomes. The program also generates synteny ribbon plots to visualize the alignment patterns between each genome and the reference before and after reorientation.

The default usage is simple, which will simply instruct anvi'o to **take care of everything *de novo***,

{{ codestart }}
%(anvi-reorient-genomes)s --fasta-txt %(fasta-txt)s \
                           --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

Another option is to tell anvi'o to **use a specific genome as reference** to to orient all others:

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --reference REF-GENOME-NAME \
                      --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

Or one can ask anvi'o to **first use the DnaA gene to orient the reference** before orienting all others:

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --use-dnaa-for-reference-orientation \
                      --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

### De novo identification of reference orientation

{:.notice}
TL;DR: Best for any set of highly similar genomes - whether circular (viral, plasmids, complete bacterial genomes) or fragmented (MAGs, draft assemblies). For circular genomes, uses secondary alignments from `minimap2` to find the most conserved position across all genomes in a data-driven manner. For fragmented genomes, orders and orients contigs based on their alignment to the reference.

If the user does not explicitly mention a reference genome, **the program will auto-select the genome with the fewest contigs as reference** (ties broken by longest total length), preferring complete circular genomes over fragmented ones. It will then rotate the reference genome to an optimal starting position that is conserved across most of the input genomes mentioned in the %(fasta-txt)s file in a *mindful* fashion. It is 'mindful' because %(anvi-reorient-genomes)s does this step by first aligning each genome to the reference **using all possible good alignments (primary + secondary alignments)** to get a true picture of conserved regions across all genomes, and then using a 1,000 bp sliding window to search for a region that is shared between all genomes in a greedy fashion (once it finds a region that works for all genomes, it stops).

{:.notice}
Please note that when a reference genome is chosen automatically, the program will likely end up rotating the automatically chosen reference genome to a more meaningful and conserved start position to maximize agreements across all genomes.

If you have a reference genome that you trust (downloaded from a reliable source or manually circularized), you can ask the program to use that as a reference. In which case the program will not be tinkering with the reference, and do its best to match every other genome to it.

### Using DnaA gene for biologically meaningful reference orientation

{:.notice}
TL;DR: Best for any highly similar set of **bacterial genomes**. If DnaA is not found in the reference, the program will issue a warning and proceed without rotating the reference.

For **bacterial genomes**, you can use the `--use-dnaa-for-reference-orientation` flag to orient the reference genome based on the DnaA gene, which typically marks the origin of replication:

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --use-dnaa-for-reference-orientation \
                      --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

This approach will use `prodigal` to call genes in the reference genome, use `hmmsearch` with the Bac_DnaA_C HMM profile from Pfam to identify the DnaA gene, rotate the reference to start at the DnaA gene position, then align all other genomes to this DnaA-oriented reference.

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

1. **Parse inputs and pick a reference**. Reads %(fasta-txt)s, and if `--reference` is not set, the program picks the genome with the fewest contigs (ties broken by longest total length). All FASTAs are sanity-checked (existence, FASTA format). The reference must be a single-contig circular genome.

2. **Determine reference orientation**. The program uses one of three strategies to orient the reference genome:
   - **DnaA-based orientation** (if `--use-dnaa-for-reference-orientation` is set): Calls genes with `prodigal`, searches for the DnaA gene using `hmmsearch` with the Bac_DnaA_C HMM profile, and rotates the reference to start at the DnaA gene position. This provides biologically meaningful orientation for bacterial genomes.
   - **De novo optimal position** (if reference is auto-selected without DnaA flag): Aligns each genome to the reference using minimap2 with `--secondary=yes -N 100 -p 0.5` to capture **all possible good alignments** (not just the best one). Builds a coverage map using 1,000 bp bins showing which positions are covered by alignments from each genome. Identifies the position with maximum coverage across all genomes (stopping early if it finds 100%% coverage). The reference is then rotated to start at this optimal position, ensuring all genomes will start at a conserved region that is genuinely shared across the dataset.
   - **User-specified reference** (if `--reference` is set): Uses the reference genome as-is without rotation.

**For circular genomes (single-contig):**

3. **Initial alignment (reference vs query)**. Runs `minimap2` (preset `asm5`, with as many threads as you asked for using `--threads`) to align each query to the reference and identifies a primary anchor near reference position 0. If the anchor is on `-` strand, the query is reverse-complemented. The query is rotated so that reference position 0 maps onto the query (first snap).

4. **Second alignment and snap**. Re-aligns the rotated query, finds the primary anchor with the smallest reference start, and rotates again to bring reference 0 onto the query (second snap).

5. **Snap-to-zero with a ref0-focused anchor**. Aligns once again, picks the primary anchor closest to reference position 0, rotates, aligns once more, and applies a final snap so that reference position 0 maps to query position 0.

6. **Iterative correction for perfect alignment**. After the final snap, the program checks if the alignment truly starts at position 0 in both the query and reference. If not (e.g., due to `minimap2` soft-clipping divergent regions), it calculates the necessary rotation, applies it, and re-aligns. This iterates up to 5 times or until the genomes are perfectly aligned at position 0.

**For fragmented genomes (multi-contig MAGs or draft assemblies):**

3. **Contig filtering**. Contigs shorter than `--min-contig-length` (default: 1,000 bp) are excluded from processing.

4. **Individual contig alignment**. Each contig is independently aligned to the reference using `minimap2`. The program detects contigs that span the circular boundary of the reference (wrap-around contigs) and splits them into separate parts for proper ordering.

5. **Contig ordering and orientation**. Contigs are sorted by their alignment position on the reference genome. Contigs aligned to the reverse strand are reverse-complemented to match the reference orientation.

6. **Output generation**. By default, contigs are written as separate sequences in the output FASTA, ordered and oriented to match the reference. If `--scaffold-fragmented` is used, contigs are concatenated into a single sequence with N-padding representing gaps based on the reference genome distances.

**For all genomes:**

7. **Write outputs**. Copies the reference %(fasta)s to the output directory (potentially rotated if auto-selected). Writes each reoriented query %(fasta)s under the same name (and using the original extension).

8. **Report per-genome stats and alignment plots**. For each genome, the program reports the orientation outcome (whether it was `TRUSTWORTHY`, `SOMEWHAT OK`, or `NOT TRUSTWORTHY` based on alignment coverage), many other statistics, and synteny ribbon plots showing the alignment patterns before and after reorientation to visualize orientation quality.

9. **Final report**. Summarizes the number of genomes in different trust categories along their output FASTA paths for the user to decide which outputs are safe for downstream analyses.


### Tips, caveats, and runtime

* **Interpreting "Start in query" and "Start in reference" (circular genomes only)**: These numbers show where `minimap2` primary alignments begin. **When both values are equal** (e.g., both 61), your genomes are correctly aligned at the same biological position. The small offset is just natural sequence variation in the first few dozen base pairs since `minimap2` soft-clips divergent regions at the beginning. In fact, Meren was so concerned about it, he triple-checked this, and discovered that despite differences in offsets across genomes, the amino acid sequences of the first genes in every aligned genome was identical even though their 'start' positions in individual genomes differed. So that's that. But **when the values differ** (e.g., query=0, reference=106), the alignment may have issues and the genome might not be well-oriented.

* **Interpreting trust labels**: For circular genomes, low coverage (either genome covered <50%%) yields `NOT TRUSTWORTHY`. For fragmented genomes, the label considers both reference coverage and the percentage of contigs that aligned. This often means the genome is too divergent or structurally different to reorient confidently. Check the alignment plots and stats before using such outputs. Genomes with `TRUSTWORTHY` labels are safe for downstream comparative analyses.

* **Understanding alignment visualizations**: The synteny ribbon plots show alignment blocks as ribbons connecting the reference (top) to the query (bottom). Green ribbons indicate forward-strand alignments, red ribbons indicate reverse-strand alignments. For fragmented genomes in the "before" plot, each contig appears as a separate segment with gaps between them. After reorientation, contigs are ordered and oriented to match the reference. Contig segments in the query genome are color-coded: white for forward-strand contigs, red for reverse-complemented contigs, and green for contigs that were split due to wrapping around the circular reference boundary.

* **Scaffolding fragmented genomes**: By default, fragmented genomes are written with contigs as separate sequences (ordered and oriented). Use `--scaffold-fragmented` to concatenate them into a single sequence with N-padding representing gaps. While this maximizes gene synteny for comparative analyses, **do not submit N-scaffolded MAGs as complete genomes** to public databases. The scaffolding is reference-based and may not represent the true genomic structure.

* **Auto-selected reference and optimal start**: When you don't specify `--reference`, the program picks the genome with the fewest contigs (ties broken by longest total length), preferring complete circular genomes over fragmented ones. It then rotates the reference to a conserved position across your dataset using secondary alignments to capture all conserved regions (not just the single best alignment). This is especially useful for sets of closely related genomes where you want them all to start at a biologically meaningful position. For bacterial genomes, consider using `--use-dnaa-for-reference-orientation` for even more consistent results based on the replication origin. If you want a specific genome or starting position, use `--reference` to override this behavior.

* **DnaA-based orientation benefits**: When working with bacterial genomes, the `--use-dnaa-for-reference-orientation` flag typically produces highly consistent alignments (e.g., all genomes starting within a few base pairs of each other) because it uses biological knowledge (the replication origin) rather than purely sequence-based heuristics. This can be particularly valuable for downstream synteny analyses or when comparing gene order across closely related strains.

* **Approximate ANI**: This is quite a garbage number at this point, and it is no one's fault. When `minimap2` emits `dv`, ANI is computed as `(1 - dv) x 100`. Otherwise it falls back to `nmatch/alen x 100`. So treat it as a quick proxy, not an accurate ANI calculation. It is not the purpose here.

* **Circular ambiguity**: Circular genomes can align equally well at different offsets. The program applies multiple snaps and iterative corrections to align reference position 0 to query position 0, but in highly repetitive cases the true biological origin may still be ambiguous.

* **Visualization options**: By default, alignment plots are generated to help you assess the quality of reorientation. Use `--skip-visualizing-alignments` to disable plotting for faster processing when you only need the FASTA files. Customize plot dimensions with `--plot-width` and `--plot-height` (note: widths below 100 characters may not display properly).

* **Runtime**: Each circular genome triggers several `minimap2` runs and `seqkit` rotations. Fragmented genomes are faster as each contig is aligned only once. The optimal start finding step (when auto-selecting reference) adds an initial survey phase that uses secondary alignments for more accurate conserved region detection. The `--use-dnaa-for-reference-orientation` flag adds gene calling and HMM search overhead (a few seconds for a typical bacterial genome). Overall, it takes no more than 30 seconds on a laptop computer to reorient 30 SAR11 genomes using the de novo approach, and slightly longer with DnaA-based orientation.
