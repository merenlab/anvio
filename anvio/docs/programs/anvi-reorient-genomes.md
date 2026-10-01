%(anvi-reorient-genomes)s aligns genomes listed in a %(fasta-txt)s to a chosen reference genome and reorients them to match the reference coordinate system. It is one of the most sophisticated individual tools in the anvi'o ecosystem, and while it is quite powerful, like most other anvi'o programs, its power depends on your understanding of what you are doing with it. Hence the very long help material below which we hope you will go through before using it.

For **circular genomes** (complete bacterial genomes, plasmids, viral genomes), %(anvi-reorient-genomes)s will rotate and/or reverse-complement them to synchronize their arbitrary circularization beginnings with the reference.

For **fragmented genomes** (SAGs, MAGs, draft assemblies of isolates), it will order and orient contigs in a query genome to match the coordinate system of a 'reference genome'. %(anvi-reorient-genomes)s will cut the occasional contig that runs past the beginning or the end of the reference and therefore has no single place to go. By doing so, it will maximize the gene order conservation across genomes, so that downstream analyses that rely on synteny can have an easier time, at the expense of literally increasing the number of contigs in your query genome up to two. It will also circularly rotate your contigs when necessary to accommodate reporting errors assemblers do (we knew a lot about the r[eporting errors assemblers make](https://doi.org/10.1038/s41587-025-02971-8), but while implementing this tool we discovered yet another one that we had never considered or observed before .. and more on this is below under '*Circularly permuted contigs*').

%(anvi-reorient-genomes)s will report (1) alignment quality between each pair of genomes (coverage between the two genomes and very crappy approximate ANI that no one should trust), (2) what actions were taken to reach a consensus (reverse-complement, rotations, contig ordering, etc), and (3) a final determination of the level of trust one should have for outcomes.

The program also generates synteny ribbon plots to visualize the alignment patterns between each genome and the reference before and after reorientation which you should certainly study carefully before moving on with your 'reoriented' genomes.

### Default usage

The default usage is simple, which will simply instruct anvi'o to **take care of everything *de novo*** (during which anvi'o will choose the longest single-contig genome as reference for you),

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

You can also ask anvi'o to **first use the DnaA gene to orient the reference** before orienting all others (which in our experience works best if you have bacterial / archaeal genomes and a good reference, and you should certainly consider doing this):

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --use-dnaa-for-reference-orientation \
                      --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

There are many other options in the `--help` menu of the program, but please keep reading since there are a few things you may want to know before start using it.

### Critical considerations of 'reference'

Well. How to treat your reference is the **single most important thing to know about this program**, and there are a few considerations to keep in mind when it comes to thinking about the reference.

Does your reference have to be a single-contig? Short answer: yes, it must be. But it doesn't necessarily be a complete, circular genome. It has to be that way *only* if you will let/ask anvi'o to **rotate** it. Here are some thoughts to clarify the situation a bit better:

* **The reference must always be a single contig.** Orienting contigs, ordering them, and scaffolding them all require a single, continuous coordinate system different contigs from a fragmented reference cannot offer. Anvi'o will stop with an error if the reference (whether you passed it explicitly or it is auto-selected by anvi'o for you) contains more than one contig. The *only way* around this error is to use `--use-auto-reference-as-is`, which comes with other caveats you must consdier (explained below).

* **The reference must *additionally* be complete and circular if anvi'o IS going to ROTATE it.** Rotation of the reference happens in two cases. If the reference is auto-selected (anvi'o rotates it to a conserved start position that all input genomes share), or when you ask anvi'o to use DnaA gene to orient the reference. But you have to keep in mind that rotating a linear sequence or a fragment of a chromosome is a biologically meaningless operation (since a fragment does not capture the entirety of a circular chromosome, so there is no coordinate to 'rotate around'). But anvi'o cannot know this since it cannot know what kind of sequences you are working with. **This requirement does not apply if nothing is being rotated**, so if you have sequences that shouldn't be rotated, you either name a reference explicitly `--reference` or tell anvi'o to use the reference it automatically choses as is with `--use-auto-reference-as-is`.

So, does your reference have to be a complete circular genome? Not necessarily. But the answer really depends on what you are planning to do with it. Here are some more ideas for you to think about before running this tool:

* **If you are working with complete bacterial or archaeal chromosomes** (or circular plasmids or viral genomes) and you care about where each genome *starts* (so that the origin of replication .. or another evolutionary meaningful and/or conserved region lines up nicely across your collection), then yes, rotation is the whole point, and both your reference and your parameters matter a great deal. Let anvi'o auto-select and rotate the reference, or reach for `--use-dnaa-for-reference-orientation`.

* **If you are working with a genomic locus**, such as a genomic region of interest, a prophage, an operon and its neighborhood, and all your other sequences are fragmented versions of that same region, then rotation is neither needed nor meaningful. All you want is for the contigs to be reverse-complemented and ordered carefully so they follow the reference as closely as possible. A single contig covering your locus is a perfectly good reference here. Just name it with `--reference` and anvi'o will leave it exactly as it is.

Finally, please note that all of this applies **only to the reference**. Every *other* entry in your %(fasta-txt)s file is free to be as fragmented as they come.

If this is confusing please reach out to us, and we will try to make it less confusing (*I will take "unfulfillable promises" for $200* (well, yes, but you can't know until we/you try)).

### De novo identification of reference orientation

{:.notice}
TL;DR: In this case anvi'o will re-set a start position for your reference genome using the secondary alignments from `minimap2` to find the *most* conserved position across ALL genomes in a data-driven manner.

By default, %(anvi-reorient-genomes)s will rotate the automatically chosen reference genome to an optimal starting position that is conserved across most of the input genomes mentioned in the %(fasta-txt)s: it will first align each genome to the reference **using all possible good alignments (primary + secondary alignments)** to get a true picture of conserved regions across all genomes, and then it will use a 1,000 bp sliding window to search for a region that is shared between all genomes in a greedy fashion (once it finds a region that works for all genomes, it stops).

When you ask the program to use a particular genome as a reference with the `--reference` parameter (because you know it to be a properly circularized complete genome, or simply becasue it is not a genome but a locus you happen to care about), the program will not be tinkering with the start position by searching for a better option that encompasses all genomes (UNLESS you also use `--use-dnaa-for-reference-orientation` flag, which will rotate any reference to the DnaA gene), and do its best to match everything else to it.

If you want anvi'o to auto-select the reference (fewest contigs, longest total length) AND use it **without any rotation**, you can use `--use-auto-reference-as-is`:

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --use-auto-reference-as-is \
                      --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

### What if none of my sequences is a single contig?

Since the reference must be a single contig, a %(fasta-txt)s file in which *every* entry is fragmented (a set of MAGs, say) will make the program stop with an error. You have two options.

**The good option** is to add a single-contig reference to your %(fasta-txt)s file and point anvi'o to it:

{{ codestart }}
anvi-reorient-genomes --fasta-txt %(fasta-txt)s \
                      --reference MY-SINGLE-CONTIG-REFERENCE \
                      --output-dir REORIENTED-FASTA-FILES/
{{ codestop }}

Remember that this reference does not have to be a complete circular genome. If you are working with whole chromosomes, a complete isolate genome from a public database is the natural choice. But if your genomes are really fragmented versions of one genomic region, then a single contig covering that region is just as good -- anvi'o will not rotate a reference you name explicitly, so circularity never enters the picture.

**The limited option** is `--use-auto-reference-as-is`, which lets anvi'o pick the least fragmented entry and use it exactly as it is. Since nothing gets rotated, nothing biologically meaningless happens to the reference, and the program will still put everything on a consistent strand. But please be clear-eyed about what you are giving up: coordinates that come from *different contigs* of a fragmented reference do not form a single continuous axis, so the contig ordering, the `--scaffold-fragmented` output, and the reported 'Start in reference' / 'Start in query' values will not be trustworthy, and the trust labels in the final report are correspondingly much weaker statements. The program will remind you of all this with a big warning when you take this route.

### Why the number of contigs may increase

{:.notice}
TL;DR: A query genome can come out of this program with **more contigs than it went in with**. This is not a bug, and it is not anvi'o being clumsy with your data :) You can always use `--keep-query-contigs-intact` to keep all your contigs intact, but you must read this section first, because if you are using that flag you are trading one problem for another.

Everything this program does (ordering contigs, orienting them, scaffolding them) depends on a single coordinate axis: the reference. In most cases, that axis has a single beginning and a single end. Your query contigs, meanwhile, are linear sequences that know nothing about where the reference happens to begin.

So what happens to a query contig that runs past the beginning or the end of that axis? Say the reference has been rotated so that it starts at the DnaA gene, and one of your query contigs looks like this:

```
        DnaA
         |
         -----------------------------------------------------------
                          reference

xxxxxxxxx---------------------           -----------------
        Contig 1                             Contig 2
```

The `xxxxxxxxx` part of `Contig 1` hangs off the left edge of the reference. Since the reference is circular, that stretch of sequence very likely aligns beautifully to the far *right* end of the reference (the region just *before* DnaA in the actual circular genome). Which means `Contig 1` has *two* places it belongs on the reference axis, at opposite ends of it, and there is no way to write it into a FASTA file in one piece that is honest about both. Written out whole, it would look like this:

```
>Contig_1
xxxxxxxxx---------------------
>Contig_2
-----------------
```

It means every downstream tool would take the first nucleotide of `Contig_1` as the leftmost thing in the genome, even though it belongs at the far right. While identities of the aligned fractions will look great, in its FASTA form, the genome will *not be aligned to the reference*. Rotating the contig here is not a solution either (Meren's note: the author will eat these words as this will be a solution for an entirely different problem we will discuss later in this document). Becasue `Contig_1` is a linear piece of DNA, not a circle, so moving its beginning to its end invents a junction between two stretches of sequence that are not adjacent in the organism.

The only honest option is to cut the contig, which is what anvi'o does by default:

```
DnaA
 |
 -----------------------------------------------------------
                     reference

 ---------------------           ----------------- xxxxxxxxx
 Contig_1_fragment_01             Contig_2         Contig_1_fragment_02
```

%(anvi-reorient-genomes)s will then align each fragment to the reference on their own merit, ordering and orienting them like any other contig, and reporting them individually.

**Anvi'o will never invent sequences and never throws any of them away**. The fragments of a contig will merge into the same contig exactly, so if you concatenate them back together you get your original contig back, nucleotide for nucleotide.

Three things trigger a cut:

* **The contig straddles the position at which the reference begins and ends**, as in the example above. This is extremely common the moment anything is rotated, but it is *not* limited to rotation: if you bring your own reference with `--reference` and it was circularized at a different position than your query genomes were, you will see exactly the same thing.

* **One end of the contig hangs off the beginning or the end of the reference.** That stretch of sequence has nowhere to go on the reference axis, and letting it ride along inside a contig that *is* placed makes it look like it was placed too. It is cut off and reported as an unaligned contig, so you can see it for what it is.

* **The contig closes on itself around the reference.** This is the sneakiest case of all, since such a contig looks perfectly well-behaved: read along it, the reference coordinate marches forward from the very beginning of the reference to its very end, and nothing looks swapped or scrambled. What gives it away is that it skips an enormous stretch of the reference on the way, and that its two ends turn out to be *neighbours* once you stop measuring their distance along the linear reference axis and measure it around the reference instead (see '*Circularly permuted contigs*' below). Written out whole, such a contig goes wherever its largest alignment block belongs, and drags the part of itself that belongs at the very *end* of the genome along with it, past every other contig. So anvi'o cuts it at the stretch it skips, and the two parts go to the two ends of the output file where they belong. Anvi'o will only do this when other contigs of the same genome cover the skipped stretch, since that is the proof that the sequence is present in the query rather than missing from it, and therefore that the junction being cut is an artifact of assembly rather than a deletion worth preserving.

The second one deserves a word of caution, because "has no alignment" is *not* by itself a reason to cut anything. Query genomes are full of variable regions that have no counterpart in the reference, and a contig that sits comfortably inside the reference with an unaligned stretch at one of its ends is just a genome being a genome. So before cutting, anvi'o asks whether there is enough reference left over before the contig's alignment, or after it to accommodate that unaligned sequence. If there is enough reference left over, it would mean that the sequence would not be hanging off anything and anvi'o would leave the poor contig with assumed hypervariability alone and intact. Anvi'o will take out its 8< only when the unaligned stretch is longer than the remaining length of the reference context. Just to put it in different words, a 20,000 nt unaligned tail on a contig that stops 400,000 nt short of the end of the reference is a variable region. But a 20,000 nt unaligned tail on a contig that stops 8 nt short of the end of the reference is dangling and may want to somewhere else.

Anvi'o will not cut a contig into a fragment shorter than `--min-contig-length` (1,000 bp by default). This is what keeps the handful of nucleotides `minimap2` routinely soft-clips off the end of a divergent alignment from being promoted into contigs of their own, and it is also the yardstick for "long compared to the reference remaining": if the part of an unaligned end that genuinely cannot fit on the reference is shorter than this, anvi'o shrugs and keeps the contig whole.

As everything that happens to your contigs/genomes, anvi'o will report every cut on your terminal, naming the contig, where it was cut, why, and where each resulting fragment ended up on the reference for your information:

```
QUERY CONTIGS CUT INTO FRAGMENTS IN BLAH
===============================================
Anvi'o had to cut 1 contig in BLAH into 2 fragments, which means
this genome will have more contigs in the output file than it had
in the input file. (...)

    - 'contig_01' (165,000 nts) was cut into 2 fragments at position 69,844,
      since it straddles the position at which the reference begins and ends:
        > 'contig_01_fragment_01' (69,844 nts) aligns to the reference at
          position 130,157 on the + strand.
        > 'contig_01_fragment_02' (95,156 nts) aligns to the reference at
          position 25 on the + strand.
```

Fragments are numbered along the original contig, from its first nucleotide to its last, so their names tell you how they were related to one another in your input assembly no matter where each of them ends up in the output file. In the example above `contig_01_fragment_02` is written *before* `contig_01_fragment_01` in the output FASTA, because that is where it belongs on the reference.

If you would rather keep your contigs intact (because a downstream tool of yours cannot cope with changing contig names or counts or something) you can always use use `--keep-query-contigs-intact`. In that case no contig will be cut, but you will most likely be doing a mistake unless you really are sure this is the right way to go with your data (in some cases, it will be the right way to go about it). Anvi'o will not pretend otherwise, either: a genome that ends up with a contig that maps to both the beginning and the end of the reference is reported as `NOT TRUSTWORTHY` regardless of how good its alignment statistics look, since its gene order no longer follows the reference (see '*Interpreting trust labels*' under '*Tips, caveats, and runtime*').

### Circularly permuted contigs, or why a contig sometimes looks scrambled

{:.notice}
TL;DR: Sometimes a query contig aligns to the reference in two big pieces whose order is *swapped*, as if half the genome had been rearranged. Usually nothing has been rearranged: the contig was cut out of a cycle in an assembly graph at an arbitrary point. Anvi'o recognizes these and rotates them back into place, and tells you when it does.

Every now and then you will have a contig that will align the reference perfectly, but in a weird way. The reference coordinate will march forward until the end of the contig, and then all of a sudden it will go back, kind of looking like this (notice the middle contig in the reference genome):

![](../../images/anvi-reorient-genomes-circular-permutated-contig-01.png)

Taken at face value, this look like a huge rearrangement (and if that was the case anvi'o would have no business touching it since rearrangements are real differences between genomes). But this one is not since its boundaries are solely determined by the boundaries of the contig. This actually **is an artifact of how the contig came out of the assembler**. Add it to the list of funny things assemblers do.

Here is what is happening here: a contig (as in a fragment of an otherwise circular chromosome) is reported from an assembly graph as a contig due to multiple complications, including unresolvable circularities that form in the construction of the cassembly graph primarily due to biology. In those cases where large circular edges are formed in the assembly graph, the assembly algorithm is literally looking at a circular genome (just sprouting from a braoder circular context like a very big balloon, which, in the example above, represents that middle contig in HIMB1485). So the algorithm cuts the circualr piece and reports a contig. But there are as many possible start positions for that pseudo-circular subgraph as there are nucleotides in that *balloon*. Depending on the assembly algorithm, the reported fragments may end up not starting from where the balloon originally conneced to the broader genomic context. And as a result, against a better reference, the linear form of that fragment as if it was *circularly permutated* which will make the the contig close on itself.

%(anvi-reorient-genomes)s will correct such contigs since they completely violate the synteny of genes in the actual genomic context. But it will require both of these to be correct before touching such a contig:

1. Read along the contig the alignment blocks jump backwards in reference coordinates **exactly once**. Since a rotation of this sort can explain one and only one jump.
2. The reference position where th e**first** block of the candidate contig begins is immediately after where its **last** block ends, within `--min-contig-length` nucleotides. So that the two ends are genuine neighbours, and rotating the contig will not invent any chimeric adjacency that is not already in the sequence (like your assembly algorithm did :') (sad lol)).

That second distance is measured **around** the reference rather than along it, since the sequences these contigs come from are circular and the position at which the reference begins and ends is an arbitrary point on that circle. A contig whose first block begins 6 nts after the beginning of the reference and whose last block ends 7 nts before the end of a 1,453,515 nt reference has its two ends 13 nts apart -- and *not* a megabase apart, which is all a linear axis can see. Getting this wrong is exactly how a contig ends up spanning the whole reference in the output file (which is the third cutting case described under '*Why the number of contigs may increase*' above).

When anvi'o rotates the contig into co-linearity, you will see messages like these in your terminal:

```
CIRCULARLY PERMUTED CONTIGS ROTATED IN HIMB1485
===============================================
Anvi'o rotated 1 contig in 'HIMB1485'. Contigs like these come out of an
assembler having been cut out of a cycle in its assembly graph (...)

    - 'HIMB1485_000000000002' (961,822 nts) was rotated by 527,750 nts, since
      its two ends are 361 nts apart on the reference and it therefore closes
      on itself. It now aligns to the reference in 3 pieces.
```

The contig will keep its name, keeps every one of its nucleotides, and does not become two contigs. Only the point at which it is considered to start has moved.

Here is the demonstration that it works using the example above where one circular rotation of a linear contig (as crazy as it sounds) collapses a contig that (in this particular case) aligned in five scattered blocks into a single one-megabase alignment, changing the reference coverage of the query from 64%% to 94.7%%:

![](../../images/anvi-reorient-genomes-circular-permutated-contig-02.png)

Again, if you would rather this did not happen, `--keep-query-contigs-intact` will stop any rotations of query contigs even if they set themselves on fire.

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

Here is an example output for two circular genomes where the query nicely aligned to the reference (the alignment plots always show before and after the alignment):

![](../../images/anvi-reorient-genomes-circular.png)

The fact that 'Start in query' and 'Start in reference' is identical is great news, and shows that the program was able to set a start position that works perfectly. The high coverage alignment between the two is also great as it shows that you brough together genomes that are quite closely related. The approximate ANI says 100%%, and it is not to be trusted. For instance, the true ANI between these two genomes in this example is about 93%%. But since we are cutting corners using the `minimap2` aligned bits between the two, we really really overestimate the ANI value. It is useful, but please don't rely on the apprximate ANI reported here for anything serious.

---

Here is another example for two circular genomes where the query did not aligned well to the reference, which means, the reorientation is highly unreliable:

![](../../images/anvi-reorient-genomes-miserable.png)

The program will give you something to read and think about at the very end:

```
(...)

FINAL REPORT
===============================================
Your genome reorientation task considered 5 genomes and 1 reference. Some
outcomes are not trustworthy (low alignment coverage can lead to unreliable
orientation). Please review the alignment plots above and the summary below to
decide which FASTA files to use downstream.

Trustworthy ..................................: 1
    - d -> /Users/meren/sandbox/miserable-genomes/reoriented/d.fa

Somewhat OK ..................................: 1
    - g -> /Users/meren/sandbox/miserable-genomes/reoriented/g.fa

Not trustworthy ..............................: 4
    - c -> /Users/meren/sandbox/miserable-genomes/reoriented/c.fa
    - b -> /Users/meren/sandbox/miserable-genomes/reoriented/b.fa
    - f -> /Users/meren/sandbox/miserable-genomes/reoriented/f.fa
    - m -> /Users/meren/sandbox/miserable-genomes/reoriented/m.fa

Failed .......................................: 0

✓ reorient_genomes.py took 0:00:15.245738
```

---

Here is a final example of re-orienting contigs in a fragmented genome using a reference:

![](../../images/anvi-reorient-genomes-scaffolds.png)

Where each contig was nicely turned around, and ordered in the final FASTA file to follow the right order that matches to the genomic context of the reference.

If you wanted, you could use the flag `--scaffold-fragmented` and make anvi'o produce a single contig FASTA file where 'missing' content in the query genome is filled with `N` bases to create your own 'complete genome' to prank your microbiologist colleagues:

![](../../images/anvi-reorient-genomes-scaffolds-scaffolded.png)

{:.warning}
Please don't do it, though -- don't 'scaffold' your MAGs based on a reference and then submit them to databases as single-contig genomes. Procedure above is excellent to have gene synteny preserved despite the missing content, but that is about it.

### What the program does

Here is a more detailed description of what is going on behind the scenes when you press ENTER:

1. **Parse inputs and pick a reference**. Reads %(fasta-txt)s, and if `--reference` is not set, the program picks the genome with the fewest contigs (ties broken by longest total length). All FASTAs are sanity-checked (existence, FASTA format). The reference must be a single contig, and this is **enforced with an error**: no matter whether the reference came from `--reference` or from auto-selection, the program refuses to continue if it has more than one contig. The single exception is `--use-auto-reference-as-is`, which proceeds with a big warning instead of an error. An auto-selected reference additionally needs to be a complete, circular genome, since it is about to be rotated -- a reference named with `--reference` does not, since it never is.

2. **Determine reference orientation**. The program uses one of three strategies to orient the reference genome:
   - **DnaA-based orientation** (if `--use-dnaa-for-reference-orientation` is set): Calls genes with `prodigal`, searches for the DnaA gene using `hmmsearch` with the Bac_DnaA_C HMM profile, and rotates the reference to start at the DnaA gene position. This provides biologically meaningful orientation for bacterial genomes. This strategy takes precedence over the two below, and applies to a user-specified reference just as much as to an auto-selected one.
   - **De novo optimal position** (if reference is auto-selected without DnaA flag): Aligns each genome to the reference using minimap2 with `--secondary=yes -N 100 -p 0.5` to capture **all possible good alignments** (not just the best one). Builds a coverage map using 1,000 bp bins showing which positions are covered by alignments from each genome. Identifies the position with maximum coverage across all genomes (stopping early if it finds 100%% coverage). The reference is then rotated to start at this optimal position, ensuring all genomes will start at a conserved region that is genuinely shared across the dataset.
   - **User-specified reference** (if `--reference` is set without the DnaA flag): Uses the reference genome as-is without rotation.
   - **Auto-selected reference used as-is** (if `--use-auto-reference-as-is` is set): Auto-selects the reference (fewest contigs, longest total length) but skips any rotation, treating the auto-chosen genome exactly like a user-specified one.

**For circular genomes (single-contig):**

3. **Initial alignment (reference vs query)**. Runs `minimap2` (preset `asm5`, with as many threads as you asked for using `--threads`) to align each query to the reference and identifies a primary anchor near reference position 0. If the anchor is on `-` strand, the query is reverse-complemented. The query is rotated so that reference position 0 maps onto the query (first snap).

4. **Second alignment and snap**. Re-aligns the rotated query, finds the primary anchor with the smallest reference start, and rotates again to bring reference 0 onto the query (second snap).

5. **Snap-to-zero with a ref0-focused anchor**. Aligns once again, picks the primary anchor closest to reference position 0, rotates, aligns once more, and applies a final snap so that reference position 0 maps to query position 0.

6. **Iterative correction for perfect alignment**. After the final snap, the program checks if the alignment truly starts at position 0 in both the query and reference. If not (e.g., due to `minimap2` soft-clipping divergent regions), it calculates the necessary rotation, applies it, and re-aligns. This iterates up to 5 times or until the genomes are perfectly aligned at position 0.

**For fragmented genomes (multi-contig MAGs or draft assemblies):**

3. **Contig filtering**. Contigs shorter than `--min-contig-length` (default: 1,000 bp) are excluded from processing.

4. **Individual contig alignment**. Each contig is independently aligned to the reference using `minimap2`.

5. **Rotating circularly permuted contigs**. Contigs that an assembler cut out of a cycle in its assembly graph -- recognizable because their two ends are neighbours on the reference -- are rotated back into co-linearity (see 'Circularly permuted contigs' above). Use `--keep-query-contigs-intact` to turn this off.

6. **Cutting contigs that run past the ends of the reference, or that close on themselves around it**. Contigs that have no single position on the reference to be placed at are cut into fragments, and each fragment is aligned to the reference on its own merit (see 'Why the number of contigs may increase' below). This happens after the rotation step, since a contig that has been made co-linear often no longer needs to be cut at all. Use `--keep-query-contigs-intact` to turn this off.

7. **Contig ordering and orientation**. Contigs are sorted by their alignment position on the reference genome. Contigs aligned to the reverse strand are reverse-complemented to match the reference orientation.

8. **Checking the layout that follows**. Once anvi'o knows where every contig goes, it checks whether the resulting order can actually be written into a FASTA file: a contig that maps to both the beginning and the end of the reference, with the sequences of other contigs belonging in between, makes the gene order in the output file unreliable no matter how good the alignments are. Anvi'o reports such a genome as `NOT TRUSTWORTHY` and names the contig responsible (see '*Interpreting trust labels*' below). With the cutting step above in place this should only happen when you use `--keep-query-contigs-intact`.

9. **Output generation**. By default, contigs are written as separate sequences in the output FASTA, ordered and oriented to match the reference. If `--scaffold-fragmented` is used, contigs are concatenated into a single sequence with N-padding representing gaps based on the reference genome distances.

**For all genomes:**

10. **Write outputs**. Copies the reference %(fasta)s to the output directory (potentially rotated if auto-selected). Writes each reoriented query %(fasta)s under the same name (and using the original extension).

11. **Report per-genome stats and alignment plots**. For each genome, the program reports the orientation outcome (whether it was `TRUSTWORTHY`, `SOMEWHAT OK`, or `NOT TRUSTWORTHY`, based on alignment coverage and, for fragmented genomes, on whether the resulting contig layout follows the reference), many other statistics, and synteny ribbon plots showing the alignment patterns before and after reorientation to visualize orientation quality.

12. **Final report**. Summarizes the number of genomes in different trust categories along their output FASTA paths for the user to decide which outputs are safe for downstream analyses.


### Tips, caveats, and runtime

* **Interpreting "Start in query" and "Start in reference" (circular genomes only)**: These numbers show where `minimap2` primary alignments begin. **When both values are equal** (e.g., both 61), your genomes are correctly aligned at the same biological position. The small offset is just natural sequence variation in the first few dozen base pairs since `minimap2` soft-clips divergent regions at the beginning. In fact, Meren was so concerned about it, he triple-checked this, and discovered that despite differences in offsets across genomes, the amino acid sequences of the first genes in every aligned genome was identical even though their 'start' positions in individual genomes differed. So that's that. But **when the values differ** (e.g., query=0, reference=106), the alignment may have issues and the genome might not be well-oriented.

* **Interpreting trust labels**: For circular genomes, low coverage (either genome covered <50%%) yields `NOT TRUSTWORTHY`. For fragmented genomes, the label considers both reference coverage and the percentage of contigs that aligned. This often means the genome is too divergent or structurally different to reorient confidently. Check the alignment plots and stats before using such outputs. Genomes with `TRUSTWORTHY` labels are safe for downstream comparative analyses.

* **A fragmented genome is also `NOT TRUSTWORTHY` when one of its contigs spans the entire reference**, no matter how good its coverage and ANI are. Coverage, ANI, and the fraction of contigs that aligned all describe how well *sequences* align, and they don't really say nothing about whether the contigs ended up in the right *order* in the output FASTA file. If the user prevented anvi'o to cuts or rotate contigs when needed (i.e., by using `--keep-query-contigs-intact`), and/or when the order is broken despite the best attempts of this program, the genome will be marked NOT TRUSTWORTHY, and anvi'o will print out various warnings in the output to name and shame those contigs:

```
CONTIGS SPANNING THE ENTIRE REFERENCE IN BLAH
===============================================
One contig in 'BLAH' maps to BOTH the beginning and the end of the reference,
and the sequences of other contigs belong in between. (...)

    - 'contig_000000000004' (217,590 nts) aligns to a 1,453,515 nt reference
      both at position 6 and at position 1,372,757, leaving 1,252,500 nts of
      reference in between that it does not cover itself, but 7 contigs of this
      genome ('BLAH_000000000001', 'BLAH_000000000008', 'BLAH_000000000003'
      and 4 more) do.
```

* **Scaffolding fragmented genomes**: By default, fragmented genomes are written with contigs as separate sequences (ordered and oriented). Use `--scaffold-fragmented` to concatenate them into a single sequence with N-padding representing gaps. While this maximizes gene synteny for comparative analyses, **do not submit N-scaffolded MAGs as complete genomes** to public databases. The scaffolding is reference-based and may not represent the true genomic structure.

* **Auto-selected reference and optimal start**: When you don't specify `--reference`, the program picks the genome with the fewest contigs (ties broken by longest total length), preferring complete circular genomes over fragmented ones. Since anvi'o is about to rotate this auto-selected reference, it should really be a complete circular genome (since if your sequences are linear loci rather than whole chromosomes, rotating them is meaningless, so name your reference explicitly with `--reference` instead). It then rotates the reference to a conserved position across your dataset using secondary alignments to capture all conserved regions (not just the single best alignment). This is especially useful for sets of closely related genomes where you want them all to start at a biologically meaningful position. For bacterial genomes, consider using `--use-dnaa-for-reference-orientation` for even more consistent results based on the replication origin. If you want a specific genome or starting position, use `--reference` to override this behavior. If you want anvi'o to auto-select the reference but skip the rotation step (i.e., keep it exactly as it appears in the input FASTA), use `--use-auto-reference-as-is`.

* **DnaA-based orientation benefits**: When working with bacterial genomes, the `--use-dnaa-for-reference-orientation` flag typically produces highly consistent alignments (e.g., all genomes starting within a few base pairs of each other) because it uses biological knowledge (the replication origin) rather than purely sequence-based heuristics. This can be particularly valuable for downstream synteny analyses or when comparing gene order across closely related strains.

* **Approximate ANI**: This is quite a garbage number at this point, and it is no one's fault. When `minimap2` emits `dv`, ANI is computed as `(1 - dv) x 100`. Otherwise it falls back to `nmatch/alen x 100`. So treat it as a quick proxy, not an accurate ANI calculation. It is not the purpose here.

* **Circular ambiguity**: Circular genomes can align equally well at different offsets. The program applies multiple snaps and iterative corrections to align reference position 0 to query position 0, but in highly repetitive cases the true biological origin may still be ambiguous.

* **Visualization options**: By default, alignment plots are generated to help you assess the quality of reorientation. Use `--skip-visualizing-alignments` to disable plotting for faster processing when you only need the FASTA files. Customize plot dimensions with `--plot-width` and `--plot-height` (note: widths below 100 characters may not display properly).

* **Runtime**: Each circular genome triggers several `minimap2` runs and `seqkit` rotations. Fragmented genomes are faster as each contig is aligned only once. The optimal start finding step (when auto-selecting reference) adds an initial survey phase that uses secondary alignments for more accurate conserved region detection. The `--use-dnaa-for-reference-orientation` flag adds gene calling and HMM search overhead (a few seconds for a typical bacterial genome). Overall, it takes no more than 30 seconds on a laptop computer to reorient 30 SAR11 genomes using the de novo approach, and slightly longer with DnaA-based orientation.
