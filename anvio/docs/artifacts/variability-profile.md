As an artifact, this describes the variability information about a single sample calculated when you ran %(anvi-profile)s. To examine variability across samples, you'll want to use this information (which is stored within your %(profile-db)s) to run %(anvi-gen-variability-profile)s.

## Details about Variability

In the context of anvi'o, variability means divergence of environmental populations from the reference used to perform metagenomic read recruitment.

Here, the term "population" describes an assemblage of co-existing microbial genomes in an environment that are similar enough to map to the context of the same reference genome.

The variability profile of a metagenome enables studies of [microbial population genetics with anvi'o](http://merenlab.org/2015/07/20/analyzing-variability/).

There are three types of variability the program %(anvi-profile)s can characterize and store: substitutions, indels, and read-edge clipping events.

### Substitutions: SNVs, SCVs, SAAVs

Anvi'o can make sense of single-nucleotide variants (SNVs), single-codon variants (SCVs), and single-amino acid variants (SAAVs). See [this article](http://merenlab.org/2015/07/20/analyzing-variability) for more information.

You can learn the name of the table in which anvi'o stores this in a given %(profile-db)s by running this command in your anvi'o environment:

``` bash
python -c 'import anvio.tables as t; print(t.variable_nts_table_name)'
```

This will tell you about its structure:

``` bash
python -c 'import anvio.tables as t; print(t.variable_nts_table_structure)'
```

### Indels: insertions and deletions

Anvi'o can also characterize insertions and deletions found within an environment based on short-read recruitment results and will store in the following table:

``` bash
python -c 'import anvio.tables as t; print(t.indels_table_name)'
```

**Notes for programmers**: The convention for the start position of an insertion is defined like so:

```
    pos_in_contig ...0123456 7890123456
    reference     ...CTACTAC TACTTCATGA...
    read              TACTAC TAC
    insertion               └──ACTG
```

In this case, the start position of the insertion in the contig is 6. The insertion _follows_ the position it is defined by. This is opposite to IGV, in which the insertion _precedes_ the position it is defined by.

For deletions, there is no such ambiguity in the start position, since the deletion starts on a reference position, not in between two reference positions.

### Clip events: soft and hard CIGAR clips at read edges

Soft (CIGAR `S`) and hard (CIGAR `H`) clips appear at the start or end of a read's alignment when the aligner could not extend the alignment further. They are not random sequencing artifacts; pile-ups of clips at the same reference position are strong evidence of structural-variant breakpoints, mobile-element insertion sites, junctions between divergent strains in metagenomes, and other recombination-like events. Anvi'o stores each clip event in the following table:

``` bash
python -c 'import anvio.tables as t; print(t.clippings_table_name)'
```

The structure:

``` bash
python -c 'import anvio.tables as t; print(t.clippings_table_structure)'
```

Each row represents one clip event with a single signature `(pos, side, type, sequence, length, partner_*)`. Multiple reads with the same signature collapse into one row whose `count` reflects the number of supporting reads. The columns are:

* `pos` / `pos_in_contig` — the breakpoint position (i.e. the read's first or last aligned reference base, depending on `side`).
* `side` — `L` if the clip is on the 5' end of this record's alignment, `R` if on the 3' end.
* `type` — `SOFT` (CIGAR `S`; bases are kept in this record's SEQ) or `HARD` (CIGAR `H`; bases live in a sibling record's SEQ, typically the primary alignment).
* `state` — one of three values that describes what the bases on the *outside* of the clip do in the read:
    * `JUNCTION` — the read continues immediately into a sibling alignment listed in the SAM `SA` tag. `sequence` is empty (the bases are in the sibling alignment, already in the contigs database). `partner_*` is populated.
    * `JUNCTION_WITH_GAP` — there is a sibling alignment in the outside direction, but with some unmapped bases between this clip and the sibling's nearest edge. `sequence` carries those gap bases. `partner_*` is populated.
    * `UNMAPPED` — there is no sibling alignment in the outside direction. The clip's outside bases simply do not map anywhere in the reference. `sequence` carries those bases. `partner_*` is empty.
* `length` — the original CIGAR clip length (independent of the gap size; `len(sequence)` can be smaller).
* `count` / `coverage` — number of reads supporting this exact event, and the per-position coverage at the breakpoint.
* `partner_contig` / `partner_junction_pos` / `partner_strand` — describe the sibling alignment that this clip joins (if any). `partner_junction_pos` is the partner's reference coordinate *adjacent to the junction* (not the SAM `SA`-tag leftmost-position): the partner's L edge for `R` clips with same-strand siblings, the partner's R edge otherwise, with the L/R flip applied for opposite-strand siblings.

**A note on what `sequence` carries.** With the three-state model above, a non-empty `sequence` always represents *novel* nucleotide content — bases that do not appear anywhere in the contigs database. For `UNMAPPED` clips this is the unmapped tail; for `JUNCTION_WITH_GAP` clips it is the gap between us and the partner. For clean `JUNCTION` clips the column is empty by design (the bases live in the partner alignment, already recorded in the contigs db). This makes the column a direct entry point for exploring microdiversity at clip events.

**A note on what's NOT stored.** Hard-clipped bases (`type='HARD'`) do not appear in the supplementary record's BAM `SEQ`. When the state requires the actual bases (e.g., `UNMAPPED` on a HARD clip), anvi'o fetches the primary alignment record by query name during profiling, extracts the bases from there, and reverse-complements them if the strands differ. Secondary alignments (BAM flag `0x100`) are not used at all for clip profiling — only primary and supplementary records contribute.

**A note on the BAM aligner.** Whether your BAM contains clip information at all depends on which aligner produced it. minimap2 (any mode), `bwa mem`, `bwa bwasw` (BWA-SW), `bowtie2 --local`, HISAT2, and STAR all emit clips by default. `bowtie2` in its default end-to-end mode and `bwa aln`/`samse`/`sampe` do not. %(anvi-profile)s inspects the `@PG` records of the BAM header and auto-skips clip profiling — with a warning — when the aligner is recognized as non-clip-emitting, so that an empty `clippings` table is never silently misleading.
