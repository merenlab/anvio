This program generates synthetic sequencing reads from a reference %(fasta)s file. It supports Illumina paired-end and single-end short reads, PacBio HiFi, and Oxford Nanopore long reads, with optional SNV injection at controlled multi-allele frequencies.

The key feature is **multi-allele SNV injection**: you can specify exactly which positions should be variable and what the base frequencies should be at each position. This is critical for benchmarking tools that rely on SNV patterns, such as %(anvi-report-dgrs)s, where real biological variability involves 3-4 different bases at a position rather than just 2.

## Quick start with presets

The easiest way to use this program is with a preset. Presets set sensible defaults for read length, insert size, error rate, and quality scores:

{{ codestart }}
anvi-script-gen-reads -f %(fasta)s \
                      -o OUTPUT_PREFIX \
                      --preset illumina-paired
{{ codestop }}

This generates `OUTPUT_PREFIX-R1.fastq` and `OUTPUT_PREFIX-R2.fastq` with 150 bp paired-end reads at 50X coverage.

Available presets:

**Short reads (Illumina):**

| Preset | Read type | Length | Insert size | Error rate | Quality |
|--------|-----------|--------|-------------|------------|---------|
| `illumina-paired` | paired-end | 150 bp | 450 bp (std 50) | 0.5% | ? (Q30) |
| `illumina-single` | single-end | 150 bp | - | 0.5% | ? (Q30) |

**Long reads (PacBio):**

| Preset | Read type | Length | Distribution | Error rate | Quality | Notes |
|--------|-----------|--------|--------------|------------|---------|-------|
| `pacbio-hifi` | long-distributed | 15 kb (std 3.5 kb) | normal | 0.1% | F (Q37) | Modern HiFi circular consensus sequencing |
| `pacbio-clr` | long-distributed | 15 kb (std 8 kb) | normal | 12% | . (Q13) | Legacy CLR (pre-HiFi), noisy but long |

**Long reads (Oxford Nanopore):**

| Preset | Read type | Length | Distribution | Error rate | Quality | Notes |
|--------|-----------|--------|--------------|------------|---------|-------|
| `ont-r9` | long-distributed | 5 kb (std 4 kb) | lognormal | 6% | 3 (Q18) | Legacy R9.4.1 chemistry |
| `ont-r10` | long-distributed | 8 kb (std 5 kb) | lognormal | 1% | = (Q28) | Modern R10.4.1 with super-accuracy basecalling |
| `ont-ultralong` | long-distributed | 50 kb (std 40 kb) | lognormal | 2% | : (Q25) | PromethION ultralong runs, very high variance |

ONT presets use a **lognormal** length distribution, which produces the right-skewed shape typical of nanopore data (mode lower than mean, with a long tail of very long reads). PacBio presets use a normal distribution, which better reflects the tighter length control of SMRT sequencing.

You can override any preset parameter individually. For example, to use the Illumina paired-end preset but with 100X coverage and 250 bp reads:

{{ codestart }}
anvi-script-gen-reads -f %(fasta)s \
                      -o OUTPUT_PREFIX \
                      --preset illumina-paired \
                      --coverage 100 \
                      --read-length 250
{{ codestop }}

## SNV injection with a mutations file

To inject SNVs at specific positions with controlled allele frequencies, provide a TAB-delimited mutations file:

{{ codestart }}
anvi-script-gen-reads -f %(fasta)s \
                      -o OUTPUT_PREFIX \
                      --preset illumina-paired \
                      --mutations-file mutations.tsv
{{ codestop }}

The mutations file must have the following columns:

```
contig_name     position    freq_A    freq_T    freq_C    freq_G
contig_1        1000        0.25      0.25      0.25      0.25
contig_1        2000        0.0       0.4       0.3       0.3
contig_2        500         0.5       0.5       0.0       0.0
```

Positions are **0-indexed** and frequencies must sum to 1.0 for each row. Each frequency represents the probability that a read covering that position will carry that base. For example, a position with `freq_A=0.25, freq_T=0.25, freq_C=0.25, freq_G=0.25` will have all four bases represented equally across reads -- the kind of multi-allele variability you see in real DGR variable regions.

## Random SNV injection

If you don't need precise control over SNV positions, you can have the program randomly place SNVs at a given density:

{{ codestart }}
anvi-script-gen-reads -f %(fasta)s \
                      -o OUTPUT_PREFIX \
                      --preset illumina-paired \
                      --snv-density 0.01 \
                      --num-alleles 3
{{ codestop }}

This places SNVs at approximately 1% of positions (1 per 100 bp), each with 3 different alleles at random frequencies. `--num-alleles` can be 2, 3, or 4.

## Custom read type (no preset)

You can specify all parameters manually instead of using a preset:

{{ codestart }}
anvi-script-gen-reads -f %(fasta)s \
                      -o OUTPUT_PREFIX \
                      --read-type paired-end \
                      --read-length 150 \
                      --insert-size 300 \
                      --insert-size-std 50 \
                      --coverage 100 \
                      --error-rate 0.005
{{ codestop }}

Available read types are `paired-end`, `single-end`, `long-fixed`, and `long-distributed`.

## Reproducibility

All runs are deterministic by default (seed = 42). To get a different random realization, change the seed:

{{ codestart }}
anvi-script-gen-reads -f %(fasta)s \
                      -o OUTPUT_PREFIX \
                      --preset illumina-paired \
                      --seed 123
{{ codestop }}

## Typical workflow

Generate reads, then run the standard anvi'o pipeline:

{{ codestart }}
# generate reads with known SNVs
anvi-script-gen-reads -f reference.fa \
                      -o sample_01 \
                      --preset illumina-paired \
                      --coverage 100 \
                      --mutations-file mutations.tsv

# build contigs database
anvi-gen-contigs-database -f reference.fa \
                          -o contigs.db

# map reads
bowtie2-build reference.fa reference
bowtie2 -x reference \
        -1 sample_01-R1.fastq \
        -2 sample_01-R2.fastq \
        -S sample_01.sam

# convert and sort
samtools view -bS sample_01.sam | samtools sort -o sample_01-raw.bam
anvi-init-bam sample_01-raw.bam -o sample_01.bam

# profile
anvi-profile -i sample_01.bam \
             -c contigs.db \
             -o sample_01_profile
{{ codestop }}
