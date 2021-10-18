A [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file that does not necessarily meet the standards of a %(contigs-fasta)s. While it is not necessary for all programs, if a given anvi'o program requires a %(contigs-fasta)s, the program %(anvi-script-reformat-fasta)s can turn a regular fasta into a %(contigs-fasta)s with the flag `--simplify-names`.

### What is a FASTA file?

A FASTA file typically contains one or more DNA, RNA, or amino acid sequences that are formatted as follows:

```
>SEQUENCE_ID VARIOUS_SEQUENCE_DATA
SEQUENCE
(...)
```

The line that starts with the character `>` is also known as the 'defline' for a given sequence. The `VARIOUS_SEQUENCE_DATA` region of the defline can be empty, or contain additional data such as the NCBI taxon ID, GI accession number, a text description of the sequence, or the start and end positions if the sequence is a portion of a larger sample. Because the FASTA file format was designed before there weren't even enough electronic calculators on the planet, there is no actual standard format to organize additional information shared in the defline.

The sequence itself is typically written in standard [IUPAC format](https://en.wikipedia.org/wiki/Nucleic_acid_notation), although you may find FASTA files with sequences that contain lower-case letter, mixed letters, no letters, or pretty much anything really. Over the years we have seen everything, and suggest you to take a careful look at your FASTA files before doing anything with them unless you generated them yourself.

You can learn more about the FASTA format on its [glorious Wikipedia page](https://en.wikipedia.org/wiki/FASTA_format).
