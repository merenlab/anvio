This program associates the tRNA reads found in your %(contigs-db)s with taxonomy information.

Once these associations are stored in your %(contigs-db)s (represented by a %(trna-taxonomy)s artifact), you'll be able to run %(anvi-estimate-trna-taxonomy)s to use the associations to estimate taxonomy on a larger scale (i.e. for a genome or metagenome).

To run this program, you'll need to have set up two things:
1. a %(trna-taxonomy-db)s, which you can set up by running %(anvi-setup-trna-taxonomy)s.
2. the 'transfer-RNAs' HMM hits in your %(contigs-db)s, which you can set up by running %(anvi-scan-trnas)s

This program will then go through the tRNA hits in your contigs database and search them against the sequences in the [GTDB](https://gtdb.ecogenomic.org/) databases that you downloaded to assign them taxonomy.

### Basic run

The following is a basic run of this program:

{{ codestart }}
anvi-run-trna-taxonomy -c %(contigs-db)s
{{ codestop }}

If you have set up the two requirements listed above, this should run smoothly.

### tRNA-seq seeds

If your %(contigs-db)s is a %(trnaseq-contigs-db)s, the sequences it holds are tRNA-seq seeds, which are labeled with the anticodons read from tRNA transcripts rather than encoded in genes. Reverse transcription reads inosine (I), a modified nucleotide, as G, and bacterial tRNA-Arg-ACG and tRNA-Leu-AAG can carry inosine at wobble position 34, so their seeds can be labeled GCG and GAG. This program therefore searches seeds labeled GCG against both GCG and ACG genes and seeds labeled GAG against both GAG and AAG genes, keeping the hits of the genes whose best alignment has the higher bitscore. The taxonomy of a seed is always stored under its labeled anticodon. In the table this program prints, a column called `as_I34` shows how many seeds were annotated with ACG or AAG genes. The percent identity of a seed with I34 counts its G at position 34 as a mismatch to the A of the genes, so a seed that differs from a gene only at that position has an identity of about 98.6%% for tRNA-Arg and 98.8%% for the longer tRNA-Leu.

### Additional Parameters

When changing these parameters, it might be a good idea to run %(anvi-estimate-trna-taxonomy)s with the `--debug` flag so that you can see what your results look like under the hood.

1. `--max-num-target-sequences`: the number of hits that this program considers for each tRNA sequence before making a final decision for the taxonomy association. The default is 100, but if you want to ensure that you have accurate data at the expense of some runtime, you can increase it.
2. `--min-percent-identity`: the minimum percent alignment needed to consider something a hit.  The default is 90, but if you're not getting any hits on a specific sequence, you can decrease it at the risk of getting some nonsense results.

Finally, this program does not usually have an output file, but if desired you can add the parameter `--all-hits-output-file` to store the list of hits that anvi'o looked at to determine the consensus hit for each sequence. The `searched_anticodon` column in this file gives the anticodon of the genes each hit comes from, which differs from the column `anticodon` only for the tRNA-seq seeds with I34 that are enumerated above.
