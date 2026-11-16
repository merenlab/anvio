This program writes a %(genes-fasta)s with one entry per gene that has at least one hit in a single %(functions)s source. You can limit the output to genes with hits to specific accession ids in that source.

It works with a %(contigs-db)s, bins in a %(collection)s, or many genomes described in %(external-genomes)s and/or %(internal-genomes)s files.

%(anvi-get-sequences-for-function-hits)s is very similar to %(anvi-get-sequences-for-hmm-hits)s, but unlike the latter, this program does not have any functionality to concatenate genes for downstream phylogenetic or phylogenomic analyses. While this was not difficult to implement, it was a mindful omission on our part. This is primarily due to the fact that genes that share a functional annotation can have extremely low sequence identity depending on the breadth of their representation in databases, so their sequences may not a reliable basis for phylogenetic or phylogenomic inferences even though it may not be immediately obvious for all functions. HMMs can also return not-so-ideal hits, but since they are built from aligned sequences by design, they are generally more appropriate to work with for such goals.

The primary purpose of this program is to give you all matching sequences annotated with a particular function, or all functions listed in a particular source.

### Learn available annotation sources

{{ codestart }}
anvi-get-sequences-for-function-hits -c %(contigs-db)s \
                                     --list-annotation-sources
{{ codestop }}

### Get sequences for all genes with hits in a source

{{ codestart }}
anvi-get-sequences-for-function-hits -c %(contigs-db)s \
                                     --annotation-source Pfam \
                                     -o %(genes-fasta)s
{{ codestop }}

### Get sequences only for genes with hits to specific accession ids

{{ codestart }}
anvi-get-sequences-for-function-hits -c %(contigs-db)s \
                                     --annotation-source Pfam \
                                     --accession-ids PF00117,PF02786 \
                                     -o %(genes-fasta)s
{{ codestop }}

Instead of a comma-separated list, `--accession-ids` also accepts a file with one accession id per line. A gene is reported if any of its hits in the source matches any of the accession ids. If an accession id matches nothing in any of your inputs, anvi'o will stop with an error rather than ignore it.

### Get sequences from bins or genomes

For every bin in a collection (use `--bin-id` or `--bin-ids-file` to focus on specific bins):

{{ codestart }}
anvi-get-sequences-for-function-hits -c %(contigs-db)s \
                                     -p %(profile-db)s \
                                     -C %(collection)s \
                                     --annotation-source COG20_FUNCTION \
                                     -o %(genes-fasta)s
{{ codestop }}

For many genomes:

{{ codestart }}
anvi-get-sequences-for-function-hits -e %(external-genomes)s \
                                     -i %(internal-genomes)s \
                                     --annotation-source KOfam \
                                     --accession-ids K00001 \
                                     -o %(genes-fasta)s
{{ codestop }}

The annotation source must be present in every input %(contigs-db)s.

### Amino acid sequences

Sequences are reported in DNA alphabet by default. Add `--get-aa-sequences` to get amino acid sequences instead.

### Deflines

The default defline looks like this:

```
{genome_name}_{gene_callers_id} source:{source}|accession:{accession}|function:{function}
```

Here `{genome_name}` is the genome name, the bin name, or the project name of the %(contigs-db)s depending on the input. If a gene has multiple hits in the source, its accessions and functions are joined with `!!!`.

You can define your own defline with `--defline-format`, and see all available variables with `--list-defline-variables`. Your format must produce a unique defline for each sequence (e.g., by including `{genome_name}` and `{gene_callers_id}`), otherwise anvi'o will complain.

{{ codestart }}
anvi-get-sequences-for-function-hits -c %(contigs-db)s \
                                     --annotation-source Pfam \
                                     --defline-format "{gene_callers_id}_{accession}" \
                                     -o %(genes-fasta)s
{{ codestop }}
