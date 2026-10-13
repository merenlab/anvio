This program processes a %(genbank-file)s, and converts it into anvi'o friendly artifacts: namely, a %(contigs-fasta)s, %(external-gene-calls)s and a %(functions-txt)s.

The %(contigs-fasta)s and %(external-gene-calls)s can be given to %(anvi-gen-contigs-database)s to create a %(contigs-db)s, and then you can use %(anvi-import-functions)s to bring the function data (in the %(functions-txt)s) into the database. Then you'll have all of the data in your %(genbank-file)s converted into a single %(contigs-db)s, which you can use for a variety of anvi'o analyses.

### Features processed by default

By default, %(anvi-script-process-genbank)s will `CDS`, `tRNA`, and `rRNA` features by default.

- `CDS` features are mapped to the anvi'o `CODING` gene call type.
- `tRNA` and `rRNA` features are mapped to the `NONCODING` gene call type.

### Handling problematic features

Genomic data often contains features that anvi'o may find difficult to process using standard workflows, such as gene calls with internal stop codons or frameshifts. This script identifies such features and handles them gracefully by reclassifying them as `NONCODING`:

1. **Pseudogenes**: Any `CDS` explicitly marked as a `/pseudogene` or having `/pseudo` in its GenBank qualifiers will be reclassified as `NONCODING`.
2. **Internal Stops and Frameshifts**: Any `CDS` with notes indicating internal stops or frameshifts (based on common NCBI PGAP terms) will also be reclassified as `NONCODING`.

This approach ensures that these features are preserved in your %(contigs-db)s without triggering translation errors during database creation.

### Notes on Output

The parameters of this program entirely deal with the outputs. Besides telling the program where to put them, you can also give the function annotation source (in the %(functions-txt)s) a custom name.

One important note about this conversion is the following: During the conversion of GenBank entries, anvi'o will assign a new gene call id to each entry, breaking the link between locus tags defined in the GenBank file and the gene entries that will later appear in the anvi'o %(contigs-db)s. One way to avoid this is to use the flag `--include-locus-tags-as-functions`, which will instruct anvi'o to add a new 'function' source for each gene in the output file for functional annotations so that the user can trace back a given gene call to the original locus tag.
