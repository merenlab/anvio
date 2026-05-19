This program **annotates genes with functions from the GlobAA gene family database**. It requires a %(globdb-data)s artifact produced by %(anvi-setup-globdb-functions)s.

GlobAA is a curated database of microbial gene families, each accompanied by gene-family-level Local Alignment Score Ratio (LASR) cutoffs. After searching your sequences against the database with DIAMOND, anvi'o uses these cutoffs to decide whether each hit is a genuine annotation or noise. Hits that pass the threshold are stored as `GlobAA` annotations in your contigs database.

### Annotate a contigs database

{{ codestart }}
anvi-run-globdb-functions -c %(contigs-db)s
{{ codestop }}

{:.warning}
Use `--num-threads` for faster DIAMOND searches on multi-core systems.

### Annotate a FASTA file of amino acid sequences

{{ codestart }}
anvi-run-globdb-functions --fasta-file my_proteins.faa \
                          --output-file my_annotations.txt
{{ codestop }}

### Custom data directory

If your %(globdb-data)s lives in a non-default location, point anvi'o to it:

{{ codestart }}
anvi-run-globdb-functions -c %(contigs-db)s \
                          --globdb-data-dir /path/to/globaa/data
{{ codestop }}

Or set the environment variable `ANVIO_GLOBAA_DATA_DIR` once and omit the flag.

### How the cutoffs work

For each DIAMOND hit, anvi'o computes a **Local Alignment Score Ratio (LASR)**, the ratio of the raw DIAMOND alignment score to the theoretical maximum self-alignment score of the query sequence computed from BLOSUM45 diagonal values as implemented by the GlobDB folk (which includes Daan Speth et al.). This is then compared against the LASR threshold (`lasr`), `selfmin`, and `selfmax` values stored in the gene family's YAML entry:

- A query whose self-alignment score falls in the expected range (`selfmin`–`selfmax`) for the family is classified as **correct_length** if the BSR passes the threshold.
- Queries shorter or longer than expected are classified as **too_short** or **too_long** respectively; these still pass and are annotated, but the classification is noted in the function description.
- Hits that do not reach the BSR threshold are labeled **below_cutoff** and are silently discarded.

This is how cutoffs work (broadly speaking), and as soon as there is a resource to cite here, we will update this information.
