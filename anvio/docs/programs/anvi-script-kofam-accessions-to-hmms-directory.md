You give this program one or more KOfam accession ids, and it generates an anvi'o compatible HMM directory %(hmm-source)s to be used with `anvi-run-hmms` by extracting them from your local %(kegg-data)s setup.

### Basic usage

You may either specify a list of KOfam accessions with `--kofam-accessions-list`:

{{ codestart }}
anvi-script-kofam-accessions-to-hmms-directory --kofam-accessions-list K00001 K00121 \
                                               -o PROFILE-NAME
{{ codestop }}

Or a file containing this list using `--kofam-accessions-file`. The file should have one accession number per line:

{{ codestart }}
anvi-script-kofam-accessions-to-hmms-directory --kofam-accessions-file file.txt \
                                               -o PROFILE-NAME
{{ codestop }}

{:.warning}
Please note that the `output` directory name will become the name for the HMM source when you use this HMM directory with `anvi-run-hmms`. So choose the output directory name accordingly, and make sure (1) it does not conflict with any existing HMM source name in your anvi'o setup, and (2) it is descriptive of the profile you are building.

If your KEGG data is not in the default location, you can specify it using `--kegg-data-dir`:

{{ codestart }}
anvi-script-kofam-accessions-to-hmms-directory --kofam-accessions-list K00001 \
                                               --kegg-data-dir /path/to/KEGGDATA \
                                               -o PROFILE-NAME
{{ codestop }}

The `PROFILE-NAME` will be a new directory created by anvi'o.
