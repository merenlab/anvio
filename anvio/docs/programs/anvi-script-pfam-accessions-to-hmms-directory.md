You give this program one or more PFAM accession ids, and it generates an anvi'o compatible HMM directory %(hmm-source)s to be used with %(anvi-run-hmms)s by downloading them from the PFAM database.

### Basic usage

You may either specify a list of PFAM accessions with `--pfam-accessions-list`:

{{ codestart }}
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list PF00001 PF00002 \
                                              -o PROFILE-NAME
{{ codestop }}

Or a file containing this list using `--pfam-accessions-file`. The file should have one accession number per line:

{{ codestart }}
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-file file.txt \
                                              -o PROFILE-NAME
{{ codestop }}

{:.warning}
Please note that the `PROFILE-NAME` will become the name for the HMM source when you use this HMM directory with %(anvi-run-hmms)s. So choose the output directory name accordingly, and make sure (1) it does not conflict with any existing HMM source name in your anvi'o setup, and (2) it is descriptive of the profile you are building.
