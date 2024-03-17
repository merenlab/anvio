You give this program one or more PFAM accession ids, and it generates an anvi'o compatible HMM directory [hmm-source](https://anvio.org/help/main/artifacts/hmm-source/) to be used with `anvi-run-hmms`.

### Basic usage

You may either specify a list of PFAM accession with `--pfam-accessions-list`:

{{ codestart }}
anvi-script-pfam-accessions-to-hmm-directory --pfam-accessions-list ACC1 ACC2 -O output
{{ codestop }}

Or a file containing this list using `--pfam-accessions-file`. The file should have one accession number per line:

{{ codestart }}
anvi-script-pfam-accessions-to-hmm-directory --pfam-accessions-file file.txt -O output
{{ codestop }}

Output folder is specified by `-O` and the folder will be created by anvi'o if it doesn't exist, otherwise anvi'o will exist.