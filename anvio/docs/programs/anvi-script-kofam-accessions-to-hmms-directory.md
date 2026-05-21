You give this program one or more KOfam accession ids, and it generates an anvi'o compatible HMM directory %(hmm-source)s to be used with `anvi-run-hmms` by extracting them from your local %(kegg-data)s setup.

### Basic usage

You may either specify a list of KOfam accessions with `--kofam-accessions-list`:

{{ codestart }}
anvi-script-kofam-accessions-to-hmms-directory --kofam-accessions-list K00001 K00121 -O output
{{ codestop }}

Or a file containing this list using `--kofam-accessions-file`. The file should have one accession number per line:

{{ codestart }}
anvi-script-kofam-accessions-to-hmms-directory --kofam-accessions-file file.txt -O output
{{ codestop }}

If your KEGG data is not in the default location, you can specify it using `--kegg-data-dir`:

{{ codestart }}
anvi-script-kofam-accessions-to-hmms-directory --kofam-accessions-list K00001 --kegg-data-dir /path/to/KEGG -O output
{{ codestop }}

The output directory is specified by `-O`. It will be created by anvi'o if it doesn't exist.
