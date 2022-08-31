The %(hmm-list)s file is a TAB-delimited file with at least three columns:

* `name`: The name of the HMM. If you are using an external HMM it MUST match the name found in the `genes.txt`
* `source`: Name of the collection of HMMs the HMM is found in e.g. Bacterial_71. If you are using an external HMM, simply put the name of the directory.
* `path`: If using an [Default HMM sources](http://127.0.0.1:4000/help/main/artifacts/hmm-source/#default-hmm-sources) simply put "INTERNAL". On the other hand, if you are using a [user-defined HMM sources](http://127.0.0.1:4000/help/main/artifacts/hmm-source/#user-defined-hmm-sources) please put the full path to the anvi'o formatted HMM directory.

Here is an example of an %(hmm-list)s using the HMM Ribosomal_L16 and Ribosomal_S2 from the internal anvi'o collection Bacteria_71:

| name          | source      | path     |
|---------------|-------------|----------|
| Ribosomal_L16 | Bacteria_71 | INTERNAL |
| Ribosomal_S2  | Bacteria_71 | INTERNAL |

You can also use external %(hmm-source)ss! An easy way to get an anvi'o ready HMM directory is to use the script %(anvi-script-pfam-accessions-to-hmms-directory)s to download a [Pfam HMM](https://pfam.xfam.org/).

{{ codestart }}
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list PF00016 -O RuBisCO_large_HMM
{{ codestop }}

Here is what the `hmm_list.txt` should look like with a combination of internal and external %(hmm-source)ss:  

| name          | source            | path                       |
|---------------|-------------------|----------------------------|
| Ribosomal_L16 | Bacteria_71       | INTERNAL                   |
| Ribosomal_S2  | Bacteria_71       | INTERNAL                   |
| RuBisCO_large | RuBisCO_large_HMM | PATH/TO/RuBisCO_large_HMM/ |