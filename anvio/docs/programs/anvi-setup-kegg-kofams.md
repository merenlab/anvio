%(anvi-setup-kegg-kofams)s downloads and organizes data from KEGG for use by other programs, namely %(anvi-run-kegg-kofams)s and %(anvi-estimate-kegg-metabolism)s. It downloads HMM profiles from the KOfams database as well as metabolism information such as that stored in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). The program generates a directory with this data (%(kegg-db)s), which by default is located at `anvio/anvio/data/misc/KEGG/`.

### Set up KEGG data

{{ codestart }}
anvi-setup-kegg-kofams
{{ codestop }}

### Set up KEGG data in non-default location

{{ codestart }}
anvi-setup-kegg-kofams --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

An important thing to note about this program is that it has rigid expectations for the format of the KEGG data that it works with. Future updates to KEGG may break things such that the data can no longer be directly obtained from KEGG or properly processed. In the event that this happens, this program still has you covered. You can provide an archived KEGG data directory to the script, and it will unpack that archive and make sure things are all in order.

### Set up from archived KEGG data

{{ codestart }}
anvi-setup-kegg-kofams --kegg-archive KEGG_archive.tar.gz
{{ codestop }}
