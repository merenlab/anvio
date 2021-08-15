A **directory of data** downloaded from the [KEGG database resource](https://www.kegg.jp/) for use in function annotation and metabolism estimation.

It is created by running the program %(anvi-setup-kegg-kofams)s. Not everything from KEGG is included in this directory, only the information relevant to downstream programs. The most critical components of this directory are KOfam HMM profiles and the %(modules-db)s which contains information on metabolic pathways as described in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html).

Programs that rely on this data directory include %(anvi-run-kegg-kofams)s and %(anvi-estimate-metabolism)s.

## Directory Location
The default location of this data is in the anvi'o folder, at `anvio/anvio/data/misc/KEGG/`. 

You can change this location when you run %(anvi-setup-kegg-kofams)s by providing a different path to the `--kegg-data-dir` parameter:

{{ codestart }}
anvi-setup-kegg-kofams --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

If you do this, you will need to provide this path to downstream programs that require this data as well.

## Directory Contents

Here is a schematic of how the %(kegg-data)s folder will look after setup:

```
KEGG
 |- MODULES.db
 |- ko_list.txt
 |- modules.keg
 |- HMMs
 |   |- Kofam.hmm
 |   |- Kofam.hmm.h3f
 |   |- (....)
 |
 |- modules
 |   |- M00001
 |   |- M00002
 |   |- (....)
 |
 |- orphan_data
     |- 01_ko_fams_with_no_threshold.txt
     |- 02_hmm_profiles_with_ko_fams_with_no_threshold.hmm

```

Typically, users will not have to work directly with any of these files, as downstream programs will interface directly with the %(modules-db)s. 

However, for the curious:
`ko_list.txt`, `modules.keg`, and all files in the `modules` subfolder are flat text files downloaded from the [KEGG website](https://www.genome.jp/kegg/). The data in these files are processed and organized into the %(modules-db)s for easier programmatic access. 

The `HMMs` subfolder contains a file of concatentated KOfam profiles (also originally downloaded from [KEGG](https://www.genome.jp/ftp/db/kofam/)), as well as the indexes for this file. Some KOfam profiles do not have a score threshold in the `ko_list.txt` file - these profiles and their corresponding entries from that file live in the `orphan_data` directory. Please note that KOs from the `orphan_data` directory will *not* be annotated in your %(contigs-db)s when you run %(anvi-run-kegg-kofams)s.
