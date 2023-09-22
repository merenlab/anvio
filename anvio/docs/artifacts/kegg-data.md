A **directory of data** downloaded from the [KEGG database resource](https://www.kegg.jp/) for use in function annotation and metabolism estimation.

It is created by running the program %(anvi-setup-kegg-data)s. Not everything from KEGG is included in this directory, only the information relevant to downstream programs. The most critical components of this directory are KOfam HMM profiles and the %(modules-db)s which contains information on metabolic pathways as described in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html), as well as functional classification hierarchies from [KEGG BRITE](https://www.genome.jp/kegg/brite.html).

Programs that rely on this data directory include %(anvi-run-kegg-kofams)s and %(anvi-estimate-metabolism)s.

## Directory Location
The default location of this data is in the anvi'o folder, at `anvio/anvio/data/misc/KEGG/`.

You can change this location when you run %(anvi-setup-kegg-data)s by providing a different path to the `--kegg-data-dir` parameter:

{{ codestart }}
anvi-setup-kegg-data --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

If you do this, you will need to provide this path to downstream programs that require this data as well.

## Directory Contents

Here is a schematic of how the %(kegg-data)s folder will look after setup:

```
KEGG
 |- MODULES.db
 |- ko_list.txt
 |- modules.keg
 |- hierarchies.json
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
 |- BRITE
 |   |- ko00001
 |   |- ko00194
 |   |- (....)
 |
 |- orphan_data
     |- 01_ko_fams_with_no_threshold.txt
     |- 02_hmm_profiles_with_ko_fams_with_no_threshold.hmm

```

### What is this data?

Typically, users will not have to work directly with any of these files, as downstream programs will interface directly with the %(modules-db)s.

However, for the curious, here is a description of each component in this data directory:
- `ko_list.txt`: a tab-delimited file from the [KEGG KOfam](https://www.genome.jp/ftp/db/kofam/) resource that describes the KOfam profile for each KEGG Ortholog (KO). It contains information like the bitscore threshold (used to differentiate between 'good' and 'bad' hits when annotating sequences), the function definition, and various data about the sequences used to generate the profile.
- The `HMMs` subfolder: contains a file of concatentated KOfam profiles (also originally downloaded from [KEGG](https://www.genome.jp/ftp/db/kofam/)), as well as the indexes for this file.
- The `orphan_data` subfolder: contains KOfam profiles for KOs that do not have a bitscore threshold in the `ko_list.txt` file (in the `.hmm` file) and their corresponding entries in from the `ko_list.txt` file (in `01_ko_fams_with_no_threshold.txt`). Please note that KOs from the `orphan_data` directory will *not* be annotated in your %(contigs-db)s when you run %(anvi-run-kegg-kofams)s. However, if you ever need to take a look at these profiles or use them in any way, here they are. :)
- `modules.keg`: a flat text file describing all metabolic modules available in the [KEGG MODULE](https://www.genome.jp/kegg/module.html) resource. This includes pathway and signature modules, but not reaction modules.
- The `modules` subfolder: contains flat text files, one for each metabolic module, downloaded using the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html). Each file describes a metabolic module's definition, classification, component orthologs, metabolic reactions, compounds, and any miscellaneous data like references and such. For an example, see the [module file for M00001](https://rest.kegg.jp/get/M00001/).
- `hierarchies.json`: a JSON-formatted file describing the available functional hierarchies in the [KEGG BRITE](https://www.genome.jp/kegg/brite.html) resource.
- The `BRITE` subfolder: contains JSON-formatted files, each one of which describes a BRITE hierarchy.
- `MODULES.db`: a SQLite database containing data parsed from the module files and BRITE hierarchies. See %(modules-db)s.

### How do we use it?

The KOfam profiles are used directly by %(anvi-run-kegg-kofams)s for annotating genes with KEGG Orthologs. The MODULE and BRITE data in the above files are processed and organized into the %(modules-db)s for easier programmatic access. %(anvi-run-kegg-kofams)s uses this database to annotate genes with BRITE categories and with the modules they participate in, when relevant. %(anvi-estimate-metabolism)s uses this database to get module information when computing completeness scores for each metabolic module.
