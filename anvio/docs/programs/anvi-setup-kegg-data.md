%(anvi-setup-kegg-data)s downloads and organizes data from KEGG for use by other programs, namely %(anvi-run-kegg-kofams)s, %(anvi-estimate-metabolism)s and %(anvi-reaction-network)s. Depending on what download mode you choose, it can download and setup one or more of the following:

- HMM profiles from the [KOfam](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907) database
  - with [an _(optional)_ step](#what-are-stray-kos-and-what-happens-when-i-include-them) to create new HMMs and estimate bit score thresholds for any 'stray KOs' without a KEGG-provided bit score threshold
- metabolic pathway information from [KEGG MODULES](https://www.genome.jp/kegg/module.html)
- functional classification information from [KEGG BRITE](https://www.genome.jp/kegg/brite.html)
- protein family information of the [KEGG Orthology database](https://www.genome.jp/kegg/ko.html)

 Typically, some processing is done following the data download to make the data work with downstream anvi'o programs. The KOfam profiles are prepared for later use by the HMMER software, and the information from MODULES and BRITE is made accessible to other anvi'o programs as a %(modules-db)s. The Orthology data is converted into a nice table that can be utilized by %(anvi-reaction-network)s. This program generates a directory with these files (%(kegg-data)s).

## Choosing a download mode

You need to pick a mode to work with this program to control which data will be downloaded from KEGG. You can see the available modes by running the following command:

{{ codestart }}
anvi-setup-kegg-data --list-modes
{{ codestop }}

You use the `--mode` parameter to tell the program which mode you want, for example:

{{ codestart }}
anvi-setup-kegg-data --mode modules
{{ codestop }}


## Default usage: downloading a KEGG snapshot

If you do not provide any arguments to this program, all KEGG data (ie, `--mode all`) will be set up in the default KEGG data directory.

{{ codestart }}
anvi-setup-kegg-data
{{ codestop }}

### How does it work?

By default, this program downloads a snapshot of the KEGG databases, already converted into an anvi'o-compatible format. The snapshot is a `.tar.gz` archive of a KEGG data directory that was (usually) generated around the time of the latest anvi'o release.

After the default KEGG archive is downloaded, it is unpacked, checked that all the expected files are present, and moved into the KEGG data directory.

### Why is this the default?

Doing it this way ensures that almost everyone uses the same version of KEGG data, which is good for reproducibility and makes it easy to share annotated datasets. The KEGG resources are updated fairly often, and we found that constantly keeping the KEGG data directory in sync with them was not ideal, because every time the data directory is updated, you have to update the KOfam annotations in all your contigs databases to keep them compatible with the current %(modules-db)s (unless you were smart enough to keep the old version of the KEGG data directory around somewhere). And of course that introduces a new nightmare as soon as you want to share datasets with your collaborators who do not have the same KEGG data directory version as you. With everyone using the same %(kegg-data)s by default, we can avoid these issues.

But the trade-off to this is that the default KEGG data version is tied to an anvi'o release, and it will not always include the most up-to-date information from KEGG. Luckily, **for those who want the most updated version of KEGG, you can still use this program to generate the KEGG data directory by downloading directly from KEGG** (see 'Getting the most up-to-date KEGG data' section below).

{:.warning}
BRITE hierarchy data is not included in the default KEGG snapshot for anvi'o `v7`. Starting from the `v7.1-dev` version of anvi'o, there is a new default KEGG snapshot including BRITE information. If you are missing this data, it can be acquired by either installing a later snapshot or by independently downloading it with this program using `--mode modules`.

{:.warning}
The data for metabolic modeling are not included in the KEGG snapshots created before anvi'o `v8`. If you are missing this data, it can be acquired by either installing a later snapshot or by independently downloading it with this program using `--mode modeling`.

### Set up KEGG data in a non-default location

You can specify a different directory in which to put this data, if you wish:

{{ codestart }}
anvi-setup-kegg-data --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

This is helpful if you don't have write access to the default directory location, or if you want to keep several different versions of the KEGG data on your computer. Just remember that when you want to use this specific KEGG data directory with later programs such as %(anvi-run-kegg-kofams)s, you will have to specify its location with the `--kegg-data-dir` flag.

### Setting up an earlier KEGG snapshot

By default, the KEGG snapshot that will be installed is the latest one, which is up-to-date with your current version of anvi'o. If, however, you want a snapshot from an earlier version, you can run something like the following to get it:

{{ codestart }}
anvi-setup-kegg-data --kegg-data-dir /path/to/directory/KEGG \
                       --kegg-snapshot v2020-04-27
{{ codestop }}

Just keep in mind that you may need to migrate the MODULES.db from these earlier versions in order to make it compatible with the current metabolism code. Anvi'o will tell you if you need to do this.

Not sure what KEGG snapshots are available for you to request? Well, you could check out the YAML file at `anvio/anvio/data/misc/KEGG-SNAPSHOTS.yaml` in your anvi'o directory, or you could just give something random to the `--kegg-snapshot` parameter and watch anvi'o freak out and tell you what is available:
{{ codestart }}
anvi-setup-kegg-data --kegg-snapshot hahaha
{{ codestop }}

Note that the latter method only shows you the date that each available snapshot was created. If you need more details about what types of data is included in each snapshot, you should look at the YAML file, which annotates each snapshot with a bit more detail. For example, the following entry does not contain metabolic modeling data (obsolete as of anvi'o `v8.0-dev`) OR models/thresholds for 'stray KOs':

```
v2023-09-18:
    url: https://figshare.com/ndownloader/files/42381873
    archive_name: KEGG_build_2023-09-18_a2b5bde358bb.tar.gz
    hash: a2b5bde358bb
    modules_db_version: 4
    no_modeling_data: True
    no_stray_KOs: True
```

## Getting the most up-to-date KEGG data: downloading directly from KEGG

This program is also capable of downloading data directly from KEGG and converting it into an anvi'o-compatible format. In fact, this is how we generate the default KEGG archive. If you want the latest KEGG data instead of the default snapshot of KEGG, try the following:

{{ codestart }}
anvi-setup-kegg-data --download-from-kegg
{{ codestop }}

Please note that this will download all the KEGG data (ie, `--mode all` is the default). If you want to independently download individual KEGG datasets, you should pick one of the other modes (the `--download-from-kegg` flag is implicitly turned on in these modes).

### How does it work?

KOfam profiles are downloadable from KEGG's [FTP site](ftp://ftp.genome.jp/pub/db/kofam/) and all other KEGG data is accessible as flat text files through their [API](https://www.kegg.jp/kegg/rest/keggapi.html). When you run this program it will first get all the files that it needs from these sources, and then it will process them by doing the following:

- determine if any KOfam profiles are missing bitscore thresholds (we call these "stray KOs"), and remove those from the standard profile location so that they are not used for annotation (if you want to see these, you will find them in the `orphan_data` folder in your KEGG data directory)
- concatenate all remaining KOfam profiles into one file and run `hmmpress` on them
- if the `--include-stray-KOs` flag is used, process the stray KO models as detailed in [this section](#what-are-stray-kos-and-what-happens-when-i-include-them)
- parse the flat text file for each KEGG module and the JSON file for each BRITE hierarchy
- store the MODULE and BRITE information in the %(modules-db)s
- parse the flat text files from KEGG Orthology and organize these into a table for metabolic modeling

An important thing to note about this option is that it has rigid expectations for the format of the KEGG data that it works with. Future updates to KEGG may break things such that the data can no longer be directly obtained from KEGG or properly processed. In the sad event that this happens, you will have to download KEGG from one of our archives instead.

### What are 'stray KOs' and what happens when I include them?

"Stray KOs" are what we call KOfam models that don't come with a bit score threshold defined by KEGG. As described in the [Aramaki et al. 2020 paper](https://doi.org/10.1093/bioinformatics/btz859), KEGG Ortholog families that have _less than 3 non-redundant sequences_ don't go through the workflow for estimating bit score thresholds. We normally avoid using these models for annotation, because without a bit score threshold, we don't have a way to distinguish between good and bad matches to the model. This is why these KOfams are removed from the set of regular profiles that is used by %(anvi-run-kegg-kofams)s.

However, in early 2024 we came up with a strategy to estimate bit score thresholds for these models so that they could be used for annotation and downstream analyses via the `--include-stray-KOs` flag. If you add this flag to your `anvi-setup-kegg-data` command, this is what will happen to each stray KO that doesn't yet have a bit score threshold:

- download its KEGG Orthology file and parse it to obtain a list of KEGG GENE sequences that belong to this KO family
- download the KEGG GENE file for each of these and parse them to get the amino acid sequences
- if more than one KEGG GENE sequence belongs to this family, we align the sequences with `muscle` and run `hmmbuild` to generate a _new_ Hidden Markov Model for this KO (named by adding `"_anvio_version"` to the end of the KO accession)
- otherwise, we use KEGG's original KOfam model for this KO
- run `hmmscan` of the KO's model (either the new one or the original) against the KEGG GENE sequences belonging to this family to get a distribution of 'positive' matches to the model
- take the minimum bit score of these 'positive' matches to use as the bit score threshold for the model

When this is complete for every KO, the program saves the estimated bit score thresholds to a tab-delimited file, located in your KEGG data directory at `orphan_data/estimated_thresholds_for_stray_kos.txt`. It also concatenates the relevant model (new or original) for each KO family into one file (at `orphan_data/anvio_hmm_profiles_for_stray_KOs.hmm`), and runs `hmmpress` on it so that these models can be used for annotation downstream (if requested).

Astute readers will notice that this process is similar to half of the workflow that KEGG uses to estimate bit scores on the larger KO families; that is, we use the 'positive' matches of the model against the gene sequences used to create it to find the minimum bit score that is required to annotate these orthologous sequences. For the sake of time and computational resources, we skip [KEGG's more rigorous process](https://doi.org/10.1093/bioinformatics/btz859) of dividing the 'positive' set into 3 random groups (which cannot be done anyway with these small families, otherwise KEGG would have done it for us already) and of including a 'negative' set of sequences that aren't orthologous to refine the bit score threshold.

Our resulting bit score thresholds are not perfect, but they are rather conservative and should mostly avoid the introduction of garbage annotations with these models downstream. Still, we always encourage caution and oversight of your data/results, especially if you decide to use these downstream (which you can do by using the `--include-stray-KOs` flag with %(anvi-run-kegg-kofams)s and %(anvi-estimate-metabolism)s). :)

{:.notice}
If you want to download a snapshot of KEGG data that includes these stray KOs processed and ready to use for annotation and beyond, you should check the `anvio/anvio/data/misc/KEGG-SNAPSHOTS.yaml` file in your anvi'o directory to see which ones have this data. Look for entries that _do not_ have `no_stray_KOs: True`.


**Why do we create new HMMs for some of these KOs?**

While testing our approach of setting the minimum 'positive' bit score to be a model's threshold, we realized that the KEGG ORTHOLOGY and KEGG GENES databases are updated far more often than the KOfam models. This causes some KOfam models to be out-of-date relative to the family of orthologs that the model is meant to represent -- that is, a KOfam profile doesn't necessarily include all of the sequences from KEGG GENES that are assigned to this KO family. This violates our assumption that every gene sequence in the family is a 'positive' match to the model. Newer sequences can still match to the model, but often with rather low bit scores, which then become the minimum score that is used to set the model's threshold. And those low thresholds can cause false positive annotations.

To circumvent this issue, we decided to remake models for any stray KO that could have sequences that aren't yet incorporated into the KO's HMM; that is, any stray KO with more than one gene sequence in the family. By doing this, we ensure that each gene sequence actually represents a 'positive' match to the model, and we can estimate thresholds that are conservative enough to avoid garbage annotations.

We clearly name the new HMMs by appending `"_anvio_version"` to the end of the KO identifier; for example, `K21692_anvio_version`. These names will be used for any annotations taken from these modified profiles. So if you don't see `"_anvio_version"` at the end of the function accession, you can be confident that the original KEGG KOfam profile was used.

### The --only-download option

The `--only-download` flag works for `KOfam` mode and `modules` mode.

Suppose you only want to download data from KEGG without processing it. For instance, perhaps you don't need a %(modules-db)s or you don't want `hmmpress` to be run on the KOfam profiles. You can instruct this program to stop after downloading by providing the `--only-download` flag:

{{ codestart }}
anvi-setup-kegg-data --mode modules \
                       --only-download \
                       --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

It's probably a good idea in this case to specify where you want this data to go using `--kegg-data-dir`, to make sure you can find it later.

{:.notice}
This option is primarily useful for developers to test `anvi-setup-kegg-data` - for instance, so that you can download the data once and run the database setup option (`--only-processing`) multiple times. However, if non-developers find another practical use-case for this flag, we'd be happy to add those ideas here. Send us a message, or feel free to edit this file and pull request your changes on the anvi'o Github repository. :)

### The --only-processing option

The `--only-processing` flag works for `KOfam` mode and `modules` mode.

Let's say you already have KEGG data on your computer that you got by running this program with the `--only-download` flag. Now you want to process the HMM files, or turn the MODULES data into a %(modules-db)s. To do that, run this program using the `--only-processing` flag and provide the location of the pre-downloaded KEGG data:

{{ codestart }}
anvi-setup-kegg-data --mode modules \
                       --only-processing \
                       --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

{:.notice}
The KEGG data that you already have on your computer has to be in the format expected by this program, or you'll run into errors. Pretty much the only reasonable way to get the data into the proper format is to run this program with the `--only-download` option. Otherwise you would have to go through a lot of manual file-changing shenanigans - possible, but not advisable.

One more note: since this flag is most often used for testing the database setup capabilities of this program, which entails running `anvi-setup-kegg-data --mode modules --only-processing` multiple times on the same KEGG data directory, there is an additional flag that may be useful in this context. To avoid having to manually delete the created modules database each time you run, you can use the `--overwrite-output-destinations` flag:

{{ codestart }}
anvi-setup-kegg-data --mode modules \
                       --only-processing \
                       --kegg-data-dir /path/to/directory/KEGG \
                       --overwrite-output-destinations
{{ codestop }}

### Avoiding files used in later versions of anvi'o

As anvi'o has expanded, new types of KEGG files have been included in the data pack. We highly recommend including these files, as this program does by default, when downloading and setting up KEGG data. However, options are available to skip these files -- perhaps you're a rebel, or adverse to change, or something is not working on your computer... whatever, all good by us.

It should make sense to you that these flags do not work when setting up from a KEGG snapshot that already includes the newer types of files.

#### BRITE hierarchies

As of anvi'o `v7.1-dev` or later, [KEGG BRITE hierarchies](https://www.genome.jp/kegg/brite.html) are added to the %(modules-db)s when running this program with `--mode modules`. These hierarchies are especially useful for classifying and making sense of KOs and other KEGG data. To avoid these files, use the `--skip-brite-hierarchies` flag:

{{ codestart }}
anvi-setup-kegg-data --mode modules --skip-brite-hierarchies
{{ codestop }}

#### Reaction network

As of anvi'o `v8.0-dev`, [KEGG binary relations files](https://www.genome.jp/brite/br08906) are included in KEGG data when running this program with `--mode modules`. These are needed for %(reaction-network)s construction by %(anvi-reaction-network)s, which is how biochemical data associated with KOs is now interpreted and analyzed. To avoid these files, use the `--skip-binary-relations` flag:

{{ codestart }}
anvi-setup-kegg-data --mode modules --skip-binary-relations
{{ codestop }}

#### Pathway visualization

As of anvi'o `v8.0-dev`, [KEGG pathway map](https://www.genome.jp/kegg/pathway.html) PNG and KGML files are included in KEGG data when running this program with `--mode modules`. These are needed by %(anvi-draw-kegg-pathways)s to visualize anvi'o data in the context of biochemical pathways. To avoid these files, use the `--skip-map-images` flag:

{{ codestart }}
anvi-setup-kegg-data --mode modules --skip-map-images
{{ codestop }}

### How do I share this data?
Suppose you have been living on the edge and annotating your contigs databases with a non-default version of %(kegg-data)s, and you share these databases with a collaborator who wants to run downstream programs like %(anvi-estimate-metabolism)s on them. Your collaborator (who has a different version of %(kegg-data)s on their computer) will likely get version errors as detailed on the %(anvi-estimate-metabolism)s help page.

In order for your collaborator to be able to work with your dataset, they need to have the same %(kegg-data)s version as you did when you ran %(anvi-run-kegg-kofams)s. If you are very lucky and KEGG has not been updated since you set up your %(kegg-data)s, they may be able to run `anvi-setup-kegg-data -D` to get it. But if not, there are a few options for you to share your version of %(kegg-data)s:

1. You could send them your KEGG data directory. First, run `tar -czvf kegg_archive.tar.gz ./KEGG` on the data directory to compress and archive it before sending it over (this command _must_ be run from its parent directory so that the archive has the expected directory structure when it is unpacked). Then your collaborator can just run `anvi-setup-kegg-data --kegg-archive kegg_archive.tar.gz --kegg-data-dir ./KEGG_ARCHIVE` and be good to go. They would just have to use `--kegg-data-dir ./KEGG_ARCHIVE` when running downstream programs. The problem here is that even the archived %(kegg-data)s is quite large, ~4-5GB, and may be unfeasible for you to send.
2. You could share with your collaborator just the %(modules-db)s. If all they want to do is to run %(anvi-estimate-metabolism)s on databases annotated by your version of the KEGG data directory, this should be all they need. They would need to pass the folder containing your %(modules-db)s to %(anvi-estimate-metabolism)s using the `--kegg-data-dir` parameter.
3. If your collaborator also wants to be able to annotate other databases with your version of %(kegg-data)s, then they need to have the KOfam profiles as well. You can send them your %(modules-db)s and have them download the KOfam profiles most similar to the ones you have from the [KOfam archives](https://www.genome.jp/ftp/db/kofam/archives/) (which are labeled by date). Then they would have to essentially construct their own KEGG data directory by copying the structure of the default one and putting the downloaded files (and the %(modules-db)s you sent them) into the correct locations. The KOfam profiles must be concatenated into a `Kofam.hmm` file and `hmmpress` must be run on that file to generate the required indices for `hmmsearch`. Your collaborator must also have the `ko_list.txt` file (which _should_ be downloaded with the profiles) in the right spot. Then they could pass their makeshift KEGG data directory to %(anvi-run-kegg-kofams)s using `--kegg-data-dir`, and they should be golden. (A word of warning: they may want to remove KOs without bitscore thresholds in the `ko_list.txt` before concatenating the profiles, otherwise they will likely get a lot of weak hits for these KOs.)

## I already have a KEGG snapshot: set up from a pre-downloaded archive file

If you have an archive (`.tar.gz`) of the KEGG data directory already on your computer (perhaps a colleague or Meren Lab developer gave you one), you can set up KEGG from this archive instead:

{{ codestart }}
anvi-setup-kegg-data --kegg-archive KEGG_archive.tar.gz
{{ codestop }}

This works the same way as the default, except that it bypasses the download step and instead uses the archive file you have provided with `--kegg-archive`.

## Info for developers: making a new KEGG snapshot available to all anvi'o users

Periodically (especially before releasing a new version of anvi'o), we want to add new KEGG database snapshots to anvi'o so that users can have more up-to-date KEGG data without having to use the `--download-from-kegg` option. In this section you will find the instructions for doing this (these instructions are also in the comments of the `anvio/data/misc/KEGG-SNAPSHOTS.yaml` file).

Available KEGG snapshots are stored in the anvi'o code repository in `anvio/data/misc/KEGG-SNAPSHOTS.yaml`. To add a new snapshot, you first need to create one by downloading and processing the data from KEGG, testing to make sure it works, and then updating this file. Here are the steps:

1. Download the latest data directly from KEGG by running `anvi-setup-kegg-data -D --kegg-data-dir ./KEGG -T 5`. This will create the new KEGG data folder with its %(modules-db)s in your current working directory. Make sure you use the exact folder name of `./KEGG`, because that is what anvi'o expects to find when it unpacks a KEGG snapshot. You may want to reduce or increase the number of threads (`-T`) according to your available compute resources.
2. Get the hash value and version info from the MODULES.db by running `anvi-db-info ./KEGG/MODULES.db`.
3. Archive the KEGG data directory by running `tar -czvf KEGG_build_YYYY-MM-DD_HASH.tar.gz ./KEGG`. Please remember to replace YYYY-MM-DD with the current date and replace HASH with the MODULES.db hash value obtained in step 2. This convention makes it easier to distinguish between KEGG snapshots by simply looking at the file name.
4. Test that setup works with this archive by running `anvi-setup-kegg-data --kegg-archive KEGG_build_YYYY-MM-DD_HASH.tar.gz --kegg-data-dir TEST_NEW_KEGG_ARCHIVE`.
5. If setup worked in the last step without errors, upload the `.tar.gz` archive to [Figshare](https://figshare.com/). If you need inspiration for filling out the keywords, categories, and description fields for the archive, you can check the previous KEGG snapshots that have been uploaded - for instance, [this one](https://figshare.com/articles/dataset/KEGG_build_2023-01-10/21862494) or [this one](https://figshare.com/articles/dataset/KEGG_build_2022-04-14/19601761). At minimum, we typically indicate the database version and hash value, and an example setup command (ie, the one from step 4), in the description of the dataset. Once the archive is published on Figshare (warning: this usually takes a while due to the large file size), you can get the download url of the archive by right-clicking on the Download button and copying the address, which should be a URL with a format similar to this example (but different numbers): `https://figshare.com/ndownloader/files/34817812`
6. Add an entry to the bottom of the `anvio/data/misc/KEGG-SNAPSHOTS.yaml` file with the Figshare download URL, archive name, and MODULES.db hash and version. If you want this to become the default snapshot (which usually only changes before the next anvi'o release), you should also update the default `self.target_snapshot` variable in `anvio/kegg.py` to be this latest version that you have added.
7. Test it by running `anvi-setup-kegg-data --kegg-data-dir TEST_NEW_KEGG`, and if it works you are done, and can push your changes to the anvi'o repository. :)

## Downloading generic KEGG data in Python

If you want to get some data from the KEGG website that is not included in our default download (or, if you only want a subset of that data without going through the whole setup process), you can use the anvi'o API to utilize our download functions. Here are some examples for using the `KeggSetup` class (for example, in the Python interpreter):

### Loading the `KeggSetup` class

`KeggSetup` is the class for downloading KEGG data (using KEGG's API). To use it in Python, you need to load the `kegg` module from anvi'o. When using it this way, we recommend skipping a variety of sanity checks using the `skip_init` parameter - this is mainly so that the class doesn't check for, remove, or complain about existing KEGG data on your computer.

```python
import anvio
import argparse
from anvio import kegg
args = argparse.Namespace(reset=False)
setup = kegg.KeggSetup(args, skip_init=True)
```

Once you have this class loaded, you can use its functions for a variety of download and processing tasks. We'll show some examples below.

### Downloading all flat files associated with a KEGG hierarchy

 The following example demonstrates the download of all KEGG COMPOUND files belonging to the BRITE hierarchy with accession `br08001`. Note that if you do not specify a download directory, the files will by default be downloaded to the current working directory.

 ```python
setup.download_kegg_files_from_hierarchy('br08001', download_dir='KEGG_COMPOUND')
 ```

 ### Downloading a hierarchical text file

 If you just want to get a KEGG `htext` file (with extension `.keg`), use the following function:

 ```python
setup.download_generic_htext('br08001', download_dir='KEGG_COMPOUND')
 ```

 ### Processing a hierarchical text file

 We have a few functions for reading KEGG's `htext` files. If all you want is a list of the accessions involved in this heirarchy (for instance, all compounds in a BRITE hierarchy for KEGG COMPOUND), use this one (the argument should be the path to the `htext` file):

 ```python
accession_list = setup.get_accessions_from_htext_file("br08001.keg")
 ```

 If you want to process the KEGG module `htext` file to get a dictionary of all modules and their names/classes/etc, use the following code. You will need to set the `kegg_module_file` attribute (of the ModulesDownload class) to point to the location of the `modules.keg` file, and the function will store the module dictionary in the `module_dict` attribute.

 ```python
modules_setup = kegg.ModulesDownload(args)
modules_setup.kegg_module_file = "modules.keg"
modules_setup.process_module_file()
modules_setup.module_dict  # this attribute now stores the module dictionary
 ```

 ### Downloading a flat file using the KEGG API

 Here is a wrapper function that will 'get' a flat file with the KEGG API. You can provide this function with the accession of the data you want (for instance, a module accession), and optionally a directory to download it into.

 ```python
setup.download_generic_flat_file('C00058',  download_dir='KEGG_COMPOUND')
 ```
