%(anvi-setup-kegg-kofams)s downloads and organizes data from KEGG for use by other programs, namely %(anvi-run-kegg-kofams)s and %(anvi-estimate-metabolism)s. It downloads HMM profiles from the [KOfam](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907) database as well as metabolism information such as that stored in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). The KOfam profiles are prepared for later use by the HMMER software, and the metabolism information is made accessible to other anvi'o programs as a %(modules-db)s. This program generates a directory with these files (%(kegg-data)s), which by default is located at `anvio/anvio/data/misc/KEGG/`.

### Default Usage

If you do not provide any arguments to this program, the KOfam profiles and KEGG information will be set up in the default KEGG data directory.

{{ codestart }}
anvi-setup-kegg-kofams
{{ codestop }}

**How does it work?**
By default, this program downloads a snapshot of the KEGG databases, already converted into an anvi'o-compatible format. The snapshot is a `.tar.gz` archive of a KEGG data directory that was generated around the time of the latest anvi'o release.  

After the default KEGG archive is downloaded, it is unpacked, checked that all the expected files are present, and moved into the KEGG data directory.

Doing it this way ensures that almost everyone uses the same version of KEGG data, which is good for reproducibility and makes it easy to share annotated datasets. The KEGG resources are updated fairly often, and we found that constantly keeping the KEGG data directory in sync with them was not ideal, because every time the data directory is updated, you have to update the KOfam annotations in all your contigs databases to keep them compatible with the current %(modules-db)s (unless you were smart enough to keep the old version of the KEGG data directory around somewhere). And of course that introduces a new nightmare as soon as you want to share datasets with your collaborators who do not have the same KEGG data directory version as you. With everyone using the same %(kegg-data)s by default, we can avoid these issues.

But the trade-off to this is that the default KEGG data version is tied to an anvi'o release, and it will not always include the most up-to-date information from KEGG. Luckily, for those who want the most updated version of KEGG, you can still use this program to generate the KEGG data directory by downloading directly from KEGG (see 'Generating an anvi'o compatible KEGG data directory from scratch' below).

### How to set up KEGG data in a non-default location

You can specify a different directory in which to put this data, if you wish:

{{ codestart }}
anvi-setup-kegg-kofams --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

This is helpful if you don't have write access to the default directory location, or if you want to keep several different versions of the KEGG data on your computer. Just remember that when you want to use this specific KEGG data directory with later programs such as %(anvi-run-kegg-kofams)s, you will have to specify its location with the `--kegg-data-dir` flag.

### Generating an anvi'o compatible KEGG data directory from scratch

This program is also capable of downloading data directly from KEGG and converting it into an anvi'o-compatible format. In fact, this is how we generate the default KEGG archive. If you want the latest KEGG data instead of the default snapshot of KEGG, try the following:

{{ codestart }}
anvi-setup-kegg-kofams --download-from-kegg
{{ codestop }}

**How does it work?**
KOfam profiles are downloadable from KEGG's [FTP site](ftp://ftp.genome.jp/pub/db/kofam/) and all other KEGG data is accessible as flat text files through their [API](https://www.kegg.jp/kegg/rest/keggapi.html). When you run this program it will first get all the files that it needs from these sources, and then it will process them by doing the following:

- determine if any KOfam profiles are missing bitscore thresholds, and remove those from the standard profile location so that they are not used for annotation (if you want to see these, you will find them in the `orphan_data` folder in your KEGG data directory)
- concatenate all remaining KOfam profiles into one file and run `hmmpress` on them
- parse the flat text file for each KEGG module and store the information into the %(modules-db)s

An important thing to note about this option is that it has rigid expectations for the format of the KEGG data that it works with. Future updates to KEGG may break things such that the data can no longer be directly obtained from KEGG or properly processed. In the sad event that this happens, you will have to download KEGG from one of our archives instead.

### Set up from archived KEGG data

If you have an archive (`.tar.gz`) of the KEGG data directory already on your computer (perhaps a colleague or Meren Lab developer gave you one), you can set up KEGG from this archive instead:

{{ codestart }}
anvi-setup-kegg-kofams --kegg-archive KEGG_archive.tar.gz
{{ codestop }}

This works the same way as the default, except that it bypasses the download step and instead uses the archive you have provided with `--kegg-archive`.
