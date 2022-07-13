This program **creates a new %(collection)s from the %(bin)ss in another collection with specific guidelines.** This is especially helpful when you wish to standardize your bin names, add project specific prefixes, and/or exclude those that do not match your criteria of completion, redundancy, and/or size estimates.

### Renaming all bins in a collection

Let's say you have a %(collection)s called `MY_COLLECTION`, which has four bins that are named poorly (which can happen due to decisions made by automatic binning tools, or after a few steps of manual refinement): `Bin_1_2_1`, `Bin_2`, `Bin_3_1_1`, and `Bin_4`. In an instance like this, running the program %(anvi-rename-bins)s the following way will standardize these bin names with a prefix specific to your project:

{{ codestart }}
anvi-rename-bins -c %(contigs-db)s \
                 -p %(profile-db)s \
                 --prefix SURFACE_OCEAN \
                 --collection-to-read MY_COLLECTION \
                 --collection-to-write SURFACE_OCEAN_SAMPLES \
                 --report-file rename.txt
{{ codestop }}

Now your %(profile-db)s will have a new collection named `SURFACE_OCEAN_SAMPLES` that will contains your four bins witht their new names `SURFACE_OCEAN_Bin_00001`, `SURFACE_OCEAN_Bin_00002`, `SURFACE_OCEAN_Bin_00003`, and `SURFACE_OCEAN_Bin_00004`. The new naming will order your bins based on their substantive completion (i.e., completion minus redunancy).

The file `rename.txt` is a TAB-delimited file that contains a summary of your renaming process. The first column has the original name of the bins that you renamed, the second has their new names, and the remaining columns contain information about those bins (like their completion, redundency, and size).

### Separating out the MAGs

You can also label your MAGs separately from your bins via the flag `--call-MAGs`:

{{ codestart }}
anvi-rename-bins -c %(contigs-db)s \
                 -p %(profile-db)s \
                 --prefix SURFACE_OCEAN \
                 --collection-to-read MY_COLLECTION \
                 --collection-to-write SURFACE_OCEAN_MAGS \
                 --report-file rename.txt \
                 --call-MAGs \
                 --min-completion-for-MAG 70
{{ codestop }}

Now, the %(collection)s `SURFACE_OCEAN_MAGS` will include  `SURFACE_OCEAN_MAG_00001`, `SURFACE_OCEAN_MAG_00002`, `SURFACE_OCEAN_MAG_00003`, and `SURFACE_OCEAN_Bin_00004`. These are exactly the same bins that the collection contained before, but now the names differenciate the wheat from the chaff.

In addition to minimum completion estimate, you can also adjust the maximum redundancy value, minimum size to call MAGs. Please see the help menu for all parameters and their descriptions. 

### Exclude bins that are not MAGs

When you use the flag `--call-MAGs`, anvi'o identifies those bins that could be considered 'MAGs' based on your specific criteria. But regardles of whether an original bin remains a bin, or tagged as a MAG, everything in your original collection will end up in your new collection. The flag `--exclude-bins` enable you to filter out those that end up not being tagged as MAGs:

{{ codestart }}
anvi-rename-bins -c %(contigs-db)s \
                 -p %(profile-db)s \
                 --prefix SURFACE_OCEAN \
                 --collection-to-read MY_COLLECTION \
                 --collection-to-write SURFACE_OCEAN_MAGS \
                 --report-file rename.txt \
                 --min-completion-for-MAG 70 \
                 --call-MAGs \
                 --exclude-bins
{{ codestop }}

With the addition of the flag `--exclude-bins` to the same command, the %(collection)s `SURFACE_OCEAN_MAGS` will no longer include %(bin)ss `SURFACE_OCEAN_Bin_00003` and `SURFACE_OCEAN_Bin_00004`.

See also the program %(anvi-delete-collection)s.

### The report file

Following is an example reporting output file anvi'o will generate at the file path declared with the parameter `--report-file`:

|**old_bin_name**|**new_bin_name**|**SCG_domain**|**completion**|**redundancy**|**size_in_Mbp**|
|:--|:--|:--|:--|:--|:--|
|Bin_2|p800_MAG_00001|eukarya|61.45|7.23|26.924911|
|Bin_1|p800_MAG_00002|bacteria|98.59|8.45|1.612349|
|Bin_3|p800_Bin_00003|blank|0.00|0.00|0.103694|
|Bin_5|p800_Bin_00004|blank|0.00|0.00|0.128382|
|Bin_4|p800_Bin_00005|bacteria|1.41|0.00|0.378418|

The column `SCG_domain` will explain which collection of single-copy core genes were used to generate these completion/redundancy estimates. The absence of any domain prediction for any given bin will be marked with the keyrowd `blank`.
