This is a 2-column TAB-delimited text file to associate a given set of items with a set of groups. Depending on the context, items here may be individual samples or genomes. The first column can have a header name `item`, `sample`, `genome` or anything else that is appropriate, and list the items that are relevant to your input data. The second column should have the header `group`, and associate each item in your data with a group.

Each item should be associated with a single group, and it is always a good idea to define groups using single words without any fancy characters. For instance, `HIGH_TEMPERATURE` or `LOW_FITNESS` are good group names. In contrast, `my group #1` or `IS-THIS-OK?`, are not quite good names for groups and may cause issues downstream depending on who uses this file.

Here is an example %(groups-txt)s file:

|item|group|
|:--|:--|
|item_01|GROUP_A|
|item_02|GROUP_B|
|item_03|GROUP_A|
|(...)|(...)|

{:.warning}
If you are passing this file to the program %(anvi-compute-metabolic-enrichment)s, the names in the `sample` column must match those in the "modules" mode output file that you provide to the program via the `--modules-txt` parameter. If you know that the sample names match but you are still getting errors, you might need to specify which column in the "modules" mode output contains those sample names using the `--sample-header` parameter.
