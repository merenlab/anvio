This is a two-column TAB-delimited file without a header that describes a %(collection)s by associating items with bin names. It can be used to import or export collections in and out of anvi'o databases, and/or transferring them between anvi'o projects. 

The first column in the file lists item names and the second column associates a given item with a bin. 

{{ codestart }}
item_01    bin_1
item_02    bin_1
item_03    bin_1
item_04    bin_2
item_05    bin_3
item_06    bin_3
{{ codestop }}

### The optinal bins info file

In addition to the essential file above, you can associate an optional TAB-delmited file with three columns with a collection to provide information about 'bins' in it, such as their source, and/or color to be used when they are displayed in %(summary)s outputs or anvi'o %(interactive)s interfaces. Here is an example:

```
bin_1	N/A	 #c9d433
bin_2	N/A	 #e86548
bin_3	N/A	 #0b8500
```

In this file format the first column is a bin name, the second column is a source, and the third column is an HTML color.

{:.notice}
The source is a free form text and can be anything. We often use `anvi-interactive` or `CONCOCT` or `anvi-refine` for our bins to track which ones were manually refined, and which ones were coming from an automated binning algorithm.

You can provide this optional file to the program %(anvi-import-collection)s with the parameter `--bins-info`.