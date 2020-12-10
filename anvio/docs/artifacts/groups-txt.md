This is a two-column, TAB-delimited file that describes which samples belong to which groups. It is a required input for computing enrichment scores of metabolic modules in %(anvi-compute-enrichment-scores)s.

Here is an example:

|sample|group|
|:--|:--|
|sample_01|GROUP-A|
|sample_02|GROUP-B|
|sample_03|GROUP-A|
|(...)|(...)|

{:.warning}
The names in the `sample` column must match those in the "modules" mode output file that you provide to the %(anvi-compute-enrichment-scores)s program. If you know that the sample names match but you are still getting errors, you might need to specify which column in the "modules" mode output contains those sample names using the `--sample-header` parameter of %(anvi-compute-enrichment-scores)s. 
