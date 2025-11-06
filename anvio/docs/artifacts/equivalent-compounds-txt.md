This file describes pairs of ModelSEED compound IDs that should be considered equivalent to each other for the purposes of predicting potential metabolic exchanges with %(anvi-predict-metabolic-exchanges)s.

## Example file

The file should be tab-delimited and have four columns. The first two columns contain the pair of ModelSEED compound IDs that you consider equivalent, and the second two columns contain their respective human-readable names (for convenience when reading the file).

|**`compound_id`**|**`equivalent_id`**|**`name`**|**`equivalent_name`**|
|:--|:--|:--|:--|
|cpd00035|cpd01003|L-Alanine|Alanine|
|cpd00039|cpd19182|L-Lysine|Lysine|
|cpd00041|cpd19181|L-Aspartate|Aspartate|

## Pro Tip: generating the file automatically

If you run %(anvi-predict-metabolic-exchanges)s with the `--use-equivalent-amino-acids` flag, it will create a file of this type containing the set of amino acid compound equivalents it automatically finds in the ModelSEED database:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --use-equivalent-amino-acids
{{ codestop }}

You can modify the resulting output file as you want to remove or add new compound equivalents and make a custom set that can be passed to %(anvi-predict-metabolic-exchanges)s with the `--custom-equivalent-compounds-file` parameter.