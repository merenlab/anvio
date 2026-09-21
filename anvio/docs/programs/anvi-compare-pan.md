This program compares two %(pan-db)s databases that were generated from the same %(genomes-storage-db)s.

It identifies gene clusters that differ in composition between the two pangenomes and classifies them as **fragmented** (a single gene cluster in one pan split into multiple in the other) or **combined** (multiple gene clusters in one pan merged into one in the other). It also classifies all gene clusters as **core**, **singleton**, or **accessory** based on genome occurrence, and adds these classifications as %(misc-data-items)s to both pan databases.

### Basic usage

{{ codestart }}
anvi-compare-pan --pan-db %(pan-db)s \
                 --compared-pan-db PAN_DB_2 \
                 --genomes-storage %(genomes-storage-db)s \
                 --output-file comparison-results.txt
{{ codestop }}

### When is this useful?

This program is most useful when you have generated two pangenomes from the same set of genomes using different parameters or methods (e.g., a sequence-based pangenome and a structure-informed pangenome) and want to understand how gene cluster assignments differ between them.

{:.notice}
Both pan databases must have been generated from the same %(genomes-storage-db)s. The program will refuse to compare pan databases with different genomes storage hashes.

### Output

The program adds comparison layers to the primary %(pan-db)s as %(misc-data-items)s, including fragmentation/combination status and gene cluster type classifications. If `--output-file` is provided, a tab-delimited text file with the comparison details is also generated.
