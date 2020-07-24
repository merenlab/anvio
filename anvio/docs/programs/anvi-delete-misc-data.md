After you've added misc-data of some kind (%(misc-data-items)s, %(misc-data-layers)s, %(misc-data-layer-orders)s, %(misc-data-nucleotides)s, or %(misc-data-amino-acids)s) using %(anvi-import-misc-data)s, you can **delete that data and remove it from the interactive interface** using this program. 

This program will release your data into the either, never to be seen again. If you would like to first export it into a text file (so that it can be seen again), you can do so with %(anvi-export-misc-data)s. 

This program only works on data that is listed as an available key (most often because it was previously imported by the user). To view available keys, call

{{ codestart }}
anvi-delete-misc-data -p %(profile-db)s \
                      --target-data-table items|layers|layer_orders \
                      --list-available-keys
{{ codestop }}

where you choose the appropriate option for the taget-data-table. 

## Data types you can delete 

### From a %(pan-db)s or %(profile-db)s, you can delete 

- items data (%(misc-data-items)s) 

{{ codestart }}
anvi-delete-misc-data -p %(profile-db)s \
                      --target-data-table items \
                      --keys-to-remove key_1
{{ codestop }}

- layers data (%(misc-data-layers)s)

{{ codestart }}
anvi-delete-misc-data -p %(pan-db)s \
                      --target-data-table layers \
                      --keys-to-remove key_2,key_3
{{ codestop }}

- layer orders data (%(misc-data-layer-orders)s)

{{ codestart }}
anvi-delete-misc-data -p %(profile-db)s \
                      --target-data-table layer_orders \
                      --keys-to-remove key_4
{{ codestop }}

#### Additional notes

If your misc-data-table is associated with a specific data group, you can provide that data group to this program with the `-D` flag. 

### From a %(contigs-db)s, you can export 

- nucleotide data (%(misc-data-nucleotides)s)

{{ codestart }}
anvi-delete-misc-data -c %(contigs-db)s \
                      --keys-to-remove key_1
{{ codestop }}

- amino acid data (%(misc-data-amino-acids)s)
