This program lets you export miscellaneous data of your choosing into a text file, which can be imported into another anvi'o project using %(anvi-import-misc-data)s. You can export the same types of data that you can import with that function. These are also listed below.

This program only works on data that is listed as an available key (most often because it was previously imported by the user). To view available keys, call

{{ codestart }}
anvi-export-misc-data -p %(profile-db)s \
                      --target-data-table items|layers|layer_orders \
                      --list-available-keys
{{ codestop }}

where you choose the appropriate option for the taget-data-table. 

## Data types you can export 

### From a %(pan-db)s or %(profile-db)s, you can export 

- items data (%(misc-data-items)s) into a %(misc-data-items-txt)s. 

{{ codestart }}
anvi-export-misc-data -p %(profile-db)s \
                      --target-data-table items \
                      --keys-to-remove key_1
{{ codestop }}

- layers data (%(misc-data-layers)s) into a %(misc-data-layers-txt)s.  

{{ codestart }}
anvi-export-misc-data -p %(pan-db)s \
                      --target-data-table layers \
                      --keys-to-remove key_2,key_3
{{ codestop }}

- layer orders data (%(misc-data-layer-orders)s) into a %(misc-data-layer-orders-txt)s. 

{{ codestart }}
anvi-export-misc-data -p %(profile-db)s \
                      --target-data-table layer_orders \
                      --keys-to-remove key_4
{{ codestop }}

#### Additional notes

If your misc-data-table is associated with a specific data group, you can provide that data group to this program with the `-D` flag. 

### From a %(contigs-db)s, you can export 

- nucleotide data (%(misc-data-nucleotides)s) into a %(misc-data-nucleotides-txt)s.

{{ codestart }}
anvi-export-misc-data -c %(contigs-db)s \
                      --keys-to-remove key_1
{{ codestop }}

- amino acid data (%(misc-data-amino-acids)s) into a %(misc-data-amino-acids-txt)s.
