This program lets you export miscellaneous data of your choosing into a text file, which can be imported into another anvi'o project using %(anvi-import-misc-data)s. You can export the same types of data that you can import with that function. These are also listed below.

To see what misc-data is available in your database, use %(anvi-show-misc-data)s. 

If your misc-data is associated with a specific data group, you can provide that data group to this program with the `-D` flag. 

## Data types you can export 

### From a pan-db or profile-db: items, layers, layer orders

**From a %(pan-db)s or %(profile-db)s, you can export**

- items data (%(misc-data-items)s) into a %(misc-data-items-txt)s. 

{{ codestart }}
anvi-export-misc-data -p %(profile-db)s \
                      --target-data-table items 
{{ codestop }}

- layers data (%(misc-data-layers)s) into a %(misc-data-layers-txt)s.  

{{ codestart }}
anvi-export-misc-data -p %(pan-db)s \
                      --target-data-table layers 
{{ codestop }}

- layer orders data (%(misc-data-layer-orders)s) into a %(misc-data-layer-orders-txt)s. 

{{ codestart }}
anvi-export-misc-data -p %(profile-db)s \
                      --target-data-table layer_orders 
{{ codestop }}

### From a contigs-db: nucleotide and amino acid information

**From a %(contigs-db)s, you can export**

- nucleotide data (%(misc-data-nucleotides)s) into a %(misc-data-nucleotides-txt)s.

{{ codestart }}
anvi-export-misc-data -c %(contigs-db)s 
                      --target-data-table nucleotides
{{ codestop }}

- amino acid data (%(misc-data-amino-acids)s) into a %(misc-data-amino-acids-txt)s.

{{ codestart }}
anvi-export-misc-data -c %(contigs-db)s 
                      --target-data-table amino_acids
{{ codestop }}
