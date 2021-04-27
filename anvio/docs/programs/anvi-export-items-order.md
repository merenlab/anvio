This program, as one might think, allows you to export a %(misc-data-items-order)s, outputing a %(misc-data-items-order-txt)s. 

You can export one of the item orders in a %(profile-db)s or %(pan-db)s as follows: 

{{ codestart }}
anvi-export-items-order -p %(profile-db)s \
                        --name cov
{{ codestop }}

The `cov` here refers to the tree that is generated using only differential coverage. Almost all anvi'o profile databases will also have available an items-order based on the tetranucleotide frequency called `tnf`, and one based on both called `tnf-cov`. 

However, to list the item orders available in this database, just don't include the name flag.  

{{ codestart }}
anvi-export-items-order -p %(pan-db)s 
{{ codestop }}

You'll get a `Config Error` that will tell you what item orders are available. 
