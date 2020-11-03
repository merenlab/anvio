This program **lists the additional data** that is stored within a %(pan-db)s, %(profile-db)s or %(contigs-db)s. This is data that can be imported with %(anvi-import-misc-data)s and is displayed in the interactive interface. 

When run, this program will output to the terminal a list of all additional data tables that are stored within the database. If you want to export a specific element of these as a text file, see %(anvi-export-misc-data)s. 

### What is displayed? 

When running on a %(profile-db)s or %(pan-db)s, the output will display the following types of data:

- %(misc-data-items)s 
- %(misc-data-layers)s
- %(misc-data-layer-orders)s (by default, this will include orders like `abundance` and `mean_coverage (newick)`)

When running on a %(contigs-db)s, the output will display the following types of data:

- %(misc-data-nucleotides)s 
- %(misc-data-amino-acids)s 

These have no default values and will only contain data that has been imported with %(anvi-import-misc-data)s. 

You also have the option to specify a specific kind of additional data table with `-t`. For example, to view only %(misc-data-items)s in a %(profile-db)s, just call

{{ codestart }}
anvi-show-misc-data -p %(profile-db)s \
                    -t items 
{{ codestop }}

Similarly to importing and exporting additional data tables, you can also focus on a specific data group with the parameter `-D`.
