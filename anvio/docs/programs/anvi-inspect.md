This lets you inspect a single split across your samples. This interface can also be opened from the %(anvi-interactive)s interface by asking for details about a specific split.

From this view, you can clearly see the coverage and detection across your split, all SNVs, and the genes identified within your split and their functional annotations. You can also  easily compare all of this data across all of the samples that this split is present in.  

To run this program, just provide a %(profile-db)s and %(contigs-db)s pair and a single split name to inspect. 

{{ codestart }}
anvi-inspect -p %(profile-db)s \
             -c %(contigs-db)s \ 
             --split-name Day17a_QCcontig9_split_00003
{{ codestop }}

You can also choose to hide SNVs marked as outliers or configure the server in various ways. 
