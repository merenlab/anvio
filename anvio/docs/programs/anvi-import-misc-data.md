This program lets you **bring additional information into your anvi'o databases** that will appear when you run %(anvi-interactive)s.   

With this, you can 
- **bring additional data about your items or layers into the anvi'o interactive interface** by putting it into a %(pan-db)s or %(profile-db)s
- **bring additional data about your nucleotides/amino acids** into a %(contigs-db)s 

## Items, Layers, and the Interactive Interface 

{:.notice}
This process, as well as the definition of an item and a layer, are described in more detail in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables). 

Basically, you can add additional information to the interactive interface by running this program on the database you want to display and a text file containing your information. You can do this with three types of data: 

1. %(misc-data-items)s by providing a %(misc-data-items-txt)s. This contains information about *each of your items in the central tree* (whether those are contigs, bins, or genes), and will appear as **additional concentric circles** when you run %(anvi-interactive)s. 

    {{ codestart }}
    anvi-import-misc-data -p %(profile-db)s \
                          -t items \
                          %(misc-data-items-txt)s 
    {{ codestop }}
        
2. %(misc-data-layers)s by providing a %(misc-data-layers-txt)s. This contains information about *each layer (or concentric circle) of the interface* (which usually correspond to your samples), and will appear as **graphs in line with your circles of data** (on the right, similar to how to the titles of each layer are displayed at the top) when you run %(anvi-interactive)s. 

    {{ codestart }}
    anvi-import-misc-data -p %(pan-db)s \
                          -t layers \
                          %(misc-data-layers-txt)s                               
    {{ codestop }}

3. %(misc-data-layer-orders)s by providing a %(misc-data-layer-orders-txt)s. This contains information about *what order you want the concentric circles to be displayed in*  (which usually correspond to your samples), and will appear as **above the misc-data-layers graphs as a tree** when you run %(anvi-interactive)s. 

    {{ codestart }}
    anvi-import-misc-data -p %(profile-db)s \
                          -t layer_orders \
                          %(misc-data-layer-orders-txt)s 
    {{ codestop }}
    
### Additional notes 

You also have the option to associate keys with only a specific data group, or transpose the input before processing (in case you misformatted it). 

If you no longer want to see data you've added with this function, you can export it as the original text file with %(anvi-export-misc-data)s and delete it from the database with %(anvi-delete-misc-data)s.

## Nucleotides, Amino Acids, and Contigs Databases


