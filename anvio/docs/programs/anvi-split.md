This program **creates smaller, self-contained anvi'o projects for each of the %(bin)ss in your project.** This is useful if you would like to share a subset of an anvi'o project. 

Simply provide either a %(contigs-db)s and %(profile-db)s pair or a %(genomes-storage-db)s and %(pan-db)s pair, as well as a %(collection)s, and it will create directories for each of your bins that contain their own databases and information. In other words, each of these directories will contain their own anvi'o projects that represent the contigs or genomes stored in that single bin. 

### An example run 

For example, let's say a %(profile-db)s has a %(collection)s with three bins, which are (very creatively) called `BIN_1`, `BIN_2`, and `BIN_3`.  

If you ran the following code: 

{{ codestart }}
anvi-split -p %(profile-db)s \
           -c %(contigs-db)s \
           -C %(collection)s \
           -o MY_PATH
{{ codestop }}

Then in the location `MY_PATH`, you would have three folders: `BIN_1`, `BIN_2`, and `BIN_3`.  Each one contains its own %(profile-db)s and %(contigs-db)s that only contains the contigs from that bin. You can then give a fellow anvi'o user just the `BIN_1` directory and they can get to work. 

Similarly, if you provide a %(genomes-storage-db)s and %(pan-db)s pair, the directories will contain their own smaller %(genomes-storage-db)s and %(pan-db)s pairs. 

### Other options 

You are also able to skip generating variability tables or compress the auxiliary data to save space. 
