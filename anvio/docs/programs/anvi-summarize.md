Anvi-summarize lets you look at a **comprehensive overview of your %(collection)s** and its many statistics that anvi'o has calculated. 

It will create a folder called `SUMMARY` that contains many different summary files, including an HTML output that conviently displays them all for you. This folder will contain anything a future user might use to import your collection, so it's useful to send to others or transfer an entire anvi'o collection and all of its data. 

In a little more detail, this program will   
* generate %(fasta)s files containing your original contigs.   
* estimate various stats about each of your bins, including competition, redundacy, and information about all of your %(hmm-hits)s    
* generate various tab-delimited matrix files with information about your bins across your samples, including various statistics.   

## Running anvi-summarize 

### Running on a profile database

A standard run of anvi-summarize on a %(profile-db)s will look something like this:

{{ codestart }}
anvi-summarize -c %(contigs-db)s \
               -p %(profile-db)s \
               -o MY_SUMMARY \
               -C %(collection)s
{{ codestop }}

This will name the output directory `MY_SUMMARY` instead of the standard `SUMMARY`. 

When running on a profile database, you also have options to 
* output very accurate (but intensely processed) coverage and detection data for each gene (using `--init-gene-coverages`)
* edit your contig names so that they contain the name of the bin that the contig is in (using `--reformat-contig-names`)
* also display the amino acid sequeunces for your gene calls.  (using `--report-aa-seqs-for-gene-calls`)

### Running on a pan database

When running on a %(pan-db)s, you'll want to instead provide the associated genomes storage database. 

{{ codestart }}
anvi-summarize -g %(genomes-storage-db)s \
               -p %(pan-db)s \
               -C %(collection)s 
{{ codestop }}

You can also choose to display DNA sequences for your gene clusters instead of amino acid sequences with the flag `--report-DNA-sequences`

### Other notes

If you're unsure what collections are in your database, you can run this program with the flag `--list-collections` or by running %(anvi-show-collections-and-bins)s.

You can also use the flag `--quick-summary` to get a less comprehensive summary with a much shorter processing time. 
