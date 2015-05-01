# User Manual

This is going to be very brief and simple, but don't hesitate to send me your questions.

--Meren

## Installation

### Dependencies

anvi'o has some dependencies, some of which will be taken care of the installer. There are two you need to make sure you have installed:

* [Chrome Web Browser](https://www.google.com/chrome/browser/desktop/). Chrome should be not only installed on your system, but unfortunately it should also be the default browser (otherwise everytime interactive interface pops up, you will need to copy-paste the address to a Chrome window). anvi'o does not support any other browser, and it will not perform optimally on others.

* [Prodigal](http://prodigal.ornl.gov/) (go to your terminal, type `prodigal` if you get an error, you need to install it). Here is a quick way to install it (the first line will not work if you don't have wget, but you can get wget installed esily typing `sudo port install wget` if you are using MacPorts system on your Mac computer):

<div style="padding-left:30px">
<pre>
wget http://prodigal.ornl.gov/zips/prodigal.v2_50.tar.gz
tar -zxvf prodigal.v2_50.tar.gz && cd prodigal.v2_50 && make
sudo cp prodigal /usr/bin/
</pre>
</div>

* [HMMer](http://hmmer.janelia.org/) (go to your terminal, type `hmmscan`, if you get an error, you need to install this one). There is an easy installer there if you follow the link. Alternatively, if you are using the port system on your Mac, you can simply type `sudo port install hmmer`.

* [SQLite](http://www.tutorialspoint.com/sqlite/sqlite_installation.htm) (go to your terminal, type sqlite3, if it gives you an error, you need to install this one). There is installation instructions on the web page if you follow the link. Or you can install it by typing `sudo port install sqlite3` if you are using the port system on your Mac.

* [MyRAST](http://blog.theseed.org/servers/) (this one is optional, but it is very useful, if you type `svr_assign_to_dna_using_figfams` and if it doesn't give you an error, you have it (press `CTRL+C` to quit)).

## Installing anvio

First get the codebase:

    git clone https://github.com/meren/anvio.git

Then go into the anvio directory, and then run the installation:

	cd anvio
	sudo python setup.py install

If you want to update your installation, you are going to need to run these commands from within the anvio directory:

    git pull
    sudo python setup.py install


## Running the mini test

"Mini test" is a minimum test set that runs almost everything in the codebase. If you have a proper installation, you shouldn't get any errors from running this test. To run it go to your anvio installation directory, and type these:

    cd tests
    ./run_mini_test.sh

All fine? Perfect!

# Running anvio like a pro

## Preparation

To run anvio you need these files:

* __A FASTA file__ of your contigs. I will call it `contigs.fa` throughout this manual.
* And __BAM files__ for your samples. Let's say, for the sake of brevity, you have two samples in your analysis, `X` and `Y`, and the BAM files for these samples are named `X-raw.bam` and `Y-raw.bam`.

Contig names in `contigs.fa` must match names found in your bam files. To make sure that is the case, do this:

    grep '>' contigs.fa | head
    anvi-profile -i X-raw.bam --list-contigs | head

Do they look identical? They better be. If they don't you need to fix that before going forward.

If you exported your FASTA file and BAM files using CLC, type these two commands, and you are going to be fine:

    sed -i 's/ .*$//g' contigs.fa
    sed -i 's/_contig_/contig/g' contigs.fa


## anvi-gen-annotation-db

Annotation database will be the essential ingredient everything you will do with anvio. It will process `contigs.fa`, and store it in a better formatted way. You can decorate your annotation database, but the one command you have to run is this:

    anvi-gen-annotation-database -f contigs.fa -o annotation.db

Each anvi-* program has a help menu that explains available parameters. Don't forget to check them for everything you run. If something is not clearly explained, let me know:

    anvi-gen-annotation-database --help

Once you have the annotation database, you can decorate it with stuff that would be very useful later.

### RAST annotation

If you have MyRAST installed, you can run these two commands to store the annotation of your contigs in the annotation database (the first line will query RAST server, which may take a while depending on the number of contigs you have), the second line will incorporate the returning info into anvio's annotation database:

    svr_assign_to_dna_using_figfams < contigs.fa > svr_assign_to_dna_using_figfams.txt 
    anvi-populate-genes-table annotation.db -p myrast_cmdline_dont_use -i svr_assign_to_dna_using_figfams.txt

I suggest you to run this, and store this file separately:

    anvi-export-genes-table annotation.db -o rast_annotation.txt

If you have to re-create the annotation database for some reason, you can use this newly generated file (`rast_annotation.txt`) to repopulate the annotation database instead of querying RAST again:

    anvi-populate-genes-table annotation.db -p default_matrix -i rast_annotation.txt

### Running HMM models

anvio can do wonders with HMM models. To decorate your annotation database with hits from standard HMM models (which, at this point, constitute published single-copy gene collections), run this command:

    anvi-populate-search-table annotation.db

### Storing CONCOCT results

If you have some CONCOCT results for your merged files (which you DON't have at this point, but you will have them later if you keep reading), you can add them into the annotation database like this:

    anvi-populate-collections-table annotation.db --parser concoct -i concoct.txt

In this case `concoct.txt` is the clustering.csv file CONCOCT generates. To see the format you can take a look at [`anvio/tests/mini_test/concoct.txt`](https://github.com/meren/anvio/blob/master/tests/mini_test/concoct.txt)

## anvi-init-bam

anvio likes BAM files when they are sorted and indexed. This is why I named your BAM file for sample X you got from CLC or Bowtie as `X-raw.bam` instead of `X.bam`. To initialize this BAM file you can run this command:

    anvi-init-bam X-raw.bam -o X

But of course it is not fun to do every BAM file you have one by one. So what to do?

A slightly better way to do would require you to do it in a `for` loop. First create a file called, say, `SAMPLE_IDs`. For your samples (`X` and `Y`) it will look like this:

    $ cat SAMPLE_IDs
    X
    Y

Then, you can run `anvi-init-bam` on all of them by typing this:

	for sample in `cat SAMPLE_IDs`; do anvi-init-bam $sample-raw.bam -o $sample; done

Good.

## anvi-profile

Profiling step makes sense of each BAM file separately by utilizing the information stored in the annotation database. The result of the profiling step is a special file that describes the run (`RUNONFO.cp`), and a profile database (`PROFILE.db`).

The minimal command to profile a BAM file looks like this:

anvi-profile -i X.bam -a annotation.db

But I encourage you to take a look at the default paramers. One of the most critical parameter is `-M` (`--min-contig-length`) parameter. The default is 10,000. Which means the profiling step will take into consideration only the contigs that are longer than 10Kb. This may be too large for your analysis. But clustering and visualization steps in anvio have some limitations, so you can't really say `-M 0` in most cases. The rule of thumb is to keep the number of contigs anvio will deal to a maximum of 20,000. How can you know how many contigs are there for a given `-M` value? Well, one thing to find that out is this:

    sqlite3 annotation.db 'select count(*) from contigs_basic_info where length > 10000;'

This command will print out the number of contigs longer than 10Kb in your dataset. You can try different values until the output is about 20,000, and use that value for `-M`. But I will not recommend you to go below 1Kb. The main reason to that is the fact that anvio relies on k-mer frequencies to better cluster contigs, and tetra-nucleotides (the default way for anvio to make sense of the sequence composotion) become very unstable very quickly.

Once you know what you `-M` is, you can, again, profile multiple samples using the similar approach we used for initializing BAM files:

    for sample in `cat SAMPLE_IDs`; do anvi-profile -i $sample.bam -M YOUR_M_VALUE -a annotation.db; done

## anvi-merge

The next step is to merge all the profiles that should be analyzed together. It goes without saying that every profiling step must have used the same parameters for analysis. If profiles have been generated with different annotation databases or with different parameters will not get merged, and you will get angry error messages from anvio.

In an ideal case, this should be enough to merge your stuff:

    anvi-merge X/RUNINFO.cp Y/RUNINFO.cp -o XY-MERGED -a annotation.db

Or alternatively you can run this:

    anvi-merge */RUNINFO.cp -o XY-MERGED -a annotation.db

When merging is done, you can export necessary files for unsupervised binning of splits using CONCOCT (to get the necessary input file for _Storing CONCOCT results_ section. For CONCOCT you need a coverages file, and splits file. You can get both running this on your merged profile:

    anvi-export-splits-and-coverages annotation.db XY-MERGED/PROFILE.db

This command will give you two files you will use for running CONCOCT (this part is going to be much better very soon).

## anvi-interactive

Once the merging is done you can finally run the interactive interface:

    anvi-interactive -r XY-MERGED/RUNINFO.cp -a annotation.db

# Other examples

## Running the visualization interface with manual input

`anvi-interactive` can also be run using text files alone, without any annotation or profile databases. In order to do that, your command line should look like this:

    anvi-interactive -f FASTA.fa -m METADATA.txt -t TREE.txt -o OUTPUT_DIRECTORY

`FASTA.fa` is a FASTA formatted file for your contigs. `METADATA.txt` is a TAB-delimited file that contains information for each entry appears in the FASTA file. Finally, the `TREE.txt` is the file where you describe the organization of your entries you have in your FASTA file. The names of the leaves of the tree, the deflines found in your FASTA file and the first column of your metadata file must match, perfectly.

Here is an example FASTA file:

<pre>
>A
ATCG
>B
GTCA
>C
GACT
>D
TACG
>E
TTCG
>F
GGCA
>G
ACAC
>H
CATC
</pre>

That can go together with this example metadata file:


| contig 	| C1   	| C2        	|
|--------	|------	|-----------	|
| A      	| 0.5  	| Cluster_1 	|
| B      	| 0.5  	| Cluster_1 	|
| C      	| 1.2  	| Cluster_0 	|
| D      	| 0.6  	| Cluster_3 	|
| E      	| 0.7  	| Cluster_3 	|
| F      	| 0.8  	| Cluster_3 	|
| G      	| 0.9  	| Cluster_4 	|
| H      	| 0.9  	| Cluster_4 	|

And finally this tree:

    ((G:0.0122349,H:0.0122349)Int11:0.450561,((C:0.0103605,D:0.0103605)Int10:0.0791402,(E:0.0262071,(F:0.00364725,(A:0,B:0)Int8:0.00364725)Int9:0.0262071)Int12:0.0791402)Int13:0.450561);

If you put these examples into cognate files, `anvi-interactive` command I showed before will work, yet it is clear this is a mock example.

Files for a bit more realistic looking example is in the [`tests/manual_interactive`](https://github.com/meren/anvio/tree/master/tests/manual_interactive) directory in the codebase. Also, `run_manual_interactive.sh` in the `tests` directory demonstrates how those files should be run.


## Generating a newick tree file from an input matrix

Please see example files stored in the [`tests/manual_interactive`](https://github.com/meren/anvio/tree/master/tests/manual_interactive) directory. If you have a matrix file (such as [`data.txt`](https://github.com/meren/anvio/blob/master/tests/manual_interactive/data.txt)), you can get a newick file for it by running this anvi'o command:

    anvi-matrix-to-newick data.txt -o tree.txt


## Exporting tables from databases

You can export the information stored in any anvio database using `anvi-get-db-table-as-matrix`.

You can see what tables are available in a given anvio dataase,

<pre>
$ anvi-get-db-table-as-matrix ANNOTATION.db --list
self
contig_sequences
kmer_contigs
kmer_splits
contigs_basic_info
splits_basic_info
hmm_hits_info
hmm_hits_in_splits
hmm_hits_in_contigs
genes_in_contigs
genes_in_splits_summary
genes_in_splits
collections_info
collections_colors
collections_of_contigs
collections_of_splits
</pre>

And export one of them:

<pre>
$ anvi-get-db-table-as-matrix ANNOTATION.db -t hmm_hits_in_contigs
Database .....................................: "ANNOTATION.db" has been initiated with its 16 tables.
Table ........................................: "hmm_hits_in_contigs" has been read with 2 entries and 8 columns.
Output .......................................: hmm_hits_in_contigs.txt
</pre>

You can also export only certain fileds from a table. To see what fields are available, use `--list` with a table:

<pre>
$ anvi-get-db-table-as-matrix ANNOTATION.db -t hmm_hits_in_contigs --list
Database .....................................: "ANNOTATION.db" has been initiated with its 16 tables.
Table ........................................: "hmm_hits_in_contigs" has been read with 2 entries and 8 columns.
Table columns ................................: "entry_id, source, contig, start, stop, gene_name, gene_id, e_value"
</pre>

Then you can specify which columns you would like to be exported from the table of interest like this:

<pre>
$ anvi-get-db-table-as-matrix ANNOTATION.db -t hmm_hits_in_contigs --fields "contig, gene_name"
Database .....................................: "ANNOTATION.db" has been initiated with its 16 tables.
Table ........................................: "hmm_hits_in_contigs" has been read with 2 entries and 8 columns.
Columns to report ............................: "contig, gene_name"
Output .......................................: hmm_hits_in_contigs.txt
</pre>

This time `hmm_hits_in_contigs.txt` will contain only two of the columns you specified. 


---

# FAQ

### I run into an issue, who's fault is it?

It is Tom's. But you can always enter a [bug report](https://github.com/meren/anvio/issues), and it would be very useful for us to remember to fix that. The codebase is young, and there will be issues.

### I have suggestions!

This is great. That's exactly why you have accesss to this codebase. Either file an [issue](https://github.com/meren/anvio/issues), or send an e-mail (a.murat.eren@gmail.com).

### Do I have permission to change the codebase?

Yes. Don't forget run `run_mini_test.sh` to make sure everything is running, and don't be worried. I will make sure you didn't break anything :)
