Stores %(hmm-hits)s for a given %(hmm-source)s in a %(contigs-db)s. In short, this is the program that will do a search for HMMs against a %(contigs-db)s and store that information into the contigs-db's %(hmm-hits)s.

This is one of the programs that users commonly run on newly generated %(contigs-db)s, along with %(anvi-scan-trnas)s, %(anvi-run-ncbi-cogs)s, %(anvi-run-scg-taxonomy)s, and so on.

### What is an HMM?

Check out the lovely vocabulary page for an example [here](http://merenlab.org/vocabulary/#hmm).

Essentially, this program will help annotate the genes in your %(contigs-db)s, using either one of the databases built into anvi'o or a custom database.

Basically, in anvi'o, Hidden Markov Models (or HMMs for short) are used to search for specific genes with known functions in a larger dataset. Nucleotide patterns for specific gene functions are contained in an %(hmm-source)s and this program uses them to search through the data in your %(contigs-db)s.

### Default Usage

To run this program with all default settings (against all default anvio %(hmm-source)s), you only need to provide a %(contigs-db)s.

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB
{{ codestop }}

### Running against a custom set of hmm-source

In order to run against your own %(hmm-source)s or a custom subset of anvi'o's hmm-sources, you have two choices.

#### Choice 1: I have my own hmm-sources on my computer

This way the source can be completely outside of anvi'o.

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB -H path_to_your_hmm_profile
{{ codestop }}

#### Choice 2: I prefer anvi'o's hmm-sources, but I don't need all of them.

By default, anvi'o will look through all of its hmm-sources when doing a search. If you only want to run against a specific one, you're in the right place. These are the currently available ones: "Bacteria_71" (type: singlecopy), "Archaea_76" (type: singlecopy), "Protista_83" (type: singlecopy), and "Ribosomal_RNAs" (type: Ribosomal_RNAs). See the page for %(hmm-source)s for more information.

For example,

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB -I Bacteria_71
{{ codestop }}

### Changing the HMMER program
Internally, this script calls an HMMER program to do its searching of sequences in your database against HMM profiles. By default, this program is `hmmscan` for amino acid HMM profiles, but you can change that to be `hmmsearch` if you want.

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB --hmmer-program hmmsearch
{{ codestop }}

This flag has no effect when your HMM profile source is for nucleotide sequences (like any of the Ribosomal RNA sources). In those cases, the program `nhmmscan` is always used.

### Saving the HMMER output

If you want to see the output from the HMMER program (eg, `hmmscan`) used to annotate your data, you can request that it be saved in a directory of your choosing. Please note that this only works when you are running on a single HMM source, as in the example below:

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB -I Bacteria_71 --hmmer-output-dir HMMER_OUTPUT_DIRNAME
{{ codestop }}

If you do this, file(s) with the prefix `hmm` will appear in that directory, with the file extension indicating the format of the output file. For example, the table output format would be called `hmm.table`.

Warning - these files are not _exactly_ the raw output of HMMER (though arguably they should be) because anvi'o does a bit of pre-processing on the raw output file(s) while jumping through some hoops to make the HMM searches multi-threaded; we hope we will get around to fixing this soon.

#### Requesting domain table output
No matter what, anvi'o will use the regular table output to annotate your contigs database. However, if you are using the --hmmer-output-dir to store the HMMER output and you also want to see the domain table output, you can request it with the `--get-domtable-output` flag.

{{ codestart }}
anvi-run-hmms -c CONTIGS_DB -I Bacteria_71 --hmmer-output-dir HMMER_OUTPUT_DIRNAME --get-domtable-output`
{{ codestop }}

Then anvi'o will run HMMER using the `--domtblout flag` to get this output for you (though again, this output won't be used to add hits to the contigs database). You should see the file `hmm.domtable` appear in the requested output directory.

This will only work with HMM profiles made for amino acid sequences, such as Bacteria_71. Profiles for nucleotide sequences, like any of the Ribosomal RNA sources (eg, Ribosomal_RNA_23S), require the use of the `nhmmscan` program, which does not have a `--domtblout flag` option.


### Other things anvi-run-hmms can do

- Add the tag `--also-scan-trnas` to basically run %(anvi-scan-trnas)s for you at the same time. It's very convenient.
- Add the tag `--just-do-it` to hide all warnings and questions in case you don't want to deal with those.
-  There are also parameters that can help speed up the runtime of this program. However, be aware of the limits of your system, especially if running on a SGE.  For example, you can increase the number of threads or switch to hmmsearch if you are scanning  a large umber of HMMs. For more information on that, check out [here](http://merenlab.org/software/anvio/vignette/#anvi-run-hmms).

### See anvi-run-hmms in action

On the [metagenomic workflow tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-hmms)!
