Stores %(hmm-hits)s for a given %(hmm-source)s in a %(contigs-db)s. In short, this is the program that will do a search for HMMs against a %(contigs-db)s and store that information into the contigs-db's %(hmm-hits)s.

This is one of the programs that users commonly run on newly generated %(contigs-db)s, along with %(anvi-scan-trnas)s, %(anvi-run-ncbi-cogs)s, %(anvi-run-scg-taxonomy)s, and so on.

### HMMs in the context of anvi'o

In a nutshell, [hidden Markov models](https://en.wikipedia.org/wiki/Hidden_Markov_model) are statistical models typically generated from known genes which enable 'searching' for similar genes in other sequence contexts.

The default anvi'o distribution includes numerous [curated HMM profiles](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm) for single-copy core genes and ribosomal RNAs, and anvi'o can work with custom HMM profiles provided by the user. In anvi'o lingo, each of these HMM profiles, whether they are built-in or user defined, is called an %(hmm-source)s.

### Default Usage

To run this program with all default settings (against all default anvi'o %(hmm-source)s), you only need to provide a %(contigs-db)s:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s
{{ codestop }}

Multithreading will dramatically improve the performance of `anvi-run-hmms`. If you have multiple CPUs or cores, you may parallelize your search:


{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              --num-threads 6
{{ codestop }}


You can also run this program on a specific built-in %(hmm-source)s:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71
{{ codestop }}

### User-defined HMMs

Running `anvi-run-hmms` with a custom model is easy. All you need to do is to create a directory with necessary files:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -H MY_HMM_PROFILE
{{ codestop }}

See the relevant section in the artifact %(hmm-source)s for details.


### Changing the HMMER program

By default, `anvi-run-hmms` will use [HMMER](http://hmmer.org/)'s `hmmscan` for amino acid HMM profiles, but you can use `hmmsearch` if you are searching a very large number of models against a relatively smaller number of sequences:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              --hmmer-program hmmsearch
{{ codestop }}

{:.notice}
This flag has no effect when your HMM profile source is for nucleotide sequences (like any of the Ribosomal RNA sources). In those cases anvi'o will use `nhmmscan` exclusively.

### Saving the HMMER output

If you want to see the output from the HMMER program (eg, `hmmscan`) used to annotate your data, you can request that it be saved in a directory of your choosing. Please note that this only works when you are running on a single HMM source, as in the example below:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71 \
              --hmmer-output-dir OUTPUT_DIR
{{ codestop }}

If you do this, file(s) with the prefix `hmm` will appear in that directory, with the file extension indicating the format of the output file. For example, the table output format would be called `hmm.table`.

{:.warning}
These resulting files are not _exactly_ the raw output of HMMER because anvi'o does quite a bit of pre-processing on the raw input and output file(s) while jumping through some hoops to make the HMM searches multi-threaded. If this is causing you a lot of headache, please let us know.

#### Requesting domain table output

No matter what, anvi'o will use the regular table output to annotate your contigs database. However, if you are using the --hmmer-output-dir to store the HMMER output, you can also request a domain table output using the flag `--get-domtable-output`.

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              -I Bacteria_71 \
              --hmmer-output-dir OUTPUT_DIR \
              --get-domtable-output`
{{ codestop }}

In this case anvi'o will run [HMMER](http://hmmer.org) using the `--domtblout` flag to generate this output file.

{:.notice}
This flag will only work with HMM profiles made for amino acid sequences. Profiles for nucleotide sequences require the use of the program `nhmmscan`, which does not have an option to store domain output.

Please note that this output **won't be used to filter hits to be added to the contigs database**. But it will give you the necessary output file to investigate the coverage of HMM hits. But you can use the program %(anvi-script-filter-hmm-hits-table)s with this file to remove weak hits from your HMM hits table later.


### Other things anvi-run-hmms can do

* Add the tag `--also-scan-trnas` to basically run %(anvi-scan-trnas)s for you at the same time. It's very convenient. (But it only works if you are not using the `-I` or `-H` flags at the same time because reasons.)

