This program **associates genes in your %(contigs-db)s with functions using the EBI's [Pfam database](https://pfam.xfam.org/).** 

Before you run this program, you'll have to set up the Pfam database on your computer with the program %(anvi-setup-pfams)s.  

The Pfam database is based on protein sequences, so anvi'o will convert your genetic information into protein sequences and then use HMMs to compare them to the database. 

{:.notice}
Unsure what an HMM is? Check out [our vocab page](http://merenlab.org/vocabulary/#hmm)

To run, you'll need to provide a %(contigs-db)s. If you stored the %(pfams-data)s that you got from running %(anvi-setup-pfams)s in a custom location, you'll need to provide that path as well. The output is a %(functions)s artifact. 

Here is a default run: 

{{ codestart }}
anvi-run-pfams -c %(contigs-db)s \
            --pfam-data-dir %(pfams-data)s 
{{ codestop }}

By default, this uses `hmmsearch` to run HMMs. You can choose to use `hmmscan` instead by running

{{ codestart }}
anvi-run-pfams -c %(contigs-db)s \
            --pfam-data-dir %(pfams-data)s \
            --hmmer-program hmmscan
{{ codestop }}

See [this article](https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/) for a discussion on the performance of the two HMMER programs. 
