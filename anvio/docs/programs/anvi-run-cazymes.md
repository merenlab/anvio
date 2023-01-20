This program **annotates genes in your %(contigs-db)s with functions using dbCAN [CAZyme HMMs](https://bcb.unl.edu/dbCAN2/download/Databases/)** 

Before you run this program, you'll have to set up the CAZyme database on your computer with the program %(anvi-setup-cazymes)s.  

The CAZyme database is based on protein sequences, so anvi'o will convert your genetic information into protein sequences and then use HMMs to compare them to the database. 

{:.notice}
Unsure what an HMM is? Check out [our vocab page](http://merenlab.org/vocabulary/#hmm)

To run, you'll need to provide a %(contigs-db)s. If you stored the %(cazyme-data)s that you got from running %(anvi-setup-cazymes)s in a custom location, you'll need to provide that path as well. The output is a %(functions)s artifact. 

Here is a default run: 

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s \
            --cazyme-data-dir %(cazyme-data)s 
{{ codestop }}

By default, this uses `hmmsearch` to run HMMs. You can choose to use `hmmscan` instead by running

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s \
            --cazyme-data-dir %(cazyme-data)s \
            --noise-cutoff-terms "-E 1" \
            --hmmer-program hmmscan
{{ codestop }}

Use the parameter `--noise-cutoff-terms` to filter out hits e.g. `--noise-cutoff-terms "-E 1"`. If you want to explore filtering options, check out the help menu of the underlying hmm program you are using e.g. `hmmsearch -h`

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s \
            --cazyme-data-dir %(cazyme-data)s \
            --noise-cutoff-terms "-E 1"
{{ codestop }}