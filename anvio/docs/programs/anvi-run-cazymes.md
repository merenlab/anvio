This program **annotates genes in your %(contigs-db)s with functions using dbCAN [CAZyme HMMs](https://bcb.unl.edu/dbCAN2/download/Databases/)** 

Before you run this program, you'll have to set up the CAZyme database on your computer with the program %(anvi-setup-cazymes)s.  

The CAZyme database is based on protein sequences, so anvi'o will convert your genetic information into protein sequences and then use HMMs to compare them to the database. 

{:.notice}
Unsure what an HMM is? Check out [our vocab page](http://merenlab.org/vocabulary/#hmm)

To run, you'll need to provide a %(contigs-db)s and the output will be a %(functions)s artifact. Here is a default run: 

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s 
{{ codestop }}

If you stored the %(cazyme-data)s that you got from running %(anvi-setup-cazymes)s in a custom location, you'll need to provide that path as well.

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s \
                 --cazyme-data-dir %(cazyme-data)s 
{{ codestop }}

By default, this uses `hmmsearch` to run HMMs. You can choose to use `hmmscan` instead by running

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s \
                 --cazyme-data-dir %(cazyme-data)s \
                 --hmmer-program hmmscan
{{ codestop }}

Use the parameter `--noise-cutoff-terms` to filter out hits. The default value is `--noise-cutoff-terms -E 1e-12`. If you want to explore filtering options, check out the help menu of the underlying hmm program you are using e.g. `hmmsearch -h`

{{ codestart }}
anvi-run-cazymes -c %(contigs-db)s \
                 --noise-cutoff-terms "-E 1e-14"
{{ codestop }}

## Exploring CAZyme annotations

The dbCAN HMM files provide only limited descriptions of CAZyme protein families and often require additional investigation following annotation. The [dbCAN3](https://bcb.unl.edu/dbCAN2/) database offers several [supplementary files](https://bcb.unl.edu/dbCAN2/download/Databases/) that can support deeper exploration, including:

- **`dbCAN_sub.hmm`** — HMM to classify proteins based on predicted substrate specificity.  
- **`FamInfo.txt.08022020.xls`** — Reference table summarizing functional details for each CAZyme family.

## Import CAZyme functions from run_dbcan

If anvi'o users are interested in importing more detailed annotations, e.g. substrate level prediction from dbCAN-sub, they should consider using [run_dbcan](https://dbcan.readthedocs.io/en/latest/) from [dbCAN3](https://bcb.unl.edu/dbCAN2/). Here are some quick steps:

Step 1. Extract amino acid sequences from your %(contigs-db)s

{{ codestart }}
anvi-get-sequences-for-gene-calls -c %(contigs-db)s --get-aa-sequences -o genes.faa
{{ codestop }}

Step 2. Run dbCAN3 via [run_dbcan](https://dbcan.readthedocs.io/en/latest/) to annotate those amino acid sequences

Step 3. Create a functions-txt with CAZyme functions

The program [run_dbcan](https://dbcan.readthedocs.io/en/latest/) has multiple [output files](https://dbcan.readthedocs.io/en/latest/user_guide/quick_start.html#understanding-the-output) which can be parsed into a %(functions-txt)s, for example, the [overview.txt](https://dbcan.readthedocs.io/en/latest/user_guide/quick_start.html#understanding-the-output). 

Step 4. Import your new CAZyme %(functions-txt)s back into your %(contigs-db)s

{{ codestart }}
anvi-import-functions -c %(contigs-db)s -i cazyme-functions.txt
{{ codestop }}