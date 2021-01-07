The input for this program is a **properly formatted** %(trnaseq-fasta)s, containing sequences from a tRNA-seq sample or split.

The program identifies tRNA among the input sequences, profiles the tRNA primary sequence and secondary structure, and filters single nucleotide variants from modified nucleotides.

The primary output of the program is a %(trnaseq-db)s. Supplemental files are also produced: an analysis summary file, a tab-separated file of unique sequences not identified as tRNA, and a tab-separated file showing the range of 5' and 3' variants trimmed from tRNA sequences.

We encourage you to read the list of options in the `anvi-trnaseq --help` menu to understand how the user can manipulate the multifaceted analyses performed by the program.

The program can generate a .ini file for tRNA feature parameterization using an alternate command, `anvi-trnaseq --default-feature-param-file <param.ini>`. The default parameterizations in the file can be modified by the user, and the file can be used as the `--feature-param-file` argument in the main mode of the program. `anvi-trnaseq --print-default-feature-params` can also be used to quickly and neatly display the defaults in the terminal.


### Create a tRNA-seq database from a FASTA file, using 16 cores

{{ codestart }}
anvi-trnaseq -f %(trnaseq-fasta)s \
             -S example_sample_name \
             -o example_empty_output_directory_path \
             -T 16
{{ codestop }}

### Create a tRNA-seq database from a sample identified as a demethylase split, overwriting the output directory if it already exists

{{ codestart }}
anvi-trnaseq -f %(trnaseq-fasta)s \
             -S example_sample_name \
             -o example_empty_output_directory_path \
             -T 16 \
             --treatment demethylase \
             -W
{{ codestop }}