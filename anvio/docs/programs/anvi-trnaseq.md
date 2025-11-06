This program **analyzes a tRNA-seq library, generating de novo predictions of tRNA sequences, structures, and modification positions**.

A FASTA file of merged paired-end tRNA-seq reads is required as input. This file is produced by the initial steps of the [trnaseq-workflow](../../workflows/trnaseq/), in which [Illumina-utils](https://github.com/merenlab/illumina-utils), merges paired-end reads and %(anvi-script-reformat-fasta)s creates anvi'o-compliant deflines in the FASTA file.

The primary output of anvi-trnaseq is a %(trnaseq-db)s. Supplemental outputs are also produced -- an analysis summary, a tabular file of unique sequences not identified as tRNA, an a tabular file of 5' and 3' extensions trimmed off mature tRNA.

The `anvi-trnaseq --help` menu provides detailed explanations of the parameters controlling the multifacted analyses performed by the program.

## Examples

*Generate a %(trnaseq-db)s from a sample using 16 cores.*

{{ codestart }}
anvi-trnaseq -f %(trnaseq-fasta)s \
             -S SAMPLE_NAME \
             -o OUTPUT_DIRECTORY \
             -T 16
{{ codestop }}

*Generate a %(trnaseq-db)s from a sample flagged as being treated with demethylase. The output directory is overwritten if it already exists.*

{{ codestart }}
anvi-trnaseq -f %(trnaseq-fasta)s \
             -S SAMPLE_NAME \
             -o OUTPUT_DIRECTORY \
             -T 16 \
             --treatment demethylase \
             --overwrite-output-destinations
{{ codestop }}

## Parameterize tRNA feature profiling

Feature profiling parameters can be modified by the user by in an optional `.ini` file. For example, the user may want a more permissive definition of a tRNA (more false positive identifications of sequences as tRNA, fewer false negative failures to identify sequences as tRNA), increasing the number of unpaired nucleotides allowed in the T stem or increasing the number of unconserved canonical nucleotides allowed in the anticodon loop. Numerous structural parameters like these can be altered.

*Write the `.ini` file to `param.ini`.*

{{ codestart }}
anvi-trnaseq --default-feature-param-file PARAM.ini
{{ codestop }}

*Nicely display the `.ini` defaults that can be written to the file in standard output.*

{{ codestart }}
anvi-trnaseq --print-default-feature-params
{{ codestop }}
