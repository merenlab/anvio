This program lets you look at the %(hmm-hits)s from a single %(hmm-source)s across multiple genomes or bins, by creating a %(hmm-hits-across-genomes-txt)s.

{:.notice}
For a simlar program that reports function hits across genomes, see %(anvi-script-gen-function-matrix-across-genomes)s.

The input of this program can be either an %(internal-genomes)s or an %(external-genomes)s.

Here are two example run on an internal-genomes:

{{ codestart }}
anvi-script-gen-hmm-hits-matrix-across-genomes -i %(internal-genomes)s \
                                               --hmm-source Bacteria_71 \
                                               -o output.txt
{{ codestop }}

To list the %(hmm-source)ss common to the datasets that you're analyzing, just add the flag `--list-hmm-sources`, as so:

{{ codestart }}
anvi-script-gen-hmm-hits-matrix-across-genomes -e %(external-genomes)s \
                                               --list-hmm-sources
{{ codestop }}
