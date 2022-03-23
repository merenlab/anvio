Generates TAB-delmited output files for %(functions)s from a single functional source across genomes.

The input genomes for this program can be provided through an %(external-genomes)s, %(internal-genomes)s, %(genomes-storage-db)s, or any combination of these sources.

This program is very similar to %(anvi-display-functions)s, and can also perform a functional enrichment analysis on-the-fly if you provide it with an optional %(groups-txt)s file. Unlike, %(anvi-display-functions)s, this program will report TAB-delmited output files for you to further analyze.

You can run the program on a set of genomes for a given annotation source:

{{ codestart }}
anvi-script-gen-function-matrix-across-genomes -e %(external-genomes)s \
                                               --annotation-source COG20_FUNCTION \
                                               --output-file-prefix MY-GENOMES
{{ codestop }}

The command above will result in two files in your work directory, both of which will be of type %(functions-across-genomes-txt)s:

* MY-GENOMES-FREQUENCY.txt
* MY-GENOMES-PRESENCE-ABSENCE.txt

{:.notice}
You can always learn about which functions are in a given %(contigs-db)s using the program %(anvi-db-info)s.

Alternatively you can run it with a %(groups-txt)s that associates sets of genomes with distinct groups,

{{ codestart }}
anvi-script-gen-function-matrix-across-genomes -i %(internal-genomes)s \
                                               --annotation-source COG20_FUNCTION \
                                               --output-file-prefix MY-GENOMES \
                                               --groups-txt groups.txt
{{ codestop }}

which would generate an additional file in your work directory of type %(functional-enrichment-txt)s:

* MY-GENOMES-FUNCTIONAL-ENRICHMENT.txt
