This program gets the coverage values from a %(bam-file)s, and puts them into a %(coverages-txt)s. 

You must provide a BAM file, but there are three ways you can choose contigs to analyze within that file: 
1. Give a contig name. Here, you can only report coverage per nucleotide position (In this example, the user is specifically asking for this anyway with the `-m` flag)

    {{ codestart }}
    anvi-script-get-coverage-from-bam -b %(bam-file)s \ 
                                     -c NAME_OF_CONTIG \ 
                                     -m pos
    {{ codestop }}

2. Give a file that contains a list of contigs (one per line; same format as the `--contigs-of-interest` tag for %(anvi-profile)s). Here, you can ask for the contig averages (as in this example) or nucleotide position coverage. 

    {{ codestart }}
    anvi-script-get-coverage-from-bam -b %(bam-file)s \ 
                                     -l NAME_OF_FILE \
                                     -m contig
    {{ codestop }}

3. Give a %(collection-txt)s file for the program to determine the coverage for all contigs in those bins. Here, you can ask for the contig averages, nucleotide position coverage or coverage per bin (as in this example). 

    {{ codestart }}
    anvi-script-get-coverage-from-bam -b %(bam-file)s \ 
                                     -C %(collection-txt)s \
                                     -m bin
    {{ codestop }}

