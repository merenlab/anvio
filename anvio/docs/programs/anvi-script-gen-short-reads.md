This program uses already assembled contigs to create a mock list of short reads. You can then use these short reads to reassemble your data in order to test alternative assembly programs or analysis methods as a positive control. 

Basically, this attempts to undo the assembly and produce a data set that could have been directly received from laboratory sequencing. While the computer's mock short reads won't be perfect, they can be used to make sure your analysis pipeline is working from step 1. 

## Example Usage

This program takes an INI file - a form of text file containing various information. For this program, the example provided in the anvi'o test suite looks like this: 

```ini
[general]
short_read_length = 10
error_rate = 0.05
coverage = 100
contig = CTGTGGTTACGCCACCTTGAGAGATATTAGTCGCGTATTGCATCCGTGCCGACAAATTGCCCAACGCATCGTTCCTTCTCCTAAGTAATTTAACATGCGT
```

Note that this file contains both the contig that you want to break down, and various information about the short reads that you want to create. To run this program, just call 

{{ codestart }}
anvi-script-gen-short-reads %(configuration-ini)s \
                            --output-file-path %(short-reads-fasta)s
{{ codestop }}
    
The resulting FASTA file with short reads will cover the `contig` with short reads that are 10 nts long at 100X coverage. There will also be an error-rate of 0.05, to mimic the sequencing errors you would get from sequencing in the wet lab. 
