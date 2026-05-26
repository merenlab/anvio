This program scans a directory for FASTA files and generates a %(fasta-txt)s file that lists each one with its name and path.

A default run looks like this:

{{ codestart }}
anvi-script-gen-fasta-txt --input-dir path/to/fasta/dir \
                          -o %(fasta-txt)s
{{ codestop }}

The program identifies all FASTA-formatted files in the given directory, uses the filename (without extension) as the sequence name, and writes a TAB-delimited %(fasta-txt)s file ready for use with other anvi'o programs.
