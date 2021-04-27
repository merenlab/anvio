This program computes the average nucleotide identity between reads in a single fasta file (using PyANI). 

To compute the ANI (or other genome distance metrics) between two genomes in different fasta files, use %(anvi-compute-genome-similarity)s. 

A default run of this program looks like this: 

{{ codestart }}
anvi-script-compute-ani-for-fasta -f %(fasta)s \ 
                                  -o path/to/output \
                                  --method ANIb
{{ codestop }}

By default, the PyANI method is ANIb (which aligns 1020 nt fragments of your sequences using BLASTN+). You can switch to ANIm, ANIblastall, or TETRA if desired. See the [PyANI documentation](https://github.com/widdowquinn/pyani) for more informaiton. 

You also have the option to change the distance metric (from the default "euclidean") or the linkage method (from the default "ward") or provide a path to a log file for debug messages. 
