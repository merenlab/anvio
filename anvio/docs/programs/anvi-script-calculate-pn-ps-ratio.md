This program **calculates the pN/pS ratio** for each gene in a %(contigs-db)s and outputs it as a %(pn-ps-data)s artifact. 

### What is the pN/pS ratio? 

The pN/pS ratio (described in [Hannigan et al. 2017](https://doi.org/10.7717/peerj.2959)) is the metagenomic equivalent of the dN/dS ratio (see [Schloissnig et al. 2012](https://doi.org/10.1038/nature11711)), both of which are measures of metagenomic variation. In short, usually the pN/pS ratio describes the ratio of non-symonmous SNPs in a gene to synonymous SNPs, normalized by the likelihood of mutations happening by chance. Anvi'o examines SCV (single codon variants) instead in order to get a more accurate representation of synonymity.  

### Neat. How do I use this program? 

Firstly, you'll need to run %(anvi-gen-variability-profile)s twice with the same parameters on the same databases. The first time, use the flag `--engine AA` to get a %(variability-profile-txt)s for SAAVs (single amino acid variants), which we'll name the `SAAVs.txt` in this example. The second time, use the flag `--engine CDN` to get a %(variability-profile-txt)s for SCVs (single codon variants), which we'll name `SCVs.txt` in this example. 

Then you can run this program like so:

{{ codestart }}
anvi-script-calculate-pn-ps-ratio -a SAAVs.txt \
                                  -b SCVs.txt \ 
                                  -c %(contigs-db)s \
                                  -o output_dir 
{{ codestop }}

This will result in a directory called `output_dir` that contains several tables that describe each of your genes. See %(pn-ps-data)s for more information. 

### Other parameters

By default, this program ignores some of the genes and variable positions in your variability profiles; you can choose to be more sensitive or ignore more positions by changing any of these three variables:  

- The minimum departure from consensus for a variable position (default: 0.10). 
- The minimum number of variants in a gene (default: 4). 
- The minimum coverage at a variable position (default: 30)
