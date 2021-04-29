This program identifies the tRNA genes in a %(contigs-db)s and stores them in an %(hmm-hits)s. 

To run, just provide a %(contigs-db)s that you want to look through. 

{{ codestart }}
anvi-scan-trnas -c %(contigs-db)s
{{ codestop }}

### Customizing the cut off score

What counts as a tRNA gene? That could be up to you. 

The default minimum score for a gene to be counted is 20.  However, you can set this cutoff to anywhere between 0-100. This value is actually used by the module tRNAScan-SE, so view [their documentation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768409/) for details. For example, to find more non-cononical tRNA genes, a user could lower the cutoff score to 10 as follows:

{{ codestart }}
anvi-scan-trnas -c %(contigs-db)s \
                --trna-cutoff-score 10
{{ codestop }}

### Other options 

- It is easy to modify where the outputs will go:

    - Use the parameter `--log-file` to provide a path for the output messages to go.
    
    - Use the parameter `--trna-hits-file` to provide a path for the raw tRNA scan data to go. 
    
- Like many anvi'o programs, you can use the tag `--just-do-it` to not have to look at questions or warnings

- You can also try to multithread whenever possible by setting the `--num-threads` parameter (it is 1 by default). This can be used to speed up runtime, but please be aware of your system and its limitations before trying this. 

### Understanding the output 

Essentially, the output of this program states the probability that each gene is a tRNA gene. See %(hmm-hits)s for more information. 
