This program runs [InteracDome](https://interacdome.princeton.edu/) on your %(contigs-db)s, which **finds the per-residue binding scores for all of your genes**. 

The full process that this program goes through is detailed in [this blog post by Evan Kiefl](https://merenlab.org/2020/07/22/interacdome/). In summary, this program runs the HMM search against all of the genes in your %(contigs-db)s, parses and filters the results, and then stores the per-residue binding frequencies for each gene into the %(contigs-db)s.  

Before running this program, you'll have to have run %(anvi-setup-interacdome)s to set up a local copy of [InteracDome's tab-separated files](https://interacdome.princeton.edu/#tab-6136-4).

### Parameters

A basic run of this program looks like this: 

{{ codestart }}
anvi-run-interacdome -c %(contigs-db)s \
{{ codestop }}

If you want to annotate potential ligand-binding positions in your sequences instead of domain-binding properties, you can choose to only use Pfams that correspond to domain-ligand interactions that had nonredundant instances across three or more distinct PDB entries and achieved a cross-validated precision of at least 0.5. 

{{ codestart }}
anvi-run-interacdome -c %(contigs-db)s \
                     --interacdome-dataset confident
{{ codestop }}

Additionally, there are three thresholds that you can set: 

1. [`--min-binding-frequency` to ignore very low frequencies](https://merenlab.org/2020/07/22/interacdome/#filtering-low-binding-frequency-scores). The InteracDome scale is from 0 (most likely not involved in binding) to 1 (most likely involved in binding). The default cutoff is 0.200000. 
2. [`--min-hit-fraction` to remove results with low detection]((https://merenlab.org/2020/07/22/interacdome/#filtering-partial-hits)). The default value is 0.5, so HMMs that are less than twice as long as the total hit length will not be considered. 
3. [`--information-content-cutoff` to ignore low-qulaity domain hits](https://merenlab.org/2020/07/22/interacdome/#filtering-bad-hits-with-information-content). The default value is 4, so for an alignment to count, the HMM sequence must be very conserved (with more than 95 percent consensus). Setting this cutoff to a very high value will keep all sequences regardless of percent consensus. 

By default, this program does not produce an output, just puts the final results into the %(misc-data-amino-acids)s (for SAAVs) or %(misc-data-nucleotides)s (for SNVs and SCVs) of your %(contigs-db)s. If you'd like a tab-delimited output, simply provide the `-O` flag and anvi'o will create a %(binding-frequencies-txt)s for your data.  
