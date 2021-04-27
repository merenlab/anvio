
This program predicts per-residue binding scores for genes in your %(contigs-db)s via the [InteracDome](https://interacdome.princeton.edu/) database.


The full process is detailed in [this blog post](https://merenlab.org/2020/07/22/interacdome/). In fact, ideally, all of that information should really be in this very document, but because the blogpost has preceded this document, it hasn't been translated over yet. So really, you should really be reading that blogpost if you want to get into the nitty gritty details. Otherwise, the quick reference herein should be sufficient.


In summary, this program runs an HMM search of the genes in your %(contigs-db)s to all the Pfam gene families that have been annotated with InteracDome binding frequencies. Then, it parses and filters results, associates binding frequencies of HMM match states to the user's genes of interest, and then stores the resulting per-residue binding frequencies for each gene into the %(contigs-db)s as %(misc-data-amino-acids)s.


Before running this program, you'll have to run %(anvi-setup-interacdome)s to set up a local copy of [InteracDome's tab-separated files](https://interacdome.princeton.edu/#tab-6136-4).



## Basic Usage

A basic run of this program looks like this:

{{ codestart }}
anvi-run-interacdome -c %(contigs-db)s -T 4
{{ codestop }}

In addition to storing per-residue binding frequencies as %(misc-data-amino-acids)s in your %(contigs-db)s, this also outputs additional files prefixed with `INTERACDOME` by default (the prefix can be changed with `-O`). These are provided as %(binding-frequencies-txt)s files named `INTERACDOME-match_state_contributors.txt` and `INTERACDOME-domain_hits.txt`. See %(binding-frequencies-txt)s for details.


## Parameters

[InteracDome](https://interacdome.princeton.edu/) offers two different binding frequency datasets that can be chosen with `--interacdome-dataset`.  Choose 'representable' to include Pfams that correspond to domain-ligand interactions that had nonredundant instances across three or more distinct PDB structures. InteracDome authors recommend using this collection to learn more about domain binding properties. Choose 'confident' to include Pfams that correspond to domain-ligand interactions that had nonredundant instances across three or more distinct PDB entries and achieved a cross-validated precision of at least 0.5. The default is 'representable', and you can change it like so:


{{ codestart }}
anvi-run-interacdome -c %(contigs-db)s \
                     --interacdome-dataset confident
{{ codestop }}

This progarm is multi-threaded, so be sure to make use of it:

{{ codestart }}
anvi-run-interacdome -c %(contigs-db)s \
                     --interacdome-dataset confident \
                     -T 8
{{ codestop }}

Additionally, there are numerous thresholds that you can set: 

1. [`--min-binding-frequency` to ignore very low frequencies](https://merenlab.org/2020/07/22/interacdome/#filtering-low-binding-frequency-scores). The InteracDome scale is from 0 (most likely not involved in binding) to 1 (most likely involved in binding). The default cutoff is 0.200000. 
2. [`--min-hit-fraction` to remove poor quality HMM hits]((https://merenlab.org/2020/07/22/interacdome/#filtering-partial-hits)). The default value is 0.5, so at least half of a profile HMM's length must align to your gene, otherwise the hit will be discarded.
3. [`--information-content-cutoff` to ignore low-qulaity domain hits](https://merenlab.org/2020/07/22/interacdome/#filtering-bad-hits-with-information-content). The default value is 4, which means every amino acid of your gene must match the consensus amino acid of the match state for each mate state with [information content](https://en.wikipedia.org/wiki/Sequence_logo) greater than 4. Decreasing this cutoff yields an increasingly stringent filter.


