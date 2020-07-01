An HMM source is a collection of one or more HMMs. 

HMMs for any set of genes can be put together and run on any anvi'o %(contigs-db)s which would yield an %(hmm-hits)s.

HMM hits in a %(contigs-db)s for a given HMM source will be accessible to anvi'o programs globally (i.e., to recover aligned or non-aligned sequences as %(fasta)s files, or display contigs that contain HMM hits in interactive interfaces, to report hits in summary outputs, etc).

Running %(anvi-db-info)s on an anvi'o contigs database will list HMM sources available in it.

### Default HMM sources

An anvi'o installation will include [multiple HMM sources](https://github.com/meren/anvio/tree/master/anvio/data/hmm) by default. These HMM sources can be run on any %(contigs-db)s with %(anvi-run-hmms)s to identify and store %(hmm-hits)s:

{{ codestart }}
%(anvi-run-hmms)s -c %(contigs-db)s
{{ codestop }}

Default HMM sources in anvi'o include,

* **Bacteria_71**: 71 single-copy core genes for domain bacteria that represent a modified version of the HMM profiles published by [Mike Lee](https://doi.org/10.1093/bioinformatics/btz188). The anvi'o collection excludes Ribosomal_S20p, PseudoU_synth_1, Exonuc_VII_S, 5-FTHF_cyc-lig, YidD and Peptidase_A8 occurred in Lee collection (as they were exceptionally redundant or rare among MAGs from various habitats), and includes Ribosomal_S3_C, Ribosomal_L5, Ribosomal_L2 to make it more compatible with [Hug et al](https://www.nature.com/articles/nmicrobiol201648)'s set of ribosomal proteins.
* **Archaea_76**: 76 single-copy core genes for domain archaea by [Mike Lee](https://doi.org/10.1093/bioinformatics/btz188).
* **Protista_83**: 83 single-copy core genes for protists (domain eukarya) by [Tom O. Delmont](http://merenlab.org/delmont-euk-scgs). 
* **Ribosomal_RNAs**: HMMs to identify ribosomal RNA genes for bacteria, archaea, and eukarya by [Torsten Seemann](https://github.com/tseemann/barrnap).

When %(anvi-run-hmms)s is run on an anvi'o %(contigs-db)s without providing any further arguments, it automatically utilizes all the default HMM sources.

{:.notice}
Similar to Ribosomal RNAs, anvi'o can also identify Transfer RNAs. Even though Transfer RNAs will also appear as an HMM source for all downstream analyses, their initial identification will require running %(anvi-scan-trnas)s program on an anvi'o %(contigs-db)s. 

### User-defined HMM sources

The user can employ additional HMM sources to identify matching genes in a given %(contigs-db)s.

Any directory with expected files in it will serve as an HMM source:

{{ codestart }}
anvi-run-hmms -c %(contigs-db)s \
              --hmm-source /PATH/TO/USER-HMM-DIRECTORY/
{{ codestop }}

Anvi'o will expect the HMM source directory to contain six files (see this for [an example directory](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Protista_83)). These files are explicitly defined as follows:

* **genes.hmm.gz**: A gzip of concatenated HMM profiles. One can (1) obtain one or more HMMs by computing them from sequence alignments or by downloading previously computed ones from online resources such as [Pfams](https://pfam.xfam.org/family/browse?browse=new), (2) concatenate all profiles into a single file called `genes.hmm`, and finally (3) compress this file using `gzip`.
* **genes.txt**: A TAB-delimited file that must contain three columns: `gene` (gene name), `accession` (gene accession number (can be anything unique)), and `hmmsource` (source of HMM profiles listed in genes.hmm.gz). The list of gene names in this file must perfectly match to the list of gene names in genes.hmm.gz.
* **kind.txt**: A flat text file which contains a single word identifying what type of profile the directory contains.
* **reference.txt**: A file containing source information for this profile to cite it properly.
* **target.txt**: A file that specifies the target *alphabet* and  *context* that defines how HMMs should be searched (this is a function of the HMM source that is used). The proper notation is 'alphabet:context'. Alphabet can be `AA`, `DNA`, or `RNA`. Context can be `GENE` or `CONTIG`. The content of this file should be any combination of one alphabet and one context term. For instance, if the content of this file is `AA:GENE`, anvi'o will search genes amino acid sequences, and so on. An exception is `AA:CONTIG`, which is an improper target since anvi'o can't translate contigs to amino acid sequences. See [this](https://github.com/meren/anvio/pull/402) for more details.

{:.warning}
All `DNA:CONTIG` targets will add new genes in the database for each hit.

* **noise_cutoff_terms.txt**: A file to specify how to deal with noise. [See this comment](https://github.com/merenlab/anvio/issues/498#issuecomment-362115921) for more information on the contents of this file.
