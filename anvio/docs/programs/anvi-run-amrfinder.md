This program **annotates genes in your %(contigs-db)s with NCBI's [AMRFinderPlus](https://github.com/ncbi/amr/wiki)**, which identifies antimicrobial resistance genes, stress response genes (such as metal, biocide, heat, and acid resistance), and virulence genes. The output is a %(functions)s artifact.

{:.notice}
Keep in mind that AMRs are difficult to annotate and AMRFinderPlus's database is mostly built from AMRs in host associated microbes. To get a comprehensive picture of AMRs in environmental populations you would need to dive much deeper than using this tool only.

Before you can run this program, you will need to install AMRFinderPlus and set up its database with %(anvi-setup-amrfinder)s.

### How does it work?

1. AMRFinderPlus requires the amino acid sequences of the genes in your %(contigs-db)s. The first thing %(anvi-run-amrfinder)s does is to **export the amino acid sequences** of all genes in your %(contigs-db)s into a temporary FASTA file while using their gene caller ids as names.
2. Anvi'o then **runs AMRFinderPlus** which compares the exported amino acid sequences to NCBI's curated reference proteins with BLASTP and to its curated protein families with HMMs. Whether a gene is identified or not depends on NCBI's curated cutoffs. Anvi'o always includes all genes in AMRFinderPlus's database, which also includes stress response genes, virulence genes, and resistance genes with uncertain effects.
3. Anvi'o then **stores the hits** in your %(contigs-db)s including six levels of information for each gene that are available in the AMRFinderPlus database (see below).

Since AMRFinderPlus works with the gene sequences that anvi'o has already identified, the program does not search the contig sequences themselves, and it does not screen for the organism-specific point mutations that the AMRFinderPlus database contains for some known pathogens.

### Running it

{{ codestart }}
anvi-run-amrfinder -c %(contigs-db)s \
                   --num-threads 4
{{ codestop }}

If you set up the %(amrfinder-data)s in a non-default location, you will need to provide that path as well:

{{ codestart }}
anvi-run-amrfinder -c %(contigs-db)s \
                   --amrfinder-data-dir path/to/directory
{{ codestop }}

### What goes into the contigs database

Each positive hit is associated with six levels of information that are contained in the AMRFinderPlus database:

|Annotation source|Accession|Function|Examples|
|:--|:--|:--|:--|
|`AMRFinder`|AMRFinderPlus hierarchy node|gene name|`bla-B3`: subclass B3 metallo-beta-lactamase|
|`AMRFinder_CLASS`|class|class|BETA-LACTAM, MERCURY|
|`AMRFinder_SUBCLASS`|subclass|subclass|CARBAPENEM, ARSENITE|
|`AMRFinder_SUBTYPE`|subtype|subtype|AMR, METAL, BIOCIDE, VIRULENCE|
|`AMRFinder_SCOPE`|scope|scope|core, plus|
|`AMRFinder_METHOD`|method|method|HMM, BLASTP, EXACTP|

AMRFinderPlus annotates each gene following an annotation hierarchy that is built like a family tree of genes, and assigns each hit to the most specific entry (node) of this tree based on the curated similarity cutoffs. Each node in this hierarchy has a unique identifier in the AMRFinderPlus database, which is more specific than gene symbols or names (several different entries share the name 'class A beta-lactamase', for instance).

Genes for which AMRFinderPlus reports no class or subclass, which includes most virulence genes, heat and acid resistance genes, and a few metal resistance genes, have no annotation at the `AMRFinder_CLASS` and `AMRFinder_SUBCLASS` level.

### E-values

Running HMMs to annotate genes generates e-values, but these are not included in the output of AMRFinderPlus because the database uses manually curated cutoff values. When you run the program in anvi'o it recovers the e-values from the hmmsearch table AMRFinderPlus leaves behind in its temporary directory, and stores them with every HMM based gene annotation. Hits without an HMM match that were identified with BLASTP only have no e-value and receive a placeholder e-value of 0. Because AMRFinderPlus calculates e-values for a fixed search space of 10,000 proteins, they are not directly comparable to the e-values of other anvi'o annotation programs.

If for some reason anvi'o cannot find the hmmsearch table (for instance, because a future version of AMRFinderPlus works differently), it will stop with an error. In that case, you can use the flag `--skip-amrfinder-e-values` to store all annotations with an e-value of 0 but still get your annotations. That does not make the annotations less reliable, the same curated cutoffs are still applied during the annotation process:

{{ codestart }}
anvi-run-amrfinder -c %(contigs-db)s \
                   --skip-amrfinder-e-values
{{ codestop }}

### The full AMRFinderPlus report

The %(contigs-db)s stores the essentials of every hit. If you would like to generate an output tab-delimited table that includes all available information for each annotation, you can use this flag:

{{ codestart }}
anvi-run-amrfinder -c %(contigs-db)s \
                   --full-amrfinder-report AMRFinder-report.txt
{{ codestop }}

In addition to what AMRFinderPlus reports for each annotated gene, anvi'o adds the HMM e-value, the HMM bit score, and NCBI's curated bit score cutoff for that HMM (`NA` where there is no HMM information).

### Re-annotating a contigs database

When you run %(anvi-run-amrfinder)s, the AMRFinderPlus database version is recorded in your %(contigs-db)s. To re-annotate the %(contigs-db)s, for example with a newer database version, you need to use the flag `--just-do-it`.

### What can you do next?

AMRFinderPlus annotations are stored as %(functions)s in your %(contigs-db)s and several other anvi'o programs can use them. Here are two examples that are useful when you are working with metagenomic assemblies rather than genomes or MAGs.

If you want to get an idea about the prevalence of specific genes in your metagenomic assemblies you can calculate the **per-population copy number** using %(anvi-gen-function-matrix)s with the `--add-per-population-copy-number` flag:

{{ codestart }}
anvi-gen-function-matrix -e %(external-genomes)s \
                         --annotation-source AMRFinder_CLASS \
                         --add-per-population-copy-number \
                         -O AMR-CLASS
{{ codestop }}

You can also look at the **coverage and detection of genes across samples** if you recruited reads to your genomes, MAGs or assemblies creating a %(profile-db)s. You need to first export the genes with the AMRFinderPlus annotations with %(anvi-export-functions)s to be able to export a file with their coverage and detection with %(anvi-export-gene-coverage-and-detection)s:

{{ codestart }}
anvi-export-functions -c %(contigs-db)s \
                      --annotation-sources AMRFinder,AMRFinder_CLASS,AMRFinder_SUBTYPE \
                      --matrix-format \
                      -o amr-genes.txt
cut -f1 amr-genes.txt | tail -n +2 > amr-gene-ids.txt
anvi-export-gene-coverage-and-detection -c %(contigs-db)s \
                                        -p %(profile-db)s \
                                        --genes-of-interest amr-gene-ids.txt \
                                        -O AMR
{{ codestop }}

Please note that genes on contigs that were too short to be profiled by %(anvi-profile)s will not have coverage values. The default cutoff is 1000 bp, but you can adjust that value when generating the %(profile-db)s.
