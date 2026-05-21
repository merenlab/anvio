%(anvi-gen-pan-representative)s generates a pangenome supplemented representative genome from a set of supplied genomes and an already computed pangenome. 

The supplied genomes should be in a fasta file that **must have simple deflines**. If your fasta file doesn’t have simple deflines, it is not a proper %(contigs-fasta)s. To fix your FASTA file use %(anvi-script-reformat-fasta)s.

Using the %(contigs-fasta)s, contigs.db need to be generated using %(anvi-gen-contigs-database)s. If you want your pangenome supplememented representative genome to include annotations, the next step is to annotate your %(contigs-db)s using any or all of the following commands:

- %(anvi-run-kegg-kofams)s
- %(anvi-run-ncbi-cogs)s
- %(anvi-run-hmms)s
- %(anvi-run-scg-taxonomomy)s

The annotated %(contigs-db)s will then all be used to generate a %(genomes-storage-db) **which stores information about your genomes, primarily for use in pangenomic analysis.**

Once you’ve generated a %(genomes-storage-db)s, you can run %(anvi-pan-genome)s, which creates a %(pan-db)s and runs various pangenomic analyses including calculating the similarities between your sequences, identifying gene clusters, and organizing your gene clusters and genomes.

Additionally, an %(external-genomes)s file is required which can be produced using %(anvi-script-gen-genomes-file)s. 

%(anvi-gen-pan-representative)s will produce a %(contigs-db)s variant `pan-genome`. The %(anvi-gen-pan-representative)s includes the whole genome which is selected as representative of the supplied genomes, supplemented with contig that includes a representative gene for each of the gene clusters in the pangenome that are not already present in the representative genome.

The default usage is simple, which will simply instruct anvi’o to select a representative with the best combination between lower contamination, higher completeness and less number of contigs. The tool will keep the contigs of the representative genome, and then select from within the remaining genomes for the genes that represent the gene clusters that haven't been covered and contatenates them in a supplementary contig. 

{{ codestart }}
%%(anvi-gen-pan-representative)s -p %(pan-db)s -g %(genomes-storage-db)s -e %(external-genomes)s \
                           -o PATH/TO/PAN-GENOME/
{{ codestop }} 

### Different options to generate the pangenome supplemented representative genome

--gap-size In the pangenome supplement contig the sequences are joined by a gap of 20 N. Use this flag to change the size of the gap.

--representative Use this flag to manually select the representative genome. If not given, the program will choose the genome representative using the alpha value.

--alpha Is a weight to define the priority for definition of representative genome. It weights contamination, completeness and number of contigs in the following formula: alpha * (completeness - contamination) + (1 - alpha) * (completeness - contamination). So an alpha value of 1, only considers completeness and contamination, and an alpha value of 0 only considers number of contigs. The default alpha value if not given is 0.8

--max-num-contigs This can be used to select a representative genome with a user defined maximal number of contigs. If not given, the program will use alpha to select a representative

--min-genomes-per-gene-cluster Allow the user to discard gene clusters to be added to the supplementary contig if a minimal number of given genomes is not represented. Defafult is 1 to keep all gene clusters. 

--keep-synteny By default, the supplementary contig is built by adding a representative gene from start codon to stop codon into the contig, joined by 20 N even if two genes are next to each other in the original contig. Keep synteny option will keep from start codon of the first gene to the stop codon of the last gene that sit together and have not been represented yet. Additionally, it will keep as much synteny as possible when selecting missing genes in the supplemental contig.

--keep-promoter Keeps the flanking regions of the genes that will be added to the supplemental contig. By definition, it also keeps synteny.