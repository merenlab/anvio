This program can work with anvi'o %(contigs-db)s, %(external-genomes)s, or %(internal-genomes)s files to return sequences for HMM hits identified through the default anvi'o %(hmm-source)ss (such as the domain-specific single-copy core genes) or user-defined %(hmm-source)ss (such as HMMs for specific antibiotic resistance gene families or any other targets).

Using it with single-copy core genes in default anvi'o HMMs make it a very versatile tool for phylogenomics as the user can define specific sets of genes to be aligned and concatenated.


### Learn available HMM sources

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                --list-hmm-sources

AVAILABLE HMM SOURCES
===============================================
* 'Bacteria_71' (type: singlecopy; num genes: 71)
* 'Archaea_76' (type: singlecopy; num genes: 76)
* 'Protista_83' (type: singlecopy; num genes: 83)
* 'Ribosomal_RNAs' (type: Ribosomal_RNAs; num genes: 12)
{{ codestop }}

### Get all sequences in a given HMM source

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                --hmm-source Bacteria_71 \
                                -o %(genes-fasta)s
{{ codestop }}

### Learn available genes in a given HMM source

Please note that the flag `--list-available-gene-names` will give you the list of genes in an **HMM collection** (for example, for `Bacteria_71` in the following use case), and it will not give you the list of genes in your genomes or metagenomes that are matching to them. You can generate a table of HMMs across your genomes or metagenomes with another program, %(anvi-script-gen-hmm-hits-matrix-across-genomes)s.

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                --hmm-source Bacteria_71 \
                                --list-available-gene-names

* Bacteria_71 [type: singlecopy]: ADK, AICARFT_IMPCHas, ATP-synt, ATP-synt_A,
Chorismate_synt, EF_TS, Exonuc_VII_L, GrpE, Ham1p_like, IPPT, OSCP, PGK,
Pept_tRNA_hydro, RBFA, RNA_pol_L, RNA_pol_Rpb6, RRF, RecO_C, Ribonuclease_P,
Ribosom_S12_S23, Ribosomal_L1, Ribosomal_L13, Ribosomal_L14, Ribosomal_L16,
Ribosomal_L17, Ribosomal_L18p, Ribosomal_L19, Ribosomal_L2, Ribosomal_L20,
Ribosomal_L21p, Ribosomal_L22, Ribosomal_L23, Ribosomal_L27, Ribosomal_L27A,
Ribosomal_L28, Ribosomal_L29, Ribosomal_L3, Ribosomal_L32p, Ribosomal_L35p,
Ribosomal_L4, Ribosomal_L5, Ribosomal_L6, Ribosomal_L9_C, Ribosomal_S10,
Ribosomal_S11, Ribosomal_S13, Ribosomal_S15, Ribosomal_S16, Ribosomal_S17,
Ribosomal_S19, Ribosomal_S2, Ribosomal_S20p, Ribosomal_S3_C, Ribosomal_S6,
Ribosomal_S7, Ribosomal_S8, Ribosomal_S9, RsfS, RuvX, SecE, SecG, SecY, SmpB,
TsaE, UPF0054, YajC, eIF-1a, ribosomal_L24, tRNA-synt_1d, tRNA_m1G_MT,
Adenylsucc_synt
{{ codestop }}

### Get sequences for some sequences in a given HMM source

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                -o %(genes-fasta)s
{{ codestop }}

### Get HMM hits in bins of a collection

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                -p %(profile-db)s \
                                -C %(collection)s
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                -o %(genes-fasta)s
{{ codestop }}

### Get amino acid sequences for HMM hits

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                -p %(profile-db)s \
                                -C %(collection)s
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --get-aa-sequences \
                                -o %(genes-fasta)s
{{ codestop }}

### Get HMM hits independently aligned and concatenated

The resulting file can be used for phylogenomics analyses via %(anvi-gen-phylogenomic-tree)s or through more sophisticated tools for curating alignments and computing trees.

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                -p %(profile-db)s \
                                -C %(collection)s
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --get-aa-sequences \
                                --concatenate-genes \
                                --return-best-hit
                                -o %(genes-fasta)s
{{ codestop }}

Please note teh presence of a new flag in this particular command line, `--return-best-hit`. This flag is most appropriate if one wishes to perform phylogenomic analyses, which ensures that for any given protein family, there will be only one gene reported from a given genome. This is necessary due to the nature of the data that goes into phylogenomic analyses, where typically multiple single-copy core genes from each genome are individually aligned and the results are concatenated into a super matrix for tree construction. This requirement will be violated *if* for a given single-copy core gene (SCG) family any given genome in the dataset has two or more genes rather than one, which can happen for a variety of technical or biological reasons. In that case, we need to pick only one of those genes, which is exactly what `--return-best-hit` flag does for us. Let's say we have two `Ribosomal_L3` gene hits in a given genome. When declared, this flag will choose the `Ribosomal_L3` gene that has the most significant hit given the hidden Markov model for `Ribosomal_L3` that was used to search for `Ribosomal_L3` genes in genomes. In cases where genome quality is sufficient and contamination is not a considerable risk, this step will choose the right hit as in many cases of multiple hits for SCGs the additional ones will have very low significance. Essentially, `--return-best-hit` makes sure you are working with the most appropriate genes for phylogenomics given the HMM modesl and significance scores for your matches in your genomes.

## Tips

### Get amino acid seqeunces for each gene in a model individually

If you are interested in recovering HMM hits for each gene in a model anvi'o knows about as a separate FASTA file, you can do it with a `for` loop easily. After learning your genes of interest, first run this to make sure your terminal environment knows about them (this is an example with a few genes from the HMM source `Bacteria_71`, but you can add as many genes as you like and use any HMM source anvi'o recognizes, of course):

``` bash
export genes="Ribosomal_L22 Ribosomal_L23 Ribosomal_L27 Ribosomal_L27A Ribosomal_L28"
export hmm_source="Bacteria_71"
```

Then, you can run the program in a loop to have your FASTA files:

``` bash
for gene in $genes
do
    anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                    --hmm-source $hmm_source \
                                    --gene-name $gene \
                                    -o ${hmm_source}-${gene}.fa
done
```

Voila!

### Exercise with the program or test scenarios

You can play with this program using the anvi'o data pack for the [infant gut data](/tutorials/infant-gut) and by replacing the parameters above with appropriate ones in the following commands.

Download the latest version of the data from here: [doi:10.6084/m9.figshare.3502445](https://doi.org/10.6084/m9.figshare.3502445)

Then, unpack it:

{{ codestart }}
tar -zxvf INFANTGUTTUTORIAL.tar.gz && cd INFANT-GUT-TUTORIAL
{{ codestop }}

Finally, import the collection `merens`:

{{ codestart }}
%(anvi-import-collection)s additional-files/collections/merens.txt \
                       -p PROFILE.db \
                       -c CONTIGS.db \
                       -C merens
{{ codestop }}

Then run the program using the `PROFILE.db`, `CONTIGS.db`, and optionally the %(collection)s `merens` to try some of the commands shown on this page.
