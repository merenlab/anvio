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

### Change the formatting of the output FASTA files

Please note that this program allows you to format the deflines of the resulting FASTA file to a great extent whenever possible. For this, it uses a set of previously-defined variables you can use to define a template of your liking. You can learn about the available variables, you can include the following flag in your command:

{{ codestart }}
anvi-get-sequences-for-hmm-hits -c %(contigs-db)s \
                                --list-defline-variables
{{ codestop }}

Which will give you an output similar to the one below:

```
WARNING
===============================================
Here are the variables you can use to provide a user-defined defline template:

* {gene_name}
* {gene_callers_id}
* {contig_name}
* {gene_unique_id}
* {bin_name}
* {source}
* {e_value}
* {start}
* {stop}
* {length}

Remember, by default, anvi'o will use the following template to format the
deflines of FASTA files it produces whenever possible.

{gene_name}___{gene_unique_id} bin_id:{bin_name}|source:{source}|e_value:{e_value}|contig:{contig_name}|gene_callers_id:{gene_callers_id}|start:{start}|stop:{stop}|length:{length}
```

With this default template, %(anvi-get-sequences-for-hmm-hits)s will provide a FASTA file with quite a comprehensive defline. The following examples use the data pack for the Infant Gut Dataset, after running the following command in the data directory:

{{ codestart }}
%(anvi-import-collection)s additional-files/collections/merens.txt \
                           -p PROFILE.db \
                           -c CONTIGS.db \
                           -C merens
{{ codestop }}

Here are a few commands that show how different deflines impact the FASTA output, starting with the default defline format:

```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -o OUTPUT.fa
```

```
>RsfS___Bacteria_71___a2cb9835d40f1ea052bc7aea2ef8c12924c76f4fa00e00b69f7dbecb bin_id:CONTIGS|source:Bacteria_71|e_value:2.6e-14|contig:Day17a_QCcontig1|gene_callers_id:22|start:17642|stop:17792|length:150
>SecE___Bacteria_71___fe4bf2883dfee62dcf2ed93acde8df5de5be421675d7231997dc5f91 bin_id:CONTIGS|source:Bacteria_71|e_value:3.6e-22|contig:Day17a_QCcontig1|gene_callers_id:105|start:95904|stop:96075|length:171
>Ribosomal_L1___Bacteria_71___9a3cc2232443c0c0e995fef5f76f34d52fb111476e0e2fa2dbe0f09e bin_id:CONTIGS|source:Bacteria_71|e_value:8.4e-59|contig:Day17a_QCcontig1|gene_callers_id:116|start:106378|stop:107068|length:690
>SecG___Bacteria_71___d23b8716a0b03543405ac9835ce127b981fc1c067d9b1e48edf9d861 bin_id:CONTIGS|source:Bacteria_71|e_value:1e-20|contig:Day17a_QCcontig1|gene_callers_id:205|start:200881|stop:201118|length:237
>SmpB___Bacteria_71___a95dd711a5d10cfb4bd9ca3d836aed077d35647161d4908a9c8252fe bin_id:CONTIGS|source:Bacteria_71|e_value:2.8e-66|contig:Day17a_QCcontig1|gene_callers_id:208|start:204418|stop:204883|length:465
```

---

```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -o OUTPUT.fa \
                                -p PROFILE.db \
                                -C merens
```

```
>RsfS___Bacteria_71___a2cb9835d40f1ea052bc7aea2ef8c12924c76f4fa00e00b69f7dbecb bin_id:E_facealis|source:Bacteria_71|e_value:2.6e-14|contig:Day17a_QCcontig1|gene_callers_id:22|start:17642|stop:17792|length:150
>SecE___Bacteria_71___fe4bf2883dfee62dcf2ed93acde8df5de5be421675d7231997dc5f91 bin_id:E_facealis|source:Bacteria_71|e_value:3.6e-22|contig:Day17a_QCcontig1|gene_callers_id:105|start:95904|stop:96075|length:171
>Ribosomal_L1___Bacteria_71___9a3cc2232443c0c0e995fef5f76f34d52fb111476e0e2fa2dbe0f09e bin_id:E_facealis|source:Bacteria_71|e_value:8.4e-59|contig:Day17a_QCcontig1|gene_callers_id:116|start:106378|stop:107068|length:690
>SecG___Bacteria_71___d23b8716a0b03543405ac9835ce127b981fc1c067d9b1e48edf9d861 bin_id:E_facealis|source:Bacteria_71|e_value:1e-20|contig:Day17a_QCcontig1|gene_callers_id:205|start:200881|stop:201118|length:237
>SmpB___Bacteria_71___a95dd711a5d10cfb4bd9ca3d836aed077d35647161d4908a9c8252fe bin_id:E_facealis|source:Bacteria_71|e_value:2.8e-66|contig:Day17a_QCcontig1|gene_callers_id:208|start:204418|stop:204883|length:465
```

---

```
anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                -c CONTIGS.db \
                                -C merens \
                                -o OUTPUT.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3
```

```
>Ribosomal_L27___Bacteria_71___477a28bece0b76ebaef1a9415cbd15b422ec581619e437a79f227f2f bin_id:E_facealis|source:Bacteria_71|e_value:4.1e-37|contig:Day17a_QCcontig2|gene_callers_id:1130|start:86504|stop:86792|length:288
>Ribosomal_L3___Bacteria_71___611453ec6ee187ebae05815acd0d20efa3be77b791a69055aa9ae776 bin_id:S_epidermidis|source:Bacteria_71|e_value:1.4e-19|contig:Day17a_QCcontig7|gene_callers_id:2339|start:19929|stop:20592|length:663
>Ribosomal_L3___Bacteria_71___3e44e2b950f135b147b3a5ca7813581228c3d91c7d2cb69c8543237d bin_id:E_facealis|source:Bacteria_71|e_value:1.5e-20|contig:Day17a_QCcontig16|gene_callers_id:3080|start:229935|stop:230565|length:630
>Ribosomal_L3___Bacteria_71___0d18b80b3c3eb77af3683e64fa4278ac08929c53c8a4593c205b8849 bin_id:P_rhinitidis|source:Bacteria_71|e_value:9.1e-17|contig:Day17a_QCcontig21|gene_callers_id:3313|start:54724|stop:55357|length:633
>Ribosomal_L28___Bacteria_71___e36828b18440c369bcf880a267cea6d55fa65cf1b0b1d66c8784901a bin_id:E_facealis|source:Bacteria_71|e_value:2.2e-24|contig:Day17a_QCcontig23|gene_callers_id:3633|start:124223|stop:124412|length:189
>Ribosomal_L27___Bacteria_71___bd93789ee7b2f24fefa4ddfc0949be06573a6f6f3cdc85a6ada54323 bin_id:S_epidermidis|source:Bacteria_71|e_value:1e-36|contig:Day17a_QCcontig29|gene_callers_id:4151|start:9097|stop:9382|length:285
>Ribosomal_L27___Bacteria_71___2b6e4f252a14a16375d107711cde0363f266f1981fa353a2421fd989 bin_id:S_aureus|source:Bacteria_71|e_value:7.2e-37|contig:Day17a_QCcontig56|gene_callers_id:6123|start:46552|stop:46837|length:285
>Ribosomal_L27___Bacteria_71___51f9c99d5530e7342657805b778105e8e75c178b1780140163e528e3 bin_id:P_rhinitidis|source:Bacteria_71|e_value:3.8e-36|contig:Day17a_QCcontig58|gene_callers_id:6253|start:92104|stop:92401|length:297
>Ribosomal_L28___Bacteria_71___4be054e14a9cb7f80dea306fbbcb391a195acda0a705853b6b010a2e bin_id:P_avidum|source:Bacteria_71|e_value:1.7e-26|contig:Day17a_QCcontig60|gene_callers_id:6345|start:76696|stop:76933|length:237
```

---

```
anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                -c CONTIGS.db \
                                -C merens \
                                -o OUTPUT.fa \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --defline-format "{"
```

```
Init .........................................: 4451 splits in 13 bin(s)


Config Error: Your f-string syntax is not working for anvi'o :/ Perhaps you forgot to open or
              close a curly bracket?
```

---

```
anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                -c CONTIGS.db \
                                -C merens \
                                -o OUTPUT.fa \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --defline-format "{lol}"
```

```
Config Error: Some of the variables in your f-string does not occur in the source dictionary
              :/ Here is the list of those that are not matching to anything: lol. In the
              meantime, these are the known keys: gene_name, gene_callers_id, contig_name,
              gene_unique_id, bin_name, source, e_value, start, stop, length.
```

---

```
anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                -c CONTIGS.db \
                                -C merens \
                                -o OUTPUT.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --return-best-hit \
                                --defline-format "{bin_name}_{gene_callers_id}"
```

```
>E_facealis_1130
>S_epidermidis_2339
>E_facealis_3080
>P_rhinitidis_3313
>E_facealis_3633
>S_epidermidis_4151
>S_aureus_6123
>P_rhinitidis_6253
```

---

```
anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                -c CONTIGS.db \
                                -C merens \
                                -o OUTPUT.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --return-best-hit \
                                --defline-format "{bin_name}_{source}_{gene_callers_id}"
```

```
>E_facealis_Bacteria_71_1130
>S_epidermidis_Bacteria_71_2339
>E_facealis_Bacteria_71_3080
>P_rhinitidis_Bacteria_71_3313
```

---

```
anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                -c CONTIGS.db \
                                -C merens \
                                -o OUTPUT.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                --return-best-hit \
                                --defline-format "{gene_name}_{gene_callers_id} source:{source}|contig:{contig_name}|start:{start}|stop:{stop}"
```

```
>Ribosomal_L27_1130 source:Bacteria_71|contig:Day17a_QCcontig2|start:86504|stop:86792
>Ribosomal_L3_2339 source:Bacteria_71|contig:Day17a_QCcontig7|start:19929|stop:20592
>Ribosomal_L3_3080 source:Bacteria_71|contig:Day17a_QCcontig16|start:229935|stop:230565
>Ribosomal_L3_3313 source:Bacteria_71|contig:Day17a_QCcontig21|start:54724|stop:55357
```

---

Please note that anvi'o will not check whether your defline format will result in FASTA entries with identical deflines.


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

### Performance optimization for aligners by passing additional params

When `--concatenate-genes` flag is used for phylogenomics applications, %(anvi-gen-phylogenomic-tree)s relies on sequence alignment using `muscle` by default. For very large number of sequences this step may fail due to various reasons, such as running out of memory, exceeding the time allocated for the job, etc. If you are having such performance issues, you may want to pass additional parameters to the aligner. For this, you can use BASH environmental variables. For instance, if you wish `muscle` to do only two iterations of alignment and stop after, you can pass that request to the anvi'o driver for `muscle` the following way by exporting a shell varaible called `MUSCLE_PARAMS`:

``` bash
# export a shell variable with additional params you learned from
# the help menu of muscle:
export MUSCLE_PARAMS="-maxiters 2"

# then run anvi-get-sequences-for-hmm-hits exactly the same way
anvi-get-sequences-for-hmm-hits (...)
```

If you are trying to make sure things are going the way you expect, feel free to turn on the debug outputs by adding `--debug` to your command line. This will allow you to see exactly what commands are running behind hte scenes, and will keep the temporary directories for each alignment so you can find the log files in them to see raw outputs from `muscle`. See the anvi'o GitHub issue [#2200](https://github.com/merenlab/anvio/issues/2200) for an extreme case and example ways to debug the process with example commands and outputs.

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
