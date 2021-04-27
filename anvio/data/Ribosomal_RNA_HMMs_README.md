# What is this?

This directory contains HMM profiles and related data for anvi'o to be able to identify bacterial, archaeal, and eukaryotic rRNA genes in contigs.

The source of these profiles is [Torsten Seemann](https://scholar.google.com/citations?user=PuH3Yp4AAAAJ)'s [Barrnap repository](https://github.com/tseemann/barrnap/tree/master/db).

# How to rebuild the genes.hmm.gz file

The HMM files in Torsten's repo is not immediately usable from within anvi'o due to a number of anvi'o- and HMMER-specific non-trivial reasons. So I had to fix a couple of things.

First checkout the repo, and make sure you are on the right time and right place for everything to work:

``` bash
cd /tmp
git clone https://github.com/tseemann/barrnap.git
cd barrnap/
git checkout acf3198a2e671f79107d8221bc32c110a6272ca6
```

Go into the databases directory. You can double check every change that will follow with `git diff`:

``` bash
cd db
```

Remove all DESC lines. HMMER converts those spaces into TABs, extending number of TAB characters arbitrarily:

``` bash
sed -i '' '/^DESC / d' *.hmm
```

Make those names look better and distinguishable:

``` bash
sed -i '' '/^NAME/ s/$/_arc/' arc.hmm
sed -i '' '/^NAME/ s/$/_bac/' bac.hmm
sed -i '' '/^NAME/ s/$/_euk/' euk.hmm
sed -i '' '/^NAME/ s/$/_mit/' mito.hmm
sed -i '' '/^ACC/ s/$/_arc/' arc.hmm
sed -i '' '/^ACC/ s/$/_bac/' bac.hmm
sed -i '' '/^ACC/ s/$/_euk/' euk.hmm
sed -i '' '/^ACC/ s/$/_mit/' mito.hmm

```

We will first split each HMM model into its own file using this Python code:

```
cat <<EOF >>split_hmms.py
import sys
with open(sys.argv[1]) as input_f:
    for section in input_f.read().split('//\n'):
        for line in section.split('\n'):
            if line.find('NAME') > -1:
                name = line.split(' ')[-1]
                with open(name + '.hmm', 'w') as output:
                    output.write(section + '//\n')
                    break
EOF
```

Then we will run it on each file:

```
for f in euk.hmm arc.hmm bac.hmm mito.hmm
do
    python split_hmms.py $f
done
```

This should result in the following files in the working directory:

```
12S_rRNA_mit.hmm  16S_rRNA_bac.hmm  18S_rRNA_euk.hmm  23S_rRNA_bac.hmm  5S_rRNA_arc.hmm  5S_rRNA_euk.hmm    5_8S_rRNA_euk.hmm
16S_rRNA_arc.hmm  16S_rRNA_mit.hmm  23S_rRNA_arc.hmm  28S_rRNA_euk.hmm  5S_rRNA_bac.hmm  5_8S_rRNA_arc.hmm
```

Now it is time to merge them into final models we will use in anvi'o:

```
cat 16S_rRNA_bac.hmm 16S_rRNA_arc.hmm 16S_rRNA_mit.hmm > Ribosomal_RNA_16S.hmm
cat 23S_rRNA_bac.hmm 23S_rRNA_arc.hmm > Ribosomal_RNA_23S.hmm
cat 18S_rRNA_euk.hmm > Ribosomal_RNA_18S.hmm
cat 28S_rRNA_euk.hmm > Ribosomal_RNA_28S.hmm
cat 5S_rRNA_arc.hmm 5S_rRNA_euk.hmm 5S_rRNA_bac.hmm 5_8S_rRNA_euk.hmm 5_8S_rRNA_arc.hmm > Ribosomal_RNA_5S.hmm
cat 12S_rRNA_mit.hmm > Ribosomal_RNA_12S.hmm

rm 12S_rRNA_mit.hmm  16S_rRNA_bac.hmm  18S_rRNA_euk.hmm  23S_rRNA_bac.hmm  5S_rRNA_arc.hmm  5S_rRNA_euk.hmm    5_8S_rRNA_euk.hmm
rm 16S_rRNA_arc.hmm  16S_rRNA_mit.hmm  23S_rRNA_arc.hmm  28S_rRNA_euk.hmm  5S_rRNA_bac.hmm  5_8S_rRNA_arc.hmm
```

At this point, one should make sure that the number of lines between these two are identical:

```
wc -l arc.hmm bac.hmm euk.hmm mito.hmm
wc -l Ribosomal_RNA_*
```


Add noise cutoffs. totally arbitrary :/ maybe at some point we will fix this behavior and the HMM framework in anvi'o will be able to use HMM profiles without any model noise cutoffs. but for now, this is it:

``` bash
sed -i '' '/CKSUM /a \
GA    750 750;\
TC    750 750;\
NC    750 750;\
' Ribosomal*.hmm
```

And ensure unique ACC ids for each HMM entry:

```
```

For each of the Ribosomal RNA classes, create individual directories (which will become anvi'o HMM data directories):

``` bash
for RNA in Ribosomal_RNA_12S  Ribosomal_RNA_16S  Ribosomal_RNA_18S  Ribosomal_RNA_23S  Ribosomal_RNA_28S  Ribosomal_RNA_5S
do
    mkdir $RNA

    mv $RNA.hmm $RNA/genes.hmm

    echo "gene accession hmmsource" > $RNA/genes.txt

    for i in `grep NAME $RNA/genes.hmm | awk '{print $2}'`
    do
        echo "$i None barrnap"
    done >> $RNA/genes.txt

    perl -p -i -e 's/ /\t/g' $RNA/genes.txt

    gzip $RNA/genes.hmm

    echo "--cut_ga" > $RNA/noise_cutoff_terms.txt
    echo $RNA > $RNA/kind.txt
    echo "Seeman T, https://github.com/tseemann/barrnap" > $RNA/reference.txt
    echo "RNA:CONTIG" > $RNA/target.txt
done
```

You are done. Now one can run `anvi-run-hmms` to test one of these:

``` bash
anvi-run-hmms -c CONTIGS.db -H Ribosomal_RNA_16S
```

And get those matchin sequences back:

```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db --hmm Ribosomal_RNA_16S
```
