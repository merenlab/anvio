# What is this?

This directory contains HMM profiles and related data for anvi'o to be able to identify bacterial, archaeal, and eukaryotic rRNA genes in contigs.

The source of these profiles is [Torsten Seemann](https://scholar.google.com/citations?user=PuH3Yp4AAAAJ)'s [Barrnap repository](https://github.com/tseemann/barrnap/tree/master/db).

# How to rebuild the genes.hmm.gz file

The HMM files in Torsten's repo is not immediately usable from within anvi'o due to a number of anvi'o- and HMMER-specific non-trivial reasons. So I had to fix a couple of things.

First checkout the repo, and make sure you are on the right time and right place for everything to work:

``` bash
cd /tmp
git clone https://github.com/tseemann/barrnap.git
git checkout ebfdc202842b4ec16ac6b3a380b17f2e00ab6b68
```

Go into the databases directory. You can double check every change that will follow with `git diff`:

``` bash
cd barrnap/db
```

Remove all DESC lines. HMMER converts those spaces into TABs, extending number of TAB characters arbitrarily:

``` bash
sed -i '' '/^DESC / d' *.hmm
```

Make those names look better and distinguishable:

``` bash
sed -i '' 's/NAME  /NAME  Bacterial_/g' bac.hmm
sed -i '' 's/NAME  /NAME  Archaeal_/g' arc.hmm
sed -i '' 's/NAME  /NAME  Mitochondrial_/g' mito.hmm
sed -i '' 's/NAME  /NAME  Eukaryotic_/g' euk.hmm
```

Add noise cutoffs. totally arbitrary :/ maybe at some point we will fix this behavior and the HMM framework in anvi'o will be able to use HMM profiles without any model noise cutoffs. but for now, this is it:

``` bash
sed -i '' '/CKSUM /a \
GA    750 750;\
TC    750 750;\
NC    750 750;\
' *.hmm
```

Create a new directory (which will become all the anvi'o data):

``` bash
mkdir Ribosomal_RNAs
```

Concatenate all into one file:

``` bash
cat *.hmm > Ribosomal_RNAs/genes.hmm
```

Create a `genes.txt` file:

```bash
echo "gene accession hmmsource" > Ribosomal_RNAs/genes.txt

for i in `grep NAME Ribosomal_RNAs/genes.hmm | awk '{print $2}'`
do
    echo "$i None barrnap"
done >> Ribosomal_RNAs/genes.txt

perl -p -i -e 's/ /\t/g' Ribosomal_RNAs/genes.txt

gzip Ribosomal_RNAs/genes.hmm
```

Create other necessary files:

``` bash
echo "Ribosomal_RNAs" > Ribosomal_RNAs/kind.txt
echo "Seeman T, https://github.com/tseemann/barrnap" > Ribosomal_RNAs/reference.txt
echo "RNA:CONTIG" > Ribosomal_RNAs/target.txt
```

You are done. Now one can run it with `anvi-run-hmms`:

``` bash
anvi-run-hmms -c CONTIGS.db -H Ribosomal_RNAs
```

And get those sequences back:

```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db --hmm Ribosomal_RNAs
```