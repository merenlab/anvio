How to re-train anvi'o domain classifier
========================================

Anvi'o uses SCGDOMAINCLASSIFIER.rf to predict domains for real time estimation
of completion and redundancy of a given genome bin. This subsystem is commonly
used by many programs including `anvi-interactive` and `anvi-estimate-genome-completeness`.

To re-train the classifier, you need to download the necessary datapack, and use
the program `anvi-script-gen-scg-domain-classifier`. To do that go to any directory:

```
cd ~/Downloads
```

Download the genomes:

```
wget https://ndownloader.figshare.com/files/15230837 -O GENOMES-TO-TRAIN-ANVIO-SCG-DOMAIN-PREDICTOR.tar.gz
tar -zxvf GENOMES-TO-TRAIN-ANVIO-SCG-DOMAIN-PREDICTOR.tar.gz
cd GENOMES-TO-TRAIN-ANVIO-SCG-DOMAIN-PREDICTOR/
```

The training data is organized into sub-directories,

```
ls
archaea  bacteria eukarya
```

Each of which contains *anvi'o contigs databases for single genomes*:

```
ls bacteria/
Acaryochloris_marina_5.db                     Aster_yellows_1606.db
Acetobacter_pasteurianus_13.db                Borrelia_afzelii_3199.db
Acetobacterium_woodii_32.db                   Calditerrivibrio_nitroreducens_4210.db
Acholeplasma_brassicae_35.db                  Candidatus_Nitrospira_defluvii_4551.db
Acidimicrobium_ferrooxidans_105.db            Candidatus_Protochlamydia_amoebophila_4582.db
Acidobacterium_capsulatum_123.db              Chlorobaculum_parvum_4919.db
Aequorivita_sublithincola_1254.db             Defluviitoga_tunisiensis_5522.db
Akkermansia_muciniphila_1366.db               Deinococcus_deserti_5540.db
Aminobacterium_colombiense_1442.db            Fusobacterium_nucleatum_9790.db
Anaerolinea_thermophila_1477.db               Isosphaera_pallida_10718.db
Aquifex_aeolicus_1539.db
```

The purpose of this is to have enough examples for each domain of life, so anvi'o
can learn about how to recognize similar genomes later based on single-copy core
genes.

To train the classifier, you can simply run the following command:

```
anvi-script-gen-scg-domain-classifier --genomes-dir . \
                                      --output /(...)/anvio/data/misc/SCGDOMAINCLASSIFIER.rf
```

How to add a new SCG collection for a new domain
================================================

In addition to the data pack directory, at this point you should have your SCG collection
for this domain in anvi'o HMM directory form. Examples to these directories can be found here:

    https://github.com/merenlab/anvio/tree/master/anvio/data/hmm

First, put your HMM directory to the proper place in the anvi'o codebase you are using.
You can learn the location by running this command:

```
DATA_DIR=`python -c 'import anvio.data.hmm as hmm_data; import os; print(os.path.dirname(hmm_data.__file__))'`

cp -r MY_HMM_DIR/ $DATA_DIR/
```

Once it is in place, you create a new directory in the data pack with the name matching to the
domain name listed in `MY_HMM_DIR/kind.txt` file, and add into this directory new contigs databases
that represent individual genomes that match to that domain.

Once you have the new directory in place, you should run `anvi-run-hmms` on all genomes in the datapack
from scratch. You can do it this way:

```
for i in `find . -name '*.db'`
do
    anvi-run-hmms -c $i -T 6
done
```

Once this is done, you are golden. Run the trainer again,

```
anvi-script-gen-scg-domain-classifier --genomes-dir . \
                                      --output /(...)/anvio/data/misc/SCGDOMAINCLASSIFIER.rf
```

And test one of your new genomes with `anvi-estimate-genome-completeness`.

If you run into any problems let us know!