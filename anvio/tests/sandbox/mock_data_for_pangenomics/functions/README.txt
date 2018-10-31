So this is how I run the EggNOGs:

    for i in aa_sequences_01 aa_sequences_02 aa_sequences_03; do time emapper.py -i $i.fa --output $i --cpu 40 --database bact --usemem --override; done ; letmeknow

Then imported them into contigs databases:

    anvi-script-run-eggnog-mapper -c 01.db --annotation aa_sequences_01.emapper.annotations --use-version 0.12.6
    anvi-script-run-eggnog-mapper -c 01.db --annotation aa_sequences_01.emapper.annotations --use-version 0.12.6
    anvi-script-run-eggnog-mapper -c 01.db --annotation aa_sequences_01.emapper.annotations --use-version 0.12.6

Then I exported them back as simple matrices:

    anvi-export-functions -c 01.db -o 01-functions.txt
    anvi-export-functions -c 02.db -o 02-functions.txt
    anvi-export-functions -c 03.db -o 03-functions.txt

Why? Because we want to include functions in our tests, but we don't
want to require people to have COGs setup for emapper or anything. This
use of simple matrices simplifies things :/
