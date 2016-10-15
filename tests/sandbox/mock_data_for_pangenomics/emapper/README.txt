this is how I run these:

for i in aa_sequences_01 aa_sequences_02 aa_sequences_03; do time emapper.py -i $i.fa --output $i --cpu 40 --database bact --usemem --override; done ; letmeknow
