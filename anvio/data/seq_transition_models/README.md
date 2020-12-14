How to regenerate the contents of this directory
================================================

## Note

For a briefing on what the purpose of this folder is, please review these pull requests:

https://github.com/merenlab/anvio/pull/1428
https://github.com/merenlab/anvio/pull/1590

In brief, these files are numpy arrays of Markov models used to identify which items (nucleotides,
amino acids, codons) are likely to follow one another in a string of sequence.

## Reproducing the AA dir

These models (the `.npy` files) were generated from a large database called
[UniRef50](https://gtdb.ecogenomic.org/). I used ~23K bacterial genomes from v89, each which had
their own contigs database. If you want to generate your own Markov model, you can follow similar
steps with your own set of contigs databases.

To begin, I created this file:

```python
#! /usr/bin/env python

import anvio
import anvio.db as db
import anvio.fastalib as u

import argparse

from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("-c", "--contigs-db-paths")
ap.add_argument("-o", "--output-fasta")

args = ap.parse_args()

# ---------------------------------------

output = u.FastaOutput(args.output_fasta)
contigs_db_paths = [x.strip() for x in open(args.contigs_db_paths, 'r').readlines()]
gc_dist = []

count = 0
for contigs_db_path in contigs_db_paths:
    print(count)

    try:
        cdb = db.DB(contigs_db_path, None, ignore_version=True)
        name = Path(contigs_db_path).stem

        basic_info = cdb.get_table_as_dataframe('contigs_basic_info')
        gc = (basic_info['gc_content'] * basic_info['length'] / basic_info['length'].sum()).sum()
        gc_dist.append(gc)

        seqs = cdb.get_table_as_dataframe('gene_amino_acid_sequences')
    except:
        count += 1
        continue

    for i, row in seqs.iterrows():
        if row['sequence'] is not '':
            output.write_id(f"{name}|{row['gene_callers_id']}|{gc:.4f}")
            output.write_seq(row['sequence'], split = False)

    count += 1

output.close()

with open('gc_distribution', 'w') as f:
    f.write('\n'.join([str(gc) for gc in gc_dist]))
```

and ran it with `python create_aa_seq_db.py -c db_paths -o aa_seq_db.fa`, where `db_paths` is the
23k db paths, one per line. The whole thing took 3 hours. As a strange artifact of my investigation,
the best Markov model is created using only the genes from these genomes that have average GC
content < 39%.  The pull request https://github.com/merenlab/anvio/pull/1590 details the evidence
for this, though it is still a mystery why. To isolate these protein sequences, we split up
`aa_seq_db.fa` with the following script:


```python
#! /usr/bin/env python

import anvio
import anvio.db as db
import anvio.fastalib as u

import argparse
import numpy as np

from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--fasta")
ap.add_argument("-b", "--bins")

args = ap.parse_args()

# ---------------------------------------

bins = np.array([float(x) for x in args.bins.split(',')])
prefix = Path(args.fasta).stem

fasta = u.SequenceSource(args.fasta)
outputs = {k: u.FastaOutput(prefix + f'_GC_{bins[k-1]:.4f}-{bins[k]:.4f}.fa') for k in range(1, len(bins))}

while next(fasta):
    b = int(np.digitize(float(fasta.id.split('|')[-1]), bins))
    outputs[b].store(fasta)

fasta.close()
for output in outputs:
    outputs[output].close()
```

And run it with `python partition_aa_seq_db.py -f aa_seq_db.fa -b 0,0.387,0.46,0.561,0.645,1.00`

Half an hour later, you end up with a couple files, most important being `aa_seq_db_GC_0.0000-0.3870.fa`.

Ok so you've got the FASTA. Take a look if you want. When you're ready, create the following script,
naming it `gen_transition_matrix.py` and putting in your current working dir:

```python
#! /usr/bin/env python

import numpy as np
import argparse
import anvio.fastalib as u
import anvio.constants as c

ap = argparse.ArgumentParser()
ap.add_argument('--fasta', '-f', required=True, help="FASTA to generate transition matrix from.")
ap.add_argument('--order', '-n', required=True, type=int, help="What order? Set -n 1 to include self. -n 2 to include nearest neighbor, etc.")
ap.add_argument('--output', '-o', required=True, help="Output matrix")
args = ap.parse_args()

fasta = u.SequenceSource(args.fasta)
aas = [c.AA_to_single_letter_code[aa] for aa in c.amino_acids if aa != 'STP']
aa_to_array_index = {aa: i for i, aa in enumerate(aas)}
num_aas = len(aas)

matrix = np.zeros(tuple([num_aas] * args.order))

while next(fasta):
    seq = fasta.seq
    for i in range(len(fasta.seq) - args.order):
        try:
            state = tuple([aa_to_array_index[aa] for aa in seq[i:i+args.order]])
            matrix[state] += 1
        except KeyError:
            # Catches bizzarre sequence characters
            pass

matrix /= matrix.sum()
np.save(args.output, matrix)
```

Then run the script to generate various orders of the model. It will take around 3 hours. If that's too long,
consider using only a portion of the fasta file (e.g. `head -n 50000 aa_seq_db_GC_0.0000-0.3870.fa >
uniref50_subset.fasta`). For example, the model currently used was generated from:

```
python gen_transition_matrix.py --fasta aa_seq_db_GC_0.0000-0.3870.fa --order 1 --output AA/MM_GC_0-39.npy
```

If you run into any problems let us know!
