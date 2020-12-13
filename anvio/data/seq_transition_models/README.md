How to regenerate the contents of this directory
================================================

FIXME FIXME

## Note

For a briefing on what the purpose of this folder is, please review this pull request:
https://github.com/merenlab/anvio/pull/1428.  In brief, these files are numpy arrays of Markov
models used to identify which items (nucleotides, amino acids, codons) are likely to follow one
another in a string of sequence.

## Reproducing the AA dir

These models (the `.npy` files) were generated from a large database called
[UniRef50](https://www.uniprot.org/help/uniref). To recreate this directory, move to the directory
of this markdown file. You'll need 15GB of space to store the database (You can delete this file
after the `.npy` files have been generated). Then...

```
mkdir -p AA
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
gunzip -d uniref50.fasta.gz
```

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

Then run the script to generate various orders of the models. **Each run will take ~5 hours**. If that's too long,
consider using only a portion of the fasta file (e.g. `head -n 50000 uniref50.fasta >
uniref50_subset.fasta`):

```
python gen_transition_matrix.py --fasta uniref50.fasta --order 1 --output AA/first_order.npy
python gen_transition_matrix.py --fasta uniref50.fasta --order 2 --output AA/second_order.npy
python gen_transition_matrix.py --fasta uniref50.fasta --order 3 --output AA/third_order.npy
python gen_transition_matrix.py --fasta uniref50.fasta --order 4 --output AA/fourth_order.npy
```

If you run into any problems let us know!
