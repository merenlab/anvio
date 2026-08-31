#!/usr/bin/env python3
"""Build per-sample variant FASTAs from a real source genome containing a known transposase.

Reads:
    source.fa                : the original genome fragment (single contig)
    transposase_coords.txt   : two tab-separated ints on one line — half-open [start, end)
    variants.txt             : TSV with `sample_name<TAB>mge_position`, where
                               `mge_position` is the insertion position in the NAKED-BACKBONE
                               coordinate frame (i.e. source with the transposase removed)

Writes:
    01_FASTA/{sample}.fa     : backbone[:pos] + mge + backbone[pos:], renamed to a stable
                               contig name so all variants share the same name (required for
                               coherent BAMs once mapped to one of them as reference).

For S01 with `mge_position == transposase_start`, the output reconstructs source.fa byte-for-byte.
"""

import os
import textwrap


CONTIG_NAME = 'c_000000000001'


def read_fasta_seq(path):
    with open(path) as f:
        parts = []
        for line in f:
            if line.startswith('>'):
                continue
            parts.append(line.strip())
    return ''.join(parts)


def write_fasta(path, name, seq, wrap=70):
    with open(path, 'w') as f:
        f.write(f">{name}\n")
        for chunk in textwrap.wrap(seq, wrap):
            f.write(chunk + "\n")


def main():
    here = os.path.dirname(os.path.abspath(__file__))

    source = read_fasta_seq(os.path.join(here, 'source.fa'))

    with open(os.path.join(here, 'transposase_coords.txt')) as f:
        start, end = (int(x) for x in f.readline().strip().split('\t'))
    assert 0 <= start < end <= len(source), \
        f"transposase coords [{start}, {end}) out of range for source of length {len(source)}"

    mge = source[start:end]
    backbone = source[:start] + source[end:]

    print(f"source length: {len(source)}; transposase [{start}, {end}) = {len(mge)} bp; naked backbone: {len(backbone)} bp")

    out_dir = os.path.join(os.getcwd(), '01_FASTA')
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(here, 'variants.txt')) as f:
        header = f.readline().rstrip('\n').split('\t')
        assert header == ['sample_name', 'mge_position'], f"unexpected variants.txt header: {header}"
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            sample_name, pos = line.split('\t')
            pos = int(pos)
            assert 0 <= pos <= len(backbone), f"insertion position {pos} for {sample_name} out of bounds for backbone of length {len(backbone)}"
            variant_seq = backbone[:pos] + mge + backbone[pos:]
            out_path = os.path.join(out_dir, f"{sample_name}.fa")
            write_fasta(out_path, CONTIG_NAME, variant_seq)
            print(f"  {sample_name}: {len(variant_seq)} bp, transposase inserted at naked-backbone position {pos}")


if __name__ == '__main__':
    main()
