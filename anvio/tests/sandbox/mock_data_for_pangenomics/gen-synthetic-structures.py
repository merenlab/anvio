#!/usr/bin/env python
"""Write synthetic protein structures for genes of a pangenome, for component tests.

Predicting structures needs MODELLER (and internet) or ColabFold (and a GPU), which a component test
can't count on. This script builds a protein-like backbone (N, CA, C, O) for the exact amino acid sequence
of a gene, so the structure files pass anvi'o's sequence integrity checks, and foldseek can compare them.

The geometry of a backbone is entirely decided by a 'fold seed': the seed draws a series of helices,
strands, and loops (as phi/psi dihedral angles), and the backbone is built from those angles. Two genes
with the same fold seed share the same geometry over their common length, so their structures are very
similar (foldseek reports a TM-score close to the length ratio of the two). Genes with different fold
seeds have unrelated structures.

Usage:

    python gen-synthetic-structures.py PAN_DB GENOMES_STORAGE SPEC_FILE OUTPUT_DIR [--all-members]

SPEC_FILE is a two-column TAB-delimited file without a header: gene cluster id, fold seed (an integer).
By default, one structure is written per gene cluster, for its representative. With --all-members, one
structure is written for each gene of each gene cluster. Next to the PDB files, OUTPUT_DIR receives

    external-structures.txt   (genome_name, gene_callers_id, path) for anvi-gen-structure-database
    structures-txt.txt        (gene_id, path) for anvi-pan-genome-structure (representatives only)
"""

import os
import sys
import math
import random
import argparse

import anvio.terminal as terminal

from anvio.dbops import PanSuperclass


THREE_LETTER = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
                'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
                'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

# ideal backbone geometry (Engh & Huber)
N_CA, CA_C, C_N, C_O = 1.458, 1.525, 1.329, 1.231
ANGLE_N_CA_C, ANGLE_CA_C_N, ANGLE_C_N_CA, ANGLE_CA_C_O = 111.2, 116.2, 121.7, 120.5
OMEGA = 180.0

HELIX, STRAND = (-57.0, -47.0), (-120.0, 130.0)
LOOP_CHOICES = [(-80.0, 150.0), (-70.0, -20.0), (60.0, 40.0), (-100.0, 10.0), (-150.0, 160.0)]


def dihedrals_for_fold(fold_seed, length):
    """Return `length` (phi, psi) pairs drawn from a series of secondary structure elements."""

    r = random.Random(fold_seed)
    angles = []
    while len(angles) < length:
        element = r.choice(['helix', 'helix', 'strand', 'loop'])
        if element == 'helix':
            angles.extend([HELIX] * r.randint(10, 18))
        elif element == 'strand':
            angles.extend([STRAND] * r.randint(5, 9))
        else:
            angles.extend([r.choice(LOOP_CHOICES) for _ in range(r.randint(3, 5))])

    return angles[:length]


def place(a, b, c, bond, angle, torsion):
    """Natural extension reference frame: the position of d, given a, b, c and the b-c-d internal coordinates."""

    angle, torsion = math.radians(angle), math.radians(torsion)

    bc = [c[i] - b[i] for i in range(3)]
    bc_len = math.sqrt(sum(x * x for x in bc))
    bc = [x / bc_len for x in bc]

    ab = [b[i] - a[i] for i in range(3)]
    n = [ab[1] * bc[2] - ab[2] * bc[1], ab[2] * bc[0] - ab[0] * bc[2], ab[0] * bc[1] - ab[1] * bc[0]]
    n_len = math.sqrt(sum(x * x for x in n))
    n = [x / n_len for x in n]

    m = [n[1] * bc[2] - n[2] * bc[1], n[2] * bc[0] - n[0] * bc[2], n[0] * bc[1] - n[1] * bc[0]]

    d2 = [-bond * math.cos(angle), bond * math.sin(angle) * math.cos(torsion), bond * math.sin(angle) * math.sin(torsion)]

    return [c[i] + d2[0] * bc[i] + d2[1] * m[i] + d2[2] * n[i] for i in range(3)]


def build_backbone(sequence, fold_seed):
    """Return a list of (residue_name, {atom_name: xyz}) for the sequence, following its fold seed."""

    dihedrals = dihedrals_for_fold(fold_seed, len(sequence))

    n = [0.0, 0.0, 0.0]
    ca = [N_CA, 0.0, 0.0]
    c = [ca[0] + CA_C * math.cos(math.radians(180 - ANGLE_N_CA_C)), CA_C * math.sin(math.radians(180 - ANGLE_N_CA_C)), 0.0]

    residues = []
    for i, aa in enumerate(sequence):
        phi, psi = dihedrals[i]

        if i > 0:
            prev_n, prev_ca, prev_c = residues[-1][1]['N'], residues[-1][1]['CA'], residues[-1][1]['C']
            prev_psi = dihedrals[i - 1][1]
            n = place(prev_n, prev_ca, prev_c, C_N, ANGLE_CA_C_N, prev_psi)
            ca = place(prev_ca, prev_c, n, N_CA, ANGLE_C_N_CA, OMEGA)
            c = place(prev_c, n, ca, CA_C, ANGLE_N_CA_C, phi)

        # the carbonyl oxygen sits opposite to the next nitrogen
        o = place(n, ca, c, C_O, ANGLE_CA_C_O, psi + 180.0)

        residues.append((THREE_LETTER[aa], {'N': n, 'CA': ca, 'C': c, 'O': o}))

    return residues


def write_pdb(sequence, fold_seed, output_path):
    lines = []
    serial = 1
    for residue_number, (residue_name, atoms) in enumerate(build_backbone(sequence, fold_seed), start=1):
        for atom_name in ['N', 'CA', 'C', 'O']:
            x, y, z = atoms[atom_name]
            lines.append("ATOM  %5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s" %
                         (serial, atom_name, residue_name, residue_number, x, y, z, 1.0, 0.0, atom_name[0]))
            serial += 1

    lines.append("TER")
    lines.append("END")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('pan_db')
    parser.add_argument('genomes_storage')
    parser.add_argument('spec_file')
    parser.add_argument('output_dir')
    parser.add_argument('--all-members', action='store_true')
    args = parser.parse_args()

    fold_seeds = {}
    for line in open(args.spec_file):
        if line.strip():
            gene_cluster_id, fold_seed = line.strip().split('\t')
            fold_seeds[gene_cluster_id] = int(fold_seed)

    pan_args = argparse.Namespace(pan_db=args.pan_db, genomes_storage=args.genomes_storage)
    pan = PanSuperclass(pan_args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
    pan.init_gene_clusters()

    if args.all_members:
        genes = [(gc, genome_name, gene_callers_id) for gc in fold_seeds
                                                    for genome_name in pan.gene_clusters[gc]
                                                    for gene_callers_id in pan.gene_clusters[gc][genome_name]]
    else:
        representatives = pan.get_gene_cluster_representative_sequences(gene_cluster_names=set(fold_seeds.keys()))
        genes = [(gc, representatives[gc]['genome_name'], representatives[gc]['gene_callers_id']) for gc in fold_seeds]

    os.makedirs(args.output_dir, exist_ok=True)
    output_dir = os.path.abspath(args.output_dir)

    external_structures = ["genome_name\tgene_callers_id\tpath"]
    structures_txt = ["gene_id\tpath"]
    for gc, genome_name, gene_callers_id in sorted(genes):
        sequence = pan.genomes_storage.get_gene_sequence(genome_name, gene_callers_id, report_DNA_sequences=False)

        unknown = set(sequence) - set(THREE_LETTER.keys())
        if unknown:
            sys.exit(f"Gene {gene_callers_id} of {genome_name} ({gc}) has non-standard amino acids ({unknown}). Pick another gene cluster.")

        path = os.path.join(output_dir, f"{genome_name}_{gene_callers_id}.pdb")
        write_pdb(sequence, fold_seeds[gc], path)

        external_structures.append(f"{genome_name}\t{gene_callers_id}\t{path}")
        structures_txt.append(f"{gc}\t{path}")

    with open(os.path.join(output_dir, 'external-structures.txt'), 'w') as f:
        f.write('\n'.join(external_structures) + '\n')

    if not args.all_members:
        with open(os.path.join(output_dir, 'structures-txt.txt'), 'w') as f:
            f.write('\n'.join(structures_txt) + '\n')

    print(f"{len(genes)} synthetic structures written to {output_dir}")


if __name__ == '__main__':
    main()
