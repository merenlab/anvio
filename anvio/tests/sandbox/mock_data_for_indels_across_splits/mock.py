#!/usr/bin/env python3
"""Mock data for the component test of deletions that cross split boundaries.

    mock.py reference          writes ref.fa (one 60 kb contig, fixed seed)
    mock.py reads CONTIGS.db   reads the contig's first split boundary B from the contigs database, and writes reads.fa:
                               reads of a genome that lacks [B - 600, B + 900), a 1,500 bp deletion across the boundary
    mock.py check indels.txt CONTIGS.db
                               the profile's indels table must hold the deletion once, whole, at B - 600, and nothing
                               starting at B (exit 1 otherwise)
"""
import csv
import random
import sqlite3
import sys

L, LEFT, RIGHT = 60000, 600, 900


def boundary(db):
    con = sqlite3.connect(db)
    starts = sorted(s for (s,) in con.execute("SELECT start FROM splits_basic_info WHERE parent = 'c_1' AND start > 0"))
    return starts[0]


def reference():
    rng = random.Random(7)
    return ''.join(rng.choice('ACGT') for _ in range(L))


if __name__ == '__main__':
    if sys.argv[1] == 'reference':
        with open('ref.fa', 'w') as fh:
            fh.write(f">c_1\n{reference()}\n")
    elif sys.argv[1] == 'reads':
        b, ref, rng = boundary(sys.argv[2]), reference(), random.Random(8)
        genome = ref[:b - LEFT] + ref[b + RIGHT:]
        with open('reads.fa', 'w') as fh:
            for i in range(400):
                n = rng.randint(3000, 8000)
                s = rng.randint(0, len(genome) - n)
                seq = ''.join(c if rng.random() > 0.002 else rng.choice('ACGT') for c in genome[s:s + n])
                fh.write(f">r{i}\n{seq}\n")
    elif sys.argv[1] == 'check':
        b = boundary(sys.argv[3])
        dels = [r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t') if r['type'] == 'DEL' and int(r['length']) >= 50]
        # the aligner may place the deletion a few bases away from where it was made (equally good alignments)
        whole = [r for r in dels if abs(int(r['pos_in_contig']) - (b - LEFT)) <= 10 and int(r['length']) == LEFT + RIGHT]
        pieces = [r for r in dels if abs(int(r['pos_in_contig']) - b) <= 2 or abs(int(r['pos_in_contig']) + int(r['length']) - b) <= 2]
        print(f"split boundary {b}; deletions >= 50 bp: {[(r['pos_in_contig'], r['length'], r['count']) for r in dels]}")
        if len(whole) == 1 and int(whole[0]['count']) > 0 and not pieces:
            print("PASS: the deletion across the split boundary is recorded once, whole")
            sys.exit(0)
        print("FAIL: the deletion across the split boundary is not recorded once, whole")
        sys.exit(1)
