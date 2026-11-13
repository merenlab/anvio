#!/usr/bin/env python3
"""Mock data for the component test of deletions that cross split boundaries.

    mock.py reference          writes ref.fa (one 60 kb contig, fixed seed)
    mock.py reads CONTIGS.db   reads the contig's first split boundary B from the contigs database, and writes reads.fa:
                               reads of a genome that lacks [B - 600, B + 900), a 1,500 bp deletion across the boundary
    mock.py check indels.txt CONTIGS.db
                               the profile's indels table must hold the deletion once, whole, at B - 600, and nothing
                               starting at B (exit 1 otherwise)
    mock.py collections CONTIGS.db
                               writes two collections of the contig's splits: cut.txt puts the splits before and after B
                               in two bins (BEFORE, AFTER), and kept.txt puts them in one bin (BOTH)
    mock.py check-split indels.txt CONTIGS.db BEFORE.txt AFTER.txt BOTH.txt
                               the indels tables of the split profiles: in BEFORE, the deletion ends at B; in AFTER, it
                               starts at the first position of the split, with the reference base of that position; in
                               BOTH, it is whole, as in the profile (exit 1 otherwise)
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


def splits(db):
    con = sqlite3.connect(db)
    return [s for (s, _) in sorted(con.execute("SELECT split, start FROM splits_basic_info WHERE parent = 'c_1'"), key=lambda x: x[1])]


def reference():
    rng = random.Random(7)
    return ''.join(rng.choice('ACGT') for _ in range(L))


def long_deletions(path):
    with open(path) as fh:
        return [r for r in csv.DictReader(fh, delimiter='\t') if r['type'] == 'DEL' and int(r['length']) >= 50]


def check_split(indels, db, before, after, both):
    b, split_names = boundary(db), splits(db)

    whole = long_deletions(indels)
    if len(whole) != 1:
        print(f"FAIL: expected one deletion in the profile, found {len(whole)}")
        return False

    whole = whole[0]
    start, end, count = int(whole['pos_in_contig']), int(whole['pos_in_contig']) + int(whole['length']), whole['count']

    expected = {
        'BEFORE': (before, {'split_name': split_names[0], 'pos_in_contig': str(start), 'length': str(b - start), 'count': count}),
        'AFTER': (after, {'split_name': split_names[1], 'pos': '0', 'pos_in_contig': str(b), 'length': str(end - b), 'count': count, 'reference': reference()[b]}),
        'BOTH': (both, {key: whole[key] for key in ['split_name', 'pos', 'pos_in_contig', 'length', 'count', 'coverage', 'reference']}),
    }

    passed = True
    for bin_name, (path, values) in expected.items():
        dels = long_deletions(path)
        found = [{key: r[key] for key in values} for r in dels]
        print(f"{bin_name}: expected {values}, found {found}")
        if found != [values]:
            passed = False

    return passed


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
    elif sys.argv[1] == 'collections':
        split_names = splits(sys.argv[2])
        with open('cut.txt', 'w') as cut, open('kept.txt', 'w') as kept:
            for i, split_name in enumerate(split_names):
                cut.write(f"{split_name}\t{['BEFORE', 'AFTER'][i] if i < 2 else 'REST'}\n")
                kept.write(f"{split_name}\t{'BOTH' if i < 2 else 'REST'}\n")
    elif sys.argv[1] == 'check-split':
        if check_split(*sys.argv[2:7]):
            print("PASS: anvi-split cuts the deletion where the bin cuts the contig, and keeps it whole otherwise")
            sys.exit(0)
        print("FAIL: anvi-split does not cut the deletion where the bin cuts the contig, or does not keep it whole otherwise")
        sys.exit(1)
