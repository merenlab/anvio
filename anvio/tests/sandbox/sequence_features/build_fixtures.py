#!/usr/bin/env python
"""Builds the small FASTA + GenBank fixtures used by the sequence-features
integration tests. Re-run this to regenerate fixtures if Biopython's emitted
format ever changes; the fixtures themselves are committed to the repo.

Coordinates in this file are zero-based half-open (Biopython native). The
written GenBank file will display the one-based inclusive form.
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation


HERE = os.path.dirname(os.path.abspath(__file__))


def _record(name, seq, description, features, topology='linear'):
    rec = SeqRecord(Seq(seq), id=name, name=name, description=description)
    rec.annotations = {'molecule_type': 'DNA', 'topology': topology, 'data_file_division': 'BCT'}
    rec.features = features
    return rec


def _write_pair(records, fa_path, gb_path):
    # write the FASTA without per-feature decoration
    fasta_records = [SeqRecord(r.seq, id=r.id, description='') for r in records]
    with open(fa_path, 'w') as out:
        SeqIO.write(fasta_records, out, 'fasta')
    with open(gb_path, 'w') as out:
        SeqIO.write(records, out, 'genbank')


def build_bacterial():
    """Two linear contigs, three gene/CDS pairs (one forward, one reverse, one with
    extra qualifiers). Used by tests 3, 7, 9, 10."""

    seq1 = 'ACGT' * 500   # 2000 bp
    seq2 = 'TGCA' * 500   # 2000 bp

    rec1 = _record('bact_c01', seq1, 'Bacterial test contig 1', [
        # forward gene / CDS pair with rich qualifiers
        SeqFeature(SimpleLocation(99, 500, strand=1),  type='gene', qualifiers={'locus_tag': ['TST_001'], 'gene': ['testA']}),
        SeqFeature(SimpleLocation(99, 500, strand=1),  type='CDS',  qualifiers={
            'locus_tag': ['TST_001'], 'gene': ['testA'],
            'codon_start': ['1'], 'transl_table': ['11'],
            'product': ['Test protein A'],
            'protein_id': ['TST_PROT_A'],
            'translation': ['MTESTPROTA'],
            'db_xref': ['REF:first', 'REF:second'],
            'note': ['this is a note with whitespace'],
        }),
        # reverse-strand gene / CDS pair
        SeqFeature(SimpleLocation(599, 900, strand=-1), type='gene', qualifiers={'locus_tag': ['TST_002'], 'gene': ['testB']}),
        SeqFeature(SimpleLocation(599, 900, strand=-1), type='CDS',  qualifiers={
            'locus_tag': ['TST_002'], 'gene': ['testB'],
            'codon_start': ['1'], 'transl_table': ['11'],
            'product': ['Test protein B'],
            'translation': ['MTESTPROTB'],
        }),
        # gene without an mRNA — CDS will link to it via derives_from
        SeqFeature(SimpleLocation(1199, 1500, strand=1), type='gene', qualifiers={'locus_tag': ['TST_003'], 'gene': ['testC']}),
        SeqFeature(SimpleLocation(1199, 1500, strand=1), type='CDS',  qualifiers={
            'locus_tag': ['TST_003'], 'gene': ['testC'],
            'codon_start': ['1'], 'transl_table': ['11'],
            'product': ['Test protein C'],
            'translation': ['MTESTPROTC'],
        }),
    ])

    rec2 = _record('bact_c02', seq2, 'Bacterial test contig 2', [
        SeqFeature(SimpleLocation(99, 700, strand=1), type='gene', qualifiers={'locus_tag': ['TST_004'], 'gene': ['testD']}),
        SeqFeature(SimpleLocation(99, 700, strand=1), type='CDS',  qualifiers={
            'locus_tag': ['TST_004'], 'gene': ['testD'],
            'codon_start': ['1'], 'transl_table': ['11'],
            'product': ['Test protein D'],
            'translation': ['MTESTPROTD'],
        }),
    ])

    _write_pair([rec1, rec2], os.path.join(HERE, 'bacterial.fa'), os.path.join(HERE, 'bacterial.gb'))


def build_eukaryotic():
    """One linear contig containing two genes with multi-exon, multi-segment CDSs:
    one forward, one reverse, both with explicit mRNA records so parent resolution
    is fully exercised. Used by test 4."""

    seq = 'ACGT' * 1500   # 6000 bp

    features = [
        # ---------- forward-strand gene with 3 exons / 3-segment CDS ----------
        SeqFeature(SimpleLocation(99, 2500, strand=1), type='gene',
                   qualifiers={'locus_tag': ['EUK_FWD_001'], 'gene': ['fwdGene']}),
        SeqFeature(CompoundLocation([SimpleLocation(99, 500, strand=1), SimpleLocation(900, 1500, strand=1), SimpleLocation(2000, 2500, strand=1)]),
                   type='mRNA',
                   qualifiers={'locus_tag': ['EUK_FWD_001'], 'gene': ['fwdGene']}),
        SeqFeature(SimpleLocation(99, 500, strand=1),  type='exon',
                   qualifiers={'locus_tag': ['EUK_FWD_001']}),
        SeqFeature(SimpleLocation(900, 1500, strand=1), type='exon',
                   qualifiers={'locus_tag': ['EUK_FWD_001']}),
        SeqFeature(SimpleLocation(2000, 2500, strand=1), type='exon',
                   qualifiers={'locus_tag': ['EUK_FWD_001']}),
        SeqFeature(CompoundLocation([SimpleLocation(99, 500, strand=1), SimpleLocation(900, 1500, strand=1), SimpleLocation(2000, 2500, strand=1)]),
                   type='CDS',
                   qualifiers={
                       'locus_tag': ['EUK_FWD_001'], 'gene': ['fwdGene'],
                       'codon_start': ['1'], 'transl_table': ['1'],
                       'product': ['Forward multi-exon protein'],
                       'translation': ['MFORWARDPROT'],
                   }),

        # ---------- reverse-strand gene with 2 exons / 2-segment CDS ----------
        # transcription order should put the GENOMIC RIGHTMOST segment first.
        SeqFeature(SimpleLocation(2999, 4500, strand=-1), type='gene',
                   qualifiers={'locus_tag': ['EUK_REV_002'], 'gene': ['revGene']}),
        SeqFeature(CompoundLocation([SimpleLocation(3999, 4500, strand=-1), SimpleLocation(2999, 3500, strand=-1)]),
                   type='mRNA',
                   qualifiers={'locus_tag': ['EUK_REV_002'], 'gene': ['revGene']}),
        SeqFeature(SimpleLocation(3999, 4500, strand=-1), type='exon',
                   qualifiers={'locus_tag': ['EUK_REV_002']}),
        SeqFeature(SimpleLocation(2999, 3500, strand=-1), type='exon',
                   qualifiers={'locus_tag': ['EUK_REV_002']}),
        SeqFeature(CompoundLocation([SimpleLocation(3999, 4500, strand=-1), SimpleLocation(2999, 3500, strand=-1)]),
                   type='CDS',
                   qualifiers={
                       'locus_tag': ['EUK_REV_002'], 'gene': ['revGene'],
                       'codon_start': ['1'], 'transl_table': ['1'],
                       'product': ['Reverse multi-exon protein'],
                       'translation': ['MREVERSEPROT'],
                   }),
    ]

    rec = _record('euk_c01', seq, 'Eukaryotic multi-exon test contig', features)
    _write_pair([rec], os.path.join(HERE, 'eukaryotic.fa'), os.path.join(HERE, 'eukaryotic.gb'))


def build_circular_origin_crossing():
    """A 1000 bp circular contig with one forward feature that crosses the origin
    via `join(900..1000, 1..200)`. Used by test 5."""

    seq = 'ACGT' * 250  # 1000 bp

    # zero-based half-open: [899, 1000) and [0, 200) — same as join(900..1000, 1..200)
    feat = SeqFeature(
        CompoundLocation([SimpleLocation(899, 1000, strand=1), SimpleLocation(0, 200, strand=1)]),
        type='gene',
        qualifiers={'locus_tag': ['ORIGIN_CROSSER_001'], 'gene': ['oriCrosser']},
    )

    rec = _record('circ_c01', seq, 'Circular contig with an origin-crossing feature', [feat], topology='circular')
    _write_pair([rec], os.path.join(HERE, 'circular.fa'), os.path.join(HERE, 'circular.gb'))


def build_malformed_linear_origin():
    """A LINEAR contig where a single feature uses a join(...) that crosses the
    origin. We expect the importer to warn and skip that feature, while still
    importing other features in the same record. Used by test 6."""

    seq = 'ACGT' * 250  # 1000 bp

    # malformed: join across a linear contig's origin
    bad = SeqFeature(
        CompoundLocation([SimpleLocation(899, 1000, strand=1), SimpleLocation(0, 200, strand=1)]),
        type='gene',
        qualifiers={'locus_tag': ['LINEAR_BAD_001'], 'gene': ['linBad']},
    )
    # one healthy feature to confirm the rest of the record still imports
    healthy = SeqFeature(SimpleLocation(399, 700, strand=1), type='gene', qualifiers={'locus_tag': ['LINEAR_OK_001']})

    rec = _record('lin_bad_c01', seq, 'Linear contig with a malformed origin-crosser', [bad, healthy], topology='linear')
    _write_pair([rec], os.path.join(HERE, 'linear_malformed.fa'), os.path.join(HERE, 'linear_malformed.gb'))


def build_unknown_type():
    """A linear contig with one builtin (gene) and one non-builtin (regulatory)
    feature. Used by test 8."""

    seq = 'ACGT' * 250
    features = [
        SeqFeature(SimpleLocation(99, 300, strand=1), type='gene',
                   qualifiers={'locus_tag': ['UNK_001'], 'gene': ['someGene']}),
        SeqFeature(SimpleLocation(400, 600, strand=1), type='regulatory',
                   qualifiers={'regulatory_class': ['promoter'], 'note': ['putative promoter region']}),
    ]
    rec = _record('unk_c01', seq, 'Contig with a non-builtin regulatory feature', features)
    _write_pair([rec], os.path.join(HERE, 'unknown_type.fa'), os.path.join(HERE, 'unknown_type.gb'))


if __name__ == '__main__':
    build_bacterial()
    build_eukaryotic()
    build_circular_origin_crossing()
    build_malformed_linear_origin()
    build_unknown_type()
    print("All fixtures rebuilt.")
