# -*- coding: utf-8 -*-
# pylint: disable=line-too-long
# v.140713
"""A very lightweight FASTA I/O library"""

import io
import sys
import gzip
import hashlib

import anvio

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class FastaOutput:
    def __init__(self, output_file_path):
        self.output_file_path = output_file_path
        self.compressed = True if self.output_file_path.endswith('.gz') else False

        if self.compressed:
            self.output_file_obj = gzip.open(output_file_path, 'wt')
        else:
            self.output_file_obj = open(output_file_path, 'w')

    def store(self, entry, split=True, store_frequencies=True):
        if entry.unique and store_frequencies:
            self.write_id('%s|%s' % (entry.id, 'frequency:%d' % len(entry.ids)))
        else:
            self.write_id(entry.id)

        self.write_seq(entry.seq, split)

    def write_id(self, id):
        self.output_file_obj.write('>%s\n' % id)

    def write_seq(self, seq, split=True):
        if split:
            seq = self.split(seq)
        self.output_file_obj.write('%s\n' % seq)

    def split(self, sequence, piece_length=80):
        ticks = list(range(0, len(sequence), piece_length)) + [len(sequence)]
        return '\n'.join([sequence[ticks[x]:ticks[x + 1]] for x in range(0, len(ticks) - 1)])

    def close(self):
        self.output_file_obj.close()


class ReadFasta:
    def __init__(self, f_name, quiet=False):
        self.ids = []
        self.sequences = []

        self.fasta = SequenceSource(f_name)

        while next(self.fasta):
            if (not quiet) and (self.fasta.pos % 1000 == 0 or self.fasta.pos == 1):
                sys.stderr.write('\r[fastalib] Reading FASTA into memory: %s' % (self.fasta.pos))
                sys.stderr.flush()
            self.ids.append(self.fasta.id)
            self.sequences.append(self.fasta.seq)

        if not quiet:
            sys.stderr.write('\n')

    def close(self):
        self.fasta.close()


class SequenceSource():
    def __init__(self, fasta_file_path, lazy_init=True, unique=False, allow_mixed_case=False):
        self.fasta_file_path = fasta_file_path
        self.name = None
        self.compressed = True if self.fasta_file_path.endswith('.gz') else False
        self.lazy_init = lazy_init
        self.allow_mixed_case = allow_mixed_case

        self.pos = 0
        self.id = None
        self.seq = None
        self.ids = []

        self.unique = unique
        self.unique_hash_dict = {}
        self.unique_hash_list = []
        self.unique_next_hash = 0

        if self.compressed:
            self.file_pointer = gzip.open(self.fasta_file_path, mode="rt")
        else:
            self.file_pointer = io.open(self.fasta_file_path, 'r', newline='')

        if not self.file_pointer.read(1) == '>':
            self.file_pointer.close()
            raise FastaLibError("File '%s' does not seem to be a FASTA file." % self.fasta_file_path)

        self.file_pointer.seek(0)

        if self.lazy_init:
            self.total_seq = None
        else:
            self.total_seq = len([l for l in self.file_pointer.readlines() if l.startswith('>')])
            self.reset()

        if self.unique:
            self.init_unique_hash()


    def init_unique_hash(self):
        while self.next_regular():
            hash = hashlib.sha1(self.seq.upper().encode('utf-8')).hexdigest()
            if hash in self.unique_hash_dict:
                self.unique_hash_dict[hash]['ids'].append(self.id)
                self.unique_hash_dict[hash]['count'] += 1
            else:
                self.unique_hash_dict[hash] = {'id': self.id,
                                               'ids': [self.id],
                                               'seq': self.seq,
                                               'count': 1}

        self.unique_hash_list = [i[1] for i in sorted([(self.unique_hash_dict[hash]['count'], hash)\
                        for hash in self.unique_hash_dict], reverse=True)]


        self.total_unique = len(self.unique_hash_dict)
        self.reset()


    def __next__(self):
        if self.unique:
            return self.next_unique()
        else:
            return self.next_regular()


    def next_unique(self):
        if self.unique:
            if self.total_unique > 0 and self.pos < self.total_unique:
                hash_entry = self.unique_hash_dict[self.unique_hash_list[self.pos]]

                self.pos += 1
                self.seq = hash_entry['seq'] if self.allow_mixed_case else hash_entry['seq'].upper()
                self.id = hash_entry['id']
                self.ids = hash_entry['ids']

                return True
            else:
                return False
        else:
            return False


    def next_regular(self):
        self.seq = None
        self.id = self.file_pointer.readline()[1:].strip()
        sequence = ''

        while True:
            line = self.file_pointer.readline()
            if not line:
                if len(sequence):
                    self.seq = sequence if self.allow_mixed_case else sequence.upper()
                    self.pos += 1
                    return True
                else:
                    return False
            if line.startswith('>'):
                self.file_pointer.seek(self.file_pointer.tell() - len(line))
                break
            sequence += line.strip()

        self.seq = sequence if self.allow_mixed_case else sequence.upper()
        self.pos += 1
        return True


    def get_seq_by_read_id(self, read_id):
        self.reset()
        while next(self):
            if self.id == read_id:
                return self.seq

        return False


    def close(self):
        self.file_pointer.close()


    def reset(self):
        self.pos = 0
        self.id = None
        self.seq = None
        self.ids = []
        self.file_pointer.seek(0)


class QualSource:
    def __init__(self, quals_file_path, lazy_init=True):
        self.quals_file_path = quals_file_path
        self.name = None
        self.lazy_init = lazy_init

        self.pos = 0
        self.id = None
        self.quals = None
        self.quals_int = None
        self.ids = []

        self.file_pointer = open(self.quals_file_path)
        self.file_pointer.seek(0)

        if self.lazy_init:
            self.total_quals = None
        else:
            self.total_quals = len([l for l in self.file_pointer.readlines() if l.startswith('>')])
            self.reset()


    def __next__(self):
        self.id = self.file_pointer.readline()[1:].strip()
        self.quals = None
        self.quals_int = None

        qualscores = ''

        while True:
            line = self.file_pointer.readline()
            if not line:
                if len(qualscores):
                    self.quals = qualscores.strip()
                    self.quals_int = [int(q) for q in self.quals.split()]
                    self.pos += 1
                    return True
                else:
                    return False
            if line.startswith('>'):
                self.file_pointer.seek(self.file_pointer.tell() - len(line))
                break
            qualscores += ' ' + line.strip()

        self.quals = qualscores.strip()
        self.quals_int = [int(q) for q in self.quals.split()]
        self.pos += 1

        return True

    def close(self):
        self.file_pointer.close()

    def reset(self):
        self.pos = 0
        self.id = None
        self.quals = None
        self.quals_int = None
        self.ids = []
        self.file_pointer.seek(0)


class FastaLibError(Exception):
    def __init__(self, e=None):
        Exception.__init__(self)
        while True:
            if e.find("  ") > -1:
                e = e.replace("  ", " ")
            else:
                break
        self.e = e
        return
    def __str__(self):
        return 'Fasta Lib Error: %s' % self.e


if __name__ == '__main__':
    pass
