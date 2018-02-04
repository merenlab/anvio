# -*- coding: utf-8 -*-
# pylint: disable=line-too-long
# v.140713
"""A very lightweight FASTA I/O library"""

import io
import sys
import gzip
import numpy
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
            self.output_file_obj = gzip.open(output_file_path, 'w')
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
            self.file_pointer = gzip.open(self.fasta_file_path)
        else:
            self.file_pointer = io.open(self.fasta_file_path, 'rU', newline='')

        if not self.file_pointer.read(1) == '>':
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

    def visualize_sequence_length_distribution(self, title, dest=None, max_seq_len=None, xtickstep=None, ytickstep=None):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        sequence_lengths = []

        self.reset()

        while next(self):
            if self.pos % 10000 == 0 or self.pos == 1:
                sys.stderr.write('\r[fastalib] Reading: %s' % (self.pos))
                sys.stderr.flush()
            sequence_lengths.append(len(self.seq))

        self.reset()

        sys.stderr.write('\n')

        if not max_seq_len:
            max_seq_len = max(sequence_lengths) + (int(max(sequence_lengths) / 100.0) or 10)

        seq_len_distribution = [0] * (max_seq_len + 1)

        for l in sequence_lengths:
            seq_len_distribution[l] += 1

        fig = plt.figure(figsize=(16, 12))
        plt.rcParams.update({'axes.linewidth': 0.9})
        plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)

        gs = gridspec.GridSpec(10, 1)

        ax1 = plt.subplot(gs[0:8])
        plt.grid(True)
        plt.subplots_adjust(left=0.05, bottom=0.03, top=0.95, right=0.98)

        plt.plot(seq_len_distribution, color='black', alpha=0.3)
        plt.fill_between(list(range(0, max_seq_len + 1)), seq_len_distribution, y2=0, color='black', alpha=0.15)
        plt.ylabel('number of sequences')
        plt.xlabel('sequence length')

        if xtickstep is None:
            xtickstep = (max_seq_len / 50) or 1

        if ytickstep is None:
            ytickstep = max(seq_len_distribution) / 20 or 1

        plt.xticks(list(range(xtickstep, max_seq_len + 1, xtickstep)), rotation=90, size='xx-small')
        plt.yticks(list(range(0, max(seq_len_distribution) + 1, ytickstep)),
                   [y for y in range(0, max(seq_len_distribution) + 1, ytickstep)],
                   size='xx-small')
        plt.xlim(xmin=0, xmax=max_seq_len)
        plt.ylim(ymin=0, ymax=max(seq_len_distribution) + (max(seq_len_distribution) / 20.0))

        plt.figtext(0.5, 0.96, '%s' % (title), weight='black', size='xx-large', ha='center')

        ax1 = plt.subplot(gs[9])
        plt.rcParams.update({'axes.edgecolor': 20})
        plt.grid(False)
        plt.yticks([])
        plt.xticks([])
        plt.text(0.02, 0.5, 'total: %s / mean: %.2f / std: %.2f / min: %s / max: %s'\
            % (len(sequence_lengths),
               numpy.mean(sequence_lengths), numpy.std(sequence_lengths),\
               min(sequence_lengths),\
               max(sequence_lengths)),\
            va='center', alpha=0.8, size='x-large')

        if dest is None:
            dest = self.fasta_file_path

        try:
            plt.savefig(dest + '.pdf')
        except:
            plt.savefig(dest + '.png')

        try:
            plt.show()
        except:
            pass

        return


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
    fasta = SequenceSource(sys.argv[1])
    fasta.visualize_sequence_length_distribution(title=sys.argv[2] if len(sys.argv) == 3 else 'None')
