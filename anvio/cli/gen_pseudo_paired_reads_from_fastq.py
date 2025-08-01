#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__requires__ = ["short-reads-fasta"]
__provides__ = ["paired-end-fastq"]
__description__ =  ("A script that takes a FASTQ file that is not paired-end (i.e., R1 alone) and converts "
                    "it into two FASTQ files that are paired-end (i.e., R1 and R2). This is a quick-and-dirty "
                    "workaround that halves each read from the original FASTQ and puts one half in "
                    "the FASTQ file for R1 and puts the reverse-complement of the second half in the "
                    "FASTQ file for R2. If you've ended up here, things have clearly not gone very well "
                    "for you, and Evan, who battled similar battles and ended up implementing this "
                    "solution wholeheartedly sympathizes")


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)



class FastqParse(object):
    def __init__(self, filepath):
        self.filepath = filepath
        filesnpaths.is_file_exists(self.filepath)

        self.set_num_reads()

        self.f = open(self.filepath, 'r')


    def set_num_reads(self):
        with open(self.filepath, 'r') as f:
            num_lines = sum(1 for line in f)

        if num_lines % 4 != 0:
            raise ConfigError("The number of lines in this 'FASTQ' file is not divisible by 4 "
                              "which makes anvi'o think this is not a FASTQ file.")

        self.num_reads = num_lines // 4


    def loop_reads(self):
        for i in range(self.num_reads):
            read = type("Read", (object,), {})()
            read.id1 = self.f.readline().strip()
            read.seq = self.f.readline().strip()
            read.id2 = self.f.readline().strip()
            read.qual = self.f.readline().strip()
            yield read


    def close(self):
        self.f.close()


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    r1_out_path = args.output_file_prefix + '_1.fastq'
    r2_out_path = args.output_file_prefix + '_2.fastq'

    filesnpaths.is_output_file_writable(r1_out_path, ok_if_exists=False)
    filesnpaths.is_output_file_writable(r2_out_path, ok_if_exists=False)

    r1_out = open(r1_out_path, 'w')
    r2_out = open(r2_out_path, 'w')

    fastq = FastqParse(args.fastq)

    progress.new('Processing reads', progress_total_items=fastq.num_reads)

    for i, read in enumerate(fastq.loop_reads()):
        if i % 100000 == 0:
            progress.increment(i)
            progress.update('%d / %d' % (i, fastq.num_reads))

        read1_len = len(read.seq)//2

        read1_seq = read.seq[:read1_len]
        read1_qual = read.qual[:read1_len]

        read2_seq = utils.rev_comp(read.seq[read1_len:])
        read2_qual = read.qual[read1_len:][::-1]

        r1_out.write('\n'.join([
            read.id1,
            read1_seq,
            read.id2,
            read1_qual,
        ]) + '\n')

        r2_out.write('\n'.join([
            read.id1,
            read2_seq,
            read.id2,
            read2_qual,
        ]) + '\n')

    r1_out.close()
    r2_out.close()

    progress.end()
    run.info_single("Done!")


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('-f', '--fastq', required=True)
    parser.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'help': \
                        'If you want final FASTQs with the format myfastq_1.fastq and myfastq_2.fastq, '
                        'then this parameter should be set to myfastq'}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
