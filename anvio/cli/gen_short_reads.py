#!/usr/bin/env python
# -*- coding: utf-8

import sys
import random
import configparser

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['configuration-ini']
__provides__ = ['short-reads-fasta']
__description__ = "Generate short reads from contigs. Useful to reconstruct mock data sets from already assembled contigs"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


class Configuration:
    def __init__(self, config):
        self.short_read_length = int(config.get('general', 'short_read_length'))
        self.coverage = int(config.get('general', 'coverage'))
        self.contig = str(config.get('general', 'contig'))

        self.SNVs = []

        for section in [s for s in config.sections() if s != 'general']:
            location = int(section)
            ratio = float(config.get(section, 'ratio'))

            if ratio > 1:
                raise ConfigError('Ratio cannot be more than 1 (error for SNV location %d)' % location)
            if ratio <= 0:
                raise ConfigError('Ratio cannot be 0 or less (error for SNV location %d)' % location)

            if location >= len(self.contig):
                raise ConfigError('SNV position at %d for a contig that is %d nts long? Really?' % (location, len(self.contig)))

            self.SNVs.append((location, ratio),)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()
    pp = terminal.pretty_print

    bases = ['A', 'T', 'C', 'G']

    if not args.output_file_path:
        raise ConfigError("You must define an output file path.")

    sequences = {}

    user_config = configparser.ConfigParser()
    user_config.read(args.configuration)
    config = Configuration(user_config)
    x = config.short_read_length
    c = config.coverage

    progress.new('Generating short reads')

    L = len(config.contig)

    av_num_short_reads_needed = int(L / x * c)

    for i in range(0, av_num_short_reads_needed):
        if (i + 1) % 100 == 0:
            progress.update('Entry %s of %s ...' % (pp(i), pp(av_num_short_reads_needed)))

        start_pos = random.randint(0, L - x)
        short_read = config.contig[start_pos:start_pos + x]

        sequences[i] = {'sequence': short_read, 'start': start_pos, 'stop': start_pos + x, 'num_SNVs': 0}

    progress.end()

    progress.new('Introducing SNVs')

    for snv_location, ratio in config.SNVs:
        progress.update('Working on location %d with ratio of %.2f' % (snv_location, ratio))

        matching_entries = []

        for entry_id in sequences:
            e = sequences[entry_id]

            if snv_location >= e['start'] and snv_location < e['stop']:
                matching_entries.append(entry_id)

        entries_to_mutate = random.sample(matching_entries, int(round(ratio * len(matching_entries))))
        current_base = config.contig[snv_location].upper()
        new_base = bases[(bases.index(current_base) + 1) % 4]

        for entry_id in entries_to_mutate:
            position_in_sequence_to_replace = snv_location - sequences[entry_id]['start']
            sequences[entry_id]['sequence'] = sequences[entry_id]['sequence'][:position_in_sequence_to_replace] + new_base + sequences[entry_id]['sequence'][position_in_sequence_to_replace + 1:]
            sequences[entry_id]['num_SNVs'] += 1

    with open(args.output_file_path, 'w') as output:
        for entry_id in sequences:
            s = sequences[entry_id]
            output.write('>%s\n' % '|'.join(['%d' % entry_id, 'start:%d' % s['start'], 'stop:%d' % s['stop'], 'num_SNVs:%d' % s['num_SNVs']]))
            output.write('%s\n' % s['sequence'])

    progress.end()

    run.info('Fasta output', args.output_file_path)


def get_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument('configuration', metavar = 'CONFIG_FILE',
                                        help = 'Configuration file')
    parser.add_argument('--output-file-path', metavar = 'FASTA_FILE',
                                        help = 'Output FASTA file path')



    return parser.get_args(parser)


if __name__ == '__main__':
    main()
