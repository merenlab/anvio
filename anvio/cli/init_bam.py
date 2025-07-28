#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import pysam

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.errors import FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__requires__ = ['raw-bam-file',]
__provides__ = ['bam-file',]
__description__ = "Sort/Index BAM files"
__resources__ = [("Another description as part of the metagenomic workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile")]


@terminal.time_program
def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    """sort and index a BAM file"""
    args = get_args()

    input_file_path = os.path.abspath(args.input_file)
    output_file_path = args.output_file
    num_threads = args.num_threads or 1

    progress = terminal.Progress()
    run = terminal.Run()

    if args.num_threads < 1:
        raise ConfigError("--num-threads should be at least 1.")

    if output_file_path:
        output_file_path = os.path.abspath(output_file_path)
        if not output_file_path.endswith('.bam'):
            raise ConfigError("The output file should end with extension '.bam'.")
    else:
        output_file_path = input_file_path + '-sorted.bam'

    filesnpaths.is_file_exists(input_file_path)
    filesnpaths.is_output_file_writable(output_file_path)

    progress.new('SORT')
    progress.update('Sorting BAM File... May take a while depending on the size.')

    # -@ is number of _additional_ threads (i.e the default is 0), so we subtract 1
    pysam.sort('-o', output_file_path, '-@', str(num_threads - 1), input_file_path)
    progress.end()

    run.info('Sorted BAM File', output_file_path)
    if not os.path.exists(output_file_path):
        raise ConfigError("Sorry. Something went wrong. Samtools thinks it generated the sorted output, yet it is not there :(")

    progress.new('INDEX')
    progress.update('Indexing BAM File...')
    pysam.index(output_file_path)
    progress.end()
    run.info('BAM File Index', output_file_path + '.bai', mc='green')


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument('input_file', metavar = 'BAM_FILE', help = 'BAM file to analyze')
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
