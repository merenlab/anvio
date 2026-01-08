#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import anvio
import numpy as np

import anvio.bamops as bamops
import anvio.terminal as terminal
import multiprocess as multiprocessing

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet']
__provides__ = []
__requires__ = ["bam-file"]
__description__ = ("This script report errors in long read assembly using read-recruitment information. "
                   "The input file should be a BAM file of long reads mapped to an assembly made from "
                   "these reads.")


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def process_contig(args, available_index_queue, output_queue, contigs_size):
    bam = bamops.BAMFileObject(args.bam_file, 'rb')

    while True:
        try:
            idx = available_index_queue.get(True)  # blocking call
            contig = list(contigs_size.keys())[idx]
            length = contigs_size[contig]

            coverage = np.zeros(length, dtype=np.uint32)
            clipping = {}

            # got through each read, compute coverage and idc (insertion, deletion, (hard - soft)clipping).
            for read in bam.fetch(contig):
                if read.flag == 4:
                    continue

                current_pos = read.reference_start
                cigar_tup = read.cigartuples
                # counting the numb of tup because position index is different for
                # hard/soft clipping if they are at the begining or at the end of cigar string
                num_tup = 0
                for tup in cigar_tup:
                    num_tup += 1
                    # if mapping, compute cov, increase current position
                    if tup[0] == 0:
                        coverage[current_pos:current_pos + tup[1]] += 1
                        current_pos += tup[1]
                    # if deletion, increase current position
                    elif tup[0] == 2:
                        current_pos += tup[1]
                    # if clipping (soft-hard: tup[0] = 4 or 5), then +1 clipping
                    elif tup[0] in [4, 5]:
                        if num_tup == 1:
                            if current_pos != 0:
                                clipping[current_pos] = clipping.get(current_pos, 0) + 1
                        elif current_pos != length:
                            clipping[current_pos - 1] = clipping.get(current_pos - 1, 0) + 1

            # bundle clipping positions with their coverage values
            clipping_with_cov = {pos: (clip_count, coverage[pos]) for pos, clip_count in clipping.items()}

            # compute zero-coverage ranges
            zero_ranges = []
            in_window = False
            window_start = 0
            for pos in range(length):
                if coverage[pos] == 0 and not in_window:
                    window_start = pos
                    in_window = True
                elif coverage[pos] > 0 and in_window:
                    zero_ranges.append((window_start, pos))
                    in_window = False
            # handle window that extends to end of contig
            if in_window:
                zero_ranges.append((window_start, length))

            output_queue.put((contig, length, clipping_with_cov, zero_ranges))

        except Exception as e:
            output_queue.put(e)


def run_program():
    args = get_args()

    if not args.just_do_it:
        raise ConfigError("This script ONLY makes sense if you are using a BAM file that was made from "
                          "mapping long read onto an assembly made with the SAME long reads. "
                          "If you are positive that you did JUST that, then re-run this program with "
                          "`--just-do-it` flag.")

    run = terminal.Run()
    progress = terminal.Progress()

    min_dist_to_end = args.min_dist_to_end
    min_clipping_ratio = args.clipping_ratio

    run.info('BAM file', args.bam_file)
    run.info('Length of contig\'s end to ignore', args.min_dist_to_end)
    run.info('Number of threads', args.num_threads)

    manager = multiprocessing.Manager()
    available_index_queue = manager.Queue()
    output_queue = manager.Queue()

    # make a dict of contigs and their size
    contigs_size = {}
    bam = bamops.BAMFileObject(args.bam_file, 'rb')
    for i, contig in enumerate(bam.references):
        contigs_size[contig] = bam.lengths[i]
        available_index_queue.put(i)
    bam.close()

    num_contigs = len(contigs_size)
    num_threads = min(num_contigs, args.num_threads)
    if num_threads != args.num_threads:
        run.info_single(f"You have {terminal.pluralize('contig', num_contigs)} which is less "
                        f"than the {args.num_threads} threads that you requested. Anvi'o will only use "
                        f"{terminal.pluralize('thread', num_contigs)}")

    progress.new(f"Computing coverage with {terminal.pluralize('thread', num_threads)} beep boop", progress_total_items = num_contigs)
    progress.update('initializing threads ...')

    processes = []
    for i in range(num_threads):
        p = multiprocessing.Process(
            target=process_contig,
            args=(args, available_index_queue, output_queue, contigs_size)
        )
        processes.append(p)
        p.start()

    # open output files for streaming writes
    clipping_ouput = args.output_prefix + '-clipping.txt'
    zero_ouput = args.output_prefix + "-zero_cov.txt"
    run.info('Output file', clipping_ouput)
    run.info('Output file', zero_ouput)

    received = 0

    with open(clipping_ouput, 'w') as clipping_file, open(zero_ouput, 'w') as zero_file:
        clipping_file.write("contig\tlength\tpos\trelative_pos\tcov\tclipping\tclipping_ratio\n")
        zero_file.write("contig\tlength\trange\trange_size\n")

        while received < num_contigs:
            try:
                result = output_queue.get()
                if isinstance(result, Exception):
                    for p in processes:
                        p.terminate()
                    raise result
                contig, contig_length, clipping_with_cov, zero_ranges = result
                received += 1
                progress.update(f"computing contigs {received}/{num_contigs}")
                progress.increment(increment_to = received)

                # write clipping results for this contig
                for pos, (clip_count, cov_at_pos) in clipping_with_cov.items():
                    clipping_ratio = clip_count/cov_at_pos
                    relative_pos = pos/contig_length
                    if clipping_ratio > min_clipping_ratio and pos > min_dist_to_end and contig_length-pos > min_dist_to_end:
                        clipping_file.write(f"{contig}\t{contig_length}\t{pos}\t{relative_pos}\t{cov_at_pos}\t{clip_count}\t{clipping_ratio}\n")

                # write zero coverage results for this contig
                for window_start, window_end in zero_ranges:
                    window_length = window_end - window_start
                    zero_file.write(f"{contig}\t{contig_length}\t{window_start}-{window_end}\t{window_length}\n")

            except KeyboardInterrupt:
                run.info_single("Received SIGINT, terminating processes...")
                for p in processes:
                    p.terminate()
                break

    for p in processes:
        p.terminate()

    progress.end()


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('REQUIRED', 'Declare your BAM file here')
    groupA.add_argument('-b', '--bam-file', required = True, help = "Sorted and indexed BAM file to analyze.")

    groupB = parser.add_argument_group('OUTPUT-PREFIX', "Choose a prefix for the two output files")
    groupB.add_argument( '-o', '--output-prefix', required=True, help = ("Output prefix for the two tab-"
            "delimited files. Will overwrite existing files."))

    groupC = parser.add_argument_group('OTHER PARAMETERS', "Various paramaters that you can change")
    groupC.add_argument('-m', '--min-dist-to-end', default=100, type=int, help = ("Sequence length on "
            "each side of a contigs that will be ignored in the clipping outut. "
            "This is to prevent aftifactual results that typically occurs at the ends of contigs."))
    groupC.add_argument( '-r', '--clipping-ratio', default=0.8, type=float,  help = ("Minimum ratio of "
            "'coverage-of-clipped-reads / total-coverage' that will be reported "
            "in the output file. Default is 0.8"))
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupC.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
