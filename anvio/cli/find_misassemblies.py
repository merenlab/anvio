#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import anvio

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

            coverage = [0] * length
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
                        for pos in range(current_pos, current_pos + tup[1]):
                            coverage[pos] += 1
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

            output_queue.put((contig, length, coverage, clipping))

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

    cov_dict = {}
    received = 0

    while received < num_contigs:
        try:
            result = output_queue.get()
            if isinstance(result, Exception):
                for p in processes:
                    p.terminate()
                raise result
            contig, length, cov, clip = result
            cov_dict[contig] = {
                'length': length,
                'cov': dict(enumerate(cov)),
                'clipping': clip
            }
            received += 1
            progress.update(f"computing contigs {received}/{num_contigs}")
            progress.increment(increment_to = received)
        except KeyboardInterrupt:
            run.info_single("Received SIGINT, terminating processes...")
            for p in processes:
                p.terminate()
            break

    for p in processes:
        p.terminate()

    progress.end()

    # writting the outputs
    clipping_ouput = args.output_prefix + '-clipping.txt'
    run.info('Output file', clipping_ouput)
    with open(clipping_ouput, 'w') as file:
        file.write("contig\tlength\tpos\trelative_pos\tcov\tclipping\tclipping_ratio\n")
        for contig, data in cov_dict.items():
            contig_length = data['length']
            for pos in data['clipping']:
                cov = cov_dict[contig]['cov'][pos]
                clipping = cov_dict[contig]['clipping'][pos]
                clipping_ratio = clipping/cov
                relative_pos = pos/contig_length
                if clipping_ratio > min_clipping_ratio and pos > min_dist_to_end and contig_length-pos > min_dist_to_end:
                    file.write(f"{contig}\t{contig_length}\t{pos}\t{relative_pos}\t{cov}\t{clipping}\t{clipping_ratio}\n")

    zero_ouput = args.output_prefix + "-zero_cov.txt"
    run.info('Output file', zero_ouput)
    with open(zero_ouput, 'w') as file:
        file.write("contig\tlength\trange\trange_size\n")
        for contig, data, in cov_dict.items():
            contig_length = data['length']
            in_window = False
            window_start = ''
            window_end = ''
            window_length = ''
            for pos in range(data['length']):
                if data['cov'][pos] == 0 and in_window == False:
                    window_start = pos
                    in_window = True
                    file.write(f"{contig}\t{contig_length}\t{window_start}-")
                elif (data['cov'][pos] > 0 and in_window == True):
                    window_end = pos
                    window_length = window_end - window_start
                    in_window = False
                    file.write(f"{window_end}\t{window_length}\n")
                # if end of contig
                if data['cov'][pos] == 0 and pos == contig_length - 1:
                    if in_window:
                        window_end = pos + 1
                        window_length = window_end - window_start
                        in_window = False
                        file.write(f"{window_end}\t{window_length}\n")
                    else:
                        window_start = pos
                        window_end = pos + 1
                        window_length = window_end - window_start
                        in_window = False
                        file.write(f"{contig}\t{contig_length}\t{window_start}-{window_end}\t{window_length}\n")
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
