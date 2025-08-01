#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.clustering as clustering

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['view-data',]
__provides__ = ['dendrogram',]
__description__ = "Takes a distance matrix, returns a newick tree"
__resources__ = [("See this program in action in the pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions")]


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()
    PL = terminal.pluralize

    if not args.output_file:
        args.output_file = args.input_matrix + '.newick'

    run.info('Input matrix file', args.input_matrix)

    try:
        # a quick check of the data to see if there are any missing
        # values
        _, _, _, vectors = utils.get_vectors_from_TAB_delim_matrix(args.input_matrix, run=run)
        num_missing_data = len([item for sublist in vectors for item in sublist if item == None])
        if num_missing_data:
            run.warning(f"Oy. Your file contains {PL('missing data item', num_missing_data)}. Anvi'o will do its "
                        f"best to accomodate for it, but you should take a careful look at your output tree (that, "
                        f"of course, if everything works after this message). Please note that anvi'o will not "
                        f"normalize your data prior to clustering and will assume that it is already normalized.")

        progress.new('Analyzing input file')
        clustering.create_newick_file_from_matrix_file(args.input_matrix, args.output_file, linkage=args.linkage,
                                        distance=args.distance, transpose=args.transpose, progress=progress,
                                        items_order_file_path=args.items_order_file, pad_with_zeros=args.pad_input_with_zeros,
                                        distance_matrix_output_path=args.distance_matrix)
        progress.end()
        run.info('Output newick', args.output_file)
        if args.items_order_file:
            run.info('Items order', args.items_order_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', 'The data you wish to cluster')
    groupA.add_argument('input_matrix', metavar='PATH', default=None,
                        help='Input matrix')

    groupB = parser.add_argument_group('OUTPUT', 'How would you like your results to be reported?')
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupB.add_argument('--items-order-file', metavar='FILE PATH', default=None,
                        help="In addition to a newick formatted output file, you can ask anvi'o to report\
                              the order of items in the resulting tree in a separate file. The content of\
                              this file will be a single-column item names the way they are ordered in the\
                              output newick dendrogram.")
    groupB.add_argument('--distance-matrix', metavar='FILE PATH', default=None,
                        help="In addition to a newick formatted output file, you can ask anvi'o to also report "
                             "a distance matrix computed from the same data that leads to the generation.")


    groupC = parser.add_argument_group('SWEETS', 'Additional options')
    groupC.add_argument(*anvio.A('transpose'), **anvio.K('transpose'))
    groupC.add_argument(*anvio.A('distance'), **anvio.K('distance'))
    groupC.add_argument(*anvio.A('linkage'), **anvio.K('linkage'))
    groupC.add_argument('--pad-input-with-zeros', default=False, action="store_true",
                        help="Do you want your vectors in the input file to be 'padded' by some zero values "
                             "so even hideous matrix files that contain entires with all empty values can be "
                             "processed with this program? Well, this flag will try to make things work for "
                             "you. But you must know that this is one of those cases where 'if you have to "
                             "use this then you have not one but two problems.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
