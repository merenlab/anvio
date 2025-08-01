#!/usr/bin/env python
# -*- coding: utf-8

import sys
import random

import anvio
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.learning import RF


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ["Tom O. Delmont"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Train a classifier for CPR prediction"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def run_program():
    args = get_args()

    filesnpaths.is_file_tab_delimited(args.matrix)
    filesnpaths.is_output_file_writable(args.output)

    cols = utils.get_columns_of_TAB_delim_file(args.matrix)
    if cols[0] != "class":
        raise ConfigError("This file doesn't seem to comply with the expected file format (the second "
                           "column is not 'class').")

    d = utils.get_TAB_delimited_file_as_dictionary(args.matrix)
    genome_to_class = dict([(g, d[g]['class']) for g in d])
    features = cols[1:]

    id_to_genome, genome_to_id, columns, vectors = utils.get_vectors_from_TAB_delim_matrix(args.matrix, cols_to_return = features)

    genomes = sorted(genome_to_class.keys())
    labels = [genome_to_class[genome] for genome in genomes]
    data = [vectors[genome_to_id[genome]] for genome in genomes]
    dim = len(data[0])

    # add some NON-CPR data points (this is garbage heuristic, and shold be improved,
    # fortunately it is pretty easy to do this part right)
    for i in range(0, int(len(genomes) * 2)):
        labels.append('NON-CPR')
        v = [1] * dim

        # add some random noise
        for j in random.sample(list(range(0, dim)), random.randint(0, int(dim / 2))):
            v[j] = 0

        data.append(v)

    rf = RF(args.output)
    rf.train(features, data, labels)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument('matrix', help = 'TAB-delimited matrix of CPR genome names, classes, and presence absence\
                                          of single-copy genes. Headers of the first two rows should be "genome", and\
                                          "class". The rest of the rows shold be single-copy genes.',metavar='MATRIX_FILE')
    parser.add_argument('-o', '--output', default = "cpr-scg.classifier", help = "Output file name for the classifier.",
                        metavar='CLASSIFIER_FILE')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
