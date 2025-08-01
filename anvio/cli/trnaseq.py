#!/usr/bin/env python
# -*- coding: utf-8
"""De novo prediction of tRNAs and their properties from a single tRNA-seq library"""

import sys

import anvio
import anvio.trnaseq as trnaseq
import anvio.trnaidentifier as trnaidentifier

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10']
__requires__ = ['trnaseq-fasta']
__provides__ = ['trnaseq-db']
__description__ = "A program to process reads from a tRNA-seq dataset to generate an anvi'o tRNA-seq database"


def main():
    args = get_args()

    try:
        if args.default_feature_param_file:
            parameterizer = trnaidentifier.TRNAFeatureParameterizer()
            parameterizer.write_param_file(args.default_feature_param_file)
        elif args.print_default_feature_params:
            parameterizer = trnaidentifier.TRNAFeatureParameterizer()
            print(parameterizer.tabulate_params())
        else:
            trnaseq_dataset = trnaseq.TRNASeqDataset(args)
            trnaseq_dataset.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1A = parser.add_argument_group('MANDATORY')
    group1A.add_argument(*anvio.A('trnaseq-fasta'), **anvio.K('trnaseq-fasta'))
    group1A.add_argument(*anvio.A('sample-name'),
                         **anvio.K('sample-name',
                                   {'help': "Unique sample name, considering all others in the experiment, "
                                            "that only includes ASCII letters and digits, without spaces"}))
    group1A.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))

    group1B = parser.add_argument_group('EXTRAS')
    group1B.add_argument(*anvio.A('treatment'), **anvio.K('treatment'))
    group1B.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))
    group1B.add_argument(*anvio.A('description'), **anvio.K('description'))

    group1C = parser.add_argument_group('ADVANCED')
    group1C.add_argument(*anvio.A('write-checkpoints'), **anvio.K('write-checkpoints'))
    group1C.add_argument(*anvio.A('load-checkpoint'), **anvio.K('load-checkpoint'))
    group1C.add_argument(*anvio.A('feature-param-file'), **anvio.K('feature-param-file'))
    group1C.add_argument(*anvio.A('threeprime-termini'), **anvio.K('threeprime-termini'))
    group1C.add_argument(*anvio.A('min-length-long-fiveprime'), **anvio.K('min-length-long-fiveprime'))
    group1C.add_argument(*anvio.A('min-trna-fragment-size'), **anvio.K('min-trna-fragment-size'))
    group1C.add_argument(*anvio.A('agglomeration-max-mismatch-freq'), **anvio.K('agglomeration-max-mismatch-freq'))
    group1C.add_argument(*anvio.A('skip-INDEL-profiling'),
                         **anvio.K('skip-INDEL-profiling',
                                   {'help': "This flag prevents the prediction of deletions in tRNA reads, which can save time."}))
    group1C.add_argument(*anvio.A('max-indel-freq'), **anvio.K('max-indel-freq'))
    group1C.add_argument(*anvio.A('left-indel-buffer'), **anvio.K('left-indel-buffer'))
    group1C.add_argument(*anvio.A('right-indel-buffer'), **anvio.K('right-indel-buffer'))

    group1D = parser.add_argument_group('PERFORMANCE')
    group1D.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    group1D.add_argument(*anvio.A('skip-fasta-check'), **anvio.K('skip-fasta-check'))
    group1D.add_argument(*anvio.A('profiling-chunk-size'), **anvio.K('profiling-chunk-size'))
    group1D.add_argument(*anvio.A('alignment-target-chunk-size'), **anvio.K('alignment-target-chunk-size'))

    group2 = parser.add_argument_group('DEFAULTS')
    group2.add_argument(*anvio.A('default-feature-param-file'), **anvio.K('default-feature-param-file'))
    group2.add_argument(*anvio.A('print-default-feature-params'), **anvio.K('print-default-feature-params'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
