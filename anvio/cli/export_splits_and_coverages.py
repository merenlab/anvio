#!/usr/bin/env python
# -*- coding: utf-8
"""A script to export a FASTA file and the coverage table from a merged database.

The purpose of this script is to export two critical files for most genome binning
software: the sequences file for tetra-nucleotide analysis, and the coverage table
that shows the coverage of each contig across samples. These output files can be
used to identify bins, and those bins can be used to populate collections table
via available parsers in anvi-populate-collections table."""

import os
import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['profile-db', 'contigs-db']
__provides__ = ['contigs-fasta', 'coverages-txt']
__description__ = ("Export split or contig sequences and coverages across samples stored in an anvi'o "
                   "profile database. This program is especially useful if you would like to 'bin' your "
                   "splits or contigs outside of anvi'o and import the binning results into anvi'o "
                   "using `anvi-import-collection` program")


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        profile_db = dbops.ProfileSuperclass(args)
        samples = profile_db.p_meta['samples']

        if args.output_dir:
            filesnpaths.gen_output_directory(args.output_dir)
        else:
            args.output_dir = os.path.dirname(os.path.abspath(args.profile_db))

        if not args.output_file_prefix:
            args.output_file_prefix = profile_db.p_meta['sample_id']

        run.info('Output directory', args.output_dir)
        run.info('Filename prefix', args.output_file_prefix)
        run.info('Splits mode', args.splits_mode)
        run.info('Report Q2Q3 coverages', args.use_Q2Q3_coverages)
        run.info('Report contigs and not splits', args.report_contigs)

        coverages_output = os.path.join(args.output_dir, args.output_file_prefix + '-COVs.txt')
        sequences_output = os.path.join(args.output_dir, args.output_file_prefix + ('-CONTIGS.fa' if args.report_contigs else '-SPLITS.fa'))

        progress.new('Bleep bloop')

        progress.update('Gathering coverage data...')
        coverages = profile_db.get_split_coverages_dict(use_Q2Q3_coverages=args.use_Q2Q3_coverages, splits_mode=args.splits_mode, report_contigs=args.report_contigs)
        utils.store_dict_as_TAB_delimited_file(coverages, coverages_output, ['contig'] + sorted(list(samples)))

        progress.update('Dealing with sequences...')
        progress.update('Gathering coverage data')
        utils.export_sequences_from_contigs_db(args.contigs_db, sequences_output, seq_names_to_export=sorted(list(coverages.keys())), splits_mode=(not args.report_contigs))

        progress.end()

        run.info('Coverages file', coverages_output, mc="green", nl_before=1)
        run.info('Sequences file', sequences_output, mc="green")
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    parser.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    parser.add_argument(*anvio.A('splits-mode'), **anvio.K('splits-mode'))
    parser.add_argument('--report-contigs', default=False, action="store_true", help="By default this program reports \
                         sequences and their coverages for 'splits'. By using this flag, you can report contig sequences\
                         and coverages instead. For obvious reasons, you can't use this flag with `--splits-mode` flag.")
    parser.add_argument('--use-Q2Q3-coverages', default=False, action="store_true", help="By default this program reports the mean \
                         coverage of a split (or contig, see --report-contigs) for each sample. By using this flag, \
                         you can report the mean Q2Q3 coverage by excluding 25 percent of the nucleotide positions with the \
                         smallest coverage values, and 25 percent of the nucleotide positions with the largest coverage values. \
                         The hope is that this removes 'outlier' positions resulting from non-specific mapping, etc. \
                         that skew the mean coverage estimate.")

    return parser.get_args(parser)

if __name__ == '__main__':
    main()
