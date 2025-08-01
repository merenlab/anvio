#!/usr/bin/env python
# -*- coding: utf-8

"""Oligotyping analysis of the linkmer reports."""

import os
import sys
from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['linkmers-txt',]
__provides__ = ['oligotypes',]
__resources__ = [("An application of the oligotyping workflow in metagenomics", "https://merenlab.org/2015/12/09/musings-over-commamox/#an-application-of-oligotyping-in-the-metagenomic-context-oligotyping-amoc")]
__description__ = "Takes an anvi'o linkmers report, generates an oligotyping output"


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
    args = get_args()
    run = terminal.Run()

    filesnpaths.is_file_tab_delimited(args.input_file)
    filesnpaths.is_output_dir_writable(args.output_dir)

    report = utils.get_TAB_delimited_file_as_dictionary(args.input_file)

    # take care of unique hashes
    hash_to_oligotype = {}
    unique_hashes = set([entry['read_unique_id'] for entry in list(report.values())])
    for unique_hash in unique_hashes:
        hash_to_oligotype[unique_hash] = []

    for entry in list(report.values()):
        hash_to_oligotype[entry['read_unique_id']].append((entry['pos_in_contig'], entry['base']),)

    for unique_hash in unique_hashes:
        hash_to_oligotype[unique_hash] = ''.join([e[1] for e in sorted(hash_to_oligotype[unique_hash])])

    d = {}

    request_ids = sorted(list(set([entry['request_id'] for entry in list(report.values())])))
    samples = sorted(list(set([entry['sample_id'] for entry in list(report.values())])))

    for request_id in request_ids:
        d[request_id] = {}
        for sample_id in samples:
            d[request_id][sample_id] = Counter()

    for entry in list(report.values()):
        request_id, sample_id, unique_hash = entry['request_id'], entry['sample_id'], entry['read_unique_id']
        oligotype = hash_to_oligotype[unique_hash]
        d[request_id][sample_id][oligotype] += 1

    # at this point our dict has inflated number of observations, and here we will normalize those observations
    # by dividing the total count to the number of nucleotide positions used for oligotyping. I have a feeling
    # that this might blow up sooner or later, but I will think of a better solution then.
    for request_id in d:
        for sample_id in d[request_id]:
            try:
                normalization_factor = len(list(list(d[request_id].values())[0].keys())[0])
            except IndexError:
                normalization_factor = 1.0

            for oligotype in d[request_id][sample_id]:
                d[request_id][sample_id][oligotype] = d[request_id][sample_id][oligotype] / normalization_factor

    run.warning('', header = "Oligotyping outputs per request", lc = 'cyan')

    for request_id in request_ids:
        output_file_path = os.path.join(args.output_dir, 'oligotype-counts-%s.txt' % request_id)
        utils.store_dict_as_TAB_delimited_file(d[request_id], output_file_path)
        run.info('Output', output_file_path)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument('-i', '--input-file', metavar = 'LINKMER_REPORT', required = True,
                        help = 'Output file of `anvi-report-linkmers`.')
    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
