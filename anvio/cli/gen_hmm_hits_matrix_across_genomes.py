#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.hmmopswrapper as hmmopswrapper

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["external-genomes", "internal-genomes", "hmm-source", "hmm-hits"]
__provides__ = ["hmm-hits-across-genomes-txt"]
__description__ = ("A simple script to generate a TAB-delimited file that reports the frequency "
                   "of HMM hits for a given HMM source across contigs databases")


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
    args, unknown = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    hmm_source = A('hmm_source') or set([])
    output_file_path = A('output_file')

    if not hmm_source:
        raise ConfigError("You need to declare an HMM source for this to work :/")

    s = hmmopswrapper.SequencesForHMMHitsWrapperForMultipleContigs(args, hmm_sources=set([hmm_source]))

    HMM_sources_common_to_all = s.get_HMM_sources_common_to_all_genomes()
    if not len(HMM_sources_common_to_all):
        raise ConfigError("There is not a single HMM source that is common to all "
                          "contigs databases you have :/")

    if hmm_source not in HMM_sources_common_to_all:
        raise ConfigError('The HMM source "%s" is not common to all of your contigs '
                          'databases you want to work with :/ Here is a list of those '
                          'that are common: %s' % (hmm_source, ', '.join(HMM_sources_common_to_all)))

    if not output_file_path:
        raise ConfigError("This will not work without an output file path.")

    s = hmmopswrapper.SequencesForHMMHitsWrapperForMultipleContigs(args, set([hmm_source]), run=terminal.Run(verbose=False))

    filesnpaths.is_output_file_writable(output_file_path)

    gene_names_in_source = [g.strip() for g in s.hmm_hits_info[hmm_source]['genes'].split(',')]

    d = {}
    progress.new('Processing contigs databases')
    for genome_name in s.genomes:
        d[genome_name] = {}
        for gene_name in gene_names_in_source:
            d[genome_name][gene_name] = 0

    for hit in s.hmm_hits.values():
        genome_hash = hit['genome_hash']
        genome_name = s.genome_hash_to_genome_name[genome_hash]
        gene_name = hit['gene_name']
        d[genome_name][gene_name] += 1

    utils.store_dict_as_TAB_delimited_file(d, output_file_path, headers=['genome_or_bin'] + sorted(gene_names_in_source))

    run.info('Output', output_file_path)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupB = parser.add_argument_group('INPUT: INTERNAL/EXTERNAL GENOMES FILE', "Yes. You need to use an internal and/or external genomes file\
                                        to tell anvi'o where your contigs databases are.")
    groupB.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupB.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))

    groupD = parser.add_argument_group('HMM STUFF', "This is where you can specify an HMM source, and/or a list of genes to filter\
                                        your results.")
    groupD.add_argument(*anvio.A('hmm-source'), **anvio.K('hmm-source'))
    groupD.add_argument(*anvio.A('list-hmm-sources'), **anvio.K('list-hmm-sources'))

    groupD = parser.add_argument_group('OUTPUTTAH')
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True }))

    return parser.parse_known_args()


if __name__ == "__main__":
    main()
