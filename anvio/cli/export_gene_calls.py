#!/usr/bin/env python
# -*- coding: utf-8

import sys
import pandas as pd

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
__authors__ = ['meren', 'ekiefl']
__requires__ = ['contigs-db']
__provides__ = ['gene-calls-txt']
__description__ = "Export gene calls from an anvi'o contigs database"


def get_aa_seq_dict(contigs_db, gene_caller_ids):
    aa_seq_dict = contigs_db.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gene_caller_ids, include_aa_sequences=True)[1]

    return {k: v['aa_sequence'] for k, v in aa_seq_dict.items()}


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_contigs_db(args.contigs_db)

        if args.list_gene_callers:
            dbops.ContigsDatabase(args.contigs_db).list_gene_caller_sources()
            sys.exit()

        if not args.output_file:
            raise ConfigError("You should provide an output file name so anvi'o does not have to make up a silly "
                               "name :/")

        if not args.gene_caller:
            raise ConfigError("You must provide a source for your gene calls. You can always use the flag "
                              "`--list-gene-callers` to see what is available in your contigs-db.")

        filesnpaths.is_output_file_writable(args.output_file)

        contigs_db = dbops.ContigsSuperclass(args)

        if not contigs_db.genes_in_contigs_dict:
            raise ConfigError("Something weird happened. The contigs database does not seem to have any gene calls :/")

        # convert gene calls initialized in ContigsSuperclass to dataframe
        gene_calls = pd.DataFrame(contigs_db.genes_in_contigs_dict).T
        available_sources = gene_calls['source'].unique()

        if not args.gene_caller:
            raise ConfigError("Please provide one or more comma-separated gene callers (--gene-caller). Use "
                              "--list-gene-callers to see your options")

        requested_sources = args.gene_caller.split(',')
        for requested_source in requested_sources:
            if requested_source not in available_sources:
                raise ConfigError("%s is not a gene caller in the contigs database. Here is a list of available "
                                 "gene callers: %s." % (requested_source, ','.join(available_sources)))

        run.info('Gene callers that will be used', ', '.join(requested_sources))

        # filter to only include requested gene callers
        gene_calls = gene_calls[gene_calls['source'].isin(requested_sources)]

        for available_source, count in gene_calls['source'].value_counts().to_dict().items():
            run.info('Num gene calls with %s' % available_source, count)
        run.info('Num gene calls total', gene_calls.shape[0])

        # also append amino acid sequences
        if not args.skip_sequence_reporting:
            gene_caller_ids_of_interest = list(gene_calls.index.values)
            gene_calls['aa_sequence'] = gene_calls.index.map(get_aa_seq_dict(contigs_db, gene_caller_ids_of_interest))

        utils.store_dataframe_as_TAB_delimited_file(gene_calls, args.output_file, include_index=True, index_label='gene_callers_id')

        run.info('Output', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SOURCE FOR GENE DATA', "You want to export gene calls from where?")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group('SOURCE FOR GENE CALLS')
    groupB.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller', {'help': "Which gene caller(s) "
                "would you like to export gene calls for? If providing multiple they should be comma-separated "
                "(no spaces). If you don't know, use --list-gene-callers.", 'default': None,}))
    groupB.add_argument(*anvio.A('list-gene-callers'), **anvio.K('list-gene-callers'))

    groupC = parser.add_argument_group('MATTERS OF REPORTING')
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupC.add_argument('--skip-sequence-reporting', action = 'store_true', help = "By default, exported gene "
                "calls have an amino acid sequences column in the output. Turn this behavior off with this flag")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
