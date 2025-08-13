#!/usr/bin/env python
# -*- coding: utf-8

"""Return short reads from BAM files mapping to gene."""

import os
import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.dbops import ContigsSuperclass
from anvio.bamops import ReadsMappingToARange
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'bam-file',]
__provides__ = ['short-reads-fasta',]
__tags__ = ["metagenomics", "profile_db", "contigs_db", "bam", "variability", "clustering"]
__description__ = ("Recover short reads from BAM files that were mapped to genes you are "
                   "interested in. It is possible to work with a single gene call, or a "
                   "bunch of them. Similarly, you can get short reads from a single "
                   "BAM file, or from many of them")


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
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    c = ContigsSuperclass(args)

    if args.gene_caller_id and args.genes_of_interest:
        raise ConfigError("You either should define a single gene caller id, or a list of "
                          "gene callers of interest. Not both")

    if not (args.genes_of_interest or args.gene_caller_id):
        raise ConfigError("You must provide either a single gene caller id, or a file with one or more gene "
                          "caller ids. You provided none of these. Well .. *starts quietly taking notes*.")

    if args.genes_of_interest:
        genes_of_interest = set([int(g) for g in utils.get_column_data_from_TAB_delim_file(args.genes_of_interest, column_indices=[0], expected_number_of_fields=1)[0] if g])
    else:
        genes_of_interest = set([int(args.gene_caller_id)])

    output_file_prefix = os.path.abspath(args.output_file_prefix or "SHORT_READS")
    filesnpaths.is_output_dir_writable(os.path.dirname(output_file_prefix))

    run.info('Num genes of interest', len(genes_of_interest))
    run.info('Leeway', args.leeway)
    run.info('Output file prefix', output_file_prefix)

    missing_gene_caller_ids = [g for g in genes_of_interest if g not in c.genes_in_contigs_dict]
    if missing_gene_caller_ids:
        msg = "Some of the unique gene caller IDs you are interested in (precisely %d of %d), are missing\
               from the contigs database you are using. " % (len(missing_gene_caller_ids), len(genes_of_interest))

        if anvio.DEBUG:
            raise ConfigError(msg + "Here are the list of gene caller IDs that's making anvi'o upset: %s" % \
                                    (", ".join([str(g) for g in missing_gene_caller_ids])))
        else:
            error_msg = msg + "Since this can be a long list in some cases, anvi'o will not show you what\
                               they are. Unless you use the `--debug` flag."
            raise ConfigError(error_msg)

    r = ReadsMappingToARange(run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))

    progress.new("Processing gene calls", progress_total_items=len(genes_of_interest))
    progress.update("...")
    for gene_caller_id in genes_of_interest:
        gene_call = c.genes_in_contigs_dict[gene_caller_id]
        mapping_range_start = gene_call['start'] - args.leeway
        mapping_range_stop = gene_call['stop'] + args.leeway
        gene_length = gene_call['stop'] - gene_call['start']

        progress.increment()
        progress.update("Working on gene call '%d' (%d nts) in contig '%s'.." % (gene_caller_id, gene_length, gene_call['contig']))

        r.process_range(args.input_files, gene_call['contig'], mapping_range_start, mapping_range_stop)
        r.report(output_file_prefix + '_GENE_CALL_%d.txt' % gene_caller_id)
    progress.end()


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES', "An anvi'o contigs database and one or more BAM files.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    groupA.add_argument('-i', '--input-files', metavar = 'INPUT_BAM(S)', nargs='+', default = None, required = True,
                        help = 'Sorted and indexed BAM files to analyze. It is essential that all BAM files must be\
                                the result of mappings against the same contigs.')

    groupB = parser.add_argument_group('GENES', "Gene calls you want to work with")
    groupB.add_argument(*anvio.A('gene-caller-id'), **anvio.K('gene-caller-id'))
    groupB.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupB.add_argument(*anvio.A('leeway'), **anvio.K('leeway'))

    groupC = parser.add_argument_group('OUTPUT', "How should results be stored.")
    groupC.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
