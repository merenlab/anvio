#!/usr/bin/env python
"""A script to export a FASTA file of contigs (or splits) from a contigs database."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


pp = terminal.pretty_print
P = terminal.pluralize


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db']
__provides__ = ['contigs-fasta']
__can_use__ = ['genes-of-interest-txt']
__description__ = "Export contigs (or splits) from an anvi'o contigs database"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        if args.contigs_of_interest and args.genes_of_interest:
            raise ConfigError("You can use either --contigs-of-interest or --genes-of-interest to limit the "
                              "sequences anvi'o will export, but not both of them at the same time. Please remove "
                              "one of these parameters from your command and try again.")

        if args.contigs_of_interest:
            filesnpaths.is_file_tab_delimited(args.contigs_of_interest, expected_number_of_fields=1)
            seq_names_to_export = [line.strip() for line in open(args.contigs_of_interest).readlines()]
        elif args.genes_of_interest:
            seq_names_to_export = get_seq_names_for_gene_caller_ids_of_interest(args, run=run)
        else:
            seq_names_to_export = None

        utils.export_sequences_from_contigs_db(args.contigs_db,
                                               args.output_file,
                                               seq_names_to_export=seq_names_to_export,
                                               splits_mode=args.splits_mode,
                                               just_do_it=args.just_do_it,
                                               truncate=(not args.no_wrap),
                                               run=run)

        run.info('Export mode', 'splits' if args.splits_mode else 'contigs')
        run.info('Output FASTA', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_seq_names_for_gene_caller_ids_of_interest(args, run=terminal.Run()):
    """Resolves a file of gene caller ids into the names of sequences that contain them.

    Depending on whether the user is in splits mode or not, the resulting names will be
    split names or contig names.

    Parameters
    ==========
    args : argparse.Namespace
        Must contain `contigs_db`, `genes_of_interest`, and `splits_mode`.

    Returns
    =======
    seq_names_to_export : list
        Sorted list of unique contig (or split) names that contain the gene calls of interest.
    """

    filesnpaths.is_file_tab_delimited(args.genes_of_interest, expected_number_of_fields=1)

    gene_caller_ids_of_interest = set([line.strip() for line in open(args.genes_of_interest).readlines() if line.strip()])

    if not len(gene_caller_ids_of_interest):
        raise ConfigError(f"The file you have provided via --genes-of-interest at '{args.genes_of_interest}' "
                          f"seems to be empty :/")

    try:
        gene_caller_ids_of_interest = set([int(g) for g in gene_caller_ids_of_interest])
    except ValueError:
        raise ConfigError(f"At least one of the entries in your genes of interest file at "
                          f"'{args.genes_of_interest}' does not look like an anvi'o gene caller id (which are "
                          f"always integers). Please make sure the file has a single column of gene caller ids and "
                          f"no column header.")

    contigs_db = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False))

    missing_gene_caller_ids = gene_caller_ids_of_interest.difference(set(contigs_db.genes_in_contigs_dict.keys()))
    if len(missing_gene_caller_ids):
        raise ConfigError(f"Anvi'o is having a hard time here: {P('gene call', len(missing_gene_caller_ids))} in your "
                          f"genes of interest file {P('is', len(missing_gene_caller_ids), alt='are')} not in the "
                          f"contigs database at '{args.contigs_db}'. Here "
                          f"{P('it is', len(missing_gene_caller_ids), alt='are some of them')}: "
                          f"{', '.join([str(g) for g in sorted(missing_gene_caller_ids)[:10]])}.")

    if args.splits_mode:
        # a single gene call may be spread across more than one split, and the `genes_in_splits`
        # table has one entry per gene call per split, so we collect every split a gene appears in
        seq_names_to_export = set([])
        for entry in contigs_db.genes_in_splits.values():
            if entry['gene_callers_id'] in gene_caller_ids_of_interest:
                seq_names_to_export.add(entry['split'])
    else:
        seq_names_to_export = set([contigs_db.genes_in_contigs_dict[g]['contig'] for g in gene_caller_ids_of_interest])

    run.info('Genes of interest', f"{pp(len(gene_caller_ids_of_interest))} gene calls found in the input file")
    run.info('Sequences that match', f"{pp(len(seq_names_to_export))} {'splits' if args.splits_mode else 'contigs'}")

    return sorted(seq_names_to_export)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('contigs-of-interest'), **anvio.K('contigs-of-interest'))
    parser.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    parser.add_argument('--splits-mode', default=False, action="store_true", help="Export split\
                        sequences instead.")
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    parser.add_argument(*anvio.A('no-wrap'), **anvio.K('no-wrap'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
