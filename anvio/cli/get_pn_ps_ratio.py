#!/usr/bin/env python
# -*- coding: utf-8
"""
Citations: doi:10.1126/science.aaz9642, doi:10.1038/nature11711
"""

import os
import sys

import numpy as np

import anvio

import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.variabilityops import VariabilityData

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__provides__ = ['pn-ps-data']
__requires__ = ['contigs-db', 'variability-profile-txt']
__description__ = ("Calculate the rates of non-synonymous and synonymous polymorphism for genes across environmetns "
                   "using the output of anvi-gen-variability-profile.")


@terminal.time_program
def main():
    args = get_args()

    try:
        calculate_pN_pS_ratio(args)
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(1)


def load_variability(args, contigs_db):
    run = terminal.Run()
    progress = terminal.Progress()

    progress.new('Loading SCVs')
    progress.update('...')

    var = VariabilityData(args)

    # filter based on user parameters
    var.filter_data(criterion='departure_from_consensus', verbose=False)
    var.filter_data(criterion='departure_from_reference', verbose=False)
    var.filter_data(criterion='coverage', verbose=False)

    # identify gene calls that are noncoding
    gene_caller_ids_for_noncoding_gene_calls = [g for g in var.data['corresponding_gene_call'].unique()
                                                if contigs_db.genes_in_contigs_dict[g]['call_type'] != constants.gene_call_types['CODING']]
    if gene_caller_ids_for_noncoding_gene_calls:
        # This should never happen with >v6.2
        var.data.drop(var.data[var.data['corresponding_gene_call'].isin(gene_caller_ids_for_noncoding_gene_calls)].index, inplace=True)

        progress.reset()
        run.warning("%d of your gene calls were 'noncoding', and were removed from downstream analyses. Here is the complete list of "
                    "gene calls that were removed: '%s'." % (len(gene_caller_ids_for_noncoding_gene_calls), ', '.join([str(g) for g in gene_caller_ids_for_noncoding_gene_calls])))

    progress.end()
    return var


def load_contigs_db(args):
    filesnpaths.is_file_exists(args.contigs_db)
    contigs_db = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
    contigs_db.init_contig_sequences()
    return contigs_db


def report(args, pNpS, pN, pS, num_SCVs, potentials):
    run = terminal.Run()

    # write it to folder
    pNpS.to_csv(os.path.join(args.output_dir, "pNpS.txt"), sep="\t", index = True)
    pN.to_csv(os.path.join(args.output_dir, "pN.txt"), sep="\t", index = True)
    pS.to_csv(os.path.join(args.output_dir, "pS.txt"), sep="\t", index = True)
    num_SCVs.to_csv(os.path.join(args.output_dir, "num_SCVs.txt"), sep="\t", index = True)
    potentials.to_csv(os.path.join(args.output_dir, "potentials.txt"), sep="\t", index = True)

    run.info_single("Done! Contents have been output to the directory '{}'.".format(args.output_dir),
                    nl_before=1,
                    nl_after=1)


def calculate_pN_pS_ratio(args):
    # gen output
    filesnpaths.check_output_directory(args.output_dir)
    filesnpaths.gen_output_directory(args.output_dir)

    # load contigs db and variability table
    contigs_db = load_contigs_db(args)

    var = load_variability(args, contigs_db)
    var.calc_pN_pS(
        contigs_db=contigs_db,
        grouping='gene',
        comparison=args.comparison,
        add_potentials=True
    )

    potentials = var.data.groupby('corresponding_gene_call')[f'nS_gene_{args.comparison}', f'nN_gene_{args.comparison}'].first()

    groupby_vars = args.groupby.split(',')

    pS_name = f'pS_gene_{args.comparison}'
    pN_name = f'pN_gene_{args.comparison}'
    pNpS_name = f'pNpS_gene_{args.comparison}'

    pN = var.data.\
        groupby(groupby_vars)\
        [pN_name].\
        sum()

    pS = var.data.\
        groupby(groupby_vars)\
        [pS_name].\
        sum()

    pNpS = pN/pS

    SCVs_per_group = var.data.\
        groupby(groupby_vars)\
        [args.comparison].\
        size()

    if args.minimum_num_variants > 1:
        pN[SCVs_per_group < args.minimum_num_variants] = np.nan
        pS[SCVs_per_group < args.minimum_num_variants] = np.nan
        pNpS[SCVs_per_group < args.minimum_num_variants] = np.nan

    pN = pN.reset_index().rename(columns={0: pN_name})
    pS = pS.reset_index().rename(columns={0: pS_name})
    pNpS = pNpS.reset_index().rename(columns={0: pNpS_name})
    SCVs_per_group = SCVs_per_group.reset_index().rename(columns={args.comparison: 'num_SCVs'})

    if args.pivot:
        pN = pN.pivot(index=groupby_vars[0], columns=groupby_vars[1], values=pN_name)
        pS = pS.pivot(index=groupby_vars[0], columns=groupby_vars[1], values=pS_name)
        pNpS = pNpS.pivot(index=groupby_vars[0], columns=groupby_vars[1], values=pNpS_name)
        SCVs_per_group = SCVs_per_group.\
            pivot(
                index=groupby_vars[0],
                columns=groupby_vars[1],
                values='num_SCVs'
            ).\
            fillna(0)

    report(args, pNpS, pN, pS, SCVs_per_group, potentials)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupV = parser.add_argument_group('VARIABILITY', 'Provide a SCV table that can be generated with anvi-gen-variability-profile.')
    groupV.add_argument('-V', '--variability-profile', help='Filepath to the SCV table.', metavar='SCV_FILE')
    groupV.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'help':'Filepath to the contigs database used \
                                                      to generate variability table.'}))

    groupE = parser.add_argument_group('TUNABLES', "Successfully tune one or more of these parameters to unlock the badge 'Advanced anvian'.")
    groupE.add_argument(*anvio.A('min-departure-from-consensus'), **anvio.K('min-departure-from-consensus', {'default': 0.0, 'help': \
                            'SCVs will be ignored if they have a departure from consensus less than this \
                            value. Note: Keep in mind you may have already supplied this parameter during anvi-gen-variability-profile.\
                            The default value is %(default).2f.'}))
    groupE.add_argument(*anvio.A('min-departure-from-reference'), **anvio.K('min-departure-from-reference', {'default': 0.0, 'help': \
                            'SCVs will be ignored if they have a departure from reference less than this \
                            value. Note: Keep in mind you may have already supplied this parameter during anvi-gen-variability-profile.\
                            The default value is %(default).2f.'}))
    groupE.add_argument('-i', '--minimum-num-variants', default=4, type=int, required=False, help='Ignore groups with less than this number\
                            of single codon variants. pN, pS, and pN/pS values will be set to NaN for these groups')
    groupE.add_argument('-m', '--min-coverage', default=30, type=int, required=False, help='If the coverage value at a codon is less than \
                            this amount, any associated SCVs will be ignored. The default is %(default)d.')
    groupE.add_argument('-x', '--comparison', default='reference', choices=['reference', 'consensus', 'popular_consensus'], help='You can determine synonymity \
                            relative to either the reference codon, or the consensus codon. The consensus codon is determined on a per-sample \
                            basis. The default is \'%(default)s.\'')

    groupO = parser.add_argument_group('OUTPUT', 'The output of this program is a folder directory with several tables.')
    groupO.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required':True}))
    groupO.add_argument('-p', '--pivot', action='store_true', help = \
                            'By default the output is in long format, however you can \
                             choose the output to be in matrix form with this flag. If you\'re not sure which one is right for you, \
                             just try one and take a look at the output--there is no cost for making a mistake :)')

    args = parser.get_args(parser)
    args.groupby = 'corresponding_gene_call,sample_id'
    args.columns_to_load = constants.codons + [
        'unique_pos_identifier',
        'corresponding_gene_call',
        'sample_id',
        'coverage',
        'departure_from_consensus',
        'departure_from_reference',
        'reference',
        'consensus',
    ]

    args.columns_to_load = list(set(args.columns_to_load))

    return args


if __name__ == '__main__':
    main()
