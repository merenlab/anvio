#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from functools import partial
from itertools import combinations, product

from anvio.dbinfo import DBInfo
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__requires__ = ['contigs-db', 'profile-db']
__provides__ = ['contigs-fasta']
__description__ = ("This script generates a FASTA file of tRNA-seq seeds "
                   "with permuted nucleotides at positions of predicted modification-induced substitutions. "
                   "The underlying nucleotide without modification is not always the most common base call. "
                   "The resulting FASTA file can be queried against a database of tRNA genes "
                   "to validate nucleotides at modified positions and find the most similar sequences.")


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
    pp = terminal.pretty_print

    contigs_db_path = args.contigs_db
    profile_db_path = args.specific_profile_db
    fasta_path = args.contigs_fasta
    min_nt_freq = args.min_nt_frequency
    max_var_positions = args.max_variable_positions
    overwrite_out_dest = args.overwrite_output_destinations

    contigs_db_info = DBInfo(contigs_db_path, expecting='contigs')
    if contigs_db_info.variant != 'trnaseq':
        raise ConfigError("The provided contigs database must be a 'trnaseq' variant. "
                          f"Yours was a '{contigs_db_info.variant}' variant.")

    profile_db_info = DBInfo(profile_db_path, expecting='profile')
    if profile_db_info.variant != 'trnaseq':
        raise ConfigError("The provided profile database must be a 'trnaseq' variant. "
                          f"Yours was a '{profile_db_info.variant}' variant.")

    filesnpaths.is_output_file_writable(fasta_path, ok_if_exists=overwrite_out_dest)


    with contigs_db_info.load_db() as contigs_db:
        contig_names_seqs = contigs_db.get_table_as_list_of_tuples('contig_sequences')

    with profile_db_info.load_db() as profile_db:
        var_nts_gb = get_var_nts_gb(profile_db)

    fasta = open(fasta_path, 'w')
    write_entries_to_fasta = partial(write_entries, fasta=fasta)
    get_var_nts_dict_with_min_nt_freq = partial(get_var_nts_dict, min_nt_freq=min_nt_freq)

    seed_count = -1
    seed_total = len(contig_names_seqs)
    permut_seq_count = 0
    max_permut_seq_count = 0
    progress.new("Writing output")
    for contig_name, contig_seq in contig_names_seqs:
        seed_count += 1
        if seed_count % 100 == 0:
            progress.update(f"{seed_count}/{seed_total} seeds")

        try:
            contig_var_nts_df = var_nts_gb.get_group(contig_name)
        except KeyError:
            write_entries_to_fasta(contig_name, contig_seq)
            continue

        var_nts_dict = get_var_nts_dict_with_min_nt_freq(contig_var_nts_df)
        if not var_nts_dict:
            write_entries_to_fasta(contig_name, contig_seq)
            continue

        permut_contigs_info = get_permut_contigs_info(var_nts_dict, max_var_positions, contig_seq)

        remove_orig_contig(permut_contigs_info, contig_seq)
        permut_seq_count += len(permut_contigs_info)
        max_permut_seq_count = max(len(permut_contigs_info), max_permut_seq_count)
        write_entries_to_fasta(contig_name, contig_seq, permut_contigs_info=permut_contigs_info)
    progress.end()

    run.info("Mean permuted seqs per contig", round(permut_seq_count / seed_total, 1))
    run.info("Max permuted seqs from a contig", pp(max_permut_seq_count))


def get_var_nts_gb(profile_db):
    var_nts_df = profile_db.get_table_as_dataframe('variable_nucleotides', columns_of_interest=['split_name', 'pos', 'A', 'C', 'G', 'T'])
    var_nts_df['contig_name'] = var_nts_df['split_name'].apply(lambda s: s.split('_split_00001')[0])
    var_nts_df = var_nts_df.drop('split_name', axis=1)
    var_nts_df = var_nts_df.set_index('contig_name')
    return var_nts_df.groupby('contig_name')


def get_var_nts_dict(contig_var_nts_df, min_nt_freq):
    var_nts_dict = {}
    for pos, contig_var_nts_df in contig_var_nts_df.set_index('pos').groupby('pos', sort=False):
        total_nts_df = contig_var_nts_df.sum(axis=0)
        total_nts_freq_series = total_nts_df / total_nts_df.sum()
        total_nts_series = total_nts_df[total_nts_freq_series > min_nt_freq]
        total_nts_series = total_nts_series.loc[~total_nts_series.index.isin([total_nts_series.idxmax()])]

        if len(total_nts_series) == 0:
            continue

        var_nts_dict[pos] = tuple(total_nts_series.index)
    return var_nts_dict


def write_entries(contig_name, contig_seq, fasta, permut_contigs_info=None):
    fasta.write(f">{contig_name}|\n{contig_seq}\n")
    if not permut_contigs_info:
        return
    for contig_seq, permut_positions, permut_nts in permut_contigs_info:
        fasta.write(f">{contig_name}|{'_'.join([str(permut_pos) + permut_nt for permut_pos, permut_nt in zip(permut_positions, permut_nts)])}\n{contig_seq}\n")


def get_permut_contigs_info(var_nts_dict, max_var_positions, contig_seq):
    permut_contigs_info = []
    for num_var_positions in range(1, min(len(var_nts_dict), max_var_positions) + 1):
        permut_combos = combinations(var_nts_dict, num_var_positions)
        for permut_positions in permut_combos:
            for permut_nts in product(*[var_nts_dict[pos] for pos in permut_positions]):
                permut_contig_seq = contig_seq
                for pos, nt in zip(permut_positions, permut_nts):
                    permut_contig_seq = permut_contig_seq[: pos] + nt + permut_contig_seq[pos + 1: ]
                permut_contigs_info.append((permut_contig_seq, tuple(permut_positions), tuple(permut_nts)))
    return permut_contigs_info


def remove_orig_contig(permut_contigs_info, contig_seq):
    orig_index = None
    for i, permut_contig_info in enumerate(permut_contigs_info):
        if permut_contig_info[0] == contig_seq:
            orig_index = i
            break
    if orig_index is None:
        return False
    permut_contigs_info.pop(orig_index)
    return True


def get_args():
    parser = ArgumentParser()

    parser.add_argument_group('MANDATORY')
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('specific-profile-db'), **anvio.K('specific-profile-db'))
    parser.add_argument(*anvio.A('contigs-fasta'),
                        **anvio.K('contigs-fasta', {'help': 'FASTA file to generate'}))

    parser.add_argument_group('OPTIONAL')
    parser.add_argument('-n', '--min-nt-frequency',
                        metavar='FLOAT', type=float,
                        default=0.05,
                        help="For a position in a contig, this is the minimum nucleotide frequency, summed across all samples, "
                             "required for the nucleotide to be substituted in permuted sequences.")
    parser.add_argument('-x', '--max-variable-positions',
                        metavar='INT', type=int,
                        default=5,
                        help="The maximum number of modified positions that can be permuted at once in a contig.")
    parser.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    return parser.parse_args()


if __name__ == '__main__':
    main()
