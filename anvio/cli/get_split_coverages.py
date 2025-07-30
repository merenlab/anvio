#!/usr/bin/env python
# -*- coding: utf-8
"""A script to export coverage values across samples of a given split in an anvi'o contigs db"""

import sys

from operator import itemgetter

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError, FilesNPathsError
from anvio.argparse import ArgumentParser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["coverages-txt",]
__requires__ = ["profile-db", "contigs-db", "collection", "bin"]
__description__ = "Export splits and the coverage table from database"
__resources__ = [("Using this program to generate split coverage visualizations across samples", "http://merenlab.org/2019/11/25/visualizing-coverages/#visualize-only-the-coverage-of-a-split-across-samples")]


@terminal.time_program
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
    run = terminal.Run()
    progress = terminal.Progress()

    utils.is_profile_db(args.profile_db)

    profile_db = dbops.ProfileSuperclass(args)
    profile_db_variant = utils.get_db_variant(args.profile_db)

    if args.list_splits:
        print('\n'.join(sorted(list(profile_db.split_names))))
        sys.exit(0)

    # further sanity checks
    if args.contigs_db:
        utils.is_contigs_db(args.contigs_db)

    # this extra flag is necessary for avoiding args.gene_caller_id == 0 which could evaluate to False
    # when gene caller id is 0 :)
    has_valid_gene_caller_id_arg = False
    if args.gene_caller_id != None:
        has_valid_gene_caller_id_arg = utils.is_gene_caller_id(args.gene_caller_id)
        if not args.contigs_db:
            raise ConfigError("If you wish to work with a gene caller id, you also need to provide "
                              "a contigs database (but you are off to a very good start).")

        if args.flank_length:
            try:
                assert(int(args.flank_length > 0))
            except:
                raise ConfigError("The flank length value should be a positive integer.")

    if args.split_name and args.contig_name:
        raise ConfigError("You have to declare either a split name or a contig name, but not both. RULES.")

    if args.flank_length and not args.gene_caller_id:
        raise ConfigError("The use of flank length is only relevant if you want to work with a gene :/")

    if not ((args.split_name or args.contig_name) or (args.collection_name or args.bin_id) or has_valid_gene_caller_id_arg):
        raise ConfigError("You must either declare a split or contig name, a bin name in a collection, or a gene caller id "
                          "for this to work.")

    if (args.collection_name or args.bin_id) and not (args.collection_name and args.bin_id):
        raise ConfigError("If you are declaring a collection name, you also need to declare a bin name, and vice "
                          "versa.")

    if args.collection_name and not args.contigs_db:
        raise ConfigError("If you want to work with a collection, you should also provide a contigs database. You "
                          "can tell it is necessary by the way it is like right there in the help menu.")

    if (args.split_name or args.contig_name) and (args.collection_name or args.bin_id):
        raise ConfigError("So you ask both for a split name and then ask for splits in a bin :/ You know who "
                          "doesn't like that? Anvi'o. Anvi'o doesn't like that.")

    if not args.output_file:
        raise ConfigError("You must declare an output file name. Because.")

    filesnpaths.is_output_file_writable(args.output_file)

    if not profile_db.auxiliary_profile_data_available:
        raise ConfigError("In order to get what you want from this program, you need the auxiliary "
                          "data file associated with this profile to be present. Now this program "
                          "will quit gracefully, and will let you imagine what might have gone "
                          "wrong.")

    splits_info_dict = {}
    target_splits_dict = {}
    if has_valid_gene_caller_id_arg:
        contigs_db = dbops.ContigsDatabase(args.contigs_db)
        matches = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_splits_table_name, where_clause=f"""gene_callers_id={args.gene_caller_id}""")

        if len(matches) > 1:
            gene_call = sorted(matches.values(), key=itemgetter('percentage_in_split'), reverse=True)[0]

            run.warning(f"Your gene call was split between {len(matches)} splits :/ So anvi'o will only use the one "
                        f"that contains the largest fraction of the gene sequence, the split {gene_call['split']}, "
                        f"which contains {gene_call['percentage_in_split']}% of the gene. We hope this doesn't ruin "
                        f"things for you (if it does, you will need to re-generate your contigs database with a larger "
                        f"split size because anvi'o programmers never imagined someone would actually get this warning).")
        else:
            gene_call = list(matches.values())[0]

        flank_length = int(args.flank_length) if args.flank_length else 0
        start = gene_call['start_in_split'] - flank_length if gene_call['start_in_split'] - flank_length > 0 else 0
        stop = gene_call['stop_in_split'] + flank_length

        if gene_call['split'] not in profile_db.split_names:
            raise ConfigError(f"The gene call '{args.gene_caller_id}' seems to be ocurring on a split ({gene_call['split']}) that is not "
                              f"in your profile database :/ This can happen since contigs that are shorter than the minimum contig length "
                              f"cutoff set during the anvi-profile step (which was {profile_db.p_meta['min_contig_length']} nts for this "
                              f"profile database) are discarded, and the coverages of the genes on them are not recovered. So it is all "
                              f"fine and it is just not a very good day for the gene call '{args.gene_caller_id}' :/")

        target_splits_dict = {gene_call['split']: {'start': start, 'stop': stop}}
    elif args.split_name:
        if args.split_name not in profile_db.split_names:
            raise ConfigError("'%s' does not seem to be a split in your profile database :/" % args.split_name)
        else:
            target_splits_dict[args.split_name] = {'start': 0, 'stop': sys.maxsize}
    elif args.contig_name:
        contigs_db = dbops.ContigsDatabase(args.contigs_db)
        splits_info_dict = contigs_db.db.get_table_as_dict(t.splits_info_table_name)

        contig_names = [tpl[0] for tpl in contigs_db.db.get_some_columns_from_table(t.contigs_info_table_name, 'contig')]

        if args.contig_name not in contig_names:
            raise ConfigError("'%s' does not seem to be a split in your profile database :/" % args.contig_name)

        split_names = [tpl[0] for tpl in contigs_db.db.get_some_columns_from_table(t.splits_info_table_name, 'split', where_clause=f'''parent="{args.contig_name}"''')]
        for split_name in split_names:
            target_splits_dict[split_name] = {'start': 0, 'stop': splits_info_dict[split_name]['length']}
    else:
        summary = summarizer.ProfileSummarizer(args, r=terminal.Run(verbose=False))
        summary.init()

        _bin = summarizer.Bin(summary, args.bin_id)
        for split_name in sorted(list(_bin.split_names)):
            target_splits_dict[split_name] = {'start': 0, 'stop': sys.maxsize}

        run.info('Collection', args.collection_name)
        run.info('Bin', args.bin_id)

    run.info('Num splits', len(target_splits_dict))

    sample_names = profile_db.p_meta['samples']
    run.info('Sample names', ', '.join(sample_names if len(sample_names) < 5 else (sample_names[:5] + ['(.. %d more ..)' % (len(sample_names) - 5)])))

    # report this stuff
    progress.new('Recovering coverages')
    entry_id = 0
    with open(args.output_file, 'w') as output_file:
        if args.contig_name:
            output_file.write('unique_entry_id\tnt_position\tcontig_name\tsample_name\tcoverage\n')
        elif has_valid_gene_caller_id_arg:
            output_file.write('unique_entry_id\tnt_position\tgene_caller_id\tsample_name\tcoverage\n')
        else:
            output_file.write('unique_entry_id\tnt_position\tsplit_name\tsample_name\tcoverage\n')

        split_names = list(target_splits_dict.keys())
        num_splits = len(split_names)

        # keep track of nucleotide positions for the contigs mode
        prior_nucleotides_in_contig = 0

        for i in range(num_splits):
            progress.update('processing split %d of %d ...' % (i, num_splits))

            split_name = split_names[i]
            start = target_splits_dict[split_name]['start']
            stop = target_splits_dict[split_name]['stop']

            if args.contig_name and i > 0:
                previous_split_name = split_names[i-1]
                prior_nucleotides_in_contig += target_splits_dict[previous_split_name]['stop']

            split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(profile_db.auxiliary_data_path, profile_db.p_meta['contigs_db_hash'], db_variant=profile_db_variant).get(split_name)

            num_nucleotides = len(list(split_coverage_values.values())[0])

            for sample_name in sorted(sample_names):
                for nt_pos in range(start, stop):
                    if nt_pos + 1 > num_nucleotides:
                        break

                    if args.contig_name:
                        output_file.write('%d\t%d\t%s\t%s\t%d\n' % (entry_id, nt_pos + prior_nucleotides_in_contig, args.contig_name, sample_name, split_coverage_values[sample_name][nt_pos]))
                    elif has_valid_gene_caller_id_arg:
                        output_file.write('%d\t%d\t%s\t%s\t%d\n' % (entry_id, nt_pos, args.gene_caller_id, sample_name, split_coverage_values[sample_name][nt_pos]))
                    else:
                        output_file.write('%d\t%d\t%s\t%s\t%d\n' % (entry_id, nt_pos, split_name, sample_name, split_coverage_values[sample_name][nt_pos]))

                    entry_id += 1

    progress.end()
    run.info('Output', "%d entries for %d splits in %d samples are stored in %s" % (entry_id, len(target_splits_dict), len(sample_names), args.output_file), mc='green')


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("ESSENTIAL ANVI'O DBs", "You need to provide a profile database, but whether you will need to provide a "
                                                               "contigs database will depend on which input option you want to go with.")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupA = parser.add_argument_group('INPUT OPTION #1: SPLIT or CONTIG NAME', "You want nothing but the coverage values in a single split .. or a contig. FINE.")
    groupA.add_argument(*anvio.A('split-name'), **anvio.K('split-name'))
    groupA.add_argument(*anvio.A('contig-name'), **anvio.K('contig-name'))

    groupB = parser.add_argument_group('INPUT OPTION #2: COLLECTION + BIN', "You want nucletide-level coverage values for all splits in a bin. FANCY.")
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupC = parser.add_argument_group('INPUT OPTION #3: GENE CALL', "You want nucletide-level coverage values for a given gene call. PRO.")
    groupC.add_argument(*anvio.A('gene-caller-id'), **anvio.K('gene-caller-id'))
    groupC.add_argument(*anvio.A('flank-length'), **anvio.K('flank-length'))

    groupD = parser.add_argument_group('BORING STUFF', "The output file and all.")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupD.add_argument(*anvio.A('list-splits'), **anvio.K('list-splits'))
    groupD.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupD.add_argument(*anvio.A('list-bins'), **anvio.K('list-bins'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
