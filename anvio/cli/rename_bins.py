#!/usr/bin/env python
# -*- coding: utf-8

import sys
import copy
import operator

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.completeness import Completeness
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.collections import TablesForCollections


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["collection", "bin"]
__requires__ = ["collection", "bin", "profile-db", "contigs-db"]
__description__ = "Rename all bins in a given collection (so they have pretty names)"


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

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    profile_db_path = A('profile_db')
    contigs_db_path = A('contigs_db')
    list_collections = A('list_collections')
    collection_name_to_read = A('collection_to_read')
    collection_name_to_write = A('collection_to_write')
    prefix = A('prefix')
    report_file_path = A('report_file')
    call_MAGs = A('call_MAGs')
    exclude_bins = A('exclude_bins')
    min_completion_for_MAG = A('min_completion_for_MAG')
    max_redundancy_for_MAG = A('max_redundancy_for_MAG')
    size_for_MAG = A('size_for_MAG')
    dry_run = A('dry_run')

    utils.is_profile_db_and_contigs_db_compatible(profile_db_path, contigs_db_path)

    collections = ccollections.Collections()
    collections.populate_collections_dict(profile_db_path)

    if list_collections:
        collections.list_collections()
        sys.exit()

    if not prefix:
        raise ConfigError("Anvi'o is having hard time believing that you called this function without "
                           "a prefix to rename bins in the collection '%s'." % collection_name_to_read)

    utils.is_this_name_OK_for_database('prefix for bins', prefix)

    if not report_file_path:
        raise ConfigError("You must provide an output file name to report file changes. It may or may not "
                           "be useful to you, but let's don't take unnecessary risks, eh? (you can use the "
                           "`--report-file` parameter)")

    filesnpaths.is_output_file_writable(report_file_path)

    if exclude_bins and not call_MAGs:
        raise ConfigError("You can use the flag `--exclude-bins` only if you are also calling MAGs "
                          "(`--call-MAGs`). If this is not clear to you why it is the case, please "
                          "take a look at the help menu (and hope for the best).")

    if min_completion_for_MAG < 0:
        raise ConfigError("Minimum completion value can't be negative because that's not how math works.")

    if min_completion_for_MAG > 100:
        raise ConfigError("A minimum completion score for a MAG that is over 100%? Anvi'o is impressed with "
                          "your expectations from your MAGs, but it cannot accept anything above 100% for "
                          "minimum completion :(")

    if max_redundancy_for_MAG < 0:
        raise ConfigError("You must choose a maximum redundancy value that is at least 0.")

    if collection_name_to_read == collection_name_to_write:
        raise ConfigError("Well :( Collection name to read from can't be identical to the collection name "
                           "you want to use to store updated bin names. You know. It kinda defeats the purpose.")

    if not collection_name_to_read:
        raise ConfigError("You must provide a collection name to read from.")

    if not collection_name_to_write:
        raise ConfigError("You must provide a new collection name to write updated bin names.")

    utils.is_this_name_OK_for_database('collection name two write', collection_name_to_write)

    if  collection_name_to_read not in collections.collections_dict:
        raise ConfigError("%s is not a valid collection name, because it doesn't exist :(. See a "
                           "list of available ones with '--list-collections' flag" % collection_name_to_read)

    if  collection_name_to_write in collections.collections_dict:
        raise ConfigError("The new collection name %s is already in the database. You must choose a new "
                           "collection name." % collection_name_to_write)

    completeness = Completeness(contigs_db_path)

    if not len(completeness.sources):
        raise ConfigError("HMM's were not run for this contigs database :/ Without that, how can this poor program can rename bins based on "
                           "their completion estimates? :(")

    run.warning('This run will use the completion and redundancy estimate '
                    'recovered from the single-copy core gene collection that provides the highest '
                    'completion estimate.')

    if not len(completeness.sources):
        raise ConfigError("It seems the contigs database does not contain HMM hits for "
                           "any of the single copy core gene collections :/ Bad news.")


    contigs_db = dbops.ContigsSuperclass(args, r = run, p = progress)

    collection_dict = collections.get_collection_dict(collection_name_to_read)
    bins_info_dict = collections.get_bins_info_dict(collection_name_to_read)

    if call_MAGs:
        warning_msg = "Please read carefully. This program is set to call any bin a 'MAG', if the bin is estimated to be\
                       less than %.0f%% redundant and more than %.0f%% complete. " % (max_redundancy_for_MAG, min_completion_for_MAG)

        if size_for_MAG:
            warning_msg += "But since you set a value for MAG size, any bin that meets the redundancy criterion and larger\
                            than %.2f Mbp in size, WILL ALSO be called a MAG regardless of its completion estimate. This is\
                            a way to not miss refined bins with low completion estimates but large sizes. If this is not what\
                            you want, set `--size-for-MAG` to 0." % (size_for_MAG)
        else:
            warning_msg += "This process will not pay attention to the bin size (Mbp). If you think it should, please see the\
                            flag `--size-for-MAG` in the help menu."

        run.warning(warning_msg, "A FRIENDLY REMINDER", lc='green')

    counter = 0
    MAGs_sorted_by_completion = []
    bins_sorted_by_completion = []
    total_num_bins = len(collection_dict)
    progress.new('Going through bins in the collection %s' % collection_name_to_read)
    for bin_name in collection_dict:
        counter += 1
        progress.update('%d in %d' % (counter, total_num_bins))
        p_completion, p_redundancy, domain, domain_probabilities, info_text, d = completeness.get_info_for_splits(set(collection_dict[bin_name]))

        if domain == 'blank':
            # if the domain is 'blank', it most likely means that our random forest classifier
            # was not certain about the domain because the genome bin wish highly incomplete
            # for a reliable prediction. in that case we will use the most reasonable C/R
            # estimate based on the highest C and set the domain to that for reporting.
            l = []
            for domain in d:
                for scg_collection_name in d[domain]:
                    l.append((d[domain][scg_collection_name]['percent_completion'], domain, scg_collection_name),)
            _, domain, scg_collection_name = sorted(l, reverse=True)[0]

            if d[domain][scg_collection_name]['percent_completion'] > 0.0:
                p_completion = d[domain][scg_collection_name]['percent_completion']
                p_redundancy = d[domain][scg_collection_name]['percent_redundancy']
            else:
                domain = 'blank'

        size_in_Mbp = sum([contigs_db.splits_basic_info[split_name]['length'] for split_name in set(collection_dict[bin_name])]) / 1000000.0
        substantive_completion = p_completion - p_redundancy

        if call_MAGs:
            if p_redundancy < max_redundancy_for_MAG:
                if p_completion >= min_completion_for_MAG:
                    MAGs_sorted_by_completion.append((bin_name, domain, substantive_completion, p_completion, p_redundancy, size_in_Mbp, 'MAG'),)
                elif size_for_MAG and size_in_Mbp >= size_for_MAG:
                    MAGs_sorted_by_completion.append((bin_name, domain, substantive_completion, p_completion, p_redundancy, size_in_Mbp, 'MAG'),)
                else:
                    bins_sorted_by_completion.append((bin_name, domain, substantive_completion, p_completion, p_redundancy, size_in_Mbp, 'Bin'),)
            else:
                bins_sorted_by_completion.append((bin_name, domain, substantive_completion, p_completion, p_redundancy, size_in_Mbp, 'Bin'),)
        else:
            bins_sorted_by_completion.append((bin_name, domain, substantive_completion, p_completion, p_redundancy, size_in_Mbp, 'Bin'),)

    MAGs_sorted_by_completion.sort(key=operator.itemgetter(1), reverse=True)
    bins_sorted_by_completion.sort(key=operator.itemgetter(1), reverse=True)

    progress.update('Finalizing the report...')

    new_collection_dict = {}
    new_bins_info_dict = {}

    counter = 0
    report = open(report_file_path, 'w')
    report.write('old_bin_name\tnew_bin_name\tSCG_domain\tcompletion\tredundancy\tsize_in_Mbp\n')

    if exclude_bins:
        list_of_bins_in_the_workings = MAGs_sorted_by_completion
    else:
        list_of_bins_in_the_workings = MAGs_sorted_by_completion + bins_sorted_by_completion

    for bin_name, domain, substantive_completion, completion, redundancy, size_in_Mbp, bin_type in list_of_bins_in_the_workings:
        counter += 1
        new_bin_name = f"{prefix}_{bin_type}_{counter:05}"
        new_collection_dict[new_bin_name] = copy.deepcopy(collection_dict[bin_name])

        if bin_name in bins_info_dict:
            new_bins_info_dict[new_bin_name] = copy.deepcopy(bins_info_dict[bin_name])

        report.write(f"{bin_name}\t{new_bin_name}\t{domain}\t{completion:.2f}\t{redundancy:.2f}\t{size_in_Mbp}\n")

    report.close()

    progress.end()

    if call_MAGs and exclude_bins:
        if len(bins_sorted_by_completion) == 0:
            run.warning(f"You have used the flag `--exclude-bins`, so anvi'o was determined to not include "
                        f"any of the bins in your previous collection '{collection_name_to_read}' that did "
                        f"not match to your criteria to call MAGs in your new collection. Although everyting "
                        f"in your original collection ended up being called a MAG, so everything ended up "
                        f"being reported anyway.", header="EXCLUDE BINS FLAG WARNING", lc="yellow")
        elif len(bins_sorted_by_completion) == 1:
            run.warning(f"You have used the flag `--exclude-bins`, and anvi'o found a single bin in your "
                        f"original collection '{collection_name_to_read}' that did not match your `--call-MAGs` "
                        f"criteria. So you will not have the bin '{bins_sorted_by_completion[0][0]}' in your "
                        f"new collection.", header="EXCLUDE BINS FLAG WARNING", lc="yellow")
        else:
            bins_no_longer_with_us = ', '.join([f"'{e[0]}'" for e in bins_sorted_by_completion])
            run.warning(f"You have used the flag `--exclude-bins` to exclude bins that were not up to your "
                        f".. standards. As a result, anvi'o excluded {len(bins_sorted_by_completion)} bins "
                        f"found in your original collection '{collection_name_to_read}', and they will not "
                        f"appear in your new collection '{collection_name_to_write}'. Here is a list of their "
                        f"names, just so you take a moment to reflect on their journey which has come to an end "
                        f"in this this project (hunger games styla): {bins_no_longer_with_us}.",
                        header="EXCLUDE BINS FLAG WARNING", lc="yellow")

    run.info('Report file', '%s' % (report_file_path))

    if dry_run:
        run.warning('This was a dry run. So nothing is updated in the profile database, but now there is '
                    'a fancy report file that shows how your bins would have been renamed. Please take a '
                    'look at it and make sure things seem just like the way you want. Once you are satisfied '
                    'you can run the program without the --dry-run flag.')

        return

    collections_table = TablesForCollections(profile_db_path, run=terminal.Run(verbose=False))
    collections_table.append(collection_name_to_write, new_collection_dict, new_bins_info_dict)

    run.info_single(f"All {counter} bins in the collection '{collection_name_to_read}' is renamed, and "
                    f"stored under the collection name '{collection_name_to_write}'.", nl_before=1,
                    nl_after=1, mc="green")


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DEFAULT INPUTS', "Standard stuff")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupA.add_argument('--collection-to-read', default = None, help = "Collection name to read from. "
                            "Anvi'o will not overwrite an existing collection, instead, it will create "
                            "a copy of your collection with new bin names.")
    groupA.add_argument('--collection-to-write', default = None, help = "The new collection name. Give "
                            "it a nice, fancy name.")

    groupB = parser.add_argument_group('OUTPUT AND TESTING', "a.k.a, sweet parameters of convenience")
    groupB.add_argument('--prefix', default = None, help = "Prefix for the bin names. Must be a "
                            "single word, composed of digits and numbers. The use of the underscore "
                            "character is OK, but that's about it (fine, the use of the dash character "
                            "is OK, too but no more!). If the prefix is 'PREFIX', each bin will be "
                            "renamed as 'PREFIX_XXX_00001, PREFIX_XXX_00002', and so on, in the order "
                            "of percent completion minus percent redundancy (what we call, 'substantive "
                            "completion'). The 'XXX' part will either be 'Bin', or 'MAG depending on "
                            "other parameters you use. Keep reading.")
    groupB.add_argument('--report-file', metavar = 'REPORT_FILE_PATH', default = None, help = "This file "
                            "will report each name change event, so you can trace back the original names "
                            "of renamed bins later.")
    groupB.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupB.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run', {'help': "When used does NOT update the "
                            "profile database, just creates the report file so you can view how things "
                            "will be renamed."}))

    groupC = parser.add_argument_group('MAG OPTIONS', "If you want to call some bins 'MAGs' because you are so cool")
    groupC.add_argument('--call-MAGs', default=False, action="store_true", help = "This program by default "
                            "rename your bins as 'PREFIX_Bin_00001', 'PREFIX_Bin_00002' and so on. If you "
                            "use this flag, it will name the ones that meet the criteria described by "
                            "MAG-related flags as 'PREFIX_MAG_00001', 'PREFIX_MAG_00002', and so on. The "
                            "ones that do not get to be named as MAGs will remain as bins.")
    groupC.add_argument('--min-completion-for-MAG', default=70, type=float, help="If --call-MAGs flag is "
                            "used, call any bin a 'MAG' if their completion estimate is above this (the "
                            "default is %(default)d), and the redundancy estimate is less than "
                            "--max-redundancy-for-MAG.")
    groupC.add_argument('--max-redundancy-for-MAG', default=10, type=float, help="If --call-MAGs flag is "
                            "used, call any bin a 'MAG' if their redundancy estimate is below this (the "
                            "default is %(default)d) and the completion estimate is above --min-completion-for-MAG.")
    groupC.add_argument('--size-for-MAG', default=0.0, type=float, metavar='MEGABASEPAIRS', help="If --call-MAGs "
                            "flag is used, call any bin a 'MAG' if their redundancy estimate is less than "
                            "--max-redundancy-for-MAG, AND THEIR SIZE IS LARGER THAN THIS VALUE REGARDLESS OF THE "
                            "COMPLETION ESTIMATE. The default behavior is to not care about this at all.")
    groupC.add_argument('--exclude-bins', default=False, action="store_true", help="Exclude bins that don't match "
                            "any of your criteria to call MAGs. By default, your new collection will include "
                            "every single item in your original collection. But if you are calling MAGs with "
                            "the `--call-MAGs` flag AND would like to keep only those that are determined to "
                            "be a MAG, you can use this flag and all the 'bins' would then be excluded from "
                            "your new collection, leaving you with only the 'good' stuff.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
