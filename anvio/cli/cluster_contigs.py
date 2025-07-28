#!/usr/bin/env python
# -*- coding: utf-8
"""A script to run automatic binning algorithms on a merged anvi'o profile"""
import os
import sys
import shutil
import random
import argparse

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers import driver_modules
from anvio.ttycolors import color_text as c
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.collections import TablesForCollections


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ["Christopher Quince"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ozcan', 'meren']
__resources__ = []
__tags__ = ["profile_db", "clustering", "collections"]
__requires__ = ['profile-db', 'contigs-db', 'collection']
__provides__ = ['collection', 'bin']
__description__ = "A program to cluster items in a merged anvi'o profile using automatic binning algorithms"


def store_clusters_in_db(contigs_db_path, profile_db_path, clusters, collection_name, source, cluster_type):

    if cluster_type == 'contig':
        contigs_db = dbops.ContigsDatabase(contigs_db_path)
        splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name, string_the_key = True)

        parent_to_split_map = {}
        for split in splits_basic_info:
            parent = splits_basic_info[split]['parent']

            if parent not in parent_to_split_map:
                parent_to_split_map[parent] = []

            parent_to_split_map[parent].append(split)

        cluster_splits = {}
        for bin_id in clusters:
            cluster_splits[bin_id] = []
            for contig_name in clusters[bin_id]:
                cluster_splits[bin_id] += parent_to_split_map[contig_name]

        clusters = cluster_splits

    bin_info_dict = {}

    for bin_id in clusters:
        bin_info_dict[bin_id] = {'html_color': '#' + ''.join(['%02X' % random.randint(50, 230) for i in range(0, 3)]), 'source': source}
                                                                # ^
                                                                #  \
                                                                #    poor man's random color generator

    c = TablesForCollections(profile_db_path, run=terminal.Run(verbose=False))
    c.append(collection_name, clusters, bin_info_dict)


def prepare_input_files(temp_path, profile_db, contigs_db):
    """This function generates files that we send to clustering algorithms.
    Returns input files inside a Namespace, which contains attributes:
        - contigs_db:
        - profile_db:
        .....
    """
    sample_names = sorted(list(profile_db.samples))
    input_files = argparse.Namespace()

    P = lambda x: os.path.join(temp_path, x)

    input_files = argparse.Namespace()
    input_files.contigs_db = os.path.realpath(contigs_db.db.db_path)
    input_files.profile_db = os.path.realpath(profile_db.db.db_path)

    input_files.contig_coverages_log_norm = P('contig_coverages_log_norm.txt')
    input_files.contig_coverages = P('contig_coverages.txt')

    input_files.split_coverages_log_norm = P('split_coverages_log_norm.txt')
    input_files.split_coverages = P('split_coverages.txt')

    input_files.splits_fasta = P('sequence_splits.fa')
    input_files.contigs_fasta = P('sequence_contigs.fa')

    splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)

    split_coverages, _ = profile_db.db.get_view_data('mean_coverage_contigs', splits_basic_info=splits_basic_info)
    split_coverages_log_norm, _ = profile_db.db.get_view_data('mean_coverage_contigs', splits_basic_info=splits_basic_info, log_norm_numeric_values=True)

    contig_coverages = {}
    contig_coverages_log_norm = {}
    contig_names = set([x['__parent__'] for x in split_coverages.values()])

    for split_name in split_coverages:
        entry = split_coverages[split_name]
        entry_log_norm = split_coverages_log_norm[split_name]
        c = entry['__parent__']
        if c in contig_names and c not in contig_coverages:
            contig_coverages[c] = entry
            contig_coverages_log_norm[c] = entry_log_norm

    # Write output files
    utils.store_dict_as_TAB_delimited_file(split_coverages, input_files.split_coverages, ['contig', *sample_names])
    utils.store_dict_as_TAB_delimited_file(split_coverages_log_norm, input_files.split_coverages_log_norm, ['contig', *sample_names])

    utils.store_dict_as_TAB_delimited_file(contig_coverages, input_files.contig_coverages, ['contig', *sample_names])
    utils.store_dict_as_TAB_delimited_file(contig_coverages_log_norm, input_files.contig_coverages_log_norm, ['contig', *sample_names])

    utils.export_sequences_from_contigs_db(contigs_db.db_path,
                                           input_files.splits_fasta,
                                           splits_mode=True)

    utils.export_sequences_from_contigs_db(contigs_db.db_path,
                                           input_files.contigs_fasta,
                                           splits_mode=False)

    return input_files


@terminal.time_program
def main():
    args, unknown, subparsers, modules = get_args()

    try:
        cluster_contigs(args, unknown, subparsers, modules)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parent_parser = argparse.ArgumentParser(description=__description__)
    parent_parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    parent_parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parent_parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name', {'required': True}))
    parent_parser.add_argument(*anvio.A('driver'), **anvio.K('driver',
        {'choices': list(driver_modules['binning'].keys()), 'type': str.lower}))
    parent_parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parent_parser.add_argument(*anvio.A('log-file'), **anvio.K('log-file'))
    parent_parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    show_help = ('--help' in sys.argv) or ('-h' in sys.argv)

    if show_help:
        parent_parser.print_help()

    modules = driver_modules['binning']
    subparsers = {}
    for name, module in modules.items():
        subparser = argparse.ArgumentParser(usage=argparse.SUPPRESS, add_help=False)

        try:
            # test if __init__() will raise ConfigError with is_program_exists
            module()

            subparser._optionals.title = " \n%s\n%s" % (c(name.upper(), "green"), ':' * 79)
            for arg_name, definiton in module.arguments.items():
                subparser.add_argument(*definiton[0], **definiton[1])

            if show_help:
                subparser.print_help()
        except ConfigError:
            if show_help:
                print("\n%s %s" % (c(name.upper(), "green"), c("[NOT FOUND]", "red")))

            continue
        finally:
            subparsers[name] = subparser

    if show_help:
        sys.exit()

    args, unknown = parent_parser.parse_known_args()

    return (args, unknown, subparsers, modules)


def cluster_contigs(args, unknown, subparsers, modules):
    run = terminal.Run()

    if not args.just_do_it:
        raise ConfigError("This is an experimental module and is not yet fully tested. Alternatively you could "
                          "use external binning algorithms with raw data or data exported from anvi'o (i.e., see program "
                          "`anvi-export-splits-and-coverages`) and import your clusters into anvi'o with `anvi-import-collection`. "
                          "If you still wish to try this workflow on your own risk, please run `anvi-cluster-contigs` with "
                          "the flag `--just-do-it`. Anvi'o developers thank you for your understanding and help.")
    else:
        run.warning("You are running an experimental workflow not every part of which may be fully and thoroughly tested :) "
                    "Please scrutinize your output carefully after analysis, and keep us posted if you see things that "
                    "surprise you.")

    collection_name = args.collection_name.strip()
    if not len(collection_name):
        raise ConfigError('Nice try. Collection name cannot be emtpy')
    try:
        utils.check_sample_id(collection_name)
    except:
        raise ConfigError('"%s" is not a proper collection name. A proper one should be a single word and not contain '
                           'ANY characters but digits, ASCII letters and underscore character(s). There should not be '
                           'any space characters, and the collection name should not start with a digit.' % collection_name)

    if args.log_file:
        filesnpaths.is_output_file_writable(args.log_file)

    utils.is_profile_db_and_contigs_db_compatible(args.profile_db, args.contigs_db)

    merged_profile_db = dbops.ProfileDatabase(args.profile_db)

    if(merged_profile_db.meta['merged'] != True):
        raise ConfigError("'%s' does not seem to be a merged profile database :/" % args.profile_db)

    contigs_db = dbops.ContigsDatabase(args.contigs_db)

    # Here we parse not sys.argv but unknown parameters
    # that is not recognized by parent parser.
    driver = args.driver
    sub_args, sub_unknown = subparsers[driver].parse_known_args(unknown)

    # --debug handled by anvio/__init__.py can be in args
    if '--debug' in sub_unknown:
        sub_unknown.remove('--debug')

    if len(sub_unknown):
        raise ConfigError("Unrecognized parameters %s" % ' '.join(sub_unknown))

    binning_module = modules[driver]()
    binning_module_name = type(binning_module).__name__

    cluster_type = getattr(binning_module, 'cluster_type', False)
    if not cluster_type or cluster_type not in ['contig', 'split']:
        raise ConfigError("Anvi'o doesn't know what type of cluster this module returns. "
                         "Please define an 'cluster_type' class variable containing either "
                         "'contig' or 'split'.")

    working_dir = filesnpaths.get_temp_directory_path()


    run.info("Contigs DB", args.contigs_db)
    run.info("Profile DB", args.profile_db)
    run.info("Binning module", binning_module_name)
    run.info("Cluster type", cluster_type)
    run.info('Working directory', working_dir)

    input_files = prepare_input_files(working_dir, merged_profile_db, contigs_db)

    if anvio.DEBUG:
        for name, path in input_files.__dict__.items():
            run.info('Input file [%s]' % name, path)

    citation = getattr(binning_module, 'citation', False)
    if citation:
        run.warning("Anvi'o is now passing all your data to the binning module '%s'. If you publish results "
                    "from this workflow, please do not forget to reference the following citation." % binning_module_name,
                     header='CITATION', lc='green')
        run.info_single("%s" % binning_module.citation, nl_after=1, mc='green')
    else:
        raise ConfigError("The module %s does not seem to have any citation information. If you are currently developing a new "
                          "module please address this in your code. If you are a user, please consider dropping a line to anvi'o "
                          "developers." % (binning_module_name))

    clusters = binning_module.cluster(input_files, sub_args, working_dir, threads=args.num_threads, log_file_path=args.log_file)

    if not anvio.DEBUG:
        shutil.rmtree(working_dir)
    else:
        run.warning("Working in debug mode, Anvi'o will not clean the temporary directories.")

    if not len(clusters):
        c = TablesForCollections(args.profile_db)
        c.delete(collection_name)

        raise ConfigError("Binning algorithm '%s' did not return any clusters :( There will not be any collection "
                          "in your profie database for these results. Please take a look at the help menu and sanity "
                          "check your parameters for this binning algorithm." % (driver))

    store_clusters_in_db(args.contigs_db, args.profile_db, clusters, collection_name, driver, cluster_type)

    run.info_single("%s formed %d clusters, which are being added to the database as a collection named %s." % \
                            (driver, len(clusters), collection_name), nl_before=1, nl_after=1, mc="green")


if __name__ == '__main__':
    main()
