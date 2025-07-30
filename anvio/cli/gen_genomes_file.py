#!/usr/bin/env python
# -*- coding: utf-8

"""A program to generate external or internal genomes files"""

import os
import sys
from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.ccollections as ccollections
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.dbinfo import FindAnvioDBs

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__authors__ = ['ivagljiva']
__provides__ = ["external-genomes", "internal-genomes"]
__requires__ = ["contigs-db", "profile-db", "collection"]
__description__ = "Generate an external genomes or internal genomes file"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(1)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    if not (args.input_dir or args.profile_db):
        raise ConfigError("You must provide either an input directory (for an external genomes file) or a "
                          "profile db/contigs db/collection name set (for an internal genomes file). Also, "
                          "just so you don't forget, we also need an output file path in either case.")

    if args.profile_db:
        if not (args.contigs_db and args.collection_name):
            raise ConfigError("Since you provided a profile db, we think you want an internal genomes file, but "
                              "anvi'o needs the corresponding contigs db and a collection name, as well.")
        if not args.output_file:
            raise ConfigError("Please provide an output file path for your internal genomes file using the -o parameter.")

        if filesnpaths.is_file_exists(args.output_file, dont_raise=True):
            raise ConfigError("The internal genomes file path that you provided already exists, and anvi'o will not overwrite "
                              "it. Please either remove the existing file or provide a different file path.")

        utils.is_pan_or_profile_db(args.profile_db)

        utils.is_profile_db_and_contigs_db_compatible(args.profile_db, args.contigs_db)

        progress.new('Accessing to the collections table')
        progress.update('...')
        collections = ccollections.Collections()
        collections.populate_collections_dict(args.profile_db)
        progress.end()

        if not collections.collections_dict:
            raise ConfigError("There are no collections in this profile database, so no internal genomes file. "
                              "Consider making some in interactive mode, or importing a collection with "
                              "`anvi-import-collection` (run `anvi-import-collection --help` for links to helpful "
                              "resources).")

        if args.collection_name not in collections.collections_dict:
            raise ConfigError(f"The collection you requested, {args.collection_name}, does not appear to be a "
                              f"valid collection in this profile database.")

        contig_db_path = os.path.abspath(args.contigs_db)
        profile_db_path = os.path.abspath(args.profile_db)

        int_genomes_dict = {}
        for bin_name in sorted(collections.collections_dict[args.collection_name]['bin_names'].split(',')):
            int_genomes_dict[bin_name] = {}
            int_genomes_dict[bin_name]['bin_id'] = bin_name
            int_genomes_dict[bin_name]['collection_id'] = args.collection_name
            int_genomes_dict[bin_name]['profile_db_path'] = profile_db_path
            int_genomes_dict[bin_name]['contigs_db_path'] = contig_db_path

        utils.store_dict_as_TAB_delimited_file(int_genomes_dict, args.output_file, key_header='name')

        run.info("Internal genomes file", args.output_file)

    elif args.input_dir:
        if not args.output_file:
            raise ConfigError("Please provide an output file path for your external genomes file using the -o parameter.")

        if not filesnpaths.is_file_exists(args.input_dir, dont_raise=True):
            raise ConfigError("The directory you provided does not exist :/")

        if filesnpaths.is_file_exists(args.output_file, dont_raise=True):
            raise ConfigError("The external genomes file path that you provided already exists, and anvi'o will not overwrite "
                              "it. Please either remove the existing file or provide a different file path.")

        anvio_dbs = FindAnvioDBs(args.input_dir, depth=(1 if not args.include_subdirs else 5)).anvio_dbs
        contigs_dbs = anvio_dbs['contigs'] if 'contigs' in anvio_dbs else []

        if not len(contigs_dbs):
            msg = "No contigs databases were found in the input directory :/ "
            if not args.include_subdirs:
                msg += ("If you would like to include subdirectories in your search as well (with a max depth of "
                        "5) you can re-run the same script with the flag `--include-subdirs`.")

            raise ConfigError(msg)
        else:
            run.info("Contigs databases", f"{len(contigs_dbs)} found.", mc="green")

        # check uniqueness of project names in contigs-db files to avoid masking files
        # see https://github.com/merenlab/anvio/issues/2239 for details
        contigs_db_project_names = [contigs_db.project_name for contigs_db in contigs_dbs]
        if len(set(contigs_db_project_names)) != len(contigs_db_project_names):
            f = Counter(contigs_db_project_names).most_common(1)
            raise ConfigError(f"Some of the contigs-db files seem to have identical project names :/ For instance, "
                              f"the project name '{f[0][0]}' occurs {f[0][1]} times among the set of {len(contigs_dbs)} "
                              f"contigs-db files this tool is currently processing. For each contigs-db, a 'project name' "
                              f"is assigned automatically based on its file name at the time of creation. A much better "
                              f"practice is to use the `--project-name` flag with `anvi-gen-contigs-database` to "
                              f"be explcit and set set a unique project name variable for each contigs-db.")

        ext_genomes_dict = {}
        for contigs_db in contigs_dbs:
            ext_genomes_dict[contigs_db.project_name] = {'contigs_db_path': os.path.abspath(contigs_db.path)}

        utils.store_dict_as_TAB_delimited_file(ext_genomes_dict, args.output_file, key_header='name')

        run.info("External genomes file", args.output_file)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("EXTERNAL GENOMES", "Provide a directory, and anvi'o will provide an external genomes "
                                                           "file containing all contigs dbs in that directory.")
    groupA.add_argument(*anvio.A('input-dir'), **anvio.K('input-dir'))
    groupA.add_argument(*anvio.A('include-subdirs'), **anvio.K('include-subdirs'))

    groupB = parser.add_argument_group("INTERNAL GENOMES", "Provide a contigs db, profile db, and collection name and anvi'o "
                                                           "will bestow upon you an internal genomes file for that collection.")
    groupB.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupB.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))


    groupC = parser.add_argument_group("OUTPUT", "Path for your internal or external genomes file")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
