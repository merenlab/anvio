#!/usr/bin/env python
# -*- coding: utf-8

import sys
import shutil

import anvio
import anvio.filesnpaths as filesnpaths

from anvio.drivers.emapper import EggNOGMapper
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Run eggnog-mapper on a contigs database, and store results"


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

    if args.annotation:
        if not args.use_version:
            raise ConfigError("If you would like to use a pre-existing annotation file to import from,\
                                you must also declare a version.")

        args.num_threads = 1
        emapper = EggNOGMapper(args, use_version=args.use_version)
        emapper.populate_annotations_dict(args.annotation)
        emapper.store_annotations_in_db(drop_previous_annotations=args.drop_previous_annotations)
    else:
        raise ConfigError("Although the codebase is ready on the anvi'o side (at least when this note was written on "
                           "December 2016, the eggnog-mapper will not run properly due to some technical issues beyond "
                           "anvi'o. Although you can still run it by hand, and then use this script with the resulting "
                           "output file usin gthe '--annotation' parameter. That way you can still import EggNOG "
                           "annotations into your contigs database. Meanwhile we will try to address this problem, and "
                           "hopoefully remove this note.")
        output_dir = filesnpaths.get_temp_directory_path()

        emapper = EggNOGMapper(args)
        emapper.process(output_dir, drop_previous_annotations=args.drop_previous_annotations)

        if not anvio.DEBUG:
            shutil.rmtree(output_dir)


def get_args():
    # quickly learn available parser versions
    from argparse import Namespace
    known_versions = EggNOGMapper(args=Namespace()).available_parsers.keys()

    # construct the help and the args
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    parser.add_argument('--drop-previous-annotations', default=False, action="store_true",
                        help='When declared, previous annotations in the database will be dropped.')
    parser.add_argument('--annotation', default=None, metavar='EMAPPER_ANNOTATION_FILE',
                        help="If you have an annotation file from a previous run, you can call this program to\
                              import the contents of that file into the database instead of a run from scratch.\
                              In that case, you must also use the `--use-version` parameter to clarify which\
                              parser version should be used to parse it.")
    parser.add_argument('--use-version', default=None, metavar="EMAPPER_VERSION",
                        help=f'The version of eggnog-mapper that generated the annotation file. Known versions \
                               of emapper include the following: {", ".join(known_versions)}')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
