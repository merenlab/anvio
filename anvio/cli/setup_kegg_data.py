#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio
import anvio.kegg as kegg

from anvio.terminal import time_program
from anvio.ttycolors import color_text as c
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'semiller10']
__provides__ = ["kegg-data", "modules-db"]
__description__ = "Download and setup various databases from KEGG"

## AVAILABLE DOWNLOAD MODES (and their parameters)
DOWNLOAD_MODES = {'KOfam': {'description': 'only KOfam annotation models (HMMs). Use this mode if '
                                           'you only want to run `anvi-run-kegg-kofams`.',
                            'arguments': {'only-download': {'flags': anvio.A('only-download'),
                                                           'definition': anvio.K('only-download')},
                                          'only-processing': {'flags': anvio.A('only-processing'),
                                                           'definition': anvio.K('only-processing')},
                                          'include-stray-KOs': {'flags': anvio.A('include-stray-KOs'),
                                                           'definition': anvio.K('include-stray-KOs')}
                                          }
                            },
                'modules': {'description': 'metabolic pathways from the KEGG MODULES database and BRITE hierarchies. Use this mode AND "KOfam" '
                                           'mode if you want to run pathway prediction with `anvi-estimate-metabolism`. This mode does the '
                                           'necessary setup that allows you to run `anvi-reaction-network` and visualize KEGG pathway maps.',
                            'arguments': {'only-download': {'flags': anvio.A('only-download'),
                                                            'definition': anvio.K('only-download')},
                                          'only-processing': {'flags': anvio.A('only-processing'),
                                                            'definition': anvio.K('only-processing')},
                                          'overwrite-output-destinations': {'flags': anvio.A('overwrite-output-destinations'),
                                                            'definition': anvio.K('overwrite-output-destinations',
                                                                {'help': "Overwrite any existing modules database "
                                                                            "in the KEGG data directory "
                                                                            "[USE WITH CAUTION]. Only relevant if you "
                                                                            "are using the --only-processing flag"})},
                                          'skip-brite-hierarchies': {'flags': anvio.A('skip-brite-hierarchies'),
                                                            'definition': anvio.K('skip-brite-hierarchies')},
                                          'skip-binary-relations': {'flags': anvio.A('skip-binary-relations'),
                                                            'definition': anvio.K('skip-binary-relations')},
                                          'skip-map-images': {'flags': anvio.A('skip-map-images'),
                                                            'definition': anvio.K('skip-map-images')}
                                        }
                            },
                'all': {'description': 'Download ALL KEGG data. This is the default mode.',
                                 'arguments': {'kegg-snapshot': {'flags': anvio.A('kegg-snapshot'),
                                                            'definition': anvio.K('kegg-snapshot')},
                                               'download-from-kegg': {'flags': anvio.A('download-from-kegg'),
                                                            'definition': anvio.K('download-from-kegg')},
                                               'kegg-archive': {'flags': anvio.A('kegg-archive'),
                                                            'definition': anvio.K('kegg-archive')}
                                            },
                }
}


@time_program
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
    args, unknown_args, subparsers = get_args()

    if args.list_modes:
        import anvio.terminal as terminal
        run = terminal.Run()
        run.warning(None, header="AVAILABLE DOWNLOAD MODES", lc="green")
        for mode, info_dict in DOWNLOAD_MODES.items():
            run.info(mode, info_dict['description'])
        sys.exit(0)

    # Here we parse mode-specific parameters that aren't recognized by the parent parser
    mode = args.mode
    mode_args, mode_unknown = subparsers[mode].parse_known_args(unknown_args)

    # if mode is 'all', we can also accept the mode-specific params
    if mode == "all":
        for sub in DOWNLOAD_MODES.keys():
            if sub != "all":
                submode_args, submode_unknown = subparsers[sub].parse_known_args(unknown_args)
                new_args_dict = vars(mode_args).copy()
                new_args_dict.update(vars(submode_args))
                mode_args = argparse.Namespace(**new_args_dict)
                for p in new_args_dict:
                    p_flag = "--" + p.replace("_", "-")
                    if p_flag in mode_unknown:
                        mode_unknown.remove(p_flag)

    args = argparse.Namespace(**vars(args), **vars(mode_args))
    # global flags are already handled by anvi'o and shouldn't be in the unknown list
    # these are coming from __init__.py and not all are relevant to this code but we catch them anyway
    global_flags_to_catch = ['--debug', '--no-progress', '--force', '--quiet', '--as-markdown',
                             '--force-overwrite', '--fix-sad-tables', '--display-db-calls', '--I-know-this-is-not-a-good-idea',
                             '--force-use-my-tree', '--debug-auto-fill-anvio-dbs']
    for gf in global_flags_to_catch:
        if gf in mode_unknown:
            mode_unknown.remove(gf)
    if len(mode_unknown):
        raise ConfigError(f"Unrecognized parameters: {' '.join(mode_unknown)}. Did you perhaps fail to specify the right mode?")

    if mode == "all" and not args.download_from_kegg:
        setup = kegg.KeggSetup(args)
        setup.setup_all_data_from_archive_or_snapshot()
    else:
        if mode == "KOfam" or mode == "all":
            args.download_from_kegg = True
            setup = kegg.KOfamDownload(args)
            setup.setup_kofams()
        if mode == "modules" or mode == "all":
            # do not reset the directory if it already happened
            if mode == "all" and args.reset:
                args.reset = False
                args.skip_init = True
            args.download_from_kegg = True
            setup = kegg.ModulesDownload(args)
            setup.setup_modules_data()


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    show_help = ('--help' in sys.argv) or ('-h' in sys.argv)

    groupM = parser.add_argument_group('MODE', "Select which data you want to download.")
    mode_help = "Depending on your choice here, this program will download and set up " + \
                "certain subsets of the data available from KEGG. Use --list-modes to see " + \
                f"a description of the options. Available modes: {', '.join(DOWNLOAD_MODES.keys())}"
    groupM.add_argument('--mode', choices=DOWNLOAD_MODES.keys(), help=mode_help, default='all')
    groupM.add_argument('--list-modes', **{'default': False, 'action': 'store_true',
                                         'help': "List the available modes and their descriptions."})

    # common arguments
    groupE = parser.add_argument_group('COMMON PARAMETERS', "These parameters apply to any mode.")
    groupE.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir'))
    groupE.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupE.add_argument(*anvio.A('reset'), **anvio.K('reset'))
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    groupA = parser.add_argument_group('MODE-SPECIFIC PARAMS', "Each section (underneath the program details) "
                                                               "below lists the parameters for one mode.")


    if show_help:
        parser.print_help()


    subparsers = {}
    for mode, info_dict in DOWNLOAD_MODES.items():
        subparser = argparse.ArgumentParser(usage=argparse.SUPPRESS, add_help=False)

        subparser._optionals.title = " \n%s\n%s" % (c(mode.upper(), "green"), ':' * 79)
        for arg_name, arg_dict in info_dict['arguments'].items():
            subparser.add_argument(*arg_dict['flags'], **arg_dict['definition'])

        if show_help:
            subparser.print_help()

        subparsers[mode] = subparser


    if show_help:
        sys.exit()

    args, unknown_args = parser.parse_known_args()

    return (args, unknown_args, subparsers)


if __name__ == '__main__':
    main()
