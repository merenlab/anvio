#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

with terminal.SuppressAllOutput():
    import anvio.data.hmm as hmm_data

available_hmm_sources = list(hmm_data.sources.keys())

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program
from anvio.tables.hmmhits import TablesForHMMHits
from anvio.tables.trnahits import TablesForTransferRNAs


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ivagljiva', 'mschecht', 'ekiefl']
__requires__ = ['contigs-db', 'hmm-source',]
__provides__ = ['hmm-hits',]
__description__ = ("This program deals with populating tables that store HMM hits in an "
                   "anvi'o contigs database")
__resources__ = [("Another description as part of the metagenomic workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile")]


@time_program
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
    P = terminal.pluralize

    # then check whether we are going to use the default HMM profiles, or run it for a new one.
    sources = {}
    if args.hmm_profile_dir and args.installed_hmm_profile:
        raise ConfigError("You must select either an installed HMM profile, or provide a path for a profile directory. "
                          "But not both :/")
    elif args.hmm_profile_dir:
        if not os.path.exists(args.hmm_profile_dir):
            raise ConfigError('No such file or directory: "%s"' % args.hmm_profile_dir)
        sources = utils.get_HMM_sources_dictionary([args.hmm_profile_dir])
        run.info('HMM profiles', '%d source%s been loaded: %s' % (len(sources),
                                                          's' if len(sources) > 1 else '',
                                                          ', '.join(['%s (%d genes)' % (s, len(sources[s]['genes']))\
                                                                                                    for s in sources])))
    elif args.installed_hmm_profile:
        args.installed_hmm_profile = [p.strip() for p in args.installed_hmm_profile.split(',') if p.strip()]

        if "Transfer_RNAs" in args.installed_hmm_profile:
            raise ConfigError("Sorry, it is not you, it is anvi'o :/ Transfer RNAs behave identically to HMMs, but "
                              "to run transfer RNA profiles on a contigs database you must use the program "
                              "`anvi-scan-trnas` and not provide it as an installed HMM source to `anvi-run-hmms`.")

        missing_installed_hmm_profiles = [p for p in args.installed_hmm_profile if p not in available_hmm_sources]
        if len(missing_installed_hmm_profiles):
            n = len(missing_installed_hmm_profiles)
            m = ', '.join([f'"{p}"' for p in missing_installed_hmm_profiles])
            raise ConfigError(f"You managed to hit {P('profile', n)} that {P('is', n, alt='are')} not actually installed "
                              f"({m}) :/ Here are the available profiles: {', '.join(available_hmm_sources)}.")

        sources = {}
        for profile in args.installed_hmm_profile:
            sources[profile] = hmm_data.sources[profile]
    else:
        # sources will be loaded from defaults.
        pass

    # storing hmmer output only works one profile at a time since the output will be overwritten with each profile
    if args.hmmer_output_dir and not len(sources.keys()) == 1:
        raise ConfigError("You requested to save the HMMER output files in a directory, but this will not work with "
                          "more than one HMM source because the output files would be overwritten with each new source. "
                          "We are sorry for this, but we humbly request that you run on each source one at a time if you "
                          "want to keep the HMMER output.")

    if args.domain_hits_table and not args.hmmer_output_dir:
        raise ConfigError("We can see that you have requested --domain-hits-table but you haven't asked us to store "
                          "this output in a directory with --hmmer-output-dir. There is no point to requesting this output "
                          "if you are never going to see it, so we figured we'd stop you right there. :)")

    search_tables = TablesForHMMHits(args.contigs_db, num_threads_to_use=args.num_threads, just_do_it=args.just_do_it,
                                     hmm_program_to_use=args.hmmer_program, hmmer_output_directory=args.hmmer_output_dir,
                                     get_domain_table_output=args.domain_hits_table, add_to_functions_table=args.add_to_functions_table)
    search_tables.populate_search_tables(sources)

    if not args.hmm_profile_dir and not args.installed_hmm_profile and args.also_scan_trnas:
        tables_for_trna_hits = TablesForTransferRNAs(args)
        tables_for_trna_hits.populate_search_tables(args.contigs_db)


def get_args():
    available_hmm_sources_pretty = '; '.join([f"'{s}' (type: {hmm_data.sources[s]['kind']})" for s in available_hmm_sources])

    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)


    groupA = parser.add_argument_group("DB", "An anvi'o contigs database to populate with HMM hits")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group("HMM OPTIONS", "If you have your own HMMs, or if you would like to run only a set of "
                                                      "default anvi'o HMM profiles rather than running them all, this is "
                                                      "your stop.")
    groupB.add_argument(*anvio.A('hmm-profile-dir'), **anvio.K('hmm-profile-dir'))
    groupB.add_argument(*anvio.A('installed-hmm-profile'), **anvio.K('installed-hmm-profile',
                            {'help': f"When you run this program without any parameter, it will run all {len(available_hmm_sources)} "
                                     f"HMM profiles installed on your system. Using this parameter, you can instruct anvi'o to run "
                                     f"only one or more of the specific profiles of you choose. You can provide a comma-separated list "
                                     f"of names for multiple profiles (but in that case don't put a space between each profile name). "
                                     f"Here is the list of installed profiles available to you: "
                                     f"{available_hmm_sources_pretty}."}))
    groupB.add_argument(*anvio.A('hmmer-output-dir'), **anvio.K('hmmer-output-dir'))
    groupB.add_argument(*anvio.A('domain-hits-table'), **anvio.K('domain-hits-table'))
    groupC = parser.add_argument_group("tRNAs", "Through this program you can also scan Transfer RNA sequences in your "
                                                "contigs database for free (instead of running `anvi-scan-trnas` later).")
    groupC.add_argument(*anvio.A('also-scan-trnas'), **anvio.K('also-scan-trnas'))

    groupD = parser.add_argument_group("PERFORMANCE", "Stuff everyone forgets to set and then get upset with how slow "
                                                      "science goes.")
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('hmmer-program'), **anvio.K('hmmer-program'))

    groupE = parser.add_argument_group("AUTHORITY", "Because you are the boss.")
    groupE.add_argument(*anvio.A('add-to-functions-table'), **anvio.K('add-to-functions-table'))
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
