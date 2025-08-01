#!/usr/bin/env python
# -*- coding: utf-8
"""Further analyze one or more bins in a collection.

   This is especially useful when there are one or more highly contaminated
   bins in a merged profile.
"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.interactive as interactive
from anvio.bottleroutes import BottleApplication

from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'blankenberg']
__resources__ = [("Refining a bin", "http://merenlab.org/2015/05/11/anvi-refine/"),
                 ("Notes on genome refinement", "http://merenlab.org/2017/05/11/anvi-refine-by-veronika/"),
                 ("A case study: Inspecting the genomic link between Archaea and Eukaryota", "http://merenlab.org/2017/01/03/loki-the-link-archaea-eukaryota/"),
                 ("As part of the metagenomic workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-refine"),
                 ("A demo", "https://www.youtube.com/watch?v=vXPKP5vKiBM")
                 ]
__requires__ = ['profile-db', 'contigs-db', 'bin',]
__provides__ = ['bin',]
__description__ = ("Start an anvi'o interactive interactive to manually curate or refine a genome, "
                   "whether it is a metagenome-assembled, single-cell, or an isolate genome")


def main():
    args = get_args()
    run = terminal.Run()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    find_from_split_name = A('find_from_split_name')
    collection_name = A('collection_name')

    try:
        if not collection_name:
            raise ConfigError("You need to provide this program with a collection name :/")

        if find_from_split_name and not collection_name:
            raise ConfigError("If you would like anvi'o to find the bin name for your split, you should "
                              "also specify a collection name.")

        if find_from_split_name:
            rows = utils.get_bin_name_from_item_name(args.profile_db, find_from_split_name, collection_name=args.collection_name)

            if not rows:
                raise ConfigError("The split name you requested was not found in collection %s :/" % collection_name)

            if not len(rows) == 1:
                raise ConfigError("There is something silly going on here. The split name is found in more "
                                  "than one collection. Which is not really possible so goodbye.")

            entry_id, collection_name, split_name, bin_name = rows[0]
            run.warning("Anvi'o found your split, and will set the bin name for your "
                        "refinement analysis to '%s'." % (bin_name),
                        header="SPLIT FOUND IN %s!" % (bin_name), lc="green")

            args.bin_id = bin_name

        if not args.bin_id or args.bin_ids_file:
            raise ConfigError("This program needs to know which bin(s) you wish to refine.")

        args.mode = 'refine'
        d = interactive.Interactive(args)
        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except DictIOError as e:
        print(e)
        sys.exit(-3)

    if args.dry_run:
        run.info_single('Dry run, eh? Fine. Bai!', nl_after=1)
        sys.exit()

    app = BottleApplication(d)
    app.run_application(args.ip_address, args.port_number)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DEFAULT INPUTS', "The interavtive interface can be started with and without\
                                                          anvi'o databases. The default use assumes you have your\
                                                          profile and contigs database, however, it is also possible\
                                                          to start the interface using ad-hoc input files. See 'MANUAL\
                                                          INPUT' section for other set of parameters that are mutually\
                                                          exclusive with datanases.")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group('REFINE-SPECIFICS', "Parameters that are essential to the refinement process.")
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupB.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))
    groupB.add_argument(*anvio.A('find-from-split-name'), **anvio.K('find-from-split-name'))

    groupC = parser.add_argument_group('ADDITIONAL STUFF', "Parameters to provide additional layers, views, or layer data.")
    groupC.add_argument(*anvio.A('tree'), **anvio.K('tree'))
    groupC.add_argument(*anvio.A('skip-hierarchical-clustering'), **anvio.K('skip-hierarchical-clustering', {'help':
                                    'Skip hierarchical clustering for the splits in the refined bin, if you skip clustering \
                                    you need to provide your own newick formatted tree using --tree parameter.'}))
    groupC.add_argument(*anvio.A('load-full-state'), **anvio.K('load-full-state'))
    groupC.add_argument(*anvio.A('additional-view'), **anvio.K('additional-view'))
    groupC.add_argument(*anvio.A('additional-layers'), **anvio.K('additional-layers'))
    groupC.add_argument(*anvio.A('annotation-source-for-per-split-summary'), **anvio.K('annotation-source-for-per-split-summary'))

    groupD = parser.add_argument_group('VISUALS RELATED', "Parameters that give access to various adjustements regarding\
                                                           the interface.")
    groupD.add_argument(*anvio.A('show-all-layers'), **anvio.K('show-all-layers'))
    groupD.add_argument(*anvio.A('split-hmm-layers'), **anvio.K('split-hmm-layers'))
    groupD.add_argument(*anvio.A('taxonomic-level'), **anvio.K('taxonomic-level'))
    groupD.add_argument(*anvio.A('hide-outlier-SNVs'), **anvio.K('hide-outlier-SNVs'))
    groupD.add_argument(*anvio.A('title'), **anvio.K('title'))
    groupD.add_argument(*anvio.A('export-svg'), **anvio.K('export-svg'))

    groupE = parser.add_argument_group('SWEET PARAMS OF CONVENIENCE', "Parameters and flags that are not quite essential (but\
                                                                       nice to have).")
    groupE.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))
    groupE.add_argument(*anvio.A('skip-init-functions'), **anvio.K('skip-init-functions'))
    groupE.add_argument(*anvio.A('skip-news'), **anvio.K('skip-news'))

    groupF = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupE.add_argument(*anvio.A('skip-auto-ordering'), **anvio.K('skip-auto-ordering'))
    groupF.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupF.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupF.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupF.add_argument(*anvio.A('read-only'), **anvio.K('read-only'))
    groupF.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupF.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
