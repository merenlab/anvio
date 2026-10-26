#!/usr/bin/env python
"""Entry point to the interactive interface.

The massage of the data is being taken care of in the interactive module,
and this file implements the bottle callbacks."""

import sys
from anvio.argparse import ArgumentParser, PopulateAnvioDBArgs

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.interactive as interactive
from anvio.bottleroutes import BottleApplication

from anvio.errors import ConfigError, FilesNPathsError, GenesDBError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ["Doğan Can Kilment", "Gökmen Göksel", "Gökmen Görgen"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ozcan', 'matthewlawrenceklein', 'isaacfink21', 'blankenberg']
__resources__ = [("A beginners tutorial on anvi'o interactive interface", "http://merenlab.org/tutorials/interactive-interface/"),
                 ("How to add more data to a display for layers and items", "http://merenlab.org/2017/12/11/additional-data-tables/"),
                 ("An overview of interactive data types", "http://merenlab.org/2016/02/27/the-anvio-interactive-interface/"),
                 ("Anvi'o 'views' demystified", "http://merenlab.org/2017/05/08/anvio-views/"),
                 ("Working with SVG files from the interactive interface", "http://merenlab.org/2016/10/27/high-resolution-figures/"),
                 ("Running remote anvi'o interactive interfaces on your local computer", "http://merenlab.org/2018/03/07/working-with-remote-interative/")]
__requires__ = ['profile-db']
__can_use__ = ['single-profile-db', 'contigs-db', 'genes-db', 'bin', 'view-data', 'dendrogram', 'phylogeny']
__provides__ = ['interactive']
__can_provide__ = ['svg', 'collection', 'bin', 'contig-inspection', 'dendrogram', 'phylogeny']
__description__ = "Start an anvi'o server for the interactive interface"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        args = fill_in_anvio_dbs(args)
        utils.is_all_npm_packages_installed()
        d = interactive.Interactive(args)
        args.port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except GenesDBError as e:
        print(e)
        sys.exit(-3)

    if args.dry_run:
        run.info_single('Dry run, eh? Fine. Bai!', nl_after=1)
        sys.exit()

    try:
        app = BottleApplication(d)
        app.run_application(args.ip_address, args.port_number)
    except ConfigError as e:
        print(e)
        sys.exit(-4)


def fill_in_anvio_dbs(args):
    """Fill in missing anvi'o databases from the work directory, unless the user is visualizing ad hoc data.

    When the user provides a tree or a view data file without a profile or contigs database, they
    may want to visualize those files alone, or may want to see them on top of anvi'o databases that
    happen to be in the work directory. We have no way of knowing which one, so if the search finds
    anvi'o databases in such cases, we ask the user to be precise rather than guessing.
    """

    # if the user explicitly asked for the manual mode, guessing databases is never the right thing to do
    if args.manual_mode:
        return args

    ad_hoc_inputs = [f for f, v in [('--tree', args.tree), ('--view-data', args.view_data)] if v]
    may_be_ad_hoc = bool(ad_hoc_inputs) and not args.profile_db and not args.contigs_db

    args = PopulateAnvioDBArgs(args, lazy_init=not anvio.DEBUG_AUTO_FILL_ANVIO_DBS).get_updated_args()

    if may_be_ad_hoc and args.profile_db:
        dbs_found = ' and '.join([f"'{db}'" for db in [args.profile_db, args.contigs_db] if db])
        raise ConfigError(f"You provided {' and '.join([f'`{f}`' for f in ad_hoc_inputs])} without a profile or contigs "
                          f"database, BUT anvi'o found {dbs_found} in your work directory, and now it is not sure "
                          f"what you want. If you would like to visualize your files ALONE, please use the flag "
                          f"`--manual-mode` along with a profile database path to store your states and collections "
                          f"(i.e., `-p MY-DISPLAY.db`, which does not need to exist). If you would like to visualize "
                          f"your files ALONG WITH the anvi'o databases anvi'o found, please provide them explicitly "
                          f"with `-p` and `-c` parameters.")

    return args


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DEFAULT INPUTS', "The interactive interface can be started with and without\
                                                          anvi'o databases. The default use assumes you have your\
                                                          profile and contigs database, however, it is also possible\
                                                          to start the interface using ad hoc input files. See 'MANUAL\
                                                          INPUT' section for required parameters.")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name', {'help':
                                    "If you have a collection in your profile database, you can use this flag to start the\
                                    interactive interface with a tree showing your bins in your collection, instead of each\
                                    split. This is very useful when you have imported your external binning results into\
                                    anvi'o, and want to see the distribution of your bins across samples. In these cases\
                                    anvi'o will cluster your bins and based on multiple metrics. Because this particular\
                                    clustering will be done on the fly within anvi'o interactive class, you get to define\
                                    a distance metric and a linkage method using --linkage and --distance parameters if\
                                    you want!"}))

    groupB = parser.add_argument_group('MANUAL INPUTS', "Input parameters to start the interactive interface without anvi'o\
                                                         databases. When you provide these without a contigs database,\
                                                         anvi'o will start the interface in manual mode automatically.\
                                                         You still need to provide a profile database path, which does\
                                                         not need to exist: anvi'o will create an empty one for you to\
                                                         store your states and collections.")
    groupB.add_argument(*anvio.A('manual-mode'), **anvio.K('manual-mode'))
    groupB.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file'))
    groupB.add_argument(*anvio.A('view-data'), **anvio.K('view-data'))
    groupB.add_argument(*anvio.A('tree'), **anvio.K('tree'))
    groupB.add_argument(*anvio.A('items-order'), **anvio.K('items-order'))

    groupC = parser.add_argument_group('ADDITIONAL STUFF', "Parameters to provide additional layers, views, or layer data.")
    groupC.add_argument(*anvio.A('additional-view'), **anvio.K('additional-view'))
    groupC.add_argument(*anvio.A('additional-layers'), **anvio.K('additional-layers'))
    groupC.add_argument(*anvio.A('annotation-source-for-per-split-summary'), **anvio.K('annotation-source-for-per-split-summary'))

    groupD = parser.add_argument_group('GENE MODE', "Gene mode related parameters.")
    groupD.add_argument(*anvio.A('gene-mode'), **anvio.K('gene-mode'))
    groupD.add_argument(*anvio.A('inseq-stats'), **anvio.K('inseq-stats'))
    groupD.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupE = parser.add_argument_group('VISUALS RELATED', "Parameters that give access to various adjustements regarding\
                                                           the interface.")
    groupE.add_argument(*anvio.A('view'), **anvio.K('view'))
    groupE.add_argument(*anvio.A('title'), **anvio.K('title'))
    groupE.add_argument(*anvio.A('taxonomic-level'), **anvio.K('taxonomic-level'))
    groupE.add_argument(*anvio.A('show-all-layers'), **anvio.K('show-all-layers'))
    groupE.add_argument(*anvio.A('split-hmm-layers'), **anvio.K('split-hmm-layers'))
    groupE.add_argument(*anvio.A('hide-outlier-SNVs'), **anvio.K('hide-outlier-SNVs'))
    groupE.add_argument(*anvio.A('state-autoload'), **anvio.K('state-autoload'))
    groupE.add_argument(*anvio.A('collection-autoload'), **anvio.K('collection-autoload'))
    groupE.add_argument(*anvio.A('export-svg'), **anvio.K('export-svg'))


    groupF = parser.add_argument_group('SWEET PARAMS OF CONVENIENCE', "Parameters and flags that are not quite essential (but\
                                                                       nice to have).")
    groupF.add_argument(*anvio.A('show-views'), **anvio.K('show-views'))
    groupF.add_argument(*anvio.A('skip-check-names'), **anvio.K('skip-check-names'))
    groupF.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupF.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))
    groupF.add_argument(*anvio.A('show-states'), **anvio.K('show-states'))
    groupF.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupF.add_argument(*anvio.A('skip-init-functions'), **anvio.K('skip-init-functions'))
    groupF.add_argument(*anvio.A('skip-auto-ordering'), **anvio.K('skip-auto-ordering'))
    groupF.add_argument(*anvio.A('skip-news'), **anvio.K('skip-news'))
    groupF.add_argument(*anvio.A('distance'), **anvio.K('distance', {'help':
                                    'The distance metric for the hierarchical clustering. Only relevant if you are running\
                                     the interactive interface in "collection" mode. The default is "%(default)s".'}))
    groupF.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'help':
                                    'The linkage method for the hierarchical clustering. Only relevant if you are running\
                                     the interactive interface in "collection" mode. The default is "%(default)s".'}))


    groupG = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupG.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupG.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupG.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupG.add_argument(*anvio.A('read-only'), **anvio.K('read-only'))
    groupG.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupG.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))
    groupG.add_argument(*anvio.A('user-server-shutdown'), **anvio.K('user-server-shutdown'))

    # we will fill in missing anvi'o databases in `main` ourselves, since we need to know what the
    # user provided explicitly before we do that (see `fill_in_anvio_dbs` for details)
    return parser.get_args(parser)


if __name__ == '__main__':
    main()
