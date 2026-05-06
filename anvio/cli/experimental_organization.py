#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering

from anvio.errors import ConfigError, FilesNPathsError
from anvio.clusteringconfuguration import ClusteringConfiguration


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['clustering-configuration']
__provides__ = ['dendrogram']
__description__ = "Create an experimental clustering dendrogram."
__resources__ = [("An example use of this program", "https://merenlab.org/2015/06/10/combining-omics-data/")]


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    split_names = []

    try:
        if not args.skip_store_in_db and not args.profile_db and not args.dry_run:
            raise ConfigError("If you are not using --skip-store-in-db or --dry-run in your call, you need to define a "
                               "profile database to be explicit about where to store the clustering result.")

        if args.profile_db:
            utils.is_profile_db(args.profile_db)
        utils.is_contigs_db(args.contigs_db)

        if not args.skip_store_in_db and not args.name and not args.dry_run:
            raise ConfigError("Since you have not used the --skip-store-in-db or --dry-run flags, this program "
                               "will perform clustering, and attempt to store the results in the profile database. "
                               "When it attempts to do it, it will need a name for this new clustering results. "
                               "Please use the --name flag to provide one.")

        if args.name:
            utils.is_this_name_OK_for_database('--name parameter', args.name)

        if args.profile_db:
            split_names = utils.get_all_item_names_from_the_database(args.profile_db)

        db_paths = {'CONTIGS.db': args.contigs_db}
        config = ClusteringConfiguration(args.config_file, args.input_directory, db_paths = db_paths, row_ids_of_interest = split_names)

        config.print_summary(run)

        if args.dry_run and not args.export_merged_matrix:
            sys.exit()

        if args.dry_run and args.export_merged_matrix:
            # the user wants the merged matrix without the costly clustering step
            clustering.scale_and_combine_matrices(config, progress=progress, run=run, debug=True)
        else:
            clustering_id, newick = clustering.order_contigs_simple(config, distance=args.distance, linkage=args.linkage, progress=progress, debug=True)

            _, distance, linkage = clustering_id.split(':')

            run.info('Distance metric used', distance, mc='green')
            run.info('Linkage method used', linkage, mc='green')

            if args.output_file:
                open(args.output_file, 'w').write(newick + '\n')
                run.info('Output', args.output_file, mc='green')

        # this here exports the final merged data matrix (i.e., the combined and scaled data from all input
        # matrices that will be used to generate dendrogram) as a TAB-delimited file, which can be useful for
        # debugging and/or for users who want to use this merged matrix for other purposes (e.g., to generate
        # a dendrogram using another software). this works with --dry-run so the user can get the matrix
        # without paying for the clustering step.
        if args.export_merged_matrix:
            header = ['items']
            for matrix in config.matrices:
                m = config.matrices_dict[matrix]
                n_cols = m['scaled_vectors'].shape[1]
                header.extend(['%s_%04d' % (matrix, j) for j in range(n_cols)])

            rows = []
            for i in range(len(config.combined_vectors)):
                rows.append([config.combined_id_to_sample[i]] + config.combined_vectors[i].tolist())

            utils.store_array_as_TAB_delimited_file(rows, args.export_merged_matrix, header)
            run.info('Merged matrix', args.export_merged_matrix, mc='green')

        if args.profile_db and not args.skip_store_in_db and not args.dry_run:
            dbops.add_items_order_to_db(anvio_db_path=args.profile_db,
                                        order_name=args.name,
                                        order_data=newick,
                                        distance=distance,
                                        linkage=linkage,
                                        make_default=False,
                                        run=run)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument('config_file', metavar = 'FILE', default = None, type=str,
                        help = 'Config file for clustering of contigs. See documentation for help.')

    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('experimental-org-input-dir'), **anvio.K('experimental-org-input-dir'))
    parser.add_argument(*anvio.A('clustering-name'), **anvio.K('clustering-name'))
    parser.add_argument(*anvio.A('distance'), **anvio.K('distance', {'default': None, 'help':
                      'The distance metric for the hierarchical clustering. If you do not use this flag,\
                       the distance metric you defined in your clustering config file will be used. If you\
                       have not defined one in your config file, then the system default will be used,\
                       which is "%s".' % constants.distance_metric_default}))
    parser.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'default': None, 'help':
                      'Same story with the `--distance`, except, the system default for this one\
                       is %s.' % constants.linkage_method_default}))
    parser.add_argument(*anvio.A('skip-store-in-db'), **anvio.K('skip-store-in-db'))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    parser.add_argument('--export-merged-matrix', metavar='FILE', default=None, help="Export the final merged "
                             "data matrix (the combined and scaled data from all input matrices that is used "
                             "to generate the dendrogram) as a TAB-delimited file.")
    parser.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

