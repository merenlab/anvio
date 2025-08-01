#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db']
__provides__ = ['gene-clusters-txt']
__description__ = ("Export gene clusters in a pan-db as a three-column, TAB-delimited file that associates each gene call "
                   "in each genome with a gene cluster")


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_pan_db(args.pan_db)

        if not args.output_file:
            raise ConfigError("You should provide an output file name so anvi'o does not have to make up a silly "
                               "name :/")

        filesnpaths.is_output_file_writable(args.output_file)

        pan_db = dbops.PanDatabase(args.pan_db)

        gene_clusters_dict = pan_db.db.get_table_as_dict(t.pan_gene_clusters_table_name)

        unique_gene_cluster_names = set()
        unique_genome_names = set()
        with open(args.output_file, 'w') as output:
            output.write("genome_name\tgene_caller_id\tgene_cluster_name\n")
            for d in gene_clusters_dict.values():
                unique_gene_cluster_names.add(d['gene_cluster_id'])
                unique_genome_names.add(d['genome_name'])
                output.write(f"{d['genome_name']}\t{d['gene_caller_id']}\t{d['gene_cluster_id']}\n")

        run.warning(f"Anvi'o has recovered {len(unique_gene_cluster_names)} gene clusters across "
                    f"{len(unique_genome_names)} genomes.", header="THINGS ARE WORKING", lc="green")

        run.info('Output', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT PAN DATABASE', "You want to export gene clusters from where?")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))

    groupB = parser.add_argument_group('MATTERS OF REPORTING')
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)

if __name__ == '__main__':
    main()
