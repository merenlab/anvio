#!/usr/bin/env python
DESCRIPTION = """This program generates a metabolic reaction network in an anvi'o contigs or pan database"""

import sys

import anvio.filesnpaths as filesnpaths

from argparse import Namespace

from anvio.argparse import ArgumentParser
from anvio.reactionnetwork import Constructor
from anvio.errors import ConfigError, FilesNPathsError
from anvio import A as A, K as K, __version__ as VERSION


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = VERSION
__authors__ = ["semiller10"]
__requires__ = ["kegg-functions"]
__can_use__ = ["contigs-db", "enzymes-txt", "reaction-ref-data", "kegg-data"]
__provides__ = ["reaction-network"]
__description__ = DESCRIPTION


def main() -> None:
    args = get_args()

    try:
        constructor = Constructor(kegg_dir=args.kegg_dir, modelseed_dir=args.modelseed_dir)

        if args.enzymes_txt:
            if not args.output_json:
                raise ConfigError(
                    "When using `--enzymes-txt`, you must also provide `--output-json` to specify "
                    "where the resulting reaction network JSON file should be written."
                )
            filesnpaths.is_output_file_writable(args.output_json)
            network = constructor.make_enzymes_txt_network(
                enzymes_txt=args.enzymes_txt,
                stats_file=args.stats_file
            )
            network.export_json(
                args.output_json,
                overwrite=args.overwrite_existing_network
            )
        elif args.contigs_db:
            constructor.make_network(
                contigs_db=args.contigs_db,
                overwrite_existing_network=args.overwrite_existing_network,
                stats_file=args.stats_file
            )
        elif args.pan_db or args.genomes_storage:
            constructor.make_network(
                pan_db=args.pan_db,
                genomes_storage_db=args.genomes_storage,
                overwrite_existing_network=args.overwrite_existing_network,
                consensus_threshold=args.consensus_threshold,
                discard_ties=args.discard_ties,
                stats_file=args.stats_file
            )
        else:
            raise ConfigError(
                "A data source is required. Provide an enzymes file (`--enzymes-txt`), a contigs "
                "database (`--contigs-db`), or a pan database (`--pan-db`) with a genomes storage "
                "database (`--genomes-storage`)."
            )
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args() -> Namespace:
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group(
        "ENZYMES FILE INPUT",
        "Generate a reaction network from a tab-delimited enzymes file, bypassing the need for a "
        "contigs database. This enables reaction networks from custom enzyme lists such as "
        "transcriptomic/proteomic data or LCA gene content predictions. The file must have columns "
        "'gene_id', 'enzyme_accession', and 'source'; only rows with source 'KOfam' are used. "
        "See 'anvi-estimate-metabolism' for the full enzymes-txt format specification."
    )
    groupA.add_argument(
        '--enzymes-txt', type=str, metavar='FILE', help=
        "Path to a tab-delimited enzymes file with required columns 'gene_id', "
        "'enzyme_accession', and 'source'. Rows where 'source' is 'KOfam' are used to build "
        "the reaction network. Use '--output-json' to specify where the resulting network file "
        "is written."
    )
    groupA.add_argument(
        '--output-json', type=str, metavar='FILE', help=
        "Path for the output reaction network JSON file. Required when using '--enzymes-txt'. "
        "The file can be loaded by 'anvi-draw-kegg-pathways' via its '--reaction-network-json' "
        "option."
    )

    groupB = parser.add_argument_group(
        "SINGLE GENOME OR METAGENOME INPUT",
        "Generate a reaction network from a contigs database and store the network in the database."
    )
    groupB.add_argument(*A('contigs-db'), **K('contigs-db', {'required': False}))

    groupC = parser.add_argument_group(
        "PANGENOME INPUT",
        "Generate a reaction network from a pan database and associated genomes storage database, "
        "and store the network in the pan database."
    )
    groupC.add_argument(*A('pan-db'), **K('pan-db', {'required': False}))
    groupC.add_argument(*A('genomes-storage'), **K('genomes-storage', {'required': False}))
    groupC.add_argument(
        '--consensus-threshold', default=None, type=float, metavar='FLOAT', help=
        "If this argument is provided, then a protein annotation must be assigned to this minimum "
        "proportion of genes in a cluster to be imputed to the cluster as a whole. By default, "
        "without this argument, the annotation assigned to the most genes becomes the annotation "
        "of the cluster (also see --discard-ties). The consensus threshold must be a number from 0 "
        "to 1."
    )
    groupC.add_argument(
        '--discard-ties', default=False, action='store_true', help=
        "By default, a gene cluster is assigned a protein annotation by finding the protein "
        "ortholog that occurs in the greatest number of genes in the cluster (see "
        "--consensus-threshold) and arbitrarily choosing one ortholog in case of a tie. With this "
        "flag, a tie instead results in an ortholog annotation not being assigned to the cluster."
    )

    groupD = parser.add_argument_group(
        "DATABASE", "KEGG and ModelSEED reference database information"
    )
    groupD.add_argument(
        '--kegg-dir', type=str, metavar='PATH', help=
        "Path to KEGG database directory. If this option is not used, the program expects a "
        "database set up in the default location used by 'anvi-setup-kegg-data'."
    )
    groupD.add_argument(
        '--modelseed-dir', type=str, metavar='PATH', help=
        "Path to ModelSEED Biochemistry database directory. If this option is not used, the "
        "program expects a database set up in the default location used by "
        "'anvi-setup-modelseed-database'."
    )

    groupE = parser.add_argument_group("OTHER OPTIONS")
    groupE.add_argument(
        '--overwrite-existing-network', default=False, action='store_true', help=
        "Overwrite an existing reaction network in the database (for contigs/pan DB modes) or "
        "overwrite the output JSON file (for enzymes-txt mode) if it already exists."
    )
    groupE.add_argument(
        '--stats-file', type=str, metavar='PATH', help=
        "Write a tab-delimited file of network overview statistics (statistics also printed to the "
        "terminal) to the output path."
    )

    args = parser.get_args(parser)
    return args


if __name__ == '__main__':
    main()
