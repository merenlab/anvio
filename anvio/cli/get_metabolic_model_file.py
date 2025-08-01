#!/usr/bin/env python
# -*- coding: utf-8
DESCRIPTION = """This program exports a metabolic reaction network to a file suitable for flux balance analysis."""

from sys import exit
from argparse import Namespace

import anvio.reactionnetwork as reactionnetwork

from anvio import A, K
from anvio.errors import ConfigError
from anvio import __version__ as VERSION
from anvio.argparse import ArgumentParser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__authors__ = ['semiller10']
__requires__ = ['contigs-db', 'reaction-network']
__provides__ = ['reaction-network-json']
__description__ = DESCRIPTION


def main() -> None:
    args = get_args()

    record_genomes = []
    for arg in args.record_genomes:
        if arg == 'cluster':
            record_genomes.append('gene cluster')
        else:
            record_genomes.append(arg)
    record_genomes = tuple(record_genomes)

    try:
        constructor = reactionnetwork.Constructor()
        if args.contigs_db:
            network: reactionnetwork.GenomicNetwork = constructor.load_network(
                contigs_db=args.contigs_db,
                check_gene_annotations=not args.ignore_changed_gene_annotations
            )
            network.export_json(
                args.output_file,
                overwrite=args.overwrite_output_destinations,
                objective='e_coli_core',
                remove_missing_objective_metabolites=args.remove_missing_objective_metabolites
            )
        elif args.pan_db or args.genomes_storage_db:
            network: reactionnetwork.PangenomicNetwork = constructor.load_network(
                pan_db=args.pan_db,
                genomes_storage_db=args.genomes_storage,
                check_gene_annotations=not args.ignore_changed_gene_annotations
            )
            network.export_json(
                args.output_file,
                overwrite=args.overwrite_output_destinations,
                objective='e_coli_core',
                remove_missing_objective_metabolites=args.remove_missing_objective_metabolites,
                record_genomes=record_genomes
            )
        else:
            raise ConfigError(
                "Either a contigs database (`--contigs-db`) OR a pan database (`--pan-db`) and genomes "
                "storage database (`--genomes-storage`) must be provided. These contain the reaction "
                "network to be loaded and exported."
            )
    except ConfigError as e:
        print(e)
        exit(-1)


def get_args() -> Namespace:
    parser = ArgumentParser(description=DESCRIPTION)

    groupA = parser.add_argument_group(
        "SINGLE GENOME OR METAGENOME",
        "Export a reaction network stored in a contigs database."
    )
    groupA.add_argument(*A('contigs-db'), **K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group(
        "PANGENOME",
        "Export a reaction network stored in a pan database, with the associated genomes storage "
        "database also required."
    )
    groupB.add_argument(*A('pan-db'), **K('pan-db', {'required': False}))
    groupB.add_argument(*A('genomes-storage'), **K('genomes-storage', {'required': False}))
    groupB.add_argument(
        '--record-genomes', nargs='*', default=['cluster', 'reaction'],
        help=(
            "This option records the genome membership of gene clusters in JSON entries. By "
            "default, genome names are recorded for gene clusters and reactions: this is "
            "equivalent to the argument 'cluster reaction'. To not record genome membership at "
            "all, use this as a flag without any arguments, i.e., '--record-genomes'. The "
            "following arguments can be provided in any combination: 'cluster', 'reaction', and "
            "'metabolite'. 'reaction' and 'metabolite' record the genomes predicted to encode "
            "enzymes associated with reactions and metabolites, respectively."
        )
    )

    groupC = parser.add_argument_group("OUTPUT")
    groupC.add_argument(*A('output-file'), **K('output-file'))
    groupC.add_argument(*A('overwrite-output-destinations'), **K('overwrite-output-destinations'))

    groupD = parser.add_argument_group("OTHER OPTIONS")
    groupD.add_argument(
        '--remove-missing-objective-metabolites', default=False, action='store_true',
        help=(
            "The metabolic model file includes the core E. coli biomass objective function (BOF), "
            "as in the COBRApy example file, 'e_coli_core.json'. With this flag, metabolites in "
            "the BOF are removed if they are not produced or consumed by enzymes in the model. "
            "Without this flag, COBRApy cannot load the output JSON file if there are missing BOF "
            "metabolites. Removed BOF metabolite data are reported to the terminal."
        )
    )
    groupD.add_argument(
        '--ignore-changed-gene-annotations', default=False, action='store_true',
        help=(
            "It is possible that the gene KO annotations used to construct the stored reaction "
            "network have since been changed in the database. By default, without using this flag, "
            "this program checks that the set of gene KO annotations that is currently stored was "
            "also that used in construction of the reaction network, and raises an error if this "
            "is not the case. Use of this flag ignores that check, permitting the set of gene "
            "annotations to have changed since creation of the reaction network."
        )
    )

    args = parser.get_args(parser)
    return args


if __name__ == '__main__':
    main()
