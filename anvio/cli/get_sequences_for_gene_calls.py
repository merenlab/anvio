#!/usr/bin/env python
# -*- coding: utf-8
"""Returns sequences for a given list of gene caller ids"""

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.genomestorage as genomestorage

from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbops import ContigsSuperclass, ContigsDatabase


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'genomes-storage-db']
__provides__ = ['genes-fasta', 'external-gene-calls']
__description__ = "A script to get back sequences for gene calls"
__resources__ = [("A tutorial on getting gene-level taxonomy for a contigs-db", "http://merenlab.org/2016/06/18/importing-taxonomy/")]


@terminal.time_program
def main():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    contigs_db_path = A('contigs_db')
    genomes_storage_db_path = A('genomes_storage')
    output_file_path = A('output_file')

    try:
        if not contigs_db_path and not genomes_storage_db_path:
            raise ConfigError("You must give this program either a contigs or a genomes storage database so it can like "
                              "export sequences for your genes? :/")

        if contigs_db_path and genomes_storage_db_path:
            raise ConfigError("You can either ask for sequences in a contigs database, or a genomes storage. But not both "
                              "(obviously).")

        if contigs_db_path:
            export_from_contigs(args)
        elif genomes_storage_db_path:
            export_from_genomes_storage(genomes_storage_db_path, output_file_path)
        else:
            raise ConfigError("o_O")
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def export_from_contigs(args):
    if args.list_annotation_sources:
        ContigsDatabase(args.contigs_db).list_function_sources()
        sys.exit()

    c = ContigsSuperclass(args)

    output_file_path = args.output_file if args.output_file else 'sequences_for_gene_calls.txt'
    filesnpaths.is_output_file_writable(output_file_path)

    gene_caller_ids = list(utils.get_gene_caller_ids_from_args(args.gene_caller_ids, args.delimiter))

    func_kwargs = dict(
        gene_caller_ids_list=gene_caller_ids,
        output_file_path=args.output_file,
        simple_headers=not args.report_extended_deflines,
        wrap=args.wrap
    )

    if args.export_gff3 and (args.list_defline_variables or args.defline_format):
        raise ConfigError("User-defined deflines for gene sequence output files is not compatible with the "
                          "GFF3 output mode. Sorry :/")

    if args.export_gff3:
        if args.get_aa_sequences:
            raise ConfigError("AA sequences can only be reported in FASTA format, please remove the --export-gff3 flag to continue.")
        if args.external_gene_calls:
            raise ConfigError("External gene calls file is not compatble with GFF format.")
        if args.flank_length:
            raise ConfigError("Due to the lazy programming, --flank-length and --export-gff3 flags are incompatible :(")
        if args.annotation_source:
            func_kwargs['gene_annotation_source'] = args.annotation_source
        if args.return_all_function_hits_for_each_gene:
            anvio.RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE = True

        c.gen_GFF3_file_of_sequences_for_gene_caller_ids(**func_kwargs)
    else:
        if args.annotation_source:
            raise ConfigError("A functional annotation source is only relevant for the GFF output format.")
        c.get_sequences_for_gene_callers_ids(**func_kwargs, report_aa_sequences=args.get_aa_sequences, flank_length=args.flank_length,
                                             output_file_path_external_gene_calls=args.external_gene_calls, defline_format=args.defline_format,
                                             list_defline_variables=args.list_defline_variables)


def export_from_genomes_storage(genomes_storage_db_path, output_file_path):
    if args.export_gff3:
        raise ConfigError("GFF output is only relevant if you are working with a contigs database :/")

    g = genomestorage.GenomeStorage(genomes_storage_db_path)
    g.gen_combined_aa_sequences_FASTA(output_file_path, report_with_genome_name=True)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('OPTION #1: GET SEQUENCES FROM A CONTIGS DB')
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    groupA.add_argument(*anvio.A('delimiter'), **anvio.K('delimiter'))
    groupA.add_argument(*anvio.A('flank-length'), **anvio.K('flank-length'))
    groupA.add_argument(*anvio.A('get-aa-sequences'), **anvio.K('get-aa-sequences'))
    groupA.add_argument(*anvio.A('external-gene-calls'), **anvio.K('external-gene-calls',
                                        {"help": ("An optional external gene calls file path that precisely "
                                                  "describes the set of gene sequences exported. Using this file "
                                                  "you can create an anvi'o contigs database from the resulting "
                                                  "genes FASTA file without having to do a gene calling from scratch.")}))

    groupB = parser.add_argument_group('OPTION #1.1: STORE SEQEUNCES AS GFF?', "This only works with OPTION #1 :/")
    groupB.add_argument(*anvio.A('export-gff3'), **anvio.K('export-gff3'))
    groupB.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source'))
    groupB.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))
    groupB.add_argument(*anvio.A('return-all-function-hits-for-each-gene'), **anvio.K('return-all-function-hits-for-each-gene'))

    groupC = parser.add_argument_group('OPTION #2: GET SEQUENCES FROM A GENOMES STORAGE')
    groupC.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))
    groupC.add_argument(*anvio.A('genomes-names'), **anvio.K('genomes-names'))

    groupD = parser.add_argument_group('FASTA OUTPUT FORMATTING OPTIONS')
    groupD.add_argument(*anvio.A('report-extended-deflines'), **anvio.K('report-extended-deflines'))
    groupD.add_argument(*anvio.A('wrap'), **anvio.K('wrap'))
    groupD.add_argument(*anvio.A('list-defline-variables'), **anvio.K('list-defline-variables'))
    groupD.add_argument(*anvio.A('defline-format'), **anvio.K('defline-format'))

    groupE = parser.add_argument_group('OPTIONS COMMON TO ALL INPUTS')
    groupE.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
