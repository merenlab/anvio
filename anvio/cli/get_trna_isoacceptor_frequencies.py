#!/usr/bin/env python
"""Get tRNA isoacceptor (anticodon-level) decoding demand from genomes, genes, and functions."""

import sys

import anvio
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.trnaisoacceptors as trnaisoacceptors
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10', 'meren']
__requires__ = ['contigs-db']
__can_use__ = ['profile-db', 'collection', 'bin', 'internal-genomes', 'external-genomes', 'hmm-hits']
__provides__ = ['trna-isoacceptor-frequencies-txt']
__description__ = (
    "Get tRNA isoacceptor (anticodon-level) decoding demand from genomes, genes, and functions")


@terminal.time_program
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
    """Prepare arguments to get tRNA isoacceptor frequencies."""

    args = get_args()
    run = terminal.Run()

    if args.get_default_decoding_weights:
        filesnpaths.is_output_file_writable(
            args.get_default_decoding_weights, ok_if_exists=False)
        trnaisoacceptors.default_decoding_weights_df.to_csv(
            args.get_default_decoding_weights, sep='\t')
        run.info("Default decoding weights file", args.get_default_decoding_weights)
        sys.exit(0)

    if (bool(args.output_file) +
        bool(args.gene_table_output) +
        bool(args.function_table_output)) != 1:
        raise ConfigError("One and only one of `--output-file`, `--gene-table-output`, or "
                          "`--function-table-output` must be provided.")
    args.return_functions = bool(args.function_table_output)

    # A generic output file is always assumed to be a per-gene table.
    if args.output_file:
        args.gene_table_output = args.output_file

    if args.gene_table_output:
        filesnpaths.is_output_file_writable(args.gene_table_output, ok_if_exists=False)
    if args.function_table_output:
        filesnpaths.is_output_file_writable(args.function_table_output, ok_if_exists=False)

    if args.function_sources is None:
        from_function_sources = False
    elif len(args.function_sources) != 1 and (args.function_accessions or args.function_names):
        raise ConfigError(
            "`--function-accessions` and `--function-names` require a single value for "
            "`--function-sources`. If select functions come from more than one source, use "
            "`--select-functions-txt`.")
    elif len(args.function_sources) == 0:
        from_function_sources = True
    else:
        from_function_sources = args.function_sources

    if (args.function_accessions or args.function_names) and args.select_functions_txt:
        raise ConfigError("`--function-accessions` and `--function-names` should not be used "
                          "with `--select-functions-txt`.")

    args.function_accessions_dict = {}
    args.function_names_dict = {}
    if isinstance(from_function_sources, list):
        function_source = from_function_sources[0]
        function_accessions = args.function_accessions if args.function_accessions else []
        for function_accession in function_accessions:
            try:
                args.function_accessions_dict[function_source].append(function_accession)
            except KeyError:
                args.function_accessions_dict[function_source] = [function_accession]
        function_names = args.function_names if args.function_names else []
        for function_name in function_names:
            try:
                args.function_names_dict[function_source].append(function_name)
            except KeyError:
                args.function_names_dict[function_source] = [function_name]

    if args.include_amino_acids and args.exclude_amino_acids:
        raise ConfigError(
            "Either `--include-amino-acids` or `--exclude-amino-acids` should be given, not both.")

    if args.include_amino_acids:
        args.exclude_amino_acids = []
        for amino_acid in constants.amino_acids:
            if amino_acid in args.include_amino_acids:
                continue
            else:
                args.exclude_amino_acids.append(amino_acid)

    args.pansequence_min_amino_acids = [int(args.pansequence_min_amino_acids[0]),
                                        float(args.pansequence_min_amino_acids[1])]

    if args.deviation_from_genome_average and (args.sum or args.average):
        raise ConfigError(
            "`--deviation-from-genome-average` cannot be combined with `--sum`/`--average`: once "
            "genes or functions are collapsed to a single row per genome, that row *is* the "
            "genome average, so every deviation value would trivially be 1.0. Please use one or "
            "the other.")

    if args.deviation_from_genome_average and not args.relative:
        raise ConfigError(
            "`--deviation-from-genome-average` requires `--relative`. Absolute weighted decoding "
            "demand scales with how many codons a gene or function has, so without `--relative`, "
            "a longer gene would appear to deviate from the genome average on every isoacceptor "
            "just by virtue of its length, not because its isoacceptor usage is actually biased. "
            "Please add `--relative`.")

    args.decoding_weights = trnaisoacceptors.parse_decoding_weights_table(args.decoding_weights_txt)

    calculator = trnaisoacceptors.TRNAIsoacceptorFrequencyCalculator(args, r=run)
    frequency_df = calculator.get_frequencies(
        from_function_sources=from_function_sources,
        return_functions=args.return_functions,
        gene_caller_ids=args.gene_caller_ids,
        function_accessions=args.function_accessions_dict,
        function_names=args.function_names_dict,
        expect_functions=args.expect_functions,
        gene_min_codons=args.gene_min_codons,
        function_min_codons=args.function_min_codons,
        min_codon_filter=args.min_codon_filter,
        drop_amino_acids=args.exclude_amino_acids,
        sequence_min_amino_acids=args.sequence_min_amino_acids,
        pansequence_min_amino_acids=args.pansequence_min_amino_acids,
        relative=args.relative,
        sum_genes=args.sum,
        average_genes=args.average,
        deviation_from_genome_average=args.deviation_from_genome_average,
        infinity_to_zero=args.infinity_to_zero)

    # Write output tables.
    if args.sum or args.average:
        table_output = args.gene_table_output if args.gene_table_output else args.function_table_output
        frequency_df.to_csv(table_output, sep='\t')
        if args.sum:
            run.info("Gene sum output", table_output)
        elif args.average:
            run.info("Gene average output", table_output)
    elif args.function_table_output:
        frequency_df.to_csv(args.function_table_output, sep='\t')
        run.info("Function table output", args.function_table_output)
    else:
        frequency_df.to_csv(args.gene_table_output, sep='\t')
        run.info("Gene table output", args.gene_table_output)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group(
        'SINGLE GENOME INPUTS',
        "Get tRNA isoacceptor frequencies from genes or functions in a single genome. A contigs "
        "database can be provided alone as an 'external' genome. An 'internal' genome (bin) also "
        "requires a profile database, collection name, and bin ID. The contigs database must have "
        "been processed by `anvi-scan-trnas` (or `anvi-run-hmms --also-scan-trnas`).")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupA.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupA.add_argument('--gene-caller-ids', type=int, nargs='+', help="Select genes by ID, space-separated.")

    groupB = parser.add_argument_group(
        'MULTIPLE GENOME INPUTS',
        "Get frequencies from genes or functions in multiple genomes by providing internal and/or "
        "external genome files listing the genomes to analyze.")
    groupB.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupB.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))

    groupC = parser.add_argument_group(
        'OUTPUT FILES',
        "This program writes tRNA isoacceptor frequency tables. When functions rather than genes "
        "are analyzed, one of two possible files can be produced: a table of frequencies per gene "
        "or per function, using `--gene-table-output` or `--function-table-output`, respectively.")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupC.add_argument(
        '--gene-table-output',
        help="A tab-delimited file of genes x isoacceptors. The index columns before frequency "
             "data contain genome names (optional, if multiple genomes are considered) and gene "
             "callers IDs.")
    groupC.add_argument(
        '--function-table-output',
        help="A tab-delimited file of functions x isoacceptors. Index columns before the "
             "frequency data contain, respectively, genome names (optional, if multiple genomes "
             "are considered), function annotation sources, accessions, and names.")
    groupC.add_argument(
        '--infinity-to-zero', default=False, action='store_true',
        help="Under `--deviation-from-genome-average`, replace any resulting +/-infinity values "
             "(a nonzero value divided by a zero genome average) with 0.0, and also replace a "
             "0/0 quotient (an isoacceptor whose value and genome average were both exactly "
             "zero) with 0.0. An isoacceptor genuinely absent from a genome is left as NaN "
             "either way. Use with caution, since infinity/NaN and 0.0 mean different things.")

    groupD = parser.add_argument_group('FREQUENCY STATISTICS',
        "Rather than absolute weighted decoding demand (default), relative composition or "
        "fold-change deviation from the genome average can be returned. Frequencies can also be "
        "summed or averaged across genes in a genome.")
    groupD.add_argument(*anvio.A('relative'), **anvio.K('relative'))
    groupD.add_argument(*anvio.A('sum'), **anvio.K('sum'))
    groupD.add_argument(*anvio.A('average'), **anvio.K('average'))
    groupD.add_argument(
        '--deviation-from-genome-average', default=False, action='store_true',
        help="Report each gene's or function's isoacceptor value as a fold-change against that "
             "genome's average isoacceptor profile. Requires `--relative`, since without it a "
             "longer gene would appear to deviate from the genome average on every isoacceptor "
             "purely due to its length. Cannot be combined with `--sum`/`--average`.")

    groupE = parser.add_argument_group('FUNCTIONS',
        "Frequencies can be calculated for functions rather than genes, summing the frequencies "
        "of the genes annotated by each function. Genes can also be subsetted to those annotated "
        "with requested functions or annotated by requested function sources.")
    groupE.add_argument(*anvio.A('function-sources'), **anvio.K('function-sources'))
    groupE.add_argument(*anvio.A('function-accessions'), **anvio.K('function-accessions'))
    groupE.add_argument(*anvio.A('function-names'), **anvio.K('function-names'))
    groupE.add_argument(*anvio.A('select-functions-txt'), **anvio.K('select-functions-txt'))
    groupE.add_argument(*anvio.A('expect-functions'), **anvio.K('expect-functions'))
    groupE.add_argument(*anvio.A('shared-function-sources'), **anvio.K('shared-function-sources'))

    groupF = parser.add_argument_group('FILTER GENES, FUNCTIONS, CODONS',
        "Genes/functions can be filtered by the number of codons they contain, and codons can be "
        "selected or excluded by amino acid, exactly as in `anvi-get-codon-frequencies`.")
    groupF.add_argument(*anvio.A('gene-min-codons'), **anvio.K('gene-min-codons'))
    groupF.add_argument(*anvio.A('function-min-codons'), **anvio.K('function-min-codons'))
    groupF.add_argument(*anvio.A('exclude-amino-acids'), **anvio.K('exclude-amino-acids'))
    groupF.add_argument(*anvio.A('include-amino-acids'), **anvio.K('include-amino-acids'))
    groupF.add_argument(*anvio.A('sequence-min-amino-acids'), **anvio.K('sequence-min-amino-acids'))
    groupF.add_argument(*anvio.A('pansequence-min-amino-acids'), **anvio.K('pansequence-min-amino-acids'))
    groupF.add_argument(*anvio.A('min-codon-filter'), **anvio.K('min-codon-filter'))

    groupG = parser.add_argument_group('tRNA MODIFICATION',
        "Wobble modification state (A34 to inosine, C34 to lysidine) can be inferred from the "
        "presence of the genes encoding the modifying enzymes, given a tRNA modification enzyme "
        "list. Without one, every isoacceptor is treated as canonically unmodified.")
    groupG.add_argument(
        '--trna-modification-enzyme-list', metavar='FILE_PATH',
        help="A tab-delimited file describing known tRNA modification enzymes, the genomic "
             "functions that encode them, and the tRNA positions they target (see the "
             "trna-modification-enzyme-list artifact for the exact format).")
    groupG.add_argument(
        '--decoding-weights-txt', metavar='FILE_PATH',
        help="A tab-delimited file of decoding efficiencies to use instead of the default "
             "bacterial wobble weights of Sabi and Tuller (2014). Use "
             "`--get-default-decoding-weights` to write a template.")
    groupG.add_argument(
        '--get-default-decoding-weights', metavar='FILE_PATH',
        help="Write the default decoding weights table to this path and exit, ignoring every "
             "other argument.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
