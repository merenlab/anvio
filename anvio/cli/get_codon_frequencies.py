#!/usr/bin/env python
# -*- coding: utf-8
"""Get codon or amino acid frequency statistics from genomes, genes, and functions."""

import sys
import pandas as pd

import anvio
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.codonusage as codonusage
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10', 'meren']
__requires__ = ['contigs-db',
                'profile-db',
                'collection',
                'bin',
                'internal-genomes',
                'external-genomes']
__provides__ = ['codon-frequencies-txt', 'aa-frequencies-txt']
__description__ = (
    "Get codon or amino acid frequency statistics from genomes, genes, and functions")


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
    """Prepare arguments to get codon frequencies."""

    args = get_args()
    run = terminal.Run()

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
        raise ConfigError("`--function-accessions` and `--function-names` should not be used with "
                          "`--select-functions-txt`.")

    # Select function sources and functions of interest.
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

    from_function_sources = parse_select_functions_table(args, from_function_sources)

    # Store the coding dictionary in `args` for use in instantiating single- and multi-genome codon
    # usage objects.
    args.codon_to_amino_acid = codonusage.get_custom_encodings(args.encodings_txt)

    if args.include_amino_acids and args.exclude_amino_acids:
        raise ConfigError(
            "Either `--include-amino-acids` or `--exclude-amino-acids` should be given, not both.")

    # Amino acids to exclude are the complement of amino acids to include.
    if args.include_amino_acids:
        args.exclude_amino_acids = []
        for amino_acid in constants.amino_acids:
            if amino_acid in args.include_amino_acids:
                continue
            else:
                args.exclude_amino_acids.append(amino_acid)

    args.pansequence_min_amino_acids = [int(args.pansequence_min_amino_acids[0]),
                                        float(args.pansequence_min_amino_acids[1])]

    # Get gene frequency tables.
    if args.internal_genomes or args.external_genomes:
        multigenome_codon_usage = codonusage.MultiGenomeCodonUsage(args, run=run)
        frequency_df = multigenome_codon_usage.get_frequencies(
            from_function_sources=from_function_sources,
            return_functions=args.return_functions,
            return_amino_acids=args.return_amino_acids,
            function_accessions=args.function_accessions_dict,
            function_names=args.function_names_dict,
            expect_functions=args.expect_functions,
            relative=args.relative,
            synonymous=args.synonymous,
            sum_genes=args.sum,
            average_genes=args.average,
            gene_min_codons=args.gene_min_codons,
            function_min_codons=args.function_min_codons,
            min_codon_filter=args.min_codon_filter,
            drop_amino_acids=args.exclude_amino_acids,
            sequence_min_amino_acids=args.sequence_min_amino_acids,
            pansequence_min_amino_acids=args.pansequence_min_amino_acids,
            label_amino_acids=args.header_amino_acids,
            infinity_to_zero=args.infinity_to_zero)
    else:
        single_genome_codon_usage = codonusage.SingleGenomeCodonUsage(args, run=run)
        frequency_df = single_genome_codon_usage.get_frequencies(
            from_function_sources=from_function_sources,
            return_functions=args.return_functions,
            return_amino_acids=args.return_amino_acids,
            gene_caller_ids=args.gene_caller_ids,
            function_accessions=args.function_accessions_dict,
            function_names=args.function_names_dict,
            expect_functions=args.expect_functions,
            relative=args.relative,
            synonymous=args.synonymous,
            sum_genes=args.sum,
            average_genes=args.average,
            gene_min_codons=args.gene_min_codons,
            function_min_codons=args.function_min_codons,
            min_codon_filter=args.min_codon_filter,
            drop_amino_acids=args.exclude_amino_acids,
            sequence_min_amino_acids=args.sequence_min_amino_acids,
            pansequence_min_amino_acids=args.pansequence_min_amino_acids,
            label_amino_acids=args.header_amino_acids,
            infinity_to_zero=args.infinity_to_zero)

    # Write output tables.
    if args.sum or args.average:
        if args.gene_table_output:
            table_output = args.gene_table_output
        else:
            table_output = args.function_table_output
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


def parse_select_functions_table(args, from_function_sources):
    """Select genes for the frequency analysis from functions specified in a tabular file."""
    if not args.select_functions_txt:
        return from_function_sources

    select_functions_df = pd.read_csv(
        args.select_functions_txt, sep='\t', header=None, names=['source', 'accession', 'name'])
    select_functions_df = select_functions_df.fillna('')

    if args.function_sources is None:
        args.function_sources = []
        from_function_sources = []
    if args.function_sources == True:
        pass
    else:
        # Add the sources from the file to the list of all function sources to check for and load in
        # the genomes.
        for source in select_functions_df['source'].unique():
            if source not in args.function_sources:
                args.function_sources.append(source)
                from_function_sources.append(source)

    for row in select_functions_df.itertuples(index=False):
        if not row.source:
            raise ConfigError(
                "Each row of the select functions table must have a source entry, as functions are "
                "allowed to come from different function annotation sources.")

        if not row.accession and not row.name:
            raise ConfigError(
                "Each row of the select functions table must have either an accession or function "
                "entry to define the function being sought in the annotation source.")

        if row.source == 'KEGG_BRITE' and not row.name:
            raise ConfigError(
                "In the select functions table, KEGG BRITE categories must be named in the third "
                "column. Categories do not have accessions in most BRITE hierarchies, and "
                "hierarchy accessions are not useful here.")

        if row.accession:
            try:
                args.function_accessions_dict[row.source].append(row.accession)
            except KeyError:
                args.function_accessions_dict[row.source] = [row.accession]

        if row.name:
            try:
                args.function_names_dict[row.source].append(row.name)
            except KeyError:
                args.function_names_dict[row.source] = [row.name]

    return from_function_sources


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group(
        'SINGLE GENOME INPUTS',
        "Get CUB from genes or functions in a single genome. A contigs database can be "
        "provided alone as an \"external\" genome. An \"internal\" genome (bin) also requires a "
        "profile database, collection name, and bin ID.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupA.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupA.add_argument(
        '--gene-caller-ids', type=int, nargs='+', help="Select genes by ID, space-separated.")

    groupB = parser.add_argument_group(
        'MULTIPLE GENOME INPUTS',
        "Get frequencies from genes or functions in multiple genomes by providing internal and/or "
        "external genome files listing the genomes to analyze.")
    groupB.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupB.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))

    groupC = parser.add_argument_group(
        'OUTPUT FILES',
        "This program writes codon or amino acid frequency tables. When functions rather than "
        "genes are analyzed, one of two possible files can be produced: a table of frequencies per "
        "gene or per function, using `--gene-table-output` or `--function-table-output`, "
        "respectively. The per-function table can be derived from the per-gene table by summing "
        "genes with the same function annotation, which are found in columns 2-4.")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupC.add_argument(
        '--gene-table-output',
        help="A tab-delimited file of genes x codons or amino acids. The index columns before "
             "frequency data contain genome names (optional, if multiple genomes are considered) "
             "and gene callers IDs.")
    groupC.add_argument(
        '--function-table-output',
        help="A tab-delimited file of functions x codons or amino acids. Index columns before the "
             "frequency data contain, respectively, genome names (optional, if multiple genomes "
             "are considered), function annotation sources, accessions, and names.")
    groupC.add_argument(
        '--header-amino-acids', default=False, action='store_true',
        help="Include the amino acid for each codon in the column header of codon output, i.e., "
             "LysAAA instead of AAA.")
    groupC.add_argument(
        '--infinity-to-zero', default=False, action='store_true',
        help="Replace NA (empty) values in output with 0.0. NA occurs with `--synonymous` when all "
             "codons for an amino acid are absent in a gene or function, resulting in 0/0, "
             "reported as NA. Use with caution, for NA and 0.0 mean different things and this will "
             "skew downstream analyses of synonymous relative frequencies, such as codon usage "
             "bias.")

    groupD = parser.add_argument_group(
        'FREQUENCY STATISTICS',
        "Rather than absolute frequencies (default), relative frequencies or synonymous (per-amino "
        "acid) relative frequencies can be returned. Absolute and relative frequencies can be "
        "calculated for the sets of codons encoding the same amino acid rather than for individual "
        "codons. Frequencies can be summed or averaged across genes in a genome.")
    groupD.add_argument(
        '--relative', default=False, action='store_true',
        help="Return relative frequencies across codons or amino acids in each gene or function.")
    groupD.add_argument(
        '--synonymous', default=False, action='store_true',
        help="Return synonymous (per-amino acid) frequencies among the codons encoding each amino "
             "acid in each gene or function.")
    groupD.add_argument(
        '--return-amino-acids', default=False, action='store_true',
        help="Return frequencies of the sets of codons encoding the same amino acid.")
    groupD.add_argument(
        '--encodings-txt',
        help="Changes to the standard genetic code can be provided in this tab-delimited file of "
             "two columns. Each entry in the first column is a codon and each entry in the second "
             "column is the three-letter code for the decoded amino acid. For example, to recode "
             "the stop codon, TGA, as Trp, 'TGA' would be placed in the first column and 'Trp' in "
             "the same row of the second. Stop/termination codons are abbreviated \"STP\". This "
             "option affects per-amino acid and synonymous codon output (when using "
             "`--return-amino-acids` or `--synonymous`).")
    groupD.add_argument(
        '--sum', default=False, action='store_true',
        help="Sum frequencies across genes, returning a single row for each genome. When functions "
             "or function sources are selected, genes are subsetted to those annotated with the "
             "requested functions or annotated by the requested sources.")
    groupD.add_argument(
        '--average', default=False, action='store_true',
        help="Average frequencies across genes in the genome, returning a single row for each "
             "genome. When functions or function sources are selected, genes are subsetted to "
             "those annotated with the requested functions or annotated by the requested sources.")

    groupE = parser.add_argument_group(
        'FUNCTIONS',
        "Frequencies can be calculated for functions rather than genes, summing the frequencies of "
        "the genes annotated by each function. Genes can also be subsetted to those annotated with "
        "requested functions or annotated by requested function sources.")
    groupE.add_argument(
        '--function-sources', nargs='*',
        help="Return frequencies for functions annotated by these sources, e.g., 'KOfam', "
             "'KEGG_BRITE', 'COG20_FUNCTION'. When used with certain other options, such as "
             "`--sum`, rather than returning statistics for each function, analyzed genes are "
             "subsetted to those annotated by the provided sources. If `--function-sources` is "
             "used as a flag without any arguments, then every source will be considered.")
    groupE.add_argument(
        '--function-accessions', nargs='+',
        help="Return frequencies for functions with these accessions from the source provided in "
             "`--function-sources`. To get accessions from multiple sources, instead use "
             "`--select-functions-txt`.")
    groupE.add_argument(
        '--function-names', nargs='+',
        help="Return frequencies for functions with these names from the source provided in "
             "`--function-sources`. To get function names from multiple sources, instead "
             "use `--select-functions-txt`.")
    groupE.add_argument(
        '--select-functions-txt',
        help="Selected functions can be listed in this tab-delimited file of three columns. The "
             "first column should contain function annotation sources, the second column "
             "accessions, and the third function names. An entry in the source column is required "
             "in every row, and either an accession or name, or both, should also be in a row. "
             "The file should not have a header of column names.")
    groupE.add_argument(
        '--expect-functions',
        default=False,
        action='store_true',
        help="By default, functions provided by `--function-accessions`, `--function-names`, and "
             "`--select-functions-txt` need not be annotated in the input genomes. With this flag, "
             "an error will be raised if any of the functions are not present in an input genome.")
    groupE.add_argument(
        '--shared-function-sources', default=False, action='store_true',
        help="Use this flag to exclude function annotation sources that were not run on every "
             "input genome.")

    groupF = parser.add_argument_group(
        'FILTER GENES, FUNCTIONS, CODONS',
        "Genes/functions can be filtered by the number of codons they contain, e.g., ignore genes "
        "shorter than 300 codons. Codons can be selected a priori, e.g., ignore Ala codons, or "
        "rarer codons can be excluded, e.g., ignore amino acids that are decoded by ≤3 codons in "
        "≥90%% of genes. Filters can improve the statistical utility of codon relative frequency "
        "data.")
    groupF.add_argument(
        '--gene-min-codons', type=int, default=0,
        help="Set the minimum number of codons required in a gene. When functions are returned "
             "rather than genes, this filter is applied to genes before grouping them as "
             "functions.")
    groupF.add_argument(
        '--function-min-codons', type=int, default=0,
        help="Set the minimum number of codons required in a function. Genes with fewer than "
             "`--gene-min-codons` are first removed, and then functional groups of the remaining "
             "genes with fewer than `--function-min-codons` are removed. This filter only applies "
             "when returning functions.")
    groupF.add_argument(
        '--exclude-amino-acids', nargs='+',
        help="Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala, "
             "and STP for stop codons). If `--synonymous`, this argument defaults to \"STP Met "
             "Trp\", and if other amino acids are excluded, for STP, Met, and Trp codons to still "
             "be excluded from the output table, they must also be explicitly provided in the "
             "argument.")
    groupF.add_argument(
        '--include-amino-acids', nargs='+',
        help="This is the complement of `--exclude-amino-acids`. Only codons for the given amino "
             "acids are analyzed and reported.")
    groupF.add_argument(
        '--sequence-min-amino-acids', type=int, default=0,
        help="Do not report codons for amino acids (and STP) that are less numerous than the given "
             "argument. For example, if the argument is 5, and a gene or function query has 4 "
             "codons encoding Asn, 2 AAT and 2 AAC, then a row for this gene in the output table "
             "will have missing values in Asn columns. This filter occurs at the end of the "
             "analysis before writing results and so does not affect prior calculations.")
    groupF.add_argument(
        '--pansequence-min-amino-acids', nargs=2, default=[0, 1.0],
        help="The first value of the argument is a positive integer representing a minimum number "
             "of codons encoding an amino acid -- 'min_amino_acids' -- and the second is a number "
             "between 0 and 1 representing a fraction of genes -- 'min_gene_fraction'. Remove "
             "codons for amino acids (and STP) that are less numerous than 'min_amino_acids' in a "
             "'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and "
             "'min_gene_fraction' is 0.9, then if there are fewer than 5 codons for an amino "
             "acid/STP in ≥90%% of genes, then the columns for these codons are dropped.")
    groupF.add_argument(
        '--min-codon-filter', choices=['length', 'remaining', 'both'], default='both',
        help="This argument arises from the ambiguity of filters that remove genes and functions "
             "by number of codons (`--gene-min-codons` and `--function-min-codons`) in relation to "
             "filters that drop codons (`--exclude/include-amino-acids` and "
             "`--pansequence-min-amino-acids`). Genes (and functions) can be filtered by "
             "their full length, e.g., genes shorter than 300 codons are ignored. They can also be "
             "filtered by the number of codons remaining after dropping codons. The codon length "
             "filter followed by dropping codons can result in genes and functions with fewer "
             "codons than the original codon threshold -- thus the option of both \"length\" and "
             "\"remaining\" filters to ensure that total codon frequencies in the output always "
             "meet the minimum codon threshold. \"both\" is needed as an option in addition to "
             "\"remaining\" so dynamic codon filtering by `--pansequence-min-amino-acids` operates "
             "on genes that passed the first length filter.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
