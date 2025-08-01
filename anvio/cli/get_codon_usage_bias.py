#!/usr/bin/env python
# -*- coding: utf-8
"""Get codon usage bias (CUB) of genes and functions."""


import os
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
__authors__ = ['semiller10']
__requires__ = ['contigs-db',
                'profile-db',
                'collection',
                'bin',
                'internal-genomes',
                'external-genomes']
__provides__ = []
__description__ = ("Get codon usage bias (CUB) statistics of genes and functions")


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
    """Prepare arguments to get codon usage bias."""

    args = get_args()
    run = terminal.Run()

    if not args.output_file:
        raise ConfigError("`--output-file` must be provided.")
    filesnpaths.is_output_file_writable(args.output_file, ok_if_exists=False)

    # Store the coding dictionary in `args` for use in instantiating single- and multi-genome codon
    # usage objects.
    args.codon_to_amino_acid = codonusage.get_custom_encodings(args.encodings_txt)

    if args.function_sources is None:
        from_function_sources = False
    elif len(args.function_sources) != 1 and (args.function_accessions or args.function_names):
        raise ConfigError(
            "`--function-accessions` and `--function-names` require a single value for "
            "`--function-sources`. If select functions come from more than one source, use "
            "`--select-functions-txt`.")
    elif len(args.function_sources) == 0:
        from_function_sources = True
    elif len(args.function_sources) >= 1:
        # This includes 1) use of `--function-sources` as a flag, 2) ≥1 argument passed to
        # `--function-sources`, and 3) no use of `--function-sources` (None).
        from_function_sources = [source for source in args.function_sources]

    if (args.function_accessions or args.function_names) and args.select_functions_txt:
        raise ConfigError("`--function-accessions` or `--function-names` should not be used with "
                          "`args.select_functions_txt`.")

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

    # Non-default reference function definitions can be provided in an tabular file of functions.
    reference_function_accessions_dict, reference_function_names_dict = parse_select_reference_functions_table(args)

    # Ensure that the reference function sources are loaded in the GenomeCodonUsage object.
    if reference_function_accessions_dict or reference_function_names_dict:
        reference_function_sources = list(
            set(reference_function_accessions_dict.keys()).union(
                set(reference_function_names_dict.keys())))
        # `args.function_sources` can assume the value None, [] when used as a flag, or [...]. When
        # used as a flag, all function sources in the genome are loaded.
        if args.function_sources is None:
            args.function_sources = reference_function_sources
        elif len(args.function_sources) > 0:
            for source in reference_function_sources:
                if source not in args.function_sources:
                    args.function_sources.append(source)

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

    output_file_root, output_file_extension = os.path.splitext(args.output_file)

    # Get CUB tables. Separate tables are produced for each combination of input genome and CUB
    # metric.
    if args.internal_genomes or args.external_genomes:
        multigenome_codon_usage = codonusage.MultiGenomeCodonUsage(args, run=run)
        output_paths = []
        # The keys are CUB metrics in the returned dict. Write files as results for each genome are
        # yielded.
        for genome_name, cub_table_dict in multigenome_codon_usage.get_codon_usage_bias(
            metrics=args.metrics,
            from_function_sources=from_function_sources,
            function_accessions=args.function_accessions_dict,
            function_names=args.function_names_dict,
            expect_functions=args.expect_functions,
            omnibias=args.omnibias,
            reference_function_accessions=reference_function_accessions_dict,
            reference_function_names=reference_function_names_dict,
            expect_reference_functions=args.expect_reference_functions,
            reference_gene_caller_ids=args.reference_gene_caller_ids,
            gene_min_codons=args.gene_min_codons,
            function_min_codons=args.function_min_codons,
            min_codon_filter=args.min_codon_filter,
            drop_amino_acids=args.exclude_amino_acids,
            sequence_min_amino_acids=args.sequence_min_amino_acids,
            pansequence_min_amino_acids=args.pansequence_min_amino_acids,
            query_min_analyzed_codons=args.query_min_analyzed_codons,
            reference_exclude_amino_acid_count=args.reference_exclude_amino_acid_count,
            reference_min_analyzed_codons=args.reference_min_analyzed_codons):

            # Write output tables for the genome.
            if args.omnibias:
                # Output an omnibias CUB table for each genome + reference metric.
                referenceless_cub_dfs = []
                for metric, cub_df in cub_table_dict.items():
                    if metric in codonusage.reference_dependent_cub_metrics:
                        output_path = (output_file_root + '-' + genome_name + '-' + metric +
                                       output_file_extension)
                        output_paths.append(output_path)
                        cub_df.to_csv(output_path, sep='\t')
                    else:
                        referenceless_cub_dfs.append(cub_df)
                if referenceless_cub_dfs:
                    referenceless_cub_df = pd.concat(referenceless_cub_dfs, axis=1)
                    output_path = 'referenceless-' + args.output_file
                    referenceless_cub_df.to_csv(output_path, sep='\t')
                    output_paths.append(output_path)
                run.info("CUB table output directory", os.path.dirname(output_path))
                run.info(
                    "CUB table output files",
                    ', '.join([os.path.basename(output_path) for output_path in output_paths]),
                    nl_after=2)
            else:
                # Output a CUB table for the genome, with a column for each metric.
                cub_df = pd.concat(cub_table_dict.values(), axis=1)
                output_path = output_file_root + '-' + genome_name + output_file_extension
                output_paths.append(output_path)
                cub_df.to_csv(output_path, sep='\t')
                run.info("CUB table output", output_path, nl_after=2)
    else:
        single_genome_codon_usage = codonusage.SingleGenomeCodonUsage(args, run=run)
        # The keys are CUB metrics in the returned dict.
        cub_table_dict = single_genome_codon_usage.get_codon_usage_bias(
            metrics=args.metrics,
            from_function_sources=from_function_sources,
            gene_caller_ids=args.gene_caller_ids,
            function_accessions=args.function_accessions_dict,
            function_names=args.function_names_dict,
            expect_functions=args.expect_functions,
            omnibias=args.omnibias,
            reference_function_accessions=reference_function_accessions_dict,
            reference_function_names=reference_function_names_dict,
            expect_reference_functions=args.expect_reference_functions,
            reference_gene_caller_ids=args.reference_gene_caller_ids,
            gene_min_codons=args.gene_min_codons,
            function_min_codons=args.function_min_codons,
            min_codon_filter=args.min_codon_filter,
            drop_amino_acids=args.exclude_amino_acids,
            sequence_min_amino_acids=args.sequence_min_amino_acids,
            pansequence_min_amino_acids=args.pansequence_min_amino_acids,
            query_min_analyzed_codons=args.query_min_analyzed_codons,
            reference_exclude_amino_acid_count=args.reference_exclude_amino_acid_count,
            reference_min_analyzed_codons=args.reference_min_analyzed_codons)

        # Single genome, one or more CUB metrics.
        if args.omnibias:
            # Output an omnibias CUB table for each metric.
            output_paths = []
            referenceless_cub_dfs = []
            for metric, cub_df in cub_table_dict.items():
                if metric in codonusage.reference_dependent_cub_metrics:
                    output_path = output_file_root + '-' + metric + output_file_extension
                    output_paths.append(output_path)
                    cub_df.to_csv(output_path, sep='\t')
                else:
                    referenceless_cub_dfs.append(cub_df)
            if referenceless_cub_dfs:
                referenceless_cub_df = pd.concat(referenceless_cub_dfs, axis=1)
                output_path = 'referenceless-' + args.output_file
                referenceless_cub_df.to_csv(output_path, sep='\t')
                output_paths.append(output_path)
            run.info("CUB table output directory", os.path.dirname(output_path))
            run.info("CUB table output files",
                     ', '.join([os.path.basename(output_path) for output_path in output_paths]))
        else:
            # Output a single CUB table with a column for each metric.
            cub_df = pd.concat(cub_table_dict.values(), axis=1)
            cub_df.to_csv(args.output_file, sep='\t')
            run.info("CUB table output", args.output_file)


def parse_select_functions_table(args, from_function_sources):
    """Select genes for the CUB analysis from functions specified in a tabular file."""
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


def parse_select_reference_functions_table(args):
    """Non-default reference genes for CUB analysis can be retrieved from a tabular file of
    functions."""
    reference_function_accessions_dict = {}
    reference_function_names_dict = {}

    if not args.select_reference_functions_txt:
        if not args.omnibias:
            if codonusage.default_reference_function_accessions:
                reference_function_accessions_dict[codonusage.default_reference_function_source] = \
                    codonusage.default_reference_function_accessions
            if codonusage.default_reference_function_names:
                reference_function_names_dict[codonusage.default_reference_function_source] = \
                    codonusage.default_reference_function_names
        return reference_function_accessions_dict, reference_function_names_dict

    select_reference_functions_df = pd.read_csv(
        args.select_reference_functions_txt,
        sep='\t',
        header=None,
        names=['source', 'accession', 'name'])
    select_reference_functions_df = select_reference_functions_df.fillna('')

    for row in select_reference_functions_df.itertuples(index=False):
        if not row.source:
            raise ConfigError(
                "Each row of the select reference functions table must have a source entry, as "
                "functions are allowed to come from different function annotation sources.")

        if not row.accession and not row.name:
            raise ConfigError(
                "Each row of the select reference functions table must have either an accession or "
                "function entry to define the function being sought in the annotation source.")

        if row.source == 'KEGG_BRITE' and not row.name:
            raise ConfigError(
                "In the select reference functions table, KEGG BRITE categories must be named in "
                "the third column. Categories do not have accessions in most BRITE hierarchies, "
                "and hierarchy accessions are not useful here.")

        if row.accession:
            try:
                reference_function_accessions_dict[row.source].append(row.accession)
            except KeyError:
                reference_function_accessions_dict[row.source] = [row.accession]

        if row.name:
            try:
                reference_function_names_dict[row.source].append(row.name)
            except KeyError:
                reference_function_names_dict[row.source] = [row.name]

    return reference_function_accessions_dict, reference_function_names_dict


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group(
        'SINGLE GENOME INPUTS',
        "Get CUB of \"query\" genes or functions in a single genome. A contigs database can be "
        "provided alone as an \"external\" genome. An \"internal\" genome (bin) also requires a "
        "profile database, collection name, and bin ID.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupA.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupA.add_argument(
        '--gene-caller-ids', type=int, nargs='+', help="Return CUB for genes with these IDs.")

    groupB = parser.add_argument_group(
        'MULTIPLE GENOME INPUTS',
        "Get CUB from genes or functions in multiple genomes by providing internal and/or "
        "external genome files listing the genomes to analyze.")
    groupB.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupB.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))

    groupC = parser.add_argument_group('OUTPUT')
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'help':
        "This program writes one or more tab-delimited files of gene/function CUB statistics. All "
        "tables have gene/function row labels. 1.a. Single genome, single CUB metric, not omnibias "
        "mode: A table with one column of CUB data is written to the provided output file path. "
        "1.b. Single genome, single CUB metric, omnibias mode: A table of gene/function x "
        "gene/function is written to the provided output file path. 2.a. Single genome, multiple "
        "CUB metrics, not omnibias mode: A table with a column of CUB data per metric is written "
        "to the provided output file path. 2.b. Single genome, multiple CUB metrics, omnibias "
        "mode: A table of gene/function x gene/function is written per reference-dependent CUB "
        "metric, and a single table with a column of CUB data per reference-independent metric is "
        "written to another file, with the provided output file path serving as a template for "
        "the new file paths. 3.a. Multiple genomes, single CUB metric, not omnibias mode: A table "
        "with one column of CUB data is written per genome, with the provided output file path "
        "serving as a template for the per-genome file paths. 3.b. Multiple genomes, single CUB "
        "metric, omnibias mode: A table of gene/function x gene/function is written per genome, "
        "with the provided output file path serving as a template for the per-genome file paths. "
        "4.a. Multiple genomes, multiple CUB metrics, not omnibias mode: A table with a column of "
        "CUB data per metric is written per genome, with the provided output file path serving as "
        "a template for the per-genome file paths. 4.b. Multiple genomes, multiple CUB metrics, "
        "omnibias mode: A table of gene/function x gene/function is written per genome x CUB "
        "metric, and a table with a column of CUB data per reference-independent metric is written "
        "per genome, with the provided output file path serving as a template for the new file "
        "paths."}))

    groupD = parser.add_argument_group('CUB_CALCULATION')
    groupD.add_argument(
        '--metrics', choices=['cai', 'delta'], nargs='*',
        help="These are the metrics used in CUB calculations, with a separate calculation and a "
             "separate section of rows in the output table for each metric provided. When used as "
             "a flag, `--metrics` returns every possible metric. Some CUB metrics depend on "
             "definition of a reference codon composition, whereas others are "
             "reference-independent. 'cai' is the reference-dependent Codon Adaptation Index of "
             "Sharp and Li (1987). 'delta' is a reference-dependent metric of Ran and Higgs (2012, "
             "Eq. 6) that involves the overall codon composition of the genome in the comparison.")
    groupD.add_argument(
        '--encodings-txt',
        help="Changes to the standard genetic code can be provided in this tab-delimited file of "
             "two columns. Each entry in the first column is a codon and each entry in the second "
             "column is the three-letter code for the decoded amino acid. For example, to recode "
             "the stop codon, TGA, as Trp, 'TGA' would be placed in the first column and 'Trp' in "
             "the same row of the second. Stop/termination codons are abbreviated \"STP\". This "
             "option affects per-amino acid and synonymous codon output (when using "
             "`--return-amino-acids` or `--synonymous`).")

    groupE = parser.add_argument_group(
        'FUNCTION SELECTION',
        "Calculate CUB for functions rather than genes. The codon frequencies of genes annotated "
        "by a function are summed, treating functions as concatenations of genes, with longer "
        "genes contributing more than shorter genes to function codon usage. Alternative methods "
        "of combining gene codon frequencies, such as averaging genes, would not treat every gene "
        "equally. Averaging means that short genes have equal importance as long genes in the CUB "
        "of a function, which could imply that short genes are rate-limiting and therefore need to "
        "be favored for translation.")
    groupE.add_argument(
        '--function-sources', nargs='*',
        help="Return CUB for functions annotated by these sources, e.g., 'KOfam', 'KEGG_BRITE', "
             "'COG20_FUNCTION'. If `--function-sources` is used as a flag without any arguments, "
             "then every source will be considered.")
    groupE.add_argument(
        '--function-accessions', nargs='+',
        help="Return CUB for functions with these accessions from the source provided in "
             "`--function-sources`. To get accessions from multiple sources, instead use "
             "`--select-functions-txt`.")
    groupE.add_argument(
        '--function-names', nargs='+',
        help="Return CUB for functions with these names from the source provided in "
             "`--function-sources`. To get function names from multiple sources, instead "
             "use `--select-functions-txt`.")
    groupE.add_argument(
        '--select-functions-txt',
        help="Functions to select can be listed in this tab-delimited file of three columns. The "
             "first column should contain function annotation sources, the second column "
             "accessions, and the third function names. An entry in the source column is required "
             "in every row, and either an accession or name, or both, should also be in a row. "
             "Column names should not be provided.")
    groupE.add_argument(
        '--expect-functions',
        default=False,
        action='store_true',
        help="By default, functions provided by `--function-accessions`, `--function-names`, and "
             "`--select-functions-txt` need not be annotated in the input genomes. With this flag, "
             "an error will be raised if any of the functions are not present.")
    groupE.add_argument(
        '--shared-function-sources', default=False, action='store_true',
        help="When analyzing multiple genomes, exclude function annotation sources that were not "
             "run on every genome.")

    groupF = parser.add_argument_group(
        'REFERENCE',
        "Certain CUB metrics rely on definition of a reference codon composition. Often, a set of "
        "housekeeping genes including ribosomal proteins is used as the reference. Reference "
        "codon frequencies are the sums of the reference gene codon frequencies (long genes "
        "contribute more to reference codon frequencies than short genes). By default, reference "
        "genes are ribosomal proteins annotated by KEGG KOfams/BRITE, requiring input contigs "
        "databases to have been annotated by `anvi-run-kegg-kofams`. Specifically, the genes "
        "must be classified within BRITE as 'Ribosome>>>Ribosomal proteins'.")
    groupF.add_argument(
        '--omnibias', default=False, action='store_true',
        help="Use every gene or function as a separate reference rather than defining a set of "
             "reference genes or functions. The resulting table of gene x gene (or function x "
             "function) CUB values is like a distance matrix of the similarity of gene codon "
             "compositions.")
    groupF.add_argument(
        '--select-reference-functions-txt',
        help="Rather than using the default reference, selected reference functions can be listed "
             "in a headerless, tab-delimited file of three columns. Each row is a selected "
             "annotation. The first column must contain function annotation sources, the second "
             "column should contain function accessions, and the third columns should contain "
             "function names. Functions can be selected by function name or accession, so a value "
             "can be provided for one or the other. In summary, every row requires a source in the "
             "first column and either an accession or name in the second and third columns, "
             "respectively.")
    groupF.add_argument(
        '--expect-reference-functions', default=False, action='store_true',
        help="By default, reference functions provided by `--select-reference-functions-txt` do "
             "not all need to be annotated in the input genomes. For example, if two KEGG BRITE "
             "categories are used as references -- 'Ribosome>>>Ribosomal proteins>>>Bacteria' and "
             "'Transcription machinery>>>Prokaryotic type>>>Bacterial type>>>RNA polymerase' -- "
             "and genes are found to be annotated by only one of the categories, then these genes "
             "will be used as the reference. If this flag is used, then an error would instead "
             "be raised since one of the reference functions is not represented in the genome.")
    groupF.add_argument(
        '--reference-gene-caller-ids', type=int, nargs='+',
        help="Include specific genes in the reference gene set, as given by their gene caller IDs "
             "in the contigs database. This option can only be used with a single genome -- a "
             "single contigs database, internal genome, or external genome.")

    groupG = parser.add_argument_group(
        'FILTER GENES, FUNCTIONS, CODONS',
        "Genes/functions can be filtered by the number of codons they contain, e.g., ignore genes "
        "shorter than 300 codons. Codons can be selected a priori, e.g., ignore Ala codons, or "
        "rarer codons can be excluded, e.g., ignore amino acids that are decoded by ≤3 codons in "
        "≥90%% of genes. Filters can improve the statistical utility of CUB data. These filters do "
        "not apply to the reference gene set.")
    groupG.add_argument(
        '--gene-min-codons', type=int, default=0,
        help="Set the minimum number of codons required in a gene forming a gene/function "
             "queries. When functions are returned rather than genes, this filter is applied to "
             "genes before grouping them as functions.")
    groupG.add_argument(
        '--function-min-codons', type=int, default=0,
        help="Set the minimum number of codons required in a function \"query\". Genes with fewer "
             "than `--gene-min-codons` are first removed, and then functional groups of the "
             "remaining genes with fewer than `--function-min-codons` are removed. This filter "
             "only applies when returning function queries.")
    groupG.add_argument(
        '--exclude-amino-acids', nargs='+', default=['STP'],
        help="Remove codons that decode the given amino acids (use three-letter codes, e.g., "
             "Ala). Amino acids encoded by a single codon (Met and Trp) are excluded perforce from "
             "CUB calculations and do not have to be made explicit here. Importantly, to continue "
             "to exclude stop codons, as per the default, make sure to include \"STP\" in the "
             "argument: exclusion of Ala and stop codons is achieved with \"Ala STP\".")
    groupG.add_argument(
        '--include-amino-acids', nargs='+',
        help="This is the complement of `--exclude-amino-acids`. Only codons for the given amino "
             "acids are analyzed.")
    groupG.add_argument(
        '--sequence-min-amino-acids', type=int, default=0,
        help="Remove codons for amino acids (and STP) that are less numerous than the given "
             "argument. For example, if the argument is 5, and a gene or function query has 4 "
             "codons encoding Asn, 2 AAT and 2 AAC, then Asn will be disregarded in the "
             "calculation of CUB for this query.")
    groupG.add_argument(
        '--pansequence-min-amino-acids', nargs=2, default=[0, 1.0],
        help="The first value of the argument is a positive integer representing a minimum number "
             "of codons encoding an amino acid -- 'min_amino_acids' -- and the second is a number "
             "between 0 and 1 representing a fraction of genes -- 'min_gene_fraction'. Remove "
             "codons for amino acids (and STP) that are less numerous than 'min_amino_acids' in a "
             "'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and "
             "'min_gene_fraction' is 0.9, then if there are fewer than 5 codons for an amino "
             "acid/STP in ≥90%% of genes, the amino acid's codons do not factor into the "
             "calculation of CUB for any query.")
    groupG.add_argument(
        '--min-codon-filter', choices=['length', 'remaining', 'both'], default='both',
        help="This argument arises from the ambiguity of filters that remove genes and functions "
             "by number of codons (`--gene-min-codons` and `--function-min-codons`) in relation to "
             "the filters that drop codons (`--exclude/include-amino-acids` and "
             "`--exclude-amino-acid-count/fraction`). Genes (and functions) can be filtered by "
             "their full length, e.g., genes shorter than 300 codons are ignored. They can also be "
             "filtered by the number of codons remaining after dropping codons. The codon length "
             "filter followed by dropping codons can result in genes and functions with fewer "
             "codons than the original codon threshold -- thus the option of both \"length\" and "
             "\"remaining\" filters to ensure that total codon frequencies always meet the minimum "
             "codon threshold. \"both\" is needed as an option in addition to \"remaining\" so "
             "dynamic codon filtering by `--exclude-amino-acid-count/fraction` operates on genes "
             "that passed the first length filter.")
    groupG.add_argument(
        '--query-min-analyzed-codons',
        type=int,
        default=codonusage.default_query_min_analyzed_codons,
        help="Only allow CUB to calculated for a query if it has at least this number of "
             "synonymous codons that will be analyzed. For reference-dependent CUB metrics, "
             "analyzed codons are those with reference compositions.")

    groupH = parser.add_argument_group(
        'FILTER REFERENCE',
        "Confidence in the reference codon composition can be increased by setting minimum "
        "frequency requirements for the set of reference genes or codons for any individual amino "
        "acid to be retained. These arguments should not be used with `--omnibias`. Queries and "
        "references are the same in omnibias mode, so references will be filtered by query options "
        "including `--exclude-amino-acid-count`, `--exclude-amino-acid-fraction`, and "
        "`--query-min-analyzed-codons`.")
    groupH.add_argument(
        '--reference-exclude-amino-acid-count',
        type=int,
        default=codonusage.default_reference_exclude_amino_acid_count,
        help="Exclude codons for amino acids with fewer than this many codons in the set of "
             "reference genes.")
    groupH.add_argument(
        '--reference-min-analyzed-codons',
        type=int,
        default=codonusage.default_reference_min_analyzed_codons,
        help="Only allow CUB to be calculated for a genome if the set of reference genes contains "
             "at least this many codons. This filter applies after excluding codons for individual "
             "amino acids using `--reference-exclude-amino-acid-count`.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
