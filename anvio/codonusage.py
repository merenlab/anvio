# -*- coding: utf-8
# pylint: disable=line-too-long
"""Codon usage analyses at the levels of genes, functional groups of genes, and genomes"""

import os
import copy
import inspect
import argparse
import numpy as np
import pandas as pd

from functools import partial
from collections import Counter

import anvio
from anvio import filesnpaths
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.ccollections as ccollections

from anvio.dbops import ContigsSuperclass
from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomedescriptions import GenomeDescriptions


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2022, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(verbose=False)

pp = terminal.pretty_print

default_codon_amino_acid_dict = {
    codon: amino_acid for codon, amino_acid in constants.codon_to_AA.items()}

# CAI is the Codon Adaptation Index of Sharp and Li (1987).
# Delta is from Ran and Higgs (2012, Eq. 6).
cub_metrics = ['cai', 'delta']
# The following CUB metrics rely upon comparison to a set of reference genes.
reference_dependent_cub_metrics = ['cai', 'delta']

# For CUB, remove single-codon amino acids (which do not affect CUB) and stop codons.
single_codon_amino_acids = ['Met', 'Trp']
ignored_cub_amino_acids = ['Met', 'Trp', 'STP']

# By default, genes annotated as ribosomal proteins by KEGG KOfams/BRITE are used as the CUB
# reference.
default_reference_function_source = 'KEGG_BRITE'
default_reference_function_accessions = []
default_reference_function_names = ['Ribosome>>>Ribosomal proteins']

default_query_min_analyzed_codons = 100
default_reference_exclude_amino_acid_count = 10
default_reference_min_analyzed_codons = 100


class SingleGenomeCodonUsage(object):
    """
    This object processes codon usage data for a single genome.

    Manipulate the raw data using the methods, `get_frequencies` and `get_codon_usage_bias`.
    """

    def __init__(self, args, r=run, rq=run_quiet, p=progress):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.contigs_db_path = A('contigs_db')

        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')

        self.gene_caller_ids_of_interest = A('gene_caller_ids')
        if self.gene_caller_ids_of_interest is None:
            self.gene_caller_ids_of_interest = set()

        self.function_sources = A('function_sources')
        self.all_brite_categories = A('all_brite_categories')

        self.ignore_start_codons = A('ignore_start_codons')
        if self.ignore_start_codons is None:
            self.ignore_start_codons = False

        self.run = r
        self.run_quiet = rq
        self.progress = p

        self._set_genetic_code()
        self._load_contigs_db_data()
        self._make_gene_codon_frequency_table()


    def _set_genetic_code(self):
        """
        Store decoding properties of the genome as object attributes.

        The dict, `args.codon_to_amino_acid`, should have keys that are codons and values that are
        three-letter amino acid codes ("STP" for stop codons). If `args.codon_to_amino_acid` is None, the
        standard genetic code is used.
        """
        if 'codon_to_amino_acid' not in self.args or self.args.codon_to_amino_acid is None:
            self.codon_amino_acid_dict = default_codon_amino_acid_dict
        else:
            check_genetic_code(self.args.codon_to_amino_acid)
            self.codon_amino_acid_dict = self.args.codon_to_amino_acid
            if self.codon_amino_acid_dict != default_codon_amino_acid_dict:
                self.run.info_single("Using a nonstandard genetic code for the genome")

        self.amino_acid_codons_dict = {}
        for codon, amino_acid in self.codon_amino_acid_dict.items():
            try:
                self.amino_acid_codons_dict[amino_acid].append(codon)
            except KeyError:
                self.amino_acid_codons_dict[amino_acid] = [codon]

        self.nonstop_amino_acid_codons_dict = {
            amino_acid: codons for amino_acid, codons in
            self.amino_acid_codons_dict.items() if amino_acid != 'STP'}

        # Ignore codons with "N" nucleotide which anvi'o records as None.
        self.ignored_cub_codons = [None]
        for amino_acid in ignored_cub_amino_acids:
            self.ignored_cub_codons += self.amino_acid_codons_dict[amino_acid]

        self.synonymous_nonstop_amino_acid_codons_dict = {
            amino_acid: codons for amino_acid, codons in self.nonstop_amino_acid_codons_dict.items()
            if amino_acid not in ignored_cub_amino_acids}


    def _load_contigs_db_data(self):
        """Load gene data from the contigs database."""
        utils.is_contigs_db(self.contigs_db_path)

        if self.profile_db_path or self.collection_name or self.bin_id:
            # Initialize the contigs superclass from the splits of the internal genome bin. An
            # exception will be raised if all three args are not provided or valid.
            self.args.split_names_of_interest = \
                ccollections.GetSplitNamesInBins(self.args).get_split_names_only()
        contigs_super = ContigsSuperclass(self.args, r=self.run_quiet)
        contigs_super.init_contig_sequences(
            gene_caller_ids_of_interest=self.gene_caller_ids_of_interest)
        self.contig_sequences_dict = contigs_super.contig_sequences

        if self.gene_caller_ids_of_interest:
            self.gene_caller_ids = list(set(self.gene_caller_ids_of_interest).intersection(
                contigs_super.gene_caller_ids_included_in_contig_sequences_initialized))
            self.genes_in_contigs_dict = {
                gene_caller_id: contigs_super.genes_in_contigs_dict[gene_caller_id]
                for gene_caller_id in self.gene_caller_ids}
        else:
            self.gene_caller_ids = \
                contigs_super.gene_caller_ids_included_in_contig_sequences_initialized
            self.genes_in_contigs_dict = contigs_super.genes_in_contigs_dict

        # Organize functional annotations into a table.
        contigs_super.init_functions(requested_sources=self.function_sources)
        if self.function_sources == []:
            self.function_sources = sorted(contigs_super.gene_function_call_sources)
            if not self.function_sources:
                self.run.warning(
                    "The value of `args.function_sources` is an empty list, indicating that all "
                    "function sources in the genome should by loaded. However, the contigs "
                    "database has not been annotated by any sources :/")

        if self.all_brite_categories and 'KEGG_BRITE' not in self.function_sources:
            raise ConfigError(
                "`all_brite_categories` can only be used with 'KEGG_BRITE' as as a function "
                "annotation source.")

        gene_function_rows = []
        for gene_caller_id, annotation_dict in contigs_super.gene_function_calls_dict.items():
            for annotation_source, annotation in annotation_dict.items():
                if annotation is None:
                    continue

                accession = annotation[0]
                name = annotation[1]
                accessions = accession.split('!!!')
                names = name.split('!!!')
                if len(accessions) == len(names):
                    # The function accession and name entries contain the same number of '!!!'
                    # separators.
                    for accession, name in zip(accessions, names):
                        if annotation_source == 'KEGG_BRITE':
                            hierarchy_accession = accession
                            categorization = name
                            if self.all_brite_categories:
                                # Include every possible depth of categorization.
                                split_categories = categorization.split('>>>')
                                for depth in range(1, len(split_categories) + 1):
                                    gene_function_rows.append(
                                        [annotation_source,
                                        hierarchy_accession,
                                        '>>>'.join(split_categories[: depth]),
                                        gene_caller_id])
                            else:
                                gene_function_rows.append(
                                    [annotation_source,
                                     hierarchy_accession,
                                     categorization,
                                     gene_caller_id])
                        else:
                            gene_function_rows.append(
                                [annotation_source, accession, name, gene_caller_id])
                else:
                    # The function accession and name entries do not contain the same number of
                    # '!!!' separators. In COG20_PATHWAY, there can be multiple accessions
                    # corresponding to the same function name.
                    if annotation_source == 'KEGG_BRITE':
                        hierarchy_accession = accession
                        categorization = name
                        if self.all_brite_categories:
                            # Include every possible depth of categorization.
                            split_categories = categorization.split('>>>')
                            for depth in range(1, len(split_categories) + 1):
                                gene_function_rows.append(
                                    [annotation_source,
                                    hierarchy_accession,
                                    '>>>'.join(split_categories[: depth]),
                                    gene_caller_id])
                        else:
                            gene_function_rows.append(
                                [annotation_source,
                                 hierarchy_accession,
                                 categorization,
                                 gene_caller_id])
                    else:
                        gene_function_rows.append(
                            [annotation_source, accession, name, gene_caller_id])
        self.gene_function_df = pd.DataFrame(
            gene_function_rows,
            columns=['function_source', 'function_accession', 'function_name', 'gene_caller_id'])


    def _make_gene_codon_frequency_table(self):
        """Generate the per-gene codon frequency DataFrame as `self.gene_codon_frequency_df`."""
        self.progress.new("Fetching codon frequency data")
        self.progress.update("...")

        gene_codon_frequencies = []
        skipped_noncoding_gene_caller_ids = []
        coding_gene_caller_ids = []
        for gene_caller_id in self.gene_caller_ids:
            # `gene_call` is a dictionary.
            gene_call = self.genes_in_contigs_dict[gene_caller_id]

            if gene_call['call_type'] != constants.gene_call_types['CODING']:
                skipped_noncoding_gene_caller_ids.append(gene_caller_id)
                continue

            coding_gene_caller_ids.append(gene_caller_id)

            gene_codons = utils.get_list_of_codons_for_gene_call(
                gene_call, self.contig_sequences_dict)
            if self.ignore_start_codons:
                gene_codons = gene_codons[1: ]
            gene_codon_frequencies.append(Counter(gene_codons))

        gene_codon_frequency_df = pd.DataFrame.from_records(gene_codon_frequencies)

        observed_codons = gene_codon_frequency_df.columns.tolist()
        for codon in constants.codon_to_AA:
            if codon not in observed_codons:
                gene_codon_frequency_df[codon] = 0

        # Drop any column named NaN for unknown codons.
        gene_codon_frequency_df = gene_codon_frequency_df[constants.codon_to_AA]

        gene_codon_frequency_df = gene_codon_frequency_df.fillna(0)
        gene_codon_frequency_df = gene_codon_frequency_df[sorted(gene_codon_frequency_df.columns)]
        gene_codon_frequency_df.index = coding_gene_caller_ids
        gene_codon_frequency_df.index.name = 'gene_caller_id'
        self.gene_codon_frequency_df = gene_codon_frequency_df

        self.progress.end()

        self.noncoding_gene_count = len(skipped_noncoding_gene_caller_ids)
        if self.noncoding_gene_count:
            self.run.warning(
                f"{pp(self.noncoding_gene_count)} of {pp(len(self.gene_caller_ids))} genes were "
                "non-coding and not added to the codon frequency table.")


    def get_frequencies(self,
                        from_function_sources=False,
                        return_functions=False,
                        return_amino_acids=False,
                        gene_caller_ids=None,
                        function_accessions=None,
                        function_names=None,
                        expect_functions=False,
                        relative=False,
                        synonymous=False,
                        sum_genes=False,
                        average_genes=False,
                        gene_min_codons=0,
                        function_min_codons=0,
                        min_codon_filter='both',
                        drop_amino_acids=None,
                        sequence_min_amino_acids=0,
                        pansequence_min_amino_acids=(0, 1.0),
                        label_amino_acids=False,
                        infinity_to_zero=False):
        """
        Get absolute (default) or relative codon or amino acid frequencies from genes or functions.

        Relative codon frequencies can be normalized per-amino acid to all synonymous codons.

        Parameters
        ==========
        from_function_sources : bool, str, or iterable of str, optional
            Select genes with functional annotations. With this argument, the first four columns of
            the returned table contain, respectively, function annotation sources, function
            accessions, function names, and gene caller IDs. There is a row for each gene/function
            combination in the table, and each row for the same gene contains the same frequency
            values. When this argument is True, use all available functional annotation sources in
            the SingleGenomeCodonUsage object. When this argument is a string, select the source
            given by the string, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. When this argument is an
            iterable, select a subset of sources given by its strings, e.g., ['KOfam',
            'COG20_FUNCTION']. By default None.
        return_functions : bool, optional
            If True (default False), output frequency tables contain function rather than gene
            results in each row. Returning per-gene results when also considering functions by using
            `--from-function-sources` facilitates analysis of per-gene variability within functions.
        return_amino_acids : bool, optional
            If True (default False), output frequency table columns are decoded amino acids (plus
            STP) rather than codons. Synonymous codon frequencies are summed to produce the amino
            acid frequencies.
        gene_caller_ids : iterable, optional
            Genes with the given IDs are selected for analysis. This parameter can be used alongside
            `function_accessions` and `function_names`. By default None.
        function_accessions : dict, optional
            Genes annotated with the given function accessions are selected for analysis. The
            argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
            parameter can be used alongside `gene_caller_ids` and `function_names`.  Note that
            'KEGG_BRITE' does not use individual function accessions but overarching hierarchy
            accessions that include multiple functions. By default None.
        function_names : dict, optional
            Genes annotated with the given function names are selected for analysis. The argument
            must be a dict keyed by function annotation source (e.g., 'KOfam', 'COG20_FUNCTION',
            'Pfam') and with values being lists of function names. Unlike `function_accessions`,
            'KEGG_BRITE' may be used as a source, with names (categorizations) given to an arbitrary
            level of the hierarchy, such as ['Ribosome>>>Ribosomal proteins', 'Transcription
            machinery>>>Prokaryotic type>>>Bacterial type>>>RNA polymerase']. All BRITE
            categorizations falling under the given level are selected. This parameter can be used
            alongside `gene_caller_ids` and `function_accessions`. By default None.
        expect_functions : bool, optional
            If True (default False), an error will be raised if any given `function_accessions` or
            `function_names` are not annotated in the input genome.
        relative : bool, optional
            If True (default False), return relative rather than absolute codon frequencies.
        synonymous : bool, optional
            If True (default False), return codon relative frequencies among synonymous codons
            decoding each amino acid (plus stop).
        sum_genes : bool, optional
            If True (default False), sum codon frequencies of genes, returning a one-row DataFrame
            of the summed frequencies. If `from_function_sources` is used, then genes are limited to
            those with the specified functional annotations. Relative/synonymous frequencies are
            calculated after summing absolute frequencies.
        average_genes : bool, optional
            If True (default False), average codon frequencies of genes, returning a one-row
            DataFrame of the averaged frequencies. If `from_function_sources` is used, then genes
            are limited to those with the specified functional annotations. Averaging occurs after
            calculation of relative/synonymous frequencies.

        Additional Parameters
        =====================
        These parameters filter genes and codons. Here is the order of all possible filters:
        gene codon frequency table ->
            drop genes on codon length ->
            drop gene/function pairs on the total length of remaining genes defining the function ->
            drop codons of defined amino acids ->
            dynamically drop codons of rarer amino acids ->
            drop genes on remaining codon frequency ->
            drop gene/function pairs on the codon sum of remaining genes defining the function ->
        filtered gene codon frequency table

        gene_min_codons : int, optional
            Ignore genes with fewer than this number of codons. By default 0.
        function_min_codons : int, optional
            Ignore gene/function pairs when the function has fewer than this number of codons. Genes
            are filtered by codon count before functions, so the total length of remaining genes
            annotated by the function is considered. By default 0.
        min_codon_filter : {"length", "remaining", "both"}, optional
            This argument arises from the ambiguity of filters that remove genes and functions by
            number of codons (`--gene-min-codons` and `--function-min-codons`) in relation to the
            filters that drop codons (`--exclude/include-amino-acids` and
            `--exclude-amino-acid-count/fraction`). Genes (and functions) can be filtered by their
            full length, e.g., genes shorter than 300 codons are ignored. They can also be
            filtered by the number of codons remaining after dropping codons. The codon length
            filter followed by dropping codons can result in genes and functions with fewer codons
            than the original codon threshold -- thus the option of both 'length' and 'remaining'
            filters to ensure that total codon frequencies in the output always meet the minimum
            codon threshold. 'both' is needed as an option in addition to 'remaining' so dynamic
            codon filtering by `--exclude-amino-acid-count/fraction` operates on genes that passed
            the first length filter. By default "both".
        drop_amino_acids : iterable, optional
            Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala, and
            STP for stop codons). If `synonymous` is True, the `drop_amino_acids` default rather
            than being None is STP plus amino acids encoded by a single codon (Met, Trp).
        sequence_min_amino_acids : int, optional
            Remove from the output rows codons for amino acids (and STP), or amino acids themselves,
            that are less numerous than `sequence_min_amino_acids`. For example, if the argument is
            5, and a gene has 4 codons encoding Asn, 2 AAT and 2 AAC, then a row for this gene in
            the output table will have missing values in Asn columns. By default 0.
        pansequence_min_amino_acids : tuple, optional
            This tuple must have two items. The first is an int >0 representing a minimum number of
            codons encoding an amino acid -- 'min_amino_acids' -- and the second is a float in the
            range (0, 1) representing a fraction of genes -- 'min_gene_fraction'. Remove codons for
            amino acids (and STP) that are less numerous than 'min_amino_acids' in a
            'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and
            'min_gene_fraction' is 0.9, then if there are fewer than 5 codons for amino acid/STP in
            ≥90% of genes, then the columns for these codons are dropped. By default (0, 1.0).
        label_amino_acids : bool, optional
            If True (default False), include the amino acid for each codon in the column header of
            the output, i.e., LysAAA instead of AAA.
        infinity_to_zero : bool, optional
            If True (default False), replace NA (empty) values in output with 0.0. NA occurs if
            `synonymous` is True and all codons for an amino acid are absent in a gene or function,
            resulting in 0/0, reported as NA. Use with caution, for NA and 0.0 mean different things
            and this will skew downstream analyses of synonymous relative frequencies, such as codon
            usage bias.

        Returns
        =======
        pandas.core.frame.DataFrame
            Frequency table of gene x codon or amino acid. If functions are not considered, then the
            Index of the returned DataFrame contains gene caller IDs. If functions are considered
            and `return_functions` is False, then each row represents a gene/function pair, and the
            same gene frequencies can be found in multiple rows with different function pairs; there
            are additional MultiIndex columns for source, accession, and name of function. If
            functions are considered and `return_functions` is True, then each row represents a
            function rather than gene/function pair. If frequencies are summed or averaged and
            functions are not considered, the returned single-row DataFrame has an Index with one
            entry, 'all'. If frequencies are summed or averaged and functions are considered, the
            returned DataFrame has a row for each function annotation source.

        Examples
        ========
        Return the amino acid frequencies of each gene annotated by KEGG KOfams.
        >>> self.get_frequencies(from_function_sources='KOfam', return_amino_acids=True)

        Return the summed amino acid frequencies of genes annotated by KEGG KOfams.
        >>> self.get_frequencies(
            from_function_sources='KOfam', return_amino_acids=True, sum_genes=True)

        Return the codon relative frequencies of each gene.
        >>> self.get_frequencies(relative=True)

        Return the relative frequencies of genes ≥300 codons in length.
        >>> self.get_frequencies(relative=True, gene_min_codons=300)

        Return the average relative frequencies of genes ≥300 codons in length.
        >>> self.get_frequencies(relative=True, average_all=True, gene_min_codons=300)

        Return the synonymous (per-amino acid) relative frequencies of genes. Remove genes <300
        codons, then remove amino acids with <5 codons in ≥90% of remaining genes, then remove genes
        with <300 codons remaining.
        >>> self.get_frequencies(
            synonymous=True, gene_min_codons=300, pansequence_min_amino_acids=(5, 0.9))

        Return the synonymous (per-amino acid) relative frequencies of KEGG KOfam and BRITE
        functions. This command removes genes <300 codons, then removes amino acids with <5 codons
        in ≥90% of remaining genes, then removes genes with <300 codons remaining, then sums gene
        frequencies in functions, then converts function codon frequencies to missing values for
        amino acids with <5 codons in the function, then calculates synonymous relative frequencies
        in functions.
        >>> self.get_frequencies(
            from_function_sources=['KOfam', 'KEGG_BRITE'],
            return_functions=True,
            synonymous=True,
            gene_min_codons=300,
            sequence_min_amino_acids=5,
            pansequence_min_amino_acids=(5, 0.9))
        """
        # CHECK ARGUMENTS AND SET UP PROCEDURE
        ######################################
        function_sources = self._establish_function_sources(from_function_sources)
        if function_sources == self.function_sources:
            gene_function_df = self.gene_function_df
        elif function_sources:
            # Subset the gene function table to the requested sources.
            gene_function_df = self.gene_function_df.set_index(
                'function_source').loc[function_sources].reset_index()

        # Check compatability of `return_amino_acids` with other arguments.
        if return_amino_acids and synonymous:
            raise ConfigError(
                "The argument `synonymous` should only be True when `return_amino_acids` is also "
                "True, as `synonymous` returns synonymous codon relative frequencies.")
        if return_amino_acids and label_amino_acids:
            # Don't bother raising an exception.
            label_amino_acids = False

        if gene_caller_ids is None:
            gene_caller_ids = []
        if function_accessions is None:
            function_accessions = {}
        if function_names is None:
            function_names = {}
        gene_subsetting = bool(gene_caller_ids or function_accessions or function_names)
        gene_codon_frequency_df = self._select_genes(
            gene_caller_ids, function_accessions, function_names, expect_functions)

        if from_function_sources and gene_subsetting:
            # Subset gene function table to genes of interest.
            gene_function_df = gene_function_df.set_index('gene_caller_id')
            gene_function_df = gene_function_df.loc[
                gene_function_df.index.intersection(gene_codon_frequency_df.index)].reset_index()

        if synonymous and not relative:
            relative = True

        if sum_genes and average_genes:
            raise ConfigError("`sum_genes` and `average_genes` cannot both be True.")

        if min_codon_filter not in ['length', 'remaining', 'both']:
            raise ConfigError("`min_codon_filter` must be one of 'length', 'remaining', or 'both'.")
        # Set `filter_gene_length`, `filter_function_length`, and `filter_remaining`.
        if min_codon_filter == 'length':
            filter_gene_length = bool(gene_min_codons)
            filter_function_length = bool(function_min_codons)
            filter_gene_remaining_codons = False
            filter_function_remaining_codons = False
        elif min_codon_filter == 'remaining':
            filter_gene_length = False
            filter_function_length = False
            filter_gene_remaining_codons = bool(gene_min_codons)
            filter_function_remaining_codons = bool(function_min_codons)
        elif min_codon_filter == 'both':
            filter_gene_length = bool(gene_min_codons)
            filter_function_length = bool(function_min_codons)
            filter_gene_remaining_codons = bool(gene_min_codons)
            filter_function_remaining_codons = bool(function_min_codons)

        # Set `drop_amino_acids`.
        if drop_amino_acids is None:
            if synonymous:
                drop_amino_acids = ignored_cub_amino_acids
            else:
                drop_amino_acids = []
        else:
            unrecognized_amino_acids = []
            for amino_acid in drop_amino_acids:
                if amino_acid not in self.amino_acid_codons_dict:
                    unrecognized_amino_acids.append(amino_acid)
            if unrecognized_amino_acids:
                raise ConfigError("The following amino acids in `drop_amino_acids` are not "
                                  f"recognized: {', '.join(unrecognized_amino_acids)}")

        if (type(pansequence_min_amino_acids[0]) != int or
            pansequence_min_amino_acids[0] < 0 or
            type(pansequence_min_amino_acids[1]) != float or
            not (0 < pansequence_min_amino_acids[1] <= 1)):
            raise ConfigError(
                "The value of `pansequence_min_amino_acids` must be a tuple with two items: the "
                "first a positive int and the second a float between 0 and 1.")

        if gene_subsetting:
            self.run.info_single(f"{pp(len(gene_codon_frequency_df))} of "
                                 f"{pp(len(self.gene_caller_ids) - self.noncoding_gene_count)} CDS "
                                 "selected from the genome")
        else:
            self.run.info_single(f"{pp(len(gene_codon_frequency_df))} CDS in the genome")

        # FILTER GENES/FUNCTIONS AND CODONS USING ADDITIONAL PARAMETERS
        ###############################################################
        # Filter genes by length.
        gene_codon_frequency_df = self._get_frequency_table(
            gene_codon_frequency_df,
            min_codons=gene_min_codons,
            filter_input_codon_count=filter_gene_length)

        # Filter functions by the total length of their genes. Gene/function pairs with fewer than
        # the required number of codons in the function are removed.
        if function_sources:
            gene_function_codon_frequency_df = gene_function_df.merge(
                gene_codon_frequency_df, how='inner', on='gene_caller_id')
            gene_function_codon_frequency_df = gene_function_codon_frequency_df.set_index(
                ['function_source', 'function_accession', 'function_name', 'gene_caller_id'])
            if filter_function_length:
                gene_function_codon_frequency_df = gene_function_codon_frequency_df.groupby(
                    ['function_source', 'function_accession', 'function_name']).filter(
                        lambda function_df: function_df.sum(axis=1).sum() >= function_min_codons)

        # Drop certain codon columns from the gene codon frequency table. Filter genes by total
        # codons remaining.
        gene_codon_frequency_df = self._get_frequency_table(
            gene_codon_frequency_df,
            drop_amino_acids=drop_amino_acids,
            pansequence_min_amino_acids=pansequence_min_amino_acids,
            min_codons=gene_min_codons,
            filter_output_codon_count=filter_gene_remaining_codons)

        need_to_filter_codons_in_gene_function_codon_frequency_table = True
        if ((drop_amino_acids or pansequence_min_amino_acids) and
            (function_min_codons and filter_function_remaining_codons)):
            # Filter functions by total codons remaining. Gene/function pairs with fewer than the
            # required number of codons in the function are removed.
            if function_sources:
                gene_function_codon_frequency_df = gene_function_df.merge(
                    gene_codon_frequency_df, how='inner', on='gene_caller_id')
                gene_function_codon_frequency_df = gene_function_codon_frequency_df.set_index(
                    ['function_source', 'function_accession', 'function_name', 'gene_caller_id'])
                if filter_function_remaining_codons:
                    gene_function_codon_frequency_df = gene_function_codon_frequency_df.groupby(
                        ['function_source', 'function_accession', 'function_name']).filter(
                            lambda function_df:
                                function_df.sum(axis=1).sum() >= function_min_codons)
                need_to_filter_codons_in_gene_function_codon_frequency_table = False

        if function_sources and need_to_filter_codons_in_gene_function_codon_frequency_table:
            gene_function_codon_frequency_df = self._get_frequency_table(
                gene_function_codon_frequency_df,
                drop_amino_acids=drop_amino_acids,
                pansequence_min_amino_acids=pansequence_min_amino_acids)

        if gene_min_codons:
            self.run.info_single(f"{pp(len(gene_codon_frequency_df))} CDS remaining after codon "
                                 "count filters")

        if pansequence_min_amino_acids[0] > 0 and pansequence_min_amino_acids[1] < 1:
            if gene_min_codons and min_codon_filter != 'remaining':
                min_gene_length_message = '≥' + pp(str(gene_min_codons)) + ' codon'
            else:
                min_gene_length_message = ''
            dynamically_dropped_amino_acids = set()
            for codon, amino_acid in self.codon_amino_acid_dict.items():
                if (codon not in gene_codon_frequency_df.columns and
                    amino_acid not in drop_amino_acids):
                    dynamically_dropped_amino_acids.add(amino_acid)
            if dynamically_dropped_amino_acids:
                self.run.warning(
                    "Codons for the following amino acids were dropped as they did not meet the "
                    f"threshold of {pp(pansequence_min_amino_acids[0])} codons in "
                    f"{pansequence_min_amino_acids[1] * 100}% of {min_gene_length_message} CDS: "
                    f"{', '.join(sorted(dynamically_dropped_amino_acids))}")

        # GET OUTPUTS WITH NO CONSIDERATION OF FUNCTIONS
        ################################################
        get_table = lambda method: method(gene_codon_frequency_df,
                                          sequence_min_amino_acids=sequence_min_amino_acids,
                                          label_amino_acids=label_amino_acids,
                                          replace_na=infinity_to_zero,
                                          output_amino_acids=return_amino_acids)
        ### Absolute frequencies
        if (not relative and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_frequency_table)
        ### Relative frequencies
        if (relative and
            not synonymous and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)
        ### Synonymous (per-amino acid) relative frequencies
        if (not return_amino_acids and
            relative and
            synonymous and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_synonymous_codon_rel_frequency_table)

        # Summed codon/amino acid results across genes:
        ### Absolute frequencies
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df),
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids,
            output_amino_acids=return_amino_acids)
        if (not relative and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_frequency_table)
        ### Relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df,
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids,
            output_amino_acids=return_amino_acids)
        if (relative and
            not synonymous and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df,
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids,
            replace_na=infinity_to_zero)
        if (not return_amino_acids and
            relative and
            synonymous and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

        # Codon results averaged across genes:
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df),
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids)
        ### Absolute frequencies
        if (not return_amino_acids and
            not relative and
            not function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df,
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids,
            replace_na=infinity_to_zero)
        if (not return_amino_acids and
            relative and
            not synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Synonymous relative frequencies
        if (not return_amino_acids and
            relative and
            synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_synonymous_codon_rel_frequency_table)

        # Amino acid results averaged across genes:
        # Gene amino acid frequencies are first calculated from codon frequencies before averaging
        # across genes.
        ### Absolute frequencies
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df),
            sequence_min_amino_acids=sequence_min_amino_acids,
            output_amino_acids=return_amino_acids)
        if (return_amino_acids and
            not relative and
            not function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (return_amino_acids and
            relative and
            not synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)

        # GET OUTPUTS WITH CONSIDERATION OF FUNCTIONS
        #############################################
        get_table = lambda method: method(
            (gene_function_codon_frequency_df.groupby(
                ['function_source', 'function_accession', 'function_name']).sum()
             if return_functions else
             gene_function_codon_frequency_df.sort_values(
                 ['function_source', 'function_accession', 'function_name'])),
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids,
            replace_na=infinity_to_zero,
            output_amino_acids=return_amino_acids)
        ### Absolute frequencies
        if (not relative and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_frequency_table)
        ### Relative frequencies
        if (relative and
            not synonymous and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)
        ### Synonymous relative frequencies
        if (not return_amino_acids and
            relative and
            synonymous and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_synonymous_codon_rel_frequency_table)

        # Remove duplicate occurrences of genes before summing or averaging frequencies of all genes
        # in the function source. A gene can have different annotations, e.g., different KOfam
        # assignments. In KEGG BRITE, a gene can be in nested categories and different hierarchies.
        gene_function_codon_frequency_df = \
            gene_function_codon_frequency_df.reset_index().groupby('function_source').apply(
                lambda source_df: source_df.drop_duplicates(
                    subset='gene_caller_id', ignore_index=True)).set_index(
                        ['function_source',
                         'function_accession',
                         'function_name',
                         'gene_caller_id'])

        # Summed codon/amino acid results across each function source:
        # If amino acid rather than codon columns are returned, then gene amino acid frequencies are
        # first calculated from the sum of codon frequencies.
        ### Absolute frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(method).droplevel(1),
            sequence_min_amino_acids=sequence_min_amino_acids,
            output_amino_acids=return_amino_acids,
            label_amino_acids=label_amino_acids)
        if (not relative and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_frequency_table)
        ### Relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(
                partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
            output_amino_acids=return_amino_acids,
            label_amino_acids=label_amino_acids)
        if (relative and
            not synonymous and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(
                partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
            label_amino_acids=label_amino_acids,
            replace_na=infinity_to_zero)
        if (not return_amino_acids and
            relative and
            synonymous and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

        # Codon results averaged across genes annotated by a function source:
        ### Absolute frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(method).droplevel(1),
            sequence_min_amino_acids=sequence_min_amino_acids,
            label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            not relative and
            function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(
                partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
            label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            relative and
            not synonymous and
            function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(
                partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
            label_amino_acids=label_amino_acids,
            replace_na=infinity_to_zero)
        if (not return_amino_acids and
            relative and
            synonymous and
            function_sources and
            average_genes):
            return get_table(self._get_average_synonymous_codon_rel_frequency_table)

        # Amino acid results averaged across genes annotated by a function source:
        ### Absolute frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(method).droplevel(1),
            sequence_min_amino_acids=sequence_min_amino_acids,
            output_amino_acids=return_amino_acids)
        if (return_amino_acids and
            not relative and
            function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('function_source').apply(
                partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
            output_amino_acids=return_amino_acids)
        if (return_amino_acids and
            relative and
            function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)

        raise ConfigError("This point should not be reached at the end of the method, "
                          "`get_frequencies`. Please contact the developers. Since you found the "
                          "end of the earth, you now get to hear a top secret mnemonic for the "
                          "rare earth elements. Scandalous Yiddish language centers praise Ned's "
                          "promise of small European garden tubs. Dinosaurs hobble erotically "
                          "thrumming yellow lutes. (scandium Sc, yttrium Y, lanthanum La, cerium "
                          "Ce, praseodymium Pr, neodymium Nd, promethium Pm, samarium Sm, europium "
                          "Eu, gadolinium Gd, terbium Tb, dysprosium Dy, holmium Ho, erbium Er, "
                          "thulium Tm, ytterbium Yb, lutetium Lu) Credit for the lanthanide series "
                          "mnemonic goes to Martyn Poliakoff: "
                          "https://www.youtube.com/watch?v=Q21clW0s0B8&ab_channel=PeriodicVideos")


    ##################################################
    # `self.get_frequencies` SETUP HELPER METHODS
    def _establish_function_sources(self, from_function_sources):
        """Gets `function_sources` for methods that take the argument, `from_function_sources`."""
        if from_function_sources == True:
            function_sources = list(self.function_sources)
            if not function_sources:
                raise ConfigError(
                    "All function sources were requested, but none exist for the object.")
        elif from_function_sources != False:
            function_sources = list(from_function_sources)
        else:
            function_sources = None

        if function_sources is not None:
            unrecognized_function_sources = []
            for function_source in function_sources:
                if function_source not in self.function_sources:
                    unrecognized_function_sources.append(function_source)
            if unrecognized_function_sources:
                raise ConfigError("The requested function annotation sources, "
                                  f"{', '.join(function_sources)}, are not among those available: "
                                  f"{', '.join(self.function_sources)}.")

        return function_sources


    def _select_genes(self,
                      gene_caller_ids,
                      function_accession_dict,
                      function_name_dict,
                      expect_functions):
        """Select genes in the gene frequency table given a list of IDs and/or functions of
        interest."""
        if not gene_caller_ids and not function_accession_dict and not function_name_dict:
            # Nothing of interest was given.
            return self.gene_codon_frequency_df

        select_gene_caller_ids = []

        unrecognized_gene_caller_ids = []
        for gene_caller_id in gene_caller_ids:
            if gene_caller_id not in self.gene_caller_ids:
                unrecognized_gene_caller_ids.append(gene_caller_id)
                continue
            select_gene_caller_ids.append(gene_caller_id)
        if unrecognized_gene_caller_ids:
            raise ConfigError(
                "The following gene caller IDs were not found in the genome: "
                f"{unrecognized_gene_caller_ids}")

        ##################################################
        # Select by function accession.
        requested_keys = []
        unrecognized_sources = []
        for function_source, function_accessions in function_accession_dict.items():
            if function_source == 'KEGG_BRITE':
                self.run.warning(
                    "Nota bene: KEGG BRITE accessions stored in anvi'o are for hierarchies as a "
                    "whole, not categories of the hierarchy. Most hierarchies do not have category "
                    "accessions. So all genes in hierarchies with the given accessions are "
                    "selected.")
            if function_source not in self.function_sources:
                unrecognized_sources.append(function_source)
            for function_accession in function_accessions:
                requested_keys.append((function_source, function_accession))
        if unrecognized_sources:
            raise ConfigError(
                "The following annotation sources in `function_accessions` were not found to have "
                f"annotated any genes in the genome: {', '.join(unrecognized_sources)}")
        requested_keys = set(requested_keys)

        gene_function_df = self.gene_function_df.set_index(
            ['function_source', 'function_accession'])
        if expect_functions:
            # All requested accessions must be found.
            available_keys = set(gene_function_df.index)
            missing_keys = requested_keys.difference(available_keys)
            select_keys = requested_keys

            if missing_keys:
                # Prepare and throw an error message.
                missing_key_dict = {}
                for missing_key in missing_keys:
                    try:
                        missing_key_dict[missing_key[0]].append(missing_key[1])
                    except KeyError:
                        missing_key_dict[missing_key[0]] = [missing_key[1]]
                missing_key_message = ''
                for source, missing_accessions in missing_key_dict.items():
                    missing_key_message += source + ': ' + ', '.join(missing_accessions) + '; '
                missing_key_message = missing_key_message[: -2]
                if missing_keys:
                    raise ConfigError(
                        "The following requested function accessions are missing from the genome: "
                        f"{missing_key_message}")
        else:
            # Not all requested accessions need be found.
            select_keys = gene_function_df.index.intersection(requested_keys)
        select_gene_caller_ids += gene_function_df.loc[select_keys]['gene_caller_id'].tolist()
        ##################################################

        ##################################################
        # Select by function name.
        requested_keys = []
        unrecognized_sources = []
        for function_source, function_names in function_name_dict.items():
            if function_source not in self.function_sources:
                unrecognized_sources.append(function_source)
            for function_name in function_names:
                requested_keys.append((function_source, function_name))
        if unrecognized_sources:
            raise ConfigError(
                "The following annotation sources in `function_names` were not found to have "
                f"annotated any genes in the genome: {', '.join(unrecognized_sources)}")
        requested_keys = set(requested_keys)

        gene_function_df = gene_function_df.reset_index()

        def select_keys(gene_function_df, requested_keys):
            gene_function_df = gene_function_df.set_index(['function_source', 'function_name'])
            if expect_functions:
                available_keys = set(gene_function_df.index)
                missing_keys = requested_keys.difference(available_keys)
                select_keys = requested_keys
            else:
                select_keys = gene_function_df.index.intersection(requested_keys)
                missing_keys = None

            return select_keys, missing_keys

        def select_keys_given_higher_brite_categories(gene_function_df, requested_keys):
            # If `self.all_brite_categories` is False, then `gene_function_df` only contains rows
            # for the most specific KEGG BRITE categories. For example, a gene can be annotated with
            # "KEGG Orthology (KO)>>>09100 Metabolism>>>09101 Carbohydrate metabolism>>>00010
            # Glycolysis / Gluconeogenesis [PATH:ko00010]", but not categories representing higher
            # levels of the hierarchy: "KEGG Orthology (KO)>>>09100 Metabolism>>>09101 Carbohydrate
            # metabolism", "KEGG Orthology (KO)>>>09100 Metabolism", or "KEGG Orthology (KO)".
            # However, `function_name_dict` can be used to select genes by higher categories. "KEGG
            # Orthology (KO)>>>09100 Metabolism", for instance, encompasses "00010 Glycolysis /
            # Gluconeogenesis [PATH:ko00010]". BRITE functions in the table are parsed to determine
            # if they are encompassed by the requested function names.
            gene_brite_df = gene_function_df[gene_function_df['function_source'] == 'KEGG_BRITE']
            gene_nonbrite_df = gene_function_df.drop(gene_brite_df.index)
            brite_requested_keys = []
            nonbrite_requested_keys = []
            for requested_key in requested_keys:
                if requested_key[0] == 'KEGG_BRITE':
                    brite_requested_keys.append(requested_key[1])
                else:
                    nonbrite_requested_keys.append(requested_key)

            gene_brite_df = gene_brite_df.set_index('function_name')
            brite_select_keys = []
            brite_found_keys = []
            for available_key in gene_brite_df.index:
                for requested_key in requested_keys:
                    try:
                        position = available_key.index(requested_key)
                    except ValueError:
                        continue
                    if position == 0 and available_key[
                        len(requested_key): len(requested_key) + 3] == '>>>':
                        brite_select_keys.append(available_key)
                        brite_found_keys.append(requested_key)
            brite_select_keys = set(
                [('KEGG_BRITE', function_name) for function_name in set(brite_select_keys)])
            brite_found_keys = set(brite_found_keys)
            brite_missing_keys = brite_found_keys.difference(brite_select_keys)
            brite_missing_keys = set(
                [('KEGG_BRITE', function_name) for function_name in set(brite_missing_keys)])

            if len(gene_nonbrite_df) > 0:
                nonbrite_select_keys, nonbrite_missing_keys = select_keys(
                    gene_nonbrite_df, nonbrite_requested_keys)

                select_keys = brite_select_keys.union(nonbrite_select_keys)
                missing_keys = brite_missing_keys.union(nonbrite_missing_keys)

            return select_keys, missing_keys

        if self.all_brite_categories:
            select_keys, missing_keys = select_keys_given_higher_brite_categories(
                gene_function_df, requested_keys)
        else:
            select_keys, missing_keys = select_keys(gene_function_df, requested_keys)

        if expect_functions and missing_keys:
            # Prepare and throw an error message.
            missing_key_dict = {}
            for missing_key in missing_keys:
                try:
                    missing_key_dict[missing_key[0]].append(missing_key[1])
                except KeyError:
                    missing_key_dict[missing_key[0]] = [missing_key[1]]
            missing_key_message = ''
            for source, missing_names in missing_key_dict.items():
                missing_key_message += source + ': ' + ', '.join(missing_names) + '; '
            missing_key_message = missing_key_message[: -2]
            if missing_keys:
                raise ConfigError(
                    "The following requested function names are missing from the genome: "
                    f"{missing_key_message}")

        gene_function_df = gene_function_df.set_index(['function_source', 'function_name'])
        select_gene_caller_ids += gene_function_df.loc[select_keys]['gene_caller_id'].tolist()
        ##################################################

        return self.gene_codon_frequency_df.loc[set(select_gene_caller_ids)]


    ##################################################
    # `self.get_frequencies` OUTPUT HELPER METHODS
    def _filter_input_codon_count_decorator(method):
        """Decorator to discard rows in the input frequency table with fewer than the minimum number
        of codons/amino acids."""
        def wrapper(*args, **kwargs):
            frequency_df = args[1]
            try:
                filter_input_codon_count = kwargs['filter_input_codon_count']
                min_codons = kwargs['min_codons']
            except KeyError:
                filter_input_codon_count = False
                min_codons = 0
            if filter_input_codon_count and min_codons:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= min_codons]
            return method(args[0], frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _drop_amino_acid_codon_columns_decorator(method):
        """Decorator to discard codon columns by amino acid from the input frequency table."""
        def wrapper(*args, **kwargs):
            codon_frequency_df = args[1]
            try:
                drop_amino_acids = kwargs['drop_amino_acids']
            except KeyError:
                drop_amino_acids = []
            if drop_amino_acids:
                drop_codons = []
                amino_acid_codons_dict = args[0].amino_acid_codons_dict
                for amino_acid in drop_amino_acids:
                    drop_codons += amino_acid_codons_dict[amino_acid]
                codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1, errors='ignore')
                if len(codon_frequency_df.columns) == 0:
                    codon_frequency_df = codon_frequency_df.drop(codon_frequency_df.index)
            return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _filter_pansequence_synonymous_codon_count_decorator(method):
        """Decorator to drop synonymous codons encoding an amino acid from the input frequency
        table based on the frequency of the amino acid across rows."""
        def wrapper(*args, **kwargs):
            codon_frequency_df = args[1]
            try:
                pansequence_min_amino_acids = kwargs['pansequence_min_amino_acids']
            except KeyError:
                pansequence_min_amino_acids = (0, 1.0)
            if pansequence_min_amino_acids[0] > 0 and 0 < pansequence_min_amino_acids[1] < 1:
                row_count = len(codon_frequency_df)
                drop_codons = []
                for synonymous_codons in args[0].amino_acid_codons_dict.values():
                    try:
                        filtered_row_count = len(codon_frequency_df[codon_frequency_df[
                            synonymous_codons].sum(axis=1) >= pansequence_min_amino_acids[0]])
                    except KeyError:
                        # This occurs when codons are missing from the frequency table.
                        continue
                    if filtered_row_count / row_count < pansequence_min_amino_acids[1]:
                        drop_codons += synonymous_codons
                codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1)
                if len(codon_frequency_df.columns) == 0:
                    codon_frequency_df = codon_frequency_df.drop(codon_frequency_df.index)
            return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _filter_sequence_synonymous_codon_count(self, codon_frequency_df, sequence_min_amino_acids):
        mask_df = pd.DataFrame()
        for synonymous_codons in self.amino_acid_codons_dict.values():
            try:
                codon_mask_series = \
                    codon_frequency_df[synonymous_codons].sum(axis=1) >= sequence_min_amino_acids
            except KeyError:
                # This occurs when codons are missing from the frequency table.
                continue
            amino_acid_mask_df = pd.DataFrame(index=codon_mask_series.index)
            for codon in synonymous_codons:
                amino_acid_mask_df[codon] = codon_mask_series
            mask_df = pd.concat([mask_df, amino_acid_mask_df], axis=1)
        codon_frequency_df = codon_frequency_df[mask_df]
        return codon_frequency_df


    def _filter_sequence_synonymous_codon_count_decorator(method):
        """Decorator to replace data with NaN for synonymous codons encoding an amino acid from the
        input frequency table based on the summed frequency of the synonymous codons in each row."""
        def wrapper(*args, **kwargs):
            codon_frequency_df = args[1]
            try:
                sequence_min_amino_acids = kwargs['sequence_min_amino_acids']
            except KeyError:
                sequence_min_amino_acids = 0
            if sequence_min_amino_acids > 0:
                mask_df = pd.DataFrame()
                for synonymous_codons in args[0].amino_acid_codons_dict.values():
                    try:
                        codon_mask_series = (codon_frequency_df[synonymous_codons].sum(axis=1) >=
                                             sequence_min_amino_acids)
                    except KeyError:
                        # This occurs when codons are missing from the frequency table.
                        continue
                    amino_acid_mask_df = pd.DataFrame(index=codon_mask_series.index)
                    for codon in synonymous_codons:
                        amino_acid_mask_df[codon] = codon_mask_series
                    mask_df = pd.concat([mask_df, amino_acid_mask_df], axis=1)
                codon_frequency_df = codon_frequency_df[mask_df]
            return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _output_amino_acids_decorator(method):
        """Decorator to output columns of amino acid rather than codon frequencies."""
        def wrapper(*args, **kwargs):
            codon_df = method(*args, **kwargs)
            try:
                output_amino_acids = kwargs['output_amino_acids']
            except KeyError:
                output_amino_acids = False
            if output_amino_acids:
                aa_df = pd.DataFrame(index=codon_df.index)
                for amino_acid, codons in args[0].amino_acid_codons_dict.items():
                    try:
                        aa_df[amino_acid] = codon_df[codons].sum(axis=1, skipna=False)
                    except KeyError:
                        # This occurs when there aren't columns for codons.
                        pass
                return aa_df
            return codon_df
        return wrapper


    def _add_amino_acid_to_header_decorator(method):
        """Decorator to add amino acid to codon column header."""
        def wrapper(*args, **kwargs):
            codon_df = method(*args, **kwargs)
            try:
                label_amino_acids = kwargs['label_amino_acids']
            except KeyError:
                label_amino_acids = False
            if label_amino_acids:
                try:
                    codon_df.columns = [
                        constants.codon_to_AA[codon] + codon for codon in codon_df.columns]
                except KeyError:
                    raise ConfigError("The columns in what should be a table of codon data are not "
                                      "recognized as codons. This is the header that is present: "
                                      f"{', '.join(list(codon_df.columns))}")
                codon_df = codon_df[sorted(codon_df.columns)]
            return codon_df
        return wrapper


    def _replace_na(method):
        """Decorator to replace NA with 0.0 in the output frequency table. This should only occur in
        synonymous relative frequency output, in which missing amino acids yield NA (not inf, though
        this is as a consequence of 0/0)."""
        def wrapper(*args, **kwargs):
            frequency_df = method(*args, **kwargs)
            try:
                replace_na = kwargs['replace_na']
            except KeyError:
                replace_na = False
            if replace_na:
                frequency_df = frequency_df.fillna(0)
            return frequency_df
        return wrapper


    def _filter_output_codon_count(method):
        """Decorator to discard rows in the output frequency table with fewer than the minimum
        number of codons/amino acids."""
        def wrapper(*args, **kwargs):
            frequency_df = method(*args, **kwargs)
            try:
                filter_output_codon_count = kwargs['filter_output_codon_count']
                min_codons = kwargs['min_codons']
            except KeyError:
                filter_output_codon_count = False
            if filter_output_codon_count:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= min_codons]
            return frequency_df
        return wrapper
    ##################################################


    # The order of decorators should not be changed (only @_output_amino_acids_decorator and
    # @_add_amino_acid_to_header_decorator, which are mutually exclusive operations, are
    # interchangeable).
    @_filter_input_codon_count_decorator
    @_drop_amino_acid_codon_columns_decorator
    @_filter_pansequence_synonymous_codon_count_decorator
    @_filter_sequence_synonymous_codon_count_decorator
    @_output_amino_acids_decorator
    @_add_amino_acid_to_header_decorator
    @_filter_output_codon_count
    def _get_frequency_table(self, codon_frequency_df, **kwargs):
        return codon_frequency_df


    # Commented decorators mean that they can theoretically be used and uncommented but are not
    # because they are not needed in the `self.get_frequencies` client.
    # @_filter_input_codon_count_decorator
    # @_drop_amino_acid_codon_columns_decorator
    # @_filter_pansequence_synonymous_codon_count_decorator
    @_output_amino_acids_decorator
    @_add_amino_acid_to_header_decorator
    # @_filter_output_codon_count
    def _get_rel_frequency_table(self, codon_frequency_df, **kwargs):
        try:
            mask_df = self._filter_sequence_synonymous_codon_count(
                codon_frequency_df, kwargs['sequence_min_amino_acids']).notna()
        except KeyError:
            mask_df = None
        codon_rel_frequency_df = codon_frequency_df.div(codon_frequency_df.sum(axis=1), axis=0)
        drop_index = codon_rel_frequency_df[codon_rel_frequency_df.isna().all(axis=1)].index
        if mask_df is not None:
            codon_rel_frequency_df = codon_rel_frequency_df[mask_df]
        # Drop rows with zero frequency.
        codon_rel_frequency_df = codon_rel_frequency_df.drop(drop_index)
        return codon_rel_frequency_df


    # @_filter_input_codon_count_decorator
    # @_drop_amino_acid_codon_columns_decorator
    # @_filter_pansequence_synonymous_codon_count_decorator
    @_filter_sequence_synonymous_codon_count_decorator
    @_add_amino_acid_to_header_decorator
    @_replace_na
    # @_filter_output_codon_count
    def _get_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        """Return the relative frequencies of codons in relation to the set of codons encoding the
        same amino acid (or stop codons). If columns for one or more codons in a synonymous set are
        missing (rather than 0), synonymous relative frequency will not be calculated for the
        remaining codons in the set."""
        try:
            mask_df = self._filter_sequence_synonymous_codon_count(
                codon_frequency_df, kwargs['sequence_min_amino_acids']).notna()
        except KeyError:
            mask_df = None
        synonymous_codon_rel_frequency_df = pd.DataFrame()
        for codons in constants.AA_to_codons.values():
            try:
                aa_codon_frequency_df = codon_frequency_df[codons]
            except KeyError:
                # This occurs when synonymous codons are missing from the frequency table.
                continue
            synonymous_codon_rel_frequency_df[codons] = aa_codon_frequency_df.div(
                aa_codon_frequency_df.sum(axis=1), axis=0)
        drop_index = synonymous_codon_rel_frequency_df[
            synonymous_codon_rel_frequency_df.isna().all(axis=1)].index
        if mask_df is not None:
            synonymous_codon_rel_frequency_df = synonymous_codon_rel_frequency_df[mask_df]
        synonymous_codon_rel_frequency_df = synonymous_codon_rel_frequency_df.drop(drop_index)
        return synonymous_codon_rel_frequency_df


    # @_filter_input_codon_count_decorator
    # @_drop_amino_acid_codon_columns_decorator
    # @_filter_pansequence_synonymous_codon_count_decorator
    # @_filter_sequence_synonymous_codon_count_decorator
    def _get_summed_frequency_table(self, frequency_df, **kwargs):
        """Return the summed frequencies across all items."""
        summed_frequency_series = frequency_df.sum()
        summed_frequency_df = \
            summed_frequency_series.to_frame('all').T.rename_axis('gene_caller_id')
        return summed_frequency_df


    def _get_summed_rel_frequency_table(self, frequency_df, **kwargs):
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in ['sequence_min_amino_acids', 'label_amino_acids']:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value

        summed_frequency_df = self._get_summed_frequency_table(frequency_df, **first_kwargs)
        summed_rel_frequency_df = self._get_rel_frequency_table(
            summed_frequency_df, **second_kwargs)
        return summed_rel_frequency_df


    def _get_summed_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in ['sequence_min_amino_acids', 'label_amino_acids', 'replace_na']:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value

        summed_codon_frequency_df = self._get_summed_frequency_table(
            codon_frequency_df, **first_kwargs)
        summed_synonymous_codon_rel_frequency_df = self._get_synonymous_codon_rel_frequency_table(
            summed_codon_frequency_df, **second_kwargs)
        return summed_synonymous_codon_rel_frequency_df


    # @_filter_input_codon_count_decorator
    # @_drop_amino_acid_codon_columns_decorator
    # @_filter_synonymous_codon_count
    def _get_average_frequency_table(self, frequency_df, **kwargs):
        """Return the average codon frequencies across all items."""
        average_codon_frequency_series = frequency_df.mean()
        average_codon_frequency_df = \
            average_codon_frequency_series.to_frame('all').T.rename_axis('gene_caller_id')
        return average_codon_frequency_df


    def _get_average_rel_frequency_table(self, frequency_df, **kwargs):
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in ['sequence_min_amino_acids', 'label_amino_acids']:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value

        average_frequency_df = self._get_average_frequency_table(frequency_df, **first_kwargs)
        average_rel_frequency_df = \
            self._get_rel_frequency_table(average_frequency_df, **second_kwargs)
        return average_rel_frequency_df


    def _get_average_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in ['sequence_min_amino_acids', 'label_amino_acids', 'replace_na']:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value

        average_codon_frequency_df = self._get_average_frequency_table(
            codon_frequency_df, **first_kwargs)
        average_synonymous_codon_rel_frequency_df = self._get_synonymous_codon_rel_frequency_table(
            average_codon_frequency_df, **second_kwargs)
        return average_synonymous_codon_rel_frequency_df


    def get_codon_usage_bias(
        self,
        metrics=None,
        from_function_sources=None,
        gene_caller_ids=None,
        function_accessions=None,
        function_names=None,
        expect_functions=False,
        omnibias=False,
        reference_function_accessions=None,
        reference_function_names=None,
        expect_reference_functions=False,
        reference_gene_caller_ids=None,
        gene_min_codons=0,
        function_min_codons=0,
        min_codon_filter='both',
        drop_amino_acids=None,
        sequence_min_amino_acids=0,
        pansequence_min_amino_acids=(0, 1.0),
        query_min_analyzed_codons=default_query_min_analyzed_codons,
        reference_exclude_amino_acid_count=default_reference_exclude_amino_acid_count,
        reference_min_analyzed_codons=default_reference_min_analyzed_codons):
        """
        Get codon usage bias (CUB) of genes or functions.

        Parameters
        ==========
        metrics : {'cai', 'delta'}, optional
            CUB metric, with valid choices being 'cai' (Codon Adaptation Index of Sharp and Li,
            1987) and 'delta' (Ran and Higgs, 2012, Eq. 6). The default of None calculates all
            available CUB metrics. For CUB metrics that require reference genes, the default
            behavior in the absence of supplied reference genes/functions or omnibias mode is to use
            the summed composition of ribosomal proteins defined by KEGG KOfams/BRITE.
        from_function_sources : bool, str, or iterable of str, optional
            Select genes with functional annotations, using their summed codon frequencies to
            calculate function CUB values. With this argument, the first three columns of the
            returned table contain, respectively, function annotation sources, function accessions,
            and function names. When this argument is True, use all available functional annotation
            sources in the SingleGenomeCodonUsage object. When this argument is a string, select the
            source given by the string, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. When this argument
            is an iterable, select a subset of sources given by its strings, e.g., ['KOfam',
            'COG20_FUNCTION']. By default None.
        gene_caller_ids : iterable, optional
            Genes with the given IDs are selected for calculation of CUB. This parameter can be
            used alongside `function_accessions` and `function_names`. By default None.
        function_accessions : dict, optional
            Genes annotated with the given function accessions are selected for calculation of CUB.
            The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
            parameter can be used alongside `gene_caller_ids` and `function_names`.  Note that
            'KEGG_BRITE' does not use individual function accessions but overarching hierarchy
            accessions that include multiple functions. By default None.
        function_names : dict, optional
            Genes annotated with the given function names are selected for calculation of CUB. The
            argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
            `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
            (categorizations) given to an arbitrary level of the hierarchy, such as
            ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
            type>>>Bacterial type>>>RNA polymerase']. All BRITE categorizations falling under the
            given level are selected. This parameter can be used alongside `gene_caller_ids` and
            `function_accessions`. By default None.
        expect_functions : bool, optional
            If True (default False), an error will be raised if any given `function_accessions` or
            `function_names` are not annotated in the input genome.
        omnibias : bool, optional
            If True (default False), use every gene or function as a separate reference rather than
            defining a set of reference genes or functions. The resulting table of gene x gene (or
            function x function) CUB values is like a distance matrix of the similarity of gene
            codon compositions.
        reference_function_accessions : dict, optional
            Genes annotated with the given function accessions are selected for the reference set.
            The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
            parameter can be used alongside `reference_gene_caller_ids` and
            `reference_function_names`. Note that 'KEGG_BRITE' does not use individual function
            accessions but overarching hierarchy accessions that include multiple functions. By
            default None.
        reference_function_names : dict, optional
            Genes annotated with the given function names are selected for the reference set. The
            argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
            `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
            (categorizations) given to an arbitrary level of the hierarchy, such as
            ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
            type>>>Bacterial type>>>RNA polymerase']. All BRITE categorizations falling under the
            given level are selected. This parameter can be used alongside
            `reference_gene_caller_ids` and `reference_function_accessions`. By default, if using
            reference-dependent metrics without omnibias mode, this argument becomes {'KEGG_BRITE':
            ['Ribosome>>>Ribosomal proteins']}.
        expect_reference_functions : bool, optional
            If True (default False), an error will be raised if any given
            `reference_function_accessions` or `reference_function_names` are not annotated in the
            input genome.
        reference_gene_caller_ids : iterable, optional
            Include specific genes in the reference gene set, as given by their gene caller IDs in
            the contigs database.

        Additional Parameters
        =====================
        These parameters filter genes and codons in CUB calculations.

        The gene codon frequency table is the input to the CUB calculation. Here is the order of
        all possible filters for \"queries\":
        gene codon frequency table ->
            drop genes on codon length ->
            drop gene/function pairs on the total length of remaining genes defining the function ->
            drop codons of defined amino acids ->
            dynamically drop codons of rarer amino acids ->
            drop genes on remaining codon frequency ->
            drop gene/function pairs on the codon sum of remaining genes defining the function ->
        filtered gene codon frequency table

        There is a filter, `query_min_analyzed_codons`, that excludes queries based on the number of
        codons participating in CUB analysis.

        There are two filters, 'reference_exclude_amino_acid_count' and
        'reference_min_analyzed_codons', that are applied to the reference gene set.

        gene_min_codons : int, optional
            Ignore genes with fewer than this number of codons in gene/function queries. By default
            0.
        function_min_codons : int, optional
            Ignore gene/function pairs when the function has fewer than this number of codons. Genes
            are filtered by codon count before functions, so the total length of remaining genes
            annotated by the function is considered. By default 0.
        min_codon_filter : {"length", "remaining", "both"}, optional
            This argument arises from the ambiguity of filters that remove genes and functions by
            number of codons (`--gene-min-codons` and `--function-min-codons`) in relation to the
            filters that drop codons (`--exclude/include-amino-acids` and
            `--exclude-amino-acid-count/fraction`). Genes (and functions) can be filtered by their
            full length, e.g., genes shorter than 300 codons are ignored. They can also be
            filtered by the number of codons remaining after dropping codons. The codon length
            filter followed by dropping codons can result in genes and functions with fewer codons
            than the original codon threshold -- thus the option of both 'length' and 'remaining'
            filters to ensure that total codon frequencies always meet the minimum codon threshold.
            'both' is needed as an option in addition to 'remaining' so dynamic codon filtering by
            `--exclude-amino-acid-count/fraction` operates on genes that passed the first length
            filter. By default "both".
        drop_amino_acids : iterable, optional
            Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala).
            Met and Trp are encoded by single codons, which perforce are excluded from CUB
            calculations. By default, stop codons are also excluded from the calculation: the
            default, None, becomes ['STP'] in the code. Importantly, to continue to exclude stop
            codons, make sure to include it in the passed value: exclusion of Ala and stop codons
            is achieved by passing ['Ala', 'STP'].
        sequence_min_amino_acids : int, optional
            Remove codons for amino acids (and STP) that are less numerous than
            `sequence_min_amino_acids`. For example, if the argument is 5, and a gene query has 4
            codons encoding Asn, 2 AAT and 2 AAC, then Asn codons will be disregarded in the
            calculation of CUB for this query. By default 0.
        pansequence_min_amino_acids : tuple, optional
            This tuple must have two items. The first is an int >0 representing a minimum number of
            codons encoding an amino acid -- 'min_amino_acids' -- and the second is a float in the
            range (0, 1) representing a fraction of genes -- 'min_gene_fraction'. Remove codons for
            amino acids (and STP) that are less numerous than 'min_amino_acids' in a
            'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and
            'min_gene_fraction' is 0.9, then if there are fewer than 5 codons for an amino acid/STP
            in ≥90% of genes, the amino acid's codons do not factor into the calculation of CUB for
            any query. By default (0, 1.0).
        query_min_analyzed_codons : int, optional
            Only allow CUB to calculated for a query if has at least this number of synonymous
            codons that will be analyzed. For reference-dependent CUB metrics, analyzed codons are
            those with reference compositions. By default 100.
        reference_exclude_amino_acid_count : int, optional
            Exclude codons for amino acids with fewer than this many codons in the set of reference
            genes. This does not apply in `omnibias` mode. By default 5.
        reference_min_analyzed_codons : int, optional
            Only allow CUB to be calculated using a reference if it has at least this number of
            codons that will be analyzed. This filter applies after excluding codons for individual
            amino acids using `reference_exclude_amino_acid_count`. By default 100.

        Returns
        =======
        dict of pandas.core.frame.DataFrame objects
            Each dict key is the name of a CUB metric (from `metrics`), and each value is a
            corresponding table of CUB data.
        """
        metrics = self._establish_cub_metrics(metrics)
        reference_metrics = []
        referenceless_metrics = []
        for metric in metrics:
            if metric in reference_dependent_cub_metrics:
                reference_metrics.append(metric)
            else:
                referenceless_metrics.append(metric)

        # Check that proper arguments were provided given the metrics.
        if reference_metrics:
            if (omnibias and
                (reference_function_accessions or
                 reference_function_names or
                 expect_reference_functions or
                 reference_gene_caller_ids)):
                raise ConfigError(
                    "Omnibias mode cannot be used when defined gene/function references are also "
                    "used. The following arguments are only relevant to defined references: "
                    "`reference_function_accessions`, `reference_function_names`, "
                    "`expect_reference_functions`, `reference_gene_caller_ids`, "
                    "`reference_exclude_amino_acid_count`, and `reference_min_analyzed_codons`. "
                    "Query filters including `sequence_min_amino_acids`, "
                    "`pansequence_min_amino_acids`, `min_gene_fraction`, and "
                    "`query_min_analyzed_codons` apply to references since queries are the same as "
                    "references in omnibias mode.")
            if omnibias and (query_min_analyzed_codons != reference_min_analyzed_codons):
                raise ConfigError(
                    "In omnibias mode, `query_min_analyzed_codons` and "
                    "`reference_min_analyzed_codons` should have the same value, since the sets of "
                    "query and reference genes/functions should be the same, and so should be "
                    "filtered in the same way.")
        elif (omnibias or
              reference_function_accessions or
              reference_function_names or
              expect_reference_functions or
              reference_gene_caller_ids):
            raise ConfigError(
                "The provided CUB metrics do not involve comparison of gene/function codon "
                "compositions. The following arguments are only relevant to "
                "reference-dependent metrics: `omnibias`, `reference_function_accessions`, "
                "`reference_function_names`, `expect_reference_functions`, "
                "`reference_gene_caller_ids`, `reference_exclude_amino_acid_count`, and "
                "`reference_min_analyzed_codons`.")

        # Get a reference codon composition when using reference-dependent metrics and not in
        # omnibias mode.
        if reference_metrics and not omnibias:
            if not reference_function_accessions and not reference_function_names:
                reference_function_accessions = {
                    default_reference_function_source: default_reference_function_accessions}
                reference_function_names = {
                    default_reference_function_source: default_reference_function_names}
            reference_codon_frequency_df = self._get_defined_reference_codon_frequencies(
                reference_metrics,
                reference_function_accessions=reference_function_accessions,
                reference_function_names=reference_function_names,
                expect_reference_functions=expect_reference_functions,
                reference_gene_caller_ids=reference_gene_caller_ids,
                reference_exclude_amino_acid_count=reference_exclude_amino_acid_count,
                reference_min_analyzed_codons=reference_min_analyzed_codons)
        else:
            reference_codon_frequency_df = None

        # Determine the encoded amino acids ignored in the analysis.
        if drop_amino_acids is None:
            drop_amino_acids = ignored_cub_amino_acids
        else:
            for amino_acid in single_codon_amino_acids:
                if amino_acid not in drop_amino_acids:
                    drop_amino_acids.append(amino_acid)

        self.run.info_single("Queries", mc='green')
        query_codon_frequency_df = self.get_frequencies(
            from_function_sources=from_function_sources,
            return_functions=bool(from_function_sources),
            gene_caller_ids=gene_caller_ids,
            function_accessions=function_accessions,
            function_names=function_names,
            expect_functions=expect_functions,
            gene_min_codons=gene_min_codons,
            function_min_codons=function_min_codons,
            min_codon_filter=min_codon_filter,
            drop_amino_acids=drop_amino_acids,
            sequence_min_amino_acids=sequence_min_amino_acids,
            pansequence_min_amino_acids=pansequence_min_amino_acids)
        if (gene_caller_ids or
            function_accessions or
            function_names or
            gene_min_codons or
            function_min_codons):
            self.run.info_single(
                f"{pp(len(query_codon_frequency_df))} "
                f"{'function' if from_function_sources else 'CDS'} queries retained")

        # Report codons missing in the query dataset.
        missing_query_codons = []
        for codon in query_codon_frequency_df.columns:
            if query_codon_frequency_df[codon].sum() == 0:
                missing_query_codons.append(codon)
        if missing_query_codons:
            missing_codon_dict = {}
            for codon in missing_query_codons:
                amino_acid = self.codon_amino_acid_dict[codon]
                try:
                    missing_codon_dict[amino_acid].append(codon)
                except KeyError:
                    missing_codon_dict[amino_acid] = [codon]
            missing_codon_message = ""
            for amino_acid, codons in sorted(missing_codon_dict.items()):
                missing_codon_message += amino_acid + ": " + ", ".join(codons) + "; "
            missing_codon_message = missing_codon_message[: -2]
            self.run.warning(
                "The following codons from synonymous decoding sets are absent in the retained "
                f"queries. {missing_codon_message}")

        if omnibias:
            reference_codon_frequency_df = query_codon_frequency_df
        else:
            # Here a defined reference codon composition is being used. Codons that are absent in
            # the query dataset obviously are not compared to the reference; check that the
            # reference has enough compared codons to meet `reference_min_analyzed_codons`.
            total_codons = reference_codon_frequency_df.drop(
                [codon for codon in missing_query_codons
                 if codon in reference_codon_frequency_df.columns], axis=1).sum(axis=1).sum()
            if total_codons < reference_min_analyzed_codons:
                # The reference does not have enough codons, so leave an empty shell of a DataFrame.
                try:
                    reference_codon_frequency_df = reference_codon_frequency_df.drop('all')
                    self.run.warning(
                        "A reference codon composition could not be established. "
                        "Reference-dependent CUB tables will be empty. The reference dataset has "
                        f"{pp(int(total_codons))} codons that will be used in the CUB analysis, "
                        "which does not meet the minimum defined by `reference_min_analyzed_codons`, "
                        f"{pp(int(reference_min_analyzed_codons))}. Note that codons are only analyzed "
                        "when present in both the query and reference datasets.")
                except KeyError:
                    # Even before removing the codons not shared by the query, the reference did not
                    # meet the minimum codon threshold.
                    pass

        if 'delta' in metrics:
            # Reference codon weights in the computation of 𝛿 involve comparison of the reference to
            # the non-reference codon set, taken from all other genes in the input dataset.
            amino_acids_dropped_from_reference = []
            for amino_acid, codons in constants.AA_to_codons.items():
                if codons[0] not in reference_codon_frequency_df.columns:
                    amino_acids_dropped_from_reference.append(amino_acid)
            default_run = self.run
            self.run = self.run_quiet
            nonreference_codon_frequency_df = self.get_frequencies(
                sum_genes=True, drop_amino_acids=amino_acids_dropped_from_reference)
            self.run = default_run
            if omnibias:
                # For perfect consistency with the "non-omnibias" mode using a defined reference
                # gene set, the codon frequencies of the reference gene/function in each omnibias
                # query-reference comparison would be subtracted from the non-reference codon
                # frequencies summed from all coding sequences in the genome, since the reference
                # composition is subtracted from the whole-genome composition to find the
                # non-reference composition in "non-omnibias" mode. For the sake of simplicity, this
                # detail is ignored.
                pass
            else:
                nonreference_codon_frequency_df = \
                    nonreference_codon_frequency_df - reference_codon_frequency_df.sum()
                # It is theoretically possible that there are fewer codons in the non-reference genes
                # than the reference genes. Ensure that the non-reference set meets the minimum codon
                # threshold of the reference set.
                total_codons = nonreference_codon_frequency_df.drop(
                    [codon for codon in missing_query_codons
                     if codon in nonreference_codon_frequency_df.columns], axis=1).sum(axis=1).sum()
                if total_codons < reference_min_analyzed_codons:
                    nonreference_codon_frequency_df = nonreference_codon_frequency_df.drop('all')
                    self.run.warning(
                        "The non-reference codon composition needed for 𝛿 could not be "
                        "established, because the non-reference genes have "
                        f"{pp(int(total_codons))} codons, fewer than the minimum of "
                        f"{pp(int(reference_min_analyzed_codons))} given by "
                        "`reference_min_analyzed_codons` that is required for the reference, and "
                        "by extension, non-reference datasets. Reference-dependent CUB tables will "
                        "be empty.")
        else:
            nonreference_codon_frequency_df = None

        cub_table_dict = {}
        for metric in reference_metrics:
            if metric == 'cai':
                cub_df = self._get_cai_table(
                    query_codon_frequency_df,
                    reference_codon_frequency_df,
                    query_min_analyzed_codons=query_min_analyzed_codons,
                    reference_min_analyzed_codons=reference_min_analyzed_codons)
                if len(cub_df.columns) > 0:
                    if omnibias:
                        cub_df.columns = cub_df.index
                    else:
                        cub_df.columns = ['CAI']
                    self.run.info_single(
                        f"CAI calculated for {pp(len(cub_df))} "
                        f"{'function' if from_function_sources else 'CDS'} queries")
            elif metric == 'delta':
                cub_df = self._get_delta_table(
                    query_codon_frequency_df,
                    reference_codon_frequency_df,
                    nonreference_codon_frequency_df,
                    query_min_analyzed_codons=query_min_analyzed_codons,
                    reference_min_analyzed_codons=reference_min_analyzed_codons)
                if len(cub_df.columns) > 0:
                    if omnibias:
                        cub_df.columns = cub_df.index
                    else:
                        cub_df.columns = ['Delta']
                    self.run.info_single(
                        f"𝛿 calculated for {pp(len(cub_df))} "
                        f"{'function' if from_function_sources else 'CDS'} queries")

            cub_table_dict[metric] = cub_df

            if len(cub_df.columns) == 0:
                continue

            # Print the number of query-reference comparisons that were thrown out.
            if omnibias:
                # This includes "self-comparisons," reported in the matrix diagonal.
                possible_comparison_count = int(len(cub_df) ** 2 / 2)
            else:
                possible_comparison_count = len(cub_df)
            unperformed_comparison_count = int(cub_df.isna().sum().sum())
            if unperformed_comparison_count:
                self.run.warning(
                    f"{pp(unperformed_comparison_count)} of {pp(possible_comparison_count)} "
                    "query-reference comparisons did not meet the minimum codon threshold in "
                    "either query or reference and so did not yield a CUB value.")

        # Calculate CUB using metrics that do not depend on a reference codon composition.
        for metric in referenceless_metrics:
            # No reference-independent CUB metrics have been programmed yet.
            pass

        return cub_table_dict


    def _establish_cub_metrics(self, metrics):
        """Establishes CUB metrics for `get_codon_usage_bias`."""
        if metrics is None:
            metrics = cub_metrics
        else:
            for metric in metrics:
                unrecognized_metrics = []
                if metric not in cub_metrics:
                    unrecognized_metrics.append(metric)
                if unrecognized_metrics:
                    raise ConfigError("The following CUB metrics are not recognized: "
                                      f"{', '.join(unrecognized_metrics)}. Here are the available "
                                      f"CUB metrics: {', '.join(cub_metrics)}")
        return metrics


    def _get_defined_reference_codon_frequencies(
        self,
        reference_metrics,
        reference_function_accessions=None,
        reference_function_names=None,
        expect_reference_functions=False,
        reference_gene_caller_ids=None,
        reference_exclude_amino_acid_count=default_reference_exclude_amino_acid_count,
        reference_min_analyzed_codons=default_reference_min_analyzed_codons):
        """
        Get a codon frequency table from a set of reference genes or functions in the genome.

        Parameters
        ==========
        reference_metrics : iterable
            CUB metrics requiring a reference codon composition.
        reference_function_accessions : dict, optional
            Genes annotated with the given function accessions are selected for the reference set.
            The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
            parameter can be used alongside `reference_gene_caller_ids` and
            `reference_function_names`. Note that 'KEGG_BRITE' does not use individual function
            accessions but overarching hierarchy accessions that include multiple functions. By
            default None.
        reference_function_names : dict, optional
            Genes annotated with the given function names are selected for the reference set. The
            argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
            `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
            (categorizations) given to an arbitrary level of the hierarchy, such as
            ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
            type>>>Bacterial type>>>RNA polymerase']. All BRITE categorizations falling under the
            given level are selected. This parameter can be used alongside
            `reference_gene_caller_ids` and `reference_function_accessions`. By default,
            {'KEGG_BRITE': ['Ribosome>>>Ribosomal proteins']}.
        expect_reference_functions : bool, optional
            If True (default False), an error will be raised if any given
            `reference_function_accessions` or `reference_function_names` are not annotated in the
            input genome.
        reference_gene_caller_ids : iterable, optional
            Include specific genes in the reference gene set, as given by their gene caller IDs in
            the contigs database.
        reference_exclude_amino_acid_count : int, optional
            Exclude codons for amino acids with fewer than this many codons in the set of reference
            genes. This does not apply in `omnibias` mode. By default 5.
        reference_min_analyzed_codons : int, optional
            Only allow CUB to be calculated using a reference if it has at least this number of
            codons that will be analyzed. This filter applies after excluding codons for individual
            amino acids using `reference_exclude_amino_acid_count`. By default 100.

        Returns
        =======
        pandas.core.frame.DataFrame
            This frequency table has a single row for the reference composition and a column per
            codon.
        """
        self.run.info_single("Reference codon composition", mc='green')

        # The default reference genes are ribosomal proteins, as annotated by KOfams/BRITE.
        if reference_function_names is None:
            if 'KEGG_BRITE' not in self.function_sources:
                raise ConfigError(
                    f"Reference-dependent metrics ({', '.join(reference_metrics)}) were requested "
                    "without defined reference genes. By default, reference genes are KEGG KOfams "
                    "classified as ribosomal proteins in BRITE. However, 'KEGG_BRITE' is not among "
                    "the function annotation sources run on the genome. This can be rectified by "
                    "rerunning `anvi-run-kegg-kofams`.")
            reference_function_names = {'KEGG_BRITE': ['Ribosome>>>Ribosomal proteins']}

        reference_codon_frequency_df = self.get_frequencies(
            gene_caller_ids=reference_gene_caller_ids,
            function_accessions=reference_function_accessions,
            function_names=reference_function_names,
            expect_functions=expect_reference_functions,
            sum_genes=True,
            drop_amino_acids=ignored_cub_amino_acids)

        # In the absence of a set of reference genes for the genome,
        # `reference_codon_frequency_df` has an index and header but no data.
        if len(reference_codon_frequency_df) == 0:
            self.run.warning(
                "A reference codon composition could not be established because none of the "
                "requested genes or functions were identified in the genome. "
                "Reference-dependent CUB tables will be empty.")
            return reference_codon_frequency_df

        # Report codons absent in the reference.
        missing_reference_codons = []
        for codon in set(constants.codons).difference(set(self.ignored_cub_codons)):
            if reference_codon_frequency_df[codon].sum() == 0:
                missing_reference_codons.append(codon)
        missing_codon_dict = {}
        for codon in missing_reference_codons:
            amino_acid = self.codon_amino_acid_dict[codon]
            try:
                missing_codon_dict[amino_acid].append(codon)
            except KeyError:
                missing_codon_dict[amino_acid] = [codon]
        missing_codon_message = ""
        for amino_acid, codons in sorted(missing_codon_dict.items()):
            missing_codon_message += amino_acid + ": " + ", ".join(codons) + "; "
        missing_codon_message = missing_codon_message[: -2]
        if missing_reference_codons:
            self.run.warning(
                "The following synonymous codons are absent in the reference. These codons will "
                "not contribute to CUB. This can skew CUB if the query contains many codons that "
                f"are absent in the reference. {missing_codon_message}")

        # Remove and report codons (columns) encoding amino acids that do not meet the minimum
        # reference amino acid count threshold, `reference_exclude_amino_acid_count`. Note that for
        # technical reasons requiring the preservation of columns for complete sets of synonymous
        # codons, columns representing complete sets of synonymous codons may be removed in this
        # method, but columns for individual missing codons are not removed.
        if reference_exclude_amino_acid_count:
            removed_amino_acids = []
            removed_codons = []
            for amino_acid, codons in self.synonymous_nonstop_amino_acid_codons_dict.items():
                present_codons = set(reference_codon_frequency_df.columns.intersection(codons))
                if (reference_codon_frequency_df[present_codons].sum().sum()
                    < reference_exclude_amino_acid_count):
                    removed_amino_acids.append(amino_acid)
                    removed_codons += list(present_codons)
            reference_codon_frequency_df = reference_codon_frequency_df.drop(removed_codons, axis=1)
            if removed_amino_acids:
                self.run.warning(
                    "The following amino acids do not meet the threshold of "
                    f"{reference_exclude_amino_acid_count} codons in the reference gene set, and "
                    "so the codons were excluded from the reference: "
                    f"{', '.join(removed_amino_acids)}")

        # Check that there are enough total codons in the reference set.
        if reference_min_analyzed_codons:
            total_codon_count = reference_codon_frequency_df.sum().sum()
            if total_codon_count < reference_min_analyzed_codons:
                reference_codon_frequency_df = reference_codon_frequency_df.drop('all')
                self.run.warning(
                    f"The reference gene set has {pp(int(total_codon_count))} codons, not meeting "
                    f"the minimum codon threshold of {pp(int(reference_min_analyzed_codons))}, so "
                    "a reference codon composition could not be established. Reference-dependent "
                    "CUB tables will be empty.")

        return reference_codon_frequency_df


    def _get_cai_table(self,
                       query_codon_frequency_df,
                       reference_codon_frequency_df,
                       query_min_analyzed_codons=default_query_min_analyzed_codons,
                       reference_min_analyzed_codons=default_reference_min_analyzed_codons):
        """
        Get a table of CAI (Sharp and Li, 1987) values for each query x reference comparison.

        Calculation of CAI:
        reference_codon_weight =
            ln(reference_codon_frequency / reference_max_synonymous_codon_frequency)
        weighted_codon_count = ∑(codon_frequency * reference_codon_weight)
        CAI = exp(weighted_codon_count / codon_count)

        CAI is maximum (1) when all of the codons in the query are the most abundant synonymous
        codons in the reference, and minimum (0-1) when all of the codons in the query are the least
        abundant synonymous codons in the reference.

        Parameters
        ==========
        query_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per item (gene/function) to be compared to each reference
            codon composition and a column per codon.
        reference_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per reference composition and a column per codon.
        query_min_analyzed_codons : int, optional
            A row of the query table must contain at least this number of codons with a reference
            codon weight for CAI to be calculated. By default 0.
        reference_min_analyzed_codons : int, optional
            A reference must contain at least this number of codons that are also present among the
            queries for it to be used in calculating CUB. By default 0.

        Returns
        =======
        pandas.core.frame.DataFrame
            This CUB table has the same row index as the input query table and a column per input
            reference. With a single reference codon composition, this table has a single column.
        """
        if len(reference_codon_frequency_df) == 0:
            return pd.DataFrame()

        # Only consider codons present in both query and reference datasets.
        nonzero_query_codon_frequency_df = query_codon_frequency_df.loc[
            :, query_codon_frequency_df.sum() > 0]
        nonzero_reference_codon_frequency_df = reference_codon_frequency_df.loc[
            :, reference_codon_frequency_df.sum() > 0]
        shared_codons = nonzero_query_codon_frequency_df.columns.intersection(
            nonzero_reference_codon_frequency_df.columns)
        query_codon_frequency_df = nonzero_query_codon_frequency_df[shared_codons]
        reference_codon_frequency_df = nonzero_reference_codon_frequency_df[shared_codons]
        # Remove queries and references that do not have enough codons remaining. When this method
        # is called from `self.get_codon_usage_bias`, a defined single-row reference will have
        # already been filtered, rendering this reference filter unneeded; on the other hand,
        # omnibias references will be the same as the query references, and
        # `query_min_analyzed_codons` was ensured to be the same as `reference_min_analyzed_codons`,
        # so the same queries and references are here removed for not meeting the same threshold.
        query_codon_frequency_df = query_codon_frequency_df[
            query_codon_frequency_df.sum(axis=1) >= query_min_analyzed_codons]
        reference_codon_frequency_df = reference_codon_frequency_df[
            reference_codon_frequency_df.sum(axis=1) >= reference_min_analyzed_codons]

        # Calculate reference codon weights.
        np.seterr(divide='ignore')
        synonymous_weight_dfs = []
        for codons in self.synonymous_nonstop_amino_acid_codons_dict.values():
            # Not every codon in the synonymous set need be present in the data.
            intermediate_df = reference_codon_frequency_df[
                reference_codon_frequency_df.columns.intersection(codons)]
            synonymous_weight_df = np.log(
                intermediate_df.div(intermediate_df.max(axis=1), axis=0))
            synonymous_weight_dfs.append(synonymous_weight_df)
        np.seterr(divide='warn')
        weight_df = pd.concat(synonymous_weight_dfs, axis=1)
        # In omnibias mode, codons may be missing from one but not all references. Zero frequencies
        # result in weights of -∞, which are replaced by NaN.
        weight_df = weight_df.replace(-np.inf, np.nan)
        weight_df = weight_df[sorted(weight_df.columns)]
        weight_array = weight_df.values
        bool_weight_array = weight_df.notna().values

        performed_comparisons = 0
        total_comparisons = len(query_codon_frequency_df) * len(weight_df)
        self.progress.new("Calculating CAI")

        cai_rows = []
        for query_codon_frequencies in query_codon_frequency_df.values:
            cai_row = []
            for weights, bool_weights in zip(weight_array, bool_weight_array):
                if performed_comparisons % 10000 == 0:
                    self.progress.update(
                        f"{performed_comparisons} / {total_comparisons} comparisons")

                query_codon_frequency_with_reference = query_codon_frequencies[bool_weights].sum()
                if query_codon_frequency_with_reference < query_min_analyzed_codons:
                    cai_row.append(np.nan)
                    performed_comparisons += 1
                    continue
                weighted_codon_count = np.nansum(query_codon_frequencies * weights)
                cai_row.append(np.exp(weighted_codon_count / query_codon_frequency_with_reference))
                performed_comparisons += 1
            cai_rows.append(cai_row)
        cai_df = pd.DataFrame(cai_rows, index=query_codon_frequency_df.index)
        cai_df = cai_df.reindex(query_codon_frequency_df.index)

        self.progress.end()

        return cai_df


    def _get_delta_table(self,
                         query_codon_frequency_df,
                         reference_codon_frequency_df,
                         nonreference_codon_frequency_df,
                         query_min_analyzed_codons=default_query_min_analyzed_codons,
                         reference_min_analyzed_codons=default_reference_min_analyzed_codons):
        """
        Get a table of 𝛿 (Ran and Higgs, 2012, Eq. 6) values for each query x reference comparison.

        Calculation of 𝛿:
        codon_weight = ln(reference_codon_synonymous_relative_frequency /
                          nonreference_codon_synonymous_relative_frequency)
        weighted_codon_count = ∑(codon_frequency * reference_codon_weight)
        𝛿 = weighted_codon_count / codon_count

        𝛿 is maximum (1) when all of the codons in the query are the most abundant synonymous codons
        in the reference, and minimum (0-1) when all of the codons in the query are the least
        abundant synonymous codons in the reference.

        𝛿 codon weights are ratios of reference to genome-wide synonymous relative frequencies. 𝛿
        ranges from (-∞, +∞), with more positive values being more similar to the reference. When
        query and reference sets are the same, 𝛿 is a likelihood ratio evaluating the distinctness
        of the set from the genome as a whole.

        Parameters
        ==========
        query_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per item (gene/function) to be compared to each reference
            codon composition and a column per codon.
        reference_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per reference composition and a column per codon.
        nonreference_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a single row of summed frequencies across all non-reference
            genes in the genome, with a column per codon.
        query_min_analyzed_codons : int, optional
            A row of the query table must contain at least this number of codons with a reference
            codon weight for 𝛿 to be calculated. By default 0.
        reference_min_analyzed_codons : int, optional
            A reference must contain at least this number of codons that are also present among the
            queries for it to be used in calculating CUB. By default 0.

        Returns
        =======
        pandas.core.frame.DataFrame
            This CUB table has the same row index as the input query table and a column per input
            reference. With a single reference codon composition, this table has a single column of
            data.
        """
        if len(reference_codon_frequency_df) == 0:
            return pd.DataFrame()

        # Calculate codon weights from the reference and non-reference gene sets.
        # First get synonymous relative frequencies from each. Do this before dropping any codons
        # missing entirely from either queries or references to ensure the accuracy of synonymous
        # relative frequencies. For example, if CysTGT is missing from the queries but not the
        # references, then the synonymous relative frequency of CysTGC in the references should
        # still be calculated as TGC / (TGC + TGT).
        reference_synonymous_codon_rel_frequency_df = \
            self._get_synonymous_codon_rel_frequency_table(reference_codon_frequency_df)
        nonreference_synonymous_codon_rel_frequency_df = \
            self._get_synonymous_codon_rel_frequency_table(nonreference_codon_frequency_df)

        # Only consider codons present in both query and reference datasets.
        nonzero_query_codon_frequency_df = query_codon_frequency_df.loc[
            :, query_codon_frequency_df.sum() > 0]
        nonzero_reference_codon_frequency_df = reference_codon_frequency_df.loc[
            :, reference_codon_frequency_df.sum() > 0]
        shared_codons = nonzero_query_codon_frequency_df.columns.intersection(
            nonzero_reference_codon_frequency_df.columns)
        query_codon_frequency_df = nonzero_query_codon_frequency_df[shared_codons]
        reference_codon_frequency_df = nonzero_reference_codon_frequency_df[shared_codons]
        reference_synonymous_codon_rel_frequency_df = \
            reference_synonymous_codon_rel_frequency_df[shared_codons]
        nonreference_synonymous_codon_rel_frequency_df = \
            nonreference_synonymous_codon_rel_frequency_df[shared_codons]
        nonreference_synonymous_codon_rel_frequencies = \
            nonreference_synonymous_codon_rel_frequency_df.squeeze().values
        # Remove queries and references that do not have enough codons remaining. When this method
        # is called from `self.get_codon_usage_bias`, a defined single-row reference will have
        # already been filtered, rendering this reference filter unneeded; on the other hand,
        # omnibias references will be the same as the query references, and
        # `query_min_analyzed_codons` was ensured to be the same as `reference_min_analyzed_codons`,
        # so the same queries and references are here removed for not meeting the same threshold.
        query_codon_frequency_df = query_codon_frequency_df[
            query_codon_frequency_df.sum(axis=1) >= query_min_analyzed_codons]
        reference_codon_frequency_df = reference_codon_frequency_df[
            reference_codon_frequency_df.sum(axis=1) >= reference_min_analyzed_codons]
        reference_synonymous_codon_rel_frequency_df = \
            reference_synonymous_codon_rel_frequency_df.loc[reference_codon_frequency_df.index]

        np.seterr(divide='ignore')
        weight_df = np.log(reference_synonymous_codon_rel_frequency_df.div(
            nonreference_synonymous_codon_rel_frequencies, axis=1))
        np.seterr(divide='warn')
        # In omnibias mode, codons may be missing from one but not all references. Zero frequencies
        # result in weights of -∞, which are replaced by NaN.
        weight_df = weight_df.replace(-np.inf, np.nan)
        weight_df = weight_df[sorted(weight_df.columns)]
        weight_array = weight_df.values
        bool_weight_array = weight_df.notna().values

        performed_comparisons = 0
        total_comparisons = len(query_codon_frequency_df) * len(weight_df)
        self.progress.new("Calculating 𝛿")

        delta_rows = []
        for query_codon_frequencies in query_codon_frequency_df.values:
            delta_row = []
            for weights, bool_weights in zip(weight_array, bool_weight_array):
                if performed_comparisons % 10000 == 0:
                    self.progress.update(
                        f"{performed_comparisons} / {total_comparisons} comparisons")

                query_codon_frequency_with_reference = query_codon_frequencies[bool_weights].sum()
                if query_codon_frequency_with_reference < query_min_analyzed_codons:
                    delta_row.append(np.nan)
                    performed_comparisons += 1
                    continue
                weighted_codon_count = np.nansum(query_codon_frequencies * weights)
                delta_row.append(weighted_codon_count / query_codon_frequency_with_reference)
                performed_comparisons += 1
            delta_rows.append(delta_row)
        delta_df = pd.DataFrame(delta_rows, index=query_codon_frequency_df.index)
        delta_df = delta_df.reindex(query_codon_frequency_df.index)

        self.progress.end()

        return delta_df


class MultiGenomeCodonUsage(object):
    """This object processes codon usage data from multiple internal and/or external genomes."""

    def __init__(self, args, r=run, rq=run_quiet, p=progress):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.contigs_db = A('contigs_db')
        self.profile_db = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_ids_path = A('bin_ids_file')

        self.internal_genomes_path = A('internal_genomes')
        self.external_genomes_path = A('external_genomes')

        self.function_sources = A('function_sources')
        self.all_brite_categories = A('all_brite_categories')
        self.use_shared_function_sources = A('shared_function_sources')
        if self.use_shared_function_sources is None:
            self.use_shared_function_sources = False

        self.preload_genomes = A('preload_genomes')
        if self.preload_genomes is None:
            self.preload_genomes = False

        self.run = r
        self.run_quiet = rq
        self.progress = p

        if self.contigs_db and (self.internal_genomes_path or self.external_genomes_path):
            raise ConfigError(
                "Internal and external genomes files cannot be provided alongside bins in a "
                "collection from a profile database. `internal_genomes` and `external_genomes` can "
                "be used by themselves or together, and `contigs_db`, `profile_db`, "
                "`collection_name`, and `bin_ids_file` can be used in various combinations.")

        if self.collection_name:
            self.bin_ids = []
            collections = ccollections.Collections()
            collections.populate_collections_dict(self.profile_db)
            if self.bin_ids_path:
                filesnpaths.is_file_plain_text(self.bin_ids_path)
                for line in open(self.bin_ids_path):
                    bin_id = line.rstrip()
                    self.bin_ids.append(bin_id)
                self.bin_ids = list(set(self.bin_ids))
                self.bin_ids.remove('')
            else:
                self.bin_ids = list(collections.get_bins_info_dict(self.collection_name))
        else:
            self.bin_ids = None
            collections = None

        # Store genome information in the following dictionary.
        self.genome_info_dict = {}
        if self.internal_genomes_path or self.external_genomes_path:
            # Store genome information starting with internal and external genome files.
            descriptions = GenomeDescriptions(args, run=self.run_quiet, progress=self.progress)
            descriptions.load_genomes_descriptions(init=False)

            if self.function_sources is None:
                pass
            elif (self.function_sources == [] and
                  self.use_shared_function_sources and
                  descriptions.function_annotation_sources_some_genomes_miss):
                # All function sources are requested (in `anvi-get-codon-usage` this corresponds to
                # `--function-sources` used as a flag), but some sources are not common to all
                # genomes and only common sources are allowed (`use_shared_function_sources`).
                # Report sources in common that will be analyzed and those not in common that will
                # be ignored.
                self.function_sources = list(descriptions.function_annotation_sources)
                self.run.info("Function sources shared across genomes: ",
                            ', '.join(self.function_sources))
                self.run.info("Function sources missing in one or more genomes: ",
                            ', '.join(descriptions.function_annotation_sources_some_genomes_miss))
            elif (len(self.function_sources) and
                  self.use_shared_function_sources and
                  descriptions.function_annotation_sources_some_genomes_miss):
                # A list of function sources is requested, but some are not common to all genomes
                # and only common sources are allowed (`use_shared_function_sources`). This causes
                # an error.
                unshared_function_sources = []
                for source in self.function_sources:
                    if source not in descriptions.function_annotation_sources_some_genomes_miss:
                        unshared_function_sources.append(source)
                if unshared_function_sources:
                    raise ConfigError(
                        "Of the requested function annotation sources, the following were not run "
                        f"on every genome: {', '.join(unshared_function_sources)}.")

            # Store information on the genomes.
            for genome_name, genome_dict in descriptions.internal_genomes_dict.items():
                self.genome_info_dict[genome_name] = genome_info = {}
                genome_info['contigs_db'] = genome_dict['contigs_db_path']
                genome_info['profile_db'] = genome_dict['profile_db_path']
                genome_info['collection_name'] = genome_dict['collection_id']
                genome_info['bin_id'] = genome_dict['bin_id']
            for genome_name, genome_dict in descriptions.external_genomes_dict.items():
                self.genome_info_dict[genome_name] = genome_info = {}
                genome_info['contigs_db'] = genome_dict['contigs_db_path']
        else:
            # Store genome information from bins in a collection for a single contigs database.
            if self.collection_name not in collections.collections_dict:
                raise ConfigError(f"The collection, '{self.collection_name}', was not found in the "
                                  f"profile database, '{self.profile_db}'.")

            missing_bin_ids = []
            for bin_id in self.bin_ids:
                try:
                    collections.is_bin_in_collection(self.collection_name, bin_id)
                except ConfigError:
                    missing_bin_ids.append(bin_id)
            if missing_bin_ids:
                raise ConfigError(
                    "The following bin IDs were not found in the collection, "
                    f"'{self.collection_name}', in the profile database, '{self.profile_db}': "
                    f"{', '.join(missing_bin_ids)}")

            for bin_id in self.bin_ids:
                self.genome_info_dict[bin_id] = genome_info = {}
                genome_info['contigs_db'] = self.contigs_db
                genome_info['profile_db'] = self.profile_db
                genome_info['collection_name'] = self.collection_name
                genome_info['bin_id'] = bin_id

        # The following information is common to all genomes regardless of input source.
        for genome_info in self.genome_info_dict.values():
            genome_info['function_sources'] = self.function_sources
            genome_info['all_brite_categories'] = self.all_brite_categories
            genome_info['codon_to_amino_acid'] = self.args.codon_to_amino_acid
            genome_info['ignore_start_codons'] = self.args.ignore_start_codons

        # There are memory- and CPU-efficient ways of setting up the object. `preload_genomes` loads
        # all of the SingleGenomeCodonUsage objects into memory, allowing methods to save time by
        # immediately accessing their frequency tables. When this argument is False,
        # SingleGenomeCodonUsage objects are loaded into memory as needed when each genome is
        # processed.
        self.genome_codon_usage_dict = None
        if self.preload_genomes:
            self.genome_codon_usage_dict = {}
            for genome_name, genome_info in self.genome_info_dict.items():
                self.genome_codon_usage_dict[genome_name] = SingleGenomeCodonUsage(
                    argparse.Namespace(**genome_info))


    def get_frequencies(self,
                        from_function_sources=False,
                        return_functions=False,
                        return_amino_acids=False,
                        function_accessions=None,
                        function_names=None,
                        expect_functions=False,
                        relative=False,
                        synonymous=False,
                        sum_genes=False,
                        average_genes=False,
                        gene_min_codons=0,
                        function_min_codons=0,
                        min_codon_filter='both',
                        drop_amino_acids=None,
                        sequence_min_amino_acids=0,
                        pansequence_min_amino_acids=(0, 1.0),
                        label_amino_acids=False,
                        infinity_to_zero=False):
        """
        Get absolute (default) or relative codon or amino acid frequencies from genes or functions
        in one or more genomes.

        Parameters
        ==========
        See the `SingleGenomeCodonUsage.get_frequencies` docstring for descriptions of each
        parameter.

        Returns
        =======
        pandas.core.frame.DataFrame
            Frequency table of gene x codon or amino acid. The first column of the MultiIndex of the
            returned DataFrame is the genome name, with blocks of results for each genome being in
            the input order of internal then external genomes. If functions are not considered, then
            the other column of the MultiIndex contains gene caller IDs. If functions are considered
            and `return_functions` is False, then each row represents a gene/function pair, and the
            same gene frequencies can be found in multiple rows with different function pairs; there
            are additional MultiIndex columns for source, accession, and name of function. If
            functions are considered and `return_functions` is True, then each row represents a
            function rather than gene/function pair. If frequencies are summed or averaged and
            functions are not considered, the returned single-row DataFrame has an Index with one
            entry, 'all'. If frequencies are summed or averaged and functions are considered, the
            returned DataFrame has a row for each function annotation source.
        """
        kwargs = {}
        arg_info = inspect.getargvalues(inspect.currentframe())
        for param in arg_info.args:
            if param == 'self':
                continue
            kwargs[param] = arg_info.locals[param]
        frequency_table_generator = self._get_genome_frequency_table(kwargs)
        frequency_dfs = []
        for genome_name in self.genome_info_dict:
            self.run.info("Genome", genome_name)
            frequency_df = next(frequency_table_generator)

            frequency_df.insert(0, 'genome_name', genome_name)
            new_index_columns = ['genome_name'] + frequency_df.index.names
            frequency_df = frequency_df.reset_index().set_index(new_index_columns)

            frequency_dfs.append(frequency_df)
        return pd.concat(frequency_dfs)


    def _get_genome_frequency_table(self, kwargs):
        """This generator yields a frequency table from each genome."""
        for genome_name, genome_info in self.genome_info_dict.items():
            if self.preload_genomes:
                genome_codon_usage = self.genome_codon_usage_dict[genome_name]
            else:
                genome_codon_usage = SingleGenomeCodonUsage(
                    argparse.Namespace(**genome_info))
            frequency_df = genome_codon_usage.get_frequencies(**kwargs)
            print()

            yield frequency_df


    def get_codon_usage_bias(
        self,
        metrics=None,
        from_function_sources=None,
        function_accessions=None,
        function_names=None,
        expect_functions=False,
        omnibias=False,
        reference_function_accessions=None,
        reference_function_names=None,
        expect_reference_functions=False,
        gene_min_codons=0,
        function_min_codons=0,
        min_codon_filter='both',
        drop_amino_acids=None,
        sequence_min_amino_acids=0,
        pansequence_min_amino_acids=(0, 1.0),
        query_min_analyzed_codons=default_query_min_analyzed_codons,
        reference_exclude_amino_acid_count=default_reference_exclude_amino_acid_count,
        reference_min_analyzed_codons=default_reference_min_analyzed_codons):
        """This generator yields a genome name and CUB table dict from each genome.

        See the `SingleGenomeCodonUsage.get_codon_usage_bias` docstring for descriptions of each
        parameter.
        """
        # Rather than individually listing a slew of arguments when calling `get_codon_usage_bias`,
        # package them in a tidy kwargs dictionary.
        kwargs = {}
        arg_info = inspect.getargvalues(inspect.currentframe())
        for param in arg_info.args:
            if param == 'self':
                continue
            kwargs[param] = arg_info.locals[param]

        for genome_name, genome_info in self.genome_info_dict.items():
            if self.preload_genomes:
                genome_codon_usage = self.genome_codon_usage_dict[genome_name]
            else:
                genome_codon_usage = SingleGenomeCodonUsage(argparse.Namespace(**genome_info))
            cub_table_dict = genome_codon_usage.get_codon_usage_bias(**kwargs)
            print()

            yield genome_name, cub_table_dict


def get_custom_encodings(encodings_txt):
    """Modify the standard genetic code given user input, as used in the scripts,
    `anvi-get-codon-frequencies` and `anvi-get-codon-usage-bias`."""
    codon_amino_acid_dict = copy.deepcopy(default_codon_amino_acid_dict)

    if not encodings_txt:
        return codon_amino_acid_dict

    filesnpaths.is_file_tab_delimited(encodings_txt)
    encodings_df = pd.read_csv(encodings_txt, sep='\t', header=None)
    codon_amino_acid_dict.update(dict(zip(encodings_df.iloc[:, 0], encodings_df.iloc[:, 1])))

    return codon_amino_acid_dict


def check_genetic_code(codon_amino_acid_dict):
    """Check that known codons and three-letter amino acid codes are used in a dict defining the
    genetic code, throwing an exception if needed."""
    unrecognized_codons = []
    unrecognized_amino_acids = []
    for codon, amino_acid in codon_amino_acid_dict.items():
        if codon not in constants.codons:
            unrecognized_codons.append(codon)
        if amino_acid not in constants.amino_acids:
            unrecognized_amino_acids.append(amino_acid)

    if unrecognized_codons:
        unrecognized_codon_message = (
            "The following codons in the provided genetic code are not recognized: "
            f"{', '.join(unrecognized_codons)}.")
    else:
        unrecognized_codon_message = ""

    if unrecognized_amino_acids:
        unrecognized_amino_acid_message = (
            "The following amino acids in the provided genetic code are not recognized: "
            f"{', '.join(unrecognized_amino_acids)}. These should be three-letter codes "
            "and \"STP\" for stop codons.")
        if unrecognized_codon_message:
            unrecognized_amino_acid_message = " " + unrecognized_amino_acid_message
    else:
        unrecognized_amino_acid_message = ""

    if unrecognized_codon_message or unrecognized_amino_acid_message:
        raise ConfigError(f"{unrecognized_codon_message}{unrecognized_amino_acid_message}")


def write_split_codon_output(codon_frequency_df,
                             output_path,
                             separate_genomes=False,
                             separate_function_sources=False,
                             output_root=None):
    """
    Store tables of codon frequency data.

    Parameters
    ----------
    codon_frequency_df : pandas.core.frame.DataFrame
        A table of codon frequencies of the format produced by
        `SingleGenomeCodonUsage.get_frequencies` or `MultiGenomeCodonUsage.get_frequencies`.
    output_path : str
        This filepath serves as a template for the paths of other files that may be generated.
        Indeed, if the arguments, `separate_genomes` and/or `separate_function_sources` are used,
        nothing is written to this path. Given the template, <root><extension>, substrings can be
        inserted in the derived filepaths, with the most elaborate such filepath being
        <root>-<genome_name>-<function_source>-<normalization_method><extension>.
    separate_genomes : bool, optional
        If True (default False), split output tables by genome. In the output paths, underscores
        replace spaces in the genome names.
    separate_function_sources : bool, optional
        If True (default False), split output tables by function source.
    output_root : str, optional
        If given (default None), then rather than splitting the output path by the extension, the
        output filepath template consists of <output_root><post_output_root>, with <genome_name>
        and <function_source> inserted between the two.
    """
    if separate_genomes:
        if 'genome' not in codon_frequency_df.index.names:
            raise ConfigError(
                "To split output by genome, as requested by `separate_genomes`, the table of "
                "codon frequencies must contain the expected index column, \"genome\".")

    if separate_function_sources:
        if 'function_source' not in codon_frequency_df.index.names:
            raise ConfigError(
                "To split output by function source, as requested by `separate_function_sources`, "
                "the table of codon frequencies must contain the expected index column, "
                "\"function_source\".")

    if output_root is not None:
        try:
            output_root_index = output_path.index(output_root)
        except ValueError:
            output_root_index = -1
        if output_root_index != 0:
            raise ConfigError(
                "The value provided for `output_root`, '{output_root}', does not but must occur at "
                "the start of `output_path`, '{output_path}'.")

    if output_root:
        output_ext = output_path[len(output_root): ]
    else:
        output_root, output_ext = os.path.splitext(output_path)
    valid_derived_output_paths = []
    invalid_derived_output_paths = []
    is_given_output_path_valid = True
    if separate_genomes and separate_function_sources:
        # Split the data by both genome and function source.
        for genome_name, genome_df in codon_frequency_df.groupby('genome_name'):
            for source, source_df in genome_df.groupby('function_source'):
                # Write a table of codon frequencies for the genome/source.
                derived_output_path = (
                    output_root + "-" +
                    genome_name.replace(" ", "_") + "-" + source +
                    output_ext)
                try:
                    filesnpaths.is_output_file_writable(derived_output_path)
                    valid_derived_output_paths.append(derived_output_path)
                except FilesNPathsError:
                    invalid_derived_output_paths.append(derived_output_path)
                    continue
                source_df.to_csv(derived_output_path, sep='\t')
    elif separate_genomes:
        # Split the data by genome.
        for genome_name, genome_df in codon_frequency_df.groupby('genome_name'):
            derived_output_path = output_root + "-" + genome_name.replace(" ", "_") + output_ext
            try:
                filesnpaths.is_output_file_writable(derived_output_path)
                valid_derived_output_paths.append(derived_output_path)
            except FilesNPathsError:
                invalid_derived_output_paths.append(derived_output_path)
                continue
            genome_df.to_csv(derived_output_path, sep='\t')
    elif separate_function_sources:
        # Split the data by function source.
        for source, source_df in codon_frequency_df.groupby('function_source'):
            derived_output_path = output_root + "-" + source + output_ext
            try:
                filesnpaths.is_output_file_writable(derived_output_path)
            except FilesNPathsError:
                invalid_derived_output_paths.append(derived_output_path)
                continue
            valid_derived_output_paths.append(derived_output_path)
            source_df.to_csv(derived_output_path, sep='\t')
    else:
        try:
            filesnpaths.is_output_file_writable(output_path)
            codon_frequency_df.to_csv(output_path, sep='\t')
            is_given_output_path_valid = 1
        except FilesNPathsError:
            is_given_output_path_valid = False

    # Report invalid output paths and, if some paths were invalid, also remove files written to
    # valid paths.
    if is_given_output_path_valid:
        template_output_path_message = ""
    else:
        template_output_path_message = ("The table of codon frequencies could not be written to "
                                        f"`output_path`, '{output_path}'.")
    if invalid_derived_output_paths:
        paths_string = ', '.join(['\'' + path + '\'' for path in invalid_derived_output_paths])
        derived_output_paths_message = (
            "Tables of codon frequencies could not be written to the following derived target "
            f"paths: {paths_string}.")
    else:
        derived_output_paths_message = ""
    if invalid_derived_output_paths or not is_given_output_path_valid:
        if is_given_output_path_valid is 1:
            # Derived output paths were invalid, and raw affinities were written to the input
            # output file.
            os.remove(output_path)
        for derived_output_path in valid_derived_output_paths:
            os.remove(derived_output_path)
        raise ConfigError(
            f"{template_output_path_message}"
            f"{' ' if invalid_derived_output_paths and not is_given_output_path_valid else ''}"
            f"{derived_output_paths_message}")
