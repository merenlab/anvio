#!/usr/bin/env python
# -*- coding: utf-8
"""
Calculate various codon frequency statistics.
"""

from __future__ import annotations

import pandas as pd

from typing import Literal
from functools import partial

import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError
from anvio.geneticcode import GeneticCode
from anvio.codonusage.genomefrequencies import GenomeCodonFrequencies

pp = terminal.pretty_print

class CodonFrequencyCalculator:
    """
    Calculates codon frequency statistics.

    Attributes
    ==========
    genome_codon_frequencies : anvio.codonusage.genomefrequencies.GenomeCodonFrequencies
        Stores raw codon frequencies of gene sequences.

    genetic_code : anvio.geneticcode.GeneticCode, GeneticCode()
        Genetic code relating codons to encoded amino acids.

    genomic_context : anvio.codonusage.genomiccontext.GenomicContext, None
        Contains genomic sequence data loaded from a contigs database.

    auto_adjust : bool, True
        If True, adjust the values of certain attributes given the values of other attributes in
        option setup. For example, if 'synonymous' is True, then 'relative' should also be True. If
        'relative' is False, then with auto adjustment turned on, 'relative' is changed to True
        during setup, whereas with auto adjustment turned off, an error is thrown. This simplifies
        parameterization.

    relative : bool, False
        If True, return frequency tables containing relative rather than absolute codon frequencies.

    synonymous : bool, False
        If True, return frequency tables containing relative frequencies among synonymous codons
        decoding each amino acid plus stop codons ("STP"), i.e., relative synonymous codon usage
        (RSCU).

    return_amino_acids : bool, False
        If True, returned frequency table columns are decoded amino acids plus stop codons ("STP")
        rather than codons. Synonymous codon frequencies are summed to produce the amino acid
        frequencies.

    label_amino_acids : bool, False
        If True, include the amino acid for each codon in the column header of returned frequency
        tables, i.e., 'LysAAA' instead of 'AAA'.

    infinity_to_zero : bool, False
        If True, replace NaN (empty) values in returned statistical tables with 0.0. NaN occurs if
        'synonymous' is True and all codons for an amino acid are absent in a gene or function,
        resulting in 0/0, reported as NaN. Use with caution, for NaN and 0.0 mean different things,
        and this will skew downstream analyses of RSCU, such as codon usage bias.

    gene_caller_ids : list[str], []
        Genes with the given IDs are selected for analysis. This parameter can be used alongside
        'function_accessions' and 'function_names'.

    function_sources : list[str], []
        Select genes with functional annotation sources in the list, e.g., ['KOfam'] or ['KOfam',
        'COG20_FUNCTION']. If the list contains items, the first four columns of returned
        statistical tables contain, respectively, function annotation sources, function accessions,
        function names, and gene caller IDs. There is a row for each gene/function combination in
        the table, and each row for the same gene with different functional annotations contains the
        same frequency values. When 'KEGG_BRITE' is in the list, in the absence of the
        'function_names' parameter narrowing the scope of the inquiry, all possible depths of each
        BRITE hierarchy in the data are returned, e.g., 'Ribosome>>>Ribosomal proteins' and the more
        general 'Ribosome' would each have rows in the output table.

    function_accessions : dict[str, list[str]], None
        Genes annotated with the given function accessions are used in analysis. Dictionary keys are
        function annotation sources (e.g., 'KOfam', 'COG20_FUNCTION') which must be available in
        'function_sources'. Values are lists of function accessions. Note that 'KEGG_BRITE'
        accessions are not individual function accessions but overarching hierarchy accessions that
        include many functions.

    function_names : dict[str, list[str]], None
        Genes annotated with the given function names are used in analysis. Dictionary keys are
        function annotation sources (e.g., 'KOfam', 'COG20_FUNCTION') which must be available in
        'function_sources'. Values are lists of function names. If 'KEGG_BRITE' is used as a source,
        names (categorizations) should be given to an arbitrary level of the hierarchy, like
        ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic type>>>Bacterial
        type>>>RNA polymerase'].

    expect_functions : bool, False
        If True, an error will be raised in getting statistical tables if any given
        'function_accessions' or 'function_names' values are not annotated in the input genome.

    return_functions : bool, False
        If True, returned statistical tables contain function rather than gene results in each row.
        Function results are aggregated from all genes annotated with the function. If False,
        per-gene results can still be returned while still considering functions with the
        'function_sources' parameter, which facilitates analysis of variability within functions.

    sum_genes : bool, False
        If True, return frequency tables containing a single row of codon frequencies summed across
        all genes. Normalization to relative/synonymous values occurs after summing absolute
        frequencies. With 'auto_adjust': 'function_sources' are ignored with this option, as
        frequencies are summed irrespective of function groups; 'return_functions' ignores this
        option, by definition summing frequencies among genes in each function.

    average_genes_in_function : bool, False
        If True, return average frequencies of genes in each function. 'relative' or 'synonymous'
        must be True when using this option, as the average only operates on relative frequencies so
        that genes are normalized for length. With 'auto_adjust': 'return_functions' is treated as
        True, returning a table of averages per function.

    std_genes_in_function : bool, False
        If True, return standard deviations of gene frequencies in each function. 'relative' or
        'synonymous' must be True when using this option, as the standard deviation only operates on
        relative frequencies so that genes are normalized for length. With 'auto_adjust':
        'return_functions' is treated as True, returning a table of standard deviations per
        function.

    gene_min_codons : int, 0
        Ignore genes with fewer than this number of codons in analyses.

    function_min_codons : int, 0
        If functions are considered, ignore gene/function pairs when the function has fewer than
        this number of codons in analyses. 'gene_min_codons' filters genes by codon count before
        functions are filtered with this option, so the total length of remaining genes annotated by
        the function is considered.

    min_codon_filter : Literal['length', 'remaining', 'both'], 'both'
        This parameter arises from the ambiguity of (1) analytical filters that remove genes and
        functions by number of codons ('gene_min_codons' and 'function_min_codons') in relation to
        (2a) the filters that drop codons from the analysis ('drop_codons', 'drop_amino_acids',
        'pansequence_min_codons', and 'pansequence_min_amino_acids') and (2b) the filters that drop
        codons from the reported table ('unreported_codons', 'unreported_amino_acids',
        'report_min_codons', and 'report_min_amino_acids').

        The possible values of this attribute are 'length', 'remaining', and 'both'.

        'length': Genes (and functions) can be filtered by their full length (concatenated length of
        annotated genes for functions) with this value, e.g., genes shorter than 300 codons are
        ignored if 'gene_min_codons' is set to 300.

        'remaining': Genes and functions can also be filtered by the number of codons remaining
        after dropping codons by providing the attribute value, 'remaining'. Applying the codon
        'length' filter followed by dropping codons can result in genes and functions with fewer
        codons than the length threshold, thus 'remaining' is a filter option in addition to
        'length' to ensure that total codon frequencies in the output always meet the minimum codon
        threshold.

        'both': Lastly, 'both' is needed as a filter option in addition to 'remaining' so that
        dynamic codon filtering by 'pansequence_min_codons' and 'pansequence_min_amino_acids' occurs
        after rather than before the first gene/function 'length' filter. See documentation in the
        'anvi-get-codon-frequencies' artifact for more information, including an example of the
        necessity of the 'both' filter.

    drop_codons : list[str], []
        Remove codons from analyses. CAUTION! Note how this affects the calculation of relative
        frequencies using 'relative' and 'synonymous' options: codon frequencies are normalized to
        the summed frequencies of retained codons. In synonymous calculations, retained codons in an
        amino acid set are considered, causing affected synonymous codon relative frequencies to
        deviate from the standard meaning of RSCU.

    drop_amino_acids : list[str], []
        Remove codons that decode the given amino acids from analyses. Specify three-letter codes,
        e.g., 'Ala', and 'STP' for stop codons. If 'synonymous' is True, the 'drop_amino_acids'
        default changes from an empty list to a list of the non-degenerate amino acids ('Met' and
        'Trp' in the standard genetic code) plus 'STP' rather than an empty list.

    pansequence_min_codons : int, (0, 0.0)
        This tuple must have two items. The first is a non-negative integer representing a minimum
        number of codons: call this 'min_codons'. The second is a float in the interval [0.0, 1.0]
        representing a fraction of genes: call this 'min_gene_fraction'. Remove codons from the
        analysis that are less numerous than 'min_codons' in a 'min_gene_fraction' of genes. For
        example, if 'min_codons' is 3 and 'min_gene_fraction' is 0.9, then if there are <3 of a
        codon in ≥90% of genes, columns for these codons are dropped. CAUTION! Note how this affects
        the calculation of relative frequencies using 'relative' and 'synonymous' options: codon
        frequencies are normalized to the summed frequency of retained codons. In synonymous
        calculations, retained codons in an amino acid set are considered, causing affected
        synonymous codon relative frequencies to deviate from the standard meaning of RSCU.

    pansequence_min_amino_acids : tuple, (0, 0.0)
        This tuple must have two items. The first is a non-negative integer representing a minimum
        number of codons encoding an amino acid: call this 'min_amino_acids'. The second is a float
        in the interval [0.0, 1.0] representing a fraction of genes: call this 'min_gene_fraction'.
        Remove codons from the analysis for amino acids (and 'STP') that are less numerous than
        'min_amino_acids' in a 'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5
        and 'min_gene_fraction' is 0.9, then if there are <5 codons for an amino acid/'STP' in ≥90%
        of genes, then the columns for these codons are dropped. CAUTION! Note how this affects the
        calculation of relative frequencies using the 'relative' option: codon frequencies are
        normalized to the summed frequency of retained codons, not all codons.

    unreported_codons : list[str], []
        Remove columns for select codons from returned frequency tables.

    unreported_amino_acids : list[str], []
        Remove columns for codons encoding select amino acids (or columns for select amino acids
        themselves using 'return_amino_acids') from returned frequency tables.

    report_min_codons : int, 0
        Replace codon values with NaN (empty) values in returned frequency tables if the codon has a
        count less than this number in a row. For example, if the attribute value is 3, and a gene
        has 2 AAT codons, then a row for this gene in the output table will have a missing value in
        the AAT column.

    report_min_amino_acids : int, 0
        Replace codon values with NaN (empty) values in returned frequency tables if the set of
        codons for an amino acid has a count less than this number in a row. This option works the
        same with 'return_amino_acids', but for amino acid values. For example, if the attribute
        value is 5, and a gene has 4 codons encoding Asn, 2 AAT and 2 AAC, then a row for this gene
        in the output table will have missing values in Asn columns.

    run : anvio.terminal.Run, anvio.terminal.Run()
        Prints run information to the terminal.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        Prints transient progress information to the terminal.

    verbose : bool, False
        Prints additional information to the terminal, including auto-adjustments made with
        'auto_adjust'.

    kwargs : dict
        All keyword arguments provided at object initialization.
    """
    def __init__(
        self,
        genome_codon_frequencies: GenomeCodonFrequencies,
        genetic_code: GeneticCode = GeneticCode(),
        auto_adjust: bool = True,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        verbose: bool = False,
        **kwargs
    ) -> None:
        """
        Parameters
        ==========
        genome_codon_frequencies : anvio.codonusage.genomefrequencies.GenomeCodonFrequencies
            This object contains raw codon frequencies of gene sequences. The 'genomic_context'
            attribute of this object is also set as the 'genomic_context' attribute of the
            initialized CodonFrequencyCalculator object. A genomic context is required for certain
            functionality, including the analysis of functional codon usage, as the genomic context
            relates genes to functional annotations.

        genetic_code : anvio.geneticcode.GeneticCode, GeneticCode()
            Genetic code relating codons to encoded amino acids and termination function.

        auto_adjust : bool, True
            If True, adjust the values of certain attributes given the values of other attributes in
            option setup. For example, if 'synonymous' is True, then 'relative' should also be True.
            If 'relative' is False, then with auto adjustment turned on, 'relative' is changed to
            True during setup, whereas with auto adjustment turned off, an error is thrown. This
            simplifies parameterization.

        run : anvio.terminal.Run, anvio.terminal.Run()
            Prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            Prints transient progress information to the terminal.

        verbose : bool, False
            Prints additional information to the terminal, including auto-adjustments made with
            'auto_adjust'.

        Keyword arguments
        =================
        Arguments that can be used in setup. Those without descriptions are documented in equivalent
        attributes of the same name.

        relative : bool

        synonymous : bool

        return_amino_acids : bool

        label_amino_acids : bool

        infinity_to_zero : bool

        gene_caller_ids : list[int]

        function_sources : list[str]
            This differs from the attribute of the same name. Without this argument, equivalent to a
            value of None, functional annotation sources are not considered in codon frequency
            calculations and do not factor into output tables: the 'function_sources' attribute
            is set to []. If this argument has a value of [], then all sources available in the
            genomic context are considered in codon frequency calculations, and the
            'function_sources' attribute is set to the list of all sources. Lastly, a subset of
            sources available in the genomic context can be selected, e.g., ['KOfam',
            'COG14_FUNCTION'], and the 'function_sources' attribute is set to this list.

        function_accessions : dict[str, list[str]]

        function_names : dict[str, list[str]]

        expect_functions : bool

        return_functions : bool

        sum_genes : bool

        average_genes_in_function : bool

        std_genes_in_function : bool

        gene_min_codons : int

        function_min_codons : int

        min_codon_filter : Literal['length', 'remaining', 'both']

        drop_codons : list[str]

        drop_amino_acids : list[str]

        pansequence_min_codons : tuple[int, float]

        pansequence_min_amino_acids : tuple[int, float]

        unreported_codons : list[str]

        unreported_amino_acids : list[str]

        report_min_codons : int

        report_min_amino_acids : int
        """
        self.genome_codon_frequencies = genome_codon_frequencies
        self.genetic_code = genetic_code
        self.auto_adjust = auto_adjust
        self.run = run
        self.progress = progress
        self.verbose = verbose
        self.kwargs = kwargs

        if self.genome_codon_frequencies.genomic_context:
            self.genomic_context = self.genome_codon_frequencies.genomic_context
        else:
            self.genomic_context = None

        A = lambda x, y: kwargs[x] if x in kwargs else y

        self.relative: bool = A('relative', False)
        self.synonymous: bool = A('synonymous', False)
        self.return_amino_acids: bool = A('return_amino_acids', False)
        self.label_amino_acids: bool = A('header_amino_acids', False)
        self.infinity_to_zero: bool = A('infinity_to_zero', False)

        self.gene_caller_ids: list[int] = A('gene_caller_ids', [])
        self.function_sources: list[str]
        if 'function_sources' in kwargs:
            # The argument should have a list value, with an empty list indicating all available
            # sources should be used.
            self.function_sources = kwargs['function_sources']
            if self.genomic_context is None:
                raise ConfigError(
                    "'function_sources' can only be used if available in the genomic context, "
                    "taken from 'genome_codon_frequencies.genomic_context'."
                )
            if not self.function_sources:
                self.function_sources = self.genomic_context.function_sources
        else:
            if self.genomic_context is None:
                self.function_sources = []
            else:
                if self.genomic_context.function_sources is None:
                    self.function_sources = []
                else:
                    self.function_sources = self.genomic_context.function_sources
        self.function_accessions: dict[str, list[str]] = A('function_accessions_dict', {})
        self.function_names: dict[str, list[str]] = A('function_names_dict', {})
        self.expect_functions: bool = A('expect_functions', False)
        self.return_functions: bool = A('return_functions', False)

        self.sum_genes: bool = A('sum', False)
        self.average_genes_in_function: bool = A('average', False)
        self.std_genes_in_function: bool = A('standard_deviation', False)

        self.gene_min_codons: int = A('gene_min_codons', 0)
        self.function_min_codons: int = A('function_min_codons', 0)
        self.min_codon_filter: Literal['length', 'remaining', 'both'] = A(
            'min_codon_filter', 'both'
        )
        self.drop_codons: list[str] = A('exclude_codons', [])
        self.drop_amino_acids: list[str] = A('exclude_amino_acids', [])
        self.pansequence_min_codons: tuple[int, float] = A('pansequence_min_codons', (0, 0.0))
        self.pansequence_min_amino_acids: tuple[int, float] = A(
            'pansequence_min_amino_acids', (0, 0.0)
        )

        self.unreported_codons: list[str] = A('unreported_codons', [])
        self.unreported_amino_acids: list[str] = A('unreported_amino_acids', [])
        self.report_min_codons: int = A('report_min_codons', 0)
        self.report_min_amino_acids: int = A('report_min_amino_acids', 0)

        self.set_up_options()

    def set_up_options(self) -> dict:
        """
        Set up attributes used in methods returning statistical tables.

        Returns
        =======
        dict
            This dictionary is empty if 'self.auto_adjust' is False. If that attribute is True, keys
            are the names of attributes that have been automatically adjusted, and values are the
            original values before adjustment.
        """
        original_values = {}
        # Hold off from setting new attribute values until the end.
        new_values = {}

        if self.synonymous and not self.relative:
            if self.auto_adjust:
                original_values['relative'] = self.relative
                new_values['relative'] = True
                if self.verbose:
                    self.run.info_single(
                        "'relative' was set to True from False due to 'synonymous' being True, "
                        "resulting in calculation of relative synonymous codon frequencies."
                    )
            else:
                raise ConfigError(
                    "'relative' cannot be False if 'synonymous' is True, as 'synonymous' "
                    "calculates relative synonymous codon frequencies."
                )

        if self.synonymous and self.return_amino_acids:
            if self.auto_adjust:
                original_values['return_amino_acids'] = self.return_amino_acids
                new_values['return_amino_acids'] = False
                if self.verbose:
                    self.run.info_single(
                        "'return_amino_acids' was set to False from True due to 'synonymous' being "
                        "True, as 'return_amino_acids' reports amino acid frequencies and "
                        "'synonymous' reports codon frequencies."
                    )
            else:
                raise ConfigError(
                    "'return_amino_acids' and 'synonymous' cannot both be True, as "
                    "'return_amino_acids' reports amino acid frequencies and 'synonymous' reports "
                    "codon frequencies."
                )

        try:
            return_amino_acids = new_values['return_amino_acids']
        except KeyError:
            return_amino_acids = self.return_amino_acids

        if return_amino_acids and self.label_amino_acids:
            if self.auto_adjust:
                original_values['label_amino_acids'] = self.label_amino_acids
                new_values['label_amino_acids'] = False
                if self.verbose:
                    self.run.info_single(
                        "'label_amino_acids' was set to False from True due to "
                        "'return_amino_acids' being True, as 'return_amino_acids' reports columns "
                        "with amino acid headers, and 'label_amino_acids' adds amino acid labels "
                        "to codon headers."
                    )
            else:
                raise ConfigError(
                    "'return_amino_acids' and 'label_amino-acids' cannot both be True, as "
                    "'return_amino_acids' reports columns with amino acid headers, and "
                    "'label_amino_acids' adds amino acid labels to codon headers."
                )

        if return_amino_acids and self.unreported_codons:
            if self.auto_adjust:
                original_values['unreported_codons'] = self.unreported_codons
                new_values['unreported_codons'] = []
                if self.verbose:
                    self.run.info_single(
                        "'unreported_codons' was set to an empty list from a list of codons due to "
                        "'return_amino_acids' being True, as 'return_amino_acids' reports columns "
                        "per amino acid, and 'unreported_codons' drops codon columns."
                    )
            else:
                raise ConfigError(
                    "'unreported_codons' does not work with 'return_amino_acids' being True, as "
                    "'return_amino_acids' reports columns per amino acid, and 'unreported_codons' "
                    "drops codon columns."
                )

        if self.expect_functions and not (self.function_accessions or self.function_names):
            if self.auto_adjust:
                original_values['expect_functions'] = self.expect_functions
                new_values['expect_functions'] = False
                if self.verbose:
                    self.run.info_single(
                        "'expect_functions' was set to False from True due to the absence of "
                        "'function_accessions' or 'function_names'."
                    )
            else:
                raise ConfigError(
                    "'expect_functions' cannot be True in the absence of 'function_accessions' or "
                    "'function_names'."
                )

        if self.sum_genes + self.average_genes_in_function + self.std_genes_in_function > 1:
            raise ConfigError(
                "'sum_genes', 'average_genes_in_function', and 'std_genes_in_function' are "
                "mutually exclusive options for calculating aggregate statistics. Only one can be "
                "True at a time."
            )

        if self.sum_genes and self.function_sources:
            if self.auto_adjust:
                original_values['function_sources'] = self.function_sources
                new_values['function_sources'] = []
                if self.verbose:
                    self.run.info_single(
                        f"'function_sources' was set to an empty list from {self.function_sources} "
                        "due to 'sum_genes' being True. With 'sum_genes', frequencies are summed "
                        "across genes in a genome, irrespective of function groups."
                    )
            else:
                raise ConfigError(
                    "'sum_genes' does not work with 'function_sources', as frequencies are summed "
                    "across genes in a genome, irrespective of function groups."
                )

        if (self.average_genes_in_function or self.std_genes_in_function) and not self.relative:
            raise ConfigError(
                "'average_genes_in_function' and 'std_genes_in_function' must be used with "
                "'relative' or 'synonymous' relative codon frequencies so that genes in a function "
                "are normalized for length."
            )

        if (
            (self.average_genes_in_function or self.std_genes_in_function)
            and not self.return_functions
        ):
            if self.average_genes_in_function:
                msg = 'average_genes_in_function'
            else:
                msg = 'std_genes_in_function'
            if self.auto_adjust:
                original_values['return_functions'] = self.return_functions
                new_values['return_functions'] = True
                if self.verbose:
                    self.run.info_single(
                        f"'return_functions' was set to True from False due to '{msg}' being True, "
                        "as the statistic is calculated among genes of a functional annotation, "
                        "returning a per-function table."
                    )
            else:
                raise ConfigError(
                    f"'{msg}' must be used with 'return_functions', as the statistic is calculated "
                    "among genes of a functional annotation, returning a per-function table."
                )

        if self.min_codon_filter == 'length':
            new_values['filter_gene_length'] = bool(self.gene_min_codons)
            new_values['filter_function_length'] = bool(self.function_min_codons)
            new_values['filter_gene_remaining_codons'] = False
            new_values['filter_function_remaining_codons'] = False
        elif self.min_codon_filter == 'remaining':
            new_values['filter_gene_length'] = False
            new_values['filter_function_length'] = False
            new_values['filter_gene_remaining_codons'] = bool(self.gene_min_codons)
            new_values['filter_function_remaining_codons'] = bool(self.function_min_codons)
        elif self.min_codon_filter == 'both':
            new_values['filter_gene_length'] = bool(self.gene_min_codons)
            new_values['filter_function_length'] = bool(self.function_min_codons)
            new_values['filter_gene_remaining_codons'] = bool(self.gene_min_codons)
            new_values['filter_function_remaining_codons'] = bool(self.function_min_codons)

        for attr_name in ('drop_codons', 'unreported_codons'):
            value = getattr(self, attr_name)
            if value is None:
                continue
            codons = list(self.genetic_code.codon_amino_acid)
            unrecognized_codons = []
            for codon in value:
                if codon not in codons:
                    unrecognized_codons.append(codon)
            if unrecognized_codons:
                msg = ', '.join([f"'{codon}'" for codon in unrecognized_codons])
                raise ConfigError(
                    f"The following codons in '{attr_name}' are not recognized: {msg}"
                )

        if self.synonymous and not self.drop_amino_acids:
            aas = list(self.genetic_code.amino_acid_codons)
            degenerate_aas = list(self.genetic_code.synonymous_amino_acid_codons)
            new_values['drop_amino_acids'] = list(set(aas).difference(set(degenerate_aas)))
            self.run.info_single(
                f"'drop_amino_acids' was set to {new_values['drop_amino_acids']} by default, since "
                "relative synonymous codon frequencies are relevant to degenerate amino acids with "
                "synonymous codons."
            )
            check_aa_attrs = ('unreported_amino_acids', )
        else:
            check_aa_attrs = ('drop_amino_acids', 'unreported_amino_acids')

        for attr_name in check_aa_attrs:
            value = getattr(self, attr_name)
            if value is None:
                continue
            aas = list(self.genetic_code.amino_acid_codons)
            unrecognized_aas = []
            for aa in value:
                if aa not in aas:
                    unrecognized_aas.append(aa)
            if unrecognized_aas:
                msg = ', '.join([f"'{aa}'" for aa in unrecognized_aas])
                raise ConfigError(
                    f"The following amino acids in '{attr_name}' are not recognized: {msg}"
                )

        if (
            type(self.pansequence_min_codons[0]) != int or
            self.pansequence_min_codons[0] < 0 or
            type(self.pansequence_min_codons[1]) != float or
            not (0 <= self.pansequence_min_codons[1] <= 1)
        ):
            raise ConfigError(
                "The value of 'pansequence_min_codons' must be a tuple with two items: the first a "
                "positive int and the second a float between 0.0 and 1.0, inclusive."
            )

        if (
            type(self.pansequence_min_amino_acids[0]) != int or
            self.pansequence_min_amino_acids[0] < 0 or
            type(self.pansequence_min_amino_acids[1]) != float or
            not (0 <= self.pansequence_min_amino_acids[1] <= 1)
        ):
            raise ConfigError(
                "The value of 'pansequence_min_amino_acids' must be a tuple with two items: the "
                "first a positive int and the second a float between 0.0 and 1.0, inclusive."
            )

        if self.function_sources:
            if not self.genomic_context:
                raise ConfigError(
                    "'function_sources' can only be used if available in the 'genomic_context'."
                )
            missing_sources = set(self.function_sources).difference(set(self.genomic_context))
            if missing_sources:
                msg = ', '.join([f'"{source}"' for source in missing_sources])
                raise ConfigError(
                    "The following functional annotation sources in 'function_sources' are not "
                    "available in the 'genomic_context'."
                )

        for attr_name, new_value in new_values.items():
            setattr(self, attr_name, new_value)

        return original_values

    def _reset_attributes(method) -> function:
        """
        Decorator to reset auto-adjusted and temporary attributes after getting statistical output.
        """
        def wrapper(*args, **kwargs):
            output = method(*args, **kwargs)

            self: CodonFrequencyCalculator = args[0]

            # Reset auto-adjusted attributes.
            if self.auto_adjust:
                for attr_name, value in self._original_values.items():
                    setattr(self, attr_name, value)
                delattr(self, '_original_values')

            # Delete temporary attributes.
            for attr_name in (
                '_dropped_codons_panseq_min_codons',
                '_dropped_codons_panseq_min_aas'
            ):
                if hasattr(self, attr_name):
                    delattr(self, attr_name)

            return output
        return wrapper

    @_reset_attributes
    def get_frequencies(self, skip_setup: bool = False) -> pd.DataFrame:
        """
        Get codon frequency statistical table.

        Parameters
        ==========
        skip_setup : bool, False
            Skip setting up attributes, including auto-adjustment.
        """
        if not skip_setup:
            self._original_values = self.set_up_options()

        if self.genome_codon_frequencies.gene_codon_frequency_df is None:
            raise ConfigError(
                "There is no gene codon frequency table available. The value of "
                "'self.genome_codon_frequencies.gene_codon_frequency_df' is None."
            )

        gene_codon_freq_df = self.select_genes()

        if self.function_sources == self.genomic_context.function_sources:
            gene_function_df = self.genomic_context.gene_function_df
        elif self.function_sources:
            gene_function_df = self.genomic_context.gene_function_df.set_index('source').loc[
                self.function_sources
            ].reset_index()
        else:
            gene_function_df = None

        gene_subsetting = self.gene_caller_ids or self.function_accessions or self.function_names

        if self.function_sources and gene_subsetting:
            gene_function_df = gene_function_df.set_index('gene_caller_id')
            gene_function_df = gene_function_df.loc[
                gene_function_df.index.intersection(gene_codon_freq_df.index)
            ].reset_index()

        if gene_subsetting:
            msg1 = pp(len(gene_codon_freq_df))
            msg2 = pp(
                len(self.genomic_context.gene_caller_ids) -
                self.genome_codon_frequencies.noncoding_gene_count
            )
            self.run.info_single(f"{msg1} of {msg2} CDS selected from the genome")
        else:
            self.run.info_single(f"{pp(len(gene_codon_freq_df))} CDS in the genome")

        self._func_cols = ['source', 'accession', 'name']
        self._gene_codon_freq_df, self._gene_func_codon_freq_df = (
            self._filter_gene_functions_codons(gene_codon_freq_df, gene_function_df)
        )

        output_df = None
        for method in (
            self._get_gene_table,
            self._get_sum_table,
            self._get_function_table,
            self._get_average_table,
            self._get_standard_deviation_table
        ):
            output_df = method()
            if output_df is not None:
                break

        # Delete temporary attributes set in this method.
        for attr_name in ('_func_cols', '_gene_codon_freq_df', '_gene_func_codon_freq_df'):
            delattr(self, attr_name)

        if output_df is not None:
            return output_df

    def select_genes(self) -> pd.DataFrame:
        """
        Select genes in the gene frequency table given a list of IDs and/or functions of interest.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Table of gene codon frequencies retaining rows for select genes.
        """
        gene_subsetting = self.gene_caller_ids or self.function_accessions or self.function_names
        if not gene_subsetting:
            return self.genome_codon_frequencies.gene_codon_frequency_df

        if self.function_sources and not self.genomic_context.gene_function_df:
            raise ConfigError(
                "There is no gene function table available. The value of "
                "'self.genomic_context.gene_function_df' is None."
            )

        select_gcids = []

        unrecognized_gcids = []
        for gcid in self.gene_caller_ids:
            if gcid not in self.gene_caller_ids:
                unrecognized_gcids.append(gcid)
                continue
            select_gcids.append(gcid)
        if unrecognized_gcids:
            raise ConfigError(
                "The following requested gene caller IDs were not found in the genome: "
                f"{unrecognized_gcids}"
            )

        # Select GCIDs annotated with select function accessions.
        index_keys = []
        unrecognized_sources = []
        for source, accessions in self.function_accessions.items():
            if source == 'KEGG_BRITE':
                self.run.warning(
                    "KEGG BRITE accessions stored in anvi'o databases are for hierarchies as a "
                    "whole, not categories of the hierarchy. Most hierarchies do not have category "
                    "accessions. So all genes classified in the selected hierarchies are being "
                    "analyzed."
                )
            if source not in self.function_sources:
                unrecognized_sources.append(source)
            for accession in accessions:
                index_keys.append((source, accession))
        if unrecognized_sources:
            raise ConfigError(
                "The following annotation sources in 'function_accessions' were not found to have "
                f"annotated the genome: {unrecognized_sources}"
            )

        gene_function_df = self.genomic_context.gene_function_df.set_index(['source', 'accession'])
        if self.expect_functions:
            index_keys = set(index_keys)
            all_keys = set(gene_function_df.index)
            missing_keys = index_keys.difference(all_keys)

            if missing_keys:
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
                        f"{missing_key_message}"
                    )
        else:
            index_keys = gene_function_df.index.intersection(index_keys)
        select_gcids += gene_function_df.loc[index_keys]['gene_caller_id'].tolist()

        # Select GCIDs annotated with select function names.
        index_keys = []
        unrecognized_sources = []
        for source, names in self.function_names.items():
            if source not in self.function_sources:
                unrecognized_sources.append(source)
            for name in names:
                index_keys.append((source, name))

        if unrecognized_sources:
            raise ConfigError(
                "The following annotation sources in 'function_names' were not found to have "
                f"annotated the genome: {unrecognized_sources}"
            )

        gene_function_df = self.genomic_context.gene_function_df.set_index(['source', 'name'])
        if self.expect_functions:
            index_keys = set(index_keys)
            all_keys = set(gene_function_df.index)
            missing_keys = index_keys.difference(all_keys)

            if missing_keys:
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
                        f"{missing_key_message}"
                    )
        else:
            index_keys = gene_function_df.index.intersection(index_keys)
        select_gcids += gene_function_df.loc[index_keys]['gene_caller_id'].tolist()

        filtered_df = self.genome_codon_frequencies.gene_codon_frequency_df.loc[set(select_gcids)]
        return filtered_df

    def _filter_gene_functions_codons(
        self,
        gene_codon_freq_df: pd.DataFrame,
        gene_function_df: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Filter codons and genes/functions by codon count.

        Parameters
        ==========
        gene_codon_freq_df : pandas.core.Frame.DataFrame
            Table of gene codon frequencies.

        gene_function_df : pandas.core.Frame.DataFrame
            Table of gene functional annotations.

        Returns
        =======
        tuple[pandas.core.Frame.DataFrame, pandas.core.Frame.DataFrame]
            Filtered tables of gene codon frequencies and gene functional annotations, respectively.
        """
        # Filter genes by length.
        gene_codon_freq_df = self._get_frequency_table(
            gene_codon_freq_df,
            min_codons=self.gene_min_codons,
            filter_input_codon_count=self.filter_gene_length
        )

        # Drop certain codon columns from the gene codon frequency table. Filter genes by total
        # codons remaining.
        gene_codon_freq_df = self._get_frequency_table(
            gene_codon_freq_df,
            drop_codons=self.drop_codons,
            drop_amino_acids=self.drop_amino_acids,
            pansequence_min_codons=self.pansequence_min_codons,
            pansequence_min_amino_acids=self.pansequence_min_amino_acids,
            min_codons=self.gene_min_codons,
            filter_output_codon_count=self.filter_gene_remaining_codons
        )

        if self.function_sources and (
            (
                self.drop_codons or
                self.drop_amino_acids or
                self.pansequence_min_codons or
                self.pansequence_min_amino_acids
            ) and
            (self.function_min_codons and self.filter_function_remaining_codons)
        ):
            # Filter functions by total codons remaining. Gene/function pairs with fewer than the
            # required number of codons in the function are removed.
            gene_func_codon_freq_df = gene_function_df.merge(
                gene_codon_freq_df, how='inner', on='gene_caller_id'
            )
            gene_func_codon_freq_df = gene_func_codon_freq_df.set_index(
                self._func_cols + ['gene_caller_id']
            )
            if self.filter_function_remaining_codons:
                gene_func_codon_freq_df = gene_func_codon_freq_df.groupby(self._func_cols).filter(
                    lambda func_df: func_df.sum(axis=1).sum() >= self.function_min_codons
                )
        elif self.function_sources:
            # Filter functions by the total length of their genes. Gene/function pairs with fewer
            # than the required number of codons in the function are removed.
            gene_func_codon_freq_df = gene_function_df.merge(
                gene_codon_freq_df, how='inner', on='gene_caller_id'
            )
            gene_func_codon_freq_df = gene_func_codon_freq_df.set_index(
                self._func_cols + ['gene_caller_id']
            )
            if self.filter_function_length:
                gene_func_codon_freq_df = gene_func_codon_freq_df.groupby(self._func_cols).filter(
                    lambda func_df: func_df.sum(axis=1).sum() >= self.function_min_codons
                )

            gene_func_codon_freq_df = self._get_frequency_table(
                gene_func_codon_freq_df,
                drop_codons=self.drop_codons,
                drop_amino_acids=self.drop_amino_acids,
                pansequence_min_codons=self.pansequence_min_codons,
                pansequence_min_amino_acids=self.pansequence_min_amino_acids
            )
        else:
            gene_func_codon_freq_df = None

        if self.gene_min_codons:
            self.run.info_single(
                f"{pp(len(gene_codon_freq_df))} CDS remaining after codon count filters"
            )

        if hasattr(self, '_dropped_codons_panseq_min_codons'):
            if self.gene_min_codons and self.min_codon_filter != 'remaining':
                min_gene_length_message = '≥' + pp(str(self.gene_min_codons)) + ' codon'
            else:
                min_gene_length_message = ''
            self.run.warning(
                "The following codons were dropped as they did not meet the threshold of "
                f"{pp(self.pansequence_min_amino_acids[0])} codons in "
                f"{self.pansequence_min_amino_acids[1] * 100}% of {min_gene_length_message} "
                f"CDS: {', '.join(sorted(self._dropped_codons_panseq_min_aas))}"
            )

        if hasattr(self, '_dropped_codons_panseq_min_aas'):
            if self.gene_min_codons and self.min_codon_filter != 'remaining':
                min_gene_length_message = '≥' + pp(str(self.gene_min_codons)) + ' codon'
            else:
                min_gene_length_message = ''
            self.run.warning(
                "Codons for the following amino acids were dropped as they did not meet the "
                f"threshold of {pp(self.pansequence_min_amino_acids[0])} codons in "
                f"{self.pansequence_min_amino_acids[1] * 100}% of {min_gene_length_message} "
                f"CDS: {', '.join(sorted(self._dropped_codons_panseq_min_aas))}"
            )

        return gene_codon_freq_df, gene_func_codon_freq_df

    def _get_gene_table(self) -> pd.DataFrame:
        """
        Get a gene-level output table with no consideration of functions.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of gene codon frequencies, which can be raw or normalized, and can be
            consolidated at the level of amino acids.
        """
        if self.function_sources or self.sum_genes:
            return

        get_table = lambda method: method(
            self._gene_codon_freq_df,
            report_min_codons=self.report_min_codons,
            report_min_amino_acids=self.report_min_amino_acids,
            unreported_codons=self.unreported_codons,
            unreported_amino_acids=self.unreported_amino_acids,
            output_amino_acids=self.return_amino_acids,
            label_amino_acids=self.label_amino_acids,
            replace_na=self.infinity_to_zero
        )

        # Absolute frequencies
        if not self.relative:
            return get_table(self._get_frequency_table)

        # Relative frequencies
        if self.relative and not self.synonymous:
            return get_table(self._get_relative_frequency_table)

        # Synonymous relative frequencies
        if not self.return_amino_acids and self.relative and self.synonymous:
            return get_table(self._get_synonymous_codon_relative_frequency_table)

    def _get_sum_table(self) -> pd.DataFrame:
        """
        Get an output table summed across genes in the genome.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of frequencies summed across genes. Frequencies can be raw or normalized,
            and can be consolidated at the level of amino acids.
        """
        if self.function_sources or not self.sum_genes:
            return

        # Absolute codon frequencies
        if not self.relative and not self.return_amino_acids:
            return self._get_frequency_table(
                self._get_summed_frequency_table(self._gene_codon_freq_df),
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids,
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids
            )

        # Absolute amino acid-level codon frequencies
        if not self.relative and self.return_amino_acids:
            return self._get_frequency_table(
                self._get_summed_frequency_table(self._gene_codon_freq_df),
                report_min_amino_acids=self.report_min_amino_acids,
                unreported_amino_acids=self.unreported_amino_acids,
                output_amino_acids=self.return_amino_acids
            )

        # Relative codon frequencies
        if self.relative and not self.synonymous and not self.return_amino_acids:
            return self._get_summed_relative_frequency_table(
                self._gene_codon_freq_df,
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids,
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids
            )

        # Relative amino acid-level codon frequencies
        if self.relative and not self.synonymous and self.return_amino_acids:
            return self._get_summed_relative_frequency_table(
                self._gene_codon_freq_df,
                report_min_amino_acids=self.report_min_amino_acids,
                unreported_amino_acids=self.unreported_amino_acids,
                output_amino_acids=self.return_amino_acids
            )

        # Synonymous relative frequencies
        if self.synonymous and self.relative:
            return self._get_summed_synonymous_codon_relative_frequency_table(
                self._gene_codon_freq_df,
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids,
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids,
                replace_na=self.infinity_to_zero
            )

    def _get_function_table(self) -> pd.DataFrame:
        """
        Get a function-level output table summing genes.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of gene codon frequencies with separate rows per gene-functional
            annotation. Frequencies can be raw or normalized, and can be consolidated at the level
            of amino acids.
        """
        if (
            not self.function_sources or
            self.average_genes_in_function or
            self.std_genes_in_function
        ):
            return

        if self.return_functions:
            input_table = self._gene_func_codon_freq_df.groupby(self._func_cols).sum()
        else:
            input_table = self._gene_func_codon_freq_df.sort_values(self._func_cols)

        get_table = lambda method: method(
            input_table,
            report_min_codons=self.report_min_codons,
            report_min_amino_acids=self.report_min_amino_acids,
            unreported_codons=self.unreported_codons,
            unreported_amino_acids=self.unreported_amino_acids,
            label_amino_acids=self.label_amino_acids
        )

        # Absolute codon frequencies
        if not self.relative and not self.return_amino_acids:
            return get_table(self._get_frequency_table)

        # Relative codon frequencies
        if self.relative and not self.synonymous and not self.return_amino_acids:
            return get_table(self._get_relative_frequency_table)

        get_table = lambda method: method(
            input_table,
            report_min_amino_acids=self.report_min_amino_acids,
            unreported_amino_acids=self.unreported_amino_acids,
            output_amino_acids=self.return_amino_acids
        )

        # Absolute amino acid-level codon frequencies
        if not self.relative and self.return_amino_acids:
            return get_table(self._get_frequency_table)

        # Relative amino acid-level codon frequencies
        if self.relative and not self.synonymous and self.return_amino_acids:
            return get_table(self._get_relative_frequency_table)

        get_table = lambda method: method(
            input_table,
            report_min_codons=self.report_min_codons,
            report_min_amino_acids=self.report_min_amino_acids,
            unreported_codons=self.unreported_codons,
            unreported_amino_acids=self.unreported_amino_acids,
            label_amino_acids=self.label_amino_acids,
            replace_na=self.infinity_to_zero
        )

        # Synonymous relative frequencies
        if self.synonymous and self.relative:
            return get_table(self._get_synonymous_codon_relative_frequency_table)

    def _get_average_table(self) -> pd.DataFrame:
        """
        Get an output table of the average relative frequencies of genes in functions.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of normalized codon frequency averages among genes in functions.
            Frequencies must be relative to avoid aggregating absolute codon counts from genes of
            different lengths. Frequencies can be consolidated at the level of amino acids.
        """
        if not self.average_genes_in_function or not self.function_sources or not self.relative:
            return

        input_table = self._gene_func_codon_freq_df.sort_values(self._func_cols)

        # Relative codon frequencies
        if not self.synonymous and not self.return_amino_acids:
            method = partial(
                self._get_average_relative_frequency_table,
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids,
            )
            return self._get_frequency_table(
                input_table.groupby(self._func_cols).apply(method).droplevel(1),
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids
            )

        # Relative amino acid-level codon frequencies
        if not self.synonymous and self.return_amino_acids:
            method = partial(
                self._get_average_relative_frequency_table,
                report_min_amino_acids=self.report_min_amino_acids,
                output_amino_acids=self.return_amino_acids
            )
            return self._get_frequency_table(
                input_table.groupby(self._func_cols).apply(method).droplevel(1),
                unreported_amino_acids=self.unreported_amino_acids
            )

        # Synonymous relative frequencies
        if self.synonymous:
            method = partial(
                self._get_average_synonymous_codon_relative_frequency_table,
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids
            )
            return self._get_frequency_table(
                input_table.groupby(self._func_cols).apply(method).droplevel(1),
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids,
                replace_na=self.infinity_to_zero
            )

    def _get_standard_deviation_table(self) -> pd.DataFrame:
        """
        Get an output table of the standard deviation of relative frequencies of genes in functions.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of normalized codon frequency standard deviations among genes in
            functions. Frequencies must be relative (including synonymous) to avoid aggregating
            absolute codon counts from genes of different lengths. Frequencies can be consolidated
            at the level of amino acids.
        """
        if not self.std_genes_in_function or not self.function_sources or not self.relative:
            return

        input_table = self._gene_func_codon_freq_df.sort_values(self._func_cols)

        # Relative codon frequencies
        if not self.synonymous and not self.return_amino_acids:
            method = partial(
                self._get_standard_deviation_relative_frequency_table,
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids
            )
            return self._get_frequency_table(
                input_table.groupby(self._func_cols).apply(method).droplevel(1),
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids
            )

        # Relative amino acid-level codon frequencies
        if not self.synonymous and self.return_amino_acids:
            method = partial(
                self._get_standard_deviation_relative_frequency_table,
                report_min_amino_acids=self.report_min_amino_acids,
                output_amino_acids=self.return_amino_acids
            )
            return self._get_frequency_table(
                input_table.groupby(self._func_cols).apply(method).droplevel(1),
                unreported_amino_acids=self.unreported_amino_acids
            )

        # Synonymous relative frequencies
        if not self.return_amino_acids and self.relative and self.synonymous:
            method = partial(
                self._get_standard_deviation_synonymous_codon_relative_frequency_table,
                report_min_codons=self.report_min_codons,
                report_min_amino_acids=self.report_min_amino_acids
            )
            return self._get_frequency_table(
                input_table.groupby(self._func_cols).apply(method).droplevel(1),
                unreported_codons=self.unreported_codons,
                unreported_amino_acids=self.unreported_amino_acids,
                label_amino_acids=self.label_amino_acids,
                replace_na=self.infinity_to_zero
            )

    def _filter_input_codon_count(method) -> function:
        """
        Decorator to discard rows in the input frequency table with fewer than the minimum number of
        codons/amino acids.
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            freq_df: pd.DataFrame = args[1]
            try:
                filter_input_codon_count: bool = kwargs['filter_input_codon_count']
                min_codons: int = kwargs['min_codons']
            except KeyError:
                filter_input_codon_count = False
                min_codons = 0
            if filter_input_codon_count and min_codons:
                freq_df = freq_df[freq_df.sum(axis=1) >= min_codons]
            return method(args[0], freq_df, *args[2: ], **kwargs)
        return wrapper

    def _drop_codon_columns(method) -> function:
        """Decorator to discard codon columns from the input frequency table."""
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = args[1]
            try:
                drop_codons: list[str] = kwargs['drop_codons']
            except KeyError:
                drop_codons: list[str] = []
            if drop_codons:
                codon_freq_df = codon_freq_df.drop(drop_codons, axis=1, errors='ignore')
                if len(codon_freq_df.columns) == 0:
                    codon_freq_df = codon_freq_df.drop(codon_freq_df.index)
            return method(args[0], codon_freq_df, *args[2: ], **kwargs)
        return wrapper

    def _drop_amino_acid_codon_columns(method) -> function:
        """Decorator to discard codon columns by amino acid from the input frequency table."""
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = args[1]
            try:
                drop_aas: list[str] = kwargs['drop_amino_acids']
            except KeyError:
                drop_aas: list[str] = []
            if drop_aas:
                drop_codons: list[str] = []
                self: CodonFrequencyCalculator = args[0]
                aa_codons = self.genetic_code.amino_acid_codons
                for aa in drop_aas:
                    drop_codons += aa_codons[aa]
                codon_freq_df = codon_freq_df.drop(drop_codons, axis=1, errors='ignore')
                if len(codon_freq_df.columns) == 0:
                    codon_freq_df = codon_freq_df.drop(codon_freq_df.index)
            return method(args[0], codon_freq_df, *args[2: ], **kwargs)
        return wrapper

    def _filter_pansequence_codon_count(method) -> function:
        """
        Decorator to drop codons from the input frequency table based on the frequency of the codon
        across rows.
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = args[1]
            try:
                panseq_min_codons: tuple[int, float] = kwargs['pansequence_min_codons']
            except KeyError:
                panseq_min_codons = (0, 0.0)
            min_codon_count = panseq_min_codons[0]
            min_row_fraction = panseq_min_codons[1]
            if min_codon_count > 0 and 0.0 < min_row_fraction <= 1.0:
                row_count = len(codon_freq_df)
                drop_codons: list[str] = []
                self: CodonFrequencyCalculator = args[0]
                for codon in self.genetic_code.codon_amino_acid:
                    if codon not in codon_freq_df.columns:
                        continue
                    filtered_row_count = sum(codon_freq_df[codon] >= min_codon_count)
                    if filtered_row_count / row_count < min_row_fraction:
                        drop_codons += codon
                if drop_codons:
                    codon_freq_df = codon_freq_df.drop(drop_codons, axis=1)
                    self._dropped_codons_panseq_min_codons = drop_codons
                if len(codon_freq_df.columns) == 0:
                    codon_freq_df = codon_freq_df.drop(codon_freq_df.index)
            return method(args[0], codon_freq_df, *args[2: ], **kwargs)
        return wrapper

    def _filter_pansequence_synonymous_codon_count(method) -> function:
        """
        Decorator to drop synonymous codons encoding an amino acid from the input frequency table
        based on the frequency of the amino acid across rows.
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = args[1]
            try:
                panseq_min_aas: tuple[int, float] = kwargs['pansequence_min_amino_acids']
            except KeyError:
                panseq_min_aas = (0, 0.0)
            min_codon_count = panseq_min_aas[0]
            min_row_fraction = panseq_min_aas[1]
            if min_codon_count > 0 and 0.0 < min_row_fraction <= 1.0:
                row_count = len(codon_freq_df)
                drop_codons: list[str] = []
                self: CodonFrequencyCalculator = args[0]
                synonymous_aa_codons = self.genetic_code.synonymous_amino_acid_codons
                for codons in synonymous_aa_codons.values():
                    present_codons = [c for c in codons if c in codon_freq_df.columns]
                    synonymous_freq_df = codon_freq_df[present_codons]
                    if len(synonymous_freq_df.columns) == 0:
                        continue
                    filtered_row_count = sum(synonymous_freq_df.sum(axis=1) >= min_codon_count)
                    if filtered_row_count / row_count < min_row_fraction:
                        drop_codons += codons
                if drop_codons:
                    codon_freq_df = codon_freq_df.drop(drop_codons, axis=1)
                    self._dropped_codons_panseq_min_aas = drop_codons
                if len(codon_freq_df.columns) == 0:
                    codon_freq_df = codon_freq_df.drop(codon_freq_df.index)
            return method(args[0], codon_freq_df, *args[2: ], **kwargs)
        return wrapper

    def _filter_sequence_codon_count(
        self,
        codon_freq_df: pd.DataFrame,
        seq_min_codons: int = 0
    ) -> pd.DataFrame:
        """
        Replace data with NaN for codons from the frequency table based on their count in the row.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item codon frequencies.

        seq_min_codons : int, 0
            Minimum codon count in the row for the codon value to be retained.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Table of item codon frequencies with NaN replacements.
        """
        if seq_min_codons == 0:
            return codon_freq_df
        filtered_df = codon_freq_df[codon_freq_df >= seq_min_codons]
        return filtered_df

    def _filter_sequence_codon_count_wrapper(method) -> function:
        """
        Decorator to replace data with NaN for codons from the frequency table based on their count
        in the row.
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            self: CodonFrequencyCalculator = args[0]
            codon_freq_df: pd.DataFrame = args[1]
            try:
                seq_min_codons: int = kwargs['report_min_codons']
            except KeyError:
                seq_min_codons = 0
            if seq_min_codons > 0:
                codon_freq_df = self._filter_sequence_codon_count(
                    codon_freq_df, seq_min_codons=seq_min_codons
                )
            return method(args[0], codon_freq_df, *args[2: ], **kwargs)
        return wrapper

    def _filter_sequence_synonymous_codon_count(
        self,
        codon_freq_df: pd.DataFrame,
        seq_min_aas: int = 0
    ) -> pd.DataFrame:
        """
        Replace data with NaN for synonymous codons encoding an amino acid from the frequency table
        based on the summed frequency of the synonymous codons in each row.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item codon frequencies.

        seq_min_aas : int, 0
            Minimum codon count for the amino acid in the row for codon values to be retained.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Table of item codon frequencies with NaN replacements.
        """
        if seq_min_aas == 0:
            return codon_freq_df
        mask_df = pd.DataFrame()
        synonymous_aa_codons = self.genetic_code.synonymous_amino_acid_codons
        for codons in synonymous_aa_codons.values():
            present_codons = [c for c in codons if c in codon_freq_df.columns]
            synonymous_freq_df = codon_freq_df[present_codons]
            if len(synonymous_freq_df.columns) == 0:
                continue
            codon_mask_series = synonymous_freq_df.sum(axis=1) >= seq_min_aas
            aa_mask_df = pd.DataFrame(index=codon_mask_series.index)
            for codon in codons:
                aa_mask_df[codon] = codon_mask_series
            mask_df = pd.concat([mask_df, aa_mask_df], axis=1)
        filtered_df = codon_freq_df[mask_df]
        return filtered_df

    def _filter_sequence_synonymous_codon_count_wrapper(method) -> function:
        """
        Decorator to replace data with NaN for synonymous codons encoding an amino acid from the
        frequency table based on the summed frequency of the synonymous codons in each row.
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            self: CodonFrequencyCalculator = args[0]
            codon_freq_df: pd.DataFrame = args[1]
            try:
                seq_min_aas: int = kwargs['report_min_amino_acids']
            except KeyError:
                seq_min_aas = 0
            if seq_min_aas > 0:
                codon_freq_df = self._filter_sequence_synonymous_codon_count(
                    codon_freq_df, seq_min_aas=seq_min_aas
                )
            return method(args[0], codon_freq_df, *args[2: ], **kwargs)
        return wrapper

    def _filter_output_codons(method) -> function:
        """Decorator to drop output columns for select codons."""
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = method(*args, **kwargs)
            try:
                unreported_codons: list[str] = kwargs['unreported_codons']
            except KeyError:
                unreported_codons: list[str] = []
            if unreported_codons:
                codon_freq_df = codon_freq_df.drop(unreported_codons, axis=1, errors='ignore')
            return codon_freq_df
        return wrapper

    def _filter_output_amino_acid_codons(method) -> function:
        """Decorator to drop output columns for synonymous codons encoding select amino acids."""
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = method(*args, **kwargs)
            try:
                unreported_aas: list[str] = kwargs['unreported_amino_acids']
            except KeyError:
                unreported_aas: list[str] = []
            if unreported_aas:
                self: CodonFrequencyCalculator = args[0]
                aa_codons = self.genetic_code.amino_acid_codons
                for aa in unreported_aas:
                    codons = aa_codons[aa]
                    codon_freq_df = codon_freq_df.drop(codons, axis=1, errors='ignore')
            return codon_freq_df
        return wrapper

    def _output_amino_acids(method) -> function:
        """Decorator to output columns of amino acid rather than codon frequencies."""
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df: pd.DataFrame = method(*args, **kwargs)
            try:
                output_aas: bool = kwargs['output_amino_acids']
            except KeyError:
                output_aas = False
            if output_aas:
                aa_df = pd.DataFrame(index=codon_freq_df.index)
                self: CodonFrequencyCalculator = args[0]
                aa_codons = self.genetic_code.amino_acid_codons
                for aa, codons in aa_codons.items():
                    present_codons = [c for c in codons if c in codon_freq_df.columns]
                    select_codon_df = codon_freq_df[present_codons]
                    if len(select_codon_df.columns) == 0:
                        continue
                    aa_df[aa] = select_codon_df[codons].sum(axis=1, skipna=False)
                    return aa_df
            return codon_freq_df
        return wrapper

    def _add_amino_acid_to_header(method) -> function:
        """Decorator to add amino acid to codon column header."""
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df = method(*args, **kwargs)
            try:
                label_amino_acids: bool = kwargs['label_amino_acids']
            except KeyError:
                label_amino_acids = False
            if label_amino_acids:
                self: CodonFrequencyCalculator = args[0]
                codon_aa = self.genetic_code.codon_amino_acid
                try:
                    codon_freq_df.columns = [codon_aa[c] + c for c in codon_freq_df.columns]
                except KeyError:
                    raise ConfigError(
                        "The columns in what should be a table of codon data are not recognized as "
                        "codons. This is the header that is present: "
                        f"{', '.join(list(codon_freq_df.columns))}"
                    )
                codon_freq_df = codon_freq_df[sorted(codon_freq_df.columns)]
            return codon_freq_df
        return wrapper

    def _replace_na(method) -> function:
        """
        Decorator to replace NA with 0.0 in the output frequency table. This should only occur in
        synonymous relative frequency output, in which missing amino acids yield NA (not inf, though
        this is as a consequence of 0/0).
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df = method(*args, **kwargs)
            try:
                replace_na = kwargs['replace_na']
            except KeyError:
                replace_na = False
            if replace_na:
                codon_freq_df = codon_freq_df.fillna(0)
            return codon_freq_df
        return wrapper

    def _filter_output_codon_count(method) -> function:
        """
        Decorator to discard rows in the output frequency table with fewer than the minimum number
        of codons/amino acids.
        """
        def wrapper(*args, **kwargs) -> pd.DataFrame:
            codon_freq_df = method(*args, **kwargs)
            try:
                filter_output_codon_count: bool = kwargs['filter_output_codon_count']
                min_codons: int = kwargs['min_codons']
            except KeyError:
                filter_output_codon_count = False
            if filter_output_codon_count:
                codon_freq_df = codon_freq_df[codon_freq_df.sum(axis=1) >= min_codons]
            return codon_freq_df
        return wrapper

    # The order of decorators should not be changed. (Only '@_output_amino_acids' and
    # '@_add_amino_acid_to_header', which are mutually exclusive operations, are interchangeable.)
    @_filter_input_codon_count
    @_drop_codon_columns
    @_drop_amino_acid_codon_columns
    @_filter_pansequence_codon_count
    @_filter_pansequence_synonymous_codon_count
    @_filter_sequence_synonymous_codon_count_wrapper
    @_filter_output_codons
    @_filter_output_amino_acid_codons
    @_output_amino_acids
    @_add_amino_acid_to_header
    @_filter_output_codon_count
    def _get_frequency_table(self, codon_freq_df: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """
        Get a table of raw codon frequencies.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of item raw codon frequencies.
        """
        return codon_freq_df

    # Commented decorators mean that they can theoretically be used and uncommented but are not
    # because they are not needed in the 'self.get_frequencies' client.
    # @_filter_input_codon_count
    # @_drop_codon_columns
    # @_drop_amino_acid_codon_columns
    # @_filter_pansequence_codon_count
    # @_filter_pansequence_synonymous_codon_count
    @_filter_output_codons
    @_filter_output_amino_acid_codons
    @_output_amino_acids
    @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of relative codon frequencies.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of item relative codon frequencies.
        """
        if 'report_min_codons' in kwargs:
            min_codon_mask_df = self._filter_sequence_codon_count(
                codon_freq_df, kwargs['report_min_codons']
            ).notna()
        else:
            min_codon_mask_df = None

        if 'report_min_amino_acids' in kwargs:
            min_aa_mask_df = self._filter_sequence_synonymous_codon_count(
                codon_freq_df, kwargs['report_min_amino_acids']
            ).notna()
        else:
            min_aa_mask_df = None

        if 'report_min_codons' in kwargs and 'report_min_amino_acids' in kwargs:
            mask_df = min_codon_mask_df | min_aa_mask_df
        elif 'report_min_codons' in kwargs:
            mask_df = min_codon_mask_df
        elif 'report_min_amino_acids' in kwargs:
            mask_df = min_aa_mask_df
        else:
            mask_df = None

        codon_rel_freq_df = codon_freq_df.div(codon_freq_df.sum(axis=1), axis=0)
        drop_index = codon_rel_freq_df[codon_rel_freq_df.isna().all(axis=1)].index
        if mask_df is not None:
            codon_rel_freq_df = codon_rel_freq_df[mask_df]
        # Drop rows with zero frequency.
        codon_rel_freq_df = codon_rel_freq_df.drop(drop_index)
        return codon_rel_freq_df

    # @_filter_input_codon_count
    # @_drop_codon_columns
    # @_drop_amino_acid_codon_columns
    # @_filter_pansequence_codon_count
    # @_filter_pansequence_synonymous_codon_count
    @_filter_output_codons
    @_filter_output_amino_acid_codons
    @_add_amino_acid_to_header
    @_replace_na
    # @_filter_output_codon_count
    def _get_synonymous_codon_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of codon relative frequencies in relation to the set of codons encoding the same
        amino acid (or stop codons). If one or more codons in a synonymous set are absent as columns
        in the input table, synonymous relative frequency is calculated for the remaining codons in
        the set.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of item synonymous relative codon frequencies.
        """
        if 'report_min_codons' in kwargs:
            min_codon_mask_df = self._filter_sequence_codon_count(
                codon_freq_df, kwargs['report_min_codons']
            ).notna()
        else:
            min_codon_mask_df = None

        if 'report_min_amino_acids' in kwargs:
            min_aa_mask_df = self._filter_sequence_synonymous_codon_count(
                codon_freq_df, kwargs['report_min_amino_acids']
            ).notna()
        else:
            min_aa_mask_df = None

        if 'report_min_codons' in kwargs and 'report_min_amino_acids' in kwargs:
            mask_df = min_codon_mask_df | min_aa_mask_df
        elif 'report_min_codons' in kwargs:
            mask_df = min_codon_mask_df
        elif 'report_min_amino_acids' in kwargs:
            mask_df = min_aa_mask_df
        else:
            mask_df = None

        synonymous_freq_df = pd.DataFrame()
        for codons in constants.AA_to_codons.values():
            present_codons = [c for c in codons if c in codon_freq_df.columns]
            aa_codon_freq_df = codon_freq_df[present_codons]
            if len(aa_codon_freq_df.columns) == 0:
                continue
            synonymous_freq_df[present_codons] = aa_codon_freq_df.div(
                aa_codon_freq_df.sum(axis=1), axis=0
            )
        drop_index = synonymous_freq_df[synonymous_freq_df.isna().all(axis=1)].index
        if mask_df is not None:
            synonymous_freq_df = synonymous_freq_df[mask_df]
        synonymous_freq_df = synonymous_freq_df.drop(drop_index)
        return synonymous_freq_df

    # @_filter_input_codon_count
    # @_drop_codon_columns
    # @_drop_amino_acid_codon_columns
    # @_filter_pansequence_codon_count
    # @_filter_pansequence_synonymous_codon_count
    # @_filter_sequence_synonymous_codon_count_wrapper
    # @_filter_output_codons
    # @_filter_output_amino_acid_codons
    # @_output_amino_acids
    # @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_summed_frequency_table(self, codon_freq_df: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """
        Get a table of codon frequencies summed across items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of raw codon frequencies summed across items.
        """
        summed_freq_series = codon_freq_df.sum()
        summed_freq_df = summed_freq_series.to_frame('all').T.rename_axis('gene_caller_id')
        return summed_freq_df

    def _get_summed_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of relative codon frequencies summed across items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of relative codon frequencies summed across items.
        """
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in [
                'report_min_codons',
                'report_min_amino_acids',
                'unreported_codons',
                'unreported_amino_acids',
                'label_amino_acids'
            ]:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value

        summed_freq_df = self._get_summed_frequency_table(codon_freq_df, **first_kwargs)
        summed_rel_freq_df = self._get_rel_freq_table(summed_freq_df, **second_kwargs)
        return summed_rel_freq_df

    def _get_summed_synonymous_codon_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of relative synonymous codon frequencies summed across items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of relative codon frequencies summed across items.
        """
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in [
                'report_min_codons',
                'report_min_amino_acids',
                'unreported_codons',
                'unreported_amino_acids',
                'label_amino_acids',
                'replace_na'
            ]:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value

        summed_codon_freq_df = self._get_summed_frequency_table(codon_freq_df, **first_kwargs)
        summed_synonymous_freq_df = self._get_synonymous_codon_relative_frequency_table(
            summed_codon_freq_df, **second_kwargs
        )
        return summed_synonymous_freq_df

    # @_filter_output_codons
    # @_filter_output_amino_acid_codons
    # @_output_amino_acids
    # @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_average_frequency_table(
        self,
        freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of average frequencies across items.

        Parameters
        ==========
        freq_df : pandas.core.Frame.DataFrame
            Table of item raw frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of frequencies averaged across items.
        """
        avg_freq_series = freq_df.mean()
        avg_freq_df = avg_freq_series.to_frame('all').T.rename_axis('gene_caller_id')
        return avg_freq_df

    def _get_average_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of relative frequencies averaged across items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of relative frequencies averaged across items.
        """
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in [
                'report_min_codons',
                'report_min_amino_acids',
                'unreported_codons',
                'unreported_amino_acids',
                'label_amino_acids'
            ]:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value
        rel_freq_df = self._get_relative_frequency_table(codon_freq_df, **first_kwargs)
        avg_rel_freq_df = self._get_average_frequency_table(rel_freq_df, **second_kwargs)
        return avg_rel_freq_df

    def _get_average_synonymous_codon_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of relative synonymous codon frequencies averaged across items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of relative synonymous codon frequencies averaged across items.
        """
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in [
                'report_min_codons',
                'report_min_amino_acids',
                'unreported_codons',
                'unreported_amino_acids',
                'label_amino_acids',
                'replace_na'
            ]:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value
        synonymous_freq_df = self._get_synonymous_codon_relative_frequency_table(
            codon_freq_df, **first_kwargs
        )
        avg_synonymous_freq_df = self._get_average_frequency_table(
            synonymous_freq_df, **second_kwargs
        )
        return avg_synonymous_freq_df

    # @_filter_output_codons
    # @_filter_output_amino_acid_codons
    # @_output_amino_acids
    # @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_standard_deviation_frequency_table(
        self,
        freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of the standard deviations of frequencies across items.

        Parameters
        ==========
        freq_df : pandas.core.Frame.DataFrame
            Table of item raw frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of the standard deviations of frequencies across items.
        """
        std_freq_series = freq_df.std()
        std_freq_df = std_freq_series.to_frame('all').T.rename_axis('gene_caller_id')
        return std_freq_df

    def _get_standard_deviation_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of the standard deviations of relative frequencies across items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of the standard deviations of relative frequencies across items.
        """
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in [
                'report_min_codons',
                'report_min_amino_acids',
                'unreported_codons',
                'unreported_amino_acids',
                'label_amino_acids'
            ]:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value
        rel_freq_df = self._get_relative_frequency_table(codon_freq_df, **first_kwargs)
        std_rel_freq_df = self._get_standard_deviation_frequency_table(rel_freq_df, **second_kwargs)
        return std_rel_freq_df

    def _get_standard_deviation_synonymous_codon_relative_frequency_table(
        self,
        codon_freq_df: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Get a table of the standard deviations of relative synonymous codon frequencies across
        items.

        Parameters
        ==========
        codon_freq_df : pandas.core.Frame.DataFrame
            Table of item raw codon frequencies.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Filtered table of the standard deviations of synonymous codon frequencies across items.
        """
        first_kwargs = {}
        second_kwargs = {}
        for key, value in kwargs.items():
            if key in [
                'report_min_codons',
                'report_min_amino_acids',
                'unreported_codons',
                'unreported_amino_acids',
                'label_amino_acids',
                'replace_na'
            ]:
                second_kwargs[key] = value
            else:
                first_kwargs[key] = value
        synonymous_freq_df = self._get_synonymous_codon_relative_frequency_table(
            codon_freq_df, **first_kwargs
        )
        std_synonymous_freq_df = self._get_standard_deviation_frequency_table(
            synonymous_freq_df, **second_kwargs
        )
        return std_synonymous_freq_df
