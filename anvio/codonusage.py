#!/usr/bin/env python
# -*- coding: utf-8
"""Module for codon usage analyses at the levels of genes, groups of genes, and genomes."""

from __future__ import annotations

import copy
import inspect
import argparse
import numpy as np
import pandas as pd

from abc import ABC
from copy import deepcopy
from functools import partial
from collections import Counter
from typing import Iterable, Union

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio import __version__ as VERSION
from anvio.dbops import ContigsSuperclass
from anvio.ccollections import GetSplitNamesInBins
from anvio.genomedescriptions import GenomeDescriptions

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"

# run = terminal.Run()
# progress = terminal.Progress()
# run_quiet = terminal.Run(verbose=False)

pp = terminal.pretty_print

class GeneticCode:
    """
    A genetic code relating codons to amino acids.

    Attributes
    ==========
    codon_amino_acid : dict[str, str]
        Mapping of codon key to amino acid value, including stop codons for "STP".

    amino_acid_codons : dict[str, list[str]]
        Mapping of amino acid key to codon values, including stop codons for "STP".

    nonstop_codon_amino_acid : dict[str, str]
        Mapping of codon key to amino acid value, excluding stop codons for "STP".

    nonstop_amino_acid_codons : dict[str, list[str]]
        Mapping of amino acid key to codon values, excluding stop codons for "STP".

    synonymous_codon_amino_acid : dict[str, str]
        Mapping of degenerate codon key to amino acid value, excluding stop codons for "STP".

    synonymous_amino_acid_codons : dict[str, list[str]]
        Mapping of degenerate amino acid key to codon values, excluding stop codons for "STP".

    run : anvio.terminal.Run, anvio.terminal.Run()
        Prints run information to the terminal.
    """
    standard_code: dict[str, str] = {codon: aa for codon, aa in constants.codon_to_AA.items()}

    def __init__(self, code: dict[str, str] = None, run: terminal.Run = terminal.Run()) -> None:
        """
        Parameters
        ==========
        code : dict[str, str]
            Maps codon keys to three-letter amino acid values ("STP" for stop codons). If None, the
            standard genetic code is used.
        """
        self.run = run
        self.set(code)

    def set(self, code: dict[str, str] = None) -> None:
        """
        Set the genetic code.

        Parameters
        ==========
        code : dict[str, str], None
            Maps codon keys to three-letter amino acid values ("STP" for stop codons). If None, the
            standard genetic code is used.
        """
        if code is None:
            self.codon_amino_acid = self.standard_code
        else:
            self.codon_amino_acid = code
            self.check()
            if self.codon_amino_acid != self.standard_code:
                self.run.info_single("Using a nonstandard genetic code for the genome")

        self.amino_acid_codons: dict[str, list[str]] = {}
        for codon, aa in self.codon_amino_acid.items():
            try:
                self.amino_acid_codons[aa].append(codon)
            except KeyError:
                self.amino_acid_codons[aa] = [codon]

        self.nonstop_codon_amino_acid: dict[str, str] = {}
        self.nonstop_amino_acid_codon: dict[str, list[str]] = {}
        for codon, aa in self.codon_amino_acid.items():
            if aa == 'STP':
                continue
            self.nonstop_codon_amino_acid[codon] = aa
            try:
                self.nonstop_amino_acid_codon[aa].append(codon)
            except KeyError:
                self.nonstop_amino_acid_codon[aa] = [codon]

        self.synonymous_codon_amino_acid: dict[str, str] = {}
        self.synonymous_amino_acid_codons: dict[str, list[str]] = {}
        for aa, codons in self.nonstop_amino_acid_codon.items():
            if len(codons) == 1:
                continue
            for codon in codons:
                self.synonymous_codon_amino_acid[codon] = aa
            self.synonymous_amino_acid_codons[aa] = codons.copy()

    def check(self) -> None:
        """
        Check that recognized codons and three-letter amino acid codes are used in the dictionary
        mapping codon to amino acid, raising an exception if untrue.
        """
        unrecognized_codons = []
        unrecognized_amino_acids = []
        for codon, amino_acid in self.codon_amino_acid.items():
            if codon not in constants.codons:
                unrecognized_codons.append(codon)
            if amino_acid not in constants.amino_acids:
                unrecognized_amino_acids.append(amino_acid)

        if unrecognized_codons:
            unrecognized_codon_message = (
                "The following codons in the genetic code are not recognized: "
                f"{', '.join(unrecognized_codons)}."
            )
        else:
            unrecognized_codon_message = ""

        if unrecognized_amino_acids:
            unrecognized_amino_acid_message = (
                "The following amino acids in the genetic code are not recognized: "
                f"{', '.join(unrecognized_amino_acids)}. These should be three-letter codes and "
                "\"STP\" for stop codons."
            )
            if unrecognized_codon_message:
                unrecognized_amino_acid_message = " " + unrecognized_amino_acid_message
        else:
            unrecognized_amino_acid_message = ""

        if unrecognized_codon_message or unrecognized_amino_acid_message:
            raise ConfigError(f"{unrecognized_codon_message}{unrecognized_amino_acid_message}")

    def modify(self, encodings: dict[str, str] = None, encodings_txt: str = None) -> None:
        """
        Modify the genetic code.

        Parameters
        ==========
        encodings : dict[str, str]
            Mapping of codon key to amino acid value. These are the changes made in the set genetic
            code.

        encodings_txt : str
            Path to a tab-delimited file with information like that passed to the 'encodings'
            option. The file should have no header, codons in the first column, and amino acids in
            the second column.
        """
        if encodings is None and encodings_txt is None:
            raise ConfigError("Either 'encodings' or 'encodings_txt' must be provided.")
        if encodings is not None and encodings_txt is not None:
            raise ConfigError("Provide either 'encodings' or 'encodings_txt', not both.")

        new_code = self.codon_amino_acid.copy()

        if encodings is not None:
            new_code.update(encodings)
            self.set(new_code)
            return

        filesnpaths.is_file_tab_delimited(encodings_txt)
        encodings_df = pd.read_csv(encodings_txt, sep='\t', header=None)
        new_code.update(dict(zip(encodings_df.iloc[:, 0], encodings_df.iloc[:, 1])))
        self.set(new_code)

class CodonUsage(ABC):
    """
    Superclass of objects that process codon usage data for a single genome or multiple genomes.

    Attributes
    ==========
    args : argparse.Namespace
        Attributes are populated from arguments.

    run : anvio.terminal.Run, anvio.terminal.Run()
        Prints run information to the terminal.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        Prints transient progress information to the terminal.

    function_sources : Union[list[str], None]
        Defines the gene annotation function sources to be made available from input databases for
        codon usage calculations, including as functions for CUB reference compositions. This can be
        a list of sources, e.g., ['KOfam', 'COG14_FUNCTION']. An empty list for
        'self.args.function_sources' causes all database annotation sources to be included. A value
        of None precludes functions from being analyzed.
    """
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()) -> None:
        self.args = args
        self._A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.run = run
        self.progress = progress

        self.genetic_code = GeneticCode()
        encodings_txt = self._A('encodings_txt')
        if encodings_txt is not None:
            self.genetic_code.modify(encodings_txt=encodings_txt)

        self.function_sources = self._A('function_sources')

class SingleGenomeCodonUsage(CodonUsage):
    """
    Codon usage data for a single genome stored in a contigs database.

    Attributes
    ==========
    contigs_db : str
        File path to a contigs database containing gene sequences.

    profile_db : str, None
        File path to a profile database associated with the contigs database, to be used with
        'collection_name' and 'bin_id'.

    collection_name : str, None
        Name of the bin collection in the profile database, to be used with 'profile_db' and
        'bin_id'.

    bin_id : str, None
        Bin name representing an internal genome in the contigs database, to be used with
        'profile_db' and 'collection_name'.

    contig_sequences_dict : dict[str, str]
        Mapping of contig name to contig nucleotide sequence from contigs database.

    genes_in_contigs_dict : dict[str, list]
        Mapping of gene caller ID to rows of the 'genes_in_contigs' table of the contigs database.

    gene_caller_ids : list[str]
        List of gene caller IDs in the contigs database.

    gene_function_df : pandas.core.frame.DataFrame
        Gene functional annotations for select function sources from the contigs database.
    """
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        super().__init__(args, run, progress)

        self.contigs_db = self._A('contigs_db')

        self.profile_db = self._A('profile_db')
        self.collection_name = self._A('collection_name')
        self.bin_id = self._A('bin_id')

        if self.contigs_db is not None:
            self._load_contigs_db_data()

        self.statistics = CodonStatistics(self, self.run, self.progress)

    def _load_contigs_db_data(self):
        """Load gene data from the contigs database."""
        utils.is_contigs_db(self.contigs_db)

        args = deepcopy(self.args)
        if self.profile_db or self.collection_name or self.bin_id:
            # Initialize the contigs superclass from the splits of the internal genome bin.
            args.split_names_of_interest = GetSplitNamesInBins(args).get_split_names_only()
        contigs_super = ContigsSuperclass(args, r=terminal.Run(verbose=False))
        contigs_super.init_contig_sequences()
        self.contig_sequences_dict = contigs_super.contig_sequences

        self.genes_in_contigs_dict = contigs_super.genes_in_contigs_dict
        self.gene_caller_ids = list(set(self.genes_in_contigs_dict))

        # Organize functional annotations into a table.
        contigs_super.init_functions(requested_sources=self.function_sources)
        if self.function_sources == []:
            self.function_sources = sorted(contigs_super.gene_function_call_sources)
            if not self.function_sources:
                self.run.warning(
                    "The value of `args.function_sources` is an empty list, indicating that all "
                    "function sources in the genome should by loaded. However, the contigs "
                    "database has not been annotated by any sources."
                )
        gene_function_rows = []
        for gcid, annotation_dict in contigs_super.gene_function_calls_dict.items():
            for source, annotation in annotation_dict.items():
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
                        if source == 'KEGG_BRITE':
                            # Include every possible depth of categorization.
                            hierarchy_accession = accession
                            categorization = name
                            split_categories = categorization.split('>>>')
                            for depth in range(1, len(split_categories) + 1):
                                depth_categorization = '>>>'.join(split_categories[: depth])
                                gene_function_rows.append(
                                    [source, hierarchy_accession, depth_categorization, gcid]
                                )
                        else:
                            gene_function_rows.append([source, accession, name, gcid])
                    continue

                # The function accession and name entries do not contain the same number of '!!!'
                # separators. In 'COG20_PATHWAY', there can be multiple accessions corresponding to
                # the same function name.
                if source == 'KEGG_BRITE':
                    # Include every possible depth of categorization.
                    hierarchy_accession = accession
                    categorization = name
                    split_categories = categorization.split('>>>')
                    for depth in range(1, len(split_categories) + 1):
                        depth_categorization = '>>>'.join(split_categories[: depth])
                        gene_function_rows.append(
                            [source, hierarchy_accession, depth_categorization, gcid]
                        )
                else:
                    gene_function_rows.append([source, accession, name, gcid])
        self.gene_function_df = pd.DataFrame(
            gene_function_rows, columns=['source', 'accession', 'name', 'gene_caller_id']
        )

class CodonStatistics:
    """
    Calculates codon statistics from sequences.

    Attributes
    ==========
    genomic_context : SingleGenomeCodonUsage, None
        Contains genomic sequence data loaded from a contigs database.

    function_sources : list[str], []
        Select genes with functional annotation sources in the list, e.g., ['KOfam'] or ['KOfam',
        'COG20_FUNCTION']. If the list contains items, the first four columns of returned
        statistical tables contain, respectively, function annotation sources, function accessions,
        function names, and gene caller IDs. There is a row for each gene/function combination in
        the table, and each row for the same gene contains the same frequency values. When
        'KEGG_BRITE' is in the list, in the absence of the 'function_names' parameter narrowing the
        scope of the inquiry, all possible depths of each BRITE hierarchy in the data are returned,
        e.g., 'Ribosome>>>Ribosomal proteins' and the more general 'Ribosome' would each have rows
        in the output table.

    return_functions : bool, False
        If True, returned statistical tables contain function rather than gene results in each row.
        If False, per-gene results can still be returned while still considering functions with the
        'function_sources' parameter, which facilitates analysis of variability within functions.

    return_amino_acids : bool, False
        If True, returned frequency table columns are decoded amino acids (plus "STP") rather than
        codons. Synonymous codon frequencies are summed to produce the amino acid frequencies.

    gene_caller_ids : list[str], []
        Genes with the given IDs are selected for analysis. This parameter can be used alongside
        'function_accessions' and 'function_names'.

    function_accessions : dict[str, list[str]], None
        Genes annotated with the given function accessions are selected for analysis. The dictionary
        keys are function annotation sources (e.g., 'KOfam', 'COG20_FUNCTION') and values are lists
        of function accessions. This parameter can be used alongside 'gene_caller_ids' and
        'function_names'. Note that 'KEGG_BRITE' accessions are not invididual function accessions
        but overarching hierarchy accessions that include many functions.

    function_names : dict[str, list[str]], None
        Genes annotated with the given function names are selected for analysis. The dictionary keys
        are function annotation sources (e.g., 'KOfam', 'COG20_FUNCTION') and values are lists of
        function names. This parameter can be used alongside 'gene_caller_ids' and
        'function_accessions'. If 'KEGG_BRITE' is used as a source, names (categorizations) should
        be given to an arbitrary level of the hierarchy, like ['Ribosome>>>Ribosomal proteins',
        'Transcription machinery>>>Prokaryotic type>>>Bacterial type>>>RNA polymerase'].

    expect_functions : bool, False
        If True, an error will be raised in getting statistical tables if any given
        'function_accessions' or 'function_names' values are not annotated in the input genome.

    relative : bool, False
        If True, return frequency tables containing relative rather than absolute codon frequencies.

    synonymous : bool, False
        If True, return frequency tables containing relative frequencies among synonymous codons
        decoding each amino acid plus "STP".

    sum_genes : bool, False
        If True, return frequency tables containing a single row per genome of codon frequencies
        summed across genes. Normalization to relative/synonymous values occurs after summing
        absolute frequencies. 'function_sources' are ignored with this option using 'auto_adjust',
        as frequencies are summed irrespective of function groups. 'return_functions' ignores this
        option using 'auto_adjust', by definition summing frequencies among genes in each function.

    average_genes_in_function : bool, False
        If True, return average frequencies of genes in each function. 'relative' or 'synonymous'
        must be True when using this option, as the average only operates on relative frequencies so
        that genes are normalized for length. 'auto_adjust' causes 'return_functions' to be treated
        as True, returning a table of averages per function.

    std_genes_in_function : bool, False
        If True, return standard deviations of gene frequencies in each function. 'relative' or
        'synonymous' must be True when using this option, as the standard deviation only operates on
        relative frequencies so that genes are normalized for length. 'auto_adjust' causes
        'return_funcitons' to be treated as True, returning a table of standard deviations per
        function.

    gene_min_codons : int, 0
        Ignore genes with fewer than this number of codons in analyses.

    function_min_codons : int, 0
        Ignore gene/function pairs when the function has fewer than this number of codons in
        analyses. Genes are filtered by codon count before functions, so the total length of
        remaining genes annotated by the function is considered.

    min_codon_filter : str, 'both'
        This parameter arises from the ambiguity of analytical filters that remove genes and
        functions by number of codons ('gene_min_filters' and 'function_min_codons') in relation to
        the filters that drop codons ('drop_codons', 'drop_amino_acids', 'report_min_amino_acids',
        and `pansequence_min_amino_acids`). Genes (and functions) can be filtered by their full
        length by providing the value, 'length', e.g., genes shorter than 300 codons are ignored.
        They can also be filtered by the number of codons remaining after dropping codons by
        providing the value, 'remaining'. Applying the codon length filter followed by dropping
        codons can result in genes and functions with fewer codons than the length threshold -- thus
        'remaining' is a filter option in addition to 'length' to ensure that total codon
        frequencies in the output always meet the minimum codon threshold. 'both' is needed as an
        option in addition to 'remaining' so dynamic codon filtering by 'report_min_amino_acids' and
        'pansequence_min_amino_acids' operates on genes that passed the first length filter.

    drop_codons : list[str], []
        Remove codons from analyses. CAUTION: Note how this affects the calculation of relative
        frequencies using 'relative' and 'synonymous' attributes, normalizing to the summed
        frequency of retained codons. In synonymous calculations, retained codons in an amino acid
        set are considered, causing affected synonymous codon relative frequencies to deviate from
        the expected meaning.

    unreported_codons : list[str], []
        Remove columns for select codons from returned frequency tables.

    drop_amino_acids : list[str], []
        Remove codons that decode the given amino acids from analyses. Specify three-letter codes,
        e.g., 'Ala' and 'STP' for stop codons. If 'synonymous' is True, the 'drop_amino_acids'
        default becomes the non-degenerate amino acids ('Met' and 'Trp' in the standard genetic
        code) plus 'STP' in setup.

    unreported_amino_acids : list[str], []
        Remove columns for codons encoding select amino acids, or columns for select amino acids
        themselves using 'return_amino_acids', from returned frequency tables.

    report_min_codons : int, 0
        Replace codon values with NaN in returned frequency tables if the codon has a count less
        than this number in a row. For example, if the attribute is 3, and a gene has 2 AAT codons,
        then a row for this gene in the output table will have a missing value in the AAT column.

    report_min_amino_acids : int, 0
        Replace codon values with NaN in returned frequency tables if the set of codons for an amino
        acid has a count less than this number in a row. This option works the same with
        'return_amino_acids', but for amino acid values. For example, if the attribute is 5, and a
        gene has 4 codons encoding Asn, 2 AAT and 2 AAC, then a row for this gene in the output
        table will have missing values in Asn columns.

    pansequence_min_codons : int, (0, 0.0)
        This tuple must have two items. The first is a positive integer representing a minimum
        number of codons -- 'min_codons' -- and the second is a float in the range [0.0, 1.0]
        representing a fraction of genes -- 'min_gene_fraction'. Remove codons from the analysis
        that are less numerous than 'min_codons' in a 'min_gene_fraction' of genes. For example, if
        'min_codons' is 3 and 'min_gene_fraction' is 0.9, then if there are <3 codons for a codon in
        ≥90% of genes, then columns for these codons are dropped. CAUTION: Note how this affects the
        calculation of relative frequencies using 'relative' and 'synonymous' attributes,
        normalizing to the summed frequency of retained codons, not all codons. In synonymous
        calculations, retained codons in an amino acid set are considered, causing affected
        synonymous codon relative frequencies to deviate from the expected meaning.

    pansequence_min_amino_acids : tuple, (0, 0.0)
        This tuple must have two items. The first is a positive integer representing a minimum
        number of codons encoding an amino acid -- 'min_amino_acids' -- and the second is a float in
        the range [0.0, 1.0] representing a fraction of genes -- 'min_gene_fraction'. Remove codons
        from the analysis for amino acids (and 'STP') that are less numerous than 'min_amino_acids'
        in a 'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and
        'min_gene_fraction' is 0.9, then if there are <5 codons for an amino acid/'STP' in ≥90% of
        genes, then the columns for these codons are dropped. CAUTION: Note how this affects the
        calculation of relative frequencies using the 'relative' attribute, normalizing to the
        summed frequency of retained codons, not all codons.

    label_amino_acids : bool, False
        If True, include the amino acid for each codon in the column header of returned frequency
        tables, i.e., 'LysAAA' instead of 'AAA'.

    infinity_to_zero : bool, False
        If True, replace NA (empty) values in returned statistical tables with 0.0. NA occurs if
        'synonymous' is True and all codons for an amino acid are absent in a gene or function,
        resulting in 0/0, reported as NA. Use with caution, for NA and 0.0 mean different things and
        this will skew downstream analyses of synonymous relative frequencies, such as codon usage
        bias.

    auto_adjust : bool, True
        If True, adjust the values of certain attributes given the values of other attributes in
        option setup. For example, if 'synonymous' is True, then 'relative' should also be True. If
        'relative' is False, then with auto adjustment turned on, 'relative' is changed to True
        during setup, and with auto adjustment turned off, an error is thrown.

    available_min_codon_filters : tuple[str]
        Available filters for filtering sequences by minimum number of codons if codons are also
        screened out.
    """
    available_min_codon_filters = ('length', 'remaining', 'both')

    def __init__(
        self,
        genomic_context: SingleGenomeCodonUsage = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        self.run = run
        self.progress = progress

        self.genomic_context = genomic_context

        if self.genomic_context is None:
            self.gene_codon_frequency_df = None
            self.function_sources = []
            self.noncoding_gene_count = 0
            args = argparse.Namespace()
        else:
            self.set_codon_frequencies_from_genome()
            if self.genomic_context.function_sources is None:
                self.function_sources = []
            else:
                self.function_sources = self.genomic_context.function_sources
            args = genomic_context.args

        A = lambda x, y: args.__dict__[x] if x in args.__dict__ else y
        self.gene_caller_ids: list[int] = A('gene_caller_ids', [])
        self.relative: bool = A('relative', False)
        self.synonymous: bool = A('synonymous', False)
        self.return_amino_acids: bool = A('return_amino_acids', False)
        self.label_amino_acids: bool = A('header_amino_acids', False)
        self.infinity_to_zero: bool = A('infinity_to_zero', False)

        self.unreported_codons: list[str] = A('unreported_codons', [])
        self.unreported_amino_acids: list[str] = A('unreported_amino_acids', [])
        self.report_min_codons: int = A('report_min_codons', 0)
        self.report_min_amino_acids: int = A('report_min_amino_acids', 0)

        self.function_accessions: dict[str, list[str]] = A('function_accessions_dict', {})
        self.function_names: dict[str, list[str]] = A('function_names_dict', {})
        self.expect_functions: bool = A('expect_functions', False)
        self.return_functions: bool = A('return_functions', False)

        self.sum_genes: bool = A('sum', False)
        self.average_genes_in_function: bool = A('average', False)
        self.std_genes_in_function: bool = A('standard_deviation', False)

        self.gene_min_codons: int = A('gene_min_codons', 0)
        self.function_min_codons: int = A('function_min_codons', 0)
        self.min_codon_filter: str = A('min_codon_filter', 'both')
        self.drop_codons: list[str] = A('exclude_codons', [])
        self.drop_amino_acids: list[str] = A('exclude_amino_acids', [])
        self.pansequence_min_codons: tuple[int, float] = A('pansequence_min_codons', (0, 0.0))
        self.pansequence_min_amino_acids: tuple[int, float] = A(
            'pansequence_min_amino_acids', (0, 0.0)
        )

        self.auto_adjust_options = A('auto_adjust_options', True)
        self.set_up_options()

    def set_codon_frequencies_from_genome(self) -> None:
        """
        Set the raw codon frequency table from genomic data.
        """
        self.progress.new("Fetching genomic codon frequency data")
        self.progress.update("...")

        gene_codon_freqs = []
        skipped_noncoding_gcids = []
        coding_gcids = []
        genome_gcids = self.genomic_context.gene_caller_ids
        for gcid in genome_gcids:
            # `gene_call` is a dictionary.
            gene_call: dict = self.genomic_context.genes_in_contigs_dict[gcid]

            if gene_call['call_type'] != constants.gene_call_types['CODING']:
                skipped_noncoding_gcids.append(gcid)
                continue

            coding_gcids.append(gcid)

            gene_codon_freqs.append(Counter(utils.get_list_of_codons_for_gene_call(
                gene_call, self.genomic_context.contig_sequences_dict
            )))

        gene_codon_freq_df = pd.DataFrame.from_records(gene_codon_freqs)

        observed_codons = gene_codon_freq_df.columns.tolist()
        for codon in constants.codon_to_AA:
            if codon not in observed_codons:
                gene_codon_freq_df[codon] = 0

        # Drop any column named NaN for unknown codons.
        gene_codon_freq_df = gene_codon_freq_df[constants.codon_to_AA]

        gene_codon_freq_df = gene_codon_freq_df.fillna(0)
        gene_codon_freq_df = gene_codon_freq_df[sorted(gene_codon_freq_df.columns)]
        gene_codon_freq_df.index = coding_gcids
        gene_codon_freq_df.index.name = 'gene_caller_id'
        self.gene_codon_frequency_df = gene_codon_freq_df

        self.progress.end()

        self.noncoding_gene_count = len(skipped_noncoding_gcids)
        if self.noncoding_gene_count:
            self.run.warning(
                f"{pp(self.noncoding_gene_count)} of {pp(len(genome_gcids))} genes were "
                "non-coding and not added to the codon frequency table."
            )

    def set_up_options(self) -> dict:
        """
        Set up attributes used in methods returning statistical tables.

        Returns
        =======
        dict
            This dictionary is always empty is 'self.auto_adjust' is False. If that attribute is
            True, keys are the names of attributes that have been automatically adjusted, and values
            are the original values before adjustment.
        """
        original_values = {}
        # Hold off from setting new attribute values until the end.
        new_values = {}

        if self.synonymous and not self.relative:
            if self.auto_adjust_options:
                original_values['relative'] = self.relative
                new_values['relative'] = True
            else:
                raise ConfigError("'relative' cannot be False if 'synonymous' is True.")

        if self.synonymous and self.return_amino_acids:
            if self.auto_adjust_options:
                original_values['return_amino_acids'] = self.return_amino_acids
                new_values['return_amino_acids'] = False
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
            if self.auto_adjust_options:
                original_values['label_amino_acids'] = self.label_amino_acids
                new_values['label_amino_acids'] = False
            else:
                raise ConfigError(
                    "'return_amino_acids' and 'label_amino-acids' cannot both be True, as "
                    "'return_amino_acids' reports columns with amino acid headers, and "
                    "'label_amino_acids' adds amino acid labels to codon headers."
                )

        if return_amino_acids and self.unreported_codons:
            if self.auto_adjust_options:
                original_values['unreported_codons'] = self.unreported_codons
                new_values['unreported_codons'] = []
            else:
                raise ConfigError(
                    "'unreported_codons' does not work with 'return_amino_acids' being True, as "
                    "'return_amino_acids' reports columns per amino acid, and 'unreported_codons' "
                    "drops codon columns."
                )

        new_values['gene_subsetting'] = bool(
            self.gene_caller_ids or self.function_accessions or self.function_names
        )

        if self.expect_functions and not (self.function_accessions or self.function_names):
            if self.auto_adjust_options:
                original_values['expect_functions'] = self.expect_functions
                new_values['expect_functions'] = False
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
            if self.auto_adjust_options:
                original_values['function_sources'] = self.function_sources
                new_values['function_sources'] = []
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
            if self.auto_adjust_options:
                original_values['return_functions'] = self.return_functions
                new_values['return_functions'] = True
            else:
                raise ConfigError(
                    "'average_genes_in_function' and 'std_genes_in_function' must be used with "
                    "'return_functions', since these aggregate statistics are calculated among "
                    "genes of a functional annotation, returning a per-function table."
                )

        if self.min_codon_filter not in self.available_min_codon_filters:
            message = ', '.join([f"'{filter}'" for filter in self.available_min_codon_filters])
            raise ConfigError(f"The value of 'min_codon_filter' must be one of {message}.")

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
            codons = list(self.genomic_context.genetic_code.codon_amino_acid)
            unrecognized_codons = []
            for codon in value:
                if codon not in codons:
                    unrecognized_codons.append(codon)
            if unrecognized_codons:
                message = ', '.join([f"'{codon}'" for codon in unrecognized_codons])
                raise ConfigError(
                    f"The following codons in '{attr_name}' are not recognized: {message}"
                )

        if self.synonymous and not self.drop_amino_acids:
            aas = list(self.genomic_context.genetic_code.amino_acid_codons)
            degenerate_aas = list(self.genomic_context.genetic_code.synonymous_amino_acid_codons)
            new_values['drop_amino_acids'] = list(set(aas).difference(set(degenerate_aas)))
            check_aa_attrs = ('unreported_amino_acids', )
        else:
            check_aa_attrs = ('drop_amino_acids', 'unreported_amino_acids')

        for attr_name in check_aa_attrs:
            value = getattr(self, attr_name)
            if value is None:
                continue
            aas = list(self.genomic_context.genetic_code.amino_acid_codons)
            unrecognized_aas = []
            for aa in value:
                if aa not in aas:
                    unrecognized_aas.append(aa)
            if unrecognized_aas:
                message = ', '.join([f"'{aa}'" for aa in unrecognized_aas])
                raise ConfigError(
                    f"The following amino acids in '{attr_name}' are not recognized: {message}"
                )

        if (
            type(self.pansequence_min_codons[0]) != int or
            self.pansequence_min_codons[0] < 0 or
            type(self.pansequence_min_codons[1]) != float or
            not (0 < self.pansequence_min_codons[1] <= 1)
        ):
            raise ConfigError(
                "The value of 'pansequence_min_codons' must be a tuple with two items: the first a "
                "positive int and the second a float between 0.0 and 1.0, inclusive."
            )

        if (
            type(self.pansequence_min_amino_acids[0]) != int or
            self.pansequence_min_amino_acids[0] < 0 or
            type(self.pansequence_min_amino_acids[1]) != float or
            not (0 < self.pansequence_min_amino_acids[1] <= 1)
        ):
            raise ConfigError(
                "The value of 'pansequence_min_amino_acids' must be a tuple with two items: the "
                "first a positive int and the second a float between 0.0 and 1.0, inclusive."
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

            self: CodonStatistics = args[0]

            # Reset auto-adjusted attributes.
            if self.auto_adjust_options:
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
        Get a codon frequency table, applying any filters defined in attributes.

        Attributes that filter genes and codons are applied to the raw gene codon frequency table
        in the following order.
            drop genes on codon length ->
            drop gene/function pairs on the total length of remaining genes defining the function ->
            drop codons for select amino acids ->
            dynamically drop codons of rarer amino acids ->
            drop genes on remaining codon frequency ->
            drop gene/function pairs on the codon sum of remaining genes defining the function

        Parameters
        ==========
        skip_setup : bool, False
            Skip the setup of filter attributes, which includes attribute cross-checks.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Codon frequency table.
        """
        if not skip_setup:
            self._original_values = self.set_up_options()

        gene_codon_freq_df = self.select_genes()

        if self.function_sources == self.genomic_context.function_sources:
            gene_function_df = self.genomic_context.gene_function_df
        elif self.function_sources:
            gene_function_df = self.genomic_context.gene_function_df.set_index('source').loc[
                self.function_sources
            ].reset_index()
        else:
            gene_function_df = None

        if self.function_sources and self.gene_subsetting:
            gene_function_df = gene_function_df.set_index('gene_caller_id')
            gene_function_df = gene_function_df.loc[
                gene_function_df.index.intersection(gene_codon_freq_df.index)
            ].reset_index()

        if self.gene_subsetting:
            self.run.info_single(
                f"{pp(len(gene_codon_freq_df))} of "
                f"{pp(len(self.genomic_context.gene_caller_ids) - self.noncoding_gene_count)} CDS "
                "selected from the genome"
            )
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

        # # Get outputs with consideration of functions.
        # output_df = self._get_nonaggregate_function_table(gene_func_codon_freq_df)
        # if output_df is not None:
        #     return output_df

        # # Remove duplicate occurrences of genes before summing or averaging frequencies of all genes
        # # in the function source. A gene can have different annotations, e.g., different KOfam
        # # assignments. In KEGG BRITE, a gene can be in nested categories and different hierarchies.
        # gene_func_codon_freq_df = gene_func_codon_freq_df.reset_index().groupby('source').apply(
        #     lambda source_df: source_df.drop_duplicates(subset='gene_caller_id', ignore_index=True)
        # ).set_index(['source', 'accession', 'name', 'gene_caller_id'])

        # for method in (
        #     self._get_sum_function_table,
        #     self._get_average_codon_function_table,
        #     self._get_average_amino_acid_function_table
        # ):
        #     output_df = method(gene_func_codon_freq_df)
        #     if output_df is not None:
        #         return output_df

        raise ConfigError(
            "This point should not be reached at the end of the 'get_frequencies' method. Please "
            "contact the developers. Here at the ends of the earth you get to hear a top secret "
            "mnemonic for the rare earth elements. \"Scandalous Yiddish language centers praise "
            "Ned's promise of small European garden tubs. Dinosaurs hobble erotically thrumming "
            "yellow lutes.\" (scandium Sc, yttrium Y, lanthanum La, cerium Ce, praseodymium Pr, "
            "neodymium Nd, promethium Pm, samarium Sm, europium Eu, gadolinium Gd, terbium Tb, "
            "dysprosium Dy, holmium Ho, erbium Er, thulium Tm, ytterbium Yb, lutetium Lu) Credit "
            "for this lanthanide series mnemonic goes to Martyn Poliakoff: "
            "https://www.youtube.com/watch?v=Q21clW0s0B8&ab_channel=PeriodicVideos"
        )

    def select_genes(self) -> pd.DataFrame:
        """
        Select genes in the gene frequency table given a list of IDs and/or functions of interest.

        Returns
        =======
        pandas.core.Frame.DataFrame
            Table of gene codon frequencies maintaining rows for select genes.
        """
        if not (self.gene_caller_ids or self.function_accessions or self.function_names):
            return self.gene_codon_frequency_df

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

        filtered_df = self.gene_codon_frequency_df.loc[set(select_gcids)]
        return filtered_df

    def _filter_gene_functions_codons(
        self,
        gene_codon_freq_df: pd.DataFrame,
        gene_function_df: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Filter codons, and genes and functions by codon count.

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

    # def _get_average_codon_gene_table(self, gene_codon_freq_df: pd.DataFrame) -> pd.DataFrame:
    #     """
    #     Get an output table of codon results averaged across genes.

    #     Parameters
    #     ==========
    #     gene_codon_freq_df : pandas.core.Frame.DataFrame
    #         Table of gene codon frequencies.

    #     Returns
    #     =======
    #     pandas.core.Frame.DataFrame
    #         Filtered table of codon frequencies averaged across genes. Frequencies can be raw or
    #         normalized.
    #     """
    #     # Absolute frequencies
    #     if (
    #         not self.return_amino_acids and
    #         not self.relative and
    #         not self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         return self._get_frequency_table(
    #             self._get_average_frequency_table(gene_codon_freq_df),
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids
    #         )

    #     # Relative frequencies
    #     if (
    #         not self.return_amino_acids and
    #         self.relative and
    #         not self.synonymous and
    #         not self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         return self._get_average_relative_frequency_table(
    #             gene_codon_freq_df,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids
    #         )

    #     # Synonymous relative frequencies
    #     if (
    #         not self.return_amino_acids and
    #         self.relative and
    #         self.synonymous and
    #         not self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         return self._get_average_synonymous_codon_relative_frequency_table(
    #             gene_codon_freq_df,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids,
    #             replace_na=self.infinity_to_zero
    #         )

    # def _get_average_amino_acid_function_table(
    #     self,
    #     gene_func_codon_freq_df: pd.DataFrame,
    # ) -> pd.DataFrame:
    #     """
    #     Get an output table of amino acid results averaged across genes in each function source.
    #     Gene amino acid frequencies are first calculated from codon frequencies before averaging
    #     across genes.

    #     Parameters
    #     ==========
    #     gene_func_codon_freq_df : pandas.core.Frame.DataFrame
    #         Table of gene codon frequencies. Functional annotation information is included for each
    #         gene. There is a separate row per gene annotation with the same codon data.

    #     Returns
    #     =======
    #     pandas.core.Frame.DataFrame
    #         Filtered table of amino acid-level frequencies averaged across genes in each function
    #         source. Frequencies can be raw or normalized.
    #     """
    #     # Absolute frequencies
    #     if (
    #         self.return_amino_acids and
    #         not self.relative and
    #         self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(
    #             self._get_average_frequency_table
    #         ).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             output_amino_acids=self.return_amino_acids
    #         )

    #     # Relative frequencies
    #     if (
    #         self.return_amino_acids and
    #         self.relative and
    #         self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         method = partial(
    #             self._get_average_relative_frequency_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids
    #         )
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(method).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             output_amino_acids=self.return_amino_acids
    #         )

    # def _get_average_amino_acid_gene_table(self, gene_codon_freq_df: pd.DataFrame) -> pd.DataFrame:
    #     """
    #     Get an output table of amino acid results averaged across genes. Gene amino acid frequencies
    #     are first calculated from codon frequencies before averaging across genes.

    #     Parameters
    #     ==========
    #     gene_codon_freq_df : pandas.core.Frame.DataFrame
    #         Table of gene codon frequencies.

    #     Returns
    #     =======
    #     pandas.core.Frame.DataFrame
    #         Filtered table of amino acid-level frequencies averaged across genes. Frequencies can be
    #         raw or normalized.
    #     """
    #     # Absolute frequencies
    #     if (
    #         self.return_amino_acids and
    #         not self.relative and
    #         not self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         return self._get_frequency_table(
    #             self._get_average_frequency_table(gene_codon_freq_df),
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             output_amino_acids=self.return_amino_acids
    #         )

    #     # Relative frequencies
    #     if (
    #         self.return_amino_acids and
    #         self.relative and
    #         not self.synonymous and
    #         not self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         return self._get_average_relative_frequency_table(
    #             gene_codon_freq_df,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             output_amino_acids=self.return_amino_acids
    #         )

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

    # def _get_sum_function_table(self, gene_func_codon_freq_df: pd.DataFrame) -> pd.DataFrame:
    #     """
    #     Get an output table summed across genes in each function source.

    #     Parameters
    #     ==========
    #     gene_func_codon_freq_df : pandas.core.Frame.DataFrame
    #         Table of gene codon frequencies. Functional annotation information is included for each
    #         gene. There is a separate row per gene annotation with the same codon data.

    #     Returns
    #     =======
    #     pandas.core.Frame.DataFrame
    #         Filtered table of codon frequencies summed across genes in each function source, raw or
    #         normalized.
    #     """
    #     # Absolute frequencies
    #     if (
    #         not self.relative and
    #         self.function_sources and
    #         self.sum_genes
    #     ):
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(
    #             self._get_summed_frequency_table
    #         ).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             output_amino_acids=self.return_amino_acids,
    #             label_amino_acids=self.label_amino_acids
    #         )

    #     # Relative frequencies
    #     if (
    #         self.relative and
    #         not self.synonymous and
    #         self.function_sources and
    #         self.sum_genes
    #     ):
    #         method = partial(
    #             self._get_summed_relative_frequency_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids
    #         )
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(method).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             output_amino_acids=self.return_amino_acids,
    #             label_amino_acids=self.label_amino_acids
    #         )

    #     # Synonymous relative frequencies
    #     if (
    #         not self.return_amino_acids and
    #         self.relative and
    #         self.synonymous and
    #         self.function_sources and
    #         self.sum_genes
    #     ):
    #         method = partial(
    #             self._get_summed_synonymous_codon_relative_frequency_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids
    #         )
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(method).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids,
    #             replace_na=self.infinity_to_zero
    #         )

    def _get_average_table(self) -> pd.DataFrame:
        """
        Get an output table of the average frequencies of genes in functions.

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
        Get an output table of the standard deviation of genes in functions.

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

    # def _get_average_codon_function_table(
    #     self,
    #     gene_func_codon_freq_df: pd.DataFrame
    # ) -> pd.DataFrame:
    #     """
    #     Get an output table of codon results averaged across genes in each function source.

    #     Parameters
    #     ==========
    #     gene_func_codon_freq_df : pandas.core.Frame.DataFrame
    #         Table of gene codon frequencies. Functional annotation information is included for each
    #         gene. There is a separate row per gene annotation with the same codon data.

    #     Returns
    #     =======
    #     pandas.core.Frame.DataFrame
    #         Filtered table of codon frequencies averaged across genes in each function source.
    #         Frequencies can be raw or normalized.
    #     """
    #     # Absolute frequencies
    #     if (
    #         not self.return_amino_acids and
    #         not self.relative and
    #         self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(
    #             self._get_average_frequency_table
    #         ).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids
    #         )

    #     # Relative frequencies
    #     if (
    #         not self.return_amino_acids and
    #         self.relative and
    #         not self.synonymous and
    #         self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         method = partial(
    #             self._get_average_relative_frequency_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids
    #         )
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(method).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids
    #         )

    #     # Synonymous relative frequencies
    #     if (
    #         not self.return_amino_acids and
    #         self.relative and
    #         self.synonymous and
    #         self.function_sources and
    #         self.average_genes_in_function
    #     ):
    #         method = partial(
    #             self._get_average_synonymous_codon_relative_frequency_table,
    #             sequence_min_codons=self.sequence_min_codons,
    #             sequence_min_amino_acids=self.sequence_min_amino_acids
    #         )
    #         input_table = gene_func_codon_freq_df.groupby('source').apply(method).droplevel(1)
    #         return self._get_frequency_table(
    #             input_table,
    #             unreported_codons=self.unreported_codons,
    #             unreported_amino_acids=self.unreported_amino_acids,
    #             label_amino_acids=self.label_amino_acids,
    #             replace_na=self.infinity_to_zero
    #         )

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
                self: CodonStatistics = args[0]
                aa_codons = self.genomic_context.genetic_code.amino_acid_codons
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
                self: CodonStatistics = args[0]
                for codon in self.genomic_context.genetic_code.codon_amino_acid:
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
                self: CodonStatistics = args[0]
                synonymous_aa_codons = (
                    self.genomic_context.genetic_code.synonymous_amino_acid_codons
                )
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
            self: CodonStatistics = args[0]
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
        synonymous_aa_codons = self.genomic_context.genetic_code.synonymous_amino_acid_codons
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
            self: CodonStatistics = args[0]
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
                self: CodonStatistics = args[0]
                aa_codons = self.genomic_context.genetic_code.amino_acid_codons
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
                self: CodonStatistics = args[0]
                aa_codons = self.genomic_context.genetic_code.amino_acid_codons
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
                self: CodonStatistics = args[0]
                codon_aa = self.genomic_context.genetic_code.codon_amino_acid
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
        codon_freq_df: pd.DataFrame,
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

    def get_codon_usage_bias(self, skip_setup: bool = False) -> dict[str, pd.DataFrame]:
        """
        Get codon usage bias (CUB) tables of genes or functions.
        """
        pass

# class SingleGenomeCodonUsage(CodonUsage):
#     """
#     Processes codon usage data from a single genome.

#     Attributes
#     ==========

#     """
#     def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
#         super().__init__(args, run, progress)

#         self.contigs_db = self._A('contigs_db')

#         self.profile_db_path = self._A('profile_db')
#         self.collection_name = self._A('collection_name')
#         self.bin_id = self._A('bin_id')

#         if self.contigs_db is not None:
#             self._load_contigs_db_data()

#         self.codon_statistics = CodonStatistics(self, self.run, self.progress)

#     def _load_contigs_db_data(self):
#         """Load gene data from the contigs database."""
#         utils.is_contigs_db(self.contigs_db)

#         args = deepcopy(self.args)
#         if self.profile_db or self.collection_name or self.bin_id:
#             # Initialize the contigs superclass from the splits of the internal genome bin.
#             args.split_names_of_interest = GetSplitNamesInBins(args).get_split_names_only()
#         contigs_super = ContigsSuperclass(args, r=terminal.Run(verbose=False))
#         contigs_super.init_contig_sequences()
#         self.contig_sequences_dict = contigs_super.contig_sequences

#         self.genes_in_contigs_dict = contigs_super.genes_in_contigs_dict
#         self.gene_caller_ids = list(set(self.genes_in_contigs_dict))

#         # Organize functional annotations into a table.
#         contigs_super.init_functions(requested_sources=self.function_sources)
#         if self.function_sources == []:
#             self.function_sources = sorted(contigs_super.gene_function_call_sources)
#             if not self.function_sources:
#                 self.run.warning(
#                     "The value of `args.function_sources` is an empty list, indicating that all "
#                     "function sources in the genome should by loaded. However, the contigs "
#                     "database has not been annotated by any sources."
#                 )
#         gene_function_rows = []
#         for gcid, annotation_dict in contigs_super.gene_function_calls_dict.items():
#             for source, annotation in annotation_dict.items():
#                 if annotation is None:
#                     continue

#                 accession = annotation[0]
#                 name = annotation[1]
#                 accessions = accession.split('!!!')
#                 names = name.split('!!!')
#                 if len(accessions) == len(names):
#                     # The function accession and name entries contain the same number of '!!!'
#                     # separators.
#                     for accession, name in zip(accessions, names):
#                         if source == 'KEGG_BRITE':
#                             # Include every possible depth of categorization.
#                             hierarchy_accession = accession
#                             categorization = name
#                             split_categories = categorization.split('>>>')
#                             for depth in range(1, len(split_categories) + 1):
#                                 depth_categorization = '>>>'.join(split_categories[: depth])
#                                 gene_function_rows.append(
#                                     [source, hierarchy_accession, depth_categorization, gcid]
#                                 )
#                         else:
#                             gene_function_rows.append([source, accession, name, gcid])
#                     continue

#                 # The function accession and name entries do not contain the same number of '!!!'
#                 # separators. In 'COG20_PATHWAY', there can be multiple accessions corresponding to
#                 # the same function name.
#                 if source == 'KEGG_BRITE':
#                     # Include every possible depth of categorization.
#                     hierarchy_accession = accession
#                     categorization = name
#                     split_categories = categorization.split('>>>')
#                     for depth in range(1, len(split_categories) + 1):
#                         depth_categorization = '>>>'.join(split_categories[: depth])
#                         gene_function_rows.append(
#                             [source, hierarchy_accession, depth_categorization, gcid]
#                         )
#                 else:
#                     gene_function_rows.append([source, accession, name, gcid])
#         self.gene_function_df = pd.DataFrame(
#             gene_function_rows, columns=['source', 'accession', 'name', 'gene_caller_id']
#         )

class CUBDefaults:
    """

    Attributes
    ==========
    metrics : tuple[str], ('cai', 'delta')
        Available CUB metrics. CAI is the Codon Adaptation Index of Sharp and Li
        (1987). Delta is from Ran and Higgs (2012, Eq. 6).

    reference_dependent_metrics : tuple[str], ('cai', 'delta')
        Available CUB metrics that rely upon comparison to a set of reference genes.

    reference_function_source : str, 'KEGG_BRITE'

    reference_function_accessions : list[str], []

    reference_function_names : list[str], ['Ribosome>>>Ribosomal proteins']
        Genes annotated by ribosomal protein KOs are used as the CUB reference.

    query_min_analyzed_codons : int, 100
        Minimum gene/function codon count to improve CUB statistical significance.

    reference_exclude_amino_acid_count : int, 10
        Minimum reference amino acid codon count to improve CUB statistical significance.

    reference_min_analyzed_codons : int, 100
        Minimum reference codon count to improve CUB statistical significance.

    genetic_code : GeneticCode, GeneticCode()
        Certain CUB defaults depend on the genetic code.

    ignored_codons : list[str]
        Codons ignored by default because they encode stop or are the only codon for an amino acid.
    """
    metrics: tuple[str] = ('cai', 'delta')
    reference_dependent_metrics: tuple[str] = ('cai', 'delta')

    reference_function_source: str = 'KEGG_BRITE'
    reference_function_accessions: list[str] = []
    reference_function_names: list[str] = ['Ribosome>>>Ribosomal proteins']

    # Minimum codon counts improve the statistical power of CUB.
    query_min_analyzed_codons: int = 100
    reference_exclude_amino_acid_count: int = 10
    reference_min_analyzed_codons: int = 100

    def __init__(self, genetic_code: GeneticCode = GeneticCode()) -> None:
        """
        Parameters
        ==========
        genetic_code : GeneticCode, GeneticCode()
            Certain CUB defaults depend on the genetic code, the standard code if not given.
        """
        self.genetic_code = genetic_code
        # Anvi'o records codons with an "N" nucleotide as None.
        self.ignored_codons = [None]
        for aa, codons in genetic_code.amino_acid_codons:
            if aa == 'STP' or len(codons) == 1:
                self.ignored_codons += codons



# class SingleGenomeCodonUsage1(CodonUsage):
#     """
#     Processes codon usage data for a single genome.

#     Attributes
#     ==========
#     contigs_db : str
#         File path to a contigs database containing gene sequences.

#     profile_db : str, None
#         File path to a profile database associated with the contigs database, to be used with
#         'collection_name' and 'bin_id'.

#     collection_name : str, None
#         Name of the bin collection in the profile database, to be used with 'profile_db' and
#         'bin_id'.

#     bin_id : str, None
#         Bin name representing an internal genome in the contigs database, to be used with
#         'profile_db' and 'collection_name'.
#     """

#     def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
#         super().__init__(args, run, progress)

#         self.contigs_db = self._A('contigs_db')

#         self.profile_db_path = self._A('profile_db')
#         self.collection_name = self._A('collection_name')
#         self.bin_id = self._A('bin_id')

#         self._load_contigs_db_data()
#         self.codon_statistics = CodonStatistics(self, self.run, self.progress)

#     # def _set_genetic_code(self):
#     #     """
#     #     Set decoding properties of the genome as attributes.

#     #     The dict, `args.codon_to_amino_acid`, should have keys that are codons and values that are
#     #     three-letter amino acid codes ("STP" for stop codons). If `args.codon_to_amino_acid` is
#     #     None, the standard genetic code is used.
#     #     """
#     #     if self.args.codon_to_amino_acid is None:
#     #         self.codon_amino_acid_dict = standard_code
#     #     else:
#     #         check_genetic_code(self.args.codon_to_amino_acid)
#     #         self.codon_amino_acid_dict = self.args.codon_to_amino_acid
#     #         if self.codon_amino_acid_dict != standard_code:
#     #             self.run.info_single("Using a nonstandard genetic code for the genome")

#     #     self.amino_acid_codons_dict = {}
#     #     for codon, amino_acid in self.codon_amino_acid_dict.items():
#     #         try:
#     #             self.amino_acid_codons_dict[amino_acid].append(codon)
#     #         except KeyError:
#     #             self.amino_acid_codons_dict[amino_acid] = [codon]

#     #     self.nonstop_amino_acid_codons_dict = {
#     #         amino_acid: codons for amino_acid, codons in
#     #         self.amino_acid_codons_dict.items() if amino_acid != 'STP'}

#     #     # Ignore codons with "N" nucleotide which anvi'o records as None.
#     #     self.ignored_cub_codons = [None]
#     #     for amino_acid in ignored_cub_amino_acids:
#     #         self.ignored_cub_codons += self.amino_acid_codons_dict[amino_acid]

#     #     self.synonymous_nonstop_amino_acid_codons_dict = {
#     #         amino_acid: codons for amino_acid, codons in self.nonstop_amino_acid_codons_dict.items()
#     #         if amino_acid not in ignored_cub_amino_acids}


#     def _load_contigs_db_data(self):
#         """Load gene data from the contigs database."""
#         utils.is_contigs_db(self.contigs_db)

#         args = deepcopy(self.args)
#         if self.profile_db or self.collection_name or self.bin_id:
#             # Initialize the contigs superclass from the splits of the internal genome bin.
#             args.split_names_of_interest = GetSplitNamesInBins(args).get_split_names_only()
#         contigs_super = ContigsSuperclass(args, r=terminal.Run(verbose=False))
#         contigs_super.init_contig_sequences()
#         self.contig_sequences_dict = contigs_super.contig_sequences

#         self.genes_in_contigs_dict = contigs_super.genes_in_contigs_dict
#         self.gene_caller_ids = list(set(self.genes_in_contigs_dict))

#         # Organize functional annotations into a table.
#         contigs_super.init_functions(requested_sources=self.function_sources)
#         if self.function_sources == []:
#             self.function_sources = sorted(contigs_super.gene_function_call_sources)
#             if not self.function_sources:
#                 self.run.warning(
#                     "The value of `args.function_sources` is an empty list, indicating that all "
#                     "function sources in the genome should by loaded. However, the contigs "
#                     "database has not been annotated by any sources."
#                 )
#         gene_function_rows = []
#         for gcid, annotation_dict in contigs_super.gene_function_calls_dict.items():
#             for source, annotation in annotation_dict.items():
#                 if annotation is None:
#                     continue

#                 accession = annotation[0]
#                 name = annotation[1]
#                 accessions = accession.split('!!!')
#                 names = name.split('!!!')
#                 if len(accessions) == len(names):
#                     # The function accession and name entries contain the same number of '!!!'
#                     # separators.
#                     for accession, name in zip(accessions, names):
#                         if source == 'KEGG_BRITE':
#                             # Include every possible depth of categorization.
#                             hierarchy_accession = accession
#                             categorization = name
#                             split_categories = categorization.split('>>>')
#                             for depth in range(1, len(split_categories) + 1):
#                                 depth_categorization = '>>>'.join(split_categories[: depth])
#                                 gene_function_rows.append(
#                                     [source, hierarchy_accession, depth_categorization, gcid]
#                                 )
#                         else:
#                             gene_function_rows.append([source, accession, name, gcid])
#                     continue

#                 # The function accession and name entries do not contain the same number of '!!!'
#                 # separators. In 'COG20_PATHWAY', there can be multiple accessions corresponding to
#                 # the same function name.
#                 if source == 'KEGG_BRITE':
#                     # Include every possible depth of categorization.
#                     hierarchy_accession = accession
#                     categorization = name
#                     split_categories = categorization.split('>>>')
#                     for depth in range(1, len(split_categories) + 1):
#                         depth_categorization = '>>>'.join(split_categories[: depth])
#                         gene_function_rows.append(
#                             [source, hierarchy_accession, depth_categorization, gcid]
#                         )
#                 else:
#                     gene_function_rows.append([source, accession, name, gcid])
#         self.gene_function_df = pd.DataFrame(
#             gene_function_rows, columns=['source', 'accession', 'name', 'gene_caller_id']
#         )


#     # def _make_gene_codon_frequency_table(self):
#     #     """Generate the per-gene codon frequency DataFrame as `self.gene_codon_frequency_df`."""

#     #     self.progress.new("Fetching codon frequency data")
#     #     self.progress.update("...")

#     #     gene_codon_frequencies = []
#     #     skipped_noncoding_gene_caller_ids = []
#     #     coding_gene_caller_ids = []
#     #     for gene_caller_id in self.gene_caller_ids:
#     #         # `gene_call` is a dictionary.
#     #         gene_call = self.genes_in_contigs_dict[gene_caller_id]

#     #         if gene_call['call_type'] != constants.gene_call_types['CODING']:
#     #             skipped_noncoding_gene_caller_ids.append(gene_caller_id)
#     #             continue

#     #         coding_gene_caller_ids.append(gene_caller_id)

#     #         gene_codon_frequencies.append(Counter(utils.get_list_of_codons_for_gene_call(
#     #             gene_call, self.contig_sequences_dict)))

#     #     gene_codon_frequency_df = pd.DataFrame.from_records(gene_codon_frequencies)

#     #     observed_codons = gene_codon_frequency_df.columns.tolist()
#     #     for codon in constants.codon_to_AA:
#     #         if codon not in observed_codons:
#     #             gene_codon_frequency_df[codon] = 0

#     #     # Drop any column named NaN for unknown codons.
#     #     gene_codon_frequency_df = gene_codon_frequency_df[constants.codon_to_AA]

#     #     gene_codon_frequency_df = gene_codon_frequency_df.fillna(0)
#     #     gene_codon_frequency_df = gene_codon_frequency_df[sorted(gene_codon_frequency_df.columns)]
#     #     gene_codon_frequency_df.index = coding_gene_caller_ids
#     #     gene_codon_frequency_df.index.name = 'gene_caller_id'
#     #     self.gene_codon_frequency_df = gene_codon_frequency_df

#     #     self.progress.end()

#     #     self.noncoding_gene_count = len(skipped_noncoding_gene_caller_ids)
#     #     if self.noncoding_gene_count:
#     #         self.run.warning(
#     #             f"{pp(self.noncoding_gene_count)} of {pp(len(self.gene_caller_ids))} genes were "
#     #             "non-coding and not added to the codon frequency table.")


#     def get_frequencies(self,
#                         from_function_sources=False,
#                         return_functions=False,
#                         return_amino_acids=False,
#                         gene_caller_ids=None,
#                         function_accessions=None,
#                         function_names=None,
#                         expect_functions=False,
#                         relative=False,
#                         synonymous=False,
#                         sum_genes=False,
#                         average_genes=False,
#                         gene_min_codons=0,
#                         function_min_codons=0,
#                         min_codon_filter='both',
#                         drop_codons: Iterable[str] = None,
#                         unreported_codons=None,
#                         drop_amino_acids=None,
#                         unreported_amino_acids=None,
#                         sequence_min_amino_acids=0,
#                         pansequence_min_amino_acids=(0, 1.0),
#                         label_amino_acids=False,
#                         infinity_to_zero=False):
#         """
#         Get absolute (default) or relative codon or amino acid frequencies from genes or functions.

#         Relative codon frequencies can be normalized per-amino acid to all synonymous codons.

#         Parameters
#         ----------
#         from_function_sources : bool, str, or iterable of str, optional
#             Select genes with functional annotations. With this argument, the first four columns of
#             the returned table contain, respectively, function annotation sources, function
#             accessions, function names, and gene caller IDs. There is a row for each gene/function
#             combination in the table, and each row for the same gene contains the same frequency
#             values. When this argument is True, use all available functional annotation sources in
#             the SingleGenomeCodonUsage object. When this argument is a string, select the source
#             given by the string, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. When this argument is an
#             iterable, select a subset of sources given by its strings, e.g., ['KOfam',
#             'COG20_FUNCTION']. When 'KEGG_BRITE' is in the argument, in the absence of
#             `function_names` narrowing the scope of the inquiry, all possible depths of each BRITE
#             hierarchy in the data are returned, e.g., 'Ribosome>>>Ribosomal proteins' and the more
#             general 'Ribosome' would each have rows in the output table. By default None.

#         return_functions : bool, optional
#             If True (default False), output frequency tables contain function rather than gene
#             results in each row. Returning per-gene results when also considering functions by using
#             `--from-function-sources` facilitates analysis of per-gene variability within functions.

#         return_amino_acids : bool, optional
#             If True (default False), output frequency table columns are decoded amino acids (plus
#             STP) rather than codons. Synonymous codon frequencies are summed to produce the amino
#             acid frequencies.

#         gene_caller_ids : iterable, optional
#             Genes with the given IDs are selected for analysis. This parameter can be used alongside
#             `function_accessions` and `function_names`. By default None.

#         function_accessions : dict, optional
#             Genes annotated with the given function accessions are selected for analysis. The
#             argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
#             parameter can be used alongside `gene_caller_ids` and `function_names`.  Note that
#             'KEGG_BRITE' does not use individual function accessions but overarching hierarchy
#             accessions that include multiple functions. By default None.

#         function_names : dict, optional
#             Genes annotated with the given function names are selected for analysis. The argument
#             must be a dict keyed by function annotation source (e.g., 'KOfam', 'COG20_FUNCTION',
#             'Pfam') and with values being lists of function names. Unlike `function_accessions`,
#             'KEGG_BRITE' may be used as a source, with names (categorizations) given to an arbitrary
#             level of the hierarchy, such as ['Ribosome>>>Ribosomal proteins', 'Transcription
#             machinery>>>Prokaryotic type>>>Bacterial type>>>RNA polymerase']. This parameter can be
#             used alongside `gene_caller_ids` and `function_accessions`. By default None.

#         expect_functions : bool, optional
#             If True (default False), an error will be raised if any given `function_accessions` or
#             `function_names` are not annotated in the input genome.

#         relative : bool, optional
#             If True (default False), return relative rather than absolute codon frequencies.

#         synonymous : bool, optional
#             If True (default False), return codon relative frequencies among synonymous codons
#             decoding each amino acid (plus stop).

#         sum_genes : bool, optional
#             If True (default False), sum codon frequencies of genes, returning a one-row DataFrame
#             of the summed frequencies. If `from_function_sources` is used, then genes are limited to
#             those with the specified functional annotations. Relative/synonymous frequencies are
#             calculated after summing absolute frequencies.

#         average_genes : bool, optional
#             If True (default False), average codon frequencies of genes, returning a one-row
#             DataFrame of the averaged frequencies. If `from_function_sources` is used, then genes
#             are limited to those with the specified functional annotations. Averaging occurs after
#             calculation of relative/synonymous frequencies.

#         Additional Parameters
#         ---------------------
#         These parameters filter genes and codons. Here is the order of all possible filters:
#         gene codon frequency table ->
#             drop genes on codon length ->
#             drop gene/function pairs on the total length of remaining genes defining the function ->
#             drop codons of defined amino acids ->
#             dynamically drop codons of rarer amino acids ->
#             drop genes on remaining codon frequency ->
#             drop gene/function pairs on the codon sum of remaining genes defining the function ->
#         filtered gene codon frequency table

#         gene_min_codons : int, optional
#             Ignore genes with fewer than this number of codons. By default 0.

#         function_min_codons : int, optional
#             Ignore gene/function pairs when the function has fewer than this number of codons. Genes
#             are filtered by codon count before functions, so the total length of remaining genes
#             annotated by the function is considered. By default 0.

#         min_codon_filter : {"length", "remaining", "both"}, optional
#             This argument arises from the ambiguity of filters that remove genes and functions by
#             number of codons (`--gene-min-codons` and `--function-min-codons`) in relation to the
#             filters that drop codons (`--exclude/include-amino-acids` and
#             `--exclude-amino-acid-count/fraction`). Genes (and functions) can be filtered by their
#             full length, e.g., genes shorter than 300 codons are ignored. They can also be
#             filtered by the number of codons remaining after dropping codons. The codon length
#             filter followed by dropping codons can result in genes and functions with fewer codons
#             than the original codon threshold -- thus the option of both 'length' and 'remaining'
#             filters to ensure that total codon frequencies in the output always meet the minimum
#             codon threshold. 'both' is needed as an option in addition to 'remaining' so dynamic
#             codon filtering by `--exclude-amino-acid-count/fraction` operates on genes that passed
#             the first length filter.

#         drop_codons : iterable, None

#         drop_amino_acids : iterable, optional
#             Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala, and
#             STP for stop codons). If `synonymous` is True, the `drop_amino_acids` default rather
#             than being None is STP plus amino acids encoded by a single codon (Met, Trp).

#         sequence_min_amino_acids : int, optional
#             Remove from the output rows codons for amino acids (and STP), or amino acids themselves,
#             that are less numerous than `sequence_min_amino_acids`. For example, if the argument is
#             5, and a gene has 4 codons encoding Asn, 2 AAT and 2 AAC, then a row for this gene in
#             the output table will have missing values in Asn columns. By default 0.

#         pansequence_min_amino_acids : tuple, optional
#             This tuple must have two items. The first is an int >0 representing a minimum number of
#             codons encoding an amino acid -- 'min_amino_acids' -- and the second is a float in the
#             range (0, 1) representing a fraction of genes -- 'min_gene_fraction'. Remove codons for
#             amino acids (and STP) that are less numerous than 'min_amino_acids' in a
#             'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and
#             'min_gene_fraction' is 0.9, then if there are fewer than 5 codons for amino acid/STP in
#             ≥90% of genes, then the columns for these codons are dropped. By default (0, 1.0).

#         label_amino_acids : bool, optional
#             If True (default False), include the amino acid for each codon in the column header of
#             the output, i.e., LysAAA instead of AAA.

#         infinity_to_zero : bool, optional
#             If True (default False), replace NA (empty) values in output with 0.0. NA occurs if
#             `synonymous` is True and all codons for an amino acid are absent in a gene or function,
#             resulting in 0/0, reported as NA. Use with caution, for NA and 0.0 mean different things
#             and this will skew downstream analyses of synonymous relative frequencies, such as codon
#             usage bias.

#         Returns
#         -------
#         pandas.core.frame.DataFrame
#             Frequency table of gene x codon or amino acid. If functions are not considered, then the
#             Index of the returned DataFrame contains gene caller IDs. If functions are considered,
#             then each row represents a gene/function pair, and the same gene frequencies can be
#             found in multiple rows with different function pairs; there are additional MultiIndex
#             columns for source, accession, and name of function. If frequencies are summed or
#             averaged and functions are not considered, the returned single-row DataFrame has an
#             Index with one entry, 'all'. If frequencies are summed or averaged and functions are
#             considered, the returned DataFrame has a row for each function annotation source.

#         Examples
#         --------
#         Return the amino acid frequencies of each gene annotated by KEGG KOfams.
#         >>> self.get_frequencies(from_function_sources='KOfam', return_amino_acids=True)

#         Return the summed amino acid frequencies of genes annotated by KEGG KOfams.
#         >>> self.get_frequencies(
#             from_function_sources='KOfam', return_amino_acids=True, sum_genes=True)

#         Return the codon relative frequencies of each gene.
#         >>> self.get_frequencies(relative=True)

#         Return the relative frequencies of genes ≥300 codons in length.
#         >>> self.get_frequencies(relative=True, gene_min_codons=300)

#         Return the average relative frequencies of genes ≥300 codons in length.
#         >>> self.get_frequencies(relative=True, average_all=True, gene_min_codons=300)

#         Return the synonymous (per-amino acid) relative frequencies of genes. Remove genes <300
#         codons, then remove amino acids with <5 codons in ≥90% of remaining genes, then remove genes
#         with <300 codons remaining.
#         >>> self.get_frequencies(
#             synonymous=True, gene_min_codons=300, pansequence_min_amino_acids=(5, 0.9))

#         Return the synonymous (per-amino acid) relative frequencies of KEGG KOfam and BRITE
#         functions. This command removes genes <300 codons, then removes amino acids with <5 codons
#         in ≥90% of remaining genes, then removes genes with <300 codons remaining, then sums gene
#         frequencies in functions, then converts function codon frequencies to missing values for
#         amino acids with <5 codons in the function, then calculates synonymous relative frequencies
#         in functions.
#         >>> self.get_frequencies(
#             from_function_sources=['KOfam', 'KEGG_BRITE'],
#             return_functions=True,
#             synonymous=True,
#             gene_min_codons=300,
#             sequence_min_amino_acids=5,
#             pansequence_min_amino_acids=(5, 0.9))
#         """
#         # CHECK ARGUMENTS AND SET UP PROCEDURE
#         ######################################
#         function_sources = self._establish_function_sources(from_function_sources)
#         if function_sources == self.function_sources:
#             gene_function_df = self.gene_function_df
#         elif function_sources:
#             # Subset the gene function table to the requested sources.
#             gene_function_df = \
#                 self.gene_function_df.set_index('source').loc[function_sources].reset_index()

#         # Check compatability of `return_amino_acids` with other arguments.
#         if return_amino_acids and synonymous:
#             raise ConfigError("The argument `synonymous` should only be True when "
#                               "`return_amino_acids` is False, as `synonymous` returns "
#                               "synonymous codon relative frequencies.")
#         if return_amino_acids and label_amino_acids:
#             # Don't bother raising an exception.
#             label_amino_acids = False

#         if gene_caller_ids is None:
#             gene_caller_ids = []
#         if function_accessions is None:
#             function_accessions = {}
#         if function_names is None:
#             function_names = {}
#         gene_subsetting = bool(gene_caller_ids or function_accessions or function_names)
#         gene_codon_frequency_df = self._select_genes(
#             gene_caller_ids, function_accessions, function_names, expect_functions)

#         if from_function_sources and gene_subsetting:
#             # Subset gene function table to genes of interest.
#             gene_function_df = gene_function_df.set_index('gene_caller_id')
#             gene_function_df = gene_function_df.loc[
#                 gene_function_df.index.intersection(gene_codon_frequency_df.index)].reset_index()

#         if synonymous and not relative:
#             relative = True

#         if sum_genes and average_genes:
#             raise ConfigError("`sum_genes` and `average_genes` cannot both be True.")

#         if min_codon_filter not in ['length', 'remaining', 'both']:
#             raise ConfigError("`min_codon_filter` must be one of 'length', 'remaining', or 'both'.")
#         # Set `filter_gene_length` `filter_function_length`, and `filter_remaining`.
#         if min_codon_filter == 'length':
#             filter_gene_length = bool(gene_min_codons)
#             filter_function_length = bool(function_min_codons)
#             filter_gene_remaining_codons = False
#             filter_function_remaining_codons = False
#         elif min_codon_filter == 'remaining':
#             filter_gene_length = False
#             filter_function_length = False
#             filter_gene_remaining_codons = bool(gene_min_codons)
#             filter_function_remaining_codons = bool(function_min_codons)
#         elif min_codon_filter == 'both':
#             filter_gene_length = bool(gene_min_codons)
#             filter_function_length = bool(function_min_codons)
#             filter_gene_remaining_codons = bool(gene_min_codons)
#             filter_function_remaining_codons = bool(function_min_codons)

#         # Set `drop_amino_acids`.
#         if drop_amino_acids is None:
#             if synonymous:
#                 drop_amino_acids = ignored_cub_amino_acids
#             else:
#                 drop_amino_acids = []
#         else:
#             unrecognized_amino_acids = []
#             for amino_acid in drop_amino_acids:
#                 if amino_acid not in self.amino_acid_codons_dict:
#                     unrecognized_amino_acids.append(amino_acid)
#             if unrecognized_amino_acids:
#                 raise ConfigError("The following amino acids in `drop_amino_acids` are not "
#                                   f"recognized: {', '.join(unrecognized_amino_acids)}")

#         if (type(pansequence_min_amino_acids[0]) != int or
#             pansequence_min_amino_acids[0] < 0 or
#             type(pansequence_min_amino_acids[1]) != float or
#             not (0 < pansequence_min_amino_acids[1] <= 1)):
#             raise ConfigError(
#                 "The value of `pansequence_min_amino_acids` must be a tuple with two items: the "
#                 "first a positive int and the second a float between 0 and 1.")

#         if gene_subsetting:
#             self.run.info_single(f"{pp(len(gene_codon_frequency_df))} of "
#                                  f"{pp(len(self.gene_caller_ids) - self.noncoding_gene_count)} CDS "
#                                  "selected from the genome")
#         else:
#             self.run.info_single(f"{pp(len(gene_codon_frequency_df))} CDS in the genome")

#         # FILTER GENES/FUNCTIONS AND CODONS USING ADDITIONAL PARAMETERS
#         ###############################################################
#         # Filter genes by length.
#         gene_codon_frequency_df = self._get_frequency_table(
#             gene_codon_frequency_df,
#             min_codons=gene_min_codons,
#             filter_input_codon_count=filter_gene_length)

#         # Filter functions by the total length of their genes. Gene/function pairs with fewer than
#         # the required number of codons in the function are removed.
#         if function_sources:
#             gene_function_codon_frequency_df = gene_function_df.merge(
#                 gene_codon_frequency_df, how='inner', on='gene_caller_id')
#             gene_function_codon_frequency_df = gene_function_codon_frequency_df.set_index(
#                 ['source', 'accession', 'name', 'gene_caller_id'])
#             if filter_function_length:
#                 gene_function_codon_frequency_df = gene_function_codon_frequency_df.groupby(
#                     ['source', 'accession', 'name']).filter(
#                         lambda function_df: function_df.sum(axis=1).sum() >= function_min_codons)

#         # Drop certain codon columns from the gene codon frequency table. Filter genes by total
#         # codons remaining.
#         gene_codon_frequency_df = self._get_frequency_table(
#             gene_codon_frequency_df,
#             drop_amino_acids=drop_amino_acids,
#             pansequence_min_amino_acids=pansequence_min_amino_acids,
#             min_codons=gene_min_codons,
#             filter_output_codon_count=filter_gene_remaining_codons)

#         need_to_filter_codons_in_gene_function_codon_frequency_table = True
#         if ((drop_amino_acids or pansequence_min_amino_acids) and
#             (function_min_codons and filter_function_remaining_codons)):
#             # Filter functions by total codons remaining. Gene/function pairs with fewer than the
#             # required number of codons in the function are removed.
#             if function_sources:
#                 gene_function_codon_frequency_df = gene_function_df.merge(
#                     gene_codon_frequency_df, how='inner', on='gene_caller_id')
#                 gene_function_codon_frequency_df = gene_function_codon_frequency_df.set_index(
#                     ['source', 'accession', 'name', 'gene_caller_id'])
#                 if filter_function_remaining_codons:
#                     gene_function_codon_frequency_df = gene_function_codon_frequency_df.groupby(
#                         ['source', 'accession', 'name']).filter(
#                             lambda function_df:
#                                 function_df.sum(axis=1).sum() >= function_min_codons)
#                 need_to_filter_codons_in_gene_function_codon_frequency_table = False

#         if function_sources and need_to_filter_codons_in_gene_function_codon_frequency_table:
#             gene_function_codon_frequency_df = self._get_frequency_table(
#                 gene_function_codon_frequency_df,
#                 drop_amino_acids=drop_amino_acids,
#                 pansequence_min_amino_acids=pansequence_min_amino_acids)

#         if gene_min_codons:
#             self.run.info_single(f"{pp(len(gene_codon_frequency_df))} CDS remaining after codon "
#                                  "count filters")

#         if pansequence_min_amino_acids[0] > 0 and pansequence_min_amino_acids[1] < 1:
#             if gene_min_codons and min_codon_filter != 'remaining':
#                 min_gene_length_message = '≥' + pp(str(gene_min_codons)) + ' codon'
#             else:
#                 min_gene_length_message = ''
#             dynamically_dropped_amino_acids = set()
#             for codon, amino_acid in self.codon_amino_acid_dict.items():
#                 if (codon not in gene_codon_frequency_df.columns and
#                     amino_acid not in drop_amino_acids):
#                     dynamically_dropped_amino_acids.add(amino_acid)
#             if dynamically_dropped_amino_acids:
#                 self.run.warning(
#                     "Codons for the following amino acids were dropped as they did not meet the "
#                     f"threshold of {pp(pansequence_min_amino_acids[0])} codons in "
#                     f"{pansequence_min_amino_acids[1] * 100}% of {min_gene_length_message} CDS: "
#                     f"{', '.join(sorted(dynamically_dropped_amino_acids))}")

#         # GET OUTPUTS WITH NO CONSIDERATION OF FUNCTIONS
#         ################################################
#         get_table = lambda method: method(gene_codon_frequency_df,
#                                           sequence_min_amino_acids=sequence_min_amino_acids,
#                                           label_amino_acids=label_amino_acids,
#                                           replace_na=infinity_to_zero,
#                                           output_amino_acids=return_amino_acids)
#         ### Absolute frequencies
#         if (not relative and
#             not function_sources and
#             not sum_genes and
#             not average_genes):
#             return get_table(self._get_frequency_table)
#         ### Relative frequencies
#         if (relative and
#             not synonymous and
#             not function_sources and
#             not sum_genes and
#             not average_genes):
#             return get_table(self._get_rel_frequency_table)
#         ### Synonymous (per-amino acid) relative frequencies
#         if (not return_amino_acids and
#             relative and
#             synonymous and
#             not function_sources and
#             not sum_genes and
#             not average_genes):
#             return get_table(self._get_synonymous_codon_rel_frequency_table)

#         # Summed codon/amino acid results across genes:
#         ### Absolute frequencies
#         get_table = lambda method: self._get_frequency_table(
#             method(gene_codon_frequency_df),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids,
#             output_amino_acids=return_amino_acids)
#         if (not relative and
#             not function_sources and
#             sum_genes):
#             return get_table(self._get_summed_frequency_table)
#         ### Relative frequencies
#         get_table = lambda method: method(
#             gene_codon_frequency_df,
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids,
#             output_amino_acids=return_amino_acids)
#         if (relative and
#             not synonymous and
#             not function_sources and
#             sum_genes):
#             return get_table(self._get_summed_rel_frequency_table)
#         ### Synonymous relative frequencies
#         get_table = lambda method: method(
#             gene_codon_frequency_df,
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids,
#             replace_na=infinity_to_zero)
#         if (not return_amino_acids and
#             relative and
#             synonymous and
#             not function_sources and
#             sum_genes):
#             return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

#         # Codon results averaged across genes:
#         get_table = lambda method: self._get_frequency_table(
#             method(gene_codon_frequency_df),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids)
#         ### Absolute frequencies
#         if (not return_amino_acids and
#             not relative and
#             not function_sources and
#             average_genes):
#             return get_table(self._get_average_frequency_table)
#         ### Relative frequencies
#         get_table = lambda method: method(
#             gene_codon_frequency_df,
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids,
#             replace_na=infinity_to_zero)
#         if (not return_amino_acids and
#             relative and
#             not synonymous and
#             not function_sources and
#             average_genes):
#             return get_table(self._get_average_rel_frequency_table)
#         ### Synonymous relative frequencies
#         if (not return_amino_acids and
#             relative and
#             synonymous and
#             not function_sources and
#             average_genes):
#             return get_table(self._get_average_synonymous_codon_rel_frequency_table)

#         # Amino acid results averaged across genes:
#         # Gene amino acid frequencies are first calculated from codon frequencies before averaging
#         # across genes.
#         ### Absolute frequencies
#         get_table = lambda method: self._get_frequency_table(
#             method(gene_codon_frequency_df),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             output_amino_acids=return_amino_acids)
#         if (return_amino_acids and
#             not relative and
#             not function_sources and
#             average_genes):
#             return get_table(self._get_average_frequency_table)
#         ### Relative frequencies
#         if (return_amino_acids and
#             relative and
#             not synonymous and
#             not function_sources and
#             average_genes):
#             return get_table(self._get_average_rel_frequency_table)

#         # GET OUTPUTS WITH CONSIDERATION OF FUNCTIONS
#         #############################################
#         get_table = lambda method: method(
#             (gene_function_codon_frequency_df.groupby(['source', 'accession', 'name']).sum()
#              if return_functions else
#              gene_function_codon_frequency_df.sort_values(['source', 'accession', 'name'])),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids,
#             replace_na=infinity_to_zero,
#             output_amino_acids=return_amino_acids)
#         ### Absolute frequencies
#         if (not relative and
#             function_sources and
#             not sum_genes and
#             not average_genes):
#             return get_table(self._get_frequency_table)
#         ### Relative frequencies
#         if (relative and
#             not synonymous and
#             function_sources and
#             not sum_genes and
#             not average_genes):
#             return get_table(self._get_rel_frequency_table)
#         ### Synonymous relative frequencies
#         if (not return_amino_acids and
#             relative and
#             synonymous and
#             function_sources and
#             not sum_genes and
#             not average_genes):
#             return get_table(self._get_synonymous_codon_rel_frequency_table)

#         # Remove duplicate occurrences of genes before summing or averaging frequencies of all genes
#         # in the function source. A gene can have different annotations, e.g., different KOfam
#         # assignments. In KEGG BRITE, a gene can be in nested categories and different hierarchies.
#         gene_function_codon_frequency_df = \
#             gene_function_codon_frequency_df.reset_index().groupby('source').apply(
#                 lambda source_df: source_df.drop_duplicates(
#                     subset='gene_caller_id', ignore_index=True)).set_index(
#                         ['source', 'accession', 'name', 'gene_caller_id'])

#         # Summed codon/amino acid results across each function source:
#         # If amino acid rather than codon columns are returned, then gene amino acid frequencies are
#         # first calculated from the sum of codon frequencies.
#         ### Absolute frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             output_amino_acids=return_amino_acids,
#             label_amino_acids=label_amino_acids)
#         if (not relative and
#             function_sources and
#             sum_genes):
#             return get_table(self._get_summed_frequency_table)
#         ### Relative frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(
#                 partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
#             output_amino_acids=return_amino_acids,
#             label_amino_acids=label_amino_acids)
#         if (relative and
#             not synonymous and
#             function_sources and
#             sum_genes):
#             return get_table(self._get_summed_rel_frequency_table)
#         ### Synonymous relative frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(
#                 partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
#             label_amino_acids=label_amino_acids,
#             replace_na=infinity_to_zero)
#         if (not return_amino_acids and
#             relative and
#             synonymous and
#             function_sources and
#             sum_genes):
#             return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

#         # Codon results averaged across genes annotated by a function source:
#         ### Absolute frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             label_amino_acids=label_amino_acids)
#         if (not return_amino_acids and
#             not relative and
#             function_sources and
#             average_genes):
#             return get_table(self._get_average_frequency_table)
#         ### Relative frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(
#                 partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
#             label_amino_acids=label_amino_acids)
#         if (not return_amino_acids and
#             relative and
#             not synonymous and
#             function_sources and
#             average_genes):
#             return get_table(self._get_average_rel_frequency_table)
#         ### Synonymous relative frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(
#                 partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
#             label_amino_acids=label_amino_acids,
#             replace_na=infinity_to_zero)
#         if (not return_amino_acids and
#             relative and
#             synonymous and
#             function_sources and
#             average_genes):
#             return get_table(self._get_average_synonymous_codon_rel_frequency_table)

#         # Amino acid results averaged across genes annotated by a function source:
#         ### Absolute frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             output_amino_acids=return_amino_acids)
#         if (return_amino_acids and
#             not relative and
#             function_sources and
#             average_genes):
#             return get_table(self._get_average_frequency_table)
#         ### Relative frequencies
#         get_table = lambda method: self._get_frequency_table(
#             gene_function_codon_frequency_df.groupby('source').apply(
#                 partial(method, sequence_min_amino_acids=sequence_min_amino_acids)).droplevel(1),
#             output_amino_acids=return_amino_acids)
#         if (return_amino_acids and
#             relative and
#             function_sources and
#             average_genes):
#             return get_table(self._get_average_rel_frequency_table)

#         raise ConfigError("This point should not be reached at the end of the method, "
#                           "`get_frequencies`. Please contact the developers. Since you found the "
#                           "end of the earth, you now get to hear a top secret mnemonic for the "
#                           "rare earth elements. Scandalous Yiddish language centers praise Ned's "
#                           "promise of small European garden tubs. Dinosaurs hobble erotically "
#                           "thrumming yellow lutes. (scandium Sc, yttrium Y, lanthanum La, cerium "
#                           "Ce, praseodymium Pr, neodymium Nd, promethium Pm, samarium Sm, europium "
#                           "Eu, gadolinium Gd, terbium Tb, dysprosium Dy, holmium Ho, erbium Er, "
#                           "thulium Tm, ytterbium Yb, lutetium Lu) Credit for the lanthanide series "
#                           "mnemonic goes to Martyn Poliakoff: "
#                           "https://www.youtube.com/watch?v=Q21clW0s0B8&ab_channel=PeriodicVideos")


#     # `self.get_frequencies` SETUP HELPER METHODS
#     #############################################
#     def _establish_function_sources(self, from_function_sources):
#         """Gets `function_sources` for methods that take the argument, `from_function_sources`."""
#         if from_function_sources == True:
#             function_sources = list(self.function_sources)
#             if not function_sources:
#                 raise ConfigError(
#                     "All function sources were requested, but none exist for the object.")
#         elif from_function_sources != False:
#             function_sources = list(from_function_sources)
#         else:
#             function_sources = None

#         if function_sources is not None:
#             unrecognized_function_sources = []
#             for function_source in function_sources:
#                 if function_source not in self.function_sources:
#                     unrecognized_function_sources.append(function_source)
#             if unrecognized_function_sources:
#                 raise ConfigError("The requested function annotation sources, "
#                                   f"{', '.join(function_sources)}, are not among those available: "
#                                   f"{', '.join(self.function_sources)}.")

#         return function_sources


#     def _select_genes(self,
#                       gene_caller_ids,
#                       function_accession_dict,
#                       function_name_dict,
#                       expect_functions):
#         """Select genes in the gene frequency table given a list of IDs and/or functions
#         of interest."""
#         if not gene_caller_ids and not function_accession_dict and not function_name_dict:
#             return self.gene_codon_frequency_df

#         select_gene_caller_ids = []

#         unrecognized_gene_caller_ids = []
#         for gene_caller_id in gene_caller_ids:
#             if gene_caller_id not in self.gene_caller_ids:
#                 unrecognized_gene_caller_ids.append(gene_caller_id)
#                 continue
#             select_gene_caller_ids.append(gene_caller_id)
#         if unrecognized_gene_caller_ids:
#             raise ConfigError("The following gene caller IDs were not found in the genome: "
#                               f"{unrecognized_gene_caller_ids}")

#         index_keys = []
#         unrecognized_sources = []
#         for function_source, function_accessions in function_accession_dict.items():
#             if function_source == 'KEGG_BRITE':
#                 self.run.warning("Nota bene: KEGG BRITE accessions stored in anvi'o are for "
#                                  "hierarchies as a whole, not categories of the hierarchy. Most "
#                                  "hierarchies do not have category accessions. So all genes in the "
#                                  "selected hierarchies are being analyzed.")
#             if function_source not in self.function_sources:
#                 unrecognized_sources.append(function_source)
#             for function_accession in function_accessions:
#                 index_keys.append((function_source, function_accession))
#         if unrecognized_sources:
#             raise ConfigError("The following annotation sources in `function_accessions` were not "
#                               f"found as having annotated the genome: {unrecognized_sources}")
#         gene_function_df = self.gene_function_df.set_index(['source', 'accession'])
#         if expect_functions:
#             index_keys = set(index_keys)
#             all_keys = set(gene_function_df.index)
#             missing_keys = index_keys.difference(all_keys)

#             if missing_keys:
#                 missing_key_dict = {}
#                 for missing_key in missing_keys:
#                     try:
#                         missing_key_dict[missing_key[0]].append(missing_key[1])
#                     except KeyError:
#                         missing_key_dict[missing_key[0]] = [missing_key[1]]
#                 missing_key_message = ''
#                 for source, missing_accessions in missing_key_dict.items():
#                     missing_key_message += source + ': ' + ', '.join(missing_accessions) + '; '
#                 missing_key_message = missing_key_message[: -2]
#                 if missing_keys:
#                     raise ConfigError("The following requested function accessions are missing "
#                                       f"from the genome: {missing_key_message}")
#         else:
#             index_keys = gene_function_df.index.intersection(index_keys)
#         select_gene_caller_ids += gene_function_df.loc[index_keys]['gene_caller_id'].tolist()

#         index_keys = []
#         unrecognized_sources = []
#         for function_source, function_names in function_name_dict.items():
#             if function_source not in self.function_sources:
#                 unrecognized_sources.append(function_source)
#             for function_name in function_names:
#                 index_keys.append((function_source, function_name))
#         if unrecognized_sources:
#             raise ConfigError("The following annotation sources in `function_names` were not "
#                               f"found as having annotated the genome: {unrecognized_sources}")
#         gene_function_df = gene_function_df.reset_index().set_index(['source', 'name'])
#         if expect_functions:
#             index_keys = set(index_keys)
#             all_keys = set(gene_function_df.index)
#             missing_keys = index_keys.difference(all_keys)

#             if missing_keys:
#                 missing_key_dict = {}
#                 for missing_key in missing_keys:
#                     try:
#                         missing_key_dict[missing_key[0]].append(missing_key[1])
#                     except KeyError:
#                         missing_key_dict[missing_key[0]] = [missing_key[1]]
#                 missing_key_message = ''
#                 for source, missing_names in missing_key_dict.items():
#                     missing_key_message += source + ': ' + ', '.join(missing_names) + '; '
#                 missing_key_message = missing_key_message[: -2]
#                 if missing_keys:
#                     raise ConfigError("The following requested function names are missing from the "
#                                       f"genome: {missing_key_message}")
#         else:
#             index_keys = gene_function_df.index.intersection(index_keys)
#         select_gene_caller_ids += gene_function_df.loc[index_keys]['gene_caller_id'].tolist()

#         return self.gene_codon_frequency_df.loc[set(select_gene_caller_ids)]


#     # `self.get_frequencies` OUTPUT HELPER METHODS
#     ##############################################
#     def _filter_input_codon_count_decorator(method):
#         """Decorator to discard rows in the input frequency table with fewer than the minimum number
#         of codons/amino acids."""
#         def wrapper(*args, **kwargs):
#             frequency_df = args[1]
#             try:
#                 filter_input_codon_count = kwargs['filter_input_codon_count']
#                 min_codons = kwargs['min_codons']
#             except KeyError:
#                 filter_input_codon_count = False
#                 min_codons = 0
#             if filter_input_codon_count and min_codons:
#                 frequency_df = frequency_df[frequency_df.sum(axis=1) >= min_codons]
#             return method(args[0], frequency_df, *args[2: ], **kwargs)
#         return wrapper


#     def _drop_amino_acid_codon_columns_decorator(method):
#         """Decorator to discard codon columns by amino acid from the input frequency table."""
#         def wrapper(*args, **kwargs):
#             codon_frequency_df = args[1]
#             try:
#                 drop_amino_acids = kwargs['drop_amino_acids']
#             except KeyError:
#                 drop_amino_acids = []
#             if drop_amino_acids:
#                 drop_codons = []
#                 amino_acid_codons_dict = args[0].amino_acid_codons_dict
#                 for amino_acid in drop_amino_acids:
#                     drop_codons += amino_acid_codons_dict[amino_acid]
#                 codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1, errors='ignore')
#                 if len(codon_frequency_df.columns) == 0:
#                     codon_frequency_df = codon_frequency_df.drop(codon_frequency_df.index)
#             return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
#         return wrapper


#     def _filter_pansequence_synonymous_codon_count_decorator(method):
#         """Decorator to drop synonymous codons encoding an amino acid from the input frequency
#         table based on the frequency of the amino acid across rows."""
#         def wrapper(*args, **kwargs):
#             codon_frequency_df = args[1]
#             try:
#                 pansequence_min_amino_acids = kwargs['pansequence_min_amino_acids']
#             except KeyError:
#                 pansequence_min_amino_acids = (0, 1.0)
#             if pansequence_min_amino_acids[0] > 0 and 0 < pansequence_min_amino_acids[1] < 1:
#                 row_count = len(codon_frequency_df)
#                 drop_codons = []
#                 for synonymous_codons in args[0].amino_acid_codons_dict.values():
#                     try:
#                         filtered_row_count = len(codon_frequency_df[codon_frequency_df[
#                             synonymous_codons].sum(axis=1) >= pansequence_min_amino_acids[0]])
#                     except KeyError:
#                         # This occurs when codons are missing from the frequency table.
#                         continue
#                     if filtered_row_count / row_count < pansequence_min_amino_acids[1]:
#                         drop_codons += synonymous_codons
#                 codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1)
#                 if len(codon_frequency_df.columns) == 0:
#                     codon_frequency_df = codon_frequency_df.drop(codon_frequency_df.index)
#             return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
#         return wrapper


#     def _filter_sequence_synonymous_codon_count(self, codon_frequency_df, sequence_min_amino_acids):
#         mask_df = pd.DataFrame()
#         for synonymous_codons in self.amino_acid_codons_dict.values():
#             try:
#                 codon_mask_series = \
#                     codon_frequency_df[synonymous_codons].sum(axis=1) >= sequence_min_amino_acids
#             except KeyError:
#                 # This occurs when codons are missing from the frequency table.
#                 continue
#             amino_acid_mask_df = pd.DataFrame(index=codon_mask_series.index)
#             for codon in synonymous_codons:
#                 amino_acid_mask_df[codon] = codon_mask_series
#             mask_df = pd.concat([mask_df, amino_acid_mask_df], axis=1)
#         codon_frequency_df = codon_frequency_df[mask_df]
#         return codon_frequency_df


#     def _filter_sequence_synonymous_codon_count_decorator(method):
#         """Decorator to replace data with NaN for synonymous codons encoding an amino acid from the
#         input frequency table based on the summed frequency of the synonymous codons in each row."""
#         def wrapper(*args, **kwargs):
#             codon_frequency_df = args[1]
#             try:
#                 sequence_min_amino_acids = kwargs['sequence_min_amino_acids']
#             except KeyError:
#                 sequence_min_amino_acids = 0
#             if sequence_min_amino_acids > 0:
#                 mask_df = pd.DataFrame()
#                 for synonymous_codons in args[0].amino_acid_codons_dict.values():
#                     try:
#                         codon_mask_series = (codon_frequency_df[synonymous_codons].sum(axis=1) >=
#                                              sequence_min_amino_acids)
#                     except KeyError:
#                         # This occurs when codons are missing from the frequency table.
#                         continue
#                     amino_acid_mask_df = pd.DataFrame(index=codon_mask_series.index)
#                     for codon in synonymous_codons:
#                         amino_acid_mask_df[codon] = codon_mask_series
#                     mask_df = pd.concat([mask_df, amino_acid_mask_df], axis=1)
#                 codon_frequency_df = codon_frequency_df[mask_df]
#             return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
#         return wrapper


#     def _output_amino_acids_decorator(method):
#         """Decorator to output columns of amino acid rather than codon frequencies."""
#         def wrapper(*args, **kwargs):
#             codon_df = method(*args, **kwargs)
#             try:
#                 output_amino_acids = kwargs['output_amino_acids']
#             except KeyError:
#                 output_amino_acids = False
#             if output_amino_acids:
#                 aa_df = pd.DataFrame(index=codon_df.index)
#                 for amino_acid, codons in args[0].amino_acid_codons_dict.items():
#                     try:
#                         aa_df[amino_acid] = codon_df[codons].sum(axis=1, skipna=False)
#                     except KeyError:
#                         # This occurs when there aren't columns for codons.
#                         pass
#                 return aa_df
#             return codon_df
#         return wrapper


#     def _add_amino_acid_to_header_decorator(method):
#         """Decorator to add amino acid to codon column header."""
#         def wrapper(*args, **kwargs):
#             codon_df = method(*args, **kwargs)
#             try:
#                 label_amino_acids = kwargs['label_amino_acids']
#             except KeyError:
#                 label_amino_acids = False
#             if label_amino_acids:
#                 try:
#                     codon_df.columns = [
#                         constants.codon_to_AA[codon] + codon for codon in codon_df.columns]
#                 except KeyError:
#                     raise ConfigError("The columns in what should be a table of codon data are not "
#                                       "recognized as codons. This is the header that is present: "
#                                       f"{', '.join(list(codon_df.columns))}")
#                 codon_df = codon_df[sorted(codon_df.columns)]
#             return codon_df
#         return wrapper


#     def _replace_na(method):
#         """Decorator to replace NA with 0.0 in the output frequency table. This should only occur in
#         synonymous relative frequency output, in which missing amino acids yield NA (not inf, though
#         this is as a consequence of 0/0)."""
#         def wrapper(*args, **kwargs):
#             frequency_df = method(*args, **kwargs)
#             try:
#                 replace_na = kwargs['replace_na']
#             except KeyError:
#                 replace_na = False
#             if replace_na:
#                 frequency_df = frequency_df.fillna(0)
#             return frequency_df
#         return wrapper


#     def _filter_output_codon_count(method):
#         """Decorator to discard rows in the output frequency table with fewer than the minimum
#         number of codons/amino acids."""
#         def wrapper(*args, **kwargs):
#             frequency_df = method(*args, **kwargs)
#             try:
#                 filter_output_codon_count = kwargs['filter_output_codon_count']
#                 min_codons = kwargs['min_codons']
#             except KeyError:
#                 filter_output_codon_count = False
#             if filter_output_codon_count:
#                 frequency_df = frequency_df[frequency_df.sum(axis=1) >= min_codons]
#             return frequency_df
#         return wrapper


#     # The order of decorators should not be changed (only @_output_amino_acids_decorator and
#     # @_add_amino_acid_to_header_decorator, which are mutually exclusive operations, are
#     # interchangeable).
#     @_filter_input_codon_count_decorator
#     @_drop_amino_acid_codon_columns_decorator
#     @_filter_pansequence_synonymous_codon_count_decorator
#     @_filter_sequence_synonymous_codon_count_decorator
#     @_output_amino_acids_decorator
#     @_add_amino_acid_to_header_decorator
#     @_filter_output_codon_count
#     def _get_frequency_table(self, codon_frequency_df, **kwargs):
#         return codon_frequency_df


#     # Commented decorators mean that they can theoretically be used and uncommented but are not
#     # because they are not needed in the `self.get_frequencies` client.
#     # @_filter_input_codon_count_decorator
#     # @_drop_amino_acid_codon_columns_decorator
#     # @_filter_pansequence_synonymous_codon_count_decorator
#     @_output_amino_acids_decorator
#     @_add_amino_acid_to_header_decorator
#     # @_filter_output_codon_count
#     def _get_rel_frequency_table(self, codon_frequency_df, **kwargs):
#         try:
#             mask_df = self._filter_sequence_synonymous_codon_count(
#                 codon_frequency_df, kwargs['sequence_min_amino_acids']).notna()
#         except KeyError:
#             mask_df = None
#         codon_rel_frequency_df = codon_frequency_df.div(codon_frequency_df.sum(axis=1), axis=0)
#         drop_index = codon_rel_frequency_df[codon_rel_frequency_df.isna().all(axis=1)].index
#         if mask_df is not None:
#             codon_rel_frequency_df = codon_rel_frequency_df[mask_df]
#         # Drop rows with zero frequency.
#         codon_rel_frequency_df = codon_rel_frequency_df.drop(drop_index)
#         return codon_rel_frequency_df


#     # @_filter_input_codon_count_decorator
#     # @_drop_amino_acid_codon_columns_decorator
#     # @_filter_pansequence_synonymous_codon_count_decorator
#     @_filter_sequence_synonymous_codon_count_decorator
#     @_add_amino_acid_to_header_decorator
#     @_replace_na
#     # @_filter_output_codon_count
#     def _get_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
#         """Return the relative frequencies of codons in relation to the set of codons encoding the
#         same amino acid (or stop codons). If columns for one or more codons in a synonymous set are
#         missing (rather than 0), synonymous relative frequency will not be calculated for the
#         remaining codons in the set."""
#         try:
#             mask_df = self._filter_sequence_synonymous_codon_count(
#                 codon_frequency_df, kwargs['sequence_min_amino_acids']).notna()
#         except KeyError:
#             mask_df = None
#         synonymous_codon_rel_frequency_df = pd.DataFrame()
#         for codons in constants.AA_to_codons.values():
#             try:
#                 aa_codon_frequency_df = codon_frequency_df[codons]
#             except KeyError:
#                 # This occurs when synonymous codons are missing from the frequency table.
#                 continue
#             synonymous_codon_rel_frequency_df[codons] = aa_codon_frequency_df.div(
#                 aa_codon_frequency_df.sum(axis=1), axis=0)
#         drop_index = synonymous_codon_rel_frequency_df[
#             synonymous_codon_rel_frequency_df.isna().all(axis=1)].index
#         if mask_df is not None:
#             synonymous_codon_rel_frequency_df = synonymous_codon_rel_frequency_df[mask_df]
#         synonymous_codon_rel_frequency_df = synonymous_codon_rel_frequency_df.drop(drop_index)
#         return synonymous_codon_rel_frequency_df


#     # @_filter_input_codon_count_decorator
#     # @_drop_amino_acid_codon_columns_decorator
#     # @_filter_pansequence_synonymous_codon_count_decorator
#     # @_filter_sequence_synonymous_codon_count_decorator
#     def _get_summed_frequency_table(self, frequency_df, **kwargs):
#         """Return the summed frequencies across all items."""
#         summed_frequency_series = frequency_df.sum()
#         summed_frequency_df = \
#             summed_frequency_series.to_frame('all').T.rename_axis('gene_caller_id')
#         return summed_frequency_df


#     def _get_summed_rel_frequency_table(self, frequency_df, **kwargs):
#         first_kwargs = {}
#         second_kwargs = {}
#         for key, value in kwargs.items():
#             if key in ['sequence_min_amino_acids', 'label_amino_acids']:
#                 second_kwargs[key] = value
#             else:
#                 first_kwargs[key] = value

#         summed_frequency_df = self._get_summed_frequency_table(frequency_df, **first_kwargs)
#         summed_rel_frequency_df = self._get_rel_frequency_table(
#             summed_frequency_df, **second_kwargs)
#         return summed_rel_frequency_df


#     def _get_summed_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
#         first_kwargs = {}
#         second_kwargs = {}
#         for key, value in kwargs.items():
#             if key in ['sequence_min_amino_acids', 'label_amino_acids', 'replace_na']:
#                 second_kwargs[key] = value
#             else:
#                 first_kwargs[key] = value

#         summed_codon_frequency_df = self._get_summed_frequency_table(
#             codon_frequency_df, **first_kwargs)
#         summed_synonymous_codon_rel_frequency_df = self._get_synonymous_codon_rel_frequency_table(
#             summed_codon_frequency_df, **second_kwargs)
#         return summed_synonymous_codon_rel_frequency_df


#     # @_filter_input_codon_count_decorator
#     # @_drop_amino_acid_codon_columns_decorator
#     # @_filter_synonymous_codon_count
#     def _get_average_frequency_table(self, frequency_df, **kwargs):
#         """Return the average codon frequencies across all items."""
#         average_codon_frequency_series = frequency_df.mean()
#         average_codon_frequency_df = \
#             average_codon_frequency_series.to_frame('all').T.rename_axis('gene_caller_id')
#         return average_codon_frequency_df


#     def _get_average_rel_frequency_table(self, frequency_df, **kwargs):
#         first_kwargs = {}
#         second_kwargs = {}
#         for key, value in kwargs.items():
#             if key in ['sequence_min_amino_acids', 'label_amino_acids']:
#                 second_kwargs[key] = value
#             else:
#                 first_kwargs[key] = value

#         average_frequency_df = self._get_average_frequency_table(frequency_df, **first_kwargs)
#         average_rel_frequency_df = \
#             self._get_rel_frequency_table(average_frequency_df, **second_kwargs)
#         return average_rel_frequency_df


#     def _get_average_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
#         first_kwargs = {}
#         second_kwargs = {}
#         for key, value in kwargs.items():
#             if key in ['sequence_min_amino_acids', 'label_amino_acids', 'replace_na']:
#                 second_kwargs[key] = value
#             else:
#                 first_kwargs[key] = value

#         average_codon_frequency_df = self._get_average_frequency_table(
#             codon_frequency_df, **first_kwargs)
#         average_synonymous_codon_rel_frequency_df = self._get_synonymous_codon_rel_frequency_table(
#             average_codon_frequency_df, **second_kwargs)
#         return average_synonymous_codon_rel_frequency_df


#     def get_codon_usage_bias(
#         self,
#         metrics=None,
#         from_function_sources=None,
#         gene_caller_ids=None,
#         function_accessions=None,
#         function_names=None,
#         expect_functions=False,
#         omnibias=False,
#         reference_function_accessions=None,
#         reference_function_names=None,
#         expect_reference_functions=False,
#         reference_gene_caller_ids=None,
#         gene_min_codons=0,
#         function_min_codons=0,
#         min_codon_filter='both',
#         drop_amino_acids=None,
#         sequence_min_amino_acids=0,
#         pansequence_min_amino_acids=(0, 1.0),
#         query_min_analyzed_codons=default_query_min_analyzed_codons,
#         reference_exclude_amino_acid_count=default_reference_exclude_amino_acid_count,
#         reference_min_analyzed_codons=default_reference_min_analyzed_codons):
#         """
#         Get codon usage bias (CUB) of genes or functions.

#         Parameters
#         ----------
#         metrics : {'cai', 'delta'}, optional
#             CUB metric, with valid choices being 'cai' (Codon Adaptation Index of Sharp and Li,
#             1987) and 'delta' (Ran and Higgs, 2012, Eq. 6). The default of None calculates all
#             available CUB metrics. For CUB metrics that require reference genes, the default
#             behavior in the absence of supplied reference genes/functions or omnibias mode is to use
#             the summed composition of ribosomal proteins defined by KEGG KOfams/BRITE.
#         from_function_sources : bool, str, or iterable of str, optional
#             Select genes with functional annotations, using their summed codon frequencies to
#             calculate function CUB values. With this argument, the first three columns of the
#             returned table contain, respectively, function annotation sources, function accessions,
#             and function names. When this argument is True, use all available functional annotation
#             sources in the SingleGenomeCodonUsage object. When this argument is a string, select the
#             source given by the string, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. When this argument
#             is an iterable, select a subset of sources given by its strings, e.g., ['KOfam',
#             'COG20_FUNCTION']. When 'KEGG_BRITE' is in the argument, in the absence of
#             `function_names` narrowing the scope of the inquiry, all possible depths of each BRITE
#             hierarchy in the data are returned, e.g., 'Ribosome>>>Ribosomal proteins' and the more
#             general 'Ribosome' would each have rows in the output table. By default None.
#         gene_caller_ids : iterable, optional
#             Genes with the given IDs are selected for calculation of CUB. This parameter can be
#             used alongside `function_accessions` and `function_names`. By default None.
#         function_accessions : dict, optional
#             Genes annotated with the given function accessions are selected for calculation of CUB.
#             The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
#             parameter can be used alongside `gene_caller_ids` and `function_names`.  Note that
#             'KEGG_BRITE' does not use individual function accessions but overarching hierarchy
#             accessions that include multiple functions. By default None.
#         function_names : dict, optional
#             Genes annotated with the given function names are selected for calculation of CUB. The
#             argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
#             `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
#             (categorizations) given to an arbitrary level of the hierarchy, such as
#             ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
#             type>>>Bacterial type>>>RNA polymerase']. This parameter can be used alongside
#             `gene_caller_ids` and `function_accessions`. By default None.
#         expect_functions : bool, optional
#             If True (default False), an error will be raised if any given `function_accessions` or
#             `function_names` are not annotated in the input genome.
#         omnibias : bool, optional
#             If True (default False), use every gene or function as a separate reference rather than
#             defining a set of reference genes or functions. The resulting table of gene x gene (or
#             function x function) CUB values is like a distance matrix of the similarity of gene
#             codon compositions.
#         reference_function_accessions : dict, optional
#             Genes annotated with the given function accessions are selected for the reference set.
#             The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
#             parameter can be used alongside `function_names`. Note that 'KEGG_BRITE' does not use
#             individual function accessions but overarching hierarchy accessions that include
#             multiple functions. By default None.
#         reference_function_names : dict, optional
#             Genes annotated with the given function names are selected for the reference set. The
#             argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
#             `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
#             (categorizations) given to an arbitrary level of the hierarchy, such as
#             ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
#             type>>>Bacterial type>>>RNA polymerase']. This parameter can be used alongside
#             `gene_caller_ids` and `function_accessions`. By default, if using reference-dependent
#             metrics without omnibias mode, this argument becomes {'KEGG_BRITE':
#             ['Ribosome>>>Ribosomal proteins']}.
#         expect_reference_functions : bool, optional
#             If True (default False), an error will be raised if any given
#             `reference_function_accessions` or `reference_function_names` are not annotated in the
#             input genome.
#         reference_gene_caller_ids : iterable, optional
#             Include specific genes in the reference gene set, as given by their gene caller IDs in
#             the contigs database.

#         Additional Parameters
#         ---------------------
#         These parameters filter genes and codons in CUB calculations.

#         The gene codon frequency table is the input to the CUB calculation. Here is the order of
#         all possible filters for \"queries\":
#         gene codon frequency table ->
#             drop genes on codon length ->
#             drop gene/function pairs on the total length of remaining genes defining the function ->
#             drop codons of defined amino acids ->
#             dynamically drop codons of rarer amino acids ->
#             drop genes on remaining codon frequency ->
#             drop gene/function pairs on the codon sum of remaining genes defining the function ->
#         filtered gene codon frequency table

#         There is a filter, `query_min_analyzed_codons`, that excludes queries based on the number of
#         codons participating in CUB analysis.

#         There are two filters, 'reference_exclude_amino_acid_count' and
#         'reference_min_analyzed_codons', that are applied to the reference gene set.

#         gene_min_codons : int, optional
#             Ignore genes with fewer than this number of codons in gene/function queries. By default
#             0.
#         function_min_codons : int, optional
#             Ignore gene/function pairs when the function has fewer than this number of codons. Genes
#             are filtered by codon count before functions, so the total length of remaining genes
#             annotated by the function is considered. By default 0.
#         min_codon_filter : {"length", "remaining", "both"}, optional
#             This argument arises from the ambiguity of filters that remove genes and functions by
#             number of codons (`--gene-min-codons` and `--function-min-codons`) in relation to the
#             filters that drop codons (`--exclude/include-amino-acids` and
#             `--exclude-amino-acid-count/fraction`). Genes (and functions) can be filtered by their
#             full length, e.g., genes shorter than 300 codons are ignored. They can also be
#             filtered by the number of codons remaining after dropping codons. The codon length
#             filter followed by dropping codons can result in genes and functions with fewer codons
#             than the original codon threshold -- thus the option of both 'length' and 'remaining'
#             filters to ensure that total codon frequencies always meet the minimum codon threshold.
#             'both' is needed as an option in addition to 'remaining' so dynamic codon filtering by
#             `--exclude-amino-acid-count/fraction` operates on genes that passed the first length
#             filter.
#         drop_amino_acids : iterable, optional
#             Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala).
#             Met and Trp are encoded by single codons, which perforce are excluded from CUB
#             calculations. By default, stop codons are also excluded from the calculation: the
#             default, None, becomes ['STP'] in the code. Importantly, to continue to exclude stop
#             codons, make sure to include it in the passed value: exclusion of Ala and stop codons
#             is achieved by passing ['Ala', 'STP'].
#         sequence_min_amino_acids : int, optional
#             Remove codons for amino acids (and STP) that are less numerous than
#             `sequence_min_amino_acids`. For example, if the argument is 5, and a gene query has 4
#             codons encoding Asn, 2 AAT and 2 AAC, then Asn codons will be disregarded in the
#             calculation of CUB for this query. By default 0.
#         pansequence_min_amino_acids : tuple, optional
#             This tuple must have two items. The first is an int >0 representing a minimum number of
#             codons encoding an amino acid -- 'min_amino_acids' -- and the second is a float in the
#             range (0, 1) representing a fraction of genes -- 'min_gene_fraction'. Remove codons for
#             amino acids (and STP) that are less numerous than 'min_amino_acids' in a
#             'min_gene_fraction' of genes. For example, if 'min_amino_acids' is 5 and
#             'min_gene_fraction' is 0.9, then if there are fewer than 5 codons for an amino acid/STP
#             in ≥90% of genes, the amino acid's codons do not factor into the calculation of CUB for
#             any query. By default (0, 1.0).
#         query_min_analyzed_codons : int, optional
#             Only allow CUB to calculated for a query if has at least this number of synonymous
#             codons that will be analyzed. For reference-dependent CUB metrics, analyzed codons are
#             those with reference compositions. By default 100.
#         reference_exclude_amino_acid_count : int, optional
#             Exclude codons for amino acids with fewer than this many codons in the set of reference
#             genes. This does not apply in `omnibias` mode. By default 5.
#         reference_min_analyzed_codons : int, optional
#             Only allow CUB to be calculated using a reference if it has at least this number of
#             codons that will be analyzed. This filter applies after excluding codons for individual
#             amino acids using `reference_exclude_amino_acid_count`. By default 100.

#         Returns
#         -------
#         dict of pandas.core.frame.DataFrame objects
#             Each dict key is the name of a CUB metric (from `metrics`), and each value is a
#             corresponding table of CUB data.
#         """
#         metrics = self._establish_cub_metrics(metrics)
#         reference_metrics = []
#         referenceless_metrics = []
#         for metric in metrics:
#             if metric in reference_dependent_cub_metrics:
#                 reference_metrics.append(metric)
#             else:
#                 referenceless_metrics.append(metric)

#         # Check that proper arguments were provided given the metrics.
#         if reference_metrics:
#             if (omnibias and
#                 (reference_function_accessions or
#                  reference_function_names or
#                  expect_reference_functions or
#                  reference_gene_caller_ids)):
#                 raise ConfigError(
#                     "Omnibias mode cannot be used when defined gene/function references are also "
#                     "used. The following arguments are only relevant to defined references: "
#                     "`reference_function_accessions`, `reference_function_names`, "
#                     "`expect_reference_functions`, `reference_gene_caller_ids`, "
#                     "`reference_exclude_amino_acid_count`, and `reference_min_analyzed_codons`. "
#                     "Query filters including `sequence_min_amino_acids`, "
#                     "`pansequence_min_amino_acids`, `min_gene_fraction`, and "
#                     "`query_min_analyzed_codons` apply to references since queries are the same as "
#                     "references in omnibias mode.")
#             if omnibias and (query_min_analyzed_codons != reference_min_analyzed_codons):
#                 raise ConfigError(
#                     "In omnibias mode, `query_min_analyzed_codons` and "
#                     "`reference_min_analyzed_codons` should have the same value, since the sets of "
#                     "query and reference genes/functions should be the same, and so should be "
#                     "filtered in the same way.")
#         elif (omnibias or
#               reference_function_accessions or
#               reference_function_names or
#               expect_reference_functions or
#               reference_gene_caller_ids):
#             raise ConfigError(
#                 "The provided CUB metrics do not involve comparison of gene/function codon "
#                 "compositions. The following arguments are only relevant to "
#                 "reference-dependent metrics: `omnibias`, `reference_function_accessions`, "
#                 "`reference_function_names`, `expect_reference_functions`, "
#                 "`reference_gene_caller_ids`, `reference_exclude_amino_acid_count`, and "
#                 "`reference_min_analyzed_codons`.")

#         # Get a reference codon composition when using reference-dependent metrics and not in
#         # omnibias mode.
#         if reference_metrics and not omnibias:
#             if not reference_function_accessions and not reference_function_names:
#                 reference_function_accessions = {
#                     default_reference_function_source: default_reference_function_accessions}
#                 reference_function_names = {
#                     default_reference_function_source: default_reference_function_names}
#             reference_codon_frequency_df = self._get_defined_reference_codon_frequencies(
#                 reference_metrics,
#                 reference_function_accessions=reference_function_accessions,
#                 reference_function_names=reference_function_names,
#                 expect_reference_functions=expect_reference_functions,
#                 reference_gene_caller_ids=reference_gene_caller_ids,
#                 reference_exclude_amino_acid_count=reference_exclude_amino_acid_count,
#                 reference_min_analyzed_codons=reference_min_analyzed_codons)
#         else:
#             reference_codon_frequency_df = None

#         # Determine the encoded amino acids ignored in the analysis.
#         if drop_amino_acids is None:
#             drop_amino_acids = ignored_cub_amino_acids
#         else:
#             for amino_acid in single_codon_amino_acids:
#                 if amino_acid not in drop_amino_acids:
#                     drop_amino_acids.append(amino_acid)

#         self.run.info_single("Queries", mc='green')
#         query_codon_frequency_df = self.get_frequencies(
#             from_function_sources=from_function_sources,
#             return_functions=bool(from_function_sources),
#             gene_caller_ids=gene_caller_ids,
#             function_accessions=function_accessions,
#             function_names=function_names,
#             expect_functions=expect_functions,
#             gene_min_codons=gene_min_codons,
#             function_min_codons=function_min_codons,
#             min_codon_filter=min_codon_filter,
#             drop_amino_acids=drop_amino_acids,
#             sequence_min_amino_acids=sequence_min_amino_acids,
#             pansequence_min_amino_acids=pansequence_min_amino_acids)
#         if (gene_caller_ids or
#             function_accessions or
#             function_names or
#             gene_min_codons or
#             function_min_codons):
#             self.run.info_single(
#                 f"{pp(len(query_codon_frequency_df))} "
#                 f"{'function' if from_function_sources else 'CDS'} queries retained")

#         # Report codons missing in the query dataset.
#         missing_query_codons = []
#         for codon in query_codon_frequency_df.columns:
#             if query_codon_frequency_df[codon].sum() == 0:
#                 missing_query_codons.append(codon)
#         if missing_query_codons:
#             missing_codon_dict = {}
#             for codon in missing_query_codons:
#                 amino_acid = self.codon_amino_acid_dict[codon]
#                 try:
#                     missing_codon_dict[amino_acid].append(codon)
#                 except KeyError:
#                     missing_codon_dict[amino_acid] = [codon]
#             missing_codon_message = ""
#             for amino_acid, codons in sorted(missing_codon_dict.items()):
#                 missing_codon_message += amino_acid + ": " + ", ".join(codons) + "; "
#             missing_codon_message = missing_codon_message[: -2]
#             self.run.warning(
#                 "The following codons from synonymous decoding sets are absent in the retained "
#                 f"queries. {missing_codon_message}")

#         if omnibias:
#             reference_codon_frequency_df = query_codon_frequency_df
#         else:
#             # Here a defined reference codon composition is being used. Codons that are absent in
#             # the query dataset obviously are not compared to the reference; check that the
#             # reference has enough compared codons to meet `reference_min_analyzed_codons`.
#             total_codons = reference_codon_frequency_df.drop(
#                 [codon for codon in missing_query_codons
#                  if codon in reference_codon_frequency_df.columns], axis=1).sum(axis=1).sum()
#             if total_codons < reference_min_analyzed_codons:
#                 # The reference does not have enough codons, so leave an empty shell of a DataFrame.
#                 try:
#                     reference_codon_frequency_df = reference_codon_frequency_df.drop('all')
#                     self.run.warning(
#                         "A reference codon composition could not be established. "
#                         "Reference-dependent CUB tables will be empty. The reference dataset has "
#                         f"{pp(int(total_codons))} codons that will be used in the CUB analysis, "
#                         "which does not meet the minimum defined by `reference_min_analyzed_codons`, "
#                         f"{pp(int(reference_min_analyzed_codons))}. Note that codons are only analyzed "
#                         "when present in both the query and reference datasets.")
#                 except KeyError:
#                     # Even before removing the codons not shared by the query, the reference did not
#                     # meet the minimum codon threshold.
#                     pass

#         if 'delta' in metrics:
#             # Reference codon weights in the computation of 𝛿 involve comparison of the reference to
#             # the non-reference codon set, taken from all other genes in the input dataset.
#             amino_acids_dropped_from_reference = []
#             for amino_acid, codons in constants.AA_to_codons.items():
#                 if codons[0] not in reference_codon_frequency_df.columns:
#                     amino_acids_dropped_from_reference.append(amino_acid)
#             default_run = self.run
#             self.run = run_quiet
#             nonreference_codon_frequency_df = self.get_frequencies(
#                 sum_genes=True, drop_amino_acids=amino_acids_dropped_from_reference)
#             self.run = default_run
#             if omnibias:
#                 # For perfect consistency with the "non-omnibias" mode using a defined reference
#                 # gene set, the codon frequencies of the reference gene/function in each omnibias
#                 # query-reference comparison would be subtracted from the non-reference codon
#                 # frequencies summed from all coding sequences in the genome, since the reference
#                 # composition is subtracted from the whole-genome composition to find the
#                 # non-reference composition in "non-omnibias" mode. For the sake of simplicity, this
#                 # detail is ignored.
#                 pass
#             else:
#                 nonreference_codon_frequency_df = \
#                     nonreference_codon_frequency_df - reference_codon_frequency_df.sum()
#                 # It is theoretically possible that there are fewer codons in the non-reference genes
#                 # than the reference genes. Ensure that the non-reference set meets the minimum codon
#                 # threshold of the reference set.
#                 total_codons = nonreference_codon_frequency_df.drop(
#                     [codon for codon in missing_query_codons
#                      if codon in nonreference_codon_frequency_df.columns], axis=1).sum(axis=1).sum()
#                 if total_codons < reference_min_analyzed_codons:
#                     nonreference_codon_frequency_df = nonreference_codon_frequency_df.drop('all')
#                     self.run.warning(
#                         "The non-reference codon composition needed for 𝛿 could not be "
#                         "established, because the non-reference genes have "
#                         f"{pp(int(total_codons))} codons, fewer than the minimum of "
#                         f"{pp(int(reference_min_analyzed_codons))} given by "
#                         "`reference_min_analyzed_codons` that is required for the reference, and "
#                         "by extension, non-reference datasets. Reference-dependent CUB tables will "
#                         "be empty.")
#         else:
#             nonreference_codon_frequency_df = None

#         cub_table_dict = {}
#         for metric in reference_metrics:
#             if metric == 'cai':
#                 cub_df = self._get_cai_table(
#                     query_codon_frequency_df,
#                     reference_codon_frequency_df,
#                     query_min_analyzed_codons=query_min_analyzed_codons,
#                     reference_min_analyzed_codons=reference_min_analyzed_codons)
#                 if len(cub_df.columns) > 0:
#                     if omnibias:
#                         cub_df.columns = cub_df.index
#                     else:
#                         cub_df.columns = ['CAI']
#                     self.run.info_single(
#                         f"CAI calculated for {pp(len(cub_df))} "
#                         f"{'function' if from_function_sources else 'CDS'} queries")
#             elif metric == 'delta':
#                 cub_df = self._get_delta_table(
#                     query_codon_frequency_df,
#                     reference_codon_frequency_df,
#                     nonreference_codon_frequency_df,
#                     query_min_analyzed_codons=query_min_analyzed_codons,
#                     reference_min_analyzed_codons=reference_min_analyzed_codons)
#                 if len(cub_df.columns) > 0:
#                     if omnibias:
#                         cub_df.columns = cub_df.index
#                     else:
#                         cub_df.columns = ['Delta']
#                     self.run.info_single(
#                         f"𝛿 calculated for {pp(len(cub_df))} "
#                         f"{'function' if from_function_sources else 'CDS'} queries")

#             cub_table_dict[metric] = cub_df

#             if len(cub_df.columns) == 0:
#                 continue

#             # Print the number of query-reference comparisons that were thrown out.
#             if omnibias:
#                 # This includes "self-comparisons," reported in the matrix diagonal.
#                 possible_comparison_count = int(len(cub_df) ** 2 / 2)
#             else:
#                 possible_comparison_count = len(cub_df)
#             unperformed_comparison_count = int(cub_df.isna().sum().sum())
#             if unperformed_comparison_count:
#                 self.run.warning(
#                     f"{pp(unperformed_comparison_count)} of {pp(possible_comparison_count)} "
#                     "query-reference comparisons did not meet the minimum codon threshold in "
#                     "either query or reference and so did not yield a CUB value.")

#         # Calculate CUB using metrics that do not depend on a reference codon composition.
#         for metric in referenceless_metrics:
#             # No reference-independent CUB metrics have been programmed yet.
#             pass

#         return cub_table_dict


#     def _establish_cub_metrics(self, metrics):
#         """Establishes CUB metrics for `get_codon_usage_bias`."""
#         if metrics is None:
#             metrics = available_cub_metrics
#         else:
#             for metric in metrics:
#                 unrecognized_metrics = []
#                 if metric not in available_cub_metrics:
#                     unrecognized_metrics.append(metric)
#                 if unrecognized_metrics:
#                     raise ConfigError("The following CUB metrics are not recognized: "
#                                       f"{', '.join(unrecognized_metrics)}. Here are the available "
#                                       f"CUB metrics: {', '.join(available_cub_metrics)}")
#         return metrics


#     def _get_defined_reference_codon_frequencies(
#         self,
#         reference_metrics,
#         reference_function_accessions=None,
#         reference_function_names=None,
#         expect_reference_functions=False,
#         reference_gene_caller_ids=None,
#         reference_exclude_amino_acid_count=default_reference_exclude_amino_acid_count,
#         reference_min_analyzed_codons=default_reference_min_analyzed_codons):
#         """
#         Get a codon frequency table from a set of reference genes or functions in the genome.

#         Parameters
#         ----------
#         reference_metrics : iterable
#             CUB metrics requiring a reference codon composition.
#         reference_function_accessions : dict, optional
#             Genes annotated with the given function accessions are selected for the reference set.
#             The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
#             parameter can be used alongside `function_names`. Note that 'KEGG_BRITE' does not use
#             individual function accessions but overarching hierarchy accessions that include
#             multiple functions. By default None.
#         reference_function_names : dict, optional
#             Genes annotated with the given function names are selected for the reference set. The
#             argument must be a dict keyed by function annotation source (e.g., 'KOfam',
#             'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
#             `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
#             (categorizations) given to an arbitrary level of the hierarchy, such as
#             ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
#             type>>>Bacterial type>>>RNA polymerase']. This parameter can be used alongside
#             `gene_caller_ids` and `function_accessions`. By default, {'KEGG_BRITE':
#             ['Ribosome>>>Ribosomal proteins']}.
#         expect_reference_functions : bool, optional
#             If True (default False), an error will be raised if any given
#             `reference_function_accessions` or `reference_function_names` are not annotated in the
#             input genome.
#         reference_gene_caller_ids : iterable, optional
#             Include specific genes in the reference gene set, as given by their gene caller IDs in
#             the contigs database.
#         reference_exclude_amino_acid_count : int, optional
#             Exclude codons for amino acids with fewer than this many codons in the set of reference
#             genes. This does not apply in `omnibias` mode. By default 5.
#         reference_min_analyzed_codons : int, optional
#             Only allow CUB to be calculated using a reference if it has at least this number of
#             codons that will be analyzed. This filter applies after excluding codons for individual
#             amino acids using `reference_exclude_amino_acid_count`. By default 100.

#         Returns
#         -------
#         pandas.core.frame.DataFrame
#             This frequency table has a single row for the reference composition and a column per
#             codon.
#         """
#         self.run.info_single("Reference codon composition", mc='green')

#         # The default reference genes are ribosomal proteins, as annotated by KOfams/BRITE.
#         if reference_function_names is None:
#             if 'KEGG_BRITE' not in self.function_sources:
#                 raise ConfigError(
#                     f"Reference-dependent metrics ({', '.join(reference_metrics)}) were "
#                     "requested without defined reference genes. By default, reference genes "
#                     "are KEGG KOfams classified as ribosomal proteins in BRITE. However, "
#                     "'KEGG_BRITE' is not among the function annotation sources run on the "
#                     "genome. This can be rectified by rerunning `anvi-run-kegg-kofams`.")
#             reference_function_names = {'KEGG_BRITE': ['Ribosome>>>Ribosomal proteins']}

#         reference_codon_frequency_df = self.get_frequencies(
#             gene_caller_ids=reference_gene_caller_ids,
#             function_accessions=reference_function_accessions,
#             function_names=reference_function_names,
#             expect_functions=expect_reference_functions,
#             sum_genes=True,
#             drop_amino_acids=ignored_cub_amino_acids)

#         # In the absence of a set of reference genes for the genome,
#         # `reference_codon_frequency_df` has an index and header but no data.
#         if len(reference_codon_frequency_df) == 0:
#             self.run.warning(
#                 "A reference codon composition could not be established because none of the "
#                 "requested genes or functions were identified in the genome. "
#                 "Reference-dependent CUB tables will be empty.")
#             return reference_codon_frequency_df

#         # Report codons absent in the reference.
#         missing_reference_codons = []
#         for codon in set(constants.codons).difference(set(self.ignored_cub_codons)):
#             if reference_codon_frequency_df[codon].sum() == 0:
#                 missing_reference_codons.append(codon)
#         missing_codon_dict = {}
#         for codon in missing_reference_codons:
#             amino_acid = self.codon_amino_acid_dict[codon]
#             try:
#                 missing_codon_dict[amino_acid].append(codon)
#             except KeyError:
#                 missing_codon_dict[amino_acid] = [codon]
#         missing_codon_message = ""
#         for amino_acid, codons in sorted(missing_codon_dict.items()):
#             missing_codon_message += amino_acid + ": " + ", ".join(codons) + "; "
#         missing_codon_message = missing_codon_message[: -2]
#         if missing_reference_codons:
#             self.run.warning(
#                 "The following synonymous codons are absent in the reference. These codons will "
#                 "not contribute to CUB. This can skew CUB if the query contains many codons that "
#                 f"are absent in the reference. {missing_codon_message}")

#         # Remove and report codons (columns) encoding amino acids that do not meet the minimum
#         # reference amino acid count threshold, `reference_exclude_amino_acid_count`. Note that for
#         # technical reasons requiring the preservation of columns for complete sets of synonymous
#         # codons, columns representing complete sets of synonymous codons may be removed in this
#         # method, but columns for individual missing codons are not removed.
#         if reference_exclude_amino_acid_count:
#             removed_amino_acids = []
#             removed_codons = []
#             for amino_acid, codons in self.synonymous_nonstop_amino_acid_codons_dict.items():
#                 present_codons = set(reference_codon_frequency_df.columns.intersection(codons))
#                 if (reference_codon_frequency_df[present_codons].sum().sum()
#                     < reference_exclude_amino_acid_count):
#                     removed_amino_acids.append(amino_acid)
#                     removed_codons += list(present_codons)
#             reference_codon_frequency_df = reference_codon_frequency_df.drop(removed_codons, axis=1)
#             if removed_amino_acids:
#                 self.run.warning(
#                     "The following amino acids do not meet the threshold of "
#                     f"{reference_exclude_amino_acid_count} codons in the reference gene set, and "
#                     "so the codons were excluded from the reference: "
#                     f"{', '.join(removed_amino_acids)}")

#         # Check that there are enough total codons in the reference set.
#         if reference_min_analyzed_codons:
#             total_codon_count = reference_codon_frequency_df.sum().sum()
#             if total_codon_count < reference_min_analyzed_codons:
#                 reference_codon_frequency_df = reference_codon_frequency_df.drop('all')
#                 self.run.warning(
#                     f"The reference gene set has {pp(int(total_codon_count))} codons, not meeting "
#                     f"the minimum codon threshold of {pp(int(reference_min_analyzed_codons))}, so "
#                     "a reference codon composition could not be established. Reference-dependent "
#                     "CUB tables will be empty.")

#         return reference_codon_frequency_df


#     def _get_cai_table(self,
#                        query_codon_frequency_df,
#                        reference_codon_frequency_df,
#                        query_min_analyzed_codons=default_query_min_analyzed_codons,
#                        reference_min_analyzed_codons=default_reference_min_analyzed_codons):
#         """
#         Get a table of CAI (Sharp and Li, 1987) values for each query x reference comparison.

#         Calculation of CAI:
#         reference_codon_weight =
#             ln(reference_codon_frequency / reference_max_synonymous_codon_frequency)
#         weighted_codon_count = ∑(codon_frequency * reference_codon_weight)
#         CAI = exp(weighted_codon_count / codon_count)

#         CAI is maximum (1) when all of the codons in the query are the most abundant synonymous
#         codons in the reference, and minimum (0-1) when all of the codons in the query are the least
#         abundant synonymous codons in the reference.

#         Parameters
#         ----------
#         query_codon_frequency_df : pandas.core.frame.DataFrame
#             This frequency table has a row per item (gene/function) to be compared to each reference
#             codon composition and a column per codon.
#         reference_codon_frequency_df : pandas.core.frame.DataFrame
#             This frequency table has a row per reference composition and a column per codon.
#         query_min_analyzed_codons : int, optional
#             A row of the query table must contain at least this number of codons with a reference
#             codon weight for CAI to be calculated. By default 0.
#         reference_min_analyzed_codons : int, optional
#             A reference must contain at least this number of codons that are also present among the
#             queries for it to be used in calculating CUB. By default 0.

#         Returns
#         -------
#         pandas.core.frame.DataFrame
#             This CUB table has the same row index as the input query table and a column per input
#             reference. With a single reference codon composition, this table has a single column.
#         """
#         if len(reference_codon_frequency_df) == 0:
#             return pd.DataFrame()

#         # Only consider codons present in both query and reference datasets.
#         nonzero_query_codon_frequency_df = query_codon_frequency_df.loc[
#             :, query_codon_frequency_df.sum() > 0]
#         nonzero_reference_codon_frequency_df = reference_codon_frequency_df.loc[
#             :, reference_codon_frequency_df.sum() > 0]
#         shared_codons = nonzero_query_codon_frequency_df.columns.intersection(
#             nonzero_reference_codon_frequency_df.columns)
#         query_codon_frequency_df = nonzero_query_codon_frequency_df[shared_codons]
#         reference_codon_frequency_df = nonzero_reference_codon_frequency_df[shared_codons]
#         # Remove queries and references that do not have enough codons remaining. When this method
#         # is called from `self.get_codon_usage_bias`, a defined single-row reference will have
#         # already been filtered, rendering this reference filter unneeded; on the other hand,
#         # omnibias references will be the same as the query references, and
#         # `query_min_analyzed_codons` was ensured to be the same as `reference_min_analyzed_codons`,
#         # so the same queries and references are here removed for not meeting the same threshold.
#         query_codon_frequency_df = query_codon_frequency_df[
#             query_codon_frequency_df.sum(axis=1) >= query_min_analyzed_codons]
#         reference_codon_frequency_df = reference_codon_frequency_df[
#             reference_codon_frequency_df.sum(axis=1) >= reference_min_analyzed_codons]

#         # Calculate reference codon weights.
#         np.seterr(divide='ignore')
#         synonymous_weight_dfs = []
#         for codons in self.synonymous_nonstop_amino_acid_codons_dict.values():
#             # Not every codon in the synonymous set need be present in the data.
#             intermediate_df = reference_codon_frequency_df[
#                 reference_codon_frequency_df.columns.intersection(codons)]
#             synonymous_weight_df = np.log(
#                 intermediate_df.div(intermediate_df.max(axis=1), axis=0))
#             synonymous_weight_dfs.append(synonymous_weight_df)
#         np.seterr(divide='warn')
#         weight_df = pd.concat(synonymous_weight_dfs, axis=1)
#         # In omnibias mode, codons may be missing from one but not all references. Zero frequencies
#         # result in weights of -∞, which are replaced by NaN.
#         weight_df = weight_df.replace(-np.inf, np.nan)
#         weight_df = weight_df[sorted(weight_df.columns)]
#         weight_array = weight_df.values
#         bool_weight_array = weight_df.notna().values

#         performed_comparisons = 0
#         total_comparisons = len(query_codon_frequency_df) * len(weight_df)
#         self.progress.new("Calculating CAI")

#         cai_rows = []
#         for query_codon_frequencies in query_codon_frequency_df.values:
#             cai_row = []
#             for weights, bool_weights in zip(weight_array, bool_weight_array):
#                 if performed_comparisons % 10000 == 0:
#                     self.progress.update(
#                         f"{performed_comparisons} / {total_comparisons} comparisons")

#                 query_codon_frequency_with_reference = query_codon_frequencies[bool_weights].sum()
#                 if query_codon_frequency_with_reference < query_min_analyzed_codons:
#                     cai_row.append(np.nan)
#                     performed_comparisons += 1
#                     continue
#                 weighted_codon_count = np.nansum(query_codon_frequencies * weights)
#                 cai_row.append(np.exp(weighted_codon_count / query_codon_frequency_with_reference))
#                 performed_comparisons += 1
#             cai_rows.append(cai_row)
#         cai_df = pd.DataFrame(cai_rows, index=query_codon_frequency_df.index)
#         cai_df = cai_df.reindex(query_codon_frequency_df.index)

#         self.progress.end()

#         return cai_df


#     def _get_delta_table(self,
#                          query_codon_frequency_df,
#                          reference_codon_frequency_df,
#                          nonreference_codon_frequency_df,
#                          query_min_analyzed_codons=default_query_min_analyzed_codons,
#                          reference_min_analyzed_codons=default_reference_min_analyzed_codons):
#         """
#         Get a table of 𝛿 (Ran and Higgs, 2012, Eq. 6) values for each query x reference comparison.

#         Calculation of 𝛿:
#         codon_weight = ln(reference_codon_synonymous_relative_frequency /
#                           nonreference_codon_synonymous_relative_frequency)
#         weighted_codon_count = ∑(codon_frequency * reference_codon_weight)
#         𝛿 = weighted_codon_count / codon_count

#         𝛿 is maximum (1) when all of the codons in the query are the most abundant synonymous codons
#         in the reference, and minimum (0-1) when all of the codons in the query are the least
#         abundant synonymous codons in the reference.

#         𝛿 codon weights are ratios of reference to genome-wide synonymous relative frequencies. 𝛿
#         ranges from (-∞, +∞), with more positive values being more similar to the reference. When
#         query and reference sets are the same, 𝛿 is a likelihood ratio evaluating the distinctness
#         of the set from the genome as a whole.

#         Parameters
#         ----------
#         query_codon_frequency_df : pandas.core.frame.DataFrame
#             This frequency table has a row per item (gene/function) to be compared to each reference
#             codon composition and a column per codon.
#         reference_codon_frequency_df : pandas.core.frame.DataFrame
#             This frequency table has a row per reference composition and a column per codon.
#         nonreference_codon_frequency_df : pandas.core.frame.DataFrame
#             This frequency table has a single row of summed frequencies across all non-reference
#             genes in the genome, with a column per codon.
#         query_min_analyzed_codons : int, optional
#             A row of the query table must contain at least this number of codons with a reference
#             codon weight for 𝛿 to be calculated. By default 0.
#         reference_min_analyzed_codons : int, optional
#             A reference must contain at least this number of codons that are also present among the
#             queries for it to be used in calculating CUB. By default 0.

#         Returns
#         -------
#         pandas.core.frame.DataFrame
#             This CUB table has the same row index as the input query table and a column per input
#             reference. With a single reference codon composition, this table has a single column of
#             data.
#         """
#         if len(reference_codon_frequency_df) == 0:
#             return pd.DataFrame()

#         # Calculate codon weights from the reference and non-reference gene sets.
#         # First get synonymous relative frequencies from each. Do this before dropping any codons
#         # missing entirely from either queries or references to ensure the accuracy of synonymous
#         # relative frequencies. For example, if CysTGT is missing from the queries but not the
#         # references, then the synonymous relative frequency of CysTGC in the references should
#         # still be calculated as TGC / (TGC + TGT).
#         reference_synonymous_codon_rel_frequency_df = \
#             self._get_synonymous_codon_rel_frequency_table(reference_codon_frequency_df)
#         nonreference_synonymous_codon_rel_frequency_df = \
#             self._get_synonymous_codon_rel_frequency_table(nonreference_codon_frequency_df)

#         # Only consider codons present in both query and reference datasets.
#         nonzero_query_codon_frequency_df = query_codon_frequency_df.loc[
#             :, query_codon_frequency_df.sum() > 0]
#         nonzero_reference_codon_frequency_df = reference_codon_frequency_df.loc[
#             :, reference_codon_frequency_df.sum() > 0]
#         shared_codons = nonzero_query_codon_frequency_df.columns.intersection(
#             nonzero_reference_codon_frequency_df.columns)
#         query_codon_frequency_df = nonzero_query_codon_frequency_df[shared_codons]
#         reference_codon_frequency_df = nonzero_reference_codon_frequency_df[shared_codons]
#         reference_synonymous_codon_rel_frequency_df = \
#             reference_synonymous_codon_rel_frequency_df[shared_codons]
#         nonreference_synonymous_codon_rel_frequency_df = \
#             nonreference_synonymous_codon_rel_frequency_df[shared_codons]
#         nonreference_synonymous_codon_rel_frequencies = \
#             nonreference_synonymous_codon_rel_frequency_df.squeeze().values
#         # Remove queries and references that do not have enough codons remaining. When this method
#         # is called from `self.get_codon_usage_bias`, a defined single-row reference will have
#         # already been filtered, rendering this reference filter unneeded; on the other hand,
#         # omnibias references will be the same as the query references, and
#         # `query_min_analyzed_codons` was ensured to be the same as `reference_min_analyzed_codons`,
#         # so the same queries and references are here removed for not meeting the same threshold.
#         query_codon_frequency_df = query_codon_frequency_df[
#             query_codon_frequency_df.sum(axis=1) >= query_min_analyzed_codons]
#         reference_codon_frequency_df = reference_codon_frequency_df[
#             reference_codon_frequency_df.sum(axis=1) >= reference_min_analyzed_codons]
#         reference_synonymous_codon_rel_frequency_df = \
#             reference_synonymous_codon_rel_frequency_df.loc[reference_codon_frequency_df.index]

#         np.seterr(divide='ignore')
#         weight_df = np.log(reference_synonymous_codon_rel_frequency_df.div(
#             nonreference_synonymous_codon_rel_frequencies, axis=1))
#         np.seterr(divide='warn')
#         # In omnibias mode, codons may be missing from one but not all references. Zero frequencies
#         # result in weights of -∞, which are replaced by NaN.
#         weight_df = weight_df.replace(-np.inf, np.nan)
#         weight_df = weight_df[sorted(weight_df.columns)]
#         weight_array = weight_df.values
#         bool_weight_array = weight_df.notna().values

#         performed_comparisons = 0
#         total_comparisons = len(query_codon_frequency_df) * len(weight_df)
#         self.progress.new("Calculating 𝛿")

#         delta_rows = []
#         for query_codon_frequencies in query_codon_frequency_df.values:
#             delta_row = []
#             for weights, bool_weights in zip(weight_array, bool_weight_array):
#                 if performed_comparisons % 10000 == 0:
#                     self.progress.update(
#                         f"{performed_comparisons} / {total_comparisons} comparisons")

#                 query_codon_frequency_with_reference = query_codon_frequencies[bool_weights].sum()
#                 if query_codon_frequency_with_reference < query_min_analyzed_codons:
#                     delta_row.append(np.nan)
#                     performed_comparisons += 1
#                     continue
#                 weighted_codon_count = np.nansum(query_codon_frequencies * weights)
#                 delta_row.append(weighted_codon_count / query_codon_frequency_with_reference)
#                 performed_comparisons += 1
#             delta_rows.append(delta_row)
#         delta_df = pd.DataFrame(delta_rows, index=query_codon_frequency_df.index)
#         delta_df = delta_df.reindex(query_codon_frequency_df.index)

#         self.progress.end()

#         return delta_df


# class MultiGenomeCodonUsage(CodonUsage):
#     """This object processes codon usage data from multiple internal and/or external genomes."""

#     def __init__(self, args, run=run, progress=progress):
#         self.args = args
#         A = lambda x: args.__dict__[x] if x in args.__dict__ else None

#         self.internal_genomes_path = A('internal_genomes')
#         self.external_genomes_path = A('external_genomes')

#         self.use_shared_function_sources = A('shared_function_sources')

#         self.preload_genomes = A('preload-genomes') or False

#         descriptions = GenomeDescriptions(args, run=run_quiet, progress=self.progress)
#         descriptions.load_genomes_descriptions(init=False)

#         if self.function_sources is None:
#             pass
#         elif (self.function_sources == [] and
#             self.use_shared_function_sources and
#             descriptions.function_annotation_sources_some_genomes_miss):
#             # All function sources are requested (in `anvi-get-codon-usage` this corresponds to
#             # `--function-sources` used as a flag), but some sources are not common to all genomes
#             # and only common sources are allowed (`use_shared_function_sources`). Report sources in
#             # common that will be analyzed and those not in common that will be ignored.
#             self.function_sources = list(descriptions.function_annotation_sources)
#             self.run.info("Function sources shared across genomes: ",
#                           ', '.join(self.function_sources))
#             self.run.info("Function sources missing in one or more genomes: ",
#                           ', '.join(descriptions.function_annotation_sources_some_genomes_miss))
#         elif (len(self.function_sources) and
#               self.use_shared_function_sources and
#               descriptions.function_annotation_sources_some_genomes_miss):
#             # A list of function sources is requested, but some are not common to all genomes and
#             # only common sources are allowed (`use_shared_function_sources`). This causes an error.
#             unshared_function_sources = []
#             for source in self.function_sources:
#                 if source not in descriptions.function_annotation_sources_some_genomes_miss:
#                     unshared_function_sources.append(source)
#             if unshared_function_sources:
#                 raise ConfigError(
#                     "Of the requested function annotation sources, the following were not run on "
#                     f"every genome: {', '.join(unshared_function_sources)}.")

#         # Store information on the genomes.
#         self.genome_info_dict = {}
#         for genome_name, genome_dict in descriptions.internal_genomes_dict.items():
#             self.genome_info_dict[genome_name] = genome_info = {}
#             genome_info['contigs_db'] = genome_dict['contigs_db_path']
#             genome_info['profile_db'] = genome_dict['profile_db_path']
#             genome_info['collection_name'] = genome_dict['collection_id']
#             genome_info['bin_id'] = genome_dict['bin_id']
#             genome_info['function_sources'] = self.function_sources
#             genome_info['codon_to_amino_acid'] = self.args.codon_to_amino_acid
#         for genome_name, genome_dict in descriptions.external_genomes_dict.items():
#             self.genome_info_dict[genome_name] = genome_info = {}
#             genome_info['contigs_db'] = genome_dict['contigs_db_path']
#             genome_info['function_sources'] = self.function_sources
#             genome_info['codon_to_amino_acid'] = self.args.codon_to_amino_acid

#         # There are memory- and CPU-efficient ways of setting up the object. `preload_genomes` loads
#         # all of the SingleGenomeCodonUsage objects into memory, allowing methods to save time by
#         # immediately accessing their frequency tables. When this argument is False,
#         # SingleGenomeCodonUsage objects are loaded into memory as needed when each genome is
#         # processed.
#         self.genome_codon_usage_dict = None
#         if self.preload_genomes:
#             self.genome_codon_usage_dict = {}
#             for genome_name, genome_info in self.genome_info_dict.items():
#                 self.genome_codon_usage_dict[genome_name] = SingleGenomeCodonUsage(
#                     argparse.Namespace(**genome_info))


#     def get_frequencies(self,
#                         from_function_sources=False,
#                         return_functions=False,
#                         return_amino_acids=False,
#                         function_accessions=None,
#                         function_names=None,
#                         expect_functions=False,
#                         relative=False,
#                         synonymous=False,
#                         sum_genes=False,
#                         average_genes=False,
#                         gene_min_codons=0,
#                         function_min_codons=0,
#                         min_codon_filter='both',
#                         drop_amino_acids=None,
#                         sequence_min_amino_acids=0,
#                         pansequence_min_amino_acids=(0, 1.0),
#                         label_amino_acids=False,
#                         infinity_to_zero=False):
#         """
#         Get absolute (default) or relative codon or amino acid frequencies from genes or functions
#         in one or more genomes.

#         See the `SingleGenomeCodonUsage.get_frequencies` docstring for descriptions of each
#         parameter.
#         """
#         kwargs = {}
#         arg_info = inspect.getargvalues(inspect.currentframe())
#         for param in arg_info.args:
#             if param == 'self':
#                 continue
#             kwargs[param] = arg_info.locals[param]
#         frequency_table_generator = self._get_genome_frequency_table(kwargs)
#         frequency_dfs = []
#         for genome_name in self.genome_info_dict:
#             self.run.info("Genome", genome_name)
#             frequency_df = next(frequency_table_generator)

#             frequency_df.insert(0, 'genome', genome_name)
#             new_index_columns = ['genome'] + frequency_df.index.names
#             frequency_df = frequency_df.reset_index().set_index(new_index_columns)

#             frequency_dfs.append(frequency_df)
#         return pd.concat(frequency_dfs)


#     def _get_genome_frequency_table(self, kwargs):
#         """This generator yields a frequency table from each genome."""
#         for genome_name, genome_info in self.genome_info_dict.items():
#             if self.preload_genomes:
#                 genome_codon_usage = self.genome_codon_usage_dict[genome_name]
#             else:
#                 genome_codon_usage = SingleGenomeCodonUsage(
#                     argparse.Namespace(**genome_info))
#             frequency_df = genome_codon_usage.get_frequencies(**kwargs)
#             print()

#             yield frequency_df


#     def get_codon_usage_bias(
#         self,
#         metrics=None,
#         from_function_sources=None,
#         function_accessions=None,
#         function_names=None,
#         expect_functions=False,
#         omnibias=False,
#         reference_function_accessions=None,
#         reference_function_names=None,
#         expect_reference_functions=False,
#         reference_gene_caller_ids=None,
#         gene_min_codons=0,
#         function_min_codons=0,
#         min_codon_filter='both',
#         drop_amino_acids=None,
#         sequence_min_amino_acids=0,
#         pansequence_min_amino_acids=(0, 1.0),
#         query_min_analyzed_codons=default_query_min_analyzed_codons,
#         reference_exclude_amino_acid_count=default_reference_exclude_amino_acid_count,
#         reference_min_analyzed_codons=default_reference_min_analyzed_codons):
#         """This generator yields a genome name and CUB table dict from each genome.

#         See the `SingleGenomeCodonUsage.get_codon_usage_bias` docstring for descriptions of each
#         parameter.
#         """
#         # Rather than individually listing a slew of arguments when calling `get_codon_usage_bias`,
#         # package them in a tidy kwargs dictionary.
#         kwargs = {}
#         arg_info = inspect.getargvalues(inspect.currentframe())
#         for param in arg_info.args:
#             if param == 'self':
#                 continue
#             kwargs[param] = arg_info.locals[param]

#         for genome_name, genome_info in self.genome_info_dict.items():
#             if self.preload_genomes:
#                 genome_codon_usage = self.genome_codon_usage_dict[genome_name]
#             else:
#                 genome_codon_usage = SingleGenomeCodonUsage(argparse.Namespace(**genome_info))
#             cub_table_dict = genome_codon_usage.get_codon_usage_bias(**kwargs)
#             print()

#             yield genome_name, cub_table_dict


# def get_custom_encodings(encodings_txt):
#     """Modify the standard genetic code given user input, as used in the scripts,
#     `anvi-get-codon-frequencies` and `anvi-get-codon-usage-bias`."""
#     codon_amino_acid_dict = copy.deepcopy(standard_code)

#     if not encodings_txt:
#         return codon_amino_acid_dict

#     encodings_df = pd.read_csv(encodings_txt, sep='\t', header=None)
#     codon_amino_acid_dict.update(dict(zip(encodings_df.iloc[:, 0], encodings_df.iloc[:, 1])))

#     return codon_amino_acid_dict


# def check_genetic_code(codon_amino_acid_dict):
#     """Check that known codons and three-letter amino acid codes are used in a dict defining the
#     genetic code, throwing an exception if needed."""
#     unrecognized_codons = []
#     unrecognized_amino_acids = []
#     for codon, amino_acid in codon_amino_acid_dict.items():
#         if codon not in constants.codons:
#             unrecognized_codons.append(codon)
#         if amino_acid not in constants.amino_acids:
#             unrecognized_amino_acids.append(amino_acid)

#     if unrecognized_codons:
#         unrecognized_codon_message = (
#             "The following codons in the provided genetic code are not recognized: "
#             f"{', '.join(unrecognized_codons)}.")
#     else:
#         unrecognized_codon_message = ""

#     if unrecognized_amino_acids:
#         unrecognized_amino_acid_message = (
#             "The following amino acids in the provided genetic code are not recognized: "
#             f"{', '.join(unrecognized_amino_acids)}. These should be three-letter codes "
#             "and \"STP\" for stop codons.")
#         if unrecognized_codon_message:
#             unrecognized_amino_acid_message = " " + unrecognized_amino_acid_message
#     else:
#         unrecognized_amino_acid_message = ""

#     if unrecognized_codon_message or unrecognized_amino_acid_message:
#         raise ConfigError(f"{unrecognized_codon_message}{unrecognized_amino_acid_message}")
