# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module for codon usage analyses at the levels of genes, groups of genes, and genomes"""


import argparse
import pandas as pd

from functools import partial
from collections import Counter, defaultdict

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.dbops import ContigsSuperclass
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

default_translation_table = constants.AA_to_codons
default_nonstop_translation_table = {
    aa: codons for aa, codons in constants.AA_to_codons.items() if aa != 'STP'}

# CAI is the Codon Adaptation Index of Sharp and Li (1987).
# Delta is the likelihood ratio of Ran and Higgs (2012).
codon_usage_bias_metrics = ['cai', 'delta']


class SingleGenomeCodonUsage(object):
    """
    This object processes codon usage data for a single genome.
    """

    def __init__(self, args=None, run=run, progress=progress, skip_init=False):
        """
        Initialize from an internal or external genome, which must have genes called.

        Raises
        ------
        ConfigError
            incomplete arguments to find internal genome
        """

        self.args = args
        self.run = run
        self.progress = progress

        # Get attributes from `args`.
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None

        self.contigs_db_path = A('contigs_db')
        self.gene_caller_ids_of_interest = list(set(A('gene_caller_ids')))

        # The following arguments are needed to a load a bin.
        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')

        num_bin_args = self.profile_db_path != None + self.collection_name != None + self.bin_id != None
        if num_bin_args == 0:
            self.is_internal_genome = False
        elif num_bin_args == 3:
            self.is_internal_genome = True
        else:
            raise ConfigError("An incomplete set of arguments was supplied to process an internal "
                              "genome. Please provide a profile database, collection, and bin. "
                              "However, if the contigs database represents a genome, provide none "
                              "of these.")

        self.function_annotation_sources = A('annotation_sources')

        self.contig_sequences_dict = None
        self.genes_in_contigs_dict = None
        self.gene_function_calls_dict = None
        self.gene_codon_frequency_df = None
        self.function_codon_frequency_dict = None
        if not skip_init:
            self._init()
            self._make_gene_codon_frequency_table()
            self._make_function_codon_frequency_tables()


    def _init(self):
        """
        Load gene data from the contigs database.

        Raises
        ------
        ConfigError
            missing some requested gene callers IDs
        """

        utils.is_contigs_db(self.contigs_db_path)

        if self.is_internal_genome:
            args = argparse.Namespace()
            args.profile_db_path = self.profile_db_path
            args.collection_name = self.collection_name
            args.bin_id = self.bin_id
            # Initialize the contigs superclass from the splits of the internal genome bin.
            args.split_names_of_interest = ccollections.GetSplitNamesInBins(args)
        contigs_super = ContigsSuperclass(args, r=run_quiet)
        self.contig_sequences_dict = contigs_super.contig_sequences

        self.genes_in_contigs_dict = contigs_super.genes_in_contigs_dict
        if self.gene_caller_ids_of_interest:
            missing_gene_caller_ids = []
            for gene_caller_id in self.gene_caller_ids_of_interest:
                if gene_caller_id not in self.gene_caller_ids_of_interest:
                    missing_gene_caller_ids.append(gene_caller_id)
            if missing_gene_caller_ids:
                raise ConfigError(f"The contigs database, {contigs_super.a_meta['project_name']}, "
                                  "does not contain the following requested gene callers IDs: "
                                  f"{', '.join(missing_gene_caller_ids)}")

        # gene_function_calls_dict[gene_caller_id][source] =
        # (function accession, function name, annotation evalue)
        if self.function_annotation_sources:
            contigs_super.init_functions(requested_sources=self.function_annotation_sources)
            self.gene_function_calls_dict = contigs_super.gene_function_calls_dict
        else:
            self.gene_function_calls_dict = {}


    def _make_gene_codon_frequency_table(self):
        """
        Generates the per-gene codon frequency DataFrame as `self.gene_codon_frequency_df`.
        """

        if self.genes_in_contigs_dict is None:
            self.run.warning("The gene codon frequency table was not generated because the "
                             f"{type(self).__name__} object's `init` method should be run.")
            return

        gene_codon_frequencies = []
        skipped_noncoding_gene_caller_ids = []
        for gene_caller_id, gene_call in self.genes_in_contigs_dict.items():
            # `gene_call` is a dictionary.
            if gene_call['call_type'] != constants.gene_call_types['CODING']:
                skipped_noncoding_gene_caller_ids.add(gene_caller_id)
                continue

            gene_codon_frequencies.append(Counter(utils.get_list_of_codons_for_gene_call(
                gene_call, self.contig_sequences_dict)))

        gene_codon_frequency_df = pd.DataFrame.from_records(gene_codon_frequencies)

        observed_codons = gene_codon_frequency_df.columns.tolist()
        for codon in constants.codon_to_AA:
            if codon not in observed_codons:
                gene_codon_frequency_df[codon] = 0

        gene_codon_frequency_df = gene_codon_frequency_df.fillna(0)
        gene_codon_frequency_df = gene_codon_frequency_df[sorted(gene_codon_frequency_df.columns)]
        gene_codon_frequency_df.index = self.genes_in_contigs_dict
        self.gene_codon_frequency_df = gene_codon_frequency_df

        if skipped_noncoding_gene_caller_ids:
            self.run.warning(f"{len(skipped_noncoding_gene_caller_ids)} of "
                             f"{len(self.genes_in_contigs_dict)} genes were non-coding and "
                             "not added to the codon frequency table.")


    def _make_function_codon_frequency_tables(self):
        """
        Generates per-function codon frequency DataFrames in `self.function_codon_frequency_dict`.

        This dictionary is keyed by function annotation source.
        """

        if not self.function_annotation_sources:
            self.run.warning("The function codon frequency table was not generated because no "
                             "function annotation sources were requested. These would be stored "
                             f"in the {type(self).__name__} object's `function_annotation_sources` "
                             "list attribute, which is set at object initialization by "
                             "`args.annotation_sources`.")
            return

        if self.gene_codon_frequency_df is None:
            self.run.warning("The function codon frequency table was not generated because the "
                             f"{type(self).__name__} object's `make_gene_codon_frequency_table` "
                             "method should be run.")
            return

        function_genes_dict = {function_source: defaultdict(list) for function_source in
                               self.function_annotation_sources}
        for gene_caller_id, function_source_dict in self.gene_function_calls_dict.items():
            for function_source, annotation_info in function_source_dict.items():
                if function_source in self.function_annotation_sources:
                    # `annotation_info` is (function accession, function name, annotation evalue)
                    function_genes_dict[function_source][
                        (annotation_info[0], annotation_info[1])].append(gene_caller_id)

        function_codon_frequency_dict = {}
        for function_source, function_source_dict in function_genes_dict.items():
            function_codon_frequency_rows = []
            annotated_gene_caller_ids = set()
            for function_id, gene_caller_ids in function_source_dict.items():
                function_frequencies = list(self.gene_codon_frequency_df.loc[gene_caller_ids].sum())
                function_codon_frequency_rows.append(list(function_id) + function_frequencies)
                annotated_gene_caller_ids.update(gene_caller_ids)

            function_codon_frequency_df = pd.DataFrame(
                function_codon_frequency_rows,
                columns=['accession', 'function'] + list(constants.codon_to_AA))
            function_codon_frequency_df = function_codon_frequency_df.set_index(
                ['accession', 'function'])
            function_codon_frequency_dict[function_source] = function_codon_frequency_df

            self.run.info(f"Genes annotated by {function_source}", len(annotated_gene_caller_ids))
        self.function_codon_frequency_dict = function_codon_frequency_dict


    def get_codon_frequencies(self,
                              as_relative_frequency=False,
                              per_amino_acid=False,
                              from_function_sources=False,
                              as_source_dict=True,
                              collapse_all=False,
                              min_length=0,
                              ignore_stop_codons=0,
                              header_amino_acid=False):
        """
        Get codon frequencies from genes or functions.

        Parameters
        ----------
        as_relative_frequency : bool, optional
            Get relative codon frequencies, by default False
        per_amino_acid : bool, optional
            Get relative codon frequencies among codons decoding each amino acid, by default False
        from_function_sources : bool, str, or list of str, optional
            Get codon frequencies from functions (one or more genes annoted by a functional
            annotation source), with True returning all functional annotation sources, a string
            returning the one source provided, and a list returning the multiple sources provided,
            by default False
        as_source_dict : bool, optional
            When codon frequencies are returned from multiple functional annotation sources, it can
            either be as a dict keyed by source (True) or as a table with an additional column
            indicating source (False), by default True
        collapse_all : bool, optional
            Sum all genes or functions in a source into a single frequency, by default False
        min_length : int, optional
            Minimum number of codons for gene or function to be reported, by default 0
        ignore_stop_codons : bool, optional
            Ignore stop codons in calculations and do not report them, by default False
        header_amino_acid : bool, optional
            Include amino acid in column header (e.g., LysAAA rather than AAA), by default False

        Returns
        -------
        pandas.core.frame.DataFrame or dict of DataFrames
            Codon table or dictionary of tables by functional annotation source

        Raises
        ------
        ConfigError
            Unrecognized functional annotation sources requested
        ConfigError
            End of method reached
        """

        if not as_relative_frequency and per_amino_acid:
            as_relative_frequency = True

        # Per-gene results:
        ### Absolute frequencies
        if (not as_relative_frequency and
            not per_amino_acid and
            not from_function_sources and
            not collapse_all):
            return self._get_gene_codon_frequency_table(min_length=min_length,
                                                        ignore_stop_codons=ignore_stop_codons,
                                                        header_amino_acid=header_amino_acid)

        per_gene_get_table = lambda method: partial(method,
                                                    self._get_gene_codon_frequency_table(),
                                                    min_length=min_length,
                                                    ignore_stop_codons=ignore_stop_codons,
                                                    header_amino_acid=header_amino_acid)

        ### Relative frequencies
        if (as_relative_frequency and
            not per_amino_acid and
            not from_function_sources and
            not collapse_all):
            return per_gene_get_table(self._get_codon_rel_frequency_table)
        ### Per-amino acid relative frequencies
        if (as_relative_frequency and
            per_amino_acid and
            not from_function_sources and
            not collapse_all):
            return per_gene_get_table(self._get_per_amino_acid_codon_rel_frequency_table)

        # Results collapsed across genes:
        ### Absolute frequencies
        if (not as_relative_frequency and
            not per_amino_acid and
            not from_function_sources and
            collapse_all):
            return per_gene_get_table(self._get_collapsed_codon_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not per_amino_acid and
            not from_function_sources and
            collapse_all):
            return per_gene_get_table(self._get_collapsed_codon_rel_frequency_table)
        ### Per-amino acid relative frequencies
        if (as_relative_frequency and
            per_amino_acid and
            not from_function_sources and
            collapse_all):
            return per_gene_get_table(self._get_collapsed_per_amino_acid_codon_rel_frequency_table)

        # Per-function results:
        if from_function_sources == False:
            raise ConfigError("This point should not be reached in `get_codon_frequencies`. "
                              "`from_function_sources` is False, but all possible avenues of "
                              "returning per-gene results appear to have been exhausted. Please "
                              "contact the developers to investigate further.")
        elif from_function_sources == True:
            function_sources = list(self.function_codon_frequency_dict)
        else:
            function_sources = list(from_function_sources)

        unrecognized_function_sources = []
        for function_source in function_sources:
            if function_source not in self.function_codon_frequency_dict:
                unrecognized_function_sources.append(function_source)
        if unrecognized_function_sources:
            raise ConfigError("The requested function annotation sources, "
                              f"{', '.join(function_sources)}, are not provided for the genome. "
                              "The following sources are available: "
                              f"{', '.join(self.function_codon_frequency_dict)}.")

        per_function_get_table = partial(self._get_function_codon_table,
                                         min_length=min_length,
                                         ignore_stop_codons=ignore_stop_codons,
                                         header_amino_acid=header_amino_acid)
        per_function_get_dict = partial(self._get_function_codon_dict,
                                        min_length=min_length,
                                        ignore_stop_codons=ignore_stop_codons,
                                        header_amino_acid=header_amino_acid)

        ### Absolute frequencies, returning table
        if (not as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            not as_source_dict and
            not collapse_all):
            return per_function_get_table(self._get_codon_frequency_table)
        ### Absolute frequencies, returning dict
        if (not as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            as_source_dict and
            not collapse_all):
            return per_function_get_dict(self._get_codon_frequency_table)
        ### Relative frequencies, returning table
        if (as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            not as_source_dict and
            not collapse_all):
            return per_function_get_table(self._get_codon_rel_frequency_table)
        ### Relative frequencies, returning dict
        if (as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            as_source_dict and
            not collapse_all):
            return per_function_get_dict(self._get_codon_rel_frequency_table)
        ### Per-amino acid relative frequencies, returning table
        if (as_relative_frequency and
            per_amino_acid and
            from_function_sources and
            not as_source_dict and
            not collapse_all):
            return per_function_get_table(self._get_per_amino_acid_codon_rel_frequency_table)
        ### Per-amino acid relative frequencies, returning dict
        if (as_relative_frequency and
            per_amino_acid and
            from_function_sources and
            as_source_dict and
            not collapse_all):
            return per_function_get_dict(self._get_per_amino_acid_codon_rel_frequency_table)

        # Results collapsed across functions in each source:
        ### Absolute frequencies, returning table
        if (not as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            not as_source_dict and
            collapse_all):
            return per_function_get_table(self._get_collapsed_codon_frequency_table)
        ### Absolute frequencies, returning dict
        if (not as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            as_source_dict and
            collapse_all):
            return per_function_get_dict(self._get_collapsed_codon_frequency_table)
        ### Relative frequencies, returning table
        if (as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            not as_source_dict and
            collapse_all):
            return per_function_get_table(self._get_collapsed_codon_rel_frequency_table)
        ### Relative frequencies, returning dict
        if (as_relative_frequency and
            not per_amino_acid and
            from_function_sources and
            as_source_dict and
            collapse_all):
            return per_function_get_dict(self._get_collapsed_codon_rel_frequency_table)
        ### Per-amino acid relative frequencies, returning table
        if (as_relative_frequency and
            per_amino_acid and
            from_function_sources and
            not as_source_dict and
            collapse_all):
            return per_function_get_table(
                self._get_collapsed_per_amino_acid_codon_rel_frequency_table)
        ### Per-amino acid relative frequencies, returning dict
        if (as_relative_frequency and
            per_amino_acid and
            from_function_sources and
            as_source_dict and
            collapse_all):
            return per_function_get_dict(
                self._get_collapsed_per_amino_acid_codon_rel_frequency_table)

        raise ConfigError("This point should not be reached at the end of the method, "
                          "`get_codon_frequencies`. Please contact the developers. Since you found "
                          "the end of the earth, you now get to hear a top secret mnemonic for the "
                          "rare earth elements. Scandalous Yiddish language centers praise Ned's "
                          "promise of small European garden tubs. Dinosaurs hobble erotically "
                          "thrumming yellow lutes. (scandium Sc, yttrium Y, lanthanum La, cerium "
                          "Ce, praseodymium Pr, neodymium Nd, promethium Pm, samarium Sm, europium "
                          "Eu, gadolinium Gd, terbium Tb, dysprosium Dy, holmium Ho, erbium Er, "
                          "thulium Tm, ytterbium Yb, lutetium Lu) Credit for the lanthanide series "
                          "mnemonic goes to Martyn Poliakoff: "
                          "https://www.youtube.com/watch?v=Q21clW0s0B8&ab_channel=PeriodicVideos")


    def get_amino_acid_frequencies(self,
                                   as_relative_frequency=False,
                                   from_function_sources=False,
                                   as_source_dict=True,
                                   collapse_all=False,
                                   min_length=0,
                                   ignore_stop_codons=0):
        """
        Get amino acid frequencies from genes or functions.

        Parameters
        ----------
        as_relative_frequency : bool, optional
            Get relative amino acid frequencies, by default False
        from_function_sources : bool, str, or list of str, optional
            Get amino acid frequencies from functions (one or more genes annoted by a functional
            annotation source), with True returning all functional annotation sources, a string
            returning the one source provided, and a list returning the multiple sources provided,
            by default False
        as_source_dict : bool, optional
            When amino acid frequencies are returned from multiple functional annotation sources, it
            can either be as a dict keyed by source (True) or as a table with an additional column
            indicating source (False), by default True
        collapse_all : bool, optional
            Sum all genes or functions in a source into a single frequency, by default False
        min_length : int, optional
            Minimum number of amino acids for gene or function to be reported, by default 0
        ignore_stop_codons : bool, optional
            Ignore stop codons in calculations and do not report them as "amino acid" STP, by
            default False

        Returns
        -------
        pandas.core.frame.DataFrame or dict of DataFrames
            Codon table or dictionary of tables by functional annotation source

        Raises
        ------
        ConfigError
            Unrecognized functional annotation sources requested
        ConfigError
            End of method reached
        """

        get_table = lambda method: partial(
            self._get_amino_acid_table,
            method(min_length=min_length, ignore_stop_codons=ignore_stop_codons))

        # Per-gene results:
        ### Absolute frequencies
        if (not as_relative_frequency and
            not from_function_sources and
            not collapse_all):
            return get_table(self._get_gene_codon_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not from_function_sources and
            not collapse_all):
            return get_table(self._get_codon_rel_frequency_table)

        # Results collapsed across genes:
        ### Absolute frequencies
        if (not as_relative_frequency and
            not from_function_sources and
            collapse_all):
            return get_table(self._get_collapsed_codon_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not from_function_sources and
            collapse_all):
            return get_table(self._get_collapsed_codon_rel_frequency_table)

        # Per-function results
        if from_function_sources == False:
            raise ConfigError("This point should not be reached in `get_amino_acid_frequencies`. "
                              "`from_function_sources` is False, but all possible avenues of "
                              "returning per-gene results appear to have been exhausted. Please "
                              "contact the developers to investigate further.")
        elif from_function_sources == True:
            function_sources = list(self.function_codon_frequency_dict)
        else:
            function_sources = list(from_function_sources)

        unrecognized_function_sources = []
        for function_source in function_sources:
            if function_source not in self.function_codon_frequency_dict:
                unrecognized_function_sources.append(function_source)
        if unrecognized_function_sources:
            raise ConfigError("The requested function annotation sources, "
                              f"{', '.join(function_sources)}, are not provided for the genome. "
                              "The following sources are available: "
                              f"{', '.join(self.function_codon_frequency_dict)}.")

        get_function_table = lambda method: partial(
            self._get_amino_acid_table,
            self._get_function_codon_table(
                method,
                function_sources=function_sources,
                min_length=min_length,
                ignore_stop_codons=ignore_stop_codons))
        get_function_dict = lambda method: partial(
            self._get_amino_acid_dict,
            self._get_function_codon_dict(
                method,
                function_sources=function_sources,
                min_length=min_length,
                ignore_stop_codons=ignore_stop_codons))

        ### Absolute frequencies, returning table
        if (not as_relative_frequency and
            from_function_sources and
            not as_source_dict and
            not collapse_all):
            return get_function_table(self._get_codon_frequency_table)
        ### Absolute frequencies, returning dict
        if (not as_relative_frequency and
            from_function_sources and
            as_source_dict and
            not collapse_all):
            return get_function_dict(self._get_codon_frequency_table)
        ### Relative frequencies, returning table
        if (as_relative_frequency and
            from_function_sources and
            not as_source_dict and
            not collapse_all):
            return get_function_table(self._get_codon_rel_frequency_table)
        ### Relative frequencies, returning dict
        if (as_relative_frequency and
            from_function_sources and
            as_source_dict and
            not collapse_all):
            return get_function_dict(self._get_codon_rel_frequency_table)

        # Results collapsed across functions in each source:
        ### Absolute frequencies, returning table
        if (not as_relative_frequency and
            from_function_sources and
            not as_source_dict and
            collapse_all):
            return get_function_table(self._get_collapsed_codon_frequency_table)
        ### Absolute frequencies, returning dict
        if (not as_relative_frequency and
            from_function_sources and
            as_source_dict and
            collapse_all):
            return get_function_dict(self._get_collapsed_codon_frequency_table)
        ### Relative frequencies, returning table
        if (as_relative_frequency and
            from_function_sources and
            not as_source_dict and
            collapse_all):
            return get_function_table(self._get_collapsed_codon_rel_frequency_table)
        ### Relative frequencies, returning dict
        if (as_relative_frequency and
            from_function_sources and
            as_source_dict and
            collapse_all):
            return get_function_table(self._get_collapsed_codon_rel_frequency_table)

        raise ConfigError("This point should not be reached at the end of the method, "
                          "`get_amino_acid_frequencies`. Please contact the developers. Since you "
                          "found the end of the earth, you now get to hear a top secret mnemonic "
                          "for the rare earth elements. Scandalous Yiddish language centers praise "
                          "Ned's promise of small European garden tubs. Dinosaurs hobble "
                          "erotically thrumming yellow lutes. (scandium Sc, yttrium Y, lanthanum "
                          "La, cerium Ce, praseodymium Pr, neodymium Nd, promethium Pm, samarium "
                          "Sm, europium Eu, gadolinium Gd, terbium Tb, dysprosium Dy, holmium Ho, "
                          "erbium Er, thulium Tm, ytterbium Yb, lutetium Lu) Credit for the "
                          "lanthanide series mnemonic goes to Martyn Poliakoff: "
                          "https://www.youtube.com/watch?v=Q21clW0s0B8&ab_channel=PeriodicVideos")


    def _filter_min_length(get_table):
        """
        Decorator to discard rows with fewer than minimum number of codons/amino acids.

        Parameters
        ----------
        get_table : function
            turns a table of row x codon or aa data into another table of row x codon or aa data
        """

        def wrapper(*args, **kwargs):
            frequency_df = args[0]
            if kwargs['min_length']:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= kwargs['min_length']]
            return get_table(frequency_df, *args[1: ], **kwargs)

        return wrapper


    def _ignore_stop_codons(get_table):
        """
        Decorator to discard stop codon columns from input codon frequency table.

        Parameters
        ----------
        get_table : function
            turns a table of row x codon or aa data into another table of row x codon or aa data
        """

        def wrapper(*args, **kwargs):
            codon_frequency_df = args[0]
            if kwargs['ignore_stop_codons']:
                codon_frequency_df = codon_frequency_df.drop(
                    default_translation_table['STP'], axis=1, errors='ignore')
            return get_table(codon_frequency_df, *args[1: ], **kwargs)

        return wrapper


    def _add_amino_acid_to_header(get_codon_table):
        """
        Decorator to add amino acid to codon column header.

        Parameters
        ----------
        get_codon_table : function
            turns a table of row x codon data into another table of row x codon data

        Raises
        ------
        KeyError
            columns are not recognized as codons
        """

        def wrapper(*args, **kwargs):
            codon_df = get_codon_table(*args, **kwargs)
            if kwargs['header_amino_acid']:
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


    @_filter_min_length
    @_ignore_stop_codons
    @_add_amino_acid_to_header
    def _get_gene_codon_frequency_table(self, **kwargs):
        return self.gene_codon_frequency_df


    @_filter_min_length
    @_ignore_stop_codons
    @_add_amino_acid_to_header
    def _get_codon_frequency_table(self, codon_frequency_df, **kwargs):
        return codon_frequency_df


    @_filter_min_length
    @_ignore_stop_codons
    @_add_amino_acid_to_header
    def _get_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        codon_relative_frequency_df = codon_frequency_df.div(codon_frequency_df.sum(axis=1))
        return codon_relative_frequency_df


    @_filter_min_length
    @_ignore_stop_codons
    @_add_amino_acid_to_header
    def _get_per_amino_acid_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        if kwargs['ignore_stop_codons']:
            aa_to_codons = default_nonstop_translation_table
        else:
            aa_to_codons = default_translation_table

        per_aa_codon_rel_frequency_df = pd.DataFrame()
        for codons in aa_to_codons.values():
            aa_codon_frequency_df = codon_frequency_df[codons]
            per_aa_codon_rel_frequency_df[codons] = aa_codon_frequency_df.div(
                aa_codon_frequency_df.sum(axis=1))

        return per_aa_codon_rel_frequency_df


    @_filter_min_length
    @_ignore_stop_codons
    @_add_amino_acid_to_header
    def _get_collapsed_codon_frequency_table(self, codon_frequency_df, **kwargs):
        collapsed_codon_frequency_df = codon_frequency_df.sum().to_frame('all').T
        return collapsed_codon_frequency_df


    def _get_collapsed_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        collapsed_codon_frequency_df = \
            self._get_collapsed_codon_frequency_table(codon_frequency_df, kwargs)

        collapsed_codon_rel_frequency_df = collapsed_codon_frequency_df.div(
            collapsed_codon_frequency_df.sum(axis=1))

        return collapsed_codon_rel_frequency_df


    def _get_collapsed_per_amino_acid_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        collapsed_codon_frequency_df = self._get_collapsed_codon_frequency_table(
            codon_frequency_df, **{'min_length': kwargs['min_length'],
                                   'ignore_stop_codons': kwargs['ignore_stop_codons']})

        collapsed_per_amino_acid_codon_rel_frequency_df = \
            self._get_per_amino_acid_codon_rel_frequency_table(
                collapsed_codon_frequency_df, **{'header_amino_acid': kwargs['header_amino_acid']})

        return collapsed_per_amino_acid_codon_rel_frequency_df


    def _get_function_codon_table(self, method, function_sources, **kwargs):
        if len(function_sources) == 1:
            function_codon_frequency_df = method(
                self.function_codon_frequency_dict[function_sources[0]], **kwargs)
        else:
            source_codon_frequency_dfs = []
            for function_source in function_sources:
                source_codon_frequency_df = method(
                    self.function_codon_frequency_dict[function_source], **kwargs)
                source_codon_frequency_df = source_codon_frequency_df.reset_index()
                source_codon_frequency_df['source'] = function_source
                source_codon_frequency_dfs.append(function_codon_frequency_df)

            function_codon_frequency_df = pd.concat(source_codon_frequency_dfs, ignore_index=True)
            function_codon_frequency_df = function_codon_frequency_df.set_index(
                ['source', 'accession', 'function'])

        return function_codon_frequency_df


    def _get_function_codon_dict(self, method, function_sources, **kwargs):
        function_codon_frequency_dict = {}
        for function_source in function_sources:
            function_codon_frequency_dict[function_source] = method(
                self.function_codon_frequency_dict[function_source], **kwargs)
        return function_codon_frequency_dict


    def _get_amino_acid_table(self, codon_df):
        aa_df = pd.DataFrame()
        for amino_acid, codons in default_translation_table.items():
            try:
                aa_df[amino_acid] = codon_df[codons].sum(axis=1)
            except KeyError:
                # This occurs when `ignore_stop_codons` was used in generating the codon table.
                pass
        return aa_df


    def _get_amino_acid_dict(self, codon_dict):
        aa_dict = {}
        for function_source, function_codon_df in codon_dict:
            aa_dict[function_source] = self._get_amino_acid_table(function_codon_df)
        return aa_dict


    def get_relative_codon_frequency_table(self, function_source):
        pass


    def get_relative_codon_frequency_table_normalized_by_amino_acid(self):
        pass


class MultiGenomeCodonUsage(GenomeDescriptions):
    """This object processes codon usage data for multiple genomes."""

    def __init__(self, args, run=run, progress=progress):
        """Initialize from internal and/or external genomes, which must contain gene calls.

        Parameters
        ----------
        args : _type_
            _description_
        run : _type_, optional
            _description_, by default run
        progress : _type_, optional
            _description_, by default progress

        Raises
        ------
        ConfigError
            _description_
        ConfigError
            _description_
        """

        self.args = args
        self.run = run
        self.progress = progress

        # Function to get attributes from `args`.
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None


        # The superclass sets the attributes `gene_caller`, `input_file_for_internal_genomes`, and
        # `input_file_for_external_genomes` from the arguments.
        GenomeDescriptions.__init__(self, args, run=self.run, progress=self.progress)
        self.load_genomes_descriptions(init=False)


        # Here are the possibilities for codon usage analysis of gene functional groups:
        # 1) No analysis because no functional annotation sources are provided.
        # 2) Annotation sources are provided and must be present in all genomes.
        # 3) All available annotation sources are requested and must be present in all genomes.
        if A('annotation_sources'):
            requested_function_annotation_sources = [
                s.strip() for s in A('annotation_sources').split(',')]
            if (len(requested_function_annotation_sources) == 1
                and requested_function_annotation_sources[0].lower() == 'all'):
                requested_function_annotation_sources = ['all']
        else:
            requested_function_annotation_sources = []
            self.run.warning("No functional annotation sources were requested, so codon usage "
                             "analysis will be conducted for individual gene calls and not for "
                             "homologous genes or functional groups of genes.")
        if (len(requested_function_annotation_sources) == 1
            and requested_function_annotation_sources[0] == 'all'):
            if self.function_annotation_sources_some_genomes_miss:
                raise ConfigError("With the 'all' argument, we assume that you want to analyze "
                                  "codon usage using functional annotation sources present in ANY "
                                  "genome, and that these sources should be present in ALL "
                                  "genomes. Unfortunately, the following sources were not present "
                                  "in every genome: "
                                  f"{', '.join(requested_function_annotation_sources)}")
        else:
            missing_sources = []
            for requested_source in requested_function_annotation_sources:
                if requested_source not in self.function_annotation_sources:
                    missing_sources.append(requested_source)
            if missing_sources:
                raise ConfigError("We assume you want the requested functional annotation "
                                  "sources to be present in every genome. Unfortunately, the "
                                  "following sources were not present in every genome: "
                                  f"{', '.join(requested_function_annotation_sources)}")


        # Initialize a `SingleGenomeCodonUsage` object for each genome.
        self.genome_codon_usage_dict = {}
        for genome_name, genome_dict in self.internal_genomes_dict.items():
            genome_codon_usage_args = {}
            genome_codon_usage_args['contigs_db_path'] = genome_dict['contigs_db_path']
            genome_codon_usage_args['profile_db_path'] = genome_dict['profile_db_path']
            genome_codon_usage_args['collection_name'] = genome_dict['collection_id']
            genome_codon_usage_args['bin_id'] = genome_dict['bin_id']
            genome_codon_usage_args['annotation_sources'] = self.function_annotation_sources
            self.genome_codon_usage_dict[genome_name] = SingleGenomeCodonUsage(genome_codon_usage_args)
        for genome_name, genome_dict in self.external_genomes_dict.items():
            genome_codon_usage_args = {}
            genome_codon_usage_args['contigs_db_path'] = genome_dict['contigs_db_path']
            genome_codon_usage_args['annotation_sources'] = self.function_annotation_sources
            self.genome_codon_usage_dict[genome_name] = SingleGenomeCodonUsage(genome_codon_usage_args)
