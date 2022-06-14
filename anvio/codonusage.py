# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module for codon usage analyses at the levels of genes, groups of genes, and genomes"""


import argparse
import numpy as np
import pandas as pd

from collections import Counter

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

translation_table = constants.AA_to_codons
nonstop_translation_table = {aa: co for aa, co in constants.AA_to_codons.items() if aa != 'STP'}

# CAI is the Codon Adaptation Index of Sharp and Li (1987).
# Delta is the likelihood ratio of Ran and Higgs (2012).
cub_metrics = ['cai', 'delta']
# The following CUB metrics rely upon comparison to a set of reference genes.
reference_dependent_cub_metrics = ['cai', 'delta']

# For CUB, remove single codon amino acids, stop codons, and codons with "N" nucleotide which were
# recorded as "None".
ignored_decodings_in_cub = ['Met', 'Trp', 'STP']
ignored_codons_in_cub = [None]
for decoding in ignored_decodings_in_cub:
    ignored_codons_in_cub += constants.AA_to_codons[decoding]
synonymous_codon_dict = {aa: co for aa, co in constants.AA_to_codons.items()
                         if aa not in ignored_decodings_in_cub}


class SingleGenomeCodonUsage(object):
    """
    This object processes codon usage data for a single genome.

    Manipulate the data using the methods, `get_codon_frequencies` and `get_codon_usage_bias`.
    """

    def __init__(self, args=None, run=run, progress=progress, skip_init=False):
        """Initialize from an internal or external genome, which must have genes called."""
        self.args = args
        self.run = run
        self.progress = progress

        # Get attributes from `args`.
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None

        self.contigs_db_path = A('contigs_db')

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

        self._load_contigs_db_data()
        self._make_gene_codon_frequency_table()


    def _load_contigs_db_data(self):
        """Load gene data from the contigs database."""
        utils.is_contigs_db(self.contigs_db_path)

        if self.is_internal_genome:
            args = argparse.Namespace()
            args.profile_db_path = self.profile_db_path
            args.collection_name = self.collection_name
            args.bin_id = self.bin_id
            # Initialize the contigs superclass from the splits of the internal genome bin.
            args.split_names_of_interest = ccollections.GetSplitNamesInBins(args)
        contigs_super = ContigsSuperclass(args, r=run_quiet)
        contigs_super.init_contig_sequences()
        self.contig_sequences_dict = contigs_super.contig_sequences

        self.genes_in_contigs_dict = contigs_super.genes_in_contigs_dict
        self.gene_caller_ids = list(set(self.genes_in_contigs_dict))

        contigs_super.init_functions()
        self.function_annotation_sources = contigs_super.list_function_sources()
        gene_function_rows = []
        for gene_caller_id, annotation_dict in contigs_super.gene_function_calls_dict.items():
            for annotation_source, annotation in annotation_dict.items():
                if annotation is None:
                    continue
                if annotation_source == 'KEGG_BRITE':
                    # Include every possible depth of categorization.
                    hierarchy_accession = annotation[0]
                    categorization = annotation[1]
                    split_categories = categorization.split('>>>')
                    for depth in range(1, len(split_categories) + 1):
                        gene_function_rows.append(
                            [gene_caller_id,
                             annotation_source,
                             hierarchy_accession,
                             '>>>'.join(split_categories[: depth])])
                else:
                    gene_function_rows.append(
                        [gene_caller_id, annotation_source, annotation[0], annotation[1]])
        self.gene_function_df = pd.DataFrame(
            gene_function_rows, columns=['gene_caller_id', 'source', 'accession', 'function'])


    def _make_gene_codon_frequency_table(self):
        """Generate the per-gene codon frequency DataFrame as `self.gene_codon_frequency_df`."""
        gene_codon_frequencies = []
        skipped_noncoding_gene_caller_ids = []
        for gene_caller_id in self.gene_caller_ids:
            # `gene_call` is a dictionary.
            gene_call = self.genes_in_contigs_dict[gene_caller_id]
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
        gene_codon_frequency_df.index = self.gene_caller_ids
        self.gene_codon_frequency_df = gene_codon_frequency_df

        if skipped_noncoding_gene_caller_ids:
            self.run.warning(f"{len(skipped_noncoding_gene_caller_ids)} of "
                             f"{len(self.gene_caller_ids)} genes were non-coding and "
                             "not added to the codon frequency table.")


    def get_codon_frequencies(self,
                              from_function_sources=False,
                              return_amino_acids=False,
                              gene_caller_ids=None,
                              function_accessions=None,
                              function_names=None,
                              as_relative_frequency=False,
                              per_amino_acid=False,
                              sum_genes=False,
                              average_genes=False,
                              gene_min_codons=0,
                              function_min_codons=0,
                              min_codon_filter='both',
                              drop_amino_acids=None,
                              min_amino_acids=0,
                              min_gene_fraction=0,
                              label_amino_acids=False):
        """
        Get absolute (default) or relative codon frequencies from genes or functions.

        Relative frequencies can be normalized per-amino acid to synonymous codons.

        Parameters
        ----------
        from_function_sources : bool, str, or iterable of str, optional
            If True (default False), return codon frequencies from functions (groups of genes)
            rather than genes. A value of True returns all available functional annotation sources.
            A string value selects the specified source, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. A
            list value selects a subset of sources, e.g., ['KOfam', 'COG20_FUNCTION']. When
            'KEGG_BRITE' is used, in the absence of `function_names` narrowing the scope of the
            inquiry, all possible depths of each BRITE hierarchy in the data is returned, e.g.,
            'Ribosome>>>Ribosomal proteins' and the more general 'Ribosome' will each have rows in
            the output table.
        return_amino_acids : bool, optional
            If True (default False), output frequency table columns are decoded amino acids (plus
            STP) rather than codons. Synonymous codon frequencies are summed to produce the amino
            acid frequencies.
        gene_caller_ids : None, optional
            Genes with the given IDs are selected for analysis. This parameter can be used alongside
            `function_accessions` and `function_names`. By default None.
        function_accessions : None, optional
            Genes annotated with the given function accessions are selected for analysis. The
            argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
            parameter can be used alongside `gene_caller_ids` and `function_names`. Note that
            'KEGG_BRITE' does not have function accessions. By default None.
        function_names : None, optional
            Genes annotated with the given function names are selected for analysis. The argument
            must be a dict keyed by function annotation source (e.g., 'KOfam', 'COG20_FUNCTION',
            'Pfam') and with values being lists of function names. Unlike `function_accessions`,
            'KEGG_BRITE' may be used as a source, with names (categorizations) given to an arbitrary
            level of the hierarchy, such as ['Ribosome>>>Ribosomal proteins', 'Transcription
            machinery>>>Prokaryotic type>>>Bacterial type>>>RNA polymerase']. This parameter can be
            used alongside `gene_caller_ids` and `function_accessions`. By default None.
        as_relative_frequency : bool, optional
            Return relative rather than absolute codon frequencies if True. By default False.
        per_amino_acid : bool, optional
            Return codon relative frequencies among (synonymous) codons decoding each amino acid
            (plus stop). By default False.
        sum_genes : bool, optional
            Sum codon frequencies of genes, returning a one-row DataFrame of the summed frequencies.
            If `from_function_sources` is used, then genes are limited to those with the functional
            annotations. Relative/per-amino acid frequencies are calculated after summing absolute
            frequencies. By default False.
        average_genes : bool, optional
            Average codon frequencies of genes, returning a one-row DataFrame of the averaged
            frequencies. If `from_function_sources` is used, then genes are limited to those with
            the functional annotations. Averaging occurs after calculation of relative/per-amino
            acid frequencies. By default False.

        Additional Parameters
        ---------------------
        These filter genes and functions and filter codons. Here is the order of all possible
        filters, with functions, not just genes, being considered:
        gene codon frequency table ->
            drop genes on codon length ->
            drop functions on their codon sum from the set of genes defining the function ->
            drop codons for defined amino acids ->
            drop codons for infrequent amino acids ->
            drop genes on remaining codon frequency ->
            drop functions on their remaining codon sum ->
        filtered gene codon frequency table

        gene_min_codons : int, optional
            Ignore genes with fewer than this number of codons. When the `function_sources` argument
            is used to return function rather than gene codon frequencies, the minimum codon filter
            is applied to genes before grouping them as functions. By default 0.
        function_min_codons : int, optional
            Ignore functions with fewer than this number of codons. Genes shorter than
            `gene_min_codons` are first removed and then functional groups of the remaining genes
            with fewer than `function_min_codons`. By default 0.
        min_codon_filter : {"length", "remaining", "both"}, optional
            This argument arises from the ambiguity of the minimum codon filters (`gene_min_codons`
            and `function_min_codons`) in relation to the filters that drop codons
            (`drop_amino_acids`, `min_amino_acids`/`min_gene_fraction`). Genes (and functions) can
            be filtered by their full LENGTH, e.g., genes shorter than 300 codons are ignored. They
            can also be filtered by the number of codons REMAINING after dropping codons defined by
            `drop_amino_acids` or determined dynamically by `min_amino_acids`/`min_gene_fraction`.
            The codon LENGTH filter followed by dropping codons can result in genes and functions
            with fewer codons than the original codon threshold -- thus the option of BOTH LENGTH
            and REMAINING filters to ensure that total codon frequencies in the output always meet
            the minimum codon threshold. The default filter type is BOTH.
        drop_amino_acids : iterable, optional
            Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala, and
            STP for stop codons). If `per_amino_acid` is True, the `drop_amino_acids` default rather
            than being None is STP plus amino acids encoded by a single codon (Met, Trp).
        min_amino_acids : int, optional
            Use this argument together with `min_gene_fraction`. Remove codons for amino acids (and
            STP) that are less numerous than `min_amino_acids` in a `min_gene_fraction` of genes.
            For example, say `min_amino_acids` is 5 and `min_gene_fraction` is 0.9. If there are
            fewer than 5 codons for an amino acid/STP in ≥90% of items, then the amino acid's codon
            columns are discarded. By default 0.
        min_gene_fraction : float, optional
            Use this argument together with `min_amino_acids`. This argument defines the fraction of
            genes in which a decoded amino acid must be found at the threshold frequency of
            `min_amino_acids` for the amino acid's codons to be retained. By default 0.
        label_amino_acids : bool, optional
            Include the amino acid for each codon in the column header of the output, i.e., LysAAA
            instead of AAA. By default False.

        Returns
        -------
        pandas.core.frame.DataFrame or dict of DataFrames
            Frequency table of gene or function x codon or amino acid. The Index of a gene frequency
            DataFrame comprises gene callers IDs. The MultiIndex of of a function frequency
            DataFrame comprises source, accession, and name of functions. The Index of a summed or
            averaged function frequency DataFrame comprises sources. Lastly, the Index of a
            single-row summed or averaged gene frequency DataFrame is 'all'.

        Examples
        --------
        Return the amino acid frequencies of each KEGG KOfam.
        >>> self.get_codon_frequencies(from_function_sources='KOfam', return_amino_acids=True)

        Return the summed amino acid frequencies of genes annotated by KEGG KOfams.
        >>> self.get_codon_frequencies(
            from_function_sources='KOfam', return_amino_acids=True, sum_genes=True)

        Return the codon relative frequencies of all genes.
        >>> self.get_codon_frequencies(as_relative_frequency=True)

        Return the relative frequencies of genes ≥300 codons in length.
        >>> self.get_codon_frequencies(as_relative_frequency=True, min_codons=300)

        Return the average relative frequencies of genes ≥300 codons in length.
        >>> self.get_codon_frequencies(as_relative_frequency=True, average_all=True, min_codons=300)

        Return the per-amino acid (synonymous) relative frequencies of genes. Remove genes <300
        codons, then remove amino acids with <5 codons in ≥90% of remaining genes, then remove genes
        with <300 codons remaining.
        >>> self.get_codon_frequencies(
            per_amino_acid=True, min_codons=300, min_amino_acids=5, min_gene_fraction=0.9)

        Return the per-amino acid (synonymous) relative frequencies of KEGG KOfams and all BRITE
        categorizations. Remove genes <300 codons, then remove amino acids with <5 codons in ≥90% of
        remaining genes, then remove genes with <300 codons remaining.
        >>> self.get_codon_frequencies(
            from_function_sources=['KOfam', 'KEGG_BRITE'],
            per_amino_acid=True,
            min_codons=300,
            min_amino_acids=5,
            min_gene_fraction=0.9)
        """
        function_sources = self._establish_function_sources(from_function_sources)
        if function_sources and from_function_sources != True:
            # Subset the gene function table to the requested sources.
            gene_function_df = \
                self.gene_function_df.set_index('source').loc[function_sources].reset_index()
        elif not function_sources:
            gene_function_df = pd.DataFrame(columns=self.gene_function_df.columns)

        # Check compatability of `return_amino_acids` with other arguments.
        if return_amino_acids and per_amino_acid:
            raise ConfigError("The argument `per_amino_acid` should only be True when "
                              "`return_amino_acids` is also True, as `per_amino_acid` returns "
                              "synonymous codon relative frequencies.")
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
            gene_caller_ids, function_accessions, function_names)

        if gene_subsetting:
            # Subset gene function table to genes of interest.
            gene_function_df = gene_function_df.set_index('gene_caller_id')
            gene_function_df = gene_function_df.loc[
                gene_function_df.index.intersection(gene_codon_frequency_df.index)].reset_index()

        if per_amino_acid and not as_relative_frequency:
            as_relative_frequency = True

        if sum_genes and average_genes:
            raise ConfigError("`sum_genes` and `average_genes` cannot both be True.")

        if min_codon_filter not in ['length', 'remaining', 'both']:
            raise ConfigError("`min_codon_filter` must be one of 'length', 'remaining', or 'both'.")
        # Set `filter_gene_length` `filter_function_length`, and `filter_remaining`.
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
            if per_amino_acid:
                drop_amino_acids = ignored_decodings_in_cub
            else:
                drop_amino_acids = []
        else:
            unrecognized_amino_acids = []
            for amino_acid in drop_amino_acids:
                if amino_acid not in translation_table:
                    unrecognized_amino_acids.append(amino_acid)
            if unrecognized_amino_acids:
                raise ConfigError("The following amino acids in `drop_amino_acids` are not "
                                  f"recognized: {', '.join(unrecognized_amino_acids)}")

        if bool(min_amino_acids) ^ bool(min_gene_fraction):
            raise ConfigError("`min_amino_acids` and `min_gene_fraction` must be used together.")

        if gene_subsetting:
            self.run.info_single(f"{len(gene_codon_frequency_df)} of {len(self.gene_caller_ids)} "
                                 "CDS selected from the genome")
        else:
            self.run.info_single(f"{len(gene_codon_frequency_df)} CDS in the genome")

        # Filter genes by length.
        gene_codon_frequency_df = self._get_frequency_table(
            gene_codon_frequency_df,
            input_min_codons=gene_min_codons,
            filter_input_codon_count=filter_gene_length)
        if gene_min_codons:
            self.run.info_single(f"{len(gene_codon_frequency_df)} CDS longer than "
                                 f"{gene_min_codons} codons")

        # Filter functions by the total number of codons in their genes.
        gene_function_frequency_df = gene_function_df.merge(
            gene_codon_frequency_df, how='inner', on='gene_caller_id')
        function_codon_frequency_df = gene_function_frequency_df.drop(
            'gene_caller_id', axis=1).groupby(['source', 'accession', 'function']).sum()
        function_codon_frequency_df = self._get_frequency_table(
            function_codon_frequency_df,
            input_min_codons=function_min_codons,
            filter_input_codon_count=filter_function_length)

        # Drop certain codons from gene codon frequency table. Filter genes by total codons
        # remaining.
        gene_codon_frequency_df = self._get_frequency_table(
            gene_codon_frequency_df,
            drop_amino_acids=drop_amino_acids,
            min_amino_acids=min_amino_acids,
            min_fraction=min_gene_fraction,
            min_codons=gene_min_codons,
            filter_output_codon_count=filter_gene_remaining_codons)
        if ((drop_amino_acids or min_amino_acids) and
            (gene_min_codons and filter_gene_remaining_codons)):
            self.run.info_single(f"{len(gene_codon_frequency_df)} CDS with ≥{gene_min_codons} "
                                 "codons remaining after dropping codons")

            # Filter functions by total codons remaining.
            gene_function_frequency_df = gene_function_df.merge(
                gene_codon_frequency_df, how='inner', on='gene_caller_id')
            function_codon_frequency_df = gene_function_frequency_df.drop(
                'gene_caller_id', axis=1).groupby(['source', 'accession', 'function']).sum()
            function_codon_frequency_df = self._get_frequency_table(
                function_codon_frequency_df,
                min_codons=function_min_codons,
                filter_output_codon_count=filter_function_remaining_codons)

        # Per-gene results:
        get_table = lambda method: method(gene_codon_frequency_df,
                                          label_amino_acids=label_amino_acids,
                                          output_amino_acids=return_amino_acids)
        ### Absolute frequencies
        if (not as_relative_frequency and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not per_amino_acid and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)
        ### Per-amino acid relative frequencies
        if (not return_amino_acids and
            as_relative_frequency and
            per_amino_acid and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_per_amino_acid_codon_rel_frequency_table)

        # Summed codon/amino acid results across genes:
        # Call `get_codon_frequency_table` on the result of the method to add amino acid names to
        # the header.
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df).to_frame('all').T,
            label_amino_acids=label_amino_acids,
            output_amino_acids=return_amino_acids)
        ### Absolute frequencies
        if (not as_relative_frequency and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not per_amino_acid and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Per-amino acid relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df, label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            per_amino_acid and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_per_amino_acid_codon_rel_frequency_table)

        # Codon results averaged across genes:
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df).to_frame('all').T,
            label_amino_acids=label_amino_acids)
        ### Absolute frequencies
        if (not return_amino_acids and
            not as_relative_frequency and
            not function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (not return_amino_acids and
            as_relative_frequency and
            not per_amino_acid and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Per-amino acid relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df, label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            per_amino_acid and
            not function_sources and
            average_genes):
            return get_table(self._get_average_per_amino_acid_codon_rel_frequency_table)

        # Amino acid results averaged across genes:
        # Gene amino acid frequencies are first calculated from codon frequencies before averaging
        # across genes.
        get_table = lambda method: method(
            self._get_frequency_table(
                gene_codon_frequency_df, output_amino_acids=return_amino_acids).to_frame('all').T)
        ### Absolute frequencies
        if (return_amino_acids and
            not as_relative_frequency and
            not function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (return_amino_acids and
            as_relative_frequency and
            not per_amino_acid and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)

        # Per-function codon/amino acid results:
        function_codon_frequency_df = function_codon_frequency_df.set_index(
            ['source', 'accession', 'function'])
        ### Absolute frequencies
        get_table = lambda method: method(function_codon_frequency_df,
                                          label_amino_acids=label_amino_acids,
                                          output_amino_acids=return_amino_acids)
        if (not as_relative_frequency and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not per_amino_acid and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)
        ### Per-amino acid relative frequencies
        if (not return_amino_acids and
            as_relative_frequency and
            per_amino_acid and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_per_amino_acid_codon_rel_frequency_table)

        # Summed codon/amino acid results across functions in each source:
        # If amino acid rather than codon columns are returned, then gene amino acid frequencies are
        # first calculated from the sum of codon frequencies.
        get_table = lambda method: self._get_frequency_table(
            self._get_frequency_table(
                function_codon_frequency_df, output_amino_acids=return_amino_acids).groupby(
                    'source').apply(method),
            label_amino_acids=label_amino_acids)
        ### Absolute frequencies
        if (not as_relative_frequency and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Per-amino acid relative frequencies
        get_table = lambda method: self._get_frequency_table(
            function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
                label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            per_amino_acid and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_per_amino_acid_codon_rel_frequency_table)

        # Codon results averaged across genes annotated by each function source:
        gene_function_frequency_df = gene_function_frequency_df.set_index(
            ['source', 'accession', 'function', 'gene_caller_id'])
        get_table = lambda method: self._get_frequency_table(
            gene_function_frequency_df.groupby('source').apply(method),
            label_amino_acids=label_amino_acids)
        ### Absolute frequencies
        if (not return_amino_acids and
            not as_relative_frequency and
            function_sources and
            sum_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (not return_amino_acids and
            as_relative_frequency and
            not per_amino_acid and
            function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Per-amino acid relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_frequency_df.groupby('source').apply(method).droplevel(1),
            label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            per_amino_acid and
            function_sources and
            average_genes):
            return get_table(self._get_average_per_amino_acid_codon_rel_frequency_table)

        # Amino acid results averaged across genes annotated by each function source:
        get_table = lambda method: self._get_frequency_table(
            gene_function_frequency_df, output_amino_acids=return_amino_acids).groupby(
                'source').apply(method)
        ### Absolute frequencies
        if (return_amino_acids and
            not as_relative_frequency and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_frequency_table)
        ### Relative frequencies
        if (return_amino_acids and
            as_relative_frequency and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)

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


    def _establish_function_sources(self, from_function_sources):
        """Establishes `function_sources` for methods that take the arguments,
        `from_function_sources` and `function_codon_frequency_dict`."""
        if from_function_sources == True:
            function_sources = list(self.function_annotation_sources)
        elif from_function_sources != False:
            function_sources = list(from_function_sources)
        else:
            function_sources = None

        if function_sources is not None:
            unrecognized_function_sources = []
            for function_source in function_sources:
                if function_source not in self.function_annotation_sources:
                    unrecognized_function_sources.append(function_source)
            if unrecognized_function_sources:
                raise ConfigError("The requested function annotation sources, "
                                  f"{', '.join(function_sources)}, are not among those available: "
                                  f"{', '.join(self.function_annotation_sources)}.")

        return function_sources


    def _select_genes(self, gene_caller_ids, function_accession_dict, function_name_dict):
        """Select gene codon frequencies based not only on a list of IDs but also functions of
        interest."""
        select_gene_caller_ids = []

        unrecognized_gene_caller_ids = []
        for gene_caller_id in gene_caller_ids:
            if gene_caller_id not in self.gene_caller_ids:
                unrecognized_gene_caller_ids.append(gene_caller_id)
                continue
            select_gene_caller_ids.append(gene_caller_id)
        if unrecognized_gene_caller_ids:
            raise ConfigError("The following gene caller IDs were not found in the genome: "
                              f"{unrecognized_gene_caller_ids}")

        index_keys = []
        unrecognized_sources = []
        for function_source, function_accessions in function_accession_dict.items():
            if function_source not in self.function_annotation_sources:
                unrecognized_sources.append(function_source)
            for function_accession in function_accessions:
                index_keys.append((function_source, function_accession))
        if unrecognized_sources:
            raise ConfigError("The following annotation sources in `function_accessions` were not "
                              f"found as having annotated the genome: {unrecognized_sources}")
        gene_function_df = self.gene_function_df.set_index(['source', 'accession'])
        select_gene_caller_ids += gene_function_df.loc[
            gene_function_df.index.intersection(index_keys)]['gene_caller_id'].tolist()

        index_keys = []
        unrecognized_sources = []
        for function_source, function_names in function_name_dict.items():
            if function_source not in self.function_annotation_sources:
                unrecognized_sources.append(function_source)
            for function_name in function_names:
                index_keys.append((function_source, function_name))
        if unrecognized_sources:
            raise ConfigError("The following annotation sources in `function_names` were not "
                              f"found as having annotated the genome: {unrecognized_sources}")
        gene_function_df = gene_function_df.reset_index.set_index(['source', 'function'])
        select_gene_caller_ids += gene_function_df.loc[
            gene_function_df.index.intersection(index_keys)]['gene_caller_id'].tolist()

        return self.gene_codon_frequency_df.loc[set(select_gene_caller_ids)]


    def _filter_input_codon_count(func):
        """Decorator to discard rows in the input frequency table with fewer than the minimum number
        of codons/amino acids."""
        def wrapper(*args, **kwargs):
            frequency_df = args[0]
            try:
                filter_input_codon_count = kwargs['filter_input_codon_count']
                input_min_codons = kwargs['input_min_codons']
            except KeyError:
                filter_input_codon_count = False
                input_min_codons = 0
            if filter_input_codon_count and input_min_codons:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= input_min_codons]
            return func(frequency_df, *args[1: ], **kwargs)
        return wrapper


    def _drop_amino_acid_codon_columns(func):
        """Decorator to discard codon columns by amino acid from the input frequency table."""

        def wrapper(*args, **kwargs):
            codon_frequency_df = args[0]
            try:
                drop_amino_acids = kwargs['drop_amino_acids']
            except KeyError:
                drop_amino_acids = []
            if drop_amino_acids:
                drop_codons = [translation_table[aa] for aa in drop_amino_acids]
                codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1, errors='ignore')
            return func(codon_frequency_df, *args[1: ], **kwargs)
        return wrapper


    def _filter_synonymous_codon_count(func):
        """Decorator to discard codon columns from the input frequency table with fewer than the
        minimum number of synonymous codons in a minimum number of rows."""
        def wrapper(*args, **kwargs):
            codon_frequency_df = args[0]
            try:
                min_amino_acids = kwargs['min_amino_acids']
                min_fraction = kwargs['min_fraction']
            except KeyError:
                min_amino_acids = 0
                min_fraction = 0
            if min_amino_acids and min_fraction:
                row_count = len(codon_frequency_df)
                drop_codons = []
                for synonymous_codons in translation_table.values():
                    try:
                        filtered_row_count = codon_frequency_df[
                            codon_frequency_df[synonymous_codons].sum() >= min_amino_acids].sum()
                    except KeyError:
                        # This occurs when codons are missing from the frequency table.
                        continue
                    if filtered_row_count / row_count < min_fraction:
                        drop_codons += synonymous_codons
                codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1)
            return func(codon_frequency_df, *args[1: ], **kwargs)
        return wrapper


    def _output_amino_acids(func):
        """Decorator to output columns of amino acid rather than codon frequencies."""
        def wrapper(*args, **kwargs):
            codon_df = func(*args, **kwargs)
            try:
                output_amino_acids = kwargs['output_amino_acids']
            except KeyError:
                output_amino_acids = False
            if output_amino_acids:
                aa_df = pd.DataFrame(index=codon_df.index)
                for amino_acid, codons in translation_table.items():
                    try:
                        aa_df[amino_acid] = codon_df[codons].sum(axis=1)
                    except KeyError:
                        # This occurs when there aren't columns for codons.
                        pass
                return aa_df
            return codon_df


    def _add_amino_acid_to_header(func):
        """Decorator to add amino acid to codon column header."""
        def wrapper(*args, **kwargs):
            codon_df = func(*args, **kwargs)
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


    def _filter_output_codon_count(func):
        """Decorator to discard rows in the output frequency table with fewer than the minimum
        number of codons/amino acids."""
        def wrapper(*args, **kwargs):
            frequency_df = func(*args, **kwargs)
            try:
                filter_output_codon_count = kwargs['filter_output_codon_count']
                min_codons = kwargs['min_codons']
            except KeyError:
                filter_output_codon_count = False
                min_codons = kwargs['min_codons']
            if filter_output_codon_count:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= min_codons]
            return frequency_df
        return wrapper


    # The order of decorators should not be changed (only @_output_amino_acids and
    # @_add_amino_acid_to_header, which are mutually exclusive operations, are interchangeable).
    @_filter_input_codon_count
    @_drop_amino_acid_codon_columns
    @_filter_synonymous_codon_count
    @_output_amino_acids
    @_add_amino_acid_to_header
    @_filter_output_codon_count
    def _get_frequency_table(self, codon_frequency_df, **kwargs):
        return codon_frequency_df


    # Commented decorators mean that they can theoretically be used but are not because they are not
    # needed in the `self.get_codon_frequencies` client.
    # @_filter_input_codon_count
    # @_drop_amino_acid_codon_columns
    # @_filter_synonymous_codon_count
    @_output_amino_acids
    @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_rel_frequency_table(self, codon_frequency_df, **kwargs):
        codon_relative_frequency_df = codon_frequency_df.div(codon_frequency_df.sum(axis=1))
        return codon_relative_frequency_df


    # @_filter_input_codon_count
    # @_drop_amino_acid_codon_columns
    # @_filter_synonymous_codon_count
    @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_per_amino_acid_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        """Return the relative frequencies of codons in relation to the set of codons encoding the
        same amino acid (or stop codons)."""
        per_aa_codon_rel_frequency_df = pd.DataFrame()
        for codons in translation_table.values():
            try:
                aa_codon_frequency_df = codon_frequency_df[codons]
            except KeyError:
                # This occurs when codons are missing from the frequency table.
                continue
            per_aa_codon_rel_frequency_df[codons] = aa_codon_frequency_df.div(
                aa_codon_frequency_df.sum(axis=1))
        return per_aa_codon_rel_frequency_df


    # @_filter_input_codon_count
    # @_drop_amino_acid_codon_columns
    # @_filter_synonymous_codon_count
    def _get_summed_frequency_table(self, frequency_df, **kwargs):
        """Return the summed frequencies across all items."""
        summed_frequency_series = frequency_df.sum()
        return summed_frequency_series


    def _get_summed_rel_frequency_table(self, frequency_df, **kwargs):
        summed_frequency_series = self._get_summed_frequency_table(frequency_df, kwargs)
        summed_rel_frequency_series = summed_frequency_series.div(summed_frequency_series.sum())
        return summed_rel_frequency_series


    @_add_amino_acid_to_header
    def _get_summed_per_amino_acid_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        kwargs_subset = {}
        for key, value in kwargs.items():
            if key == 'label_amino_acids':
                continue
            kwargs_subset[key] = value
        summed_codon_frequency_series = self._get_summed_frequency_table(
            codon_frequency_df, **kwargs_subset)
        summed_codon_frequency_df = summed_codon_frequency_series.to_frame('all').T

        summed_per_amino_acid_codon_rel_frequency_df = \
            self._get_per_amino_acid_codon_rel_frequency_table(summed_codon_frequency_df)

        return summed_per_amino_acid_codon_rel_frequency_df


    # @_filter_input_codon_count
    # @_drop_amino_acid_codon_columns
    # @_filter_synonymous_codon_count
    def _get_average_frequency_table(self, frequency_df, **kwargs):
        """Return the average codon frequencies across all items."""
        average_codon_frequency_series = frequency_df.mean()
        return average_codon_frequency_series


    def _get_average_rel_frequency_table(self, frequency_df, **kwargs):
        average_frequency_series = self._get_average_frequency_table(frequency_df, **kwargs)
        average_rel_frequency_series = average_frequency_series.div(average_frequency_series.sum())
        return average_rel_frequency_series


    @_add_amino_acid_to_header
    def _get_average_per_amino_acid_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        kwargs_subset = {}
        for key, value in kwargs.items():
            if key == 'label_amino_acids':
                continue
            kwargs_subset[key] = value
        average_codon_frequency_series = self._get_average_frequency_table(
            codon_frequency_df, **kwargs)
        average_codon_frequency_df = average_codon_frequency_series.to_frame('all').T

        average_per_amino_acid_codon_rel_frequency_df = \
            self._get_per_amino_acid_codon_rel_frequency_table(average_codon_frequency_df)

        return average_per_amino_acid_codon_rel_frequency_df

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
