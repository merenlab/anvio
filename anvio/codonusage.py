# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module for codon usage analyses at the levels of genes, groups of genes, and genomes"""


import inspect
import argparse
import numpy as np
import pandas as pd

from copy import deepcopy
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

    Manipulate the raw data using the methods, `get_frequencies` and `get_codon_usage_bias`.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.contigs_db_path = A('contigs_db')

        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')

        self.function_sources = A('function_sources')

        self.run = run
        self.progress = progress

        self._load_contigs_db_data()
        self._make_gene_codon_frequency_table()


    def _load_contigs_db_data(self):
        """Load gene data from the contigs database."""
        utils.is_contigs_db(self.contigs_db_path)

        if self.profile_db_path or self.collection_name or self.bin_id:
            # Initialize the contigs superclass from the splits of the internal genome bin.
            self.args.split_names_of_interest = \
                ccollections.GetSplitNamesInBins(self.args).get_split_names_only()
        contigs_super = ContigsSuperclass(self.args, r=run_quiet)
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
                    "database has not been annotated by any sources :/")
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
                            # Include every possible depth of categorization.
                            hierarchy_accession = accession
                            categorization = name
                            split_categories = categorization.split('>>>')
                            for depth in range(1, len(split_categories) + 1):
                                gene_function_rows.append(
                                    [annotation_source,
                                     hierarchy_accession,
                                     '>>>'.join(split_categories[: depth]),
                                     gene_caller_id])
                        else:
                            gene_function_rows.append(
                                [annotation_source, accession, name, gene_caller_id])
                else:
                    # The function accession and name entries do not contain the same number of
                    # '!!!' separators. In COG20_PATHWAY, there can be multiple accessions
                    # corresponding to the same function name.
                    if annotation_source == 'KEGG_BRITE':
                        # Include every possible depth of categorization.
                        hierarchy_accession = accession
                        categorization = name
                        split_categories = categorization.split('>>>')
                        for depth in range(1, len(split_categories) + 1):
                            gene_function_rows.append(
                                [annotation_source,
                                 hierarchy_accession,
                                 '>>>'.join(split_categories[: depth]),
                                 gene_caller_id])
                    else:
                        gene_function_rows.append(
                            [annotation_source, accession, name, gene_caller_id])
        self.gene_function_df = pd.DataFrame(
            gene_function_rows, columns=['source', 'accession', 'name', 'gene_caller_id'])


    def _make_gene_codon_frequency_table(self):
        """Generate the per-gene codon frequency DataFrame as `self.gene_codon_frequency_df`."""
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

            gene_codon_frequencies.append(Counter(utils.get_list_of_codons_for_gene_call(
                gene_call, self.contig_sequences_dict)))

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

        if skipped_noncoding_gene_caller_ids:
            self.run.warning(
                f"{len(skipped_noncoding_gene_caller_ids)} of {len(self.gene_caller_ids)} genes "
                "were non-coding and not added to the codon frequency table.")


    def get_frequencies(self,
                        from_function_sources=False,
                        return_functions=False,
                        return_amino_acids=False,
                        gene_caller_ids=None,
                        function_accessions=None,
                        function_names=None,
                        as_relative_frequency=False,
                        as_synonymous=False,
                        sum_genes=False,
                        average_genes=False,
                        gene_min_codons=0,
                        function_min_codons=0,
                        min_codon_filter='both',
                        drop_amino_acids=None,
                        min_amino_acids=0,
                        min_gene_fraction=1,
                        label_amino_acids=False):
        """
        Get absolute (default) or relative codon or amino acid frequencies from genes or functions.

        Relative codon frequencies can be normalized per-amino acid to all synonymous codons.

        Parameters
        ----------
        from_function_sources : bool, str, or iterable of str, optional
            Select genes with functional annotations. With this argument, the first four columns of
            the returned table contain, respectively, gene caller IDs, function annotation sources,
            function accessions, and function names. There is a row for each gene/function
            combination in the table, and each row for the same gene contains the same frequency
            values. When this argument is True, use all available functional annotation sources in
            the SingleGenomeCodonUsage object. When this argument is a string, select the source
            given by the string, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. When this argument is an
            iterable, select a subset of sources given by its strings, e.g., ['KOfam',
            'COG20_FUNCTION']. When 'KEGG_BRITE' is in the argument, in the absence of
            `function_names` narrowing the scope of the inquiry, all possible depths of each BRITE
            hierarchy in the data are returned, e.g., 'Ribosome>>>Ribosomal proteins' and the more
            general 'Ribosome' would each have rows in the output table.
        return_functions : bool, optional
            If True (default False), output frequency tables contain function rather than gene
            results in each row. Returning per-gene results when also considering functions by using
            `--from-function-sources` facilitates analysis of per-gene variability within functions.
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
            If True (default False), return relative rather than absolute codon frequencies.
        as_synonymous : bool, optional
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
        ---------------------
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
            This argument arises from the ambiguity of the minimum codon filters (`gene_min_codons`
            and `function_min_codons`) in relation to the filters that drop codons
            (`drop_amino_acids`, `min_amino_acids`/`min_gene_fraction`). Genes and functions can be
            filtered by their full LENGTH, e.g., genes shorter than 300 codons are ignored. They can
            also be filtered by the number of codons REMAINING after dropping codons set by
            `drop_amino_acids` or determined dynamically by `min_amino_acids`/`min_gene_fraction`.
            The codon LENGTH filter followed by dropping codons can result in genes and functions
            with fewer codons than the original codon threshold -- thus the option of BOTH LENGTH
            and REMAINING filters to ensure that total codon frequencies in the output always meet
            the minimum codon threshold. The default filter type is BOTH.
        drop_amino_acids : iterable, optional
            Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala, and
            STP for stop codons). If `as_synonymous` is True, the `drop_amino_acids` default rather
            than being None is STP plus amino acids encoded by a single codon (Met, Trp).
        min_amino_acids : int, optional
            Use this argument together with `min_gene_fraction`. Remove codons for amino acids (and
            STP) that are less numerous than `min_amino_acids` in a `min_gene_fraction` of genes.
            For example, say `min_amino_acids` is 5 and `min_gene_fraction` is 0.9. If there are
            fewer than 5 codons for an amino acid/STP in ≥90% of genes, then the amino acid's codons
            are removed. By default 0.
        min_gene_fraction : float, optional
            Use this argument together with `min_amino_acids`. This argument defines the fraction of
            genes in which a decoded amino acid must be found at the threshold frequency of
            `min_amino_acids` for the amino acid's codons to be removed. By default 1.
        label_amino_acids : bool, optional
            If True (default False), include the amino acid for each codon in the column header of
            the output, i.e., LysAAA instead of AAA.

        Returns
        -------
        pandas.core.frame.DataFrame
            Frequency table of gene x codon or amino acid. If functions are not considered, then the
            Index of the returned DataFrame contains gene caller IDs. If functions are considered,
            then each row represents a gene/function pair, and the same gene frequencies can be
            found in multiple rows with different function pairs; there are additional MultiIndex
            columns for source, accession, and name of function. If frequencies are summed or
            averaged and functions are not considered, the returned single-row DataFrame has an
            Index with one entry, 'all'. If frequencies are summed or averaged and functions are
            considered, the returned DataFrame has a row for each function annotation source.

        Examples
        --------
        Return the amino acid frequencies of each gene annotated by KEGG KOfams.
        >>> self.get_frequencies(from_function_sources='KOfam', return_amino_acids=True)

        Return the summed amino acid frequencies of genes annotated by KEGG KOfams.
        >>> self.get_frequencies(
            from_function_sources='KOfam', return_amino_acids=True, sum_genes=True)

        Return the codon relative frequencies of each gene.
        >>> self.get_frequencies(as_relative_frequency=True)

        Return the relative frequencies of genes ≥300 codons in length.
        >>> self.get_frequencies(as_relative_frequency=True, gene_min_codons=300)

        Return the average relative frequencies of genes ≥300 codons in length.
        >>> self.get_frequencies(as_relative_frequency=True, average_all=True, gene_min_codons=300)

        Return the synonymous (per-amino acid) relative frequencies of genes. Remove genes <300
        codons, then remove amino acids with <5 codons in ≥90% of remaining genes, then remove genes
        with <300 codons remaining.
        >>> self.get_frequencies(
            as_synonymous=True, gene_min_codons=300, min_amino_acids=5, min_gene_fraction=0.9)

        Return the synonymous (per-amino acid) relative frequencies of KEGG KOfam and BRITE
        functions. This command removes genes <300 codons, then removes amino acids with <5 codons
        in ≥90% of remaining genes, then removes genes with <300 codons remaining, then sums gene
        frequencies in functions, then calculates synonymous relative frequencies in functions.
        >>> self.get_frequencies(
            from_function_sources=['KOfam', 'KEGG_BRITE'],
            return_functions=True,
            as_synonymous=True,
            gene_min_codons=300,
            min_amino_acids=5,
            min_gene_fraction=0.9)
        """
        # CHECK ARGUMENTS AND SET UP PROCEDURE
        ######################################
        function_sources = self._establish_function_sources(from_function_sources)
        if function_sources == self.function_sources:
            gene_function_df = self.gene_function_df
        elif function_sources:
            # Subset the gene function table to the requested sources.
            gene_function_df = \
                self.gene_function_df.set_index('source').loc[function_sources].reset_index()

        # Check compatability of `return_amino_acids` with other arguments.
        if return_amino_acids and as_synonymous:
            raise ConfigError("The argument `as_synonymous` should only be True when "
                              "`return_amino_acids` is also True, as `as_synonymous` returns "
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

        if from_function_sources and gene_subsetting:
            # Subset gene function table to genes of interest.
            gene_function_df = gene_function_df.set_index('gene_caller_id')
            gene_function_df = gene_function_df.loc[
                gene_function_df.index.intersection(gene_codon_frequency_df.index)].reset_index()

        if as_synonymous and not as_relative_frequency:
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
            if as_synonymous:
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

        if min_amino_acids > 0 and min_gene_fraction == 1:
            raise ConfigError(
                "A positive value of `min_amino_acids`, rather than the default value of 0, "
                "indicates that codons should be dynamically filtered by their occurrence across "
                "genes, but `min_gene_fraction` has the default value of 1, so no filtering would "
                "be done.")
        if min_amino_acids == 0 and min_gene_fraction < 1:
            raise ConfigError(
                "A value of `min_gene_fraction` of less than 1, rather than the default value of "
                "1, indicates that codons should be dynamically filtered by their occurrence "
                "across genes, but `min_amino_acids` has the default value of 0, so no filtering "
                "would be done.")

        if gene_subsetting:
            self.run.info_single(f"{len(gene_codon_frequency_df)} of {len(self.gene_caller_ids)} "
                                 "CDS selected from the genome",
                                 nl_before=1,
                                 nl_after=1)
        else:
            self.run.info_single(f"{len(gene_codon_frequency_df)} CDS in the genome",
                                 nl_before=1,
                                 nl_after=1)

        # FILTER GENES/FUNCTIONS AND CODONS USING ADDITIONAL PARAMETERS
        ###############################################################
        # Filter genes by length.
        gene_codon_frequency_df = self._get_frequency_table(
            gene_codon_frequency_df,
            input_min_codons=gene_min_codons,
            filter_input_codon_count=filter_gene_length)

        # Filter functions by the total length of their genes. Gene/function pairs with fewer than
        # the required number of codons in the function are removed.
        if function_sources:
            gene_function_codon_frequency_df = gene_function_df.merge(
                gene_codon_frequency_df, how='inner', on='gene_caller_id')
            gene_function_codon_frequency_df = gene_function_codon_frequency_df.set_index(
                ['source', 'accession', 'name', 'gene_caller_id'])
            if filter_function_length:
                gene_function_codon_frequency_df = gene_function_codon_frequency_df.groupby(
                    ['source', 'accession', 'name']).filter(
                        lambda function_df: function_df.sum(axis=1).sum() >= function_min_codons)

        # Drop certain codons from gene codon frequency table. Filter genes by total codons
        # remaining.
        gene_codon_frequency_df = self._get_frequency_table(
            gene_codon_frequency_df,
            drop_amino_acids=drop_amino_acids,
            min_amino_acids=min_amino_acids,
            min_fraction=min_gene_fraction,
            min_codons=gene_min_codons,
            filter_output_codon_count=filter_gene_remaining_codons)

        need_to_filter_codons_in_gene_function_codon_frequency_table = True
        if ((drop_amino_acids or min_amino_acids) and
            (gene_min_codons and filter_gene_remaining_codons)):
            # Filter functions by total codons remaining. Gene/function pairs with fewer than the
            # required number of codons in the function are removed.
            if function_sources:
                gene_function_codon_frequency_df = gene_function_df.merge(
                    gene_codon_frequency_df, how='inner', on='gene_caller_id')
                gene_function_codon_frequency_df = gene_function_codon_frequency_df.set_index(
                    ['source', 'accession', 'name', 'gene_caller_id'])
                if filter_function_remaining_codons:
                    gene_function_codon_frequency_df = gene_function_codon_frequency_df.groupby(
                        ['source', 'accession', 'name']).filter(
                            lambda function_df:
                                function_df.sum(axis=1).sum() >= function_min_codons)
                need_to_filter_codons_in_gene_function_codon_frequency_table = False

        if function_sources and need_to_filter_codons_in_gene_function_codon_frequency_table:
            gene_function_codon_frequency_df = self._get_frequency_table(
                gene_function_codon_frequency_df,
                drop_amino_acids=drop_amino_acids,
                min_amino_acids=min_amino_acids,
                min_fraction=min_gene_fraction)

        if gene_min_codons:
            self.run.info_single(f"{len(gene_codon_frequency_df)} CDS remaining after codon count filters")

        # GET OUTPUTS WITH NO CONSIDERATION OF FUNCTIONS
        ################################################
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
            not as_synonymous and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)
        ### Synonymous (per-amino acid) relative frequencies
        if (not return_amino_acids and
            as_relative_frequency and
            as_synonymous and
            not function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_synonymous_codon_rel_frequency_table)

        # Summed codon/amino acid results across genes:
        # Call `_get_codon_frequency_table` on the result of the method to add amino acid names to
        # the header or change codon to amino acid frequencies (change columns). Note in this and
        # further `get_table` functions that `_get_codon_frequency_table` is used for the purpose of
        # its decorator functions.
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
            not as_synonymous and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df, label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            as_synonymous and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

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
            not as_synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df, label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            as_synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_synonymous_codon_rel_frequency_table)

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
            not as_synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)

        # GET OUTPUTS WITH CONSIDERATION OF FUNCTIONS
        #############################################
        get_table = lambda method: method(
            (gene_function_codon_frequency_df.groupby(['source', 'accession', 'name']).sum()
             if return_functions else
             gene_function_codon_frequency_df.sort_values(['source', 'accession', 'name'])),
            label_amino_acids=label_amino_acids,
            output_amino_acids=return_amino_acids)
        ### Absolute frequencies
        if (not as_relative_frequency and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_frequency_table)
        ### Relative frequencies
        if (as_relative_frequency and
            not as_synonymous and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_rel_frequency_table)
        ### Synonymous relative frequencies
        if (not return_amino_acids and
            as_relative_frequency and
            as_synonymous and
            function_sources and
            not sum_genes and
            not average_genes):
            return get_table(self._get_synonymous_codon_rel_frequency_table)

        # Remove duplicate occurrences of genes before summing or averaging frequencies of all genes
        # in the function source. A gene can have different annotations, e.g., different KOfam
        # assignments. In KEGG BRITE, a gene can be in nested categories and different hierarchies.
        gene_function_codon_frequency_df = \
            gene_function_codon_frequency_df.reset_index().drop_duplicates(
                subset='gene_caller_id', ignore_index=True).set_index(
                    ['source', 'accession', 'name', 'gene_caller_id'])

        # Summed codon/amino acid results across each function source:
        # If amino acid rather than codon columns are returned, then gene amino acid frequencies are
        # first calculated from the sum of codon frequencies.
        get_table = lambda method: self._get_frequency_table(
            self._get_frequency_table(
                gene_function_codon_frequency_df, output_amino_acids=return_amino_acids).groupby(
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
        ### Synonymous relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
                label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            as_synonymous and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

        # Codon results averaged across genes annotated by a function source:
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('source').apply(method),
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
            not as_synonymous and
            function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
            label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            as_relative_frequency and
            as_synonymous and
            function_sources and
            average_genes):
            return get_table(self._get_average_synonymous_codon_rel_frequency_table)

        # Amino acid results averaged across genes annotated by a function source:
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df, output_amino_acids=return_amino_acids).groupby(
                'source').apply(method)
        ### Absolute frequencies
        if (return_amino_acids and
            not as_relative_frequency and
            function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (return_amino_acids and
            as_relative_frequency and
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


    # `self.get_frequencies` SETUP HELPER METHODS
    #############################################
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


    def _select_genes(self, gene_caller_ids, function_accession_dict, function_name_dict):
        """Select genes in the gene frequency table given not only a list of IDs but also functions
        of interest."""
        if not gene_caller_ids and not function_accession_dict and not function_name_dict:
            return self.gene_codon_frequency_df

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
            if function_source not in self.function_sources:
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
            if function_source not in self.function_sources:
                unrecognized_sources.append(function_source)
            for function_name in function_names:
                index_keys.append((function_source, function_name))
        if unrecognized_sources:
            raise ConfigError("The following annotation sources in `function_names` were not "
                              f"found as having annotated the genome: {unrecognized_sources}")
        gene_function_df = gene_function_df.reset_index().set_index(['source', 'name'])
        select_gene_caller_ids += gene_function_df.loc[
            gene_function_df.index.intersection(index_keys)]['gene_caller_id'].tolist()

        return self.gene_codon_frequency_df.loc[set(select_gene_caller_ids)]


    # `self.get_frequencies` OUTPUT HELPER METHODS
    ##############################################
    def _filter_input_codon_count(method):
        """Decorator to discard rows in the input frequency table with fewer than the minimum number
        of codons/amino acids."""
        def wrapper(*args, **kwargs):
            frequency_df = args[1]
            try:
                filter_input_codon_count = kwargs['filter_input_codon_count']
                input_min_codons = kwargs['input_min_codons']
            except KeyError:
                filter_input_codon_count = False
                input_min_codons = 0
            if filter_input_codon_count and input_min_codons:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= input_min_codons]
            return method(args[0], frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _drop_amino_acid_codon_columns(method):
        """Decorator to discard codon columns by amino acid from the input frequency table."""
        def wrapper(*args, **kwargs):
            codon_frequency_df = args[1]
            try:
                drop_amino_acids = kwargs['drop_amino_acids']
            except KeyError:
                drop_amino_acids = []
            if drop_amino_acids:
                drop_codons = []
                for aa in drop_amino_acids:
                    drop_codons += translation_table[aa]
                codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1, errors='ignore')
            return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _filter_synonymous_codon_count(method):
        """Decorator to discard codon columns from the input frequency table with fewer than the
        minimum number of synonymous codons in a minimum number of rows."""
        def wrapper(*args, **kwargs):
            codon_frequency_df = args[1]
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
                        filtered_row_count = len(codon_frequency_df[codon_frequency_df[
                            synonymous_codons].sum(axis=1) >= min_amino_acids])
                    except KeyError:
                        # This occurs when codons are missing from the frequency table.
                        continue
                    if filtered_row_count / row_count < min_fraction:
                        drop_codons += synonymous_codons
                codon_frequency_df = codon_frequency_df.drop(drop_codons, axis=1)
            return method(args[0], codon_frequency_df, *args[2: ], **kwargs)
        return wrapper


    def _output_amino_acids(method):
        """Decorator to output columns of amino acid rather than codon frequencies."""
        def wrapper(*args, **kwargs):
            codon_df = method(*args, **kwargs)
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
        return wrapper


    def _add_amino_acid_to_header(method):
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


    def _filter_output_codon_count(method):
        """Decorator to discard rows in the output frequency table with fewer than the minimum
        number of codons/amino acids."""
        def wrapper(*args, **kwargs):
            frequency_df = method(*args, **kwargs)
            try:
                filter_output_codon_count = kwargs['filter_output_codon_count']
                min_codons = kwargs['input_min_codons']
            except KeyError:
                filter_output_codon_count = False
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


    # Commented decorators mean that they can theoretically be used and uncommented but are not
    # because they are not needed in the `self.get_frequencies` client.
    # @_filter_input_codon_count
    # @_drop_amino_acid_codon_columns
    # @_filter_synonymous_codon_count
    @_output_amino_acids
    @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_rel_frequency_table(self, codon_frequency_df, **kwargs):
        codon_relative_frequency_df = codon_frequency_df.div(codon_frequency_df.sum(axis=1), axis=0)
        # Drop rows with zero frequency.
        codon_relative_frequency_df = codon_relative_frequency_df.dropna()
        return codon_relative_frequency_df


    # @_filter_input_codon_count
    # @_drop_amino_acid_codon_columns
    # @_filter_synonymous_codon_count
    @_add_amino_acid_to_header
    # @_filter_output_codon_count
    def _get_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        """Return the relative frequencies of codons in relation to the set of codons encoding the
        same amino acid (or stop codons)."""
        synonymous_codon_rel_frequency_df = pd.DataFrame()
        for codons in translation_table.values():
            try:
                aa_codon_frequency_df = codon_frequency_df[codons]
            except KeyError:
                # This occurs when codons are missing from the frequency table.
                continue
            synonymous_codon_rel_frequency_df[codons] = aa_codon_frequency_df.div(
                aa_codon_frequency_df.sum(axis=1), axis=0)
        synonymous_codon_rel_frequency_df = synonymous_codon_rel_frequency_df.dropna(how='all')
        return synonymous_codon_rel_frequency_df


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
    def _get_summed_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        kwargs_subset = {}
        for key, value in kwargs.items():
            if key == 'label_amino_acids':
                continue
            kwargs_subset[key] = value
        summed_codon_frequency_series = self._get_summed_frequency_table(
            codon_frequency_df, **kwargs_subset)
        summed_codon_frequency_df = summed_codon_frequency_series.to_frame('all').T

        summed_synonymous_codon_rel_frequency_df = \
            self._get_synonymous_codon_rel_frequency_table(summed_codon_frequency_df)

        return summed_synonymous_codon_rel_frequency_df


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
    def _get_average_synonymous_codon_rel_frequency_table(self, codon_frequency_df, **kwargs):
        kwargs_subset = {}
        for key, value in kwargs.items():
            if key == 'label_amino_acids':
                continue
            kwargs_subset[key] = value
        average_codon_frequency_series = self._get_average_frequency_table(
            codon_frequency_df, **kwargs)
        average_codon_frequency_df = average_codon_frequency_series.to_frame('all').T

        average_synonymous_codon_rel_frequency_df = \
            self._get_synonymous_codon_rel_frequency_table(average_codon_frequency_df)

        return average_synonymous_codon_rel_frequency_df

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


class MultiGenomeCodonUsage(object):
    """This object processes codon usage data from multiple internal and/or external genomes."""

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.internal_genomes_path = A('internal_genomes')
        self.external_genomes_path = A('external_genomes')

        gene_output_path = A('gene_table_output')
        function_output_path = A('function_table_output')
        if not gene_output_path and not function_output_path:
            raise ConfigError(
                "Neither a gene nor function output table path was given. One or both is required.")

        self.function_sources = A('function_sources')
        self.use_shared_function_sources = A('shared_function_sources')

        self.preload_genomes = A('preload-genomes') or False

        self.run = run
        self.progress = progress

        descriptions = GenomeDescriptions(args, run=self.run, progress=self.progress)
        descriptions.load_genomes_descriptions(init=False)

        if (self.function_sources and
            self.use_shared_function_sources and
            descriptions.function_annotation_sources_some_genomes_miss):
            raise ConfigError(
                "Of the requested function annotation sources, the following were not run on every "
                f"genome: {', '.join(descriptions.function_annotation_sources_some_genomes_miss)}")

        # Store information on accessing the genomes.
        self.genome_info_dict = {}
        for genome_name, genome_dict in descriptions.internal_genomes_dict.items():
            self.genome_info_dict[genome_name] = genome_info = {}
            genome_info['contigs_db'] = genome_dict['contigs_db_path']
            genome_info['profile_db'] = genome_dict['profile_db_path']
            genome_info['collection_name'] = genome_dict['collection_id']
            genome_info['bin_id'] = genome_dict['bin_id']
            genome_info['function_sources'] = self.function_sources
        for genome_name, genome_dict in descriptions.external_genomes_dict.items():
            self.genome_info_dict[genome_name] = genome_info = {}
            genome_info['contigs_db'] = genome_dict['contigs_db_path']
            genome_info['function_sources'] = self.function_sources

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
                        as_relative_frequency=False,
                        as_synonymous=False,
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
        Get absolute (default) or relative codon or amino acid frequencies from genes or functions
        in one or more genomes.

        See the SingleGenomeCodonUsage docstring for descriptions of each parameter.
        """
        args = {}
        arg_info = inspect.getargvalues(inspect.currentframe())
        for param in arg_info.args:
            if param == 'self':
                continue
            args[param] = arg_info.locals[param]
        frequency_table_generator = self._get_genome_frequency_table(args)
        frequency_dfs = []
        for genome_name in self.genome_info_dict:
            frequency_df = next(frequency_table_generator)

            frequency_df.insert(0, 'genome', genome_name)
            new_index_columns = ['genome'] + list(frequency_df.index.names)
            frequency_df = frequency_df.reset_index().set_index(new_index_columns)

            frequency_dfs.append(frequency_df)
        return pd.concat(frequency_dfs)


    def _get_genome_frequency_table(self, args):
        """This generator yields a frequency table from each genome."""
        for genome_name, genome_info in self.genome_info_dict.items():
            if self.preload_genomes:
                genome_codon_usage = self.genome_codon_usage_dict[genome_name]
            else:
                genome_codon_usage = SingleGenomeCodonUsage(
                    argparse.Namespace(**genome_info))
            frequency_df = genome_codon_usage.get_frequencies(**args)
            yield frequency_df
