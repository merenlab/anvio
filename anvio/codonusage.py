# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module for codon usage analyses at the levels of genes, groups of genes, and genomes"""


import inspect
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

pp = terminal.pretty_print

translation_table = constants.AA_to_codons
nonstop_translation_table = {aa: co for aa, co in translation_table.items() if aa != 'STP'}
decoding_table = constants.codon_to_AA

# CAI is the Codon Adaptation Index of Sharp and Li (1987).
# Delta is from Ran and Higgs (2012, Eq. 6).
cub_metrics = ['cai', 'delta']
# The following CUB metrics rely upon comparison to a set of reference genes.
reference_dependent_cub_metrics = ['cai', 'delta']

# For CUB, remove single codon amino acids, stop codons, and codons with "N" nucleotide which were
# recorded as "None".
ignored_decodings_in_cub = ['Met', 'Trp', 'STP']
ignored_codons_in_cub = [None]
for decoding in ignored_decodings_in_cub:
    ignored_codons_in_cub += translation_table[decoding]
synonymous_codon_dict = {aa: co for aa, co in translation_table.items()
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
                f"{pp(len(skipped_noncoding_gene_caller_ids))} of {pp(len(self.gene_caller_ids))} "
                "genes were non-coding and not added to the codon frequency table.")


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
            the returned table contain, respectively, function annotation sources, function
            accessions, function names, and gene caller IDs. There is a row for each gene/function
            combination in the table, and each row for the same gene contains the same frequency
            values. When this argument is True, use all available functional annotation sources in
            the SingleGenomeCodonUsage object. When this argument is a string, select the source
            given by the string, e.g., 'KOfam', 'COG20_FUNCTION', 'Pfam'. When this argument is an
            iterable, select a subset of sources given by its strings, e.g., ['KOfam',
            'COG20_FUNCTION']. When 'KEGG_BRITE' is in the argument, in the absence of
            `function_names` narrowing the scope of the inquiry, all possible depths of each BRITE
            hierarchy in the data are returned, e.g., 'Ribosome>>>Ribosomal proteins' and the more
            general 'Ribosome' would each have rows in the output table. By default None.
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
            machinery>>>Prokaryotic type>>>Bacterial type>>>RNA polymerase']. This parameter can be
            used alongside `gene_caller_ids` and `function_accessions`. By default None.
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
            the first length filter.
        drop_amino_acids : iterable, optional
            Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala, and
            STP for stop codons). If `synonymous` is True, the `drop_amino_acids` default rather
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
        >>> self.get_frequencies(relative=True)

        Return the relative frequencies of genes ≥300 codons in length.
        >>> self.get_frequencies(relative=True, gene_min_codons=300)

        Return the average relative frequencies of genes ≥300 codons in length.
        >>> self.get_frequencies(relative=True, average_all=True, gene_min_codons=300)

        Return the synonymous (per-amino acid) relative frequencies of genes. Remove genes <300
        codons, then remove amino acids with <5 codons in ≥90% of remaining genes, then remove genes
        with <300 codons remaining.
        >>> self.get_frequencies(
            synonymous=True, gene_min_codons=300, min_amino_acids=5, min_gene_fraction=0.9)

        Return the synonymous (per-amino acid) relative frequencies of KEGG KOfam and BRITE
        functions. This command removes genes <300 codons, then removes amino acids with <5 codons
        in ≥90% of remaining genes, then removes genes with <300 codons remaining, then sums gene
        frequencies in functions, then calculates synonymous relative frequencies in functions.
        >>> self.get_frequencies(
            from_function_sources=['KOfam', 'KEGG_BRITE'],
            return_functions=True,
            synonymous=True,
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
        if return_amino_acids and synonymous:
            raise ConfigError("The argument `synonymous` should only be True when "
                              "`return_amino_acids` is also True, as `synonymous` returns "
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
            if synonymous:
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
            self.run.info_single(f"{pp(len(gene_codon_frequency_df))} of "
                                 f"{pp(len(self.gene_caller_ids))} CDS selected from the genome",
                                 nl_after=1)
        else:
            self.run.info_single(f"{pp(len(gene_codon_frequency_df))} CDS in the genome",
                                 nl_after=1)

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
            self.run.info_single(f"{pp(len(gene_codon_frequency_df))} CDS remaining after codon "
                                 "count filters",
                                 nl_after=1)

        if min_amino_acids:
            if gene_min_codons and min_codon_filter != 'remaining':
                min_gene_length_message = '≥' + pp(str(gene_min_codons)) + ' codon'
            else:
                min_gene_length_message = ''
            dynamically_dropped_amino_acids = set()
            for codon, amino_acid in decoding_table.items():
                if codon not in gene_codon_frequency_df.columns:
                    dynamically_dropped_amino_acids.add(amino_acid)
            self.run.warning(
                "Codons for the following amino acids were dropped as they did not meet the "
                f"threshold of {pp(min_amino_acids)} codons in {min_gene_fraction * 100}% of "
                f"{min_gene_length_message} CDS: {', '.join(dynamically_dropped_amino_acids)}")

        # GET OUTPUTS WITH NO CONSIDERATION OF FUNCTIONS
        ################################################
        get_table = lambda method: method(gene_codon_frequency_df,
                                          label_amino_acids=label_amino_acids,
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
        # Call `_get_codon_frequency_table` on the result of the method to add amino acid names to
        # the header or change codon to amino acid frequencies (change columns). Note in this and
        # further `get_table` functions that `_get_codon_frequency_table` is used for the purpose of
        # its decorator functions.
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df).to_frame('all').T.rename_axis('gene_caller_ids'),
            label_amino_acids=label_amino_acids,
            output_amino_acids=return_amino_acids)
        ### Absolute frequencies
        if (not relative and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_frequency_table)
        ### Relative frequencies
        if (relative and
            not synonymous and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df, label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            relative and
            synonymous and
            not function_sources and
            sum_genes):
            return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

        # Codon results averaged across genes:
        get_table = lambda method: self._get_frequency_table(
            method(gene_codon_frequency_df).to_frame('all').T.rename_axis('gene_caller_ids'),
            label_amino_acids=label_amino_acids)
        ### Absolute frequencies
        if (not return_amino_acids and
            not relative and
            not function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (not return_amino_acids and
            relative and
            not synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: method(
            gene_codon_frequency_df, label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            relative and
            synonymous and
            not function_sources and
            average_genes):
            return get_table(self._get_average_synonymous_codon_rel_frequency_table)

        # Amino acid results averaged across genes:
        # Gene amino acid frequencies are first calculated from codon frequencies before averaging
        # across genes.
        get_table = lambda method: method(
            gene_codon_frequency_df, output_amino_acids=return_amino_acids).to_frame(
                'all').T.rename_axis('gene_caller_ids')
        ### Absolute frequencies
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
            (gene_function_codon_frequency_df.groupby(['source', 'accession', 'name']).sum()
             if return_functions else
             gene_function_codon_frequency_df.sort_values(['source', 'accession', 'name'])),
            label_amino_acids=label_amino_acids,
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
            gene_function_codon_frequency_df.reset_index().groupby('source').apply(
                lambda source_df: source_df.drop_duplicates(
                    subset='gene_caller_id', ignore_index=True)).set_index(
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
        if (not relative and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_frequency_table)
        ### Relative frequencies
        if (relative and
            not synonymous and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
                label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            relative and
            synonymous and
            function_sources and
            sum_genes):
            return get_table(self._get_summed_synonymous_codon_rel_frequency_table)

        # Codon results averaged across genes annotated by a function source:
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('source').apply(method),
            label_amino_acids=label_amino_acids)
        ### Absolute frequencies
        if (not return_amino_acids and
            not relative and
            function_sources and
            sum_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
        if (not return_amino_acids and
            relative and
            not synonymous and
            function_sources and
            average_genes):
            return get_table(self._get_average_rel_frequency_table)
        ### Synonymous relative frequencies
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df.groupby('source').apply(method).droplevel(1),
            label_amino_acids=label_amino_acids)
        if (not return_amino_acids and
            relative and
            synonymous and
            function_sources and
            average_genes):
            return get_table(self._get_average_synonymous_codon_rel_frequency_table)

        # Amino acid results averaged across genes annotated by a function source:
        get_table = lambda method: self._get_frequency_table(
            gene_function_codon_frequency_df, output_amino_acids=return_amino_acids).groupby(
                'source').apply(method)
        ### Absolute frequencies
        if (return_amino_acids and
            not relative and
            function_sources and
            average_genes):
            return get_table(self._get_average_frequency_table)
        ### Relative frequencies
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


    def _select_genes(
        self, gene_caller_ids, function_accession_dict, function_name_dict, expect_functions):
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
            if function_source == 'KEGG_BRITE':
                self.run.warning("Nota bene: KEGG BRITE accessions stored in anvi'o are for "
                                 "hierarchies as a whole, not categories of the hierarchy. Most "
                                 "hierarchies do not have category accessions. So all genes in the "
                                 "selected hierarchies are being analyzed.")
            if function_source not in self.function_sources:
                unrecognized_sources.append(function_source)
            for function_accession in function_accessions:
                index_keys.append((function_source, function_accession))
        if unrecognized_sources:
            raise ConfigError("The following annotation sources in `function_accessions` were not "
                              f"found as having annotated the genome: {unrecognized_sources}")
        gene_function_df = self.gene_function_df.set_index(['source', 'accession'])
        if expect_functions:
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
                    raise ConfigError("The following requested function accessions are missing "
                                      f"from the genome: {missing_key_message}")
        else:
            index_keys = gene_function_df.index.intersection(index_keys)
        select_gene_caller_ids += gene_function_df.loc[index_keys]['gene_caller_id'].tolist()

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
        if expect_functions:
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
                    raise ConfigError("The following requested function names are missing from the "
                                      f"genome: {missing_key_message}")
        else:
            index_keys = gene_function_df.index.intersection(index_keys)
        select_gene_caller_ids += gene_function_df.loc[index_keys]['gene_caller_id'].tolist()

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
                min_codons = kwargs['min_codons']
            except KeyError:
                filter_input_codon_count = False
                min_codons = 0
            if filter_input_codon_count and min_codons:
                frequency_df = frequency_df[frequency_df.sum(axis=1) >= min_codons]
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
                if len(codon_frequency_df.columns) == 0:
                    codon_frequency_df = codon_frequency_df.drop(codon_frequency_df.index)
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
                if len(codon_frequency_df.columns) == 0:
                    codon_frequency_df = codon_frequency_df.drop(codon_frequency_df.index)
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
                min_codons = kwargs['min_codons']
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
        summed_frequency_series = self._get_summed_frequency_table(frequency_df, **kwargs)
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
        summed_codon_frequency_df = \
            summed_codon_frequency_series.to_frame('all').T.rename_axis('gene_caller_ids')

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
        average_codon_frequency_df = average_codon_frequency_series.to_frame(
            'all').T.rename_axis('gene_caller_ids')

        average_synonymous_codon_rel_frequency_df = \
            self._get_synonymous_codon_rel_frequency_table(average_codon_frequency_df)

        return average_synonymous_codon_rel_frequency_df


    def get_codon_usage_bias(self,
                             metrics=None,
                             from_function_sources=None,
                             gene_caller_ids=None,
                             function_accessions=None,
                             function_names=None,
                             expect_functions=False,
                             omnibias=False,
                             reference_gene_caller_ids=None,
                             reference_function_accessions=None,
                             reference_function_names=None,
                             expect_reference_functions=False,
                             gene_min_codons=0,
                             function_min_codons=0,
                             min_codon_filter='both',
                             drop_amino_acids=['STP'],
                             min_amino_acids=0,
                             min_gene_fraction=1,
                             query_min_analyzed_codons=100,
                             reference_exclude_amino_acid_count=5,
                             reference_min_codons=100):
        """
        Get codon usage bias (CUB) of genes or functions.

        Parameters
        ----------
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
            'COG20_FUNCTION']. When 'KEGG_BRITE' is in the argument, in the absence of
            `function_names` narrowing the scope of the inquiry, all possible depths of each BRITE
            hierarchy in the data are returned, e.g., 'Ribosome>>>Ribosomal proteins' and the more
            general 'Ribosome' would each have rows in the output table. By default None.
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
            type>>>Bacterial type>>>RNA polymerase']. This parameter can be used alongside
            `gene_caller_ids` and `function_accessions`. By default None.
        expect_functions : bool, optional
            If True (default False), an error will be raised if any given `function_accessions` or
            `function_names` are not annotated in the input genome.
        omnibias : bool, optional
            If True (default False), use every gene or function as a separate reference rather than
            defining a set of reference genes or functions. The resulting table of gene x gene (or
            function x function) CUB values is like a distance matrix of the similarity of gene
            codon compositions.
        reference_gene_caller_ids : iterable, optional
            Genes with the given IDs are selected for the reference set. The parameter can be used
            alongside `reference_function_accessions` and `reference_function_names`. By default
            None.
        reference_function_accessions : dict, optional
            Genes annotated with the given function accessions are selected for the reference set.
            The argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function accessions. This
            parameter can be used alongside `reference_gene_caller_ids` and `function_names`. Note
            that 'KEGG_BRITE' does not use individual function accessions but overarching hierarchy
            accessions that include multiple functions. By default None.
        reference_function_names : dict, optional
            Genes annotated with the given function names are selected for the reference set. The
            argument must be a dict keyed by function annotation source (e.g., 'KOfam',
            'COG20_FUNCTION', 'Pfam') and with values being lists of function names. Unlike
            `function_accessions`, 'KEGG_BRITE' may be used as a source, with names
            (categorizations) given to an arbitrary level of the hierarchy, such as
            ['Ribosome>>>Ribosomal proteins', 'Transcription machinery>>>Prokaryotic
            type>>>Bacterial type>>>RNA polymerase']. This parameter can be used alongside
            `gene_caller_ids` and `function_accessions`. By default, {'KEGG_BRITE':
            ['Ribosome>>>Ribosomal proteins']}.
        expect_reference_functions : bool, optional
            If True (default False), an error will be raised if any given
            `reference_function_accessions` or `reference_function_names` are not annotated in the
            input genome.

        Additional Parameters
        ---------------------
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

        There are two filters, 'reference_exclude_amino_acid_count' and 'reference_min_codons', that
        are applied to the reference gene set.

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
            filters to ensure that total codon frequencies always meet the minimum codon threshold.
            'both' is needed as an option in addition to 'remaining' so dynamic codon filtering by
            `--exclude-amino-acid-count/fraction` operates on genes that passed the first length
            filter.
        drop_amino_acids : iterable, optional
            Remove codons that decode the given amino acids (use three-letter codes, e.g., Ala). By
            default, the value, ['STP'], excludes stop codons from CUB calculations.
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
        query_min_analyzed_codons : int, optional
            Only allow CUB to calculated for a query if has at least this number of synonymous
            codons that will be analyzed. For reference-dependent CUB metrics, analyzed codons are
            those with reference compositions. By default 100.
        reference_exclude_amino_acid_count : int, optional
            Exclude codons for amino acids with fewer than this many codons in the set of reference
            genes. This does not apply in `omnibias` mode. By default 5.
        reference_min_codons : int, optional
            Only allow CUB to be calculated using a reference if it has at least this number of
            codons. This filter applies after excluding codons for individual amino acids using
            `reference_exclude_amino_acid_count`. By default 100.

        Returns
        -------
        dict of pandas.core.frame.DataFrame objects
            Each dict key is the name of a CUB metric (from `metrics`), and each value is a
            corresponding table of CUB data.
        """
        metrics = self._establish_cub_metrics(metrics)
        reference_metrics = []
        nonreference_metrics = []
        for metric in metrics:
            if metric in reference_dependent_cub_metrics:
                reference_metrics.append(metric)
            else:
                nonreference_metrics.append(metric)

        # Check that proper arguments were provided given the metrics.
        if reference_metrics:
            if (omnibias and
                (reference_gene_caller_ids or
                 reference_function_accessions or
                 reference_function_names or
                 expect_reference_functions or
                 reference_exclude_amino_acid_count or
                 reference_min_codons)):
                raise ConfigError(
                    "Omnibias mode cannot be used when defined gene/function references are also "
                    "used. The following arguments are only relevant to defined references: "
                    "`reference_gene_caller_ids`, `reference_function_accessions`, "
                    "`reference_function_names`, `expect_reference_functions`, "
                    "`reference_exclude_amino_acid_count`, and `reference_min_codons`.")
        else:
            if (omnibias or
                reference_gene_caller_ids or
                reference_function_accessions or
                reference_function_names or
                expect_reference_functions or
                reference_exclude_amino_acid_count or
                reference_min_codons):
                raise ConfigError(
                    "The provided CUB metrics do not involve comparison of gene/function codon "
                    "compositions. The following arguments are only relevant to "
                    "reference-dependent metrics: `omnibias`, `reference_gene_caller_ids`, "
                    "`reference_function_accessions`, `reference_function_names`, "
                    "`expect_reference_functions`, `reference_exclude_amino_acid_count`, and "
                    "`reference_min_codons`.")

        # Get a reference codon composition when using reference-dependent metrics and not in
        # omnibias mode.
        if not reference_metrics or omnibias:
            reference_codon_frequency_df = None
        else:
            # The default reference genes are KOfams classified as ribosomal proteins in BRITE.
            if reference_function_names is None:
                if 'KEGG_BRITE' not in self.function_sources:
                    raise ConfigError(
                        f"Reference-dependent metrics ({', '.join(reference_metrics)}) were "
                        "requested without defined reference genes. By default, reference genes "
                        "are KEGG KOfams classified as ribosomal proteins in BRITE. However, "
                        "'KEGG_BRITE' is not among the function annotation sources run on the "
                        "genome. This can be rectified by rerunning `anvi-run-kegg-kofams`.")
                reference_function_names = {'KEGG_BRITE': ['Ribosome>>>Ribosomal proteins']}

            reference_codon_frequency_df = self.get_frequencies(
                gene_caller_ids=reference_gene_caller_ids,
                function_accessions=reference_function_accessions,
                function_names=reference_function_names,
                expect_functions=expect_reference_functions,
                sum_genes=True,
                drop_amino_acids=ignored_decodings_in_cub)

            # In the absence of a set of reference genes for the genome, generate a DataFrame with
            # an index and header but no data.
            if len(reference_codon_frequency_df) == 0:
                found_reference_genes = False
            else:
                found_reference_genes = True

            # Remove codons for rarer amino acids (columns) from the reference set.
            if found_reference_genes and reference_exclude_amino_acid_count:
                reference_drop_amino_acids = []
                reference_drop_codons = []
                for amino_acid, codons in translation_table.items():
                    if (reference_codon_frequency_df[codons].sum().sum()
                        < reference_exclude_amino_acid_count):
                        reference_drop_amino_acids.append(amino_acid)
                        reference_drop_codons += codons
                reference_codon_frequency_df = \
                    reference_codon_frequency_df.drop(reference_drop_codons, axis=1)

                self.run.warning(
                    "The following amino acids do not meet the threshold of "
                    f"{reference_exclude_amino_acid_count} codons in the reference gene set, "
                    "and so the codons were excluded from the reference: "
                    f"{', '.join(reference_drop_codons)}")

            # If there aren't enough codons total remaining in the reference set, remove the data
            # from the DataFrame, leaving just an index and header.
            if found_reference_genes and reference_min_codons:
                total_codon_count = reference_codon_frequency_df.sum().sum()
                if total_codon_count < reference_min_codons:
                    reference_codon_frequency_df = reference_codon_frequency_df.drop('all')

                    self.run.warning(
                        f"The reference gene set has {pp(total_codon_count)} codons, not "
                        f"meeting the minimum codon threshold of {pp(reference_min_codons)}.")

            absent_amino_acids = []
            for amino_acid, codons in translation_table.items():
                if len(codons) == 1:
                    continue
                try:
                    if reference_codon_frequency_df[codons].sum().sum() == 0:
                        absent_amino_acids.append(amino_acid)
                except KeyError:
                    continue
                reference_codon_frequency_df = reference_codon_frequency_df.drop(codons, axis=1)
            if absent_amino_acids:
                self.run.warning("Codons for the following multi-codon amino acids were not "
                                 f"detected in the reference set: {', '.join(absent_amino_acids)}")

            # If the absolute frequency table has no data, then the synonymous relative frequency
            # table will also not have data.
            reference_synonymous_codon_rel_frequency_series = \
                self._get_synonymous_codon_rel_frequency_table(
                    reference_codon_frequency_df).squeeze()

            if len(reference_synonymous_codon_rel_frequency_series) == 0:
                self.run.warning("No reference codon composition was established.")

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
            min_amino_acids=min_amino_acids,
            min_gene_fraction=min_gene_fraction)

        if omnibias:
            reference_codon_frequency_df = query_codon_frequency_df

        if 'delta' in metrics:
            # Reference codon weights in the computation of delta involve comparison of the
            # reference codon set to the overall codon set, which here is taken from all genes in
            # the input dataset.
            overall_codon_frequency_df = self.get_frequencies(sum_genes=True)
        else:
            overall_codon_frequency_df = None

        cub_table_dict = {}

        for metric in reference_metrics:
            if metric == 'cai':
                cub_df = self._get_cai_table(
                    query_codon_frequency_df,
                    reference_codon_frequency_df,
                    query_min_codons=query_min_analyzed_codons,
                    reference_min_codons=reference_min_codons)
                if omnibias:
                    cub_df.columns = reference_codon_frequency_df.index
                else:
                    cub_df.columns = ['CAI']
            elif metric == 'delta':
                cub_df = self._get_delta_table(
                    query_codon_frequency_df,
                    reference_codon_frequency_df,
                    overall_codon_frequency_df,
                    query_min_codons=query_min_analyzed_codons,
                    reference_min_codons=reference_min_codons)
                if omnibias:
                    cub_df.columns = reference_codon_frequency_df.index
                else:
                    cub_df.columns = ['Delta']

            cub_table_dict[metric] = cub_df

            # Print the number of query-reference comparisons that were thrown out.
            if omnibias:
                possible_comparison_count = \
                    len(query_codon_frequency_df) * len(reference_codon_frequency_df) / 2
            else:
                possible_comparison_count = len(query_codon_frequency_df)
            unperformed_comparison_count = cub_df.isna().sum().sum()
            self.run.warning(f"{pp(unperformed_comparison_count)} of "
                             f"{pp(possible_comparison_count)} query-reference comparisons did "
                             "not meet the minimum codon threshold in either query or reference "
                             "and so did not yield a CUB value.")

        # Calculate CUB using metrics that do not depend on a reference codon composition.
        for metric in nonreference_metrics:
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


    def _get_cai_table(self,
                       query_codon_frequency_df,
                       reference_codon_frequency_df,
                       query_min_codons=0,
                       reference_min_codons=0):
        """Get a table of CAI (Sharp and Li, 1987) values for each query x reference comparison.

        Calculation of CAI:
        reference_codon_weight =
            ln(reference_codon_frequency / reference_max_synonymous_codon_frequency)
        weighted_codon_count = ∑(codon_frequency * reference_codon_weight)
        CAI = exp(weighted_codon_count / codon_count)

        CAI is maximum (1) when all of the codons in the query are the most abundant synonymous
        codons in the reference, and minimum (0-1) when all of the codons in the query are the least
        abundant synonymous codons in the reference.

        Parameters
        ----------
        query_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per item (gene/function) to be compared to each reference
            codon composition and a column per codon.
        reference_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per reference composition and a column per codon.
        query_min_codons : int, optional
            A row of the query table must contain at least this number of codons with a reference
            codon weight for CAI to be calculated. By default 0.
        reference_min_codons : int, optional
            A row of the reference table must contain at least this number of synonymous codons to
            be used as a reference. By default 0.

        Returns
        -------
        pandas.core.frame.DataFrame
            This CUB table has the same row index as the input query table and a column per input
            reference. With a single reference codon composition, this table has a single column.
        """
        reference_weight_series_list = self._get_cai_reference_weights(
            reference_codon_frequency_df, reference_min_codons=reference_min_codons)

        cai_rows = []
        for _, query_codon_frequency_series in query_codon_frequency_df.iterrows():
            # Check that the query row meets the minimum codon threshold.
            if query_codon_frequency_series.sum() < query_min_codons:
                cai_rows.append([np.nan] * len(reference_weight_series_list))
                continue

            cai_row = []
            for reference_weight_series in reference_weight_series_list:
                query_codons_with_reference_frequency_series = query_codon_frequency_series[
                    reference_weight_series.index]

                # Check that the number of query codons with a reference meets the minimum codon
                # threshold.
                if query_codons_with_reference_frequency_series.sum() < query_min_codons:
                    cai_row.append(np.nan)
                    continue

                weighted_codon_count = query_codons_with_reference_frequency_series.dot(
                    reference_weight_series)
                cai_row.append(np.exp(
                    weighted_codon_count / query_codons_with_reference_frequency_series.sum()))
            cai_rows.append(cai_row)

        # The returned table has a column per reference and default integer column names.
        return pd.DataFrame(cai_rows, index=query_codon_frequency_df.index)


    def _get_cai_reference_weights(self, codon_frequency_df, reference_min_codons=100):
        """Get reference codon weights needed for calculation of CAI.

        A CAI reference weight is calculated as the natural log of the frequency of a codon divided
        by the maximum frequency of synonymous codons.

        Parameters
        ----------
        codon_frequency_df : pandas.core.frame.DataFrame
            This is the reference codon frequency table.
        reference_min_codons : int, optional
            If a row of `codon_frequency_df` contains fewer than this number of synonymous codons,
            then the returned series is empty. By default 0.

        Returns
        -------
        list of pandas.core.series.Series
            This is a list of reference weight series for each row of the input table. Each series
            has a codon index ordered like the input table. If there are fewer reference codons than
            `reference_min_codons`, then the series will be empty. If a group of synonymous codons
            has zero frequency, then these codons are dropped in the series. Hence the returned
            series can be of unequal length.
        """
        reference_weight_series_list = []
        for _, codon_frequency_series in codon_frequency_df.iterrows():
            reference_weight_dict = {}
            retained_frequency = 0
            for codons in synonymous_codon_dict.values():
                try:
                    synonymous_codon_frequency_series = codon_frequency_series[codons]
                except KeyError:
                    # Treat one missing codon in the table as the absence of all synonymous codons.
                    continue
                max_frequency = max(synonymous_codon_frequency_series)

                retained_frequency += sum(codon_frequency_series[codons])

                for codon in codons:
                    codon_frequency = codon_frequency_series[codon]
                    if codon_frequency == 0:
                        # Do not record reference weights for codons with zero frequency, as the log
                        # is undefined.
                        continue
                    reference_weight_dict[codon] = np.log(codon_frequency / max_frequency)

            # Check for an invalid input row with fewer codons than allowed in a reference.
            if retained_frequency < reference_min_codons:
                reference_weight_dict = {codon: np.nan for codon in reference_weight_dict}

            # Codons without reference weights are absent from the series index.
            reference_weight_series = pd.Series(reference_weight_dict)
            reference_weight_series = reference_weight_series[reference_weight_series.index.intersection(codon_frequency_df.columns)]
            reference_weight_series_list.append(reference_weight_series)

        return reference_weight_series_list


    def _get_delta_table(self,
                         query_codon_frequency_df,
                         reference_codon_frequency_df,
                         overall_codon_frequency_df,
                         query_min_codons=0,
                         reference_min_codons=0):
        """
        Get a table of delta (Ran and Higgs, 2012, Eq. 6) values for each query x reference
        comparison.

        Calculation of delta:
        reference_codon_weight = ln(
            reference_codon_synonymous_relative_frequency /
            overall_codon_synonymous_relative_frequency)
        weighted_codon_count = ∑(codon_frequency * reference_codon_weight)
        delta = weighted_codon_count / codon_count

        Delta and CAI differ in the calculation of reference codon weights. Whereas CAI compares
        reference codon frequency to the maximum frequency of synonymous reference codons, delta
        compares reference codon synonymous frequency to the overall synonymous frequency of the
        codon in the genome. When the query and reference sets are the same, delta is a likelihood
        ratio evaluating the distinctness of the set from the genome as a whole. Delta ranges from
        (-∞, +∞), with more positive values being more similar to the reference.

        Parameters
        ----------
        query_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per item (gene/function) to be compared to each reference
            codon composition and a column per codon.
        reference_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a row per reference composition and a column per codon.
        overall_codon_frequency_df : pandas.core.frame.DataFrame
            This frequency table has a single row of summed frequencies across all genes in the
            genome, with a column per codon.
        query_min_codons : int, optional
            A row of the query table must contain at least this number of codons with a reference
            codon weight for delta to be calculated. By default 0.
        reference_min_codons : int, optional
            A row of the reference table must contain at least this number of synonymous codons to
            be used as a reference. By default 0.

        Returns
        -------
        pandas.core.frame.DataFrame
            This CUB table has the same row index as the input query table and a column per input
            reference. With a single reference codon composition, this table has a single column.
        """
        reference_weight_series_list = self._get_delta_reference_weights(
            reference_codon_frequency_df,
            overall_codon_frequency_df,
            reference_min_codons=reference_min_codons)

        delta_rows = []
        for _, query_codon_frequency_series in query_codon_frequency_df.iterrows():
            # Check that the query row meets the minimum codon threshold.
            if query_codon_frequency_series.sum() < query_min_codons:
                delta_rows.append([np.nan] * len(reference_weight_series_list))
                continue

            delta_row = []
            for reference_weight_series in reference_weight_series_list:
                query_codons_with_reference_frequency_series = query_codon_frequency_series[
                    reference_weight_series.index]

                # Check that the number of query codons with a reference meets the minimum codon
                # threshold.
                if query_codons_with_reference_frequency_series.sum() < query_min_codons:
                    delta_row.append(np.nan)
                    continue

                weighted_codon_count = query_codons_with_reference_frequency_series.dot(
                    reference_weight_series)
                delta_row.append(
                    weighted_codon_count / query_codons_with_reference_frequency_series.sum())
            delta_rows.append(delta_row)

        # The returned table has a column per reference and default integer column names.
        return pd.DataFrame(delta_rows, index=query_codon_frequency_df.index)


    def _get_delta_reference_weights(self,
                                     reference_codon_frequency_df,
                                     overall_codon_frequency_df,
                                     reference_min_codons=100):
        """Get reference codon weights needed for calculation of delta.

        A delta reference weight is calculated as the natural log of the synonymous relative
        frequency of a codon in the reference set divided by that of a codon in the genome as a
        whole.

        Parameters
        ----------
        reference_codon_frequency_df : pandas.core.frame.DataFrame
            This is the reference codon frequency table.
        overall_codon_frequency_df : pandas.core.frame.DataFrame
            This is the summed frequency table, with one row, of the genome as a whole.
        reference_min_codons : int, optional
            If a row of `reference_codon_frequency_df` contains fewer than this number of synonymous
            codons, then the returned series is empty. By default 0.

        Returns
        -------
        list of pandas.core.series.Series
            This is a list of reference weight series for each row of the reference table. Each
            series has a codon index ordered like the input table. If there are fewer reference
            codons than `reference_min_codons`, then the series will be empty. If a group of
            synonymous codons has zero frequency, then these codons are dropped in the series. Hence
            the returned series can be of unequal length.
        """
        reference_synonymous_codon_rel_frequency_df = \
            self._get_synonymous_codon_rel_frequency_table(reference_codon_frequency_df)
        overall_synonymous_codon_rel_frequency_series = \
            self._get_synonymous_codon_rel_frequency_table(overall_codon_frequency_df).squeeze()

        reference_weight_series_list = []
        reference_codon_frequency_rows = \
            reference_codon_frequency_df.iterrows()
        reference_synonymous_codon_rel_frequency_rows = \
            reference_synonymous_codon_rel_frequency_df.iterrows()
        while reference_codon_frequency_rows:
            _, reference_codon_frequency_series = next(reference_codon_frequency_rows)
            _, reference_synonymous_codon_rel_frequency_series = next(
                reference_synonymous_codon_rel_frequency_rows)

            # Divide reference by overall synonymous relative frequencies. No codons should be
            # present in the overall but not the reference series, but drop such codons yielding NaN
            # just in case.
            normalized_synonymous_codon_rel_frequency_series = \
                reference_synonymous_codon_rel_frequency_series.div(
                    overall_synonymous_codon_rel_frequency_series).dropna()
            # Do not record reference weights for codons with zero frequency, as the log is
            # undefined.
            normalized_synonymous_codon_rel_frequency_series = \
                normalized_synonymous_codon_rel_frequency_series[
                    normalized_synonymous_codon_rel_frequency_series > 0]

            # Ignore references with fewer than the threshold number of codons.
            retained_frequency = reference_codon_frequency_series[
                normalized_synonymous_codon_rel_frequency_series.index].sum()
            if retained_frequency < reference_min_codons:
                reference_weight_series = pd.Series()
            else:
                reference_weight_series = np.log(normalized_synonymous_codon_rel_frequency_series)

            reference_weight_series_list.append(reference_weight_series)

        return reference_weight_series_list


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
                        expect_functions=False,
                        relative=False,
                        synonymous=False,
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
        Get absolute (default) or relative codon or amino acid frequencies from genes or functions
        in one or more genomes.

        See the SingleGenomeCodonUsage.get_frequencies docstring for descriptions of each parameter.
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
            self.run.info("Genome", genome_name)
            frequency_df = next(frequency_table_generator)

            frequency_df.insert(0, 'genome', genome_name)
            new_index_columns = ['genome'] + frequency_df.index.names
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


    def get_codon_usage_bias(self,
                             metrics=None,
                             from_function_sources=None,
                             function_accessions=None,
                             function_names=None,
                             expect_functions=False,
                             omnibias=False,
                             reference_function_source=None,
                             reference_function_accessions=None,
                             reference_function_names=None,
                             expect_reference_functions=False,
                             gene_min_codons=0,
                             function_min_codons=0,
                             min_codon_filter='both',
                             drop_amino_acids=None,
                             min_amino_acids=0,
                             min_gene_fraction=1,
                             query_min_analyzed_codons=100,
                             reference_exclude_amino_acid_count=5,
                             reference_min_codons=100):
        """
        Get CUB from genes or functions in one or more genomes.

        See the SingleGenomeCodonUsage.get_frequencies docstring for descriptions of each parameter.
        """
        args = {}
        arg_info = inspect.getargvalues(inspect.currentframe())
        for param in arg_info.args:
            if param == 'self':
                continue
            args[param] = arg_info.locals[param]
        cub_table_dict_generator = self._get_cub_table_dict(args)
        cub_table_master_dict = {}
        for genome_name in self.genome_info_dict:
            self.run.info("Genome", genome_name)
            cub_table_master_dict[genome_name] = next(cub_table_dict_generator)

        return cub_table_master_dict


    def _get_genome_frequency_table(self, args):
        """This generator yields a CUB table dict from each genome."""
        for genome_name, genome_info in self.genome_info_dict.items():
            if self.preload_genomes:
                genome_codon_usage = self.genome_codon_usage_dict[genome_name]
            else:
                genome_codon_usage = SingleGenomeCodonUsage(
                    argparse.Namespace(**genome_info))
            cub_table_dict = genome_codon_usage.get_codon_usage_bias(**args)

            yield cub_table_dict
