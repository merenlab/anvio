#!/usr/bin/env python
# -*- coding: utf-8
"""
Get raw codon frequencies.
"""

import pandas as pd

from collections import Counter

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.codonusage.genomiccontext import GenomicContext

pp = terminal.pretty_print

class GenomeCodonFrequencies:
    """
    Stores raw codon frequencies of gene sequences.

    The codon frequency table can be automatically set up from a genomic context.

    Without a genomic context, the codon frequency table can be provided manually, which can be
    useful, for example, in handling codon frequencies of arbitrary sequences.

    Attributes
    ==========
    genomic_context : anvio.codonusage.genomiccontext.GenomicContext, None
        Source of genes for which codon frequencies are calculated.

    gene_codon_frequency_df : pandas.core.Frame.DataFrame, None
        Table of codon frequencies. The index column is 'gene_caller_id'. Columns have codon
        headers, with a column for each of the 64 possible codons.

    noncoding_gene_count : int, 0
        Number of non-coding genes in the genomic context.

    run : anvio.terminal.Run, anvio.terminal.Run()
        Prints run information to the terminal.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        Prints transient progress information to the terminal.
    """
    def __init__(
        self,
        genomic_context: GenomicContext = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Parameters
        ==========
        genomic_context : anvio.codonusage.genomiccontext.GenomicContext, None
            Used as source of genes for which codon frequencies are calculated.

        run : anvio.terminal.Run, anvio.terminal.Run()
            Prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            Prints transient progress information to the terminal.
        """
        self.genomic_context = genomic_context
        self.run = run
        self.progress = progress

        if self.genomic_context:
            self.set_from_genomic_context()
        else:
            self.gene_codon_frequency_df: pd.DataFrame = None
            self.noncoding_gene_count: int = 0

    def set_from_genomic_context(self):
        """Set the codon frequency table from genomic context."""
        if self.genomic_context is None:
            return

        self.progress.new("Fetching genomic codon frequency data")
        self.progress.update("...")

        gene_codon_freqs = []
        skipped_noncoding_gcids = []
        coding_gcids = []
        genome_gcids = self.genomic_context.gene_caller_ids
        for gcid in genome_gcids:
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
            self.run.info_single(
                f"{pp(self.noncoding_gene_count)} of {pp(len(genome_gcids))} genes were "
                "non-coding and not added to the codon frequency table."
            )
