#!/usr/bin/env python
# -*- coding: utf-8
"""
A genomic context for codon usage analyses includes gene sequences and their available functional
annotations.
"""

import pandas as pd

from typing import Any
from copy import deepcopy

import anvio.utils as utils
import anvio.terminal as terminal

from argparse import Namespace
from anvio.dbops import ContigsSuperclass
from anvio.ccollections import GetSplitNamesInBins

class GenomicContext:
    """
    Stores information on gene sequences and their functional annotations.

    The genomic context can be automatically set up from a contigs database (external genome) or a
    bin (internal genome) from the contigs database. It is also possible to minimally initialize
    this object without a contigs database and fill out the attributes of a genomic context.

    Attributes
    ==========
    contigs_db : str, None
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

    function_sources : Union[list[str], None]
        Defines the gene annotation function sources that are presented by the genomic context,
        e.g., ['KOfam', 'COG14_FUNCTION']. Note that this can be subset of all available sources in
        the contigs database. A value of None means no sources are made available.

    contig_sequences_dict : dict[str, dict[str, Any]], {}
        Keys are contig names. Values are contig information. Each dictionary value should contain a
        key, 'sequence', that maps to the contig's nucleotide sequence string. This attribute
        represents the 'contig_sequences' table of the contigs database.

    genes_in_contigs_dict : dict[int, dict[str, Any]], {}
        Keys are gene caller IDs. Values are information on the gene call in the context of contigs.
        Each dictionary value should contain attributes, 'contig', 'start', 'stop', 'direction', and
        'call_type', for the name of the contig containing the gene call, the start and stop
        positions of the gene in the contig, the orientation of the gene ('f' or 'r'), and whether
        the gene is coding (1) or not (not 1).

    gene_caller_ids : list[int], []
        List of gene caller IDs in the contigs database corresponding to the keys of
        'genes_in_contigs_dict'.

    gene_function_df : pandas.core.Frame.DataFrame, None
        Gene functional annotations for select 'function_sources'. The table should have columns,
        'source', 'accession', 'name', and 'gene_caller_id'.

    run : anvio.terminal.Run, anvio.terminal.Run()
        Prints run information to the terminal.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        Prints transient progress information to the terminal.

    kwargs : dict[str, Any]
        Keyword arguments provided at initialization.
    """
    def __init__(
        self,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        **kwargs
    ) -> None:
        """
        Parameters
        ==========
        run : anvio.terminal.Run, anvio.terminal.Run()
            Prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            Prints transient progress information to the terminal.

        Keyword arguments
        =================
        Arguments that can be used in setup. Those without descriptions are documented in equivalent
        attributes of the same name.

        contigs_db : str

        profile_db : str

        collection_name : str

        bin_id : str

        function_sources : list[str]
            By default, load gene functional annotations from all sources used in annotation of the
            contigs database. A subset of sources can be made available, e.g., by providing the
            argument, ['KOfam', 'COG14_FUNCTION']. To make no sources available, provide an empty
            list as the argument.
        """
        self.run = run
        self.progress = progress
        self.kwargs = kwargs

        A = lambda x: kwargs[x] if x in kwargs else None

        self.contigs_db = A('contigs_db')

        self.profile_db = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')

        if self.contigs_db:
            self.function_sources = A('function_sources')
            self._load_contigs_db_data()
        else:
            self.contig_sequences_dict: dict[str, dict[str, Any]] = {}
            self.genes_in_contigs_dict: dict[int, dict[str, Any]] = {}
            self.gene_caller_ids: list[int] = []
            self.function_sources: list[str] = []
            self.gene_function_df: pd.DataFrame = None

    def _load_contigs_db_data(self):
        """Load gene data from the contigs database."""
        utils.is_contigs_db(self.contigs_db)

        args = Namespace(**self.kwargs)
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
        if self.function_sources is None:
            self.function_sources = sorted(contigs_super.gene_function_call_sources)
            if not self.function_sources:
                self.run.warning(
                    "No 'function_sources' keyword argument was provided, indicating that all gene "
                    "function annotation sources in the contigs database should by loaded. "
                    "However, the contigs database has not been annotated by any sources."
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
