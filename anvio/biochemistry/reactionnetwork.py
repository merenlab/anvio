# -*- coding: utf-8
# pylint: disable=line-too-long
"""Generate a metabolic reaction network from gene annotations."""

import os
import re
import glob
import json
import math
import time
import shutil
import hashlib
import tarfile
import zipfile
import argparse
import fractions
import functools
import collections
import numpy as np
import pandas as pd
import multiprocessing as mp

from argparse import Namespace
from typing import Dict, List, Set, Tuple

import anvio.utils as utils
import anvio.dbinfo as dbinfo
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio import DEBUG, __file__ as ANVIO_PATH, __version__ as VERSION
from anvio.dbops import ContigsDatabase, PanDatabase, ContigsSuperclass, PanSuperclass


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


run_quiet = terminal.Run(verbose=False)


class ModelSEEDCompound:
    """Representation of a chemical in the network, with properties given by the ModelSEED Biochemistry database."""
    def __init__(self) -> None:
        self.modelseed_id: str = None
        self.modelseed_name: str = None
        self.kegg_aliases: Tuple[str] = None
        self.charge: int = None
        self.formula: str = None

class ModelSEEDReaction:
    """Representation of a reaction in the network, with properties given by the ModelSEED Biochemistry database."""
    def __init__(self) -> None:
        self.modelseed_id: str = None
        self.modelseed_name: str = None
        self.kegg_aliases: Tuple[str] = None
        self.ec_number_aliases: Tuple[str] = None
        # Compounds, coefficients, and compartments have corresponding elements.
        self.compounds: Tuple[ModelSEEDCompound] = None
        self.coefficients: Tuple[int] = None
        self.compartments: Tuple[str] = None
        self.reversibility: bool = None

class KO:
    """Representation of a KEGG Ortholog in the network."""
    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        # Map *ModelSEED reaction ID* to *ModelSEED reaction object or reaction aliases* in the
        # following 3 dictionaries.
        self.reactions: Dict[str, ModelSEEDReaction] = {}
        # Record the KEGG REACTION IDs *encoded by the KO* that are aliases of the ModelSEED
        # reaction ID. These could be a subset of the KEGG reaction aliases of the ModelSEED
        # reaction. The same is true of EC numbers.
        self.kegg_reaction_aliases: Dict[str, List[str]] = {}
        self.ec_number_aliases: Dict[str, List[str]] = {}

class Gene:
    """Representation of a gene in the metabolic network."""
    def __init__(self) -> None:
        self.gcid: int = None
        # KOs matching the gene
        self.kos: List[KO] = []
        # Record the strength of each KO match.
        self.e_values: List[float] = []

class GeneCluster:
    """Representation of a gene cluster in the metabolic network."""
    def __init__(self) -> None:
        self.gene_cluster_id: int = None
        # Consensus KO among the genes in the cluster. The KO annotations of genes in genomes that
        # underlie each consensus annotation are not tracked. This would require storing more data
        # on the consensus annotations in `dbops.PanSuperclass.get_gene_cluster_function_summary`.
        self.ko: KO = None

class Bin:
    """Representation of a bin of genes or gene clusters."""
    pass

class GeneBin(Bin):
    """Representation of a bin of genes."""
    def __init__(self) -> None:
        self.genes: List[Gene] = []

class GeneClusterBin(Bin):
    """Representation of a bin of gene clusters."""
    def __init__(self) -> None:
        self.gene_clusters: List[GeneCluster] = []

class BinCollection:
    """Representation of a collection of bins."""
    def __init__(self) -> None:
        self.bins: List[Bin] = []

class ReactionNetwork:
    """A reaction network predicted from KEGG KO and ModelSEED annotations."""
    def __init__(self) -> None:
        # Map KO ID to KO object.
        self.kos: Dict[str, KO] = {}
        # Map ModelSEED reaction ID to reaction object.
        self.reactions: Dict[str, ModelSEEDReaction] = {}
        # Map ModelSEED compound ID to compound object.
        self.metabolites: Dict[str, ModelSEEDCompound] = {}
        # The following dictionaries map reaction aliases in the network: as in, not all known
        # aliases, but only those sourced from KOs and contributing ModelSEEDReaction objects.
        self.kegg_modelseed_aliases: Dict[str, List[str]] = {}
        self.ec_number_modelseed_aliases: Dict[str, List[str]] = {}
        self.modelseed_kegg_aliases: Dict[str, List[str]] = {}
        self.modelseed_ec_number_aliases: Dict[str, List[str]] = {}

class GenomicNetwork(ReactionNetwork):
    """
    A reaction network predicted from KEGG KO and ModelSEED annotations of genes.

    Here is an example of the information used to create a genomic network.
    gene 1 ---> KO 1 ---> KEGG reaction 1 ---> ModelSEED reaction 1 ---> ModelSEED metabolites 1, 2, ...
    |      |         |
    |      |         ---> EC number 1 -------> ModelSEED reaction 1 ---> ModelSEED metabolites 1, 2, ...
    |      |         |                |
    |      |         |                -------> ModelSEED reaction 2 ---> ...
    |      |         |
    |      |         ---> EC number 2 -------> ...
    |      |
    |      ---> KO 2 ---> ...
    |
    gene 2 ---> ...
    """
    def __init__(self) -> None:
        # map gene caller ID to gene object
        super().__init__()
        self.genes: Dict[int, Gene] = {}
        self.bins: Dict[str, GeneBin] = {}
        self.collection: BinCollection = None

    def remove_metabolites_without_formula(self, path: str = None) -> None:
        """
        Remove metabolites without a formula in the ModelSEED database from the network.

        Reactions that involve the metabolite are also removed from the network, as are KOs that are
        predicted to only catalyze these reactions and genes that are only associated with these KOs.

        Additional metabolites that are exclusive to removed reactions can also be removed from the
        network. In the output table of removed metabolites, these additional metabolites have a formula.

        path : str, None
            If not None, write 3 tab-delimited files of metabolites, reactions, and genes removed
            from the network to file locations based on the path. For example, if the argument,
            'removed.tsv', is provided, then the following files will be written:
            'removed-metabolites.tsv', 'removed-reactions.tsv', and 'removed-genes.tsv'.
        """
        if path:
            filesnpaths.is_output_file_writable(path)
            path_basename, path_extension = os.path.splitext(path)
            metabolite_path = f"{path_basename}-metabolites{path_extension}"
            reaction_path = f"{path_basename}-reactions{path_extension}"
            gene_path = f"{path_basename}-genes{path_extension}"
        else:
            metabolite_path = None
            reaction_path = None
            gene_path = None

        metabolites_to_remove = []
        for modelseed_compound_id, metabolite in self.metabolites.items():
            # ModelSEED compounds without a formula have a formula value of None in the network object.
            if metabolite.formula is None:
                metabolites_to_remove.append(modelseed_compound_id)
        removed = self.purge_metabolites(metabolites_to_remove)

        if not path:
            return

        metabolite_table = []
        for metabolite in removed['metabolite']:
            metabolite: ModelSEEDCompound
            metabolite_table.append((metabolite.modelseed_id, metabolite.modelseed_name, metabolite.formula))

        reaction_table = []
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            row = []
            row.append(reaction.modelseed_id)
            row.append(reaction.modelseed_name)
            removed_metabolites: List[str] = []
            for metabolite in reaction.compounds:
                if metabolite.modelseed_id in metabolites_to_remove:
                    removed_metabolites.append(metabolite.modelseed_id)
            row.append(', '.join(removed_metabolites))
            reaction_table.append(row)

        gene_table = []
        removed_reactions: Set[str] = set([row[0] for row in reaction_table])
        for gene in removed['gene']:
            gene: Gene
            row = []
            row.append(gene.gcid)
            gene_reactions = []
            for ko in gene.kos:
                for modelseed_reaction_id in ko.reactions:
                    gene_reactions.append(modelseed_reaction_id)
            gene_reactions = set(gene_reactions)
            row.append(', '.join(gene_reactions.intersection(removed_reactions)))
            row.append(', '.join(ko.id for ko in gene.kos))
            row.append('!!!'.join(ko.name for ko in gene.kos))
            gene_table.append(row)

        pd.DataFrame(
            metabolite_table,
            columns=['modelseed_id', 'modelseed_name', 'formula']
        ).to_csv(metabolite_path, sep='\t', index=False)
        pd.DataFrame(
            reaction_table,
            columns=['modelseed_id', 'modelseed_name', 'removed_metabolite_modelseed_ids']
        ).to_csv(reaction_path, sep='\t', index=False)
        pd.DataFrame(
            gene_table,
            columns=['gene_callers_id', 'removed_reaction_modelseed_ids', 'ko_ids', 'ko_names']
        ).to_csv(gene_path, sep='\t', index=False)

    def purge_metabolites(self, metabolites_to_remove: List[str]) -> Dict[str, List]:
        """
        Remove any trace of the given metabolites from the network.

        Reactions involving the metabolite are also purged from the network. KOs that were only
        associated with removed reactions are purged; genes that were only associated with removed
        KOs are purged.

        Removal of reactions involving the metabolite can also result in other metabolites being
        being removed from the network, those that exclusively participate in these reactions.

        Parameters
        ==========
        metabolites_to_remove : list
            list of ModelSEED compound IDs identifying metabolites to remove

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, 'purge_reactions', then the dictionary
            will look like the following:
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'ko': [<removed KO objects>],
                'gene': [<removed Gene objects>]
            }

            If this method is called from the method, 'purge_reactions', then the dictionary will
            only contain one entry:
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }
        """
        removed_metabolites = []
        for modelseed_compound_id in metabolites_to_remove:
            try:
                removed_metabolites.append(self.metabolites.pop(modelseed_compound_id))
            except KeyError:
                # This occurs when removing other "unintended" metabolites from the network.
                # 'purge_metabolites' was first called with metabolites of interest, then
                # 'purge_reactions' was called from within the method to remove reactions involving
                # the metabolites of interest, and then 'purge_metabolites' was called AGAIN from
                # within 'purge_reactions' to remove other metabolites that were exclusively found
                # in the removed reactions. In this last call of 'purge_metabolites', the
                # 'metabolites_to_remove' also include the metabolites of interest that were already
                # removed from 'self.metabolites' in the original 'purge_metabolites' call. This
                # KeyError occurs when trying to remove those already-removed metabolites. A bit complicated!

                # Alternatively, this can occur if the metabolite from 'metabolites_to_remove' is
                # not in the network.
                pass

        if not removed_metabolites:
            return {'metabolite': [], 'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'ko': [], 'gene': []}

        reactions_to_remove = []
        for modelseed_reaction_id, reaction in self.reactions.items():
            for compound in reaction.compounds:
                if compound.modelseed_id in metabolites_to_remove:
                    reactions_to_remove.append(modelseed_reaction_id)
                    break

        removed = {'metabolite': removed_metabolites}
        if reactions_to_remove:
            removed_cascading_up = self.purge_reactions(reactions_to_remove)
            # There may be other metabolites exclusively involved in the removed reactions; these
            # metabolites were therefore also removed.
            removed['metabolite'] = removed_metabolites + removed_cascading_up.pop('metabolite')
        else:
            # This method must have been called from the method, 'purge_reactions', because the
            # reactions containing the metabolites were already removed from the network.
            removed_cascading_up = {'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'ko': [], 'gene': []}
        removed.update(removed_cascading_up)
        return removed

    def purge_reactions(self, reactions_to_remove: List[str]) -> Dict[str, List]:
        """
        Remove any trace of the given reactions from the network.

        Metabolites that exclusively participate in removed reactions are purged. KOs that were only
        associated with removed reactions are purged; genes that were only associated with removed
        KOs are purged.

        Parameters
        ==========
        reactions_to_remove : list
            list of ModelSEED reaction IDs identifying reactions to remove

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, 'purge_metabolites', or the method,
            'purge_kos', then the dictionary will look like the following:
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>],
                'ko': [<removed KO objects>],
                'gene': [<removed Gene objects>]
            }

            If this method is called from the method, 'purge_metabolites', then the dictionary will
            look like the following:
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [],
                'ko': [<removed KO objects>],
                'gene': [<removed Gene objects>]
            }

            If this method is called from the method, 'purge_kos', then the dictionary will look
            like the following:
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects],
                'ko': [],
                'gene': []
            }
        """
        removed_reactions: List[ModelSEEDReaction] = []
        for modelseed_reaction_id in reactions_to_remove:
            try:
                removed_reactions.append(self.reactions.pop(modelseed_reaction_id))
            except KeyError:
                # This occurs when the original method called is 'purge_reactions', followed by
                # 'purge_kos', followed by this method again -- 'removed_reactions' will be empty.
                # Alternatively, this occurs if the reaction in 'reactions_to_remove' is not in the network.
                continue
            try:
                self.modelseed_kegg_aliases.pop(modelseed_reaction_id)
            except KeyError:
                # The reaction has no KO KEGG REACTION aliases.
                pass
            try:
                self.modelseed_ec_number_aliases.pop(modelseed_reaction_id)
            except KeyError:
                # The reaction has no KO EC number aliases.
                pass

        if not removed_reactions:
            return {'metabolite': [], 'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'ko': [], 'gene': []}

        kegg_reactions_to_remove = []
        for kegg_reaction_id, modelseed_reaction_ids in self.kegg_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the KEGG reaction were removed, so remove the KEGG reaction as well.
                kegg_reactions_to_remove.append(kegg_reaction_id)
                continue
            for idx in sorted(aliases_to_remove, reverse=True):
                # Only some of the ModelSEED reactions aliased by the KO KEGG reaction were removed.
                modelseed_reaction_ids.pop(idx)
        removed_kegg_reactions = []
        for kegg_reaction_id in kegg_reactions_to_remove:
            removed_kegg_reactions.append(self.kegg_modelseed_aliases.pop(kegg_reaction_id))

        ec_numbers_to_remove = []
        for ec_number, modelseed_reaction_ids in self.ec_number_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the EC number were removed, so remove the EC number as well.
                ec_numbers_to_remove.append(ec_number)
                continue
            for idx in sorted(aliases_to_remove, reverse=True):
                # Only some of the ModelSEED reactions aliased by the EC number were removed.
                modelseed_reaction_ids.pop(idx)
        removed_ec_numbers = []
        for ec_number in ec_numbers_to_remove:
            removed_ec_numbers.append(self.ec_number_modelseed_aliases.pop(ec_number))

        metabolites_to_remove: List[str] = []
        for reaction in removed_reactions:
            for metabolite in reaction.compounds:
                metabolites_to_remove.append(metabolite.modelseed_id)
        metabolites_to_remove = list(set(metabolites_to_remove))
        for reaction in self.reactions.values():
            metabolites_to_spare: List[int] = []
            for metabolite in reaction.compounds:
                for idx, modelseed_compound_id in enumerate(metabolites_to_remove):
                    if modelseed_compound_id == metabolite.modelseed_id:
                        # Do not remove the metabolite, because it participates in a retained reaction.
                        metabolites_to_spare.append(idx)
            for idx in sorted(metabolites_to_spare, reverse=True):
                metabolites_to_remove.pop(idx)
        if metabolites_to_remove:
            removed_cascading_down = self.purge_metabolites(metabolites_to_remove)
            removed_cascading_down.pop('reaction')
            removed_cascading_down.pop('kegg_reaction')
            removed_cascading_down.pop('ec_number')
        else:
            # No metabolites were exclusive to the removed reactions. This cannot happen if this
            # method was called from within 'purge_metabolites'.
            removed_cascading_down = {'metabolite': []}

        kos_to_remove = []
        for ko_id, ko in self.kos.items():
            ko_reactions_to_remove = []
            for modelseed_reaction_id in ko.reactions:
                if modelseed_reaction_id in reactions_to_remove:
                    ko_reactions_to_remove.append(modelseed_reaction_id)
            if len(ko_reactions_to_remove) == len(ko.reactions):
                # All reactions associated with the KO were removed, so remove the KO as well.
                kos_to_remove.append(ko_id)
                continue
            for modelseed_reaction_id in ko_reactions_to_remove:
                # Only some of the reactions associated with the KO are invalid.
                ko.reactions.pop(modelseed_reaction_id)
                try:
                    ko.kegg_reaction_aliases.pop(modelseed_reaction_id)
                except KeyError:
                    # The reaction has no KO KEGG REACTION aliases.
                    pass
                try:
                    ko.ec_number_aliases.pop(modelseed_reaction_id)
                except KeyError:
                    # The reaction has no KO EC number aliases.
                    pass
        # If this method was called from 'purge_kos' then the KOs that are only associated with
        # reactions removed here were already removed from the network, and 'kos_to_remove' would be
        # empty. In contrast, if this method was not called from 'purge_kos', but zero KOs are
        # exclusively associated with reactions removed here, then 'kos_to_remove' would likewise be empty.
        if kos_to_remove:
            removed_cascading_up = self.purge_kos(kos_to_remove)
            removed_cascading_up.pop('reaction')
            removed_cascading_up.pop('kegg_reaction')
            removed_cascading_up.pop('ec_number')
        else:
            removed_cascading_up = {'ko': [], 'gene': []}

        removed = {'reaction': removed_reactions, 'kegg_reaction': removed_kegg_reactions, 'ec_number': removed_ec_numbers}
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def purge_kos(self, kos_to_remove: List[str]) -> Dict[str, List]:
        """
        Remove any trace of the given KOs from the network.

        Reactions and metabolites that were only associated with removed KOs are purged. Genes that
        were only associated with removed KOs are purged.

        Parameters
        ==========
        kos_to_remove : list
            list of KO IDs identifying KOs to remove

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, 'purge_reactions', or the method,
            'purge_genes', then the dictionary will look like the following:
            {
                'ko': [<removed KO objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>],
                'gene': [<removed Gene objects>]
            }

            If this method is called from the method, 'purge_reactions', then the dictionary will
            look like the following:
            {
                'ko': [<removed KO objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': [],
                'gene': [<removed Gene objects>]
            }

            If this method is called from the method, 'purge_genes', then the dictionary will
            look like the following:
            {
                'ko': [<removed KO objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects],
                'genes': []
            }
        """
        removed_kos: List[KO] = []
        for ko_id in kos_to_remove:
            try:
                removed_kos.append(self.kos.pop(ko_id))
            except KeyError:
                # This occurs when the original method called is 'purge_kos', followed by
                # 'purge_genes', followed by this method again -- 'removed_kos' will be empty.
                # Alternatively, this occurs if the KO in 'kos_to_remove' is not in the network.
                pass

        if not removed_kos:
            return {'metabolite': [], 'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'ko': [], 'gene': []}

        reactions_to_remove: List[str] = []
        for ko in removed_kos:
            for modelseed_reaction_id in ko.reactions:
                reactions_to_remove.append(modelseed_reaction_id)
        reactions_to_remove = list(set(reactions_to_remove))
        for ko in self.kos.values():
            reactions_to_spare: List[int] = []
            for modelseed_reaction_id in ko.reactions:
                for idx, modelseed_reaction_id_to_remove in enumerate(reactions_to_remove):
                    if modelseed_reaction_id == modelseed_reaction_id_to_remove:
                        # The reaction is associated with a retained KO, so do not remove the reaction.
                        reactions_to_spare.append(idx)
            for idx in sorted(reactions_to_spare, reverse=True):
                reactions_to_remove.pop(idx)
        if reactions_to_remove:
            removed_cascading_down = self.purge_reactions(reactions_to_remove)
            removed_cascading_down.pop('ko')
        else:
            # This method must have been called from the method, 'purge_reactions', because the
            # reactions that are only associated with the removed KOs were already removed from the network.
            removed_cascading_down = {'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'metabolite': []}

        genes_to_remove: List[str] = []
        for gcid, gene in self.genes.items():
            gene_kos_to_remove: List[int] = []
            for idx, ko in enumerate(gene.kos):
                if ko.id in kos_to_remove:
                    gene_kos_to_remove.append(idx)
            if len(gene_kos_to_remove) == len(gene.kos):
                # All KOs matching the gene were removed, so remove it.
                genes_to_remove.append(gcid)
                continue
            for idx in sorted(gene_kos_to_remove, reverse=True):
                gene.kos.pop(idx)
                gene.e_values.pop(idx)
        # If this method was called from 'purge_genes' then the genes that are only associated with
        # KOs removed here were already removed from the network, and 'genes_to_remove' would be
        # empty. In contrast, if this method was not called from 'purge_genes', but zero genes are
        # only associated with KOs removed here, then 'genes_to_remove' would likewise be empty.
        if genes_to_remove:
            removed_cascading_up = self.purge_genes(genes_to_remove)
            removed_cascading_up.pop('ko')
        else:
            removed_cascading_up = {'gene': []}

        removed = {'ko': removed_kos}
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def purge_genes(self, genes_to_remove: List[str]) -> Dict[str, List]:
        """
        Remove any trace of the given genes from the network.

        KOs, reactions, and metabolites that were only associated with removed genes are purged.

        Parameters
        ==========
        genes_to_remove : list
            list of gene callers IDs identifying genes to remove

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, 'purge_kos', then the dictionary will
            look like the following:
            {
                'gene': [<removed Gene objects>],
                'ko': [<removed KO objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>]
            }

            If this method is called from the method, 'purge_kos', then the dictionary will look
            like the following:
            {
                'gene': [<removed Gene objects>],
                'ko': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }
        """
        removed_genes: List[Gene] = []
        for gcid in genes_to_remove:
            try:
                removed_genes.append(self.genes.pop(gcid))
            except KeyError:
                # This occurs if the gene in 'genes_to_remove' is not in the network.
                pass

        if not removed_genes:
            return {'metabolite': [], 'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'ko': [], 'gene': []}

        kos_to_remove: List[str] = []
        for gene in removed_genes:
            for ko in gene.kos:
                kos_to_remove.append(ko.id)
        kos_to_remove = list(set(kos_to_remove))
        for gene in self.genes.values():
            kos_to_spare: List[int] = []
            for ko in gene.kos:
                for idx, ko_id in enumerate(kos_to_remove):
                    if ko.id == ko_id:
                        # The KO is associated with a retained gene, so do not remove the KO.
                        kos_to_spare.append(idx)
            for idx in sorted(kos_to_spare, reverse=True):
                kos_to_remove.pop(idx)
        if kos_to_remove:
            removed_cascading_down = self.purge_kos(kos_to_remove)
            removed_cascading_down.pop('gene')
        else:
            # This method must have been called from the method, 'purge_kos', because the
            # KOs that are only associated with the removed genes were already removed from the network.
            removed_cascading_down = {'ko': [], 'reaction': [], 'kegg_reaction': [], 'ec_number': [], 'metabolite': []}

        # TODO: remove genes from self.bins

        removed = {'gene': removed_genes}
        removed.update(removed_cascading_down)
        return removed

    def export_json(
        self,
        path: str,
        overwrite: bool = False,
        objective: str = None,
        remove_missing_objective_metabolites: bool = False,
        # record_bins: tuple = ('gene', ),
        indent: int = 2,
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Export the network to a metabolic model file in JSON format.

        All information from the network is included in the JSON so that the file can by imported by
        anvi'o as a GenomicNetwork object containing the same information.

        Parameters
        ==========
        path : str
            output JSON file path

        overwrite : bool, False
            Overwrite the JSON file if it already exists.

        objective : str, None
            An objective to use in the model, stored as the first entry in the JSON 'reactions'
            array. Currently, the only valid options are None and 'e_coli_core'.

            None does not add an objective to the JSON, meaning that FBA cannot be performed on the model.

            'e_coli_core' is the biomass objective from the COBRApy example JSON file of E. coli
            "core" metabolism, 'e_coli_core.json'.

        remove_missing_objective_metabolites : bool, False
            If True, remove metabolites from the JSON objective that are not produced or consumed in
            the reaction network. FBA fails with metabolites outside the network.

        record_bins : tuple, ('gene', )
            Record bin membership in JSON entries, if a collection of bins is present in the
            reaction network. By default, bin membership is only recorded for genes with the
            argument, ('gene', ). 'reaction' and 'metabolite' can also be provided in a tuple
            argument (e.g., ('reaction', ) or ('metabolite', 'reaction', 'gene')) to likewise record
            in which bins the reaction and metabolite entries occur. Pass an argument of None to not
            record bin membership in the JSON at all.

        indent : int, 2
            spaces of indentation per nesting level in JSON file

        progress : terminal.Progress, terminal.Progress()
        """
        progress.new("Constructing JSON")
        progress.update("Setting up")
        filesnpaths.is_output_file_writable(path, ok_if_exists=overwrite)
        json_dict = JSONStructure.get()
        json_genes: List[Dict] = json_dict['genes']
        json_reactions: List[Dict] = json_dict['reactions']
        json_metabolites: List[Dict] = json_dict['metabolites']
        if objective == 'e_coli_core':
            objective_dict = JSONStructure.get_e_coli_core_objective()
            if remove_missing_objective_metabolites:
                objective_metabolites: Dict = objective_dict['metabolites']
                objective_original_metabolites: Dict = objective_dict['notes']['original_metabolite_ids']
                missing_metabolite_ids = []
                # The original objective had metabolite BiGG IDs, which were replaced with KEGG COMPOUND IDs.
                missing_original_metabolite_ids = []
                for metabolite_id, original_metabolite_id in zip(objective_metabolites, objective_original_metabolites):
                    if metabolite_id[:-2] not in self.metabolites:
                        # The metabolite (removing localization substring) is not in the network.
                        missing_metabolite_ids.append(metabolite_id)
                        missing_original_metabolite_ids.append(original_metabolite_id)
                for metabolite_id in missing_metabolite_ids:
                    objective_metabolites.pop(metabolite_id)
                for original_metabolite_id in missing_original_metabolite_ids:
                    objective_original_metabolites.pop(original_metabolite_id)
            json_reactions.append(objective_dict)
        elif objective != None:
            raise ConfigError(f"Anvi'o does not recognize an objective with the name, '{objective}'.")

        progress.update("Genes")
        reaction_genes: Dict[str, List[str]] = {}
        reaction_kos: Dict[str, List[KO]] = {}
        for gcid, gene in self.genes.items():
            gene_entry = JSONStructure.get_gene_entry()
            gcid_str = str(gcid)
            gene_entry['id'] = gcid_str
            notes = gene_entry['notes']
            # Record KO IDs and annotation e-values in the notes section of the gene entry.
            notes['ko'] = notes_kos = {}
            for ko, e_value in zip(gene.kos, gene.e_values):
                notes_kos[ko.id] = str(e_value)
                for modelseed_reaction_id in ko.reactions:
                    try:
                        reaction_genes[modelseed_reaction_id].append(gcid_str)
                    except KeyError:
                        reaction_genes[modelseed_reaction_id] = [gcid_str]
                    try:
                        reaction_kos[modelseed_reaction_id].append(ko)
                    except KeyError:
                        reaction_kos[modelseed_reaction_id] = [ko]
            json_genes.append(gene_entry)

        progress.update("Reactions")
        compound_compartments: Dict[str, Set[str]] = {}
        for modelseed_reaction_id, reaction in self.reactions.items():
            reaction_entry = JSONStructure.get_reaction_entry()
            reaction_entry['id'] = modelseed_reaction_id
            reaction_entry['name'] = reaction.modelseed_name
            metabolites = reaction_entry['metabolites']
            for compound, compartment, coefficient in zip(reaction.compounds, reaction.compartments, reaction.coefficients):
                modelseed_compound_id = compound.modelseed_id
                metabolites[f"{modelseed_compound_id}_{compartment}"] = coefficient
                try:
                    compound_compartments[modelseed_compound_id].add(compartment)
                except KeyError:
                    compound_compartments[modelseed_compound_id] = set(compartment)
            if not reaction.reversibility:
                # By default, the reaction entry was set up to be reversible; here make it irreversible.
                reaction_entry['lower_bound'] = 0.0
            reaction_entry['gene_reaction_rule'] = " or ".join([gcid for gcid in reaction_genes[modelseed_reaction_id]])
            notes = reaction_entry['notes']
            # Record gene KO annotations which aliased the reaction via KEGG REACTION or EC number.
            notes['ko'] = ko_notes = {}
            ko_kegg_aliases = []
            ko_ec_number_aliases = []
            for ko in reaction_kos[modelseed_reaction_id]:
                try:
                    kegg_aliases = ko.kegg_reaction_aliases[modelseed_reaction_id]
                except KeyError:
                    kegg_aliases = []
                try:
                    ec_number_aliases = ko.ec_number_aliases[modelseed_reaction_id]
                except KeyError:
                    ec_number_aliases = []
                ko_notes[ko.id] = {'kegg.reaction': kegg_aliases, 'ec-code': ec_number_aliases}
                ko_kegg_aliases += kegg_aliases
                ko_ec_number_aliases += ec_number_aliases
            ko_kegg_aliases = set(ko_kegg_aliases)
            ko_ec_number_aliases = set(ko_ec_number_aliases)
            # Record other KEGG REACTION or EC number aliases of the reaction in the ModelSEED
            # database that did not happen to be associated with KO annotations.
            notes['other_aliases'] = {
                'kegg.reaction': list(set(reaction.kegg_aliases).difference(ko_kegg_aliases)),
                'ec-code': list(set(reaction.ec_number_aliases).difference(ko_ec_number_aliases))
            }
            json_reactions.append(reaction_entry)

        progress.update("Metabolites")
        for modelseed_compound_id, metabolite in self.metabolites.items():
            modelseed_compound_name = metabolite.modelseed_name
            charge = metabolite.charge
            formula = metabolite.formula
            kegg_compound_aliases = list(metabolite.kegg_aliases)
            for compartment in compound_compartments[modelseed_compound_id]:
                metabolite_entry = JSONStructure.get_metabolite_entry()
                metabolite_entry['id'] = f"{modelseed_compound_id}_{compartment}"
                metabolite_entry['name'] = modelseed_compound_name
                metabolite_entry['compartment'] = compartment
                # Compounds without a formula have a nominal charge of 10000000 in the ModelSEED
                # compounds database, which is replaced by None in the reaction network and 0 in the JSON.
                metabolite_entry['charge'] = charge if charge is not None else 0
                metabolite_entry['formula'] = formula if formula is not None else ""
                metabolite_entry['annotation']['kegg.compound'] = kegg_compound_aliases
                json_metabolites.append(metabolite_entry)

        progress.update("Saving")
        with open(path, 'w') as f:
            json.dump(json_dict, f, indent=indent)
        progress.end()

# class PangenomicNetwork(ReactionNetwork):
#     """A reaction network predicted from KEGG KO and ModelSEED annotations of pangenomic gene clusters."""
#     def __init__(self) -> None:
#         # map gene cluster ID to gene cluster object
#         self.gene_clusters: Dict[str, GeneCluster] = {}
#         self.bins: Dict[str, GeneClusterBin] = {}
#         self.collection: BinCollection = None
#         super().__init__()

#     def export_json(
#         self,
#         path: str,
#         annotate_genes: tuple = ('genome', 'bin', ),
#         annotate_reactions: tuple = ('genome', 'bin', 'kegg_reaction', 'ec_number'),
#         annotate_metabolites: tuple = ('genome', 'bin', 'kegg_compound'),
#         run: terminal.Run = terminal.Run(),
#         progress: terminal.Progress = terminal.Progress()
#     ) -> None:
#         """
#         Export the network to a metabolic model file in JSON format. *Gene entries in this file
#         represent gene clusters.* Optionally, gene, reaction, and metabolite entries in this file
#         are annotated with the names of genomes and names of gene cluster bins in which they occur.

#         Parameters
#         ==========
#         path : str
#             output JSON file path

#         annotate_genes : tuple, ('genome', 'bin', )
#             Annotate gene (cluster) entries in the JSON file with additional data, selecting
#             from the following:

#             'genome' : genomes in which the genes of the cluster occur

#             'bin' : bins in which the gene cluster occurs

#             'all_ko' : all KOs associated with genes in the cluster, sorted in descending order of
#                 the number of genes in the cluster that were associated with each KO and then mean
#                 e-value of gene-KO assignments

#             'ko' : KOs associated with the gene cluster that yielded reactions in the network,
#                 sorted in descending order of the number of genes in the cluster that were
#                 associated with each KO and then mean e-value of gene-KO assignments

#             'ko_count' : number of genes in the cluster that were associated with each KO; if
#                 'all_ko' is provided, then each value corresponds to a KO in 'all_ko', whereas if
#                 only 'ko' is provided, then each value corresponds to a KO in 'ko'

#             'e_value' : mean scores of KO associations with genes in the cluster; if 'all_ko' is
#                 provided, then each value corresponds to a KO in 'all_ko', whereas if only 'ko' is
#                 provided, then each value corresponds to a KO in 'ko'

#         annotate_reactions : tuple, ('genome', 'bin', 'kegg_reaction', 'ec_number')
#             Annotate reaction entries in the JSON file with additional data, selecting from the following:

#             'genome' : genomes in which the reaction occurs

#             'bin' : bins in which the reaction occurs

#             'kegg_reaction' : KO-associated KEGG reaction IDs yielding the ModelSEED reaction

#             'ec_number' : KO-associated EC numbers yielding the ModelSEED reaction

#             'ko' : KOs yielding the ModelSEED reaction

#         annotate_metabolites : tuple, ('genome', 'bin', 'kegg_compound')
#             Annotate metabolite entries in the JSON file with additional data, selecting from the following:

#             'genome' : genomes in which the metabolite occurs

#             'bin' : bins in which the metabolite occurs

#             'kegg_compound' : KEGG compound aliases of the ModelSEED compound

#             'ko' : KOs yielding the ModelSEED compound

#         run : terminal.Run, terminal.Run()

#         progress : terminal.Progress, terminal.Progress()
#         """
#         pass

class JSONStructure:
    """JSON structure of metabolic model file."""
    def get() -> Dict:
        """Top-level file framework."""
        return {
            'metabolites': [],
            'reactions': [],
            'genes': [],
            'id': '',
            'compartments': {
                'c': 'cytosol',
                'e': 'extracellular space'
            },
            'version': '1'
        }

    def get_metabolite_entry() -> Dict:
        """"Format of each object in the 'metabolites' array."""
        return {
            'id': '',
            'name': '',
            'compartment': '',
            'charge': 0, # placeholder: uncharged
            'formula': '',
            'notes': {},
            'annotation': {}
        }

    def get_reaction_entry() -> Dict:
        """Format of each object in the 'reactions' array."""
        return {
            'id': '',
            'name': '',
            'metabolites': {},
            # By default, make the reaction perfectly reversible.
            'lower_bound': -1000.0,
            'upper_bound': 1000.0,
            'gene_reaction_rule': '',
            'subsystem': '',
            'notes': {},
            'annotation': {}
        }

    def get_gene_entry() -> Dict:
        """Format of each object in the 'genes' array."""
        return {
            'id': '',
            'name': '',
            'notes': {},
            'annotation': {}
        }

    def get_e_coli_core_objective() -> Dict:
        """
        Biomass objective from the 'reactions' array in the COBRApy example JSON file,
        'e_coli_core.json', with KBase/ModelSEED compound IDs replacing BiGG metabolite IDs.
        """
        return {
            'id': 'BIOMASS_Ecoli_core_w_GAM',
            'name': 'Biomass Objective Function with GAM',
            'metabolites': {
                'cpd00169_c': -1.496,
                'cpd00022_c': -3.7478,
                'cpd00008_c': 59.81,
                'cpd00024_c': 4.1182,
                'cpd00002_c': -59.81,
                'cpd00010_c': 3.7478,
                'cpd00236_c': -0.361,
                'cpd00072_c': -0.0709,
                'cpd00102_c': -0.129,
                'cpd00079_c': -0.205,
                'cpd00053_c': -0.2557,
                'cpd00023_c': -4.9414,
                'cpd00001_c': -59.81,
                'cpd00067_c': 59.81,
                'cpd00003_c': -3.547,
                'cpd00004_c': 3.547,
                'cpd00006_c': 13.0279,
                'cpd00005_c': -13.0279,
                'cpd00032_c': -1.7867,
                'cpd00061_c': -0.5191,
                'cpd00009_c': 59.81,
                'cpd00020_c': -2.8328,
                'cpd00101_c': -0.8977
            },
            'lower_bound': 0.0,
            'upper_bound': 1000.0,
            'gene_reaction_rule': '',
            'objective_coefficient': 1.0,
            'subsystem': 'Biomass and maintenance functions',
            'notes': {
                'original_bigg_ids': [
                    'Biomass_Ecoli_core_w_GAM'
                ],
                'original_metabolite_ids': {
                    '3pg_c': -1.496,
                    'accoa_c': -3.7478,
                    'adp_c': 59.81,
                    'akg_c': 4.1182,
                    'atp_c': -59.81,
                    'coa_c': 3.7478,
                    'e4p_c': -0.361,
                    'f6p_c': -0.0709,
                    'g3p_c': -0.129,
                    'g6p_c': -0.205,
                    'gln__L_c': -0.2557,
                    'glu__L_c': -4.9414,
                    'h2o_c': -59.81,
                    'h_c': 59.81,
                    'nad_c': -3.547,
                    'nadh_c': 3.547,
                    'nadp_c': 13.0279,
                    'nadph_c': -13.0279,
                    'oaa_c': -1.7867,
                    'pep_c': -0.5191,
                    'pi_c': 59.81,
                    'pyr_c': -2.8328,
                    'r5p_c': -0.8977
                }
            },
            'annotation': {
                'bigg.reaction': [
                    'BIOMASS_Ecoli_core_w_GAM'
                ],
                'metanetx.reaction': [
                    'MNXR96280'
                ],
                'sbo': 'SBO:0000629'
            }
        }

class KODatabase:
    """
    Representation of the KEGG KO database used in the construction of reaction networks.

    Unless an alternative directory is provided, the database is downloaded and set up in a
    default anvi'o data directory, and loaded from this directory in network construction.
    """
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/MISC/KEGG/KO_REACTION_NETWORK')
    expected_files = ['ko_info.txt', 'ko_data.tsv']

    def __init__(self, ko_dir: str = None) -> None:
        """
        Load the table derived from downloaded KEGG KO entry files that relates KOs to KEGG
        reactions and EC numbers.

        Parameters
        ==========
        ko_dir : str, None
            Directory containing KO data table to use rather than the default.
        """
        if ko_dir:
            if not os.path.isdir(ko_dir):
                raise ConfigError(f"There is no such directory, '{ko_dir}'.")
        else:
            ko_dir = self.default_dir

        for expected_file in self.expected_files:
            if not os.path.isfile(os.path.join(ko_dir, expected_file)):
                raise ConfigError(f"No required file named '{expected_file}' was found in the KO directory, '{ko_dir}'.")

        f = open(os.path.join(ko_dir, 'ko_info.txt'))
        f.readline()
        self.release = ' '.join(f.readline().strip().split()[1:])
        f.close()

        self.ko_table = pd.read_csv(os.path.join(ko_dir, 'ko_data.tsv'), sep='\t', header=0, index_col=0, low_memory=False)

    def set_up(
        num_threads: int = 1,
        dir: str = None,
        reset: bool = False,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Download KEGG KO entry files and parse these files to construct a tab-delimited file
        relating KOs to KEGG reactions and EC numbers.

        Parameters
        ==========
        num_threads : int, 1
            Number of threads to use in parallelizing the download of KO files.

        dir : str, None
            Directory in which to create a subdirectory called `KO_REACTION_NETWORK`,
            in which files are downloaded and set up. This argument overrides
            the default directory.

        reset : bool, False
            If True, remove any existing 'KO_REACTION_NETWORK' database directory and the files
            therein. If False, an exception is raised if there are files in this directory.

        run : anvio.terminal.Run, None

        progress : anvio.terminal.Progress, None
        """
        if dir:
            if os.path.isdir(dir):
                ko_dir = os.path.join(dir, 'KO_REACTION_NETWORK')
            else:
                raise ConfigError(f"There is no such directory, '{dir}'. You should create it "
                                   "first if you want to use it.")
        else:
            ko_dir = KODatabase.default_dir
            parent_dir = os.path.dirname(ko_dir)
            if not os.path.exists(parent_dir):
                os.mkdir(parent_dir)
        if os.path.exists(ko_dir):
            if reset:
                shutil.rmtree(ko_dir)
            else:
                raise ConfigError(
                    f"The KO database directory, '{ko_dir}', already exists. 'reset' can be used "
                    "to remove the database at this location and set it up again."
                )
        os.mkdir(ko_dir)

        if num_threads == 1:
            run.warning(
                "Only 1 thread will be used to download KO files. It is advisable to set a higher "
                "number of threads to download faster."
            )
        assert type(num_threads) is int and num_threads > 0

        # Download a file for each entry in a KEGG database.
        download_root = 'https://rest.kegg.jp/'
        while True:
            # Break out of this loop upon confirming that the KEGG release didn't change in the
            # middle of downloading KO files.
            progress.new(f"Downloading KEGG KO files")
            # Get the database version before download.
            progress.update("Database info")
            info_before_path = os.path.join(ko_dir, 'ko_info_before.txt')
            utils.download_file(f'{download_root}info/ko', info_before_path)
            f = open(info_before_path)
            f.readline()
            release_before = ' '.join(f.readline().strip().split()[1:])
            f.close()

            # Get a list of all KO IDs.
            progress.update("KO list")
            list_path = os.path.join(ko_dir, 'ko_list.txt')
            utils.download_file(f'{download_root}list/ko', list_path)
            ko_ids = []
            f = open(list_path)
            for line in f:
                line.split()[0]
                ko_ids.append(line.split('\t')[0])
            f.close()

            # Download KO entry files.
            manager = mp.Manager()
            input_queue = manager.Queue()
            output_queue = manager.Queue()
            for ko_id in ko_ids:
                input_queue.put((f'{download_root}get/{ko_id}', os.path.join(ko_dir, f'{ko_id}.txt')))
            workers: List[mp.Process] = []
            for _ in range(num_threads):
                worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
                workers.append(worker)
                worker.start()
            downloaded_count = 0
            undownloaded_count = 0
            total = len(ko_ids)
            undownloaded = []
            while downloaded_count + undownloaded_count < total:
                output = output_queue.get()
                if output is True:
                    downloaded_count += 1
                    progress.update(f"{downloaded_count} / {total} KO files")
                else:
                    undownloaded_count += 1
                    undownloaded.append(os.path.splitext(os.path.basename(output))[0])
            for worker in workers:
                worker.terminate()
            if undownloaded:
                raise ConfigError(
                    "Unfortunately, files for the following KOs failed to download despite multiple attempts, "
                    f"and so the database needs to be set up again: {', '.join(undownloaded)}"
                )

            # Get the database version after download.
            progress.update("Database info (again)")
            info_after_path = os.path.join(ko_dir, 'ko_info.txt')
            utils.download_file(f'{download_root}info/ko', info_after_path)
            f = open(info_after_path)
            f.readline()
            release_after = ' '.join(f.readline().strip().split()[1:])
            f.close()

            # Check that the database had the same version before and after download.
            progress.end()
            if release_before == release_after:
                # Retain one of the info files and delete the other.
                info_path = info_after_path
                os.remove(info_before_path)
                break
            else:
                run.warning(
                    "It's your lucky day! The version of KEGG appears to have changed from "
                    f"'{release_before}' to '{release_after}' while anvi'o was downloading files "
                    "from the KO database. Anvi'o will now attempt to redownload all of the files. "
                )
        run.info(f"Total number of KOs/entry files", total)
        run.info("KEGG KO database version", release_after)
        run.info("KEGG KO list", list_path)
        run.info("KEGG KO info", info_path)

        progress.new("Processing KEGG KO database")
        # Make a tab-delimited file relating KO IDs and names to KEGG reactions and EC numbers.
        kos_data = {}
        paths = glob.glob(os.path.join(ko_dir, 'K*.txt'))
        for num_processed, path in enumerate(paths):
            progress.update(f"{num_processed} / {total} KO files")
            # Parse the KO file.
            ko_data = {}
            section = None
            # Unfortunately, a non-unicode character can crop up.
            f = open(path, 'rb')
            for line in f.read().decode(errors='replace').split('\n'):
                if line[0] == ' ':
                    pass
                else:
                    section = line.split()[0]
                if section == 'NAME':
                    # The name value follows 'NAME' at the beginning of the line.
                    ko_data['name'] = line[4:].strip()
                    # EC numbers associated with the KO are recorded at the end of the name value.
                    ec_string = re.search('\[EC:.*\]', line)
                    if ec_string:
                        ko_data['ec_numbers'] = ec_string[0][4:-1]
                elif section == 'DBLINKS':
                    # There is a row for each linked database in this section. There can be a row
                    # for KEGG REACTION database entries. The first line of the section starts with
                    # 'DBLINKS' and is followed by a value for a linked database. Values from the
                    # linked database are separated by ': ' from the name of the database, e.g.,
                    # 'RN: R00001'.
                    split_line = line.split()
                    try:
                        rn_index = split_line.index('RN:')
                    except ValueError:
                        continue
                    ko_data['reactions'] = ' '.join(split_line[rn_index + 1:])
                elif section == 'GENES':
                    # This is the section after DBLINKS.
                    break
            f.close()
            ko_id = os.path.splitext(os.path.basename(path))[0]
            kos_data[ko_id] = ko_data
        progress.update("Making a table mapping KOs to KEGG reactions and EC numbers")
        columns = {h: [] for h in ['name', 'reactions', 'ec_numbers']}
        for ko_data in kos_data.values():
            for h, column in columns.items():
                try:
                    value = ko_data[h]
                except KeyError:
                    value = None
                column.append(value)
        table: pd.DataFrame = pd.DataFrame.from_dict(columns)
        table.index = kos_data
        table = table.sort_index()
        table_path = os.path.join(ko_dir, 'ko_data.tsv')
        table.to_csv(table_path, sep='\t')
        progress.end()
        run.info("Table of select KEGG KO data", table_path)

        # Tarball the KO entry files.
        progress.new("Compressing downloaded KEGG KO entry files")
        progress.update("...")
        ko_entries_dir = os.path.join(ko_dir, 'ko_entries')
        os.mkdir(ko_entries_dir)
        for path in paths:
            shutil.move(path, ko_entries_dir)
        tar_path = os.path.join(ko_dir, 'ko_entries.tar.gz')
        with tarfile.open(tar_path, mode='w:gz') as tar:
            tar.add(ko_entries_dir, arcname='.')
        progress.end()
        shutil.rmtree(ko_entries_dir)
        run.info("Archived KEGG KO entry files", tar_path)

def _download_worker(
    input_queue: mp.Queue,
    output_queue: mp.Queue,
    max_num_tries: int = 100,
    wait_secs: float = 10.0
) -> None:
    """Multiprocessing download worker."""
    while True:
        url, path = input_queue.get()
        num_tries = 0
        while True:
            try:
                utils.download_file(url, path)
                output = True
                break
            except (ConfigError, ConnectionResetError) as e:
                num_tries += 1
                if num_tries > max_num_tries:
                    output = path
                    break
                time.sleep(wait_secs)
        output_queue.put(output)

class ModelSEEDDatabase:
    """
    The ModelSEED Biochemistry database set up by anvi'o.

    By default, the database is loaded from a default directory of ModelSEED files unless an
    alternative directory is provided.
    """
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/MISC/ModelSEED')

    # Compounds are identified as cytosolic or extracellular in ModelSEED reactions.
    compartment_ids = {0: 'c', 1: 'e'}

    def __init__(self, modelseed_dir: str = None) -> None:
        """
        Load and set up reorganized tables of reactions and compounds from the ModelSEED directory.

        Parameters
        ==========
        modelseed_dir : str, None
            Directory of ModelSEED files to use instead of the default.
        """
        if modelseed_dir:
            if not os.path.isdir(modelseed_dir):
                raise ConfigError(f"There is no such directory, '{modelseed_dir}'.")
        else:
            modelseed_dir = self.default_dir
        sha_path = os.path.join(modelseed_dir, 'sha.txt')
        if not os.path.isfile(sha_path):
            raise ConfigError(
                f"No required file named 'sha.txt' was found in the ModelSEED directory, '{modelseed_dir}'."
            )
        reactions_path = os.path.join(modelseed_dir, 'reactions.tsv')
        if not os.path.isfile(reactions_path):
            raise ConfigError(
                f"No required file named 'reactions.tsv' was found in the ModelSEED directory, '{modelseed_dir}'."
            )
        compounds_path = os.path.join(modelseed_dir, 'compounds.tsv')
        if not os.path.isfile(compounds_path):
            raise ConfigError(
                f"No required file named 'compounds.tsv' was found in the ModelSEED directory, '{modelseed_dir}'."
            )

        with open(sha_path) as f:
            self.sha = f.read().strip()
        reactions_table = pd.read_csv(reactions_path, sep='\t', header=0, low_memory=False)
        self.compounds_table = pd.read_csv(compounds_path, sep='\t', header=0, index_col='id', low_memory=False)

        # Facilitate lookup of reaction data by KEGG REACTION ID via a reorganized reactions table.
        # Remove reactions without KEGG aliases.
        reactions_table_without_na = reactions_table.dropna(subset=['KEGG'])
        expanded = []
        ko_id_col = []
        for ko_ids, row in zip(
            reactions_table_without_na['KEGG'],
            reactions_table_without_na.itertuples(index=False)
        ):
            ko_ids: str
            # A ModelSEED reaction can have multiple KEGG aliases.
            for ko_id in ko_ids.split('; '):
                ko_id_col.append(ko_id)
                expanded.append(row)
        kegg_reactions_table = pd.DataFrame(expanded)
        kegg_reactions_table['KEGG_REACTION_ID'] = ko_id_col
        self.kegg_reactions_table = kegg_reactions_table

        # Facilitate lookup of reaction data by EC number via a reorganized reactions table.
        # Remove reactions without EC number aliases.
        reactions_table_without_na = reactions_table.dropna(subset=['ec_numbers'])
        expanded = []
        ec_number_col = []
        for ec_numbers, row in zip(
            reactions_table_without_na['ec_numbers'],
            reactions_table_without_na.itertuples(index=False)
        ):
            ec_numbers: str
            # A ModelSEED reaction can have multiple EC number aliases.
            for ec_number in ec_numbers.split('|'):
                ec_number_col.append(ec_number)
                expanded.append(row)
        ec_reactions_table = pd.DataFrame(expanded)
        ec_reactions_table['EC_number'] = ec_number_col
        self.ec_reactions_table = ec_reactions_table

    def set_up(
        dir: str = None,
        reset: bool = False,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Download the ModelSEED Biochemistry database, which consists of two tables of reaction and
        metabolite data, and reorganize the tables.

        Parameters
        ==========
        dir : str, None
            Directory in which to create a new subdirectory called 'MODELSEED', in which files are
            downloaded and set up. This argument overrides the default directory.

        reset : bool, False
            If True, remove any existing 'MODELSEED' database directory and the files therein. If
            False, an exception is raised if there are files in this directory.

        run : anvio.terminal.Run, None

        progress : anvio.terminal.Progress, None
        """
        if dir:
            if os.path.isdir(dir):
                modelseed_dir = os.path.join(dir, 'MODELSEED')
            else:
                raise ConfigError(f"There is no such directory, '{dir}'.")
        else:
            modelseed_dir = ModelSEEDDatabase.default_dir
            parent_dir = os.path.dirname(modelseed_dir)
            if not os.path.exists(parent_dir):
                os.mkdir(parent_dir)
        if os.path.exists(modelseed_dir):
            if reset:
                shutil.rmtree(modelseed_dir)
            else:
                raise ConfigError(
                    f"The ModelSEED database directory, '{modelseed_dir}', already exists. 'reset' "
                    "can be used to remove the database at this location and set it up again."
                )
        os.mkdir(modelseed_dir)

        def download(url, path):
            max_num_tries = 100
            wait_secs = 10.0
            num_tries = 0
            while True:
                try:
                    utils.download_file(url, path, progress=progress)
                    break
                except ConnectionResetError:
                    num_tries += 1
                    if num_tries > max_num_tries:
                        raise ConnectionResetError(
                            f"The connection was reset by the peer more than {max_num_tries} times, "
                            "the maximum number of attempts. Try setting up the ModelSEED database again."
                        )
                    time.sleep(wait_secs)
        # The commit SHA taken from the following file is stored in a text file to track the version
        # of the ModelSEED database.
        json_url = 'https://api.github.com/repos/ModelSEED/ModelSEEDDatabase/commits'
        json_path = os.path.join(modelseed_dir, 'commits.json')
        download(json_url, json_path)
        with open(json_path) as f:
            sha = json.load(f)[0]['sha']
        zip_url = f'https://github.com/ModelSEED/ModelSEEDDatabase/archive/{sha}.zip'
        zip_path = os.path.join(modelseed_dir, f'ModelSEEDDatabase-{sha}.zip')
        download(zip_url, zip_path)

        progress.new("Setting up ModelSEED files")
        progress.update("Extracting")
        with zipfile.ZipFile(zip_path, 'r') as f:
            f.extractall(modelseed_dir)
        reactions_path = os.path.join(modelseed_dir, f'ModelSEEDDatabase-{sha}', 'Biochemistry', 'reactions.tsv')
        compounds_path = os.path.join(modelseed_dir, f'ModelSEEDDatabase-{sha}', 'Biochemistry', 'compounds.tsv')
        shutil.move(reactions_path, modelseed_dir)
        shutil.move(compounds_path, modelseed_dir)
        reactions_path = os.path.join(modelseed_dir, 'reactions.tsv')
        compounds_path = os.path.join(modelseed_dir, 'compounds.tsv')
        sha_path = os.path.join(modelseed_dir, 'sha.txt')
        with open(sha_path, 'w') as f:
            f.write(sha)
        os.remove(json_path)
        os.remove(zip_path)
        shutil.rmtree(os.path.join(modelseed_dir, f'ModelSEEDDatabase-{sha}'))

        progress.update("Loading")
        reactions_table = pd.read_csv(reactions_path, sep='\t', header=0, low_memory=False)
        compounds_table = pd.read_csv(compounds_path, sep='\t', header=0, low_memory=False)

        progress.update("Reorganizing tables")
        # Reorganize the downloaded tables, storing in the same locations. The tables each have a
        # column of aliases, or IDs for the same reaction or compound from various databases. Split
        # these IDs into separate columns added to the end of the table, dropping the alias column.
        def expand_aliases(table: pd.DataFrame) -> pd.DataFrame:
            new_rows = []
            for aliases in table.aliases:
                aliases: str
                new_row = {}
                if pd.isna(aliases):
                    new_rows.append(new_row)
                    continue
                split_aliases = aliases.split('|')
                for alias in split_aliases:
                    sep_index = alias.index(': ')
                    alias_key = alias[: sep_index]
                    alias_value = alias[sep_index + 2:].lstrip()
                    new_row[alias_key] = alias_value
                new_rows.append(new_row)
            alias_df = pd.DataFrame(new_rows)
            alias_df.fillna('')
            new_table = pd.concat([table.drop('aliases', axis=1), alias_df], axis=1)
            return new_table
        reactions_table = expand_aliases(reactions_table)
        compounds_table = expand_aliases(compounds_table)

        progress.update("Saving reorganized tables")
        reactions_table.to_csv(reactions_path, sep='\t', index=None)
        compounds_table.to_csv(compounds_path, sep='\t', index=None)
        progress.end()

        run.info("ModelSEED database version (git commit hash)", sha)
        run.info("Reorganized ModelSEED reactions table", reactions_path)
        run.info("Reorganized ModelSEED compounds table", compounds_path)

class Constructor:
    """Construct metabolic reaction network objects."""
    def __init__(
        self,
        ko_dir: str = None,
        modelseed_dir: str = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        self.ko_dir = ko_dir
        self.modelseed_dir = modelseed_dir
        self.run = run
        self.progress = progress

    # def import_network(self, json: str) -> ReactionNetwork:
    #     """Import a metabolic model JSON file as a reaction network object."""
    #     pass

    def load_network(
        self,
        contigs_db: str = None,
        genomes_storage_db: str = None,
        pan_db: str = None
    ) -> ReactionNetwork:
        """
        Load a reaction network stored in a database as a reaction network object.

        Parameters
        ==========
        contigs_db : str, None
            Path to a contigs database in which a reaction network is stored.

        genomes_storage_db: str, None
            Path to a genomes storage database in which KO annotations are stored. 'pan_db' is also
            required.

        pan_db: str, None
            Path to a pan database in which a reaction network is stored. 'genomes_storage_db' is
            also required.

        Returns
        =======
        ReactionNetwork
            Reaction network loaded from the input database.
        """
        # Check that the reaction network stored in a database is derived from the current gene KO
        # annotations in the database.
        if contigs_db:
            network = self.load_contigs_database_network(contigs_db)
        elif genomes_storage_db or pan_db:
            network = self.load_pangenomic_network(genomes_storage_db=genomes_storage_db, pan_db=pan_db)
        else:
            raise ConfigError(
                "A reaction network must be loaded from a database source. "
                "Either a contigs database or a genomes storage database and pan database are required."
            )
        return network

    def load_contigs_database_network(
        self,
        contigs_db: str,
        check_gene_annotations: bool = True
    ) -> GenomicNetwork:
        """
        Load reaction network data stored in a contigs database as a reaction network object.

        Parameters
        ==========
        contigs_db : str
            Path to a contigs database in which a reaction network is stored.

        check_gene_annotations : bool, True
            If True, check that the reaction network stored in the contigs database was made from
            the same set of KEGG KO gene annotations as currently in the database, and throw an
            error if this is not the case. If False, allow the stored reaction network to have been
            made from a different set of KO gene annotations than is currently stored in the
            database. This can result in different genes being associated with KOs in the returned
            GenomicNetwork than in the original network that was stored.

        Returns
        =======
        GenomicNetwork
            Reaction network loaded from the contigs database.
        """
        # Load the contigs database.
        utils.is_contigs_db(contigs_db)
        args = argparse.Namespace()
        args.contigs_db = contigs_db
        contigs_super = ContigsSuperclass(args, r=run_quiet)
        contigs_super.init_functions(requested_sources=['KOfam'])

        # Check that the network stored in the contigs database was made from the same set of KEGG
        # KO gene annotations as currently in the database.
        stored_hash = contigs_super.a_meta['reaction_network_ko_annotations_hash']
        current_hash = self.hash_ko_annotations(contigs_super.gene_function_calls_dict)
        if check_gene_annotations:
            if stored_hash != current_hash:
                ConfigError(
                    "The reaction network stored in the contigs database was made from a different set "
                    "of KEGG KO gene annotations than is currently in the database. There are two "
                    "solutions to this problem. First, 'anvi-reaction-network' can be run again to "
                    "overwrite the existing network stored in the database with a new network from the "
                    "new KO gene annotations. Second, 'check_gene_annotations' can be made False rather "
                    "than True, allowing the stored network to have been made from a different set of KO "
                    "gene annotations than is currently stored in the database. This can result in "
                    "different genes being associated with KOs in the returned GenomicNetwork than in "
                    "the original network that was stored. The available version of the KO database is "
                    "used to fill in data for KOs in the network that are not current gene annotations."
                )
                # The KO database is needed if KOs in the stored network are among the current gene annotations.
                ko_db = KODatabase(ko_dir=self.ko_dir)
        else:
            if stored_hash != current_hash:
                self.run.warning(
                    "The reaction network stored in the contigs database was made from a different set "
                    "of KEGG KO gene annotations than is currently in the database. This will be "
                    "ignored since 'check_gene_annotations' is False. This can result in different genes "
                    "being associated with KOs in the returned GenomicNetwork than in the original "
                    "network that was stored."
                )

        network = GenomicNetwork()

        cdb = ContigsDatabase(contigs_db)
        functions_table = cdb.db.get_table_as_dataframe('gene_functions', where_clause='source = "KOfam"')
        reactions_table = cdb.db.get_table_as_dataframe('gene_function_reactions')
        metabolites_table = cdb.db.get_table_as_dataframe('gene_function_metabolites')
        cdb.disconnect()

        # Make gene objects for all genes with KO annotations in the contigs database, including
        # genes that are not in the network, which are later removed from the network.
        for gcid, ko_id, ko_name, e_value in zip(functions_table['gene_callers_id'], functions_table['accession'], functions_table['function'], functions_table['e_value']):
            try:
                # This is not the first annotation involving the gene, so an object for it already exists.
                gene = network.genes[gcid]
            except KeyError:
                gene = Gene()
                gene.gcid = gcid
                network.genes[gcid] = gene
            try:
                # This is not the first annotation involving the KO, so an object for it already exists.
                ko = network.kos[ko_id]
            except KeyError:
                ko = KO()
                ko.id = ko_id
                ko.name = ko_name
                network.kos[ko_id] = ko
            gene.kos.append(ko)
            gene.e_values.append(e_value)

        # Make ModelSEED reaction objects. Relate them to KOs through KEGG REACTION and EC number aliases.
        for row in reactions_table.itertuples():
            # Each row of the table contains information on a different ModelSEED reaction.
            reaction = ModelSEEDReaction()
            modelseed_reaction_id: str = row.modelseed_reaction_id
            reaction.modelseed_id = modelseed_reaction_id
            reaction.modelseed_name = row.modelseed_reaction_name
            network.reactions[modelseed_reaction_id] = reaction

            modelseed_compound_ids: str = row.metabolite_modelseed_ids
            reaction.compounds = []
            for modelseed_compound_id in modelseed_compound_ids.split(', '):
                try:
                    compound = network.metabolites[modelseed_compound_id]
                except KeyError:
                    # This is not the first reaction involving the compound, so an object for it already exists.
                    compound = ModelSEEDCompound()
                    compound.modelseed_id = modelseed_compound_id
                    reaction.compounds.append(compound)
                    network.metabolites[modelseed_compound_id] = compound
            reaction.compounds = tuple(reaction.compounds)

            stoichiometry: str = row.stoichiometry
            reaction.coefficients = tuple(int(coeff) for coeff in stoichiometry.split(', '))
            compartments: str = row.compartments
            reaction.compartments = tuple(compartments.split(', '))
            reversibility: int = row.reversibility
            reaction.reversibility = bool(reversibility)

            # Map KEGG reaction aliases of the ModelSEED reaction to all KOs that were associated with the KEGG reaction.
            kegg_reaction_ko_ids: Dict[str, List[str]] = {}
            kegg_reaction_sources: str = row.ko_kegg_reaction_source
            for kegg_reaction_item in kegg_reaction_sources.split('; '):
                if not kegg_reaction_item:
                    # The ModelSEED reaction was not sourced from KEGG reactions.
                    continue
                kegg_reaction_id, ko_ids = kegg_reaction_item.split(': (')
                ko_ids = ko_ids[:-1].split(', ')
                kegg_reaction_ko_ids[kegg_reaction_id] = ko_ids
            # Record *all* KEGG reaction aliases of the ModelSEED reaction, including those not associated with KO annotations.
            other_kegg_reaction_ids: str = row.other_kegg_reaction_ids
            reaction.kegg_aliases = list(kegg_reaction_ko_ids)
            if other_kegg_reaction_ids:
                reaction.kegg_aliases += other_kegg_reaction_ids.split(', ')
            reaction.kegg_aliases = tuple(reaction.kegg_aliases)

            network.modelseed_kegg_aliases[modelseed_reaction_id] = modelseed_kegg_aliases = []
            orphan_ko_ids = []
            reaction_added_to_ko = False
            for kegg_reaction_id, ko_ids in kegg_reaction_ko_ids.items():
                # Record the ModelSEED reaction as one of the aliases of the KEGG reaction in the network.
                try:
                    network.kegg_modelseed_aliases[kegg_reaction_id].append(modelseed_reaction_id)
                except KeyError:
                    network.kegg_modelseed_aliases[kegg_reaction_id] = [modelseed_reaction_id]
                modelseed_kegg_aliases.append(kegg_reaction_id)
                for ko_id in ko_ids:
                    try:
                        ko = network.kos[ko_id]
                    except KeyError:
                        # This only happens when the current set of gene KO annotations in the contigs
                        # database does not match the set from which the reaction network was originally
                        # made, and the KO in the network is no longer a gene annotation in the database.
                        ko = KO()
                        ko.ko_id = ko_id
                        # The KO name is unknown from the contigs database, so take it from the KO database.
                        ko.ko_name = ko_db.ko_table.loc[ko_id, 'name']
                        network.kos[ko_id] = ko
                        orphan_ko_ids.append(ko_id)
                    if not reaction_added_to_ko:
                        # This is the first encounter with the reaction for the KO.
                        ko.reactions[modelseed_reaction_id] = reaction
                        reaction_added_to_ko = True
                    try:
                        ko.kegg_reaction_aliases[modelseed_reaction_id].append(kegg_reaction_id)
                    except KeyError:
                        ko.kegg_reaction_aliases[modelseed_reaction_id] = [kegg_reaction_id]

            # Map EC number aliases of the ModelSEED reaction to all KOs that were associated with the EC number.
            ec_number_ko_ids: Dict[str, List[str]] = {}
            ec_number_sources: str = row.ko_ec_number_source
            for ec_number_item in ec_number_sources.split('; '):
                if not ec_number_item:
                    # The ModelSEED reaction was not sourced from EC numbers.
                    continue
                ec_number, ko_ids = ec_number_item.split(': (')
                ko_ids = ko_ids[:-1].split(', ')
                ec_number_ko_ids[ec_number] = ko_ids
            # Record *all* EC number aliases of the ModelSEED reaction, including those not associated with KO annotations.
            other_ec_numbers: str = row.other_ec_numbers
            reaction.ec_number_aliases = list(ec_number_ko_ids)
            if other_ec_numbers:
                reaction.ec_number_aliases += other_ec_numbers.split(', ')
            reaction.ec_number_aliases = tuple(reaction.ec_number_aliases)

            network.modelseed_ec_number_aliases[modelseed_reaction_id] = modelseed_ec_number_aliases = []
            for ec_number, ko_ids in ec_number_ko_ids.items():
                # Record the ModelSEED reaction as one of the aliases of the EC number in the network.
                try:
                    network.ec_number_modelseed_aliases[ec_number].append(modelseed_reaction_id)
                except KeyError:
                    network.ec_number_modelseed_aliases[ec_number] = [modelseed_reaction_id]
                modelseed_ec_number_aliases.append(ec_number)
                for ko_id in ko_ids:
                    try:
                        ko = network.kos[ko_id]
                    except KeyError:
                        # This only happens when the current set of gene KO annotations in the contigs
                        # database does not match the set from which the reaction network was originally
                        # made, and the KO in the network is no longer a gene annotation in the database.
                        ko = KO()
                        ko.ko_id = ko_id
                        # The KO name is unknown from the contigs database, so take it from the KO database.
                        ko.ko_name = ko_db.ko_table.loc[ko_id, 'name']
                        network.kos[ko_id] = ko
                        orphan_ko_ids.append(ko_id)
                    if not reaction_added_to_ko:
                        # This is the first encounter with the reaction for the KO.
                        ko.reactions[modelseed_reaction_id] = reaction
                        reaction_added_to_ko = True
                    try:
                        ko.ec_number_aliases[modelseed_reaction_id].append(ec_number)
                    except KeyError:
                        ko.ec_number_aliases[modelseed_reaction_id] = [ec_number]

            if DEBUG:
                self.run.info_single(
                    "The following KOs are in the stored reaction network in the contigs database "
                    "but are not among the gene KO annotations; the available version of the KO database "
                    f"was used to fill in the function names of these KOs: {', '.join(orphan_ko_ids)}"
                )

        # Assign ModelSEED compound object attributes.
        for row in metabolites_table.itertuples():
            modelseed_compound_id = row.modelseed_compound_id
            compound = network.metabolites[modelseed_compound_id]
            modelseed_compound_name: str = row.modelseed_compound_name
            compound.modelseed_name = modelseed_compound_name
            kegg_aliases: str = row.kegg_aliases
            compound.kegg_aliases = tuple(kegg_aliases.split(', '))
            # Compounds without a formula, recorded here as None, have a nominal charge of 10000000
            # in the ModelSEED compounds database. This is replaced by NaN in the table and here as
            # None in the reaction network.
            formula: str = row.formula
            compound.formula = formula
            charge: int = row.charge
            compound.charge = charge if not np.isnan(charge) else None

        # Remove any trace of genes that do not contribute to the reaction network. Also remove
        # unnetworked KO links to genes.
        unnetworked_gcids = []
        for gcid, gene in network.genes.items():
            gene_in_network = False
            unnetworked_ko_indices: List[int] = []
            idx = 0
            for ko, e_value in zip(gene.kos, gene.e_values):
                if ko.reactions:
                    gene_in_network = True
                else:
                    unnetworked_ko_indices.append(idx)
                idx += 1
            if gene_in_network:
                for unnetworked_ko_index in reversed(unnetworked_ko_indices):
                    gene.kos.pop(unnetworked_ko_index)
                    gene.e_values.pop(unnetworked_ko_index)
            else:
                unnetworked_gcids.append(gcid)
        for gcid in unnetworked_gcids:
            network.genes.pop(gcid)

        # Remove any trace of KOs that do not contribute to the reaction network.
        unnetworked_ko_ids = []
        for ko_id, ko in network.kos.items():
            if not ko.reactions:
                unnetworked_ko_ids.append(ko_id)
        for ko_id in unnetworked_ko_ids:
            network.kos.pop(ko_id)

        # Remove entries showing KO KEGG reaction aliases of ModelSEED reactions where there are no aliases.
        modelseed_reaction_ids = []
        for modelseed_reaction_id, kegg_reaction_ids in network.modelseed_kegg_aliases.items():
            if not kegg_reaction_ids:
                modelseed_reaction_ids.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids:
            network.modelseed_kegg_aliases.pop(modelseed_reaction_id)

        # Remove entries showing EC number aliases of ModelSEED reactions where there are no aliases.
        modelseed_reaction_ids = []
        for modelseed_reaction_id, ec_numbers in network.modelseed_ec_number_aliases.items():
            if not ec_numbers:
                modelseed_reaction_ids.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids:
            network.modelseed_ec_number_aliases.pop(modelseed_reaction_id)

        return network

    def make_network(
        self,
        contigs_db: str = None,
        genomes_storage_db: str = None,
        pan_db: str = None,
        store: bool = True,
        overwrite_existing_network: bool = False
    ) -> ReactionNetwork:
        """
        Make a metabolic reaction network from KEGG Orthologs stored in an anvi'o database,
        associated KEGG annotations, and the ModelSEED Biochemistry database.

        Parameters
        ==========
        contigs_db : str, None
            Path to a contigs database. The database can represent different types of samples,
            including a single genome, metagenome, or transcriptome. The network is derived from
            gene KO annotations stored in the database. If 'store' is True, the network is saved in
            the database.

        genomes_storage_db : str, None
            Path to a genomes storage database. The pangenomic network is derived from gene KO
            annotations stored in the database. 'pan_db' is also required.

        pan_db : str, None
            Path to a pan database. The pangenomic network is determined for gene clusters stored in
            the database. If 'store' is True, the network is saved in the database.
            'genomes_storage_db' is also required.

        store : bool, True
            Save the network. A network constructed from a contigs database is stored in that
            database. A pangenomic network constructed from a genomes stroage database and pan
            database is stored in the pan database.

        overwrite_existing_network : bool, False
            Overwrite an existing network stored in the contigs or pan database. 'store' is also
            required.
        """
        if contigs_db:
            self.run.info_single(
                "A reaction network will be made from protein orthology annotations in the contigs database."
            )
            network = self.make_contigs_database_network(contigs_db, store=store, overwrite_existing_network=overwrite_existing_network)
        elif genomes_storage_db or pan_db:
            self.run.info_single(
                "A pangenomic reaction network will be made from protein orthology annotations "
                "in the genomes storage database and gene clusters in the pan database."
            )
            network = self.make_pangenomic_network(genomes_storage_db, pan_db, store=store, overwrite_existing_network=overwrite_existing_network)
        else:
            raise ConfigError(
                "A reaction network cannot be made without a database source. "
                "Either a contigs database or a genomes storage database and pan database are required."
            )
        return network

    def make_contigs_database_network(
        self,
        contigs_db: str,
        store: bool = True,
        overwrite_existing_network: bool = False
    ) -> GenomicNetwork:
        """
        Make a metabolic reaction network from KEGG Orthologs stored in a contigs database.

        Parameters
        ==========
        contigs_db : str
            Path to a contigs database. The database can represent different types of samples,
            including a single genome, metagenome, or transcriptome. The network is derived from
            gene KO annotations stored in the database.

        store : bool, True
            Save the network to the contigs database.

        overwrite_existing_network : bool, False
            Overwrite an existing network stored in the contigs database. 'store' is also required.

        Returns
        =======
        GenomicNetwork
            The network derived from the contigs database.
        """
        # Load the contigs database.
        self.run.info("Contigs database", contigs_db)
        utils.is_contigs_db(contigs_db)
        args = argparse.Namespace()
        args.contigs_db = contigs_db
        contigs_super = ContigsSuperclass(args, r=run_quiet)
        if store and contigs_super.a_meta['reaction_network_ko_annotations_hash'] and not overwrite_existing_network:
            raise ConfigError(
                "The existing reaction network in the contigs database must be explicitly overwritten."
            )
        contigs_super.init_functions(requested_sources=['KOfam'])

        self.progress.new("Building reaction network")
        self.progress.update("Loading reference databases")

        ko_db = KODatabase(self.ko_dir)
        modelseed_db = ModelSEEDDatabase(self.modelseed_dir)

        network = GenomicNetwork()

        modelseed_kegg_reactions_table = modelseed_db.kegg_reactions_table
        modelseed_ec_reactions_table = modelseed_db.ec_reactions_table
        modelseed_compounds_table = modelseed_db.compounds_table

        # List KOs that annotated genes in the contigs database but for some reason are not found in
        # the KO database.
        undefined_ko_ids = []

        # Parse gene-KO matches recorded in the contigs database.
        gene_function_calls_dict: Dict = contigs_super.gene_function_calls_dict
        total_ko_matches = len(gene_function_calls_dict)
        num_ko_matches_parsed = -1
        for gcid, gene_dict in gene_function_calls_dict.items():
            num_ko_matches_parsed += 1
            self.progress.update(f"Gene-KO matches parsed: {num_ko_matches_parsed} / {total_ko_matches}")

            if gcid in network.genes:
                # An object representing the gene was already added to the network.
                gene = network.genes[gcid]
            else:
                gene = Gene()
                gene.gcid = gcid
                # Add the gene to the network, regardless of whether it yields reactions. Genes not
                # contributing to the reaction network are removed later.
                network.genes[gcid] = gene

            ko_data = gene_dict['KOfam']
            gene.e_values.append(float(ko_data[2]))
            ko_id = ko_data[0]
            if ko_id in network.kos:
                # The KO was associated with an already encountered gene and added to the network.
                # Objects representing ModelSEED reactions and metabolites and other data associated
                # with the KO were added to the network as well.
                gene.kos.append(network.kos[ko_id])
                continue
            else:
                ko = KO()
                ko.id = ko_id
                ko.name = ko_data[1]
                gene.kos.append(ko)
                # Add the KO to the network, regardless of whether it yields reactions. KOs not
                # contributing to the network are removed later.
                network.kos[ko_id] = ko

            # Find KEGG reactions and EC numbers associated with the newly encountered KO.
            try:
                ko_info = ko_db.ko_table.loc[ko.id]
            except KeyError:
                undefined_ko_ids.append(ko_id)
                continue
            ko_kegg_reaction_info: str = ko_info.loc['reactions']
            if pd.isna(ko_kegg_reaction_info):
                # The KO is not associated with KEGG reactions.
                ko_kegg_reaction_ids = []
            else:
                ko_kegg_reaction_ids = ko_kegg_reaction_info.split()
            ko_ec_number_info: str = ko_info.loc['ec_numbers']
            if pd.isna(ko_ec_number_info):
                # The KO is not associated with EC numbers.
                ko_ec_numbers = []
            else:
                ko_ec_numbers = ko_ec_number_info.split()

            if not (ko_kegg_reaction_ids or ko_ec_numbers):
                # The KO is not associated with any KEGG reactions or EC numbers, and thereby cannot
                # be associated with ModelSEED reactions.
                continue

            # If a KEGG reaction has already been encountered, then aliased ModelSEED reactions have
            # been processed and added as ModelSEEDReaction objects to the network. Therefore, KEGG
            # reactions that have already been encountered are treated differently than KEGG
            # reactions encountered for the first time.
            new_kegg_reaction_ids = []
            for kegg_reaction_id in ko_kegg_reaction_ids:
                try:
                    # The KEGG reaction has already been encountered. Retrieve ModelSEED reactions
                    # aliased by the KEGG reaction.
                    modelseed_reaction_ids = network.kegg_modelseed_aliases[kegg_reaction_id]
                except KeyError:
                    new_kegg_reaction_ids.append(kegg_reaction_id)
                    # The following list of ModelSEED reaction IDs associated with the KEGG reaction
                    # is filled in later. If no ModelSEED reactions are associated with the KEGG
                    # reaction, the entry in the dictionary will be removed.
                    network.kegg_modelseed_aliases[kegg_reaction_id] = []
                    continue
                for modelseed_reaction_id in modelseed_reaction_ids:
                    try:
                        # Retrieve the existing ModelSEEDReaction object.
                        reaction = network.reactions[modelseed_reaction_id]
                    except KeyError:
                        # The ModelSEED reaction associated with the EC number did not have valid
                        # data: for example, when the 'stoichiometry' field is empty.
                        continue
                    # Associate the ModelSEED reaction with the newly encountered KO.
                    ko.reactions[modelseed_reaction_id] = reaction
                    # Record which KEGG REACTION IDs and EC numbers from the KO yield the ModelSEED reaction.
                    ko.kegg_reaction_aliases[modelseed_reaction_id] = list(
                        set(ko_kegg_reaction_ids).intersection(set(reaction.kegg_aliases))
                    )
                    ko.ec_number_aliases[modelseed_reaction_id] = list(
                        set(ko_ec_numbers).intersection(set(reaction.ec_number_aliases))
                    )

            # As above with KEGG reactions, if an EC number has already been encountered, then
            # aliased ModelSEED reactions have been processed and added as ModelSEEDReaction objects
            # to the network. Therefore, EC numbers that have already been encountered are treated
            # differently than EC numbers encountered for the first time.
            new_ec_numbers = []
            for ec_number in ko_ec_numbers:
                try:
                    # The EC number has already been encountered. Retrieve ModelSEED reactions
                    # aliased by the EC number.
                    modelseed_reaction_ids = network.ec_number_modelseed_aliases[ec_number]
                except KeyError:
                    new_ec_numbers.append(ec_number)
                    # The following list of ModelSEED reaction IDs associated with the EC number is
                    # filled in later. If no ModelSEED reactions are associated with the EC number,
                    # the entry in the dictionary will be removed.
                    network.ec_number_modelseed_aliases[ec_number] = []
                    continue
                for modelseed_reaction_id in modelseed_reaction_ids:
                    try:
                        # Retrieve the existing ModelSEEDReaction object.
                        reaction = network.reactions[modelseed_reaction_id]
                    except KeyError:
                        # The ModelSEED reaction associated with the EC number did not have valid
                        # data: for example, when the 'stoichiometry' field is empty.
                        continue
                    if modelseed_reaction_id in reaction.ec_number_aliases:
                        # A KEGG reaction associated with the newly encountered KO was also
                        # associated with the ModelSEED reaction. KO EC number aliases were
                        # previously recorded along with KO KEGG reaction aliases. Redundant work
                        # can be avoided here linking the ModelSEED reaction to the KO in the network.
                        continue
                    ko.reactions[modelseed_reaction_id] = reaction
                    ko.kegg_reaction_aliases[modelseed_reaction_id] = list(
                        set(ko_kegg_reaction_ids).intersection(set(reaction.kegg_aliases))
                    )
                    ko.ec_number_aliases[modelseed_reaction_id] = list(
                        set(ko_ec_numbers).intersection(set(reaction.ec_number_aliases))
                    )

            if not (new_kegg_reaction_ids or new_ec_numbers):
                # All of the KEGG reactions and EC numbers associated with the KO have already been
                # encountered in previously processed KOs and added to the network, so proceed to
                # the next gene KO annotation.
                continue

            # Get data on ModelSEED reactions aliased by newly encountered KEGG REACTION IDs and EC numbers.
            modelseed_reactions_data = {}
            if new_kegg_reaction_ids:
                # Each row of the table represents a unique KEGG reaction -> ModelSEED reaction mapping.
                modelseed_kegg_reactions_dict: Dict[str, Dict] = modelseed_kegg_reactions_table[
                    modelseed_kegg_reactions_table['KEGG_REACTION_ID'].isin(new_kegg_reaction_ids)
                ].to_dict(orient='index')
                for modelseed_reaction_data in modelseed_kegg_reactions_dict.values():
                    kegg_reaction_id = modelseed_reaction_data['KEGG_REACTION_ID']
                    modelseed_reaction_id = modelseed_reaction_data['id']
                    # Record the association between the KEGG reaction and ModelSEED reaction in the
                    # network, and vice versa.
                    network.kegg_modelseed_aliases[kegg_reaction_id].append(modelseed_reaction_id)
                    try:
                        network.modelseed_kegg_aliases[modelseed_reaction_id].append(kegg_reaction_id)
                    except KeyError:
                        # This is the first time the ModelSEED reaction has been encountered.
                        network.modelseed_kegg_aliases[modelseed_reaction_id] = [kegg_reaction_id]
                        network.modelseed_ec_number_aliases[modelseed_reaction_id] = []
                    if modelseed_reaction_id in modelseed_reactions_data:
                        # One of the other newly encountered KEGG reactions also mapped to this
                        # ModelSEED reaction, so do not record redundant ModelSEED reaction data.
                        continue
                    modelseed_reactions_data[modelseed_reaction_id] = modelseed_reaction_data
            if new_ec_numbers:
                # Each row of the table represents a unique EC number -> ModelSEED reaction mapping.
                modelseed_ec_reactions_dict: Dict[str, Dict] = modelseed_ec_reactions_table[
                    modelseed_ec_reactions_table['EC_number'].isin(new_ec_numbers)
                ].to_dict(orient='index')
                for modelseed_reaction_data in modelseed_ec_reactions_dict.values():
                    ec_number = modelseed_reaction_data['EC_number']
                    modelseed_reaction_id = modelseed_reaction_data['id']
                    # Record the association between the EC number and ModelSEED reaction in the
                    # network, and vice versa.
                    network.ec_number_modelseed_aliases[ec_number].append(modelseed_reaction_id)
                    try:
                        network.modelseed_ec_number_aliases[modelseed_reaction_id].append(ec_number)
                    except KeyError:
                        # This is the first time the ModelSEED reaction has been encountered.
                        network.modelseed_ec_number_aliases[modelseed_reaction_id] = [ec_number]
                        network.modelseed_kegg_aliases[modelseed_reaction_id] = []
                    if modelseed_reaction_id in modelseed_reactions_data:
                        # One of the other newly encountered KEGG reactions or EC numbers also
                        # mapped to this ModelSEED reaction, so do not record redundant ModelSEED reaction data.
                        continue
                    modelseed_reactions_data[modelseed_reaction_id] = modelseed_reaction_data
            if not modelseed_reactions_data:
                # The newly encountered KEGG REACTION IDs and EC numbers do not map to ModelSEED
                # reactions (are not in the table).
                continue

            # Process the ModelSEED reactions aliased by newly encountered KEGG reactions and EC numbers.
            for modelseed_reaction_id, modelseed_reaction_data in modelseed_reactions_data.items():
                if modelseed_reaction_id in network.reactions:
                    # The ModelSEED reaction is aliased by previously encountered KEGG reactions and
                    # EC numbers, and so has already been added to the network.
                    continue
                # Make a new reaction object for the ModelSEED ID. This object does not yet have
                # metabolite objects (for the ModelSEED compound IDs) added to it yet.
                reaction, modelseed_compound_ids = self._get_modelseed_reaction(modelseed_reaction_data)
                if reaction is None:
                    # For some reason, the reaction does not have a equation in the ModelSEED
                    # database. Associations between such reactions without equations and sourcing
                    # KEGG reactions and EC numbers are later removed from the network attributes,
                    # 'kegg_modelseed_aliases', 'ec_number_modelseed_aliases',
                    # 'modelseed_kegg_aliases', and 'modelseed_ec_number_aliases'.
                    continue
                ko.reactions[modelseed_reaction_id] = reaction
                # Record which KEGG REACTION IDs and EC numbers from the KO yield the ModelSEED reaction.
                ko.kegg_reaction_aliases[modelseed_reaction_id] = list(
                    set(new_kegg_reaction_ids).intersection(set(reaction.kegg_aliases))
                )
                ko.ec_number_aliases[modelseed_reaction_id] = list(
                    set(new_ec_numbers).intersection(set(reaction.ec_number_aliases))
                )
                network.reactions[modelseed_reaction_id] = reaction

                # If the ModelSEED compound ID has been encountered in previously processed
                # reactions, then there is already a ModelSEEDCompound object for it.
                new_modelseed_compound_ids = []
                reaction_compounds = []
                for modelseed_compound_id in modelseed_compound_ids:
                    if modelseed_compound_id in network.metabolites:
                        reaction_compounds.append(network.metabolites[modelseed_compound_id])
                    else:
                        new_modelseed_compound_ids.append(modelseed_compound_id)

                # Generate new metabolite objects in the network
                for modelseed_compound_id in new_modelseed_compound_ids:
                    try:
                        modelseed_compound_series: pd.Series = modelseed_compounds_table.loc[modelseed_compound_id]
                    except KeyError:
                        raise ConfigError(
                            f"A row for the ModelSEED compound ID, '{modelseed_compound_id}', was expected "
                            "but not found in the ModelSEED compounds table. This ID was found in the equation "
                            f"for the ModelSEED reaction, '{modelseed_reaction_id}'."
                        )
                    modelseed_compound_data = modelseed_compound_series.to_dict()
                    modelseed_compound_data['id'] = modelseed_compound_id
                    compound = self._get_modelseed_compound(modelseed_compound_data)
                    reaction_compounds.append(compound)
                    network.metabolites[modelseed_compound_id] = compound
                reaction.compounds = tuple(reaction_compounds)

        # List genes that do not contribute to the reaction network. Remove any trace of these genes
        # from the network.
        unnetworked_gcids = []
        for gcid, gene in network.genes.items():
            for ko in gene.kos:
                if ko.reactions:
                    break
            else:
                unnetworked_gcids.append(gcid)
        for gcid in unnetworked_gcids:
            network.genes.pop(gcid)

        # List KOs that do not contribute to the reaction network. Remove any trace of these KOs
        # from the network.
        unnetworked_ko_ids = []
        for ko_id, ko in network.kos.items():
            if not ko.reactions:
                unnetworked_ko_ids.append(ko_id)
        for ko_id in unnetworked_ko_ids:
            network.kos.pop(ko_id)

        # List KO KEGG reactions that do not map to ModelSEED reactions. Remove any trace of these
        # KEGG reactions from the network.
        unnetworked_kegg_reaction_ids = []
        for kegg_reaction_id, modelseed_reaction_ids in network.kegg_modelseed_aliases.items():
            if not modelseed_reaction_ids:
                unnetworked_kegg_reaction_ids.append(kegg_reaction_id)
        for kegg_reaction_id in unnetworked_kegg_reaction_ids:
            network.kegg_modelseed_aliases.pop(kegg_reaction_id)

        # List KO EC numbers that do not map to ModelSEED reactions. Remove any trace of these EC
        # numbers from the network.
        unnetworked_ec_numbers = []
        for ec_number, modelseed_reaction_ids in network.ec_number_modelseed_aliases.items():
            if not modelseed_reaction_ids:
                unnetworked_ec_numbers.append(ec_number)
        for ec_number in unnetworked_ec_numbers:
            network.ec_number_modelseed_aliases.pop(ec_number)

        # List aliased ModelSEED reactions that did not yield a ModelSEEDReaction object due to the
        # lack of an equation for the reaction in the ModelSEED database. Remove any trace of these
        # reactions from the network.
        undefined_modelseed_reaction_ids = list(
            set(network.modelseed_kegg_aliases).difference(set(network.reactions))
        )
        for modelseed_reaction_id in undefined_modelseed_reaction_ids:
            network.modelseed_kegg_aliases.pop(modelseed_reaction_id)
            network.modelseed_ec_number_aliases.pop(modelseed_reaction_id)
        self.progress.end()

        if DEBUG:
            self.run.info_single(
                "The following ModelSEED reactions would have been added to the reaction network "
                "had there been a chemical equation in the ModelSEED database; perhaps it is worth "
                "investigating the ModelSEED reactions table to understand why this is not the case: "
                f"{', '.join(undefined_modelseed_reaction_ids)}"
            )

        if undefined_ko_ids:
            self.run.info_single(
                "Certain genes matched KOs that were not found in the KO database. It could be that the "
                "KOfams used to annotate genes were not from the same KEGG database version as the KO files. "
                f"Here are the unrecognized KO IDs from the contigs database: {', '.join(undefined_ko_ids)}"
            )

        if store:
            if contigs_super.a_meta['reaction_network_ko_annotations_hash']:
                self.run.warning("Deleting existing reaction network from contigs database")
                cdb = ContigsDatabase(contigs_db)
                cdb.db._exec(f'''DELETE from {tables.gene_function_reactions_table_name}''')
                cdb.db._exec(f'''DELETE from {tables.gene_function_metabolites_table_name}''')
                cdb.disconnect()
                self.run.info_single("Deleted data in gene function reactions and metabolites tables", nl_after=1)

            self.progress.new("Saving reaction network to contigs database")
            self.progress.update("Reactions table")
            self._store_contigs_database_reactions(network, contigs_db)
            self.progress.update("Metabolites table")
            self._store_contigs_database_metabolites(network, contigs_db)

            self.progress.update("Metadata")
            ko_annotations_hash = self.hash_ko_annotations(gene_function_calls_dict)
            cdb = ContigsDatabase(contigs_db)
            cdb.db.set_meta_value('reaction_network_ko_annotations_hash', ko_annotations_hash)
            cdb.db.set_meta_value('reaction_network_kegg_database_release', ko_db.release)
            cdb.db.set_meta_value('reaction_network_modelseed_database_sha', modelseed_db.sha)
            cdb.disconnect()
            self.progress.end()

        stats = {}
        self.run.info_single("METABOLIC REACTION NETWORK STATISTICS", mc='green', nl_after=1)

        self.progress.new("Counting genes and KOs")
        self.progress.update("...")

        cdb = ContigsDatabase(contigs_db)
        stats['Total gene calls in genome'] = gene_call_count = cdb.db.get_row_counts_from_table('genes_in_contigs')
        cdb.disconnect()
        stats['Genes annotated with protein KOs'] = ko_annotated_gene_count = len(network.genes) + len(unnetworked_gcids)
        stats['Genes in network'] = networked_gene_count = len(network.genes)
        stats['Protein KOs annotating genes'] = annotating_ko_count = len(network.kos) + len(unnetworked_ko_ids)
        stats['KOs in network'] = networked_ko_count = len(network.kos)
        self.progress.end()

        self.run.info_single("Gene calls and KEGG Ortholog (KO) annotations")
        self.run.info("Total gene calls in genome", gene_call_count)
        self.run.info("Genes annotated with protein KOs", ko_annotated_gene_count)
        self.run.info("Genes in network", networked_gene_count)
        self.run.info("Protein KOs annotating genes", annotating_ko_count)
        self.run.info("KOs in network", networked_ko_count, nl_after=1)

        self.progress.new("Counting reactions and KO sources")
        self.progress.update("...")

        stats['Reactions in network'] = reaction_count = len(network.reactions)
        reaction_counts = []
        for ko in network.kos.values():
            reaction_counts.append(len(ko.reactions))
        stats['Mean reactions per KO'] = mean_reactions_per_ko = round(np.mean(reaction_counts), 1)
        stats['Stdev reactions per KO'] = std_reactions_per_ko = round(np.std(reaction_counts), 1)
        stats['Max reactions per KO'] = max_reactions_per_ko = max(reaction_counts)
        self.progress.end()

        self.run.info_single("ModelSEED reactions in network and KO sources")
        self.run.info("Reactions in network", reaction_count)
        self.run.info("Mean reactions per KO", mean_reactions_per_ko)
        self.run.info("Stdev reactions per KO", std_reactions_per_ko)
        self.run.info("Max reactions per KO", max_reactions_per_ko, nl_after=1)

        self.progress.new("Counting reactions from each alias source")
        self.progress.update("...")

        kegg_reaction_source_count = 0
        ec_number_source_count = 0
        both_source_count = 0
        for modelseed_reaction_id, kegg_reaction_ids in network.modelseed_kegg_aliases.items():
            ec_numbers = network.modelseed_ec_number_aliases[modelseed_reaction_id]
            if kegg_reaction_ids:
                kegg_reaction_source_count += 1
            if ec_numbers:
                ec_number_source_count += 1
            if kegg_reaction_ids and ec_numbers:
                both_source_count += 1
        stats['Reactions aliased by KEGG reaction'] = kegg_reaction_source_count
        stats['Reactions aliased by EC number'] = ec_number_source_count
        stats['Rxns aliased by both KEGG rxn & EC number'] = both_source_count
        stats['Reactions aliased only by KEGG reaction'] = only_kegg_reaction_source_count = kegg_reaction_source_count - both_source_count
        stats['Reactions aliased only by EC number'] = only_ec_number_source_count = ec_number_source_count - both_source_count

        stats['KEGG reactions contributing to network'] = kegg_reaction_count = len(network.kegg_modelseed_aliases)
        reaction_counts = []
        for kegg_reaction_id, modelseed_reaction_ids in network.kegg_modelseed_aliases.items():
            reaction_counts.append(len(modelseed_reaction_ids))
        stats['Mean reactions per KEGG reaction'] = mean_reactions_per_kegg_reaction = round(np.mean(reaction_counts), 1)
        stats['Stdev reactions per KEGG reaction'] = std_reactions_per_kegg_reaction = round(np.std(reaction_counts), 1)
        stats['Max reactions per KEGG reaction'] = max_reactions_per_kegg_reaction = max(reaction_counts)

        stats['EC numbers contributing to network'] = ec_number_count = len(network.ec_number_modelseed_aliases)
        reaction_counts = []
        for ec_number, modelseed_reaction_ids in network.ec_number_modelseed_aliases.items():
            reaction_counts.append(len(modelseed_reaction_ids))
        stats['Mean reactions per EC number'] = mean_reactions_per_ec_number = round(np.mean(reaction_counts), 1)
        stats['Stdev reactions per EC number'] = std_reactions_per_ec_number = round(np.std(reaction_counts), 1)
        stats['Max reactions per EC number'] = max_reactions_per_ec_number = max(reaction_counts)
        self.progress.end()

        self.run.info_single("Reaction alias source comparison")
        self.run.info("Reactions aliased by KEGG reaction", kegg_reaction_source_count)
        self.run.info("Reactions aliased by EC number", ec_number_source_count)
        self.run.info("Rxns aliased by both KEGG rxn & EC number", both_source_count)
        self.run.info("Reactions aliased only by KEGG reaction", only_kegg_reaction_source_count)
        self.run.info("Reactions aliased only by EC number", only_ec_number_source_count)
        self.run.info("KEGG reactions contributing to network", kegg_reaction_count)
        self.run.info("Mean reactions per KEGG reaction", mean_reactions_per_kegg_reaction)
        self.run.info("Stdev reactions per KEGG reaction", std_reactions_per_kegg_reaction)
        self.run.info("Max reactions per KEGG reaction", max_reactions_per_kegg_reaction)
        self.run.info("EC numbers contributing to network", ec_number_count)
        self.run.info("Mean reactions per EC number", mean_reactions_per_ec_number)
        self.run.info("Stdev reactions per EC number", std_reactions_per_ec_number)
        self.run.info("Max reactions per EC number", max_reactions_per_ec_number, nl_after=1)

        self.progress.new("Counting reactions and metabolites by property")
        self.progress.update("...")

        reversible_count = 0
        irreversible_count = 0
        cytoplasmic_compound_ids = []
        extracellular_compound_ids = []
        consumed_compound_ids = []
        produced_compound_ids = []
        compound_reaction_counts = {}
        for reaction in network.reactions.values():
            if reaction.reversibility:
                reversible_count += 1
            else:
                irreversible_count += 1
            encountered_compound_ids = []
            for compartment, coefficient, compound in zip(reaction.compartments, reaction.coefficients, reaction.compounds):
                compound_id = compound.modelseed_id
                if compartment == 'c':
                    cytoplasmic_compound_ids.append(compound_id)
                else:
                    extracellular_compound_ids.append(compound_id)
                if reaction.reversibility:
                    consumed_compound_ids.append(compound_id)
                    produced_compound_ids.append(compound_id)
                elif coefficient < 0:
                    consumed_compound_ids.append(compound_id)
                else:
                    produced_compound_ids.append(compound_id)
                if compound_id not in encountered_compound_ids:
                    try:
                        compound_reaction_counts[compound_id] += 1
                    except KeyError:
                        compound_reaction_counts[compound_id] = 1
        stats['Reversible reactions'] = reversible_count
        stats['Irreversible reactions'] = irreversible_count
        cytoplasmic_compound_ids = set(cytoplasmic_compound_ids)
        extracellular_compound_ids = set(extracellular_compound_ids)
        stats['Metabolites in network'] = metabolite_count = len(network.metabolites)
        stats['Cytoplasmic metabolites'] = cytoplasmic_count = len(cytoplasmic_compound_ids)
        stats['Extracellular metabolites'] = extracellular_count = len(extracellular_compound_ids)
        stats['Exclusively cytoplasmic metabolites'] = exclusively_cytoplasmic_count = len(cytoplasmic_compound_ids.difference(extracellular_compound_ids))
        stats['Exclusively extracellular metabolites'] = exclusively_extracellular_count = len(extracellular_compound_ids.difference(cytoplasmic_compound_ids))
        stats['Cytoplasmic/extracellular metabolites'] = cytoplasmic_plus_extracellular_count = len(cytoplasmic_compound_ids.intersection(extracellular_compound_ids))
        consumed_compound_ids = set(consumed_compound_ids)
        produced_compound_ids = set(produced_compound_ids)
        stats['Consumed metabolites'] = consumed_count = len(consumed_compound_ids)
        stats['Produced metabolites'] = produced_count = len(produced_compound_ids)
        stats['Both consumed & produced metabolites'] = consumed_plus_produced_count = len(consumed_compound_ids.intersection(produced_compound_ids))
        stats['Exclusively consumed metabolites'] = exclusively_consumed_count = len(consumed_compound_ids.difference(produced_compound_ids))
        stats['Exclusively produced metabolites'] = exclusively_produced_count = len(produced_compound_ids.difference(consumed_compound_ids))
        metabolite_reaction_counts = collections.Counter(compound_reaction_counts.values())
        stats['Metabolites consumed or produced by 1 rxns'] = one_reaction_count = metabolite_reaction_counts[1]
        stats['Metabolites consumed or produced by 2 rxns'] = two_reactions_count = metabolite_reaction_counts[2]
        stats['Metabolites consumed or produced by 3+ rxns'] = three_plus_reactions_count = metabolite_count - one_reaction_count - two_reactions_count
        self.progress.end()

        self.run.info_single("Reaction reversibility")
        self.run.info("Reversible reactions", reversible_count)
        self.run.info("Irreversible reactions", irreversible_count, nl_after=1)

        self.run.info_single("Metabolites and localization")
        self.run.info("Metabolites in network", metabolite_count)
        self.run.info("Cytoplasmic metabolites", cytoplasmic_count)
        self.run.info("Extracellular metabolites", extracellular_count)
        self.run.info("Exclusively cytoplasmic metabolites", exclusively_cytoplasmic_count)
        self.run.info("Exclusively extracellular metabolites", exclusively_extracellular_count)
        self.run.info("Cytoplasmic/extracellular metabolites", cytoplasmic_plus_extracellular_count, nl_after=1)

        self.run.info_single("Metabolite consumption and production")
        self.run.info("Consumed metabolites", consumed_count)
        self.run.info("Produced metabolites", produced_count)
        self.run.info("Both consumed & produced metabolites", consumed_plus_produced_count)
        self.run.info("Exclusively consumed metabolites", exclusively_consumed_count)
        self.run.info("Exclusively produced metabolites", exclusively_produced_count)
        self.run.info("Metabolites consumed or produced by 1 rxn", one_reaction_count)
        self.run.info("Metabolites consumed or produced by 2 rxns", two_reactions_count)
        self.run.info("Metabolites consumed or produced by 3+ rxns", three_plus_reactions_count)

        return network

    def _get_modelseed_reaction(self, modelseed_reaction_data: Dict) -> Tuple[ModelSEEDReaction, List[str]]:
        """
        Generate a ModelSEED reaction object and list of associated ModelSEED compound IDs from the
        ModelSEED reaction table entry. The reaction object is not populated with metabolite objects
        from the list of associated compound IDs.

        Parameters
        ==========
        modelseed_reaction_data : Dict
            A dictionary representation of a row for a reaction in the ModelSEED reaction table set
            up by anvi'o.

        Returns
        =======
        ModelSEEDReaction
            An object representation of the ModelSEED reaction.

        List[str]
            ModelSEED compound IDs of reactants and products.
        """
        stoichiometry: str = modelseed_reaction_data['stoichiometry']
        if pd.isna(stoichiometry):
            # ignore any reaction lacking a chemical equation for some reason
            return None, None

        reaction = ModelSEEDReaction()

        modelseed_id = modelseed_reaction_data['id']
        if pd.isna(modelseed_id):
            raise ConfigError(
                "The row for the reaction in the ModelSEED table does not but should have an ID. "
                f"Here is the data in the row: '{modelseed_reaction_data}'"
            )
        reaction.modelseed_id = modelseed_id

        modelseed_name = modelseed_reaction_data['name']
        if pd.isna(modelseed_name):
            reaction.modelseed_name = None
        else:
            reaction.modelseed_name = modelseed_name

        kegg_reaction_ids: str = modelseed_reaction_data['KEGG']
        if pd.isna(kegg_reaction_ids):
            reaction.kegg_aliases = tuple()
        else:
            reaction.kegg_aliases = tuple(kegg_reaction_ids.split('; '))

        ec_numbers: str = modelseed_reaction_data['ec_numbers']
        if pd.isna(ec_numbers):
            reaction.ec_number_aliases = []
        else:
            reaction.ec_number_aliases = ec_numbers.split('|')

        reversibility = modelseed_reaction_data['reversibility']
        if pd.isna(reversibility):
            raise ConfigError(
                "The row for the reaction in the ModelSEED table was expected to have a 'reversibility' value. "
                f"Here is the data in the row: '{modelseed_reaction_data}'"
            )
        if reversibility == '=' or reversibility == '?':
            # Assume that reactions lacking data ('?') are reversible.
            reaction.reversibility = True
        else:
            reaction.reversibility = False

        decimal_reaction_coefficients = []
        split_stoichiometry = stoichiometry.split(';')
        modelseed_compound_ids = []
        compartments = []
        for entry in split_stoichiometry:
            split_entry = entry.split(':')
            decimal_reaction_coefficients.append(split_entry[0])
            modelseed_compound_ids.append(split_entry[1])
            compartments.append(ModelSEEDDatabase.compartment_ids[int(split_entry[2])])
        reaction.compartments = tuple(compartments)
        reaction_coefficients = self._to_lcm_denominator(decimal_reaction_coefficients)
        direction = modelseed_reaction_data['direction']
        if pd.isna(direction):
            raise ConfigError(
                "The row for the reaction in the ModelSEED table was expected to have a 'direction' value. "
                f"Here is the data in the row: '{modelseed_reaction_data}'"
            )
        if (direction == '>' and reversibility == '<') or (direction == '<' and reversibility == '>'):
            # The way the reaction is written is the opposite of the way the reaction proceeds.
            reaction_coefficients = [-c for c in reaction_coefficients]
        reaction.coefficients = tuple(reaction_coefficients)

        return reaction, modelseed_compound_ids

    def _to_lcm_denominator(self, floats: List[float]) -> Tuple[int]:
        """
        Convert a list of numbers to their lowest common integer multiples.

        Parameters
        ==========
        floats : List[float]

        Returns
        =======
        List[int]
        """
        def lcm(a, b):
            return a * b // math.gcd(a, b)
        rationals = [fractions.Fraction(f).limit_denominator() for f in floats]
        lcm_denom = functools.reduce(lcm, [r.denominator for r in rationals])
        return list(int(r.numerator * lcm_denom / r.denominator) for r in rationals)

    def _get_modelseed_compound(self, modelseed_compound_data: Dict) -> ModelSEEDCompound:
        """
        Generate a ModelSEED compound object from its entry in the ModelSEED table.

        Parameters
        ==========
        modelseed_compound_data : Dict
            A dictionary representation of a row for a compound in the ModelSEED compound table set
            up by anvi'o.

        Returns
        =======
        ModelSEEDCompound
            An object representation of the ModelSEED compound.
        """
        compound = ModelSEEDCompound()
        compound.modelseed_id = modelseed_compound_data['id']

        modelseed_name = modelseed_compound_data['name']
        if pd.isna(modelseed_name):
            compound.modelseed_name = None
        else:
            compound.modelseed_name = modelseed_name

        kegg_aliases: str = modelseed_compound_data['KEGG']
        if pd.isna(kegg_aliases):
            compound.kegg_aliases = tuple()
        else:
            compound.kegg_aliases = tuple(kegg_aliases.split('; '))

        formula = modelseed_compound_data['formula']
        if pd.isna(formula):
            compound.formula = None
            # compounds without formulas have a nominal charge of 10000000 in compounds.tsv
            compound.charge = None
        else:
            compound.formula = formula
            charge = modelseed_compound_data['charge']
            if pd.isna(charge):
                raise ConfigError(
                    f"The charge of a ModelSEED compound, '{compound.modelseed_id}', was not recorded "
                    "in 'compounds.tsv' but is expected to be present as an integer. Here is the data "
                    f"in the row for the compound: '{modelseed_compound_data}'"
                )
            compound.charge = charge

        return compound

    def _store_contigs_database_reactions(self, network: GenomicNetwork, contigs_db: str) -> None:
        """
        Store reaction data in the relevant contigs database table.

        Parameters
        ==========
        network : GenomicNetwork
            The reaction network generated from gene KO annotations in the contigs database.

        contigs_db : str
            The path to the contigs database from which the reaction network was generated.
        """
        # Transfer data from reaction objects to dictionaries mapping to table entries.
        reactions_data: Dict[str, Dict] = {}
        for modelseed_reaction_id, reaction in network.reactions.items():
            reaction_data = {}
            reaction_data['modelseed_reaction_id'] = modelseed_reaction_id
            reaction_data['modelseed_reaction_name'] = reaction.modelseed_name
            reaction_data['metabolite_modelseed_ids'] = ', '.join([c.modelseed_id for c in reaction.compounds])
            reaction_data['stoichiometry'] = ', '.join([str(c) for c in reaction.coefficients])
            reaction_data['compartments'] = ', '.join(reaction.compartments)
            reaction_data['reversibility'] = reaction.reversibility
            # Record KEGG REACTION IDs and EC numbers that are aliases of ModelSEED reactions but
            # are *NOT* associated with gene KO annotations; associated aliases are recorded later.
            reaction_data['other_kegg_reaction_ids'] = ', '.join(
                set(reaction.kegg_aliases).difference(set(network.modelseed_kegg_aliases[modelseed_reaction_id]))
            )
            reaction_data['other_ec_numbers'] = ', '.join(
                set(reaction.ec_number_aliases).difference(set(network.modelseed_ec_number_aliases[modelseed_reaction_id]))
            )
            reactions_data[modelseed_reaction_id] = reaction_data

        # Get *KO* KEGG REACTION ID and EC number aliases of each ModelSEED reaction. These are not
        # all possible aliases, but only those associated with KOs that matched genes. Structure
        # alias data as follows:
        # <ModelSEED reaction ID>: {
        #   <KEGG REACTION ID 1>: [<KO IDs associated with KEGG REACTION ID 1>],
        #   <KEGG REACTION ID 2>: [<KO IDs associated with KEGG REACTION ID 2>],
        #   ...
        # }
        # <ModelSEED reaction ID>: {
        #   <EC number 1>: [<KO IDs associated with EC number 1>],
        #   <EC number 2>: [<KO IDs associated with EC number 2>],
        # ...
        # }
        ko_reaction_aliases: Dict[str, Tuple[Dict[str, List[str]], Dict[str, List[str]]]] = {
            modelseed_reaction_id: ({}, {}) for modelseed_reaction_id in reactions_data
        }
        for ko_id, ko in network.kos.items():
            for modelseed_reaction_id, reaction in ko.reactions.items():
                aliases = ko_reaction_aliases[modelseed_reaction_id]

                kegg_reaction_aliases = aliases[0]
                kegg_reaction_ids = ko.kegg_reaction_aliases[modelseed_reaction_id]
                for kegg_reaction_id in kegg_reaction_ids:
                    try:
                        ko_ids: List = kegg_reaction_aliases[kegg_reaction_id]
                    except KeyError:
                        kegg_reaction_aliases[kegg_reaction_id] = ko_ids = []
                    ko_ids.append(ko_id)

                ec_number_aliases = aliases[1]
                ec_numbers = ko.ec_number_aliases[modelseed_reaction_id]
                for ec_number in ec_numbers:
                    try:
                        ko_ids: List = ec_number_aliases[ec_number]
                    except KeyError:
                        ec_number_aliases[ec_number] = ko_ids = []
                    ko_ids.append(ko_id)
        for modelseed_reaction_id, aliases in ko_reaction_aliases.items():
            reaction_data = reactions_data[modelseed_reaction_id]

            # Make the entry for KO KEGG REACTION aliases, which looks akin to the following arbitrary example:
            # 'R00001: (K00010, K00100); R01234: (K54321)'
            kegg_reaction_aliases = aliases[0]
            entry = []
            for kegg_reaction_id, ko_ids in kegg_reaction_aliases.items():
                entry.append(f'{kegg_reaction_id}: ({", ".join(sorted(ko_ids))})')
            reaction_data['ko_kegg_reaction_source'] = '; '.join(sorted(entry))

            # Make the entry for KO EC number aliases, which looks akin to the following arbitrary example:
            # '1.1.1.1: (K00010, K00100); 1.2.3.4: (K12345); 6.7.8.99: (K65432)
            ec_number_aliases = aliases[1]
            entry = []
            for ec_number, ko_ids in ec_number_aliases.items():
                entry.append(f'{ec_number}: ({", ".join(sorted(ko_ids))})')
            reaction_data['ko_ec_number_source'] = '; '.join(sorted(entry))

        reactions_table = pd.DataFrame.from_dict(reactions_data, orient='index').reset_index(drop=True).sort_values('modelseed_reaction_id')
        reactions_table = reactions_table[tables.gene_function_reactions_table_structure]

        cdb = ContigsDatabase(contigs_db)
        cdb.db._exec_many(
            f'''INSERT INTO {tables.gene_function_reactions_table_name} VALUES ({','.join('?' * len(tables.gene_function_reactions_table_structure))})''',
            reactions_table.values
        )
        cdb.disconnect()

    def _store_contigs_database_metabolites(self, network: GenomicNetwork, contigs_db: str) -> None:
        """
        Store metabolite data in the relevant contigs database table.

        Parameters
        ==========
        network : GenomicNetwork
            The reaction network generated from gene KO annotations in the contigs database.

        contigs_db : str
            The path to the contigs database from which the reaction network was generated.

        Returns
        =======
        None
        """
        # Transfer data from metabolite objects to dictionaries mapping to table entries.
        metabolites_data = {}
        for modelseed_compound_id, compound in network.metabolites.items():
            metabolite_data = {}
            metabolite_data['modelseed_compound_id'] = modelseed_compound_id
            metabolite_data['modelseed_compound_name'] = compound.modelseed_name
            metabolite_data['kegg_aliases'] = ', '.join(compound.kegg_aliases)
            metabolite_data['formula'] = compound.formula
            metabolite_data['charge'] = compound.charge
            metabolites_data[modelseed_compound_id] = metabolite_data

        metabolites_table = pd.DataFrame.from_dict(metabolites_data, orient='index').reset_index(drop=True).sort_values('modelseed_compound_id')
        metabolites_table = metabolites_table[tables.gene_function_metabolites_table_structure]

        cdb = ContigsDatabase(contigs_db)
        cdb.db._exec_many(
            f'''INSERT INTO {tables.gene_function_metabolites_table_name} VALUES ({','.join('?' * len(tables.gene_function_metabolites_table_structure))})''',
            metabolites_table.values
        )
        cdb.disconnect()

    def hash_ko_annotations(self, gene_function_calls_dict: Dict) -> str:
        """
        Hash gene KO annotations in a contigs database to concisely represent the data used to
        construct a reaction network.

        Parameters
        ==========
        gene_function_calls_dict : str
            Dictionary containing gene KO annotations loaded by a contigs superclass.

        Returns
        =======
        str
            Hash representation of all gene KO annotations.
        """
        ko_annotations = []
        for gcid, gene_dict in gene_function_calls_dict.items():
            ko_data = gene_dict['KOfam']
            ko_id = ko_data[0]
            ko_name = ko_data[1]
            e_value = ko_data[2]
            ko_annotations.append((str(gcid), ko_id, ko_name, str(e_value)))
        # Sort KO annotation data in ascending order of gene caller ID and KO accession.
        ko_annotations = sorted(ko_annotations, key=lambda x: (x[0], x[1]))

        ko_annotations_string = ''
        for ko_annotation in ko_annotations:
            ko_annotations_string += ''.join(ko_annotation)

        hashed_ko_annotations = hashlib.sha1(ko_annotations_string.encode('utf-8')).hexdigest()
        return hashed_ko_annotations

    # def make_pangenomic_network(
    #     self,
    #     genomes_storage_db: str,
    #     pan_db: str,
    #     store: bool = True,
    #     overwrite_existing_network: bool = False
    # ) -> PangenomicNetwork:
    #     """
    #     Make a pangenomic metabolic reaction network from KEGG Orthologs stored a genomes storage
    #     database and gene clusters stored in a pan database.

    #     Parameters
    #     ==========
    #     genomes_storage_db : str
    #         Path to a genomes storage database. The pangenomic network is derived from gene KO
    #         annotations stored in the database.

    #     pan_db : str
    #         Path to a pan database. The pangenomic network is determined for gene clusters stored in
    #         the database.

    #     Returns
    #     =======
    #     PangenomicNetwork
    #         The network derived from the pangenomic databases.
    #     """
    #     return
