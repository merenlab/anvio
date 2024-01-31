# -*- coding: utf-8
# pylint: disable=line-too-long
"""Generate, manipulate, and export metabolic reaction networks from gene annotations."""

from __future__ import annotations

import os
import re
import glob
import json
import math
import time
import random
import shutil
import hashlib
import tarfile
import zipfile
import argparse
import tempfile
import fractions
import functools
import collections
import numpy as np
import pandas as pd
import multiprocessing as mp

from copy import deepcopy
from argparse import Namespace
from dataclasses import dataclass, field
from typing import Any, Dict, List, Set, Tuple, Union, Iterable

import anvio.utils as utils
import anvio.dbinfo as dbinfo
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.db import DB
from anvio.errors import ConfigError
from anvio import DEBUG, __file__ as ANVIO_PATH, __version__ as VERSION
from anvio.dbops import ContigsDatabase, PanDatabase, PanSuperclass, ProfileDatabase


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


run_quiet = terminal.Run(verbose=False)

# Network statistics are stored in a dictionary of dictionaries. Keys in the outer dictionary are
# "classes" of network statistics. Keys in the inner dictionary are statistics themselves.
GenomicNetworkStats = Dict[str, Dict[str, Any]]
PangenomicNetworkStats = Dict[str, Dict[str, Any]]

RANDOM_SEED = 1066


@dataclass
class ModelSEEDCompound:
    """
    Representation of a chemical (a compound, element, or ions thereof) or a class of chemicals
    (either abstract, like 'Cofactors' and 'Biomass', or defined, like 'Carboxylic acid' and
    'Polynucleotides'), with properties given by the ModelSEED Biochemistry database.

    Attributes
    ==========
    modelseed_id : str, None
        The ModelSEED compound ID, formatted 'cpdXXXXX', where each X is a digit, e.g., 'cpd00001'.

    modelseed_name : str, None
        Name of the ModelSEED compound, e.g., 'cpd00001' has the name, 'H2O'. When absent in the
        database, assumes a value of None.

    kegg_aliases : Tuple[str], None
        The KEGG COMPOUND IDs that are known to possibly alias the ModelSEED compound, according to
        the ModelSEED database, e.g., 'cpd00001' has the aliases, ('C00001', 'C01328'). A KEGG
        COMPOUND ID is formatted 'CXXXXX', where each X is a digit, e.g., 'C00001'.

    charge : int, None
        The electrical charge of the ModelSEED compound, e.g., 'cpd00001' has charge 0. ModelSEED
        compounds without a formula have a nominal charge of 10000000 in the database.

    formula : str, None
        The formula of the ModelSEED compound, e.g., 'cpd00001' has the formula, 'H2O'. When absent
        in the database, assumes a value of None.

    abundances : Dict[str, float], dict()
        Abundance profile data (from metabolomics, for instance) with each key being a sample name
        and each value being the abundance of the ModelSEED compound in that sample.
    """
    modelseed_id: str = None
    modelseed_name: str = None
    kegg_aliases: Tuple[str] = None
    charge: int = None
    formula: str = None
    abundances: Dict[str, float] = field(default_factory=dict)

@dataclass
class ModelSEEDReaction:
    """
    Representation of a reaction, with properties given by the ModelSEED Biochemistry database.

    Attributes
    ==========
    modelseed_id : str, None
        The ModelSEED reaction ID, formatted 'rxnXXXXX', where each X is a digit, e.g.,
        'rxn00001'.

    modelseed_name : str, None
        Name of the reaction, e.g., 'rxn00001' has the name, 'diphosphate phosphohydrolase'. When
        absent in the database, assumes a value of None.

    kegg_aliases : Tuple[str], None
        The KEGG REACTION IDs that are known to possibly alias the ModelSEED reaction, according to
        the ModelSEED database, e.g., 'rxn00001' has the aliases, ('R00004'). A KEGG REACTION ID is
        formatted 'RXXXXX', where each X is a digit, e.g., 'R00001'.

    ec_number_aliases : Tuple[str], None
        The EC numbers that are known to possibly alias the ModelSEED reaction, according to the
        ModelSEED database, e.g., 'rxn00001' has the aliases, ('3.6.1.1').

    compounds : Tuple[ModelSEEDCompound], None
        ModelSEED compound IDs of reactants and products involved in the reaction, e.g., 'rxn00001'
        involves the compounds, ('cpd00001', 'cpd00012', 'cpd00009', 'cpd00067'). A compound ID is
        formatted 'cpdXXXXX', where each X is a digit, e.g., 'cpd00001'. Each compound item has a
        corresponding stoichiometric reaction coefficient in the attribute, 'coefficients', and a
        corresponding cellular compartment in the attribute, 'compartments'.

    coefficients : Tuple[int], None
        Integer stoichiometric reaction coefficients of reactants and products, with negative
        coefficients indicating reactants and positive coefficients indicating products, e.g.,
        'rxn00001' has the coefficients, (-1, -1, 2, 1). Each coefficient item has a corresponding
        ModelSEED compound ID in the attribute, 'compounds', and a corresponding cellular
        compartment in the attribute, 'compartments'.

    compartments : Tuple[str], None
        Cellular compartments of reactants and products, with valid values being 'c' for 'cytosolic'
        and 'e' for 'extracellular', e.g., 'rxn00001' involves the compartments, ('c', 'c', 'c',
        'c'). Each compartment item has a corresponding ModelSEED compound ID in the attribute,
        'compounds', and a corresponding stoichiometric reaction coefficient in the attribute,
        'coefficients'.

    reversibility : bool, None
        Reaction reversibility, with True indicating the reaction is reversible and False indicating
        the reaction is irreversible given the equation encoded in the attributes, 'compounds',
        'coefficients', and 'compartments'. For example, 'rxn00001' has a value of False.
    """
    modelseed_id: str = None
    modelseed_name: str = None
    kegg_aliases: Tuple[str] = None
    ec_number_aliases: Tuple[str] = None
    compounds: Tuple[ModelSEEDCompound] = None
    coefficients: Tuple[int] = None
    compartments: Tuple[str] = None
    reversibility: bool = None

@dataclass
class KO:
    """
    Representation of a KEGG Ortholog (KO).

    Attributes
    ==========
    id : str, None
        KEGG ORTHOLOGY ID in the format, 'KXXXXX', where X is a digit, e.g., 'K00001'.

    name : str, None
        Name of the KO, e.g., 'K00001' has the name, 'alcohol dehydrogenase [EC:1.1.1.1]'.

    modules : Dict[str, KEGGModule], dict
        KEGG modules containing the KO, with keys being module IDs.

    hierarchies : Dict[str, Dict[Tuple[str], Tuple[BRITECategory]]], dict
        Membership of the KO in BRITE hierarchies. Keys are hierarchy IDs. Values are dictionary
        representations of categorizations in the hierarchy. For example, 'K00844', hexokinase, is
        classified multiple ways in the 'KEGG Orthology (KO)' hierarchy, 'ko00001', including '09100
        Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis
        [PATH:ko00010]' and '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00051 Fructose
        and mannose metabolism [PATH:ko00051]'. This hierarchy and these classifications would be
        represented as follows: {'ko00001': {('09100 Metabolism', '09101 Carbohydrate metabolism',
        '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'): (<BRITECategory for '09100 ...'>,
        <BRITECategory for '09101 ...'>, <BRITECategory for '00010 ...'>), ('09100 Metabolism',
        '09101 Carbohydrate metabolism', '00051 Fructose and mannose metabolism [PATH:ko00051]'):
        (<BRITECategory for '09100 ...'>, <BRITECategory for '09101 ...'>, <BRITECategory for '00051
        ...'>)}}

    pathways : Dict[str, KEGGPathway], dict
        KEGG pathways containing the KO, with keys being pathway IDs.

    reactions : Dict[str, ModelSEEDReaction], dict()
        ModelSEED reactions associated with the KO via KO KEGG reaction and EC number annotations.
        Keys are ModelSEED reaction IDs and values are 'ModelSEEDReaction' objects. A ModelSEED
        reaction ID is formatted 'rxnXXXXX', where each X is a digit, e.g., 'rxn00001'.

    kegg_reaction_aliases : Dict[str, List[str]], dict()
        KEGG reaction annotations of the KO that alias ModelSEED reactions. A KEGG REACTION ID is
        formatted 'RXXXXX', where each X is a digit, e.g., 'R00001'. For example, KO 'K00003' has
        two KEGG reaction annotations, both of which are associated with ModelSEED reactions via the
        ModelSEED database: {'R01773': ['rxn01301', 'rxn27933'], 'R01775': ['rxn01302',
        'rxn27932']}. Note that a ModelSEED reaction may have more KEGG reaction aliases than those
        annotating the KO: all known KEGG reaction aliases of the ModelSEED reaction in the
        ModelSEED database are recorded in the 'kegg_aliases' attribute of a 'ModelSEEDReaction'
        object.

    ec_number_aliases : Dict[str, List[str]], dict()
        EC number annotations of the KO that alias ModelSEED reactions. For example, KO 'K00003' has
        one EC number annotation, which is associated with ModelSEED reactions via the ModelSEED
        database: {'1.1.1.3': ['rxn01301', 'rxn01302', 'rxn19904', 'rxn27931', 'rxn27932',
        'rxn27933', 'rxn33957']}. Note that a ModelSEED reaction may have more EC number aliases
        than those annotating the KO: all known EC number aliases of the ModelSEED reaction in the
        ModelSEED database are recorded in the 'ec_number_aliases' attribute of a
        'ModelSEEDReaction' object.
    """
    id: str = None
    name: str = None
    modules: Dict[str, KEGGModule] = field(default_factory=dict)
    hierarchies: Dict[str, Dict[Tuple[str], Tuple[BRITECategory]]] = field(default_factory=dict)
    pathways: Dict[str, KEGGPathway] = field(default_factory=dict)
    reactions: Dict[str, ModelSEEDReaction] = field(default_factory=dict)
    kegg_reaction_aliases: Dict[str, List[str]] = field(default_factory=dict)
    ec_number_aliases: Dict[str, List[str]] = field(default_factory=dict)

@dataclass
class KEGGModule:
    """
    Representation of a KEGG module with KOs in the reaction network.

    Attributes
    ==========
    id : str, None
        KEGG MODULE ID in the format, 'MXXXXX', where X is a digit, e.g., 'M00001'.

    name : str, None
        Name of the module, e.g., 'M00001' has the name, 'Glycolysis (Embden-Meyerhof pathway),
        glucose => pyruvate'.

    kos : Dict[str, KO], dict()
        Reaction network KOs that are in the module, with keys being KO IDs.

    pathways : Dict[str, KEGGPathway], dict()
        KEGG pathways containing the module, with keys being KEGG pathway IDs.
    """
    id: str = None
    name: str = None
    kos: Dict[str, KO] = field(default_factory=dict)
    pathways: Dict[str, KEGGPathway] = field(default_factory=dict)

@dataclass
class KEGGPathway:
    """
    Representation of a KEGG pathway with KOs in the reaction network.

    Attributes
    ==========
    id : str, None
        The KEGG PATHWAY ID in the format, 'XXXXX', where X is a digit, e.g., '00010' represents
        'Glycolysis / Gluconeogensis', and corresponds to the reference pathway map, 'map00010'.

    name : str, None
        Name of the pathway, e.g., 'Glycolysis / Gluconeogenesis' for pathway ID '00010'.

    category : BRITECategory, None
        Certain pathways are equivalent to categories in KEGG BRITE hierarchies, such as bottom-most
        categories in the hierarchy, 'ko00001', e.g., '00010 Glycolysis / Gluconeogenesis
        [PATH:ko00010]'.

    kos : Dict[str, KO], dict
        Reaction network KOs in the pathway. To reiterate, this does not include KOs in the pathway
        that are not in the reaction network. Keys are KO IDs.

    modules : Dict[str, KEGGModule], dict
        Modules in the pathway that contain reaction network KOs. To reiterate, this does not
        include modules in the pathway that are not in the reaction network. Keys are module IDs.
    """
    id: str = None
    name: str = None
    category: BRITECategory = None
    kos: Dict[str, KO] = field(default_factory=dict)
    modules: Dict[str, KEGGModule] = field(default_factory=dict)

@dataclass
class BRITEHierarchy:
    """
    Representation of a KEGG BRITE hierarchy with KOs in the reaction network.

    Attributes
    ==========
    id : str, None
        BRITE hierarchy ID in the format, 'koXXXXX' or 'brXXXXX', where X is a digit, e.g.,
        'ko00001'.

    name : str, None
        Name of the hierarchy, e.g., 'ko00001' has the name, 'KEGG Orthology (KO)'.

    categories : Dict[Tuple[str], Tuple[BRITECategory]], dict
        Categorizations of reaction network KOs in the hierarchy. Categories at each level receive
        their own entries. For example, 'K00844', hexokinase, is classified multiple ways in the
        'KEGG Orthology (KO)' hierarchy, 'ko00001', including '09100 Metabolism >>> 09101
        Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis [PATH:00010]' and '09100
        Metabolism >>> 09101 Carbohydrate metabolism >>> 00051 Fructose and mannose metabolism
        [PATH:00051]'. These categorizations would yield entries like the following: {('09100
        Metabolism', ): (<BRITECategory for '09100 ...'>, ), ('09100 Metabolism', '09101
        Carbohydrate metabolism'): (<BRITECategory for '09100 ...'>, <BRITECategory for '09101
        ...'>), ('09100 Metabolism', '09101 Carbohydrate metabolism', '00010 Glycolysis /
        Gluconeogenesis [PATH:00010]'): (<BRITECategory for '09100 ...'>, <BRITECategory for '09101
        ...'>, <BRITECategory for '00010 ...'>), ('09100 Metabolism', '09101 Carbohydrate
        metabolism', '00051 Fructose and mannose metabolism [PATH:00051]'): (<BRITECategory for
        '09100 ...'>, <BRITECategory for '09101 ...'>, <BRITECategory for '00051 ...'>)}

    kos : Dict[str, KO], dict
        Reaction network KOs in the hierarchy. To reiterate, this does not include KOs in the
        hierarchy that are not in the reaction network. Keys are KO IDs.
    """
    id: str = None
    name: str = None
    categories: Dict[Tuple[str], Tuple[BRITECategory]] = field(default_factory=dict)
    kos: Dict[str, KO] = field(default_factory=dict)

@dataclass
class BRITECategory:
    """
    Representation of a KEGG BRITE hierarchy category with KOs in the reaction network.

    Attributes
    ==========
    id : str
        Unique ID for the category comprising the hierarchy ID and the hierarchical categorization.
        The following example demonstrates ID format. In the 'KEGG Orthology (KO)' hierarchy,
        'ko00001', there is a category, '09100 Metabolism >>> 09101 Carbohydrate metabolism >>>
        00010 Glycolysis / Gluconeogenesis [PATH:00010]'. This yields the ID, 'ko00001: 09100
        Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis
        [PATH:00010]'.

    name : str
        Name of the category. These need not be unique in a hierarchy. For example, there are
        multiple categories called 'Small subunit' and 'Large subunit' in the 'Ribosome' hierarchy.

    hierarchy : BRITEHierarchy, None
        The BRITE hierarchy containing the category.

    supercategory : BRITECategory, None
        The encompassing category. This is None if there is no category higher in the hierarchy. For
        example, the supercategory immediately above the category, '00010 Glycolysis /
        Gluconeogenesis' in the hierarchy, 'ko00001', is '09101 Carbohydrate metabolism'.

    subcategories : List[BRITECategory], list
        The encompassed categories containing KOs in the reaction network. This is an empty list if
        there are no categories lower in the hierarchy. For example, the category, 'Polyketide
        synthase (PKS) >>> Modular type I PKS' in the hierarchy, 'ko01008' encompasses the
        categories, 'cis-AT PKS' and 'trans-AT PKS'.

    pathway : KEGGPathway, None
        Certain categories are equivalent to KEGG pathways, such as bottom-most categories in the
        hierarchy, 'ko00001', e.g., '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'.

    kos : Dict[str, KO], dict
        Reaction network KOs in the category (and all subcategories). To reiterate, this does not
        include KOs in the category that are not in the reaction network. Keys are KO IDs.
    """
    id: str = None
    name: str = None
    hierarchy: BRITEHierarchy = None
    supercategory: BRITECategory = None
    subcategories: List[BRITECategory] = field(default_factory=list)
    pathway: KEGGPathway = None
    kos: Dict[str, KO] = field(default_factory=dict)

@dataclass
class Gene:
    """
    Representation of a gene.

    Attributes
    ==========
    gcid : int, None
        The gene callers ID, or unique anvi'o identifier, of the gene: a non-negative integer.

    kos : Dict[str, KO], dict()
        KEGG Orthologs (KOs) annotating the gene. Keys are KO IDs, which are formatted as 'KXXXXX',
        where each X is a digit, e.g., 'K00001'. Values are 'KO' objects.

    e_values : Dict[str, float], dict()
        E-values express the strength of KO-gene associations. Keys are KO IDs; values are
        non-negative numbers.

    protein : Protein, None
        This object is used for storing abundance data on the protein expressed by the gene (from
        proteomics, for instance).
    """
    gcid: int = None
    kos: Dict[str, KO] = field(default_factory=dict)
    e_values: Dict[str, float] = field(default_factory=dict)
    protein: Protein = None

@dataclass
class Protein:
    """
    This object stores protein abundance data (from proteomics, for instance).

    Attributes
    ==========
    id : int, None
        The unique anvi'o ID for the protein: a non-negative integer.

    genes : Dict[int, Gene], dict()
        Genes that can express the protein. Keys are gene callers IDs; values are 'Gene' objects.

    abundances : Dict[str, float], dict()
        Protein abundance profile data with each key being a sample name and each value being the
        abundance of the protein expressed by the gene in that sample.
    """
    id: int = None
    genes: Dict[int, Gene] = field(default_factory=dict)
    abundances: Dict[str, float] = field(default_factory=dict)

@dataclass
class GeneCluster:
    """
    Representation of a gene cluster.

    Attributes
    ==========
    gene_cluster_id : int, None
        The unique anvi'o ID for the gene cluster: a non-negative integer.

    genomes : List[str], []
        The names of the genomes contributing the genes in the cluster.

    ko : KO, None
        The consensus KO among the genes in the cluster. (Consensus KOs can be found from a
        pangenome by the anvi'o method, 'dbops.PanSuperclass.get_gene_cluster_function_summary'.)
        Note that the individual gene KO annotations underlying the consensus annotation are not
        tracked.
    """
    gene_cluster_id: int = None
    genomes: List[str] = field(default_factory=list)
    ko: KO = None

class ReactionNetwork:
    """
    A reaction network predicted from KEGG KO and ModelSEED annotations.

    A reaction network need not be fully connected: it is not guaranteed that there exists a path
    through the network from one arbitrary reaction to another.

    Attributes
    ==========
    kos : Dict[str, KO], dict()
        KOs in the network, with keys being KO IDs.

    modules : Dict[str, KEGGModule], dict()
        KEGG modules containing KOs in the network, with keys being module IDs.

    pathways : Dict[str, KEGGPathway], dict()
        KEGG pathways containing KOs in the network, with keys being pathway IDs.

    hierarchies : Dict[str, BRITEHierarchy], dict()
        KEGG BRITE hierarchies containing KOs in the network, with keys being hierarchy IDs.

    categories : Dict[str, Dict[Tuple[str], Tuple[BRITECategory]]], dict()
        KEGG BRITE hierarchy categories containing KOs in the network. Keys are hierarchy IDs.
        Values are dictionary representations of categorizations in the hierarchy. Categories at
        each level receive their own entries. For example, 'K00844', hexokinase, is classified
        multiple ways in the 'KEGG Orthology (KO)' hierarchy, 'ko00001', including '09100
        Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis
        [PATH:00010]' and '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00051 Fructose
        and mannose metabolism [PATH:00051]'. These categorizations would yield entries like the
        following: {'ko00001': {('09100 Metabolism', ): (<BRITECategory for '09100 ...'>, ), ('09100
        Metabolism', '09101 Carbohydrate metabolism'): (<BRITECategory for '09100 ...'>,
        <BRITECategory for '09101 ...'>), ('09100 Metabolism', '09101 Carbohydrate metabolism',
        '00010 Glycolysis / Gluconeogenesis [PATH:00010]'): (<BRITECategory for '09100 ...'>,
        <BRITECategory for '09101 ...'>, <BRITECategory for '00010 ...'>), ('09100 Metabolism',
        '09101 Carbohydrate metabolism', '00051 Fructose and mannose metabolism [PATH:00051]'):
        (<BRITECategory for '09100 ...'>, <BRITECategory for '09101 ...'>, <BRITECategory for '00051
        ...'>)}}

    reactions : Dict[str, ModelSEEDReaction], dict()
        ModelSEED reactions in the network, with keys being reaction IDs.

    metabolites : Dict[str, ModelSEEDCompound], dict()
        ModelSEED compounds in the network, with keys being metabolite IDs.

    kegg_modelseed_aliases : Dict[str, List[str]], dict()
        This maps KEGG REACTION IDs associated with KOs in the network to ModelSEED reactions
        aliased by the KEGG reaction. KO-associated KEGG reactions that do not alias ModelSEED
        reactions are not included.

    ec_number_modelseed_aliases : Dict[str, List[str]], dict()
        This maps EC numbers associated with KOs in the network to ModelSEED reactions aliased by
        the EC number. KO-associated EC numbers that do not alias ModelSEED reactions are not
        included.

    modelseed_kegg_aliases : Dict[str, List[str]], dict()
        This maps the IDs of ModelSEED reactions in the network to lists of KEGG REACTION IDs that
        are associated with KOs in the network and alias the ModelSEED reaction.

    modelseed_ec_number_aliases : Dict[str, List[str]], dict()
        This maps the IDs of ModelSEED reactions in the network to lists of EC numbers that are
        associated with KOs in the network and alias the ModelSEED reaction.

    run : anvio.terminal.Run, anvio.terminal.Run()
        This object prints run information to the terminal. This attribute is assigned the argument
        of the same name upon initialization.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        This object prints transient progress information to the terminal. This attribute is
        assigned the argument of the same name upon initialization.

    verbose : bool, True
        Report more information to the terminal if True.
    """
    def __init__(
        self,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        verbose: bool = True
    ) -> None:
        """
        Parameters
        ==========
        run : anvio.terminal.Run, anvio.terminal.Run()
            This object sets the 'run' attribute, which prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object sets the 'progress' attribute, which prints transient progress information
            to the terminal.

        verbose : bool, True
            This sets the 'verbose' attribute, causing more information to be reported to the
            terminal if True.

        Returns
        =======
        None
        """
        self.kos: Dict[str, KO] = {}
        self.modules: Dict[str, KEGGModule] = {}
        self.pathways: Dict[str, KEGGPathway] = {}
        self.hierarchies: Dict[str, BRITEHierarchy] = {}
        self.categories: Dict[str, Dict[Tuple[str], Tuple[BRITECategory]]] = {}
        self.reactions: Dict[str, ModelSEEDReaction] = {}
        self.metabolites: Dict[str, ModelSEEDCompound] = {}
        # The following dictionaries map reaction aliases in the network: as in, not all known
        # aliases, but only those sourced from KOs and contributing ModelSEEDReaction objects.
        self.kegg_modelseed_aliases: Dict[str, List[str]] = {}
        self.ec_number_modelseed_aliases: Dict[str, List[str]] = {}
        self.modelseed_kegg_aliases: Dict[str, List[str]] = {}
        self.modelseed_ec_number_aliases: Dict[str, List[str]] = {}

        self.run = run
        self.progress = progress
        self.verbose = verbose

    def _copy(self, copied_network: ReactionNetwork) -> None:
        """
        This method is used in the process of copying a reaction network, and contains steps in
        common to different types of network: copy the attributes of the network BESIDES genes (in a
        GenomicNetwork) / gene clusters (in a PangenomicNetwork) and protein abundances (which can
        only be stored in a GenomicNetwork).

        Parameters
        ==========
        copied_network : ReactionNetwork
            The network copy under construction.
        """
        # Copy metabolites.
        for modelseed_id, metabolite in self.metabolites.items():
            copied_metabolite = ModelSEEDCompound()
            copied_metabolite.modelseed_id = modelseed_id
            copied_metabolite.modelseed_name = metabolite.modelseed_name
            copied_metabolite.kegg_aliases = metabolite.kegg_aliases
            copied_metabolite.charge = metabolite.charge
            copied_metabolite.formula = metabolite.formula
            abundances: Dict[str, float] = getattr(metabolite, 'abundances', None)
            if abundances is None:
                # Metabolite sample abundances cannot be recorded in a pangenomic network.
                continue
            copied_metabolite.abundances = abundances.copy()

            copied_network.metabolites[modelseed_id] = copied_metabolite

        # Copy reactions.
        for modelseed_id, reaction in self.reactions.items():
            copied_reaction = ModelSEEDReaction()
            copied_reaction.modelseed_id = modelseed_id
            copied_reaction.modelseed_name = reaction.modelseed_name
            copied_reaction.kegg_aliases = reaction.kegg_aliases
            copied_reaction.ec_number_aliases = reaction.ec_number_aliases
            metabolites = []
            for metabolite in reaction.compounds:
                metabolites.append(copied_network.metabolites[metabolite.modelseed_id])
            copied_reaction.compounds = tuple(metabolites)
            copied_reaction.coefficients = reaction.coefficients
            copied_reaction.compartments = reaction.compartments
            copied_reaction.reversibility = reaction.reversibility

            copied_network.reactions[modelseed_id] = copied_reaction

        # Copy KOs.
        for ko_id, ko in self.kos.items():
            copied_ko = KO()
            copied_ko.id = ko_id
            copied_ko.name = ko.name
            for modelseed_id in ko.reactions:
                copied_ko.reactions[modelseed_id] = copied_network.reactions[modelseed_id]
            for modelseed_id, kegg_ids in ko.kegg_reaction_aliases.items():
                copied_ko.kegg_reaction_aliases[modelseed_id] = kegg_ids.copy()
            for modelseed_id, ec_numbers in ko.ec_number_aliases.items():
                copied_ko.ec_number_aliases[modelseed_id] = ec_numbers.copy()

            copied_network.kos[ko_id] = copied_ko

        # Copy reaction alias mappings.
        for kegg_id, modelseed_ids in self.kegg_modelseed_aliases.items():
            copied_network.kegg_modelseed_aliases[kegg_id] = modelseed_ids.copy()

        for ec_number, modelseed_ids in self.ec_number_modelseed_aliases.items():
            copied_network.ec_number_modelseed_aliases[ec_number] = modelseed_ids.copy()

        for modelseed_id, kegg_ids in self.modelseed_kegg_aliases.items():
            copied_network.modelseed_kegg_aliases[modelseed_id] = kegg_ids.copy()

        for modelseed_id, ec_numbers in self.modelseed_ec_number_aliases.items():
            copied_network.modelseed_ec_number_aliases[modelseed_id] = ec_numbers.copy()

        # Copy KEGG modules.
        for module_id, module in self.modules.items():
            copied_module = KEGGModule()
            copied_module.id = module_id
            copied_module.name = module.name
            for ko_id in module.kos:
                copied_module.kos[ko_id] = copied_ko = copied_network.kos[ko_id]
                copied_ko.modules[module_id] = copied_module

            copied_network.modules[module_id] = copied_module

        # Copy KEGG hierarchies.
        for hierarchy_id, hierarchy in self.hierarchies.items():
            copied_hierarchy = BRITEHierarchy()
            copied_hierarchy.id = hierarchy_id
            copied_hierarchy.name = hierarchy.name
            for ko_id in hierarchy.kos:
                copied_hierarchy.kos[ko_id] = copied_network.kos[ko_id]

            copied_network.hierarchies[hierarchy_id] = copied_hierarchy

        # Copy KEGG hierarchy categories.
        copied_categories: Dict[str, BRITECategory] = {}
        for hierarchy_id, categorizations in self.categories.items():
            hierarchy = self.hierarchies[hierarchy_id]

            for categorization in categorizations.values():
                category = categorization[-1]
                category_id = category.id
                assert category_id not in copied_categories

                copied_category = BRITECategory()
                copied_category.id = category_id
                copied_category.name = category.name
                copied_category.hierarchy = hierarchy

                supercategory = category.supercategory
                if supercategory is None:
                    copied_category.supercategory = None
                else:
                    try:
                        copied_supercategory = copied_categories[supercategory.id]
                    except KeyError:
                        copied_supercategory = None
                    copied_category.supercategory = copied_supercategory
                    if copied_supercategory is None:
                        # Replace the subcategory placeholder ID in the copied supercategory with
                        # the copied category.
                        copied_supercategory[
                            copied_supercategory.subcategories.index(category_id)
                        ] = copied_category

                subcategories = category.subcategories
                if not subcategories:
                    copied_category.subcategories = []
                else:
                    for subcategory in subcategories:
                        try:
                            copied_subcategory = copied_categories[subcategory.id]
                        except KeyError:
                            copied_subcategory = None
                        if copied_subcategory is None:
                            # Use the subcategory ID as a placeholder to be replaced by the copied
                            # subcategory object when it is created.
                            copied_category.subcategories.append(subcategory.id)
                        else:
                            copied_category.subcategories.append(copied_subcategory)
                            # The subcategory had already been copied, but the current category (its
                            # supercategory) had not, so the copied subcategory had a supercategory
                            # attribute of None. This is here replaced by the copied category.
                            copied_subcategory.supercategory = copied_category

                for ko_id in category.kos:
                    copied_category.kos[ko_id] = copied_network.kos[ko_id]

                copied_categories[category_id] = copied_category

        for hierarchy_id, categorizations in self.categories.items():
            copied_network_hierarchy: Dict[Tuple[str], Tuple[BRITECategory]] = {}
            copied_network.categories[hierarchy_id] = copied_network_hierarchy
            for category_key, categorization in categorizations.items():
                copied_network_hierarchy[category_key] = tuple(
                    [copied_categories[category.id] for category in categorization]
                )

        # Copy KEGG pathways. For any BRITE categories equivalent to a pathway, reference the copied
        # pathway from the copied category.
        for pathway_id, pathway in self.pathways.items():
            copied_pathway = KEGGPathway()
            copied_pathway.id = pathway_id
            copied_pathway.name = pathway.name
            for ko_id in pathway.kos:
                copied_pathway.kos[ko_id] = copied_ko = copied_network.kos[ko_id]
                copied_ko.pathways[pathway_id] = copied_pathway
            if pathway.category is not None:
                category = pathway.category
                copied_category = copied_categories[category.id]
                copied_pathway.category = copied_category
                assert copied_category.pathway is None
                copied_category.pathway = copied_pathway

            copied_network.pathways[pathway_id] = copied_pathway

        # Record BRITE hierarchical categorizations of each KO.
        for ko_id, ko in self.kos.items():
            copied_ko = copied_network.kos[ko_id]
            for hierarchy_id, categorizations in ko.hierarchies.items():
                copied_ko.hierarchies[hierarchy_id] = copied_categorizations = {}
                for category_key, categorization in categorizations.items():
                    copied_categorizations[category_key] = tuple(
                        [copied_categories[category.id] for category in categorization]
                    )

        # Record the pathway membership of modules.
        for pathway_id, pathway in self.pathways.items():
            copied_pathway = copied_network.pathways[pathway_id]
            for module_id in pathway.modules:
                copied_pathway.modules[module_id] = copied_network.modules[module_id]

        for module_id, module in self.modules.items():
            copied_module = copied_network.modules[module_id]
            for pathway_id in module.pathways:
                copied_module.pathways[pathway_id] = copied_network.pathways[pathway_id]

        # Record categories of each BRITE hierarchy.
        for hierarchy_id, hierarchy in self.hierarchies.items():
            copied_hierarchy = copied_network.hierarchies[hierarchy_id]
            copied_categorizations = copied_hierarchy.categories
            for category_key, categorization in hierarchy.categories.items():
                copied_categorizations[category_key] = tuple(
                    [copied_categories[category.id] for category in categorization]
                )

    def remove_missing_objective_metabolites(self, objective_dict: Dict) -> None:
        """
        Remove metabolites from a biomass objective dictionary that are not produced or consumed by
        any reactions in the network.

        Parameters
        ==========
        objective_dict : dict
            Biomass objective in COBRApy JSON format, like that returned by the method,
            'JSONStructure.get_e_coli_core_objective'.

        Returns
        =======
        None
        """
        objective_metabolites: Dict = objective_dict['metabolites']
        missing_metabolite_ids = []
        if 'original_metabolite_ids' in objective_dict['notes']:
            # The E. coli objective had metabolite BiGG IDs, which were replaced with KEGG COMPOUND
            # IDs, and the original BiGG IDs were recorded in the 'notes' section of the objective.
            missing_original_metabolite_ids = []
            objective_original_metabolites: Dict = objective_dict['notes'][
                'original_metabolite_ids'
            ]
            for metabolite_id, original_metabolite_id in zip(
                objective_metabolites, objective_original_metabolites
            ):
                if metabolite_id[:-2] not in self.metabolites:
                    # The metabolite (removing localization substring) is not in the network.
                    missing_metabolite_ids.append(metabolite_id)
                    missing_original_metabolite_ids.append(original_metabolite_id)
            for original_metabolite_id in missing_original_metabolite_ids:
                objective_original_metabolites.pop(original_metabolite_id)
        else:
            for metabolite_id in objective_metabolites:
                if metabolite_id[:-2] not in self.metabolites:
                    # The metabolite (removing localization substring) is not in the network.
                    missing_metabolite_ids.append(metabolite_id)
        for metabolite_id in missing_metabolite_ids:
            objective_metabolites.pop(metabolite_id)

        if not self.verbose:
            return

        if 'original_metabolite_ids' in objective_dict['notes']:
            id_string = ""
            for original_id, modelseed_id in zip(
                missing_original_metabolite_ids, missing_metabolite_ids
            ):
                id_string += f"{original_id} ({modelseed_id}), "
            id_string = id_string[:-2]
            self.run.info_single(
                f"""\
                The following metabolites were removed from the biomass objective, with the original
                IDs aliasing the ModelSEED compound IDs in parentheses: {id_string}\
                """
            )
        else:
            self.run.info_single(
                f"""\
                The following metabolites, given by their ModelSEED compound IDs, were removed from
                the biomass objective: {', '.join(missing_metabolite_ids)}\
                """
            )

    def _write_remove_metabolites_without_formula_output(
        self,
        output_path : str,
        removed: Dict[str, List]
    ) -> None:
        """
        Parameters
        ==========
        output_path : str
            Write tab-delimited files of metabolites, reactions, KOs, KEGG modules, KEGG pathways,
            KEGG BRITE hierarchies, and KEGG BRITE hierarchy categories removed from the network to
            file locations based on the provided path. For example, if the argument, 'removed.tsv',
            is provided, then the following files will be written: 'removed-metabolites.tsv',
            'removed-reactions.tsv', 'removed-kos.tsv', 'removed-modules.tsv',
            'removed-pathways.tsv', 'removed-hierarchies.tsv', and 'removed-categories.tsv'.

        removed : Dict[str, List]
            Data removed from the network. The dictionary looks like the following for a genomic
            network. (For a pangenomic network, the last gene entry is replaced by a gene cluster
            entry, 'gene_cluster': [<removed GeneCluster objects>].)
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'gene': [<removed Gene objects>]
            }
        """
        # Record the reactions removed as a consequence of involving formulaless metabolites, and
        # record the formulaless metabolites involved in removed reactions.
        removed_metabolite_ids: List[str] = [metabolite.id for metabolite in removed['metabolite']]
        metabolite_removed_reactions: Dict[str, List[str]] = {}
        reaction_removed_metabolites: Dict[str, List[str]] = {}
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            reaction_removed_metabolites[reaction.modelseed_id] = metabolite_ids = []
            for metabolite in reaction.compounds:
                if metabolite.modelseed_id in removed_metabolite_ids:
                    try:
                        metabolite_removed_reactions[metabolite.modelseed_id].append(
                            reaction.modelseed_id
                        )
                    except KeyError:
                        metabolite_removed_reactions[metabolite.modelseed_id] = [
                            reaction.modelseed_id
                        ]
                    metabolite_ids.append(metabolite.modelseed_id)

        metabolite_table = []
        for metabolite in removed['metabolite']:
            metabolite: ModelSEEDCompound
            row = []
            row.append(metabolite.modelseed_id)
            row.append(metabolite.modelseed_name)
            row.append(metabolite.formula)
            try:
                # The metabolite did not have a formula.
                removed_reaction_ids = metabolite_removed_reactions[metabolite.modelseed_id]
            except KeyError:
                # The metabolite had a formula but was removed as a consequence of all the reactions
                # involving the metabolite being removed due to them containing formulaless
                # metabolites: the metabolite did not cause any reactions to be removed.
                row.append("")
                continue
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus the reaction is recorded multiple times.
            row.append(", ".join(sorted(set(removed_reaction_ids))))

        reaction_table = []
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            row = []
            row.append(reaction.modelseed_id)
            row.append(reaction.modelseed_name)
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus is recorded multiple times.
            row.append(
                ", ".join(set(reaction_removed_metabolites[reaction.modelseed_id]))
            )
            row.append(", ".join([metabolite.modelseed_id for metabolite in reaction.compounds]))
            row.append(get_chemical_equation(reaction))
            reaction_table.append(row)

        ko_table = []
        for ko in removed['ko']:
            ko: KO
            row = []
            row.append(ko.id)
            row.append(ko.name)
            row.append(", ".join(ko.reactions))
            ko_table.append(row)

        module_table = []
        for module in removed['module']:
            module: KEGGModule
            row = []
            row.append(module.id)
            row.append(module.name)
            row.append(", ".join(module.kos))
            module_table.append(row)

        pathway_table = []
        for pathway in removed['pathway']:
            pathway: KEGGPathway
            row = []
            row.append(pathway.id)
            row.append(pathway.name)
            row.append(", ".join(pathway.kos))
            row.append(", ".join(pathway.modules))
            pathway_table.append(row)

        hierarchy_table = []
        for hierarchy in removed['hierarchy']:
            hierarchy: BRITEHierarchy
            row = []
            row.append(hierarchy.id)
            row.append(hierarchy.name)
            row.append(", ".join(hierarchy.kos))
            pathway_table.append(row)

        category_table = []
        for category in removed['category']:
            category: BRITECategory
            row = []
            row.append(category.hierarchy.id)
            row.append(category.hierarchy.name)
            row.append(category.id[len(category.hierarchy.id) + 2:])
            row.append(", ".join(category.kos))
            category_table.append(row)

        path_basename, path_extension = os.path.splitext(output_path)
        metabolite_path = f"{path_basename}-metabolites{path_extension}"
        reaction_path = f"{path_basename}-reactions{path_extension}"
        ko_path = f"{path_basename}-kos{path_extension}"
        module_path = f"{path_basename}-modules{path_extension}"
        pathway_path = f"{path_basename}-pathways{path_extension}"
        hierarchy_path = f"{path_basename}-hierarchies{path_extension}"
        category_path = f"{path_basename}-categories{path_extension}"

        pd.DataFrame(
            metabolite_table,
            columns=[
                "ModelSEED compound ID",
                "ModelSEED compound name",
                "Formula",
                "Removed reaction ModelSEED IDs"
            ]
        ).to_csv(metabolite_path, sep='\t', index=False)

        pd.DataFrame(
            reaction_table,
            columns=[
                "ModelSEED reaction ID",
                "ModelSEED reaction name",
                "Removed ModelSEED compound IDs",
                "Reaction ModelSEED compound IDs",
                "Equation"
            ]
        ).to_csv(reaction_path, sep='\t', index=False)

        pd.DataFrame(
            ko_table,
            columns=[
                "KO ID",
                "KO name",
                "KO ModelSEED reaction IDs"
            ]
        ).to_csv(ko_path, sep='\t', index=False)

        pd.DataFrame(
            module_table,
            columns=[
                "KEGG module ID",
                "KEGG module name",
                "Module KOs"
            ]
        ).to_csv(module_path, sep='\t', index=False)

        pd.DataFrame(
            pathway_table,
            columns=[
                "KEGG pathway ID",
                "KEGG pathway name",
                "Pathway KOs",
                "Pathway modules"
            ]
        ).to_csv(pathway_path, sep='\t', index=False)

        pd.DataFrame(
            hierarchy_table,
            columns=[
                "KEGG BRITE hierarchy ID",
                "KEGG BRITE hierarchy name",
                "Hierarchy KOs"
            ]
        ).to_csv(hierarchy_path, sep='\t', index=False)

        pd.DataFrame(
            category_table,
            columns=[
                "KEGG BRITE hierarchy ID",
                "KEGG BRITE hierarchy name",
                "KEGG BRITE hierarchy categorization",
                "Category KOs"
            ]
        ).to_csv(category_path, sep='\t', index=False)

    def _purge_metabolites(self, metabolites_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given metabolites from the network.

        Reactions involving the metabolite are also purged from the network. KOs that are only
        associated with removed reactions are purged. In genomic networks, genes that are only
        associated with removed KOs are purged. In pangenomic networks, gene clusters assigned
        removed KOs are purged. KEGG modules, pathways, BRITE hierarchies, and BRITE hierarchy
        categories only associated with purged KOs are removed.

        Removal of reactions involving the metabolite can also result in other metabolites being
        being removed from the network, those that exclusively participate in these reactions.

        Parameters
        ==========
        metabolites_to_remove : Iterable[str]
            ModelSEED compound IDs to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            The dictionary examples below are for a genomic network. For a pangenomic network, the
            gene entry is replaced by the gene cluster entry, 'gene_cluster': [<removed GeneCluster
            objects>] or 'gene_cluster': []. The examples show protein entries as if the genomic
            network has been annotated with protein abundances; these are absent for genomic
            networks lacking protein annotations and for pangenomic networks.

            If this method is NOT called from the method, '_purge_reactions', then the dictionary
            will look like the following.
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>]
            }

            If this method is called from the method, '_purge_reactions', then the dictionary will
            look like the following.
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene': [],
                'protein': []
            }

            If no metabolites are removed from the network, then the dictionary will look like the
            following regardless of calling method.
            {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene': [],
                'protein': []
            }
        """
        metabolites_to_remove = set(metabolites_to_remove)
        removed_metabolites: List[ModelSEEDCompound] = []
        for modelseed_compound_id in metabolites_to_remove:
            try:
                removed_metabolites.append(self.metabolites.pop(modelseed_compound_id))
            except KeyError:
                # This can occur for two reasons. First, the metabolite from 'metabolites_to_remove'
                # could not be in the network.

                # Second, this can occur when removing other "unintended" metabolites from the
                # network. '_purge_metabolites' was first called with metabolites of interest, then
                # '_purge_reactions' was called from within the method the remove reactions
                # involving the metabolites of interest, and then '_purge_metabolites' was called
                # again from within '_purge_reactions' to remove other metabolites exclusively found
                # in the removed reactions. In this last call of '_purge_metabolites', the
                # 'metabolites_to_remove' also include the metabolites of interest that were already
                # removed from 'self.metabolites' in the original '_purge_metabolites' call. This
                # KeyError occurs when trying to remove those already-removed metabolites.
                pass
        if not removed_metabolites:
            removed = {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': []
            }
            if isinstance(self, GenomicNetwork):
                removed['gene'] = []
                if self.proteins:
                    removed['protein'] = []
            elif isinstance(self, PangenomicNetwork):
                removed['gene_cluster'] = []
            else:
                raise AssertionError
            return removed

        # Purge reactions from the record that involve removed metabolites.
        reactions_to_remove = []
        for modelseed_reaction_id, reaction in self.reactions.items():
            for compound in reaction.compounds:
                if compound.modelseed_id in metabolites_to_remove:
                    reactions_to_remove.append(modelseed_reaction_id)
                    break

        removed = {'metabolite': removed_metabolites}
        if reactions_to_remove:
            removed_cascading_up = self._purge_reactions(reactions_to_remove)
            # There may be other metabolites exclusively involved in the removed reactions; these
            # metabolites were therefore also removed.
            removed['metabolite'] = removed_metabolites + removed_cascading_up.pop('metabolite')
        else:
            # This method must have been called from the method, '_purge_reactions', because the
            # reactions containing the metabolites were already removed from the network.
            removed_cascading_up = {
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': []
            }
            if isinstance(self, GenomicNetwork):
                removed_cascading_up['gene'] = []
                if self.proteins:
                    removed_cascading_up['protein'] = []
            elif isinstance(self, PangenomicNetwork):
                removed_cascading_up['gene_cluster'] = []
            else:
                raise AssertionError
        removed.update(removed_cascading_up)
        return removed

    def _purge_reactions(self, reactions_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given reactions from the network.

        Metabolites that exclusively participate in removed reactions are purged. KOs that are only
        associated with removed reactions are purged. In genomic networks, genes that are only
        associated with removed KOs are purged. In pangenomic networks, gene clusters assigned
        removed KOs are purged. KEGG modules, pathways, BRITE hierarchies, and BRITE hierarchy
        categories only associated with purged KOs are removed.

        Parameters
        ==========
        reactions_to_remove : Iterable[str]
            ModelSEED reaction IDs to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            The dictionary examples below are for a genomic network. For a pangenomic network, the
            gene entry is replaced by the gene cluster entry, 'gene_cluster': [<removed GeneCluster
            objects>] or 'gene_cluster': []. The examples show protein entries as if the genomic
            network has been annotated with protein abundances; these are absent for genomic
            networks lacking protein annotations and for pangenomic networks.

            If this method is NOT called from the method, '_purge_metabolites', or the method,
            '_purge_kos', then the dictionary will look like the following.
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>]
            }

            If this method is called from the method, '_purge_metabolites', then the dictionary will
            look like the following.
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>]
            }

            If this method is called from the method, '_purge_kos', then the dictionary will look
            like the following.
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene': [],
                'protein': []
            }

            If no reactions are removed from the network, then the dictionary will look like the
            following regardless of calling method.
            {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene': [],
                'protein': []
            }
        """
        reactions_to_remove = set(reactions_to_remove)
        removed_reactions: List[ModelSEEDReaction] = []
        for modelseed_reaction_id in reactions_to_remove:
            try:
                removed_reactions.append(self.reactions.pop(modelseed_reaction_id))
            except KeyError:
                # This occurs when the original method called is '_purge_reactions', followed by
                # '_purge_kos', followed by this method again -- 'removed_reactions' will be empty.
                # Alternatively, this occurs if the reaction in 'reactions_to_remove' is not in the
                # network.
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
            removed = {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': []
            }
            if isinstance(self, GenomicNetwork):
                removed['gene'] = []
                if self.proteins:
                    removed['protein'] = []
            elif isinstance(self, PangenomicNetwork):
                removed['gene_cluster'] = []
            else:
                raise AssertionError
            return removed

        # Remove KEGG reaction aliases of ModelSEED reactions from the network.
        kegg_reactions_to_remove = []
        for kegg_reaction_id, modelseed_reaction_ids in self.kegg_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the KEGG reaction were removed, so remove the
                # KEGG reaction as well.
                kegg_reactions_to_remove.append(kegg_reaction_id)
                continue
            for idx in sorted(aliases_to_remove, reverse=True):
                # Only some of the ModelSEED reactions aliased by the KO KEGG reaction were removed.
                modelseed_reaction_ids.pop(idx)
        removed_kegg_reactions = []
        for kegg_reaction_id in kegg_reactions_to_remove:
            removed_kegg_reactions.append(self.kegg_modelseed_aliases.pop(kegg_reaction_id))

        # Remove EC number aliases of ModelSEED reactions from the network.
        ec_numbers_to_remove = []
        for ec_number, modelseed_reaction_ids in self.ec_number_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the EC number were removed, so remove the EC
                # number as well.
                ec_numbers_to_remove.append(ec_number)
                continue
            for idx in sorted(aliases_to_remove, reverse=True):
                # Only some of the ModelSEED reactions aliased by the EC number were removed.
                modelseed_reaction_ids.pop(idx)
        removed_ec_numbers = []
        for ec_number in ec_numbers_to_remove:
            removed_ec_numbers.append(self.ec_number_modelseed_aliases.pop(ec_number))

        # Purge metabolites from the network that are exclusive to removed reactions.
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
                        # Do not remove the metabolite, because it participates in a retained
                        # reaction.
                        metabolites_to_spare.append(idx)
            for idx in sorted(metabolites_to_spare, reverse=True):
                metabolites_to_remove.pop(idx)
        if metabolites_to_remove:
            removed_cascading_down = self._purge_metabolites(metabolites_to_remove)
            for key in [
                'reaction',
                'kegg_reaction',
                'ec_number',
                'ko',
                'module',
                'pathway',
                'hierarchy',
                'category'
            ]:
                removed_cascading_down.pop(key)
            if isinstance(self, GenomicNetwork):
                removed_cascading_down.pop('gene')
                if self.proteins:
                    removed_cascading_down.pop('protein')
            elif isinstance(self, PangenomicNetwork):
                removed_cascading_down.pop('gene_cluster')
            else:
                raise AssertionError
        else:
            # No metabolites were exclusive to the removed reactions. (This point cannot be reached
            # if this method were called from within '_purge_metabolites'.)
            removed_cascading_down = {'metabolite': []}

        # Purge KOs from the network that are only associated with removed reactions.
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
        # If this method was called from '_purge_kos' then the KOs that are only associated with
        # reactions removed here were already removed from the network, and 'kos_to_remove' would be
        # empty. In contrast, if this method was not called from '_purge_kos', but zero KOs are
        # exclusively associated with reactions removed here, then 'kos_to_remove' would likewise be
        # empty.
        if kos_to_remove:
            removed_cascading_up = self._purge_kos(kos_to_remove)
            for key in ['reaction', 'kegg_reaction', 'ec_number', 'metabolite']:
                removed_cascading_up.pop(key)
        else:
            removed_cascading_up = {
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': []
            }
            if isinstance(self, GenomicNetwork):
                removed_cascading_up['gene'] = []
                if self.proteins:
                    removed_cascading_up['protein'] = []
            elif isinstance(self, PangenomicNetwork):
                removed_cascading_up['gene_cluster'] = []
            else:
                raise AssertionError

        removed = {
            'reaction': removed_reactions,
            'kegg_reaction': removed_kegg_reactions,
            'ec_number': removed_ec_numbers
        }
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def _purge_kos(
        self,
        kos_to_remove: Iterable[str] = None,
        modules_to_remove: Iterable[str] = None,
        pathways_to_remove: Iterable[str] = None,
        hierarchies_to_remove: Iterable[str] = None,
        categories_to_remove: Dict[str, List[Tuple[str]]] = None
    ) -> Dict[str, List]:
        """
        Remove any trace of the given KOs from the network.

        If KEGG modules, pathways, BRITE hierarchies, or BRITE hierarchy categories are provided,
        remove any trace of the KOs that are in these from the network.

        Reactions and metabolites that are only associated with removed KOs are purged. In genomic
        networks, genes that are only associated with removed KOs are purged. In pangenomic
        networks, gene clusters assigned removed KOs are purged. KEGG modules, pathways, BRITE
        hierarchies, and BRITE hierarchy categories only associated with purged KOs are removed.

        Parameters
        ==========
        kos_to_remove : Iterable[str], None
            KO IDs to remove.

        modules_to_remove : Iterable[str], None
            KEGG module IDs to remove, equivalent to giving the KOs in the modules to the argument,
            'kos_to_remove'.

        pathways_to_remove : Iterable[str], None
            KEGG pathway IDs to remove, equivalent to giving the KOs in the pathways to the
            argument, 'kos_to_remove'.

        hierarchies_to_remove : Iterable[str], None
            KEGG BRITE hierarchy IDs to remove, equivalent to giving the KOs in the hierarchies to
            the argument, 'kos_to_remove'.

        categories_to_remove : Dict[str, List[Tuple[str]]], None
            KEGG BRITE hierarchy categories to remove, equivalent to giving the KOs in the
            categories to the argument, 'kos_to_remove'. The dictionary argument is keyed by BRITE
            hierarchy ID and has values that list category tuples. For example, to purge KOs from
            the network contained in the 'ko00001' 'KEGG Orthology (KO)' hierarchy categories,
            '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis /
            Gluconeogenesis [PATH:ko00010]' and '09100 Metabolism >>> 09101 Carbohydrate
            metabolism >>> 00051 Fructose and mannose metabolism [PATH:ko00051]', the dictionary
            argument would need to look like the following: {'ko00001': [('09100 Metabolism', '09101
            Carbohydrate metabolism', '00010 Glycolysis / Gluconeogenesis'), ('09100 Metabolism',
            '09101 Carbohydrate metabolism', '00051 Fructose and mannose metabolism
            [PATH:ko00051]')]}

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            The dictionary examples below are for a genomic network. For a pangenomic network, the
            gene entry is replaced by the gene cluster entry, 'gene_cluster': [<removed GeneCluster
            objects>] or 'gene_cluster': []. The examples show protein entries as if the genomic
            network has been annotated with protein abundances; these are absent for genomic
            networks lacking protein annotations and for pangenomic networks.

            If this method is NOT called from the method, '_purge_reactions', or the method,
            '_purge_genes', then the dictionary will look like the following.
            {
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>],
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>]
            }

            If this method is called from the method, '_purge_reactions', then the dictionary will
            look like the following.
            {
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': [],
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>]
            }

            If this method is called from the method, '_purge_genes', then the dictionary will look
            like the following.
            {
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>],
                'gene': [],
                'protein': []
            }

            If no KOs are removed from the network, then the dictionary will look like the following
            regardless of calling method.
            {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene': [],
                'protein': []
            }
        """
        assert (
            kos_to_remove or
            modules_to_remove or
            pathways_to_remove or
            hierarchies_to_remove or
            categories_to_remove
        )

        if kos_to_remove is None:
            kos_to_remove: List[str] = []
        if modules_to_remove is None:
            modules_to_remove: List[str] = []
        if pathways_to_remove is None:
            pathways_to_remove: List[str] = []
        if hierarchies_to_remove is None:
            hierarchies_to_remove: List[str] = []
        if categories_to_remove is None:
            categories_to_remove: Dict[str, List[Tuple[str]]] = {}

        # Get the KOs to remove from requested modules, pathways, hierarchies, and hierarchy
        # categories.
        for module_id in modules_to_remove:
            try:
                module = self.modules[module_id]
            except KeyError:
                # The requested module is not in the network.
                continue
            kos_to_remove += module.kos
        for pathway_id in pathways_to_remove:
            try:
                pathway = self.pathways[pathway_id]
            except KeyError:
                # The requested pathway is not in the network.
                continue
            kos_to_remove += pathway.kos
        for hierarchy_id in hierarchies_to_remove:
            try:
                hierarchy = self.hierarchies[hierarchy_id]
            except KeyError:
                # The requested hierarchy is not in the network.
                continue
            kos_to_remove += hierarchy.kos
        for hierarchy_id, category_keys in categories_to_remove.items():
            try:
                categorization = self.categories[hierarchy_id][category_keys]
            except KeyError:
                # The requested category is not in the network.
                continue
            category = categorization[-1]
            kos_to_remove += category.kos

        # Remove requested KOs from the network.
        kos_to_remove = set(kos_to_remove)
        removed_kos: List[KO] = []
        for ko_id in kos_to_remove:
            try:
                removed_kos.append(self.kos.pop(ko_id))
            except KeyError:
                # This occurs when the original method called is '_purge_kos', followed by
                # '_purge_genes' or '_purge_gene_clusters', followed by this method again --
                # 'removed_kos' will be empty. Alternatively, this occurs if the KO in
                # 'kos_to_remove' is not in the network.
                pass

        if not removed_kos:
            removed = {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': []
            }
            if isinstance(self, GenomicNetwork):
                removed['gene'] = []
            elif isinstance(self, PangenomicNetwork):
                removed['gene_cluster'] = []
            else:
                raise AssertionError
            return removed

        # Remove modules from the network that exclusively contain removed KOs.
        modules_to_remove.clear()
        for removed_ko in removed_kos:
            for module_id, module in removed_ko.modules.items():
                if module_id in modules_to_remove:
                    # The module has already been considered via another removed KO.
                    continue
                removed_kos_in_module: List[str] = []
                for ko_id, ko in module.kos.items():
                    if ko_id in kos_to_remove:
                        removed_kos_in_module.append(ko_id)
                if len(removed_kos_in_module) == len(module.kos):
                    # Remove the module since all of its KOs were removed from the network.
                    modules_to_remove.append(module_id)
                    continue
                for ko_id in removed_kos_in_module:
                    # Remove obsolete KO references from the retained module.
                    module.kos.pop(ko_id)
        removed_modules: List[KEGGModule] = []
        for module_id in modules_to_remove:
            module = self.modules.pop(module_id)
            # Remove obsolete module references from pathways.
            for pathway in module.pathways.values():
                pathway.modules.pop(module_id)
            removed_modules.append(module)

        # Remove pathways from the network that exclusively contain removed KOs.
        pathways_to_remove.clear()
        for removed_ko in removed_kos:
            for pathway_id, pathway in removed_ko.pathways.items():
                if pathway_id in pathways_to_remove:
                    # The pathway has already been considered via another removed KO.
                    continue
                removed_kos_in_pathway: List[str] = []
                for ko_id, ko in pathway.kos.items():
                    if ko_id in kos_to_remove:
                        removed_kos_in_pathway.append(ko_id)
                if len(removed_kos_in_pathway) == len(pathway.kos):
                    # Remove the pathway since all of its KOs were removed from the network.
                    pathways_to_remove.append(pathway_id)
                    continue
                for ko_id in removed_kos_in_pathway:
                    # Remove obsolete KO references from the retained pathway.
                    pathway.kos.pop(ko_id)
        removed_pathways: List[KEGGPathway] = []
        for pathway_id in pathways_to_remove:
            pathway = self.pathways.pop(pathway_id)
            # Remove obsolete pathway references from modules.
            for module in pathway.modules.values():
                module.pathways.pop(pathway_id)
            removed_pathways.append(pathway)

        # For simplicity's sake, remove hierarchies and hierarchy categories from the network
        # separately.
        # Remove hierarchies from the network that exclusively contain removed KOs.
        hierarchies_to_remove.clear()
        for removed_ko in removed_kos:
            for hierarchy_id in removed_ko.hierarchies:
                if hierarchy_id in hierarchies_to_remove:
                    # The hierarchy has already been considered via another removed KO.
                    continue
                removed_kos_in_hierarchy: List[str] = []
                for ko_id, ko in hierarchy.kos.items():
                    if ko_id in kos_to_remove:
                        removed_kos_in_hierarchy.append(ko_id)
                if len(removed_kos_in_hierarchy) == len(hierarchy.kos):
                    # Remove the hierarchy since all of its KOs were removed from the network.
                    hierarchies_to_remove.append(hierarchy_id)
                    continue
                for ko_id in removed_kos_in_pathway:
                    # Remove obsolete KO references from the retained hierarchy.
                    hierarchy.kos.pop(ko_id)
        removed_hierarchies: List[BRITEHierarchy] = []
        removed_categories: List[BRITECategory] = []
        for hierarchy_id in hierarchies_to_remove:
            hierarchy = self.hierarchies.pop(hierarchy_id)
            # All of the categories of the hierarchy are consequently removed. Pathways with
            # equivalent categories that are here removed would already have been removed.
            categorizations = self.categories.pop(hierarchy_id)
            for categorization in categorizations.values():
                category = categorization[-1]
                removed_categories.append(category)
            removed_hierarchies.append(hierarchy)

        # Remove categories from the network that exclusively contain removed KOs.
        categories_to_remove: List[Tuple[str, Tuple[str]]] = []
        for removed_ko in removed_kos:
            for hierarchy_id, categorizations in removed_ko.hierarchies.items():
                if hierarchy_id in hierarchies_to_remove:
                    # The hierarchy, and thus all of its categories, was already removed.
                    continue
                removed_kos_in_category: List[str] = []
                # Consider the (most specific) category and supercategories containing the KO.
                for category_key, categorization in categorizations.items():
                    is_supercategory_removed = False
                    for rank, category in enumerate(categorization, 1):
                        key = category_key[:rank]
                        if is_supercategory_removed:
                            # This category should be removed since it was already established that
                            # the supercategory in the categorization should be removed.
                            categories_to_remove.append((hierarchy_id, key))
                            continue
                        if (hierarchy_id, key) in categories_to_remove:
                            # The category has already been considered via another removed KO.
                            break
                        for ko_id, ko in category.kos.items():
                            if ko_id in kos_to_remove:
                                removed_kos_in_category.append(ko_id)
                        if len(removed_kos_in_category) == len(category.kos):
                            # Remove the category since all of its KOs were removed from the
                            # network.
                            categories_to_remove.append((hierarchy_id, key))
                            is_supercategory_removed = True
                            continue
                        for ko_id in removed_kos_in_category:
                            # Remove obsolete KO references from the retained category.
                            category.kos.pop(ko_id)
        for hierarchy_id, category_key in categories_to_remove:
            categorization = self.categories[hierarchy_id].pop(category_key)
            category = categorization[-1]
            if category.pathway is not None:
                # The equivalent pathway should already have been removed from the network.
                assert category.pathway.id not in self.pathways
            # Remove the obsolete category reference from the hierarchy.
            category.hierarchy.categories.pop(category_key)
            removed_categories.append(category)

        # Purge reactions from the network that are exclusive to removed KOs.
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
                        # The reaction is associated with a retained KO, so do not remove the
                        # reaction.
                        reactions_to_spare.append(idx)
            for idx in sorted(reactions_to_spare, reverse=True):
                reactions_to_remove.pop(idx)
        if reactions_to_remove:
            removed_cascading_down = self._purge_reactions(reactions_to_remove)
            for key in ['ko', 'module', 'pathway', 'hierarchy', 'category']:
                removed_cascading_down.pop(key)
            if isinstance(self, GenomicNetwork):
                removed_cascading_down.pop('gene')
            elif isinstance(self, PangenomicNetwork):
                removed_cascading_down.pop('gene_cluster')
            else:
                raise AssertionError
        else:
            # This method must have been called "cascading up" from the method, '_purge_reactions',
            # because the reactions that are only associated with the removed KOs were already
            # removed from the network.
            removed_cascading_down = {
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        if isinstance(self, GenomicNetwork):
            # Purge genes from the network that are are only associated with removed KOs.
            genes_to_remove: List[str] = []
            for gcid, gene in self.genes.items():
                gene_kos_to_remove: List[str] = []
                for ko_id in gene.kos:
                    if ko_id in kos_to_remove:
                        gene_kos_to_remove.append(ko_id)
                if len(gene_kos_to_remove) == len(gene.kos):
                    # All KOs matching the gene were removed, so remove it as well.
                    genes_to_remove.append(gcid)
                    continue
                for ko_id in gene_kos_to_remove:
                    gene.kos.pop(ko_id)
                    gene.e_values.pop(ko_id)
            # If this method was called from '_purge_genes', then the genes that are only associated
            # with KOs removed here were already removed from the network, and 'genes_to_remove'
            # would be empty. In contrast, if this method was not called from '_purge_genes', but
            # zero genes are only associated with KOs removed here, then 'genes_to_remove' would
            # likewise be empty.
            if genes_to_remove:
                removed_cascading_up = self._purge_genes(genes_to_remove)
                for key in [
                    'ko',
                    'module',
                    'pathway',
                    'reaction',
                    'kegg_reaction',
                    'ec_number',
                    'metabolite'
                ]:
                    removed_cascading_up.pop(key)
            else:
                removed_cascading_up = {'gene': []}
        elif isinstance(self, PangenomicNetwork):
            gene_clusters_to_remove: List[str] = []
            for gene_cluster_id, cluster in self.gene_clusters.items():
                if cluster.ko.id in kos_to_remove:
                    gene_clusters_to_remove.append(gene_cluster_id)
            # If this method was called from 'purge_gene_clusters' then the gene clusters that are only
            # associated with KOs removed here were already removed from the network, and
            # 'gene_clusters_to_remove' would be empty.
            if gene_clusters_to_remove:
                removed_cascading_up = self._purge_gene_clusters(gene_clusters_to_remove)
                removed_cascading_up.pop('ko')
            else:
                removed_cascading_up = {'gene_cluster': []}
        else:
            raise AssertionError

        removed = {
            'ko': removed_kos,
            'module': removed_modules,
            'pathway': removed_pathways,
            'hierarchy': removed_hierarchies,
            'category': removed_categories
        }
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def _subset_network_by_kos(
        self,
        ko_ids: Iterable[str],
        subnetwork: ReactionNetwork = None
    ) -> ReactionNetwork:
        """
        Subset the network by KOs.

        Parameters
        ==========
        ko_ids : Iterable[str]
            KO IDs to subset.

        subnetwork : ReactionNetwork, None
            This network under construction is provided when the KOs being added to the network
            annotate already subsetted genes or gene clusters.

        Returns
        =======
        ReactionNetwork
            If a 'subnetwork' argument is provided, then that network is returned after
            modification. Otherwise, a new subsetted reaction network is returned.
        """
        if isinstance(self, GenomicNetwork):
            if subnetwork is None:
                subnetwork = GenomicNetwork()
                # Signify that genes annotated by subsetted KOs are to be added to the network.
                subset_referencing_genes = True
            else:
                assert isinstance(subnetwork, GenomicNetwork)
                # Signify that the KOs being added to the network annotate subsetted genes that were
                # already added to the network.
                subset_referencing_genes = False
        elif isinstance(self, PangenomicNetwork):
            if subnetwork is None:
                subnetwork = PangenomicNetwork()
                # Signify that gene clusters annotated by subsetted KOs are to be added to the
                # network.
                subset_referencing_gene_clusters = True
            else:
                assert isinstance(subnetwork, PangenomicNetwork)
                # Signify that the KOs being added to the network annotate subsetted gene clusters
                # that were already added to the network.
                subset_referencing_gene_clusters = False
        else:
            raise AssertionError

        for ko_id in ko_ids:
            try:
                ko = self.kos[ko_id]
            except KeyError:
                # This occurs if the requested KO ID is not in the source network.
                continue

            subsetted_ko = KO()
            subsetted_ko.id = ko.id
            subsetted_ko.name = ko.name
            subsetted_ko.kegg_reaction_aliases = deepcopy(ko.kegg_reaction_aliases)
            subsetted_ko.ec_number_aliases = deepcopy(ko.ec_number_aliases)

            # Add reactions annotating the KO to the subsetted network as new objects, and then
            # reference these objects in the KO object.
            reaction_ids = [reaction_id for reaction_id in ko.reactions]
            subnetwork = self._subset_network_by_reactions(reaction_ids, subnetwork=subnetwork)
            subsetted_ko.reactions = {
                reaction_id: subnetwork.reactions[reaction_id] for reaction_id in reaction_ids
            }

            subnetwork.kos[ko_id] = subsetted_ko

        # Add modules, pathways, and hierarchies/categories to the subnetwork after new KO objects
        # have already been added so that references to the new KO objects can be made from the new
        # module, pathway, hierarchy, and category objects.
        # Add modules to the subnetwork.
        for ko_id, subsetted_ko in subnetwork.kos.items():
            ko = self.kos[ko_id]
            for module_id, module in ko.modules.items():
                if module_id in subnetwork.modules:
                    # The module was already added to the subnetwork via another KO.
                    continue

                subnetwork_module = KEGGModule()
                subnetwork_module.id = module_id
                subnetwork_module.name = module.name

                for module_ko_id in module.kos:
                    try:
                        subsetted_module_ko = subnetwork.kos[module_ko_id]
                    except KeyError:
                        continue
                    subnetwork_module.kos[module_ko_id] = subsetted_module_ko
                    subsetted_module_ko.modules[module_id] = subnetwork_module

        # Add pathways to the subnetwork. Add pathway references to subnetwork modules.
        for ko_id, subsetted_ko in subnetwork.kos.items():
            for pathway_id, pathway in ko.pathways.items():
                if pathway_id in subnetwork.pathways:
                    # The pathway was already added to the subnetwork via another KO.
                    continue

                subnetwork_pathway = KEGGPathway()
                subnetwork_pathway.id = pathway_id
                subnetwork_pathway.name = pathway.name

                for pathway_ko_id in pathway.kos:
                    try:
                        subsetted_pathway_ko = subnetwork.kos[pathway_ko_id]
                    except KeyError:
                        continue
                    subnetwork_pathway.kos[pathway_ko_id] = subsetted_pathway_ko
                    subsetted_pathway_ko.pathways[pathway_ko_id] = subnetwork_pathway

                for module_id in pathway.modules:
                    try:
                        subnetwork_module = subnetwork.modules[module_id]
                    except KeyError:
                        continue
                    subnetwork_pathway.modules[module_id] = subnetwork_module
                    subnetwork_module.pathways[pathway_id] = subnetwork_pathway

        # Add hierarchies and categories to the subnetwork.
        for ko_id, subsetted_ko in subnetwork.kos.items():
            for hierarchy_id, categorizations in ko.hierarchies.items():
                if hierarchy_id in subnetwork.hierarchies:
                    # The hierarchy was already added to the subnetwork via another KO.
                    continue

                # Add the hierarchy to the subnetwork.
                hierarchy = self.hierarchies[hierarchy_id]
                subnetwork_hierarchy = BRITEHierarchy()
                subnetwork_hierarchy.id = hierarchy_id
                subnetwork_hierarchy.name = hierarchy.name

                for hierarchy_ko_id in hierarchy.kos:
                    try:
                        subsetted_hierarchy_ko = subnetwork.kos[hierarchy_ko_id]
                    except KeyError:
                        continue
                    subnetwork_hierarchy.kos[hierarchy_ko_id] = subsetted_hierarchy_ko
                    subsetted_hierarchy_ko.hierarchies[hierarchy_id] = {}

                # Create new category objects for the subnetwork.
                subnetwork_hierarchy_categories = subnetwork.categories[hierarchy_id]
                for category_key, categorization in categorizations.items():
                    # Also consider each supercategory of the (most specific) category referenced by
                    # the KO.
                    subnetwork_categorization: List[BRITECategory] = []
                    for rank, category in enumerate(categorization, 1):
                        key = category_key[:rank]
                        category_id = category.id
                        try:
                            # The category was already encountered via another KO.
                            subnetwork_category = subnetwork_hierarchy_categories[key]
                            is_new_category = False
                        except KeyError:
                            subnetwork_category = BRITECategory()
                            is_new_category = True
                        subnetwork_categorization.append(subnetwork_category)

                        if not is_new_category:
                            continue

                        subnetwork_category.id = category_id
                        subnetwork_category.name = category.name
                        subnetwork_category.hierarchy = subnetwork_hierarchy

                        # Relate the category to an equivalent pathway in the subnetwork.
                        pathway = category.pathway
                        if pathway is not None:
                            subnetwork_pathway = subnetwork.pathways[pathway.id]
                            subnetwork_pathway.category = subnetwork_category
                            subnetwork_category.pathway = subnetwork_pathway

                        # Use subcategory IDs as placeholders to later fill in subcategory
                        # object references.
                        for subcategory in category.subcategories:
                            subnetwork_category.subcategories.append(subcategory.id)

                        subnetwork_hierarchy.categories[key] = subnetwork_category

                        subnetwork_hierarchy_categories[subnetwork_category][
                            key
                        ] = subnetwork_category

                        # Add the new category object to the hierarchies attribute of each
                        # subsetted KO in the category.
                        for category_ko_id in category.kos:
                            subsetted_category_ko = subnetwork.kos[category_ko_id]
                            subnetwork_category.kos[category_ko_id] = subsetted_category_ko
                            subsetted_category_ko.hierarchies[hierarchy_id][key] = tuple(
                                subnetwork_categorization
                            )

                        # Fill in subnetwork supercategory/subcategory relationships.
                        if len(subnetwork_categorization) > 1:
                            subnetwork_supercategory = subnetwork_categorization[-2]
                            subnetwork_category.supercategory = subnetwork_supercategory
                            # Replace the placeholder subcategory ID with the object.
                            subnetwork_supercategory.subcategories[
                                subnetwork_supercategory.subcategories.index(category_id)
                            ] = subnetwork_category

        if isinstance(self, GenomicNetwork):
            if subset_referencing_genes:
                # Add genes that are annotated by the subsetted KOs to the network.
                self._subset_genes_via_kos(subnetwork)
        elif isinstance(self, PangenomicNetwork):
            if subset_referencing_gene_clusters:
                # Add gene clusters that are annotated by the subsetted KOs to the network.
                self._subset_gene_clusters_via_kos(subnetwork)

        return subnetwork

    def _subset_network_by_reactions(
        self,
        reaction_ids: Iterable[str],
        subnetwork: ReactionNetwork = None
    ) -> ReactionNetwork:
        """
        Subset the network by ModelSEED reactions.

        Parameters
        ==========
        reaction_ids : Iterable[str]
            ModelSEED reaction IDs to subset.

        subnetwork : ReactionNetwork, None
            This network under construction is provided when the reactions being added to the
            network annotate already subsetted KOs.

        Returns
        =======
        ReactionNetwork
            If a 'subnetwork' argument is provided, then that network is returned after
            modification. Otherwise, a new subsetted reaction network is returned.
        """
        if isinstance(self, GenomicNetwork):
            if subnetwork is None:
                subnetwork = GenomicNetwork()
                # Signify that KOs annotated by subsetted reactions are to be added to the network.
                subset_referencing_kos = True
            else:
                assert isinstance(subnetwork, GenomicNetwork)
                # Signify that the reactions being added to the network annotate subsetted KOs that
                # were already added to the network.
                subset_referencing_kos = False
        elif isinstance(self, PangenomicNetwork):
            if subnetwork is None:
                subnetwork = PangenomicNetwork()
                # Signify that KOs annotated by subsetted reactions are to be added to the network.
                subset_referencing_kos = True
            else:
                assert isinstance(subnetwork, PangenomicNetwork)
                # Signify that the reactions being added to the network annotate subsetted KOs that
                # were already added to the network.
                subset_referencing_kos = False
        else:
            raise AssertionError

        # Copy the network attributes mapping reaction aliases.
        kegg_modelseed_aliases: Dict[str, List[str]] = {}
        ec_number_modelseed_aliases: Dict[str, List[str]] = {}

        for reaction_id in reaction_ids:
            try:
                reaction = self.reactions[reaction_id]
            except KeyError:
                # This occurs if the requested reaction is not in the source network.
                continue

            # Copy the reaction object, including referenced metabolite objects, from the source
            # network.
            subsetted_reaction: ModelSEEDReaction = deepcopy(reaction)
            subnetwork.reactions[reaction_id] = subsetted_reaction
            # Record the metabolites involved in the reaction, and add them to the network.
            for metabolite in subsetted_reaction.compounds:
                compound_id = metabolite.modelseed_id
                subnetwork.metabolites[compound_id] = metabolite

            try:
                subnetwork.modelseed_kegg_aliases[reaction_id] += list(reaction.kegg_aliases)
            except KeyError:
                subnetwork.modelseed_kegg_aliases[reaction_id] = list(reaction.kegg_aliases)

            try:
                subnetwork.modelseed_ec_number_aliases[reaction_id] += list(
                    reaction.ec_number_aliases
                )
            except KeyError:
                subnetwork.modelseed_ec_number_aliases[reaction_id] = list(
                    reaction.ec_number_aliases
                )

            for kegg_id in reaction.kegg_aliases:
                try:
                    kegg_modelseed_aliases[kegg_id].append(reaction_id)
                except KeyError:
                    kegg_modelseed_aliases[kegg_id] = [reaction_id]

            for ec_number in reaction.ec_number_aliases:
                try:
                    ec_number_modelseed_aliases[ec_number].append(reaction_id)
                except KeyError:
                    ec_number_modelseed_aliases[ec_number] = [reaction_id]

        if subnetwork.kegg_modelseed_aliases:
            for kegg_id, modelseed_ids in kegg_modelseed_aliases.items():
                try:
                    subnetwork.kegg_modelseed_aliases[kegg_id] += modelseed_ids
                except KeyError:
                    subnetwork.kegg_modelseed_aliases[kegg_id] = modelseed_ids
        else:
            subnetwork.kegg_modelseed_aliases = kegg_modelseed_aliases

        if subnetwork.ec_number_modelseed_aliases:
            for ec_number, modelseed_ids in ec_number_modelseed_aliases.items():
                try:
                    subnetwork.ec_number_modelseed_aliases[ec_number] += modelseed_ids
                except KeyError:
                    subnetwork.ec_number_modelseed_aliases[ec_number] = modelseed_ids
        else:
            subnetwork.ec_number_modelseed_aliases = ec_number_modelseed_aliases

        if subset_referencing_kos:
            # Add KOs that are annotated by the subsetted reactions to the network.
            self._subset_kos_via_reactions(subnetwork)

        return subnetwork

    def _subset_kos_via_reactions(self, subnetwork: ReactionNetwork) -> None:
        """
        Add KOs that are annotated with subsetted reactions to the subsetted network.

        Then add genes that are annotated with these added KOs to the subsetted network.

        Parameters
        ==========
        subnetwork : ReactionNetwork
            The subsetted reaction network under construction.

        Returns
        =======
        None
        """
        if isinstance(self, GenomicNetwork):
            assert isinstance(subnetwork, GenomicNetwork)
        elif isinstance(self, PangenomicNetwork):
            assert isinstance(subnetwork, PangenomicNetwork)
        else:
            raise AssertionError

        subsetted_reaction_ids = list(subnetwork.reactions)
        for ko_id, ko in self.kos.items():
            # Check all KOs in the source network for subsetted reactions.
            subsetted_ko = None
            for reaction_id in ko.reactions:
                if reaction_id not in subsetted_reaction_ids:
                    # The KO is not annotated by the subsetted reaction.
                    continue

                if not subsetted_ko:
                    # Create a new KO object for the subsetted KO. The subsetted KO object would
                    # already have been created had another subsetted reaction been among the
                    # reactions annotating the KO.
                    subsetted_ko = KO()
                    subsetted_ko.id = ko_id
                    subsetted_ko.name = ko.name
                subsetted_ko.reactions[reaction_id] = subnetwork.reactions[reaction_id]
                subsetted_ko.kegg_reaction_aliases = deepcopy(ko.kegg_reaction_aliases)
                subsetted_ko.ec_number_aliases = deepcopy(ko.ec_number_aliases)

            if subsetted_ko:
                subnetwork.kos[ko_id] = subsetted_ko

        if isinstance(self, GenomicNetwork):
            # Add genes that are annotated with the added KOs to the subsetted network.
            self._subset_genes_via_kos(subnetwork)
        elif isinstance(self, PangenomicNetwork):
            # Add gene clusters that are annotated with the added KOs to the subsetted network.
            self._subset_gene_clusters_via_kos(subnetwork)

    def _merge_network(self, network: ReactionNetwork, merged_network: ReactionNetwork) -> None:
        """
        This method is used in the process of merging the network with another network to produce a
        merged network, and contains steps common to different types of network: merge the
        attributes of the networks BESIDES genes (in a GenomicNetwork) / gene clusters (in a
        PangenomicNetwork) and protein abundances (which can only be stored in a GenomicNetwork).

        Parameters
        ==========
        network : ReactionNetwork
            The other reaction network being merged.

        merged_network : ReactionNetwork
            The merged reaction network under construction.

        Returns
        =======
        None
        """
        if isinstance(network, GenomicNetwork):
            assert isinstance(merged_network, GenomicNetwork)
        else:
            assert isinstance(merged_network, PangenomicNetwork)

        # Add metabolites to the merged network, starting with metabolites in the first network and
        # continuing with metabolites exclusive to the second network. Assume objects representing
        # the same metabolites in both networks properly have identical attributes.
        for metabolite_id, first_metabolite in self.metabolites.items():
            merged_metabolite = ModelSEEDCompound()
            merged_metabolite.modelseed_id = metabolite_id
            merged_metabolite.modelseed_name = first_metabolite.modelseed_name
            merged_metabolite.kegg_aliases = first_metabolite.kegg_aliases
            merged_metabolite.charge = first_metabolite.charge
            merged_metabolite.formula = first_metabolite.formula
            abundances: Dict[str, float] = getattr(first_metabolite, 'abundances', None)
            if abundances is None:
                continue
            merged_metabolite.abundances = abundances.copy()

            merged_network.metabolites[metabolite_id] = merged_metabolite

        for metabolite_id in set(network.metabolites).difference(self.metabolites):
            second_metabolite = network.metabolites[metabolite_id]

            merged_metabolite = ModelSEEDCompound()
            merged_metabolite.modelseed_id = metabolite_id
            merged_metabolite.modelseed_name = second_metabolite.modelseed_name
            merged_metabolite.kegg_aliases = second_metabolite.kegg_aliases
            merged_metabolite.charge = second_metabolite.charge
            merged_metabolite.formula = second_metabolite.formula
            abundances: Dict[str, float] = getattr(second_metabolite, 'abundances', None)
            if abundances is None:
                continue
            merged_metabolite.abundances = abundances.copy()

            merged_network.metabolites[metabolite_id] = merged_metabolite

        # Add reactions to the merged network, starting with reactions in the first network and
        # continuing with reactions exclusive to the second network. Assume objects representing the
        # same reactions in both networks properly have identical attributes.

        # Determine network attributes mapping reaction aliases.
        kegg_modelseed_aliases: Dict[str, List[str]] = {}
        ec_number_modelseed_aliases: Dict[str, List[str]] = {}

        for reaction_id, first_reaction in self.reactions.items():
            merged_reaction = ModelSEEDReaction()
            merged_reaction.modelseed_id = reaction_id
            merged_reaction.modelseed_name = first_reaction.modelseed_name
            merged_reaction.kegg_aliases = first_reaction.kegg_aliases
            merged_reaction.ec_number_aliases = first_reaction.ec_number_aliases
            metabolites = []
            for metabolite in first_reaction.compounds:
                metabolites.append(merged_network.metabolites[metabolite.modelseed_id])
            merged_reaction.compounds = tuple(metabolites)
            merged_reaction.coefficients = first_reaction.coefficients
            merged_reaction.compartments = first_reaction.compartments
            merged_reaction.reversibility = first_reaction.reversibility

            merged_network.reactions[reaction_id] = merged_reaction

            try:
                merged_network.modelseed_kegg_aliases[reaction_id] += list(
                    first_reaction.kegg_aliases
                )
            except KeyError:
                merged_network.modelseed_kegg_aliases[reaction_id] = list(
                    first_reaction.kegg_aliases
                )

            try:
                merged_network.modelseed_ec_number_aliases[reaction_id] += list(
                    first_reaction.ec_number_aliases
                )
            except KeyError:
                merged_network.modelseed_ec_number_aliases[reaction_id] = list(
                    first_reaction.ec_number_aliases
                )

            for kegg_id in first_reaction.kegg_aliases:
                try:
                    kegg_modelseed_aliases[kegg_id].append(reaction_id)
                except KeyError:
                    kegg_modelseed_aliases[kegg_id] = [reaction_id]

            for ec_number in first_reaction.ec_number_aliases:
                try:
                    ec_number_modelseed_aliases[ec_number].append(reaction_id)
                except KeyError:
                    ec_number_modelseed_aliases[ec_number] = [reaction_id]

        for reaction_id in set(network.reactions).difference(self.reactions):
            second_reaction = network.reactions[reaction_id]

            merged_reaction = ModelSEEDReaction()
            merged_reaction.modelseed_id = reaction_id
            merged_reaction.modelseed_name = second_reaction.modelseed_name
            merged_reaction.kegg_aliases = second_reaction.kegg_aliases
            merged_reaction.ec_number_aliases = second_reaction.ec_number_aliases
            metabolites = []
            for metabolite in second_reaction.compounds:
                metabolites.append(merged_network.metabolites[metabolite.modelseed_id])
            merged_reaction.compounds = tuple(metabolites)
            merged_reaction.coefficients = second_reaction.coefficients
            merged_reaction.compartments = second_reaction.compartments
            merged_reaction.reversibility = second_reaction.reversibility

            merged_network.reactions[reaction_id] = merged_reaction

            try:
                merged_network.modelseed_kegg_aliases[reaction_id] += list(
                    second_reaction.kegg_aliases
                )
            except KeyError:
                merged_network.modelseed_kegg_aliases[reaction_id] = list(
                    second_reaction.kegg_aliases
                )

            try:
                merged_network.modelseed_ec_number_aliases[reaction_id] += list(
                    second_reaction.ec_number_aliases
                )
            except KeyError:
                merged_network.modelseed_ec_number_aliases[reaction_id] = list(
                    second_reaction.ec_number_aliases
                )

            for kegg_id in second_reaction.kegg_aliases:
                try:
                    kegg_modelseed_aliases[kegg_id].append(reaction_id)
                except KeyError:
                    kegg_modelseed_aliases[kegg_id] = [reaction_id]

            for ec_number in second_reaction.ec_number_aliases:
                try:
                    ec_number_modelseed_aliases[ec_number].append(reaction_id)
                except KeyError:
                    ec_number_modelseed_aliases[ec_number] = [reaction_id]

        if merged_network.kegg_modelseed_aliases:
            for kegg_id, modelseed_ids in kegg_modelseed_aliases.items():
                try:
                    merged_network.kegg_modelseed_aliases[kegg_id] += modelseed_ids
                except KeyError:
                    merged_network.kegg_modelseed_aliases[kegg_id] = modelseed_ids
        else:
            merged_network.kegg_modelseed_aliases = kegg_modelseed_aliases

        if merged_network.ec_number_modelseed_aliases:
            for ec_number, modelseed_ids in ec_number_modelseed_aliases.items():
                try:
                    merged_network.ec_number_modelseed_aliases[ec_number] += modelseed_ids
                except KeyError:
                    merged_network.ec_number_modelseed_aliases[ec_number] = modelseed_ids
        else:
            merged_network.ec_number_modelseed_aliases = ec_number_modelseed_aliases

        # Add KOs to the merged network, first adding KOs present in both source networks, and then
        # adding KOs present exclusively in each source network.
        first_ko_ids = set(self.kos)
        second_ko_ids = set(network.kos)

        for ko_id in first_ko_ids.intersection(second_ko_ids):
            first_ko = self.kos[ko_id]
            second_ko = network.kos[ko_id]

            merged_ko = KO()
            merged_ko.id = ko_id
            merged_ko.name = first_ko.name

            # The new object representing the KO in the merged network should have all reaction
            # annotations from both source KO objects, as these objects can have different reaction
            # references.
            reaction_ids = set(first_ko.reactions).union(set(second_ko.reactions))
            merged_ko.reactions = {
                reaction_id: merged_network.reactions[reaction_id] for reaction_id in reaction_ids
            }
            for reaction_id in reaction_ids:
                try:
                    merged_ko.kegg_reaction_aliases[reaction_id] = first_ko.kegg_reaction_aliases[
                        reaction_id
                    ]
                except KeyError:
                    # The reaction has no KO KEGG REACTION aliases.
                    pass
                try:
                    merged_ko.ec_number_aliases[reaction_id] = first_ko.ec_number_aliases[
                        reaction_id
                    ]
                except KeyError:
                    # The reaction has no KO KEGG REACTION aliases.
                    pass

            merged_network.kos[ko_id] = merged_ko

        for ko_id in first_ko_ids.difference(second_ko_ids):
            first_ko = self.kos[ko_id]

            ko = KO()
            ko.id = ko_id
            ko.name = first_ko.name
            ko.reactions = {
                reaction_id: merged_network.reactions[reaction_id]
                for reaction_id in first_ko.reactions
            }
            ko.kegg_reaction_aliases = deepcopy(first_ko.kegg_reaction_aliases)
            ko.ec_number_aliases = deepcopy(first_ko.ec_number_aliases)

            merged_network.kos[ko_id] = ko

        for ko_id in second_ko_ids.difference(first_ko_ids):
            second_ko = network.kos[ko_id]

            ko = KO()
            ko.id = ko_id
            ko.name = second_ko.name
            ko.reactions = {
                reaction_id: merged_network.reactions[reaction_id]
                for reaction_id in second_ko.reactions
            }
            ko.kegg_reaction_aliases = deepcopy(second_ko.kegg_reaction_aliases)
            ko.ec_number_aliases = deepcopy(second_ko.ec_number_aliases)

            merged_network.kos[ko_id] = ko

        # Add modules to the merged network, first adding modules present in both source networks,
        # and then adding modules present exclusively in each source network.
        first_module_ids = set(self.modules)
        second_module_ids = set(network.modules)

        for module_id in first_module_ids.intersection(second_module_ids):
            first_module = self.modules[module_id]
            second_module = network.modules[module_id]

            merged_module = KEGGModule()
            merged_module.id = module_id
            merged_module.name = first_module.name
            for ko_id in set(first_module.kos).union(set(second_module.kos)):
                merged_module.kos[ko_id] = merged_ko = merged_network.kos[ko_id]
                merged_ko.modules[module_id] = merged_module

            merged_network.modules[module_id] = merged_module
            # TODO: Add merged pathway references to modules.

        for module_id in first_module_ids.difference(second_module_ids):
            first_module = self.modules[module_id]

            module = KEGGModule()
            module.id = module_id
            module.name = first_module.name
            for ko_id in first_module.kos:
                module.kos[ko_id] = merged_ko = merged_network.kos[ko_id]
                merged_ko.modules[module_id] = module

            merged_network.modules[module_id] = module

        for module_id in second_module_ids.difference(first_module_ids):
            second_module = network.modules[module_id]

            module = KEGGModule()
            module.id = module_id
            module.name = second_module.name
            for ko_id in second_module.kos:
                module.kos[ko_id] = merged_ko = merged_network.kos[ko_id]
                merged_ko.modules[module_id] = module

            merged_network.modules[module_id] = module

        # Add hierarchies to the merged network, first adding hierarchies present in the first
        # network, and then adding hierarchies present exclusively in the second network.
        for hierarchy_id, first_hierarchy in self.hierarchies.items():
            hierarchy = BRITEHierarchy()
            hierarchy.id = hierarchy_id
            hierarchy.name = first_hierarchy.name

            merged_network.hierarchies[hierarchy_id] = hierarchy

        exclusive_second_hierarchy_ids = set(network.hierarchies).difference(set(self.hierarchies))
        for hierarchy_id in exclusive_second_hierarchy_ids:
            second_hierarchy = network.hierarchies[hierarchy_id]

            hierarchy = BRITEHierarchy()
            hierarchy.id = hierarchy_id
            hierarchy.name = second_hierarchy.name

            merged_network.hierarchies[hierarchy_id] = hierarchy

        # Add categories to the merged network, first adding categories present in both source
        # networks, and then adding categories present exclusively in each source network.
        first_category_ids = set(self.categories)
        second_category_ids = set(network.categories)

        for category_id in first_category_ids.intersection(second_category_ids):
            first_category = self.categories[category_id]
            second_category = network.categories[category_id]

            merged_category = BRITECategory()
            merged_category.id = first_category.id
            merged_category.name = first_category.name
            merged_category.hierarchy = merged_network.hierarchies[first_category.hierarchy.id]
            for ko_id in set(first_category.kos).union(set(second_category.kos)):
                merged_category.kos[ko_id] = merged_network.kos[ko_id]

            merged_network.categories[category_id] = merged_category

        for category_id in first_category_ids.difference(second_category_ids):
            first_category = self.categories[category_id]

            category = BRITECategory()
            category.id = category_id
            category.name = first_category.name
            category.hierarchy = merged_network.hierarchies[first_category.hierarchy.id]
            for ko_id in first_category.kos:
                category.kos[ko_id] = merged_network.kos[ko_id]

            merged_network.categories[category_id] = category

        for category_id in second_category_ids.difference(first_category_ids):
            second_category = network.categories[category_id]

            category = BRITECategory()
            category.id = category_id
            category.name = second_category.name
            category.hierarchy = merged_network.hierarchies[second_category.hierarchy.id]
            for ko_id in second_category.kos:
                category.kos[ko_id] = merged_network.kos[ko_id]

            merged_network.categories[category_id] = category

        # Add pathways to the merged network, first adding pathways present in both source networks,
        # and then adding pathways present exclusively in each source network.
        first_pathway_ids = set(self.pathways)
        second_pathway_ids = set(network.pathways)

        for pathway_id in first_pathway_ids.intersection(second_pathway_ids):
            first_pathway = self.pathways[pathway_id]
            second_pathway = network.pathways[pathway_id]

            merged_pathway = KEGGPathway()
            merged_pathway.id = pathway_id
            merged_pathway.name = first_pathway.name
            for ko_id in set(first_pathway.kos).intersection(set(second_pathway.kos)):
                merged_pathway.kos[ko_id] = merged_ko = merged_network.kos[ko_id]
                merged_ko.pathways[pathway_id] = merged_pathway
            if first_pathway.category is not None:
                merged_pathway.category = merged_network.categories[first_pathway.category.id]

            merged_network.pathways[pathway_id] = merged_pathway

        for pathway_id in first_pathway_ids.difference(second_pathway_ids):
            first_pathway = self.pathways[pathway_id]

            pathway = KEGGPathway()
            pathway.id = pathway_id
            pathway.name = first_pathway.name
            for ko_id in first_pathway.kos:
                pathway.kos[ko_id] = merged_ko = merged_network.kos[ko_id]
                merged_ko.pathways[pathway_id] = pathway
            if first_pathway.category is not None:
                pathway.category = merged_network.categories[first_pathway.category.id]

            merged_network.pathways[pathway_id] = pathway

        for pathway_id in second_pathway_ids.difference(first_pathway_ids):
            second_pathway = network.pathways[pathway_id]

            pathway = KEGGPathway()
            pathway.id = pathway_id
            pathway.name = second_pathway.name
            for ko_id in second_pathway.kos:
                pathway.kos[ko_id] = merged_ko = merged_network.kos[ko_id]
                merged_ko.pathways[pathway_id] = pathway
            if second_pathway.category is not None:
                pathway.category = merged_network.categories[second_pathway.category.id]

            merged_network.pathways[pathway_id] = pathway

        # Record BRITE hierarchical categorizations of each KO.
        for ko_id, first_ko in self.kos.items():
            merged_ko = merged_network.kos[ko_id]
            for hierarchy_id, categories in first_ko.hierarchies.items():
                merged_ko.hierarchies[hierarchy_id] = merged_categories = []
                for categorization in categories:
                    merged_categorization = []
                    for category in categorization:
                        merged_categorization.append(merged_network.categories[category.id])
                    merged_categories.append(tuple(merged_categorization))

        for ko_id in second_ko_ids.difference(first_ko_ids):
            second_ko = network.kos[ko_id]
            merged_ko = merged_network.kos[ko_id]
            for hierarchy_id, categories in second_ko.hierarchies.items():
                merged_ko.hierarchies[hierarchy_id] = merged_categories = []
                for categorization in categories:
                    merged_categorization = []
                    for category in categorization:
                        merged_categorization.append(merged_network.categories[category.id])
                    merged_categories.append(tuple(merged_categorization))

        # Record the pathway membership of modules.
        for pathway_id, pathway in self.pathways.items():
            merged_pathway = merged_network.pathways[pathway_id]
            for module_id in pathway.modules:
                merged_pathway.modules[module_id] = merged_network.modules[module_id]

        for pathway_id in second_pathway_ids.difference(first_pathway_ids):
            second_pathway = network.pathways[pathway_id]
            merged_pathway = merged_network.pathways[pathway_id]
            for module_id in second_pathway.modules:
                merged_pathway.modules[module_id] = merged_network.modules[module_id]

        for module_id, module in self.modules.items():
            merged_module = merged_network.modules[module_id]
            for pathway_id in module.pathways:
                merged_module.pathways[pathway_id] = merged_network.pathways[pathway_id]

        for module_id in second_module_ids.difference(first_module_ids):
            second_module = network.modules[module_id]
            merged_module = merged_network.modules[module_id]
            for pathway_id in module.pathways:
                merged_module.pathways[pathway_id] = merged_network.pathways[pathway_id]

        # Record categories of each BRITE hierarchy.
        for hierarchy_id, hierarchy in self.hierarchies.items():
            merged_hierarchy = merged_network.hierarchies[hierarchy_id]
            for categorization in hierarchy.categories:
                merged_categorization = []
                for category in categorization:
                    merged_categorization.append(merged_network.categories[category.id])
                merged_hierarchy.categories.append(tuple(merged_categorization))

        for hierarchy_id in exclusive_second_hierarchy_ids:
            second_hierarchy = network.hierarchies[hierarchy_id]
            merged_hierarchy = merged_network.hierarchies[hierarchy_id]
            for categorization in second_hierarchy.categories:
                merged_categorization = []
                for category in categorization:
                    merged_categorization.append(merged_network.categories[category.id])
                merged_hierarchy.categories.append(tuple(merged_categorization))

        # Finish filling out attributes of merged BRITE categories.
        for category_id, category in self.categories.items():
            merged_category = merged_network.categories[category_id]
            # Record any supercategory of the category.
            supercategory = category.supercategory
            if supercategory is not None:
                merged_category.supercategory = merged_network.categories[supercategory.id]
            # Record any subcategories of the category.
            for subcategory in category.subcategories:
                merged_category.subcategories.append(merged_network.categories[subcategory.id])
            # Record any pathway equivalent to the category.
            pathway = category.pathway
            if pathway is not None:
                merged_category.pathway = merged_network.pathways[pathway.id]

        for category_id in second_category_ids.difference(first_category_ids):
            second_category = network.categories[category_id]
            merged_category = merged_network.categories[category_id]
            # Record any supercategory of the category.
            supercategory = second_category.supercategory
            if supercategory is not None:
                merged_category.supercategory = merged_network.categories[supercategory.id]
            # Record any subcategories of the category.
            for subcategory in second_category.subcategories:
                merged_category.subcategories.append(merged_network.categories[subcategory.id])
            # Record any pathway equivalent to the category.
            pathway = second_category.pathway
            if pathway is not None:
                merged_category.pathway = merged_network.pathways[pathway.id]

    def _get_common_overview_statistics(
        self,
        stats: Union[GenomicNetworkStats, PangenomicNetworkStats]
    ) -> None:
        """
        Calculate overview statistics that are found the same way for both genomic and pangenomic
        networks.

        Parameters
        ==========
        stats : Union[GenomicNetworkStats, PangenomicNetworkStats]
            Network statistics are stored in a dictionary of dictionaries. Keys in the outer
            dictionary are "classes" of network statistics. Keys in the inner dictionary are
            the names of the statistics themselves.

        Returns
        =======
        None
        """
        self.progress.new("Counting reactions and KO sources")
        self.progress.update("...")
        stats['Reactions and KO sources'] = stats_group = {}

        stats_group['Reactions in network'] = len(self.reactions)
        reaction_counts = []
        for ko in self.kos.values():
            reaction_counts.append(len(ko.reactions))
        stats_group['Mean reactions per KO'] = round(np.mean(reaction_counts), 1)
        stats_group['Stdev reactions per KO'] = round(np.std(reaction_counts), 1)
        stats_group['Max reactions per KO'] = max(reaction_counts)

        self.progress.end()

        self.progress.new("Counting reactions from each alias source")
        self.progress.update("...")
        stats['Reaction alias sources'] = stats_group = {}

        kegg_aliased_modelseed_reaction_ids = []
        for modelseed_reaction_id, kegg_reaction_ids in self.modelseed_kegg_aliases.items():
            if len(kegg_reaction_ids) > 0:
                kegg_aliased_modelseed_reaction_ids.append(modelseed_reaction_id)
        ec_number_aliased_modelseed_reaction_ids = []
        for modelseed_reaction_id, ec_numbers in self.modelseed_ec_number_aliases.items():
            if len(ec_numbers) > 0:
                ec_number_aliased_modelseed_reaction_ids.append(modelseed_reaction_id)
        kegg_reaction_source_count = len(kegg_aliased_modelseed_reaction_ids)
        ec_number_source_count = len(ec_number_aliased_modelseed_reaction_ids)
        both_source_count = len(
            set(kegg_aliased_modelseed_reaction_ids).intersection(
                set(ec_number_aliased_modelseed_reaction_ids)
            )
        )
        stats_group['Reactions aliased by KEGG reaction'] = kegg_reaction_source_count
        stats_group['Reactions aliased by EC number'] = ec_number_source_count
        stats_group['Rxns aliased by both KEGG rxn & EC number'] = both_source_count
        stats_group['Reactions aliased only by KEGG reaction'] = (
            kegg_reaction_source_count - both_source_count
        )
        stats_group['Reactions aliased only by EC number'] = (
            ec_number_source_count - both_source_count
        )

        stats_group['KEGG reactions contributing to network'] = len(self.kegg_modelseed_aliases)
        reaction_counts = []
        for modelseed_reaction_ids in self.kegg_modelseed_aliases.values():
            reaction_counts.append(len(modelseed_reaction_ids))
        stats_group['Mean reactions per KEGG reaction'] = round(np.mean(reaction_counts), 1)
        stats_group['Stdev reactions per KEGG reaction'] = round(np.std(reaction_counts), 1)
        stats_group['Max reactions per KEGG reaction'] = max(reaction_counts)

        stats_group['EC numbers contributing to network'] = len(self.ec_number_modelseed_aliases)
        reaction_counts = []
        for modelseed_reaction_ids in self.ec_number_modelseed_aliases.values():
            reaction_counts.append(len(modelseed_reaction_ids))
        stats_group['Mean reactions per EC number'] = round(np.mean(reaction_counts), 1)
        stats_group['Stdev reactions per EC number'] = round(np.std(reaction_counts), 1)
        stats_group['Max reactions per EC number'] = max(reaction_counts)

        self.progress.end()

        self.progress.new("Counting reactions and metabolites by property")
        self.progress.update("...")
        stats['Reaction and metabolite properties'] = stats_group = {}

        reversible_count = 0
        irreversible_count = 0
        cytoplasmic_compound_ids = []
        extracellular_compound_ids = []
        consumed_compound_ids = []
        produced_compound_ids = []
        compound_reaction_counts = {}
        for reaction in self.reactions.values():
            if reaction.reversibility:
                reversible_count += 1
            else:
                irreversible_count += 1
            encountered_compound_ids = []
            for compartment, coefficient, compound in zip(
                reaction.compartments, reaction.coefficients, reaction.compounds
            ):
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
        stats_group['Reversible reactions'] = reversible_count
        stats_group['Irreversible reactions'] = irreversible_count
        cytoplasmic_compound_ids = set(cytoplasmic_compound_ids)
        extracellular_compound_ids = set(extracellular_compound_ids)
        stats_group['Metabolites in network'] = metabolite_count = len(self.metabolites)
        stats_group['Cytoplasmic metabolites'] = len(cytoplasmic_compound_ids)
        stats_group['Extracellular metabolites'] = len(extracellular_compound_ids)
        stats_group['Exclusively cytoplasmic metabolites'] = len(
            cytoplasmic_compound_ids.difference(extracellular_compound_ids)
        )
        stats_group['Exclusively extracellular metabolites'] = len(
            extracellular_compound_ids.difference(cytoplasmic_compound_ids)
        )
        stats_group['Cytoplasmic/extracellular metabolites'] = len(
            cytoplasmic_compound_ids.intersection(extracellular_compound_ids)
        )
        consumed_compound_ids = set(consumed_compound_ids)
        produced_compound_ids = set(produced_compound_ids)
        stats_group['Consumed metabolites'] = len(consumed_compound_ids)
        stats_group['Produced metabolites'] = len(produced_compound_ids)
        stats_group['Both consumed & produced metabolites'] = len(
            consumed_compound_ids.intersection(produced_compound_ids)
        )
        stats_group['Exclusively consumed metabolites'] = len(
            consumed_compound_ids.difference(produced_compound_ids)
        )
        stats_group['Exclusively produced metabolites'] = len(
            produced_compound_ids.difference(consumed_compound_ids)
        )
        metabolite_reaction_counts = collections.Counter(compound_reaction_counts.values())
        one_reaction_count = metabolite_reaction_counts[1]
        stats_group['Metabolites consumed or produced by 1 rxn'] = one_reaction_count
        two_reactions_count = metabolite_reaction_counts[2]
        stats_group['Metabolites consumed or produced by 2 rxns'] = two_reactions_count
        three_plus_reactions_count = metabolite_count - one_reaction_count - two_reactions_count
        stats_group['Metabolites consumed or produced by 3+ rxns'] = three_plus_reactions_count

        self.progress.end()

    def _print_common_overview_statistics(
        self,
        stats: Union[GenomicNetworkStats, PangenomicNetworkStats]
    ) -> None:
        """
        Print overview statistics that are the same for both genomic and pangenomic networks.

        Parameters
        ==========
        stats : Union[GenomicNetworkStats, PangenomicNetworkStats]
            Network statistics are stored in a dictionary of dictionaries. Keys in the outer
            dictionary are "classes" of network statistics. Keys in the inner dictionary are
            the names of the statistics themselves.

        Returns
        =======
        None
        """
        self.run.info_single("ModelSEED reactions in network and KO sources")
        stats_group = stats['Reactions and KO sources']
        for key in (
            'Reactions in network',
            'Mean reactions per KO',
            'Stdev reactions per KO',
            'Max reactions per KO'
        ):
            self.run.info(key, stats_group[key])

        self.run.info_single("Reaction alias source comparison", nl_before=1)
        stats_group = stats['Reaction alias sources']
        for key in (
            'Reactions aliased by KEGG reaction',
            'Reactions aliased by EC number',
            'Rxns aliased by both KEGG rxn & EC number',
            'Reactions aliased only by KEGG reaction',
            'Reactions aliased only by EC number',
            'KEGG reactions contributing to network',
            'Mean reactions per KEGG reaction',
            'Stdev reactions per KEGG reaction',
            'Max reactions per KEGG reaction',
            'EC numbers contributing to network',
            'Mean reactions per EC number',
            'Stdev reactions per EC number',
            'Max reactions per EC number'
        ):
            self.run.info(key, stats_group[key])

        stats_group = stats['Reaction and metabolite properties']
        self.run.info_single("Reaction reversibility", nl_before=1)
        for key in (
            'Reversible reactions',
            'Irreversible reactions'
        ):
            self.run.info(key, stats_group[key])

        self.run.info_single("Metabolites and localization", nl_before=1)
        for key in (
            'Metabolites in network',
            'Cytoplasmic metabolites',
            'Extracellular metabolites',
            'Exclusively cytoplasmic metabolites',
            'Exclusively extracellular metabolites',
            'Cytoplasmic/extracellular metabolites'
        ):
            self.run.info(key, stats_group[key])

        self.run.info_single("Metabolite consumption and production", nl_before=1)
        for key in (
            'Consumed metabolites',
            'Produced metabolites',
            'Both consumed & produced metabolites',
            'Exclusively consumed metabolites',
            'Exclusively produced metabolites',
            'Metabolites consumed or produced by 1 rxn',
            'Metabolites consumed or produced by 2 rxns',
            'Metabolites consumed or produced by 3+ rxns'
        ):
            self.run.info(key, stats_group[key])
        print()

    def write_overview_statistics(
        self,
        stats_file: str,
        stats: Union[GenomicNetworkStats, PangenomicNetworkStats] = None
    ) -> None:
        """
        Write a tab-delimited file of overview statistics for the metabolic network.

        Parameters
        ==========
        stats_file : str
            Path to output tab-delimited file of overview statistics.

        stats : Union[GenomicNetworkStats, PangenomicNetworkStats], None
            With the default value of None, network statistics will be calculated and written to
            file. Alternatively, provided network statistics will be written to file without
            calculating anew.

        Returns
        =======
        None
        """
        if not stats:
            # Subclasses must have a method, 'get_overview_statistics'.
            stats = self.get_overview_statistics()

        filesnpaths.is_output_file_writable(stats_file)

        table = []
        for stats_group_name, stats_group in stats.items():
            for stat_name, stat_value in stats_group.items():
                table.append([stats_group_name, stat_name, stat_value])
        pd.DataFrame(table, columns=['Group', 'Statistic', 'Value']).to_csv(
            stats_file, sep='\t', index=False
        )

        self.run.info("Metabolic network statistics output file", stats_file)

class GenomicNetwork(ReactionNetwork):
    """
    A reaction network predicted from KEGG Ortholog annotations of genes and ModelSEED data.

    Attributes
    ==========
    kos : Dict[str, KO], dict()
        This dictionary maps the IDs of KOs in the network to object representations of the KOs.

    reactions : Dict[str, ModelSEEDReaction], dict()
        This maps the IDs of ModelSEED reactions in the network to object representations of the
        reactions.

    metabolites : Dict[str, ModelSEEDCompound], dict()
        This maps the IDs of ModelSEED metabolites in the network to object representations of the
        metabolites.

    kegg_modelseed_aliases : Dict[str, List[str]], dict()
        This maps KEGG REACTION IDs associated with KOs in the network to ModelSEED reactions
        aliased by the KEGG reaction. KO-associated KEGG reactions that do not alias ModelSEED
        reactions are not included.

    ec_number_modelseed_aliases : Dict[str, List[str]], dict()
        This maps EC numbers associated with KOs in the network to ModelSEED reactions aliased by
        the EC number. KO-associated EC numbers that do not alias ModelSEED reactions are not
        included.

    modelseed_kegg_aliases : Dict[str, List[str]], dict()
        This maps the IDs of ModelSEED reactions in the network to lists of KEGG REACTION IDs that
        are associated with KOs in the network and alias the ModelSEED reaction.

    modelseed_ec_number_aliases : Dict[str, List[str]], dict()
        This maps the IDs of ModelSEED reactions in the network to lists of EC numbers that are
        associated with KOs in the network and alias the ModelSEED reaction.

    contigs_db_source_path : str, None
        Path to the contigs database from which the network was built.

    profile_db_source_path : str, None
        Path to the profile database from which protein and metabolite abundance data was loaded.

    genes : Dict[int, Gene], dict()
        This maps gene callers IDs to object representations of genes in the network.

    bins : Dict[str, GeneBin], dict()

    collection : BinCollection, None

    proteins : Dict[int, Protein], dict()
        This maps protein IDs to object representations of proteins with abundance data in the
        network.
    """
    def __init__(
        self,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        verbose: bool = True
    ) -> None:
        """
        Parameters
        ==========
        run : anvio.terminal.Run, anvio.terminal.Run()
            This object sets the 'run' attribute, which prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object sets the 'progress' attribute, which prints transient progress information
            to the terminal.

        verbose : bool, True
            This sets the 'verbose' attribute, causing more information to be reported to the
            terminal if True.

        Returns
        =======
        None
        """
        super().__init__(run=run, progress=progress, verbose=verbose)
        self.contigs_db_source_path: str = None
        self.profile_db_source_path: str = None
        self.genes: Dict[int, Gene] = {}
        self.bins: Dict[str, GeneBin] = {}
        self.collection: BinCollection = None
        self.proteins: Dict[int, Protein] = {}

    def copy(self) -> GenomicNetwork:
        """
        Create a deep copy of the reaction network.

        Returns
        =======
        GenomicNetwork
            Deep copy of the reaction network.
        """
        copied_network = GenomicNetwork()

        self._copy(copied_network)

        for gcid, gene in self.genes.items():
            copied_gene = Gene()
            copied_gene.gcid = gcid
            for ko_id, ko in gene.kos.items():
                copied_gene.kos[ko_id] = ko
            for ko_id, e_value in gene.e_values.items():
                copied_gene.e_values[ko_id] = e_value

            copied_network.genes[gcid] = copied_gene

        if self.proteins:
            for protein_id, protein in self.proteins.items():
                copied_protein = Protein()
                copied_protein.id = protein_id
                for gcid, gene in protein.genes.items():
                    copied_protein.genes[gcid] = gene = copied_network.genes[gcid]
                    gene.protein = copied_protein
                copied_protein.abundances = protein.abundances.copy()

                copied_network.proteins[protein_id] = copied_protein

        return copied_network

    def remove_metabolites_without_formula(self, output_path: str = None) -> None:
        """
        Remove metabolites without a formula in the ModelSEED database from the network.

        Other items can be removed from the network by association: reactions that involve a
        formulaless metabolite; other metabolites with formulas that are exclusive to such
        reactions; KOs predicted to exclusively catalyze such reactions; and genes exclusively
        annotated with such KOs. Removed metabolites with a formula are reported alongside
        formulaless metabolites to the optional output table of removed metabolites.

        output_path : str, None
            If not None, write four tab-delimited files of metabolites, reactions, KEGG Orthologs,
            and genes removed from the network to file locations based on the provided path. For
            example, if the argument, 'removed.tsv', is provided, then the following files will be
            written: 'removed-metabolites.tsv', 'removed-reactions.tsv', 'removed-kos.tsv', and
            'removed-genes.tsv'.
        """
        if self.verbose:
            self.progress.new("Removing metabolites without a formula in the network")
            self.progress.update("...")

        if output_path:
            path_basename, path_extension = os.path.splitext(output_path)
            metabolite_path = f"{path_basename}-metabolites{path_extension}"
            reaction_path = f"{path_basename}-reactions{path_extension}"
            ko_path = f"{path_basename}-kos{path_extension}"
            gene_path = f"{path_basename}-genes{path_extension}"
            for path in (metabolite_path, reaction_path, ko_path, gene_path):
                filesnpaths.is_output_file_writable(path)

        metabolites_to_remove = []
        for modelseed_compound_id, metabolite in self.metabolites.items():
            # ModelSEED compounds without a formula have a formula value of None in the network
            # object.
            if metabolite.formula is None:
                metabolites_to_remove.append(modelseed_compound_id)
        removed = self.purge_metabolites(metabolites_to_remove)

        if self.verbose:
            self.progress.end()
            self.run.info("Removed metabolites", len(removed['metabolite']))
            self.run.info("Removed reactions", len(removed['reaction']))
            self.run.info("Removed KOs", len(removed['ko']))
            self.run.info("Removed genes", len(removed['gene']))

        if not output_path:
            return

        if self.verbose:
            self.progress.new("Writing output files of removed network items")
            self.progress.update("...")

        # Record the reactions removed as a consequence of involving formulaless metabolites, and
        # record the formulaless metabolites involved in removed reactions.
        metabolite_removed_reactions: Dict[str, List[str]] = {}
        reaction_removed_metabolites: Dict[str, List[str]] = {}
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            reaction_removed_metabolites[reaction.modelseed_id] = metabolite_ids = []
            for metabolite in reaction.compounds:
                if metabolite.modelseed_id in metabolites_to_remove:
                    try:
                        metabolite_removed_reactions[metabolite.modelseed_id].append(
                            reaction.modelseed_id
                        )
                    except KeyError:
                        metabolite_removed_reactions[metabolite.modelseed_id] = [
                            reaction.modelseed_id
                        ]
                    metabolite_ids.append(metabolite.modelseed_id)

        metabolite_table = []
        for metabolite in removed['metabolite']:
            metabolite: ModelSEEDCompound
            row = []
            row.append(metabolite.modelseed_id)
            row.append(metabolite.modelseed_name)
            row.append(metabolite.formula)
            try:
                # The metabolite did not have a formula.
                removed_reaction_ids = metabolite_removed_reactions[metabolite.modelseed_id]
            except KeyError:
                # The metabolite had a formula but was removed as a consequence of all the reactions
                # involving the metabolite being removed due to them containing formulaless
                # metabolites: the metabolite did not cause any reactions to be removed.
                row.append("")
                continue
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus the reaction is recorded multiple times.
            row.append(", ".join(sorted(set(removed_reaction_ids))))

        reaction_table = []
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            row = []
            row.append(reaction.modelseed_id)
            row.append(reaction.modelseed_name)
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus is recorded multiple times.
            row.append(
                ", ".join(set(reaction_removed_metabolites[reaction.modelseed_id]))
            )
            row.append(", ".join([metabolite.modelseed_id for metabolite in reaction.compounds]))
            row.append(get_chemical_equation(reaction))
            reaction_table.append(row)

        ko_table = []
        for ko in removed['ko']:
            ko: KO
            row = []
            row.append(ko.id)
            row.append(ko.name)
            row.append(", ".join(ko.reactions))
            ko_table.append(row)

        gene_table = []
        for gene in removed['gene']:
            gene: Gene
            row = []
            row.append(gene.gcid)
            row.append(", ".join(gene.kos))
            gene_table.append(row)

        pd.DataFrame(
            metabolite_table,
            columns=[
                "ModelSEED compound ID",
                "ModelSEED compound name",
                "Formula",
                "Removed reaction ModelSEED IDs"
            ]
        ).to_csv(metabolite_path, sep='\t', index=False)
        pd.DataFrame(
            reaction_table,
            columns=[
                "ModelSEED reaction ID",
                "ModelSEED reaction name",
                "Removed ModelSEED compound IDs",
                "Reaction ModelSEED compound IDs",
                "Equation"
            ]
        ).to_csv(reaction_path, sep='\t', index=False)
        pd.DataFrame(
            ko_table,
            columns=[
                "KO ID",
                "KO name",
                "KO ModelSEED reaction IDs"
            ]
        ).to_csv(ko_path, sep='\t', index=False)
        pd.DataFrame(
            gene_table,
            columns=[
                "Gene callers ID",
                "KO IDs"
            ]
        ).to_csv(gene_path, sep='\t', index=False)

        if self.verbose:
            self.progress.end()
            self.run.info("Table of removed metabolites", metabolite_path)
            self.run.info("Table of removed reactions", reaction_path)
            self.run.info("Table of removed KOs", ko_path)
            self.run.info("Table of removed genes", gene_path)

    def purge_metabolites(self, metabolites_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given metabolites from the network.

        Reactions involving the metabolite are also purged from the network. KOs that were only
        associated with removed reactions are purged; genes that were only associated with removed
        KOs are purged.

        Removal of reactions involving the metabolite can also result in other metabolites being
        being removed from the network, those that exclusively participate in these reactions.

        Parameters
        ==========
        metabolites_to_remove : Iterable[str]
            Metabolites to remove by ModelSEED compound ID.

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
            only contain one significant entry:
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }
        """
        removed_metabolites: List[ModelSEEDCompound] = []
        for modelseed_compound_id in metabolites_to_remove:
            try:
                removed_metabolites.append(self.metabolites.pop(modelseed_compound_id))
            except KeyError:
                # This can occur for two reasons. First, the metabolite from 'metabolites_to_remove'
                # could not be in the network.

                # Second, this can occur when removing other "unintended" metabolites from the
                # network. 'purge_metabolites' was first called with metabolites of interest, then
                # 'purge_reactions' was called from within the method the remove reactions involving
                # the metabolites of interest, and then 'purge_metabolites' was called again from
                # within 'purge_reactions' to remove other metabolites exclusively found in the
                # removed reactions. In this last call of 'purge_metabolites', the
                # 'metabolites_to_remove' also include the metabolites of interest that were already
                # removed from 'self.metabolites' in the original 'purge_metabolites' call. This
                # KeyError occurs when trying to remove those already-removed metabolites.
                pass
        if not removed_metabolites:
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }

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
            removed_cascading_up = {
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }
        removed.update(removed_cascading_up)
        return removed

    def purge_reactions(self, reactions_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given reactions from the network.

        Metabolites that exclusively participate in removed reactions are purged. KOs that were only
        associated with removed reactions are purged; genes that were only associated with removed
        KOs are purged.

        Parameters
        ==========
        reactions_to_remove : Iterable[str]
            Reactions to remove by ModelSEED reaction ID.

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
                # Alternatively, this occurs if the reaction in 'reactions_to_remove' is not in the
                # network.
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
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }

        # Remove KEGG reaction aliases of ModelSEED reactions in the network.
        kegg_reactions_to_remove = []
        for kegg_reaction_id, modelseed_reaction_ids in self.kegg_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the KEGG reaction were removed, so remove the
                # KEGG reaction as well.
                kegg_reactions_to_remove.append(kegg_reaction_id)
                continue
            for idx in sorted(aliases_to_remove, reverse=True):
                # Only some of the ModelSEED reactions aliased by the KO KEGG reaction were removed.
                modelseed_reaction_ids.pop(idx)
        removed_kegg_reactions = []
        for kegg_reaction_id in kegg_reactions_to_remove:
            removed_kegg_reactions.append(self.kegg_modelseed_aliases.pop(kegg_reaction_id))

        # Remove EC number aliases of ModelSEED reactions in the network.
        ec_numbers_to_remove = []
        for ec_number, modelseed_reaction_ids in self.ec_number_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the EC number were removed, so remove the EC
                # number as well.
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
                        # Do not remove the metabolite, because it participates in a retained
                        # reaction.
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
        # exclusively associated with reactions removed here, then 'kos_to_remove' would likewise be
        # empty.
        if kos_to_remove:
            removed_cascading_up = self.purge_kos(kos_to_remove)
            removed_cascading_up.pop('reaction')
            removed_cascading_up.pop('kegg_reaction')
            removed_cascading_up.pop('ec_number')
        else:
            removed_cascading_up = {'ko': [], 'gene': []}

        removed = {
            'reaction': removed_reactions,
            'kegg_reaction': removed_kegg_reactions,
            'ec_number': removed_ec_numbers
        }
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def purge_kos(self, kos_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given KOs from the network.

        Reactions and metabolites that were only associated with removed KOs are purged. Genes that
        were only associated with removed KOs are purged.

        Parameters
        ==========
        kos_to_remove : Iterable[str]
            KOs to remove by ID.

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
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }

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
                        # The reaction is associated with a retained KO, so do not remove the
                        # reaction.
                        reactions_to_spare.append(idx)
            for idx in sorted(reactions_to_spare, reverse=True):
                reactions_to_remove.pop(idx)
        if reactions_to_remove:
            removed_cascading_down = self.purge_reactions(reactions_to_remove)
            removed_cascading_down.pop('ko')
        else:
            # This method must have been called from the method, 'purge_reactions', because the
            # reactions that are only associated with the removed KOs were already removed from the
            # network.
            removed_cascading_down = {
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        genes_to_remove: List[str] = []
        for gcid, gene in self.genes.items():
            gene_kos_to_remove: List[str] = []
            for ko_id in gene.kos:
                if ko_id in kos_to_remove:
                    gene_kos_to_remove.append(ko_id)
            if len(gene_kos_to_remove) == len(gene.kos):
                # All KOs matching the gene were removed, so remove it as well.
                genes_to_remove.append(gcid)
                continue
            for ko_id in gene_kos_to_remove:
                gene.kos.pop(ko_id)
                gene.e_values.pop(ko_id)
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

    def purge_genes(self, genes_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given genes from the network.

        KOs, reactions, and metabolites that were only associated with removed genes are purged.

        Parameters
        ==========
        genes_to_remove : Iterable[str]
            Genes to remove by gene callers ID.

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
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene': []
            }

        kos_to_remove: List[str] = []
        for gene in removed_genes:
            for ko_id in gene.kos:
                kos_to_remove.append(ko_id)
        kos_to_remove = list(set(kos_to_remove))
        for gene in self.genes.values():
            kos_to_spare: List[str] = []
            for ko_id in gene.kos:
                if ko_id in kos_to_remove:
                    # The KO is associated with a retained gene, so do not remove the KO.
                    kos_to_spare.append(ko_id)
            for ko_id in kos_to_spare:
                kos_to_remove.remove(ko_id)
        if kos_to_remove:
            removed_cascading_down = self.purge_kos(kos_to_remove)
            removed_cascading_down.pop('gene')
        else:
            # This method must have been called from the method, 'purge_kos', because the KOs that
            # are only associated with the removed genes were already removed from the network.
            removed_cascading_down = {
                'ko': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        # TODO: remove genes from self.bins

        removed = {'gene': removed_genes}
        removed.update(removed_cascading_down)
        return removed

    def subset_network(
        self,
        kegg_modules_to_subset: Iterable[str] = None,
        brite_categories_to_subset: Iterable[str] = None,
        genes_to_subset: Iterable[int] = None,
        kos_to_subset: Iterable[str] = None,
        reactions_to_subset: Iterable[str] = None,
        metabolites_to_subset: Iterable[str] = None
    ) -> GenomicNetwork:
        """
        Subset a smaller network from the metabolic network.

        If requested KEGG modules, BRITE categories, genes, KOs, reactions, or metabolites are not
        present in the network, no error is raised.

        Subsetted items are not represented by the same objects as in the source network, i.e., new
        gene, KO, reaction, and metabolite objects are created and added to the subsetted network.

        Network items (i.e., genes, KOs, reactions, and metabolites) that reference requested items
        (e.g., genes in the network referencing requested KOs; KOs in the network referencing
        requested reactions) are added to the subsetted network. A gene added to the subsetted
        network due to references to requested KOs will be missing references to any other
        "unrequested" KOs annotating the gene in the source network. Likewise, genes and KOs that
        are added to the subsetted network due to references to requested reactions will be missing
        references to any other unrequested reactions. In other words, certain KO and reaction
        annotations can be selected to the exclusion of others, e.g., a KO encoding two reactions
        can be "redefined" or "pruned" to encode one requested reaction in the subsetted network; a
        KO encoding multiple reactions can be pruned to encode only those reactions involving
        requested metabolites.

        If the 'verbose' attribute of the source 'GenomicNetwork' object is True, then report to the
        terminal the identities of requested KEGG modules, BRITE categories, genes, KOs, reactions,
        and metabolites that are not present in the network.

        Parameters
        ==========
        kegg_modules_to_subset : Iterable[str], None
            KEGG module IDs to subset.

        brite_categories_to_subset : Iterable[str], None
            BRITE categories to subset.

        genes_to_subset : Iterable[int], None
            Gene callers IDs to subset.

        kos_to_subset : Iterable[str], None
            KO IDs to subset.

        reactions_to_subset : Iterable[str], None
            ModelSEED reaction IDs to subset.

        metabolites_to_subset : Iterable[str], None
            ModelSEED metabolite IDs to subset.

        Returns
        =======
        GenomicNetwork
            New subsetted reaction network.
        """
        # Sequentially subset the network for each type of request. Upon generating two subsetted
        # networks from two types of request, merge the networks into a single subsetted network;
        # repeat.
        first_subnetwork = None
        for items_to_subset, subset_network_method in (
            (kegg_modules_to_subset, self._subset_network_by_modules),
            (brite_categories_to_subset, self._subset_network_by_brite),
            (genes_to_subset, self._subset_network_by_genes),
            (kos_to_subset, self._subset_network_by_kos),
            (reactions_to_subset, self._subset_network_by_reactions),
            (metabolites_to_subset, self._subset_network_by_metabolites)
        ):
            if not items_to_subset:
                continue

            second_subnetwork = subset_network_method(items_to_subset)

            if first_subnetwork is None:
                first_subnetwork = second_subnetwork
            else:
                first_subnetwork = first_subnetwork.merge_network(second_subnetwork)

        return first_subnetwork

    def _subset_network_by_modules(self, kegg_modules: Iterable[str]) -> GenomicNetwork:
        """
        Subset the network by KOs in requested KEGG modules.

        Parameters
        ==========
        kegg_modules : Iterable[str]
            KEGG modules (of KOs) to subset by ID.

        Returns
        =======
        GenomicNetwork
            New subsetted reaction network.
        """
        pass

    def _subset_network_by_brite(self, brite_categories: Iterable[str]) -> GenomicNetwork:
        """
        Subset the network by KOs in requested KEGG BRITE hierarchy categories.

        Parameters
        ==========
        brite_categories : Iterable[str]
            KEGG BRITE hierarchy categories (of KOs) to subset.

        Returns
        =======
        GenomicNetwork
            New subsetted reaction network.
        """
        pass

    def _subset_network_by_genes(self, gcids: Iterable[int]) -> GenomicNetwork:
        """
        Subset the network by genes with requested gene callers IDs.

        Parameters
        ==========
        gcids : Iterable[int]
            Genes to subset by gene callers ID.

        Returns
        =======
        GenomicNetwork
            New subsetted reaction network.
        """
        subnetwork = GenomicNetwork()

        for gcid in gcids:
            try:
                gene = self.genes[gcid]
            except KeyError:
                # This occurs if the requested gene callers ID is not in the source network.
                continue

            subsetted_gene = Gene()
            subsetted_gene.gcid = gcid
            subsetted_gene.e_values = gene.e_values.copy()

            # Add KOs annotating the gene to the subsetted network as new objects, and then
            # reference these objects in the gene object.
            ko_ids = list(gene.kos)
            subnetwork = self._subset_network_by_kos(ko_ids, subnetwork=subnetwork)
            subsetted_gene.kos = {ko_id: subnetwork.kos[ko_id] for ko_id in ko_ids}

            subnetwork.genes[gcid] = subsetted_gene
        self._subset_proteins(subnetwork)

        return subnetwork

    def _subset_proteins(self, subnetwork: GenomicNetwork) -> None:
        """
        Add protein abundance data to the subsetted network.

        Parameters
        ==========
        subnetwork : GenomicNetwork
            The subsetted reaction network under construction.

        Returns
        =======
        None
        """
        if not self.proteins:
            # Protein abundance profile data is not present in the source network.
            return subnetwork

        # Parse each protein with abundance data.
        for protein_id, protein in self.proteins.items():
            subsetted_gcids: List[int] = []
            for gcid in protein.genes:
                if gcid in subnetwork.genes:
                    # A subsetted gene encodes the protein.
                    subsetted_gcids.append(gcid)
            if not subsetted_gcids:
                # No genes expressing the protein were subsetted, so the protein data is not added.
                continue

            subsetted_protein = Protein()
            subsetted_protein.id = protein_id
            subsetted_protein.abundances = protein.abundances.copy()
            for gcid in subsetted_gcids:
                subsetted_gene = subnetwork.genes[gcid]
                subsetted_gene.protein = subsetted_protein
                subsetted_protein.genes[gcid] = subsetted_gene

            subnetwork.proteins[protein_id] = subsetted_protein

    def _subset_network_by_kos(
        self,
        ko_ids: Iterable[str],
        subnetwork: GenomicNetwork = None
    ) -> GenomicNetwork:
        """
        Subset the network by KOs with requested KO IDs.

        Parameters
        ==========
        ko_ids : Iterable[str]
            KOs to subset by KO ID.

        subnetwork : GenomicNetwork, None
            This network under construction is provided when the KOs being added to the network
            annotate already subsetted genes.

        Returns
        =======
        GenomicNetwork
            If a 'subnetwork' argument is provided, then that network is returned after
            modification. Otherwise, a new subsetted reaction network is returned.
        """
        if subnetwork is None:
            subnetwork = GenomicNetwork()
            # Signify that genes annotated by subsetted KOs are to be added to the network.
            subset_referencing_genes = True
        else:
            assert isinstance(subnetwork, GenomicNetwork)
            # Signify that the KOs being added to the network annotate subsetted genes that were
            # already added to the network.
            subset_referencing_genes = False

        for ko_id in ko_ids:
            try:
                ko = self.kos[ko_id]
            except KeyError:
                # This occurs if the requested KO ID is not in the source network.
                continue

            subsetted_ko = KO()
            subsetted_ko.id = ko.id
            subsetted_ko.name = ko.name
            subsetted_ko.kegg_reaction_aliases = deepcopy(ko.kegg_reaction_aliases)
            subsetted_ko.ec_number_aliases = deepcopy(ko.ec_number_aliases)

            # Add reactions annotating the KO to the subsetted network as new objects, and then
            # reference these objects in the KO object.
            reaction_ids = [reaction_id for reaction_id in ko.reactions]
            subnetwork = self._subset_network_by_reactions(reaction_ids, subnetwork=subnetwork)
            subsetted_ko.reactions = {
                reaction_id: subnetwork.reactions[reaction_id] for reaction_id in reaction_ids
            }

            subnetwork.kos[ko_id] = subsetted_ko

        if subset_referencing_genes:
            # Add genes that are annotated by the subsetted KOs to the network.
            self._subset_genes_via_kos(subnetwork)

        return subnetwork

    def _subset_network_by_reactions(
        self,
        reaction_ids: Iterable[str],
        subnetwork: GenomicNetwork = None
    ) -> GenomicNetwork:
        """
        Subset the network by reactions with ModelSEED reaction IDs.

        Parameters
        ==========
        reaction_ids : Iterable[str]
            Reactions to subset by ModelSEED reaction ID.

        subnetwork : GenomicNetwork, None
            This network under construction is provided when the reactions being added to the
            network annotate already subsetted KOs.

        Returns
        =======
        GenomicNetwork
            If a 'subnetwork' argument is provided, then that network is returned after
            modification. Otherwise, a new subsetted reaction network is returned.
        """
        if subnetwork is None:
            subnetwork = GenomicNetwork()
            # Signify that KOs annotated by subsetted reactions are to be added to the network.
            subset_referencing_kos = True
        else:
            assert isinstance(subnetwork, GenomicNetwork)
            # Signify that the reactions being added to the network annotate subsetted KOs that were
            # already added to the network.
            subset_referencing_kos = False

        # Copy the network attributes mapping reaction aliases.
        kegg_modelseed_aliases: Dict[str, List[str]] = {}
        ec_number_modelseed_aliases: Dict[str, List[str]] = {}

        for reaction_id in reaction_ids:
            try:
                reaction = self.reactions[reaction_id]
            except KeyError:
                # This occurs if the requested reaction is not in the source network.
                continue

            # Copy the reaction object, including referenced metabolite objects, from the source
            # network.
            subsetted_reaction: ModelSEEDReaction = deepcopy(reaction)
            subnetwork.reactions[reaction_id] = subsetted_reaction
            # Record the metabolites involved in the reaction, and add them to the network.
            for metabolite in subsetted_reaction.compounds:
                compound_id = metabolite.modelseed_id
                subnetwork.metabolites[compound_id] = metabolite

            try:
                subnetwork.modelseed_kegg_aliases[reaction_id] += list(reaction.kegg_aliases)
            except KeyError:
                subnetwork.modelseed_kegg_aliases[reaction_id] = list(reaction.kegg_aliases)

            try:
                subnetwork.modelseed_ec_number_aliases[reaction_id] += list(
                    reaction.ec_number_aliases
                )
            except KeyError:
                subnetwork.modelseed_ec_number_aliases[reaction_id] = list(
                    reaction.ec_number_aliases
                )

            for kegg_id in reaction.kegg_aliases:
                try:
                    kegg_modelseed_aliases[kegg_id].append(reaction_id)
                except KeyError:
                    kegg_modelseed_aliases[kegg_id] = [reaction_id]

            for ec_number in reaction.ec_number_aliases:
                try:
                    ec_number_modelseed_aliases[ec_number].append(reaction_id)
                except KeyError:
                    ec_number_modelseed_aliases[ec_number] = [reaction_id]

        if subnetwork.kegg_modelseed_aliases:
            for kegg_id, modelseed_ids in kegg_modelseed_aliases.items():
                try:
                    subnetwork.kegg_modelseed_aliases[kegg_id] += modelseed_ids
                except KeyError:
                    subnetwork.kegg_modelseed_aliases[kegg_id] = modelseed_ids
        else:
            subnetwork.kegg_modelseed_aliases = kegg_modelseed_aliases

        if subnetwork.ec_number_modelseed_aliases:
            for ec_number, modelseed_ids in ec_number_modelseed_aliases.items():
                try:
                    subnetwork.ec_number_modelseed_aliases[ec_number] += modelseed_ids
                except KeyError:
                    subnetwork.ec_number_modelseed_aliases[ec_number] = modelseed_ids
        else:
            subnetwork.ec_number_modelseed_aliases = ec_number_modelseed_aliases

        if subset_referencing_kos:
            # Add KOs that are annotated by the subsetted reactions to the network.
            self._subset_kos_via_reactions(subnetwork)

        return subnetwork

    def _subset_genes_via_kos(self, subnetwork: GenomicNetwork) -> None:
        """
        Add genes that are annotated with subsetted KOs to the subsetted network.

        These gene objects only reference subsetted KOs and not other KOs that also annotate the
        gene but are not subsetted.

        Parameters
        ==========
        subnetwork : GenomicNetwork
            The subsetted reaction network under construction.

        Returns
        =======
        None
        """
        subsetted_ko_ids = list(subnetwork.kos)
        for gcid, gene in self.genes.items():
            # Check all genes in the source network for subsetted KOs.
            subsetted_gene = None
            for ko_id in gene.kos:
                if ko_id not in subsetted_ko_ids:
                    # The gene is not annotated by the subsetted KO.
                    continue

                if not subsetted_gene:
                    # Create a new gene object for the subsetted gene. The gene object would already
                    # have been created had another subsetted KO been among the KOs annotating the
                    # gene.
                    subsetted_gene = Gene()
                    subsetted_gene.gcid = gcid
                subsetted_gene.kos[ko_id] = subnetwork.kos[ko_id]
                subsetted_gene.e_values[ko_id] = gene.e_values[ko_id]

            if subsetted_gene:
                subnetwork.genes[gcid] = subsetted_gene

    def _subset_kos_via_reactions(self, subnetwork: GenomicNetwork) -> None:
        """
        Add KOs that are annotated with subsetted reactions to the subsetted network.

        Then add genes that are annotated with these added KOs to the subsetted network.

        Parameters
        ==========
        subnetwork : GenomicNetwork
            The subsetted reaction network under construction.

        Returns
        =======
        None
        """
        subsetted_reaction_ids = list(subnetwork.reactions)
        for ko_id, ko in self.kos.items():
            # Check all KOs in the source network for subsetted reactions.
            subsetted_ko = None
            for reaction_id in ko.reactions:
                if reaction_id not in subsetted_reaction_ids:
                    # The KO is not annotated by the subsetted reaction.
                    continue

                if not subsetted_ko:
                    # Create a new KO object for the subsetted KO. The subsetted KO object would
                    # already have been created had another subsetted reaction been among the
                    # reactions annotating the KO.
                    subsetted_ko = KO()
                    subsetted_ko.id = ko_id
                    subsetted_ko.name = ko.name
                subsetted_ko.reactions[reaction_id] = subnetwork.reactions[reaction_id]
                subsetted_ko.kegg_reaction_aliases = deepcopy(ko.kegg_reaction_aliases)
                subsetted_ko.ec_number_aliases = deepcopy(ko.ec_number_aliases)

            if subsetted_ko:
                subnetwork.kos[ko_id] = subsetted_ko

        # Add genes that are annotated with the added KOs to the subsetted network.
        self._subset_genes_via_kos(subnetwork)

    def _subset_network_by_metabolites(self, compound_ids: Iterable[str]) -> GenomicNetwork:
        """
        Subset the network by metabolites with ModelSEED compound IDs.

        Parameters
        ==========
        compound_ids : Iterable[str]
            Metabolites to subset by ModelSEED compound ID.

        Returns
        =======
        GenomicNetwork
            New subsetted reaction network.
        """
        subnetwork = GenomicNetwork()

        for reaction_id, reaction in self.reactions.items():
            # Check all reactions in the source network for subsetted metabolites.
            for metabolite in reaction.compounds:
                if metabolite.modelseed_id in compound_ids:
                    break
            else:
                # The reaction does not involve any of the requested metabolites.
                continue

            # Copy the reaction object, including referenced metabolite objects, from the source
            # network.
            subsetted_reaction: ModelSEEDReaction = deepcopy(reaction)
            subnetwork.reactions[reaction_id] = subsetted_reaction

            # Add the metabolites involved in the reaction to the subsetted network. (There can be
            # unavoidable redundancy here in readding previously encountered metabolites.)
            for subsetted_metabolite in subsetted_reaction.compounds:
                subnetwork.metabolites[subsetted_metabolite.modelseed_id] = subsetted_metabolite

        # Add KOs that are annotated with the added reactions to the subsetted network, and then add
        # genes annotated with the added KOs to the subsetted network.
        self._subset_kos_via_reactions(subnetwork)

        return subnetwork

    def merge_network(self, network: GenomicNetwork) -> GenomicNetwork:
        """
        Merge the genomic reaction network with another genomic reaction network.

        Each network can contain different genes, KOs, and reactions/metabolites. Merging
        nonredundantly incorporates all of this data as new objects in the new network.

        Objects representing genes or KOs in both networks can have different sets of references:
        genes can be annotated by different KOs, and KOs can be annotated by different reactions.

        Otherwise, object attributes should be consistent between the networks. For instance, the
        same ModelSEED reactions and metabolites in both networks should have identical attributes.
        If applicable, both networks should have been annotated with the same protein and metabolite
        abundance data.

        The purpose of this method is to combine different, but potentially overlapping, subnetworks
        from the same pangenome.

        Parameters
        ==========
        network : GenomicNetwork
            The other genomic reaction network being merged.

        Returns
        =======
        GenomicNetwork
            The merged genomic reaction network.
        """
        assert not (
            (self.proteins is None and network.proteins is not None) and
            (self.proteins is not None and network.proteins is None)
        )

        merged_network = GenomicNetwork()

        self._merge_network(network, merged_network)

        # Add genes to the merged network, first adding genes present in both source networks, and
        # then adding genes present exclusively in each source network.
        first_gcids = set(self.genes)
        second_gcids = set(network.genes)

        for gcid in first_gcids.intersection(second_gcids):
            first_gene = self.genes[gcid]
            second_gene = network.genes[gcid]

            # The new object representing the gene in the merged network should have all KO
            # annotations from each source gene object, as these objects can have different KO
            # references.
            merged_gene = Gene()
            merged_gene.gcid = gcid
            ko_ids = set(first_gene.kos).union(set(second_gene.kos))
            for ko_id in ko_ids:
                merged_gene.kos[ko_id] = merged_network.kos[ko_id]
            first_ko_ids = set(first_gene.kos)
            second_ko_ids = set(second_gene.kos).difference(set(first_gene.kos))
            for ko_id in first_ko_ids:
                merged_gene.e_values[ko_id] = first_gene.e_values[ko_id]
            for ko_id in second_ko_ids:
                merged_gene.e_values[ko_id] = second_gene.e_values[ko_id]

            merged_network.genes[gcid] = merged_gene

        for gcid in first_gcids.difference(second_gcids):
            first_gene = self.genes[gcid]

            gene = Gene()
            gene.gcid = gcid
            gene.kos = {ko_id: merged_network.kos[ko_id] for ko_id in first_gene.kos}
            gene.e_values = first_gene.e_values.copy()

            merged_network.genes[gcid] = gene

        for gcid in second_gcids.difference(first_gcids):
            second_gene = network.genes[gcid]

            gene = Gene()
            gene.gcid = gcid
            gene.kos = {ko_id: merged_network.kos[ko_id] for ko_id in second_gene.kos}
            gene.e_values = second_gene.e_values.copy()

            merged_network.genes[gcid] = gene

        if not self.proteins and not network.proteins:
            # No protein abundance data is present and needs to be added to the merged network.
            return merged_network

        # Add protein abundance data to the merged network, first adding proteins annotating genes
        # in both source networks, and then adding proteins annotating genes exclusively in each
        # source network.
        first_protein_ids = set(self.proteins)
        second_protein_ids = set(network.proteins)

        # Assume that each source network was annotated with the same protein annotation data, so
        # that the same gene in each network should have the same protein abundance profile.
        for protein_id in first_protein_ids.intersection(second_protein_ids):
            first_protein = self.proteins[protein_id]
            second_protein = network.proteins[protein_id]

            merged_protein = Protein()
            merged_protein.id = protein_id
            for gcid in first_protein.genes:
                merged_protein.genes[gcid] = merged_network.genes[gcid]
            for gcid in set(second_protein.genes).difference(set(first_protein.genes)):
                merged_protein.genes[gcid] = merged_network.genes[gcid]
            merged_protein.abundances = first_protein.abundances.copy()

            merged_network.proteins[protein_id] = merged_protein

        for protein_id in first_protein_ids.difference(second_protein_ids):
            first_protein = self.proteins[protein_id]

            protein = Protein()
            protein.id = protein_id
            protein.genes = {gcid: merged_network.genes[gcid] for gcid in first_protein.genes}
            protein.abundances = first_protein.abundances.copy()

            merged_network.proteins[protein_id] = protein

        for protein_id in second_protein_ids.difference(first_protein_ids):
            second_protein = network.proteins[protein_id]

            protein = Protein()
            protein.id = protein_id
            protein.genes = {gcid: merged_network.genes[gcid] for gcid in second_protein.genes}
            protein.abundances = second_protein.abundances.copy()

            merged_network.proteins[protein_id] = protein

        return merged_network

    # def export_table(
    #     self,
    #     path: str
    # ) -> None:
    #     d = {}
    #     for ko_id, ko in self.kos.items():
    #         d[ko_id] = e = {}
    #         e['ko_id'] = ko_id
    #         e['ko_name'] = ko.name
    #         e['gene_count'] = 0
    #         e['gene_ids'] = ''
    #         e['protein_id'] = ''
    #         e['protein_abundance'] = 0
    #         e['reaction_ids'] = ''
    #         for m in self.metabolites.values():
    #             if m.abundances:
    #                 e[m.modelseed_name] = np.nan
    #     for gcid, gene in self.genes.items():
    #         if gene.protein:
    #             protein_id = str(gene.protein.id)
    #             mean_protein_abund = np.mean(list(gene.protein.abundances.values()))
    #         else:
    #             protein_id = ''
    #             mean_protein_abund = 0
    #         for ko_id, ko in gene.kos.items():
    #             e = d[ko_id]
    #             e['gene_count'] += 1
    #             e['gene_ids'] += f', {gcid}'
    #             e['protein_id'] = protein_id
    #             e['protein_abundance'] = mean_protein_abund
    #     for ko_id, ko in self.kos.items():
    #         e = d[ko_id]
    #         for reaction_id, reaction in ko.reactions.items():
    #             e['reaction_ids'] += f'{reaction_id}, '
    #             for metabolite in reaction.compounds:
    #                 if metabolite.abundances:
    #                     e[metabolite.modelseed_name] = np.mean(list(metabolite.abundances.values()))
    #     table = []
    #     for e in d.values():
    #         e: dict
    #         e['reaction_ids'] = e['reaction_ids'].strip(' ').strip(',')
    #         e['gene_ids'] = e['gene_ids'].strip(' ').strip(',')
    #         table.append([v for v in e.values()])
    #     pd.DataFrame(table, columns=[key for key in e]).to_csv(path, sep='\t', index=False)

    def get_overview_statistics(
        self,
        precomputed_counts: Dict[str, int] = None
    ) -> GenomicNetworkStats:
        """
        Calculate overview statistics for the genomic metabolic network.

        Parameters
        ==========
        precomputed_counts : Dict[str, int], None
            To spare additional computations that involve loading and parsing the contigs database,
            this dictionary can contain two pieces of precomputed data: the value for the key,
            'total_genes', should be the number of genes in the genome; the value for the key,
            'genes_assigned_kos', should be the number of genes in the genome assigned KOs; the
            value for the key, 'kos_assigned_genes', should be the number of unique KOs assigned to
            genes in the genome.

        Returns
        =======
        GenomicNetworkStats
            Network statistics are stored in a dictionary of dictionaries. Keys in the outer
            dictionary are "classes" of network statistics. Keys in the inner dictionary are
            statistics themselves.
        """
        if (
            precomputed_counts is not None and
            sorted(precomputed_counts) != [
                    'genes_assigned_kos', 'kos_assigned_genes', 'total_genes'
            ]
        ):
            raise ConfigError(
                "The 'precomputed_counts' argument must be a dictionary only containing the keys, "
                "'total_genes', 'genes_assigned_kos', and 'kos_assigned_genes'."
            )

        stats: GenomicNetworkStats = {}

        self.progress.new("Counting genes and KOs")
        self.progress.update("...")
        stats['Gene and KO counts'] = stats_group = {}

        if precomputed_counts:
            assert (
                type(precomputed_counts['total_genes']) is int and
                precomputed_counts['total_genes'] >= 0
            )
            gene_count = precomputed_counts['total_genes']
            assert (
                type(precomputed_counts['genes_assigned_kos']) is int and
                precomputed_counts['genes_assigned_kos'] >= 0
            )
            ko_annotated_gene_count = precomputed_counts['genes_assigned_kos']
            assert (
                type(precomputed_counts['kos_assigned_genes']) is int and
                precomputed_counts['kos_assigned_genes'] >= 0
            )
            annotating_ko_count = precomputed_counts['kos_assigned_genes']
        else:
            if self.contigs_db_source_path:
                cdb = ContigsDatabase(self.contigs_db_source_path)
                gene_count = cdb.db.get_row_counts_from_table('genes_in_contigs')
                gene_ko_id_table = cdb.db.get_table_as_dataframe(
                    'gene_functions',
                    where_clause='source = "KOfam"',
                    columns_of_interest=['gene_callers_id', 'source']
                )
                ko_annotated_gene_count = gene_ko_id_table['gene_callers_id'].nunique()
                annotating_ko_count = gene_ko_id_table['KOfam'].nunique()
                cdb.disconnect()
            else:
                gene_count = None
                ko_annotated_gene_count = None
                annotating_ko_count = None

        if gene_count is not None:
            stats_group['Total gene calls in genome'] = gene_count
        if ko_annotated_gene_count is not None:
            stats_group['Genes annotated with protein KOs'] = ko_annotated_gene_count
        stats_group['Genes in network'] = len(self.genes)
        if annotating_ko_count is not None:
            stats_group['Protein KOs annotating genes'] = annotating_ko_count
        stats_group['KOs in network'] = len(self.kos)
        self.progress.end()

        self._get_common_overview_statistics(stats)

        if precomputed_counts:
            return stats

        if not self.contigs_db_source_path:
            self.run.info_single(
                f"""\
                Since the genomic network was not associated with a contigs database, the following
                statistics could not be calculated and were not reported to the output file:
                'Total gene calls in genome', 'Genes annotated with protein KOs', and 'Protein KOs
                annotating genes'.\
                """
            )

        return stats

    def print_overview_statistics(self, stats: GenomicNetworkStats = None) -> None:
        """
        Print overview statistics for the genomic metabolic network.

        Parameters
        ==========
        stats : GenomicNetworkStats, None
            With the default value of None, network statistics will be calculated and printed.
            Alternatively, provided network statistics will be printed without calculating anew.

        Returns
        =======
        None
        """
        if not stats:
            stats = self.get_overview_statistics()

        self.run.info_single("METABOLIC REACTION NETWORK STATISTICS", mc='green', nl_after=1)

        self.run.info_single("Gene calls and KEGG Ortholog (KO) annotations")
        stats_group = stats['Gene and KO counts']
        self.run.info("Total gene calls in genome", stats_group['Total gene calls in genome'])
        self.run.info(
            "Genes annotated with protein KOs", stats_group['Genes annotated with protein KOs']
        )
        self.run.info("Genes in network", stats_group['Genes in network'])
        self.run.info("Protein KOs annotating genes", stats_group['Protein KOs annotating genes'])
        self.run.info("KOs in network", stats_group['KOs in network'], nl_after=1)

        self._print_common_overview_statistics(stats)

    def export_json(
        self,
        path: str,
        overwrite: bool = False,
        objective: str = None,
        remove_missing_objective_metabolites: bool = False,
        # record_bins: Tuple[str] = ('gene', ),
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

            None means that no objective is added to the JSON, meaning that FBA cannot be performed
            on the model.

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
            in which bins the reaction and metabolite entries occur. To not record bins at all, pass
            either an empty tuple or None.

        indent : int, 2
            spaces of indentation per nesting level in JSON file

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
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
                self.remove_missing_objective_metabolites(objective_dict)
            json_reactions.append(objective_dict)
        elif objective != None:
            raise ConfigError(f"Anvi'o does not recognize an objective with the name, '{objective}'.")

        progress.update("Genes")
        reaction_genes: Dict[str, List[str]] = {}
        reaction_kos: Dict[str, List[KO]] = {}
        for gcid, gene in self.genes.items():
            gene_entry = JSONStructure.get_gene_entry()
            json_genes.append(gene_entry)
            gcid_str = str(gcid)
            gene_entry['id'] = gcid_str
            # Record KO IDs and annotation e-values in the annotation section of the gene entry.
            annotation = gene_entry['annotation']
            annotation['ko'] = annotation_kos = {}
            for ko_id, ko in gene.kos.items():
                annotation_kos[ko_id] = str(gene.e_values[ko_id])
                for modelseed_reaction_id in ko.reactions:
                    try:
                        reaction_genes[modelseed_reaction_id].append(gcid_str)
                    except KeyError:
                        reaction_genes[modelseed_reaction_id] = [gcid_str]
                    try:
                        reaction_kos[modelseed_reaction_id].append(ko)
                    except KeyError:
                        reaction_kos[modelseed_reaction_id] = [ko]

        progress.update("Reactions")
        compound_compartments: Dict[str, Set[str]] = {}
        for modelseed_reaction_id, reaction in self.reactions.items():
            reaction_entry = JSONStructure.get_reaction_entry()
            json_reactions.append(reaction_entry)
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

        progress.update("Metabolites")
        for modelseed_compound_id, metabolite in self.metabolites.items():
            modelseed_compound_name = metabolite.modelseed_name
            charge = metabolite.charge
            formula = metabolite.formula
            kegg_compound_aliases = list(metabolite.kegg_aliases)
            for compartment in compound_compartments[modelseed_compound_id]:
                metabolite_entry = JSONStructure.get_metabolite_entry()
                json_metabolites.append(metabolite_entry)
                metabolite_entry['id'] = f"{modelseed_compound_id}_{compartment}"
                metabolite_entry['name'] = modelseed_compound_name
                metabolite_entry['compartment'] = compartment
                # Compounds without a formula have a nominal charge of 10000000 in the ModelSEED
                # compounds database, which is replaced by None in the reaction network and 0 in the JSON.
                metabolite_entry['charge'] = charge if charge is not None else 0
                metabolite_entry['formula'] = formula if formula is not None else ""
                metabolite_entry['annotation']['kegg.compound'] = kegg_compound_aliases

        progress.update("Saving")
        with open(path, 'w') as f:
            json.dump(json_dict, f, indent=indent)
        progress.end()

class PangenomicNetwork(ReactionNetwork):
    """
    A reaction network predicted from KEGG KO and ModelSEED annotations of pangenomic gene clusters.

    Attributes
    ==========
    kos : Dict[str, KO], dict()
        This dictionary maps the IDs of KOs in the network to object representations of the KOs.

    reactions : Dict[str, ModelSEEDReaction], dict()
        This maps the IDs of ModelSEED reactions in the network to object representations of the
        reactions.

    metabolites : Dict[str, ModelSEEDCompound], dict()
        This maps the IDs of ModelSEED metabolites in the network to object representations of the
        metabolites.

    kegg_modelseed_aliases : Dict[str, List[str]], dict()
        This maps KEGG REACTION IDs associated with KOs in the network to ModelSEED reactions
        aliased by the KEGG reaction. KO-associated KEGG reactions that do not alias ModelSEED
        reactions are not included.

    ec_number_modelseed_aliases : Dict[str, List[str]], dict()
        This maps EC numbers associated with KOs in the network to ModelSEED reactions aliased by
        the EC number. KO-associated EC numbers that do not alias ModelSEED reactions are not
        included.

    modelseed_kegg_aliases : Dict[str, List[str]], dict()
        This maps the IDs of ModelSEED reactions in the network to lists of KEGG REACTION IDs that
        are associated with KOs in the network and alias the ModelSEED reaction.

    modelseed_ec_number_aliases : Dict[str, List[str]], dict()
        This maps the IDs of ModelSEED reactions in the network to lists of EC numbers that are
        associated with KOs in the network and alias the ModelSEED reaction.

    pan_db_source_path : str, None
        Path to the pan database from which the network was built.

    genomes_storage_db_source_path : str, None
        Path to the genomes storage database from which the network was built.

    consensus_threshold : float, None
        A parameter used in the selection of the gene cluster consensus KOs from which the network
        was built.

    discard_ties : bool, None
        A parameter used in the selection of the gene cluster consensus KOs from which the network
        was built.

    consistent_annotations : bool, None
        A loaded network may be based on a set of gene KO annotations in the genomes storage
        database that has since changed, in which case this attribute would be False.

    gene_clusters : Dict[str, GeneCluster], dict()
        This maps the IDs of gene clusters in the network to object representations of the clusters.

    bins : Dict[str, GeneClusterBin], dict()

    collection : BinCollection, None
    """
    def __init__(
        self,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        verbose: bool = True
    ) -> None:
        """
        Parameters
        ==========
        run : anvio.terminal.Run, anvio.terminal.Run()
            This object sets the 'run' attribute, which prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object sets the 'progress' attribute, which prints transient progress information
            to the terminal.

        verbose : bool, True
            This sets the 'verbose' attribute, causing more information to be reported to the
            terminal if True.

        Returns
        =======
        None
        """
        super().__init__(run=run, progress=progress, verbose=verbose)
        self.pan_db_source_path: str = None
        self.genomes_storage_db_source_path: str = None
        self.consensus_threshold: float = None
        self.discard_ties: bool = None
        self.consistent_annotations: bool = None
        self.gene_clusters: Dict[str, GeneCluster] = {}
        self.bins: Dict[str, GeneClusterBin] = {}
        self.collection: BinCollection = None

    def copy(self) -> PangenomicNetwork:
        """
        Create a deep copy of the reaction network.

        Returns
        =======
        PangenomicNetwork
            Deep copy of the reaction network.
        """
        copied_network = PangenomicNetwork()

        self._copy(copied_network)

        for cluster_id, cluster in self.gene_clusters.items():
            copied_cluster = GeneCluster()
            copied_cluster.gene_cluster_id = cluster_id
            copied_cluster.genomes = cluster.genomes.copy()
            copied_cluster.ko = copied_network.kos[cluster.ko.id]

            copied_network.gene_clusters[cluster_id] = copied_cluster

        return copied_network

    def remove_metabolites_without_formula(self, output_path: str = None) -> None:
        """
        Remove metabolites without a formula in the ModelSEED database from the network.

        Other items can be removed from the network by association: reactions that involve a
        formulaless metabolite; other metabolites with formulas that are exclusive to such
        reactions; KOs predicted to exclusively catalyze such reactions; and gene clusters annotated
        with such KOs. Removed metabolites with a formula are reported alongside formulaless
        metabolites to the output table of removed metabolites.

        output_path : str, None
            If not None, write four tab-delimited files of metabolites, reactions, KEGG Orthologs,
            and gene clusters removed from the network to file locations based on the provided path.
            For example, if the argument, 'removed.tsv', is provided, then the following files will
            be written: 'removed-metabolites.tsv', 'removed-reactions.tsv', 'removed-kos.tsv', and
            'removed-gene-clusters.tsv'.
        """
        if output_path:
            path_basename, path_extension = os.path.splitext(output_path)
            metabolite_path = f"{path_basename}-metabolites{path_extension}"
            reaction_path = f"{path_basename}-reactions{path_extension}"
            ko_path = f"{path_basename}-kos{path_extension}"
            gene_cluster_path = f"{path_basename}-gene-clusters{path_extension}"
            for path in (metabolite_path, reaction_path, ko_path, gene_cluster_path):
                filesnpaths.is_output_file_writable(path)

        metabolites_to_remove = []
        for modelseed_compound_id, metabolite in self.metabolites.items():
            # ModelSEED compounds without a formula have a formula value of None in the network
            # object.
            if metabolite.formula is None:
                metabolites_to_remove.append(modelseed_compound_id)
        removed = self.purge_metabolites(metabolites_to_remove)

        if self.verbose:
            self.run.info("Removed metabolites", len(removed['metabolite']))
            self.run.info("Removed reactions", len(removed['reaction']))
            self.run.info("Removed KOs", len(removed['ko']))
            self.run.info("Removed gene clusters", len(removed['gene_cluster']))

        if not output_path:
            return

        # Record the reactions removed as a consequence of involving formulaless metabolites, and
        # record the formulaless metabolites involved in removed reactions.
        metabolite_removed_reactions: Dict[str, List[str]] = {}
        reaction_removed_metabolites: Dict[str, List[str]] = {}
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            reaction_removed_metabolites[reaction.modelseed_id] = metabolite_ids = []
            for metabolite in reaction.compounds:
                if metabolite.modelseed_id in metabolites_to_remove:
                    try:
                        metabolite_removed_reactions[metabolite.modelseed_id].append(
                            reaction.modelseed_id
                        )
                    except KeyError:
                        metabolite_removed_reactions[metabolite.modelseed_id] = [
                            reaction.modelseed_id
                        ]
                    metabolite_ids.append(metabolite.modelseed_id)

        metabolite_table = []
        for metabolite in removed['metabolite']:
            metabolite: ModelSEEDCompound
            row = []
            row.append(metabolite.modelseed_id)
            row.append(metabolite.modelseed_name)
            row.append(metabolite.formula)
            try:
                # The metabolite did not have a formula.
                removed_reaction_ids = metabolite_removed_reactions[metabolite.modelseed_id]
            except KeyError:
                # The metabolite had a formula but was removed as a consequence of all the reactions
                # involving the metabolite being removed due to them containing formulaless
                # metabolites: the metabolite did not cause any reactions to be removed.
                row.append("")
                continue
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus the reaction is recorded multiple times.
            row.append(", ".join(sorted(set(removed_reaction_ids))))

        reaction_table = []
        for reaction in removed['reaction']:
            reaction: ModelSEEDReaction
            row = []
            row.append(reaction.modelseed_id)
            row.append(reaction.modelseed_name)
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus is recorded multiple times.
            row.append(
                ", ".join(set(reaction_removed_metabolites[reaction.modelseed_id]))
            )
            row.append(", ".join([metabolite.modelseed_id for metabolite in reaction.compounds]))
            row.append(get_chemical_equation(reaction))
            reaction_table.append(row)

        ko_table = []
        for ko in removed['ko']:
            ko: KO
            row = []
            row.append(ko.id)
            row.append(ko.name)
            row.append(", ".join(ko.reactions))
            ko_table.append(row)

        gene_cluster_table = []
        for cluster in removed['gene_cluster']:
            cluster: GeneCluster
            row = []
            row.append(cluster.gene_cluster_id)
            row.append(cluster.ko.id)
            row.append(", ".join(cluster.genomes))
            gene_cluster_table.append(row)

        pd.DataFrame(
            metabolite_table,
            columns=[
                "ModelSEED compound ID",
                "ModelSEED compound name",
                "Formula",
                "Removed reaction ModelSEED IDs"
            ]
        ).to_csv(metabolite_path, sep='\t', index=False)
        pd.DataFrame(
            reaction_table,
            columns=[
                "ModelSEED reaction ID",
                "ModelSEED reaction name",
                "Removed ModelSEED compound IDs",
                "Reaction ModelSEED compound IDs",
                "Equation"
            ]
        ).to_csv(reaction_path, sep='\t', index=False)
        pd.DataFrame(
            ko_table,
            columns=[
                "KO ID",
                "KO name",
                "KO ModelSEED reaction IDs"
            ]
        ).to_csv(ko_path, sep='\t', index=False)
        pd.DataFrame(
            gene_cluster_table,
            columns=[
                "Gene cluster ID",
                "KO ID",
                "Gene cluster genomes"
            ]
        ).to_csv(gene_cluster_path, sep='\t', index=False)

        if self.verbose:
            self.run.info("Table of removed metabolites", metabolite_path)
            self.run.info("Table of removed reactions", reaction_path)
            self.run.info("Table of removed KOs", ko_path)
            self.run.info("Table of removed gene clusters", gene_cluster_path)

    def purge_metabolites(self, metabolites_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given metabolites from the network.

        Reactions involving the metabolite are also purged from the network. KOs that were only
        associated with removed reactions are purged; gene clusters that were only associated with
        removed KOs are purged.

        Removal of reactions involving the metabolite can also result in other metabolites being
        being removed from the network, those that exclusively participate in these reactions.

        Parameters
        ==========
        metabolites_to_remove : Iterable[str]
            ModelSEED compound IDs identifying metabolites to remove.

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
                'gene_cluster': [<removed GeneCluster objects>]
            }

            If this method is called from the method, 'purge_reactions', then the dictionary will
            only contain one significant entry:
            {
                'metabolite': [<removed ModelSEEDCompound objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene_cluster': []
            }
        """
        removed_metabolites: List[ModelSEEDCompound] = []
        for compound_id in metabolites_to_remove:
            try:
                removed_metabolites.append(self.metabolites.pop(compound_id))
            except KeyError:
                # This can occur for two reasons. First, the metabolite from 'metabolites_to_remove'
                # could not be in the network.

                # Second, this can occur when removing other "unintended" metabolites from the
                # network. 'purge_metabolites' was first called with metabolites of interest, then
                # 'purge_reactions' was called from within the method the remove reactions involving
                # the metabolites of interest, and then 'purge_metabolites' was called again from
                # within 'purge_reactions' to remove other metabolites exclusively found in the
                # removed reactions. In this last call of 'purge_metabolites', the
                # 'metabolites_to_remove' also include the metabolites of interest that were already
                # removed from 'self.metabolites' in the original 'purge_metabolites' call. This
                # KeyError occurs when trying to remove those already-removed metabolites.
                pass
        if not removed_metabolites:
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene_cluster': []
            }

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
            removed_cascading_up = {
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene_cluster': []
            }
        removed.update(removed_cascading_up)
        return removed

    def purge_reactions(self, reactions_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given reactions from the network.

        Metabolites that exclusively participate in removed reactions are purged. KOs that were only
        associated with removed reactions are purged; gene clusters that were only associated with
        removed KOs are purged.

        Parameters
        ==========
        reactions_to_remove : Iterable[str]
            ModelSEED reaction IDs identifying reactions to remove.

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
                'gene_cluster': [<removed GeneCluster objects>]
            }

            If this method is called from the method, 'purge_metabolites', then the dictionary will
            look like the following:
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [],
                'ko': [<removed KO objects>],
                'gene_cluster': [<removed GeneCluster objects>]
            }

            If this method is called from the method, 'purge_kos', then the dictionary will look
            like the following:
            {
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects],
                'ko': [],
                'gene_cluster': []
            }
        """
        removed_reactions: List[ModelSEEDReaction] = []
        for modelseed_reaction_id in reactions_to_remove:
            try:
                removed_reactions.append(self.reactions.pop(modelseed_reaction_id))
            except KeyError:
                # This occurs when the original method called is 'purge_reactions', followed by
                # 'purge_kos', followed by this method again -- 'removed_reactions' will be empty.
                # Alternatively, this occurs if the reaction in 'reactions_to_remove' is not in the
                # network.
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
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene_cluster': []
            }

        # Remove KEGG reaction aliases of ModelSEED reactions in the network.
        kegg_reactions_to_remove = []
        for kegg_reaction_id, modelseed_reaction_ids in self.kegg_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the KEGG reaction were removed, so remove the
                # KEGG reaction as well.
                kegg_reactions_to_remove.append(kegg_reaction_id)
                continue
            for idx in sorted(aliases_to_remove, reverse=True):
                # Only some of the ModelSEED reactions aliased by the KO KEGG reaction were removed.
                modelseed_reaction_ids.pop(idx)
        removed_kegg_reactions = []
        for kegg_reaction_id in kegg_reactions_to_remove:
            removed_kegg_reactions.append(self.kegg_modelseed_aliases.pop(kegg_reaction_id))

        # Remove EC number aliases of ModelSEED reactions in the network.
        ec_numbers_to_remove = []
        for ec_number, modelseed_reaction_ids in self.ec_number_modelseed_aliases.items():
            aliases_to_remove: List[int] = []
            for idx, modelseed_reaction_id in enumerate(modelseed_reaction_ids):
                if modelseed_reaction_id in reactions_to_remove:
                    aliases_to_remove.append(idx)
            if len(aliases_to_remove) == len(modelseed_reaction_ids):
                # All ModelSEED reactions aliased by the EC number were removed, so remove the EC
                # number as well.
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
                        # Do not remove the metabolite, because it participates in a retained
                        # reaction.
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
        # exclusively associated with reactions removed here, then 'kos_to_remove' would likewise be
        # empty.
        if kos_to_remove:
            removed_cascading_up = self.purge_kos(kos_to_remove)
            removed_cascading_up.pop('reaction')
            removed_cascading_up.pop('kegg_reaction')
            removed_cascading_up.pop('ec_number')
        else:
            removed_cascading_up = {'ko': [], 'gene_cluster': []}

        removed = {
            'reaction': removed_reactions,
            'kegg_reaction': removed_kegg_reactions,
            'ec_number': removed_ec_numbers
        }
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def purge_kos(self, kos_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given KOs from the network.

        Reactions and metabolites that were only associated with removed KOs are purged. Genes that
        were only associated with removed KOs are purged.

        Parameters
        ==========
        kos_to_remove : Iterable[str]
            KO IDs identifying KOs to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, 'purge_reactions', or the method,
            'purge_gene_clusters', then the dictionary will look like the following:
            {
                'ko': [<removed KO objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>],
                'gene_cluster': [<removed GeneCluster objects>]
            }

            If this method is called from the method, 'purge_reactions', then the dictionary will
            look like the following:
            {
                'ko': [<removed KO objects>],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': [],
                'gene_cluster': [<removed GeneCluster objects>]
            }

            If this method is called from the method, 'purge_gene_clusters', then the dictionary
            will look like the following:
            {
                'ko': [<removed KO objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects],
                'gene_clusters': []
            }
        """
        removed_kos: List[KO] = []
        for ko_id in kos_to_remove:
            try:
                removed_kos.append(self.kos.pop(ko_id))
            except KeyError:
                # This occurs when the original method called is 'purge_kos', followed by
                # 'purge_gene_clusters', followed by this method again -- 'removed_kos' will be
                # empty. Alternatively, this occurs if the KO in 'kos_to_remove' is not in the
                # network.
                pass

        if not removed_kos:
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene_cluster': []
            }

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
                        # The reaction is associated with a retained KO, so do not remove the
                        # reaction.
                        reactions_to_spare.append(idx)
            for idx in sorted(reactions_to_spare, reverse=True):
                reactions_to_remove.pop(idx)
        if reactions_to_remove:
            removed_cascading_down = self.purge_reactions(reactions_to_remove)
            removed_cascading_down.pop('ko')
        else:
            # This method must have been called from the method, 'purge_reactions', because the
            # reactions that are only associated with the removed KOs were already removed from the
            # network.
            removed_cascading_down = {
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        gene_clusters_to_remove: List[str] = []
        for gene_cluster_id, cluster in self.gene_clusters.items():
            if cluster.ko.id in kos_to_remove:
                gene_clusters_to_remove.append(gene_cluster_id)
        # If this method was called from 'purge_gene_clusters' then the gene clusters that are only
        # associated with KOs removed here were already removed from the network, and
        # 'gene_clusters_to_remove' would be empty.
        if gene_clusters_to_remove:
            removed_cascading_up = self.purge_gene_clusters(gene_clusters_to_remove)
            removed_cascading_up.pop('ko')
        else:
            removed_cascading_up = {'gene_cluster': []}

        removed = {'ko': removed_kos}
        removed.update(removed_cascading_down)
        removed.update(removed_cascading_up)
        return removed

    def purge_gene_clusters(self, gene_clusters_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given gene clusters from the network.

        KOs, reactions, and metabolites that were only associated with removed gene clusters are
        purged.

        Parameters
        ==========
        gene_clusters_to_remove : Iterable[str]
            Gene cluster IDs identifying clusters to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, 'purge_kos', then the dictionary will
            look like the following:
            {
                'gene_cluster': [<removed GeneCluster objects>],
                'ko': [<removed KO objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>]
            }

            If this method is called from the method, 'purge_kos', then the dictionary will look
            like the following:
            {
                'gene_cluster': [<removed GeneCluster objects>],
                'ko': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }
        """
        removed_gene_clusters: List[GeneCluster] = []
        for gene_cluster_id in gene_clusters_to_remove:
            try:
                removed_gene_clusters.append(self.gene_clusters.pop(gene_cluster_id))
            except KeyError:
                # This occurs if the cluster in 'gene_clusters_to_remove' is not in the network.
                pass

        if not removed_gene_clusters:
            return {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'gene_cluster': []
            }

        kos_to_remove: List[str] = []
        for cluster in removed_gene_clusters:
            kos_to_remove.append(cluster.ko.id)
        kos_to_remove = list(set(kos_to_remove))
        for gene_cluster in self.gene_clusters.values():
            kos_to_spare: List[str] = []
            if gene_cluster.ko.id in kos_to_remove:
                # The KO is associated with a retained gene cluster, so do not remove the KO.
                kos_to_spare.append(gene_cluster.ko.id)
            for ko_id in kos_to_spare:
                kos_to_remove.remove(ko_id)
        if kos_to_remove:
            removed_cascading_down = self.purge_kos(kos_to_remove)
            removed_cascading_down.pop('gene_cluster')
        else:
            # This method must have been called from the method, 'purge_kos', because the KOs that
            # are only associated with the removed gene clusters were already removed from the
            # network.
            removed_cascading_down = {
                'ko': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        # TODO: remove gene clusters from self.bins

        removed = {'gene_cluster': removed_gene_clusters}
        removed.update(removed_cascading_down)
        return removed

    def subset_network(
        self,
        kegg_modules_to_subset: Iterable[str] = None,
        brite_categories_to_subset: Iterable[str] = None,
        gene_clusters_to_subset: Iterable[int] = None,
        kos_to_subset: Iterable[str] = None,
        reactions_to_subset: Iterable[str] = None,
        metabolites_to_subset: Iterable[str] = None
    ) -> PangenomicNetwork:
        """
        Subset a smaller network from the metabolic network.

        If requested KEGG modules, BRITE categories, gene clusters, KOs, reactions, or metabolites
        are not present in the network, no error is raised.

        Subsetted items are not represented by the same objects as in the source network, i.e., new
        gene cluster, KO, reaction, and metabolite objects are created and added to the subsetted
        network.

        Network items (i.e., gene clusters, KOs, reactions, and metabolites) that reference
        requested items (e.g., gene clusters in the network referencing requested KOs; KOs in the
        network referencing requested reactions) are added to the subsetted network. KOs (and by
        extension, gene clusters referencing KOs) that are added to the subsetted network due to
        references to requested reactions will be missing references to any other unrequested
        reactions. In other words, certain reaction annotations can be selected to the exclusion of
        others, e.g., a KO encoding two reactions can be "redefined" or "pruned" to encode one
        requested reaction in the subsetted network; a KO encoding multiple reactions can be pruned
        to encode only those reactions involving requested metabolites.

        If the 'verbose' attribute of the source 'PangenomicNetwork' object is True, then report to
        the terminal the identities of requested KEGG modules, BRITE categories, gene clusters, KOs,
        reactions, and metabolites that are not present in the network.

        Parameters
        ==========
        kegg_modules_to_subset : List[str], None
            KEGG modules (of KOs) to subset by ID.

        brite_categories_to_subset : List[str], None
            KEGG BRITE hierarchy categories (of KOs) to subset.

        gene_clusters_to_subset : List[int], None
            Gene clusters to subset by ID.

        kos_to_subset : List[str], None
            KOs to subset by ID.

        reactions_to_subset : List[str], None
            ModelSEED reactions to subset by ID.

        metabolites_to_subset : List[str], None
            ModelSEED metabolites to subset by ID.

        Returns
        =======
        PangenomicNetwork
            New subsetted reaction network.
        """
        # Sequentially subset the network for each type of request. Upon generating two subsetted
        # networks from two types of request, merge the networks into a single subsetted network;
        # repeat.
        first_subnetwork = None
        for items_to_subset, subset_network_method in (
            (kegg_modules_to_subset, self._subset_network_by_modules),
            (brite_categories_to_subset, self._subset_network_by_brite),
            (gene_clusters_to_subset, self._subset_network_by_gene_clusters),
            (kos_to_subset, self._subset_network_by_kos),
            (reactions_to_subset, self._subset_network_by_reactions),
            (metabolites_to_subset, self._subset_network_by_metabolites)
        ):
            if not items_to_subset:
                continue

            second_subnetwork = subset_network_method(items_to_subset)

            if first_subnetwork is None:
                first_subnetwork = second_subnetwork
            else:
                first_subnetwork = first_subnetwork.merge_network(second_subnetwork)

        return first_subnetwork

    def _subset_network_by_modules(self, kegg_modules: Iterable[str]) -> PangenomicNetwork:
        """
        Subset the network by KOs in requested KEGG modules.

        Parameters
        ==========
        kegg_modules : Iterable[str]
            KEGG modules (of KOs) to subset by ID.

        Returns
        =======
        PangenomicNetwork
            New subsetted reaction network.
        """
        pass

    def _subset_network_by_brite(self, brite_categories: Iterable[str]) -> PangenomicNetwork:
        """
        Subset the network by KOs in requested KEGG BRITE hierarchy categories.

        Parameters
        ==========
        brite_categories : Iterable[str]
            KEGG BRITE hierarchy categories (of KOs) to subset.

        Returns
        =======
        PangenomicNetwork
            New subsetted reaction network.
        """
        pass

    def _subset_network_by_gene_clusters(
        self,
        gene_cluster_ids: Iterable[int]
    ) -> PangenomicNetwork:
        """
        Subset the network by gene clusters with requested IDs.

        Parameters
        ==========
        gene_cluster_ids : Iterable[int]
            Gene clusters to subset by ID.

        Returns
        =======
        PangenomicNetwork
            New subsetted reaction network.
        """
        subnetwork = PangenomicNetwork()

        for gene_cluster_id in gene_cluster_ids:
            try:
                cluster = self.gene_clusters[gene_cluster_id]
            except KeyError:
                # This occurs if the requested gene cluster ID is not in the source network.
                continue

            subsetted_cluster = GeneCluster()
            subsetted_cluster.gene_cluster_id = gene_cluster_id
            subsetted_cluster.genomes = cluster.genomes.copy()

            # Add KOs annotating the gene cluster to the subsetted network as new objects, and then
            # reference these objects in the cluster object.
            ko_id = cluster.ko.id
            subnetwork = self._subset_network_by_kos([ko_id], subnetwork=subnetwork)
            subsetted_cluster.ko = subnetwork.kos[ko_id]

            subnetwork.gene_clusters[gene_cluster_id] = subsetted_cluster

        return subnetwork

    def _subset_network_by_kos(
        self,
        ko_ids: Iterable[str],
        subnetwork: PangenomicNetwork = None
    ) -> PangenomicNetwork:
        """
        Subset the network by KOs with requested KO IDs.

        Parameters
        ==========
        ko_ids : Iterable[str]
            KOs to subset by ID.

        subnetwork : PangenomicNetwork, None
            This network under construction is provided when the KOs being added to the network
            annotate already subsetted gene clusters.

        Returns
        =======
        PangenomicNetwork
            If a 'subnetwork' argument is provided, then that network is returned after
            modification. Otherwise, a new subsetted reaction network is returned.
        """
        if subnetwork is None:
            subnetwork = PangenomicNetwork()
            # Signify that gene clusters annotated by subsetted KOs are to be added to the network.
            subset_referencing_gene_clusters = True
        else:
            assert isinstance(subnetwork, PangenomicNetwork)
            # Signify that the KOs being added to the network annotate subsetted gene clusters that
            # were already added to the network.
            subset_referencing_gene_clusters = False

        for ko_id in ko_ids:
            try:
                ko = self.kos[ko_id]
            except KeyError:
                # This occurs if the requested KO ID is not in the source network.
                continue

            subsetted_ko = KO()
            subsetted_ko.id = ko.id
            subsetted_ko.name = ko.name
            subsetted_ko.kegg_reaction_aliases = deepcopy(ko.kegg_reaction_aliases)
            subsetted_ko.ec_number_aliases = deepcopy(ko.ec_number_aliases)

            # Add reactions annotating the KO to the subsetted network as new objects, and then
            # reference these objects in the KO object.
            reaction_ids = [reaction_id for reaction_id in ko.reactions]
            subnetwork = self._subset_network_by_reactions(reaction_ids, subnetwork=subnetwork)
            subsetted_ko.reactions = {
                reaction_id: subnetwork.reactions[reaction_id] for reaction_id in reaction_ids
            }

            subnetwork.kos[ko_id] = subsetted_ko

        if subset_referencing_gene_clusters:
            # Add gene clusters that are annotated by the subsetted KOs to the network.
            self._subset_gene_clusters_via_kos(subnetwork)

        return subnetwork

    def _subset_network_by_reactions(
        self,
        reaction_ids: Iterable[str],
        subnetwork: PangenomicNetwork = None
    ) -> PangenomicNetwork:
        """
        Subset the network by reactions with ModelSEED reaction IDs.

        Parameters
        ==========
        reaction_ids : Iterable[str]
            Reactions to subset by ModelSEED reaction ID.

        subnetwork : PangenomicNetwork, None
            This network under construction is provided when the reactions being added to the
            network annotate already subsetted KOs.

        Returns
        =======
        PangenomicNetwork
            If a 'subnetwork' argument is provided, then that network is returned after
            modification. Otherwise, a new subsetted reaction network is returned.
        """
        if subnetwork is None:
            subnetwork = PangenomicNetwork()
            # Signify that KOs annotated by subsetted reactions are to be added to the network.
            subset_referencing_kos = True
        else:
            assert isinstance(subnetwork, PangenomicNetwork)
            # Signify that the reactions being added to the network annotate subsetted KOs that were
            # already added to the network.
            subset_referencing_kos = False

        # Copy the network attributes mapping reaction aliases.
        kegg_modelseed_aliases: Dict[str, List[str]] = {}
        ec_number_modelseed_aliases: Dict[str, List[str]] = {}

        for reaction_id in reaction_ids:
            try:
                reaction = self.reactions[reaction_id]
            except KeyError:
                # This occurs if the requested reaction is not in the source network.
                continue

            # Copy the reaction object, including referenced metabolite objects, from the source
            # network.
            subsetted_reaction: ModelSEEDReaction = deepcopy(reaction)
            subnetwork.reactions[reaction_id] = subsetted_reaction
            # Record the metabolites involved in the reaction, and add them to the network.
            for metabolite in subsetted_reaction.compounds:
                compound_id = metabolite.modelseed_id
                subnetwork.metabolites[compound_id] = metabolite

            try:
                subnetwork.modelseed_kegg_aliases[reaction_id] += list(reaction.kegg_aliases)
            except KeyError:
                subnetwork.modelseed_kegg_aliases[reaction_id] = list(reaction.kegg_aliases)

            try:
                subnetwork.modelseed_ec_number_aliases[reaction_id] += list(
                    reaction.ec_number_aliases
                )
            except KeyError:
                subnetwork.modelseed_ec_number_aliases[reaction_id] = list(
                    reaction.ec_number_aliases
                )

            for kegg_id in reaction.kegg_aliases:
                try:
                    kegg_modelseed_aliases[kegg_id].append(reaction_id)
                except KeyError:
                    kegg_modelseed_aliases[kegg_id] = [reaction_id]

            for ec_number in reaction.ec_number_aliases:
                try:
                    ec_number_modelseed_aliases[ec_number].append(reaction_id)
                except KeyError:
                    ec_number_modelseed_aliases[ec_number] = [reaction_id]

        if subnetwork.kegg_modelseed_aliases:
            for kegg_id, modelseed_ids in kegg_modelseed_aliases.items():
                try:
                    subnetwork.kegg_modelseed_aliases[kegg_id] += modelseed_ids
                except KeyError:
                    subnetwork.kegg_modelseed_aliases[kegg_id] = modelseed_ids
        else:
            subnetwork.kegg_modelseed_aliases = kegg_modelseed_aliases

        if subnetwork.ec_number_modelseed_aliases:
            for ec_number, modelseed_ids in ec_number_modelseed_aliases.items():
                try:
                    subnetwork.ec_number_modelseed_aliases[ec_number] += modelseed_ids
                except KeyError:
                    subnetwork.ec_number_modelseed_aliases[ec_number] = modelseed_ids
        else:
            subnetwork.ec_number_modelseed_aliases = ec_number_modelseed_aliases

        if subset_referencing_kos:
            # Add KOs that are annotated by the subsetted reactions to the network.
            self._subset_kos_via_reactions(subnetwork)

        return subnetwork

    def _subset_gene_clusters_via_kos(self, subnetwork: PangenomicNetwork) -> None:
        """
        Add gene clusters that are annotated with subsetted KOs to the subsetted network.

        Parameters
        ==========
        subnetwork : PangenomicNetwork
            The subsetted reaction network under construction.

        Returns
        =======
        None
        """
        subsetted_ko_ids = list(subnetwork.kos)
        for gene_cluster_id, cluster in self.gene_clusters.items():
            # Check all gene clusters in the source network for subsetted KOs.
            if cluster.ko.id in subsetted_ko_ids:
                # Create a new gene cluster object for the subsetted cluster.
                subsetted_cluster = GeneCluster()
                subsetted_cluster.gene_cluster_id = gene_cluster_id
                subsetted_cluster.genomes = cluster.genomes.copy()
                subsetted_cluster.ko = subnetwork.kos[cluster.ko.id]
                subnetwork.gene_clusters[gene_cluster_id] = subsetted_cluster

    def _subset_kos_via_reactions(self, subnetwork: PangenomicNetwork) -> None:
        """
        Add KOs that are annotated with subsetted reactions to the subsetted network.

        Then add gene clusters that are annotated with these added KOs to the subsetted network.

        Parameters
        ==========
        subnetwork : PangenomicNetwork
            The subsetted reaction network under construction.

        Returns
        =======
        None
        """
        subsetted_reaction_ids = list(subnetwork.reactions)
        for ko_id, ko in self.kos.items():
            # Check all KOs in the source network for subsetted reactions.
            subsetted_ko = None
            for reaction_id in ko.reactions:
                if reaction_id not in subsetted_reaction_ids:
                    # The KO is not annotated by the subsetted reaction.
                    continue

                if not subsetted_ko:
                    # Create a new KO object for the subsetted KO. The subsetted KO object would
                    # already have been created had another subsetted reaction been among the
                    # reactions annotating the KO.
                    subsetted_ko = KO()
                    subsetted_ko.id = ko_id
                    subsetted_ko.name = ko.name
                subsetted_ko.reactions[reaction_id] = subnetwork.reactions[reaction_id]
                subsetted_ko.kegg_reaction_aliases = deepcopy(ko.kegg_reaction_aliases)
                subsetted_ko.ec_number_aliases = deepcopy(ko.ec_number_aliases)

            if subsetted_ko:
                subnetwork.kos[ko_id] = subsetted_ko

        # Add gene clusters that are annotated with the added KOs to the subsetted network.
        self._subset_gene_clusters_via_kos(subnetwork)

    def _subset_network_by_metabolites(self, compound_ids: List[str]) -> PangenomicNetwork:
        """
        Subset the network by metabolites with ModelSEED compound IDs.

        Parameters
        ==========
        compound_ids : List[str]
            List of metabolites to subset by ModelSEED compound ID.

        Returns
        =======
        PangenomicNetwork
            New subsetted reaction network.
        """
        subnetwork = PangenomicNetwork()

        for reaction_id, reaction in self.reactions.items():
            # Check all reactions in the source network for subsetted metabolites.
            for metabolite in reaction.compounds:
                if metabolite.modelseed_id in compound_ids:
                    break
            else:
                # The reaction does not involve any of the requested metabolites.
                continue

            # Copy the reaction object, including referenced metabolite objects, from the source
            # network.
            subsetted_reaction: ModelSEEDReaction = deepcopy(reaction)
            subnetwork.reactions[reaction_id] = subsetted_reaction

            # Add the metabolites involved in the reaction to the subsetted network. (There can be
            # unavoidable redundancy here in readding previously encountered metabolites.)
            for subsetted_metabolite in subsetted_reaction.compounds:
                subnetwork.metabolites[subsetted_metabolite.modelseed_id] = subsetted_metabolite

        # Add KOs that are annotated with the added reactions to the subsetted network, and then add
        # gene_clusters annotated with the added KOs to the subsetted network.
        self._subset_kos_via_reactions(subnetwork)

        return subnetwork

    def merge_network(self, network: PangenomicNetwork) -> PangenomicNetwork:
        """
        Merge the pangenomic reaction network with another pangenomic reaction network.

        Each network can contain different gene clusters, KOs, and reactions/metabolites. Merging
        nonredundantly incorporates all of this data as new objects in the new network.

        Objects representing KOs in both networks can have different sets of references: KOs can be
        annotated by different reactions.

        However, the same gene cluster in each network should have the same consensus KO annotation:
        the method of annotation should have been the same in each network. This reflects the
        purpose of the method, which is to combine different, but potentially overlapping,
        subnetworks from the same pangenome.

        ModelSEED reactions and metabolites in both networks should have identical attributes.

        Parameters
        ==========
        network : PangenomicNetwork
            The other pangenomic reaction network being merged.

        Returns
        =======
        PangenomicNetwork
            The merged pangenomic reaction network.
        """
        merged_network = PangenomicNetwork()

        self._merge_network(network, merged_network)

        # Add gene clusters to the merged network, first adding clusters present in both source
        # networks, and then adding clusters present exclusively in each source network.
        first_gene_cluster_ids = set(self.gene_clusters)
        second_gene_cluster_ids = set(network.gene_clusters)

        for gene_cluster_id in first_gene_cluster_ids.intersection(second_gene_cluster_ids):
            first_cluster = self.gene_clusters[gene_cluster_id]
            second_cluster = network.gene_clusters[gene_cluster_id]

            # The new object representing the gene cluster in the merged network should have all KO
            # annotations from each source cluster object, as these objects can have different KO
            # references.
            merged_cluster = GeneCluster()
            merged_cluster.gene_cluster_id = gene_cluster_id
            assert first_cluster.genomes == second_cluster.genomes
            merged_cluster.genomes = first_cluster.genomes.copy()
            assert first_cluster.ko.id == second_cluster.ko.id
            merged_cluster.ko = merged_network.kos[first_cluster.ko.id]

            merged_network.gene_clusters[gene_cluster_id] = merged_cluster

        for gene_cluster_id in first_gene_cluster_ids.difference(second_gene_cluster_ids):
            first_cluster = self.gene_clusters[gene_cluster_id]

            cluster = GeneCluster()
            cluster.gene_cluster_id = gene_cluster_id
            cluster.genomes = first_cluster.genomes.copy()
            cluster.ko = merged_network.kos[first_cluster.ko.id]

            merged_network.gene_clusters[gene_cluster_id] = cluster

        for gene_cluster_id in second_gene_cluster_ids.difference(first_gene_cluster_ids):
            second_cluster = network.gene_clusters[gene_cluster_id]

            cluster = GeneCluster()
            cluster.gene_cluster_id = gene_cluster_id
            cluster.genomes = second_cluster.genomes.copy()
            cluster.ko = merged_network.kos[second_cluster.ko.id]

            merged_network.gene_clusters[gene_cluster_id] = cluster

        return merged_network

    def get_overview_statistics(
        self,
        precomputed_counts: Dict[str, int] = None
    ) -> PangenomicNetworkStats:
        """
        Calculate overview statistics for the pangenomic metabolic network.

        Parameters
        ==========
        precomputed_counts : Dict[str, int], None
            To spare additional computations that involve loading and parsing databases, this
            dictionary can contain three pieces of precomputed data: the value for the key,
            'total_gene_clusters', should be the number of gene clusters in the pangenome; the value
            for the key, 'gene_clusters_assigned_ko', should be the number of gene clusters in the
            pangenome assigned a consensus KO (or None if 'self.consistent_annotations' is False);
            the value for the key, 'kos_assigned_gene_clusters', should be the number of consensus
            KOs assigned to gene clusters in the pangenome (or None if 'self.consistent_annotations'
            is False).

        Returns
        =======
        PangenomicNetworkStats
            Network statistics are stored in a dictionary of dictionaries. Keys in the outer
            dictionary are "classes" of network statistics. Keys in the inner dictionary are
            statistics themselves.
        """
        if (
            precomputed_counts is not None and
            sorted(precomputed_counts) != [
                'gene_clusters_assigned_ko', 'kos_assigned_gene_clusters', 'total_gene_clusters'
            ]
        ):
            raise ConfigError(
                "The 'precomputed_counts' argument must be a dictionary only containing the keys, "
                "'total_gene_clusters', 'gene_clusters_assigned_ko', and "
                "'kos_assigned_gene_clusters'."
            )

        stats: PangenomicNetworkStats = {}

        self.progress.new("Counting gene clusters and KOs")
        self.progress.update("...")
        stats['Gene cluster and KO counts'] = stats_group = {}

        if precomputed_counts:
            assert (
                type(precomputed_counts['total_gene_clusters']) is int and
                precomputed_counts['total_gene_clusters'] >= 0
            )
            gene_cluster_count = precomputed_counts['total_gene_clusters']
            assert (
                precomputed_counts['gene_clusters_assigned_ko'] is None or
                (
                    type(precomputed_counts['gene_clusters_assigned_ko']) is int and
                    precomputed_counts['gene_clusters_assigned_ko'] >= 0
                )
            )
            ko_annotated_gene_cluster_count = precomputed_counts['gene_clusters_assigned_ko']
            assert (
                precomputed_counts['kos_assigned_gene_clusters'] is None or
                (
                    type(precomputed_counts['kos_assigned_gene_clusters']) is int and
                    precomputed_counts['kos_assigned_gene_clusters'] >= 0
                )
            )
            annotating_ko_count = precomputed_counts['kos_assigned_gene_clusters']
            assert not (
                (ko_annotated_gene_cluster_count is None and annotating_ko_count is not None) or
                (ko_annotated_gene_cluster_count is not None and annotating_ko_count is None)
            )
        else:
            # One database cannot be available without the other.
            assert not (
                (
                    self.pan_db_source_path is None and
                    self.genomes_storage_db_source_path is not None
                ) or
                (
                    self.pan_db_source_path is not None and
                    self.genomes_storage_db_source_path is None
                )
            )

            if self.pan_db_source_path and self.genomes_storage_db_source_path:
                pdb = PanDatabase(self.pan_db_source_path)
                gene_cluster_count = pdb.meta['num_gene_clusters']
                pdb.disconnect()
            else:
                gene_cluster_count = None

            if (
                self.pan_db_source_path and
                self.genomes_storage_db_source_path and
                self.consistent_annotations is False
            ):
                args = argparse.Namespace()
                args.genomes_storage = self.genomes_storage_db_source_path
                args.consensus_threshold = self.consensus_threshold
                args.discard_ties = self.discard_ties
                pan_super = PanSuperclass(args, r=run_quiet)
                pan_super.init_gene_clusters()
                pan_super.init_gene_clusters_functions()
                pan_super.init_gene_clusters_functions_summary_dict()
                gene_clusters_functions_summary_dict: Dict = (
                    pan_super.gene_clusters_functions_summary_dict
                )
                ko_annotated_gene_cluster_count = 0
                ko_ids = []
                for gene_cluster_functions_data in gene_clusters_functions_summary_dict.values():
                    gene_cluster_ko_data = gene_cluster_functions_data['KOfam']
                    if gene_cluster_ko_data != {'function': None, 'accession': None}:
                        # A KO was assigned to the cluster.
                        ko_annotated_gene_cluster_count += 1
                        ko_ids.append(gene_cluster_ko_data['accession'])
                annotating_ko_count = len(set(ko_ids))
            else:
                ko_annotated_gene_cluster_count = None
                annotating_ko_count = None

        if gene_cluster_count is not None:
            stats_group['Total gene clusters in pangenome'] = gene_cluster_count
        if ko_annotated_gene_cluster_count is not None:
            stats_group['Gene clusters assigned protein KO'] = ko_annotated_gene_cluster_count
        stats_group['Gene clusters in network'] = len(self.gene_clusters)
        if annotating_ko_count is not None:
            stats_group['Protein KOs assigned to gene clusters'] = annotating_ko_count
        stats_group['KOs in network'] = len(self.kos)
        self.progress.end()

        self._get_common_overview_statistics(stats)

        if precomputed_counts:
            return stats

        if not (self.pan_db_source_path and self.genomes_storage_db_source_path):
            self.run.info_single(
                f"""\
                Since the pangenomic network was not associated with a pan database and genomes
                storage database, the following statistics could not be calculated and were not
                reported to the output file: 'Total gene clusters in pangenome', 'Gene clusters
                assigned protein KOs', and 'Protein KOs assigned to gene clusters'.\
                """
            )
        elif self.consistent_annotations is False:
            self.run.info_single(
                f"""\
                The network attribute, 'consistent_annotations', is False, which indicates that the
                reaction network stored in the pan database was made from a different set of KO gene
                annotations than is currently in the genomes storage database. Therefore, the
                following statistics were not calculated and reported to the output file to avoid
                potential inaccuracies: 'Gene clusters assigned protein KO' and 'Protein KOs
                assigned to gene clusters'.\
                """
            )

        return stats

    def print_overview_statistics(self, stats: GenomicNetworkStats = None) -> None:
        """
        Print overview statistics for the genomic metabolic network.

        Parameters
        ==========
        stats : GenomicNetworkStats, None
            With the default value of None, network statistics will be calculated and printed.
            Alternatively, provided network statistics will be printed without calculating anew.

        Returns
        =======
        None
        """
        if not stats:
            stats = self.get_overview_statistics()

        self.run.info_single("METABOLIC REACTION NETWORK STATISTICS", mc='green', nl_after=1)

        self.run.info_single("Gene clusters and KEGG Ortholog (KO) annotations")
        stats_group = stats['Gene cluster and KO counts']
        self.run.info(
            "Total gene clusters in pangenome", stats_group['Total gene clusters in pangenome']
        )
        self.run.info(
            "Gene clusters annotated with protein KO",
            stats_group['Gene clusters assigned protein KO']
        )
        self.run.info("Gene clusters in network", stats_group['Gene clusters in network'])
        self.run.info(
            "Protein KOs assigned to gene clusters",
            stats_group['Protein KOs assigned to gene clusters']
        )
        self.run.info("KOs in network", stats_group['KOs in network'], nl_after=1)

        self._print_common_overview_statistics(stats)

    def export_json(
        self,
        path: str,
        overwrite: bool = False,
        objective: str = None,
        remove_missing_objective_metabolites: bool = False,
        record_genomes: Tuple[str] = ('gene', 'reaction'),
        # record_bins: Tuple[str] = ('gene', 'reaction'),
        indent: int = 2,
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Export the network to a metabolic model file in JSON format. Entries in the "gene" section
        of this file represent gene clusters.

        All information from the network is included in the JSON so that the file can by imported by
        anvi'o as a PangenomicNetwork object containing the same information.

        Parameters
        ==========
        path : str
            output JSON file path

        overwrite : bool, False
            Overwrite the JSON file if it already exists.

        objective : str, None
            An objective to use in the model, stored as the first entry in the JSON 'reactions'
            array. Currently, the only valid options are None and 'e_coli_core'.

            None means that no objective is added to the JSON, meaning that FBA cannot be performed
            on the model.

            'e_coli_core' is the biomass objective from the COBRApy example JSON file of E. coli
            "core" metabolism, 'e_coli_core.json'.

        remove_missing_objective_metabolites : bool, False
            If True, remove metabolites from the JSON objective that are not produced or consumed in
            the reaction network. FBA fails with metabolites outside the network.

        record_genomes : tuple, ('gene cluster', 'reaction')
            Record the genome membership of gene clusters in JSON entries. By default, genome names
            are recorded for gene clusters and reactions with the argument, ('gene cluster',
            'reaction'). To not record genomes at all, pass either an empty tuple or None. The
            following valid strings can be provided in a tuple in any combination: 'gene cluster',
            'reaction', and 'metabolite'. 'reaction' and 'metabolite' record the genomes predicted
            to encode enzymes associated with reactions and metabolites, respectively.

        indent : int, 2
            spaces of indentation per nesting level in JSON file

        progress : terminal.Progress, terminal.Progress()
        """
        if record_genomes is None:
            record_genomes = ()
        valid_items = ('gene cluster', 'reaction', 'metabolite')
        invalid_items = []
        for item in record_genomes:
            if item not in valid_items:
                invalid_items.append(item)
        if invalid_items:
            raise ConfigError(
                f"The following items in the 'record_genomes' argument are invalid: {', '.join(invalid_items)}"
            )

        progress.new("Constructing JSON")
        progress.update("Setting up")
        filesnpaths.is_output_file_writable(path, ok_if_exists=overwrite)
        json_dict = JSONStructure.get()
        json_gene_clusters: List[Dict] = json_dict['genes']
        json_reactions: List[Dict] = json_dict['reactions']
        json_metabolites: List[Dict] = json_dict['metabolites']
        if objective == 'e_coli_core':
            objective_dict = JSONStructure.get_e_coli_core_objective()
            if remove_missing_objective_metabolites:
                self.remove_missing_objective_metabolites(objective_dict)
            json_reactions.append(objective_dict)
        elif objective != None:
            raise ConfigError(f"Anvi'o does not recognize an objective with the name, '{objective}'.")

        progress.update("Gene clusters")
        reaction_gene_clusters: Dict[str, List[str]] = {}
        reaction_kos: Dict[str, List[KO]] = {}
        # The following two dictionaries are only needed for recording the occurrence of reactions
        # and metabolites in genomes.
        reaction_genomes: Dict[str, List[str]] = {}
        metabolite_genomes: Dict[str, List[str]] = {}
        for cluster_id, gene_cluster in self.gene_clusters.items():
            gene_cluster_entry = JSONStructure.get_gene_entry()
            json_gene_clusters.append(gene_cluster_entry)
            cluster_id_str = str(cluster_id)
            gene_cluster_entry['id'] = cluster_id_str
            # Record KO IDs in the annotation section of the gene cluster entry. In a JSON file produced
            # from a 'GenomicNetwork', KO IDs are paired with their gene annotation e-values, which
            # can't be done with consensus KOs for gene clusters. Therefore, where the e-value would
            # be, put an empty string.
            annotation = gene_cluster_entry['annotation']
            annotation['ko'] = annotation_kos = {}
            ko = gene_cluster.ko
            annotation_kos[ko.id] = ""
            for modelseed_reaction_id in ko.reactions:
                try:
                    reaction_gene_clusters[modelseed_reaction_id].append(cluster_id_str)
                except KeyError:
                    reaction_gene_clusters[modelseed_reaction_id] = [cluster_id_str]
                try:
                    reaction_kos[modelseed_reaction_id].append(ko)
                except KeyError:
                    reaction_kos[modelseed_reaction_id] = [ko]
            if not record_genomes:
                continue
            genome_names = gene_cluster.genomes
            if 'gene cluster' in record_genomes:
                # Record the names of the genomes contributing to the gene cluster in the notes section
                # of the gene cluster entry.
                gene_cluster_entry['notes']['genomes'] = genome_names
            if 'reaction' in record_genomes:
                for modelseed_reaction_id in ko.reactions:
                    try:
                        reaction_genomes[modelseed_reaction_id] += genome_names
                    except KeyError:
                        reaction_genomes[modelseed_reaction_id] = genome_names
            if 'metabolite' in record_genomes:
                for reaction in ko.reactions.values():
                    for compartment, metabolite in zip(reaction.compartments, reaction.compounds):
                        entry_id = f"{metabolite.modelseed_id}_{compartment}"
                        try:
                            metabolite_genomes[entry_id] += genome_names
                        except KeyError:
                            metabolite_genomes[entry_id] = genome_names

        progress.update("Reactions")
        compound_compartments: Dict[str, Set[str]] = {}
        for modelseed_reaction_id, reaction in self.reactions.items():
            reaction_entry = JSONStructure.get_reaction_entry()
            json_reactions.append(reaction_entry)
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
            reaction_entry['gene_reaction_rule'] = " or ".join([gcid for gcid in reaction_gene_clusters[modelseed_reaction_id]])
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
            if 'reaction' not in record_genomes:
                continue
            notes['genomes'] = sorted(set(reaction_genomes[modelseed_reaction_id]))

        progress.update("Metabolites")
        for modelseed_compound_id, metabolite in self.metabolites.items():
            modelseed_compound_name = metabolite.modelseed_name
            charge = metabolite.charge
            formula = metabolite.formula
            kegg_compound_aliases = list(metabolite.kegg_aliases)
            for compartment in compound_compartments[modelseed_compound_id]:
                metabolite_entry = JSONStructure.get_metabolite_entry()
                json_metabolites.append(metabolite_entry)
                entry_id = f"{modelseed_compound_id}_{compartment}"
                metabolite_entry['id'] = entry_id
                metabolite_entry['name'] = modelseed_compound_name
                metabolite_entry['compartment'] = compartment
                # Compounds without a formula have a nominal charge of 10000000 in the ModelSEED
                # compounds database, which is replaced by None in the reaction network and 0 in the JSON.
                metabolite_entry['charge'] = charge if charge is not None else 0
                metabolite_entry['formula'] = formula if formula is not None else ""
                metabolite_entry['annotation']['kegg.compound'] = kegg_compound_aliases
                if 'metabolite' not in record_genomes:
                    continue
                notes['genomes'] = sorted(set(metabolite_genomes[entry_id]))

        progress.update("Saving")
        with open(path, 'w') as f:
            json.dump(json_dict, f, indent=indent)
        progress.end()

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
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/misc/KEGG/KO_REACTION_NETWORK')
    expected_files = ['ko_info.txt', 'ko_data.tsv']

    def __init__(self, ko_dir: str = None) -> None:
        """
        Load the table derived from downloaded KEGG KO entry files that relates KOs to KEGG
        reactions and EC numbers.

        Parameters
        ==========
        ko_dir : str, None
            The directory containing reference KEGG Orthology (KO) tables set up by anvi'o. The
            default argument of None expects KO data to be set up in the default anvi'o directory
            used by the program `anvi-setup-kegg-data`.
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

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
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
                os.makedirs(parent_dir)
        if os.path.exists(ko_dir):
            if reset:
                shutil.rmtree(ko_dir)
            else:
                raise ConfigError(
                    f"The KO database directory, '{ko_dir}', already exists. 'reset' can be used "
                    "to remove the database at this location and set it up again."
                )
        os.makedirs(ko_dir)

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

class ModelSEEDDatabase:
    """
    The ModelSEED Biochemistry database set up by anvi'o.

    By default, the database is loaded from a default directory of ModelSEED files unless an
    alternative directory is provided.
    """
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/misc/MODELSEED')

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
                f"""\
                No required file named 'sha.txt' was found in the ModelSEED directory,\
                '{modelseed_dir}'.\
                """
            )
        reactions_path = os.path.join(modelseed_dir, 'reactions.tsv')
        if not os.path.isfile(reactions_path):
            raise ConfigError(
                f"""\
                No required file named 'reactions.tsv' was found in the ModelSEED directory,\
                '{modelseed_dir}'.\
                """
            )
        compounds_path = os.path.join(modelseed_dir, 'compounds.tsv')
        if not os.path.isfile(compounds_path):
            raise ConfigError(
                f"""\
                No required file named 'compounds.tsv' was found in the ModelSEED directory,\
                '{modelseed_dir}'.\
                """
            )

        with open(sha_path) as f:
            self.sha = f.read().strip()
        reactions_table = pd.read_csv(reactions_path, sep='\t', header=0, low_memory=False)
        self.compounds_table: pd.DataFrame = pd.read_csv(
            compounds_path,
            sep='\t',
            header=0,
            index_col='id',
            low_memory=False
        )

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

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
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
                    f"""\
                    The ModelSEED database directory, '{modelseed_dir}', already exists. 'reset'\
                    can be used to remove the database at this location and set it up again.\
                    """
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
                            f"""\
                            The connection was reset by the peer more than {max_num_tries} times,\
                            the maximum number of attempts. Try setting up the ModelSEED database\
                            again.\
                            """
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
        reactions_path = os.path.join(
            modelseed_dir, f'ModelSEEDDatabase-{sha}', 'Biochemistry', 'reactions.tsv'
        )
        compounds_path = os.path.join(
            modelseed_dir, f'ModelSEEDDatabase-{sha}', 'Biochemistry', 'compounds.tsv'
        )
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
    """Make, store, and load metabolic reaction networks."""
    def __init__(
        self,
        ko_dir: str = None,
        modelseed_dir: str = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Parameters
        ==========
        ko_dir : str, None
            The directory containing reference KEGG Orthology (KO) tables set up by anvi'o. The
            default argument of None expects KO data to be set up in the default anvi'o directory
            used by the program `anvi-setup-kegg-data`.

        modelseed_dir : str, None
            The directory containing reference ModelSEED Biochemistry tables set up by anvi'o. The
            default argument of None expects ModelSEED data to be set up in the default anvi'o
            directory used by the program `anvi-setup-modelseed-database`.

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
        """
        self.ko_dir = ko_dir
        self.modelseed_dir = modelseed_dir
        self.run = run
        self.progress = progress

    def load_network(
        self,
        contigs_db: str = None,
        pan_db: str = None,
        genomes_storage_db: str = None,
        check_gene_annotations: bool = True,
        load_protein_abundances: bool = False,
        load_metabolite_abundances: bool = False,
        profile_db: str = None,
        quiet: bool = False,
        stats_file: str = None
    ) -> ReactionNetwork:
        """
        Load a reaction network stored in a database as a reaction network object.

        Parameters
        ==========
        contigs_db : str, None
            Path to a contigs database in which a reaction network is stored.

        pan_db : str, None
            Path to a pan database in which a reaction network is stored. 'genomes_storage_db' is
            also required.

        genomes_storage_db : str, None
            Path to a genomes storage database in which KO annotations are stored. 'pan_db' is also
            required.

        check_gene_annotations : bool, True
            If True, as by default, check that the stored reaction network was made from the set of
            gene KO annotations that is currently stored. An exception is raised if this is not the
            case. If False, allow the stored reaction network to have been made from a different set
            of gene KO annotations than is currently stored. This can result in different KOs in the
            returned ReactionNetwork than in the original network that was stored.

        load_protein_abundances : bool, False
            If loading the network from a contigs database, also load abundance measurements of
            proteins that can be expressed by genes in the network. 'profile_db' is also required,
            as abundance profile data is stored there.

        load_metabolite_abundances : bool, False
            If loading the network from a contigs database, also load stored abundance measurements
            of metabolites found in the network. 'profile_db' is also required, as abundance profile
            data is stored there.

        profile_db : str, None
            If loading protein or metabolite abundance data, this database is required, as abundance
            profile data is stored there.

        quiet : bool, False
            Print network overview statistics to the terminal if False.

        stats_file : str, None
            Write network overview statistics to a tab-delimited file at this output path.

        Returns
        =======
        ReactionNetwork
            Reaction network loaded from the input database.
        """
        # Check that the reaction network stored in a database is derived from the current gene KO
        # annotations in the database.
        if contigs_db:
            network = self.load_contigs_database_network(
                contigs_db,
                check_gene_annotations=check_gene_annotations,
                load_protein_abundances=load_protein_abundances,
                load_metabolite_abundances=load_metabolite_abundances,
                profile_db=profile_db,
                quiet=quiet,
                stats_file=stats_file
            )
        elif genomes_storage_db or pan_db:
            network = self.load_pan_database_network(
                genomes_storage_db=genomes_storage_db,
                pan_db=pan_db,
                check_gene_annotations=check_gene_annotations,
                quiet=quiet,
                stats_file=stats_file
            )
        else:
            raise ConfigError(
                f"""\
                A reaction network must be loaded from a database source. Either a contigs database\
                or a genomes storage database and pan database are required.\
                """
            )
        return network

    def load_contigs_database_network(
        self,
        contigs_db: str,
        check_gene_annotations: bool = True,
        load_protein_abundances: bool = False,
        load_metabolite_abundances: bool = False,
        profile_db: str = None,
        quiet: bool = False,
        stats_file: str = None
    ) -> GenomicNetwork:
        """
        Load reaction network data stored in a contigs database as a reaction network object.

        Parameters
        ==========
        contigs_db : str
            Path to a contigs database in which a reaction network is stored.

        check_gene_annotations : bool, True
            If True, as by default, check that the reaction network stored in the contigs database
            was made from the same set of gene KO annotations as currently in the database, and
            throw an error if this is not the case. If False, allow the stored reaction network to
            have been made from a different set of gene KO annotations than is currently stored in
            the database. This can result in different KO assignments to genes in the returned
            GenomicNetwork than in the original network that was stored.

        load_protein_abundances : bool, False
            Load stored abundance measurements of proteins that can be expressed by genes in the
            network. 'profile_db' is also required, as abundance profile data is stored there.

        load_metabolite_abundances : bool, False
            Load stored abundance measurements of metabolites found in the network. 'profile_db' is
            also required, as abundance profile data is stored there.

        profile_db : str, None
            If loading protein or metabolite abundance data, this database is required, as abundance
            profile data is stored there.

        quiet : bool, False
            Print network overview statistics to the terminal if False.

        stats_file : str, None
            Write network overview statistics to a tab-delimited file at this output path.

        Returns
        =======
        GenomicNetwork
            Reaction network loaded from the contigs database.
        """
        # Preemptively check the statistics file path.
        if stats_file is not None:
            filesnpaths.is_output_file_writable(stats_file)

        # Load the contigs database.
        utils.is_contigs_db(contigs_db)
        cdb = ContigsDatabase(contigs_db)
        cdb_db: DB = cdb.db
        sources: List[str] = cdb.meta['gene_function_sources']
        if not sources or not 'KOfam' in sources:
            raise ConfigError(
                f"""\
                The contigs database indicates that genes were never annotated with KOs. This is\
                especially strange since to load a reaction network means that a network had to be\
                constructed from gene KO annotations in the database.\
                """
            )

        # Check that the network stored in the contigs database was made from the same set of KO
        # gene annotations as currently in the database.
        stored_hash = cdb_db.get_meta_value('reaction_network_ko_annotations_hash')
        gene_ko_hits_table = cdb_db.get_table_as_dataframe(
            'gene_functions',
            where_clause='source = "KOfam"',
            columns_of_interest=['gene_callers_id', 'accession', 'function', 'e_value']
        )
        current_hash = self.hash_contigs_db_ko_hits(gene_ko_hits_table)
        if stored_hash != current_hash:
            if check_gene_annotations:
                raise ConfigError(
                    f"""\
                    The reaction network stored in the contigs database was made from a different\
                    set of KEGG KO gene annotations than is currently in the database. There are\
                    two solutions to this problem. First, 'anvi-reaction-network' can be run again\
                    to overwrite the existing network stored in the database with a new network\
                    from the new KO gene annotations. Second, 'check_gene_annotations' can be made\
                    False rather than True, allowing the stored network to have been made from a\
                    different set of KO gene annotations than is currently stored in the database.\
                    This can result in different genes being associated with KOs in the returned\
                    GenomicNetwork than in the original network that was stored. The available\
                    version of the KO database that has been set up by anvi'o is used to fill in\
                    data for KOs in the network that are not current gene annotations.\
                    """
                )
            self.run.warning(
                f"""\
                The reaction network stored in the contigs database was made from a different set of
                KEGG KO gene annotations than is currently in the database. This will be ignored
                since 'check_gene_annotations' is False. This can result in different genes being
                associated with KOs in the returned GenomicNetwork than in the original network that
                was stored.\
                """
            )

        network = GenomicNetwork(run=self.run, progress=self.progress)
        network.contigs_db_source_path = os.path.abspath(contigs_db)

        # Make objects representing all genes with KO annotations in the contigs database, including
        # genes that are not in the network, which are later removed from the network.
        for row in gene_ko_hits_table.itertuples(index=False):
            gcid = int(row.gene_callers_id)
            ko_id = row.accession
            ko_name = row.function
            e_value = float(row.e_value)

            try:
                # This is not the first annotation involving the gene, so an object for it already
                # exists.
                gene = network.genes[gcid]
            except KeyError:
                gene = Gene()
                gene.gcid = gcid
                network.genes[gcid] = gene
            try:
                # This is not the first annotation involving the KO, so an object for it already
                # exists.
                ko = network.kos[ko_id]
            except KeyError:
                ko = KO()
                ko.id = ko_id
                ko.name = ko_name
                network.kos[ko_id] = ko
            gene.kos[ko_id] = ko
            gene.e_values[ko_id] = e_value

        self._load_modelseed_reactions(cdb, network)
        self._load_modelseed_compounds(cdb, network)

        # Remove any trace of genes that do not contribute to the reaction network. Also remove
        # unnetworked KO links to genes.
        unnetworked_gcids = []
        for gcid, gene in network.genes.items():
            gene_in_network = False
            unnetworked_kos: List[str] = []
            for ko_id, ko in gene.kos.items():
                if ko.reactions:
                    gene_in_network = True
                else:
                    unnetworked_kos.append(ko_id)
            if gene_in_network:
                for unnetworked_ko_id in unnetworked_kos:
                    gene.kos.pop(unnetworked_ko_id)
                    gene.e_values.pop(unnetworked_ko_id)
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

        # Remove entries in the network attribute mapping ModelSEED reaction IDs to KO KEGG
        # REACTION ID aliases if no such aliases were found to exist.
        modelseed_reaction_ids = []
        for modelseed_reaction_id, kegg_reaction_ids in network.modelseed_kegg_aliases.items():
            if not kegg_reaction_ids:
                modelseed_reaction_ids.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids:
            network.modelseed_kegg_aliases.pop(modelseed_reaction_id)

        # Remove entries in the network attribute mapping ModelSEED reaction IDs to KO EC number
        # aliases if no such aliases were found to exist.
        modelseed_reaction_ids = []
        for modelseed_reaction_id, ec_numbers in network.modelseed_ec_number_aliases.items():
            if not ec_numbers:
                modelseed_reaction_ids.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids:
            network.modelseed_ec_number_aliases.pop(modelseed_reaction_id)

        if load_protein_abundances or load_metabolite_abundances:
            network.profile_db_source_path = os.path.abspath(profile_db)
            pdb = ProfileDatabase(profile_db)
            if load_protein_abundances:
                self._load_protein_abundances(pdb, cdb, network)
            if load_metabolite_abundances:
                self._load_metabolite_abundances(pdb, network)
            pdb.disconnect()

        if quiet and not stats_file:
            return network

        precomputed_counts = {
            'total_genes': cdb_db.get_row_counts_from_table('genes_in_contigs'),
            'genes_assigned_kos': len(network.genes) + len(unnetworked_gcids),
            'kos_assigned_genes': len(network.kos) + len(unnetworked_ko_ids)
        }
        cdb.disconnect()
        stats = network.get_overview_statistics(precomputed_counts=precomputed_counts)
        if not quiet:
            network.print_overview_statistics(stats=stats)
        if stats_file:
            network.write_overview_statistics(stats_file, stats=stats)

        return network

    def _load_protein_abundances(
        self,
        profile_database: ProfileDatabase,
        contigs_database: ContigsDatabase,
        network: GenomicNetwork
    ) -> None:
        """
        Load abundance data for proteins that can be expressed by genes in the metabolic network.

        Protein isoforms are not supported.

        Parameters
        ==========
        profile_database : ProfileDatabase
            The database storing protein measurements.

        contigs_database : ContigsDatabase
            The database storing associations between genes and proteins.

        network : GenomicNetwork
            The genomic network under construction.

        Returns
        =======
        None
        """
        protein_abundances_table = profile_database.db.get_table_as_dataframe(
            tables.protein_abundances_table_name
        )
        if len(protein_abundances_table) == 0:
            return

        gene_functions_table = contigs_database.db.get_table_as_dataframe(
            'gene_functions', columns_of_interest=['gene_callers_id', 'source', 'accession']
        )
        gene_functions_table = gene_functions_table[
            gene_functions_table['gene_callers_id'].isin(network.genes)
        ]
        gene_functions_table = gene_functions_table.rename(
            {'source': 'reference_source', 'accession': 'reference_id'}, axis=1
        )

        protein_abundances_table = protein_abundances_table.merge(
            gene_functions_table, how='inner', on=['reference_source', 'reference_id']
        )

        multiprotein_genes: Dict[int, List[int]] = {}
        for key, protein_table in protein_abundances_table.groupby(
            ['protein_id', 'reference_source', 'reference_id']
        ):
            protein_id = key[0]
            protein = Protein()
            network.proteins[protein_id] = protein
            protein.id = protein_id
            for gcid in protein_table['gene_callers_id'].unique():
                gene = network.genes[gcid]
                protein.genes[gcid] = gene
                if gene.protein:
                    try:
                        multiprotein_genes[gene].append(protein_id)
                    except KeyError:
                        multiprotein_genes[gene] = [protein_id]
                else:
                    gene.protein = protein
            for row in protein_table.itertuples():
                protein.abundances[row.sample_name] = row.abundance_value

        if multiprotein_genes:
            s = ""
            for gcid, protein_ids in multiprotein_genes.items():
                s += f"{gcid}: {', '.join(protein_ids)}; "
            s = s[: -1]
            raise ConfigError(
                f"""\
                Certain genes were unexpectedly associated with multiple proteins with abundance\
                data. Unfortunately, multiple protein products are not currently allowed in anvi'o,\
                so the protein abundance data must be edited down in the profile database to permit\
                use with the reaction network. These are as follows, with the gene callers ID\
                separated by a comma-separated\ list of protein IDs. {s}\
                """
            )

    def _load_metabolite_abundances(
        self,
        profile_database: ProfileDatabase,
        network: GenomicNetwork
    ) -> None:
        """
        Load abundance data for metabolites represented in the metabolic network.

        Parameters
        ==========
        profile_database : ProfileDatabase
            The database storing protein measurement data that is loaded into the genomic network.

        network : GenomicNetwork
            The genomic network under construction.

        Returns
        =======
        None
        """
        metabolite_abundances_table = profile_database.db.get_table_as_dataframe(
            tables.metabolite_abundances_table_name
        )
        metabolite_abundances_table = metabolite_abundances_table[
            metabolite_abundances_table['reference_id'].isin(network.metabolites)
        ]
        if len(metabolite_abundances_table) == 0:
            return

        for modelseed_compound_id, metabolite_table in metabolite_abundances_table.groupby(
            'reference_id'
        ):
            metabolite = network.metabolites[modelseed_compound_id]
            for row in metabolite_table.itertuples():
                metabolite.abundances[row.sample_name] = row.abundance_value

    def load_pan_database_network(
        self,
        pan_db: str,
        genomes_storage_db: str,
        check_gene_annotations: bool = True,
        quiet: bool = False,
        stats_file: str = None
    ) -> PangenomicNetwork:
        """
        Load reaction network data stored in a pan database as a reaction network object.

        Parameters
        ==========
        pan_db : str
            Path to the pan database in which the reaction network is stored.

        genomes_storage_db : str
            Path to the genomes storage database associated with the pan database.

        check_annotations : bool, True
            If True, as by default, check that the reaction network stored in the pan database was
            made from the set of gene KO annotations currently stored in the associated genomes
            storage database. An exception is raised if this is not the case. If False, allow the
            stored reaction network to have been made from a different set of gene KO annotations
            than is currently stored in the genomes storage database. This can result in different
            consensus KOs assigned to gene clusters in the returned PangenomicNetwork than in the
            original network that was stored.

        quiet : bool, False
            Print network overview statistics to the terminal if False.

        stats_file : str, None
            Write network overview statistics to a tab-delimited file at this output path.

        Returns
        =======
        PangenomicNetwork
            The network derived from the pangenomic databases.
        """
        # Preemptively check the statistics file path.
        if stats_file is not None:
            filesnpaths.is_output_file_writable(stats_file)

        # Load the pan database.
        pan_db_info = dbinfo.PanDBInfo(pan_db)
        self_table = pan_db_info.get_self_table()
        # No consensus threshold may have been used in network construction, in which case the value
        # of the parameter is None.
        consensus_threshold = self_table['reaction_network_consensus_threshold']
        if consensus_threshold is not None:
            consensus_threshold = float(consensus_threshold)
        discard_ties = bool(int(self_table['reaction_network_discard_ties']))
        args = argparse.Namespace()
        args.pan_db = pan_db
        args.genomes_storage = genomes_storage_db
        args.consensus_threshold = consensus_threshold
        args.discard_ties = discard_ties
        pan_super = PanSuperclass(args, r=run_quiet)
        pan_super.init_gene_clusters()
        pan_super.init_gene_clusters_functions()
        pan_super.init_gene_clusters_functions_summary_dict()
        gene_clusters_functions_summary_dict: Dict = pan_super.gene_clusters_functions_summary_dict

        # Check that the network stored in the pan database was made from the same set of KO gene
        # annotations currently in the associated genomes storage database.
        stored_hash = self_table['reaction_network_ko_annotations_hash']
        current_hash = self.hash_pan_db_ko_annotations(
            genomes_storage_db,
            gene_clusters_functions_summary_dict,
            consensus_threshold,
            discard_ties
        )
        if stored_hash != current_hash:
            if check_gene_annotations:
                # Note that another unstated possible cause of the error could be due to manual
                # meddling with the metavariables, 'consensus_threshold' and 'discard_ties', in the
                # database. Assume that the user was not engaged in mischief.
                raise ConfigError(
                    f"""\
                    The reaction network stored in the pan database was made from a different set\
                    of KO gene annotations than is currently in the associated genomes storage\
                    database. There are two solutions to this problem. First, the program,\
                    'anvi-reaction-network', can be run again to overwrite the existing network\
                    stored in the pan database with a new network from the new KO gene annotations.\
                    Second, 'check_gene_annotations' can be given an argument of False instead of\
                    True, preventing this exception from being raised if the stored network was\
                    made from a different set of KO gene annotations than is currently in the\
                    genomes storage database. This can result in different consensus KOs assigned\
                    to gene clusters in the returned PangenomicNetwork than in the original network\
                    that was stored. The available version of the KO database that has been set up\
                    by anvi'o is used to fill in data for any KOs in the network that are not\
                    current gene annotations in the genomes storage database.\
                    """
                )
            self.run.warning(
                f"""\
                The reaction network stored in the pan database was made from a different set of KO
                gene annotations than is currently in the genomes storage database. This will be
                ignored since 'check_gene_annotations' is False. This can result in different
                consensus KO assignments to gene clusters in the returned PangenomicNetwork than in
                the original network that was stored.\
                """
            )

        # Create the reaction network object.
        network = PangenomicNetwork(run=self.run, progress=self.progress)
        network.pan_db_source_path = os.path.abspath(pan_db)
        network.genomes_storage_db_source_path = os.path.abspath(genomes_storage_db)
        network.consensus_threshold = consensus_threshold
        network.discard_ties = discard_ties
        if stored_hash == current_hash:
            network.consistent_annotations = True
        else:
            network.consistent_annotations = False

        # Make objects representing all gene clusters with consensus KO annotations.
        num_gene_clusters_assigned_ko = 0
        for cluster_id, gene_cluster_functions_data in gene_clusters_functions_summary_dict.items():
            # Retrieve the consensus KO across genes in the cluster. Parameterization of the method
            # used to select consensus KOs occurred in pan super initialization. Parameter values
            # were loaded from pan database metavariables.
            gene_cluster_ko_data = gene_cluster_functions_data['KOfam']
            if gene_cluster_ko_data == {'function': None, 'accession': None}:
                # No KO was assigned to the cluster.
                continue
            ko_id = gene_cluster_ko_data['accession']
            num_gene_clusters_assigned_ko += 1

            gene_cluster = GeneCluster()
            gene_cluster.gene_cluster_id = cluster_id
            gene_cluster.genomes = list(pan_super.gene_clusters[cluster_id])
            # Add the gene cluster to the network, regardless of whether it yields reactions. Gene
            # clusters not contributing to the reaction network are removed later.
            network.gene_clusters[cluster_id] = gene_cluster

            try:
                # This is not the first gene cluster that has been encountered with the KO assigned
                # to it, so an object for the KO already exists.
                ko = network.kos[ko_id]
            except KeyError:
                ko = KO()
                ko.id = ko_id
                ko.name = gene_cluster_ko_data['function']
                network.kos[ko_id] = ko
            gene_cluster.ko = ko

        pdb = PanDatabase(pan_db)
        self._load_modelseed_reactions(pdb, network)
        self._load_modelseed_compounds(pdb, network)

        # Remove any trace of gene clusters that do not contribute to the reaction network.
        unnetworked_cluster_ids = []
        for cluster_id, gene_cluster in network.gene_clusters.items():
            if gene_cluster.ko.reactions:
                continue
            unnetworked_cluster_ids.append(cluster_id)
        for cluster_id in unnetworked_cluster_ids:
            network.gene_clusters.pop(cluster_id)

        # Remove any trace of KOs that do not contribute to the reaction network.
        unnetworked_ko_ids = []
        for ko_id, ko in network.kos.items():
            if not ko.reactions:
                unnetworked_ko_ids.append(ko_id)
        for ko_id in unnetworked_ko_ids:
            network.kos.pop(ko_id)

        # Remove entries in the network attribute mapping ModelSEED reaction IDs to KO KEGG
        # REACTION ID aliases if no such aliases were found to exist.
        modelseed_reaction_ids = []
        for modelseed_reaction_id, kegg_reaction_ids in network.modelseed_kegg_aliases.items():
            if not kegg_reaction_ids:
                modelseed_reaction_ids.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids:
            network.modelseed_kegg_aliases.pop(modelseed_reaction_id)

        # Remove entries in the network attribute mapping ModelSEED reaction IDs to KO EC number
        # aliases if no such aliases were found to exist.
        modelseed_reaction_ids = []
        for modelseed_reaction_id, ec_numbers in network.modelseed_ec_number_aliases.items():
            if not ec_numbers:
                modelseed_reaction_ids.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids:
            network.modelseed_ec_number_aliases.pop(modelseed_reaction_id)

        if quiet and not stats_file:
            return network

        if network.consistent_annotations:
            precomputed_counts = {
                'total_gene_clusters': pdb.meta['num_gene_clusters'],
                'gene_clusters_assigned_ko': num_gene_clusters_assigned_ko,
                'kos_assigned_gene_clusters': len(network.kos) + len(unnetworked_ko_ids)
            }
        else:
            precomputed_counts = {
                'total_gene_clusters': pdb.meta['num_gene_clusters'],
                'gene_clusters_assigned_ko': None,
                'kos_assigned_gene_clusters': None
            }
        pdb.disconnect()
        stats = network.get_overview_statistics(precomputed_counts=precomputed_counts)
        if not quiet:
            network.print_overview_statistics(stats=stats)
        if stats_file:
            network.write_overview_statistics(stats_file, stats=stats)

        return network

    def _load_modelseed_reactions(
        self,
        database: Union[ContigsDatabase, PanDatabase],
        network: Union[GenomicNetwork, PangenomicNetwork]
    ) -> None:
        """
        Add ModelSEED reactions to the network being loaded from the database.

        ModelSEED reaction objects are related to KOs through KEGG REACTION and EC number aliases.

        Parameters
        ==========
        database : ContigsDatabase or PanDatabase
            The database storing a reaction network. In loading a genomic network, provide a contigs
            database; in loading a pangenomic network, provide a pan database.

        network : GenomicNetwork or PangenomicNetwork
            The reaction network under construction.

        Returns
        =======
        None
        """
        # Load the table of reactions data.
        if type(database) is ContigsDatabase:
            reactions_table = database.db.get_table_as_dataframe('gene_function_reactions')
            if type(network) is not GenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'ContigsDatabase', so the provided "
                    "'network' must be of type 'GenomicNetwork'. Instead, the reaction network "
                    f"argument was of type '{type(network)}'."
                )
        elif type(database) is PanDatabase:
            reactions_table = database.db.get_table_as_dataframe('gene_cluster_function_reactions')
            if type(network) is not PangenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'PanDatabase', so the provided 'network' "
                    "must be of type 'PangenomicNetwork'. Instead, the reaction network argument "
                    f"was of type '{type(network)}'."
                )
        else:
            raise ConfigError(
                "The provided 'database' must be of type 'ContigsDatabase' or 'PanDatabase'. "
                f"Instead, the argument was of type '{type(database)}'."
            )

        # The KO database is needed if KOs in the stored network aren't among the current gene
        # annotations.
        try:
            ko_db = KODatabase(ko_dir=self.ko_dir)
        except ConfigError as e:
            raise ConfigError(
                f"{e} Please set up the KO database in the default directory with the program, "
                "'anvi-reaction-network'."
            )

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
                    # This is not the first reaction involving the compound, so an object for it
                    # already exists.
                    compound = network.metabolites[modelseed_compound_id]
                except KeyError:
                    compound = ModelSEEDCompound()
                    compound.modelseed_id = modelseed_compound_id
                    network.metabolites[modelseed_compound_id] = compound
                reaction.compounds.append(compound)
            reaction.compounds = tuple(reaction.compounds)

            stoichiometry: str = row.stoichiometry
            reaction.coefficients = tuple(int(coeff) for coeff in stoichiometry.split(', '))
            compartments: str = row.compartments
            reaction.compartments = tuple(compartments.split(', '))
            reversibility: int = row.reversibility
            reaction.reversibility = bool(reversibility)

            # Map KEGG reaction aliases of the ModelSEED reaction to all KOs that were associated
            # with the KEGG reaction.
            kegg_reaction_ko_ids: Dict[str, List[str]] = {}
            kegg_reaction_sources: str = row.ko_kegg_reaction_source
            for kegg_reaction_item in kegg_reaction_sources.split('; '):
                if not kegg_reaction_item:
                    # The ModelSEED reaction was not sourced from KEGG reactions.
                    continue
                kegg_reaction_id, ko_ids = kegg_reaction_item.split(': (')
                ko_ids = ko_ids[:-1].split(', ')
                kegg_reaction_ko_ids[kegg_reaction_id] = ko_ids
            # Record *all* KEGG reaction aliases of the ModelSEED reaction, including those not
            # associated with KO annotations.
            other_kegg_reaction_ids: str = row.other_kegg_reaction_ids
            reaction.kegg_aliases = list(kegg_reaction_ko_ids)
            if other_kegg_reaction_ids:
                reaction.kegg_aliases += other_kegg_reaction_ids.split(', ')
            reaction.kegg_aliases = tuple(reaction.kegg_aliases)

            network.modelseed_kegg_aliases[modelseed_reaction_id] = modelseed_kegg_aliases = []
            orphan_ko_ids = []
            for kegg_reaction_id, ko_ids in kegg_reaction_ko_ids.items():
                # Record the ModelSEED reaction as one of the aliases of the KEGG reaction in the
                # network.
                try:
                    network.kegg_modelseed_aliases[kegg_reaction_id].append(modelseed_reaction_id)
                except KeyError:
                    network.kegg_modelseed_aliases[kegg_reaction_id] = [modelseed_reaction_id]
                modelseed_kegg_aliases.append(kegg_reaction_id)
                for ko_id in ko_ids:
                    try:
                        ko = network.kos[ko_id]
                    except KeyError:
                        # In the case of a genomic network, this error arises when the current set
                        # of gene KO annotations in the contigs database does not match the set from
                        # which the reaction network was originally made, and the KO under
                        # consideration in the network is no longer a gene annotation in the
                        # database. In the case of a pangenomic network, this error arises when the
                        # current set of gene cluster consensus KO annotations does not match the
                        # set from which the reaction network was originally made and the consensus
                        # KO under consideration in the network no longer annotates a gene cluster
                        # in the pan database. (The current set of gene cluster consensus KO
                        # annotations is derived from the pan and genomes storage databases using
                        # the parameters, 'consensus_threshold' and 'discard_ties'.)
                        ko = KO()
                        ko.ko_id = ko_id
                        # The KO name is unknown from the database, so take it from the KO database.
                        ko.ko_name = ko_db.ko_table.loc[ko_id, 'name']
                        network.kos[ko_id] = ko
                        orphan_ko_ids.append(ko_id)
                    if not modelseed_reaction_id in ko.reactions:
                        # This is the first time encountering the reaction as a reference of the KO.
                        ko.reactions[modelseed_reaction_id] = reaction
                    try:
                        ko.kegg_reaction_aliases[modelseed_reaction_id].append(kegg_reaction_id)
                    except KeyError:
                        ko.kegg_reaction_aliases[modelseed_reaction_id] = [kegg_reaction_id]

            # Map EC number aliases of the ModelSEED reaction to all KOs that were associated with
            # the EC number.
            ec_number_ko_ids: Dict[str, List[str]] = {}
            ec_number_sources: str = row.ko_ec_number_source
            for ec_number_item in ec_number_sources.split('; '):
                if not ec_number_item:
                    # The ModelSEED reaction was not sourced from EC numbers.
                    continue
                ec_number, ko_ids = ec_number_item.split(': (')
                ko_ids = ko_ids[:-1].split(', ')
                ec_number_ko_ids[ec_number] = ko_ids
            # Record *all* EC number aliases of the ModelSEED reaction, including those not
            # associated with KO annotations.
            other_ec_numbers: str = row.other_ec_numbers
            reaction.ec_number_aliases = list(ec_number_ko_ids)
            if other_ec_numbers:
                reaction.ec_number_aliases += other_ec_numbers.split(', ')
            reaction.ec_number_aliases = tuple(reaction.ec_number_aliases)

            modelseed_ec_number_aliases = []
            network.modelseed_ec_number_aliases[modelseed_reaction_id] = modelseed_ec_number_aliases
            for ec_number, ko_ids in ec_number_ko_ids.items():
                # Record the ModelSEED reaction as one of the aliases of the EC number in the
                # network.
                try:
                    network.ec_number_modelseed_aliases[ec_number].append(modelseed_reaction_id)
                except KeyError:
                    network.ec_number_modelseed_aliases[ec_number] = [modelseed_reaction_id]
                modelseed_ec_number_aliases.append(ec_number)
                for ko_id in ko_ids:
                    try:
                        ko = network.kos[ko_id]
                    except KeyError:
                        # This error arises for the same reason as before (in processing KEGG
                        # reactions).
                        ko = KO()
                        ko.ko_id = ko_id
                        # The KO name is unknown from the database, so take it from the KO database.
                        ko.ko_name = ko_db.ko_table.loc[ko_id, 'name']
                        network.kos[ko_id] = ko
                        orphan_ko_ids.append(ko_id)
                    if not modelseed_reaction_id in ko.reactions:
                        # This is the first time encountering the reaction as a reference of the KO.
                        ko.reactions[modelseed_reaction_id] = reaction
                    try:
                        ko.ec_number_aliases[modelseed_reaction_id].append(ec_number)
                    except KeyError:
                        ko.ec_number_aliases[modelseed_reaction_id] = [ec_number]

            if DEBUG:
                # "Orphan" KOs can only arise when 'check_gene_annotations' is False in the calling
                # method, 'load_contigs_database_network' or 'load_pan_database_network'.
                if type(network) is GenomicNetwork:
                    self.run.info_single(
                        "The following KOs are found in the stored reaction network in the contigs "
                        "database, but they are not found among the current gene KO annotations in "
                        "the contigs database. The available version of the KO database set up by "
                        "anvi'o was used to retrieve the function 'names' of these KOs: "
                        f"{', '.join(orphan_ko_ids)}"
                    )
                elif type(network) is PangenomicNetwork:
                    self.run.info_single(
                        "The following KOs are found in the stored reaction network in the pan "
                        "database, but they are not found among the current gene KO annotations in "
                        "the genomes storage database. The available version of the KO database "
                        "set up by anvi'o was used to retrieve the function 'names' of these KOs: "
                        f"{', '.join(orphan_ko_ids)}"
                    )

    def _load_modelseed_compounds(
        self,
        database: Union[ContigsDatabase, PanDatabase],
        network: Union[GenomicNetwork, PangenomicNetwork]
    ) -> None:
        """
        Add ModelSEED compounds to the network being loaded from the database.

        Parameters
        ==========
        database : ContigsDatabase or PanDatabase
            The database storign a reaction network. In loading a genomic network, provide a contigs
            database; in loading a pangenomic network, provide a pan database.

        network : GenomicNetwork or PangenomicNetwork
            The reaction network under construction.

        Returns
        =======
        None
        """
        # Load the table of compounds data.
        if type(database) is ContigsDatabase:
            metabolites_table = database.db.get_table_as_dataframe('gene_function_metabolites')
            if type(network) is not GenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'ContigsDatabase', so the provided "
                    "'network' must be of type 'GenomicNetwork'. Instead, the reaction network "
                    f"argument was of type '{type(network)}'."
                )
        elif type(database) is PanDatabase:
            metabolites_table = database.db.get_table_as_dataframe('gene_cluster_function_metabolites')
            if type(network) is not PangenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'PanDatabase', so the provided 'network' "
                    "must be of type 'PangenomicNetwork'. Instead, the reaction network argument "
                    f"was of type '{type(database)}'."
                )
        else:
            raise ConfigError(
                "The provided 'database' must be of type 'ContigsDatabase' or 'PanDatabase'. "
                f"Instead, the argument was of type '{type(database)}'."
            )

        for row in metabolites_table.itertuples():
            # Each row of the table contains information on a different ModelSEED compound.
            modelseed_compound_id = row.modelseed_compound_id
            compound = network.metabolites[modelseed_compound_id]
            modelseed_compound_name: str = row.modelseed_compound_name
            compound.modelseed_name = modelseed_compound_name
            kegg_aliases: str = row.kegg_aliases
            if kegg_aliases:
                compound.kegg_aliases = tuple(kegg_aliases.split(', '))
            else:
                compound.kegg_aliases = tuple()
            # Compounds without a formula, recorded here as None, have a nominal charge of 10000000
            # in the ModelSEED compounds database. This is replaced by NaN in the table and here as
            # None in the reaction network.
            formula: str = row.formula
            compound.formula = formula
            charge: int = row.charge
            compound.charge = None if np.isnan(charge) else int(charge)

    def make_network(
        self,
        contigs_db: str = None,
        pan_db: str = None,
        genomes_storage_db: str = None,
        store: bool = True,
        overwrite_existing_network: bool = False,
        consensus_threshold: float = None,
        discard_ties: bool = False,
        stats_file: str = None
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

        pan_db : str, None
            Path to a pan database. The pangenomic network is determined for gene clusters stored in
            the database. If 'store' is True, the network is saved in the database.
            An argument for the paired 'genomes_storage_db' is also required.

        genomes_storage_db : str, None
            Path to a genomes storage database. The pangenomic network is derived from gene KO
            annotations stored in the database. An argument for the paired 'pan_db' is also
            required.

        store : bool, True
            Save the network. A network constructed from a contigs database is stored in that
            database. A pangenomic network constructed from a genomes stroage database and pan
            database is stored in the pan database.

        overwrite_existing_network : bool, False
            Overwrite an existing network stored in the contigs or pan database. 'store' is also
            required.

        consensus_threshold : float, None
            This parameter applies to pangenomes. With the default of None, the protein annotation
            most frequent among genes in a cluster is assigned to the cluster itself. If a
            non-default argument is provided (a value on [0, 1]), at least this proportion of genes
            in the cluster must have the most frequent annotation for the cluster to be annotated.

        discard_ties : bool, False
            This parameter applies to pangenomes. If multiple protein annotations are most frequent
            among genes in a cluster, then do not assign an annotation to the cluster itself when
            this argument is True. By default, this argument is False, so one of the most frequent
            annotations would be arbitrarily chosen.

        stats_file : str, None
            Write network overview statistics to a tab-delimited file at this output path.

        Returns
        =======
        ReactionNetwork
            Reaction network loaded from the input database.
        """
        if contigs_db and (pan_db or genomes_storage_db):
            raise ConfigError(
                "Either a contigs database OR both a pan database and genomes storage database are required "
                "to make either a (meta)genomic reaction network or a pangenomic reaction network, respectively."
            )
        elif contigs_db:
            self.run.info_single(
                "A reaction network will be made from protein orthology annotations in the contigs database."
            )
            network = self.make_contigs_database_network(
                contigs_db,
                store=store,
                overwrite_existing_network=overwrite_existing_network,
                stats_file=stats_file
            )
        elif genomes_storage_db or pan_db:
            self.run.info_single(
                "A pangenomic reaction network will be made from protein orthology annotations "
                "in the genomes storage database and gene clusters in the pan database."
            )
            network = self.make_pangenomic_network(
                pan_db,
                genomes_storage_db,
                store=store,
                overwrite_existing_network=overwrite_existing_network,
                consensus_threshold=consensus_threshold,
                discard_ties=discard_ties,
                stats_file=stats_file
            )
        else:
            raise ConfigError(
                "A reaction network cannot be made without a database source. Either a contigs database OR "
                "a pan database and genomes storage database are required to make either a (meta)genomic "
                "reaction network or a pangenomic reaction network, respectively."
            )
        return network

    def make_contigs_database_network(
        self,
        contigs_db: str,
        store: bool = True,
        overwrite_existing_network: bool = False,
        stats_file: str = None
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

        stats_file : str, None
            Write network overview statistics to a tab-delimited file at this output path.

        Returns
        =======
        GenomicNetwork
            The network derived from the contigs database.
        """
        # Here is an example of the information used to create a genomic network.
        # gene 1 ---> KO 1 ---> KEGG rxn 1 ---> ModelSEED rxn 1 ---> ModelSEED metabs 1, 2, ...
        # |      |         |
        # |      |         ---> EC number 1 --> ModelSEED rxn 1 ---> ModelSEED metabs 1, 2, ...
        # |      |         |                |
        # |      |         |                --> ModelSEED rxn 2 ---> ...
        # |      |         |
        # |      |         ---> EC number 2 --> ...
        # |      |
        # |      ---> KO 2 ---> ...
        # |
        # gene 2 ---> ...

        # Preemptively check the statistics file path.
        if stats_file is not None:
            filesnpaths.is_output_file_writable(stats_file)

        # Load the contigs database.
        self.run.info("Contigs database", contigs_db)
        utils.is_contigs_db(contigs_db)
        cdb = ContigsDatabase(contigs_db)
        cdb_db: DB = cdb.db
        sources: List[str] = cdb.meta['gene_function_sources']
        if not sources or not 'KOfam' in sources:
            raise ConfigError(
                f"""\
                The contigs database indicates that genes were never annotated with KOs, which is\
                required to build a reaction network. This can be solved by running \
                'anvi-run-kegg-kofams' on the contigs database.\
                """
            )
        if (
            store and
            cdb_db.get_meta_value('reaction_network_ko_annotations_hash') and
            not overwrite_existing_network
        ):
            raise ConfigError(
                f"""\
                The existing reaction network in the contigs database must be explicitly\
                overwritten.\
                """
            )

        self.progress.new("Building reaction network")
        self.progress.update("Loading reference databases")

        network = GenomicNetwork(run=self.run, progress=self.progress)
        network.contigs_db_source_path = os.path.abspath(contigs_db)

        # Load reference databases.
        ko_db = KODatabase(self.ko_dir)
        modelseed_db = ModelSEEDDatabase(self.modelseed_dir)
        modelseed_kegg_reactions_table = modelseed_db.kegg_reactions_table
        modelseed_ec_reactions_table = modelseed_db.ec_reactions_table
        modelseed_compounds_table = modelseed_db.compounds_table

        # Record KOs that annotated genes in the contigs database but for some reason aren't found
        # in the KEGG KO database set up by anvi'o.
        undefined_ko_ids: List[str] = []

        # Record ModelSEED reactions that would have been added to the reaction network if the
        # reaction had a chemical equation in the ModelSEED Biochemistry database.
        undefined_modelseed_reaction_ids: List[str] = []

        # Parse gene-KO matches recorded in the contigs database.
        gene_ko_hits_table = cdb_db.get_table_as_dataframe(
            'gene_functions',
            where_clause='source = "KOfam"',
            columns_of_interest=['gene_callers_id', 'accession', 'function', 'e_value']
        )
        total_ko_matches = len(gene_ko_hits_table)
        num_ko_matches_parsed = -1
        for row in gene_ko_hits_table.itertuples(index=False):
            num_ko_matches_parsed += 1
            self.progress.update(
                f"Gene-KO matches parsed: {num_ko_matches_parsed} / {total_ko_matches}"
            )

            # Get data on the gene-KO match.
            gcid = int(row.gene_callers_id)
            ko_id = row.accession
            ko_name = row.function
            e_value = float(row.e_value)

            # Represent the gene as an object.
            if gcid in network.genes:
                # The gene has already been added to the network.
                gene = network.genes[gcid]
                is_new_gene = False
            else:
                gene = Gene()
                gene.gcid = gcid
                is_new_gene = True

            try:
                # The KO and its associated reactions and metabolites have already been added to the
                # network.
                ko = network.kos[ko_id]
                is_new_ko = False
            except KeyError:
                is_new_ko = True
            if not is_new_ko:
                # Have the gene reference the KO.
                gene.kos[ko_id] = ko
                gene.e_values[ko_id] = e_value
                if is_new_gene:
                    # Add the newly encountered gene to the network.
                    network.genes[gcid] = gene
                # Proceed to the next gene-KO match.
                continue

            # Get KEGG REACTION IDs and EC numbers associated with the KO.
            try:
                ko_info = ko_db.ko_table.loc[ko_id]
            except KeyError:
                # For some reason the KO annotated a gene in the contigs database but is not found
                # in the KEGG KO database set up by anvi'o. Do not add the alien KO or the gene, if
                # newly encountered, to the network.
                undefined_ko_ids.append(ko_id)
                continue
            ko_kegg_reaction_ids = self._get_ko_kegg_reaction_ids(ko_info)
            ko_ec_numbers = self._get_ko_ec_numbers(ko_info)

            if not ko_kegg_reaction_ids and not ko_ec_numbers:
                # The KO is not associated with any KEGG reactions or EC numbers, and therefore
                # anvi'o can't relate the KO to ModelSEED reactions. Do not add the unsystematizable
                # KO or the gene, if newly encountered, to the network.
                continue

            # Check if KEGG reactions and EC numbers associated with the KO have already been added
            # to the network in processing other gene-KO matches. To have been added to the network,
            # KEGG reactions and EC numbers must have aliased ModelSEED reactions.
            old_kegg_reaction_ids, new_kegg_reaction_ids = self._find_kegg_reaction_ids(
                ko_kegg_reaction_ids, network
            )
            old_ec_numbers, new_ec_numbers = self._find_ec_numbers(
                ko_ec_numbers, network
            )

            # Retrieve data on ModelSEED reactions aliasing KEGG reactions that haven't been added
            # to the network. Each row of the reference table represents a unique mapping of KEGG
            # reaction to ModelSEED reaction.
            modelseed_kegg_reactions_dict: Dict[int, Dict] = modelseed_kegg_reactions_table[
                modelseed_kegg_reactions_table['KEGG_REACTION_ID'].isin(new_kegg_reaction_ids)
            ].to_dict(orient='index')

            # Retrieve data on ModelSEED reactions aliasing EC numbers that haven't been added to
            # the network. Each row of the reference table represents a unique mapping of EC number
            # to ModelSEED reaction.
            modelseed_ec_reactions_dict: Dict[int, Dict] = modelseed_ec_reactions_table[
                modelseed_ec_reactions_table['EC_number'].isin(new_ec_numbers)
            ].to_dict(orient='index')

            if not (
                old_kegg_reaction_ids or
                old_ec_numbers or
                modelseed_kegg_reactions_dict or
                modelseed_ec_reactions_dict
            ):
                # None of the KEGG REACTION IDs and EC numbers associated with the KO map to
                # ModelSEED reactions (none are in the ModelSEED Biochemistry reactions table).
                continue

            # Find "undefined" ModelSEED reactions without an equation.
            undefined_modelseed_reaction_ids += self._remove_undefined_reactions(
                ko_kegg_reaction_ids,
                new_kegg_reaction_ids,
                ko_ec_numbers,
                new_ec_numbers,
                modelseed_kegg_reactions_dict,
                modelseed_ec_reactions_dict
            )

            if not (
                old_kegg_reaction_ids or
                new_kegg_reaction_ids or
                old_ec_numbers or
                new_ec_numbers
            ):
                # The KO is not added to the network as none of the reactions have equations.
                continue
            # The newly encountered KO is now known to be associated with KEGG reactions or EC
            # numbers that alias ModelSEED reactions with an equation; this new data can be added to
            # the network.

            if is_new_gene:
                # Add an object representing the newly encountered gene to the network.
                network.genes[gcid] = gene

            # Add an object representing the newly encountered KO to the network.
            ko = KO()
            ko.id = ko_id
            ko.name = ko_name
            network.kos[ko_id] = ko
            gene.kos[ko_id] = ko
            gene.e_values[ko_id] = e_value

            # Associate ModelSEED reactions that have been previously added to the network under
            # construction with the newly encountered KO.
            self._process_added_reactions(
                old_kegg_reaction_ids,
                old_ec_numbers,
                network,
                ko,
                ko_kegg_reaction_ids,
                ko_ec_numbers
            )

            # Add ModelSEED reactions aliasing newly encountered KEGG reactions and EC numbers to
            # the network.
            self._add_reactions(
                modelseed_kegg_reactions_dict,
                modelseed_ec_reactions_dict,
                network,
                modelseed_compounds_table,
                ko,
                old_kegg_reaction_ids,
                new_kegg_reaction_ids,
                old_ec_numbers,
                new_ec_numbers
            )

        if DEBUG:
            for gene in network.genes.values():
                for ko in gene.kos.values():
                    assert ko.reactions

            assert sorted(network.modelseed_kegg_aliases) == sorted(network.reactions)
            assert sorted(network.modelseed_ec_number_aliases) == sorted(network.reactions)

        undefined_modelseed_reaction_ids = set(undefined_modelseed_reaction_ids)
        undefined_ko_ids = set(undefined_ko_ids)

        self.progress.end()

        if DEBUG:
            self.run.info_single(
                f"""\
                The following ModelSEED reactions would have been added to the reaction network had
                there been a chemical equation in the ModelSEED database; perhaps it is worth
                investigating the ModelSEED reactions table to understand why this is not the case:
                {', '.join(undefined_modelseed_reaction_ids)}\
                """
            )

        if undefined_ko_ids:
            self.run.info_single(
                f"""\
                Certain genes matched KOs that were not found in the reference KO database. These
                KOs will not be used in network construction. It could be that the KOfams used to
                annotate genes were not from the same KEGG database version as the reference KO
                files. Here are the unrecognized KO IDs from the contigs database:
                {','.join(undefined_ko_ids)}\
                """
            )

        ko_dir = KODatabase.default_dir if self.ko_dir is None else self.ko_dir
        if self.modelseed_dir is None:
            modelseed_dir = ModelSEEDDatabase.default_dir
        else:
            modelseed_dir = self.modelseed_dir
        self.run.info("Reference KEGG KO database directory", ko_dir, nl_before=1)
        self.run.info("Reference ModelSEED database directory", modelseed_dir, nl_after=1)

        if store:
            if cdb_db.get_meta_value('reaction_network_ko_annotations_hash'):
                self.run.warning("Deleting existing reaction network from contigs database")
                cdb_db._exec(f'''DELETE from {tables.gene_function_reactions_table_name}''')
                cdb_db._exec(f'''DELETE from {tables.gene_function_metabolites_table_name}''')
                self.run.info_single(
                    "Deleted data in gene function reactions and metabolites tables", nl_after=1
                )

            self.progress.new("Saving reaction network to contigs database")
            self.progress.update("Reactions table")
            reactions_table = self._get_database_reactions_table(network)
            sql_statement = (
                f"INSERT INTO {tables.gene_function_reactions_table_name} VALUES "
                f"({','.join('?' * len(tables.gene_function_reactions_table_structure))})"
            )
            cdb_db._exec_many(sql_statement, reactions_table.values)
            self.progress.update("Metabolites table")
            metabolites_table = self._get_database_metabolites_table(network)
            sql_statement = (
                f"INSERT INTO {tables.gene_function_metabolites_table_name} VALUES "
                f"({','.join('?' * len(tables.gene_function_metabolites_table_structure))})"
            )
            cdb_db._exec_many(sql_statement, metabolites_table.values)

            self.progress.update("Metadata")
            ko_annotations_hash = self.hash_contigs_db_ko_hits(gene_ko_hits_table)
            cdb_db.set_meta_value('reaction_network_ko_annotations_hash', ko_annotations_hash)
            cdb_db.set_meta_value('reaction_network_kegg_database_release', ko_db.release)
            cdb_db.set_meta_value('reaction_network_modelseed_database_sha', modelseed_db.sha)
            self.progress.end()

        precomputed_counts = {
            'total_genes': cdb_db.get_row_counts_from_table('genes_in_contigs'),
            'genes_assigned_kos': gene_ko_hits_table['gene_callers_id'].nunique(),
            'kos_assigned_genes': gene_ko_hits_table['accession'].nunique()
        }
        cdb.disconnect()
        stats = network.get_overview_statistics(precomputed_counts=precomputed_counts)
        network.print_overview_statistics(stats=stats)
        if stats_file:
            network.write_overview_statistics(stats_file, stats=stats)

        return network

    def make_pangenomic_network(
        self,
        pan_db: str,
        genomes_storage_db: str,
        store: bool = True,
        overwrite_existing_network: bool = False,
        consensus_threshold: float = None,
        discard_ties: bool = False,
        stats_file: str = None
    ) -> PangenomicNetwork:
        """
        Make a pangenomic metabolic reaction network from KEGG Orthologs stored a genomes storage
        database and gene clusters stored in a pan database.

        Parameters
        ==========
        pan_db : str
            Path to a pan database. The pangenomic network is determined for gene clusters stored in
            the database.

        genomes_storage_db : str
            Path to a genomes storage database. The pangenomic network is derived from gene KO
            annotations stored in the database.

        store : bool, True
            Save the network to the pan database.

        overwrite_existing_network : bool, False
            Overwrite an existing network stored in the pan database. 'store' is also required.

        consensus_threshold : float, None
            With the default of None, the protein annotation most frequent among genes in a cluster
            is assigned to the cluster itself. If a non-default argument is provided (a value on [0,
            1]), at least this proportion of genes in the cluster must have the most frequent
            annotation for the cluster to be annotated.

        discard_ties : bool, False
            If multiple protein annotations are most frequent among genes in a cluster, then do not
            assign an annotation to the cluster itself when this argument is True. By default, this
            argument is False, so one of the most frequent annotations would be arbitrarily chosen.

        stats_file : str, None
            Write network overview statistics to a tab-delimited file at this output path.

        Returns
        =======
        PangenomicNetwork
            The network derived from the pangenomic databases.
        """
        # Preemptively check the statistics file path.
        if stats_file is not None:
            filesnpaths.is_output_file_writable(stats_file)

        # Load the pan database.
        args = Namespace()
        args.pan_db = pan_db
        args.genomes_storage = genomes_storage_db
        args.discard_ties = discard_ties
        args.consensus_threshold = consensus_threshold
        pan_super = PanSuperclass(args, r=run_quiet)

        if (
            store and
            pan_super.p_meta['reaction_network_ko_annotations_hash'] and
            not overwrite_existing_network
        ):
            raise ConfigError(
                "The existing reaction network in the pan database must be explicitly overwritten."
            )

        # Check that genome contigs databases were annotated with KOs before building the pan
        # database. Unlike in contigs super, the initialization of functions by a method of pan
        # super does not allow specification of particular functional annotation sources, with
        # concomitant checks for their existence.
        gs_info = dbinfo.GenomeStorageDBInfo(genomes_storage_db)
        gs_sources: str = gs_info.get_self_table()['gene_function_sources']
        if 'KOfam' not in [source.strip() for source in gs_sources.split(',')]:
            raise ConfigError(
                f"""\
                The genomes of the pangenome were not annotated with KOs, which can be rectified\
                by running `anvi-run-kegg-kofams` on the genome contigs databases and remaking\
                the pangenome.\
                """
            )
        pan_super.init_gene_clusters()
        pan_super.init_gene_clusters_functions()
        pan_super.init_gene_clusters_functions_summary_dict()

        self.progress.new("Building reaction network")
        self.progress.update("Loading reference databases")

        # Create the reaction network object.
        network = PangenomicNetwork(run=self.run, progress=self.progress)
        network.pan_db_source_path = os.path.abspath(pan_db)
        network.genomes_storage_db_source_path = os.path.abspath(genomes_storage_db)
        network.consensus_threshold = consensus_threshold
        network.discard_ties = discard_ties
        network.consistent_annotations = True

        # Load reference databases.
        ko_db = KODatabase(self.ko_dir)
        modelseed_db = ModelSEEDDatabase(self.modelseed_dir)
        modelseed_kegg_reactions_table = modelseed_db.kegg_reactions_table
        modelseed_ec_reactions_table = modelseed_db.ec_reactions_table
        modelseed_compounds_table = modelseed_db.compounds_table

        # Record KOs that annotated gene clusters in the pan database but for some reason aren't
        # found in the KEGG KO database set up by anvi'o.
        undefined_ko_ids: List[str] = []

        # Record ModelSEED reactions that would have been added to the reaction network if the
        # reaction had a chemical equation in the ModelSEED Biochemistry database.
        undefined_modelseed_reaction_ids: List[str] = []

        # Parse gene clusters.
        gene_clusters_functions_summary_dict: Dict = pan_super.gene_clusters_functions_summary_dict
        total_gene_clusters = len(pan_super.gene_clusters)
        num_gene_clusters_parsed = -1
        num_gene_clusters_assigned_ko = 0
        ko_ids_assigned_gene_cluster: List[str] = []
        for cluster_id, gene_cluster_functions_data in gene_clusters_functions_summary_dict.items():
            num_gene_clusters_parsed += 1
            self.progress.update(
                f"Gene clusters parsed: {num_gene_clusters_parsed} / {total_gene_clusters}"
            )

            # Retrieve the consensus KO across genes in the cluster. Parameterization of the method
            # used to select consensus KOs occurred in pan super initialization.
            gene_cluster_ko_data = gene_cluster_functions_data['KOfam']
            if gene_cluster_ko_data == {'function': None, 'accession': None}:
                # No KO was assigned to the cluster.
                continue
            ko_id = gene_cluster_ko_data['accession']
            num_gene_clusters_assigned_ko += 1
            ko_ids_assigned_gene_cluster.append(ko_id)

            # Represent the gene cluster as an object.
            gene_cluster = GeneCluster()
            gene_cluster.gene_cluster_id = cluster_id
            gene_cluster.genomes = list(pan_super.gene_clusters[cluster_id])

            try:
                # The KO and its associated reactions and metabolites have already been added to the
                # network.
                ko = network.kos[ko_id]
                is_new_ko = False
            except KeyError:
                is_new_ko = True
            if not is_new_ko:
                # Have the gene cluster reference the KO.
                gene_cluster.ko = ko
                # Add the newly encountered gene cluster to the network.
                network.gene_clusters[cluster_id] = gene_cluster
                # Proceed to the next gene cluster.
                continue

            # Get KEGG REACTION IDs and EC numbers associated with the KO.
            try:
                ko_info = ko_db.ko_table.loc[ko_id]
            except KeyError:
                # For some reason the KO annotated genes in the pangenomic databases but is not
                # found in the KEGG KO database set up by anvi'o. Do not add the alien KO or the
                # gene cluster to the network.
                undefined_ko_ids.append(ko_id)
                continue
            ko_kegg_reaction_ids = self._get_ko_kegg_reaction_ids(ko_info)
            ko_ec_numbers = self._get_ko_ec_numbers(ko_info)

            if not ko_kegg_reaction_ids and not ko_ec_numbers:
                # The KO is not associated with any KEGG reactions or EC numbers, and therefore
                # anvi'o can't relate the KO to ModelSEED reactions. Do not add the unsystematizable
                # KO or the gene cluster to the network.
                continue

            # Check if KEGG reactions and EC numbers associated with the KO have already been added
            # to the network in processing other gene clusters. To have been added to the network,
            # KEGG reactions and EC numbers must have aliased ModelSEED reactions.
            old_kegg_reaction_ids, new_kegg_reaction_ids = self._find_kegg_reaction_ids(
                ko_kegg_reaction_ids, network
            )
            old_ec_numbers, new_ec_numbers = self._find_ec_numbers(
                ko_ec_numbers, network
            )

            # Retrieve data on ModelSEED reactions aliasing KEGG reactions that haven't been added
            # to the network. Each row of the reference table represents a unique mapping of KEGG
            # reaction to ModelSEED reaction.
            modelseed_kegg_reactions_dict: Dict[int, Dict] = modelseed_kegg_reactions_table[
                modelseed_kegg_reactions_table['KEGG_REACTION_ID'].isin(new_kegg_reaction_ids)
            ].to_dict(orient='index')

            # Retrieve data on ModelSEED reactions aliasing EC numbers that haven't been added to
            # the network. Each row of the reference table represents a unique mapping of EC number
            # to ModelSEED reaction.
            modelseed_ec_reactions_dict: Dict[int, Dict] = modelseed_ec_reactions_table[
                modelseed_ec_reactions_table['EC_number'].isin(new_ec_numbers)
            ].to_dict(orient='index')

            if not (
                old_kegg_reaction_ids or
                old_ec_numbers or
                modelseed_kegg_reactions_dict or
                modelseed_ec_reactions_dict
            ):
                # None of the KEGG REACTION IDs and EC numbers associated with the KO map to
                # ModelSEED reactions (none are in the ModelSEED Biochemistry reactions table).
                continue

            # Find "undefined" ModelSEED reactions without an equation.
            undefined_modelseed_reaction_ids += self._remove_undefined_reactions(
                ko_kegg_reaction_ids,
                new_kegg_reaction_ids,
                ko_ec_numbers,
                new_ec_numbers,
                modelseed_kegg_reactions_dict,
                modelseed_ec_reactions_dict
            )

            if not (
                old_kegg_reaction_ids or
                new_kegg_reaction_ids or
                old_ec_numbers or
                new_ec_numbers
            ):
                # The KO is not added to the network as none of the reactions have equations.
                continue
            # The newly encountered KO is now known to be associated with KEGG reactions or EC
            # numbers that alias ModelSEED reactions with an equation; this new data can be added to
            # the network.

            # Add an object representing the gene cluster to the network.
            network.gene_clusters[cluster_id] = gene_cluster

            # Add an object representing the newly encountered KO to the network.
            ko = KO()
            ko.id = ko_id
            ko.name = gene_cluster_ko_data['function']
            network.kos[ko_id] = ko
            gene_cluster.ko = ko

            # Associate ModelSEED reactions that have been previously added to the network under
            # construction with the newly encountered KO.
            self._process_added_reactions(
                old_kegg_reaction_ids,
                old_ec_numbers,
                network,
                ko,
                ko_kegg_reaction_ids,
                ko_ec_numbers
            )

            # Add ModelSEED reactions aliasing newly encountered KEGG reactions and EC numbers to
            # the network.
            self._add_reactions(
                modelseed_kegg_reactions_dict,
                modelseed_ec_reactions_dict,
                network,
                modelseed_compounds_table,
                ko,
                old_kegg_reaction_ids,
                new_kegg_reaction_ids,
                old_ec_numbers,
                new_ec_numbers
            )

        if DEBUG:
            for gene_cluster in network.gene_clusters.values():
                assert gene_cluster.ko.reactions

            assert sorted(network.modelseed_kegg_aliases) == sorted(network.reactions)
            assert sorted(network.modelseed_ec_number_aliases) == sorted(network.reactions)

        undefined_modelseed_reaction_ids = set(undefined_modelseed_reaction_ids)
        undefined_ko_ids = set(undefined_ko_ids)

        self.progress.end()

        if DEBUG:
            self.run.info_single(
                f"""\
                The following ModelSEED reactions would have been added to the reaction network had
                there been a chemical equation in the ModelSEED database; perhaps it is worth
                investigating the ModelSEED reactions table to understand why this is not the case:
                {', '.join(undefined_modelseed_reaction_ids)}\
                """
            )

        if undefined_ko_ids:
            self.run.info_single(
                f"""\
                Certain gene clusters were assigned KOs that were not found in the reference KO
                database. These KOs will not be used in network construction. It could be that the
                KOfams used to annotate genes were not from the same KEGG database version as the
                reference KO files. Here are the unrecognized KO IDs from the pangenomic databases:
                {','.join(undefined_ko_ids)}\
                """
            )

        ko_dir = KODatabase.default_dir if self.ko_dir is None else self.ko_dir
        if self.modelseed_dir is None:
            modelseed_dir = ModelSEEDDatabase.default_dir
        else:
            modelseed_dir = self.modelseed_dir
        self.run.info("Reference KEGG KO database directory", ko_dir, nl_before=1)
        self.run.info("Reference ModelSEED database directory", modelseed_dir, nl_after=1)

        pdb = PanDatabase(pan_db)
        if store:
            if pan_super.p_meta['reaction_network_ko_annotations_hash']:
                self.run.warning("Deleting existing reaction network from pan database")
                pdb.db._exec(
                    f'''DELETE from {tables.pan_gene_cluster_function_reactions_table_name}'''
                )
                pdb.db._exec(
                    f'''DELETE from {tables.pan_gene_cluster_function_metabolites_table_name}'''
                )
                self.run.info_single(
                    "Deleted data in gene cluster function reactions and metabolites tables",
                    nl_after=1
                )

            self.progress.new("Saving reaction network to pan database")
            self.progress.update("Reactions table")
            reactions_table = self._get_database_reactions_table(network)
            table_name = tables.pan_gene_cluster_function_reactions_table_name
            table_structure = tables.pan_gene_cluster_function_reactions_table_structure
            pdb.db._exec_many(
                f'''INSERT INTO {table_name} VALUES ({','.join('?' * len(table_structure))})''',
                reactions_table.values
            )
            self.progress.update("Metabolites table")
            metabolites_table = self._get_database_metabolites_table(network)
            table_name = tables.pan_gene_cluster_function_metabolites_table_name
            table_structure = tables.gene_function_metabolites_table_structure
            pdb.db._exec_many(
                f'''INSERT INTO {table_name} VALUES ({','.join('?' * len(table_structure))})''',
                metabolites_table.values
            )

            self.progress.update("Metadata")
            ko_annotations_hash = self.hash_pan_db_ko_annotations(
                genomes_storage_db,
                gene_clusters_functions_summary_dict,
                consensus_threshold=consensus_threshold,
                discard_ties=discard_ties
            )
            pdb.db.set_meta_value('reaction_network_ko_annotations_hash', ko_annotations_hash)
            pdb.db.set_meta_value('reaction_network_kegg_database_release', ko_db.release)
            pdb.db.set_meta_value('reaction_network_modelseed_database_sha', modelseed_db.sha)
            pdb.db.set_meta_value('reaction_network_consensus_threshold', consensus_threshold)
            pdb.db.set_meta_value('reaction_network_discard_ties', int(discard_ties))
            self.progress.end()

        ko_ids_assigned_gene_cluster = set(ko_ids_assigned_gene_cluster)
        precomputed_counts = {
            'total_gene_clusters': pdb.meta['num_gene_clusters'],
            'gene_clusters_assigned_ko': num_gene_clusters_assigned_ko,
            'kos_assigned_gene_clusters': len(ko_ids_assigned_gene_cluster)
        }
        pdb.disconnect()
        stats = network.get_overview_statistics(precomputed_counts=precomputed_counts)
        network.print_overview_statistics(stats=stats)
        if stats_file:
            network.write_overview_statistics(stats_file, stats=stats)

        return network

    def _get_ko_kegg_reaction_ids(self, ko_info: pd.Series) -> Set[str]:
        """
        Get the set of KEGG REACTION IDs associated with the KO under consideration in network
        construction.

        Parameters
        ==========
        ko_info : pd.Series
            The row for the KO in the KO data table set up by anvi'o.

        Returns
        =======
        Set[str]
            KEGG REACTION IDs associated with the KO.
        """
        ko_kegg_reaction_info: str = ko_info.loc['reactions']
        if pd.isna(ko_kegg_reaction_info):
            # The KO is not associated with KEGG reactions.
            ko_kegg_reaction_ids: Set[str] = set()
        else:
            ko_kegg_reaction_ids = set(ko_kegg_reaction_info.split())
        return ko_kegg_reaction_ids

    def _get_ko_ec_numbers(self, ko_info: pd.Series) -> Set[str]:
        """
        Get the set of EC numbers associated with the KO under consideration in network
        construction.

        Parameters
        ==========
        ko_info : pd.Series
            The row for the KO in the KO data table set up by anvi'o.

        Returns
        =======
        Set[str]
            EC numbers associated with the KO.
        """
        ko_ec_number_info: str = ko_info.loc['ec_numbers']
        if pd.isna(ko_ec_number_info):
            # The KO is not associated with EC numbers.
            ko_ec_numbers: Set[str] = set()
        else:
            ko_ec_numbers = set(ko_ec_number_info.split())
        return ko_ec_numbers

    def _find_kegg_reaction_ids(
        self,
        ko_kegg_reaction_ids: Set[str],
        network: ReactionNetwork
    ) -> Tuple[List[str], List[str]]:
        """
        Find which KEGG reactions associated with a KO under consideration in network construction
        are already in the network.

        Parameters
        ==========
        ko_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs associated with the KO.

        Returns
        =======
        Tuple[Set[str], Set[str]]
            A set of KEGG REACTION IDs in the network and a set of those not in the network.
        """
        old_kegg_reaction_ids: List[str] = []
        new_kegg_reaction_ids: List[str] = []
        for kegg_reaction_id in ko_kegg_reaction_ids:
            if kegg_reaction_id in network.kegg_modelseed_aliases:
                old_kegg_reaction_ids.append(kegg_reaction_id)
            else:
                new_kegg_reaction_ids.append(kegg_reaction_id)
        old_kegg_reaction_ids: Set[str] = set(old_kegg_reaction_ids)
        new_kegg_reaction_ids: Set[str] = set(new_kegg_reaction_ids)
        return old_kegg_reaction_ids, new_kegg_reaction_ids

    def _find_ec_numbers(
        self,
        ko_ec_numbers: Set[str],
        network: ReactionNetwork
    ) -> Tuple[Set[str], Set[str]]:
        """
        Find which EC numbers associated with a KO under consideration in network construction are
        already in the network.

        Parameters
        ==========
        ko_ec_numbers : Set[str]

        Returns
        =======
        Tuple[List[str], List[str]]
            A set of EC numbers in the network and a set of those not in the network.
        """
        old_ec_numbers: List[str] = []
        new_ec_numbers: List[str] = []
        for ec_number in ko_ec_numbers:
            if ec_number in network.ec_number_modelseed_aliases:
                old_ec_numbers.append(ec_number)
            else:
                new_ec_numbers.append(ec_number)
        old_ec_numbers: Set[str] = set(old_ec_numbers)
        new_ec_numbers: Set[str] = set(new_ec_numbers)
        return old_ec_numbers, new_ec_numbers

    def _remove_undefined_reactions(
        self,
        ko_kegg_reaction_ids: Set[str],
        new_kegg_reaction_ids: Set[str],
        ko_ec_numbers: Set[str],
        new_ec_numbers: Set[str],
        modelseed_kegg_reactions_dict: Dict[int, Dict],
        modelseed_ec_reactions_dict: Dict[int, Dict]
    ) -> List[str]:
        """
        Find "undefined" ModelSEED reactions lacking a chemical formula, and remove these from
        further consideration in network construction.

        Parameters
        ==========
        ko_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs associated with the KO under consideration in network construction.

        new_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs associated with the KO that are not already in the network.

        ko_ec_numbers : Set[str]
            EC numbers associated with the KO under consideration in network construction.

        new_ec_numbers : Set[str]
            EC numbers associated with the KO that are not already in the network.

        modelseed_kegg_reactions_dict : Dict[int, Dict]
            Data on ModelSEED reactions aliasing the newly encountered KEGG reactions.

        modelseed_ec_reactions_dict : Dict[int, Dict]
            Data on ModelSEED reactions aliasing the newly encountered EC numbers.

        Returns
        =======
        List[str]
            Undefined ModelSEED reaction IDs.
        """
        undefined_modelseed_reaction_ids: List[str] = []

        defined_kegg_reactions: Dict[str, bool] = {}.fromkeys(new_kegg_reaction_ids, False)
        undefined_kegg_indices: List[int] = []
        for idx, modelseed_reaction_data in modelseed_kegg_reactions_dict.items():
            if pd.isna(modelseed_reaction_data['stoichiometry']):
                undefined_modelseed_reaction_ids.append(modelseed_reaction_data['id'])
                undefined_kegg_indices.append(idx)
            else:
                defined_kegg_reactions[modelseed_reaction_data['KEGG_REACTION_ID']] = True

        defined_ec_numbers: Dict[str, bool] = {}.fromkeys(new_ec_numbers, False)
        undefined_ec_indices: List[int] = []
        for idx, modelseed_reaction_data in modelseed_ec_reactions_dict.items():
            if pd.isna(modelseed_reaction_data['stoichiometry']):
                undefined_modelseed_reaction_ids.append(modelseed_reaction_data['id'])
                undefined_ec_indices.append(idx)
            else:
                defined_ec_numbers[modelseed_reaction_data['EC_number']] = True

        for kegg_reaction_id, is_defined in defined_kegg_reactions.items():
            if not is_defined:
                ko_kegg_reaction_ids.remove(kegg_reaction_id)
                new_kegg_reaction_ids.remove(kegg_reaction_id)
        for idx in undefined_kegg_indices:
            modelseed_kegg_reactions_dict.pop(idx)

        for ec_number, is_defined in defined_ec_numbers.items():
            if not is_defined:
                ko_ec_numbers.remove(ec_number)
                new_ec_numbers.remove(ec_number)
        for idx in undefined_ec_indices:
            modelseed_ec_reactions_dict.pop(idx)

        return undefined_modelseed_reaction_ids

    def _add_reactions(
        self,
        modelseed_kegg_reactions_dict: Dict[int, Dict],
        modelseed_ec_reactions_dict: Dict[int, Dict],
        network: ReactionNetwork,
        modelseed_compounds_table: pd.DataFrame,
        ko: KO,
        old_kegg_reaction_ids: Set[str],
        new_kegg_reaction_ids: Set[str],
        old_ec_numbers: Set[str],
        new_ec_numbers: Set[str]
    ) -> None:
        """
        Add ModelSEED reactions aliasing newly encountered KEGG reactions and EC numbers, which are
        referenced by the KO under consideration, to the network under construction.

        Parameters
        ==========
        modelseed_kegg_reactions_dict : Dict[int, Dict]
            Data on ModelSEED reactions aliasing the newly encountered KEGG reactions.

        modelseed_ec_reactions_dict : Dict[int, Dict]
            Data on ModelSEED reactions aliasing the newly encountered EC numbers.

        network : ReactionNetwork
            The reaction network under construction.

        modelseed_compounds_table : DataFrame
            Loaded compounds table of ModelSEED Biochemistry database set up by anvi'o.

        ko : KO
            The representation of the KO being processed.

        old_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs referenced by the KO that are in the network.

        new_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs referenced by the KO that are not in the network.

        old_ec_numbers : Set[str]
            EC numbers referenced by the KO that are in the network.

        new_ec_numbers : Set[str]
            EC numbers referenced by the KO that are not in the network.

        Returns
        =======
        None
        """
        # Add ModelSEED reactions aliasing newly encountered KEGG reactions to the network.

        # The following dictionary maps KEGG REACTION IDs referenced by the KO to aliasing ModelSEED
        # reaction IDs.
        kegg_modelseed_alias_dict: Dict[str, List[str]] = {}
        # The following dictionary maps ModelSEED reaction ID to aliasing KEGG REACTION IDs
        # referenced by the KO.
        modelseed_kegg_alias_dict: Dict[str, Tuple[ModelSEEDReaction, List[str]]] = {}
        for modelseed_reaction_data in modelseed_kegg_reactions_dict.values():
            # Each entry in the dictionary is unique to a KEGG reaction aliasing a ModelSEED
            # reaction.
            kegg_reaction_id = modelseed_reaction_data['KEGG_REACTION_ID']
            modelseed_reaction_id = modelseed_reaction_data['id']

            if kegg_reaction_id in kegg_modelseed_alias_dict:
                if DEBUG:
                    assert modelseed_reaction_id not in kegg_modelseed_alias_dict[kegg_reaction_id]
                kegg_modelseed_alias_dict[kegg_reaction_id].append(modelseed_reaction_id)
            else:
                kegg_modelseed_alias_dict[kegg_reaction_id] = [modelseed_reaction_id]

            if modelseed_reaction_id in modelseed_kegg_alias_dict:
                # The ModelSEED reaction was already added to the network, aliased by another KEGG
                # reaction referenced by the KO.
                reaction = modelseed_kegg_alias_dict[modelseed_reaction_id][0]
                is_added = True
            else:
                try:
                    # The ModelSEED reaction was already added to the network through another KO.
                    reaction = network.reactions[modelseed_reaction_id]
                    is_added = True
                except KeyError:
                    # Generate a new ModelSEED reaction object.
                    reaction = self._get_modelseed_reaction(
                        modelseed_reaction_data,
                        modelseed_compounds_table,
                        network=network
                    )
                    is_added = False
                if DEBUG:
                    # No reactions lacking an equation should make it into the network or be under
                    # consideration when this method is called during network construction.
                    assert reaction.coefficients

            if not is_added:
                # Add the new reaction to the network.
                network.reactions[modelseed_reaction_id] = reaction
                for metabolite in reaction.compounds:
                    if metabolite.modelseed_id not in network.metabolites:
                        network.metabolites[metabolite.modelseed_id] = metabolite

            if DEBUG:
                if is_added and (modelseed_reaction_id not in modelseed_kegg_alias_dict):
                    # Previously processed KO(s) must have referenced KEGG REACTION ID(s) and/or EC
                    # number(s) that aliased the ModelSEED reaction.
                    try:
                        other_kegg_reaction_ids = network.modelseed_kegg_aliases[
                            modelseed_reaction_id
                        ]
                        assert not set(other_kegg_reaction_ids).intersection(old_kegg_reaction_ids)
                        assert not set(other_kegg_reaction_ids).intersection(new_kegg_reaction_ids)
                    except KeyError:
                        pass
                    try:
                        other_ec_numbers = network.modelseed_ec_number_aliases[
                            modelseed_reaction_id
                        ]
                        assert not set(other_ec_numbers).intersection(old_ec_numbers)
                        assert not set(other_ec_numbers).intersection(new_ec_numbers)
                    except KeyError:
                        pass

            # Associate the reaction with the KO.
            ko.reactions[modelseed_reaction_id] = reaction

            try:
                modelseed_kegg_alias_tuple = modelseed_kegg_alias_dict[modelseed_reaction_id]
                modelseed_kegg_alias_tuple[1].append(kegg_reaction_id)
            except KeyError:
                modelseed_kegg_alias_dict[modelseed_reaction_id] = (reaction, [kegg_reaction_id])
                continue

        # Record KEGG reaction aliases in the network and KO.
        for kegg_reaction_id, modelseed_reaction_ids in kegg_modelseed_alias_dict.items():
            if DEBUG:
                assert kegg_reaction_id not in network.kegg_modelseed_aliases
            network.kegg_modelseed_aliases[kegg_reaction_id] = modelseed_reaction_ids
        for modelseed_reaction_id, modelseed_kegg_alias_tuple in modelseed_kegg_alias_dict.items():
            reaction = modelseed_kegg_alias_tuple[0]
            if DEBUG:
                assert not old_kegg_reaction_ids.intersection(set(reaction.kegg_aliases))
                assert not old_ec_numbers.intersection(set(reaction.ec_number_aliases))

            kegg_reaction_ids = modelseed_kegg_alias_tuple[1]
            try:
                kegg_aliases = network.modelseed_kegg_aliases[modelseed_reaction_id]
            except KeyError:
                network.modelseed_kegg_aliases[modelseed_reaction_id] = kegg_aliases = []
            if DEBUG:
                assert not set(kegg_reaction_ids).intersection(set(kegg_aliases))
            kegg_aliases += kegg_reaction_ids
            if modelseed_reaction_id not in network.modelseed_ec_number_aliases:
                # No previously processed KO(s) referenced EC number(s) that aliased the ModelSEED
                # reaction.
                network.modelseed_ec_number_aliases[modelseed_reaction_id] = []

            if DEBUG:
                assert not modelseed_reaction_id in ko.kegg_reaction_aliases
            ko.kegg_reaction_aliases[modelseed_reaction_id] = kegg_reaction_ids
            if DEBUG:
                assert not modelseed_reaction_id in ko.ec_number_aliases
            ko.ec_number_aliases[modelseed_reaction_id] = []

        # Add ModelSEED reactions aliasing newly encountered EC numbers to the network.

        # The following dictionary maps EC numbers referenced by the KO to aliasing ModelSEED
        # reaction IDs.
        ec_modelseed_alias_dict: Dict[str, List[str]] = {}
        # The following dictionary maps ModelSEED reaction ID to aliasing EC numbers referenced by
        # the KO.
        modelseed_ec_alias_dict: Dict[str, Tuple[ModelSEEDReaction, List[str]]] = {}
        for modelseed_reaction_data in modelseed_ec_reactions_dict.values():
            # Each entry in the dictionary is unique to a EC number aliasing a ModelSEED reaction.
            ec_number = modelseed_reaction_data['EC_number']
            modelseed_reaction_id = modelseed_reaction_data['id']

            if ec_number in ec_modelseed_alias_dict:
                if DEBUG:
                    assert modelseed_reaction_id not in ec_modelseed_alias_dict[ec_number]
                ec_modelseed_alias_dict[ec_number].append(modelseed_reaction_id)
            else:
                ec_modelseed_alias_dict[ec_number] = [modelseed_reaction_id]

            if modelseed_reaction_id in modelseed_kegg_alias_dict:
                # The ModelSEED reaction was aliased by a KEGG reaction referenced by the KO and so
                # was already added to the network.
                reaction = modelseed_kegg_alias_dict[modelseed_reaction_id][0]
                if modelseed_reaction_id in modelseed_ec_alias_dict:
                    ec_aliases = modelseed_ec_alias_dict[modelseed_reaction_id][1]
                    if DEBUG:
                        assert ec_number not in ec_aliases
                    ec_aliases.append(ec_number)
                else:
                    modelseed_ec_alias_dict[modelseed_reaction_id] = (reaction, [ec_number])
                continue

            if modelseed_reaction_id in modelseed_ec_alias_dict:
                # The ModelSEED reaction was already added to the network, aliased by another EC
                # number referenced by the KO.
                reaction = modelseed_ec_alias_dict[modelseed_reaction_id][0]
                is_added = True
            else:
                try:
                    # The ModelSEED reaction was already added to the network through another KO.
                    reaction = network.reactions[modelseed_reaction_id]
                    is_added = True
                except KeyError:
                    # Generate a new ModelSEED reaction object.
                    reaction = self._get_modelseed_reaction(
                        modelseed_reaction_data,
                        modelseed_compounds_table,
                        network=network
                    )
                    is_added = False
                if DEBUG:
                    # No reactions lacking an equation should make it into the network or be under
                    # consideration when this method is called during network construction.
                    assert reaction.coefficients

            if not is_added:
                # Add the new reaction to the network.
                network.reactions[modelseed_reaction_id] = reaction
                for metabolite in reaction.compounds:
                    if metabolite.modelseed_id not in network.metabolites:
                        network.metabolites[metabolite.modelseed_id] = metabolite

            if DEBUG:
                if is_added and (modelseed_reaction_id not in modelseed_ec_alias_dict):
                    # Previously processed KO(s) must have referenced KEGG REACTION ID(s) and/or EC
                    # number(s) that aliased the ModelSEED reaction.
                    try:
                        other_kegg_reaction_ids = network.modelseed_kegg_aliases[
                            modelseed_reaction_id
                        ]
                        assert not set(other_kegg_reaction_ids).intersection(old_kegg_reaction_ids)
                        assert not set(other_kegg_reaction_ids).intersection(new_kegg_reaction_ids)
                    except KeyError:
                        pass
                    try:
                        other_ec_numbers = network.modelseed_ec_number_aliases[
                            modelseed_reaction_id
                        ]
                        assert not set(other_ec_numbers).intersection(old_ec_numbers)
                        assert not set(other_ec_numbers).intersection(new_ec_numbers)
                    except KeyError:
                        pass

            # Associate the reaction with the KO.
            ko.reactions[modelseed_reaction_id] = reaction

            try:
                modelseed_ec_alias_tuple = modelseed_ec_alias_dict[modelseed_reaction_id]
                modelseed_ec_alias_tuple[1].append(ec_number)
            except KeyError:
                modelseed_ec_alias_dict[modelseed_reaction_id] = (reaction, [ec_number])

        # Record EC number aliases in the network and KO.
        for ec_number, modelseed_reaction_ids in ec_modelseed_alias_dict.items():
            if DEBUG:
                assert ec_number not in network.ec_number_modelseed_aliases
            network.ec_number_modelseed_aliases[ec_number] = modelseed_reaction_ids
        for modelseed_reaction_id, modelseed_ec_alias_tuple in modelseed_ec_alias_dict.items():
            reaction = modelseed_ec_alias_tuple[0]
            if DEBUG:
                assert not old_kegg_reaction_ids.intersection(set(reaction.kegg_aliases))
                assert not old_ec_numbers.intersection(set(reaction.ec_number_aliases))

            ec_numbers = modelseed_ec_alias_tuple[1]
            try:
                ec_aliases = network.modelseed_ec_number_aliases[modelseed_reaction_id]
            except KeyError:
                network.modelseed_ec_number_aliases[modelseed_reaction_id] = ec_aliases = []
            if DEBUG:
                assert not set(ec_numbers).intersection(set(ec_aliases))
            ec_aliases += ec_numbers
            if modelseed_reaction_id not in network.modelseed_kegg_aliases:
                # Neither this KO nor any previously processed KO(s) referenced KEGG reaction(s)
                # that aliased the ModelSEED reaction.
                network.modelseed_kegg_aliases[modelseed_reaction_id] = []

            if modelseed_reaction_id in ko.ec_number_aliases:
                # The ModelSEED reaction aliased KEGG reaction(s) refereced by the KO. An empty
                # list was added for the ModelSEED reaction in the following attribute.
                if DEBUG:
                    assert not ko.ec_number_aliases[modelseed_reaction_id]
                ko.ec_number_aliases[modelseed_reaction_id] += ec_numbers
            else:
                # The ModelSEED reaction did not alias any KEGG reactions referenced by the KO.
                ko.ec_number_aliases[modelseed_reaction_id] = ec_numbers
            if modelseed_reaction_id not in ko.kegg_reaction_aliases:
                ko.kegg_reaction_aliases[modelseed_reaction_id] = []

    def _get_modelseed_reaction(
        self,
        modelseed_reaction_data: Dict,
        modelseed_compounds_table: pd.DataFrame,
        network: ReactionNetwork = None
    ) -> ModelSEEDReaction:
        """
        Get an object representation of the ModelSEED reaction.

        The reaction object references objects representing associated ModelSEED compounds.

        Parameters
        ==========
        modelseed_reaction_data : Dict
            The dictionary representation of a row of the ModelSEED reaction table set up by anvi'o,
            containing data on the reaction.

        modelseed_compounds_table : pd.DataFrame
            Loaded ModelSEED Biochemistry compounds database.

        network : ReactionNetwork, None
            The reaction network under construction, with reaction compound objects drawn from the
            network, if possible, rather than created anew, as is done when a network is not
            provided. New reaction and compound objects are not added to the network by this method.

        Returns
        =======
        ModelSEEDReaction
            The representation of the reaction with data sourced from ModelSEED Biochemistry.
        """
        reaction = ModelSEEDReaction()

        modelseed_reaction_id = modelseed_reaction_data['id']
        if DEBUG:
            assert pd.notna(modelseed_reaction_id)
        reaction.modelseed_id = modelseed_reaction_id

        modelseed_name: str = modelseed_reaction_data['name']
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
            reaction.ec_number_aliases = tuple()
        else:
            reaction.ec_number_aliases = tuple(ec_numbers.split('|'))

        reversibility: str = modelseed_reaction_data['reversibility']
        if pd.isna(reversibility):
            reaction.reversibility = None
        elif reversibility == '=' or reversibility == '?':
            # Assume that reactions lacking data ('?') are reversible.
            reaction.reversibility = True
        else:
            reaction.reversibility = False

        stoichiometry: str = modelseed_reaction_data['stoichiometry']
        modelseed_compound_ids: List[str] = []
        if pd.isna(stoichiometry):
            if DEBUG:
                assert pd.isna(modelseed_reaction_data['reversibility'])
                assert pd.isna(modelseed_reaction_data['direction'])
            reaction.compartments = None
            reaction.coefficients = None
        else:
            if DEBUG:
                assert pd.notna(modelseed_reaction_data['reversibility'])
                assert pd.notna(modelseed_reaction_data['direction'])
            decimal_reaction_coefficients: List[float] = []
            split_stoichiometry = stoichiometry.split(';')
            compartments: List[str] = []
            for entry in split_stoichiometry:
                split_entry = entry.split(':')
                if DEBUG:
                    assert len(split_entry) > 3
                decimal_reaction_coefficients.append(float(split_entry[0]))
                modelseed_compound_ids.append(split_entry[1])
                compartments.append(ModelSEEDDatabase.compartment_ids[int(split_entry[2])])
            reaction.compartments = tuple(compartments)
            reaction_coefficients = to_lcm_denominator(decimal_reaction_coefficients)
            direction = modelseed_reaction_data['direction']
            if (
                (direction == '>' and reversibility == '<') or
                (direction == '<' and reversibility == '>')
            ):
                # The way the reaction is written is the opposite of the way the reaction proceeds.
                reaction_coefficients = [-c for c in reaction_coefficients]
            reaction.coefficients = tuple(reaction_coefficients)

        if not modelseed_compound_ids:
            reaction.compounds = None
            return reaction

        reaction_compounds: List[ModelSEEDCompound] = []
        for modelseed_compound_id in modelseed_compound_ids:
            if network:
                try:
                    # The ModelSEED compound ID has been encountered in previously processed
                    # reactions, so there is already a 'ModelSEEDCompound' object for it.
                    reaction_compounds.append(network.metabolites[modelseed_compound_id])
                    continue
                except KeyError:
                    pass

            # Generate a new metabolite object.
            try:
                modelseed_compound_series: pd.Series = modelseed_compounds_table.loc[
                    modelseed_compound_id
                ]
            except KeyError:
                raise ConfigError(
                    f"""\
                    A row for the ModelSEED compound ID, '{modelseed_compound_id}', was expected\
                    but not found in the ModelSEED compounds table. This ID was found in the\
                    equation for the ModelSEED reaction, '{modelseed_reaction_id}'.\
                    """
                )
            modelseed_compound_data = modelseed_compound_series.to_dict()
            modelseed_compound_data['id'] = modelseed_compound_id
            compound = self._get_modelseed_compound(modelseed_compound_data)
            reaction_compounds.append(compound)
        reaction.compounds = tuple(reaction_compounds)

        return reaction

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

    def _process_added_reactions(
        self,
        old_kegg_reaction_ids: Set[str],
        old_ec_numbers: Set[str],
        network: ReactionNetwork,
        ko: KO,
        ko_kegg_reaction_ids: Set[str],
        ko_ec_numbers: Set[str]
    ) -> None:
        """
        Associate ModelSEED reactions that have been previously added to the network under
        construction with the newly encountered KO.

        Parameters
        ==========
        old_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs previously added to the network.

        old_ec_numbers : Set[str]
            EC numbers previously added to the network.

        network : ReactionNetwork
            The reaction network under construction.

        ko : KO
            The representation of the KO being processed.

        ko_kegg_reaction_ids : Set[str]
            KEGG REACTION IDs associated with the KO under consideration in network construction.

        ko_ec_numbers : Set[str]
            EC numbers associated with the KO under consideration in network construction.

        Returns
        =======
        None
        """
        # Associate reactions aliasing KEGG reactions with the KO.
        for kegg_reaction_id in old_kegg_reaction_ids:
            for modelseed_reaction_id in network.kegg_modelseed_aliases[kegg_reaction_id]:
                reaction = network.reactions[modelseed_reaction_id]
                ko.reactions[modelseed_reaction_id] = reaction
                # Record which KEGG REACTION IDs and EC numbers associated with the KO alias the
                # ModelSEED reaction.
                ko.kegg_reaction_aliases[modelseed_reaction_id] = list(
                    ko_kegg_reaction_ids.intersection(set(reaction.kegg_aliases))
                )
                ko.ec_number_aliases[modelseed_reaction_id] = list(
                    ko_ec_numbers.intersection(set(reaction.ec_number_aliases))
                )

        # Associate reactions aliasing EC numbers with the KO.
        for ec_number in old_ec_numbers:
            for modelseed_reaction_id in network.ec_number_modelseed_aliases[ec_number]:
                if modelseed_reaction_id in ko.reactions:
                    # The ModelSEED reaction has already been associated with the KO, as it was
                    # aliased by a KEGG reaction referenced by the KO -- addressed above -- as well
                    # as the EC number.
                    continue
                reaction = network.reactions[modelseed_reaction_id]
                ko.reactions[modelseed_reaction_id] = reaction
                if DEBUG:
                    assert not ko_kegg_reaction_ids.intersection(set(reaction.kegg_aliases))
                ko.kegg_reaction_aliases[modelseed_reaction_id] = []
                ko.ec_number_aliases[modelseed_reaction_id] = list(
                    ko_ec_numbers.intersection(set(reaction.ec_number_aliases))
                )

    def _get_database_reactions_table(self, network: ReactionNetwork) -> pd.DataFrame:
        """
        Make a reactions table that can be stored in either a contigs or pan database, as the tables
        have the same structure. A `ReactionNetwork` can be reconstructed with the same data from
        the reactions and metabolites tables of the database.

        Parameters
        ==========
        network : ReactionNetwork
            The reaction network generated from gene or gene cluster KO annotations.

        Returns
        =======
        pd.DataFrame
            The table of reactions data to be stored in the contigs or pan database.
        """
        if DEBUG:
            assert (
                tables.gene_function_reactions_table_structure ==
                tables.pan_gene_cluster_function_reactions_table_structure
            )
            assert (
                tables.gene_function_reactions_table_types ==
                tables.pan_gene_cluster_function_reactions_table_types
            )

        # Transfer data from reaction objects to dictionaries mapping to table entries.
        reactions_data: Dict[str, Dict] = {}
        for modelseed_reaction_id, reaction in network.reactions.items():
            reaction_data = {}
            reaction_data['modelseed_reaction_id'] = modelseed_reaction_id
            reaction_data['modelseed_reaction_name'] = reaction.modelseed_name
            reaction_data['metabolite_modelseed_ids'] = ', '.join(
                [c.modelseed_id for c in reaction.compounds]
            )
            reaction_data['stoichiometry'] = ', '.join([str(c) for c in reaction.coefficients])
            reaction_data['compartments'] = ', '.join(reaction.compartments)
            reaction_data['reversibility'] = reaction.reversibility
            # Record KEGG REACTION IDs and EC numbers that are aliases of ModelSEED reactions but
            # are *NOT* associated with gene KO annotations; associated aliases are recorded later.
            reaction_data['other_kegg_reaction_ids'] = ', '.join(
                set(reaction.kegg_aliases).difference(
                    set(network.modelseed_kegg_aliases[modelseed_reaction_id])
                )
            )
            reaction_data['other_ec_numbers'] = ', '.join(
                set(reaction.ec_number_aliases).difference(
                    set(network.modelseed_ec_number_aliases[modelseed_reaction_id])
                )
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

            # Make the entry for KO KEGG REACTION aliases, which looks akin to the following
            # arbitrary example: 'R00001: (K00010, K00100); R01234: (K54321)'
            kegg_reaction_aliases = aliases[0]
            entry = []
            for kegg_reaction_id, ko_ids in kegg_reaction_aliases.items():
                entry.append(f'{kegg_reaction_id}: ({", ".join(sorted(ko_ids))})')
            reaction_data['ko_kegg_reaction_source'] = '; '.join(sorted(entry))

            # Make the entry for KO EC number aliases, which looks akin to the following arbitrary
            # example: '1.1.1.1: (K00010, K00100); 1.2.3.4: (K12345); 6.7.8.99: (K65432)
            ec_number_aliases = aliases[1]
            entry = []
            for ec_number, ko_ids in ec_number_aliases.items():
                entry.append(f'{ec_number}: ({", ".join(sorted(ko_ids))})')
            reaction_data['ko_ec_number_source'] = '; '.join(sorted(entry))

        reactions_table = pd.DataFrame.from_dict(
            reactions_data, orient='index'
        ).reset_index(drop=True).sort_values('modelseed_reaction_id')
        reactions_table = reactions_table[tables.gene_function_reactions_table_structure]
        return reactions_table

    def _get_database_metabolites_table(self, network: ReactionNetwork) -> pd.DataFrame:
        """
        Make a metabolites table that can be stored in either a contigs or pan database, as the
        tables have the same structure. A `ReactionNetwork` can be reconstructed with the same data
        from the reactions and metabolites tables of the database.

        Parameters
        ==========
        network : ReactionNetwork
            The reaction network generated from gene or gene cluster KO annotations.

        Returns
        =======
        pd.DataFrame
            The table of metabolites data to be stored in the contigs or pan database.
        """
        if DEBUG:
            assert (
                tables.gene_function_metabolites_table_structure ==
                tables.pan_gene_cluster_function_metabolites_table_structure
            )
            assert (
                tables.gene_function_metabolites_table_types ==
                tables.pan_gene_cluster_function_metabolites_table_types
            )

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

        metabolites_table = pd.DataFrame.from_dict(
            metabolites_data, orient='index'
        ).reset_index(drop=True).sort_values('modelseed_compound_id')
        metabolites_table = metabolites_table[tables.gene_function_metabolites_table_structure]
        return metabolites_table

    def hash_contigs_db_ko_hits(self, gene_ko_hits_table: pd.DataFrame) -> str:
        """
        To concisely represent the data underlying a reaction network, hash all gene KO annotations
        in the contigs database.

        Parameters
        ==========
        gene_ko_hits_table : pd.DataFrame
            This table contains gene KO hit data from the contigs database 'gene_functions' table.

        Returns
        =======
        str
            Hash representation of all gene KO annotations.
        """
        gene_ko_hits_table = gene_ko_hits_table.sort_values(['gene_callers_id', 'accession'])

        gene_ko_hits_string = ''
        for row in gene_ko_hits_table.itertuples(index=False):
            gene_ko_hits_string += str(row.gene_callers_id)
            gene_ko_hits_string += row.accession
            gene_ko_hits_string += row.function
            gene_ko_hits_string += str(row.e_value)

        hashed_gene_ko_hits = hashlib.sha1(gene_ko_hits_string.encode('utf-8')).hexdigest()
        return hashed_gene_ko_hits

    def hash_pan_db_ko_annotations(
        self,
        genomes_storage_db: str,
        gene_clusters_functions_summary_dict: Dict,
        consensus_threshold: float,
        discard_ties: bool
    ) -> str:
        """
        To concisely represent the data underlying a reaction network, hash all gene KO annotations
        in the constituent genomes, all consensus KO annotations of the gene clusters, and
        parameters used to select consensus KOs.

        Parameters
        ==========
        genomes_storage_db : str
            This is the path to a genomes storage database with the underlying gene KO annotations.

        gene_clusters_functions_summary_dict : dict
            This dictionary is loaded by a pan superclass and contains gene cluster KO annotations.

        consensus_threshold : float, None
            This parameter was used in setting consensus KO annotations of gene clusters.

        discard_ties : bool, False
            This parameter was used in setting consensus KO annotations of gene clusters.

        Returns
        =======
        str
            Hash representation of all gene cluster consensus KO annotations and the parameters used
            to select consensus KOs.
        """
        gsdb = dbinfo.GenomeStorageDBInfo(genomes_storage_db).load_db()
        functions_table = gsdb.get_table_as_dataframe('gene_function_calls', where_clause='source = "KOfam"')
        gsdb.disconnect()
        ko_annotations = []
        for row in functions_table.itertuples(index=False):
            ko_annotations.append((row.genome_name, str(row.gene_callers_id), row.accession, row.function, str(row.e_value)))
        ko_annotations = sorted(ko_annotations, key=lambda x: (x[0], x[1], x[2]))

        ko_annotations = []
        for cluster_id, gene_cluster_dict in gene_clusters_functions_summary_dict.items():
            ko_data = gene_cluster_dict['KOfam']
            ko_id = ko_data['accession']
            ko_name = ko_data['function']
            # When the KO ID and name are None, convert them into 'None'.
            ko_annotations.append((str(cluster_id), str(ko_id), str(ko_name)))
        ko_annotations = sorted(ko_annotations, key=lambda x: x[0])

        ko_annotations_string = f'{consensus_threshold}_{int(discard_ties)}_'
        for ko_annotation in ko_annotations:
            ko_annotations_string += ''.join(ko_annotation)

        hashed_ko_annotations = hashlib.sha1(ko_annotations_string.encode('utf-8')).hexdigest()
        return hashed_ko_annotations

class Tester:
    """
    This class tests reaction network construction and operations.

    Attributes
    ==========
    ko_dir : str, None
        The directory containing reference KEGG Orthology (KO) tables set up by anvi'o. This
        attribute is assigned the argument of the same name upon initialization.

    modelseed_dir : str, None
        The directory containing reference ModelSEED Biochemistry tables set up by anvi'o. This
        attribute is assigned the argument of the same name upon initialization.

    test_dir : str, None
        The directory storing test files, including copied input files and output files. With the
        default value of None, temporary directories are created and deleted as needed by methods.
        None of the test files in a provided directory, in contrast, are deleted. This attribute is
        assigned the argument of the same name upon initialization.

    run : anvio.terminal.Run, anvio.terminal.Run()
        This object prints run information to the terminal. This attribute is assigned the argument
        of the same name upon initialization.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        This object prints transient progress information to the terminal. This attribute is
        assigned the argument of the same name upon initialization.
    """
    def __init__(
        self,
        ko_dir: str = None,
        modelseed_dir: str = None,
        test_dir: str = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Parameters
        ==========
        ko_dir : str, None
            The directory containing reference KEGG Orthology (KO) tables set up by anvi'o. The
            default argument of None expects KO data to be set up in the default anvi'o directory
            used by the program `anvi-setup-kegg-data`.

        modelseed_dir : str, None
            The directory containing reference ModelSEED Biochemistry tables set up by anvi'o. The
            default argument of None expects ModelSEED data to be set up in the default anvi'o
            directory used by the program `anvi-setup-modelseed-database`.

        test_dir : str, None
            The directory storing test files. With the default value of None, temporary test
            directories are created and deleted by Tester methods; these methods operate on copies
            of input files in the test directories. In contrast, a provided directory will not be
            deleted, which can be useful for further work on output files.

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
        """
        self.ko_dir = ko_dir
        self.modelseed_dir = modelseed_dir
        self.test_dir = test_dir
        self.run = run
        self.progress = progress

    def test_contigs_database_network(self, contigs_db: str, copy_db: bool = True) -> None:
        """
        Test the construction of a reaction network from a contigs database, and test that network
        methods are able to run and do not fail certain basic (by no means comprehensive) tests.

        Parameters
        ==========
        contigs_db : str
            Path to a contigs database. The database can represent different types of samples,
            including a single genome, metagenome, or transcriptome. The network is derived from
            gene KO annotations stored in the database.

        copy_db : bool, True
            If True, as by default, store the generated reaction network in a copy of the input
            contigs database. If a test directory has been set, the database copy is placed there
            with a derived filename, e.g., "my-CONTIGS.db" is copied to a file like
            "TEST/my-CONTIGS-k2z9jxjd.db". If False, store the reaction network in the input contigs
            database, overwriting any that is already stored.

        Returns
        =======
        None
        """
        if self.test_dir is None:
            test_dir = filesnpaths.get_temp_directory_path()
        else:
            test_dir = self.test_dir
        self.run.info("Test directory", test_dir, nl_after=1)

        self.run.info_single("NETWORK CONSTRUCTION:", mc='magenta', level=0)
        utils.is_contigs_db(contigs_db)

        if copy_db:
            # Operations are performed on a copy of the contigs database in the (provided or
            # temporary) test directory.
            basename = os.path.basename(contigs_db)
            prefix, suffix = os.path.splitext(basename)
            contigs_db_target = tempfile.NamedTemporaryFile(
                prefix=f"{prefix}-", suffix=suffix, dir=test_dir
            ).name
            shutil.copy(contigs_db, contigs_db_target)
        else:
            contigs_db_target = contigs_db

        con = Constructor(
            ko_dir=self.ko_dir,
            modelseed_dir=self.modelseed_dir,
            run=self.run,
            progress=self.progress
        )

        make_stats_file_target = os.path.join(test_dir, "make_contigs_db_network_stats.tsv")
        network = con.make_contigs_database_network(
            contigs_db=contigs_db_target,
            overwrite_existing_network=True,
            stats_file=make_stats_file_target
        )

        self.run.info_single("NETWORK LOADING:", mc='magenta', level=0)
        load_stats_file_target = os.path.join(test_dir, "load_contigs_db_network_stats.tsv")
        con.load_contigs_database_network(contigs_db_target, stats_file=load_stats_file_target)

        # Check that the statistics on the network constructed and saved in the contigs database are
        # the same as the statistics on the same network loaded back into memory from the contigs
        # database.
        make_stats_table = pd.read_csv(
            make_stats_file_target,
            sep='\t',
            header=0,
            index_col='Statistic',
            usecols=['Statistic', 'Value']
        )
        make_stats_table = make_stats_table.rename({'Value': 'make'}, axis=1)
        load_stats_table = pd.read_csv(
            load_stats_file_target,
            sep='\t',
            header=0,
            index_col='Statistic',
            usecols=['Statistic', 'Value']
        )
        load_stats_table = load_stats_table.rename({'Value': 'load'}, axis=1)
        stats_table = pd.merge(
            make_stats_table, load_stats_table, left_index=True, right_index=True
        )
        inconsistent_stats: Dict[str, Tuple[float, float]] = {}
        for row in stats_table.itertuples():
            if row.make != row.load:
                inconsistent_stats[row.Index] = (row.make, row.load)
        if inconsistent_stats:
            s = ""
            for stat, stat_tuple in inconsistent_stats.items():
                s += f"{stat}: {stat_tuple[0]}, {stat_tuple[1]}; "
            s = s[:-2]
            raise AssertionError(
                f"""\
                Statistics on the network constructed and saved to the contigs database differ from\
                what should be the same statistics on the same network loaded from the contigs\
                database. Here are the different statistics, with the value from network\
                construction before the value from network loading: {s}\
                """
            )

        self.run.info_single(
            "PURGE OF METABOLITES WITHOUT FORMULA:", mc='magenta', nl_before=1, level=0
        )
        network.copy().remove_metabolites_without_formula(
            output_path=os.path.join(test_dir, "removed.tsv")
        )
        print()

        self.progress.new("Testing network purge methods")
        self.progress.update("...")
        # Network pruning tests use a random sample of half the network items (nodes) of each type.
        random.seed(RANDOM_SEED)
        metabolite_sample = set(random.sample(
            list(network.metabolites), math.ceil(len(network.metabolites) / 2)
        ))
        random.seed(RANDOM_SEED)
        reaction_sample = set(random.sample(
            list(network.reactions), math.ceil(len(network.reactions) / 2)
        ))
        random.seed(RANDOM_SEED)
        ko_sample = set(random.sample(list(network.kos), math.ceil(len(network.kos) / 2)))
        random.seed(RANDOM_SEED)
        gene_sample = set(random.sample(list(network.genes), math.ceil(len(network.genes) / 2)))

        copied_network = network.copy()
        # The basic tests of the copy method check that the network-level attributes appear to
        # contain the same items. What remains untested is that all of the references between nodes
        # are identical, e.g., the reactions referenced by each KO.
        assert list(network.metabolites) == list(copied_network.metabolites)
        assert list(network.reactions) == list(copied_network.reactions)
        assert list(network.kos) == list(copied_network.kos)
        assert list(network.genes) == list(copied_network.genes)
        assert list(network.proteins) == list(copied_network.proteins)
        assert network.kegg_modelseed_aliases == copied_network.kegg_modelseed_aliases
        assert network.modelseed_kegg_aliases == copied_network.modelseed_kegg_aliases
        assert network.ec_number_modelseed_aliases == copied_network.ec_number_modelseed_aliases
        assert network.modelseed_ec_number_aliases == copied_network.modelseed_ec_number_aliases
        removed = copied_network.purge_metabolites(metabolite_sample)
        # The most basic test of the purge (pruning) method is that the network no longer contains
        # the items that were requested to be removed. What remains untested, and would require a
        # curated test dataset, is the removal of certain other "upstream" and "downstream" nodes
        # associated with the nodes requested to be removed, e.g., KOs and genes upstream and
        # metabolites downstream of requested reactions.
        assert metabolite_sample.difference(set(copied_network.metabolites)) == metabolite_sample
        assert not metabolite_sample.difference(
            set([metabolite.modelseed_id for metabolite in removed['metabolite']])
        )

        copied_network = network.copy()
        removed = copied_network.purge_reactions(reaction_sample)
        assert reaction_sample.difference(set(copied_network.reactions)) == reaction_sample
        assert not reaction_sample.difference(
            set([reaction.modelseed_id for reaction in removed['reaction']])
        )

        copied_network = network.copy()
        removed = copied_network.purge_kos(ko_sample)
        assert ko_sample.difference(set(copied_network.kos)) == ko_sample
        assert not ko_sample.difference(set([ko.id for ko in removed['ko']]))

        copied_network = network.copy()
        removed = copied_network.purge_genes(gene_sample)
        assert gene_sample.difference(set(copied_network.genes)) == gene_sample
        assert not gene_sample.difference(set([gene.gcid for gene in removed['gene']]))
        self.progress.end()

        self.progress.new("Testing network subset methods")
        self.progress.update("...")
        subnetwork = network.subset_network(metabolites_to_subset=metabolite_sample)
        # The most basic test of the subset method is that the new network contains the requested
        # items. What remains untested, and would require a curated test dataset, is the inclusion
        # of certain other "upstream" and "downstream" nodes associated with the nodes requested to
        # be removed, e.g., KOs and genes upstream and metabolites downstream of requested
        # reactions.
        assert not metabolite_sample.difference(set(subnetwork.metabolites))

        subnetwork = network.subset_network(reactions_to_subset=reaction_sample)
        assert not reaction_sample.difference(set(subnetwork.reactions))

        subnetwork = network.subset_network(kos_to_subset=ko_sample)
        assert not ko_sample.difference(set(subnetwork.kos))

        subnetwork = network.subset_network(genes_to_subset=gene_sample)
        assert not gene_sample.difference(set(subnetwork.genes))

        # Network merging functionality is tested within the following command.
        subnetwork = network.subset_network(
            genes_to_subset=gene_sample,
            kos_to_subset=ko_sample,
            reactions_to_subset=reaction_sample,
            metabolites_to_subset=metabolite_sample
        )
        assert not metabolite_sample.difference(set(subnetwork.metabolites))
        assert not reaction_sample.difference(set(subnetwork.reactions))
        assert not ko_sample.difference(set(subnetwork.kos))
        assert not gene_sample.difference(set(subnetwork.genes))
        self.progress.end()

        if self.test_dir is None:
            shutil.rmtree(test_dir)

        self.run.info_single(
            "All tests passed for the contigs database reaction network",
            mc='magenta',
            nl_before=1,
            level=0
        )
        self.run.info_single("Network construction and storage in the contigs database")
        self.run.info_single("Purge metabolites without formula")
        self.run.info_single("Purge select metabolites")
        self.run.info_single("Purge select reactions")
        self.run.info_single("Purge select KOs")
        self.run.info_single("Purge select genes")
        self.run.info_single("Subset select metabolites")
        self.run.info_single("Subset select reactions")
        self.run.info_single("Subset select KOs")
        self.run.info_single("Subset select genes")
        self.run.info_single("Subset select metabolites, reactions, KOs, and genes", nl_after=1)

    def test_pan_database_network(
        self,
        pan_db: str,
        genomes_storage_db: str,
        copy_db: bool = True,
        consensus_threshold: float = None,
        discard_ties: bool = False
    ) -> None:
        """
        Test the construction of a reaction network from a pan database, and test that network
        methods are able to run and do not fail certain basic (by no means comprehensive) tests.

        Parameters
        ==========
        pan_db : str
            Path to a pan database. The pangenomic network is determined for gene clusters stored in
            the database.

        genomes_storage_db : str
            Path to a genomes storage database. The pangenomic network is derived from gene KO
            annotations stored in the database.

        copy_db : bool, True
            If True, as by default, store the generated reaction network in a copy of the input pan
            database. If a test directory has been set, the database copy is placed there with a
            derived filename, e.g., "my-PAN.db" is copied to a file like "TEST/my-PAN-spiba5e7.db".
            If False, store the reaction network in the input pan database, overwriting any that is
            already stored.

        consensus_threshold : float, None
            With the default of None, the protein annotation most frequent among genes in a cluster
            is assigned to the cluster itself. If a non-default argument is provided (a value on [0,
            1]), at least this proportion of genes in the cluster must have the most frequent
            annotation for the cluster to be annotated.

        discard_ties : bool, False
            If multiple protein annotations are most frequent among genes in a cluster, then do not
            assign an annotation to the cluster itself when this argument is True. By default, this
            argument is False, so one of the most frequent annotations would be arbitrarily chosen.

        Returns
        =======
        None
        """
        if self.test_dir is None:
            test_dir = filesnpaths.get_temp_directory_path()
        else:
            test_dir = self.test_dir
        self.run.info("Test directory", test_dir, nl_after=1)

        self.run.info_single("NETWORK CONSTRUCTION:", mc='magenta', level=0)
        utils.is_pan_db(pan_db)
        utils.is_genome_storage(genomes_storage_db)

        if copy_db:
            # Operations are performed on a copy of the pan database in the (provided or temporary)
            # test directory.
            basename = os.path.basename(pan_db)
            prefix, suffix = os.path.splitext(basename)
            pan_db_target = tempfile.NamedTemporaryFile(
                prefix=f"{prefix}-", suffix=suffix, dir=test_dir
            ).name
            shutil.copy(pan_db, pan_db_target)
        else:
            pan_db_target = pan_db

        con = Constructor(
            ko_dir=self.ko_dir,
            modelseed_dir=self.modelseed_dir,
            run=self.run,
            progress=self.progress
        )

        make_stats_file_target = os.path.join(test_dir, "make_pan_db_network_stats.tsv")
        network = con.make_pangenomic_network(
            pan_db=pan_db_target,
            genomes_storage_db=genomes_storage_db,
            overwrite_existing_network=True,
            consensus_threshold=consensus_threshold,
            discard_ties=discard_ties,
            stats_file=make_stats_file_target
        )

        self.run.info_single("NETWORK LOADING:", mc='magenta', level=0)
        load_stats_file_target = os.path.join(test_dir, "load_pan_db_network_stats.tsv")
        con.load_pan_database_network(
            pan_db_target, genomes_storage_db, stats_file=load_stats_file_target
        )

        # Check that the statistics on the network constructed and saved in the pan database are the
        # same as the statistics on the same network loaded back into memory from the pan database.
        # At the moment, this is used as the best substitute for an equality method that compares
        # every item in two networks.
        make_stats_table = pd.read_csv(
            make_stats_file_target,
            sep='\t',
            header=0,
            index_col='Statistic',
            usecols=['Statistic', 'Value']
        )
        make_stats_table = make_stats_table.rename({'Value': 'make'}, axis=1)
        load_stats_table = pd.read_csv(
            load_stats_file_target,
            sep='\t',
            header=0,
            index_col='Statistic',
            usecols=['Statistic', 'Value']
        )
        load_stats_table = load_stats_table.rename({'Value': 'load'}, axis=1)
        stats_table = pd.merge(
            make_stats_table, load_stats_table, left_index=True, right_index=True
        )
        inconsistent_stats: Dict[str, Tuple[float, float]] = {}
        for row in stats_table.itertuples():
            if row.make != row.load:
                inconsistent_stats[row.Index] = (row.make, row.load)
        if inconsistent_stats:
            s = ""
            for stat, stat_tuple in inconsistent_stats.items():
                s += f"{stat}: {stat_tuple[0]}, {stat_tuple[1]}; "
            s = s[:-2]
            raise AssertionError(
                f"""\
                Statistics on the network constructed and saved to the pan database differ from\
                what should be the same statistics on the same network loaded from the pan\
                database. Here are the different statistics, with the value from network\
                construction before the value from network loading: {s}\
                """
            )

        self.run.info_single(
            "PURGE OF METABOLITES WITHOUT FORMULA:", mc='magenta', nl_before=1, level=0
        )
        network.copy().remove_metabolites_without_formula(
            output_path=os.path.join(test_dir, "removed.tsv")
        )
        print()

        self.progress.new("Testing network purge methods")
        self.progress.update("...")
        # Network pruning tests use a random sample of half the network items (nodes) of each type.
        random.seed(RANDOM_SEED)
        metabolite_sample = set(random.sample(
            list(network.metabolites), math.ceil(len(network.metabolites) / 2)
        ))
        random.seed(RANDOM_SEED)
        reaction_sample = set(random.sample(
            list(network.reactions), math.ceil(len(network.reactions) / 2)
        ))
        random.seed(RANDOM_SEED)
        ko_sample = set(random.sample(list(network.kos), math.ceil(len(network.kos) / 2)))
        random.seed(RANDOM_SEED)
        gene_cluster_sample = set(random.sample(
            list(network.gene_clusters), math.ceil(len(network.gene_clusters) / 2)
        ))

        copied_network = network.copy()
        # The basic tests of the copy method check that the network-level attributes appear to
        # contain the same items. What remains untested is that all of the references between nodes
        # are identical, e.g., the reactions referenced by each KO.
        assert list(network.metabolites) == list(copied_network.metabolites)
        assert list(network.reactions) == list(copied_network.reactions)
        assert list(network.kos) == list(copied_network.kos)
        assert list(network.gene_clusters) == list(copied_network.gene_clusters)
        assert network.kegg_modelseed_aliases == copied_network.kegg_modelseed_aliases
        assert network.modelseed_kegg_aliases == copied_network.modelseed_kegg_aliases
        assert network.ec_number_modelseed_aliases == copied_network.ec_number_modelseed_aliases
        assert network.modelseed_ec_number_aliases == copied_network.modelseed_ec_number_aliases
        removed = copied_network.purge_metabolites(metabolite_sample)
        # The most basic test of the purge (pruning) method is that the network no longer contains
        # the items that were requested to be removed. What remains untested, and would require a
        # curated test dataset, is the removal of certain other "upstream" and "downstream" nodes
        # associated with the nodes requested to be removed, e.g., KOs and gene clusters upstream
        # and metabolites downstream of requested reactions.
        assert metabolite_sample.difference(set(copied_network.metabolites)) == metabolite_sample
        assert not metabolite_sample.difference(
            set([metabolite.modelseed_id for metabolite in removed['metabolite']])
        )

        copied_network = network.copy()
        removed = copied_network.purge_reactions(reaction_sample)
        assert reaction_sample.difference(set(copied_network.reactions)) == reaction_sample
        assert not reaction_sample.difference(
            set([reaction.modelseed_id for reaction in removed['reaction']])
        )

        copied_network = network.copy()
        removed = copied_network.purge_kos(ko_sample)
        assert ko_sample.difference(set(copied_network.kos)) == ko_sample
        assert not ko_sample.difference(set([ko.id for ko in removed['ko']]))

        copied_network = network.copy()
        removed = copied_network.purge_gene_clusters(gene_cluster_sample)
        assert (
            gene_cluster_sample.difference(set(copied_network.gene_clusters)) ==
            gene_cluster_sample
        )
        assert not gene_cluster_sample.difference(
            set([gene_cluster.gene_cluster_id for gene_cluster in removed['gene_cluster']])
        )
        self.progress.end()

        self.progress.new("Testing network subset methods")
        self.progress.update("...")
        subnetwork = network.subset_network(metabolites_to_subset=metabolite_sample)
        # The most basic test of the subset method is that the new network contains the requested
        # items. What remains untested, and would require a curated test dataset, is the inclusion
        # of certain other "upstream" and "downstream" nodes associated with the nodes requested to
        # be removed, e.g., KOs and gene clusters upstream and metabolites downstream of requested
        # reactions.
        assert not metabolite_sample.difference(set(subnetwork.metabolites))

        subnetwork = network.subset_network(reactions_to_subset=reaction_sample)
        assert not reaction_sample.difference(set(subnetwork.reactions))

        subnetwork = network.subset_network(kos_to_subset=ko_sample)
        assert not ko_sample.difference(set(subnetwork.kos))

        subnetwork = network.subset_network(gene_clusters_to_subset=gene_cluster_sample)
        assert not gene_cluster_sample.difference(set(subnetwork.gene_clusters))

        # Network merging functionality is tested within the following command.
        subnetwork = network.subset_network(
            gene_clusters_to_subset=gene_cluster_sample,
            kos_to_subset=ko_sample,
            reactions_to_subset=reaction_sample,
            metabolites_to_subset=metabolite_sample
        )
        assert not metabolite_sample.difference(set(subnetwork.metabolites))
        assert not reaction_sample.difference(set(subnetwork.reactions))
        assert not ko_sample.difference(set(subnetwork.kos))
        assert not gene_cluster_sample.difference(set(subnetwork.gene_clusters))
        self.progress.end()

        if self.test_dir is None:
            shutil.rmtree(test_dir)

        self.run.info_single(
            "All tests passed for the pan database reaction network",
            mc='magenta',
            nl_before=1,
            level=0
        )
        self.run.info_single("Network construction and storage in the pan database")
        self.run.info_single("Purge metabolites without formula")
        self.run.info_single("Purge select metabolites")
        self.run.info_single("Purge select reactions")
        self.run.info_single("Purge select KOs")
        self.run.info_single("Purge select gene clusters")
        self.run.info_single("Subset select metabolites")
        self.run.info_single("Subset select reactions")
        self.run.info_single("Subset select KOs")
        self.run.info_single("Subset select gene clusters")
        self.run.info_single(
            "Subset select metabolites, reactions, KOs, and gene clusters", nl_after=1
        )

def get_chemical_equation(reaction: ModelSEEDReaction) -> str:
    """
    Get a decent-looking chemical equation.

    Parameters
    ==========
    reaction : ModelSEEDReaction
        The representation of the reaction with data sourced from ModelSEED Biochemistry.

    Returns
    =======
    str
        The stoichiometric equation has integer coefficients; reactants and products are represented
        by ModelSEED Biochemistry compound names and compartment symbols "(c)" if cytosolic and
        "(e)" if extracellular; and a unidirectional arrow, "->", if irreversible and bidirectional
        arrow, "<->", if reversible.
    """
    equation = ""
    leftside = True
    for coefficient, metabolite, compartment in zip(
        reaction.coefficients, reaction.compounds, reaction.compartments
    ):
        if leftside and coefficient > 0:
            leftside = False
            if reaction.reversibility:
                equation += "<-> "
            else:
                equation += "-> "

        if leftside:
            coeff = -coefficient
        else:
            coeff = coefficient
        equation += f"{coeff} {metabolite.modelseed_name} [{compartment}] + "

    return equation.rstrip('+ ')

def to_lcm_denominator(floats: Iterable[float]) -> Tuple[int]:
    """
    Convert a list of floats to a list of integers, with a list containing fractional numbers
    transformed to a list of lowest common integer multiples.

    Parameters
    ==========
    floats : Iterable[float]
        Numbers to convert.

    Returns
    =======
    List[int]
        List of integers transformed from the input list.
    """
    def lcm(a, b):
        return a * b // math.gcd(a, b)

    rationals = [fractions.Fraction(f).limit_denominator() for f in floats]
    lcm_denom = functools.reduce(lcm, [r.denominator for r in rationals])

    return list(int(r.numerator * lcm_denom / r.denominator) for r in rationals)

def _download_worker(
    input_queue: mp.Queue,
    output_queue: mp.Queue,
    max_num_tries: int = 100,
    wait_secs: float = 10.0
) -> None:
    """
    Multiprocessing worker to download files from a queue.

    Parameters
    ==========
    input_queue : multiprocessing.Queue
        Queue of length-two iterables of the URL and local path for each file to download.

    output_queue : multiprocessing.Queue
        Queue in which the success of each download operation is recorded, with True put in the
        output queue if the download succeeded and the local path from the input queue put in the
        output queue if the download failed (after exceeding the maximum number of tries).

    max_num_tries : int, 100
        The maximum number of times to try downloading a file (in case of a connection reset).

    wait_secs : float, 10.0
        The number of seconds to wait between each file download attempt.

    Returns
    =======
    None
    """
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
