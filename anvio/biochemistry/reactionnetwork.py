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

import anvio.kegg as kegg
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
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
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

    Objects of this class are stored in the 'metabolites' attribute of a ReactionNetwork instance.

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

    Objects of this class are stored in the 'reactions' attribute of a ReactionNetwork instance.

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

    compound_ids : Tuple[str], None
        ModelSEED IDs of reactants and products involved in the reaction. For example, 'rxn00001'
        involves the ModelSEED compounds, 'cpd00001', 'cpd00012', 'cpd00009', and 'cpd00067'. A
        compound ID is formatted 'cpdXXXXX', where each X is a digit, e.g., 'cpd00001'. IDs can be
        used to look up metabolite objects in the 'metabolites' attribute of the ReactionNetwork
        containing the reaction. Each metabolite object has a corresponding stoichiometric reaction
        coefficient in the reaction attribute, 'coefficients', and a corresponding cellular
        compartment in the reaction attribute, 'compartments'.

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
    compound_ids: Tuple[str] = None
    coefficients: Tuple[int] = None
    compartments: Tuple[str] = None
    reversibility: bool = None

@dataclass
class KO:
    """
    Representation of a KEGG Ortholog (KO) in a reaction network.

    Objects of this class are stored in the 'kos' attribute of a ReactionNetwork instance.

    Attributes
    ==========
    id : str, None
        KEGG ORTHOLOGY ID in the format, 'KXXXXX', where X is a digit, e.g., 'K00001'.

    name : str, None
        Name of the KO, e.g., 'K00001' has the name, 'alcohol dehydrogenase [EC:1.1.1.1]'.

    module_ids : List[str], list()
        IDs of KEGG modules containing the KO, which can be used to look up module objects in the
        'pathways' attribute of the ReactionNetwork containing the KO.

    hierarchies : Dict[str, List[Tuple[str]]], dict()
        Membership of the KO in BRITE hierarchies. Keys are hierarchy IDs. Values are dictionary
        representations of categorizations in the hierarchy. For example, 'K00844', hexokinase, is
        classified multiple ways in the 'KEGG Orthology (KO)' hierarchy, 'ko00001', including '09100
        Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis
        [PATH:ko00010]' and '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00051 Fructose
        and mannose metabolism [PATH:ko00051]'. This hierarchy and these classifications would be
        represented as follows: {'ko00001': [('09100 Metabolism', '09101 Carbohydrate metabolism',
        '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'), ('09100 Metabolism', '09101
        Carbohydrate metabolism', '00051 Fructose and mannose metabolism [PATH:ko00051]'), ...],
        ...} Hierarchy IDs and categorization tuples can be used to look up category objects in the
        'categories' attribute of the ReactionNetwork containing the KO.

    pathway_ids : List[str], list()
        IDs of KEGG pathways containing the KO, which can be used to look up pathway objects in the
        'pathways' attribute of the ReactionNetwork containing the KO.

    reaction_ids : List[str], list()
        IDs of ModelSEED reactions associated with the KO via KEGG reaction and EC number
        annotations of the KO. A ModelSEED reaction ID is formatted 'rxnXXXXX', where each X is a
        digit, e.g., 'rxn00001'. ModelSEED reaction IDs can be used to look up reaction objects in
        the 'reactions' attribute of the ReactionNetwork containing the KO.

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
    module_ids: List[str] = field(default_factory=list)
    hierarchies: Dict[str, List[Tuple[str]]] = field(default_factory=dict)
    pathway_ids: List[str] = field(default_factory=list)
    reaction_ids: List[str] = field(default_factory=list)
    kegg_reaction_aliases: Dict[str, List[str]] = field(default_factory=dict)
    ec_number_aliases: Dict[str, List[str]] = field(default_factory=dict)

@dataclass
class KEGGModule:
    """
    Representation of a KEGG module with KOs in a reaction network.

    Objects of this class are stored in the 'modules' attribute of a ReactionNetwork instance.

    Attributes
    ==========
    id : str, None
        KEGG MODULE ID in the format, 'MXXXXX', where X is a digit, e.g., 'M00001'.

    name : str, None
        Name of the module, e.g., 'M00001' has the name, 'Glycolysis (Embden-Meyerhof pathway),
        glucose => pyruvate'.

    ko_ids : List[str], list()
        IDs of reaction network KOs that are in the module, which can be used to look up KO objects
        in the 'kos' attribute of the ReactionNetwork containing the module. To reiterate, this does
        not include KOs in the module that are not in the reaction network.

    pathway_ids : List[str], list()
        IDs of KEGG pathways containing the module, which can be used to look up KEGG pathway
        objects in the 'pathways' attribute of the ReactionNetwork containing the module.
    """
    id: str = None
    name: str = None
    ko_ids: List[str] = field(default_factory=list)
    pathway_ids: List[str] = field(default_factory=list)

@dataclass
class KEGGPathway:
    """
    Representation of a KEGG pathway with KOs in a reaction network.

    Objects of this class are stored in the 'pathways' attribute of a ReactionNetwork instance.

    Attributes
    ==========
    id : str, None
        The KEGG PATHWAY ID in the format, 'XXXXX', where X is a digit, e.g., '00010' represents
        'Glycolysis / Gluconeogensis', and corresponds to the reference pathway map, 'map00010'.

    name : str, None
        Name of the pathway, e.g., 'Glycolysis / Gluconeogenesis' for pathway ID '00010'.

    categorization : Tuple[str], None
        Certain pathways are equivalent to bottommost categories in the KEGG BRITE hierarchy,
        'ko00001', e.g., '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]', which is represented
        in the categorization tuple, ('09100 Metabolism', '09101 Carbohydrate metabolism', '00010
        Glycolysis / Gluconeogenesis [PATH:ko00010]'). The categorization tuple can be used to
        retrieve the category object from the 'categories' attribute of the ReactionNetwork
        containing the pathway, e.g., `category = network.categories['ko00001'][('09100 Metabolism',
        '09101 Carbohydrate metabolism', '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]')]`

    ko_ids : List[str], list()
        IDs of reaction network KOs that are in the pathway, which can be used to look up KO objects in the
        'kos' attribute of the ReactionNetwork containing the pathway. To reiterate, this does not
        include KOs in the pathway that are not in the reaction network.

    module_ids : List[str], list()
        IDs of modules in the pathway that contain reaction network KOs. Module IDs can be used to
        look up module objects in the 'modules' attribute of the ReactionNetwork containing the
        pathway. To reiterate, this does not include modules in the pathway that do not contain KOs
        in the reaction network.
    """
    id: str = None
    name: str = None
    categorization: Tuple[str] = None
    ko_ids: List[str] = field(default_factory=list)
    module_ids: List[str] = field(default_factory=list)

@dataclass
class BRITEHierarchy:
    """
    Representation of a KEGG BRITE hierarchy with KOs in a reaction network.

    Objects of this class are stored in the 'hierarchies' attribute of a ReactionNetwork instance.

    Attributes
    ==========
    id : str, None
        BRITE hierarchy ID in the format, 'koXXXXX' or 'brXXXXX', where X is a digit, e.g.,
        'ko00001'. Currently, given the anvi'o KEGG data setup, the ReactionNetwork will only
        contain hierarchies with IDs in the 'koXXXXX' format.

    name : str, None
        Name of the hierarchy, e.g., 'ko00001' has the name, 'KEGG Orthology (KO)'.

    categorizations : List[Tuple[str]], list()
        Categorizations of reaction network KOs in the hierarchy. To reiterate, this does not
        include categories that do not contain KOs in the reaction network. Categories at each level
        receive their own entries. For example, 'K00844', hexokinase, is classified multiple ways in
        the 'KEGG Orthology (KO)' hierarchy, 'ko00001', including '09100 Metabolism >>> 09101
        Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis [PATH:00010]' and '09100
        Metabolism >>> 09101 Carbohydrate metabolism >>> 00051 Fructose and mannose metabolism
        [PATH:00051]'. These categorizations would yield four entries like the following: [('09100
        Metabolism', ), ('09100 Metabolism', '09101 Carbohydrate metabolism'), ('09100 Metabolism',
        '09101 Carbohydrate metabolism', '00010 Glycolysis / Gluconeogenesis [PATH:00010]'), ('09100
        Metabolism', '09101 Carbohydrate metabolism', '00051 Fructose and mannose metabolism
        [PATH:00051]')]. Each categorization tuple can be used to retrieve the corresponding
        category object from the 'categories' attribute of the ReactionNetwork containing the
        hierarchy, e.g., `category = network.categories['ko00001'][('09100 Metabolism', '09101
        Carbohydrate metabolism')]`

    ko_ids : List[str], list()
        IDs of reaction network KOs in the hierarchy, which can be used to look up KO objects in the
        'kos' attribute of the ReactionNetwork containing the hierarchy. To reiterate, this does not
        include KOs in the hierarchy that are not in the reaction network.
    """
    id: str = None
    name: str = None
    categorizations: List[Tuple[str]] = field(default_factory=list)
    ko_ids: List[str] = field(default_factory=list)

@dataclass
class BRITECategory:
    """
    Representation of a KEGG BRITE hierarchy category with KOs in a reaction network.

    Objects of this class are stored in the 'categories' attribute of a ReactionNetwork instance.

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

    hierarchy_id : str, None
        ID of the BRITE hierarchy containing the category, which can be used to look up the
        hierarchy object in the 'hierarchies' attribute of the ReactionNetwork containing the
        category.

    subcategory_names : List[str], list()
        The names of encompassed categories containing KOs in the reaction network. This is an empty
        list if there are no categories lower in the hierarchy. For example, the category,
        'Polyketide synthase (PKS) >>> Modular type I PKS' in the hierarchy, 'ko01008' encompasses
        the categories, 'cis-AT PKS' and 'trans-AT PKS'. Objects representing these subcategories
        can be looked up in the 'categories' attribute of the ReactionNetwork containing the
        category, e.g., `cis_category = network.categories['ko01008'][('Polyketide synthase (PKS)',
        'Modular type I PKS', 'cis-AT PKS')]` and `trans_category = network.categories['ko01008']
        [('Polyketide synthase (PKS)', 'Modular type I PKS', 'trans-AT PKS')]`

    pathway_id : str, None
        Certain bottommost categories in the hierarchy, 'ko00001', are equivalent to KEGG pathways,
        e.g., '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'. This attribute encodes any
        equivalent pathway ID, which can be used to look up the pathway object using the 'pathways'
        attribute of the ReactionNetwork containing the category.

    ko_ids : List[str], list()
        IDs of Reaction network KOs in the category (and all subcategories), which can be used to
        look up KO objects in the 'kos' attribute of the ReactionNetwork containing the category. To
        reiterate, this does not include KOs in the category that are not in the reaction network.
    """
    id: str = None
    name: str = None
    hierarchy_id: str = None
    subcategory_names: List[str] = field(default_factory=list)
    pathway_id: str = None
    ko_ids: List[str] = field(default_factory=list)

@dataclass
class Gene:
    """
    Representation of a gene in a genomic reaction network.

    Objects of this class are stored in the 'categories' attribute of a GenomicNetwork instance.

    Attributes
    ==========
    gcid : int, None
        The gene callers ID, or unique anvi'o identifier, of the gene: a non-negative integer.

    ko_ids : List[str], list()
        IDs of KOs annotating the gene, which can be used to look up KO objects in the 'kos'
        attribute of the GenomicNetwork containing the gene.

    e_values : Dict[str, float], dict()
        E-values express the strength of KO-gene associations. Keys are KO IDs; values are
        non-negative numbers.

    protein_id : Protein, None
        ID of the protein expressed by the gene. The protein is used for storing abundance data,
        from proteomics, for instance.
    """
    gcid: int = None
    ko_ids: List[str] = field(default_factory=list)
    e_values: Dict[str, float] = field(default_factory=dict)
    protein_id: str = None

@dataclass
class Protein:
    """
    This object stores protein abundance data (from proteomics, for instance) in a reaction network.

    Objects of this class are stored in the 'proteins' attribute of a GenomicNetwork instance.

    Attributes
    ==========
    id : int, None
        The unique anvi'o ID for the protein: a non-negative integer.

    gcids : List[int], list()
        Anvi'o gene callers IDs of genes that can express the protein. These can be used to look up
        gene objects in the 'genes' attribute of the GenomicNetwork containing the gene.

    abundances : Dict[str, float], dict()
        Protein abundance profile data with each key being a sample name and each value being the
        abundance of the protein expressed by the gene in that sample.
    """
    id: int = None
    gcids: List[int] = field(default_factory=list)
    abundances: Dict[str, float] = field(default_factory=dict)

@dataclass
class GeneCluster:
    """
    Representation of a gene cluster in a pangenomic reaction network.

    Objects of this class are stored in the 'gene_clusters' attribute of a PangenomicNetwork
    instance.

    Attributes
    ==========
    gene_cluster_id : int, None
        The unique anvi'o ID for the gene cluster: a non-negative integer.

    genomes : List[str], []
        The names of the genomes contributing the genes in the cluster.

    ko_id : str, None
        ID of the consensus KO among the genes in the cluster, which can be used to look up the KO
        object in the 'kos' attribute of the PangenomicNetwork containing the gene cluster.
        (Consensus KOs can be found from a pangenome by the anvi'o method,
        'dbops.PanSuperclass.get_gene_cluster_function_summary'.) Note that the individual gene KO
        annotations underlying the consensus annotation are not tracked.
    """
    gene_cluster_id: int = None
    genomes: List[str] = field(default_factory=list)
    ko_id: str = None

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
                "The following metabolites were removed from the biomass objective, with the "
                f"original IDs aliasing the ModelSEED compound IDs in parentheses: {id_string}"
            )
        else:
            self.run.info_single(
                "The following metabolites, given by their ModelSEED compound IDs, were removed "
                f"from the biomass objective: {', '.join(missing_metabolite_ids)}"
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
        removed_metabolites: List[ModelSEEDCompound] = removed['metabolite']
        removed_metabolite_ids: List[str] = [
            metabolite.modelseed_id for metabolite in removed_metabolites
        ]
        metabolite_removed_reactions: Dict[str, List[str]] = {}
        reaction_removed_metabolites: Dict[str, List[str]] = {}
        removed_reactions: List[ModelSEEDReaction] = removed['reaction']
        for reaction in removed_reactions:
            reaction_removed_metabolites[reaction.modelseed_id] = metabolite_ids = []
            for compound_id in reaction.compound_ids:
                if compound_id in removed_metabolite_ids:
                    try:
                        metabolite_removed_reactions[compound_id].append(reaction.modelseed_id)
                    except KeyError:
                        metabolite_removed_reactions[compound_id] = [reaction.modelseed_id]
                    metabolite_ids.append(compound_id)

        metabolite_table = []
        for metabolite in removed_metabolites:
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
        for reaction in removed_reactions:
            row = []
            row.append(reaction.modelseed_id)
            row.append(reaction.modelseed_name)
            # The set accounts for the theoretical possibility that a compound is present on both
            # sides of the reaction equation and thus is recorded multiple times.
            row.append(
                ", ".join(set(reaction_removed_metabolites[reaction.modelseed_id]))
            )
            row.append(", ".join(reaction.compound_ids))
            row.append(get_chemical_equation(reaction))
            reaction_table.append(row)

        ko_table = []
        removed_kos: List[KO] = removed['ko']
        for ko in removed_kos:
            row = []
            row.append(ko.id)
            row.append(ko.name)
            row.append(", ".join(ko.reaction_ids))
            ko_table.append(row)

        module_table = []
        removed_modules: List[KEGGModule] = removed['module']
        for module in removed_modules:
            row = []
            row.append(module.id)
            row.append(module.name)
            row.append(", ".join(module.ko_ids))
            module_table.append(row)

        pathway_table = []
        removed_pathways: List[KEGGPathway] = removed['pathway']
        for pathway in removed_pathways:
            row = []
            row.append(pathway.id)
            row.append(pathway.name)
            row.append(", ".join(pathway.ko_ids))
            row.append(", ".join(pathway.module_ids))
            pathway_table.append(row)

        hierarchy_table = []
        removed_hierarchies: List[BRITEHierarchy] = removed['hierarchy']
        removed_hierarchy_names: Dict[str, str] = {}
        for hierarchy in removed_hierarchies:
            row = []
            row.append(hierarchy.id)
            row.append(hierarchy.name)
            row.append(", ".join(hierarchy.ko_ids))
            pathway_table.append(row)
            removed_hierarchy_names[hierarchy.id] = hierarchy.name

        category_table = []
        removed_categories: List[BRITECategory] = removed['category']
        for category in removed_categories:
            row = []
            row.append(category.hierarchy_id)
            try:
                row.append(self.hierarchies[category.hierarchy_id].name)
            except KeyError:
                row.append(removed_hierarchy_names[category.hierarchy_id])
            row.append(category.id[len(category.hierarchy_id) + 2:])
            row.append(", ".join(category.ko_ids))
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
        for compound_id in metabolites_to_remove:
            try:
                removed_metabolites.append(self.metabolites.pop(compound_id))
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
        removed_metabolite_ids = [metabolite.modelseed_id for metabolite in removed_metabolites]

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
        reactions_to_remove: List[str] = []
        for reaction_id, reaction in self.reactions.items():
            for compound_id in reaction.compound_ids:
                if compound_id in removed_metabolite_ids:
                    reactions_to_remove.append(reaction_id)
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
        for reaction_id in reactions_to_remove:
            try:
                removed_reactions.append(self.reactions.pop(reaction_id))
            except KeyError:
                # This occurs when the original method called is '_purge_reactions', followed by
                # '_purge_kos', followed by this method again -- 'removed_reactions' will be empty.
                # Alternatively, this occurs if the reaction in 'reactions_to_remove' is not in the
                # network.
                continue
            try:
                self.modelseed_kegg_aliases.pop(reaction_id)
            except KeyError:
                # The reaction has no KO KEGG REACTION aliases.
                pass
            try:
                self.modelseed_ec_number_aliases.pop(reaction_id)
            except KeyError:
                # The reaction has no KO EC number aliases.
                pass
        removed_reaction_ids = [reaction.modelseed_id for reaction in removed_reactions]

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
            for idx, reaction_id in enumerate(modelseed_reaction_ids):
                if reaction_id in removed_reaction_ids:
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
            for idx, reaction_id in enumerate(modelseed_reaction_ids):
                if reaction_id in removed_reaction_ids:
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
            for compound_id in reaction.compound_ids:
                metabolites_to_remove.append(compound_id)
        metabolites_to_remove = list(set(metabolites_to_remove))
        metabolites_to_spare: List[int] = []
        for idx, compound_id in enumerate(metabolites_to_remove):
            for reaction in self.reactions.values():
                for reaction_compound_id in reaction.compound_ids:
                    if compound_id == reaction_compound_id:
                        # Do not remove the metabolite, as it participates in a retained reaction.
                        metabolites_to_spare.append(idx)
                        break
                else:
                    continue
                break
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
            for reaction_id in ko.reaction_ids:
                if reaction_id in removed_reaction_ids:
                    ko_reactions_to_remove.append(reaction_id)
            if len(ko_reactions_to_remove) == len(ko.reaction_ids):
                # All reactions associated with the KO were removed, so remove the KO as well.
                kos_to_remove.append(ko_id)
                continue
            for reaction_id in ko_reactions_to_remove:
                # Only some of the reactions associated with the KO are invalid.
                ko.reaction_ids.remove(reaction_id)
                try:
                    ko.kegg_reaction_aliases.pop(reaction_id)
                except KeyError:
                    # The reaction has no KO KEGG REACTION aliases.
                    pass
                try:
                    ko.ec_number_aliases.pop(reaction_id)
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
        else:
            kos_to_remove = list(set(kos_to_remove))
        if modules_to_remove is None:
            modules_to_remove: List[str] = []
        else:
            modules_to_remove = list(set(modules_to_remove))
        if pathways_to_remove is None:
            pathways_to_remove: List[str] = []
        else:
            pathways_to_remove = list(set(pathways_to_remove))
        if hierarchies_to_remove is None:
            hierarchies_to_remove: List[str] = []
        else:
            hierarchies_to_remove = list(set(hierarchies_to_remove))
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
            kos_to_remove += module.ko_ids
        for pathway_id in pathways_to_remove:
            try:
                pathway = self.pathways[pathway_id]
            except KeyError:
                # The requested pathway is not in the network.
                continue
            kos_to_remove += pathway.ko_ids
        for hierarchy_id in hierarchies_to_remove:
            try:
                hierarchy = self.hierarchies[hierarchy_id]
            except KeyError:
                # The requested hierarchy is not in the network.
                continue
            kos_to_remove += hierarchy.ko_ids
        for hierarchy_id, categorizations in categories_to_remove.items():
            try:
                hierarchy_categorizations = self.categories[hierarchy_id]
            except KeyError:
                # The requested hierarchy is not in the network.
                continue
            for categorization in categorizations:
                try:
                    categories = hierarchy_categorizations[categorization]
                except KeyError:
                    # The requested category is not in the network.
                    continue
                category = categories[-1]
                kos_to_remove += category.ko_ids

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
        removed_ko_ids = [ko.id for ko in removed_kos]

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

        # Remove modules and pathways from the network that exclusively contain removed KOs.
        removed_modules: List[KEGGModule] = []
        removed_pathways: List[KEGGPathway] = []
        removed_hierarchies: List[BRITEHierarchy] = []
        removed_categories: List[BRITECategory] = []
        for removed_ko in removed_kos:
            removed_ko_id = removed_ko.id

            for module_id in removed_ko.module_ids:
                module = self.modules[module_id]
                module.ko_ids.remove(removed_ko_id)
                if not module.ko_ids:
                    self.modules.pop(module_id)
                    # Remove obsolete module references from pathways.
                    for pathway_id in module.pathway_ids:
                        pathway = self.pathways[pathway_id]
                        pathway.module_ids.remove(module_id)
                    removed_modules.append(module)

            for pathway_id in removed_ko.pathway_ids:
                pathway = self.pathways[pathway_id]
                pathway.ko_ids.remove(removed_ko_id)
                if not pathway.ko_ids:
                    self.pathways.pop(pathway_id)
                    # Remove obsolete pathway references from modules.
                    for module_id in pathway.module_ids:
                        module = self.modules[module_id]
                        module.pathway_ids.remove(pathway_id)
                    removed_pathways.append(pathway)

            for hierarchy_id, categorizations in removed_ko.hierarchies.items():
                hierarchy = self.hierarchies[hierarchy_id]
                hierarchy.ko_ids.remove(removed_ko_id)
                if not hierarchy.ko_ids:
                    self.hierarchies.pop(hierarchy_id)
                    network_categorizations = self.categories.pop(hierarchy_id)
                    removed_hierarchies.append(hierarchy)
                    # Record categories removed along with the hierarchy.
                    for categorization in hierarchy.categorizations:
                        removed_categories.append(network_categorizations[categorization][-1])
                    continue

                network_categorizations = self.categories[hierarchy_id]
                for categorization in categorizations:
                    categories = network_categorizations[categorization]
                    supercategory = None
                    for depth, category in enumerate(categories, 1):
                        try:
                            category.ko_ids.remove(removed_ko_id)
                        except ValueError:
                            # The KO has already been removed from the supercategory, as it was
                            # already encountered in another subcategory.
                            supercategory = category
                            continue
                        if not category.ko_ids:
                            focus_categorization = categorization[:depth]
                            # Remove obsolete category references from the hierarchy.
                            hierarchy.categorizations.remove(focus_categorization)
                            # Remove obsolete categories from the network.
                            network_categorizations.pop(focus_categorization)
                            # Remove obsolete category references from the supercategory if the
                            # supercategory also hasn't been removed.
                            if supercategory:
                                supercategory.subcategory_names.remove(category.name)
                            supercategory = None
                            removed_categories.append(category)
                            continue
                        supercategory = category

        if 'ko00001' in self.hierarchies:
            network_categorizations = self.categories['ko00001']
            for pathway in removed_pathways:
                if pathway.categorization is not None:
                    assert pathway.categorization not in network_categorizations
        for category in removed_categories:
            if category.pathway_id is not None:
                assert category.pathway_id not in self.pathways

        # Purge reactions from the network that are exclusive to removed KOs.
        reactions_to_remove: List[str] = []
        for ko in removed_kos:
            for reaction_id in ko.reaction_ids:
                reactions_to_remove.append(reaction_id)
        reactions_to_remove = list(set(reactions_to_remove))
        for ko in self.kos.values():
            reactions_to_spare: List[int] = []
            for reaction_id in ko.reaction_ids:
                for idx, reaction_id_to_remove in enumerate(reactions_to_remove):
                    if reaction_id == reaction_id_to_remove:
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
                for ko_id in gene.ko_ids:
                    if ko_id in removed_ko_ids:
                        gene_kos_to_remove.append(ko_id)
                if len(gene_kos_to_remove) == len(gene.ko_ids):
                    # All KOs matching the gene were removed, so remove it as well.
                    genes_to_remove.append(gcid)
                    continue
                for ko_id in gene_kos_to_remove:
                    gene.ko_ids.remove(ko_id)
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
                    'hierarchy',
                    'category',
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
                if cluster.ko_id in removed_ko_ids:
                    gene_clusters_to_remove.append(gene_cluster_id)
            # If this method was called from 'purge_gene_clusters' then the gene clusters that are only
            # associated with KOs removed here were already removed from the network, and
            # 'gene_clusters_to_remove' would be empty.
            if gene_clusters_to_remove:
                removed_cascading_up = self._purge_gene_clusters(gene_clusters_to_remove)
                for key in [
                    'ko',
                    'module',
                    'pathway',
                    'hierarchy',
                    'category',
                    'reaction',
                    'kegg_reaction',
                    'ec_number',
                    'metabolite'
                ]:
                    removed_cascading_up.pop(key)
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
        subnetwork: ReactionNetwork = None,
        inclusive: bool = False
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

        inclusive : bool, False
            This option applies to genomic and not pangenomic networks. If True, "inclusive"
            subsetting applies a "Midas touch" where all items in the network that are however
            associated with requested KOs are "turned to gold" and included in the subsetted
            network. In default "exclusive" subsetting, a gene added to the subsetted network due to
            references to requested KOs will be missing its references to any other unrequested KOs
            in the source network.

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

            # Copy reactions annotating the KO to the subsetted network.
            self._subset_network_by_reactions(ko.reaction_ids, subnetwork=subnetwork)

            subnetwork.kos[ko_id] = deepcopy(ko)

            # Copy modules annotating the KO to the subsetted network.
            for module_id in ko.module_ids:
                if module_id in subnetwork.modules:
                    # The module was already added to the subnetwork via another KO.
                    continue

                module = self.modules[module_id]
                subnetwork.modules[module_id] = deepcopy(module)

            # Copy pathways annotating the KO to the subsetted network.
            for pathway_id in ko.pathway_ids:
                if pathway_id in subnetwork.pathways:
                    # The pathway was already added to the subnetwork via another KO.
                    continue

                pathway = self.pathways[pathway_id]
                subnetwork.pathways[pathway_id] = deepcopy(pathway)

            # Copy hierarchies annotating the KO to the subsetted network.
            for hierarchy_id in ko.hierarchies:
                if hierarchy_id in subnetwork.hierarchies:
                    # The hierarchy and all categories were already added to the subnetwork via
                    # another KO.
                    continue

                hierarchy = self.hierarchies[hierarchy_id]
                subnetwork.hierarchies[hierarchy_id] = deepcopy(hierarchy)
                subnetwork_hierarchy_categorizations: Dict[Tuple[str], Tuple[BRITECategory]] = {}
                subnetwork.categories[hierarchy_id] = subnetwork_hierarchy_categorizations

                # Copy all categories in the hierarchy to the subsetted network.
                hierarchy_categorizations = self.categories[hierarchy_id]
                for categorization in hierarchy.categorizations:
                    if categorization in subnetwork_hierarchy_categorizations:
                        # The category must have been a supercategory of another category already
                        # copied into the subnetwork along with all of its other supercategories.
                        continue
                    categories = hierarchy_categorizations[categorization]
                    subnetwork_categories = []
                    for depth, category in enumerate(categories, 1):
                        focus_categorization = categorization[:depth]
                        try:
                            # The supercategory must have been a supercategory of another category
                            # already copied into the subnetwork.
                            subnetwork_category = subnetwork_hierarchy_categorizations[
                                focus_categorization
                            ]
                            is_category_added = True
                        except KeyError:
                            is_category_added = False
                        if not is_category_added:
                            subnetwork_category = deepcopy(category)
                        subnetwork_categories.append(subnetwork_category)
                        if is_category_added:
                            continue
                        subnetwork_hierarchy_categorizations[focus_categorization] = tuple(
                            subnetwork_categories
                        )

        if isinstance(self, GenomicNetwork):
            if subset_referencing_genes:
                # Add genes that are annotated by the subsetted KOs to the network.
                self._subset_genes_via_kos(subnetwork, inclusive=inclusive)
        elif isinstance(self, PangenomicNetwork):
            if subset_referencing_gene_clusters:
                # Add gene clusters that are annotated by the subsetted KOs to the network.
                self._subset_gene_clusters_via_kos(subnetwork)

        return subnetwork

    def _subset_network_by_reactions(
        self,
        reaction_ids: Iterable[str],
        subnetwork: ReactionNetwork = None,
        inclusive: bool = False
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

        inclusive : bool, False
            If True, "inclusive" subsetting applies a "Midas touch" where all items in the network
            that are however associated with requested reactions are "turned to gold" and included
            in the subsetted network. In default "exclusive" subsetting, KOs and genes or gene
            clusters that are added to the subsetted network due to references to requested
            reactions will be missing references to any other unrequested reactions.

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

        for reaction_id in reaction_ids:
            try:
                reaction = self.reactions[reaction_id]
            except KeyError:
                # This occurs if the requested reaction is not in the source network.
                continue
            self._subset_reaction(subnetwork, reaction)

        if subset_referencing_kos:
            # Add KOs that are annotated by the subsetted reactions to the network.
            self._subset_kos_via_reactions(subnetwork, inclusive=inclusive)

        return subnetwork

    def _subset_reaction(self, subnetwork: ReactionNetwork, reaction: ModelSEEDReaction) -> None:
        """
        Add a reaction to a subsetted network along with metabolites involved in the reaction.

        Parameters
        ==========
        subnetwork : ReactionNetwork

        reaction : ModelSEEDReaction

        Returns
        =======
        None
        """
        reaction_id = reaction.modelseed_id
        subnetwork.reactions[reaction_id] = deepcopy(reaction)

        # Copy metabolites involved in the reaction to the subnetwork.
        for compound_id in reaction.compound_ids:
            if compound_id in subnetwork.metabolites:
                continue
            metabolite = self.metabolites[compound_id]
            subnetwork.metabolites[compound_id] = deepcopy(metabolite)

        # Add KEGG reaction and EC number aliases of the reaction to the subsetted network.
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
                subnetwork.kegg_modelseed_aliases[kegg_id].append(reaction_id)
            except KeyError:
                subnetwork.kegg_modelseed_aliases[kegg_id] = [reaction_id]

        for ec_number in reaction.ec_number_aliases:
            try:
                subnetwork.ec_number_modelseed_aliases[ec_number].append(reaction_id)
            except KeyError:
                subnetwork.ec_number_modelseed_aliases[ec_number] = [reaction_id]

    def _subset_kos_via_reactions(
        self,
        subnetwork: ReactionNetwork,
        inclusive: bool = False
    ) -> None:
        """
        Add KOs that are annotated with subsetted reactions to the subsetted network.

        Then add genes that are annotated with these added KOs to the subsetted network.

        Parameters
        ==========
        subnetwork : ReactionNetwork
            The subsetted reaction network under construction.

        inclusive : bool, False
            This option applies to genomic and not pangenomic networks. If True, "inclusive"
            subsetting applies a "Midas touch" where all items in the network that are however
            associated with requested KOs are "turned to gold" and included in the subsetted
            network. In default "exclusive" subsetting, a gene added to the subsetted network due to
            references to requested KOs will be missing its references to any other unrequested KOs
            in the source network.

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

        for ko_id, ko in self.kos.items():
            # Check all KOs in the source network for subsetted reactions.
            subsetted_reaction_ids: List[str] = []
            for reaction_id in ko.reaction_ids:
                if reaction_id in subnetwork.reactions:
                    # The KO is annotated by the subsetted reaction.
                    subsetted_reaction_ids.append(reaction_id)
            if not subsetted_reaction_ids:
                # The KO is not annotated by any subsetted reactions.
                continue

            if inclusive:
                # Copy the KO, including all its references, to the subsetted network.
                subnetwork.kos[ko_id] = deepcopy(ko)

                # Add "unrequested" reactions associated with the KO to the subsetted network if not
                # already added.
                for reaction_id in ko.reaction_ids:
                    if reaction_id in subsetted_reaction_ids:
                        continue
                    self._subset_reaction(subnetwork, self.reactions[reaction_id])
                continue

            # Subsetting is exclusive, not inclusive. Add the KO only with references to subsetted
            # reactions.
            subnetwork_ko = KO(
                id=ko_id,
                name=ko.name,
                module_ids=ko.module_ids.copy(),
                hierarchies=deepcopy(ko.hierarchies),
                pathway_ids=ko.pathway_ids.copy()
            )

            for reaction_id in subsetted_reaction_ids:
                subnetwork_ko.reaction_ids.append(reaction_id)

            for kegg_id, modelseed_reaction_ids in ko.kegg_reaction_aliases.items():
                for reaction_id in modelseed_reaction_ids:
                    if reaction_id not in subsetted_reaction_ids:
                        continue
                    try:
                        subnetwork_ko.kegg_reaction_aliases[kegg_id].append(reaction_id)
                    except KeyError:
                        subnetwork_ko.kegg_reaction_aliases[kegg_id] = [reaction_id]

            for ec_number, modelseed_reaction_ids in ko.ec_number_aliases.items():
                for reaction_id in modelseed_reaction_ids:
                    if reaction_id not in subsetted_reaction_ids:
                        continue
                    try:
                        subnetwork_ko.ec_number_aliases[ec_number].append(reaction_id)
                    except KeyError:
                        subnetwork_ko.ec_number_aliases[ec_number] = [reaction_id]

            subnetwork.kos[ko_id] = subnetwork_ko

        if isinstance(self, GenomicNetwork):
            # Copy genes that are annotated with the added KOs to the subsetted network.
            self._subset_genes_via_kos(subnetwork, inclusive=inclusive)
        elif isinstance(self, PangenomicNetwork):
            # Copy gene clusters that are annotated with the added KOs to the subsetted network.
            self._subset_gene_clusters_via_kos(subnetwork)

    def _subset_network_by_metabolites(
        self,
        compound_ids: Iterable[str],
        inclusive: bool = False
    ) -> ReactionNetwork:
        """
        Subset the network by metabolites.

        Parameters
        ==========
        compound_ids : Iterable[str]
            ModelSEED compound IDs to subset.

        inclusive : bool, False
            If True, "inclusive" subsetting applies a "Midas touch" where all items in the network
            that are however associated with requested metabolites are "turned to gold" and included
            in the subsetted network. In default "exclusive" subsetting, KOs and genes or gene
            clusters that are added to the subsetted network due to references to reactions
            involving requested metabolites will be missing references to any other reactions not
            involving any requested metabolites.

        Returns
        =======
        ReactionNetwork
            New subsetted reaction network.
        """
        if isinstance(self, GenomicNetwork):
            subnetwork = GenomicNetwork()
        elif isinstance(self, PangenomicNetwork):
            subnetwork = PangenomicNetwork()
        else:
            raise AssertionError

        for reaction in self.reactions.values():
            # Check all reactions in the source network for subsetted metabolites.
            for compound_id in reaction.compound_ids:
                if compound_id in compound_ids:
                    break
            else:
                # The reaction does not involve any of the requested metabolites.
                continue
            self._subset_reaction(subnetwork, reaction)

        # Add KOs that are annotated with the added reactions to the subsetted network, and then add
        # genes or gene clusters annotated with the added KOs to the subsetted network.
        self._subset_kos_via_reactions(subnetwork, inclusive=inclusive)

        return subnetwork

    def _merge_network(self, network: ReactionNetwork) -> ReactionNetwork:
        """
        This method is used in the process of merging the network with another network to produce a
        merged network, and contains steps common to different types of network: merge the
        attributes of the networks BESIDES genes (in a GenomicNetwork) / gene clusters (in a
        PangenomicNetwork) and protein abundances (which can only be stored in a GenomicNetwork).

        Parameters
        ==========
        network : ReactionNetwork
            The other reaction network being merged.

        Returns
        =======
        merged_network : ReactionNetwork
            A merged reaction network to be completed in the calling method.
        """
        if isinstance(network, GenomicNetwork):
            merged_network = GenomicNetwork()
        elif isinstance(network, PangenomicNetwork):
            merged_network = PangenomicNetwork()
        else:
            raise AssertionError

        merged_network.metabolites = deepcopy(self.metabolites)
        merged_network.reactions = deepcopy(self.reactions)
        merged_network.kos = deepcopy(self.kos)
        merged_network.modules = deepcopy(self.modules)
        merged_network.pathways = deepcopy(self.pathways)
        merged_network.hierarchies = deepcopy(self.hierarchies)
        merged_network.categories = deepcopy(self.categories)
        merged_network.kegg_modelseed_aliases = deepcopy(self.kegg_modelseed_aliases)
        merged_network.ec_number_modelseed_aliases = deepcopy(self.ec_number_modelseed_aliases)
        merged_network.modelseed_kegg_aliases = deepcopy(self.modelseed_kegg_aliases)
        merged_network.modelseed_ec_number_aliases = deepcopy(self.modelseed_ec_number_aliases)

        # Copy unique metabolites from the second network. Assume objects representing the same
        # metabolite in both networks have identical attributes.
        for compound_id, metabolite in network.metabolites.items():
            if compound_id in merged_network.metabolites:
                continue
            merged_network.metabolites[compound_id] = deepcopy(metabolite)

        # Copy unique reactions from the second network. Assume objects representing the same
        # reaction in both networks have identical attributes.
        for reaction_id, reaction in network.reactions.items():
            if reaction_id in merged_network.reactions:
                continue
            merged_network.reactions[reaction_id] = deepcopy(reaction)

        # Reconcile reaction ID aliases, which can differ between the networks depending on the KO
        # sources of the reactions.
        for kegg_reaction_id, modelseed_reaction_ids in network.kegg_modelseed_aliases.items():
            try:
                merged_modelseed_reaction_ids = merged_network.kegg_modelseed_aliases[
                    kegg_reaction_id
                ]
            except KeyError:
                merged_network.kegg_modelseed_aliases[
                    kegg_reaction_id
                ] = modelseed_reaction_ids.copy()
                continue
            merged_network.kegg_modelseed_aliases[kegg_reaction_id] = sorted(
                set(modelseed_reaction_ids + merged_modelseed_reaction_ids)
            )

        for ec_number, modelseed_reaction_ids in network.ec_number_modelseed_aliases.items():
            try:
                merged_modelseed_reaction_ids = merged_network.ec_number_modelseed_aliases[
                    ec_number
                ]
            except KeyError:
                merged_network.ec_number_modelseed_aliases[
                    ec_number
                ] = modelseed_reaction_ids.copy()
                continue
            merged_network.ec_number_modelseed_aliases[ec_number] = sorted(
                set(modelseed_reaction_ids + merged_modelseed_reaction_ids)
            )

        for modelseed_reaction_id, kegg_reaction_ids in network.modelseed_kegg_aliases.items():
            try:
                merged_kegg_reaction_ids = merged_network.modelseed_kegg_aliases[
                    modelseed_reaction_id
                ]
            except KeyError:
                merged_network.modelseed_kegg_aliases[
                    modelseed_reaction_id
                ] = kegg_reaction_ids.copy()
                continue
            merged_network.modelseed_kegg_aliases[modelseed_reaction_id] = sorted(
                set(kegg_reaction_ids + merged_kegg_reaction_ids)
            )

        for modelseed_reaction_id, ec_numbers in network.modelseed_ec_number_aliases.items():
            try:
                merged_ec_numbers = merged_network.modelseed_ec_number_aliases[
                    modelseed_reaction_id
                ]
            except KeyError:
                merged_network.modelseed_ec_number_aliases[
                    modelseed_reaction_id
                ] = ec_numbers.copy()
                continue
            merged_network.modelseed_ec_number_aliases[modelseed_reaction_id] = sorted(
                set(ec_numbers + merged_ec_numbers)
            )

        # Copy KOs from the second network. These can have different reaction annotations, so take
        # the union of the reactions associated with the same KO. Assume KOs with the same ID have
        # the same name.
        for ko_id, ko in network.kos.items():
            try:
                merged_ko = merged_network.kos[ko_id]
            except KeyError:
                merged_network.kos[ko_id] = deepcopy(ko)
                continue

            # The KOs should be classified in the same modules, pathways, and hierarchy categories,
            # unless the networks are derived from different reference database versions.
            merged_ko.module_ids = sorted(set(ko.module_ids + merged_ko.module_ids))
            merged_ko.pathway_ids = sorted(set(ko.pathway_ids + merged_ko.pathway_ids))
            for hierarchy_id, categorizations in ko.hierarchies.items():
                try:
                    merged_hierarchy_categorizations = merged_ko.hierarchies[hierarchy_id]
                except KeyError:
                    merged_ko.hierarchies[hierarchy_id] = categorizations.copy()
                    continue
                merged_ko.hierarchies[hierarchy_id] = sorted(
                    set(categorizations + merged_hierarchy_categorizations)
                )

            merged_ko.reaction_ids = sorted(set(ko.reaction_ids + merged_ko.reaction_ids))

            for kegg_reaction_id, modelseed_reaction_ids in merged_ko.kegg_reaction_aliases.items():
                try:
                    merged_modelseed_reaction_ids = merged_ko.kegg_reaction_aliases[
                        kegg_reaction_id
                    ]
                except KeyError:
                    merged_ko.kegg_reaction_aliases = modelseed_reaction_ids.copy()
                    continue
                merged_ko.kegg_reaction_aliases[kegg_reaction_id] = sorted(
                    set(merged_modelseed_reaction_ids + modelseed_reaction_ids)
                )

            for ec_number, modelseed_reaction_ids in merged_ko.ec_number_aliases.items():
                try:
                    merged_modelseed_reaction_ids = merged_ko.ec_number_aliases[ec_number]
                except KeyError:
                    merged_ko.ec_number_aliases = modelseed_reaction_ids.copy()
                    continue
                merged_ko.ec_number_aliases[ec_number] = sorted(
                    set(merged_modelseed_reaction_ids + modelseed_reaction_ids)
                )

        # Copy modules from the second network. Modules from the two networks can contain different
        # KOs.
        module_pathway_ids: List[str] = []
        for module_id, module in network.modules.items():
            try:
                merged_module = merged_network.modules[module_id]
            except KeyError:
                merged_network.modules[module_id] = deepcopy(module)
                continue

            merged_module.ko_ids = sorted(set(module.ko_ids + merged_module.ko_ids))
            merged_module.pathway_ids = sorted(set(module.pathway_ids + merged_module.pathway_ids))
            module_pathway_ids += merged_module.pathway_ids

        # Copy pathways from the second network. Pathways from the two networks can contain
        # different KOs.
        pathway_module_ids: List[str] = []
        for pathway_id, pathway in network.pathways.items():
            try:
                merged_pathway = merged_network.pathways[pathway_id]
            except KeyError:
                merged_network.pathways[pathway_id] = deepcopy(pathway)
                continue

            merged_pathway.ko_ids = sorted(set(pathway.ko_ids + merged_pathway.ko_ids))
            merged_pathway.module_ids = sorted(set(pathway.module_ids + merged_pathway.module_ids))
            pathway_module_ids += merged_pathway.module_ids
            assert pathway.categorization == merged_pathway.categorization

        for pathway_id in module_pathway_ids:
            assert pathway_id in merged_network.pathways
        for module_id in pathway_module_ids:
            assert module_id in merged_network.modules

        # Copy hierarchies from the second network. Hierarchies from the two networks can contain
        # different categories due to different KOs. Assume hierarchies with the same ID have the
        # same name.
        for hierarchy_id, hierarchy in network.hierarchies.items():
            try:
                merged_hierarchy = merged_network.hierarchies[hierarchy_id]
            except KeyError:
                merged_network.hierarchies[hierarchy_id] = deepcopy(hierarchy)
                continue

            merged_hierarchy.categorizations = sorted(
                set(hierarchy.categorizations + merged_hierarchy.categorizations)
            )

            merged_hierarchy.ko_ids = sorted(set(hierarchy.ko_ids + merged_hierarchy.ko_ids))

        # Copy hierarchy categories from the second network.
        for hierarchy_id, categorizations in network.categories.items():
            try:
                merged_hierarchy_categorizations = merged_network.categories[hierarchy_id]
            except KeyError:
                merged_network.categories[hierarchy_id] = deepcopy(categorizations)
                continue

            for categorization, categories in categorizations.items():
                try:
                    merged_categories = merged_hierarchy_categorizations[categorization]
                    is_category_copied = True
                except KeyError:
                    is_category_copied = False

                if is_category_copied:
                    # Reconcile the subcategories and KOs contained in the category from the two
                    # networks.
                    category = categories[-1]
                    merged_category = merged_categories[-1]

                    merged_category.subcategory_names = sorted(
                        set(category.subcategory_names + merged_category.subcategory_names)
                    )

                    merged_category.ko_ids = sorted(set(category.ko_ids + merged_category.ko_ids))

                    continue

                # Copy the category and supercategories that have not already been copied.
                copied_categories = []
                for depth, category in enumerate(categories, 1):
                    focus_categorization = categorization[:depth]
                    try:
                        # The supercategory has already been copied.
                        copied_category = merged_hierarchy_categorizations[focus_categorization]
                        is_focus_category_copied = True
                    except KeyError:
                        is_focus_category_copied = False
                    if not is_focus_category_copied:
                        copied_category = deepcopy(category)
                    copied_categories.append(copied_category)
                    if is_focus_category_copied:
                        continue
                    merged_hierarchy_categorizations[focus_categorization] = tuple(
                        copied_categories
                    )

        return merged_network

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
        self.progress.new("Counting KO biological classifications")
        self.progress.update("...")
        stats['KO biological classification'] = stats_group = {}

        ko_in_module_count = 0
        ko_in_pathway_count = 0
        for ko in self.kos.values():
            if ko.module_ids:
                ko_in_module_count += 1
            if ko.pathway_ids:
                ko_in_pathway_count += 1

        module_ko_counts = []
        for module in self.modules.values():
            module_ko_counts.append(len(module.ko_ids))

        pathway_ko_counts = []
        for pathway in self.pathways.values():
            pathway_ko_counts.append(len(pathway.ko_ids))

        all_level_category_count = 0
        low_level_category_count = 0
        for categorizations in self.categories.values():
            for categories in categorizations.values():
                all_level_category_count += 1
                category = categories[-1]
                if not category.subcategory_names:
                    low_level_category_count += 1

        stats_group['KEGG modules in network'] = len(self.modules)
        stats_group['KOs in modules'] = ko_in_module_count
        stats_group['Mean KOs per module'] = round(np.mean(module_ko_counts), 1)
        stats_group['Max KOs per module'] = max(module_ko_counts)
        stats_group['KEGG pathways in network'] = len(self.pathways)
        stats_group['KOs in pathways'] = ko_in_pathway_count
        stats_group['Mean KOs per pathway'] = round(np.mean(pathway_ko_counts), 1)
        stats_group['Max KOs per pathway'] = max(pathway_ko_counts)
        stats_group['KEGG BRITE hierarchies in network'] = len(self.hierarchies)
        stats_group['All-level hierarchy categories in network'] = all_level_category_count
        stats_group['Low-level hierarchy categories in network'] = low_level_category_count

        self.progress.end()

        self.progress.new("Counting reactions and KO sources")
        self.progress.update("...")
        stats['Reactions and KO sources'] = stats_group = {}

        stats_group['Reactions in network'] = len(self.reactions)
        reaction_counts = []
        for ko in self.kos.values():
            reaction_counts.append(len(ko.reaction_ids))
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
            for compartment, coefficient, compound_id in zip(
                reaction.compartments, reaction.coefficients, reaction.compound_ids
            ):
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
        self.run.info_single("KO biological classification")
        stats_group = stats['KO biological classification']
        for key in (
            'KEGG modules in network',
            'KOs in modules',
            'Mean KOs per module',
            'Max KOs per module',
            'KEGG pathways in network',
            'KOs in pathways',
            'Mean KOs per pathway',
            'Max KOs per pathway',
            'KEGG BRITE hierarchies in network',
            'All-level hierarchy categories in network',
            'Low-level hierarchy categories in network'
        ):
            self.run.info(key, stats_group[key])

        self.run.info_single("ModelSEED reactions in network and KO sources", nl_before=1)
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

    contigs_db_source_path : str, None
        Path to the contigs database from which the network was built.

    profile_db_source_path : str, None
        Path to the profile database from which protein and metabolite abundance data was loaded.

    genes : Dict[int, Gene], dict()
        This maps gene callers IDs to object representations of genes in the network.

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
        self.proteins: Dict[int, Protein] = {}

    def remove_metabolites_without_formula(self, output_path: str = None) -> None:
        """
        Remove metabolites without a formula in the ModelSEED database from the network.

        Other items can be removed from the network by association: reactions that involve a
        formulaless metabolite; other metabolites with formulas that are exclusive to such
        reactions; KOs predicted to exclusively catalyze such reactions; and genes exclusively
        annotated with such KOs. Removed metabolites with a formula are reported alongside
        formulaless metabolites to the optional output table of removed metabolites.

        output_path : str, None
            If not None, write tab-delimited files of metabolites, reactions, KOs, KEGG modules,
            KEGG pathways, KEGG BRITE hierarchies, KEGG BRITE hierarchy categories, and genes
            removed from the network to file locations based on the provided path. For example, if
            the argument, 'removed.tsv', is provided, then the following files will be written:
            'removed-metabolites.tsv', 'removed-reactions.tsv', 'removed-kos.tsv',
            'removed-modules.tsv', 'removed-pathways.tsv', 'removed-hierarchies.tsv',
            'removed-categories.tsv', and 'removed-genes.tsv'.
        """
        if self.verbose:
            self.progress.new("Removing metabolites without a formula in the network")
            self.progress.update("...")

        if output_path:
            path_basename, path_extension = os.path.splitext(output_path)
            metabolite_path = f"{path_basename}-metabolites{path_extension}"
            reaction_path = f"{path_basename}-reactions{path_extension}"
            ko_path = f"{path_basename}-kos{path_extension}"
            module_path = f"{path_basename}-modules{path_extension}"
            pathway_path = f"{path_basename}-pathways{path_extension}"
            hierarchy_path = f"{path_basename}-hierarchies{path_extension}"
            category_path = f"{path_basename}-categories{path_extension}"
            gene_path = f"{path_basename}-genes{path_extension}"
            for path in (
                metabolite_path,
                reaction_path,
                ko_path,
                module_path,
                pathway_path,
                hierarchy_path,
                category_path,
                gene_path
            ):
                filesnpaths.is_output_file_writable(path)

        metabolites_to_remove: List[str] = []
        for compound_id, metabolite in self.metabolites.items():
            # ModelSEED compounds without a formula have a formula value of None in the network
            # object.
            if metabolite.formula is None:
                metabolites_to_remove.append(compound_id)
        removed = self._purge_metabolites(metabolites_to_remove)

        if self.verbose:
            self.progress.end()
            self.run.info("Removed metabolites", len(removed['metabolite']))
            self.run.info("Removed reactions", len(removed['reaction']))
            self.run.info("Removed KOs", len(removed['ko']))
            self.run.info("Removed KEGG modules", len(removed['module']))
            self.run.info("Removed KEGG pathways", len(removed['pathway']))
            self.run.info("Removed KEGG BRITE hierarchies", len(removed['hierarchy']))
            self.run.info("Removed KEGG BRITE hierarchy categories", len(removed['category']))
            self.run.info("Removed genes", len(removed['gene']))

        if not output_path:
            return

        if self.verbose:
            self.progress.new("Writing output files of removed network items")
            self.progress.update("...")

        gene_table = []
        for gene in removed['gene']:
            gene: Gene
            row = []
            row.append(gene.gcid)
            row.append(", ".join(gene.ko_ids))
            gene_table.append(row)

        self._write_remove_metabolites_without_formula_output(output_path, removed)

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
            self.run.info("Table of removed KEGG modules", module_path)
            self.run.info("Table of removed KEGG pathways", pathway_path)
            self.run.info("Table of removed KEGG BRITE hierarchies", hierarchy_path)
            self.run.info("Table of removed KEGG BRITE hierarchy categories", category_path)
            self.run.info("Table of removed genes", gene_path)

    def prune(
        self,
        genes_to_remove: Union[str, Iterable[int]] = None,
        proteins_to_remove: Union[str, Iterable[str]] = None,
        kos_to_remove: Union[str, Iterable[str]] = None,
        modules_to_remove: Union[str, Iterable[str]] = None,
        pathways_to_remove: Union[str, Iterable[str]] = None,
        hierarchies_to_remove: Union[str, Iterable[str]] = None,
        categories_to_remove: Dict[str, List[Tuple[str]]] = None,
        reactions_to_remove: Union[str, Iterable[str]] = None,
        metabolites_to_remove: Union[str, Iterable[str]] = None
    ) -> Dict[str, List]:
        """
        Prune items from the metabolic network.

        Pruning modifies the network in situ: use the network 'copy' method as needed to create a
        backup of the network.

        If requested genes, proteins, KOs, KEGG modules, KEGG pathways, KEGG BRITE hierarchies, KEGG
        BRITE hierarchy categories, reactions, or metabolites are not present in the network, no
        error is raised.

        Network items (e.g., genes, KOs, reactions, and metabolites) that are exclusively associated
        with requested items are also removed from the network. Example: Consider a KO that is
        requested to be removed from the network. The KO is associated with two reactions. The first
        reaction is exclusive to the KO and thus is also removed, whereas the second reaction is
        also associated with another retained KO and thus is retained in the network. The first
        reaction involves four metabolites, and two are exclusive to the reaction: these are also
        removed from the network. The KO is the only KO annotating a certain gene but one of two KOs
        annotating another gene: the former but not the latter gene is removed from the network
        along with the KO. Note that KO annotations of genes and reaction annotations of KOs can be
        selected to the exclusion of others. In the example, the latter gene is left with one KO.

        Parameters
        ==========
        genes_to_remove : Union[int, Iterable[int]], None
            Gene callers ID(s) to remove.

        proteins_to_remove : Union[str, Iterable[str]], None
            Protein ID(s) to remove if the network has been annotated with protein abundances,
            equivalent to giving the genes encoding the protein(s) to the argument,
            'genes_to_remove'.

        kos_to_remove : Union[str, Iterable[str]], None
            KO ID(s) to remove.

        modules_to_remove : Union[str, Iterable[str]], None
            KEGG module ID(s) to remove, with the effect of giving the KOs in the module(s) to the
            argument, 'kos_to_remove'. This does not remove other module annotations of these KOs
            from the network that also annotate other KOs.

        pathways_to_remove : Union[str, Iterable[str]], None
            KEGG pathway ID(s) to remove, with the effect of giving the KOs in the pathway(s) to the
            argument, 'kos_to_remove'. This does not remove other pathway annotations of these KOs
            from the network that also annotate other KOs.

        hierarchies_to_remove : Union[str, Iterable[str]], None
            KEGG BRITE hierarchy (or hierarchies) to remove, with the effect of giving the KOs in
            the hierarchy to the argument, 'kos_to_remove'. This does not remove other hierarchy
            annotations of these KOs from the network that also annotate other KOs.

        categories_to_remove : Dict[str, List[Tuple[str]]], None
            KEGG BRITE hierarchy categories to remove, with the effect of giving the KOs in the
            categories to the argument, 'kos_to_remove'. This does not remove other category
            annotations of these KOs from the network that also annotate other KOs. The dictionary
            argument is keyed by BRITE hierarchy ID and has values that list category tuples. For
            example, to remove KOs from the network contained in the 'ko00001' 'KEGG Orthology (KO)'
            hierarchy categories, '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00010
            Glycolysis / Gluconeogenesis [PATH:ko00010]' and '09100 Metabolism >>> 09101
            Carbohydrate metabolism >>> 00051 Fructose and mannose metabolism [PATH:ko00051]', the
            dictionary argument would need to look like the following: {'ko00001': [('09100
            Metabolism', '09101 Carbohydrate metabolism', '00010 Glycolysis / Gluconeogenesis'),
            ('09100 Metabolism', '09101 Carbohydrate metabolism', '00051 Fructose and mannose
            metabolism [PATH:ko00051]')]}

        reactions_to_remove : Union[str, Iterable[str]], None
            ModelSEED reaction ID(s) to remove.

        metabolites_to_remove : Union[str, Iterable[str]], None
            ModelSEED compound ID(s) to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            The dictionary has the following format. It shows protein entries as if the network has
            been annotated with protein abundances; these are absent for genomic networks lacking
            protein annotations.
            {
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>]
            }
        """
        assert (
            genes_to_remove or
            proteins_to_remove or
            kos_to_remove or
            modules_to_remove or
            pathways_to_remove or
            hierarchies_to_remove or
            categories_to_remove or
            reactions_to_remove or
            metabolites_to_remove
        )

        removed: Dict[str, List] = {
            'gene': [],
            'protein': [],
            'ko': [],
            'module': [],
            'pathway': [],
            'hierarchy': [],
            'category': [],
            'reaction': [],
            'kegg_reaction': [],
            'ec_number': [],
            'metabolite': []
        }
        if not self.proteins:
            removed.pop('protein')

        if genes_to_remove or proteins_to_remove:
            for item_type, removed_items in self._purge_genes(
                genes_to_remove=genes_to_remove,
                proteins_to_remove=proteins_to_remove
            ).items():
                removed[item_type] += removed_items

        if (
            kos_to_remove or
            modules_to_remove or
            pathways_to_remove or
            hierarchies_to_remove or
            categories_to_remove
        ):
            for item_type, removed_items in self._purge_kos(
                kos_to_remove=kos_to_remove,
                modules_to_remove=modules_to_remove,
                pathways_to_remove=pathways_to_remove,
                hierarchies_to_remove=hierarchies_to_remove,
                categories_to_remove=categories_to_remove
            ).items():
                removed[item_type] += removed_items

        if reactions_to_remove:
            for item_type, removed_items in self._purge_reactions(reactions_to_remove).items():
                removed[item_type] += removed_items

        if metabolites_to_remove:
            for item_type, removed_items in self._purge_metabolites(metabolites_to_remove).items():
                removed[item_type] += removed_items

        return removed

    def _purge_genes(
        self,
        genes_to_remove: Iterable[int] = None,
        proteins_to_remove: Iterable[str] = None
    ) -> Dict[str, List]:
        """
        Remove any trace of the given genes from the network.

        If proteins are provided, remove any trace of the genes that encode these from the network.

        KOs, reactions, and metabolites that are only associated with removed genes are purged. KEGG
        modules, pathways, BRITE hierarchies, and BRITE hierarchy categories only associated with
        purged KOs are removed.

        Parameters
        ==========
        genes_to_remove : Iterable[int], None
            Gene callers IDs to remove.

        proteins_to_remove : Iterable[str], None
            Protein IDs to remove if the network has been annotated with protein abundances,
            equivalent to giving the genes encoding the proteins to the argument, 'genes_to_remove'.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            The dictionary examples below show protein entries as if the network has been annotated
            with protein abundances; these are absent for genomic networks lacking protein
            annotations.

            If this method is NOT called from the method, '_purge_kos', then the dictionary will
            look like the following.
            {
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>]
            }

            If this method is called from the method, '_purge_kos', then the dictionary will look
            like the following.
            {
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

            If no genes are removed from the network, then the dictionary will look like the
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
        assert genes_to_remove or proteins_to_remove

        if genes_to_remove is None:
            genes_to_remove: List[int] = []
        if proteins_to_remove is None:
            proteins_to_remove: List[str] = []

        # Get genes to remove from requested proteins.
        for protein_id in proteins_to_remove:
            try:
                protein = self.proteins[protein_id]
            except KeyError:
                # The requested protein ID is not in the network.
                continue
            genes_to_remove += protein.gcids

        # Remove requested genes from the network.
        genes_to_remove = set(genes_to_remove)
        removed_genes: List[Gene] = []
        for gcid in genes_to_remove:
            try:
                removed_genes.append(self.genes.pop(gcid))
            except KeyError:
                # This occurs if the gene in 'genes_to_remove' is not in the network.
                pass

        if not removed_genes:
            removed = {
                'metabolite': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene': []
            }
            if self.proteins:
                removed['protein'] = []
            return removed

        # Purge KOs from the network that are exclusive to removed genes.
        kos_to_remove: List[str] = []
        for gene in removed_genes:
            for ko_id in gene.ko_ids:
                kos_to_remove.append(ko_id)
        kos_to_remove = list(set(kos_to_remove))
        for gene in self.genes.values():
            kos_to_spare: List[str] = []
            for ko_id in gene.ko_ids:
                if ko_id in kos_to_remove:
                    # The KO is not removed because it is associated with a retained gene.
                    kos_to_spare.append(ko_id)
            for ko_id in kos_to_spare:
                kos_to_remove.remove(ko_id)
        if kos_to_remove:
            removed_cascading_down = self._purge_kos(kos_to_remove)
            removed_cascading_down.pop('gene')
        else:
            # This method must have been called from the method, '_purge_kos', because the KOs that
            # are only associated with the removed genes were already removed from the network.
            removed_cascading_down = {
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        if not self.proteins:
            removed = {'gene': removed_genes}
            removed.update(removed_cascading_down)
            return removed

        # Purge protein abundance annotations from the network that are exclusive to removed genes.
        removed_gcids: List[str] = []
        proteins_to_remove: List[str] = []
        for protein_id, protein in self.proteins.items():
            for gcid in protein.gcids:
                if gcid not in removed_gcids:
                    # The protein is not removed because it is associated with a retained gene.
                    break
            else:
                proteins_to_remove.append(protein_id)
        removed_proteins: List[Protein] = []
        for protein_id in proteins_to_remove:
            removed_proteins.append(self.proteins.pop(protein_id))

        removed = {'gene': removed_genes, 'protein': removed_proteins}
        removed.update(removed_cascading_down)
        return removed

    def subset_network(
        self,
        genes_to_subset: Union[int, Iterable[int]] = None,
        proteins_to_subset: Union[str, Iterable[str]] = None,
        kos_to_subset: Union[str, Iterable[str]] = None,
        modules_to_subset: Union[str, Iterable[str]] = None,
        pathways_to_subset: Union[str, Iterable[str]] = None,
        hierarchies_to_subset: Union[str, Iterable[str]] = None,
        categories_to_subset: Dict[str, List[Tuple[str]]] = None,
        reactions_to_subset: Union[str, Iterable[str]] = None,
        metabolites_to_subset: Union[str, Iterable[str]] = None,
        inclusive: bool = False
    ) -> GenomicNetwork:
        """
        Subset a smaller network from the metabolic network.

        If requested genes, proteins, KOs, KEGG modules, KEGG pathways, KEGG BRITE hierarchies, KEGG
        BRITE hierarchy categories, reactions, or metabolites are not present in the network, no
        error is raised.

        Subsetted items are not represented by the same objects as in the source network, i.e., new
        gene, KO, reaction, metabolite, and other objects are created and added to the subsetted
        network.

        Network items (e.g., genes, KOs, reactions, and metabolites) that are associated with
        requested items (e.g., genes in the network that reference requested KOs; metabolites
        referenced by requested reactions) are added to the subsetted network.

        The choice of "inclusive" or, by default, "exclusive" subsetting determines which associated
        items are included in the subsetted network. In exclusive subsetting, a gene added to the
        subsetted network due to references to requested KOs will be missing its references to any
        other unrequested KOs in the source network. Likewise, genes and KOs that are added to the
        subsetted network due to references to requested reactions will be missing references to any
        other unrequested reactions. In other words, certain KO and reaction annotations can be
        selected to the exclusion of others, e.g., a KO encoding two reactions can be restricted to
        encode one requested reaction in the subsetted network; a KO encoding multiple reactions can
        be restricted to encode only those reactions involving requested metabolites.

        "Inclusive" subsetting applies a "Midas touch" where all items in the network that are
        however associated with requested KOs, reactions, and metabolites are "turned to gold" and
        included in the subsetted network. A gene added to the subsetted network due to references
        to requested KOs will still include all of its other references to unrequested KOs in the
        source network. Likewise, KOs and genes that are added to the subsetted network due to
        references to requested reactions and metabolites will include all their other references to
        unrequested reactions and metabolites. Inclusive subsetting precludes the emendation of gene
        KO annotations and KO reaction annotations.

        Parameters
        ==========
        genes_to_subset : Union[int, Iterable[int]], None
            Gene callers ID(s) to subset.

        proteins_to_subset : Union[str, Iterable[str]], None
            Protein ID(s) to subset if the network has been annotated with protein abundances,
            equivalent to giving the genes encoding the protein(s) to the argument,
            'genes_to_subset'.

        kos_to_subset : Union[str, Iterable[str]], None
            KO ID(s) to subset.

        modules_to_subset : Union[str, Iterable[str]], None
            KEGG module ID(s) to subset, with the effect of giving the KOs in the module(s) to the
            argument, 'kos_to_subset'. This does not exclude other module annotations of these KOs
            from the network.

        pathways_to_subset : Union[str, Iterable[str]], None
            KEGG pathway ID(s) to subset, with the effect of giving the KOs in the pathway(s) to the
            argument, 'kos_to_subset'. This does not exclude other pathway annotations of these KOs
            from the network.

        hierarchies_to_subset : Union[str, Iterable[str]], None
            KEGG BRITE hierarchy (or hierarchies) to subset, with the effect of giving the KOs in
            the hierarchy to the argument, 'kos_to_subset'. This does not exclude other hierarchy
            annotations of these KOs from the network.

        categories_to_subset : Dict[str, List[Tuple[str]]], None
            KEGG BRITE hierarchy categories to subset, with the effect of giving the KOs in the
            categories to the argument, 'kos_to_subset'. This does not exclude other category
            annotations of these KOs from the network. The dictionary argument is keyed by BRITE
            hierarchy ID and has values that list category tuples. For example, to subset KOs from
            the network contained in the 'ko00001' 'KEGG Orthology (KO)' hierarchy categories,
            '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis /
            Gluconeogenesis [PATH:ko00010]' and '09100 Metabolism >>> 09101 Carbohydrate
            metabolism >>> 00051 Fructose and mannose metabolism [PATH:ko00051]', the dictionary
            argument would need to look like the following: {'ko00001': [('09100 Metabolism', '09101
            Carbohydrate metabolism', '00010 Glycolysis / Gluconeogenesis'), ('09100 Metabolism',
            '09101 Carbohydrate metabolism', '00051 Fructose and mannose metabolism
            [PATH:ko00051]')]}

        reactions_to_subset : Union[str, Iterable[str]], None
            ModelSEED reaction ID(s) to subset.

        metabolites_to_subset : Union[str, Iterable[str]], None
            ModelSEED compound ID(s) to subset.

        inclusive : bool, False
            If True, "inclusive" subsetting applies a "Midas touch" where all items in the network
            that are however associated with requested KOs, reactions, and metabolites are "turned
            to gold" and included in the subsetted network. In default "exclusive" subsetting, a
            gene added to the subsetted network due to references to requested KOs will be missing
            its references to any other unrequested KOs in the source network; KOs and genes that
            are added to the subsetted network due to references to requested reactions and
            metabolites will be missing references to any other unrequested reactions and
            metabolites.

        Returns
        =======
        GenomicNetwork
            New subsetted reaction network.
        """
        assert (
            genes_to_subset or
            proteins_to_subset or
            kos_to_subset or
            modules_to_subset or
            pathways_to_subset or
            hierarchies_to_subset or
            categories_to_subset or
            reactions_to_subset or
            metabolites_to_subset
        )

        if genes_to_subset is None:
            genes_to_subset: List[int] = []
        if proteins_to_subset is None:
            proteins_to_subset: List[str] = []

        # Get genes to subset from requested proteins.
        for protein_id in proteins_to_subset:
            try:
                protein = self.proteins[protein_id]
            except KeyError:
                # The requested protein ID is not in the network.
                continue
            genes_to_subset += protein.gcids

        if kos_to_subset is None:
            kos_to_subset: List[str] = []
        else:
            kos_to_subset = list(kos_to_subset)
        if modules_to_subset is None:
            modules_to_subset: List[str] = []
        if pathways_to_subset is None:
            pathways_to_subset: List[str] = []
        if hierarchies_to_subset is None:
            hierarchies_to_subset: List[str] = []
        if categories_to_subset is None:
            categories_to_subset: Dict[str, List[Tuple[str]]] = {}

        # Get KOs to subset from requested modules, pathways, hierarchies, and hierarchy categories.
        for module_id in modules_to_subset:
            try:
                module = self.modules[module_id]
            except KeyError:
                # The requested module is not in the network.
                continue
            kos_to_subset += module.ko_ids
        for pathway_id in pathways_to_subset:
            try:
                pathway = self.pathways[pathway_id]
            except KeyError:
                # The requested pathway is not in the network.
                continue
            kos_to_subset += pathway.ko_ids
        for hierarchy_id in hierarchies_to_subset:
            try:
                hierarchy = self.hierarchies[hierarchy_id]
            except KeyError:
                # The requested hierarchy is not in the network.
                continue
            kos_to_subset += hierarchy.ko_ids
        for hierarchy_id, categorizations in categories_to_subset.items():
            hierarchy_categorizations = self.categories[hierarchy_id]
            for categorization in categorizations:
                try:
                    categories = hierarchy_categorizations[categorization]
                except KeyError:
                    # The requested category is not in the network.
                    continue
                category = categories[-1]
                kos_to_subset += category.ko_ids
        kos_to_subset = set(kos_to_subset)

        # Sequentially subset the network for each type of request. Upon generating two subsetted
        # networks from two types of request, merge the networks into a single subsetted network;
        # repeat.
        first_subnetwork = None
        for items_to_subset, subset_network_method in (
            (genes_to_subset, self._subset_network_by_genes),
            (kos_to_subset, functools.partial(self._subset_network_by_kos, inclusive=inclusive)),
            (reactions_to_subset, functools.partial(
                self._subset_network_by_reactions, inclusive=inclusive
            )),
            (metabolites_to_subset, functools.partial(
                self._subset_network_by_metabolites, inclusive=inclusive
            ))
        ):
            if not items_to_subset:
                continue

            second_subnetwork = subset_network_method(items_to_subset)

            if first_subnetwork is None:
                first_subnetwork = second_subnetwork
            else:
                first_subnetwork = first_subnetwork.merge_network(second_subnetwork)

        return first_subnetwork

    def _subset_network_by_genes(self, gcids: Iterable[int]) -> GenomicNetwork:
        """
        Subset the network by genes.

        Parameters
        ==========
        gcids : Iterable[int]
            Gene callers IDs to subset.

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

            # Subset KOs annotating the gene.
            self._subset_network_by_kos(gene.ko_ids, subnetwork=subnetwork)

            subnetwork.genes[gcid] = deepcopy(gene)

            # Only include protein abundances of subsetted genes, ignoring references to unsubsetted
            # genes not encoding the protein.
            if gene.protein_id is not None:
                try:
                    subnetwork_protein = subnetwork.proteins[gene.protein_id]
                except KeyError:
                    protein = self.proteins[gene.protein_id]
                    subnetwork_protein = Protein(id=protein.id, abundances=protein.abundances)
                    subnetwork.proteins[gene.protein_id] = subnetwork_protein
                subnetwork_protein.gcids.append(gcid)

        return subnetwork

    def _subset_genes_via_kos(
        self,
        subnetwork: GenomicNetwork,
        inclusive: bool = False
    ) -> None:
        """
        Add genes that are annotated with subsetted KOs to the subsetted network.

        These gene objects only reference subsetted KOs and not other KOs that also annotate the
        gene but are not subsetted.

        Parameters
        ==========
        subnetwork : GenomicNetwork
            The subsetted reaction network under construction.

        inclusive : bool, False
            If True, "inclusive" subsetting applies a "Midas touch" where all items in the network
            that are however associated with requested KOs are "turned to gold" and included in the
            subsetted network. In default "exclusive" subsetting, a gene added to the subsetted
            network due to references to requested KOs will be missing its references to any other
            unrequested KOs in the source network.

        Returns
        =======
        None
        """
        for gcid, gene in self.genes.items():
            # Check all genes in the source network for subsetted KOs.
            subsetted_ko_ids: List[str] = []
            for ko_id in gene.ko_ids:
                if ko_id in subnetwork.kos:
                    # The gene is annotated by the subsetted KO.
                    subsetted_ko_ids.append(ko_id)
            if not subsetted_ko_ids:
                # The gene is not annotated by any subsetted KOs.
                continue

            if inclusive:
                # Copy the gene, including all its references, to the subsetted network.
                subnetwork.genes[gcid] = deepcopy(gene)

                # Only include protein abundances of subsetted genes, ignoring references to
                # unsubsetted genes not encoding the protein.
                if gene.protein_id is not None:
                    try:
                        subnetwork_protein = subnetwork.proteins[gene.protein_id]
                    except KeyError:
                        protein = self.proteins[gene.protein_id]
                        subnetwork_protein = Protein(id=protein.id, abundances=protein.abundances)
                        subnetwork.proteins[gene.protein_id] = subnetwork_protein
                    subnetwork_protein.gcids.append(gcid)

                # Add "unrequested" KOs associated with the gene to the subsetted network if not
                # already added.
                for ko_id in gene.ko_ids:
                    if ko_id in subsetted_ko_ids:
                        continue
                    ko = self.kos[ko_id]
                    subnetwork.kos[ko_id] = deepcopy(ko)

                    # Add reactions associated with the unrequested KO to the subsetted network if
                    # not already added.
                    for reaction_id in ko.reaction_ids:
                        if reaction_id in subnetwork.reactions:
                            continue
                        reaction = self.reactions[reaction_id]
                        self._subset_reaction(subnetwork, reaction)
                continue

            # Subsetting is exclusive, not inclusive. Add the gene only with references to subsetted
            # KOs.
            subnetwork_gene = Gene(gcid=gcid, protein_id=gene.protein_id)

            if gene.protein_id is not None:
                try:
                    subnetwork_protein = subnetwork.proteins[gene.protein_id]
                except KeyError:
                    protein = self.proteins[gene.protein_id]
                    subnetwork_protein = Protein(id=protein.id, abundances=protein.abundances)
                    subnetwork.proteins[gene.protein_id] = subnetwork_protein
                subnetwork_protein.gcids.append(gcid)

            for ko_id in subsetted_ko_ids:
                subnetwork_gene.ko_ids.append(ko_id)
                subnetwork_gene.e_values[ko_id] = gene.e_values[ko_id]

    def merge_network(self, network: GenomicNetwork) -> GenomicNetwork:
        """
        Merge the genomic reaction network with another genomic reaction network derived from the
        same contigs database.

        The purpose of the network is to combine different, but potentially overlapping, subnetworks
        from the same genome.

        Each network can contain different genes, KOs, and reactions/metabolites. Merging
        nonredundantly incorporates all of this data as new objects in the new network.

        Objects representing genes, KOs, KEGG modules, KEGG pathways, BRITE hierarchies, and
        hierarchy categories in both networks can have different sets of references: genes can be
        annotated by different KOs; KOs can be annotated by different reactions; depending on the
        KOs in each network, different modules, pathways, and hierarchies/categories can be present.

        Other object attributes should be consistent between the networks. For instance, the same
        ModelSEED reactions and metabolites in both networks should have identical attributes. The
        same gene-KO should have the same e-values. If applicable, both networks should have been
        annotated with the same protein and metabolite abundance data.

        Parameters
        ==========
        network : GenomicNetwork
            The other genomic reaction network being merged.

        Returns
        =======
        GenomicNetwork
            The merged genomic reaction network.
        """
        merged_network: GenomicNetwork = self._merge_network(network)

        merged_network.genes = deepcopy(self.genes)
        merged_network.proteins = deepcopy(self.proteins)

        # Copy genes from the second network. Assume the same gene KO references have the same
        # e-values. Assume identical protein abundance assignments.
        for gcid, gene in network.genes.items():
            try:
                merged_gene = merged_network.genes[gcid]
            except KeyError:
                merged_network.genes[gcid] = deepcopy(gene)
                continue

            merged_gene.ko_ids = sorted(set(gene.ko_ids + merged_gene.ko_ids))

            for ko_id, e_value in gene.e_values.items():
                if ko_id in merged_gene.e_values:
                    continue

                merged_gene.e_values[ko_id] = e_value

        # Copy proteins from the second network. Assume identical abundance profiles.
        merged_proteins = merged_network.proteins
        for protein_id, protein in network.proteins.items():
            try:
                merged_protein = merged_proteins[protein_id]
            except KeyError:
                merged_proteins[protein_id] = deepcopy(protein)
                continue

            merged_protein.gcids = sorted(set(protein.gcids + merged_protein.gcids))

        return merged_network

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
            this dictionary must contain certain precomputed data: the key, 'total_genes', should
            have a value of the number of genes in the genome; the key, 'kos_assigned_genes', should
            have a value of the number of genes in the genome that are assigned KOs; the key,
            'kos_assigned_genes', should have a value of the number of unique KOs assigned to genes
            in the genome.

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
                "Since the genomic network was not associated with a contigs database, the "
                "following statistics could not be calculated and were not reported to the output "
                "file: 'Total gene calls in genome', 'Genes annotated with protein KOs', and "
                "'Protein KOs annotating genes'."
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
        indent: int = 2,
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Export the network to a metabolic model file in JSON format.

        All information from the network is included in the JSON so that the file can be loaded as a
        GenomicNetwork object containing the same information.

        Parameters
        ==========
        path : str
            Output JSON file path.

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

        indent : int, 2
            Spaces of indentation per nesting level in JSON file.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            Prints transient progress information to the terminal.
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
            raise ConfigError(
                f"Anvi'o does not recognize an objective with the name, '{objective}'."
            )

        progress.update("Genes")
        reaction_genes: Dict[str, List[str]] = {}
        reaction_kos: Dict[str, List[KO]] = {}
        for gcid, gene in self.genes.items():
            gene_entry = JSONStructure.get_gene_entry()
            json_genes.append(gene_entry)
            gcid_str = str(gcid)
            gene_entry['id'] = gcid_str

            # Record KO IDs, annotation e-values, and KO classifications in the annotation section
            # of the gene entry.
            annotation = gene_entry['annotation']
            annotation['ko'] = annotation_kos = {}
            for ko_id in gene.ko_ids:
                annotation_kos[ko_id] = annotation_ko = {
                    'e_value': str(gene.e_values[ko_id]),
                    'modules': {},
                    'pathways': {},
                    'hierarchies': {}
                }

                # Record KEGG modules containing the KO.
                ko = self.kos[ko_id]
                annotation_ko_modules = annotation_ko['modules']
                for module_id in ko.module_ids:
                    module = self.modules[module_id]
                    module_annotation = module.name
                    if not module.pathway_ids:
                        annotation_ko_modules[module_id] = module_annotation
                        continue
                    # Cross-reference KEGG pathways containing the module.
                    module_annotation += "[pathways:"
                    for pathway_id in module.pathway_ids:
                        module_annotation += f" {pathway_id}"
                    module_annotation += "]"
                    annotation_ko_modules[module_id] = module_annotation

                # Record KEGG pathways containing the KO.
                annotation_ko_pathways = annotation_ko['pathways']
                for pathway_id in ko.pathway_ids:
                    pathway = self.pathways[pathway_id]
                    annotation_ko_pathways[pathway_id] = pathway.name

                # Record membership of the KO in KEGG BRITE hierarchies.
                annotation_ko_hierarchies: Dict[str, List[str]] = annotation_ko['hierarchies']
                for hierarchy_id, categorizations in ko.hierarchies.items():
                    hierarchy_name = self.hierarchies[hierarchy_id].name
                    annotation_ko_hierarchies[
                        f"{hierarchy_id}: {hierarchy_name}"
                    ] = annotation_ko_categories = []
                    hierarchy_categorizations = self.categories[hierarchy_id]
                    for categorization in categorizations:
                        categories = hierarchy_categorizations[categorization]
                        category = categories[-1]
                        category_id = category.id
                        annotation_ko_categories.append(category_id[len(hierarchy_id) + 2:])

                # Set up dictionaries needed to fill out reaction entries.
                for reaction_id in ko.reaction_ids:
                    try:
                        reaction_genes[reaction_id].append(gcid_str)
                    except KeyError:
                        reaction_genes[reaction_id] = [gcid_str]
                    try:
                        reaction_kos[reaction_id].append(ko)
                    except KeyError:
                        reaction_kos[reaction_id] = [ko]

                if not self.proteins:
                    continue

                # A protein section is added if the network has been annotated with protein
                # abundances.
                annotation['protein'] = annotation_protein = {
                    'id': None,
                    'abundances': {}
                }
                if gene.protein_id:
                    # Record abundances of the protein encoded by the gene.
                    protein_id = gene.protein_id
                    protein = self.proteins[protein_id]
                    annotation_protein['id'] = protein_id
                    annotation_protein_abundances = annotation_protein['abundances']
                    for sample_name, abundance_value in protein.abundances.items():
                        annotation_protein_abundances[sample_name] = abundance_value

        progress.update("Reactions")
        compound_compartments: Dict[str, Set[str]] = {}
        for reaction_id, reaction in self.reactions.items():
            reaction_entry = JSONStructure.get_reaction_entry()
            json_reactions.append(reaction_entry)
            reaction_entry['id'] = reaction_id
            reaction_entry['name'] = reaction.modelseed_name
            metabolites = reaction_entry['metabolites']
            for compound_id, compartment, coefficient in zip(
                reaction.compound_ids, reaction.compartments, reaction.coefficients
            ):
                metabolites[f"{compound_id}_{compartment}"] = coefficient
                try:
                    compound_compartments[compound_id].add(compartment)
                except KeyError:
                    compound_compartments[compound_id] = set(compartment)
            if not reaction.reversibility:
                # By default, the reaction entry was set up to be reversible; here make it
                # irreversible.
                reaction_entry['lower_bound'] = 0.0
            reaction_entry['gene_reaction_rule'] = " or ".join(
                [gcid for gcid in reaction_genes[reaction_id]]
            )

            notes = reaction_entry['notes']
            # Record gene KO annotations which aliased the reaction via KEGG REACTION or EC number.
            notes['ko'] = ko_notes = {}
            ko_kegg_aliases = []
            ko_ec_number_aliases = []
            for ko in reaction_kos[reaction_id]:
                try:
                    kegg_aliases = ko.kegg_reaction_aliases[reaction_id]
                except KeyError:
                    kegg_aliases = []
                try:
                    ec_number_aliases = ko.ec_number_aliases[reaction_id]
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
        for compound_id, metabolite in self.metabolites.items():
            modelseed_compound_name = metabolite.modelseed_name
            charge = metabolite.charge
            formula = metabolite.formula
            kegg_compound_aliases = list(metabolite.kegg_aliases)
            for compartment in compound_compartments[compound_id]:
                metabolite_entry = JSONStructure.get_metabolite_entry()
                json_metabolites.append(metabolite_entry)
                metabolite_entry['id'] = f"{compound_id}_{compartment}"
                metabolite_entry['name'] = modelseed_compound_name
                metabolite_entry['compartment'] = compartment
                # Compounds without a formula have a nominal charge of 10000000 in the ModelSEED
                # compounds database, which is replaced by None in the reaction network and 0 in the
                # JSON.
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

    def remove_metabolites_without_formula(self, output_path: str = None) -> None:
        """
        Remove metabolites without a formula in the ModelSEED database from the network.

        Other items can be removed from the network by association: reactions that involve a
        formulaless metabolite; other metabolites with formulas that are exclusive to such
        reactions; KOs predicted to exclusively catalyze such reactions; and gene clusters annotated
        with such KOs. Removed metabolites with a formula are reported alongside formulaless
        metabolites to the output table of removed metabolites.

        output_path : str, None
            If not None, write tab-delimited files of metabolites, reactions, KOs, KEGG modules,
            KEGG pathways, KEGG BRITE hierarchies, KEGG BRITE hierarchy categories, and gene
            clusters removed from the network to file locations based on the provided path. For
            example, if the argument, 'removed.tsv', is provided, then the following files will be
            written: 'removed-metabolites.tsv', 'removed-reactions.tsv', 'removed-kos.tsv',
            'removed-modules.tsv', 'removed-pathways.tsv', 'removed-hierarchies.tsv',
            'removed-categories.tsv', and 'removed-gene-clusters.tsv'.
        """
        if self.verbose:
            self.progress.new("Removing metabolites without a formula in the network")
            self.progress.update("...")

        if output_path:
            path_basename, path_extension = os.path.splitext(output_path)
            metabolite_path = f"{path_basename}-metabolites{path_extension}"
            reaction_path = f"{path_basename}-reactions{path_extension}"
            ko_path = f"{path_basename}-kos{path_extension}"
            module_path = f"{path_basename}-modules{path_extension}"
            pathway_path = f"{path_basename}-pathways{path_extension}"
            hierarchy_path = f"{path_basename}-hierarchies{path_extension}"
            category_path = f"{path_basename}-categories{path_extension}"
            gene_cluster_path = f"{path_basename}-gene-clusters{path_extension}"
            for path in (
                metabolite_path,
                reaction_path,
                ko_path,
                module_path,
                pathway_path,
                hierarchy_path,
                category_path,
                gene_cluster_path
            ):
                filesnpaths.is_output_file_writable(path)

        metabolites_to_remove = []
        for compound_id, metabolite in self.metabolites.items():
            # ModelSEED compounds without a formula have a formula value of None in the network
            # object.
            if metabolite.formula is None:
                metabolites_to_remove.append(compound_id)
        removed = self._purge_metabolites(metabolites_to_remove)

        if self.verbose:
            self.run.info("Removed metabolites", len(removed['metabolite']))
            self.run.info("Removed reactions", len(removed['reaction']))
            self.run.info("Removed KOs", len(removed['ko']))
            self.run.info("Removed KEGG modules", len(removed['module']))
            self.run.info("Removed KEGG pathways", len(removed['pathway']))
            self.run.info("Removed KEGG BRITE hierarchies", len(removed['hierarchy']))
            self.run.info("Removed KEGG BRITE hierarchy categories", len(removed['category']))
            self.run.info("Removed gene clusters", len(removed['gene_cluster']))

        if not output_path:
            return

        if self.verbose:
            self.progress.new("Writing output files of removed network items")
            self.progress.update("...")

        gene_cluster_table = []
        for cluster in removed['gene_cluster']:
            cluster: GeneCluster
            row = []
            row.append(cluster.gene_cluster_id)
            row.append(cluster.ko_id)
            row.append(", ".join(cluster.genomes))
            gene_cluster_table.append(row)

        self._write_remove_metabolites_without_formula_output(output_path, removed)

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
            self.run.info("Table of removed KEGG modules", module_path)
            self.run.info("Table of removed KEGG pathways", pathway_path)
            self.run.info("Table of removed KEGG BRITE hierarchies", hierarchy_path)
            self.run.info("Table of removed KEGG BRITE hierarchy categories", category_path)
            self.run.info("Table of removed gene clusters", gene_cluster_path)

    def prune(
        self,
        gene_clusters_to_remove: Union[int, Iterable[int]] = None,
        kos_to_remove: Union[str, Iterable[str]] = None,
        modules_to_remove: Union[str, Iterable[str]] = None,
        pathways_to_remove: Union[str, Iterable[str]] = None,
        hierarchies_to_remove: Union[str, Iterable[str]] = None,
        categories_to_remove: Dict[str, List[Tuple[str]]] = None,
        reactions_to_remove: Union[str, Iterable[str]] = None,
        metabolites_to_remove: Union[str, Iterable[str]] = None
    ) -> Dict[str, List]:
        """
        Prune items from the metabolic network.

        Pruning modifies the network in situ: use the network 'copy' method as needed to create a
        backup of the network.

        If requested gene clusters, KOs, KEGG modules, KEGG pathways, KEGG BRITE hierarchies, KEGG
        BRITE hierarchy categories, reactions, or metabolites are not present in the network, no
        error is raised.

        Network items (e.g., gene clusters, KOs, reactions, and metabolites) that are exclusively
        associated with requested items are also removed from the network. Example: Consider a KO
        that is requested to be removed from the network. The KO is associated with two reactions.
        The first reaction is exclusive to the KO and thus is also removed, whereas the second
        reaction is also associated with another retained KO and thus is retained in the network.
        The first reaction involves four metabolites, and two are exclusive to the reaction: these
        are also removed from the network. Each gene cluster has a single consensus KO annotation,
        so any gene clusters assigned this KO are removed from the network. Note that reaction
        annotations of KOs can be selected to the exclusion of others. In the example, the latter
        gene cluster is left with one KO.

        Parameters
        ==========
        gene_clusters_to_remove : Union[str, Iterable[int]], None
            Gene cluster ID(s) to remove.

        kos_to_remove : Union[str, Iterable[str]], None
            KO ID(s) to remove.

        modules_to_remove : Union[str, Iterable[str]], None
            KEGG module ID(s) to remove, with the effect of giving the KOs in the module(s) to the
            argument, 'kos_to_remove'. This does not remove other module annotations of these KOs
            that also annotate other KOs.

        pathways_to_remove : Union[str, Iterable[str]], None
            KEGG pathway ID(s) to remove, with the effect of giving the KOs in the pathway(s) to the
            argument, 'kos_to_remove'. This does not remove other pathway annotations of these KOs
            that also annotate other KOs.

        hierarchies_to_remove : Union[str, Iterable[str]], None
            KEGG BRITE hierarchy (or hierarchies) to remove, with the effect of giving the KOs in
            the hierarchy to the argument, 'kos_to_remove'. This does not remove other hierarchy
            annotations of these KOs that also annotate other KOs.

        categories_to_remove : Dict[str, List[Tuple[str]]], None
            KEGG BRITE hierarchy categories to remove, with the effect of giving the KOs in the
            categories to the argument, 'kos_to_remove'. This does not remove other category
            annotations of these KOs that also annotate other KOs. The dictionary argument is keyed
            by BRITE hierarchy ID and has values that list category tuples. For example, to remove
            KOs contained in the 'ko00001' 'KEGG Orthology (KO)' hierarchy categories, '09100
            Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis / Gluconeogenesis
            [PATH:ko00010]' and '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00051
            Fructose and mannose metabolism [PATH:ko00051]', the dictionary argument would need to
            look like the following: {'ko00001': [('09100 Metabolism', '09101 Carbohydrate
            metabolism', '00010 Glycolysis / Gluconeogenesis'), ('09100 Metabolism', '09101
            Carbohydrate metabolism', '00051 Fructose and mannose metabolism [PATH:ko00051]')]}

        reactions_to_remove : Union[str, Iterable[str]], None
            ModelSEED reaction ID(s) to remove.

        metabolites_to_remove : Union[str, Iterable[str]], None
            ModelSEED compound ID(s) to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            The dictionary has the following format. It shows protein entries as if the network has
            been annotated with protein abundances; these are absent for genomic networks lacking
            protein annotations.
            {
                'gene': [<removed Gene objects>],
                'protein': [<removed Protein objects>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>]
            }
        """
        assert (
            gene_clusters_to_remove or
            kos_to_remove or
            modules_to_remove or
            pathways_to_remove or
            hierarchies_to_remove or
            categories_to_remove or
            reactions_to_remove or
            metabolites_to_remove
        )

        removed: Dict[str, List] = {
            'gene_cluster': [],
            'ko': [],
            'module': [],
            'pathway': [],
            'hierarchy': [],
            'category': [],
            'reaction': [],
            'kegg_reaction': [],
            'ec_number': [],
            'metabolite': []
        }

        if gene_clusters_to_remove:
            for item_type, removed_items in self._purge_gene_clusters(
                gene_clusters_to_remove=gene_clusters_to_remove
            ).items():
                removed[item_type] += removed_items

        if (
            kos_to_remove or
            modules_to_remove or
            pathways_to_remove or
            hierarchies_to_remove or
            categories_to_remove
        ):
            for item_type, removed_items in self._purge_kos(
                kos_to_remove=kos_to_remove,
                modules_to_remove=modules_to_remove,
                pathways_to_remove=pathways_to_remove,
                hierarchies_to_remove=hierarchies_to_remove,
                categories_to_remove=categories_to_remove
            ).items():
                removed[item_type] += removed_items

        if reactions_to_remove:
            for item_type, removed_items in self._purge_reactions(reactions_to_remove).items():
                removed[item_type] += removed_items

        if metabolites_to_remove:
            for item_type, removed_items in self._purge_metabolites(metabolites_to_remove).items():
                removed[item_type] += removed_items

        return removed

    def _purge_gene_clusters(self, gene_clusters_to_remove: Iterable[str]) -> Dict[str, List]:
        """
        Remove any trace of the given gene clusters from the network.

        KOs, reactions, and metabolites that are only associated with removed gene clusters are
        purged. KEGG modules, pathways, BRITE hierarchies, and BRITE hierarchy categories only
        associated with purged KOs are removed.

        Parameters
        ==========
        gene_clusters_to_remove : Iterable[str]
            Gene cluster IDs to remove.

        Returns
        =======
        dict
            This dictionary contains data removed from the network.

            If this method is NOT called from the method, '_purge_kos', then the dictionary will
            look like the following.
            {
                'gene_cluster': [<removed GeneCluster objects>],
                'ko': [<removed KO objects>],
                'module': [<removed KEGGModule objects>],
                'pathway': [<removed KEGGPathway objects>],
                'hierarchy': [<removed BRITEHierarchy objects>],
                'category': [<removed BRITECategory objects>],
                'reaction': [<removed ModelSEEDReaction objects>],
                'kegg_reaction': [<removed KEGG REACTION IDs>],
                'ec_number': [<removed EC numbers>],
                'metabolite': [<removed ModelSEEDCompound objects>]
            }

            If this method is called from the method, '_purge_kos', then the dictionary will look
            like the following.
            {
                'gene_cluster': [<removed GeneCluster objects>],
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

            If no gene clusters are removed from the network, then the dictionary will look like the
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
                'gene_cluster': []
            }
        """
        gene_clusters_to_remove = set(gene_clusters_to_remove)
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
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'gene_cluster': []
            }

        # Purge KOs from the network that are exclusively assigned to removed gene clusters.
        kos_to_remove: List[str] = []
        for cluster in removed_gene_clusters:
            kos_to_remove.append(cluster.ko_id)
        kos_to_remove = list(set(kos_to_remove))
        for gene_cluster in self.gene_clusters.values():
            kos_to_spare: List[str] = []
            if gene_cluster.ko_id in kos_to_remove:
                # The KO is associated with a retained gene cluster, so do not remove the KO.
                kos_to_spare.append(gene_cluster.ko_id)
            for ko_id in kos_to_spare:
                kos_to_remove.remove(ko_id)
        if kos_to_remove:
            removed_cascading_down = self._purge_kos(kos_to_remove)
            removed_cascading_down.pop('gene_cluster')
        else:
            # This method must have been called from the method, '_purge_kos', because the KOs that
            # are only associated with the removed gene clusters were already removed from the
            # network.
            removed_cascading_down = {
                'ko': [],
                'module': [],
                'pathway': [],
                'hierarchy': [],
                'category': [],
                'reaction': [],
                'kegg_reaction': [],
                'ec_number': [],
                'metabolite': []
            }

        removed = {'gene_cluster': removed_gene_clusters}
        removed.update(removed_cascading_down)
        return removed

    def subset_network(
        self,
        gene_clusters_to_subset: Union[int, Iterable[int]] = None,
        kos_to_subset: Union[str, Iterable[str]] = None,
        modules_to_subset: Union[str, Iterable[str]] = None,
        pathways_to_subset: Union[str, Iterable[str]] = None,
        hierarchies_to_subset: Union[str, Iterable[str]] = None,
        categories_to_subset: Dict[str, List[Tuple[str]]] = None,
        reactions_to_subset: Union[str, Iterable[str]] = None,
        metabolites_to_subset: Union[str, Iterable[str]] = None,
        inclusive: bool = False
    ) -> PangenomicNetwork:
        """
        Subset a smaller network from the metabolic network.

        If requested gene clusters, KOs, KEGG modules, KEGG pathways, KEGG BRITE hierarchies, KEGG
        BRITE hierarchy categories, reactions, or metabolites are not present in the network, no
        error is raised.

        Subsetted items are not represented by the same objects as in the source network, i.e., new
        gene cluster, KO, reaction, metabolite, and other objects are created and added to the
        subsetted network.

        Network items (e.g., gene clusters, KOs, reactions, and metabolites) that are associated with
        requested items (e.g., gene clusters in the network that reference requested KOs; metabolites
        referenced by requested reactions) are added to the subsetted network.

        The choice of "inclusive" or, by default, "exclusive" subsetting determines which associated
        items are included in the subsetted network. In exclusive subsetting, KOs (and by extension,
        gene clusters assigned the KOs) that are added to the subsetted network due to references to
        requested reactions will be missing references to any other unrequested reactions. In other
        words, certain reaction annotations can be selected to the exclusion of others, e.g., a KO
        encoding two reactions can be restricted to encode one requested reaction in the subsetted
        network; a KO encoding multiple reactions can be restricted to encode only those reactions
        involving requested metabolites.

        "Inclusive" subsetting applies a "Midas touch" where all items in the network that are
        however associated with requested reactions and metabolites are "turned to gold" and
        included in the subsetted network. KOs and gene clusters that are added to the subsetted
        network due to references to requested reactions and metabolites will include all their
        other references to unrequested reactions and metabolites. Inclusive subsetting precludes
        the emendation of KO reaction annotations.

        Parameters
        ==========
        gene_clusters_to_subset : Union[int, Iterable[int]], None
            Gene cluster ID(s) to subset.

        kos_to_subset : Union[str, Iterable[str]], None
            KO ID(s) to subset.

        modules_to_subset : List[str], None
            KEGG module ID(s) to subset, with the effect of giving the KOs in the module(s) to the
            argument, 'kos_to_subset'. This does not exclude other module annotations of these KOs
            from the network.

        pathways_to_subset : Union[str, Iterable[str]], None
            KEGG pathway ID(s) to subset, with the effect of giving the KOs in the pathway(s) to the
            argument, 'kos_to_subset'. This does not exclude other pathway annotations of these KOs
            from the network.

        hierarchies_to_subset : Union[str, Iterable[str]], None
            KEGG BRITE hierarchy (or hierarchies) to subset, with the effect of giving the KOs in
            the hierarchy to the argument, 'kos_to_subset'. This does not exclude other hierarchy
            annotations of these KOs from the network.

        categories_to_subset : Dict[str, List[Tuple[str]]], None
            KEGG BRITE hierarchy categories to subset, with the effect of giving the KOs in the
            categories to the argument, 'kos_to_subset'. This does not exclude other category
            annotations of these KOs from the network. The dictionary argument is keyed by BRITE
            hierarchy ID and has values that list category tuples. For example, to subset KOs from
            the network contained in the 'ko00001' 'KEGG Orthology (KO)' hierarchy categories,
            '09100 Metabolism >>> 09101 Carbohydrate metabolism >>> 00010 Glycolysis /
            Gluconeogenesis [PATH:ko00010]' and '09100 Metabolism >>> 09101 Carbohydrate
            metabolism >>> 00051 Fructose and mannose metabolism [PATH:ko00051]', the dictionary
            argument would need to look like the following: {'ko00001': [('09100 Metabolism', '09101
            Carbohydrate metabolism', '00010 Glycolysis / Gluconeogenesis'), ('09100 Metabolism',
            '09101 Carbohydrate metabolism', '00051 Fructose and mannose metabolism
            [PATH:ko00051]')]}

        reactions_to_subset : Union[str, Iterable[str]], None
            ModelSEED reaction ID(s) to subset.

        metabolites_to_subset : Union[str, Iterable[str]], None
            ModelSEED compound ID(s) to subset.

        inclusive : bool, False
            If True, "inclusive" subsetting applies a "Midas touch" where all items in the network
            that are however associated with requested reactions and metabolites are "turned to
            gold" and included in the subsetted network. In default "exclusive" subsetting, KOs and
            gene clusters that are added to the subsetted network due to references to requested
            reactions and metabolites will be missing references to any other unrequested reactions
            and metabolites.

        Returns
        =======
        PangenomicNetwork
            New subsetted reaction network.
        """
        assert (
            gene_clusters_to_subset or
            kos_to_subset or
            modules_to_subset or
            pathways_to_subset or
            hierarchies_to_subset or
            categories_to_subset or
            reactions_to_subset or
            metabolites_to_subset
        )

        if kos_to_subset is None:
            kos_to_subset: List[str] = []
        else:
            kos_to_subset = list(kos_to_subset)
        if modules_to_subset is None:
            modules_to_subset: List[str] = []
        if pathways_to_subset is None:
            pathways_to_subset: List[str] = []
        if hierarchies_to_subset is None:
            hierarchies_to_subset: List[str] = []
        if categories_to_subset is None:
            categories_to_subset: Dict[str, List[Tuple[str]]] = {}

        # Get KOs to subset from requested modules, pathways, hierarchies, and hierarchy categories.
        for module_id in modules_to_subset:
            try:
                module = self.modules[module_id]
            except KeyError:
                # The requested module is not in the network.
                continue
            kos_to_subset += module.ko_ids
        for pathway_id in pathways_to_subset:
            try:
                pathway = self.pathways[pathway_id]
            except KeyError:
                # The requested pathway is not in the network.
                continue
            kos_to_subset += pathway.ko_ids
        for hierarchy_id in hierarchies_to_subset:
            try:
                hierarchy = self.hierarchies[hierarchy_id]
            except KeyError:
                # The requested hierarchy is not in the network.
                continue
            kos_to_subset += hierarchy.ko_ids
        for hierarchy_id, categorizations in categories_to_subset.items():
            try:
                hierarchy_categorizations = self.categories[hierarchy_id]
            except KeyError:
                # The requested hierarchy is not in the network.
                continue
            for categorization in categorizations:
                try:
                    categories = hierarchy_categorizations[categorization]
                except KeyError:
                    # The requested category is not in the network.
                    continue
                category = categories[-1]
                kos_to_subset += category.ko_ids
        kos_to_subset = set(kos_to_subset)

        # Sequentially subset the network for each type of request. Upon generating two subsetted
        # networks from two types of request, merge the networks into a single subsetted network;
        # repeat.
        first_subnetwork = None
        for items_to_subset, subset_network_method in (
            (gene_clusters_to_subset, self._subset_network_by_gene_clusters),
            (kos_to_subset, functools.partial(self._subset_network_by_kos, inclusive=inclusive)),
            (reactions_to_subset, functools.partial(
                self._subset_network_by_reactions, inclusive=inclusive
            )),
            (metabolites_to_subset, functools.partial(
                self._subset_network_by_metabolites, inclusive=inclusive
            ))
        ):
            if not items_to_subset:
                continue

            second_subnetwork = subset_network_method(items_to_subset)

            if first_subnetwork is None:
                first_subnetwork = second_subnetwork
            else:
                first_subnetwork = first_subnetwork.merge_network(second_subnetwork)

        return first_subnetwork

    def _subset_network_by_gene_clusters(
        self,
        gene_cluster_ids: Iterable[int]
    ) -> PangenomicNetwork:
        """
        Subset the network by gene clusters with requested IDs.

        Parameters
        ==========
        gene_cluster_ids : Iterable[int]
            Gene cluster IDs to subset.

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

            # Subset the consensus KO annotating the gene cluster.
            self._subset_network_by_kos([cluster.ko_id], subnetwork=subnetwork)

            subnetwork.gene_clusters[gene_cluster_id] = deepcopy(cluster)

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
            if cluster.ko_id in subsetted_ko_ids:
                # A gene cluster is annotated by the subsetted KO.
                subnetwork.gene_clusters[gene_cluster_id] = deepcopy(cluster)

    def merge_network(self, network: PangenomicNetwork) -> PangenomicNetwork:
        """
        Merge the pangenomic reaction network with another pangenomic reaction network derived from
        the same pan database.

        The purpose of the network is to combine different, but potentially overlapping, subnetworks
        from the same pangenome.

        Each network can contain different gene clusters, KOs, and reactions/metabolites. Merging
        nonredundantly incorporates all of this data as new objects in the new network.

        Objects representing KOs in both networks can have different sets of references: KOs can be
        annotated by different reactions. However, the same gene cluster in each network should have
        the same consensus KO annotation. ModelSEED reactions and metabolites in both networks
        should have identical attributes.

        Parameters
        ==========
        network : PangenomicNetwork
            The other pangenomic reaction network being merged.

        Returns
        =======
        PangenomicNetwork
            The merged pangenomic reaction network.
        """
        merged_network: PangenomicNetwork = self._merge_network(network)

        merged_network.gene_clusters = deepcopy(self.gene_clusters)

        # Copy gene clusters from the second network. Assume they have the same consensus KOs.
        merged_gene_clusters = merged_network.gene_clusters
        for gene_cluster_id, cluster in network.gene_clusters.items():
            if gene_cluster_id not in merged_gene_clusters:
                merged_gene_clusters[gene_cluster_id] = deepcopy(cluster)

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
            dictionary must contain certain precomputed data: the key, 'total_gene_clusters', should
            have a value of the number of gene clusters in the pangenome; the key,
            'gene_clusters_assigned_ko', should have a value of the number of gene clusters in the
            pangenome assigned a consensus KO (or None if 'self.consistent_annotations' is False);
            the key, 'kos_assigned_gene_clusters', should have a value of the number of consensus
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
                "Since the pangenomic network was not associated with a pan database and genomes "
                "storage database, the following statistics could not be calculated and were not "
                "reported to the output file: 'Total gene clusters in pangenome', 'Gene clusters "
                "assigned protein KOs', and 'Protein KOs assigned to gene clusters'."
            )
        elif self.consistent_annotations is False:
            self.run.info_single(
                "The network attribute, 'consistent_annotations', is False, which indicates that "
                "the reaction network stored in the pan database was made from a different set of "
                "KO gene annotations than is currently in the genomes storage database. Therefore, "
                "the following statistics were not calculated and reported to the output file to "
                "avoid potential inaccuracies: 'Gene clusters assigned protein KO' and 'Protein "
                "KOs assigned to gene clusters'."
            )

        return stats

    def print_overview_statistics(self, stats: PangenomicNetworkStats = None) -> None:
        """
        Print overview statistics for the genomic metabolic network.

        Parameters
        ==========
        stats : PangenomicNetworkStats, None
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
        indent: int = 2,
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Export the network to a metabolic model file in JSON format.

        Entries in the "gene" section of this file represent gene clusters.

        All information from the network is included in the JSON so that the file can be loaded as a
        PangenomicNetwork object containing the same information.

        Parameters
        ==========
        path : str
            Output JSON file path.

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
            Spaces of indentation per nesting level in JSON file.

        progress : terminal.Progress, terminal.Progress()
            Prints transient progress information to the terminal.
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
                f"The following items in the 'record_genomes' argument are invalid: "
                f"{', '.join(invalid_items)}"
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
            raise ConfigError(
                f"Anvi'o does not recognize an objective with the name, '{objective}'."
            )

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

            # Record the consensus KO ID and classifications in the annotation section of the gene
            # cluster entry.
            annotation = gene_cluster_entry['annotation']
            ko_id = gene_cluster.ko_id
            annotation['ko'] = annotation_ko = {
                'id': ko_id,
                'modules': {},
                'pathways': {},
                'hierarchies': {}
            }

            # Record KEGG modules containing the KO.
            ko = self.kos[ko_id]
            annotation_ko_modules = annotation_ko['modules']
            for module_id in ko.module_ids:
                module = self.modules[module_id]
                module_annotation = module.name
                if not module.pathway_ids:
                    annotation_ko_modules[module_id] = module_annotation
                    continue
                # Cross-reference KEGG pathways containing the module.
                module_annotation += "[pathways:"
                for pathway_id in module.pathway_ids:
                    module_annotation += f" {pathway_id}"
                module_annotation += "]"
                annotation_ko_modules[module_id] = module_annotation

            # Record KEGG pathways containing the KO.
            annotation_ko_pathways = annotation_ko['pathways']
            for pathway_id in ko.pathway_ids:
                pathway = self.pathways[pathway_id]
                annotation_ko_pathways[pathway_id] = pathway.name

            # Record membership of the KO in KEGG BRITE hierarchies.
            annotation_ko_hierarchies: Dict[str, List[str]] = annotation_ko['hierarchies']
            for hierarchy_id, categorizations in ko.hierarchies.items():
                hierarchy_name = self.hierarchies[hierarchy_id].name
                annotation_ko_hierarchies[
                    f"{hierarchy_id}: {hierarchy_name}"
                ] = annotation_ko_categories = []
                hierarchy_categorizations = self.categories[hierarchy_id]
                for categorization in categorizations:
                    categories = hierarchy_categorizations[categorization]
                    category = categories[-1]
                    category_id = category.id
                    annotation_ko_categories.append(category_id[len(hierarchy_id) + 2:])

            # Set up dictionaries needed to fill out reaction entries.
            for reaction_id in ko.reaction_ids:
                try:
                    reaction_gene_clusters[reaction_id].append(cluster_id_str)
                except KeyError:
                    reaction_gene_clusters[reaction_id] = [cluster_id_str]
                try:
                    reaction_kos[reaction_id].append(ko)
                except KeyError:
                    reaction_kos[reaction_id] = [ko]

            if not record_genomes:
                continue

            genome_names = gene_cluster.genomes
            if 'gene cluster' in record_genomes:
                # Record the names of the genomes contributing to the gene cluster in the notes
                # section of the gene cluster entry.
                gene_cluster_entry['notes']['genomes'] = genome_names
            if 'reaction' in record_genomes:
                for reaction_id in ko.reaction_ids:
                    try:
                        reaction_genomes[reaction_id] += genome_names
                    except KeyError:
                        reaction_genomes[reaction_id] = genome_names
            if 'metabolite' in record_genomes:
                for reaction_id in ko.reaction_ids:
                    reaction = self.reactions[reaction_id]
                    for compartment, compound_id in zip(
                        reaction.compartments, reaction.compound_ids
                    ):
                        entry_id = f"{compound_id}_{compartment}"
                        try:
                            metabolite_genomes[entry_id] += genome_names
                        except KeyError:
                            metabolite_genomes[entry_id] = genome_names

        progress.update("Reactions")
        compound_compartments: Dict[str, Set[str]] = {}
        for reaction_id, reaction in self.reactions.items():
            reaction_entry = JSONStructure.get_reaction_entry()
            json_reactions.append(reaction_entry)
            reaction_entry['id'] = reaction_id
            reaction_entry['name'] = reaction.modelseed_name
            metabolites = reaction_entry['metabolites']
            for compound_id, compartment, coefficient in zip(
                reaction.compound_ids, reaction.compartments, reaction.coefficients
            ):
                metabolites[f"{compound_id}_{compartment}"] = coefficient
                try:
                    compound_compartments[compound_id].add(compartment)
                except KeyError:
                    compound_compartments[compound_id] = set(compartment)
            if not reaction.reversibility:
                # By default, the reaction entry was set up to be reversible; here make it
                # irreversible.
                reaction_entry['lower_bound'] = 0.0
            reaction_entry['gene_reaction_rule'] = " or ".join(
                [gcid for gcid in reaction_gene_clusters[reaction_id]]
            )

            notes = reaction_entry['notes']
            # Record gene KO annotations which aliased the reaction via KEGG REACTION or EC number.
            notes['ko'] = ko_notes = {}
            ko_kegg_aliases = []
            ko_ec_number_aliases = []
            for ko in reaction_kos[reaction_id]:
                try:
                    kegg_aliases = ko.kegg_reaction_aliases[reaction_id]
                except KeyError:
                    kegg_aliases = []
                try:
                    ec_number_aliases = ko.ec_number_aliases[reaction_id]
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
            notes['genomes'] = sorted(set(reaction_genomes[reaction_id]))

        progress.update("Metabolites")
        for compound_id, metabolite in self.metabolites.items():
            modelseed_compound_name = metabolite.modelseed_name
            charge = metabolite.charge
            formula = metabolite.formula
            kegg_compound_aliases = list(metabolite.kegg_aliases)
            for compartment in compound_compartments[compound_id]:
                metabolite_entry = JSONStructure.get_metabolite_entry()
                json_metabolites.append(metabolite_entry)
                entry_id = f"{compound_id}_{compartment}"
                metabolite_entry['id'] = entry_id
                metabolite_entry['name'] = modelseed_compound_name
                metabolite_entry['compartment'] = compartment
                # Compounds without a formula have a nominal charge of 10000000 in the ModelSEED
                # compounds database, which is replaced by None in the reaction network and 0 in the
                # JSON.
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

class KEGGData:
    """
    This object handles KEGG reference data used in reaction networks.

    Attributes
    ==========
    kegg_context : anvio.kegg.KeggContext
        This contains anvi'o KEGG database attributes, such as filepaths.

    modules_db : anvio.kegg.ModulesDatabase
        The anvi'o modules database from which KEGG data is loaded.

    modules_db_hash : str
        The unique identifier of the anvi'o modules database, unique to its contents.

    ko_data : Dict[str, Dict[str, Any]]
        This dictionary relates KO IDs to various data, as shown in the following schematic.
        ko_data = {
            <KO 1 ID>:
                {
                    'EC':   (<EC numbers>),
                    'RN':   (<KEGG reaction IDs>),
                    'MOD':  (<KEGG module IDs>),
                    'PTH':  (<KEGG pathway IDs>),
                    'HIE':  {
                        <BRITE hierarchy ID A>: (
                            (<hierarchical category names X, Y, Z, ... classifying KO 1>),
                            (<hierarchical category names A, B, C, ... classifying KO 1>),
                            ...
                        ),
                        <BRITE hierarchy ID B>: (...),
                        ...
                    }
                },
            <KO 2>: {...},
            ...
        }

    module_data : Dict[str, Dict[str, Any]]
        This dictionary relates module IDs to module names and pathways, as shown in the following
        schematic.
        module_data = {
            <module 1 ID>:
                {
                    'NAME': <module 1 name>,
                    'PTH':  (<KEGG PATHWAY IDS>)
                },
            <module 2 ID>: {...},
            ...
        }

    pathway_data : Dict[str, Dict[str, Any]]
        This dictionary relates pathway IDs to pathway names and equivalent categories in the BRITE
        hierarchy, 'ko00001', as shown in the following schematic.
        pathway_data = {
            <pathway 1 ID>:
                {
                    'NAME': <pathway 1 name>,
                    'CAT':  (<hierarchical category names X, Y, Z equivalent to pathway 1>)
                },
            <pathway 2 ID>: {...},
            ...
        }

    hierarchy_data : Dict[str, str]
        This dictionary relates BRITE hierarchy IDs to hierarchy names.
    """
    def __init__(self, kegg_dir: str = None) -> None:
        """
        Set up the KEGG data in attributes designed for reaction network construction.

        Parameters
        ==========
        kegg_dir : str, None
            Directory containing an anvi'o KEGG database. The default argument of None expects the
            KEGG database to be set up in the default directory used by the program
            `anvi-setup-kegg-data`.
        """
        args = argparse.Namespace()
        args.kegg_data_dir = kegg_dir
        self.kegg_context = kegg.KeggContext(args)

        missing_paths = self.check_for_binary_relation_files()
        if missing_paths:
            raise ConfigError(
                "Unfortunately, the KEGG database needs to be reinstalled, because the following "
                "KEGG binary relation files used in reaction networks were not found: "
                f"{', '.join(missing_paths)}"
            )

        utils.is_kegg_modules_db(self.kegg_context.kegg_modules_db_path)

        self.modules_db = kegg.ModulesDatabase(
            self.kegg_context.kegg_modules_db_path, argparse.Namespace(quiet=True)
        )
        self.modules_db_hash = self.modules_db.db.get_meta_value('hash')

        self.ko_data: Dict[str, Dict[str, Any]] = {}
        self.module_data: Dict[str, Dict[str, Any]] = {}
        self.pathway_data: Dict[str, Dict[str, Any]] = {}
        self.hierarchy_data: Dict[str, str] = {}

        self._load_ko_binary_relation_data()
        self._load_module_data()
        self._load_ko_module_data()
        self._load_ko_pathway_data()
        self._load_ko_hierarchy_data()
        self._load_pathway_data()
        self._load_hierarchy_data()

        # Add placeholder entries for missing KO data.
        for ko_dict in self.ko_data.values():
            if 'EC' not in ko_dict:
                ko_dict['EC'] = tuple()
            if 'RN' not in ko_dict:
                ko_dict['RN'] = tuple()
            if 'MOD' not in ko_dict:
                ko_dict['MOD'] = tuple()
            if 'PTH' not in ko_dict:
                ko_dict['PTH'] = tuple()
            if 'HIE' not in ko_dict:
                ko_dict['HIE'] = {}

        self.modules_db.disconnect()

    def check_for_binary_relation_files(self) -> List[str]:
        """
        Check for expected binary relation files.

        Returns
        =======
        List[str]
            The paths of missing binary relation files not found at the expected locations.
        """
        missing_paths = []
        for file in self.kegg_context.kegg_binary_relation_files.values():
            path = os.path.join(self.kegg_context.binary_relation_data_dir, file)
            if not os.path.exists(path):
                missing_paths.append(path)

        return missing_paths

    def _load_ko_binary_relation_data(self) -> None:
        """
        Load KO binary relations to EC numbers and KEGG reactions into 'ko_data'.

        Returns
        =======
        None
        """
        ko_data = self.ko_data
        for binary_relation, file in self.kegg_context.kegg_binary_relation_files.items():
            binary_relation_path = os.path.join(self.kegg_context.binary_relation_data_dir, file)
            binary_relation_df = pd.read_csv(binary_relation_path, sep='\t', header=0).rename(
                {'#KO': 'KO'}, axis=1
            )

            entry_label = binary_relation[1]
            if entry_label == 'EC':
                col_name = 'EC number'
            elif entry_label == 'RN':
                col_name = 'Reaction'
            else:
                raise AssertionError
            binary_relation_df = binary_relation_df.rename({col_name: entry_label}, axis=1)
            col_name = entry_label
            assert binary_relation_df.columns.tolist() == ['KO', col_name]

            prefix_length = len(f'[{entry_label}:')
            for line in binary_relation_df.itertuples():
                ko_id: str = line.KO
                try:
                    ko_data_dict = ko_data[ko_id]
                except KeyError:
                    ko_data[ko_id] = ko_data_dict = {}
                entry: str = getattr(line, col_name)
                # An entry has a format like '[EC:1.1.1.18 1.1.1.369]' or '[RN:R00842]'.
                ko_data_dict[col_name] = tuple(entry[prefix_length:-1].split())

    def _load_module_data(self) -> None:
        """
        Load module name and pathway membership into 'module_data'.

        Returns
        =======
        None
        """
        modules_names_table = self.modules_db.db.get_table_as_dataframe(
            'modules',
            where_clause='data_name="NAME"',
            columns_of_interest=['module', 'data_value']
        ).rename({'module': 'module_id', 'data_value': 'module_name'}, axis=1)

        modules_pathways_table = self.modules_db.db.get_table_as_dataframe(
            'modules',
            where_clause='data_name="PATHWAY"',
            columns_of_interest=['module', 'data_value']
        ).rename({'module': 'module_id', 'data_value': 'pathway_id'}, axis=1)

        module_data = self.module_data
        for key, module_table in modules_names_table.merge(
            modules_pathways_table, how='left', on='module_id'
        ).groupby(['module_id', 'module_name']):
            module_id = key[0]
            module_name = key[1]
            module_data[module_id] = module_dict = {}
            module_dict['NAME'] = module_name
            pathway_id_tuple = tuple(sorted(module_table['pathway_id']))
            if pd.isna(pathway_id_tuple[0]):
                # The module is not part of any pathway.
                module_dict['PTH'] = tuple()
            else:
                module_dict['PTH'] = pathway_id_tuple

    def _load_ko_module_data(self) -> None:
        """
        Load KO classification within modules into 'ko_data'.

        Returns
        =======
        None
        """
        kos_modules_table = self.modules_db.db.get_table_as_dataframe(
            'modules',
            where_clause='data_name="ORTHOLOGY"',
            columns_of_interest=['module', 'data_value', 'data_definition']
        ).rename({'module': 'module_id', 'data_value': 'ko_id'}, axis=1)

        ko_id_pattern = re.compile('K\d{5}')
        for orthology_entry, ko_modules_table in kos_modules_table.groupby('ko_id'):
            if not re.match(ko_id_pattern, orthology_entry):
                # Screen for "orthology" entries that are validly formatted KO IDs. There are also
                # orthology entires for modules that are part of other modules.
                continue
            ko_id = orthology_entry

            try:
                ko_dict = self.ko_data[ko_id]
            except KeyError:
                self.ko_data[ko_id] = ko_dict = {}
            ko_dict['MOD'] = tuple(sorted(ko_modules_table['module_id']))

    def _load_ko_pathway_data(self) -> None:
        """
        Load KO classification within pathways into 'ko_data'.

        Only pathways that are equivalent to categories in the KO BRITE hierarchy, 'ko00001', are
        considered, because the complete KO membership of these pathways (maps) is accessible in the
        modules database given how BRITE hierarchy files are processed by 'anvi-setup-kegg-data'.
        This excludes global and overview metabolism maps, such as the global 'Metabolic pathways'
        and overview 'Degradation of aromatic compounds' maps. This also excludes maps corresponding
        to hierarchies with IDs starting 'br' rather than 'ko'. 'br' maps involve other data besides
        KOs, such as drugs, diseases, reactions, and compounds.

        Returns
        =======
        None
        """
        hierarchy_table = self.modules_db.db.get_table_as_dataframe(
            'brite_hierarchies',
            where_clause='hierarchy_accession="ko00001"',
            columns_of_interest=['ortholog_accession', 'categorization']
        ).rename({'hierarchy_accession': 'hierarchy_id', 'ortholog_accession': 'ko_id'}, axis=1)

        for ko_id, ko_table in hierarchy_table.groupby('ko_id'):
            try:
                ko_dict = self.ko_data[ko_id]
            except KeyError:
                self.ko_data[ko_id] = ko_dict = {}

            pathway_ids: List[str] = []
            for categorization in ko_table['categorization']:
                categorization: str
                category_name = categorization.split('>>>')[-1]
                # 'ko00001' categories that are equivalent to pathways have names formatted like
                # '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'.
                if category_name[-15:-8] != ' [PATH:' and category_name[-1] != ']':
                    continue
                pathway_ids.append(f'map{category_name[-6:-1]}')
            ko_dict['PTH'] = tuple(pathway_ids)

    def _load_ko_hierarchy_data(self) -> None:
        """
        Load KO BRITE hierarchy membership into 'ko_data'.

        Only KO hierarchies with IDs starting 'ko' are considered given modules database setup.

        Returns
        =======
        None
        """
        hierarchies_table = self.modules_db.db.get_table_as_dataframe(
            'brite_hierarchies',
            columns_of_interest=['hierarchy_accession', 'ortholog_accession', 'categorization']
        ).rename({'hierarchy_accession': 'hierarchy_id', 'ortholog_accession': 'ko_id'}, axis=1)

        for ko_id, ko_table in hierarchies_table.groupby('ko_id'):
            try:
                ko_dict = self.ko_data[ko_id]
            except KeyError:
                self.ko_data[ko_id] = ko_dict = {}

            ko_hierarchies_dict: Dict[str, Tuple[Tuple[str]]] = {}
            for hierarchy_id, ko_hierarchy_table in ko_table.groupby('hierarchy_id'):
                categorizations: List[Tuple[str]] = []
                for categorization in ko_hierarchy_table['categorization']:
                    categorization: str
                    categorizations.append(tuple(categorization.split('>>>')))
                ko_hierarchies_dict[hierarchy_id] = tuple(categorizations)
            ko_dict['HIE'] = ko_hierarchies_dict

    def _load_pathway_data(self) -> None:
        """
        Load pathway name and BRITE categorization into 'pathway_data'.

        Only pathways that are equivalent to categories in the KO BRITE hierarchy, 'ko00001', are
        considered.

        Returns
        =======
        None
        """
        categorizations = self.modules_db.db.get_single_column_from_table(
            'brite_hierarchies',
            'categorization',
            unique=True,
            where_clause='hierarchy_accession="ko00001"'
        )

        for categorization in categorizations:
            categorization: str
            categorization_tuple = tuple(categorization.split('>>>'))
            category_name = categorization_tuple[-1]
            # 'ko00001' categories that are equivalent to pathways have names formatted like '00010
            # Glycolysis / Gluconeogenesis [PATH:ko00010]'.
            if category_name[-15:-8] != ' [PATH:' and category_name[-1] != ']':
                continue
            pathway_id = f'map{category_name[-6:-1]}'
            assert category_name[:6] == f'{pathway_id[3:]} '
            pathway_dict: Dict[str, Any] = {}
            pathway_name = f'{category_name[6:-14]}'
            pathway_dict['NAME'] = pathway_name
            pathway_dict['CAT'] = categorization_tuple
            self.pathway_data[pathway_id] = pathway_dict

    def _load_hierarchy_data(self) -> None:
        """
        Load hierarchy names into 'hierarchy_data'.

        Returns
        =======
        None
        """
        hierarchies_table = self.modules_db.db.get_table_as_dataframe(
            'brite_hierarchies', columns_of_interest=['hierarchy_accession', 'hierarchy_name']
        ).rename({'hierarchy_accession': 'hierarchy_id'}, axis=1).drop_duplicates()

        for row in hierarchies_table.itertuples():
            hierarchy_id = row.hierarchy_id
            hierarchy_name = row.hierarchy_name
            self.hierarchy_data[hierarchy_id] = hierarchy_name

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
                raise ConfigError(
                    f"No required file named '{expected_file}' was found in the KO directory, "
                    f"'{ko_dir}'."
                )

        f = open(os.path.join(ko_dir, 'ko_info.txt'))
        f.readline()
        self.release = ' '.join(f.readline().strip().split()[1:])
        f.close()

        self.ko_table = pd.read_csv(
            os.path.join(ko_dir, 'ko_data.tsv'), sep='\t', header=0, index_col=0, low_memory=False
        )

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
                raise ConfigError(
                    f"There is no such directory, '{dir}'. You should create it first if you want "
                    "to use it."
                )
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
                input_queue.put(
                    (f'{download_root}get/{ko_id}', os.path.join(ko_dir, f'{ko_id}.txt'))
                )
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
                    "Unfortunately, files for the following KOs failed to download despite "
                    "multiple attempts, and so the database needs to be set up again: "
                    f"{', '.join(undownloaded)}"
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
                    "from the KO database. Anvi'o will now attempt to redownload all of the files."
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

    Attributes
    ==========
    sha : str
        The git commit SHA used to track the version of the downloaded ModelSEED database.

    compounds_table : pandas.core.frame.DataFrame
        The ModelSEED Biochemistry Compound database.

    kegg_reactions_table : pandas.core.frame.DataFrame
        The ModelSEED Biochemistry Reaction database, retaining ModelSEED reactions with KEGG reaction aliases
        and keying by KEGG reaction ID.

    ec_reactions_table : pandas.core.frame.DataFrame
        The ModelSEED Biochemistry Reaction database, retaining ModelSEED reactions with EC number
        aliases and keying by EC number.
    """
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/misc/MODELSEED')

    # Compounds are identified as cytosolic or extracellular in ModelSEED reactions.
    compartment_ids = {0: 'c', 1: 'e'}

    def __init__(self, modelseed_dir: str = None) -> None:
        """
        Load and set up reorganized tables of reactions and compounds from the ModelSEED directory
        to facilitate reaction network construction.

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
                "No required file named 'sha.txt' was found in the ModelSEED directory, "
                f"'{modelseed_dir}'."
            )
        reactions_path = os.path.join(modelseed_dir, 'reactions.tsv')
        if not os.path.isfile(reactions_path):
            raise ConfigError(
                "No required file named 'reactions.tsv' was found in the ModelSEED directory, "
                f"'{modelseed_dir}'."
            )
        compounds_path = os.path.join(modelseed_dir, 'compounds.tsv')
        if not os.path.isfile(compounds_path):
            raise ConfigError(
                "No required file named 'compounds.tsv' was found in the ModelSEED directory, "
                f"'{modelseed_dir}'."
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
                            "the maximum number of attempts. Try setting up the ModelSEED database "
                            "again."
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
        kegg_dir: str = None,
        modelseed_dir: str = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Parameters
        ==========
        kegg_dir : str, None
            The directory containing the anvi'o KEGG database. The default argument of None expects
            the KEGG database to be set up in the default directory used by the program
            `anvi-setup-kegg-data`.

        modelseed_dir : str, None
            The directory containing reference ModelSEED Biochemistry tables set up by anvi'o. The
            default argument of None expects ModelSEED data to be set up in the default anvi'o
            directory used by the program `anvi-setup-modelseed-database`.

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
        """
        self.kegg_dir = kegg_dir
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
                "A reaction network must be loaded from a database source. Either a contigs "
                "database or a genomes storage database and pan database are required."
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
        if stats_file is not None:
            filesnpaths.is_output_file_writable(stats_file)

        # Load the contigs database.
        utils.is_contigs_db(contigs_db)
        cdb = ContigsDatabase(contigs_db)
        cdb_db: DB = cdb.db
        sources: List[str] = cdb.meta['gene_function_sources']
        if not sources or not 'KOfam' in sources:
            raise ConfigError(
                "The contigs database indicates that genes were never annotated with KOs. This is "
                "especially strange since to load a reaction network means that a network had to "
                "be constructed from gene KO annotations in the database. "
            )

        # Check that the network stored in the contigs database was made from the same set of KO
        # gene annotations as is in the database.
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
                    "The reaction network stored in the contigs database was made from a different "
                    "two solutions to this problem. First, 'anvi-reaction-network' can be run "
                    "again to overwrite the existing network stored in the database with a new "
                    "network from the new KO gene annotations. Second, 'check_gene_annotations' "
                    "can be made False rather than True, allowing the stored network to have been "
                    "made from a different set of KO gene annotations than is currently stored in "
                    "the database. This can result in different genes being associated with KOs in "
                    "the returned GenomicNetwork than in the original network that was stored. The "
                    "available version of the KO database that has been set up by anvi'o is used "
                    "to fill in data for KOs in the network that are not current gene annotations."
                )
            self.run.warning(
                "The reaction network stored in the contigs database was made from a different set "
                "of KEGG KO gene annotations than is currently in the database. This will be "
                "ignored since 'check_gene_annotations' is False. This can result in different "
                "KO assignments to genes in the returned GenomicNetwork than in the original "
                "network that was stored."
            )

        network = GenomicNetwork(run=self.run, progress=self.progress)
        network.contigs_db_source_path = os.path.abspath(contigs_db)

        # Identify KOs present in the reaction network as it was stored that have disappeared from
        # among the gene-KO annotations.
        ko_id_pattern = re.compile('K\d{5}')
        reaction_network_ko_ids: Set[str] = set([
            kegg_id for kegg_id in
            set(cdb_db.get_single_column_from_table('reaction_network_kegg', 'kegg_id'))
            if re.fullmatch(ko_id_pattern, kegg_id)
        ])
        contigs_db_ko_ids = set(gene_ko_hits_table['accession'])
        missing_ko_ids = reaction_network_ko_ids.difference(contigs_db_ko_ids)
        if missing_ko_ids:
            self.run.warning(
                "The following KOs present in the reaction network as it was originally stored are "
                "not present among the current gene-KO hits in the contigs database, indicating "
                f"that genes were reannotated: {', '.join(missing_ko_ids)}"
            )

        # Count the genes with KO hits, a summary statistic.
        num_genes_assigned_kos = gene_ko_hits_table['gene_callers_id'].nunique()

        # Make objects representing genes with KO annotations in the stored reaction network. Make
        # objects representing the KOs, initially only assigning their ID attribute.
        gene_ko_hits_table: pd.DataFrame = gene_ko_hits_table.set_index('accession').loc[
            reaction_network_ko_ids.intersection(contigs_db_ko_ids)
        ]
        gene_ko_hits_table = gene_ko_hits_table.reset_index().set_index('gene_callers_id')
        for row in gene_ko_hits_table.itertuples():
            gcid = row.Index
            ko_id = row.accession
            ko_name = row.function
            e_value = float(row.e_value)

            try:
                # This is not the first annotation involving the gene, so an object for it already
                # exists.
                gene = network.genes[gcid]
            except KeyError:
                gene = Gene(gcid=gcid)
                network.genes[gcid] = gene

            try:
                # This is not the first annotation involving the KO, so an object for it already
                # exists.
                ko = network.kos[ko_id]
            except KeyError:
                ko = KO(id=ko_id, name=ko_name)
                network.kos[ko_id] = ko
            gene.ko_ids.append(ko_id)
            gene.e_values[ko_id] = e_value

        self._load_modelseed_reactions(cdb, network)
        self._load_modelseed_compounds(cdb, network)
        self._load_ko_classifications(cdb, network)

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
            'genes_assigned_kos': num_genes_assigned_kos,
            'kos_assigned_genes': len(contigs_db_ko_ids)
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
            Database storing protein measurements.

        contigs_database : ContigsDatabase
            Database storing associations between genes and proteins.

        network : GenomicNetwork
            Genomic network under construction.

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
                protein.gcids.append(gcid)
                if gene.protein_id:
                    try:
                        multiprotein_genes[gcid].append(protein_id)
                    except KeyError:
                        multiprotein_genes[gcid] = [protein_id]
                else:
                    gene.protein_id = protein_id
            for row in protein_table.itertuples():
                protein.abundances[row.sample_name] = row.abundance_value

        if multiprotein_genes:
            msg = ""
            for gcid, protein_ids in multiprotein_genes.items():
                msg += f"{gcid}: {', '.join(protein_ids)}; "
            msg = msg[: -1]
            raise ConfigError(
                "Certain genes were unexpectedly associated with multiple proteins with abundance "
                "data. Unfortunately, multiple protein products are not currently allowed in "
                "anvi'o, so the protein abundance data must be edited down in the profile database "
                "to permit use with the reaction network. These are as follows, with the gene "
                f"callers ID separated by a comma-separated\ list of protein IDs. {msg}"
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
            Database storing protein measurement data that is loaded into the genomic network.

        network : GenomicNetwork
            Genomic network under construction.

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

        for compound_id, metabolite_table in metabolite_abundances_table.groupby('reference_id'):
            metabolite = network.metabolites[compound_id]
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
            Path to a pan database in which a reaction network is stored.

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
            Reaction network loaded from the pangenomic databases.
        """
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
                    "The reaction network stored in the pan database was made from a different set "
                    "of KO gene annotations than is currently in the associated genomes storage "
                    "database. There are two solutions to this problem. First, the program, "
                    "'anvi-reaction-network', can be run again to overwrite the existing network "
                    "stored in the pan database with a new network from the new KO gene "
                    "annotations. Second, 'check_gene_annotations' can be given an argument of "
                    "False instead of True, preventing this exception from being raised if the "
                    "stored network was made from a different set of KO gene annotations than is "
                    "currently in the genomes storage database. This can result in different "
                    "consensus KOs assigned to gene clusters in the returned PangenomicNetwork "
                    "than in the original network that was stored. The available version of the KO "
                    "database that has been set up by anvi'o is used to fill in data for any KOs "
                    "in the network that are not current gene annotations in the genomes storage "
                    "database."
                )
            self.run.warning(
                "The reaction network stored in the pan database was made from a different set of "
                "KO gene annotations than is currently in the genomes storage database. This will "
                "be ignored since 'check_gene_annotations' is False. This can result in different "
                "consensus KO assignments to gene clusters in the returned PangenomicNetwork than "
                "in the original network that was stored."
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

        # Find gene clusters with consensus KO annotations. Make objects representing gene clusters
        # with KO annotations in the stored reaction network. Make objects representing the KOs,
        # initially only assigning their ID attribute.

        pdb = PanDatabase(pan_db)
        pdb_db: DB = pdb.db
        ko_id_pattern = re.compile('K\d{5}')
        reaction_network_ko_ids: List[str] = [
            kegg_id for kegg_id in
            set(pdb_db.get_single_column_from_table('reaction_network_kegg', 'kegg_id'))
            if re.fullmatch(ko_id_pattern, kegg_id)
        ]

        num_gene_clusters_assigned_ko = 0
        ko_ids_assigned_gene_cluster = []
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
            ko_ids_assigned_gene_cluster.append(ko_id)

            if ko_id not in reaction_network_ko_ids:
                # The KO is not in the stored reaction network, indicating that it is a newer
                # annotation.
                continue

            gene_cluster = GeneCluster(gene_cluster_id=cluster_id, ko_id=ko_id)
            gene_cluster.genomes = list(pan_super.gene_clusters[cluster_id])
            network.gene_clusters[cluster_id] = gene_cluster

            try:
                # This is not the first gene cluster that has been encountered with the KO assigned
                # to it, so an object for the KO already exists.
                ko = network.kos[ko_id]
            except KeyError:
                ko = KO(id=ko_id, name=gene_cluster_ko_data['function'])
                network.kos[ko_id] = ko
        ko_ids_assigned_gene_cluster = set(ko_ids_assigned_gene_cluster)

        missing_ko_ids = set(reaction_network_ko_ids).difference(set(network.kos))
        if missing_ko_ids:
            self.run.warning(
                "The following KOs present in the reaction network as it was originally stored are "
                "not present among the current gene cluster consensus KOs derived from the "
                "pangenomic databases, suggesting that genes underlying the gene clusters were "
                "reannotated."
            )

        self._load_modelseed_reactions(pdb, network)
        self._load_modelseed_compounds(pdb, network)
        self._load_ko_classifications(pdb, network)

        if quiet and not stats_file:
            return network

        if network.consistent_annotations:
            precomputed_counts = {
                'total_gene_clusters': pdb.meta['num_gene_clusters'],
                'gene_clusters_assigned_ko': num_gene_clusters_assigned_ko,
                'kos_assigned_gene_clusters': len(ko_ids_assigned_gene_cluster)
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
        network: ReactionNetwork
    ) -> None:
        """
        Add ModelSEED reactions to the network being loaded from the database.

        ModelSEED reaction objects are related to KOs through KEGG REACTION and EC number aliases.

        Parameters
        ==========
        database : ContigsDatabase or PanDatabase
            Database storing a reaction network.

        network : ReactionNetwork
            Network under construction. A GenomicNetwork is loaded from a ContigsDatabase; a
            PangenomicNetwork is loaded from a PanDatabase.

        Returns
        =======
        None
        """
        # Load the table of reactions data.
        if type(database) is ContigsDatabase:
            reactions_table = database.db.get_table_as_dataframe(
                tables.reaction_network_reactions_table_name
            )
            if type(network) is not GenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'ContigsDatabase', so the provided "
                    "'network' must be of type 'GenomicNetwork'. Instead, the reaction network "
                    f"argument was of type '{type(network)}'."
                )
        elif type(database) is PanDatabase:
            reactions_table = database.db.get_table_as_dataframe(
                tables.pan_reaction_network_reactions_table_name
            )
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

        # Each row of the table contains information on a different ModelSEED reaction.
        for row in reactions_table.itertuples():
            # Check that the reaction is associated with a KO that matches a gene in the contigs
            # database (or is assigned to a gene cluster in the pan database). If KO gene
            # annotations have been updated from those used to create the stored reaction network,
            # then certain KOs and everything inferred from the KOs may no longer annotate genes
            # and gene clusters, and are therefore not loaded.
            is_networked = False

            # Map KEGG reaction aliases of the ModelSEED reaction to all KOs that were associated
            # with the KEGG reaction.
            kegg_reaction_sources: str = row.ko_kegg_reaction_source
            kegg_reaction_kos: Dict[str, List[KO]] = {}
            for kegg_reaction_item in kegg_reaction_sources.split('; '):
                if not kegg_reaction_item:
                    # The ModelSEED reaction was not sourced from KEGG reactions.
                    continue
                kegg_reaction_id, ko_ids = kegg_reaction_item.split(': (')
                ko_ids = ko_ids[:-1].split(', ')
                kegg_reaction_kos[kegg_reaction_id] = kos = []
                for ko_id in ko_ids:
                    try:
                        kos.append(network.kos[ko_id])
                    except KeyError:
                        continue
                if kos:
                    is_networked = True

            # Map EC number aliases of the ModelSEED reaction to all KOs that were associated with
            # the EC number.
            ec_number_sources: str = row.ko_ec_number_source
            ec_number_kos: Dict[str, List[KO]] = {}
            for ec_number_item in ec_number_sources.split('; '):
                if not ec_number_item:
                    # The ModelSEED reaction was not sourced from EC numbers.
                    continue
                ec_number, ko_ids = ec_number_item.split(': (')
                ko_ids = ko_ids[:-1].split(', ')
                ec_number_kos[ec_number] = kos = []
                for ko_id in ko_ids:
                    try:
                        kos.append(network.kos[ko_id])
                    except KeyError:
                        continue
                if kos:
                    is_networked = True

            if not is_networked:
                continue

            reaction = ModelSEEDReaction()
            reaction_id: str = row.modelseed_reaction_id
            reaction.modelseed_id = reaction_id
            reaction.modelseed_name = row.modelseed_reaction_name
            network.reactions[reaction_id] = reaction

            modelseed_compound_ids: str = row.metabolite_modelseed_ids
            reaction_compound_ids = []
            for compound_id in modelseed_compound_ids.split(', '):
                reaction_compound_ids.append(compound_id)
                if compound_id not in network.metabolites:
                    metabolite = ModelSEEDCompound()
                    metabolite.modelseed_id = compound_id
                    network.metabolites[compound_id] = metabolite
            reaction.compound_ids = tuple(reaction_compound_ids)

            stoichiometry: str = row.stoichiometry
            reaction.coefficients = tuple(int(coeff) for coeff in stoichiometry.split(', '))
            compartments: str = row.compartments
            reaction.compartments = tuple(compartments.split(', '))
            reversibility: int = row.reversibility
            reaction.reversibility = bool(reversibility)

            # Record *all* KEGG reaction aliases of the ModelSEED reaction, including those not
            # associated with KO annotations.
            other_kegg_reaction_ids: str = row.other_kegg_reaction_ids
            reaction.kegg_aliases = list(kegg_reaction_kos)
            if other_kegg_reaction_ids:
                reaction.kegg_aliases += other_kegg_reaction_ids.split(', ')
            reaction.kegg_aliases = tuple(reaction.kegg_aliases)

            network.modelseed_kegg_aliases[reaction_id] = modelseed_kegg_aliases = []
            for kegg_reaction_id, kos in kegg_reaction_kos.items():
                # Record the ModelSEED reaction as one of the aliases of the KEGG reaction in the
                # network.
                try:
                    network.kegg_modelseed_aliases[kegg_reaction_id].append(reaction_id)
                except KeyError:
                    network.kegg_modelseed_aliases[kegg_reaction_id] = [reaction_id]
                modelseed_kegg_aliases.append(kegg_reaction_id)
                for ko in kos:
                    if not reaction_id in ko.reaction_ids:
                        # This is the first time encountering the reaction as a reference of the KO.
                        ko.reaction_ids.append(reaction_id)
                    try:
                        ko.kegg_reaction_aliases[reaction_id].append(kegg_reaction_id)
                    except KeyError:
                        ko.kegg_reaction_aliases[reaction_id] = [kegg_reaction_id]

            # Record *all* EC number aliases of the ModelSEED reaction, including those not
            # associated with KO annotations.
            other_ec_numbers: str = row.other_ec_numbers
            reaction.ec_number_aliases = list(ec_number_kos)
            if other_ec_numbers:
                reaction.ec_number_aliases += other_ec_numbers.split(', ')
            reaction.ec_number_aliases = tuple(reaction.ec_number_aliases)

            modelseed_ec_number_aliases = []
            network.modelseed_ec_number_aliases[reaction_id] = modelseed_ec_number_aliases
            for ec_number, kos in ec_number_kos.items():
                # Record the ModelSEED reaction as one of the aliases of the EC number in the
                # network.
                try:
                    network.ec_number_modelseed_aliases[ec_number].append(reaction_id)
                except KeyError:
                    network.ec_number_modelseed_aliases[ec_number] = [reaction_id]
                modelseed_ec_number_aliases.append(ec_number)
                for ko in kos:
                    if not reaction_id in ko.reaction_ids:
                        # This is the first time encountering the reaction as a reference of the KO.
                        ko.reaction_ids.append(reaction_id)
                    try:
                        ko.ec_number_aliases[reaction_id].append(ec_number)
                    except KeyError:
                        ko.ec_number_aliases[reaction_id] = [ec_number]

    def _load_modelseed_compounds(
        self,
        database: Union[ContigsDatabase, PanDatabase],
        network: ReactionNetwork
    ) -> None:
        """
        Add ModelSEED compounds to the network being loaded from the database.

        Parameters
        ==========
        database : ContigsDatabase or PanDatabase
            Database storing a reaction network.

        network : GenomicNetwork or PangenomicNetwork
            Network under construction. A GenomicNetwork is loaded from a ContigsDatabase; a
            PangenomicNetwork is loaded from a PanDatabase.

        Returns
        =======
        None
        """
        # Load the table of compounds data.
        if type(database) is ContigsDatabase:
            metabolites_table = database.db.get_table_as_dataframe(
                tables.reaction_network_metabolites_table_name
            )
            if type(network) is not GenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'ContigsDatabase', so the provided "
                    "'network' must be of type 'GenomicNetwork'. Instead, the reaction network "
                    f"argument was of type '{type(network)}'."
                )
        elif type(database) is PanDatabase:
            metabolites_table = database.db.get_table_as_dataframe(
                tables.pan_reaction_network_metabolites_table_name
            )
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

        # Each row of the table contains information on a different ModelSEED compound.
        for row in metabolites_table.itertuples():
            compound_id = row.modelseed_compound_id
            try:
                metabolite = network.metabolites[compound_id]
            except KeyError:
                # The metabolite in the stored network is not loaded. The metabolite only
                # participates in reactions that also were not loaded. The reactions are only
                # associated with KOs that no longer match genes in the contigs database (or are no
                # longer assigned to gene clusters in the pan database). Gene KO annotations must
                # have been updated from those used to create the stored reaction network.
                continue
            modelseed_compound_name: str = row.modelseed_compound_name
            metabolite.modelseed_name = modelseed_compound_name
            kegg_aliases: str = row.kegg_aliases
            if kegg_aliases:
                metabolite.kegg_aliases = tuple(kegg_aliases.split(', '))
            else:
                metabolite.kegg_aliases = tuple()
            # Compounds without a formula, recorded here as None, have a nominal charge of 10000000
            # in the ModelSEED compounds database. This is replaced by NaN in the table and here as
            # None in the reaction network.
            formula: str = row.formula
            metabolite.formula = formula
            charge: int = row.charge
            metabolite.charge = None if np.isnan(charge) else int(charge)


    def _load_ko_classifications(
        self,
        database: Union[ContigsDatabase, PanDatabase],
        network: ReactionNetwork
    ) -> None:
        """
        Add information on KEGG module, pathway, and BRITE hierarchy membership to the network being
        loaded from the database.

        Parameters
        ==========
        database : ContigsDatabase or PanDatabase
            Database storing a reaction network.

        network : ReactionNetwork
            Network under construction. A GenomicNetwork is loaded from a ContigsDatabase; a
            PangenomicNetwork is loaded from a PanDatabase.

        Returns
        =======
        None
        """
        # Load the table of compounds data.
        if type(database) is ContigsDatabase:
            kegg_table = database.db.get_table_as_dataframe(tables.reaction_network_kegg_table_name)
            if type(network) is not GenomicNetwork:
                raise ConfigError(
                    "The provided 'database' was of type 'ContigsDatabase', so the provided "
                    "'network' must be of type 'GenomicNetwork'. Instead, the reaction network "
                    f"argument was of type '{type(network)}'."
                )
        elif type(database) is PanDatabase:
            kegg_table = database.db.get_table_as_dataframe(
                tables.pan_reaction_network_kegg_table_name
            )
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

        # Create a separate table for each type of KEGG information.
        kegg_table = kegg_table.fillna('')
        ko_id_pattern = re.compile('K\d{5}')
        kos_table = kegg_table[kegg_table['kegg_id'].apply(
            lambda ko_id: True if re.fullmatch(ko_id_pattern, ko_id) else False
        )]
        module_id_pattern = re.compile('M\d{5}')
        modules_table = kegg_table[kegg_table['kegg_id'].apply(
            lambda module_id: True if re.fullmatch(module_id_pattern, module_id) else False
        )]
        pathway_id_pattern = re.compile('map\d{5}')
        pathways_table = kegg_table[kegg_table['kegg_id'].apply(
            lambda pathway_id: True if re.fullmatch(pathway_id_pattern, pathway_id) else False
        )]
        hierarchy_id_pattern = re.compile('ko\d{5}')
        hierarchies_table = kegg_table[kegg_table['kegg_id'].apply(
            lambda hierarchy_id: True if re.fullmatch(hierarchy_id_pattern, hierarchy_id) else False
        )]

        # Remove entries in the KO table for stored network KOs that no longer annotate genes in the
        # contigs database. (For pangenomics, instead read this and following comments as indicating
        # that the KOs no longer are consensus annotations of gene clusters from the pan database.)
        reaction_network_ko_ids: Set[str] = set(kos_table['kegg_id'])
        kos_table: pd.DataFrame = kos_table.set_index('kegg_id').loc[
            set(network.kos).intersection(reaction_network_ko_ids)
        ]

        # Fill out KEGG classification attributes of KO objects in the loaded network.
        # Record KOs that have names that differ between (newer) gene annotations in the database
        # and (older) information in the stored network due to a KEGG database update.
        inconsistent_kos: Dict[str, Tuple[str, str]] = {}
        for row in kos_table.itertuples():
            ko_id: str = row.Index
            ko = network.kos[ko_id]

            ko_name: str = row.name
            if ko.name != ko_name:
                inconsistent_kos[ko_id] = (ko.name, ko_name)

            module_ids_str: str = row.modules
            if module_ids_str:
                module_ids = module_ids_str.split(', ')
            else:
                # The KO is not classified in a module.
                module_ids = []
            for module_id in module_ids:
                ko.module_ids.append(module_id)
                try:
                    # Another KO was loaded that is classified in the module.
                    module = network.modules[module_id]
                except KeyError:
                    # Create an object in the network for the newly encountered module.
                    network.modules[module_id] = module = KEGGModule(id=module_id)
                module.ko_ids.append(ko_id)

            pathway_ids_str: str = row.pathways
            if pathway_ids_str:
                pathway_ids = pathway_ids_str.split(', ')
            else:
                # The KO is not classified in a pathway.
                pathway_ids = []
            for pathway_id in pathway_ids:
                ko.pathway_ids.append(pathway_id)
                try:
                    # Another KO was loaded that is classified in the pathway.
                    pathway = network.pathways[pathway_id]
                except KeyError:
                    # Create an object in the network for the newly encountered pathway.
                    network.pathways[pathway_id] = pathway = KEGGPathway(id=pathway_id)
                pathway.ko_ids.append(ko_id)

            categorizations_str: str = row.brite_categorization
            if categorizations_str:
                categorization_strs = categorizations_str.split(' !!! ')
            else:
                # The KO is not classified in a BRITE category.
                categorization_strs = []
            loaded_categories: Dict[str, List[Tuple[str]]] = {}
            if categorization_strs:
                for categorization_str in categorization_strs:
                    full_categorization = categorization_str.split(' >>> ')
                    hierarchy_id = full_categorization[0]
                    assert re.fullmatch(hierarchy_id_pattern, hierarchy_id)
                    try:
                        categorizations = loaded_categories[hierarchy_id]
                    except KeyError:
                        loaded_categories[hierarchy_id] = categorizations = []
                    categorization = tuple(full_categorization[1:])
                    categorizations.append(categorization)

            for hierarchy_id, categorizations in loaded_categories.items():
                ko.hierarchies[hierarchy_id] = categorizations

                try:
                    # Another KO was loaded that is classified in the hierarchy.
                    hierarchy = network.hierarchies[hierarchy_id]
                except KeyError:
                    # Create an object in the network for the newly encountered hierarchy.
                    network.hierarchies[hierarchy_id] = hierarchy = BRITEHierarchy(id=hierarchy_id)
                hierarchy.ko_ids.append(ko_id)

                try:
                    network_categorizations = network.categories[hierarchy_id]
                except KeyError:
                    network.categories[hierarchy_id] = network_categorizations = {}

                for categorization in categorizations:
                    if categorization in hierarchy.categorizations:
                        # Another KO was loaded that is classified in the category.
                        categories = network_categorizations[categorization]
                        for category in categories:
                            if ko_id not in category.ko_ids:
                                category.ko_ids.append(ko_id)
                        continue

                    hierarchy.categorizations.append(categorization)
                    categories: List[BRITECategory] = []

                    # Add the category and unencountered supercategories to the network.
                    for depth, category_name in enumerate(categorization, 1):
                        focus_categorization = categorization[:depth]
                        try:
                            # Another KO was loaded that is classified in the supercategory.
                            category = network_categorizations[focus_categorization][-1]
                            is_added = True
                        except KeyError:
                            is_added = False

                        if is_added:
                            if ko_id not in category.ko_ids:
                                category.ko_ids.append(ko_id)
                            categories.append(category)
                            continue

                        if depth > 1:
                            # The unencountered category is a subcategory of its supercategory.
                            categories[-1].subcategory_names.append(category_name)

                        category = BRITECategory()
                        category.id = f'{hierarchy_id}: {" >>> ".join(focus_categorization)}'
                        category.name = category_name
                        category.hierarchy_id = hierarchy_id
                        category.ko_ids.append(ko_id)
                        categories.append(category)
                        network_categorizations[focus_categorization] = tuple(categories)

        if inconsistent_kos:
            msg = ''
            for ko_id, ko_names in inconsistent_kos.items():
                gene_ko_name, network_ko_name = ko_names
                msg += f"{ko_id}: '{gene_ko_name}' ||| '{network_ko_name}', "
            self.run.warning(
                "KO names differ between certain records in the stored reaction network and the "
                "current gene-KO hits, indicating that genes were reannotated with KOs from a "
                "newer version of the KEGG database in which certain names have changed. Hopefully "
                "the names are similar enough to indicate that the underlying KO ID that unites "
                "them is from the same ortholog. The first name following the KO ID is from the "
                "gene annotation and the second after the '|||' separator is from the stored "
                f"network: {msg}"
            )

        # Fill out module and pathway attributes.
        for row in modules_table.itertuples():
            module_id: str = row.kegg_id
            try:
                module = network.modules[module_id]
            except KeyError:
                # The module is not loaded because it no longer contains any KOs that annotate genes
                # in the database.
                continue

            module_name: str = row.name
            module.name = module_name

            pathway_ids: str = row.pathways
            if not pathway_ids:
                continue
            for pathway_id in pathway_ids.split(', '):
                module.pathway_ids.append(pathway_id)
                pathway = network.pathways[pathway_id]
                pathway.module_ids.append(module_id)

        # Fill out pathway and equivalent BRITE category attributes.
        for row in pathways_table.itertuples():
            pathway_id: str = row.kegg_id
            try:
                pathway = network.pathways[pathway_id]
            except KeyError:
                # The pathway is not loaded because it no longer contains any KOs that annotate
                # genes in the database.
                continue

            pathway_name: str = row.name
            pathway.name = pathway_name

            categorization_str: str = row.brite_categorization
            full_categorization = categorization_str.split(' >>> ')
            hierarchy_id = full_categorization[0]
            assert hierarchy_id == 'ko00001'
            categorization = tuple(full_categorization[1:])
            pathway.categorization = categorization
            category = network.categories[hierarchy_id][categorization][-1]
            category.pathway_id = pathway_id

        # Fill out hierarchy names.
        for row in hierarchies_table.itertuples():
            hierarchy_id: str = row.kegg_id
            try:
                hierarchy = network.hierarchies[hierarchy_id]
            except KeyError:
                # The hierarchy is not loaded because it no longer contains any KOs that annotate
                # genes in the database.
                continue

            hierarchy_name: str = row.name
            hierarchy.name = hierarchy_name

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
                "Either a contigs database OR both a pan database and genomes storage database are "
                "required to make either a (meta)genomic reaction network or a pangenomic reaction "
                "network, respectively."
            )
        elif contigs_db:
            self.run.info_single(
                "A reaction network will be made from protein orthology annotations in the contigs "
                "database."
            )
            network = self.make_contigs_database_network(
                contigs_db,
                store=store,
                overwrite_existing_network=overwrite_existing_network,
                stats_file=stats_file
            )
        elif genomes_storage_db or pan_db:
            self.run.info_single(
                "A pangenomic reaction network will be made from protein orthology annotations in "
                "the genomes storage database and gene clusters in the pan database."
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
                "A reaction network cannot be made without a database source. Either a contigs "
                "database OR a pan database and genomes storage database are required to make "
                "either a (meta)genomic reaction network or a pangenomic reaction network, "
                "respectively."
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
                "The contigs database indicates that genes were never annotated with KOs, which is "
                "required to build a reaction network. This can be solved by running "
                "'anvi-run-kegg-kofams' on the contigs database."
            )
        if (
            store and
            cdb_db.get_meta_value('reaction_network_ko_annotations_hash') and
            not overwrite_existing_network
        ):
            raise ConfigError(
                "The existing reaction network in the contigs database must be explicitly "
                "overwritten."
            )

        self.progress.new("Building reaction network")

        network = GenomicNetwork(run=self.run, progress=self.progress)
        network.contigs_db_source_path = os.path.abspath(contigs_db)

        # Load reference databases.
        self.progress.update("Loading KEGG reference database")
        kegg_db = KEGGData(kegg_dir=self.kegg_dir)
        kegg_kos_data = kegg_db.ko_data
        kegg_modules_data = kegg_db.module_data
        kegg_pathways_data = kegg_db.pathway_data
        kegg_hierarchies_data = kegg_db.hierarchy_data

        self.progress.update("Loading ModelSEED Biochemistry reference database")
        modelseed_db = ModelSEEDDatabase(modelseed_dir=self.modelseed_dir)
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
                gene.ko_ids.append(ko_id)
                gene.e_values[ko_id] = e_value
                if is_new_gene:
                    # Add the unadded gene to the network.
                    network.genes[gcid] = gene
                # Proceed to the next gene-KO match.
                continue

            # Get KEGG REACTION IDs and EC numbers associated with the KO.
            try:
                ko_info = kegg_kos_data[ko_id]
            except KeyError:
                # For some reason the KO annotated a gene in the contigs database but is not found
                # in the KEGG database set up by anvi'o. Do not add the alien KO or the gene, if
                # unadded, to the network.
                undefined_ko_ids.append(ko_id)
                continue
            ko_kegg_reaction_ids = self._get_ko_kegg_reaction_ids(ko_info)
            ko_ec_numbers = self._get_ko_ec_numbers(ko_info)

            if not ko_kegg_reaction_ids and not ko_ec_numbers:
                # The KO is not associated with any KEGG reactions or EC numbers, and therefore
                # anvi'o can't relate the KO to ModelSEED reactions. Do not add the unsystematizable
                # KO or the gene, if unadded, to the network.
                continue

            # Check if KEGG reactions and EC numbers associated with the KO have already been added
            # to the network in processing other gene-KO matches. To have been added to the network,
            # KEGG reactions and EC numbers must have aliased ModelSEED reactions.
            old_kegg_reaction_ids, new_kegg_reaction_ids = self._find_kegg_reaction_ids(
                ko_kegg_reaction_ids, network
            )
            old_ec_numbers, new_ec_numbers = self._find_ec_numbers(ko_ec_numbers, network)

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
                # Add an object representing the unadded gene to the network.
                network.genes[gcid] = gene

            # Add an object representing the newly encountered KO to the network.
            ko = KO()
            ko.id = ko_id
            ko.name = ko_name
            network.kos[ko_id] = ko
            gene.ko_ids.append(ko_id)
            gene.e_values[ko_id] = e_value

            # Associate ModelSEED reactions that have previously been added to the network under
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

            # Add KEGG classifications of the KO (modules, pathways, and BRITE hierarchies) to the
            # network.
            self._add_ko_classification(
                ko,
                network,
                ko_info,
                kegg_modules_data,
                kegg_pathways_data,
                kegg_hierarchies_data
            )

        self._relate_modules_pathways(network, kegg_modules_data)

        if DEBUG:
            for gene in network.genes.values():
                for ko_id in gene.ko_ids:
                    ko = network.kos[ko_id]
                    assert ko.reaction_ids

            assert sorted(network.modelseed_kegg_aliases) == sorted(network.reactions)
            assert sorted(network.modelseed_ec_number_aliases) == sorted(network.reactions)

        undefined_modelseed_reaction_ids = set(undefined_modelseed_reaction_ids)
        undefined_ko_ids = set(undefined_ko_ids)

        self.progress.end()

        if DEBUG:
            self.run.info_single(
                "The following ModelSEED reactions would have been added to the reaction network "
                "had there been a chemical equation in the ModelSEED database; perhaps it is worth "
                "investigating the ModelSEED reactions table to understand why this is not the "
                f"case: {', '.join(undefined_modelseed_reaction_ids)}"
            )

        if undefined_ko_ids:
            self.run.info_single(
                "Certain genes matched KOs that were not found in the KEGG reference database. "
                "These KOs will not be used in network construction. It could be that the KOfams "
                "used to annotate genes were not from the same KEGG database version as the "
                "reference. Here are the unrecognized KO IDs from the contigs database: "
                f"{','.join(undefined_ko_ids)}"
            )

        kegg_dir = kegg_db.kegg_context.kegg_data_dir
        if self.modelseed_dir is None:
            modelseed_dir = ModelSEEDDatabase.default_dir
        else:
            modelseed_dir = self.modelseed_dir
        self.run.info("Reference KEGG database directory", kegg_dir, nl_before=1)
        self.run.info("Reference ModelSEED database directory", modelseed_dir, nl_after=1)

        if store:
            if cdb_db.get_meta_value('reaction_network_ko_annotations_hash'):
                self.run.warning("Deleting existing reaction network from contigs database")
                cdb_db._exec(f'''DELETE from {tables.reaction_network_reactions_table_name}''')
                cdb_db._exec(f'''DELETE from {tables.reaction_network_metabolites_table_name}''')
                cdb_db._exec(f'''DELETE from {tables.reaction_network_kegg_table_name}''')
                self.run.info_single(
                    "Deleted data in gene function reactions and metabolites tables", nl_after=1
                )

            self.progress.new("Saving reaction network to contigs database")
            self.progress.update("Reactions table")
            reactions_table = self._get_database_reactions_table(network)
            sql_statement = (
                f"INSERT INTO {tables.reaction_network_reactions_table_name} VALUES "
                f"({','.join('?' * len(tables.reaction_network_reactions_table_structure))})"
            )
            cdb_db._exec_many(sql_statement, reactions_table.values)
            self.progress.update("Metabolites table")
            metabolites_table = self._get_database_metabolites_table(network)
            sql_statement = (
                f"INSERT INTO {tables.reaction_network_metabolites_table_name} VALUES "
                f"({','.join('?' * len(tables.reaction_network_metabolites_table_structure))})"
            )
            cdb_db._exec_many(sql_statement, metabolites_table.values)
            self.progress.update("KEGG KO information table")
            kegg_table = self._get_database_kegg_table(network)
            sql_statement = (
                f"INSERT INTO {tables.reaction_network_kegg_table_name} VALUES "
                f"({','.join('?' * len(tables.reaction_network_kegg_table_structure))})"
            )
            cdb_db._exec_many(sql_statement, kegg_table.values)

            self.progress.update("Metadata")
            ko_annotations_hash = self.hash_contigs_db_ko_hits(gene_ko_hits_table)
            cdb_db.set_meta_value('reaction_network_ko_annotations_hash', ko_annotations_hash)
            # The KEGG database release is now not explicitly that, but instead the hash of the
            # anvi'o modules database.
            cdb_db.set_meta_value('reaction_network_kegg_database_release', kegg_db.modules_db_hash)
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
                "The genomes of the pangenome were not annotated with KOs, which can be rectified "
                "by running `anvi-run-kegg-kofams` on the genome contigs databases and remaking "
                "the pangenome."
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
        kegg_db = KEGGData(kegg_dir=self.kegg_dir)
        kegg_kos_data = kegg_db.ko_data
        kegg_modules_data = kegg_db.module_data
        kegg_pathways_data = kegg_db.pathway_data
        kegg_hierarchies_data = kegg_db.hierarchy_data

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
                gene_cluster.ko_id = ko_id
                # Add the newly encountered gene cluster to the network.
                network.gene_clusters[cluster_id] = gene_cluster
                # Proceed to the next gene cluster.
                continue

            # Get KEGG REACTION IDs and EC numbers associated with the KO.
            try:
                ko_info = kegg_kos_data[ko_id]
            except KeyError:
                # For some reason the KO annotated genes in the pangenomic databases but is not
                # found in the KEGG database set up by anvi'o. Do not add the alien KO or the gene
                # cluster to the network.
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
            old_ec_numbers, new_ec_numbers = self._find_ec_numbers(ko_ec_numbers, network)

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
            gene_cluster.ko_id = ko_id

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

            # Add KEGG classifications of the KO (modules, pathways, and BRITE hierarchies) to the
            # network.
            self._add_ko_classification(
                ko,
                network,
                ko_info,
                kegg_modules_data,
                kegg_pathways_data,
                kegg_hierarchies_data
            )

        self._relate_modules_pathways(network, kegg_modules_data)

        if DEBUG:
            for gene_cluster in network.gene_clusters.values():
                ko = network.kos[gene_cluster.ko_id]
                assert ko.reaction_ids

            assert sorted(network.modelseed_kegg_aliases) == sorted(network.reactions)
            assert sorted(network.modelseed_ec_number_aliases) == sorted(network.reactions)

        undefined_modelseed_reaction_ids = set(undefined_modelseed_reaction_ids)
        undefined_ko_ids = set(undefined_ko_ids)

        self.progress.end()

        if DEBUG:
            self.run.info_single(
                "The following ModelSEED reactions would have been added to the reaction network "
                "had there been a chemical equation in the ModelSEED database; perhaps it is worth "
                "investigating the ModelSEED reactions table to understand why this is not the "
                f"case: {', '.join(undefined_modelseed_reaction_ids)}"
            )

        if undefined_ko_ids:
            self.run.info_single(
                "Certain gene clusters were assigned KOs that were not found in the KEGG reference "
                "database. These KOs will not be used in network construction. It could be that "
                "the KOfams used to annotate genes were not from the same KEGG database version as "
                "the reference. Here are the unrecognized KO IDs from the pangenomic databases: "
                f"{','.join(undefined_ko_ids)}"
            )

        ko_dir = kegg_db.kegg_context.kegg_data_dir
        if self.modelseed_dir is None:
            modelseed_dir = ModelSEEDDatabase.default_dir
        else:
            modelseed_dir = self.modelseed_dir
        self.run.info("Reference KEGG database directory", ko_dir, nl_before=1)
        self.run.info("Reference ModelSEED database directory", modelseed_dir, nl_after=1)

        pdb = PanDatabase(pan_db)
        if store:
            if pan_super.p_meta['reaction_network_ko_annotations_hash']:
                self.run.warning("Deleting existing reaction network from pan database")
                pdb.db._exec(
                    f'''DELETE from {tables.pan_reaction_network_reactions_table_name}'''
                )
                pdb.db._exec(
                    f'''DELETE from {tables.pan_reaction_network_metabolites_table_name}'''
                )
                self.run.info_single(
                    "Deleted data in gene cluster function reactions and metabolites tables",
                    nl_after=1
                )

            self.progress.new("Saving reaction network to pan database")
            self.progress.update("Reactions table")
            reactions_table = self._get_database_reactions_table(network)
            table_name = tables.pan_reaction_network_reactions_table_name
            table_structure = tables.pan_reaction_network_reactions_table_structure
            pdb.db._exec_many(
                f'''INSERT INTO {table_name} VALUES ({','.join('?' * len(table_structure))})''',
                reactions_table.values
            )
            self.progress.update("Metabolites table")
            metabolites_table = self._get_database_metabolites_table(network)
            table_name = tables.pan_reaction_network_metabolites_table_name
            table_structure = tables.pan_reaction_network_metabolites_table_structure
            pdb.db._exec_many(
                f'''INSERT INTO {table_name} VALUES ({','.join('?' * len(table_structure))})''',
                metabolites_table.values
            )
            self.progress.update("KEGG KO information table")
            kegg_table = self._get_database_kegg_table(network)
            table_name = tables.pan_reaction_network_kegg_table_name
            table_structure = tables.pan_reaction_network_kegg_table_structure
            pdb.db._exec_many(
                f'''INSERT INTO {table_name} VALUES ({','.join('?' * len(table_structure))})''',
                kegg_table.values
            )

            self.progress.update("Metadata")
            ko_annotations_hash = self.hash_pan_db_ko_annotations(
                genomes_storage_db,
                gene_clusters_functions_summary_dict,
                consensus_threshold=consensus_threshold,
                discard_ties=discard_ties
            )
            pdb.db.set_meta_value('reaction_network_ko_annotations_hash', ko_annotations_hash)
            # The KEGG database release is now not explicitly that, but instead the hash of the
            # anvi'o modules database.
            pdb.db.set_meta_value('reaction_network_kegg_database_release', kegg_db.modules_db_hash)
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

    def _get_ko_kegg_reaction_ids(self, ko_info: Dict[str, Any]) -> Set[str]:
        """
        Get the set of KEGG REACTION IDs associated with the KO under consideration in network
        construction.

        Parameters
        ==========
        ko_info : Dict[str, Any]
            Information on the KO loaded from the anvi'o KEGG database.

        Returns
        =======
        Set[str]
            KEGG REACTION IDs associated with the KO.
        """
        try:
            ko_kegg_reaction_ids = set(ko_info['RN'])
        except KeyError:
            # The KO is not associated with KEGG reactions.
            ko_kegg_reaction_ids = set()
        ko_kegg_reaction_ids: Set[str]

        return ko_kegg_reaction_ids

    def _get_ko_ec_numbers(self, ko_info: Dict[str, Any]) -> Set[str]:
        """
        Get the set of EC numbers associated with the KO under consideration in network
        construction.

        Parameters
        ==========
        ko_info : Dict[str, Any]
            Information on the KO loaded from the anvi'o KEGG database.

        Returns
        =======
        Set[str]
            EC numbers associated with the KO.
        """
        try:
            ko_ec_numbers = set(ko_info['EC'])
        except KeyError:
            # The KO is not associated with EC numbers.
            ko_ec_numbers = set()
        ko_ec_numbers: Set[str]

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
        old_kegg_reaction_ids = []
        new_kegg_reaction_ids = []
        for kegg_reaction_id in ko_kegg_reaction_ids:
            if kegg_reaction_id in network.kegg_modelseed_aliases:
                old_kegg_reaction_ids.append(kegg_reaction_id)
            else:
                new_kegg_reaction_ids.append(kegg_reaction_id)
        old_kegg_reaction_ids = set(old_kegg_reaction_ids)
        new_kegg_reaction_ids = set(new_kegg_reaction_ids)

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
            EC numbers associated with the KO.

        Returns
        =======
        Tuple[List[str], List[str]]
            A set of EC numbers in the network and a set of those not in the network.
        """
        old_ec_numbers = []
        new_ec_numbers = []
        for ec_number in ko_ec_numbers:
            if ec_number in network.ec_number_modelseed_aliases:
                old_ec_numbers.append(ec_number)
            else:
                new_ec_numbers.append(ec_number)
        old_ec_numbers = set(old_ec_numbers)
        new_ec_numbers = set(new_ec_numbers)

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
            Reaction network under construction.

        modelseed_compounds_table : pandas.core.frame.DataFrame
            Loaded compounds table of ModelSEED Biochemistry database set up by anvi'o.

        ko : KO
            KO being added to the network.

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
                    reaction, metabolites = self._get_modelseed_reaction(
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
                for metabolite in metabolites:
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
            ko.reaction_ids.append(modelseed_reaction_id)

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
                    reaction, reaction_metabolites = self._get_modelseed_reaction(
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
                for metabolite in reaction_metabolites:
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
            ko.reaction_ids.append(modelseed_reaction_id)

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
    ) -> Tuple[ModelSEEDReaction, List[ModelSEEDCompound]]:
        """
        Get an object representation of the ModelSEED reaction and object representations of the
        associated ModelSEED compounds involved in the reaction.

        Parameters
        ==========
        modelseed_reaction_data : Dict
            Dictionary representation of a row of the ModelSEED reaction table set up by anvi'o,
            containing data on the reaction.

        modelseed_compounds_table : pandas.core.frame.DataFrame
            Loaded ModelSEED Biochemistry compounds database.

        network : ReactionNetwork, None
            Reaction network under construction, with reaction compound objects drawn from the
            network, if possible, rather than created anew, as is done when a network is not
            provided. New reaction and compound objects are not added to the network by this method.

        Returns
        =======
        ModelSEEDReaction
            Representation of the reaction with data sourced from ModelSEED Biochemistry.

        List[ModelSEEDCompound]
            Representations of metabolites involved in the reaction, with data sourced from
            ModelSEED Biochemistry.
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
            reaction.compound_ids = modelseed_compound_ids

        if not modelseed_compound_ids:
            return reaction, []

        reaction_metabolites: List[ModelSEEDCompound] = []
        for compound_id in modelseed_compound_ids:
            if network:
                try:
                    # The ModelSEED compound ID has been encountered in previously processed
                    # reactions, so there is already a 'ModelSEEDCompound' object for it.
                    reaction_metabolites.append(network.metabolites[compound_id])
                    continue
                except KeyError:
                    pass

            # Generate a new metabolite object.
            try:
                modelseed_compound_series: pd.Series = modelseed_compounds_table.loc[compound_id]
            except KeyError:
                raise ConfigError(
                    f"A row for the ModelSEED compound ID, '{compound_id}', was expected but not "
                    "found in the ModelSEED compounds table. This ID was found in the equation for "
                    f"the ModelSEED reaction, '{modelseed_reaction_id}'."
                )
            modelseed_compound_data = modelseed_compound_series.to_dict()
            modelseed_compound_data['id'] = compound_id
            metabolite = self._get_modelseed_compound(modelseed_compound_data)
            reaction_metabolites.append(metabolite)

        return reaction, reaction_metabolites

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
            Reaction network under construction.

        ko : KO
            KO being added to the network.

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
                ko.reaction_ids.append(modelseed_reaction_id)
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
                if modelseed_reaction_id in ko.reaction_ids:
                    # The ModelSEED reaction has already been associated with the KO, as it was
                    # aliased by a KEGG reaction referenced by the KO -- addressed above -- as well
                    # as the EC number.
                    continue
                reaction = network.reactions[modelseed_reaction_id]
                ko.reaction_ids.append(modelseed_reaction_id)
                if DEBUG:
                    assert not ko_kegg_reaction_ids.intersection(set(reaction.kegg_aliases))
                ko.kegg_reaction_aliases[modelseed_reaction_id] = []
                ko.ec_number_aliases[modelseed_reaction_id] = list(
                    ko_ec_numbers.intersection(set(reaction.ec_number_aliases))
                )

    def _add_ko_classification(
        self,
        ko: KO,
        network: ReactionNetwork,
        ko_info: Dict[str, Any],
        kegg_modules_data: Dict[str, Dict[str, Any]],
        kegg_pathways_data: Dict[str, Dict[str, Any]],
        kegg_hierarchies_data: Dict[str, str]
    ) -> None:
        """
        Add KEGG classifications of the KO (modules, pathways, and BRITE hierarchies) to the
        network under construction.

        Parameters
        ==========
        ko : KO
            KO being added to the network.

        network : ReactionNetwork
            Reaction network under construction.

        ko_info : Dict[str, Any]
            Information on the KO loaded from the anvi'o KEGG database.

        kegg_modules_data : Dict[str, Dict[str, Any]]
            This dictionary of KEGG reference data relates module IDs to module names and pathways.

        kegg_pathways_data : Dict[str, Dict[str, Any]]
            This dictionary of KEGG reference data relates pathway IDs to pathway names and
            equivalent categories in the BRITE hierarchy, 'ko00001'.

        kegg_hierarchies_data : Dict[str, str]
            This dictionary of KEGG reference data relates BRITE hierarchy IDs to hierarchy names.

        Returns
        =======
        None
        """
        ko_id = ko.id

        # Reference module IDs in the KO.
        ko_info_mod: Tuple[str] = ko_info['MOD']
        for module_id in ko_info_mod:
            ko.module_ids.append(module_id)

        # Reference pathway IDs in the KO.
        ko_info_pth: Tuple[str] = ko_info['PTH']
        for pathway_id in ko_info_pth:
            ko.pathway_ids.append(pathway_id)

        # Reference BRITE hierarchy categorizations in the KO.
        ko_info_hie: Dict[str, Tuple[Tuple[str]]] = ko_info['HIE']
        for hierarchy_id, categorizations in ko_info_hie.items():
            ko_hierarchy_categorizations: List[Tuple[str]] = []
            for categorization in categorizations:
                ko_hierarchy_categorizations.append(categorization)
            ko.hierarchies[hierarchy_id] = ko_hierarchy_categorizations

        # Fill out module objects in the network.
        # Track module membership of pathways to facilitate pathway object creation.
        pathway_modules: Dict[str, List[str]] = {}
        for module_id in ko.module_ids:
            try:
                # The module has already been added to the network via another KO.
                module = network.modules[module_id]
            except KeyError:
                # Create a module object and add it to the network.
                module_info = kegg_modules_data[module_id]
                module = KEGGModule(id=module_id)
                module.name = module_info['NAME']
                for pathway_id in module_info['PTH']:
                    if pathway_id not in kegg_pathways_data:
                        # Only pathways that are equivalent to categories in the KO BRITE hierarchy,
                        # 'ko00001', are considered.
                        continue
                    module.pathway_ids.append(pathway_id)
                    try:
                        pathway_modules[pathway_id].append(module_id)
                    except KeyError:
                        pathway_modules[pathway_id] = [module_id]
                network.modules[module_id] = module
            module.ko_ids.append(ko_id)

        # Fill out pathway objects in the network.
        # Track BRITE categories that are equivalent to pathways to facilitate category object
        # creation.
        category_pathways: Dict[Tuple[str], str] = {}
        for pathway_id in ko.pathway_ids:
            try:
                # The pathway has already been added to the network via another KO.
                pathway = network.pathways[pathway_id]
            except KeyError:
                # Create a pathway object and add it to the network.
                pathway_info = kegg_pathways_data[pathway_id]
                pathway = KEGGPathway(id=pathway_id)
                pathway.name = pathway_info['NAME']
                categorization = pathway_info['CAT']
                pathway.categorization = categorization
                category_pathways[categorization] = pathway_id
                network.pathways[pathway_id] = pathway

            try:
                # The KO was associated with newly encountered modules, which need to be referenced
                # by the pathway.
                module_ids = pathway_modules[pathway_id]
                for module_id in module_ids:
                    pathway.module_ids.append(module_id)
            except KeyError:
                pass

            pathway.ko_ids.append(ko_id)

        # Fill out hierarchy objects in the network.
        for hierarchy_id, categorizations in ko.hierarchies.items():
            try:
                # The hierarchy has already been added to the network via another KO.
                hierarchy = network.hierarchies[hierarchy_id]
                network_hierarchy_categories = network.categories[hierarchy_id]
            except KeyError:
                # Create a new hierarchy object and add it to the network.
                hierarchy_name = kegg_hierarchies_data[hierarchy_id]
                hierarchy = BRITEHierarchy(id=hierarchy_id)
                hierarchy.name = hierarchy_name
                network.hierarchies[hierarchy_id] = hierarchy
                network_hierarchy_categories: Dict[Tuple[str], Tuple[BRITECategory]] = {}
                network.categories[hierarchy_id] = network_hierarchy_categories

            hierarchy.ko_ids.append(ko_id)

            # Fill out category objects in the network.
            for categorization in categorizations:
                try:
                    # The category has already been added to the network via another KO.
                    categories = network_hierarchy_categories[categorization]
                except KeyError:
                    # Add a category object to the network for each level of the categorization.
                    categories: List[BRITECategory] = []
                    for depth, focus_category_name in enumerate(categorization, 1):
                        focus_categorization = categorization[:depth]
                        try:
                            # The supercategory has already been added to the network.
                            focus_categories = network_hierarchy_categories[focus_categorization]
                            categories = list(focus_categories)
                            continue
                        except KeyError:
                            pass

                        # Add the previously unencountered category to the network.
                        category = BRITECategory()
                        category.id = f'{hierarchy_id}: {">>>".join(focus_categorization)}'
                        category.name = focus_category_name
                        category.hierarchy_id = hierarchy_id
                        if depth == len(categorization) and hierarchy_id == 'ko00001':
                            try:
                                category.pathway_id = category_pathways[categorization]
                            except KeyError:
                                pass
                        categories.append(category)

                        if len(categories) > 1:
                            # Consider the supercategory of the newly encountered category. Add the
                            # category name as a subcategory reference of the supercategory.
                            categories[-2].subcategory_names.append(focus_category_name)

                        hierarchy.categorizations.append(focus_categorization)
                        network_hierarchy_categories[focus_categorization] = tuple(categories)

    def _relate_modules_pathways(
        self,
        network: ReactionNetwork,
        kegg_modules_data: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Link modules and pathways.

        Certain KOs but not others in a module can be in a pathway. Only KOs in the network are
        relevant. A module is only linked to pathways via KOs in the network, so relationships
        between modules and pathways in the network are only resolved here after all KOs have been
        added to the network.

        Parameters
        ==========
        network : ReactionNetwork
            Reaction network under construction.

        kegg_modules_data : Dict[str, Dict[str, Any]]
            This dictionary of KEGG reference data relates module IDs to module names and pathways.

        Returns
        =======
        None
        """
        for module_id, module in network.modules.items():
            module_info = kegg_modules_data[module_id]
            for pathway_id in module_info['PTH']:
                try:
                    pathway = network.pathways[pathway_id]
                except KeyError:
                    continue
                module.pathway_ids.append(pathway_id)
                pathway.module_ids.append(module_id)

    def _get_database_reactions_table(self, network: ReactionNetwork) -> pd.DataFrame:
        """
        Make a reactions table that can be stored in either a contigs or pan database, as the tables
        have the same structure. A ReactionNetwork can be reconstructed with the same data from the
        reactions, metabolites, and KEGG tables of the database.

        Parameters
        ==========
        network : ReactionNetwork
            Network generated from gene or gene cluster KO annotations.

        Returns
        =======
        pd.DataFrame
            Table of reactions data to be stored in the contigs or pan database.
        """
        if DEBUG:
            assert (
                tables.reaction_network_reactions_table_structure ==
                tables.pan_reaction_network_reactions_table_structure
            )
            assert (
                tables.reaction_network_reactions_table_types ==
                tables.pan_reaction_network_reactions_table_types
            )

        # Transfer data from reaction objects to dictionaries mapping to table entries.
        reactions_data: Dict[str, Dict] = {}
        for modelseed_reaction_id, reaction in network.reactions.items():
            reaction_data = {}
            reaction_data['modelseed_reaction_id'] = modelseed_reaction_id
            reaction_data['modelseed_reaction_name'] = reaction.modelseed_name
            reaction_data['metabolite_modelseed_ids'] = ', '.join(reaction.compound_ids)
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
            for modelseed_reaction_id in ko.reaction_ids:
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
        reactions_table = reactions_table[tables.reaction_network_reactions_table_structure]

        return reactions_table

    def _get_database_metabolites_table(self, network: ReactionNetwork) -> pd.DataFrame:
        """
        Make a metabolites table that can be stored in either a contigs or pan database, as the
        tables have the same structure. A ReactionNetwork can be reconstructed with the same data
        from the reactions, metabolites, and KEGG tables of the database.

        Parameters
        ==========
        network : ReactionNetwork
            Network generated from gene or gene cluster KO annotations.

        Returns
        =======
        pd.DataFrame
            Table of metabolites data to be stored in the contigs or pan database.
        """
        if DEBUG:
            assert (
                tables.reaction_network_metabolites_table_structure ==
                tables.pan_reaction_network_metabolites_table_structure
            )
            assert (
                tables.reaction_network_metabolites_table_types ==
                tables.pan_reaction_network_metabolites_table_types
            )

        # Transfer data from metabolite objects to dictionaries mapping to table entries.
        metabolites_data = {}
        for modelseed_compound_id, metabolite in network.metabolites.items():
            metabolite_data = {}
            metabolite_data['modelseed_compound_id'] = modelseed_compound_id
            metabolite_data['modelseed_compound_name'] = metabolite.modelseed_name
            metabolite_data['kegg_aliases'] = ', '.join(metabolite.kegg_aliases)
            metabolite_data['formula'] = metabolite.formula
            metabolite_data['charge'] = metabolite.charge
            metabolites_data[modelseed_compound_id] = metabolite_data

        metabolites_table = pd.DataFrame.from_dict(
            metabolites_data, orient='index'
        ).reset_index(drop=True).sort_values('modelseed_compound_id')
        metabolites_table = metabolites_table[tables.reaction_network_metabolites_table_structure]

        return metabolites_table

    def _get_database_kegg_table(self, network: ReactionNetwork) -> pd.DataFrame:
        """
        Make a table recording the relationships between KEGG KOs, modules, pathways, and BRITE
        hierarchies in the reaction network that can be stored in either a contigs or a pan
        database, as tables have the same structure. A ReactionNetwork can be reconstructed with the
        same data from the reaction, metabolites, and KEGG tables of the database.

        Parameters
        ==========
        network : ReactionNetwork
            Network generated from gene or gene cluster KO annotations.

        Returns
        =======
        pd.DataFrame
            Table of KEGG information to be stored.
        """
        if DEBUG:
            assert (
                tables.reaction_network_kegg_table_structure ==
                tables.pan_reaction_network_kegg_table_structure
            )
            assert (
                tables.reaction_network_kegg_table_types ==
                tables.pan_reaction_network_kegg_table_types
            )

        # Transfer data from KEGG objects to dictionaries mapping to table entries.
        kegg_data = {}

        # The first rows in the table are for KOs.
        ko_id_pattern = re.compile('K\d{5}')
        for ko_id, ko in network.kos.items():
            ko_data = {}
            assert re.fullmatch(ko_id_pattern, ko_id)
            ko_data['id'] = ko_id
            ko_data['name'] = ko.name
            ko_data['modules'] = ', '.join(ko.module_ids)
            ko_data['pathways'] = ', '.join(ko.pathway_ids)
            brite_categorizations = []
            for hierarchy_id, categorizations in ko.hierarchies.items():
                for categorization in categorizations:
                    brite_categorizations.append(
                        f'{hierarchy_id} >>> {" >>> ".join(categorization)}'
                    )
            ko_data['brite_categorization'] = ' !!! '.join(brite_categorizations)
            kegg_data[f'1{ko_id}'] = ko_data

        # Modules are second in the table.
        module_id_pattern = re.compile('M\d{5}')
        for module_id, module in network.modules.items():
            module_data = {}
            assert re.fullmatch(module_id_pattern, module_id)
            module_data['id'] = module_id
            module_data['name'] = module.name
            module_data['modules'] = ''
            module_data['pathways'] = ', '.join(module.pathway_ids)
            # TODO: There is a BRITE hierarchy of modules, 'ko00002'. This hierarchy is not
            # currently downloaded or handled by kegg.py, but the classification of modules, like
            # that of pathways as categories in the hierarchy, 'ko00001', could be useful.
            module_data['brite_categorization'] = ''
            kegg_data[f'2{module_id}'] = module_data

        # Pathways are third in the table.
        pathway_id_pattern = re.compile('map\d{5}')
        for pathway_id, pathway in network.pathways.items():
            pathway_data = {}
            assert re.fullmatch(pathway_id_pattern, pathway_id)
            pathway_data['id'] = pathway_id
            pathway_data['name'] = pathway.name
            pathway_data['modules'] = ''
            pathway_data['pathways'] = ''
            pathway_data[
                'brite_categorization'
            ] = f'ko00001 >>> {" >>> ".join(pathway.categorization)}'
            kegg_data[f'3{pathway_id}'] = pathway_data

        # Hierarchies are fourth in the table.
        hierarchy_id_pattern = re.compile('ko\d{5}')
        for hierarchy_id, hierarchy in network.hierarchies.items():
            hierarchy_data = {}
            # Only hierarchies of KOs should be in consideration. Hierarchies of other KEGG items
            # that do not resolve to KOs, such as reactions and drugs, have IDs that start with 'br'
            # rather than 'ko'.
            assert re.fullmatch(hierarchy_id_pattern, hierarchy_id)
            hierarchy_data['id'] = hierarchy_id
            hierarchy_data['name'] = hierarchy.name
            hierarchy_data['modules'] = ''
            hierarchy_data['pathways'] = ''
            hierarchy_data['brite_categorization'] = ''
            kegg_data[f'4{hierarchy_id}'] = hierarchy_data

        kegg_table = pd.DataFrame.from_dict(
            kegg_data, orient='index'
        ).sort_index().reset_index(drop=True)
        kegg_table = kegg_table[tables.reaction_network_kegg_table_structure]

        return kegg_table

    def hash_contigs_db_ko_hits(self, gene_ko_hits_table: pd.DataFrame) -> str:
        """
        To concisely represent the data underlying a reaction network, hash all gene KO annotations
        in the contigs database.

        Parameters
        ==========
        gene_ko_hits_table : pandas.core.frame.DataFrame
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
        functions_table = gsdb.get_table_as_dataframe(
            'gene_function_calls', where_clause='source = "KOfam"'
        )
        gsdb.disconnect()
        ko_annotations = []
        for row in functions_table.itertuples(index=False):
            ko_annotations.append((
                row.genome_name,
                str(row.gene_callers_id),
                row.accession,
                row.function,
                str(row.e_value)
            ))
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
    kegg_dir : str, None
        Directory containing an anvi'o KEGG database. This attribute is assigned the argument of the
        same name upon initialization.

    modelseed_dir : str, None
        Directory containing reference ModelSEED Biochemistry tables set up by anvi'o. This
        attribute is assigned the argument of the same name upon initialization.

    test_dir : str, None
        Directory storing test files, including copied input and output files. With the default
        value of None, temporary directories are created and deleted as needed by methods. In
        contrast, if a directory is provided, it and its contents will not be deleted. This
        attribute is assigned the argument of the same name upon initialization.

    run : anvio.terminal.Run, anvio.terminal.Run()
        This object prints run information to the terminal. This attribute is assigned the argument
        of the same name upon initialization.

    progress : anvio.terminal.Progress, anvio.terminal.Progress()
        This object prints transient progress information to the terminal. This attribute is
        assigned the argument of the same name upon initialization.
    """
    def __init__(
        self,
        kegg_dir: str = None,
        modelseed_dir: str = None,
        test_dir: str = None,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Parameters
        ==========
        kegg_dir : str, None
            Directory containing an anvi'o KEGG database. The default argument of None expects KEGG
            data to be set up in the default anvi'o directory used by the program,
            `anvi-setup-kegg-data`.

        modelseed_dir : str, None
            Directory containing reference ModelSEED Biochemistry tables set up by anvi'o. The
            default argument of None expects ModelSEED data to be set up in the default anvi'o
            directory used by the program, `anvi-setup-modelseed-database`.

        test_dir : str, None
            Directory storing test files. With the default value of None, temporary test directories
            are created and deleted by Tester methods; these methods operate on copies of input
            files in the test directories. In contrast, if a directory is provided, it and its
            contents will not be deleted.

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
        """
        self.kegg_dir = kegg_dir
        self.modelseed_dir = modelseed_dir
        self.test_dir = test_dir
        self.run = run
        self.progress = progress

    def test_contigs_database_network(self, contigs_db: str, copy_db: bool = True) -> None:
        """
        Test the construction of a reaction network from a contigs database, and test that network
        methods are able to run and do not fail certain basic tests.

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
            database, overwriting any network that is already stored.

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
            kegg_dir=self.kegg_dir,
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

        # Check that the statistics for the network constructed and saved in the contigs database
        # are the same as the statistics for the same network loaded back into memory from the
        # contigs database.
        inconsistent_stats = self._get_inconsistent_statistics(
            make_stats_file_target, load_stats_file_target
        )
        if inconsistent_stats:
            msg = ""
            for stat, stat_tuple in inconsistent_stats.items():
                msg += f"{stat}: {stat_tuple[0]}, {stat_tuple[1]}; "
            msg = msg[:-2]
            raise AssertionError(
                "Statistics on the network constructed and saved to the contigs database differ "
                "from what should be the same statistics on the same network loaded from the "
                "contigs database. Here are the different statistics, with the value from network "
                f"construction before the value from network loading: {msg}"
            )

        self.run.info_single(
            "PURGE OF METABOLITES WITHOUT FORMULA:", mc='magenta', nl_before=1, level=0
        )
        deepcopy(network).remove_metabolites_without_formula(
            output_path=os.path.join(test_dir, "removed.tsv")
        )
        print()

        self.progress.new("Testing network purge methods")
        self.progress.update("...")
        # Network pruning tests use a random sample of half the network items (nodes) of each type.
        sample_proportion = 0.5
        sample_seed = RANDOM_SEED
        samples = self._get_common_item_samples(
            network, proportion=sample_proportion, seed=sample_seed
        )

        random.seed(sample_seed)
        gene_sample = set(random.sample(
            list(network.genes), sample_proportion * math.ceil(len(network.genes) / 2)
        ))

        self._test_common_prune(network, samples)

        copied_network = deepcopy(network)
        removed = copied_network.prune(genes_to_remove=gene_sample)
        assert gene_sample.difference(set(copied_network.genes)) == gene_sample
        assert not gene_sample.difference(set([gene.gcid for gene in removed['gene']]))
        self.progress.end()

        self.progress.new("Testing network subset methods")
        self.progress.update("...")
        self._test_common_subset(network, samples)

        subnetwork = network.subset_network(genes_to_subset=gene_sample)
        assert not gene_sample.difference(set(subnetwork.genes))

        # Test network merging functionality by subsetting samples of items of all types.
        metabolite_sample: Set[str] = samples['metabolite']
        reaction_sample: Set[str] = samples['reaction']
        ko_sample: Set[str] = samples['ko']
        module_sample: Set[str] = samples['module']
        pathway_sample: Set[str] = samples['pathway']
        hierarchy_sample: Set[str] = samples['hierarchy']
        category_sample_dict: Dict[str, List[Tuple[str]]] = samples['category_dict']
        subnetwork = network.subset_network(
            genes_to_subset=gene_sample,
            kos_to_subset=ko_sample,
            modules_to_subset=module_sample,
            pathways_to_subset=pathway_sample,
            hierarchies_to_subset=hierarchy_sample,
            categories_to_subset=category_sample_dict,
            reactions_to_subset=reaction_sample,
            metabolites_to_subset=metabolite_sample
        )
        assert not metabolite_sample.difference(set(subnetwork.metabolites))
        assert not reaction_sample.difference(set(subnetwork.reactions))
        assert not ko_sample.difference(set(subnetwork.kos))
        assert not module_sample.difference(set(subnetwork.modules))
        assert not pathway_sample.difference(set(subnetwork.pathways))
        assert not hierarchy_sample.difference(set(subnetwork.hierarchies))
        remaining_category_ids: List[str] = []
        for categorizations in subnetwork.categories.values():
            for categories in categorizations.values():
                remaining_category_ids.append(categories[-1].id)
        category_sample: Set[str] = samples['category']
        assert not category_sample.difference(set(remaining_category_ids))
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
        self.run.info_single("Purge KOs in select KEGG modules")
        self.run.info_single("Purge KOs in select KEGG pathways")
        self.run.info_single("Purge KOs in select KEGG BRITE hierarchies")
        self.run.info_single("Purge KOs in select KEGG BRITE hierarchy categories")
        self.run.info_single("Purge select genes")
        self.run.info_single("Subset select metabolites")
        self.run.info_single("Subset select reactions")
        self.run.info_single("Subset select KOs")
        self.run.info_single("Subset KOs in select KEGG modules")
        self.run.info_single("Subset KOs in select KEGG pathways")
        self.run.info_single("Subset KOs in select KEGG BRITE hierarchies")
        self.run.info_single("Subset KOs in select KEGG BRITE hierarchy categories")
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
        methods are able to run and do not fail certain basic tests.

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
            If False, store the reaction network in the input pan database, overwriting any network
            that is already stored.

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
            kegg_dir=self.kegg_dir,
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

        # Check that the statistics for the network constructed and saved in the pan database are
        # the same as the statistics for the same network loaded back into memory from the pan
        # database.
        inconsistent_stats = self._get_inconsistent_statistics(
            make_stats_file_target, load_stats_file_target
        )
        if inconsistent_stats:
            msg = ""
            for stat, stat_tuple in inconsistent_stats.items():
                msg += f"{stat}: {stat_tuple[0]}, {stat_tuple[1]}; "
            msg = msg[:-2]
            raise AssertionError(
                "Statistics on the network constructed and saved to the pan database differ from "
                "what should be the same statistics on the same network loaded from the pan "
                "database. Here are the different statistics, with the value from network "
                f"construction before the value from network loading: {msg}"
            )

        self.run.info_single(
            "PURGE OF METABOLITES WITHOUT FORMULA:", mc='magenta', nl_before=1, level=0
        )
        deepcopy(network).remove_metabolites_without_formula(
            output_path=os.path.join(test_dir, "removed.tsv")
        )
        print()

        self.progress.new("Testing network purge methods")
        self.progress.update("...")
        # Network pruning tests use a random sample of half the network items (nodes) of each type.
        sample_proportion = 0.5
        sample_seed = RANDOM_SEED
        samples = self._get_common_item_samples(
            network, proportion=sample_proportion, seed=sample_seed
        )

        random.seed(sample_seed)
        gene_cluster_sample = set(random.sample(
            list(network.gene_clusters), sample_proportion * math.ceil(len(network.gene_clusters))
        ))

        self._test_common_prune(network, samples)

        copied_network = deepcopy(network)
        removed = copied_network.prune(gene_clusters_to_remove=gene_cluster_sample)
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
        self._test_common_subset(network, samples)

        subnetwork = network.subset_network(gene_clusters_to_subset=gene_cluster_sample)
        assert not gene_cluster_sample.difference(set(subnetwork.gene_clusters))

        # Test network merging functionality by subsetting samples of items of all types.
        metabolite_sample: Set[str] = samples['metabolite']
        reaction_sample: Set[str] = samples['reaction']
        ko_sample: Set[str] = samples['ko']
        module_sample: Set[str] = samples['module']
        pathway_sample: Set[str] = samples['pathway']
        hierarchy_sample: Set[str] = samples['hierarchy']
        category_sample_dict: Dict[str, List[Tuple[str]]] = samples['category_dict']
        subnetwork = network.subset_network(
            gene_clusters_to_subset=gene_cluster_sample,
            kos_to_subset=ko_sample,
            modules_to_subset=module_sample,
            pathways_to_subset=pathway_sample,
            hierarchies_to_subset=hierarchy_sample,
            categories_to_subset=category_sample_dict,
            reactions_to_subset=reaction_sample,
            metabolites_to_subset=metabolite_sample
        )
        assert not metabolite_sample.difference(set(subnetwork.metabolites))
        assert not reaction_sample.difference(set(subnetwork.reactions))
        assert not ko_sample.difference(set(subnetwork.kos))
        assert not module_sample.difference(set(subnetwork.modules))
        assert not pathway_sample.difference(set(subnetwork.pathways))
        assert not hierarchy_sample.difference(set(subnetwork.hierarchies))
        remaining_category_ids: List[str] = []
        for categorizations in subnetwork.categories.values():
            for categories in categorizations.values():
                remaining_category_ids.append(categories[-1].id)
        category_sample: Set[str] = samples['category']
        assert not category_sample.difference(set(remaining_category_ids))
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
        self.run.info_single("Purge KOs in select KEGG modules")
        self.run.info_single("Purge KOs in select KEGG pathways")
        self.run.info_single("Purge KOs in select KEGG BRITE hierarchies")
        self.run.info_single("Purge KOs in select KEGG BRITE hierarchy categories")
        self.run.info_single("Purge select gene clusters")
        self.run.info_single("Subset select metabolites")
        self.run.info_single("Subset select reactions")
        self.run.info_single("Subset select KOs")
        self.run.info_single("Subset KOs in select KEGG modules")
        self.run.info_single("Subset KOs in select KEGG pathways")
        self.run.info_single("Subset KOs in select KEGG BRITE hierarchies")
        self.run.info_single("Subset KOs in select KEGG BRITE hierarchy categories")
        self.run.info_single("Subset select gene clusters")
        self.run.info_single(
            "Subset select metabolites, reactions, KOs, and gene clusters", nl_after=1
        )

    def _get_inconsistent_statistics(
        self,
        make_stats_file_target: str,
        load_stats_file_target: str
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compare statistics for the network constructed and saved in an anvi'o database with
        statistics for what should be the same network loaded back into memory from the database,
        returning any inconsistent statistics.

        Parameters
        ==========
        make_stats_file_target : str
            Path to file of statistics for network construction.

        load_stats_file_target : str
            Path to file of statistics for network loading.

        Returns
        =======
        Dict[str, Tuple[float, float]]
            Inconsistent statistics, with keys being names of the statistics and values being pairs
            of the statistic from network construction and loading, respectively.
        """
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

        return inconsistent_stats

    def _get_common_item_samples(
        self,
        network: ReactionNetwork,
        proportion: float = 0.5,
        seed: int = 0
    ) -> Dict:
        """
        Get a random sample of half of the network items (nodes) of each type shared in common
        between genomic and pangenomic networks.

        Parameters
        ==========
        network : ReactionNetwork
            Network generated from gene or gene cluster KO annotations.

        proportion : float, 0.5
            Proportion to be randomly sampled of network items of each type.

        seed : int, 0
            Seed for random number generation.

        Returns
        =======
        dict
            Dictionary with keys indicating item type and values being samples of item IDs.
        """
        samples = {}

        random.seed(seed)
        samples['metabolite'] = set(random.sample(
            list(network.metabolites), math.ceil(proportion * len(network.metabolites))
        ))

        random.seed(seed)
        samples['reaction'] = set(
            random.sample(list(network.reactions), math.ceil(proportion * len(network.reactions)))
        )

        random.seed(seed)
        samples['ko'] = set(
            random.sample(list(network.kos), math.ceil(proportion * len(network.kos)))
        )

        random.seed(seed)
        samples['module'] = set(
            random.sample(list(network.modules), math.ceil(proportion * len(network.modules)))
        )

        random.seed(seed)
        samples['pathway'] = set(
            random.sample(list(network.pathways), math.ceil(proportion * len(network.pathways)))
        )

        random.seed(seed)
        samples['pathway'] = set(random.sample(
            list(network.hierarchies), math.ceil(proportion * len(network.hierarchies))
        ))

        random.seed(seed)
        all_category_ids: List[str] = []
        for categorizations in network.categories.values():
            for categories in categorizations.values():
                all_category_ids.append(categories[-1].id)
        samples['category'] = category_sample = set(
            random.sample(all_category_ids, math.ceil(proportion * len(all_category_ids)))
        )
        # Reformat the category IDs into an argument for pruning and subsetting.
        category_sample_dict: Dict[str, List[Tuple[str]]] = {}
        hierarchy_id_pattern = re.compile('ko\d{5}')
        for category_id in category_sample:
            hierarchy_id = category_id.split(':')[0]
            assert re.fullmatch(hierarchy_id_pattern, hierarchy_id)
            try:
                categorizations = category_sample_dict[hierarchy_id]
            except KeyError:
                category_sample_dict[hierarchy_id] = categorizations = []
            categorizations.append(tuple(category_id[len(hierarchy_id) + 2:].split(' >>> ')))
        samples['category_dict'] = category_sample_dict

        return samples

    def _test_common_prune(
        self,
        network: Union[GenomicNetwork, PangenomicNetwork],
        samples: Dict
    ) -> None:
        """
        Test the prune method of the reaction network, purging items of types in common to genomic
        and pangenomic networks.

        Parameters
        ==========
        network : Union[GenomicNetwork, PangenomicNetwork]
            Network generated from gene or gene cluster KO annotations.

        samples : Dict
            Dictionary with keys indicating item type and values being samples of item IDs.

        Returns
        =======
        None
        """
        copied_network = deepcopy(network)
        metabolite_sample: Set[str] = samples['metabolite']
        removed = copied_network.prune(metabolites_to_remove=metabolite_sample)
        # The most basic test of the purge (pruning) method is that the network no longer contains
        # the items that were requested to be removed. What remains untested, and would require a
        # curated test dataset, is the removal of certain other "upstream" and "downstream" nodes
        # associated with the nodes requested to be removed, e.g., KOs and genes or gene clusters
        # upstream and metabolites downstream of requested reactions.
        assert metabolite_sample.difference(set(copied_network.metabolites)) == metabolite_sample
        assert not metabolite_sample.difference(
            set([metabolite.modelseed_id for metabolite in removed['metabolite']])
        )

        copied_network = deepcopy(network)
        reaction_sample: Set[str] = samples['reaction']
        removed = copied_network.prune(reactions_to_remove=reaction_sample)
        assert reaction_sample.difference(set(copied_network.reactions)) == reaction_sample
        assert not reaction_sample.difference(
            set([reaction.modelseed_id for reaction in removed['reaction']])
        )

        copied_network = deepcopy(network)
        ko_sample: Set[str] = samples['ko']
        removed = copied_network.prune(kos_to_remove=ko_sample)
        assert ko_sample.difference(set(copied_network.kos)) == ko_sample
        assert not ko_sample.difference(set([ko.id for ko in removed['ko']]))

        copied_network = deepcopy(network)
        module_sample: Set[str] = samples['module']
        removed = copied_network.prune(modules_to_remove=module_sample)
        assert module_sample.difference(set(copied_network.modules)) == module_sample
        assert not module_sample.difference(set([module.id for module in removed['module']]))

        copied_network = deepcopy(network)
        pathway_sample: Set[str] = samples['pathway']
        removed = copied_network.prune(pathways_to_remove=pathway_sample)
        assert pathway_sample.difference(set(copied_network.pathways)) == pathway_sample
        assert not pathway_sample.difference(set([pathway.id for pathway in removed['pathway']]))

        copied_network = deepcopy(network)
        hierarchy_sample: Set[str] = samples['hierarchy']
        removed = copied_network.prune(hierarchies_to_remove=hierarchy_sample)
        assert hierarchy_sample.difference(set(copied_network.hierarchies)) == hierarchy_sample
        assert not hierarchy_sample.difference(
            set([hierarchy.id for hierarchy in removed['hierarchies']])
        )

        copied_network = deepcopy(network)
        category_sample_dict: Dict[str, List[Tuple[str]]] = samples['category_dict']
        removed = copied_network.prune(categories_to_remove=category_sample_dict)
        remaining_category_ids: List[str] = []
        for categorizations in copied_network.categories.values():
            for categories in categorizations.values():
                remaining_category_ids.append(categories[-1].id)
        category_sample: Set[str] = samples['category']
        assert category_sample.difference(set(remaining_category_ids)) == category_sample
        removed_category_ids: List[str] = []
        for category in removed['category']:
            category: BRITECategory
            removed_category_ids.append(category.id)
        assert not category_sample.difference(set(removed_category_ids))

    def _test_common_subset(
        self,
        network: Union[GenomicNetwork, PangenomicNetwork],
        samples: Dict
    ) -> None:
        """
        Test the subset method of the reaction network, subsetting items of types in common to
        genomic and pangenomic networks.

        Parameters
        ==========
        network : Union[GenomicNetwork, PangenomicNetwork]
            Network generated from gene or gene cluster KO annotations.

        samples : Dict
            Dictionary with keys indicating item type and values being samples of item IDs.

        Returns
        =======
        None
        """
        metabolite_sample: Set[str] = samples['metabolite']
        subnetwork = network.subset_network(metabolites_to_subset=metabolite_sample)
        # The most basic test of the subset method is that the new network contains the requested
        # items. What remains untested, and would require a curated test dataset, is the inclusion
        # of certain other "upstream" and "downstream" nodes associated with the nodes requested to
        # be removed, e.g., KOs and gene clusters upstream and metabolites downstream of requested
        # reactions.
        # Assert that all of the items requested to be subsetted were subsetted.
        assert not metabolite_sample.difference(set(subnetwork.metabolites))

        reaction_sample: Set[str] = samples['reaction']
        subnetwork = network.subset_network(reactions_to_subset=reaction_sample)
        assert not reaction_sample.difference(set(subnetwork.reactions))

        ko_sample: Set[str] = samples['ko']
        subnetwork = network.subset_network(kos_to_subset=ko_sample)
        assert not ko_sample.difference(set(subnetwork.kos))

        module_sample: Set[str] = samples['module']
        subnetwork = network.subset_network(modules_to_subset=module_sample)
        assert not module_sample.difference(set(subnetwork.modules))

        pathway_sample: Set[str] = samples['pathway']
        subnetwork = network.subset_network(pathways_to_subset=pathway_sample)
        assert not pathway_sample.difference(set(subnetwork.pathways))

        hierarchy_sample: Set[str] = samples['hierarchy']
        subnetwork = network.subset_network(hierarchies_to_subset=hierarchy_sample)
        assert not hierarchy_sample.difference(set(subnetwork.hierarchies))

        category_sample_dict: Dict[str, List[Tuple[str]]] = samples['category_dict']
        subnetwork = network.subset_network(categories_to_subset=category_sample_dict)
        remaining_category_ids: List[str] = []
        for categorizations in subnetwork.categories.values():
            for categories in categorizations.values():
                remaining_category_ids.append(categories[-1].id)
        category_sample: Set[str] = samples['category']
        assert not category_sample.difference(set(remaining_category_ids))

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
    for coefficient, modelseed_compound_id, compartment in zip(
        reaction.coefficients, reaction.compound_ids, reaction.compartments
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
        equation += f"{coeff} {modelseed_compound_id} [{compartment}] + "

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
