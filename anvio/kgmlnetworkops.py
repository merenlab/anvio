import os
import sys
import pandas as pd

from copy import deepcopy
from argparse import Namespace
from typing import Any, Literal, Union
from dataclasses import dataclass, field

import anvio.kegg as kegg
import anvio.kgml as kgml
import anvio.terminal as terminal
import anvio.reactionnetwork as rn

from anvio.dbops import ContigsDatabase
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

@dataclass
class Chain:
    """
    Chain of compounds linked by reactions occurring in a KGML representation of a KEGG pathway. The
    chain may be contextualized in a reaction network, with information added to the attributes
    'gaps', 'aliased_modelseed_compounds', 'aliased_modelseed_reactions', and 'network_kos'.

    Attributes
    ==========
    kgml_compound_entries : list[kgml.Entry], []
        KGML compound entries linked by kgml_reactions oriented in kgml_reaction_directions.

    is_consumed : bool, None
        True for a "consumption chain" in which the first through the penultimate compounds in the
        chain are consumed by the reactions in the chain. False for a "production chain" in which
        the first through penultimate compounds are produced. A chain is one or the other.

    kgml_reactions : list[kgml.Reaction], []
        KGML reactions linking kgml_compound_entries. Given is_consumed and
        kgml_reaction_directions, one compound entry is a substrate and one is a product of the
        reaction.

    kgml_reaction_directions : list[bool], []
        Items correspond to kgml_reactions, with True indicating the reaction runs as given from
        substrates to products, and False indicating the reaction runs in reverse from products to
        substrates.

    gaps : list[bool], []
        Items correspond to kgml_reactions, with False indicating the reaction is genomically
        encoded by a KO in the reaction network, and True indicating the reaction is not encoded and
        thus is a gap.

    aliased_modelseed_compounds : list[tuple[anvio.reactionnetwork.ModelSEEDCompound]], []
        Items correspond to kgml_compound_entries. Each item is the tuple of ModelSEED compounds
        aliased by KEGG compound IDs in the name attribute of the KGML compound entry.

    network_kos : list[tuple[anvio.reactionnetwork.KO]], []
        Items correspond to kgml_reactions. Each item is the tuple of KOs in the reaction network
        encoding KEGG reaction IDs in the name attribute of the KGML reaction.

    aliased_modelseed_reactions : list[tuple[anvio.reactionnetwork.ModelSEEDReaction]], []
        Items correspond to kgml_reactions. Each item is the tuple of ModelSEED reactions aliased by
        KEGG reaction IDs in the name attribute of the KGML reaction.

    is_consumption_terminus : bool, None
        This value holds for the first consumed KGML compound in the chain, which is the first
        compound of kgml_compound_entries if is_consumed is True and the last if is_consumed is
        False. KGMLNetworkWalker.check_kgml_compound_entry_consumption_terminus returns True if the
        compound is a consumption terminus and False if not.

    is_production_terminus : bool, None
        This value holds for the last produced KGML compound in the chain, which is the first
        compound of kgml_compound_entries if is_consumed is False and the last if is_consumed is
        True. KGMLNetworkWalker.check_kgml_compound_entry_production_terminus returns True if the
        compound is a production terminus and False if not.

    consumption_reversibility_range : tuple[int, int], None
        Pythonic start and stop indices corresponding to a range of the first consumed KGML
        compounds in the chain, which are the first compounds of kgml_compound_entries if
        is_consumed is True and the last if is_consumed is False. This range indicates compounds
        that are linked by KGML reactions with a type attribute of "reversible". For example, if the
        compound at the beginning of a consumption chain (is_consumed True) is linked to the next
        compound by an irreversible reaction, then the value is (0, 1).

    production_reversibility_range : tuple[int, int], None
        Pythonic start and stop indices corresponding to a range of the last produced KGML compounds
        in the chain, which are the first compounds of kgml_compound_entries if is_consumed is False
        and the last if is_consumed is True. This range indicates compounds that are linked by KGML
        reactions with a type attribute of "reversible". For example, if the compound at the
        beginning of a production chain (is_consumed False) is linked to the next compound by an
        irreversible reaction, then the value is (0, 1).

    cyclic_branch_index : int, None
        Non-negative integer if the chain is "partly cyclic," meaning that the chain loops a cycle
        via the entry or exit of an uncycled compound in a "branching" reaction that is part of the
        cycle, such as acetyl-CoA entering the Krebs cycle by combining with oxaloacetate to form
        citrate using citrate synthase. The branching reaction occurs twice in the partly cyclic
        chain, once at the end. The non-negative integer is the index of the first occurrence of the
        branching KGML reaction that also ends the chain. -1 if the chain is not partly cyclic.
    """
    kgml_compound_entries: list[kgml.Entry] = field(default_factory=list)
    is_consumed: bool = None
    kgml_reactions: list[kgml.Reaction] = field(default_factory=list)
    kgml_reaction_directions: list[bool] = field(default_factory=list)
    gaps: list[bool] = field(default_factory=list)
    aliased_modelseed_compounds: list[tuple[rn.ModelSEEDCompound]] = field(default_factory=list)
    aliased_modelseed_reactions: list[tuple[rn.ModelSEEDReaction]] = field(default_factory=list)
    network_kos: list[tuple[rn.KO]] = field(default_factory=list)
    is_consumption_terminus: bool = None
    is_production_terminus: bool = None
    consumption_reversibility_range: tuple[int, int] = None
    production_reversibility_range: tuple[int, int] = None
    cyclic_branch_index: int = None

class KGMLNetworkWalker:
    """
    Walk chains of compounds linked by reactions in a KGML representation of a KEGG pathway.

    Attributes
    ==========
    kegg_pathway_number : str
        Numerical ID of the pathway to walk. The pathway must have a KEGG reaction (RN) type KGML
        file available. Valid pathways are in categories 1.0 - 1.11 as of the March 3, 2025 release
        of KEGG (see https://www.genome.jp/kegg/pathway.html).

    kegg_data : anvio.reactionnetwork.KEGGData
        Contains information on an anvi'o KEGG installation, which is assumed to be at the default
        location. A reaction network stored in the network attribute should be constructed from that
        version of the database.

    kgml_rn_pathway : anvio.kgml.Pathway
        Loaded KEGG reaction (RN) type KGML file for the ID number.

    kgml_ko_pathway : anvio.kgml.Pathway
        Loaded KEGG KO type KGML file for the ID number.

    kgml_ec_pathway : anvio.kgml.Pathway
        Loaded KEGG EC number type KGML file for the ID number.

    rn_pathway_kgml_compound_id_to_kgml_reactions : dict[str, list[kgml.Reaction]]
        From the RN type KGML pathway, map KGML compound IDs to the IDs of KGML reactions involving
        the compound as substrate or product.

    pathway_nonenzymatic_kgml_reaction_ids : list[str]
        The RN, KO, and EC type KGML pathways are compared to find nonenzymatic reactions in the RN
        file that are not in the KO or EC files, the KGML IDs of which are recorded in this list.

    rn_pathway_kgml_compound_id_to_kgml_compound_entry : dict[str, list[kgml.Entry]]
        Map KGML compound IDs to the corresponding entries in the RN type pathway.

    ko_pathway_kgml_ortholog_id_to_kgml_ortholog_entry : dict[str, kgml.Entry]
        Map KGML ortholog entry IDs to the corresponding entries in the KO type pathway.

    rn_pathway_kgml_reaction_id_to_ko_ids : dict[str, list[str]]
        Map KGML reaction IDs to KO IDs. The reactions are from the RN type KGML pathway and the KO
        IDs are from the name attribute of a corresponding ortholog entry, if one exists, in the KO
        type KGML pathway.

    rn_pathway_keggcpd_ids_in_kgml_reactions : list[str]
        KEGG compound IDs of KGML compounds that participate in pathway KGML reactions.

    contigs_db_path : str, None
        Path to contigs database containing reaction network.

    network : anvio.reactionnetwork.GenomicNetwork, None
        Reaction network that can either be independent of a contigs database (contigs_db_path value
        of None) or associated with a contigs database.

    network_keggrn_id_to_modelseed_reactions : dict[str, list[rn.ModelSEEDReaction]], {}
        Map the IDs of KEGG reactions (not KGML reaction IDs) to aliased ModelSEED reactions in the
        reaction network.

    network_keggcpd_id_to_modelseed_compounds : dict[str, list[rn.ModelSEEDCompound]], {}
        Map the IDs of KEGG compounds (not KGML compound IDs) to aliased ModelSEED compounds in the
        reaction network.

    network_keggcpd_ids_in_pathway : list[str]
        KEGG compound IDs in the pathway from the reaction network. KEGG compounds in the reaction
        network are selected to ensure that they are linked to KO annotations by ModelSEED reactions
        with KEGG reaction aliases, not EC number aliases, since EC numbers associated with KOs can
        alias a large number of reactions of questionable validity for the enzyme.

    compound_fate : Literal['consume', 'produce', 'both'], 'consume'
        Seek chains that consume or produce compounds in the network. If 'consume' or 'produce',
        only consumption or production chains are sought, respectively. If 'both', both consumption
        and production chains are sought, yielding pairs of chains, one consumption chain ordered
        from the first reactant to the last product, and one production chain ordered from the last
        product to the first reactant. Therefore, chains that only contain reversible reactions are
        found in both directions with 'both', as the first reactant can be the last product, and the
        last product can be the first reactant, yielding four chains traversing identical compounds
        and reactions.

    max_reactions : int, None
        Truncate chains at this number of reactions. If None, chains can be continued to
        indeterminate length.

    keep_intermediate_chains : bool, False
        Chains starting from compounds in the network can be subchains of chains from other
        compounds in the network. If False, such intermediate chains are ignored.

    max_gaps : int, 0
        Chains can contain up to this number of reactions not genomically encoded in the reaction
        network.

    allow_terminal_gaps : bool, False
        Chains can start or end with reactions not found in the reaction network if True.

    allow_alternative_reaction_gaps : bool, False
        If a chain links two compounds by a reaction in the reaction network, and there are other
        "parallel" KGML reactions not in the reaction network that also link the compounds, then
        treat these parallel reactions as gaps with a value of True. With a value of False, ignore
        parallel reaction gaps.

    run : anvio.terminal.Run, anvio.terminal.Run()
        This object prints run information to the terminal.

    verbose : bool, False
        Print additional runtime information to the terminal, such as reaction network summary
        statistics.
    """
    def __init__(self, args: Namespace):
        """
        Parameters
        ==========
        args : argparse.Namespace
            Contains arguments, listed below. See the class docstring for more information on
            arguments set as attributes, including default values. The only required argument is
            kegg_pathway_number.

        'kegg_pathway_number' : str

        'contigs_db_path' : str

        'network' : anvio.reactionnetwork.GenomicNetwork

        'compound_fate' : Literal['consume', 'produce', 'both']

        'max_reactions' : int

        'keep_intermediate_chains' : bool

        'max_gaps' : int

        'allow_terminal_gaps' : bool

        'allow_alternative_reaction_gaps' : bool

        'run' : anvio.terminal.Run

        'verbose' : bool
        """
        A = lambda x, y: args.__dict__[x] if x in args.__dict__ else y

        self.kegg_pathway_number: str = args.kegg_pathway_number

        self.contigs_db_path: str = A('contigs_db_path', None)
        self.network: rn.GenomicNetwork = A('network', None)
        self.verbose = A('verbose', False)
        if self.contigs_db_path is not None and self.network is None:
            constructor = rn.Constructor()
            self.network = constructor.load_contigs_database_network(
                self.contigs_db_path, quiet=not self.verbose
            )
        self.compound_fate: str = A('compound_fate', 'consume')
        self.max_reactions: int = A('max_reactions', None)
        self.keep_intermediate_chains: bool = A('keep_intermediate_chains', False)
        self.max_gaps: int = A('max_gaps', 0)
        self.allow_terminal_gaps: bool = A('allow_terminal_gaps', False)
        self.allow_alternative_reaction_gaps: bool = A('allow_alternative_reaction_gaps', False)

        # Assume that the reaction network was constructed with the KEGG and ModelSEED databases
        # found at the default anvi'o location.
        self.kegg_data = rn.KEGGData()

        self.run: terminal.Run = A('run', terminal.Run())

        self.sanity_check()

        # Load the set of KGML files for the pathway.
        kgml_rn_dir = self.kegg_data.kegg_context.kgml_1x_rn_dir
        kgml_ko_dir = self.kegg_data.kegg_context.kgml_1x_ko_dir
        kgml_ec_dir = self.kegg_data.kegg_context.kgml_1x_ec_dir
        xml_ops = kgml.XMLOps()
        kgml_rn_path = os.path.join(kgml_rn_dir, f'rn{self.kegg_pathway_number}.xml')
        kgml_ko_path = os.path.join(kgml_ko_dir, f'ko{self.kegg_pathway_number}.xml')
        kgml_ec_path = os.path.join(kgml_ec_dir, f'ec{self.kegg_pathway_number}.xml')
        self.kgml_rn_pathway = xml_ops.load(kgml_rn_path)
        self.kgml_ko_pathway = xml_ops.load(kgml_ko_path)
        self.kgml_ec_pathway = xml_ops.load(kgml_ec_path)

        # Make attributes storing key KGML pathway data.
        ko_pathway_kgml_reaction_ids: list[str] = []
        for reaction_uuid in self.kgml_ko_pathway.children['reaction']:
            kgml_reaction: kgml.Reaction = self.kgml_ko_pathway.uuid_element_lookup[reaction_uuid]
            ko_pathway_kgml_reaction_ids.append(kgml_reaction.id)

        ec_pathway_kgml_reaction_ids: list[str] = []
        for reaction_uuid in self.kgml_ec_pathway.children['reaction']:
            kgml_reaction: kgml.Reaction = self.kgml_ec_pathway.uuid_element_lookup[reaction_uuid]
            ec_pathway_kgml_reaction_ids.append(kgml_reaction.id)

        self.rn_pathway_kgml_compound_id_to_kgml_reactions: dict[str, list[kgml.Reaction]] = {}
        self.pathway_nonenzymatic_kgml_reaction_ids: list[str] = []
        for reaction_uuid in self.kgml_rn_pathway.children['reaction']:
            kgml_reaction: kgml.Reaction = self.kgml_rn_pathway.uuid_element_lookup[reaction_uuid]

            for substrate_uuid in kgml_reaction.children['substrate']:
                substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[substrate_uuid]
                try:
                    self.rn_pathway_kgml_compound_id_to_kgml_reactions[substrate.id].append(
                        kgml_reaction
                    )
                except KeyError:
                    self.rn_pathway_kgml_compound_id_to_kgml_reactions[substrate.id] = [
                        kgml_reaction
                    ]
            for product_uuid in kgml_reaction.children['product']:
                product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[product_uuid]
                try:
                    self.rn_pathway_kgml_compound_id_to_kgml_reactions[product.id].append(
                        kgml_reaction
                    )
                except KeyError:
                    self.rn_pathway_kgml_compound_id_to_kgml_reactions[product.id] = [kgml_reaction]

            if kgml_reaction.id in ko_pathway_kgml_reaction_ids:
                continue
            if kgml_reaction.id in ec_pathway_kgml_reaction_ids:
                continue
            self.pathway_nonenzymatic_kgml_reaction_ids.append(kgml_reaction.id)

        self.rn_pathway_kgml_compound_id_to_kgml_compound_entry: dict[str, list[kgml.Entry]] = {}
        for entry_uuid in self.kgml_rn_pathway.children['entry']:
            kgml_entry: kgml.Entry = self.kgml_rn_pathway.uuid_element_lookup[entry_uuid]
            if kgml_entry.type != 'compound':
                continue
            self.rn_pathway_kgml_compound_id_to_kgml_compound_entry[kgml_entry.id] = kgml_entry

        self.ko_pathway_kgml_ortholog_id_to_kgml_ortholog_entry: dict[str, kgml.Entry] = {}
        for entry_uuid in self.kgml_ko_pathway.children['entry']:
            kgml_entry: kgml.Entry = self.kgml_ko_pathway.uuid_element_lookup[entry_uuid]
            if kgml_entry.type != 'ortholog':
                continue
            self.ko_pathway_kgml_ortholog_id_to_kgml_ortholog_entry[kgml_entry.id] = kgml_entry

        self.rn_pathway_kgml_reaction_id_to_ko_ids: dict[str, list[str]] = {}
        for reaction_uuid in self.kgml_rn_pathway.children['reaction']:
            kgml_reaction: kgml.Reaction = self.kgml_rn_pathway.uuid_element_lookup[reaction_uuid]
            try:
                kgml_ortholog_entry = self.ko_pathway_kgml_ortholog_id_to_kgml_ortholog_entry[
                    kgml_reaction.id
                ]
            except KeyError:
                continue
            ko_ids: list[str] = []
            for candidate_ko_id in kgml_ortholog_entry.name.split():
                if candidate_ko_id[:3] != 'ko:':
                    continue
                ko_ids.append(candidate_ko_id[3:])
            self.rn_pathway_kgml_reaction_id_to_ko_ids[kgml_reaction.id] = ko_ids

        self.rn_pathway_keggcpd_ids_in_kgml_reactions = \
            self._get_kgml_reaction_keggcpd_ids_in_pathway()

        # Make attributes storing key reaction network data.
        if self.network:
            self.network_keggrn_id_to_modelseed_reactions: dict[
                str, list[rn.ModelSEEDReaction]
            ] = {}
            for modelseed_reaction in self.network.reactions.values():
                for keggrn_id in modelseed_reaction.kegg_aliases:
                    try:
                        self.network_keggrn_id_to_modelseed_reactions[keggrn_id].append(
                            modelseed_reaction
                        )
                    except KeyError:
                        self.network_keggrn_id_to_modelseed_reactions[keggrn_id] = [
                            modelseed_reaction
                        ]

            self.network_keggcpd_id_to_modelseed_compounds: dict[
                str, list[rn.ModelSEEDCompound]
            ] = {}
            for modelseed_compound in self.network.metabolites.values():
                for keggcpd_id in modelseed_compound.kegg_aliases:
                    try:
                        self.network_keggcpd_id_to_modelseed_compounds[keggcpd_id].append(
                            modelseed_compound
                        )
                    except KeyError:
                        self.network_keggcpd_id_to_modelseed_compounds[keggcpd_id] = [
                            modelseed_compound
                        ]

            self.network_keggcpd_ids_in_pathway = \
                self._get_reaction_network_keggcpd_ids_in_pathway()
        else:
            self.network_keggrn_id_to_modelseed_reactions = None
            self.network_keggcpd_id_to_modelseed_compounds = None
            self.network_keggcpd_ids_in_pathway = None

    @staticmethod
    def check_pathway_number(kegg_pathway_number: str, kegg_data: rn.KEGGData) -> bool:
        """
        Check that the KEGG pathway number has an RN type KGML file available in the anvi'o KEGG
        installation.

        Parameters
        ==========
        kegg_pathway_number : str
            Numerical ID of a pathway.

        kegg_data : anvio.reactionnetwork.KEGGData
            Contains information on an anvi'o KEGG installation.

        Returns
        =======
        bool
            RN type KGML file exists if True.
        """
        kgml_rn_path = os.path.join(
            kegg_data.kegg_context.kgml_1x_rn_dir, f'rn{kegg_pathway_number}.xml'
        )
        if os.path.exists(kgml_rn_path):
            return True
        return False

    def sanity_check(self) -> None:
        """Check the validity of various attributes."""
        if not self.check_pathway_number(
            kegg_pathway_number=self.kegg_pathway_number, kegg_data=self.kegg_data
        ):
            raise ConfigError(
                "The KEGG pathway must have a reaction (RN) type KGML file available in the anvi'o "
                "KEGG installation. Valid pathways are in categories 1.0 - 1.11 as of the March 3, "
                "2025 release of KEGG (see https://www.genome.jp/kegg/pathway.html)."
            )

        if not (
            self.max_reactions is None or
            isinstance(self.max_reactions, int) and self.max_reactions > 0
        ):
            raise ConfigError("'max_reactions' must have a value of None or a positive int.")

        if not (isinstance(self.max_gaps, int) and self.max_gaps >= 0):
            raise ConfigError("'max_gaps' must have a non-negative int value.")

        if self.compound_fate not in ('consume', 'produce', 'both'):
            raise ConfigError(
                "'compound_fate' must have a value of 'consume', 'produce', or 'both'."
            )

    def _get_kgml_reaction_keggcpd_ids_in_pathway(self) -> list[str]:
        """
        Get KEGG compound IDs of KGML compounds that participate in pathway KGML reactions.

        Returns
        =======
        list[str]
            KEGG compound IDs.
        """
        rn_pathway_keggcpd_ids_in_kgml_reactions = []
        for reaction_uuid in self.kgml_rn_pathway.children['reaction']:
            kgml_reaction: kgml.Reaction = self.kgml_rn_pathway.uuid_element_lookup[reaction_uuid]
            for uuid in kgml_reaction.children['substrate']:
                substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[uuid]
                kgml_compound_entry = self.rn_pathway_kgml_compound_id_to_kgml_compound_entry[
                    substrate.id
                ]
                for candidate_keggcpd_id in kgml_compound_entry.name.split():
                    if candidate_keggcpd_id[:4] != 'cpd:':
                        continue
                    keggcpd_id = candidate_keggcpd_id[4:]
                    rn_pathway_keggcpd_ids_in_kgml_reactions.append(keggcpd_id)
        rn_pathway_keggcpd_ids_in_kgml_reactions = list(
            set(rn_pathway_keggcpd_ids_in_kgml_reactions)
        )

        return rn_pathway_keggcpd_ids_in_kgml_reactions

    def _get_reaction_network_keggcpd_ids_in_pathway(self) -> list[str]:
        """
        Get KEGG compound IDs in the pathway from the reaction network.

        KEGG compounds in the reaction network are selected to ensure that they are linked to KO
        annotations by ModelSEED reactions with KEGG reaction aliases, not EC number aliases, since
        EC numbers associated with KOs can alias a large number of reactions of questionable
        validity for the enzyme.

        Returns
        =======
        list[str]
            KEGG compound IDs.
        """
        keggcpd_ids = []
        pathway_id = f'map{self.kegg_pathway_number}'
        for ko in self.network.kos.values():
            if pathway_id not in ko.pathway_ids:
                continue
            for modelseed_reaction_id in ko.reaction_ids:
                modelseed_reaction = self.network.reactions[modelseed_reaction_id]
                if not modelseed_reaction.kegg_aliases:
                    continue
                for modelseed_compound_id in modelseed_reaction.compound_ids:
                    modelseed_compound = self.network.metabolites[modelseed_compound_id]
                    keggcpd_id_aliases = list(modelseed_compound.kegg_aliases)
                    keggcpd_ids += keggcpd_id_aliases
        keggcpd_ids = sorted(set(keggcpd_ids))

        return keggcpd_ids

    def get_chains(
        self,
        keggcpd_ids: Union[str, list[str]] = None,
        modelseed_compound_ids: Union[str, list[str]] = None
    ) -> dict[str, list[Chain]]:
        """
        Get chains in the pathway starting from compounds. Select compounds can be requested. If
        none are requested and a reaction network is available, chains are sought from all compounds
        in the network and pathway; if no network is available, chains are sought for all compounds
        in the pathway.

        Parameters
        ==========
        keggcpd_ids : Union[str, list[str]], None
            KEGG compound ID(s) to seek in the KGML file.

        modelseed_compound_ids : Union[str, list[str]], None
            ModelSEED compound ID(s), which are mapped to KEGG compound IDs via a reaction
            network.

        Returns
        =======
        dict[str, list[Chain]]
            Keys are KEGG compound IDs or, if requested, ModelSEED compound IDs. Values are lists of
            chains starting from the corresponding KGML compounds. If certain requested compounds
            are not found in the pathway, empty lists are returned for them.
        """
        if keggcpd_ids is not None and modelseed_compound_ids is not None:
            raise ConfigError(
                "Chains can be sought from either KEGG compound IDs or ModelSEED compound IDs, but "
                "not both."
            )

        if keggcpd_ids is None and modelseed_compound_ids is None:
            if self.network is None:
                keggcpd_ids = self.rn_pathway_keggcpd_ids_in_kgml_reactions
            else:
                keggcpd_ids = self.network_keggcpd_ids_in_pathway

        if isinstance(keggcpd_ids, str):
            keggcpd_ids = [keggcpd_ids]
        if isinstance(keggcpd_ids, list):
            compound_id_chains = self._get_chains_from_kegg_compound_ids(keggcpd_ids)
            return compound_id_chains

        if modelseed_compound_ids is not None and self.network is None:
            raise ConfigError(
                "A reaction network is required to get chains from ModelSEED compound IDs, but "
                "none happen to be stored in the 'network' attribute."
            )
        if isinstance(modelseed_compound_ids, str):
            modelseed_compound_ids = [modelseed_compound_ids]
        if isinstance(modelseed_compound_ids, list):
            run_verbosity = self.run.verbose
            if not modelseed_compound_ids:
                self.run.verbose = False
                modelseed_compound_ids = list(self.network.metabolites)
            else:
                self.run.verbose = True
            compound_id_chains = self._get_chains_from_modelseed_compound_ids(modelseed_compound_ids)
            self.run.verbose = run_verbosity

        return compound_id_chains

    def _get_chains_from_kegg_compound_ids(self, keggcpd_ids: list[str]) -> dict[str, list[Chain]]:
        """
        Get chains in the pathway starting from select KEGG compounds.

        Parameters
        ==========
        keggcpd_ids : list[str]
            KEGG compound IDs to seek in the KGML file.

        Returns
        =======
        dict[str, list[Chain]]
            Keys are the select KEGG compound IDs. Values are lists of chains starting from the
            corresponding KGML compounds. If certain of the KEGG compounds are not found in the
            pathway, empty lists are returned for them.
        """
        keggcpd_id_chains: dict[str, list[Chain]] = {}
        for keggcpd_id in keggcpd_ids:
            chains = self._get_chains_from_kegg_compound_id(keggcpd_id)
            keggcpd_id_chains[keggcpd_id] = chains

        if not self.keep_intermediate_chains:
            keggcpd_id_chains = self.remove_intermediate_chains(keggcpd_id_chains)

        return keggcpd_id_chains

    def _get_chains_from_kegg_compound_id(self, keggcpd_id: str) -> list[Chain]:
        """
        Get chains in the pathway starting from a KEGG compound.

        Parameters
        ==========
        keggcpd_id : str
            KEGG compound ID to seek in the KGML file.

        Returns
        =======
        list[Chain]
            List of chains starting from corresponding KGML compounds. If the KEGG compound is not
            found in the pathway, an empty list is returned.
        """
        kgml_compound_entries = self.kgml_rn_pathway.get_entries(kegg_ids=[keggcpd_id])

        chains: list[Chain] = []
        for kgml_compound_entry in kgml_compound_entries:
            kgml_compound_entry_chains: list[Chain] = self._get_chains_from_kgml_compound_entry(
                kgml_compound_entry
            )

            new_chains: list[Chain] = []
            for kgml_compound_entry_chain in kgml_compound_entry_chains:
                candidate_kgml_compound_entry_ids = [
                    c.id for c in kgml_compound_entry_chain.kgml_compound_entries
                ]
                for chain in chains:
                    kgml_compound_entry_ids = [c.id for c in chain.kgml_compound_entries]
                    if candidate_kgml_compound_entry_ids == kgml_compound_entry_ids:
                        break
                else:
                    new_chains.append(kgml_compound_entry_chain)
            chains = chains + new_chains

        return chains

    def _get_chains_from_modelseed_compound_ids(
        self,
        modelseed_compound_ids: list[str]
    ) -> dict[str, list[Chain]]:
        """
        Get chains in the pathway starting from select ModelSEED compounds.

        Parameters
        ==========
        modelseed_compound_ids : list[str]
            ModelSEED compound IDs, which are mapped to KEGG compound IDs via a reaction network and
            are sought in the KGML file.

        Returns
        =======
        dict[str, list[Chain]]
            Keys are the select ModelSEED compound IDs. Values are lists of chains starting from the
            corresponding KGML compounds. If certain of the ModelSEED compounds are not found in the
            pathway, empty lists are returned for them.
        """
        modelseed_compound_id_chains: dict[str, list[Chain]] = {}
        self.run.verbose = True if self.verbose else False
        for modelseed_compound_id in modelseed_compound_ids:
            chains = self._get_chains_from_modelseed_compound_id(modelseed_compound_id)
            modelseed_compound_id_chains[modelseed_compound_id] = chains

        if not self.keep_intermediate_chains:
            modelseed_compound_id_chains = self.remove_intermediate_chains(
                modelseed_compound_id_chains
            )

        return modelseed_compound_id_chains

    def _get_chains_from_modelseed_compound_id(self, modelseed_compound_id: str) -> list[Chain]:
        """
        Get chains in the pathway starting from a ModelSEED compound.

        Parameters
        ==========
        modelseed_compound_id : str
            ModelSEED compound ID, which is mapped to KEGG compound IDs and sought in the KGML file.

        Returns
        =======
        list[Chain]
            List of chains starting from corresponding KGML compounds. If the ModelSEED compound is
            not found in the pathway, an empty list is returned.
        """
        try:
            compound = self.network.metabolites[modelseed_compound_id]
        except KeyError:
            self.run.info_single(
                "There is no ModelSEED compound in the network with the ID, "
                f"'{modelseed_compound_id}'."
            )
            return []

        keggcpd_ids = compound.kegg_aliases
        if len(keggcpd_ids) == 0:
            self.run.info_single(
                f"ModelSEED compound '{modelseed_compound_id}' has no KEGG compound aliases."
            )
            return []
        keggcpd_ids = sorted(keggcpd_ids)

        chains: list[Chain] = []
        for keggcpd_id in keggcpd_ids:
            keggcpd_chains = self._get_chains_from_kegg_compound_id(keggcpd_id)

            new_chains: list[Chain] = []
            for keggcpd_chain in keggcpd_chains:
                candidate_kgml_compound_entry_ids = [
                    c.id for c in keggcpd_chain.kgml_compound_entries
                ]
                for chain in chains:
                    kgml_compound_entry_ids = [c.id for c in chain.kgml_compound_entries]
                    if candidate_kgml_compound_entry_ids == kgml_compound_entry_ids:
                        break
                else:
                    new_chains.append(keggcpd_chain)
            chains = chains + new_chains

        return chains

    def _get_chains_from_kgml_compound_entry(
        self,
        kgml_compound_entry: kgml.Entry,
        current_chain: Chain = None,
        terminal_chains: list[Chain] = None
    ) -> Union[Chain, list[Chain]]:
        """
        Recursive method to get chains in the pathway starting from a KGML compound entry (a circle
        on the map).

        Parameters
        ==========
        kgml_compound_entry : anvio.kgml.Entry
            KGML Entry of type "compound", which, in recursion, is at the end of a chain under
            construction.

        current_chain : Chain, None
            Chain under construction, which should be None in the initial method call.

        terminal_chains : list[Chain], None
            List of finished chains from the initial compound that fulfill
            criteria set in class attributes. This should be None in the initial method call.

        Returns
        =======
        Union[Chain, list[Chain]]
            Return the current chain under construction from a recursive method call. Return the
            list of terminal chains from the initial method call.
        """
        if terminal_chains is None:
            # This condition should only occur in the initial method call.
            terminal_chains = []

        # Check that the target compound KGML ID is involved in any reactions at all in the
        # pathway.
        try:
            kgml_reactions = self.rn_pathway_kgml_compound_id_to_kgml_reactions[
                kgml_compound_entry.id
            ]
        except KeyError:
            if terminal_chains:
                raise AssertionError(
                    "The target KGML compound ID has been found not to participate in KGML "
                    "reactions. This is only expected to occur for the initial method call, as "
                    "recursive calls target a compound known to participate in a reaction. "
                    "However, this error occurred because chains have been found, which should not "
                    "be possible at this early point in the initial call."
                )
            return terminal_chains

        if self.network:
            # Find aliased ModelSEED compounds.
            modelseed_compounds = {}
            for candidate_keggcpd_id in kgml_compound_entry.name.split():
                if candidate_keggcpd_id[:4] != 'cpd:':
                    continue
                keggcpd_id = candidate_keggcpd_id[4:]
                try:
                    keggcpd_modelseed_compounds = self.network_keggcpd_id_to_modelseed_compounds[
                        keggcpd_id
                    ]
                except KeyError:
                    continue
                for modelseed_compound in keggcpd_modelseed_compounds:
                    modelseed_compounds[modelseed_compound.modelseed_id] = modelseed_compound
            modelseed_compounds = tuple(
                sorted(modelseed_compounds.values(), key=lambda c: c.modelseed_id)
            )

        if current_chain is None:
            # This occurs in the initial method call. Seed a chain with the target compound.
            current_chain = Chain(
                kgml_compound_entries=[kgml_compound_entry],
                aliased_modelseed_compounds=[modelseed_compounds] if self.network else []
            )

            if self.compound_fate == 'consume':
                consumption_options = [True]
            elif self.compound_fate == 'produce':
                consumption_options = [False]
            elif self.compound_fate == 'both':
                consumption_options = [True, False]
            else:
                raise AssertionError("'compound_fate' does not have an accepted value.")

            # At the end of the initial method call, return terminal chains identified after
            # recursion.
            return_terminal_chains = True
        else:
            # This occurs in a recursive method call. Add the target compound to the chain under
            # construction.
            current_chain.kgml_compound_entries.append(kgml_compound_entry)
            if self.network:
                current_chain.aliased_modelseed_compounds.append(modelseed_compounds)

            if self.max_reactions is not None:
                if self.max_reactions == len(current_chain.kgml_reactions):
                    return current_chain

            # A cycle has been encountered, terminating the chain.
            if kgml_compound_entry.id in [c.id for c in current_chain.kgml_compound_entries[:-1]]:
                return current_chain

            consumption_options = [current_chain.is_consumed]

            # At the end of a recursive method call, return the completed chain, a candidate
            # terminal chain.
            return_terminal_chains = False

        for is_explore_consumption in consumption_options:
            # Record each KGML reaction involving the compound that is to be explored. If comparing
            # to a reaction network, also record whether the reaction represents a gap (if not in
            # the network), and the network KOs encoding the reaction.
            kgml_reaction_info: list[tuple[kgml.Reaction, bool, list[rn.KO]]] = []
            for kgml_reaction in kgml_reactions:
                if current_chain.kgml_reactions:
                    if kgml_reaction.id == current_chain.kgml_reactions[-1].id:
                        # Do not record the same reaction as the previous one in the chain.
                        continue

                if not self.network:
                    kgml_reaction_info.append((kgml_reaction, None, None))
                    continue

                if kgml_reaction.id in self.pathway_nonenzymatic_kgml_reaction_ids:
                    kgml_reaction_info.append((kgml_reaction, False, []))
                    continue

                try:
                    ko_ids = self.rn_pathway_kgml_reaction_id_to_ko_ids[kgml_reaction.id]
                except KeyError:
                    ko_ids = []

                network_kos: list[rn.KO] = []
                for ko_id in ko_ids:
                    try:
                        network_kos.append(self.network.kos[ko_id])
                    except KeyError:
                        pass

                if network_kos:
                    # The reaction is associated with one or more KOs in the reaction network.
                    kgml_reaction_info.append((kgml_reaction, False, network_kos))
                    continue

                # The reaction would be a gap in the reaction network.
                if sum(current_chain.gaps) >= self.max_gaps:
                    continue

                if not self.allow_terminal_gaps:
                    # Gaps are not allowed at the beginning or end of the chain.
                    if not current_chain.kgml_reactions:
                        # The reaction chain would start with a gap.
                        continue

                kgml_reaction_info.append((kgml_reaction, True, network_kos))

            for kgml_reaction, is_gap, network_kos in kgml_reaction_info:
                if not self.allow_alternative_reaction_gaps and is_gap:
                    # Ignore KGML reaction gaps that transform the same KGML compounds as other
                    # "parallel" KGML reactions that are in the reaction network.
                    substrate_ids: list[str] = []
                    for uuid in kgml_reaction.children['substrate']:
                        substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[uuid]
                        substrate_ids.append(substrate.id)
                    product_ids: list[str] = []
                    for uuid in kgml_reaction.children['product']:
                        product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[uuid]
                        product_ids.append(product.id)

                    for other_kgml_reaction, other_is_gap, _ in kgml_reaction_info:
                        if kgml_reaction.id == other_kgml_reaction.id:
                            continue
                        if other_is_gap:
                            continue
                        other_substrate_ids: list[str] = []
                        for uuid in other_kgml_reaction.children['substrate']:
                            substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[uuid]
                            other_substrate_ids.append(substrate.id)
                        other_product_ids: list[str] = []
                        for uuid in other_kgml_reaction.children['product']:
                            product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[uuid]
                            other_product_ids.append(product.id)
                        if substrate_ids == other_substrate_ids and product_ids == other_product_ids:
                            is_alternative_reaction = True
                            break
                    else:
                        is_alternative_reaction = False
                    if is_alternative_reaction:
                        continue

                # The compound under consideration is either recorded as a substrate or product of
                # the KGML reaction.
                for substrate_uuid in kgml_reaction.children['substrate']:
                    # Search for the compound under consideration among KGML reaction substrates.
                    kgml_substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[
                        substrate_uuid
                    ]
                    if kgml_compound_entry.id != kgml_substrate.id:
                        continue

                    if is_explore_consumption:
                        is_forward = True
                    else:
                        if kgml_reaction.type == 'reversible':
                            # Although the compound under consideration is recorded as a substrate
                            # of the KGML reaction and it needs to be a product, the reaction is
                            # reversible.
                            is_forward = False
                        else:
                            is_forward = None
                            break

                    next_kgml_compound_ids = []
                    for product_uuid in kgml_reaction.children['product']:
                        kgml_product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[
                            product_uuid
                        ]
                        next_kgml_compound_ids.append(kgml_product.id)
                    break
                else:
                    for product_uuid in kgml_reaction.children['product']:
                        # Search for the compound under consideration among KGML reaction products.
                        kgml_product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[
                            product_uuid
                        ]
                        if kgml_compound_entry.id != kgml_product.id:
                            continue

                        if not is_explore_consumption:
                            is_forward = True
                        else:
                            if kgml_reaction.type == 'reversible':
                                # Although the compound under consideration is recorded as a product
                                # of the KGML reaction and it needs to be a substrate, the reaction
                                # is reversible.
                                is_forward = False
                            else:
                                is_forward = None
                                break

                        next_kgml_compound_ids = []
                        for substrate_uuid in kgml_reaction.children['substrate']:
                            kgml_substrate: kgml.Substrate = \
                                self.kgml_rn_pathway.uuid_element_lookup[substrate_uuid]
                            next_kgml_compound_ids.append(kgml_substrate.id)
                        break
                    else:
                        raise AssertionError

                if is_forward is None:
                    continue

                if self.network:
                    # Get the ModelSEED reactions aliasing the KEGG reactions underlying the KGML
                    # reaction.
                    network_kos = tuple(sorted(network_kos, key=lambda ko: ko.id))

                    modelseed_reaction_dict: dict[str, rn.ModelSEEDReaction] = {}
                    for candidate_keggrn_id in kgml_reaction.name.split():
                        if candidate_keggrn_id[:3] != 'rn:':
                            continue
                        keggrn_id = candidate_keggrn_id[3:]
                        try:
                            modelseed_reactions = self.network_keggrn_id_to_modelseed_reactions[
                                keggrn_id
                            ]
                        except KeyError:
                            continue
                        for modelseed_reaction in modelseed_reactions:
                            modelseed_reaction_dict[
                                modelseed_reaction.modelseed_id
                            ] = modelseed_reaction
                    modelseed_reactions = [
                        modelseed_reaction_dict[modelseed_reaction_id]
                        for modelseed_reaction_id in sorted(modelseed_reaction_dict)
                    ]

                # Recurse on each KGML compound on the other side of the reaction.
                for next_kgml_compound_id in next_kgml_compound_ids:
                    if (
                        len(current_chain.kgml_compound_entries) >= 2 and
                        next_kgml_compound_id == current_chain.kgml_compound_entries[-2].id
                    ):
                        # Avoid backtracking to the previous KGML compound in the chain. This causes
                        # the shortest cycles with a single intermediate, such as the cyclic reuse
                        # of thiamine diphosphate (ThPP) in the Krebs cycle (00020), to not close in
                        # chains.
                        continue

                    new_chain = Chain(
                        kgml_compound_entries=current_chain.kgml_compound_entries.copy(),
                        is_consumed=is_explore_consumption,
                        kgml_reactions=current_chain.kgml_reactions + [kgml_reaction],
                        kgml_reaction_directions=\
                            current_chain.kgml_reaction_directions + [is_forward],
                        gaps=current_chain.gaps + [is_gap] if self.network else [],
                        aliased_modelseed_compounds=\
                            current_chain.aliased_modelseed_compounds.copy() \
                                if self.network else [],
                        network_kos=current_chain.network_kos + [tuple(network_kos)] \
                            if self.network else [],
                        aliased_modelseed_reactions=\
                            current_chain.aliased_modelseed_reactions + \
                                [tuple(modelseed_reactions)] if self.network else []
                    )
                    candidate_terminal_chain: Chain = \
                        self._get_chains_from_kgml_compound_entry(
                            self.rn_pathway_kgml_compound_id_to_kgml_compound_entry[
                                next_kgml_compound_id
                            ],
                            current_chain=new_chain,
                            terminal_chains=terminal_chains
                        )

                    if self.network:
                        if not self.allow_terminal_gaps:
                            # Gaps are not allowed at the beginning or end of the chain.
                            if candidate_terminal_chain.gaps[-1]:
                                # The reaction chain would end with a gap.
                                continue

                    # A subchain to each compound in a terminal chain is also generated as a
                    # "candidate terminal chain" and is ignored.
                    if terminal_chains:
                        candidate_kgml_compound_ids = [
                            c.id for c in candidate_terminal_chain.kgml_compound_entries
                        ]
                        terminal_kgml_compound_ids = [
                            c.id for c in terminal_chains[-1].kgml_compound_entries
                        ]
                        candidate_kgml_reaction_ids = [
                            r.id for r in candidate_terminal_chain.kgml_reactions
                        ]
                        terminal_kgml_reaction_ids = [
                            r.id for r in terminal_chains[-1].kgml_reactions
                        ]
                        if (
                            candidate_kgml_compound_ids == terminal_kgml_compound_ids[
                                : len(candidate_kgml_compound_ids)
                            ] and candidate_kgml_reaction_ids == terminal_kgml_reaction_ids[
                                : len(candidate_kgml_reaction_ids)
                            ]
                        ):
                            if len(candidate_kgml_compound_ids) == len(terminal_kgml_compound_ids):
                                raise AssertionError(
                                    "Only shorter subchains of the last found terminal chain were "
                                    "expected here, but a candidate terminal chain with the same "
                                    "KGML reactions and compounds was found."
                                )
                            continue

                    if candidate_terminal_chain.is_consumed:
                        candidate_terminal_chain.is_consumption_terminus = \
                            self.check_kgml_compound_entry_consumption_terminus(
                                candidate_terminal_chain.kgml_compound_entries[0]
                            )
                        candidate_terminal_chain.is_production_terminus = \
                            self.check_kgml_compound_entry_production_terminus(
                                candidate_terminal_chain.kgml_compound_entries[-1]
                            )
                    else:
                        candidate_terminal_chain.is_consumption_terminus = \
                            self.check_kgml_compound_entry_consumption_terminus(
                                candidate_terminal_chain.kgml_compound_entries[-1]
                            )
                        candidate_terminal_chain.is_production_terminus = \
                            self.check_kgml_compound_entry_production_terminus(
                                candidate_terminal_chain.kgml_compound_entries[0]
                            )

                    i = 0
                    for candidate_kgml_reaction in candidate_terminal_chain.kgml_reactions:
                        if candidate_kgml_reaction.type == 'irreversible':
                            break
                        i += 1
                    if candidate_terminal_chain.is_consumed:
                        candidate_terminal_chain.consumption_reversibility_range = (0, i + 1)
                    else:
                        candidate_terminal_chain.production_reversibility_range = (0, i + 1)

                    i = 0
                    for candidate_kgml_reaction in candidate_terminal_chain.kgml_reactions[::-1]:
                        if candidate_kgml_reaction.type == 'irreversible':
                            break
                        i += 1
                    if candidate_terminal_chain.is_consumed:
                        candidate_terminal_chain.production_reversibility_range = (
                            len(candidate_terminal_chain.kgml_reactions) - i,
                            len(candidate_terminal_chain.kgml_reactions) + 1
                        )
                    else:
                        candidate_terminal_chain.consumption_reversibility_range = (
                            len(candidate_terminal_chain.kgml_reactions) - i,
                            len(candidate_terminal_chain.kgml_reactions) + 1
                        )

                    candidate_terminal_chain.cyclic_branch_index = self.get_cyclic_branch_index(
                        candidate_terminal_chain
                    )

                    terminal_chains.append(candidate_terminal_chain)

        if return_terminal_chains:
            return terminal_chains
        else:
            return current_chain

    def check_kgml_compound_entry_consumption_terminus(
        self,
        kgml_compound_entry: kgml.Entry
    ) -> bool:
        """
        Check if a KGML Entry of type "compound" is a consumption terminus in a pathway map. A
        consumption terminus is defined in a negative sense as follows.

        The compound is not a consumption terminus if an irreversible reaction produces the
        compound. Also, the compound is not a consumption terminus if a) there are multiple
        reactions that consume the compound, b) any one of these reactions is reversible, and c) not
        all of the reactions have the same KGML compound products.

        Examples of consumption termini are acetyl-CoA in Fatty Acid Biosynthesis (00061) and
        L-valine, L-leucine, and L-isoleucine in Valine, Leucine, and Isoleucine Degradation
        (00280).

        Parameters
        ==========
        kgml_compound_entry : anvio.kgml.Entry
            Compound to check.

        Returns
        =======
        bool
            True if consumption terminus, False otherwise.
        """
        try:
            kgml_reactions = self.rn_pathway_kgml_compound_id_to_kgml_reactions[
                kgml_compound_entry.id
            ]
        except KeyError:
            return

        # Inspect each reaction involving the compound.
        reaction_reversibility_data: list[bool] = []
        reaction_compound_data: list[tuple[str]] = []
        for kgml_reaction in kgml_reactions:
            is_reversible = kgml_reaction.type == 'reversible'

            kgml_product_ids = []
            for product_uuid in kgml_reaction.children['product']:
                kgml_product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[product_uuid]
                if kgml_compound_entry.id == kgml_product.id and not is_reversible:
                    # The compound is the product of an irreversible reaction.
                    return False
                kgml_product_ids.append(kgml_product.id)

            kgml_substrate_ids = []
            for substrate_uuid in kgml_reaction.children['substrate']:
                kgml_substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[
                    substrate_uuid
                ]
                kgml_substrate_ids.append(kgml_substrate.id)

            reaction_reversibility_data.append(is_reversible)
            if kgml_compound_entry.id in kgml_substrate_ids:
                reaction_compound_data.append(tuple(kgml_product_ids))
            elif kgml_compound_entry.id in kgml_product_ids:
                reaction_compound_data.append(tuple(kgml_substrate_ids))
            else:
                raise AssertionError(
                    "The compound should either be a substrate or product of the reaction. This "
                    "point should not be reached."
                )

        if any(reaction_reversibility_data) and len(set(reaction_compound_data)) > 1:
            # The compound is involved in at least one reversible reaction. The compound is also
            # involved in distinct reactions with different sets of chemical species.
            return False

        return True

    def check_kgml_compound_entry_production_terminus(
        self,
        kgml_compound_entry: kgml.Entry
    ) -> bool:
        """
        Check if a KGML Entry of type "compound" is a production terminus in a pathway map. A
        production terminus is defined in a negative sense as follows.

        The compound is not a production terminus if an irreversible reaction consumes the
        compound. Also, the compound is not a production terminus if a) there are multiple
        reactions that produce the compound, b) any one of these reactions is reversible, and c) not
        all of the reactions have the same KGML compound reactants.

        Examples of production termini are acetyl-CoA in Fatty Acid Degradation (00071) and
        L-valine, L-leucine, and L-isoleucine in Valine, Leucine, and Isoleucine Biosynthesis
        (00290).

        Parameters
        ==========
        kgml_compound_entry : anvio.kgml.Entry
            Compound to check.

        Returns
        =======
        bool
            True if production terminus, False otherwise.
        """
        try:
            kgml_reactions = self.rn_pathway_kgml_compound_id_to_kgml_reactions[
                kgml_compound_entry.id
            ]
        except KeyError:
            return

        reaction_reversibility_data: list[bool] = []
        reaction_compound_data: list[tuple[str]] = []
        for kgml_reaction in kgml_reactions:
            is_reversible = kgml_reaction.type == 'reversible'

            kgml_substrate_ids = []
            for substrate_uuid in kgml_reaction.children['substrate']:
                kgml_substrate: kgml.Substrate = self.kgml_rn_pathway.uuid_element_lookup[
                    substrate_uuid
                ]
                if kgml_compound_entry.id == kgml_substrate.id and not is_reversible:
                    return False
                kgml_substrate_ids.append(kgml_substrate.id)

            kgml_product_ids = []
            for product_uuid in kgml_reaction.children['product']:
                kgml_product: kgml.Product = self.kgml_rn_pathway.uuid_element_lookup[product_uuid]
                kgml_product_ids.append(kgml_product.id)

            reaction_reversibility_data.append(is_reversible)
            if kgml_compound_entry.id in kgml_substrate_ids:
                reaction_compound_data.append(tuple(kgml_product_ids))
            elif kgml_compound_entry.id in kgml_product_ids:
                reaction_compound_data.append(tuple(kgml_substrate_ids))
            else:
                raise AssertionError

        if any(reaction_reversibility_data) and len(set(reaction_compound_data)) > 1:
            return False

        return True

    @staticmethod
    def get_cyclic_branch_index(chain: Chain) -> int:
        """
        To identify a "partly cyclic" chain looping a cycle via the entry or exit of an uncycled
        compound through a reaction in the cycle, check if the last reaction in the chain is
        traversed earler in the chain.

        This method avoids cyclic chains that do not branch off and only include cycled compounds.

        Stringent criteria are imposed in the identification of a partly cyclic chain beyond the
        last KGML reaction occurring a second time in the chain, although that might be sufficient
        given experience. The branch KGML reaction must occur in the same direction and produce or
        consume the same KGML compound, which is the cycle-closing compound last in the chain.

        Parameters
        ==========
        chain : Chain
            Chain to be checked.

        Returns
        =======
        int
            If the chain is partly cyclic, a non-negative int representing the index of the first
            occurrence of the branching KGML reaction that also ends the chain. -1 if the chain is
            not partly cyclic.
        """
        last_kgml_reaction_id = chain.kgml_reactions[-1].id
        last_direction = chain.kgml_reaction_directions[-1]
        candidate_cycled_kgml_compound_id = chain.kgml_compound_entries[-1].id

        for i, (kgml_reaction, direction) in enumerate(zip(
            chain.kgml_reactions[: -1], chain.kgml_reaction_directions[: -1]
        )):
            if kgml_reaction.id == last_kgml_reaction_id and direction == last_direction:
                break
        else:
            return -1

        if (
            chain.is_consumed and
            candidate_cycled_kgml_compound_id == chain.kgml_compound_entries[i + 1].id
        ):
            # The partly cyclic consumption chain traversed the same KGML reaction in the same
            # direction producing the same compound as before.
            return i
        elif (
            not chain.is_consumed and
            candidate_cycled_kgml_compound_id == chain.kgml_compound_entries[i].id
        ):
            # The partly cyclic production chain traversed the same KGML reaction in the same
            # direction consuming the same reactant as before.
            return i

        return -1

    def remove_intermediate_chains(
        self,
        compound_id_chains: dict[str, list[Chain]]
    ) -> dict[str, list[Chain]]:
        """
        Remove subchains, or chains running in the same direction within other chains.

        Assume that there are no subchains of other chains starting from the same compound, since
        these are not found by methods that get chains (via _get_chains_from_kgml_compound_entry):
        subchains can only start from an intermediate compound in another chain.

        Parameters
        ==========
        compound_id_chains : dict[str, list[Chain]]
            Keys are IDs of compounds from which chains start. Values are lists of chains from the
            compounds.

        Returns
        =======
        dict[str, list[Chain]]
            Dictionary like compound_id_chains but with subchains removed.
        """
        # Determine whether each chain is a subchain.
        derep_compound_id_chains: dict[str, list[Chain]] = {}
        for compound_id, chains in compound_id_chains.items():
            derep_chains: list[Chain] = []

            for chain in chains:
                kgml_compound_ids: list[str] = [c.id for c in chain.kgml_compound_entries]
                kgml_reaction_ids: list[str] = [r.id for r in chain.kgml_reactions]

                for other_compound_id, other_chains in compound_id_chains.items():
                    if compound_id == other_compound_id:
                        # Assume that there are no subchains starting from the same compound as
                        # other chains.
                        continue

                    for other_chain in other_chains:
                        if chain.is_consumed != other_chain.is_consumed:
                            # Subchains run in the same direction as the other chain.
                            continue

                        other_kgml_compound_ids: list[str] = [
                            c.id for c in other_chain.kgml_compound_entries
                        ]
                        delta_length = len(other_kgml_compound_ids) - len(kgml_compound_ids)
                        if delta_length < 1:
                            # Subchains must be shorter than the other chain.
                            continue

                        try:
                            i = other_kgml_compound_ids.index(kgml_compound_ids[0])
                        except ValueError:
                            # The chain isn't a subchain as its starting compound isn't in the other
                            # chain.
                            continue

                        other_kgml_reaction_ids: list[str] = [
                            r.id for r in other_chain.kgml_reactions
                        ]
                        if (
                            kgml_compound_ids == other_kgml_compound_ids[
                                i: i + len(kgml_compound_ids)
                            ] and kgml_reaction_ids == other_kgml_reaction_ids[
                                i: i + len(kgml_reaction_ids)
                            ]
                        ):
                            # The chain is a subchain because it shares the same compounds and
                            # reactions in the same order as the other chain.
                            break
                    else:
                        continue
                    # Subchain was found.
                    break
                else:
                    # The chain was not a subchain, so record it.
                    derep_chains.append(chain)

            if derep_chains:
                derep_compound_id_chains[compound_id] = derep_chains
        return derep_compound_id_chains

@dataclass
class GapChainRelations:
    """
    Records how chains with fewer (including 0) genomic gaps, "ungappy" chains, relate to a chain
    with more gaps, the "gappy" chain.

    Attributes
    ==========
    gappy_chain : Chain, None
        Chain with more gaps.

    ungappy_chains : list[Chain], []
        Chains with fewer (including 0) gaps.

    overlaps : list[tuple[tuple[int, int]]], []
        This list has a tuple item for each ungappy chain. There is an inner tuple for each reaction
        step shared between the gappy and ungappy chain, with the first item of the inner tuple
        being the index of the reaction in the gappy chain and the second item being the index of
        the reaction in the ungappy chain. Overlap involves both chains consuming the same reactant
        and producing the same product in the reaction step.

    is_subchain : list[bool], []
        This list has an item for each ungappy chain. An ungappy chain has a value of True if the
        chain is a subchain of the gappy chain, meaning that none of the gaps in the gappy chain
        occur in the ungappy subchain. Otherwise the ungappy chain has a value of False.
    """
    gappy_chain: Chain = None
    ungappy_chains: list[Chain] = field(default_factory=list)
    overlaps: list[tuple[tuple[int, int]]] = field(default_factory=list)
    is_subchain: list[bool] = field(default_factory=list)

@dataclass
class SharedGaps:
    """
    Information associated with the set of genomic gaps that exists in one or more gappy chains:
    different chains may share the same set of gaps.

    Attributes
    ==========
    gap_kgml_reactions : list[anvio.kgml.Reaction], []
        The reaction for each gap.

    gap_chain_relations : list[GapChainRelations], []
        This list contains an item per gappy chain that has the set of gaps represented here.
    """
    gap_kgml_reactions: list[kgml.Reaction] = field(default_factory=list)
    gap_chain_relations: list[GapChainRelations] = field(default_factory=list)

class GapAnalyzer:
    """
    Analyze chains of KGML compounds linked by reactions, with some reactions designated as genomic
    gaps. Compare two sets of chains found from the same KGML source but with the set of gappy
    chains permitting more gaps than ungappy chains.

    Attributes
    ==========
    gappy_chains : list[Chain]
        Chains with more gaps permitted than ungappy chains.

    ungappy_chains : list[Chain]
        Chains with fewer gaps permitted than gappy chains.

    gap_relations : dict[tuple[str], SharedGaps]
        Information associated with sets of gaps found in one or more chains, including
        relationships between ungappy and gappy chains.
    """
    def __init__(self, gappy_chains: list[Chain], ungappy_chains: list[Chain]):
        """
        Parameters
        ==========
        gappy_chains : list[Chain]
            Set as attribute of same name.

        ungappy_chains : list[Chain]
            Set as attribute of same name.
        """
        self.gappy_chains = gappy_chains
        self.ungappy_chains = ungappy_chains
        self.gap_relations = self._get_gap_relations()

    def _get_gap_relations(self) -> dict[tuple[str], SharedGaps]:
        """
        Get information associated with each set of gaps that exists in one or more gappy chains,
        including relationships between gappy chains and overlapping ungappy chains.

        Returns
        =======
        dict[tuple[str], SharedGaps]
            Keys are the KGML reaction IDs of sets of gaps in gappy chains. Values are information
            associated with each set of gaps.
        """
        def sort_ungappy_chains(
            ungappy_chains: list[Chain],
            overlaps: list[tuple[tuple[int, int]]]
        ) -> tuple[list[Chain], list[tuple[tuple[int, int]]]]:
            """
            Sort ungappy chains overlapping with a particular gappy chain by the index of the first
            overlapping reaction step in the gappy chain.

            Parameters
            ==========
            ungappy_chains : list[Chain]
                Ungappy chains overlapping with a gappy chain.

            overlaps : list[tuple[tuple[int, int]]]
                Each outer tuple corresponds to each ungappy chain. Inner tuples represent KGML
                reactions shared, in order, between the gappy and ungappy chains. The items of the
                tuple are the reaction indices in the gappy and ungappy chains, respectively.

            Returns
            =======
            list[Chain]
                Sorted ungappy chains.

            list[tuple[int, int]]
                Sorted overlaps, corresponding to sorted ungappy chains.
            """
            ungappy_keys: dict[int, int] = {}
            for i, (ungappy_chain, overlap) in enumerate(zip(ungappy_chains, overlaps)):
                ungappy_keys[i] = overlap[0][0]
            sorted_ungappy_chain_indices = sorted(
                ungappy_keys, key=lambda ungappy_chain_index: ungappy_keys[ungappy_chain_index]
            )
            sorted_ungappy_chains = []
            sorted_overlaps = []
            for i in sorted_ungappy_chain_indices:
                sorted_ungappy_chains.append(ungappy_chains[i])
                sorted_overlaps.append(overlaps[i])

            return sorted_ungappy_chains, sorted_overlaps

        gap_relations = {}
        for gappy_chain in self.gappy_chains:
            if not any(gappy_chain.gaps):
                # The chain has no gaps. (Chains with and without gaps can be returned when seeking
                # chains permitting gaps.)
                continue

            gap_kgml_reactions = [
                kgml_reaction for is_gap, kgml_reaction in
                zip(gappy_chain.gaps, gappy_chain.kgml_reactions) if is_gap
            ]
            gap_kgml_reaction_ids = tuple(
                [kgml_reaction.id for kgml_reaction in gap_kgml_reactions]
            )
            try:
                # Another gappy chain had the same gaps as the current chain.
                shared_gaps = gap_relations[gap_kgml_reaction_ids]
            except KeyError:
                gap_relations[gap_kgml_reaction_ids] = shared_gaps = SharedGaps()
            shared_gaps.gap_kgml_reactions = gap_kgml_reactions
            gap_chain_relations = GapChainRelations(gappy_chain=gappy_chain)
            shared_gaps.gap_chain_relations.append(gap_chain_relations)

            gappy_kgml_reaction_ids = [kr.id for kr in gappy_chain.kgml_reactions]
            gappy_kgml_compound_ids = [kc.id for kc in gappy_chain.kgml_compound_entries]

            # Record information on ungappy chains overlapping with the gappy chain.
            overlaps: list[tuple[tuple[int, int]]] = []
            ungappy_chains: list[Chain] = []
            for ungappy_chain in self.ungappy_chains:
                if gappy_chain.is_consumed != ungappy_chain.is_consumed:
                    # Ignore ungappy chains with reactions running in the opposite direction to the
                    # gappy chain.
                    continue

                ungappy_kgml_reaction_ids = [kr.id for kr in ungappy_chain.kgml_reactions]

                if gappy_kgml_reaction_ids == ungappy_kgml_reaction_ids:
                    # Ignore gappy and ungappy chains with identical reactions. The same chains can
                    # be returned when seeking chains allowing for more and fewer gaps. Chains can
                    # also have the same sequence of reactions but differ in which one of multiple
                    # initial reactants or final products of the reaction sequence are included in
                    # each chain.
                    continue

                ungappy_kgml_compound_ids = [kc.id for kc in ungappy_chain.kgml_compound_entries]

                overlap: list[tuple[int, int]] = []
                prev_ungappy_overlap_index = -1
                for gappy_index, gappy_kgml_reaction_id in enumerate(gappy_kgml_reaction_ids):
                    for i, ungappy_kgml_reaction_id in enumerate(
                        ungappy_kgml_reaction_ids[prev_ungappy_overlap_index + 1: ]
                    ):
                        ungappy_index = prev_ungappy_overlap_index + i + 1
                        if (
                            gappy_kgml_reaction_id == ungappy_kgml_reaction_id and
                            gappy_kgml_compound_ids[gappy_index] ==
                            ungappy_kgml_compound_ids[ungappy_index] and
                            gappy_kgml_compound_ids[gappy_index + 1] ==
                            ungappy_kgml_compound_ids[ungappy_index + 1]
                        ):
                            # Record the index of the reaction in the gappy and ungappy chains,
                            # respectively.
                            overlap.append((gappy_index, ungappy_index))
                            # Backtracking to previous ungappy chain reactions is not allowed.
                            prev_ungappy_overlap_index = ungappy_index
                            break

                if not overlap:
                    # Ignore ungappy chains that do not overlap the gappy chain.
                    continue
                ungappy_chains.append(ungappy_chain)
                # Represent each ungappy chain's overlap of the gappy chain as a tuple of tuples.
                overlaps.append(tuple(overlap))

            # Sort ungappy chains associated with the gappy chain.
            sorted_ungappy_chains, sorted_overlaps = sort_ungappy_chains(ungappy_chains, overlaps)
            gap_chain_relations.ungappy_chains = sorted_ungappy_chains
            gap_chain_relations.overlaps = sorted_overlaps

            # Record whether the ungappy chain is a subchain of the gappy chain.
            for ungappy_chain, overlap in zip(
                gap_chain_relations.ungappy_chains, gap_chain_relations.overlaps
            ):
                if len(overlap) == len(ungappy_chain.kgml_reactions):
                    gap_chain_relations.is_subchain.append(True)
                else:
                    gap_chain_relations.is_subchain.append(False)
        return gap_relations

    def rank_gaps(self) -> list[tuple[str]]:
        """
        Rank sets of gaps that are in the gappy chains but not the ungappy chains. This algorithm
        can be used to find candidate reactions to "gap-fill" in pathways. In essence, the algorithm
        assigns higher ranks to sets of gaps that occur in the middle of longer chains than toward
        the edges of shorter chains. If the gappy chains are allowed at most one gap and the ungappy
        chains are allowed no gaps, then the ranked "sets of gaps" are single reaction gaps.

        Here are the steps in the algorithm:

        1. Loop through each set of gaps (one or more gaps) in the gappy chains. (Gappy chains may
        share gap KGML reactions that branch to multiple KGML compounds, yielding multiple chains
        sharing the same gaps.)

        2. Inner loop through each gappy chain containing the set of gaps.

        3. Ignore gappy chains that, beside gap reactions, contain the same reactions as an ungappy
        chain. Gaps in the gappy chain create "shortcuts" in what is otherwise the same chain.

        4. Find the longest ungappy chains that encompass the reactions of the gappy chain beside
        its gaps not in ungappy chains. Call these unique segments.

        5. Sort unique segments in ascending order of length, with ties broken by index position in
        the chain. This is the last step in the inner loop.

        6. Sort gappy chains sharing the same gap reactions (gappy chains that only differ due to
        one or more gap reactions involving multiple substrates or products that are in the
        different chains). Gappy chains are sorted in descending order of unique segment length,
        first considering the shortest unique segment from each chain, then the next shortest to
        break ties, etc. The top-ranking chain has the longest of the shortest unique segments.

        7. Select the top-ranking gappy chain to represent the set of gaps. This is the last step of
        the outer loop.

        8. Sort the gappy chains selected in step (7) to represent each set of gaps. Use the
        ranking procedure of step (6).

        Returns
        =======
        list[tuple[str]]
            Ranked sets of gaps in gappy chains, with higher ranks first. Each tuple of KGML
            reaction IDs corresponds to a key of the gap_relations attribute.
        """
        def is_subsequence(t1: tuple[int], t2: tuple[int]) -> bool:
            """
            Check if the first sequence is a subsequence of the second.

            Parameters
            ==========
            t1 : tuple[int]
                First integer sequence.

            t2 : tuple[int]
                Second integer sequence.

            Returns
            =======
            bool
                True if a subsequence was found, and False otherwise.
            """
            if len(t1) >= len(t2):
                return False

            for i in range(len(t2) - len(t1) + 1):
                if t2[i: i + len(t1)] == t1:
                    return True
            return False

        def rank_gappy_chains_by_segment_lengths(
            gappy_chain_segments: dict[Any, list[int]]
        ) -> list[Any]:
            """
            Sort gappy chains in descending order of unique segment length, first considering the
            shortest segment from each chain, then the next shortest to break ties, etc. The
            top-ranking chain has the longest of the shortest segments.

            Parameters
            ==========
            gappy_chain_segments : dict[Any, list[int]]
                Keys identify a gappy chain. Values are lengths, in ascending order, of the unique
                segments in the gappy chain.

            Returns
            =======
            list[Any]
                Ranked keys of gappy_chain_segments.
            """
            # Find the maximum number of segments in the gappy chains.
            max_segment_count = 0
            for segments in gappy_chain_segments.values():
                if len(segments) > max_segment_count:
                    max_segment_count = len(segments)

            # Gappy chains with fewer segments than the maximum have their lists of segments lengths
            # padded at the end with zero lengths.
            gappy_chain_segment_lengths: dict[int, list[int]] = {}
            for gappy_chain_index, segments in gappy_chain_segments.items():
                padded_segments = segments + [
                    tuple() for i in range(max_segment_count - len(segments))
                ]
                gappy_chain_segment_lengths[gappy_chain_index] = [
                    len(segment) for segment in padded_segments
                ]

            ranked_gappy_chain_ids = [item[0] for item in sorted(
                gappy_chain_segment_lengths.items(),
                key=lambda item: tuple(-segment_length for segment_length in item[1])
            )]
            return ranked_gappy_chain_ids

        gap_unique_segments: dict[tuple[str], list[tuple[int]]] = {}
        for gap_kgml_reaction_ids, shared_gaps in self.gap_relations.items():
            gappy_chain_unique_segments: dict[int, list[tuple[int]]] = {}

            for gappy_chain_index, gap_chain_relations in enumerate(
                shared_gaps.gap_chain_relations
            ):
                gappy_chain = gap_chain_relations.gappy_chain

                # Ignore gappy chains that, beside gap reactions, consist of a smaller subset of
                # ungappy chain reactions. This is caused by gaps in the gappy chain creating
                # "shortcuts" in the ungappy chain.
                kgml_reaction_ids_absent_gaps = [
                    kgml_reaction.id for kgml_reaction in gappy_chain.kgml_reactions
                    if kgml_reaction.id not in gap_kgml_reaction_ids
                ]
                for ungappy_chain, overlap in zip(
                    gap_chain_relations.ungappy_chains, gap_chain_relations.overlaps
                ):
                    if len(overlap) == len(kgml_reaction_ids_absent_gaps):
                        assert (
                            len(ungappy_chain.kgml_reactions) > len(kgml_reaction_ids_absent_gaps)
                        )
                        is_shortcut = True
                        break
                else:
                    is_shortcut = False
                if is_shortcut:
                    continue

                # Find "unique segments", the longest ungappy chains that between them contain the
                # compounds and reactions of the gappy chain, besides the gap reactions not in
                # ungappy chains.

                # The case of cyclic gappy chains with the entry or exit of an uncycled compound
                # requires special treatment and creates an exception to the compound overlap
                # requirement. Consider a cycle of just two reactions. Reaction 1 incorporates A, a
                # compound produced by the cycle, A + B -> C. Reaction 2 emits D, a compound
                # produced by the cycle, C -> B + D.

                # First consider reaction 2 as the gap in the gappy chain, G, and reaction 1 as the
                # only reaction in the two related ungappy chains, U and V. Chain U contains
                # compounds A and C linked by reaction 1, and chain V contains compounds B and C
                # linked by reaction 1. Gappy chain G contains compounds A, C, B, and then C again:
                # reaction 1 reacts A to C, gap reaction 2 reacts C to B, and reaction 1 occurs
                # again, reacting B to C. The cycle is complete encountering C again. The gappy
                # chain is a consumption chain, found by the KGML network walker parameterized with
                # a compound fate value of 'consume' or 'both', but not 'produce'.

                # Both chains U and V fulfill the "unique segment" criterion of containing the
                # reactions that are not unique gaps in the gappy chain: both contain reaction 1.
                # However, neither fulfills the other criterion of containing the gappy chain
                # compounds, as U contains A but not B, and V contains B but not A. It does not make
                # sense to treat both as unique segments to fulfill the criterion, as this would
                # doubly trace reaction 1. Instead, select the ungappy chain that contains the
                # uncycled compound (here A). Chain U instead of V is chosen as the single unique segment
                # covering chain G.

                # For the sake of completeness, here is a complementary example in which reaction 1
                # is now the gap. Chain U contains compounds C and D linked by reaction 2, and chain
                # V contains compounds C and B linked by reaction 2. Gappy chain G contains
                # compounds C, B, C, and then D: reaction 2 reacts C to B, gap reaction 1 reacts B
                # to C, and reaction 2 occurs again, reacting C to D. The gappy chain is a
                # production chain, found by the KGML network walker parameterized with a compound
                # fate value of 'produce' or 'both', but not 'consume': therefore the compounds and
                # reactions are recorded in the opposite order just given, as the chain was found
                # starting from the final product, D, and completed encountered C a second time.

                # Both chains U and V fulfill the "unique segment" criterion of containing the
                # reactions that are not unique gaps in the gappy chain: both contain reaction 2.
                # However, neither fulfills the other criterion of containing the gappy chain
                # compounds, as U contains D but not B, and V contains B but not D. Follow the rule
                # of selecting the ungappy chain that contains the uncycled compound (here D). Chain
                # U instead of V is chosen as the single unique segment covering chain G.

                # To identify a chain looping a cycle via the entry or exit of an uncycled compound,
                # check if the last reaction in the chain is traversed earlier. In the first example
                # above, the last reaction of the gappy consumption chain, reaction 1, produces
                # compound C; reaction 1 is encountered earlier entering the cycle, also producing
                # C. In the second example, the last reaction of the gappy production chain,
                # reaction 2, consumes compound C; reaction 2 is encountered earlier exiting the
                # cycle, also consuming C. This method avoids cyclic chains that do not branch off
                # and only include cycled compounds. Unique segments can be found for such purely
                # cyclic gappy chains using the standard method.
                branch_index = gappy_chain.cyclic_branch_index
                if branch_index == -1:
                    # The gappy chain is not partly cyclic.
                    segments = [
                        tuple([overlap_indices[0] for overlap_indices in overlap])
                        for overlap in gap_chain_relations.overlaps
                    ]
                else:
                    # The gappy chain is partly cyclic. The branch index records the position in the
                    # chain of the particular occurrence of the reaction that enters or exits the
                    # cycle.
                    gappy_kgml_compound_id = gappy_chain.kgml_compound_entries[branch_index].id
                    segments: list[tuple[int]] = []
                    for ungappy_chain, overlap in zip(
                        gap_chain_relations.ungappy_chains, gap_chain_relations.overlaps
                    ):

                        gappy_overlap_indices = [overlap_indices[0] for overlap_indices in overlap]
                        if branch_index not in gappy_overlap_indices:
                            # The ungappy chain does not overlap with the reaction step that enters
                            # or exits the cycle, so ignore it.
                            break
                        segments.append(tuple(gappy_overlap_indices))

                # Find unique segments: the longest ungappy chains that encompass the reactions
                # of the gappy chain, besides gap reactions not in ungappy chains, and the
                # involved compounds.
                unique_segments: list[tuple[int]] = []
                for i, segment in enumerate(segments):
                    for j, other_segment in enumerate(segments):
                        if i == j:
                            continue
                        if is_subsequence(segment, other_segment):
                            break
                    else:
                        unique_segments.append(segment)
                # Sort unique segments in ascending order of length, with ties broken by index
                # position in the chain.
                unique_segments = sorted(set(unique_segments), key=lambda segment: len(segment))
                gappy_chain_unique_segments[gappy_chain_index] = unique_segments

                if not gappy_chain_unique_segments:
                    raise AssertionError(
                        "The gappy chain would only lack unique segments if the gappy and ungappy "
                        "chains do not have the proper relationship, in which they were found "
                        "identically except that more gaps were allowed in the gappy chains."
                    )

            if not gappy_chain_unique_segments:
                # All of the gappy chains containing the set of gaps were shortcuts in ungappy
                # chains. The set of gaps is thus ignored.
                continue

            # Sort gappy chains sharing the same gap reactions.
            ranked_gappy_chain_indices: list[int] = rank_gappy_chains_by_segment_lengths(
                gappy_chain_unique_segments
            )
            # Select the top-ranking gappy chain to represent the set of gaps.
            gap_unique_segments[gap_kgml_reaction_ids] = gappy_chain_unique_segments[
                ranked_gappy_chain_indices[0]
            ]

        # Sort the gappy chains selected to represent each set of gaps.
        ranked_gap_kgml_reaction_ids = rank_gappy_chains_by_segment_lengths(gap_unique_segments)
        return ranked_gap_kgml_reaction_ids

class GapFiller:
    """
    Analyze how genomic gaps in a KEGG pathway can be filled.

    Attributes
    ==========
    kegg_pathway_number : str
        Numerical ID of the pathway to analyze. The pathway must have a KEGG reaction (RN) type KGML
        file available. Valid pathways are in categories 1.0 - 1.11 as of the March 3, 2025 release
        of KEGG (see https://www.genome.jp/kegg/pathway.html).

    contigs_db_path : str
        Path to contigs database containing reaction network.

    all_ko_hits_path : str
        Path to table of all KO hits to genes in the contigs database. This file is generated by
        anvi-run-kegg-kofams with --log-bitscores.

    ko_cog_path : str
        Path to KEGG binary relations file mapping KO to COG IDs.

    compound_fate : Literal['consume', 'produce'], 'consume'
        Fill gaps in chains that consume or produce compounds in the network. The two options simply
        consider the same chains in reverse, with 'consume' treating chains ordered from reactant to
        product, and 'produce' treating chains ordered from product to reactant.

    max_reactions : int, None
        Truncate chains at this number of reactions. If None, chains can be continued to
        indeterminate length.

    allow_alternative_reaction_gaps : bool, False
        If a chain links two compounds by a reaction in the reaction network, and there are other
        "parallel" KGML reactions not in the reaction network that also link the compounds, then
        treat these parallel reactions as gaps with a value of True. With a value of False, ignore
        parallel reaction gaps.

    walker : KGMLNetworkWalker
        Walks chains of compounds linked by reactions in a KGML representation of a KEGG pathway.

    contigs_db : anvio.dbops.ContigsDatabase
        Contigs database containing reaction network.

    cog_function_source : str
        Most recent COG functional annotation source available in the contigs database. Gene hits to
        this source are considered.

    genes_in_contigs_df : pandas.core.frame.DataFrame
        Table of information on gene calls in contigs.

    gene_kos_df : pandas.core.frame.DataFrame
        Table of gene-KO hits, including all hits from the file at 'all_ko_hits_path', not just top
        hits recorded in the contigs database.

    top_gene_kos_df : pandas.core.frame.DataFrame
        Table of top-ranking gene-KO hits, those stored in the contigs database.

    gene_in_pathway : dict[str, bool]
        Records whether genes are linked via top-ranking KO hits to the pathway. Keys are gene
        caller IDs of genes with top-ranking KO hits recorded in the contigs database. Values are
        True if a KO is found in the pathway, and False if not.

    gene_cogs_df : pandas.core.frame.DataFrame
        Table of gene COG annotations, which is expected to only include the top hits.

    ko_cogs : dict[str, list[str]]
        Mapping of KO to COG IDs loaded from ko_cog_path.

    cog_kos : dict[str, list[str]]
        Mapping of COG to KO IDs loaded from ko_cog_path.

    ko_list_df : pandas.core.frame.DataFrame
        Table of all KOs and their names from the anvi'o KEGG installation.

    ungapped_chains : dict[str, list[Chain]]
        Chains from compounds in the KGML representation of the pathway, allowing zero genomic gaps
        in the chains.

    gapped_chains : dict[str, list[Chain]]
        Chains from compounds in the KGML representation of the pathway, allowing up to one genomic
        gap per chain.

    gap_analyzer : GapAnalyzer
        Compares the ungapped and gapped chains.

    ranked_gaps : list[str]
        KGML reaction ID of each gap in gapped_chains as ranked by gap_analyzer.rank_gaps. In
        essence, higher ranking (first) gaps occur toward the middle of longer chains and lower
        ranking (last) gaps occur toward the edges of shorter chains.

    kgml_reaction_id_ko_definitions : dict[str, list[tuple[str, str]]]
        KOs associated with KGML reactions in gapped chains. Keys are the KGML reaction IDs of gaps
        in the pathway. Value lists have an entry for each KO associated with the KGML reaction. The
        tuple consists of the KO ID and name.

    self.gapped_chains_ko_ids : dict[str, list[list[tuple[str]]]]
        KOs linking genes to KGML reactions in gapped chains. Keys are the KGML reaction IDs of gaps
        in the pathway. Value outer lists have an entry for each chain containing the gap. The inner
        list for each chain has an entry for each KGML reaction in the chain, in order. The tuple
        value for the reaction contains KO IDs in the contigs database reaction network that are
        top-ranking gene hits.

    pathway_ortholog_entries : list[Entry]
        Ortholog entries in the KGML KO type pathway.

    gap_syntenous_regions : dict[str, list[list[int]]]
        Syntenous regions around genes in gapped chains. Keys are the KGML reaction IDs of gaps in
        the pathway. Values are lists of unique syntenous regions around genes in chains containing
        the gap. Each inner list identifies a syntenous region as the row indices of full-length
        gene calls in 'genes_in_contigs_df'.

    gap_chain_gcids : dict[str, list[int]]
        Genes annotated by KOs associated with KGML reactions in gapped chains. Keys are the KGML
        reaction IDs of gaps in the pathway. Values are lists of gene caller IDs linked to reactions
        in the chains containing the gap.
    """
    def __init__(self, args: Namespace):
        """
        Parameters
        ==========
        args : argparse.Namespace
            Contains arguments, listed below. See the class docstring for more information on
            arguments set as attributes, including default values. The required arguments are
            kegg_pathway_number, contigs_db_path, all_ko_hits_path, and ko_cog_path.

        'kegg_pathway_number' : str

        'contigs_db_path' : str

        'all_ko_hits_path' : str

        'ko_cog_path' : str

        'compound_fate' : Literal['consume', 'produce']

        'max_reactions' : int

        'allow_alternative_reaction_gaps' : bool
        """
        A = lambda x, y: args.__dict__[x] if x in args.__dict__ else y

        self.kegg_pathway_number: str = A('kegg_pathway_number', None)
        self.contigs_db_path: str = A('contigs_db_path', None)
        if self.kegg_pathway_number is None or self.contigs_db_path is None:
            raise ConfigError(
                "A KEGG pathway number (args.kegg_pathway_number) and the path to a contigs "
                "database containing a reaction network (args.contigs_db_path) should be provided "
                "for initialization."
            )

        self.all_ko_hits_path: str = A('all_ko_hits_path', None)
        if self.all_ko_hits_path is None:
            raise ConfigError(
                "The path to the table of all KO hits to genes in the contigs "
                "(args.all_ko_hits_path) should be provided for initialization."
            )

        self.ko_cog_path: str = A('ko_cog_path', None)
        if self.ko_cog_path is None:
            raise ConfigError(
                "The path to the KEGG binary relations file mapping KO to COG IDs "
                "(args.ko_cog_path) should be provided for initialization."
            )

        self.compound_fate: str = A('compound_fate', 'consume')
        self.max_reactions: Union[int, None] = A('max_reactions', None)
        self.allow_alternative_reaction_gaps: bool = A('allow_alternative_reaction_gaps', False)
        walker_args = Namespace(
            kegg_pathway_number=self.kegg_pathway_number,
            contigs_db_path=self.contigs_db_path,
            compound_fate=self.compound_fate,
            max_reactions=self.max_reactions,
            allow_alternative_reaction_gaps=self.allow_alternative_reaction_gaps
        )
        self.walker = KGMLNetworkWalker(walker_args)

        self.contigs_db = ContigsDatabase(self.contigs_db_path)
        if 'KOfam' not in self.contigs_db.meta['gene_function_sources']:
            raise ConfigError(
                "The contigs database has no 'KOfam' functional annotation source. This is "
                "puzzling since a reaction network was loaded from the contigs database, and the "
                "network should contain gene KO hits that were initially recorded in the database."
            )
        for version in ['24', '20', '14']:
            # Take the most recent version available starting with COG 2014.
            if f'COG{version}_FUNCTION' in self.contigs_db.meta['gene_function_sources']:
                break
        else:
            raise ConfigError(
                "The contigs database has no COG functional annotation source ('COG24_FUNCTION', "
                "'COG20_FUNCTION', or 'COG14_FUNCTION'). The gap filler requires COG annotation."
            )
        self.cog_function_source = f'COG{version}_FUNCTION'

        # Load the table containing information on gene calls and where they are located in contigs.
        # Sort gene rows by contig start position.
        self.genes_in_contigs_df = self.contigs_db.db.get_table_as_dataframe('genes_in_contigs')
        self.genes_in_contigs_df = self.genes_in_contigs_df.sort_values(
            ['contig', 'start']
        ).reset_index()

        # Load a table containing information on all KO hits to genes in the contigs database,
        # including lower-ranking hits.
        gene_kos_df = pd.read_csv(
            self.all_ko_hits_path, sep='\t',
            usecols=['gene_callers_id', 'gene_hmm_id', 'gene_name', 'e_value']
        )[['gene_callers_id', 'gene_hmm_id', 'gene_name', 'e_value']]
        gene_kos_df = gene_kos_df.rename(
            {'gene_hmm_id': 'accession', 'gene_name': 'function'}, axis=1
        )
        self.gene_kos_df = gene_kos_df.sort_values(['gene_callers_id', 'e_value', 'accession'])

        # Load a table of the top gene-KO hits stored in the contigs database. Top hits to genes
        # from the table of all hits can be excluded from the table in the contigs database due to
        # not meeting a score threshold.
        gene_functions_df = self.contigs_db.db.get_table_as_dataframe('gene_functions')
        self.top_gene_kos_df = gene_functions_df[gene_functions_df['source'] == 'KOfam']

        # Record whether genes are linked to the KEGG pathway via top KO hits.
        self.gene_in_pathway: dict[str, bool] = {}
        ko_id_pathway_ids = lambda ko_id: self.walker.kegg_data.ko_data[ko_id]['PTH']
        for gcid, gene_df in self.top_gene_kos_df.groupby('gene_callers_id'):
            for ko_id in gene_df['accession']:
                try:
                    pathway_ids = ko_id_pathway_ids(ko_id)
                except KeyError:
                    # The KO is not in a pathway, at least as recorded in the anvi'o KEGG
                    # installation (modules database).
                    continue
                if self.kegg_pathway_number in pathway_ids:
                    self.gene_in_pathway[gcid] = True
                    break
            else:
                self.gene_in_pathway[gcid] = False

        # Load a table containing information on COG hits, splitting the table into one line per
        # hit, as in the gene KO hit table.
        expanded_gene_hits = []
        for row in gene_functions_df[
            gene_functions_df['source'] == self.cog_function_source
        ].itertuples(index=False):
            if '!!!' not in row.accession:
                expanded_gene_hits.append(row)
                continue
            for cog_id, cog_name in zip(row.accession.split('!!!'), row.function.split('!!!')):
                new_row = row._replace(accession=cog_id, function=cog_name)
                expanded_gene_hits.append(new_row)
        self.gene_cogs_df = pd.DataFrame(expanded_gene_hits, columns=gene_functions_df.columns)

        # Load the KEGG binary relations file mapping KO to COG IDs, and make dictionaries mapping
        # KOs and COGs to each other.
        self.ko_cogs: dict[str, list[str]] = {}
        self.cog_kos: dict[str, list[str]] = {}
        ko_cog_df = pd.read_csv(self.ko_cog_path, sep='\t')
        ko_cog_df.columns = ['ko', 'cog']
        for row in ko_cog_df.itertuples():
            cog_ids = row.cog[5: -1].split()
            self.ko_cogs[row.ko] = cog_ids
            for cog_id in cog_ids:
                try:
                    self.cog_kos[cog_id].append(row.ko)
                except KeyError:
                    self.cog_kos[cog_id] = [row.ko]

        # Load the table of all KO definitions in the anvi'o KEGG installation, needed to retrieve
        # the names of KOs associated with KGML reactions and gene COG hits.
        self.ko_list_df = pd.read_csv(
            self.walker.kegg_data.kegg_context.ko_list_file_path, sep='\t',
            usecols=['knum', 'definition']
        )
        self.ko_list_df = self.ko_list_df.rename(
            {'knum': 'accession', 'definition': 'function'}, axis=1
        )

        # Find chains in the pathway with zero gaps and up to one gap.
        self.ungapped_chains: list[Chain] = []
        for chains in self.walker.get_chains().values():
            self.ungapped_chains += chains
        self.walker.max_gaps = 1
        self.gapped_chains: list[Chain] = []
        for chains in self.walker.get_chains().values():
            self.gapped_chains += chains

        # Compare gapped to ungapped chains. Rank gaps to prioritize for gap filling.
        self.gap_analyzer = GapAnalyzer(self.gapped_chains, self.ungapped_chains)
        self.ranked_gaps = [
            kgml_reaction_ids[0] for kgml_reaction_ids in self.gap_analyzer.rank_gaps()
        ]

        # For each KGML reaction in a gapped chain, record the associated KOs (ID and name).
        self.kgml_reaction_id_ko_definitions: dict[str, list[tuple[str, str]]] = {}
        ko_list_df = self.ko_list_df.set_index('accession')
        for shared_gaps in self.gap_analyzer.gap_relations.values():
            for relations in shared_gaps.gap_chain_relations:
                gapped_chain = relations.gappy_chain
                for kgml_reaction in gapped_chain.kgml_reactions:
                    kgml_reaction_id = kgml_reaction.id
                    if kgml_reaction_id in self.kgml_reaction_id_ko_definitions:
                        continue
                    ko_ids = self.walker.rn_pathway_kgml_reaction_id_to_ko_ids[kgml_reaction_id]
                    ko_definitions = []
                    for ko_id in ko_ids:
                        ko_name = ko_list_df.loc[ko_id]['function']
                        ko_definitions.append((ko_id, ko_name))
                    self.kgml_reaction_id_ko_definitions[kgml_reaction_id] = ko_definitions

        # Record which KOs link genes to reactions in gapped chains.
        self.gapped_chains_ko_ids: dict[str, list[list[tuple[str]]]] = {}
        for key, shared_gaps in self.gap_analyzer.gap_relations.items():
            kgml_reaction_id = key[0]
            gapped_chains = [relations.gappy_chain for relations in shared_gaps.gap_chain_relations]
            self.gapped_chains_ko_ids[kgml_reaction_id] = [
                [tuple([ko.id for ko in kos]) for kos in gapped_chain.network_kos]
                for gapped_chain in gapped_chains
            ]

        # Find the syntenous regions around genes in gapped chains. Map gap reaction IDs to
        # syntenous regions.
        self.pathway_ortholog_entries = self.walker.kgml_ko_pathway.get_entries(
            entry_type='ortholog'
        )
        self.gap_syntenous_regions: dict[str, list[list[int]]] = {}
        self.gap_chain_gcids: dict[str, list[int]] = {}
        for kgml_reaction_id in self.ranked_gaps:
            syntenous_regions, gcids = self.get_gap_kgml_reaction_syntenous_regions(
                kgml_reaction_id
            )
            self.gap_syntenous_regions[kgml_reaction_id] = syntenous_regions
            self.gap_chain_gcids[kgml_reaction_id] = gcids

    def get_gap_kgml_reaction_syntenous_regions(
        self,
        kgml_reaction_id: str
    ) -> tuple[list[list[int]], list[int]]:
        """
        Find "syntenous regions" of genes in the same orientation around genes in chains with the
        reaction gap.

        Parameters
        ==========
        kgml_reaction_id : str
            Gap reaction ID that occurs in one or more chains.

        Returns
        =======
        list[list[int]]
            List of unique syntenous regions. An empty list is returned if the reaction ID does not
            represent a gap.

        list[int]
            Gene caller IDs of genes annotated by KOs associated with reactions in chains containing
            the gap.
        """
        key = (kgml_reaction_id, )
        try:
            shared_gaps = self.gap_analyzer.gap_relations[key]
        except KeyError:
            return []

        gcids: list[int] = []
        for gap_chain_relations in shared_gaps.gap_chain_relations:
            gapped_chain = gap_chain_relations.gappy_chain
            chain_gcids = self.get_chain_gcids(gapped_chain)
            gcids += chain_gcids
        gcids = sorted(set(gcids))
        syntenous_regions = self.get_syntenous_regions(gcids)

        return syntenous_regions, gcids

    def get_chain_gcids(self, chain: Chain) -> list[int]:
        """
        Get the gene caller IDs of genes annotated by KOs associated with reactions in the chain.

        Parameters
        ==========
        chain : Chain
            Chain in the pathway map.

        Returns
        =======
        list[int]
            Gene caller IDs.
        """
        chain_gcids: list[int] = []
        for is_gap, kgml_reaction in chain.kgml_reactions:
            if is_gap:
                continue

            # Find the KGML ortholog entry for the reaction.
            for entry in self.pathway_ortholog_entries:
                if entry.id == kgml_reaction.id:
                    kgml_entry = entry
                    break
            else:
                raise AssertionError

            # The entry corresponds to one or more KOs. Record IDs of genes annotated by the KOs.
            ko_ids: list[str] = []
            for name_substring in kgml_entry.name.split():
                if name_substring[:3] == 'ko:':
                    ko_ids.append(name_substring[3:])
            chain_gcids += self.top_gene_kos_df[
                self.top_gene_kos_df['accession'].isin(ko_ids)
            ]['gene_callers_id'].unique().tolist()
        chain_gcids = sorted(set(chain_gcids))

        return chain_gcids

    def get_syntenous_regions(self, gcids: list[int]) -> list[list[int]]:
        """
        Find "syntenous regions" of genes in the same orientation around any of the given genes.

        Parameters
        ==========
        gcids : list[int]
            Gene caller IDs of genes to search around.

        Returns
        =======
        list[list[int]]
            Unique list of syntenous regions, avoiding duplicate regions around multiple input
            genes.
        """
        syntenous_regions: list[list[int]] = []
        for gcid in gcids:
            syntenous_region = self.get_syntenous_region(gcid)
            if syntenous_region in syntenous_regions:
                # The syntenous region has already been encountered in considering another gene.
                continue
            syntenous_regions.append(syntenous_region)

        return syntenous_regions

    def get_syntenous_region(self, gcid: int) -> list[int]:
        """
        Find the "syntenous region" of genes in the same orientation around the given gene.

        Parameters
        ==========
        gcid : int
            Gene caller ID of gene to search around.

        Returns
        =======
        list[int]
            The syntenous region is recorded as a list of row indices of full-length gene calls in
            the 'genes_in_contigs_df' attribute, which is the 'genes_in_contigs' table of the
            contigs database sorted by contig ID and, within each contig, start position of the gene
            call. The input gene is included in the syntenous region.
        """
        # The search around each gene stops in either direction when a gene in the opposite
        # orientation is found or the first or last gene in the contig is reached.
        row = self.genes_in_contigs_df[
            self.genes_in_contigs_df['gene_callers_id'] == gcid
        ].squeeze()
        row_index = row.name
        direction = row['direction']

        syntenous_region: list[int] = []
        # Search preceding genes.
        for row_index in range(row_index, -1, -1):
            lower_row = self.genes_in_contigs_df.loc[row_index]
            if lower_row['direction'] != direction:
                break
            if lower_row['partial'] == 1:
                # Ignore partial gene calls.
                continue
            syntenous_region.append(row_index)
        # Search succeeding genes.
        for row_index in range(row_index + 1, len(self.genes_in_contigs_df)):
            higher_row = self.genes_in_contigs_df.loc[row_index]
            if higher_row['direction'] != direction:
                break
            if higher_row['partial'] == 1:
                continue
            syntenous_region.append(row_index)
        syntenous_region.sort()

        return syntenous_region

    def eval_pathway(self) -> dict:
        """
        Evaluate genomic evidence to fill reaction gaps in the pathway map.

        The returned dictionary contains a dictionary for each gap. These are ordered, roughly
        speaking, with those closer to the center of longer chains coming before those closer to the
        edge of shorter chains. Each dictionary contains a list of candidate genes to fill the gap,
        which is empty if no candidates were found.

        Returns
        =======
        dict
            JSON-formatted dictionary of evidence for how gaps can be filled.
        """
        json_pathway = {}
        json_pathway['pathway_number'] = self.kegg_pathway_number
        json_pathway['pathway_name'] = self.walker.kgml_rn_pathway.title
        json_gaps: list[dict] = []
        json_pathway['gaps'] = json_gaps
        for kgml_reaction_id in self.ranked_gaps:
            json_gap = self.eval_gap_kgml_reaction(kgml_reaction_id)
            json_gaps.append(json_gap)

        return json_pathway

    def eval_gap_kgml_reaction(self, kgml_reaction_id : str) -> Union[dict, None]:
        """
        Evaluate genomic evidence to fill a particular reaction gap in the pathway map.

        Evidence supports hypotheses for how particular genes can fill the gap. Evidence is based on
        alternative KO annotations of genes and synteny with other genes in the pathway. Here are
        lines of evidence to fill a gap.

        1. A gene's top-scoring KO annotation can be replaced by a lower-scoring KO annotation
        that corresponds to the gap reaction.

        2. A gene's top-scoring COG annotation maps to a KO that corresponds to the gap reaction.
        If the gene does not have a KO annotation, or has a different KO annotation, it can be
        replaced by the KO suggested by the COG.

        3. A gene is adjacent to other syntenous genes (in the same orientation) in the pathway,
        or stronger yet, in the chain with the gap. This gene may also have an alternative KO
        annotation to fill the gap (see 1 and 2), or may be unannotated.

        Here are practical examples of how reported evidence can be used in gap-filling.

        1. A gap can be filled by reannotation of a gene with a lower-ranking KO for a bifunctional
        protein. For instance, there is a gap for the ArgA reaction in arginine biosynthesis, but a
        gene is annotated with the ArgB KO for the subsequent reaction. The gene has a lower-scoring
        ArgAB bifunctional KO which fills the gap. The gene is adjacent to syntenous genes in the
        gapped biosynthetic chain, and there are no adjacent unannotated genes that may represent
        ArgA. ArgAB would represent a gene fusion, which can be confirmed by separate protein domain
        hits to ArgA and ArgB in the amino acid sequence.

        2. A gap can represent an alternative reaction that is missing from the map, the
        identification of which allows the chain to be completed. For instance, some bacteria such
        as Bacteroides fragilis have an arginine biosynthetic pathway variant that produces arginine
        through succinylated derivatives of glutamate rather than typical acetylated derivatives.
        There is a KO for N-succinyl-L-ornithine carbamoyltransferase which is similar to the KO for
        N-acetyl-L-ornithine carbamoyltransferase. However, the KEGG map only contains reactions
        involving acetylated intermediates. Therefore, genomes with the succinyl version of the
        pathway have a map that looks complete except for the carbamoyltransferase gap. However, the
        gene in question can have a lower-ranking KO for the N-acetyl variant, which fills the gap,
        and this is especially obvious if the gene is adjacent to syntenous arginine biosynthetic
        genes in the gapped chain. In the context of metabolic modeling, identification of this
        variant indicates that an automatically annotated acetyl version of the pathway might need
        to be replaced by a succinyl version of the pathway, with potentially important implications
        for the metabolic network given the alternative origins of succinyl-CoA and acetyl-CoA.

        Parameters
        ==========
        kgml_reaction_id : str
            ID of KGML reaction gap.

        Returns
        =======
        Union[dict, None]
            JSON-formatted dictionary of gap-filling evidence. The dictionary is empty if no KOs
            associated with the KGML reaction hit genes in the contigs database. None is returned if
            the requested KGML reaction ID is not found in the pathway or is not a gap.

            The following shows the template of the returned dictionary.
            {
                "pathway_number": <pathway number>,
                "pathway_name": <pathway name>,
                "gaps": [
                    {
                        "kgml_reaction_id": <KGML reaction ID>,
                        "candidate_genes": [
                            {
                                "gene_callers_id": <GCID>,
                                "gap_ko_hits": [
                                    {
                                        "ko_id": <KO ID>,
                                        "ko_name": <KO name>,
                                        "e_value": <E value>
                                    },
                                    ...
                                ],
                                "other_ko_hits": [
                                    {
                                        "ko_id": <KO ID>,
                                        "ko_name": <KO name>,
                                        "e_value": <E value>
                                    },
                                    ...
                                ],
                                "gap_kos_via_cog_top_hits": [
                                    {
                                        "ko_id": <KO ID>,
                                        "ko_name": <KO name>,
                                        "cog_ids": [<COG ID>, ...]
                                    },
                                    ...
                                ],
                                "cog_top_hits": [
                                    {
                                        "cog_id": <COG ID>,
                                        "cog_name": <COG name>,
                                        "e_value": <E value>
                                    },
                                    ...
                                ],
                                "syntenous_region": { or null
                                    "id": <ID>,
                                    "gene_index": <gene index>,
                                    "is_gene_elsewhere_in_gapped_chain": <true or false>,
                                    "is_gene_elsewhere_in_pathway": <true or false>,
                                    "region_gene_count": <count of full genes in region>,
                                    "region_pathway_gene_count":
                                        <count of genes in region and pathway>,
                                    "region_gapped_chain_gene_count":
                                        <count of genes in region and gapped chains>,
                                    "adjacent_pathway_genes": [
                                        {
                                            "gene_index": <gene index>,
                                            "is_gene_in_gapped_chain": <true or false>
                                        },
                                        ...
                                    ]
                                }
                            },
                            ...
                        ],
                        "syntenous_regions": [
                            {
                                "id": <id>,
                                "contig_id": <contig ID>,
                                "gene_orientation": <"f" or "r">,
                                "full_genes": [
                                    {
                                        "gene_index": <gene index>,
                                        "is_gap_candidate": <true or false>,
                                        "start_position": <start position>,
                                        "stop_position": <stop position>,
                                        "gene_callers_id": <GCID>,
                                        "gene_call_source": <gene call source>,
                                        "ko_top_hits": [
                                            {
                                                "ko_id": <KO ID>,
                                                "ko_name": <KO name>,
                                                "e_value": <E value>
                                            },
                                            ...
                                        ],
                                        "cog_top_hits": [
                                            {
                                                "cog_id": <COG ID>,
                                                "cog_name": <COG name>,
                                                "e_value": <E value>
                                            },
                                            ...
                                        ],
                                        "is_in_pathway": <true or false>,
                                        "gapped_chains": [
                                            {
                                                "gapped_chain_index": <gapped chain index>,
                                                "reaction_count": <reaction count>,
                                                "reaction_indices_in_chain":
                                                    [<reaction_index_in_chain>, ...]
                                            },
                                            ...
                                        ]
                                    }
                                ]
                            }
                        ],
                        "gapped_chains": [
                            {
                                "gapped_chain_index": <gapped chain index>,
                                "reactions": [
                                    {
                                        "reaction_index_in_chain": <reaction index in chain>,
                                        "is_gap": <true or false>,
                                        "kgml_reaction_id": <KGML reaction ID>,
                                        "kos": [
                                            {
                                                "ko_id": <KO ID>,
                                                "ko_name": <KO name>
                                            },
                                            ...
                                        ]
                                    },
                                    ...
                                ]
                            },
                            ...
                        ]
                    },
                    ...
                ]
            }
        """
        def get_json_candidate_gene_sort_key(
            json_candidate_gene: dict
        ) -> tuple[int, float, float, int, int, int]:
            """
            Get a sort key for the JSON-formatted dictionary of information on a gap-filling
            candidate gene.

            The sort first depends on the presence or absence of a particular type of gap-filling
            evidence and then, to break ties within presence/absence tiers, numerical values from
            each tier.

            The following are the types of presence/absence evidence for a candidate gene.
            A. The gene has a top-scoring COG hit that is associated with a KO linked to the gap
            KGML reaction.
            B. The gene has lower-ranking KO hit linked to the gap KGML reaction.
            C. The gene is syntenous with other genes linked to chains containing the gap.

            The following are the presence/absence tiers indicated by their identifying values (0-5)
            for the first level of the sort.
            0: The gene has evidence types A + B + C.
            1: A + C
            2: B + C
            3: A
            4: B
            5: C

            The following are the numerical values for the second through sixth levels of the sort.
            The second and third levels are sorted in ascending order, and the fourth, fifth, and
            sixth levels are descending.
            1. E value of the best-scoring, lower-ranking KO hit linked to the gap reaction. In the
            absence of a lower-ranking KO hit, infinity is used as a placeholder value.
            2. E value of the best-scoring COG associated with a KO linked to the reaction. In the
            absence of such a COG, infinity is used as a placeholder value.
            3. Count of adjacent pathway genes in the syntenous region that contains genes linked to
            gapped chains. If the candidate gene is not in a syntenous region flanked by genes in
            gapped chains, then 0 is used as a placeholder value.
            4. Count of genes (not necessarily adjacent to the candidate gene) linked to gapped
            chains in the syntenous region. If the candidate gene is not in a syntenous region
            containing genes in gapped chains, then 0 is used as a placeholder value.
            5. Count of pathway genes (not necessarily adjacent to the candidate gene) in the
            syntenous region. If the candidate gene is not in a syntenous region containing genes in
            the pathway, then 0 is used as a placeholder value.

            Further ties are unlikely and not systematically resolved.

            Parameters
            ==========
            json_candidate_gene : list[dict]
                JSON-formatted dictionary of information on a gap-filling candidate gene.

            Returns
            =======
            tuple[int, float, float, int, int, int]
                Sort key, with the first level being presence/absence tier of evidence (0-5), and
                the second, third, and fourth levels being a numerical value for each type of
                evidence.
            """
            evidence = []
            has_cog_ko_evidence = len(json_candidate_gene['gap_kos_via_cog_top_hits']) > 0
            evidence.append(has_cog_ko_evidence)
            has_lower_ranking_ko_evidence = len(json_candidate_gene['gap_ko_hits']) > 0
            evidence.append(has_lower_ranking_ko_evidence)
            has_synteny_evidence = json_candidate_gene['syntenous_region'] is not None
            evidence.append(has_synteny_evidence)
            evidence = tuple(evidence)
            if evidence == (True, True, True):
                evidence_tier = 0
            elif evidence == (True, False, True):
                evidence_tier = 1
            elif evidence == (False, True, True):
                evidence_tier = 2
            elif evidence == (True, False, False):
                evidence_tier = 3
            elif evidence == (False, True, False):
                evidence_tier = 4
            elif evidence == (False, False, False):
                evidence_tier = 5
            sort_key = [evidence_tier]

            ko_e_values = [
                gap_ko_hit['e_value'] for gap_ko_hit in json_candidate_gene['gap_ko_hits']
            ]
            min_ko_e_value = min(
                ko_e_values
            ) if len(json_candidate_gene['gap_ko_hits']) > 0 else float('inf')
            sort_key.append(min_ko_e_value)

            cog_ids: list[str] = []
            for gap_ko_via_cog_top_hits in json_candidate_gene['gap_kos_via_cog_top_hits']:
                for cog_id in gap_ko_via_cog_top_hits['cog_ids']:
                    cog_ids.append(cog_id)
            cog_e_values = [
                cog_top_hit['e_value'] for cog_top_hit in json_candidate_gene['cog_top_hits']
                if cog_top_hit['cog_id'] in set(cog_ids)
            ]
            min_cog_e_value = min(cog_e_values) if cog_e_values else float('inf')
            sort_key.append(min_cog_e_value)

            json_syntenous_region = json_candidate_gene['syntenous_region']
            if json_syntenous_region is None:
                sort_key.extend([0, 0, 0])
            else:
                sort_key.append(len(json_syntenous_region['adjacent_pathway_genes']))
                sort_key.append(json_syntenous_region['region_gapped_chain_gene_count'])
                sort_key.append(json_syntenous_region['region_pathway_gene_count'])

            return sort_key

        try:
            # Get KO IDs associated with the KGML reaction.
            ko_ids = self.walker.rn_pathway_kgml_reaction_id_to_ko_ids[kgml_reaction_id]
        except KeyError:
            # The KGML reaction ID is not found in the pathway.
            return None

        # Get chains containing the reaction gap. Make template JSON dictionaries used in every JSON
        # candidate gene entry: these contain information on the position of the gap and the number
        # of reactions in the chains. Make JSON dictionaries, added later to the gap JSON dictionary
        # after other attributes, that record information on the gapped chains as a whole.
        try:
            shared_gaps = self.gap_analyzer.gap_relations[(kgml_reaction_id, )]
        except KeyError:
            # The KGML reaction ID is not a gap in a chain.
            return None
        gapped_chains = [relations.gappy_chain for relations in shared_gaps.gap_chain_relations]
        json_full_gene_gapped_chains_template: list[dict] = []
        json_gapped_chains: list[dict] = []
        for gapped_chain_index, gapped_chain in enumerate(gapped_chains):
            json_full_gene_gapped_chain_template = {}
            json_full_gene_gapped_chains_template.append(json_full_gene_gapped_chain_template)
            json_full_gene_gapped_chain_template['gapped_chain_index'] = gapped_chain_index
            json_full_gene_gapped_chain_template['reaction_count'] = len(
                gapped_chain.kgml_reactions
            )

            json_gapped_chain = {}
            json_gapped_chains.append(json_gapped_chain)
            json_gapped_chain['gapped_chain_index'] = gapped_chain_index
            json_gapped_chain_reactions: list[dict] = []
            json_gapped_chain['reactions'] = json_gapped_chain_reactions

            for kgml_reaction_index, kgml_reaction in enumerate(gapped_chain.kgml_reactions):
                json_gapped_chain_reaction = {}
                json_gapped_chain_reactions.append(json_gapped_chain_reaction)

                json_gapped_chain_reaction['reaction_index_in_chain'] = kgml_reaction_index

                if kgml_reaction.id == kgml_reaction_id:
                    json_full_gene_gapped_chain_template[
                        'reaction_indices_in_chain'
                    ] = kgml_reaction_index
                    json_gapped_chain_reaction['is_gap'] = True
                else:
                    json_gapped_chain_reaction['is_gap'] = False

                json_gapped_chain_reaction['kgml_reaction_id'] = kgml_reaction.id

                json_kos: list[dict] = []
                json_gapped_chain_reaction['kos'] = json_kos
                for ko_definitions in self.kgml_reaction_id_ko_definitions[kgml_reaction.id]:
                    json_ko = {}
                    json_kos.append(json_ko)
                    json_ko['ko_id'], json_ko['ko_name'] = ko_definitions

        # Get the GCIDs linked to the gapped chains.
        gap_chain_gcids = self.gap_chain_gcids[kgml_reaction_id]

        # Get information on the syntenous regions around genes in the chains containing the gap.
        gap_syntenous_regions = self.gap_syntenous_regions[kgml_reaction_id]
        json_syntenous_regions = self.get_gap_synteny_info(kgml_reaction_id)

        # Get information regarding gap KOs associated with genes, recording at the gene level.
        gcid_records: dict[int, list[dict]] = {}
        for ko_id in ko_ids:
            ko_info = self.get_ko_info(ko_id)
            if not ko_info:
                continue
            for json_candidate_gene in ko_info:
                gcid = json_candidate_gene['gene_callers_id']
                try:
                    gcid_records[gcid].append(json_candidate_gene)
                except KeyError:
                    gcid_records[gcid] = [json_candidate_gene]

        # Create a JSON dictionary for the gap. Copy data into new gene records, adding KO-related
        # records for each gene. Then add information on the syntenous region of each gene.
        json_gap = {}
        json_gap['kgml_reaction_id'] = kgml_reaction_id
        json_candidate_genes: list[dict] = []
        json_gap['candidate_genes'] = json_candidate_genes
        for gcid, records in gcid_records.items():
            json_candidate_gene = {}
            json_candidate_genes.append(json_candidate_gene)
            json_candidate_gene['gene_callers_id'] = gcid

            # Compile gap KO hits to the gene.
            json_gap_ko_hits: list[dict] = []
            json_candidate_gene['gap_ko_hits'] = json_gap_ko_hits
            gap_ko_ids: list[str] = []
            for record in records:
                assert len(record['gap_ko_hits']) == 1
                json_gap_ko_hit = record['gap_ko_hits'][0]
                json_gap_ko_hits.append(json_gap_ko_hit)
                gap_ko_ids.append(json_gap_ko_hit['ko_id'])

            # Get other KO hits to the gene.
            json_other_ko_hits: list[dict] = []
            json_candidate_gene['other_ko_hits'] = json_other_ko_hits
            record = records[0]
            assert len(record['other_ko_hits']) == 1
            for json_other_ko_hit in record['other_ko_hits']:
                if json_other_ko_hit['ko_id'] not in gap_ko_ids:
                    json_other_ko_hits.append(json_other_ko_hit)

            # Compile gap KOs associated with COG hits to the gene.
            json_gap_kos_via_cog_top_hits: list[dict] = []
            json_candidate_gene['gap_kos_via_cog_top_hits'] = json_gap_kos_via_cog_top_hits
            for record in records:
                assert len(record['gap_kos_via_cog_top_hits']) == 1
                json_gap_ko_via_cog_top_hits = record['gap_kos_via_cog_top_hits'][0]
                json_gap_kos_via_cog_top_hits.append(json_gap_ko_via_cog_top_hits)

            # Get COG top hits to the gene.
            json_candidate_gene['cog_top_hits'] = records[0]['cog_top_hits']

            # Get information on the syntenous region around the gene candidate suggested by KO.
            syntenous_region = self.get_syntenous_region(gcid)
            for gap_syntenous_region_index, gap_syntenous_region in enumerate(gap_syntenous_regions):
                if syntenous_region == gap_syntenous_region:
                    # The region is the same as one of the regions found around genes in the chains
                    # containing the gap.
                    json_syntenous_region = json_syntenous_regions[gap_syntenous_region_index]
                    break
            else:
                # The gene candidate to fill the gap does not occur in any of the syntenous regions
                # around genes in the chains containing the gap. Record this region.
                json_syntenous_region = self.make_json_syntenous_region(syntenous_region)
                json_syntenous_regions.append(json_syntenous_region)
                json_syntenous_region['id'] = gap_syntenous_region_index + 1
            # Fill in information on the chains containing the gap, including the number of KGML
            # reactions in the chain and the index of the gap reaction in the chain.
            for gap_gene_index, json_full_gene in enumerate(json_syntenous_region['full_genes']):
                if json_full_gene['gene_callers_id'] == gcid:
                    break
            else:
                raise AssertionError(
                    "The gene candidate to fill the gap should have been recorded in the syntenous "
                    "region."
                )
            json_full_gene['is_gap_candidate'] = True
            json_full_gene['gapped_chains'] = deepcopy(json_full_gene_gapped_chains_template)

            json_candidate_gene_syntenous_region = {}
            json_candidate_gene_syntenous_region['id'] = json_syntenous_region['id']
            json_candidate_gene_syntenous_region['gene_index'] = gap_gene_index
            # The gap gene candidate can be a gene already in a chain containing the gap, e.g., the
            # gene is reannotated as a bifunctional gene, such as ArgB being reannotated as ArgAB to
            # fill the gap for ArgA. Record whether the gene is elsewhere in a gapped chain.
            is_gene_elsewhere_in_gapped_chain = gcid in gap_chain_gcids
            json_candidate_gene_syntenous_region[
                'is_gene_elsewhere_in_gapped_chain'
            ] = is_gene_elsewhere_in_gapped_chain
            # Likewise the gap gene candidate may be found elsewhere in the pathway map.
            json_candidate_gene_syntenous_region[
                'is_gene_elsewhere_in_pathway'
            ] = True if is_gene_elsewhere_in_gapped_chain else self.gene_in_pathway[gcid]
            json_candidate_gene_syntenous_region['region_gene_count'] = len(
                json_syntenous_region['full_genes']
            )
            json_candidate_gene_syntenous_region['region_pathway_gene_count'] = sum(
                [jg['is_in_pathway'] for jg in json_syntenous_region['full_genes']]
            )
            json_candidate_gene_syntenous_region['region_gapped_chain_gene_count'] = sum(
                [1 if jg['gapped_chains'] else 0 for jg in json_syntenous_region['full_genes']]
            )
            # Record genes that directly surround the gap candidate gene in the syntenous region and
            # that are also in the pathway.
            json_adjacent_pathway_genes = self.make_json_adjacent_pathway_genes(
                json_syntenous_region, gap_gene_index, gap_chain_gcids
            )
            json_candidate_gene_syntenous_region[
                'adjacent_pathway_genes'
            ] = json_adjacent_pathway_genes

        # All syntenous regions related to the gap have been found at this point.
        json_gap['syntenous_regions'] = json_syntenous_regions
        # Record genes in syntenous regions that are not gap-filling candidates.
        for json_syntenous_region in json_syntenous_regions:
            for json_full_gene in json_syntenous_region['full_genes']:
                if json_full_gene['is_gap_candidate'] is None:
                    json_full_gene['is_gap_candidate'] = False

        # Search syntenous regions for KO-unannotated genes that are adjacent to genes in a gapped
        # chain. Record adjacent unannotated genes as gap-filling candidates.
        for json_syntenous_region in json_syntenous_regions[: len(gap_syntenous_regions)]:
            for gene_index, json_full_gene in enumerate(json_syntenous_region['full_genes']):
                if json_full_gene['ko_top_hits'] or json_full_gene['is_gap_candidate']:
                    # The gene is annotated with a KO or is a candidate to fill the gap.
                    continue

                # Check on both sides for a gene in a gapped chain.
                is_adjacent = False
                if gene_index - 1 >= 0:
                    other_json_full_gene = json_syntenous_region['full_genes'][gene_index - 1]
                    if other_json_full_gene['gapped_chains']:
                        is_adjacent = True
                if not is_adjacent and gene_index + 1 < len(json_syntenous_region['full_genes']):
                    other_json_full_gene = json_syntenous_region['full_genes'][gene_index + 1]
                    if other_json_full_gene['gapped_chains']:
                        is_adjacent = True
                if not is_adjacent:
                    # There is no flanking gene in a gapped chain.
                    continue

                # Create a record for the gap-filling candidate gene.
                json_candidate_gene = {}
                json_candidate_genes.append(json_candidate_gene)
                json_candidate_gene['gene_callers_id'] = json_full_gene['gene_callers_id']
                # The gene doesn't have any KO hits.
                json_candidate_gene['gap_ko_hits'] = []
                json_candidate_gene['other_ko_hits'] = []
                # The gene doesn't have COG hits aliasing KOs that can fill the gap.
                json_candidate_gene['gap_kos_via_cog_top_hits'] = []
                # The gene may have COG hits.
                json_candidate_gene['cog_top_hits'] = deepcopy(json_full_gene['cog_top_hits'])
                json_candidate_gene_syntenous_region = {}
                json_candidate_gene['syntenous_region'] = json_candidate_gene_syntenous_region
                json_candidate_gene_syntenous_region['id'] = json_syntenous_region['id']
                json_candidate_gene_syntenous_region['gene_index'] = gene_index
                # Without KO hits, the gene can't be linked to reactions in the gapped chain, nor
                # can it be in the pathway.
                json_candidate_gene_syntenous_region['is_gene_elsewhere_in_gapped_chain'] = False
                json_candidate_gene_syntenous_region['is_gene_elsewhere_in_pathway'] = False
                json_candidate_gene_syntenous_region['region_gene_count'] = len(
                    json_syntenous_region['full_genes']
                )
                json_candidate_gene_syntenous_region['region_pathway_gene_count'] = sum(
                    [jg['is_in_pathway'] for jg in json_syntenous_region['full_genes']]
                )
                json_candidate_gene_syntenous_region['region_gapped_chain_gene_count'] = sum(
                    [1 if jg['gapped_chains'] else 0 for jg in json_syntenous_region['full_genes']]
                )
                # Record genes that directly surround the gap candidate gene in the syntenous region
                # and that are also in the pathway.
                json_adjacent_pathway_genes = self.make_json_adjacent_pathway_genes(
                    json_syntenous_region, gene_index, gap_chain_gcids
                )
                json_candidate_gene_syntenous_region[
                    'adjacent_pathway_genes'
                ] = json_adjacent_pathway_genes

        # Add information on the gapped chains to the gap record.
        json_gap['gapped_chains'] = json_gapped_chains

        # Sort gap-filling candidate genes.
        json_gap['candidate_genes'] = sorted(
            json_gap['candidate_genes'], key=get_json_candidate_gene_sort_key
        )

        return json_gap

    def get_gap_synteny_info(self, kgml_reaction_id: str) -> Union[list[dict], None]:
        """
        Get information on all syntenous regions around genes linked to chains containing the gap.

        Parameters
        ==========
        kgml_reaction_id : str
            ID of KGML reaction gap.

        Returns
        =======
        list[dict]
            JSON-formatted dictionary for each syntenous region. None is returned if the requested
            KGML reaction ID is not a gap in a chain.
        """
        # Find which KOs link genes to reactions in each chain.
        try:
            gapped_chains_ko_ids = self.gapped_chains_ko_ids[kgml_reaction_id]
        except KeyError:
            # The KGML reaction ID does not correspond to a gap in a chain.
            return None

        # Get syntenous regions around the genes in the chains containing the reaction gap. Create a
        # JSON-formatted dictionary for each syntenous region.
        json_syntenous_regions: list[dict] = []
        for syntenous_region_index, syntenous_region in enumerate(
            self.gap_syntenous_regions[kgml_reaction_id]
        ):
            json_syntenous_region = self.make_json_syntenous_region(syntenous_region)
            json_syntenous_regions.append(json_syntenous_region)
            json_syntenous_region['id'] = syntenous_region_index

            # Record the gapped chains that contain the gene, including the number of KGML reactions
            # in the chain, and the indices of KGML reactions linked to the gene via top KO hits.
            for json_full_genes in json_syntenous_region['full_genes']:
                json_ko_top_hits: list[dict] = json_full_genes['ko_top_hits']
                json_gapped_chains: list[dict] = json_full_genes['gapped_chains']
                for gapped_chain_index, gapped_chain_ko_ids in enumerate(gapped_chains_ko_ids):
                    reaction_indices: list[int] = []
                    for json_ko_top_hit in json_ko_top_hits:
                        ko_id = json_ko_top_hit['ko_id']
                        for reaction_index, chain_gene_ko_ids in enumerate(gapped_chain_ko_ids):
                            if ko_id in chain_gene_ko_ids:
                                reaction_indices.append(reaction_index)
                    if not reaction_indices:
                        # The gene is not linked to any reactions in the chain.
                        continue
                    json_gapped_chain = {}
                    json_gapped_chains.append(json_gapped_chain)
                    json_gapped_chain['gapped_chain_index'] = gapped_chain_index
                    json_gapped_chain['reaction_count'] = len(gapped_chain_ko_ids)
                    json_gapped_chain['reaction_indices_in_chain'] = sorted(set(reaction_indices))

        return json_syntenous_regions

    def make_json_syntenous_region(self, syntenous_region: list[int]) -> dict:
        """
        Make a JSON-formatted dictionary for a syntenous region.

        Various attributes are left to be filled out.

        Parameters
        ==========
        syntenous_region : list[int]
            The syntenous region is recorded as a list of row indices of full-length gene calls in
            the 'genes_in_contigs_df' attribute.

        Returns
        =======
        dict
            JSON-formatted dictionary of information on the syntenous region.
        """
        # Get information on the genes in the syntenous region.
        gene_rows = self.genes_in_contigs_df.loc[syntenous_region]
        first_gene_row = gene_rows.iloc[0]

        # Record information regarding the whole syntenous region. The unique ID of the region is
        # left to be filled out. IDs are meant to identify regions related to a particular gap.
        json_syntenous_region = {}
        json_syntenous_region['id'] = None
        json_syntenous_region['contig_id'] = first_gene_row.contig
        json_syntenous_region['gene_orientation'] = first_gene_row.direction

        # Record information on each gene in the region.
        json_full_genes: list[dict] = []
        json_syntenous_region['full_genes'] = json_full_genes
        for gene_index, gene_row in enumerate(gene_rows.itertuples()):
            json_full_gene = {}
            json_full_genes.append(json_full_gene)

            # Record the identity and position of the gene. Whether the gene is a candidate to fill
            # a particular gap is left to be filled out.
            json_full_gene['gene_index'] = gene_index
            json_full_gene['is_gap_candidate'] = None
            json_full_gene['start_position'] = gene_row.start
            json_full_gene['stop_position'] = gene_row.stop
            json_full_gene['gene_callers_id'] = gcid = gene_row.gene_callers_id
            json_full_gene['gene_call_source'] = (
                f"{gene_row.source}"
                f"{'_' + gene_row.version if gene_row.version != 'unknown' else ''}"
            )

            # Record the top KO hits to the gene.
            json_ko_top_hits: list[dict] = []
            json_full_gene['ko_top_hits'] = json_ko_top_hits
            ko_df = self.top_gene_kos_df[self.top_gene_kos_df['gene_callers_id'] == gcid]
            for ko_row in ko_df.itertuples():
                json_ko_top_hits['ko_id'] = ko_id = ko_row.accession
                json_ko_top_hits['ko_name'] = ko_row.function
                json_ko_top_hits['e_value'] = ko_row.e_value

            # Record the top COG hits to the gene.
            json_cog_top_hits: list[dict] = []
            json_full_gene['cog_top_hits'] = json_cog_top_hits
            cog_df = self.gene_cogs_df[self.gene_cogs_df['gene_callers_id'] == gcid]
            for cog_row in cog_df.itertuples():
                json_cog_top_hits['cog_id'] = cog_row.accession
                json_cog_top_hits['cog_name'] = cog_row.function
                json_cog_top_hits['e_value'] = cog_row.e_value

            # Record whether the gene is linked to the pathway through any of its top KO hits.
            json_full_gene['is_in_pathway'] = self.gene_in_pathway[gcid]

            # It is left to record the gapped chains that contain the gene, including the number of KGML
            # reactions in the chain, and the indices of KGML reactions linked to the gene via
            # top KO hits.
            json_full_gene['gapped_chains'] = []

        return json_syntenous_region

    def get_ko_info(self, ko_id: str) -> list[dict]:
        """
        Get information on KO annotations associated with genes.

        Parameters
        ==========
        ko_id : str
            KO ID.

        Returns
        =======
        list[dict]
            JSON-formatted dictionary for each gene associated with the KO.
        """
        # Find the genes with lower-ranked hits to the KO and the genes with hits to COGs that map
        # to the KO.
        json_ko_gene_hits = self.get_lower_scoring_ko_hit_info(ko_id)
        json_cog_gene_hits = self.get_ko_via_cog_hit_info(ko_id)

        if not json_ko_gene_hits and not json_cog_gene_hits:
            # No genes are associated with the KO.
            return {}

        # Combine gene records for the two types of KO associations. First map GCID to records.
        gcids: set[int] = set(
            [json_ko_hit['gene_callers_id'] for json_ko_hit in json_ko_gene_hits] +
            [json_cog_hit['gene_callers_id'] for json_cog_hit in json_cog_gene_hits]
        )
        gcid_records = {gcid: {'ko': None, 'cog': None} for gcid in gcids}
        for json_ko_gene_hit in json_ko_gene_hits:
            gcid_records[json_ko_gene_hit['gene_callers_id']]['ko'] = json_ko_gene_hit
        for json_gene_cog_hit in json_cog_gene_hits:
            gcid_records[json_gene_cog_hit['gene_callers_id']]['cog'] = json_gene_cog_hit

        # Make new gene records, copying data from the two record types.
        json_candidate_genes: list[dict] = []
        for gcid, records in gcid_records.items():
            json_candidate_gene = {'gene_callers_id': gcid}

            # Record information on KO hits.
            json_ko_gene_hit = records['ko']
            if json_ko_gene_hit:
                # The gene has lower-ranking hits to the KO.
                json_candidate_gene['gap_ko_hits'] = json_ko_gene_hit['gap_ko_hits']
                json_candidate_gene['other_ko_hits'] = json_ko_gene_hit['other_ko_hits']
            else:
                # The gene does not have lower-ranking hits to the KO, but is associated with the KO
                # by means of COG hits.
                json_candidate_gene['gap_ko_hits'] = []
                # Find and record other KO hits to the gene.
                other_gene_hits_df = self.gene_kos_df[self.gene_kos_df['gene_callers_id'] == gcid]
                assert ko_id not in other_gene_hits_df['accession']
                json_other_ko_hits: list[dict] = []
                json_candidate_gene['other_ko_hits'] = json_other_ko_hits
                for other_gene_hit_row in other_gene_hits_df.itertuples():
                    json_other_ko_hit = {}
                    json_other_ko_hits.append(json_other_ko_hit)
                    json_other_ko_hit['ko_id'] = other_gene_hit_row.accession
                    json_other_ko_hit['ko_name'] = other_gene_hit_row.function
                    json_other_ko_hit['e_value'] = other_gene_hit_row.e_value

            # Record information on associations to the gap KO via COG hits.
            json_cog_gene_hit = records['cog']
            if json_cog_gene_hit:
                json_candidate_gene['gap_kos_via_cog_top_hits'] = json_cog_gene_hit[
                    'gap_kos_via_cog_top_hits'
                ]
            else:
                json_candidate_gene['gap_kos_via_cog_top_hits'] = []

            # Find and record COG hits.
            cog_hits_df = self.gene_cogs_df[self.gene_cogs_df['gene_callers_id'] == gcid]
            json_cog_top_hits: list[dict] = []
            json_candidate_gene['cog_top_hits'] = json_cog_top_hits
            for cog_hit_row in cog_hits_df.itertuples():
                json_cog_top_hit = {}
                json_cog_top_hits.append(json_cog_top_hit)
                json_cog_top_hit['cog_id'] = cog_hit_row.accession
                json_cog_top_hit['cog_name'] = cog_hit_row.function
                json_cog_top_hit['e_value'] = cog_hit_row.e_value

        return json_candidate_genes

    def get_lower_scoring_ko_hit_info(self, ko_id: str) -> list[dict]:
        """
        Get information on all gene hits to the KO.

        Parameters
        ==========
        ko_id : str
            KO ID.

        Returns
        =======
        list[dict]
            JSON-formatted dictionary for each gene associated with the KO, formatted for use in the
            dictionary of gap-filling evidence.
        """
        # Get all gene hits to the KO.
        gene_hits_df = self.gene_kos_df[self.gene_kos_df['accession'] == ko_id]

        if not len(gene_hits_df):
            # No genes hit the KO.
            return {}

        # Check proper data formatting. Get the KO name.
        assert gene_hits_df['gene_callers_id'].nunique() == len(gene_hits_df)
        assert gene_hits_df['ko_name'].nunique() == 1
        ko_name = gene_hits_df.iloc[0]['function']

        # Record information on each gene.
        json_candidate_genes: list[dict] = []
        for gene_hit_row in gene_hits_df.itertuples():
            json_candidate_gene = {}
            json_candidate_genes.append(json_candidate_gene)
            json_candidate_gene['gene_callers_id'] = gcid = gene_hit_row.gene_callers_id

            json_gap_ko_hits: list[dict] = []
            json_candidate_gene['gap_ko_hits'] = json_gap_ko_hits
            json_gap_ko_hit = {}
            json_gap_ko_hits.append(json_gap_ko_hit)
            json_gap_ko_hit['ko_id'] = ko_id
            json_gap_ko_hit['ko_name'] = ko_name
            json_gap_ko_hit['e_value'] = gene_hit_row.e_value

            # Record information on other KOs that hit the gene.
            other_gene_hits_df = self.gene_kos_df[self.gene_kos_df['gene_callers_id'] == gcid]
            other_gene_hits_df = other_gene_hits_df[other_gene_hits_df['accession'] != ko_id]

            json_other_ko_hits: list[dict] = []
            for other_gene_hit_row in other_gene_hits_df.itertuples():
                json_other_ko_hit = {}
                json_other_ko_hits.append(json_other_ko_hit)
                json_other_ko_hit['ko_id'] = other_gene_hit_row.accession
                json_other_ko_hit['ko_name'] = other_gene_hit_row.function
                json_other_ko_hit['e_value'] = other_gene_hit_row.e_value

        return json_candidate_genes

    def get_ko_via_cog_hit_info(self, ko_id: str) -> list[dict]:
        """
        Get information on genes with COG annotations that map to the KO via the KEGG binary
        relations table.

        Parameters
        ==========
        ko_id : str
            KO ID.

        Returns
        =======
        list[dict]
            JSON-formatted dictionary for each gene associated with the KO via COG, formatted for
            use in the dictionary of gap-filling evidence. If the KO does not map to any COGs, or if
            no genes hit the COGs, then an empty dictionary is returned.
        """
        # Find COGs that KEGG deems equivalent to the KO, and find gene hits to the COGs.
        try:
            equivalent_cog_ids = self.ko_cogs[ko_id]
        except KeyError:
            # The KO does not have equivalent COGs.
            return {}
        gene_hits_df = self.gene_cogs_df[self.gene_cogs_df['accession'].isin(equivalent_cog_ids)]

        if not len(gene_hits_df):
            # No genes are annotated by equivalent COGs.
            return {}

        # Get the KO name, checking proper data formatting.
        ko_name_df = self.ko_list_df[self.ko_list_df['accession'] == ko_id]
        assert len(ko_name_df) == 1
        ko_name = ko_name_df.iloc[0]['function']

        # Record information on each gene.
        json_candidate_genes: list[dict] = []
        for gcid, hits_df in gene_hits_df.groupby('gene_callers_id'):
            json_candidate_gene = {}
            json_candidate_genes.append(json_candidate_gene)
            json_candidate_gene['gene_callers_id'] = gcid

            json_gap_kos_via_cog_top_hits: list[dict] = []
            json_candidate_gene['gap_kos_via_cog_top_hits'] = json_gap_kos_via_cog_top_hits
            json_gap_ko_via_cog_top_hits = {}
            json_gap_kos_via_cog_top_hits.append(json_gap_ko_via_cog_top_hits)
            json_gap_ko_via_cog_top_hits['ko_id'] = ko_id
            json_gap_ko_via_cog_top_hits['ko_name'] = ko_name
            json_associated_cog_hits: list[dict] = []
            json_gap_ko_via_cog_top_hits['cog_ids'] = hits_df['accession'].tolist()

        return json_candidate_genes

    def make_json_adjacent_pathway_genes(
        self,
        json_syntenous_region: dict,
        gap_gene_index: int,
        gap_chain_gcids: list[int]
    ) -> list[dict]:
        """
        Make JSON-formatted dictionaries containing information on the genes in a syntenous region
        that directly surround a candidate gap-filling gene and that have top KO annotations found
        in the pathway.

        Dictionaries for adjacent genes are ordered along the contig, from genes before to genes
        after the gap gene.

        Parameters
        ==========
        json_syntenous_region : dict
            JSON-formatted dictionary of information on the syntenous region.

        gap_gene_index : int
            Index of the candidate gap-filling gene in the syntenous region.

        gap_chain_gcids : list[int]
            Gene caller IDs of all genes linked to KGML reactions in the chain by top KO
            annotations.

        Returns
        =======
        list[dict]
            JSON-formatted dictionaries containing information on adjacent genes in the pathway.
        """
        json_adjacent_pathway_genes: list[dict] = []
        # Investigate adjacent genes moving from the gap gene toward the start of the contig.
        for relative_gene_index, json_full_gene in enumerate(
            json_syntenous_region['full_genes'][: gap_gene_index][::-1], 1
        ):
            is_gene_in_pathway = json_full_gene['is_in_pathway']
            if not is_gene_in_pathway:
                # The gene's top KO annotations are not in the pathway, so stop the traversal.
                break
            json_adjacent_pathway_gene = {}
            json_adjacent_pathway_genes.append(json_adjacent_pathway_gene)
            json_adjacent_pathway_gene['gene_index'] = gap_gene_index - relative_gene_index
            json_adjacent_pathway_gene[
                'is_gene_in_gapped_chain'
            ] = json_full_gene['gene_callers_id'] in gap_chain_gcids

        # Investigate adjacent genes moving from the gap gene toward the end of the contig.
        for relative_gene_index, json_full_gene in enumerate(
            json_syntenous_region['full_genes'][gap_gene_index + 1: ], 1
        ):
            is_gene_in_pathway = json_full_gene['is_in_pathway']
            if not is_gene_in_pathway:
                # The gene's top KO annotations are not in the pathway, so stop the traversal.
                break
            json_adjacent_pathway_gene = {}
            json_adjacent_pathway_genes.append(json_adjacent_pathway_gene)
            json_adjacent_pathway_gene['gene_index'] = gap_gene_index + relative_gene_index
            json_adjacent_pathway_gene[
                'is_gene_in_gapped_chain'
            ] = json_full_gene['gene_callers_id'] in gap_chain_gcids

        return json_adjacent_pathway_genes
