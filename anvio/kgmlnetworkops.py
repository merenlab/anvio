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
    'gaps', 'aliased_modelseed_compounds', 'network_kos', 'aliased_modelseed_reactions',
    'network_kos', and 'aliased_modelseed_reactions'.

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
    """
    kgml_compound_entries: list[kgml.Entry] = field(default_factory=list)
    is_consumed: bool = None
    kgml_reactions: list[kgml.Reaction] = field(default_factory=list)
    kgml_reaction_directions: list[bool] = field(default_factory=list)
    gaps: list[bool] = field(default_factory=list)
    aliased_modelseed_compounds: list[tuple[rn.ModelSEEDCompound]] = field(default_factory=list)
    network_kos: list[tuple[rn.KO]] = field(default_factory=list)
    aliased_modelseed_reactions: list[tuple[rn.ModelSEEDReaction]] = field(default_factory=list)
    is_consumption_terminus: bool = None
    is_production_terminus: bool = None
    consumption_reversibility_range: tuple[int, int] = None
    production_reversibility_range: tuple[int, int] = None

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

    network_keggcpd_id_to_modelseed_compounds : dict[str, list[rn.ModelSEEDCompound]], {}
        Map the IDs of KEGG compounds (not KGML compound IDs) in the reaction network to aliased
        ModelSEED compounds.

    network_keggcpd_ids_in_pathway : list[str]
        KEGG compound IDs in the pathway from the reaction network. KEGG compounds in the reaction
        network are selected to ensure that they are linked to KO annotations by ModelSEED reactions
        with KEGG reaction aliases, not EC number aliases, since EC numbers associated with KOs can
        alias a large number of reactions of questionable validity for the enzyme.

    compound_fate : Literal['consume', 'produce', 'both'], 'both'
        Seek chains that consume or produce compounds in the network. If 'consume' or 'produce',
        only consumption or production chains are sought, respectively. If 'both', both consumption
        and production chains are sought. Chains that only contain reversible reactions are
        therefore found in both directions with 'both'.

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
            Contains arguments. See the class docstring for more information on arguments set as
            attributes, including default values. The only required argument is kegg_pathway_number.

            kegg_pathway_number : str

            contigs_db_path : str

            network : anvio.reactionnetwork.GenomicNetwork

            compound_fate : Literal['consume', 'produce', 'both']

            max_reactions : int

            keep_intermediate_chains : bool

            max_gaps : int

            allow_terminal_gaps : bool

            allow_alternative_reaction_gaps : bool

            run : anvio.terminal.Run

            verbose : bool
        """
        A = lambda x, y: args.__dict__[x] if x in args.__dict__ else y

        self.kegg_pathway_number: str = args.kegg_pathway_number

        self.contigs_db_path: str = A(args.contigs_db, None)
        self.network: rn.GenomicNetwork = A(args.network, None)
        self.verbose = A('verbose', False)
        if self.contigs_db_path is not None and not self.network:
            constructor = rn.Constructor()
            self.network = constructor.load_contigs_database_network(
                self.contigs_db_path, quiet=not self.verbose
            )

        self.compound_fate: str = A('compound_fate', 'both')
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

            if self.compound_fate == 'both':
                consumption_options = [True, False]
            elif self.compound_fate == 'consume':
                consumption_options = [True]
            elif self.compound_fate == 'produce':
                consumption_options = [False]
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

                    keggrn_ids: list[str] = []
                    for candidate_keggrn_id in kgml_reaction.name.split():
                        if candidate_keggrn_id[:3] != 'rn:':
                            continue
                        keggrn_id = candidate_keggrn_id[3:]

                    modelseed_reaction_ids: list[str] = []
                    for ko in network_kos:
                        for keggrn_id in keggrn_ids:
                            try:
                                modelseed_reaction_ids += ko.kegg_reaction_aliases[keggrn_id]
                            except KeyError:
                                continue

                    modelseed_reactions: list[rn.ModelSEEDReaction] = []
                    for modelseed_reaction_id in sorted(set(modelseed_reaction_ids)):
                        modelseed_reactions.append(self.network.reactions[modelseed_reaction_id])

                # Recurse on each KGML compound on the other side of the reaction.
                for next_kgml_compound_id in next_kgml_compound_ids:
                    if len(current_chain.kgml_compound_entries) > 2:
                        if next_kgml_compound_id == current_chain.kgml_compound_entries[-2].id:
                            # Avoid backtracking to the previous KGML compound in the chain.
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

                    # For each compound in a terminal chain, a subchain that extends to the terminus
                    # is also generated as a "candidate terminal chain" and is ignored.
                    if terminal_chains:
                        if [c.id for c in candidate_terminal_chain.kgml_compound_entries] == [
                            c.id for c in terminal_chains[-1].kgml_compound_entries[
                                :len(candidate_terminal_chain.kgml_compound_entries)
                            ]
                        ]:
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
                    for kgml_reaction in candidate_terminal_chain.kgml_reactions:
                        if kgml_reaction.type == 'irreversible':
                            break
                        i += 1
                    if candidate_terminal_chain.is_consumed:
                        candidate_terminal_chain.consumption_reversibility_range = (0, i + 1)
                    else:
                        candidate_terminal_chain.production_reversibility_range = (0, i + 1)

                    i = 0
                    for kgml_reaction in candidate_terminal_chain.kgml_reactions[::-1]:
                        if kgml_reaction.type == 'irreversible':
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
    Records how chains with fewer genomic gaps relate to a chain with more gaps.

    Attributes
    ==========
    gappy_chain : Chain, None
        Chain with more gaps.

    ungappy_chains : list[Chain], []
        Chains with fewer gaps.

    overlaps : list[tuple[tuple[int, int]]], []
        This list has a tuple item for each ungappy chain. There is an inner tuple for each reaction
        shared between the gappy and ungappy chain, with the first item of the inner tuple being the
        index of the reaction in the gappy chain and the second item being the index of the reaction
        in the ungappy chain.

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
        Each of the gaps is a reaction.

    gap_chain_relations : list[GapChainRelations], []
        This list contains an item per gappy chain that has the set of gaps represented here.
    """
    gap_kgml_reactions: list[kgml.Reaction] = field(default_factory=list)
    gap_chain_relations: list[GapChainRelations] = field(default_factory=list)

class GapAnalyzer:
    """
    Analyze chains of KGML compounds linked by reactions, some of which are designated as gaps.
    Compare two sets of chains found from the same KGML source but with the set of gappy chains
    permitting more gaps than ungappy chains.

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
        self.gap_relations = self.get_gap_relations()

    def get_gap_relations(self) -> dict[tuple[str], SharedGaps]:
        """
        Get information associated with each set of gaps that exists in one or more gappy chains,
        including relationships between gappy chains and overlapping ungappy chains.

        Returns
        =======
        dict[tuple[str], SharedGaps]
            Keys are the KGML reaction IDs of sets of gaps in gappy chains. Values are information
            associated with each set of gaps.
        """
        gap_relations = {}
        for gappy_chain in self.gappy_chains:
            if not any(gappy_chain.gaps):
                # The chain has no gaps. (Gapless along with gapped chains can be returned when
                # seeking chains allowing for gaps.)
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

            gappy_chain_kgml_reaction_ids = [
                kgml_reaction.id for kgml_reaction in gappy_chain.kgml_reactions
            ]

            # Record information on ungappy chains overlapping with the gappy chain.
            overlaps: list[tuple[tuple[int, int]]] = []
            ungappy_chains: list[Chain] = []
            for ungappy_chain in self.ungappy_chains:
                if gappy_chain.is_consumed != ungappy_chain.is_consumed:
                    # Ignore ungappy chains with reactions running in the opposite direction to the
                    # gappy chain.
                    continue

                ungappy_chain_kgml_reaction_ids = [
                    kgml_reaction.id for kgml_reaction in ungappy_chain.kgml_reactions
                ]

                if gappy_chain_kgml_reaction_ids == ungappy_chain_kgml_reaction_ids:
                    # Ignore gappy and ungappy chains with identical reactions. (The same chains can
                    # be returned when seeking chains allowing for more and fewer gaps.)
                    continue

                overlap: list[tuple[int, int]] = []
                for i, gappy_chain_kgml_reaction_id in enumerate(gappy_chain_kgml_reaction_ids):
                    for j, ungappy_chain_kgml_reaction_id in enumerate(
                        ungappy_chain_kgml_reaction_ids
                    ):
                        if gappy_chain_kgml_reaction_id == ungappy_chain_kgml_reaction_id:
                            # Record the index of the reaction in the gappy and ungappy chains,
                            # respectively.
                            overlap.append((i, j))
                if not overlap:
                    # Ignore ungappy chains that do not overlap with the gappy chain.
                    continue
                ungappy_chains.append(ungappy_chain)
                # Each ungappy chain's overlap with the gappy chain is represented by a tuple of
                # tuples.
                overlaps.append(tuple(overlap))

            # Sort ungappy chains associated with the gappy chain by index of first overlapping
            # reaction in the gappy chain.
            sorted_overlaps = sorted(overlaps, key=lambda overlap: overlap[0][0])
            sorted_ungappy_chains: list[Chain] = []
            for overlap in sorted_overlaps:
                overlap: tuple[tuple[int]]
                sorted_ungappy_chains.append(ungappy_chains[overlaps.index(overlap)])
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

                # Ignore gappy chains that, beside gap reactions, contain the same reactions as an
                # ungappy chain. Gaps in the gappy chain create "shortcuts" in what is otherwise the
                # same chain.
                kgml_reaction_ids_absent_gaps = set([
                    kgml_reaction.id for kgml_reaction in gappy_chain.kgml_reactions
                    if kgml_reaction.id not in gap_kgml_reaction_ids
                ])
                for ungappy_chain in gap_chain_relations.ungappy_chains:
                    if not kgml_reaction_ids_absent_gaps.difference(
                        set([kgml_reaction.id for kgml_reaction in ungappy_chain.kgml_reactions])
                    ):
                        is_difference_gaps = True
                        break
                else:
                    is_difference_gaps = False
                if is_difference_gaps:
                    continue

                segments = [
                    tuple([indices[0] for indices in overlap])
                    for overlap in gap_chain_relations.overlaps
                ]

                # Find unique segments, or the longest ungappy chains that encompass the reactions
                # of the gappy chain beside its gaps not in ungappy chains.
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
    def __init__(self, args: Namespace):
        A = lambda x, y: args.__dict__[x] if x in args.__dict__ else y

        self.walker: KGMLNetworkWalker = A(args.walker, None)
        if self.walker is None:
            raise ConfigError("A KGML network walker is required as the 'walker' argument.")
        if self.walker.contigs_db_path is None:
            raise ConfigError(
                "The KGML network walker must be associated with a contigs database."
            )
        if self.walker.network is None:
            raise ConfigError(
                "The KGML network walker must be associated with a reaction network."
            )

        self.gap_analyzer: GapAnalyzer = A(args.gap_analyzer, None)
        # if self.gap_analyzer is None:
        #     raise ConfigError(
        #         "The gap analyzer used to compare gapped and ungapped chains is required as the "
        #         "'gap_analyzer' argument."
        #     )

        self.ko_cog_path: str = A(args.ko_cog, None)

        self.pathway_ortholog_entries = self.walker.kgml_ko_pathway.get_entries(
            entry_type='ortholog'
        )

        self.contigs_db = ContigsDatabase(self.contigs_db_path)
        function_sources = (
            self.contigs_db.meta['gene_function_sources']
            if self.meta['gene_function_sources'] else []
        )
        if 'KOfam' not in function_sources:
            raise ConfigError(
                "Genes of the contigs database should have been annotated with KOs, but there is "
                "no 'KOfam' function source."
            )
        # if 'COG20_FUNCTION' not in function_sources and self.ko_cog_path:
        #     raise ConfigError(
        #         "KO and COG20 annotations are to be compared, so genes of the contigs database "
        #         "should have been annotated with COG20 functions, but there is no 'COG20_FUNCTION' "
        #         "source."
        #     )
        self.genes_in_contigs_df = self.contigs_db.db.get_table_as_dataframe('genes_in_contigs')
        # Sort genes by contig start position.
        self.genes_in_contigs_df = self.genes_in_contigs_df.sort_values(
            ['contig', 'start']
        ).reset_index()
        gene_functions_df = self.contigs_db.db.get_table_as_dataframe('gene_functions')
        self.gene_kos_df = gene_functions_df[gene_functions_df['source'] == 'KOfam']
        if 'COG20_FUNCTION' not in function_sources:
            self.gene_cogs_df = None
        else:
            self.gene_cogs_df = gene_functions_df[gene_functions_df['source'] == 'COG20_FUNCTION']

        if self.gap_analyzer:
            self.gap_key_ranks = self.gap_analyzer.rank_gaps()

        if self.ko_cog_path is None:
            self.cog_kos = None
        else:
            cog_kos: dict[str, list[str]] = {}
            ko_cog_df = pd.read_csv(self.ko_cog_path, sep='\t')
            ko_cog_df.columns = ['ko', 'cog']
            for row in ko_cog_df.itertuples():
                for cog_id in row.cog[5: -1].split():
                    try:
                        cog_kos[cog_id].append(row.ko)
                    except KeyError:
                        cog_kos[cog_id] = [row.ko]
            self.cog_kos = cog_kos

    def eval_gap_ko(self, ko_id: str) -> Union[dict, None]:
        gene_hits_df = self.gene_kos_df[self.gene_kos_df['accession'] == ko_id]

        if not len(gene_hits_df):
            return None

        other_gene_hits_df = self.gene_kos_df[
            self.gene_kos_df['gene_callers_id'].isin(gene_hits_df['gene_callers_id'])
        ]
        other_gene_hits_df = other_gene_hits_df[other_gene_hits_df['accession'] != ko_id]

        json_obj = {}
        if len(gene_hits_df):
            json_obj['ko_id'] = ko_id
            json_obj['ko_name'] = gene_hits_df.iloc[0]['function']

        if gcid_cog_df is None:
            gcid_cog_df = None
            equivalent_cog_ids = None
        else:
            gcid_cog_df: dict[int, pd.DataFrame] = {}
            for gcid in gene_hits_df['gene_callers_id']:
                gcid_cog_df[gcid] = self.get_cog_hits(gcid)

            if self.cog_kos is None:
                equivalent_cog_ids = None
            else:
                equivalent_cog_ids: list[str] = []
                for gcid, cog_df in gcid_cog_df.items():
                    for cog_id in cog_df['accession'].unique():
                        equivalent_ko_ids = self.cog_kos[cog_id]
                        for equivalent_ko_id in equivalent_ko_ids:
                            if equivalent_ko_id != ko_id:
                                continue
                            equivalent_cog_ids.append(cog_id)

        json_gene_hits: list[dict] = []
        json_obj['annotated_genes'] = json_gene_hits
        for gene_hit_row in gene_hits_df.itertuples():
            json_gene_hit = {}
            json_gene_hits.append(json_gene_hit)

            json_gene_hit['gene_callers_id'] = gcid = gene_hit_row.gene_callers_id
            json_gene_hit['e_value'] = gene_hit_row.e_value

            json_other_kos: list[dict] = []
            json_gene_hit['other_kos'] = json_other_kos
            for other_gene_hit_row in other_gene_hits_df[
                other_gene_hits_df['gene_callers_id'] == gcid
            ].itertuples():
                json_other_ko = {}
                json_other_kos.append(json_other_ko)
                json_other_ko['other_ko_id'] = other_gene_hit_row.accession
                json_other_ko['other_ko_name'] = other_gene_hit_row.function
                json_other_ko['other_e_value'] = other_gene_hit_row.e_value

            if gcid_cog_df is None:
                json_cog_annotations = None
            else:
                json_cog_annotations: list[dict] = []
                cog_df = gcid_cog_df[gcid]
                for cog_row in cog_df.itertuples():
                    json_cog_annotation = {}
                    json_cog_annotation['cog20_id'] = cog_id = cog_row.accession
                    json_cog_annotation['cog20_name'] = cog_row.function
                    if self.cog_kos is None:
                        json_cog_annotation['is_equivalent_to_ko'] = None
                    elif cog_id in equivalent_cog_ids:
                        json_cog_annotation['is_equivalent_to_ko'] = True
                    else:
                        json_cog_annotation['is_equivalent_to_ko'] = False
            json_gene_hit['cog20_annotations'] = json_cog_annotations

        return json_obj

    def get_cog_hits(self, gcid: int) -> Union[pd.DataFrame, None]:
        if self.gene_cogs_df is None:
            return None

        gene_hits_df = self.gene_cogs_df[self.gene_cogs_df['gene_callers_id'] == gcid]

        expanded_gene_hits = []
        for row in gene_hits_df.itertuples(index=False):
            if '!!!' not in row.accession:
                expanded_gene_hits.append(row)
                continue
            for cog_id, cog_name in zip(row.accession.split('!!!'), row.function.split('!!!')):
                new_row = row._replace(accession=cog_id, function=cog_name)
                expanded_gene_hits.append(new_row)
        expanded_gene_hits_df = pd.DataFrame(expanded_gene_hits, columns=gene_hits_df.columns)

        return expanded_gene_hits_df

    def eval_kgml_reaction(self, kgml_reaction_id):
        pass

    def eval_pathway(self):
        pass

    def eval_chain(self, chains):
        pass

    def eval_gap_in_chains(self, kgml_reaction_id, chains):
        pass

    def get_chain_gcids(self, chain: Chain) -> list[int]:
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

            # The entry corresponds to one or more KOs: look for KOs that annotate genes.
            network_kos: list[rn.KO] = []
            for name_substring in kgml_entry.name.split():
                if name_substring[:3] != 'ko:':
                    continue

                ko_id = name_substring[3:]
                try:
                    network_kos.append(self.walker.network.kos[ko_id])
                except KeyError:
                    # The KO does not annotate a gene.
                    continue

            # Record IDs of genes annotated by the KOs.
            for ko in network_kos:
                for gene in self.walker.network.genes.values():
                    if ko.id in gene.ko_ids:
                        chain_gcids.append(gene.gcid)
        return chain_gcids

    def get_syntenous_regions(self, gcids: list[int]) -> list[list[int]]:
        syntenous_regions: list[list[int]] = []
        for gcid in gcids:
            # Search around the gene. The search stops in either direction when a gene in the
            # opposite orientation is found or the first or last gene in the contig is reached.
            row = self.genes_in_contigs_df[
                self.genes_in_contigs_df['gene_callers_id'] == gcid
            ].squeeze()
            row_index = row.name
            direction = row['direction']

            syntenous_region: list[int] = []
            # Search preceding genes.
            for row_index in range(row_index, -1, -1):
                prior_row = self.genes_in_contigs_df.loc[row_index]
                if prior_row['direction'] != direction:
                    break
                if prior_row['partial'] == 1:
                    # Ignore partial gene calls.
                    continue
                syntenous_region.append(row_index)
            # Search succeeding genes.
            for row_index in range(row_index + 1, len(self.genes_in_contigs_df)):
                next_row = self.genes_in_contigs_df.loc[row_index]
                if next_row['direction'] != direction:
                    break
                if next_row['partial'] == 1:
                    continue
                syntenous_region.append(row_index)

            syntenous_region = sorted(syntenous_region)
            if syntenous_region in syntenous_regions:
                # The syntenous region has already been encountered in considering another gene.
                continue
            syntenous_regions.append(syntenous_region)
        return syntenous_regions

    """
    {
        "pathway_number": <pathway number>, # Class to evaluate each gap
        "pathway_name": <pathway name>,
        "gaps": [
            {
                "kgml_reaction_id": <KGML reaction ID>, # have a function to return result from this point
                "ko_hits": [
                    {
                        "ko_id": <KO ID>, # have a function to return result from this point
                        "ko_name": <KO name>,
                        "gene_hits": [
                            {
                                "gene_callers_id": <GCID>,
                                "e_value": <E value>,
                                "other_kos": [
                                    {
                                        "other_ko_id": <KO ID>,
                                        "other_ko_name": <KO name>,
                                        "other_e_value": <E value>
                                    },
                                    ...
                                ],
                                "cog20_annotations": [ or null
                                    {
                                        "cog20_id": <COG20 ID>,
                                        "cog20_name": <COG20 name>,
                                        "is_equivalent_to_ko": <true or false> or null,
                                    },
                                    ...
                                ],
                                "gene_index_in_syntenous_region": <gene index in syntenous region>,
                                "syntenous_region_in_pathway": {
                                    "contig_id": <contig ID>,
                                    "gene_orientation": <"f" or "r">,
                                    "gene_count": <gene count>,
                                    "surrounding_syntenous_genes": [
                                        {
                                            "gene_callers_id": <GCID>,
                                            "start_position": <start position>,
                                            "stop_position": <stop position>,
                                            "gene_call_source": <gene call source>,
                                            "ko_id": <KO ID or null>,
                                            "ko_name": <KO name or null>,
                                            "cog20_id": <COG20 ID or null>,
                                            "cog20_name": <COG20 ID or null>,
                                            "is_in_pathway": <true or false>,
                                            "chains": [
                                                {
                                                    "chain_id": <chain ID>,
                                                    "reaction_count": <reaction count>,
                                                    "reaction_indices_in_chain": [<reaction index in chain>, ...]
                                                },
                                                ...
                                            ]
                                        },
                                        ...
                                    ]
                                }
                            },
                            ...
                        ]
                    },
                    ...
                ]
            },
            ...
        ],
        "chains": [
            {
                "id": <ID>,
            },
            ...
        ]
    }
    """

    # def evaluate(self):
    #     json = {}
    #     json['kegg_pathway_number'] = self.walker.kegg_pathway_number
    #     json['kegg_pathway_name'] = self.walker.kgml_ko_pathway.name
    #     json['gaps'] = []
    #     # have separate evaluate gap function for arbitrary gap
    #     for gap_kgml_reaction_ids in self.gap_key_ranks:
    #         json_gap = {}
    #         shared_gaps = self.gap_analyzer.gap_relations[gap_kgml_reaction_ids]

    #         for gap_index in range(len(gap_kgml_reaction_ids)):
    #             gap_kgml_reaction = shared_gaps.gap_kgml_reactions[gap_index]

    #             for entry in self.pathway_ortholog_entries:
    #                 if gap_kgml_reaction.id == entry.id:
    #                     gap_entry = entry
    #                     break
    #             else:
    #                 raise AssertionError(
    #                     f"Every KGML reaction in the KO type file should correspond to an ortholog "
    #                     f"type entry. The reaction with KGML ID '{gap_kgml_reaction.id}' does not "
    #                     "have a corresponding entry with the same ID in pathway "
    #                     f"'{self.walker.kegg_pathway_number}'."
    #                 )

    #             gap_ko_gene_hits: dict[str, pd.DataFrame] = {}
    #             for candidate_ko_id in gap_entry.name.split():
    #                 if candidate_ko_id[:3] != 'ko:':
    #                     continue
    #                 ko_id = candidate_ko_id[3:]
    #                 gene_hits = self.gene_kos_df[self.gene_kos_df['accession'] == ko_id]
    #                 if len(gene_hits):
    #                     gap_ko_gene_hits[ko_id] = gene_hits

    #             if gap_ko_gene_hits:
    #                 pass

    #         for gap_chain_relations in shared_gaps.gap_chain_relations:
    #             gappy_chain = gap_chain_relations.gappy_chain
    #             # Find genes annotated by KOs in the chain.
    #             chain_gcids = get_chain_gcids(gappy_chain)
    #             chain_syntenous_regions = get_syntenous_regions(chain_gcids)

    #         json['gaps'] = json_gap

# class GapFiller:
#     def __init__(self, args: Namespace) -> None:
#         A = lambda x, y: args.__dict__[x] if x in args.__dict__ else y
#         self.gapped_chains: list[Chain] = A(args.gapped_chains, None)
#         self.ungapped_chains: list[Chain] = A(args.ungapped_chains, None)
#         self.walker: KGMLNetworkWalker = A(args.walker, None)
#         self.contigs_db_path: str = A(args.contigs_db, None)
#         self.ko_cog_path: str = A(args.ko_cog, None)

#         self.gap_analyzer = GapAnalyzer(self.gapped_chains, self.ungapped_chains)
#         self.chain_gaps = self.gap_analyzer.rank_gaps()

#         self.pathway_ortholog_entries = self.walker.kgml_ko_pathway.get_entries(
#             entry_type='ortholog'
#         )

#         self.contigs_db = ContigsDatabase(self.contigs_db_path)
#         function_sources = (
#             self.contigs_db.meta['gene_function_sources']
#             if self.meta['gene_function_sources'] else []
#         )
#         assert 'KOfam' in function_sources
#         assert 'COG20_FUNCTION' in function_sources
#         self.genes_in_contigs_df = self.contigs_db.db.get_table_as_dataframe('genes_in_contigs')
#         self.genes_in_contigs_df = self.genes_in_contigs_df.sort_values('start').reset_index()
#         self.gene_functions_df = self.contigs_db.db.get_table_as_dataframe('gene_functions')

#         cog_kos: dict[str, list[str]] = {}
#         ko_cog_df = pd.read_csv(self.ko_cog_path, sep='\t')
#         ko_cog_df.columns = ['ko', 'cog']
#         for row in ko_cog_df.itertuples():
#             for cog_id in row.cog[5: -1].split():
#                 try:
#                     cog_kos[cog_id].append(row.ko)
#                 except KeyError:
#                     cog_kos[cog_id] = [row.ko]
#         self.cog_kos = cog_kos

#     def get_syntenous_regions(self, chain_gcids: list[int]) -> list[list[int]]:
#         # Find syntenous regions containing genes in the chain.
#         syntenous_regions: list[list[int]] = []
#         for gcid in chain_gcids:
#             # Search for the syntenous region around the GCID.
#             row = self.genes_in_contigs_df[
#                 self.genes_in_contigs_df['gene_callers_id'] == gcid
#             ].squeeze()
#             row_index = row.name
#             direction = row['direction']

#             syntenous_region: list[int] = []
#             # Search for preceding syntenous genes. The search stops when a gene in the opposite
#             # direction is found or the first gene is reached.
#             for candidate_row_index in range(row_index, -1, -1):
#                 prior_row = self.genes_in_contigs_df.loc[candidate_row_index]
#                 if prior_row['direction'] != direction:
#                     break
#                 if prior_row['partial'] == 1:
#                     # Ignore partial gene calls.
#                     continue
#                 syntenous_region.append(candidate_row_index)
#             # Search for succeeding syntenous genes.
#             for candidate_row_index in range(row_index + 1, len(self.genes_in_contigs_df)):
#                 next_row = self.genes_in_contigs_df.loc[candidate_row_index]
#                 if next_row['direction'] != direction:
#                     break
#                 if next_row['partial'] == 1:
#                     continue
#                 syntenous_region.append(candidate_row_index)

#             syntenous_region = sorted(syntenous_region)
#             if syntenous_region in syntenous_regions:
#                 # The syntenous region has already been encountered in considering another gene in
#                 # the chain.
#                 continue
#             syntenous_regions.append(syntenous_region)
#         return syntenous_regions

#     def run(self):
#         for chain in self.gapped_chains:
#             if True not in chain.gaps:
#                 continue

#             # Find genes annotated by KOs in the chain.
#             chain_gcids: list[int] = []
#             for is_gap, kgml_reaction in chain.kgml_reactions:
#                 if is_gap:
#                     continue

#                 # Find the KGML ortholog entry for the reaction.
#                 for entry in self.ortholog_entries:
#                     if entry.id == kgml_reaction.id:
#                         kgml_entry = entry
#                         break
#                 else:
#                     raise AssertionError

#                 # The entry corresponds to one or more KOs. Look for the KOs in the reaction
#                 # network.
#                 network_kos: list[rn.KO] = []
#                 for name_substring in kgml_entry.name.split():
#                     if name_substring[:3] != 'ko:':
#                         continue

#                     ko_id = name_substring[3:]
#                     try:
#                         network_kos.append(self.walker.network.kos[ko_id])
#                     except KeyError:
#                         continue

#                 # Get the annotated genes.
#                 for ko in network_kos:
#                     for gene in self.walker.network.genes.values():
#                         if ko.id in gene.ko_ids:
#                             chain_gcids.append(gene.gcid)

#             syntenous_regions = self.get_syntenous_regions(chain_gcids)



#         for gap_kgml_reaction_ids in self.chain_gaps:
#             for gap_kgml_reaction_id in gap_kgml_reaction_ids:
#                 pass
