# -*- coding: utf-8
# pylint: disable=line-too-long
"""Generate a metabolic reaction network from gene annotations."""

import os
import numpy as np
import pandas as pd

from math import gcd
from hashlib import sha1
from functools import reduce
from argparse import Namespace
from fractions import Fraction
from collections import Counter
from typing import Dict, List, Tuple

import anvio.tables as tables
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.utils import is_contigs_db
from anvio.dbops import ContigsDatabase
from anvio.dbops import ContigsSuperclass
from anvio import __file__ as ANVIO_PATH, __version__ as VERSION


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
        self.kegg_id_aliases: Tuple[str] = None
        self.charge: int = None
        self.formula: str = None

class ModelSEEDReaction:
    """Representation of a reaction in the network, with properties given by the ModelSEED Biochemistry database."""
    def __init__(self) -> None:
        self.modelseed_id: str = None
        self.modelseed_name: str = None
        self.kegg_id_aliases: Tuple[str] = None
        self.ec_number_aliases: Tuple[str] = None
        # compounds, coefficients, and compartments have corresponding elements
        self.compounds: Tuple[ModelSEEDCompound] = None
        self.coefficients: Tuple[int] = None
        self.compartments: Tuple[str] = None
        self.reversibility: bool = None

class KO:
    """Representation of a KEGG Ortholog in the network."""
    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        # map *ModelSEED reaction ID* to *ModelSEED reaction object or reaction aliases* in the
        # following dictionaries
        self.reactions: Dict[str, ModelSEEDReaction] = {}
        # Record the KEGG REACTION IDs *encoded by the KO* that are aliases of the ModelSEED
        # reaction ID. These could be a subset of the KEGG reaction aliases of the ModelSEED
        # reaction. The same is true of EC numbers.
        self.kegg_reaction_aliases: Dict[str, Tuple[str]] = {}
        self.ec_number_aliases: Dict[str, Tuple[str]] = {}

class Gene:
    """Representation of a gene in the metabolic network."""
    def __init__(self) -> None:
        self.gcid: int = None
        # KOs matching the gene
        self.kos: List[KO] = []
        # record the strength of each KO match
        self.e_values: List[float] = []

class GeneCluster:
    """Representation of a gene cluster."""
    def __init__(self) -> None:
        # genes in the gene cluster
        self.genes: List[Gene] = []

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
        # map KO ID to KO object
        self.kos: Dict[int, KO] = {}
        # map ModelSEED reaction ID to reaction object
        self.reactions: Dict[str, ModelSEEDReaction] = {}
        # map ModelSEED compound ID to compound object
        self.metabolites: Dict[str, ModelSEEDCompound] = {}
        # map KEGG REACTION ID to ModelSEED reaction IDs
        self.kegg_reactions: Dict[str, List[str]] = {}
        # map EC number to ModelSEED reaction IDs
        self.ec_numbers: Dict[str, List[str]] = {}

class GenomicNetwork(ReactionNetwork):
    """A reaction network predicted from KEGG KO and ModelSEED annotations of genes."""
    def __init__(self) -> None:
        # map gene caller ID to gene object
        self.genes: Dict[int, Gene] = {}
        self.bins: Dict[str, GeneBin] = {}
        self.collection: BinCollection = None
        super().__init__()

    def export_json(
        self,
        path: str,
        annotate_genes: tuple = ('bin', ),
        annotate_reactions: tuple = ('bin', 'kegg_reaction', 'ec_number'),
        annotate_metabolites: tuple = ('bin', 'kegg_compound'),
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Export the network to a metabolic model file in JSON format.

        Parameters
        ==========
        path : str
            output JSON file path

        annotate_genes : tuple, ('bin', )
            Annotate gene entries in the JSON file with additional data, selecting from the following:

            'bin' : bins in which the gene occurs

            'all_ko' : all KOs associated with the gene

            'ko' : KOs associated with the gene that yielded reactions in the network

            'ko_e_value' : scores of KO associations with the gene; if 'all_ko' is provided, then
                each value corresponds to a KO in 'all_ko', whereas if only 'ko' is provided, then each
                value corresponds to a KO in 'ko'

        annotate_reactions : tuple, ('bin', 'kegg_reaction', 'ec_number')
            Annotate reaction entries in the JSON file with additional data, selecting from the following:

            'bin' : bins in which the reaction occurs

            'kegg_reaction' : KO-associated KEGG reaction IDs yielding the ModelSEED reaction

            'ec_number' : KO-associated EC numbers yielding the ModelSEED reaction

            'ko' : KOs yielding the ModelSEED reaction

        annotate_metabolites : tuple, ('bin', 'kegg_compound')
            Annotate metabolite entries in the JSON file with additional data, selecting from the following:

            'bin' : bins in which the metabolite occurs

            'kegg_compound' : KEGG compound aliases of the ModelSEED compound

            'ko' : KOs yielding the ModelSEED compound

        run : terminal.Run, terminal.Run()

        progress : terminal.Progress, terminal.Progress()
        """
        pass

class PangenomicNetwork(ReactionNetwork):
    """A reaction network predicted from KEGG KO and ModelSEED annotations of pangenomic gene clusters."""
    def __init__(self) -> None:
        # map gene cluster ID to gene cluster object
        self.gene_clusters: Dict[str, GeneCluster] = {}
        self.bins: Dict[str, GeneClusterBin] = {}
        self.collection: BinCollection = None
        super().__init__()

    def export_json(
        self,
        path: str,
        annotate_genes: tuple = ('genome', 'bin', ),
        annotate_reactions: tuple = ('genome', 'bin', 'kegg_reaction', 'ec_number'),
        annotate_metabolites: tuple = ('genome', 'bin', 'kegg_compound'),
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        """
        Export the network to a metabolic model file in JSON format. *Gene entries in this file
        represent gene clusters.* Optionally, gene, reaction, and metabolite entries in this file
        are annotated with the names of genomes and names of gene cluster bins in which they occur.

        Parameters
        ==========
        Parameters
        ==========
        path : str
            output JSON file path

        annotate_genes : tuple, ('genome', 'bin', )
            Annotate gene (cluster) entries in the JSON file with additional data, selecting
            from the following:

            'genome' : genomes in which the genes of the cluster occur

            'bin' : bins in which the gene cluster occurs

            'all_ko' : all KOs associated with genes in the cluster, sorted in descending order of
                the number of genes in the cluster that were associated with each KO and then mean
                e-value of gene-KO assignments

            'ko' : KOs associated with the gene cluster that yielded reactions in the network,
                sorted in descending order of the number of genes in the cluster that were
                associated with each KO and then mean e-value of gene-KO assignments

            'ko_count' : number of genes in the cluster that were associated with each KO; if
                'all_ko' is provided, then each value corresponds to a KO in 'all_ko', whereas if
                only 'ko' is provided, then each value corresponds to a KO in 'ko'

            'e_value' : mean scores of KO associations with genes in the cluster; if 'all_ko' is
                provided, then each value corresponds to a KO in 'all_ko', whereas if only 'ko' is
                provided, then each value corresponds to a KO in 'ko'

        annotate_reactions : tuple, ('genome', 'bin', 'kegg_reaction', 'ec_number')
            Annotate reaction entries in the JSON file with additional data, selecting from the following:

            'genome' : genomes in which the reaction occurs

            'bin' : bins in which the reaction occurs

            'kegg_reaction' : KO-associated KEGG reaction IDs yielding the ModelSEED reaction

            'ec_number' : KO-associated EC numbers yielding the ModelSEED reaction

            'ko' : KOs yielding the ModelSEED reaction

        annotate_metabolites : tuple, ('genome', 'bin', 'kegg_compound')
            Annotate metabolite entries in the JSON file with additional data, selecting from the following:

            'genome' : genomes in which the metabolite occurs

            'bin' : bins in which the metabolite occurs

            'kegg_compound' : KEGG compound aliases of the ModelSEED compound

            'ko' : KOs yielding the ModelSEED compound

        run : terminal.Run, terminal.Run()

        progress : terminal.Progress, terminal.Progress()
        """
        pass

class KEGGDatabase:
    """
    The KEGG KO and REACTION databases set up by anvi'o.

    By default, the databases are loaded from a default directory of KEGG KO and REACTION files
    unless an alternative directory is provided.
    """
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/MISC/PROTEIN_DATA/kegg')

    def __init__(self, db_dir: str = None) -> None:
        """
        Load and set up KO and reaction tables derived from downloaded definition files. These
        tables facilitate the lookup of KO IDs and names, EC numbers, and KEGG REACTION IDs.

        Parameters
        ==========
        db_dir : str, None
            Directory of KEGG KO and REACTION files to use instead of the default.
        """
        self.ko_table: pd.DataFrame = None
        self.reaction_table: pd.DataFrame = None
        # The KEGG release is stored in the KO and REACTION info files downloaded at the time of the
        # other database files.
        self.release: str = None

        if db_dir:
            if not os.path.isdir(db_dir):
                raise ConfigError(f"The provided KEGG database directory, '{db_dir}', was not recognized as a directory.")
        else:
            db_dir = self.default_dir
        ko_data_path = os.path.join(db_dir, 'ko_data.tsv')
        if not os.path.isfile(ko_data_path):
            raise ConfigError(f"The KO data table, 'ko_data.tsv', was not found in the database directory, '{db_dir}'.")
        reaction_data_path = os.path.join(db_dir, 'reaction_data.tsv')
        if not os.path.isfile(reaction_data_path):
            raise ConfigError(f"The KEGG REACTION data table, 'reaction_data.tsv', was not found in the database directory, '{db_dir}'.")
        ko_info_path = os.path.join(db_dir, 'ko_info.txt')
        if not os.path.isfile(ko_info_path):
            raise ConfigError(f"The KO info file, 'ko_info.txt', was not found in the database directory, '{db_dir}'.")
        reaction_info_path = os.path.join(db_dir, 'reaction_info.txt')
        if not os.path.isfile(reaction_info_path):
            raise ConfigError(f"The KEGG REACTION info file, 'reaction_info.txt', was not found in the database directory, '{db_dir}'.")

        self.ko_table = pd.read_csv(ko_data_path, sep='\t', header=0, index_col=0, low_memory=False)
        self.reaction_table = pd.read_csv(reaction_data_path, sep='\t', header=0, index_col=0, low_memory=False)

        ko_release = None
        with open(ko_info_path) as f:
            for line in f:
                if line[:2] == 'ko':
                    ko_release = line[2:].strip()
                    break
        assert ko_release is not None
        reaction_release = None
        with open(reaction_info_path) as f:
            for line in f:
                if line[:2] == 'rn':
                    reaction_release = line[2:].strip()
                    break
        assert reaction_release is not None
        assert ko_release == reaction_release
        self.release = ko_release

class ModelSEEDDatabase:
    """
    The ModelSEED Biochemistry database set up by anvi'o.

    By default, the database is loaded from a default directory of ModelSEED files unless an
    alternative directory is provided.
    """
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/MISC/PROTEIN_DATA/modelseed')

    # Compounds are identified as cytosolic or extracellular in ModelSEED reactions.
    compartment_ids = {0: 'c', 1: 'e'}

    def __init__(self, db_dir: str = None) -> None:
        """
        Load and set up reaction and compound tables from the data directory.

        Parameters
        ==========
        db_dir : str, None
            Directory of ModelSEED files to use instead of the default.
        """
        # The KEGG and EC tables are rearrangements of the ModelSEED reactions table facilitating
        # lookup of reaction data by KEGG REACTION ID or EC number rather than ModelSEED reaction ID.
        self.kegg_reactions_table: pd.DataFrame = None
        self.ec_reactions_table: pd.DataFrame = None
        self.compounds_table: pd.DataFrame = None
        # The SHA is the git commit hash of the ModelSEED repository, stored in a file when
        # anvi'o downloaded and set up the ModelSEED files.
        self.sha: str = None

        if db_dir:
            if not os.path.isdir(db_dir):
                raise ConfigError(f"The provided ModelSEED database directory, '{db_dir}', was not recognized as a directory.")
        else:
            db_dir = self.default_dir
        reactions_path = os.path.join(db_dir, 'reactions.tsv')
        if not os.path.isfile(reactions_path):
            raise ConfigError(f"The ModelSEED reactions table, 'reactions.tsv', was not found in the database directory, '{db_dir}'.")
        compounds_path = os.path.join(db_dir, 'compounds.tsv')
        if not os.path.isfile(compounds_path):
            raise ConfigError(f"The ModelSEED compounds table, 'compounds.tsv', was not found in the database directory, '{db_dir}'.")
        sha_path = os.path.join(db_dir, 'sha.txt')
        if not os.path.isfile(sha_path):
            raise ConfigError(f"The file, 'sha.txt', which contains the git commit hash of the ModelSEED database, was not found in the database directory, '{db_dir}'.")

        reactions_table = pd.read_csv(reactions_path, sep='\t', header=0, low_memory=False)
        self.compounds_table = pd.read_csv(compounds_path, sep='\t', header=0, index_col='id', low_memory=False)
        with open(sha_path) as f:
            self.sha = f.readline().strip()

        # Reorganize the reactions table to facilitate lookup of reaction data by KEGG REACTION ID.
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

        # Reorganize the reactions table to facilitate lookup of reaction data by EC number.
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

class Constructor:
    """Construct metabolic reaction network objects."""
    def __init__(
        self,
        kegg_dir: str,
        modelseed_dir: str,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        self.kegg_dir = kegg_dir
        self.modelseed_dir = modelseed_dir
        self.run = run
        self.progress = progress

        self.kegg_db = KEGGDatabase(self.kegg_dir)
        self.modelseed_db = ModelSEEDDatabase(self.modelseed_dir)

    def import_network(self, json: str) -> ReactionNetwork:
        """Import a reaction network from a metabolic model JSON file."""
        pass

    def load_network(
        self,
        contigs_db: str = None,
        genomes_storage_db: str = None,
        pan_db: str = None
    ) -> ReactionNetwork:
        self.check_network(contigs_db=contigs_db, genomes_storage_db=genomes_storage_db, pan_db=pan_db)

    def check_network(
        self,
        contigs_db: str = None,
        genomes_storage_db: str = None,
        pan_db: str = None
    ) -> bool:
        """
        Check that the reaction network stored in a database is derived from the current gene KO
        annotations in the database.

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
        bool
            True if the reaction network is derived from the current gene KO annotations in the
            database, else False.
        """
        return


        self.progress.new("Building reaction network")
        self.progress.update("...")

        network = SingleGenomeNetwork()

        modelseed_kegg_reactions_table = self.modelseed_db.kegg_reactions_table
        modelseed_ec_reactions_table = self.modelseed_db.ec_reactions_table
        modelseed_compounds_table = self.modelseed_db.compounds_table

        # the following dictionaries help track objects already added to the network
        # KEGG REACTION ID -> ModelSEED reaction objects
        added_kegg_reaction_dict: Dict[str, List[ModelSEEDReaction]] = {}
        # EC number -> ModelSEED reaction objects
        added_ec_number_dict: Dict[str, List[ModelSEEDReaction]] = {}

        # record KOs that annotated genes but for some reason are not found in the KO database
        undefined_ko_ids = []

        gene_function_calls_dict: Dict = contigs_super.gene_function_calls_dict
        total_ko_matches = len(gene_function_calls_dict)
        num_ko_matches_parsed = -1
        # loop through gene-KO matches recorded in the contigs database
        for gcid, gene_dict in gene_function_calls_dict.items():
            num_ko_matches_parsed += 1
            self.progress.update(f"Gene-KO matches parsed: {num_ko_matches_parsed} / {total_ko_matches}")

            if gcid in network.genes:
                # an object representing the gene was already added to the network
                gene = network.genes[gcid]
            else:
                gene = Gene()
                gene.gcid = gcid
                # add the gene to the network
                network.genes[gcid] = gene

            ko_data = gene_dict['KOfam']
            gene.e_values.append(float(ko_data[2]))
            ko_id = ko_data[0]
            if ko_id in network.kos:
                # An object representing the KO was already added to the network. Therefore, objects
                # representing associated reactions and metabolites were already added as well.
                gene.kos.append(network.kos[ko_id])
                continue
            else:
                ko = KO()
                ko.id = ko_id
                ko.name = ko_data[1]
                # add the KO to the gene
                gene.kos.append(ko)
                # add the KO to the network
                network.kos[ko_id] = ko

            # find KEGG reactions and EC numbers associated with the KO
            try:
                ko_info = self.kegg_db.ko_table.loc[ko.id]
            except KeyError:
                undefined_ko_ids.append(ko_id)
                continue
            ko_kegg_reaction_info: str = ko_info.loc['reactions']
            if pd.isna(ko_kegg_reaction_info):
                # the KO is not associated with KEGG reactions
                kegg_reaction_ids = []
            else:
                kegg_reaction_ids = ko_kegg_reaction_info.split()
            ko_ec_number_info: str = ko_info.loc['ec_numbers']
            if pd.isna(ko_ec_number_info):
                # the KO is not associated with EC numbers
                ec_numbers = []
            else:
                ec_numbers = ko_ec_number_info.split()

            if not (kegg_reaction_ids or ec_numbers):
                # ModelSEED reaction objects cannot be defined for the KO in the absence of
                # either KEGG reactions or EC numbers associated with the KO
                continue

            # separate KEGG REACTION IDs that have previously been used to create ModelSEEDReaction
            # objects from those that have not
            unadded_kegg_reaction_ids = []
            for kegg_reaction_id in kegg_reaction_ids:
                if not kegg_reaction_id in added_kegg_reaction_dict:
                    unadded_kegg_reaction_ids.append(kegg_reaction_id)
                    continue
                # retrieve the one or more ModelSEED reactions aliased by the KEGG REACTION ID
                reactions = added_kegg_reaction_dict[kegg_reaction_id]
                for reaction in reactions:
                    modelseed_id = reaction.modelseed_id
                    ko.reactions[modelseed_id] = reaction
                    ko.kegg_reaction_aliases[modelseed_id] = tuple(
                        set(kegg_reaction_ids).intersection(set(reaction.kegg_id_aliases))
                    )
                    ko.ec_number_aliases[modelseed_id] = tuple(
                        set(ec_numbers).intersection(set(reaction.ec_number_aliases))
                    )

            # separate EC numbers that have previously been used to create ModelSEEDReaction objects
            # from from those that have not
            unadded_ec_numbers = []
            for ec_number in ec_numbers:
                if not ec_number in added_ec_number_dict:
                    unadded_ec_numbers.append(ec_number)
                    continue
                # retrieve the one or more ModelSEED reactions aliased by the EC number
                reactions = added_ec_number_dict[ec_number]
                for reaction in reactions:
                    modelseed_id = reaction.modelseed_id
                    if reaction.modelseed_id in reaction.ec_number_aliases:
                        # a KEGG REACTION ID associated with the KO was already linked to the
                        # ModelSEED ID, so it would be redundant to create a new reaction object
                        # given the EC number association
                        continue
                    # an EC number but no KEGG REACTION ID associated with the KO was linked to
                    # the ModelSEED ID
                    ko.reactions[modelseed_id] = reaction
                    ko.kegg_reaction_aliases[modelseed_id] = tuple(
                        set(kegg_reaction_ids).intersection(set(reaction.kegg_id_aliases))
                    )
                    assert len(ko.kegg_reaction_aliases) == 0
                    ko.ec_number_aliases[modelseed_id] = tuple(
                        set(ec_numbers).intersection(set(reaction.ec_number_aliases))
                    )

            if not (unadded_kegg_reaction_ids or unadded_ec_numbers):
                # all of the KEGG reactions and EC numbers associated with the KO have already been
                # encountered in previously processed KOs, so proceed to the next gene KO annotation
                continue

            # get data on unencountered ModelSEED reactions aliased by KEGG REACTION IDs and EC
            # numbers
            modelseed_reactions_data = {}
            if unadded_kegg_reaction_ids:
                modelseed_kegg_reactions_dict: Dict[str, Dict] = modelseed_kegg_reactions_table[
                    modelseed_kegg_reactions_table['KEGG_REACTION_ID'].isin(unadded_kegg_reaction_ids)
                ].to_dict(orient='index')
                for modelseed_reaction_data in modelseed_kegg_reactions_dict.values():
                    # data on the ModelSEED reaction has already been added to the list using a
                    # different alias
                    modelseed_reaction_id = modelseed_reaction_data['id']
                    if modelseed_reaction_id not in modelseed_reactions_data:
                        modelseed_reactions_data[modelseed_reaction_id] = modelseed_reaction_data
            if unadded_ec_numbers:
                modelseed_ec_reactions_dict: Dict[str, Dict] = modelseed_ec_reactions_table[
                    modelseed_ec_reactions_table['EC_number'].isin(unadded_ec_numbers)
                ].to_dict(orient='index')
                for modelseed_reaction_data in modelseed_ec_reactions_dict.values():
                    # data on the ModelSEED reaction has already been added to the list using a
                    # different alias
                    modelseed_reaction_id = modelseed_reaction_data['id']
                    if modelseed_reaction_id not in modelseed_reaction_data:
                        modelseed_reactions_data[modelseed_reaction_id] = modelseed_reaction_data
            if not modelseed_reactions_data:
                # the unadded KEGG REACTION IDs and EC numbers are not in the ModelSEED reactions table
                continue

            # generate new reaction objects in the network
            for modelseed_reaction_id, modelseed_reaction_data in modelseed_reactions_data.items():
                reaction, modelseed_compound_ids = self._get_modelseed_reaction(modelseed_reaction_data)
                if reaction is None:
                    # the reaction does not have an equation in ModelSEED for some reason
                    continue
                ko.reactions[modelseed_reaction_id] = reaction
                ko.kegg_reaction_aliases[modelseed_reaction_id] = tuple(
                    set(unadded_kegg_reaction_ids).intersection(set(reaction.kegg_id_aliases))
                )
                ko.ec_number_aliases[modelseed_reaction_id] = tuple(
                    set(unadded_ec_numbers).intersection(set(reaction.ec_number_aliases))
                )
                network.reactions[modelseed_reaction_id] = reaction

                # separate ModelSEED compund IDs that have previously been used to create
                # ModelSEEDCompound objects from those that have not
                unadded_modelseed_compound_ids = []
                reaction_compounds = []
                for modelseed_compound_id in modelseed_compound_ids:
                    if modelseed_compound_id in network.metabolites:
                        # the ModelSEED compound ID has been encountered in previously processed reactions
                        reaction_compounds.append(network.metabolites[modelseed_compound_id])
                    else:
                        unadded_modelseed_compound_ids.append(modelseed_compound_id)

                # generate new metabolite objects in the network
                for modelseed_compound_id in unadded_modelseed_compound_ids:
                    try:
                        modelseed_compound_series: pd.Series = modelseed_compounds_table.loc[modelseed_compound_id]
                    except KeyError:
                        raise ConfigError(
                            f"A row for the ModelSEED compound ID, '{modelseed_compound_id}', was expected "
                            "but not found in the ModelSEED compounds table. This ID was found in the equation "
                            f"for ModelSEED reaction, '{modelseed_reaction_id}'."
                        )
                    modelseed_compound_data = modelseed_compound_series.to_dict()
                    modelseed_compound_data['id'] = modelseed_compound_id
                    compound = self._get_modelseed_compound(modelseed_compound_data)
                    reaction_compounds.append(compound)
                    network.metabolites[modelseed_compound_id] = compound
                reaction.compounds = tuple(reaction_compounds)

        if undefined_ko_ids:
            self.run.info_single(
                "Certain genes matched KOs that were not found in the KO database. "
                "It could be that the KOfams used to annotate genes were not from the same KEGG "
                "database version as the KO definition files. Here are the unrecognized KO IDs "
                f"matching genes in the contigs database: {', '.join(undefined_ko_ids)}"
            )

        self.progress.update("Writing network to contigs database")
        self._store_reactions_in_contigs_db(network, contigs_db_path)
        self._store_metabolites_in_contigs_db(network, contigs_db_path)

        contigs_db = ContigsDatabase(contigs_db_path)
        gene_call_count = contigs_db.db.get_row_counts_from_table('genes_in_contigs')
        contigs_db.disconnect()

        self.progress.update("Calculating network statistics")

        contributing_gcids = []
        contributing_ko_ids = []
        for gcid, gene in network.genes.items():
            for ko in gene.kos:
                if ko.reactions:
                    contributing_gcids.append(gcid)
                    contributing_ko_ids.append(ko.id)
                    break
        contributing_gcids = set(contributing_gcids)
        contributing_ko_ids = set(contributing_ko_ids)

        kegg_plus_ec_alias_count = 0
        only_kegg_alias_count = 0
        only_ec_alias_count = 0
        contributing_kegg_reaction_ids = []
        contributing_ec_numbers = []
        for reaction in network.reactions.values():
            if reaction.kegg_id_aliases and reaction.ec_number_aliases:
                kegg_plus_ec_alias_count += 1
            elif reaction.kegg_id_aliases:
                only_kegg_alias_count += 1
            elif reaction.ec_number_aliases:
                only_ec_alias_count += 1
            else:
                raise ConfigError("It should not be possible for a reaction to have neither a KEGG REACTION ID nor EC number alias!")
            contributing_kegg_reaction_ids += list(reaction.kegg_id_aliases)
            contributing_ec_numbers += list(reaction.ec_number_aliases)
        contributing_kegg_reaction_ids = set(contributing_kegg_reaction_ids)
        contributing_ec_numbers = set(contributing_ec_numbers)

        self.progress.end()

        self.run.info("Total gene calls", gene_call_count)
        self.run.info("Genes annotated with KOs", len(network.genes))
        self.run.info("Genes contributing to network", len(contributing_gcids))
        self.run.info("KOs annotating genes", len(network.kos))
        self.run.info("KOs contributing to network", len(contributing_ko_ids))
        self.run.info("Reactions (ModelSEED) in network", len(network.reactions))
        self.run.info("Reactions aliased by KEGG rxn & EC number", kegg_plus_ec_alias_count)
        self.run.info("Reactions aliased only by KEGG reaction", only_kegg_alias_count)
        self.run.info("Reactions aliased only by EC number", only_ec_alias_count)
        self.run.info("KEGG reactions contributing to network", len(contributing_kegg_reaction_ids))
        self.run.info("EC numbers contributing to network", len(contributing_ec_numbers))
        self.run.info("Metabolites in network", len(network.metabolites))

    def _store_reactions_in_contigs_db(self, network: SingleGenomeNetwork, contigs_db_path: str) -> None:
        """Store reaction data in the contigs database table."""
        reactions_data = {}
        # transfer data from reaction objects
        for modelseed_reaction_id, reaction in network.reactions.items():
            reaction_data = {}
            reaction_data['modelseed_reaction_id'] = modelseed_reaction_id
            reaction_data['modelseed_reaction_name'] = reaction.modelseed_name
            reaction_data['metabolite_modelseed_ids'] = ', '.join([c.modelseed_id for c in reaction.compounds])
            reaction_data['stoichiometry'] = ', '.join([str(c) for c in reaction.coefficients])
            reaction_data['compartments'] = ', '.join(reaction.compartments)
            reaction_data['reversibility'] = reaction.reversibility
            reactions_data[modelseed_reaction_id] = reaction_data

        # Get *KO* KEGG REACTION ID and EC number aliases of each ModelSEED reaction. These are not
        # all possible aliases, but only those associated with KOs that matched genes. Structure
        # alias data as follows:
        # <ModelSEED reaction ID>: {
        #   <KEGG REACTION ID 1>: [<KO IDs associated with KEGG REACTION ID 1>],
        #   <KEGG REACTION ID 2>: [<KO IDs associated with  KEGG REACTION ID 2>], ...}
        # <ModelSEED reaction ID>: {
        #   <EC number 1>: [<KO IDs associated with EC number 1>],
        #   <EC number 2>: [<KO IDs associated with EC number 2>], ...}
        ko_reaction_aliases = {modelseed_reaction_id: ({}, {}) for modelseed_reaction_id in reactions_data}
        for ko_id, ko in network.kos.items():
            for modelseed_reaction_id, reaction in ko.reactions.items():
                aliases = ko_reaction_aliases[modelseed_reaction_id]

                kegg_reaction_aliases = aliases[0]
                kegg_reaction_ids = ko.kegg_reaction_aliases[modelseed_reaction_id]
                for kegg_reaction_id in kegg_reaction_ids:
                    try:
                        ko_ids: list = kegg_reaction_aliases[kegg_reaction_id]
                    except KeyError:
                        kegg_reaction_aliases[kegg_reaction_id] = ko_ids = []
                    ko_ids.append(ko_id)

                ec_number_aliases = aliases[1]
                ec_numbers = ko.ec_number_aliases[modelseed_reaction_id]
                for ec_number in ec_numbers:
                    try:
                        ko_ids: list = ec_number_aliases[ec_number]
                    except KeyError:
                        ec_number_aliases[ec_number] = ko_ids = []
                    ko_ids.append(ko_id)
        for modelseed_reaction_id, aliases in ko_reaction_aliases.items():
            reaction_data = reactions_data[modelseed_reaction_id]

            kegg_reaction_aliases = aliases[0]
            entry = []
            for kegg_reaction_id, ko_ids in kegg_reaction_aliases.items():
                entry.append(f'{kegg_reaction_id}: ({", ".join(sorted(ko_ids))})')
            reaction_data['ko_kegg_reaction_source'] = '; '.join(sorted(entry))

            ec_number_aliases = aliases[1]
            entry = []
            for ec_number, ko_ids in ec_number_aliases.items():
                entry.append(f'{ec_number}: ({", ".join(sorted(ko_ids))})')
            reaction_data['ko_ec_number_source'] = '; '.join(sorted(entry))

        reactions_table = pd.DataFrame.from_dict(reactions_data, orient='index').reset_index(drop=True).sort_values('modelseed_reaction_id')
        reactions_table = reactions_table[tables.gene_function_reactions_table_structure]

        contigs_db = ContigsDatabase(contigs_db_path)
        contigs_db.db._exec_many(
            f'''INSERT INTO {tables.gene_function_reactions_table_name} VALUES ({','.join('?' * len(tables.gene_function_reactions_table_structure))})''',
            reactions_table.values
        )
        contigs_db.disconnect()

    def _store_metabolites_in_contigs_db(self, network: SingleGenomeNetwork, contigs_db_path: str) -> None:
        """Store metabolite data in the contigs database table."""
        metabolites_data = {}
        for modelseed_compound_id, compound in network.metabolites.items():
            metabolite_data = {}
            metabolite_data['modelseed_compound_id'] = modelseed_compound_id
            metabolite_data['modelseed_compound_name'] = compound.modelseed_name
            metabolite_data['formula'] = compound.formula
            metabolite_data['charge'] = compound.charge
            metabolites_data[modelseed_compound_id] = metabolite_data
        metabolites_table = pd.DataFrame.from_dict(metabolites_data, orient='index').reset_index(drop=True).sort_values('modelseed_compound_id')
        metabolites_table = metabolites_table[tables.gene_function_metabolites_table_structure]

        contigs_db = ContigsDatabase(contigs_db_path)
        contigs_db.db._exec_many(
            f'''INSERT INTO {tables.gene_function_metabolites_table_name} VALUES ({','.join('?' * len(tables.gene_function_metabolites_table_structure))})''',
            metabolites_table.values
        )
        contigs_db.disconnect()

    def _load_contigs_db(self, contigs_db_path: str) -> ContigsSuperclass:
        is_contigs_db(contigs_db_path)
        args = Namespace()
        args.contigs_db = contigs_db_path
        contigs_super = ContigsSuperclass(args, r=run_quiet)
        contigs_super.init_functions(requested_sources=['KOfam'])
        return contigs_super

    def _get_modelseed_reaction(self, modelseed_reaction_data: Dict) -> Tuple[ModelSEEDReaction, List[str]]:
        """
        Create a ModelSEED reaction object from its entry in the ModelSEED table.

        Do not populate the reaction object with metabolite objects. Return both the new reaction
        object and a list of associated ModelSEED compound IDs.
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
            reaction.kegg_id_aliases = []
        else:
            reaction.kegg_id_aliases = kegg_reaction_ids.split('; ')

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
            compartments.append(self.compartment_ids[int(split_entry[2])])
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

    def _to_lcm_denominator(self, floats) -> Tuple[int]:
        def lcm(a, b):
            return a * b // gcd(a, b)
        rationals = [Fraction(f).limit_denominator() for f in floats]
        lcm_denom = reduce(lcm, [r.denominator for r in rationals])
        return tuple(int(r.numerator * lcm_denom / r.denominator) for r in rationals)

    def _get_modelseed_compound(self, modelseed_compound_data: Dict) -> ModelSEEDCompound:
        compound = ModelSEEDCompound()
        compound.modelseed_id = modelseed_compound_data['id']

        modelseed_name = modelseed_compound_data['name']
        if pd.isna(modelseed_name):
            compound.modelseed_name = None
        else:
            compound.modelseed_name = modelseed_name

        kegg_id_aliases: str = modelseed_compound_data['KEGG']
        if pd.isna(kegg_id_aliases):
            compound.kegg_id_aliases = []
        else:
            compound.kegg_id_aliases = kegg_id_aliases.split('; ')

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
