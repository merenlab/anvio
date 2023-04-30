# -*- coding: utf-8
# pylint: disable=line-too-long
"""Representation of gene protein homologs."""

from __future__ import annotations

import pandas as pd

from typing import Dict, List, Tuple

import anvio.terminal as terminal
import anvio.proteinorthology.refdbs as refdbs

from anvio.errors import ConfigError
from anvio import __version__ as VERSION


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


run_quiet = terminal.Run(verbose=False)

class Chemical:
    """A chemical."""
    def __init__(self) -> None:
        # The chemical may have one or more IDs (values) per reference database (key).
        self.reference_ids: Dict[str, List[str]] = {}
        self.charge: int = None
        self.formula: str = None
        self.smiles_string: str = None

    @property
    def id(self) -> str:
        """The first recorded reference database/ID. If none is recorded, return None."""
        try:
            return next(iter(self.reference_ids.values()))[0]
        except StopIteration:
            return None

class Reaction:
    """A chemical reaction."""
    def __init__(self) -> None:
        # The reaction may have one or more IDs (values) per reference database (key).
        self.reference_ids: Dict[str, List[str]] = {}
        self.chemicals: List[Chemical] = []
        self.coefficients: List[float] = []
        self.compartments: List[str] = []
        self.reversibility: bool = None

    @property
    def id(self) -> str:
        """The first recorded reference database/ID. If none is recorded, return None."""
        try:
            return next(iter(self.reference_ids.values()))[0]
        except StopIteration:
            return None

class OrthologAnnotation:
    """Data regarding a group of protein functional orthologs."""
    def __init__(self) -> None:
        self.source: str = None # e.g., 'KEGG' or 'Pfam' databases
        self.accession: str = None # e.g., KEGG ortholog 'K00001'

class AnvioOrthologAnnotation(OrthologAnnotation):
    """
    An ortholog annotation of a gene/gene cluster as stored in an anvi'o database.

    Parameters
    ==========
    entry : Dict
        This dictionary contains information on an ortholog annotation. For a gene, this
        corresponds to a row of the gene functions table of a contigs database.
    """
    def __init__(self, entry: Dict) -> None:
        if 'gene_callers_id' in entry:
            self.id = int(entry['gene_callers_id'])
        else:
            self.id = entry['gene_cluster_id']
        source = entry['source']
        if source == 'KEGG' or source == 'KOfam':
            # In anvi'o databases, the source of KEGG orthologs is recorded as 'KOfam'. This
            # suggests that KOs represented here must derived from a KOfam HMM, which need not be
            # the case.
            self.source = 'KEGG'
        else:
            self.source = source # e.g., 'Pfam'
        self.accession = entry['accession'] # e.g., KEGG ortholog 'K00001'
        self.function = entry['function'] # description of gene function
        self.e_value = float(entry['evalue']) # annotation confidence score

    def __str__(self) -> str:
        if isinstance(self.id, int):
            return f"Gene '{self.id}' {self.source} ortholog accession, '{self.accession}'"
        else:
            return f"Gene cluster '{self.id}' {self.source} ortholog accession, '{self.accession}'"

class AnvioKOAnnotation(AnvioOrthologAnnotation):
    """A KEGG ortholog (KO) of a gene as stored in an anvi'o contigs database."""
    def __init__(self, entry: Dict) -> None:
        super().__init__(entry)

    def get_reactions(
        self,
        kegg_db: refdbs.KEGGDatabase,
        cross_reference_dbs: Tuple[refdbs.ProteinReferenceDatabase] = None
    ) -> List[Reaction]:
        """
        Get reaction data for the ortholog.

        Parameters
        ==========
        kegg_db : anvio.proteinorthology.refdbs.KEGGDatabase
            KEGG reference database.
        cross_reference_dbs : tuple
            Protein reference databases ('anvio.proteinorthology.refdbs.ProteinReferenceDatabase')
            with which KOs are cross-referenced. For now, a ModelSEED database must be supplied as
            the sole cross-reference database.

        Returns
        =======
        List[Reaction]
        """
        for db in (kegg_db, ) + cross_reference_dbs:
            db._check_reference_database_initialization()
        cross_ref_db_names = tuple(db.db_name for db in cross_reference_dbs)
        if cross_ref_db_names == ('modelseed', ):
            reactions = self._get_reactions_from_kegg_and_modelseed(
                kegg_db, cross_reference_dbs[0]
            )
        else:
            raise ConfigError(
                "For now, a ModelSEED database must be supplied as the sole cross-reference "
                "database."
            )
        return reactions

    def _get_reactions_from_kegg_and_modelseed(
        self,
        kegg_db: refdbs.KEGGDatabase,
        modelseed_db: refdbs.ModelSEEDDatabase
    ) -> List[Reaction]:
        """
        Get reaction data for the ortholog from KEGG and ModelSEED reference databases.

        A KO may be associated with KEGG Reactions and EC numbers. This method first considers KEGG
        Reactions and, absent these, EC numbers. Reaction IDs and EC numbers are cross-referenced to
        the ModelSEED database. Reaction data is not returned for Reaction IDs and EC numbers that
        are not found in the ModelSEED database!

        Parameters
        ==========
        kegg_db : anvio.proteinorthology.refdbs.KEGGDatabase
            KEGG reference database.
        modelseed_db : anvio.proteinorthology.refdbs.ModelSEEDDatabase
            ModelSEED reference database. The ModelSEED Biochemistry Database has harmonized
            reaction and compound data with KEGG and other databases.

        Returns
        =======
        List[Reaction]
        """
        # Retrieve or generate ModelSEED reactions tables for looking up reactions by KEGG Reaction
        # ID and EC number, respectively.
        if hasattr(modelseed_db, 'reaction_lookup_tables'):
            lookup_dict_existed = True
            if (
                'KEGG' in modelseed_db.reaction_lookup_tables and
                'ec_numbers' in modelseed_db.reaction_lookup_tables
            ):
                lookup_tables_existed = True
            else:
                lookup_tables_existed = False
        else:
            lookup_dict_existed = False
            lookup_tables_existed = False
        if not lookup_tables_existed:
            modelseed_db._check_reference_database_initialization()
            modelseed_db._set_reaction_lookup_table('KEGG')
            modelseed_db._set_reaction_lookup_table('ec_numbers')
        kegg_reactions_table = modelseed_db.reaction_lookup_tables['KEGG']
        ec_numbers_reactions_table = modelseed_db.reaction_lookup_tables['ec_numbers']
        # Restore the passed ModelSEED database object to its original state by removing KEGG or EC
        # lookup tables added by the current method.
        if not lookup_dict_existed:
            delattr(modelseed_db, 'reaction_lookup_tables')
        elif not lookup_tables_existed:
            modelseed_db.reaction_lookup_tables.pop('KEGG')
            modelseed_db.reaction_lookup_tables.pop('ec_numbers')
        # Use any KEGG Reaction IDs associated with the KO to find cross-referenced ModelSEED
        # reactions.
        kegg_series = kegg_db.ko_data.loc[self.accession]
        reactions = []
        reaction_accessions = kegg_series.loc['reactions']
        if pd.notna(reaction_accessions):
            reaction_accessions: str
            for reaction_data in kegg_reactions_table[
                kegg_reactions_table['KEGG'].isin(reaction_accessions.split())
            ].to_dict(orient='index').values():
                reaction = modelseed_db.get_reaction(reaction_data)
                if reaction is None:
                    continue
                reactions.append(reaction)
            if reactions:
                return reactions
        # Reaching this point, either no KEGG Reaction IDs were associated with the KO, or Reaction
        # IDs for the KO were not cross-referenced with any ModelSEED reactions. Next use any EC
        # numbers associated with the KO to find cross-referenced ModelSEED reactions. KO Reactions
        # are preferred as a more precise representation of the reactions that can be catalyzed than
        # EC numbers.
        ec_numbers = kegg_series.loc['ec_numbers']
        if pd.notna(ec_numbers):
            ec_numbers: str
            for reaction_data in ec_numbers_reactions_table[
                ec_numbers_reactions_table['ec_numbers'].isin(ec_numbers.split())
            ].to_dict(orient='index').values():
                reaction = modelseed_db.get_reaction(reaction_data)
                if reaction is None:
                    continue
                reactions.append(reaction)
            if reactions:
                return reactions
        # No reaction data could be recovered for the KO.
        return reactions
