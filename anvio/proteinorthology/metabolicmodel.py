# -*- coding: utf-8
# pylint: disable=line-too-long
"""Metabolic model generation tools."""

from __future__ import annotations

import os
import json
import pandas as pd

from argparse import Namespace
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Set, Tuple

import anvio.terminal as terminal
import anvio.proteinorthology.refdbs as refdbs
import anvio.proteinorthology.protein as protein

from anvio.errors import ConfigError
from anvio.terminal import Run, Progress
from anvio.ccollections import Collections
from anvio.filesnpaths import is_output_file_writable
from anvio.dbops import ContigsSuperclass, PanSuperclass
from anvio import TABULATE, QUIET, __version__ as VERSION
from anvio.utils import is_contigs_db, is_genome_storage, is_pan_db


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


run_quiet = terminal.Run(verbose=False)

class COBRApyJSONStructure:
    """COBRApy JSON input file structure."""

    @staticmethod
    def get() -> Dict[str, Any]:
        """JSON format."""
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

    @staticmethod
    def get_metabolite_entry() -> Dict[str, Any]:
        """Format of each object in the JSON 'metabolites' array."""
        return {
            'id': '',
            'name': '',
            'compartment': '',
            'charge': 0, # placeholder: uncharged
            'formula': '',
            'notes': {},
            'annotation': {}
        }

    @staticmethod
    def get_reaction_entry() -> Dict[str, Any]:
        """Format of each reaction object in the JSON 'reactions' array."""
        return {
            'id': '',
            'name': '',
            'metabolites': {},
            'lower_bound': -1000.0, # placeholder: reversible reaction
            'upper_bound': 1000.0,
            'gene_reaction_rule': '',
            'subsystem': '',
            'notes': {},
            'annotation': {}
        }

    @staticmethod
    def get_gene_entry() -> Dict[str, Any]:
        """Format of each object in the JSON 'genes' array."""
        return {
            'id': '',
            'name': '',
            'notes': {},
            'annotation': {}
        }

    @staticmethod
    def get_ecoli_objective() -> Dict[str, Any]:
        """Biomass objective from JSON 'reactions' array in the COBRApy example file,
        'e_coli_core.json'. BiGG metabolite IDs have been replaced with KBase/ModelSEED compound
        IDs."""
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
                'original_metabolite_names': {
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

class ModelInput(ABC):
    """Metabolic model input data."""
    def set_protein_annotation_source(
        self,
        source: str = 'KEGG',
        cross_references: Tuple[str] = ('ModelSEED', ),
        db_superdir: str = refdbs.ProteinReferenceDatabase.default_superdir
    ) -> None:
        """
        Set the source of protein annotations used in the model and prepare for analysis of
        orthologs.

        Parameters
        ==========
        source : str
            The source of protein ortholog annotations. For now, only KEGG orthologs can be
            processed. Default, 'KEGG'.
        cross_references : tuple
            Protein reference databases with which orthologs are cross-referenced. For now,
            'ModelSEED' must be supplied as the sole cross-reference with source 'KEGG'. Default,
            ('ModelSEED', )
        db_superdir : str
            The directory containing protein reference database subdirectories, e.g., 'kegg',
            'modelseed'. The source and cross-references correspond to these databases.
        """
        if source != 'KEGG':
            raise ConfigError(
                "For now, the only protein annotation source that can be processed is 'KEGG', not "
                f"'{source}'."
            )
        if source == 'KEGG':
            if cross_references != ('ModelSEED', ):
                raise ConfigError(
                    "For now, the 'ModelSEED' database is the sole required cross-reference for "
                    "source, 'KEGG', and not "
                    f"'{', '.join(r for r in cross_references if r != 'ModelSEED')}'."
                )
        self.protein_annotation_source = source
        self.cross_references = cross_references
        for db in (self.protein_annotation_source, ) + self.cross_references:
            db: str
            if not os.path.isdir(os.path.join(db_superdir, db.lower())):
                raise ConfigError(
                    f"'{db}' database files should be located in a directory named '{db.lower()}' "
                )
        if self.protein_annotation_source == 'KEGG':
            self.source_db = refdbs.KEGGDatabase(db_superdir=db_superdir)
        self.source_db.load()
        cross_ref_dbs = []
        for r in self.cross_references:
            if r == 'ModelSEED':
                cross_ref_db = refdbs.ModelSEEDDatabase(db_superdir=db_superdir)
                cross_ref_db.load()
                cross_ref_dbs.append(cross_ref_db)
        self.cross_ref_dbs = tuple(sorted(cross_ref_dbs, key=lambda db: type(db).__name__))

    def write_cobrapy_json(self, path: str, indent: int = 2) -> None:
        """
        Write a COBRApy JSON file from the input data.

        Parameters
        ==========
        path : str
            Output file path.
        indent : int
            Number of spaces of indentation per nesting level in JSON file. Default, 2.
        """
        is_output_file_writable(path)
        if self.protein_annotation_source == 'KEGG' and self.cross_references == ('ModelSEED', ):
            # When obtaining ortholog data from KEGG cross-referenced with ModelSEED, setting
            # reaction lookup tables indexed by KEGG ID and EC number greatly speeds up the process.
            modelseed_db: refdbs.ModelSEEDDatabase = self.cross_ref_dbs[0]
            modelseed_db._set_reaction_lookup_table('KEGG')
            modelseed_db._set_reaction_lookup_table('ec_numbers')
        cobrapy_dict = self._get_cobrapy_json_dict()
        self._find_missing_objective_metabolites(
            cobrapy_dict, remove=self.remove_missing_objective_metabolites
        )
        with open(path, 'w') as f:
            json.dump(cobrapy_dict, f, indent=indent)
        run: Run = self.run
        run.info("Metabolic model file", path, nl_before=1)

    @abstractmethod
    def _get_cobrapy_json_dict(self) -> Dict[str, Any]:
        """Get a dictionary representation of a COBRApy JSON file from the input data."""
        raise NotImplementedError

    @abstractmethod
    def _get_cobrapy_gene_dict(self) -> Dict[str, Any]:
        """Get a dictionary representation of a gene object in a COBRApy JSON file genes array."""
        raise NotImplementedError

    @abstractmethod
    def _get_cobrapy_reaction_dict(self) -> Tuple[Dict[str, Any], bool]:
        """
        Get a dictionary representation of a reaction object in a COBRApy JSON file reactions array.

        Returns
        =======
        Dict[str, Any]
            Reaction entry
        bool
            True if the reaction had already been recorded in the COBRApy JSON dictionary
        """
        raise NotImplementedError

    @abstractmethod
    def _get_cobrapy_metabolite_dict(self) -> Tuple[Dict[str, Any], bool]:
        """
        Get a dictionary representation of a metabolite object in a COBRApy JSON file metabolites
        array.

        Returns
        =======
        Dict[str, Any]
            Metabolite entry
        bool
            True if the metabolite had already been recorded in the COBRApy JSON dictionary
        """
        raise NotImplementedError

    def _find_missing_objective_metabolites(
        self,
        cobrapy_dict: Dict[str, Any],
        remove: bool = False
    ) -> None:
        """
        Find metabolites in the biomass objective function that are not produced or consumed by
        enzymes in the model.

        Parameters
        ==========
        cobrapy_dict : Dict
            Dictionary representation of a COBRApy JSON input file
        remove : bool
            If True (default False), remove metabolites from the biomass objective function that are
            not produced or consumed by enzymes in the model.
        """
        objective: Dict[str, Any] = cobrapy_dict['reactions'][0]
        try:
            objective_metabolites: Dict[str, float] = objective['metabolites']
        except KeyError:
            raise ConfigError(
                "The objective function of the metabolic model does not have a 'metabolites' "
                "section characteristic of an expected biomass objective function."
            )
        reaction_metabolites: List[Dict[str, Any]] = cobrapy_dict['metabolites']
        reaction_metabolite_ids = [m['id'] for m in reaction_metabolites]
        missing_objective_metabolite_data = []
        for metabolite_id, stoichiometric_coefficient in objective_metabolites.items():
            if metabolite_id in reaction_metabolite_ids:
                continue
            missing_objective_metabolite_data.append((metabolite_id, stoichiometric_coefficient))
        if not missing_objective_metabolite_data:
            return
        if remove:
            for metabolite_id, stoichiometric_coefficient in missing_objective_metabolite_data:
                objective_metabolites.pop(metabolite_id)
        run: Run = self.run
        if remove:
            run.warning(
                f"Missing metabolite data removed from '{objective['name']}'", nl_after=0
            )
        else:
            run.warning(
                f"Metabolites in '{objective['name']}' missing from the model", nl_after=0
            )
        if QUIET:
            return objective
        # Print removed metabolites.
        TABULATE(
            pd.DataFrame(
                missing_objective_metabolite_data,
                index=range(1, len(missing_objective_metabolite_data) + 1)
            ),
            ('Metabolite', 'Stoichiometric coefficient')
        )
        return objective

class Pangenome(ModelInput):
    """
    Metabolic model input data from a pangenome.

    Parameters
    ==========
    genomes_storage_path : str
        Anvi'o database of each genome in the pangenome.
    pan_db_path : str
        Anvi'o database of pangenome gene clusters.
    collection_name : str
        With bin ID, use to select gene clusters in the pangenome. Default, None.
    bin_id : str
        With collection name, use to select gene clusters in the pangenome. Default, None.
    protein_annotation_source : str
        The source of protein ortholog annotations. For now, only KEGG orthologs can be processed.
        Default, 'KEGG'.
    cross_references : tuple
        Protein reference databases with which orthologs are cross-referenced. For now, 'ModelSEED'
        must be supplied as the sole cross-reference with source 'KEGG'. Default, ('ModelSEED', )
    db_superdir : str
        The directory containing protein reference database subdirectories, e.g., 'kegg',
        'modelseed'. The source and cross-references correspond to these databases.
    discard_ties : bool
        If multiple protein annotations are most frequent among genes in a cluster, then do not
        assign an annotation to the cluster itself when this argument is True. By default, this
        argument is False, so one of the most frequent annotations would be arbitrarily chosen.
    consensus_threshold : float
        Without this argument (default None), the protein annotation most frequent among genes
        in a cluster is assigned to the cluster itself. With this argument (a value on [0, 1]),
        at least this proportion of genes in the cluster must have the most frequent annotation
        for the cluster to be annotated.
    remove_missing_objective_metabolites : bool
        Metabolites in the biomass objective function are removed if they are not recorded as being
        produced or consumed by enzymes in the model.
    """
    def __init__(
        self,
        genomes_storage_path: str,
        pan_db_path: str,
        collection_name: str = None,
        bin_id: str = None,
        protein_annotation_source: str = 'KEGG',
        cross_references: Tuple[str] = ('ModelSEED', ),
        db_superdir: str = refdbs.ProteinReferenceDatabase.default_superdir,
        discard_ties: bool = False,
        consensus_threshold: float = None,
        remove_missing_objective_metabolites: bool = False,
        run: Run = Run(),
        progress: Progress = Progress()
    ) -> None:
        self.genomes_storage_path = genomes_storage_path
        is_genome_storage(self.genomes_storage_path)
        self.pan_db_path = pan_db_path
        is_pan_db(self.pan_db_path)
        self.collection_name = collection_name
        self.bin_id = bin_id
        self._init_bin()
        self.set_protein_annotation_source(
            source=protein_annotation_source,
            cross_references=cross_references,
            db_superdir=db_superdir
        )
        self.discard_ties = discard_ties
        self.consensus_threshold = consensus_threshold
        self.remove_missing_objective_metabolites = remove_missing_objective_metabolites
        self.run = run
        self.progress = progress
        self._init_pan_super()

    def _init_bin(self) -> None:
        """Select a bin of gene clusters for consideration."""
        if self.collection_name or self.bin_id:
            self.collections = Collections(r=run_quiet)
            self.collections.populate_collections_dict(self.pan_db_path)
            self.collections.is_bin_in_collection(self.collection_name, self.bin_id)
            self.select_gene_cluster_ids = set(
                self.collections.get_collection_dict(self.collection_name)[self.bin_id]
            )
        else:
            self.collections = None
            self.select_gene_cluster_ids = set()

    def _init_pan_super(self) -> None:
        """Set up the pangenomic data."""
        args = Namespace()
        args.pan_db = self.pan_db_path
        args.genomes_storage = self.genomes_storage_path
        args.discard_ties = self.discard_ties
        args.consensus_threshold = self.consensus_threshold
        self.pan_super = PanSuperclass(args, r=run_quiet)
        self.pan_super.init_gene_clusters(gene_cluster_ids_to_focus=self.select_gene_cluster_ids)
        self.pan_super.init_gene_clusters_functions()
        self.pan_super.init_gene_clusters_functions_summary_dict()

    def _get_cobrapy_json_dict(self) -> Dict[str, Any]:
        cobrapy_dict = COBRApyJSONStructure.get()
        self.cobrapy_reactions: List[Dict] = cobrapy_dict['reactions']
        self.cobrapy_reactions.append(COBRApyJSONStructure.get_ecoli_objective())
        self.cobrapy_metabolites: List[Dict] = cobrapy_dict['metabolites']
        self.recorded_reactions: Dict[str, Dict] = {}
        self.recorded_metabolites: Dict[str, Dict] = {}
        self.progress.new("Analyzing gene clusters")
        num_analyzed = 0
        total = len(self.pan_super.gene_clusters)
        for gene_cluster_id, genome_gcids in self.pan_super.gene_clusters.items():
            self.progress.update(f"{num_analyzed} / {total}")
            num_analyzed += 1
            genome_gcids: Dict[str, List[str]]
            # Find consensus orthologs across genes in the cluster.
            orthologs_data = self.pan_super.gene_clusters_functions_summary_dict[gene_cluster_id]
            if self.protein_annotation_source == 'KEGG':
                anvio_source = 'KOfam'
            else:
                anvio_source = self.protein_annotation_source
            source_data = orthologs_data[anvio_source]
            if source_data['accession'] is None:
                # No ortholog from the source was assigned to the cluster.
                continue
            # Generate an ortholog object containing protein reference data.
            entry = {
                'gene_cluster_id': gene_cluster_id,
                'source': self.protein_annotation_source,
                'accession': source_data['accession'],
                'function': source_data['function'],
                'evalue': 0.0 # arbitrary, since annotation in majority of genes is used
            }
            if self.protein_annotation_source == 'KEGG':
                ortholog = protein.AnvioKOAnnotation(entry)
            else:
                raise ConfigError(
                    f"The protein annotation source, '{self.protein_annotation_source}', is not "
                    "recognized."
                )
            reactions = ortholog.get_reactions(self.source_db, cross_reference_dbs=self.cross_ref_dbs)
            if not reactions:
                # No reference reaction data could be assigned to the ortholog.
                continue
            cobrapy_gene_dict = self._get_cobrapy_gene_dict(gene_cluster_id, genome_gcids, reactions)
            cobrapy_genes: List = cobrapy_dict['genes']
            cobrapy_genes.append(cobrapy_gene_dict)
        self.progress.end()
        # Delete convenience attributes.
        delattr(self, 'cobrapy_reactions')
        delattr(self, 'cobrapy_metabolites')
        delattr(self, 'recorded_reactions')
        delattr(self, 'recorded_metabolites')
        return cobrapy_dict

    def _get_cobrapy_gene_dict(
        self,
        gene_cluster_id: str,
        genome_gcids: Dict[str, List[str]],
        reactions: List[protein.Reaction]
    ) -> Dict[str, Any]:
        cobrapy_gene_dict = COBRApyJSONStructure.get_gene_entry()
        cobrapy_gene_dict['id'] = gene_cluster_id
        for reaction in reactions:
            cobrapy_reaction_dict, already_recorded_reaction = self._get_cobrapy_reaction_dict(
                reaction, genome_gcids
            )
            if not already_recorded_reaction:
                self.cobrapy_reactions.append(cobrapy_reaction_dict)
                self.recorded_reactions[reaction.id] = cobrapy_reaction_dict
        return cobrapy_gene_dict

    def _get_cobrapy_reaction_dict(
        self,
        reaction: protein.Reaction,
        genome_gcids: Dict[str, List[str]]
    ) -> Tuple[Dict[str, Any], bool]:
        # Find the genomes encoding the reaction.
        genome_ids = set()
        for genome_id, gcids in genome_gcids.items():
            if gcids:
                genome_ids.add(genome_id)
        try:
            # There is already a JSON object for the reaction.
            cobrapy_reaction_dict = self.recorded_reactions[reaction.id]
            already_recorded_reaction = True
        except KeyError:
            cobrapy_reaction_dict = COBRApyJSONStructure.get_reaction_entry()
            already_recorded_reaction = False
        # List the genomes encoding the reaction in the JSON reaction object.
        cobrapy_reaction_annotation = cobrapy_reaction_dict['annotation']
        if already_recorded_reaction:
            recorded_genomes: List = cobrapy_reaction_annotation['genomes']
            cobrapy_reaction_annotation['genomes'] = sorted(
                genome_ids.union(set(recorded_genomes))
            )
            return cobrapy_reaction_dict, already_recorded_reaction
        cobrapy_reaction_dict['id'] = reaction.id
        cobrapy_reaction_dict['name'] = reaction.name
        reversibility = reaction.reversibility
        if not reversibility:
            cobrapy_reaction_dict['lower_bound'] = 0.0
        cobrapy_reaction_annotation['genomes'] = sorted(genome_ids)
        for chemical, coefficient, compartment in zip(
            reaction.chemicals, reaction.coefficients, reaction.compartments
        ):
            cobrapy_metabolite_dict, already_recorded_metabolite = self._get_cobrapy_metabolite_dict(
                chemical, coefficient, compartment, reversibility, genome_ids
            )
            metabolite_id = cobrapy_metabolite_dict['id']
            cobrapy_reaction_dict['metabolites'][metabolite_id] = float(coefficient)
            if not already_recorded_metabolite:
                self.cobrapy_metabolites.append(cobrapy_metabolite_dict)
                self.recorded_metabolites[metabolite_id] = cobrapy_metabolite_dict
        return cobrapy_reaction_dict, already_recorded_reaction

    def _get_cobrapy_metabolite_dict(
        self,
        chemical: protein.Chemical,
        coefficient: float,
        compartment: str,
        reversibility: bool,
        genome_ids: Set[str]
    ) -> Tuple[Dict[str, Any], bool]:
        if pd.isna(chemical.select_bigg_id):
            compound_id = chemical.modelseed_compound_id
        else:
            compound_id = chemical.select_bigg_id
        metabolite_id = f'{compound_id}_{compartment}'
        try:
            cobrapy_metabolite_dict = self.recorded_metabolites[metabolite_id]
            already_recorded_metabolite = True
        except KeyError:
            cobrapy_metabolite_dict = COBRApyJSONStructure.get_metabolite_entry()
            already_recorded_metabolite = False
        cobrapy_metabolite_annotation = cobrapy_metabolite_dict['annotation']
        if already_recorded_metabolite:
            if reversibility:
                recorded_consuming_genome_ids: List = cobrapy_metabolite_annotation['consuming_genomes']
                recorded_producing_genome_ids: List = cobrapy_metabolite_annotation['producing_genomes']
                cobrapy_metabolite_annotation['consuming_genomes'] = sorted(
                    genome_ids.union(set(recorded_consuming_genome_ids))
                )
                cobrapy_metabolite_annotation['producing_genomes'] = sorted(
                    genome_ids.union(set(recorded_producing_genome_ids))
                )
            elif coefficient < 0:
                recorded_consuming_genome_ids: List = cobrapy_metabolite_annotation['consuming_genomes']
                cobrapy_metabolite_annotation['consuming_genomes'] = sorted(
                    genome_ids.union(set(recorded_consuming_genome_ids))
                )
            elif coefficient > 0:
                recorded_producing_genome_ids: List = cobrapy_metabolite_annotation['producing_genomes']
                cobrapy_metabolite_annotation['producing_genomes'] = sorted(
                    genome_ids.union(set(recorded_producing_genome_ids))
                )
            return cobrapy_metabolite_dict, already_recorded_metabolite
        cobrapy_metabolite_dict['id'] = metabolite_id
        cobrapy_metabolite_dict['name'] = chemical.name if chemical.name else ""
        cobrapy_metabolite_dict['compartment'] = compartment
        cobrapy_metabolite_dict['charge'] = chemical.charge if chemical.charge else 0
        cobrapy_metabolite_dict['formula'] = chemical.formula if chemical.formula else ""
        if reversibility:
            cobrapy_metabolite_annotation['consuming_genomes'] = sorted(genome_ids)
            cobrapy_metabolite_annotation['producing_genomes'] = sorted(genome_ids)
        elif coefficient < 0:
            cobrapy_metabolite_annotation['consuming_genomes'] = sorted(genome_ids)
            cobrapy_metabolite_annotation['producing_genomes'] = []
        elif coefficient > 0:
            cobrapy_metabolite_annotation['consuming_genomes'] = []
            cobrapy_metabolite_annotation['producing_genomes'] = sorted(genome_ids)
        return cobrapy_metabolite_dict, already_recorded_metabolite

class ExternalGenome(ModelInput):
    """
    Metabolic model input data from an external genome.

    Parameters
    ==========
    contigs_db_path : str
        Stores data on the genome.
    protein_annotation_source : str
        The source of protein ortholog annotations. For now, only KEGG orthologs can be processed.
        Default, 'KEGG'.
    remove_missing_objective_metabolites : bool
        Metabolites in the biomass objective function are removed if they are not recorded as being
        produced or consumed by enzymes in the model.
    """
    def __init__(
        self,
        contigs_db_path: str,
        protein_annotation_source: str = 'KEGG',
        db_superdir: str = refdbs.ProteinReferenceDatabase.default_superdir,
        remove_missing_objective_metabolites: bool = False,
        run: Run = Run(),
        progress: Progress = Progress()
    ) -> None:
        self.contigs_db_path = contigs_db_path
        is_contigs_db(self.contigs_db_path)
        self.remove_missing_objective_metabolites = remove_missing_objective_metabolites
        self.run = run
        self.progress = progress
        self.set_protein_annotation_source(source=protein_annotation_source, db_superdir=db_superdir)
        self._init_contigs_super()

    def _init_contigs_super(self) -> None:
        args = Namespace()
        args.contigs_db = self.contigs_db_path
        self.contigs_super = ContigsSuperclass(args, r=run_quiet)
        if self.protein_annotation_source == 'KEGG':
            anvio_source = 'KOfam'
        else:
            anvio_source = self.protein_annotation_source
        self.contigs_super.init_functions(requested_sources=[anvio_source])

    def _get_cobrapy_json_dict(self) -> Dict[str, Any]:
        cobrapy_dict = COBRApyJSONStructure.get()
        self.cobrapy_reactions: List[Dict] = cobrapy_dict['reactions']
        self.cobrapy_reactions.append(COBRApyJSONStructure.get_ecoli_objective())
        self.cobrapy_metabolites: List[Dict] = cobrapy_dict['metabolites']
        self.recorded_reactions: Dict[str, Dict] = {}
        self.recorded_metabolites: Dict[str, Dict] = {}
        self.progress.new("Analyzing genes")
        num_analyzed = 0
        total = len(self.contigs_super.gene_function_calls_dict)
        for gcid, gene_dict in self.contigs_super.gene_function_calls_dict.items():
            self.progress.update(f"{num_analyzed} / {total}")
            num_analyzed += 1
            if self.protein_annotation_source == 'KEGG':
                anvio_source = 'KOfam'
            else:
                anvio_source = self.protein_annotation_source
            source_data = gene_dict[anvio_source]
            if not source_data:
                # No ortholog was assigned to the gene.
                continue
            # Generate an ortholog object containing protein reference data.
            entry = {
                'gene_callers_id': gcid,
                'source': self.protein_annotation_source,
                'accession': source_data[0],
                'function': source_data[1],
                'evalue': source_data[2]
            }
            if self.protein_annotation_source == 'KEGG':
                ortholog = protein.AnvioKOAnnotation(entry)
            else:
                raise ConfigError(
                    f"The protein annotation source, '{self.protein_annotation_source}', is not "
                    "recognized."
                )
            reactions = ortholog.get_protein_data(
                self.source_db, cross_reference_dbs=self.cross_ref_dbs
            ).reactions
            if not reactions:
                # No reference reaction data could be assigned to the ortholog.
                continue
            cobrapy_gene_dict = self._get_cobrapy_gene_dict(gcid, reactions)
            cobrapy_genes: List = cobrapy_dict['genes']
            cobrapy_genes.append(cobrapy_gene_dict)
        self.progress.end()
        # Delete convenience attributes.
        delattr(self, 'cobrapy_reactions')
        delattr(self, 'cobrapy_metabolites')
        delattr(self, 'recorded_reactions')
        delattr(self, 'recorded_metabolites')
        return cobrapy_dict

    def _get_cobrapy_gene_dict(self, gcid: int, reactions: List[protein.Reaction]) -> Dict[str, Any]:
        cobrapy_gene_dict = COBRApyJSONStructure.get_gene_entry()
        cobrapy_gene_dict['id'] = gcid
        for reaction in reactions:
            cobrapy_reaction_dict, already_recorded_reaction = self._get_cobrapy_reaction_dict(reaction)
            if not already_recorded_reaction:
                self.cobrapy_reactions.append(cobrapy_reaction_dict)
                self.recorded_reactions[reaction.id] = cobrapy_reaction_dict
        return cobrapy_gene_dict

    def _get_cobrapy_reaction_dict(self, reaction: protein.Reaction) -> Tuple[Dict[str, Any], bool]:
        try:
            # There is already a JSON object for the reaction.
            cobrapy_reaction_dict = self.recorded_reactions[reaction.id]
            return cobrapy_reaction_dict, True
        except KeyError:
            cobrapy_reaction_dict = COBRApyJSONStructure.get_reaction_entry()
        # List the genomes encoding the reaction in the JSON reaction object.
        cobrapy_reaction_dict['id'] = reaction.id
        cobrapy_reaction_dict['name'] = reaction.name
        reversibility = reaction.reversibility
        if not reversibility:
            cobrapy_reaction_dict['lower_bound'] = 0.0
        for chemical, coefficient, compartment in zip(
            reaction.chemicals, reaction.coefficients, reaction.compartments
        ):
            cobrapy_metabolite_dict, already_recorded_metabolite = self._get_cobrapy_metabolite_dict(
                chemical, compartment
            )
            metabolite_id = cobrapy_metabolite_dict['id']
            cobrapy_reaction_dict['metabolites'][metabolite_id] = float(coefficient)
            if not already_recorded_metabolite:
                self.cobrapy_metabolites.append(cobrapy_metabolite_dict)
                self.recorded_metabolites[metabolite_id] = cobrapy_metabolite_dict
        return cobrapy_reaction_dict, False

    def _get_cobrapy_metabolite_dict(
        self,
        chemical: protein.Chemical,
        compartment: str
    ) -> Tuple[Dict[str, Any], bool]:
        if pd.isna(chemical.select_bigg_id):
            compound_id = chemical.modelseed_compound_id
        else:
            compound_id = chemical.select_bigg_id
        metabolite_id = f'{compound_id}_{compartment}'
        try:
            cobrapy_metabolite_dict = self.recorded_metabolites[metabolite_id]
            return cobrapy_metabolite_dict, True
        except KeyError:
            cobrapy_metabolite_dict = COBRApyJSONStructure.get_metabolite_entry()
        cobrapy_metabolite_dict['id'] = metabolite_id
        cobrapy_metabolite_dict['name'] = chemical.name if chemical.name else ""
        cobrapy_metabolite_dict['compartment'] = compartment
        cobrapy_metabolite_dict['charge'] = chemical.charge if chemical.charge else 0
        cobrapy_metabolite_dict['formula'] = chemical.formula if chemical.formula else ""
        return cobrapy_metabolite_dict, False
