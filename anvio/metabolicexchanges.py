#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains classes for predicting metabolic exchanges via the reaction network and KGML processing subsystems."""

from copy import deepcopy
from argparse import Namespace
from collections import defaultdict

import anvio
import anvio.kgml as kgml
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.kgmlnetworkops as nw
import anvio.reactionnetwork as rn

from anvio.dbops import ContigsDatabase
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Iva Veseli"
__email__ = "iva.veseli@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)

MAPS_TO_EXCLUDE = set(["00470", # D-amino acid biosynthesis. This map mostly connects other maps.
                       "00195", # Photosynthesis. This map describes photosystems and complexes rather than reactions.
                       "00190", # Oxidative Phosphorylation. This map describes enzyme complexes rather than reactions.
                       "00543", # Exopolysaccharide biosynthesis. Does not have a RN-type KGML file.
])

class ExchangePredictorArgs():
    def __init__(self, args, format_args_for_single_estimator=False, run=run, progress=progress):
        """A base class to assign arguments to attributes for ExchangePredictor classes.
        
        Parameters
        ==========
        format_args_for_single_estimator: bool
            This is a special case where an args instance is generated to be passed to the
            'single' predictor subclass from within the 'multi' predictor subclass. Setting 
            it to True ensures that any args specific to the multi subclass are removed to 
            avoid sanity check issues.
        """

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_prefix = A('output_file_prefix')
        self.modelseed_data_dir = A('modelseed_data_dir')
        self.kegg_data_dir = A('kegg_data_dir')
        self.num_threads = A('num_threads')
        self.use_equivalent_amino_acids = A('use_equivalent_amino_acids')
        self.custom_equivalent_compounds_file = A('custom_equivalent_compounds_file')
        self.maximum_gaps = A('maximum_gaps')
        self.add_reactions_to_output = A('add_reactions_to_output')
        self.no_pathway_walk = A('no_pathway_walk')
        self.pathway_walk_only = A('pathway_walk_only')

        self.sanity_check_args()

        # PRINT INFO for arguments common to subclasses
        self.run.info("Predicting exchanges from KEGG Pathway Map walks", not args.no_pathway_walk)
        self.run.info("Predicting exchanges from merged Reaction Network", not args.pathway_walk_only)

    def sanity_check_args(self):
        """Here we sanity check all the common arguments to make sure they are sensibly set."""

        if self.use_equivalent_amino_acids and self.custom_equivalent_compounds_file:
            raise ConfigError("You can only provide one of `--use-equivalent-amino-acids` and `--custom-equivalent-compounds-file`. If you "
                        "want to equate L- and non-stereo-specific amino acid IDs, you should include them in your custom equivalents file. "
                        "(Pro tip: the easiest way to do that is to run this program once with just `--use-equivalent-amino-acids` to get the "
                        "AA equivalents file, and modify from there)")
        if self.no_pathway_walk and self.pathway_walk_only:
            raise ConfigError("The parameters --no-pathway-walk and --pathway-walk-only are mutually exclusive.")
        

class ExchangePredictorSingle(ExchangePredictorArgs):
    """Class for predicting exchanges between a single pair of genomes.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-predict-metabolic-exchanges
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_1 = A('contigs_db_1')
        self.contigs_db_2 = A('contigs_db_2')

        # INIT BASE CLASS to format common arguments
        ExchangePredictorArgs.__init__(self, self.args)

        # INPUT SANITY CHECKS (for anything that was not already checked by the base class)
        if not self.contigs_db_1 or not self.contigs_db_2:
            raise ConfigError("The ExchangePredictorSingle class needs two contigs databases to work with.")

    def predict_exchanges(self):
        """This is the driver function to predict metabolic exchanges between two genomes."""

        # LOAD DATA
        self.eq_compounds = {}
        self.run.info("Using equivalent L- and non-stereo-specific amino acid compound IDs", self.use_equivalent_amino_acids)
        if self.use_equivalent_amino_acids:
            self.eq_compounds = self.find_equivalent_amino_acids()
        elif self.custom_equivalent_compounds_file:
            self.run.info("Custom equivalent compounds file", self.custom_equivalent_compounds_file)
            self.eq_compounds = self.load_equivalent_compounds()

        self.load_reaction_networks()
        
        # MERGE NETWORKS
        self.merged, db_names = self.genomes_to_compare['A']['network']._merge_two_genome_networks(self.genomes_to_compare['B']['network'])
        self.genomes_to_compare['A']['name'] = db_names[0]
        self.genomes_to_compare['B']['name'] = db_names[1]
        self.run.info_single(f"Created a merged network with {len(self.merged.genes)} genes.")
        self.run.info("Number of metabolites in merged network to process", len(self.merged.metabolites))

        # SETUP DICTIONARIES
        self.map_kegg_ids_to_modelseed_ids()
        self.map_kegg_ids_to_compound_names()
        self.all_pathway_maps = self.merged._get_pathway_map_set(map_ids_to_exclude=MAPS_TO_EXCLUDE, id_selection_prefix = "00")
        # this will store the output of the KEGG Pathway Map walks
        # Dictionary structure: {compound_id (modelseed ID): {pathway_id: {organism_id: {fate: [chains]}}}}
        self.compound_to_pathway_walk_chains = defaultdict(dict)

        if not self.no_pathway_walk:
            # STEP 1: PATHWAY MAP WALKS (TODO multithread: split pathway maps between workers)
            self.walk_all_pathway_maps()
        

    def find_equivalent_amino_acids(self, print_to_file=True, output_file_name="equivalent_amino_acids.txt"):
        """Looks through the ModelSEED compound table to identify L-amino acids and their non-stereo-specific counterparts.

        Returns a dictionary in which keys are compound IDs matched to their names and the equivalent compound ID. Each match is
        in the dictionary twice to enable O(1) lookups with either compound ID.

        Parameters
        ==========
        print_to_file : Boolean
            if True, the dictionary of equivalent compounds will be written to a tab-delimited output file
        output_file_name : string
            the name of the output file to print to (if desired)
        """

        MS_db = rn.ModelSEEDDatabase(modelseed_dir=self.modelseed_data_dir)
        comps = MS_db.compounds_table.dropna(subset="name")

        aa_list = anvio.constants.amino_acids_long + ["Selenocysteine", "Pyrrolysine"]
        equivalent_AAs = {}
        for a in aa_list:
            eqs = comps[(comps.name == a) | (comps.name == f"L-{a}")]['name']
            id_list = eqs.index.to_list()
            if not len(id_list) == 2:
                if anvio.DEBUG:
                    compound_list = [f"{i} ({eqs.loc[i]})" for i in id_list]
                    run.warning(f"While looking for equivalent compound IDs for the amino acid {a}, we didn't find the "
                                f"expected 2 matching compounds '{a}' and 'L-{a}'. Instead, here is what we found: "
                                f"{', '.join(compound_list)}. No equivalencies will be stored for this amino acid.", header='DEBUG', lc='yellow')
                continue
            equivalent_AAs[id_list[0]] = {'equivalent_id': id_list[1], 'name': eqs.loc[id_list[0]], 'equivalent_name': eqs.loc[id_list[1]]}
            equivalent_AAs[id_list[1]] = {'equivalent_id': id_list[0], 'name': eqs.loc[id_list[1]], 'equivalent_name': eqs.loc[id_list[0]]}

        if print_to_file:
            utils.store_dict_as_TAB_delimited_file(equivalent_AAs, output_file_name, key_header='compound_id',
                        headers = ['compound_id', 'equivalent_id', 'name', 'equivalent_name'])
            self.run.info("File of equivalent amino acid compound IDs", output_file_name)

        return equivalent_AAs

    def load_equivalent_compounds(self):
        """Loads the provided table to equivalent compounds and returns it as a dictionary"""

        eq_comp_dict = utils.get_TAB_delimited_file_as_dictionary(self.custom_equivalent_compounds_file, expected_fields=['compound_id','equivalent_id'])
        # make sure each match is in there twice for easy lookups from either direction
        pairs_to_add = {}
        for c, match in eq_comp_dict.items():
            c2 = match['equivalent_id']
            if c2 not in eq_comp_dict:
                pairs_to_add[c2] = {'equivalent_id': c}
            elif eq_comp_dict[c2]['equivalent_id'] != c:
                raise ConfigError(f"While parsing your file of custom equivalent compound IDs, we found a pair of equivalents that don't match: "
                                    f"compound {c} is paired with {c2}, while {c2} is paired with {eq_comp_dict[c2]['equivalent_id']}. You should probably "
                                    f"fix that.")
        if pairs_to_add:
            eq_comp_dict.update(pairs_to_add)
            if anvio.DEBUG:
                self.run.warning(f"Found {len(pairs_to_add)} inverse pairs of equivalent compounds that were not in the equivalent compounds file. "
                                    f"These have been added to the dictionary of equivalent compounds. Here is the set of missing pairs: {pairs_to_add}. ",
                                    header="DEBUG", lc="yellow")

        return eq_comp_dict

    def load_reaction_networks(self):
        """Loads the reaction network for each genome and establishes the self.genomes_to_compare attribute of genome information."""

        self.genomes_to_compare = {'A': {'contigs_db_path': self.contigs_db_1},
                       'B': {'contigs_db_path': self.contigs_db_2}}
        constructor = rn.Constructor()
        for g in self.genomes_to_compare:
            self.run.info("Loading reaction network from database", self.genomes_to_compare[g]['contigs_db_path'])
            self.genomes_to_compare[g]['network'] = constructor.load_network(contigs_db=self.genomes_to_compare[g]['contigs_db_path'], quiet=True)

    def map_kegg_ids_to_modelseed_ids(self):
        """Creates the self.kegg_id_to_modelseed_id dictionary. 
        The keys are KEGG compound IDs and the values are ModelSEED IDs for all compounds in the provided network.
        """

        self.kegg_id_to_modelseed_id = {}
        for mid, compound in self.merged.metabolites.items():
            for kid in compound.kegg_aliases:
                self.kegg_id_to_modelseed_id[kid] = mid
        
    def map_kegg_ids_to_compound_names(self):
        """Creates the self.id_to_name_dict dictionary.
        The keys are KEGG compound IDs and the values are compound names for all compounds in the provided network.
        """

        self.kegg_id_to_compound_name = {}
        for mid, compound in self.merged.metabolites.items():
            for kid in compound.kegg_aliases:
                self.kegg_id_to_compound_name[kid] = compound.modelseed_name

    def get_compound_name_from_kegg_id(self, kid):
        """A safer dictionary access that returns None if a given compound is not in self.kegg_id_to_compound_name"""

        if kid in self.kegg_id_to_compound_name:
            return self.kegg_id_to_compound_name[kid]
        else:
            return "None"

    def get_reaction_equation(self, reaction_value):
        """Looks up all compound names in the merged network and returns an equation for the chemical reaction with those names."""

        name_list = [self.merged.metabolites[c].modelseed_name for c in reaction_value.compound_ids]
        return rn.get_chemical_equation(reaction_value, use_compound_names=name_list, ignore_compartments = True)

    def get_args_for_pathway_walker(self, net, pathway_map, fate, gaps):
        """Returns a Namespace with arguments for KGMLNetworkWalker"""

        walker_args = Namespace()
        walker_args.network = net
        walker_args.kegg_pathway_number = pathway_map
        walker_args.compound_fate = fate
        walker_args.max_gaps = gaps
        walker_args.keep_intermediate_chains = True
        walker_args.verbose = False
        return walker_args

    def walk_all_pathway_maps(self):
        """Loops over all KEGG Pathway Maps in the merged network and fills self.compound_to_pathway_walk_chains."""

        num_pms_to_process = len(self.all_pathway_maps)
        processed_count = 0
        self.progress.new('Walking through KEGG Pathway Maps', progress_total_items=num_pms_to_process)
        for pm in self.all_pathway_maps:
            for g in self.genomes_to_compare: 
                wargs = self.get_args_for_pathway_walker(self.genomes_to_compare[g]['network'], pm, fate='produce', gaps=self.maximum_gaps)
                walker = nw.KGMLNetworkWalker(wargs)
                production_chains = walker.get_chains()
                walker.compound_fate = 'consume'
                consumption_chains = walker.get_chains()
                for compound in set(production_chains.keys()).union(set(consumption_chains.keys())):
                    if compound not in self.kegg_id_to_modelseed_id:
                        raise ConfigError(f"We didn't find a modelseed compound associated with {compound} in pathway map {pm}")
                    modelseed_id = self.kegg_id_to_modelseed_id[compound]
                    
                    if pm not in self.compound_to_pathway_walk_chains[modelseed_id]:
                        self.compound_to_pathway_walk_chains[modelseed_id][pm] = {}
                    self.compound_to_pathway_walk_chains[modelseed_id][pm][g] = {'produce': production_chains[compound] if compound in production_chains else None, 
                                                                            'consume': consumption_chains[compound] if compound in consumption_chains else None}
            processed_count += 1
            self.progress.update(f"{processed_count} / {num_pms_to_process} Pathway Maps")
            self.progress.increment(increment_to=processed_count)
        self.progress.end()