#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains classes for predicting metabolic exchanges via the reaction network and KGML processing subsystems."""

import os
from copy import deepcopy
from argparse import Namespace
from collections import defaultdict

import anvio
import anvio.kgml as kgml
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.kgmlnetworkops as nw
import anvio.reactionnetwork as rn
import anvio.filesnpaths as filesnpaths

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
    def __init__(self, args, format_args_for_single_estimator=False, run=run_quiet):
        """A base class to assign arguments to attributes for ExchangePredictor classes.
        
        PARAMETERS
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
        run.info("Predicting exchanges from KEGG Pathway Map walks", not self.no_pathway_walk)
        run.info("Predicting exchanges from merged Reaction Network", not self.pathway_walk_only)

        # establish expected outputs based on provided arguments
        self.output_types = ['potentially-exchanged-compounds', 'unique-compounds']
        if not self.no_pathway_walk:
            self.output_types.append('evidence')

        # to fool a single estimator into passing sanity checks, nullify multi estimator args here
        if format_args_for_single_estimator:
            self.databases = None
            self.external_genomes_file = None
            self.internal_genomes_file = None

    def sanity_check_args(self):
        """Here we sanity check all the common arguments to make sure they are sensibly set."""

        if self.use_equivalent_amino_acids and self.custom_equivalent_compounds_file:
            raise ConfigError("You can only provide one of `--use-equivalent-amino-acids` and `--custom-equivalent-compounds-file`. If you "
                        "want to equate L- and non-stereo-specific amino acid IDs, you should include them in your custom equivalents file. "
                        "(Pro tip: the easiest way to do that is to run this program once with just `--use-equivalent-amino-acids` to get the "
                        "AA equivalents file, and modify from there)")
        if self.no_pathway_walk and self.pathway_walk_only:
            raise ConfigError("The parameters --no-pathway-walk and --pathway-walk-only are mutually exclusive.")

    def setup_output_for_appending(self):
        """Initializes and returns a dictionary of AppendableFile objects, one for each output type"""

        output_dict = {}
        for typ in self.output_types:
            output_path = self.output_file_prefix + "-" + typ + ".txt"
            if filesnpaths.is_file_exists(output_path, dont_raise=True) and not anvio.FORCE_OVERWRITE:
                raise ConfigError(f"It seems like output files with your requested prefix already exist, for "
                                  f"example: {output_path}. Please delete the existing files or provide a "
                                  f"different output prefix, OR add the --force-overwrite flag to your command "
                                  f"(in which case, we will delete the existing files before starting a new one "
                                  f"to append to).")
            if filesnpaths.is_file_exists(output_path, dont_raise=True) and anvio.FORCE_OVERWRITE:
                os.remove(output_path)
                self.run.info_single(f"Removed existing output file at {output_path}")
            output_file_for_type = filesnpaths.AppendableFile(output_path, append_type=dict, fail_if_file_exists=True)
            output_dict[typ] = output_file_for_type

        return output_dict
        

class ExchangePredictorSingle(ExchangePredictorArgs):
    """Class for predicting exchanges between a single pair of genomes.
    
    PARAMETERS
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
        ExchangePredictorArgs.__init__(self, self.args, run=self.run)

        # INPUT SANITY CHECKS (for anything that was not already checked by the base class)
        if not self.contigs_db_1 or not self.contigs_db_2:
            raise ConfigError("The ExchangePredictorSingle class needs two contigs databases to work with.")

    def predict_exchanges(self, output_files_dictionary=None, return_data_dicts=False):
        """This is the driver function to predict metabolic exchanges between two genomes.
        
        PARAMETERS
        ==========
        output_files_dictionary : dictionary of output type (str), AppendableFile object pairs
            contains an initialized AppendableFile object to append output to for each output type
            (used in multi-mode to direct all output from several estimators to the same files)
        return_data_dicts : boolean
            If True, this function will return the prediction result dictionaries instead of 
            appending them to the output files
        """

        # LOAD DATA
        self.eq_compounds = {}
        self.run.info("Using equivalent L- and non-stereo-specific amino acid compound IDs", self.use_equivalent_amino_acids)
        if self.use_equivalent_amino_acids:
            self.eq_compounds = self.find_equivalent_amino_acids()
        elif self.custom_equivalent_compounds_file:
            self.run.info("Custom equivalent compounds file", self.custom_equivalent_compounds_file)
            self.eq_compounds = self.load_equivalent_compounds()

        if output_files_dictionary:
            self.output_file_dict = output_files_dictionary
        else:
            self.output_file_dict = self.setup_output_for_appending()
        
        self.load_reaction_networks()
        
        # MERGE NETWORKS
        self.merged, db_names = self.genomes_to_compare['A']['network']._merge_two_genome_networks(self.genomes_to_compare['B']['network'])
        self.genomes_to_compare['A']['name'] = db_names[0]
        self.genomes_to_compare['B']['name'] = db_names[1]
        self.run.info_single(f"Created a merged network with {len(self.merged.genes)} genes.")
        self.run.info("Number of metabolites in merged network to process", len(self.merged.metabolites))

        # SET UP DICTIONARIES
        self.map_kegg_ids_to_modelseed_ids()
        self.map_kegg_ids_to_compound_names()
        self.all_pathway_maps = self.merged._get_pathway_map_set(map_ids_to_exclude=MAPS_TO_EXCLUDE, id_selection_prefix = "00")
        # this will store the output of the KEGG Pathway Map walks
        # Dictionary structure: {compound_id (modelseed ID): {pathway_id: {organism_id: {fate: [chains]}}}}
        self.compound_to_pathway_walk_chains = {}

        # dictionaries to store output
        data_dicts = {t: {} for t in self.output_types}
        compounds_not_in_maps = set()
        self.processed_compound_ids = set() # this is how we'll make sure we don't process equivalent compounds twice

        if not self.no_pathway_walk:
            # STEP 1: PATHWAY MAP WALKS (TODO multithread: split pathway maps between workers)
            self.walk_all_pathway_maps()
            compounds_not_in_maps = set(self.merged.metabolites.keys()).difference(set(self.compound_to_pathway_walk_chains.keys()))

            # STEP 2: PREDICT EXCHANGES using pathway walk evidence
            pot_exchanged_pw, uniq_pw, data_dicts['evidence'] = self.predict_from_pathway_walks()
            data_dicts['potentially-exchanged-compounds'].update(pot_exchanged_pw)
            data_dicts['unique-compounds'].update(uniq_pw)
            
        else:
            compounds_not_in_maps = set(self.merged.metabolites.keys())
        self.run.info("Number of compounds from merged network not in Pathway Maps", len(compounds_not_in_maps))

        if not self.pathway_walk_only:
            # STEP 3: PREDICT EXCHANGES using reaction network
            pot_exchanged_rn, uniq_rn = self.predict_from_reaction_network(compounds_not_in_maps)
            compounds_in_both_dicts = set(data_dicts['potentially-exchanged-compounds'].keys()).intersection(pot_exchanged_rn.keys())
            if compounds_in_both_dicts:
                raise ConfigError(f"We found the following compound IDs that were found to be potentially-exchanged "
                            f"both by Pathway Walk and Reaction Network. Normally this shouldn't happen, and it means "
                            f"there is something weird with the set of compounds we gave to the reaction network. Regardless, "
                            f"we don't want to overwrite any of the Pathway Walk results, so we are stopping the show. Here "
                            f"are the affected compounds for debugging purposes: {compounds_in_both_dicts}")
            data_dicts['potentially-exchanged-compounds'].update(pot_exchanged_rn)
            compounds_in_both_dicts = set(data_dicts['unique-compounds'].keys()).intersection(uniq_rn.keys())
            if compounds_in_both_dicts:
                raise ConfigError(f"We found the following compound IDs that were found to be unique "
                            f"both by Pathway Walk and Reaction Network. Normally this shouldn't happen, and it means "
                            f"there is something weird with the set of compounds we gave to the reaction network. Regardless, "
                            f"we don't want to overwrite any of the Pathway Walk results, so we are stopping the show. Here "
                            f"are the affected compounds for debugging purposes: {compounds_in_both_dicts}")
            data_dicts['unique-compounds'].update(uniq_rn)

        self.run.warning(f"Identified {len(data_dicts['potentially-exchanged-compounds'])} potentially exchanged compounds and "
                         f"{len(data_dicts['unique-compounds'])} compounds unique to one genome.", header='OVERALL RESULTS', lc='green')
        
        # STEP 4: OUTPUT or RETURN
        if return_data_dicts:
            return data_dicts
        else:
            self.append_output_from_dicts(data_dicts)
            for typ, file_object in self.output_file_dict.items():
                self.run.info(f"Output with {typ}", file_object.path)
                file_object.close()       

    def find_equivalent_amino_acids(self, print_to_file=True, output_file_name="equivalent_amino_acids.txt"):
        """Looks through the ModelSEED compound table to identify L-amino acids and their non-stereo-specific counterparts.

        Returns a dictionary in which keys are compound IDs matched to their names and the equivalent compound ID. Each match is
        in the dictionary twice to enable O(1) lookups with either compound ID.

        PARAMETERS
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
        constructor = rn.Constructor(run=self.run, progress=self.progress)
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
                all_compounds_in_map = list(set(production_chains.keys()).union(set(consumption_chains.keys())))
                for compound in all_compounds_in_map:
                    if compound not in self.kegg_id_to_modelseed_id:
                        raise ConfigError(f"The merged reaction network doesn't contain a modelseed compound associated with {compound} in pathway map {pm}")
                    modelseed_id = self.kegg_id_to_modelseed_id[compound]
                    if modelseed_id in self.compound_to_pathway_walk_chains:
                        if anvio.DEBUG:
                            self.progress.reset()
                            all_kegg_compounds_with_same_modelseed_id = [x for x in self.kegg_id_to_modelseed_id if self.kegg_id_to_modelseed_id[x] == modelseed_id]
                            if len(all_kegg_compounds_with_same_modelseed_id) > 1:
                                self.run.warning(f"While processing KEGG compound {compound} in Pathway Map {pm}, we found that the associated "
                                                f"ModelSEED compound ID ({modelseed_id}) was already in the pathway walk dictionary. These are "
                                                f"all the KEGG compounds with that same ModelSEED ID: {','.join(all_kegg_compounds_with_same_modelseed_id)}")
                    else:
                        self.compound_to_pathway_walk_chains[modelseed_id] = {}
                    if pm not in self.compound_to_pathway_walk_chains[modelseed_id]:
                        self.compound_to_pathway_walk_chains[modelseed_id][pm] = {}
                    if g in self.compound_to_pathway_walk_chains[modelseed_id][pm]:
                        self.compound_to_pathway_walk_chains[modelseed_id][pm][g]['produce'] += production_chains[compound]
                        self.compound_to_pathway_walk_chains[modelseed_id][pm][g]['consume'] += consumption_chains[compound]
                    else:
                        self.compound_to_pathway_walk_chains[modelseed_id][pm][g] = {'produce': production_chains[compound] if compound in production_chains else None, 
                                                                            'consume': consumption_chains[compound] if compound in consumption_chains else None}
            processed_count += 1
            self.progress.update(f"{processed_count} / {num_pms_to_process} Pathway Maps")
            self.progress.increment(increment_to=processed_count)
        self.progress.end()
        self.run.info("Number of Pathway Maps processed", len(self.all_pathway_maps))
        self.run.info("Number of compounds processed from Pathway Maps", len(self.compound_to_pathway_walk_chains.keys()))

    def producer_consumer_decision_tree(self, genomes_can_produce, genomes_can_consume):
        """Given two genomes, decides which genome is the predicted producer/consumer of the current compound.
        If the compound is neither potentially-exchanged or unique to one genome, returns None so we can move on.

        Note that having a 'primary' producer genome doesn't mean the other genome cannot produce the compound 
        (and likewise for consumption), because we accept "2 producers, 1 consumer" or "1 producer, 2 consumers"
        as potential exchange scenarios. And in these cases, the genome that is alone in its role gets assigned
        as primary, and the other genome by default takes the other role even though there are technically 2 of them.

        PARAMETERS
        ==========
        genomes_can_produce : set of str
            which genome(s) have the capacity to produce the compound (how we figured that out doesn't matter at this point)
        genomes_can_consume : set of str
            which genome(s) have the capacity to consume the compound (again, how we know doesn't matter)

        RETURNS
        =======
        producer : str
            the 'primary' producer of the compound, if it is unique or potentially-exchanged. is None otherwise
        consumer : str
            the 'primary' consumer of the compound, if it is unique or potentially-exchanged. is None otherwise
        """

        producer = None
        consumer = None
        if len(genomes_can_produce) == 1 or len(genomes_can_consume) == 1:
                # when the same genome produces and/or consumes, compound is unique to one genome
                if genomes_can_produce == genomes_can_consume:
                    producer = genomes_can_produce.pop()
                    consumer = genomes_can_consume.pop()
                elif len(genomes_can_produce) == 0: # only one can consume
                    consumer = genomes_can_consume.pop()
                elif len(genomes_can_consume) == 0: # only one can produce
                    producer = genomes_can_produce.pop()
                else: # potentially-exchanged
                    if len(genomes_can_produce) == 2:
                        consumer = genomes_can_consume.pop()
                        genomes_can_produce.remove(consumer)
                        producer = genomes_can_produce.pop()
                    elif len(genomes_can_consume) == 2:
                        producer = genomes_can_produce.pop()
                        genomes_can_consume.remove(producer)
                        consumer = genomes_can_consume.pop()
                    else:
                        producer = genomes_can_produce.pop()
                        consumer = genomes_can_consume.pop()

        return producer, consumer

    def predict_exchange_from_pathway_walk(self, chains_for_all_equivalent_compounds):
        """Determines whether a compound is unique or potentially-exchanged, and if so returns the internal IDs of 
        the producer/consumer genome (or None, otherwise). Only works for 2 genomes.

        PARAMETERS
        ==========
        chains_for_all_equivalent_compounds : dictionary
            A subset of self.compound_to_pathway_walk_chains. Reminder of dict structure:
            {compound_id (modelseed ID): {pathway_id: {organism_id: {fate: [chains]}}}}

        RETURNS
        =======
        prod : str
            ID of the 'primary' producer. See producer_consumer_decision_tree().
        cons : str
            ID of the 'primary' consumer. See producer_consumer_decision_tree().
        """

        genomes_produce = set()
        genomes_consume = set()
        for cid, map_dict in chains_for_all_equivalent_compounds.items():
            for pm, genome_dict in map_dict.items():
                for g in genome_dict:
                    if genome_dict[g]['produce']: # if this is not None, there is at least one production chain
                        genomes_produce.add(g)
                    if genome_dict[g]['consume']: # if this is not None, there is at least one production chain
                        genomes_consume.add(g)

        prod, cons = self.producer_consumer_decision_tree(genomes_produce, genomes_consume)

        return prod, cons

    def get_longest_chains(self, chain_list):
        """Loops over all chains in provided list, computes length, and returns the maximum length and longest chain.
        
        PARAMETERS
        ==========
        chain_list : list of Chain objects (a class from kgmlnetworkops)
            Reaction chains from Pathway Map walk within a given reaction network.

        RETURNS
        =======
        max_length : int
            Length of the longest Chain
        longest : list of Chain objects
            all Chains with length max_length in the provided list
        """

        max_length = None
        longest = []

        if chain_list:
            chain_lengths = [len(chain.kgml_reactions) for chain in chain_list]
            max_length = max(chain_lengths)
            longest = [chain_list[i] for i in range(len(chain_list)) if chain_lengths[i] == max_length]

        return max_length, longest

    def get_max_overlap(self, reference_chains, comparison_chains):
        """Computes the longest reaction overlap between two lists of reaction chains.
        Also returns the reference chain that exhibits the maximum overlap with the comparison chains. Note 
        that if there are multiple solutions with the same max overlap, we only report the first 'longest' 
        chain with that overlap value.

        PARAMETERS
        ==========
        reference_chains : list of Chain objects (a class from kgmlnetworkops)
            The Chains we are using as our reference.
        comparison_chains : list of Chain objects
            The Chains we are comparing to the reference.

        RETURNS
        =======
        max_overlap : int
            Length of the longest overlap (highest number of shared reactions) between a reference chain 
            and a comparison chain
        longest_with_max_overlap : Chain object
            the first reference Chain (from the reference_chains list) that exhibits the maximum overlap with 
            a comparison chain.
        """

        max_overlap = None
        longest_with_max_overlap = None
        for sc in reference_chains:
            reaction_chain = [r.name for r in sc.kgml_reactions]
            for c in comparison_chains:
                compare_chain = [r.name for r in c.kgml_reactions]
                intersection = 0
                i = 0
                while (i < len(reaction_chain)) and (i < len(compare_chain)) and reaction_chain[i] == compare_chain[i]:
                    intersection += 1
                    i += 1
                if not max_overlap:
                    max_overlap = intersection
                    longest_with_max_overlap = sc
                elif intersection > max_overlap:
                    max_overlap = intersection
                    longest_with_max_overlap = sc

        return max_overlap, longest_with_max_overlap

    def get_pathway_walk_evidence_for_compound_in_map(self, production_chains_in_producer, production_chains_in_consumer, 
        consumption_chains_in_producer, consumption_chains_in_consumer):
        """Compares the pathway walk output to compute max reaction chain length, overlap, etc.
        
        PARAMETERS
        ==========
        production_chains_in_producer : list of Chain objects (a class from kgmlnetworkops)
            Reaction Chains to produce a compound in one Pathway Map in the primary producer
        production_chains_in_consumer : list of Chain objects
            Reaction Chains to produce a compound in one Pathway Map in the primary consumer
        consumption_chains_in_producer : list of Chain objects
            Reaction Chains to consume a compound in one Pathway Map in the primary producer
        consumption_chains_in_consumer : list of Chain objects
            Reaction Chains to consume a compound in one Pathway Map in the primary consumer

        RETURNS
        =======
        results : dict
            Contains evidence supporting the potential exchange of the compound, including length of
            reaction chains and length of overlap between the producer and consumer
        """

        max_production_length, longest_producer_chains = self.get_longest_chains(production_chains_in_producer)
        max_production_overlap, longest_with_max_overlap = self.get_max_overlap(longest_producer_chains, production_chains_in_consumer)

        if longest_with_max_overlap:
            longest_producer_chain_strings = {"reactions": [r.name for r in longest_with_max_overlap.kgml_reactions],
                                            "compounds": [c.name[4:] for c in longest_with_max_overlap.kgml_compound_entries]}
        else:
            longest_producer_chain_strings = {"reactions": [], "compounds": []}

        max_consumption_length, longest_consumer_chains = self.get_longest_chains(consumption_chains_in_consumer)
        max_consumption_overlap, longest_with_max_overlap = self.get_max_overlap(longest_consumer_chains, consumption_chains_in_producer)

        if longest_with_max_overlap:
            longest_consumer_chain_strings = {"reactions": [r.name for r in longest_with_max_overlap.kgml_reactions],
                                            "compounds": [c.name[4:] for c in longest_with_max_overlap.kgml_compound_entries]}
        else:
            longest_consumer_chain_strings = {"reactions": [], "compounds": []}

        results = {"max_production_length": max_production_length, # in producer network
                "max_consumption_length": max_consumption_length, # in consumer network
                "max_production_overlap": max_production_overlap, # between longest chain in producer network and all in consumer
                "max_consumption_overlap": max_consumption_overlap, # between longest chain in consumer network and all in producer
                "prop_production_overlap": max_production_overlap / max_production_length if max_production_length and max_production_overlap else None,
                "prop_consumption_overlap": max_consumption_overlap / max_consumption_length if max_consumption_length and max_consumption_overlap else None,
                "longest_chain_producer_strings": longest_producer_chain_strings,
                "longest_chain_consumer_strings": longest_consumer_chain_strings
        }

        return results

    def get_pathway_walk_evidence(self, reaction_chains_for_compound, producer, consumer):
        """Loops over all equivalent compounds and all pathway walk output to obtain per-map evidence.
        
        PARAMETERS
        ==========
        reaction_chains_for_compound : dict
            A subset of self.compound_to_pathway_walk_chains. Reminder of dict structure:
            {compound_id (modelseed ID): {pathway_id: {organism_id: {fate: [chains]}}}}
        producer : str
            ID of the organism predicted to be the 'primary' producer in a potential exchange
        consumer : str
            ID of the organism predicted to be the 'primary' consumer in a potential exchange
        """

        all_evidence_for_compound = {}
        for cid, map_dict in reaction_chains_for_compound.items():
            for pm, organism_dict in map_dict.items():
                prod_chains_in_producer = organism_dict[producer]['produce'] if producer in organism_dict else []
                prod_chains_in_consumer = organism_dict[consumer]['produce'] if consumer in organism_dict else []
                cons_chains_in_producer = organism_dict[producer]['consume'] if producer in organism_dict else []
                cons_chains_in_consumer = organism_dict[consumer]['consume'] if consumer in organism_dict else []
                if pm in all_evidence_for_compound:
                    raise ConfigError(f"While processing pathway walk compound {cid}, we noticed that the pathway map {pm} was already "
                                    f"processed, which means that an equivalent compound was also present in this pathway map. Not "
                                    f"how to handle this yet, so instead we say bye-bye :(")
                all_evidence_for_compound[pm] = self.get_pathway_walk_evidence_for_compound_in_map(prod_chains_in_producer, 
                                                            prod_chains_in_consumer, cons_chains_in_producer, cons_chains_in_consumer)

        return all_evidence_for_compound

    def predict_from_pathway_walks(self):
        """Loops over each compound in KEGG Pathway Maps and predicts whether it is potentially-exchanged or unique.

        RETURNS
        =======
        3 dictionaries: 
            potentially-exchanged compounds
            unique compounds
            evidence from Pathway Walk for potentially-exchanged compounds
        """

        potentially_exchanged_compounds = {}
        unique_compounds = {}
        pathway_walk_evidence = {}
        pathway_walk_dict_key = 0

        num_compounds_to_process = len(self.compound_to_pathway_walk_chains)
        processed_count = 0
        self.progress.new('Processing compounds in KEGG Pathway Maps', progress_total_items=num_compounds_to_process)
        # TODO multithread: split compounds between workers
        for compound_id in self.compound_to_pathway_walk_chains:
            if compound_id in self.processed_compound_ids:
                continue
            compound_reaction_chains = {compound_id: self.compound_to_pathway_walk_chains[compound_id]}
            if compound_id in self.eq_compounds:
                eq_comp = self.eq_compounds[compound_id]['equivalent_id']
                if eq_comp in self.compound_to_pathway_walk_chains:
                    compound_reaction_chains[eq_comp] = self.compound_to_pathway_walk_chains[eq_comp]
                self.processed_compound_ids.add(eq_comp)
            
            producer,consumer = self.predict_exchange_from_pathway_walk(compound_reaction_chains)
            if producer or consumer:
                compound_name = self.merged.metabolites[compound_id].modelseed_name
                producer_name = self.genomes_to_compare[producer]['name'] if producer else None
                consumer_name = self.genomes_to_compare[consumer]['name'] if consumer else None
                if producer == consumer or (not producer) or (not consumer): # unique to one genome
                    unique_compounds[compound_id] = {'compound_name': compound_name, 
                                                    'genomes': ",".join([x for x in set([producer_name,consumer_name]) if x]),
                                                    'produced_by': producer_name,
                                                    'consumed_by': consumer_name,
                                                    'prediction_method': 'Pathway_Map_Walk',
                                                    }
                else: # potentially-exchanged
                    potentially_exchanged_compounds[compound_id] = {'compound_name': compound_name,
                                                                    'genomes': ",".join([producer_name,consumer_name]),
                                                                    'produced_by': producer_name,
                                                                    'consumed_by': consumer_name,
                                                                    'prediction_method': 'Pathway_Map_Walk',
                                                                    }

                    per_map_evidence_for_compound = self.get_pathway_walk_evidence(compound_reaction_chains, producer, consumer)
                    # set up some variables to find the longest chain of reactions to use as the summary evidence for an exchange
                    # if there are multiple 'longest chains', we'll report the one with not-None or smallest overlap
                    overall_max_prior = None
                    overall_max_posterior = None
                    overall_overlap_prior = None
                    overall_overlap_posterior = None
                    prop_overlap_prior = None
                    prop_overlap_posterior = None
                    reported_map_prior = None
                    reported_map_posterior = None
                    for map_id, map_evidence in per_map_evidence_for_compound.items():
                        pathway_walk_evidence[pathway_walk_dict_key] = {'compound': compound_id,
                                                                        'compound_name': compound_name,
                                                                        'type': "production",
                                                                        'organism': producer_name,
                                                                        'pathway_map': map_id,
                                                                        'longest_reaction_chain_length': map_evidence["max_production_length"],
                                                                        'maximum_overlap': map_evidence["max_production_overlap"],
                                                                        'proportion_overlap': map_evidence["prop_production_overlap"],
                                                                        'longest_chain_reactions': ",".join(map_evidence["longest_chain_producer_strings"]["reactions"]),
                                                                        'longest_chain_compounds': ",".join(map_evidence["longest_chain_producer_strings"]["compounds"]),
                                                                        'longest_chain_compound_names': ",".join([self.get_compound_name_from_kegg_id(c) for c in map_evidence["longest_chain_producer_strings"]["compounds"]])}
                        pathway_walk_dict_key += 1
                        pathway_walk_evidence[pathway_walk_dict_key] = {'compound': compound_id,
                                                                        'compound_name': compound_name,
                                                                        'type': "consumption",
                                                                        'organism': consumer_name,
                                                                        'pathway_map': map_id,
                                                                        'longest_reaction_chain_length': map_evidence["max_consumption_length"],
                                                                        'maximum_overlap': map_evidence["max_consumption_overlap"],
                                                                        'proportion_overlap': map_evidence["prop_consumption_overlap"],
                                                                        'longest_chain_reactions': ",".join(map_evidence["longest_chain_consumer_strings"]["reactions"]),
                                                                        'longest_chain_compounds': ",".join(map_evidence["longest_chain_consumer_strings"]["compounds"]),
                                                                        'longest_chain_compound_names': ",".join([self.get_compound_name_from_kegg_id(c) for c in map_evidence["longest_chain_consumer_strings"]["compounds"]])}
                        pathway_walk_dict_key += 1
                        
                        def update_reported_pathway_evidence_for_prior():
                            """Updates variables like overall_max_prior with values from the current pathway map as the new 'best' 
                            evidence for an exchange."""
                            return map_evidence["max_production_length"], map_evidence["max_production_overlap"], \
                                map_evidence["prop_production_overlap"], map_id
                        def update_reported_pathway_evidence_for_posterior():
                            """Updates variables like overall_max_posterior with values from the current pathway map as the new 'best' 
                            evidence for an exchange."""
                            return map_evidence["max_consumption_length"], map_evidence["max_consumption_overlap"], map_evidence["prop_consumption_overlap"], map_id

                        # we want to find the longest chain of production reactions + the longest chain of consumption reactions
                        if (not overall_max_prior and map_evidence["max_production_length"]) or \
                            (overall_max_prior and map_evidence["max_production_length"] and overall_max_prior < map_evidence["max_production_length"]):
                            overall_max_prior, overall_overlap_prior, prop_overlap_prior, reported_map_prior = update_reported_pathway_evidence_for_prior()
                        # if we found a production chain of the same max length, report the one with smaller (or not None) overlap proportion
                        elif overall_max_prior and map_evidence["max_production_length"] and overall_max_prior == map_evidence["max_production_length"]:
                            if overall_overlap_prior is None and map_evidence["max_production_overlap"] is not None or \
                              prop_overlap_prior is None and map_evidence["prop_production_overlap"] is not None or \
                             (prop_overlap_prior and map_evidence["prop_production_overlap"] and prop_overlap_prior > map_evidence["prop_production_overlap"]):
                                overall_max_prior, overall_overlap_prior, prop_overlap_prior, reported_map_prior = update_reported_pathway_evidence_for_prior()
                        if (not overall_max_posterior and map_evidence["max_consumption_length"]) or \
                            (overall_max_posterior and map_evidence["max_consumption_length"] and overall_max_posterior < map_evidence["max_consumption_length"]):
                            overall_max_posterior, overall_overlap_posterior, prop_overlap_posterior, reported_map_posterior = update_reported_pathway_evidence_for_posterior()
                        # if we found a consumption chain of the same max length, report the one with smaller (or not None) overlap proportion
                        elif overall_max_posterior and map_evidence["max_consumption_length"] and overall_max_posterior == map_evidence["max_consumption_length"]:
                            if overall_overlap_posterior is None and map_evidence["max_consumption_overlap"] is not None or \
                             prop_overlap_posterior is None and map_evidence["prop_consumption_overlap"] is not None or \
                             (prop_overlap_posterior and map_evidence["prop_consumption_overlap"] and prop_overlap_posterior > map_evidence["prop_consumption_overlap"]):
                                overall_max_posterior, overall_overlap_posterior, prop_overlap_posterior, reported_map_posterior = update_reported_pathway_evidence_for_posterior()

                    longest_overall_chain = None
                    if overall_max_prior and overall_max_posterior:
                        longest_overall_chain = overall_max_prior + overall_max_posterior
                    elif overall_max_prior:
                        longest_overall_chain = overall_max_prior
                    elif overall_max_posterior:
                        longest_overall_chain = overall_max_posterior
                    potentially_exchanged_compounds[compound_id]['max_reaction_chain_length'] = longest_overall_chain
                    potentially_exchanged_compounds[compound_id]['max_production_chain_length'] = overall_max_prior
                    potentially_exchanged_compounds[compound_id]['max_consumption_chain_length'] = overall_max_posterior
                    potentially_exchanged_compounds[compound_id]['production_overlap_length'] = overall_overlap_prior
                    potentially_exchanged_compounds[compound_id]['consumption_overlap_length'] = overall_overlap_posterior
                    potentially_exchanged_compounds[compound_id]['production_overlap_proportion'] = prop_overlap_prior
                    potentially_exchanged_compounds[compound_id]['consumption_overlap_proportion'] = prop_overlap_posterior
                    potentially_exchanged_compounds[compound_id]['production_chain_pathway_map'] = reported_map_prior
                    potentially_exchanged_compounds[compound_id]['consumption_chain_pathway_map'] = reported_map_posterior
            self.processed_compound_ids.add(compound_id)
            processed_count += 1
            self.progress.update(f"{processed_count} / {num_compounds_to_process} compounds processsed")
            self.progress.increment(increment_to=processed_count)
        self.progress.end()
        self.run.info("Number of exchanged compounds predicted from KEGG Pathway Map walks", len(potentially_exchanged_compounds))
        self.run.info("Number of unique compounds predicted from KEGG Pathway Map walks", len(unique_compounds))

        return potentially_exchanged_compounds, unique_compounds, pathway_walk_evidence

    def predict_from_reaction_network(self, compound_set):
        """For each compound in the provided set, predicts whether it is potentially-exchanged or unique from the reaction network.

        RETURNS
        =======
        2 dictionaries: 
            potentially-exchanged compounds
            unique compounds
        """

        unique_compounds = {}
        potentially_exchanged_compounds = {}
        num_compounds_to_process = len(compound_set)
        processed_count=0
        self.progress.new('Processing compounds in reaction network', progress_total_items=num_compounds_to_process)
        for compound_id in compound_set:
            if compound_id in self.processed_compound_ids:
                continue
            compound_name = self.merged.metabolites[compound_id].modelseed_name
            if anvio.DEBUG:
                self.progress.reset()
                run.info_single(f"Working on compound {compound_id} ({compound_name})")
            sub_ids = [compound_id]
            if compound_id in self.eq_compounds:
                eq_comp = self.eq_compounds[compound_id]['equivalent_id']
                sub_ids.append(self.eq_compounds[compound_id]['equivalent_id'])
                self.processed_compound_ids.add(self.eq_compounds[compound_id]['equivalent_id'])
            sub_network = self.merged.subset_network(metabolites_to_subset=sub_ids)

            genomes_produce = set([])
            genomes_consume = set([])
            production_reactions = {g: {} for g in self.genomes_to_compare}
            consumption_reactions = {g: {} for g in self.genomes_to_compare}
            for rid, reaction in sub_network.reactions.items():
                # replace all equivalent IDs with the compound ID we are currently working on in the network
                if compound_id in self.eq_compounds:
                    eq_comp =self.eq_compounds[compound_id]['equivalent_id']
                    reaction.compound_ids = [compound_id if x == eq_comp else x for x in reaction.compound_ids]

                idx = reaction.compound_ids.index(compound_id)
                #TODO: make a dict of compounds with transport reactions and use as evidence in final output
                if reaction.compound_ids.count(compound_id) > 1: # likely a transport reaction, ignore
                    if anvio.DEBUG:
                        run.warning(f"Found {compound_id} more than once in {rid}. We are skipping this reaction",
                                        header="DEBUG", lc="yellow")
                    continue
                # which genomes produce this compound?
                elif reaction.coefficients[idx] > 0:
                    for g, genome_info in self.genomes_to_compare.items():
                        if genome_info['name'] in reaction.genomes_of_origin:
                            genomes_produce.add(g)
                            production_reactions[g][rid] = self.get_reaction_equation(reaction)

                # which genomes utilize this compound?
                elif reaction.coefficients[idx] < 0:
                    for g, genome_info in self.genomes_to_compare.items():
                        if genome_info['name'] in reaction.genomes_of_origin:
                            genomes_consume.add(g)
                            consumption_reactions[g][rid] = self.get_reaction_equation(reaction)

            def add_reactions_to_dict_for_compound(compound_dict):
                """Modifies the compound dictionary in place to add production and consumption reaction output"""
                for g in self.genomes_to_compare:
                    prod_rxn_ids = sorted(production_reactions[g].keys())
                    cons_rxn_ids = sorted(consumption_reactions[g].keys())
                    prod_rxn_eqs = [production_reactions[g][rid] for rid in prod_rxn_ids]
                    cons_rxn_eqs = [consumption_reactions[g][rid] for rid in cons_rxn_ids]
                    compound_dict[f"production_rxn_ids_{g}"] = " / ".join(prod_rxn_ids) if len(prod_rxn_ids) else None
                    compound_dict[f"consumption_rxn_ids_{g}"] = " / ".join(cons_rxn_ids) if len(cons_rxn_ids) else None
                    compound_dict[f"production_rxn_eqs_{g}"] = " / ".join(prod_rxn_eqs) if len(prod_rxn_eqs) else None
                    compound_dict[f"consumption_rxn_eqs_{g}"] = " / ".join(cons_rxn_eqs) if len(cons_rxn_eqs) else None
            
            producer,consumer = self.producer_consumer_decision_tree(genomes_produce, genomes_consume)
            if producer or consumer:
                producer_name = self.genomes_to_compare[producer]['name'] if producer else None
                consumer_name = self.genomes_to_compare[consumer]['name'] if consumer else None
                if producer == consumer or (not producer) or (not consumer): # unique to one genome
                    unique_compounds[compound_id] = {'compound_name': compound_name, 
                                                    'genomes': ",".join([x for x in set([producer_name,consumer_name])if x]),
                                                    'produced_by': producer_name,
                                                    'consumed_by': consumer_name,
                                                    'prediction_method': 'Reaction_Network_Subset',
                                                    }
                    add_reactions_to_dict_for_compound(unique_compounds[compound_id])
                else: # potentially-exchanged
                    potentially_exchanged_compounds[compound_id] = {'compound_name': compound_name,
                                                                    'genomes': ",".join([producer_name,consumer_name]),
                                                                    'produced_by': producer_name,
                                                                    'consumed_by': consumer_name,
                                                                    'prediction_method': 'Reaction_Network_Subset',
                                                                    }
                    add_reactions_to_dict_for_compound(potentially_exchanged_compounds[compound_id])

                    # fill these in with blanks to avoid issues later
                    potentially_exchanged_compounds[compound_id]['max_reaction_chain_length'] = None
                    potentially_exchanged_compounds[compound_id]['max_production_chain_length'] = None
                    potentially_exchanged_compounds[compound_id]['max_consumption_chain_length'] = None
                    potentially_exchanged_compounds[compound_id]['production_overlap_length'] = None
                    potentially_exchanged_compounds[compound_id]['consumption_overlap_length'] = None
                    potentially_exchanged_compounds[compound_id]['production_overlap_proportion'] = None
                    potentially_exchanged_compounds[compound_id]['consumption_overlap_proportion'] = None
                    potentially_exchanged_compounds[compound_id]['production_chain_pathway_map'] = None
                    potentially_exchanged_compounds[compound_id]['consumption_chain_pathway_map'] = None

            self.processed_compound_ids.add(compound_id)
            processed_count += 1
            self.progress.update(f"{processed_count} / {num_compounds_to_process} compounds processsed")
            self.progress.increment(increment_to=processed_count)

        self.progress.end()
        self.run.info("Number of exchanged compounds predicted from Reaction Network subset approach", len(potentially_exchanged_compounds))
        self.run.info("Number of unique compounds predicted from Reaction Network subset approach", len(unique_compounds))

        return potentially_exchanged_compounds, unique_compounds

    def append_output_from_dicts(self, output_dicts):
        """This function appends the output dictionaries to initialized AppendableFile objects in self.output_file_dict.
        
        PARAMETERS
        ==========
        output_dicts : dictionary of dictionaries
            Key is output type, and value is the data dictionary associated with that output type
        """

        output_header = ['compound_id', 'compound_name', 'genomes', 'produced_by', 'consumed_by', 'prediction_method']
        if self.add_reactions_to_output:
            output_header += [f"production_rxn_ids_{g}" for g in db_names] + [f"consumption_rxn_ids_{g}" for g in db_names] + \
            [f"production_rxn_eqs_{g}" for g in db_names] + [f"consumption_rxn_eqs_{g}" for g in db_names]
        exchange_header = deepcopy(output_header)
        if not self.no_pathway_walk:
            exchange_header += ['max_reaction_chain_length', 
                                'max_production_chain_length', 'production_overlap_length', 'production_overlap_proportion', 
                                'production_chain_pathway_map',
                                'max_consumption_chain_length', 'consumption_overlap_length', 'consumption_overlap_proportion', 
                                'consumption_chain_pathway_map'
                                ]

        for mode, file_obj in self.output_file_dict.items():  
            if mode not in output_dicts:
                raise ConfigError(f"Uh oh. You've requested to generate output of type '{mode}' but we don't "
                                  f"have a data dictionary associated with that type.")

            header_list = output_header
            if mode == 'potentially-exchanged-compounds':
                header_list = exchange_header

            if mode == 'evidence':
                file_obj.append(output_dicts[mode], do_not_write_key_column=True, none_value="None")
            else:
                file_obj.append(output_dicts[mode], headers=header_list, key_header='compound_id', none_value="None")