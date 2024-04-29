# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Common functions and classes for SCG/TRNA taxonomy.
"""

import os
import sys
import gzip
import copy
import hashlib
import argparse
import numpy as np
import pandas as pd

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from collections import OrderedDict, Counter

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.drivers.diamond import Diamond
from anvio.tables.scgtaxonomy import TableForSCGTaxonomy
from anvio.tables.trnataxonomy import TableForTRNATaxonomy
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.dbops import ContigsSuperclass, ContigsDatabase, ProfileSuperclass, ProfileDatabase

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print

HASH = lambda d: str(hashlib.sha224(''.join([str(d[level]) for level in constants.levels_of_taxonomy]).encode('utf-8')).hexdigest()[0:8])

class TerminologyHelper(object):
    def __init__(self):
        proper_foci = ['trnas', 'scgs']

        if not hasattr(self.ctx, 'focus'):
            raise ConfigError("You are lost :/ Any class initializes from shared taxonomy classes "
                              "must have a `self.ctx.focus` variable with either of these values: %s." % ', '.join(proper_foci))

        if self.ctx.focus not in proper_foci:
            raise ConfigError("Unknown focus %s :/" % (self.ctx.focus))

        self.scgs_focus = self.ctx.focus == 'scgs'
        self.trna_focus = self.ctx.focus == 'trnas'

        C = lambda x, y: x if self.scgs_focus else y
        self._DELIVERABLE = C("SCG taxonomy", "tRNA taxonomy")
        self._ITEM = C("SCG", 'anticodon')
        self._ITEMS = C("SCGs", 'anticodons')
        self._SUPPORTING_ITEMS = C("supporting_scgs", 'supporting_anticodons')
        self._TOTAL_ITEMS = C("total_scgs", 'total_anticodons')
        self._SOURCE_DATA = C("single-copy core gene", "tRNA gene")
        self._SETUP_PROGRAM = C("anvi-setup-scg-taxonomy", "anvi-setup-trna-taxonomy")
        self._COMPUTE_COVS_FLAG = C("--compute-scg-coverages", "--compute-anticodon-coverages")
        self._ITEM_FOR_METAGENOME_MODE_PARAM = C("--scg-name-for-metagenome-mode", "--anticodon-for-metagenome-mode")
        self._VARIABLE_NAME_IN_TABLE = C("gene_name", "anticodon")


class AccessionIdToTaxonomy(object):
    """A base classs that populates `self.accession_to_taxonomy_dict`"""

    def __init__(self):
        self.accession_to_taxonomy_dict = {}

        if not os.path.exists(self.accession_to_taxonomy_file_path):
            return None

        letter_to_level = dict([(l.split('_')[1][0], l) for l in self.levels_of_taxonomy])

        self.progress.new("Reading the accession to taxonomy file")
        self.progress.update('...')

        with gzip.open(self.accession_to_taxonomy_file_path, 'rb') as taxonomy_file:
            for line in taxonomy_file.readlines():
                line = line.decode('utf-8')

                if line.startswith('#'):
                    continue

                accession, taxonomy_text = line.strip('\n').split('\t')
                # taxonomy_text kinda looks like these:
                #
                #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Alcaligenes;s__Alcaligenes faecalis_C
                #    d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus_B;s__Enterococcus_B faecalis
                #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter sp1
                #    d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A cereus_AU
                #    d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Tissierellales;f__Helcococcaceae;g__Finegoldia;s__Finegoldia magna_H

                d = {}
                for letter, taxon in [e.split('__', 1) for e in taxonomy_text.split(';')]:
                    if letter in letter_to_level:
                        d[letter_to_level[letter]] = taxon
                    else:
                        self.run.warning("Some weird letter found in '%s' :(" % taxonomy_text)

                self.accession_to_taxonomy_dict[accession] = d

        # let's add one more accession for all those missing accessions
        self.accession_to_taxonomy_dict['unknown_accession'] = dict([(taxon, None) for taxon in self.levels_of_taxonomy])

        self.progress.end()


class TaxonomyEstimatorSingle(TerminologyHelper):
    def __init__(self, skip_init=False):
        TerminologyHelper.__init__(self)

        # update variables
        # FIXME: there is a better way to do it.
        self.ctx.items = self.ctx.SCGs if self.scgs_focus else self.ctx.anticodons
        self.report_item_frequencies_path = self.report_scg_frequencies_path if self.scgs_focus else self.report_anticodon_frequencies_path
        self.compute_item_coverages = self.compute_scg_coverages if self.scgs_focus else self.compute_anticodon_coverages
        self.item_name_for_metagenome_mode = self.scg_name_for_metagenome_mode if self.scgs_focus else self.anticodon_for_metagenome_mode
        self.per_item_output_file = self.per_scg_output_file if self.scgs_focus else self.per_anticodon_output_file

        # these dictionaries that will be initiated later
        self.contigs_db_project_name = "Unknown"
        self.item_name_to_gene_caller_id_dict = {}
        self.frequency_of_items_with_taxonomy = {}
        self.gene_callers_id_to_item_taxonomy_dict = {}
        self.split_name_to_gene_caller_ids_dict = {}
        self.gene_callers_id_to_split_name_dict = {}
        self.sample_names_in_profile_db = None

        self.initialized = False

        self.run.info('Contigs DB', self.contigs_db_path)
        if self.profile_db_path:
            self.run.info('Profile DB', self.profile_db_path, mc="green")
        self.run.info('Metagenome mode', self.metagenome_mode)
        if self.metagenome_mode:
            self.run.info(f'{self._ITEM} for metagenome', self.item_name_for_metagenome_mode)

        if not skip_init:
            self.init()


    def init(self):
        self.init_items_data()

        if self.report_item_frequencies_path:
            with open(self.report_item_frequencies_path, 'w') as output:
                for item_name, frequency in self.frequency_of_items_with_taxonomy.items():
                    output.write("%s\t%d\n" % (item_name, frequency))

            self.run.info(f'{self._ITEM} frequencies in contigs db', self.report_item_frequencies_path, nl_before=1)
            sys.exit()

        if self.profile_db_path:
            self.sample_names_in_profile_db = ProfileDatabase(self.profile_db_path).samples

        self.initialized = True


    def init_items_data(self):
        """Initialize {self._ITEM} taxonomy for the entire contigs database"""

        if not self.contigs_db_path:
            return None

        self.progress.new('Initializing')
        self.progress.update(f'{self._ITEM} taxonomy dictionary')

        for item_name in self.ctx.items:
            self.item_name_to_gene_caller_id_dict[item_name] = set([])

        if self.scgs_focus:
            anvio_taxonomy_table_name = t.scg_taxonomy_table_name
        elif self.trna_focus:
            anvio_taxonomy_table_name = t.trna_taxonomy_table_name
        else:
            anvio_taxonomy_table_name = None

        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        self.contigs_db_project_name = contigs_db.meta['project_name']
        anvio_taxonomy_table = contigs_db.db.get_table_as_dict(anvio_taxonomy_table_name)
        genes_in_splits = contigs_db.db.get_some_columns_from_table(t.genes_in_splits_table_name, "split, gene_callers_id")
        min_contig_length_in_contigs_db = contigs_db.db.get_max_value_in_column(t.contigs_info_table_name, "length", return_min_instead=True)
        contigs_db.disconnect()

        # this is important. before we begin, we need to filter out gene caller ids and splits from main dictionaries if
        # they shouldn't be there. read the warning below to see the utility of this step.
        if self.profile_db_path and self.compute_item_coverages:
            split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.profile_db_path))
            split_names_in_contigs_db = set([tpl[0] for tpl in genes_in_splits])
            splits_missing_in_profile_db = split_names_in_contigs_db.difference(split_names_in_profile_db)

            min_contig_length_in_profile_db = ProfileDatabase(self.profile_db_path).meta['min_contig_length']

            if len(splits_missing_in_profile_db):
                self.progress.reset()
                self.run.warning(f"Please note that anvi'o found {pp(len(split_names_in_contigs_db))} splits in your contigs database. "
                                 f"But only {pp(len(split_names_in_profile_db))} of them appeared ot be in the profile database. As a "
                                 f"result, anvi'o will now remove the {pp(len(splits_missing_in_profile_db))} splits that occur only in "
                                 f"the contigs db from all downstream analyses here (if you didn't use the flag `{self._COMPUTE_COVS_FLAG}` "
                                 f"this wouldn't have been necessary, but with the current settings this is really the best for everyone). "
                                 f"Where is this difference coming from though? Well. This is often the case because the 'minimum contig "
                                 f"length parameter' set during the `anvi-profile` step can exclude many contigs from downstream analyses "
                                 f"(often for good reasons, too). For instance, in your case the minimum contig length goes as low as "
                                 f"{pp(min_contig_length_in_contigs_db)} nts in your contigs database. Yet, the minimum contig length set "
                                 f"in the profile databaes is {pp(min_contig_length_in_profile_db)} nts. Hence the difference. Anvi'o "
                                 f"hopes that this explaines some things.")

                self.progress.update("Removing %s splits missing form the profile db" % pp(len(splits_missing_in_profile_db)))
                genes_in_splits = [tpl for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db]

                # so now we know the final list of split names and gene caller ids as they are stored in the updated
                # `genes_in_splits` variable. time to clean up the `anvio_taxonomy_table` dictionary as well.
                final_set_of_gene_caller_ids = set([tpl[1] for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db])
                entry_ids_to_remove = [entry for entry in anvio_taxonomy_table if anvio_taxonomy_table[entry]['gene_callers_id'] not in final_set_of_gene_caller_ids]
                [anvio_taxonomy_table.pop(e) for e in entry_ids_to_remove]

        # NOTE: This will modify the taxonomy strings read from the contigs database. see the
        # function header for `trim_taxonomy_dict_entry` for more information.
        if self.simplify_taxonomy_information:
            self.progress.update(f'{self._DELIVERABLE} dicts ... trimming main')
            for key in anvio_taxonomy_table:
                anvio_taxonomy_table[key] = self.trim_taxonomy_dict_entry(anvio_taxonomy_table[key])

        self.progress.update(f'{self._DELIVERABLE} dicts ... building g->tax')
        for entry in anvio_taxonomy_table.values():
            gene_callers_id = entry['gene_callers_id']
            self.gene_callers_id_to_item_taxonomy_dict[gene_callers_id] = entry

        self.progress.update(f'{self._DELIVERABLE} dicts ... building tax->g')
        for entry in self.gene_callers_id_to_item_taxonomy_dict.values():
            item_gene_name = entry[self._VARIABLE_NAME_IN_TABLE]
            gene_callers_id = entry['gene_callers_id']
            self.item_name_to_gene_caller_id_dict[item_gene_name].add(gene_callers_id)

        self.progress.update(f'{self._DELIVERABLE} dicts ... building s->g + g->s')
        for split_name, gene_callers_id in genes_in_splits:
            if gene_callers_id not in self.gene_callers_id_to_item_taxonomy_dict:
                continue

            if split_name not in self.split_name_to_gene_caller_ids_dict:
                self.split_name_to_gene_caller_ids_dict[split_name] = set()

            self.split_name_to_gene_caller_ids_dict[split_name].add(gene_callers_id)
            self.gene_callers_id_to_split_name_dict[gene_callers_id] = split_name

        self.progress.end()

        self.frequency_of_items_with_taxonomy = OrderedDict(sorted([(g, len(self.item_name_to_gene_caller_id_dict[g])) for g in self.item_name_to_gene_caller_id_dict], key = lambda x: x[1], reverse=True))

        if self.metagenome_mode or anvio.DEBUG:
            self.run.info_single(f"A total of %s {self._SOURCE_DATA}s with taxonomic affiliations were successfully initialized "
                                 f"from the contigs database 🎉 Following shows the frequency of these {self._ITEMS}: %s." % \
                                            (pp(len(self.gene_callers_id_to_item_taxonomy_dict)),
                                             ', '.join(["%s (%d)" % (g, self.frequency_of_items_with_taxonomy[g]) \
                                                                for g in self.frequency_of_items_with_taxonomy])), nl_before=1)


    def trim_taxonomy_dict_entry(self, taxonomy_dict_entry):
        """ Remove excess information from taxonomy information.

        The purpose of this is to give an option to the user to simplify GTDB names, that
        will have a text information for every level of taxonomy depending on what branches
        genomes fit, but it is not always helpful to the user. Such as this one:

             t_domain Bacteria
             t_phylum Firmicutes
             t_class Clostridia
             t_order Monoglobales
             t_family UBA1381
             t_genus CAG-41
             t_species CAG-41 sp900066215

         in this case the user may want to get this instead:

             t_domain Bacteria
             t_phylum Firmicutes
             t_class Clostridia
             t_order Monoglobales
             t_family None
             t_genus None
             t_species None

         So this function will take a taxonomy dict entry , and will return a simplified
         version of it if trimming is applicable.

        Paremeters
        ==========
        taxonomy_dict_entry: dict
            a dictionary that contains keys for all taxon names. such as this one:
                {[...],
                 't_domain': 'Bacteria',
                 't_phylum': 'Firmicutes',
                 't_class': 'Clostridia',
                 't_order': 'Oscillospirales',
                 't_family': 'Acutalibacteraceae',
                 't_genus': 'Ruminococcus',
                 't_species': 'Ruminococcus sp002491825'
                }
         """

        # for optimization, these letters should all have three characters. if that
        # behavior needs to change, the code down below must be updates. the purpose
        # of this is not to have a comprehensive list of EVERY single GTDB-specific
        # clade designations, but to make sure teh vast majority of names are covered.
        GTDB_specific_clade_prefixes = ['CAG', 'GCA', 'UBA', 'FUL', 'PAL', '2-0', 'Fen', 'RF3', 'TAN']

        taxonomic_levels_to_nullify = []

        for taxonomic_level in self.ctx.levels_of_taxonomy[::-1]:
            if not taxonomy_dict_entry[taxonomic_level]:
                continue

            if taxonomic_level == 't_species':
                species_name = taxonomy_dict_entry[taxonomic_level].split(' ')[1]
                try:
                    int(species_name[2])
                    taxonomic_levels_to_nullify.append(taxonomic_level)
                except:
                    None
            else:
                if taxonomy_dict_entry[taxonomic_level][0:3] in GTDB_specific_clade_prefixes:
                    taxonomic_levels_to_nullify.append(taxonomic_level)

        # this is the best way to make sure we are not going to nullify order, but leave behind a family name.
        if taxonomic_levels_to_nullify:
            level_below_which_to_nullify = min([self.ctx.levels_of_taxonomy.index(l) for l in taxonomic_levels_to_nullify])
            for taxonomic_level in self.ctx.levels_of_taxonomy[level_below_which_to_nullify:]:
                taxonomy_dict_entry[taxonomic_level] = None

        return taxonomy_dict_entry


    def get_blank_hit_template_dict(self):
        hit = {}

        for level in self.ctx.levels_of_taxonomy[::-1]:
            hit[level] = None

        return hit


    def get_consensus_taxonomy(self, items_taxonomy_dict):
        """Takes in a items_taxonomy_dict, returns a final taxonomic string that summarize all"""

        if not len(items_taxonomy_dict):
            return dict([(l, None) for l in self.ctx.levels_of_taxonomy])

        pd.set_option('mode.chained_assignment', None)

        item_hits = list([v for v in items_taxonomy_dict.values() if v['t_domain']])

        if not len(item_hits):
            return self.get_blank_hit_template_dict()

        df = pd.DataFrame.from_records(item_hits)

        # we have already stored a unique hash for taxonomy strings. here we will figure out most frequent
        # hash values in the df
        tax_hash_counts = df['tax_hash'].value_counts()
        tax_hash_df = tax_hash_counts.rename_axis('tax_hash').reset_index(name='frequency')
        max_frequency = tax_hash_df.frequency.max()
        tax_hash_df_most_frequent = tax_hash_df[tax_hash_df.frequency == max_frequency]

        if len(tax_hash_df_most_frequent.index) == 1:
            # if there is only a single winner, we're golden
            winner_tax_hash = tax_hash_df_most_frequent.tax_hash[0]

            # get the consensus hit based on the winner hash
            consensus_hit = df[df.tax_hash == winner_tax_hash].head(1)

            # turn it into a Python dict before returning
            return consensus_hit.to_dict('records')[0]
        else:
            # if there are competing hashes, we need to be more careful to decide
            # which taxonomic level should we use to cut things off.
            consensus_hit = self.get_blank_hit_template_dict()
            for level in self.ctx.levels_of_taxonomy[::-1]:
                if len(df[level].unique()) == 1:
                    consensus_hit[level] = df[level].unique()[0]

            return consensus_hit


    def print_taxonomy_hits_in_splits(self, hits, bin_name=None):
        self.progress.reset()
        self.run.warning(None, header='Hits for %s' % (bin_name if bin_name else "a bunch of splits"), lc="green")

        target = 'gene_name' if self.scgs_focus else 'anticodon'

        if len(hits) == 1 and hits[0][target] == 'CONSENSUS':
            self.run.info_single("No hits :/")
        else:
            if self.scgs_focus:
                header = [self._ITEM, 'gene', 'pct id', 'taxonomy']
                field_names = [self._VARIABLE_NAME_IN_TABLE, 'percent_identity', 'gene_callers_id']
            else:
                header = [self._ITEM, 'amino_acid', 'gene', 'pct id', 'taxonomy']
                field_names = [self._VARIABLE_NAME_IN_TABLE, 'amino_acid', 'percent_identity', 'gene_callers_id']

            table = []

            for hit in hits:
                taxon_text = ' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy])

                # if the hit we are working on sent here as 'consensus', we will color it up a bit so it shows up
                # more clearly in the debug output.
                if hit[self._VARIABLE_NAME_IN_TABLE] == 'CONSENSUS':
                    taxon_text = terminal.c(taxon_text, color='red')

                    for field_name in field_names:
                        if field_name in hit:
                            hit[field_name] = terminal.c(hit[field_name], color='red')
                        else:
                            hit[field_name] = terminal.c('--', color='red')

                if self.scgs_focus:
                    table.append([hit[self._VARIABLE_NAME_IN_TABLE], str(hit['gene_callers_id']), str(hit['percent_identity']), taxon_text])
                else:
                    table.append([hit[self._VARIABLE_NAME_IN_TABLE], str(hit['amino_acid']), str(hit['gene_callers_id']), str(hit['percent_identity']), taxon_text])

            anvio.TABULATE(table, header)


    def get_taxonomy_dict(self, gene_caller_ids, bin_name=None):
        items_taxonomy_dict = {}

        improper_gene_caller_ids = [g for g in gene_caller_ids if g not in self.gene_callers_id_to_item_taxonomy_dict]
        if improper_gene_caller_ids:
            raise ConfigError("Something weird is going on. Somehow anvi'o has a bunch of gene caller ids for which it is "
                              "supposed to estimate taxonomy. However, %d of them do not occur in a key dictionary. The code "
                              "here does not know what to suggest :( Apologies." % len(improper_gene_caller_ids))

        for gene_callers_id in gene_caller_ids:
            items_taxonomy_dict[gene_callers_id] = self.gene_callers_id_to_item_taxonomy_dict[gene_callers_id]
            items_taxonomy_dict[gene_callers_id]["tax_hash"] = HASH(self.gene_callers_id_to_item_taxonomy_dict[gene_callers_id])

        return items_taxonomy_dict


    def estimate_for_list_of_splits(self, split_names=None, bin_name=None):
        """Estimate {self._ITEM} taxonomy for a bunch of splits that belong to a single population.

           The purpose of this function is to to do critical things: identify genes we use for taxonomy in `split_names`,
           and generate a consensus taxonomy with the assumption that these are coming from splits that represents a
           single population.

           It will return a dictionary with multiple items, including a dictionary that contains the final consensus\
           taxonomy, another one that includes every {self._ITEM} and their raw associations with taxon names (from which the\
           consensus taxonomy was computed), as well as information about how many {self._ITEMS} were analyzed and supported the\
           consesnus.
        """

        if self.metagenome_mode:
            raise ConfigError("Someone is attempting to estimate taxonomy for a set of splits using a class inherited in "
                              "`metagenome mode`. If you are a programmer please note that it is best to use the member "
                              "function `estimate` directly.")

        consensus_taxonomy = None

        gene_caller_ids_of_interest = self.get_gene_caller_ids_for_splits(split_names)
        items_taxonomy_dict = self.get_taxonomy_dict(gene_caller_ids_of_interest)

        try:
            consensus_taxonomy = self.get_consensus_taxonomy(items_taxonomy_dict)
            consensus_taxonomy[self._VARIABLE_NAME_IN_TABLE] = 'CONSENSUS'
            consensus_taxonomy['percent_identity'] = '--'
            consensus_taxonomy['gene_callers_id'] = '--'

        except Exception as e:
            self.print_taxonomy_hits_in_splits(list(items_taxonomy_dict.values()))

            raise ConfigError(f"While trying to sort out the consensus taxonomy for %s anvi'o failed :( The list of {self._ITEM} taxon hits that "
                              f"caused the failure is printed in your terminal. But the actual error message that came from the depths "
                              f"of the codebase was this: '%s'." % (('the bin "%s"' % bin_name) if bin_name else 'a bunch of splits', e))

        if anvio.DEBUG:
            self.print_taxonomy_hits_in_splits(list(items_taxonomy_dict.values()) + [consensus_taxonomy], bin_name)

        # set some useful information. `total_items` is the number of items with taxonomy found in the collection of splits. the
        # `supporting_items` shows how many of them supports the consensus taxonomy fully
        total_items = len(items_taxonomy_dict)
        supporting_items = 0

        consensus_taxonomy_levels_occupied = [level for level in self.ctx.levels_of_taxonomy if consensus_taxonomy[level]]
        consensus_taxonomy_str = ' / '.join([consensus_taxonomy[level] for level in consensus_taxonomy_levels_occupied])

        for item_taxonomy_hit in items_taxonomy_dict.values():
            item_taxonomy_hit_str = ' / '.join([str(item_taxonomy_hit[level]) for level in consensus_taxonomy_levels_occupied])

            if item_taxonomy_hit_str == consensus_taxonomy_str:
                item_taxonomy_hit['supporting_consensus'] = True
                supporting_items += 1
            else:
                item_taxonomy_hit['supporting_consensus'] = False

        return {'consensus_taxonomy': consensus_taxonomy,
                self._ITEMS.lower(): items_taxonomy_dict,
                self._TOTAL_ITEMS: total_items,
                self._SUPPORTING_ITEMS: supporting_items,
                'metagenome_mode': False}


    def estimate_for_bins_in_collection(self):
        bins_taxonomy_dict = {}

        bin_name_to_split_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
        self.run.info_single("%s split names associated with %s bins of in collection '%s' have been "
                             "successfully recovered 🎊" % (pp(sum([len(v) for v in bin_name_to_split_names_dict.values()])),
                                                           pp(len(bin_name_to_split_names_dict)),
                                                           self.collection_name), nl_before=1)

        for bin_name in bin_name_to_split_names_dict:
            split_names = bin_name_to_split_names_dict[bin_name]
            bins_taxonomy_dict[bin_name] = self.estimate_for_list_of_splits(split_names, bin_name)

        return bins_taxonomy_dict


    def estimate_for_contigs_db_for_genome(self):
        contigs_db_taxonomy_dict = {}

        item_frequencies = self.frequency_of_items_with_taxonomy.values()
        if len([sf for sf in item_frequencies if sf > 1]) * 100 / len(item_frequencies) > 20:
            if self.scgs_focus:
                if self.just_do_it:
                    self.run.warning("Because you asked anvi'o to just do it, it will do it, but you seem to have too much contamination "
                                     "in this contigs database for it to represent a genome. So probably taxonomy estimations are all "
                                     "garbage, but hey, at least it runs?")
                else:
                    raise ConfigError("Because you haven't used the `--metagenome-mode` flag, anvi'o was trying to treat your contigs "
                                      "database as a genome. But there seems to be too much redundancy of single-copy core genes in this "
                                      "contigs database to assign taxonomy with any confidence :/ A more proper way to do this is to use the "
                                      "`--metagenome-mode` flag. Or you can also tell anvi'o to `--just-do-it`. It is your computer after "
                                      "all :( But you should still be aware that in that case you would likely get a completely irrelevant "
                                      "answer from this program.")
            elif self.trna_focus:
                reminder = ("Please note that since you haven't used the `--metagenome-mode` flag, anvi'o will treat this contigs "
                            "database as a genome and not a metagenome and will try to report a 'consensus' taxonomy based on all "
                            "the tRNA seqeunces and their taxonomic affiliations. ")

                if not self.per_item_output_file:
                    reminder += ("Since this will collapse all anticodons into a single result, you may want to use the flag "
                                 "`--per-anticodon-output-file` to your command line if you would like to see individual "
                                 "anticodons and their taxonomy for your genome.")

                self.run.warning(reminder)

        splits_in_contigs_database = self.split_name_to_gene_caller_ids_dict.keys()
        contigs_db_taxonomy_dict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(split_names=splits_in_contigs_database,
                                                                                                  bin_name=self.contigs_db_project_name)
        return contigs_db_taxonomy_dict


    def estimate_for_contigs_db_for_metagenome(self):
        """Treat a given contigs database as a metagenome.

           This function deserves some attention. It relies on a single SCG or anticodon to estimate the composition of a metagenome.
           For instance, its sister function, `estimate_for_contigs_db_for_genome`, works with a list of splits that are
           assumed to belong to the same genome. In which case a consensus taxonomy learned from all {self._ITEMS} is most
           appropriate. In this case, however, we don't know which split will go together, hence, we can't pull together
           {self._ITEMS} to learn a consensus taxonomy for independent populations in the metagenome. The best we can do is to stick
           with a single {self._ITEM} with the hope that (1) it will cut through as many populations as possible and (2) will have
           reasonable power to resolve taxonomy all by itself. These independent assumptions will both work in some cases
           and both fail in others.
        """

        # we first need to decide which {self._ITEM} we should use to survey taxonomy
        most_frequent_item = next(iter(self.frequency_of_items_with_taxonomy))
        if self.item_name_for_metagenome_mode:
            frequency_of_user_chosen_item = self.frequency_of_items_with_taxonomy[self.item_name_for_metagenome_mode]
            frequency_of_most_frequent_item = self.frequency_of_items_with_taxonomy[most_frequent_item]

            if frequency_of_user_chosen_item < frequency_of_most_frequent_item:
                additional_note = f" And just so you know, there is another {self._ITEM} that was observed more times (i.e., \
                                    {most_frequent_item}; {frequency_of_most_frequent_item} times) in this metagenome compared \
                                    to yours (i.e., {frequency_of_most_frequent_item} times). You're the boss, of course."
            else:
                additional_note = ""

            self.run.warning(f"As per your request anvi'o set '{self.item_name_for_metagenome_mode}' to be THE {self._SOURCE_DATA} \
                               to survey your metagenome for its taxonomic composition.{additional_note}")
        else:
            self.item_name_for_metagenome_mode = most_frequent_item

            self.run.warning(f"Anvi'o automatically set '{self.item_name_for_metagenome_mode}' to be THE {self._SOURCE_DATA} to "
                             f"survey your metagenome for its taxonomic composition. If you are not happy with that, you could "
                             f"change it with the parameter `{self._ITEM_FOR_METAGENOME_MODE_PARAM}`.")

        gene_caller_ids_of_interest = self.item_name_to_gene_caller_id_dict[self.item_name_for_metagenome_mode]
        items_taxonomy_dict = self.get_taxonomy_dict(gene_caller_ids=gene_caller_ids_of_interest,
                                                     bin_name=self.contigs_db_project_name)

        return {self.contigs_db_project_name: {self._ITEMS.lower(): items_taxonomy_dict,
                                               'metagenome_mode': True}}


    def get_items_taxonomy_super_dict(self):
        """Function that returns the `items_taxonomy_super_dict` for SCGs or anticodons.

           `items_taxonomy_super_dict` contains a wealth of information regarding samples, genes,
           gene taxonomic affiliations, consensus taxonomy, and coverages of genes across samples.
        """
        items_taxonomy_super_dict = {}

        if not self.initialized:
            self.init()

        if self.profile_db_path and not self.metagenome_mode:
            items_taxonomy_super_dict['taxonomy'] = self.estimate_for_bins_in_collection()
        elif not self.profile_db_path and not self.metagenome_mode:
            items_taxonomy_super_dict['taxonomy'] = self.estimate_for_contigs_db_for_genome()
        elif self.metagenome_mode:
            items_taxonomy_super_dict['taxonomy'] = self.estimate_for_contigs_db_for_metagenome()
        else:
            raise ConfigError("This class doesn't know how to deal with that yet :/")

        if self.compute_item_coverages and self.metagenome_mode:
            items_taxonomy_super_dict['coverages'] = self.get_item_coverages_across_samples_dict_in_metagenome_mode(items_taxonomy_super_dict)
        elif self.compute_item_coverages and not self.metagenome_mode:
            items_taxonomy_super_dict['coverages'] = self.get_item_coverages_across_samples_dict_in_genome_mode(items_taxonomy_super_dict)
        else:
            items_taxonomy_super_dict['coverages'] = None

        return items_taxonomy_super_dict


    def estimate(self):
        items_taxonomy_super_dict = self.get_items_taxonomy_super_dict()

        if self.update_profile_db_with_taxonomy:
            self.add_taxonomy_as_additional_layer_data(items_taxonomy_super_dict)

        self.print_items_taxonomy_super_dict(items_taxonomy_super_dict)

        if self.output_file_path:
            self.store_items_taxonomy_super_dict(items_taxonomy_super_dict)

        if self.per_item_output_file:
            self.store_taxonomy_per_item(items_taxonomy_super_dict)

        if self.sequences_file_path_prefix:
            self.store_sequences_for_items(items_taxonomy_super_dict)


    def print_items_taxonomy_super_dict(self, items_taxonomy_super_dict):
        if anvio.QUIET:
            return

        self.progress.reset()

        if self.collection_name:
            self.run.warning(None, header='Estimated taxonomy for collection "%s"' % self.collection_name, lc="green")
        elif self.metagenome_mode:
            self.run.warning(None, header='Taxa in metagenome "%s"' % self.contigs_db_project_name, lc="green")
        else:
            self.run.warning(None, header='Estimated taxonomy for "%s"' % self.contigs_db_project_name, lc="green")

        d = self.get_print_friendly_items_taxonomy_super_dict(items_taxonomy_super_dict)

        ordered_bin_names = sorted(list(d.keys()))

        if self.metagenome_mode:
            header = ['percent_identity', 'taxonomy']
        else:
            header = ['', self._TOTAL_ITEMS, self._SUPPORTING_ITEMS, 'taxonomy']

        # if we are in `--compute-xxx-coverages` mode, and more than 5 sample names, we are in trouble since they will\
        # unlikely fit into the display while printing them. so here we will cut it to make sure things look OK.
        samples_not_shown = 0
        sample_names_to_display = None
        if self.compute_item_coverages:
            sample_names_to_display = sorted(self.sample_names_in_profile_db)[0:5]
            samples_not_shown = sorted(self.sample_names_in_profile_db)[5:]

            header += sample_names_to_display

            if samples_not_shown:
                header += ['... %d more' % len(samples_not_shown)]

            # since we know coverages and sample names, we have a chance here to order the output
            # based on coverage. so let's do that.
            if self.metagenome_mode:
                sorted_bin_coverage_tuples = sorted([(bin_name, sum([d[bin_name]['coverages'][sample_name] for sample_name in self.sample_names_in_profile_db])) for bin_name in d], key=lambda x: x[1], reverse=True)
            else:
                sorted_bin_coverage_tuples = sorted([(bin_name, sum([(d[bin_name]['coverages'][sample_name] if d[bin_name][self._SUPPORTING_ITEMS] else 0) for sample_name in self.sample_names_in_profile_db])) for bin_name in d], key=lambda x: x[1], reverse=True)
            ordered_bin_names = [tpl[0] for tpl in sorted_bin_coverage_tuples]


        table = []
        for bin_name in ordered_bin_names:
            bin_data = d[bin_name]

            # set the taxonomy text depending on how much room we have. if there are sample coverages, keep it simple,
            # otherwise show the entire taxonomy text.
            if self.compute_item_coverages:
                taxon_text_l = ['(%s) %s' % (l.split('_')[1][0], bin_data[l]) for l in self.ctx.levels_of_taxonomy[::-1] if bin_data[l]]
                taxon_text = taxon_text_l[0] if taxon_text_l else '(NA) NA'
            else:
                taxon_text = ' / '.join([bin_data[l] if bin_data[l] else '' for l in self.ctx.levels_of_taxonomy])

            # setting up the table columns here.
            if self.metagenome_mode:
                row = [bin_name, str(bin_data['percent_identity']), taxon_text]
            else:
                row = [bin_name, str(bin_data[self._TOTAL_ITEMS]), str(bin_data[self._SUPPORTING_ITEMS]), taxon_text]

            # if there are coverages, add samples to the display too
            if self.compute_item_coverages:
                row += [d[bin_name]['coverages'][sample_name] for sample_name in sample_names_to_display]

            if samples_not_shown:
                row += ['... %d more' % len(samples_not_shown)]

            table.append(row)

        # if we are not in metagenome mode let's sort the output table based on total and
        # supporting items
        if not self.metagenome_mode:
            table = sorted(table, key=lambda x: (int(x[1]), int(x[2])), reverse=True)

        anvio.TABULATE(table, header)


    def store_items_taxonomy_super_dict(self, items_taxonomy_super_dict):
        d = self.get_print_friendly_items_taxonomy_super_dict(items_taxonomy_super_dict)

        if self.metagenome_mode:
            if self.scgs_focus:
                headers = ['scg_name', 'percent_identity']
            else:
                headers = ['anticodon', 'amino_acid', 'percent_identity']
        else:
            headers = ['bin_name', self._TOTAL_ITEMS, self._SUPPORTING_ITEMS]

        headers += self.ctx.levels_of_taxonomy

        if self.compute_item_coverages:
            headers_for_samples = sorted(self.sample_names_in_profile_db)
        else:
            headers_for_samples = []

        with open(self.output_file_path, 'w') as output:
            output.write('\t'.join(headers + headers_for_samples) + '\n')
            for item in d:
                line = [item] + [d[item][h] for h in headers[1:]]

                if self.compute_item_coverages:
                    for sample_name in headers_for_samples:
                        line.append(d[item]['coverages'][sample_name])

                output.write('\t'.join([str(f) for f in line]) + '\n')

        self.run.info("Output file", self.output_file_path, nl_before=1)


    def store_sequences_for_items(self, items_taxonomy_super_dict):
        """Report sequences for items if possible"""

        if self.ctx.focus != 'scgs':
            raise ConfigError("This function is only tested in SCGs mode. If you need to report "
                              "sequences for taxonomy items reported in other foci, please get in "
                              "touch with anvi'o developers.")

        if not self.scg_name_for_metagenome_mode:
            raise ConfigError("You can't ask anvi'o to store seqeunces for SCGs unless you are "
                              "working with a specific SCG name :(")

        d = self.get_print_friendly_items_taxonomy_super_dict(items_taxonomy_super_dict)

        c = ContigsSuperclass(self.args, r=run_quiet)

        gene_caller_ids = [v['gene_callers_id'] for v in d.values()]

        gene_caller_ids_list, sequences_dict = c.get_sequences_for_gene_callers_ids(gene_caller_ids, include_aa_sequences=True)
        if not len(gene_caller_ids_list):
            raise ConfigError("Something that should have never happened, happened :/ Please re-run the same command with "
                              "`--debug` and send the Traceback to an anvi'o developer.")

        dna_sequences_output_file_path = self.sequences_file_path_prefix + '_DNA.fa'
        amino_acid_sequences_output_file_path = self.sequences_file_path_prefix + '_AA.fa'

        contigs_db_name = c.a_meta['project_name']

        with open(amino_acid_sequences_output_file_path, 'w') as aa_sequences_output, open(dna_sequences_output_file_path, 'w') as dna_sequences_output:
            for entry in d.values():
                header = f"{contigs_db_name}_{entry['gene_name']}_{entry['gene_callers_id']}"
                dna_sequence = sequences_dict[entry['gene_callers_id']]['sequence']
                amino_acid_sequence = sequences_dict[entry['gene_callers_id']]['aa_sequence']

                aa_sequences_output.write(f">{header}\n{amino_acid_sequence}\n")
                dna_sequences_output.write(f">{header}\n{dna_sequence}\n")

        self.run.info("DNA sequences for SCGs", dna_sequences_output_file_path, nl_before=1)
        self.run.info("AA sequences for SCGs", amino_acid_sequences_output_file_path)


    def store_taxonomy_per_item(self, items_taxonomy_super_dict):
        if self.scgs_focus:
            headers = ['bin_name', 'gene_name', 'percent_identity']
        else:
            headers = ['bin_name','anticodon', 'amino_acid', 'percent_identity']

        headers += self.ctx.levels_of_taxonomy

        if self.compute_item_coverages:
            headers_for_samples = sorted(self.sample_names_in_profile_db)
        else:
            headers_for_samples = []

        with open(self.per_item_output_file, 'w') as output:
            output.write('\t'.join(headers + headers_for_samples) + '\n')
            for bin_name in items_taxonomy_super_dict["taxonomy"]:
                for item_entry in items_taxonomy_super_dict["taxonomy"][bin_name][self._ITEMS.lower()].values():
                    line = [bin_name] + [item_entry[h] for h in headers[1:]]

                    if self.compute_item_coverages:
                        for sample_name in headers_for_samples:
                            line.append(item_entry['coverages'][sample_name])

                    output.write('\t'.join([str(f) for f in line]) + '\n')

        self.run.info(f"A detailed {self._ITEM} output file", self.per_item_output_file)


    def get_print_friendly_items_taxonomy_super_dict(self, items_taxonomy_super_dict):
        contigs_db_name = anvio.dbops.ContigsDatabase(self.contigs_db_path).meta['project_name_str']

        d = {}

        if self.metagenome_mode:
            for item_hit in items_taxonomy_super_dict['taxonomy'][self.contigs_db_project_name][self._ITEMS.lower()].values():
                item_hit_name = '%s_%s_%d' % (contigs_db_name, item_hit[self._VARIABLE_NAME_IN_TABLE], item_hit['gene_callers_id'])

                d[item_hit_name] = item_hit

                if self.compute_item_coverages:
                    d[item_hit_name]['coverages'] = items_taxonomy_super_dict['coverages'][item_hit['gene_callers_id']]
        else:
            for bin_name in items_taxonomy_super_dict['taxonomy']:
                d[bin_name] = items_taxonomy_super_dict['taxonomy'][bin_name]['consensus_taxonomy']
                d[bin_name][self._TOTAL_ITEMS] = items_taxonomy_super_dict['taxonomy'][bin_name][self._TOTAL_ITEMS]
                d[bin_name][self._SUPPORTING_ITEMS] = items_taxonomy_super_dict['taxonomy'][bin_name][self._SUPPORTING_ITEMS]

                if self.compute_item_coverages:
                    d[bin_name]['coverages'] = items_taxonomy_super_dict['coverages'][bin_name]

        return d


    def add_taxonomy_as_additional_layer_data(self, items_taxonomy_super_dict):
        """A function that adds taxonomy to additional data tables of a given profile
           database. This will only work in metagenome mode."""

        if not self.metagenome_mode or not self.compute_item_coverages:
            return

        self.progress.new("Adding summary taxonomy for samples")
        self.progress.update('...')

        items_dict = list(items_taxonomy_super_dict['taxonomy'].values())[0][self._ITEMS.lower()]

        # at this stage each items_dict entry will look like this, and most critically will
        # have the same items
        #
        # "7660": {
        #       "gene_callers_id": 7660,
        #       "gene_name": "Ribosomal_S6",
        #       "accession": "CONSENSUS",
        #       "percent_identity": "98.9",
        #       "t_domain": "Bacteria",
        #       "t_phylum": "Firmicutes",
        #       "t_class": "Bacilli",
        #       "t_order": "Staphylococcales",
        #       "t_family": "Staphylococcaceae",
        #       "t_genus": "Staphylococcus",
        #       "t_species": null,
        #       "tax_hash": "b310c392"
        #     },
        #
        # this will enable us to learn the which item has been used to calculate coverage
        # information by only looking at a single entry:
        item_name = list(items_dict.values())[0][self._VARIABLE_NAME_IN_TABLE]

        # the might for loop to go through all taxonomic levels one by one
        for level in self.ctx.levels_of_taxonomy[::-1]:
            # setting the data group early on:
            data_group = '%s_%s' % (item_name, level[2:])
            self.progress.update('Working on %s-level data' % level)
            data_dict = {}
            data_keys_list = set([])
            for sample_name in self.sample_names_in_profile_db:
                data_dict[sample_name] = Counter()
                for gene_callers_id in items_dict:
                    # starting with a tiny hack to fill in missing values. here we first find
                    # the most highly resolved level of taxonomy that is not null for this
                    # particular item taxonomy
                    i = 0
                    for i in range(self.ctx.levels_of_taxonomy.index(level), 0, -1):
                        if items_dict[gene_callers_id][self.ctx.levels_of_taxonomy[i]]:
                            break

                    # just some abbreviations
                    l = self.ctx.levels_of_taxonomy[i][2:]
                    m = items_dict[gene_callers_id][self.ctx.levels_of_taxonomy[i]]

                    # if the best level we found in the previous step is matching to the level
                    # set by the main for loop, we're good to go with that name:
                    if level == self.ctx.levels_of_taxonomy[i]:
                        taxon_name = m
                    # otherwise we will try to replace that None name with something that is more
                    # sensible:
                    else:
                        taxon_name = "Unknown_%s_%s_%d" % (l, m, gene_callers_id)

                    # a key that will turn these data into stacked bar charts in the interface once
                    # they are added to the database:
                    key = '%s!%s' % (data_group, taxon_name)

                    # step where we add up all the values for each identical taxon names as we build
                    # the data dictionary:
                    data_dict[sample_name][key] += items_taxonomy_super_dict['coverages'][gene_callers_id][sample_name]
                    data_keys_list.add(key)

            # next few lines demonstrate the power of anvi'o quite nicely:
            self.progress.update("Updating additional data tables...")
            args = argparse.Namespace(profile_db=self.profile_db_path, target_data_group=data_group, just_do_it=True)
            T = TableForLayerAdditionalData(args, r=run_quiet, p=progress_quiet)
            T.add(data_dict, list(data_keys_list))

            self.progress.reset()

            self.run.info_single("%s level taxonomy is added to the profile database." % (level.capitalize()))

        self.progress.end()


    def get_gene_caller_ids_for_splits(self, split_names_list):
        """Returns gene caller ids found in a list of splits"""

        gene_caller_ids_for_splits = set([])
        for split_name in split_names_list:
            if split_name in self.split_name_to_gene_caller_ids_dict:
                gene_caller_ids_for_splits.update(self.split_name_to_gene_caller_ids_dict[split_name])

        return gene_caller_ids_for_splits


    def get_split_names_for_items_taxonomy_super_dict(self, items_taxonomy_super_dict):
        """Returns a list of split names associated with items found in a items_taxonomy_super_dict."""

        if self._ITEMS.lower() not in list(items_taxonomy_super_dict['taxonomy'].values())[0]:
            raise ConfigError("Someone called this function with something that doesn't look like the kind "
                              "of input data it was expecting (sorry for the vagueness of the message, but "
                              "anvi'o hopes that will be able to find out why it is happening).")

        split_names = set([])

        for entry_name in items_taxonomy_super_dict['taxonomy']:
            for gene_callers_id in items_taxonomy_super_dict['taxonomy'][entry_name][self._ITEMS.lower()]:
                split_names.add(self.gene_callers_id_to_split_name_dict[gene_callers_id])

        return split_names


    def get_item_coverages_across_samples_dict_in_genome_mode(self, items_taxonomy_super_dict):
        self.progress.reset()
        self.run.info_single(f"Anvi'o will now attempt to recover {self._ITEM} coverages in GENOME MODE from the profile "
                             "database, which contains %d samples." % (len(self.sample_names_in_profile_db)), nl_before=1, nl_after=1)

        item_coverages_across_samples_dict = self.get_item_coverages_across_samples_dict(items_taxonomy_super_dict)

        bin_avg_coverages_across_samples_dict = {}
        for bin_name in items_taxonomy_super_dict['taxonomy']:
            bin_avg_coverages_across_samples_dict[bin_name] = dict([(sample_name, None) for sample_name in self.sample_names_in_profile_db])
            for sample_name in self.sample_names_in_profile_db:
                average_coverage_across_samples = [item_coverages_across_samples_dict[gene_callers_id][sample_name] for gene_callers_id in items_taxonomy_super_dict['taxonomy'][bin_name][self._ITEMS.lower()]]
                if average_coverage_across_samples:
                    bin_avg_coverages_across_samples_dict[bin_name][sample_name] = np.mean(average_coverage_across_samples)

        self.run.warning(f"Anvi'o has just finished recovering {self._ITEM} coverages from the profile database to estimate "
                         f"the average coverage of your bins across your samples. Please note that anvi'o {self._DELIVERABLE} "
                         f"framework is using only {len(self.ctx.items)} {self._ITEMS} to estimate taxonomy. Which means, even a highly complete bin "
                         f"may be missing all of them. In which case, the coverage of that bin will be `None` across all "
                         f"your samples. The best way to prevent any misleading insights is take these results with a "
                         f"huge grain of salt, and use the `anvi-summarize` output for critical applications.",
                         header="FRIENDLY REMINDER", lc="blue")

        return bin_avg_coverages_across_samples_dict


    def get_item_coverages_across_samples_dict_in_metagenome_mode(self, items_taxonomy_super_dict):
        """Get item coverages in metagenome mode."""

        if not self.metagenome_mode:
            raise ConfigError("You're calling the wrong function. Your class is not in metagenome mode.")

        self.progress.reset()
        self.run.info_single(f"Anvi'o will now attempt to recover {self._ITEM} coverages from the profile database, which "
                             f"contains {len(self.sample_names_in_profile_db)} samples.", nl_before=1, nl_after=1)

        return self.get_item_coverages_across_samples_dict(items_taxonomy_super_dict)


    def get_item_coverages_across_samples_dict(self, items_taxonomy_super_dict):
        """Get item coverages"""
        item_coverages_across_samples_dict = {}

        self.progress.new('Recovering coverages')
        self.progress.update(f'Learning all split names affiliated with {self._ITEMS} ..')
        split_names_of_interest = self.get_split_names_for_items_taxonomy_super_dict(items_taxonomy_super_dict)
        self.progress.end()

        # initialize split coverages for splits that have anything to do with our items
        args = copy.deepcopy(self.args)
        args.split_names_of_interest = split_names_of_interest
        args.collection_name = None
        profile_db = ProfileSuperclass(args, p=self.progress, r=run_quiet)
        profile_db.init_split_coverage_values_per_nt_dict()

        # recover all gene caller ids that occur in our taxonomy estimation dictionary
        # and ge their coverage stats from the profile super
        gene_caller_ids_of_interest = set([])
        for bin_name in items_taxonomy_super_dict['taxonomy']:
            for gene_callers_id in items_taxonomy_super_dict['taxonomy'][bin_name][self._ITEMS.lower()]:
                gene_caller_ids_of_interest.add(gene_callers_id)

        # at this point we have everything. splits of interest are loaded in memory in `profile_db`, and we know
        # which gene caller ids we are interested in recovering coverages for. the way to access to gene coverages
        # is a bit convoluted in the dbops for historical reasons, but it is quite straightforward. the most
        # weird part is that we need a copy of a contigs super. so we will start with that:
        self.progress.new(f"Recovering {self._ITEM} coverages")
        self.progress.update("Initiating the contigs super class")
        contigs_db = ContigsSuperclass(self.args, r=run_quiet, p=progress_quiet)

        for split_name in split_names_of_interest:
            self.progress.update("Working with %s" % split_name)
            # note for the curious: yes, here we are sending the same gene caller ids of interest over and over to
            # the `get_gene_level_coverage_stats` for each split, but that function is smart enough to not spend any
            # time on those gene caller ids that do not occur in the split name we are interested in.
            all_item_stats_in_split, failed_gene_calls = profile_db.get_gene_level_coverage_stats(split_name, contigs_db, gene_caller_ids_of_interest=gene_caller_ids_of_interest)

            for item_stats in all_item_stats_in_split.values():
                for entry in item_stats.values():
                    gene_callers_id = int(entry['gene_callers_id'])
                    sample_name = entry['sample_name']
                    coverage = entry['non_outlier_mean_coverage']

                    if gene_callers_id not in item_coverages_across_samples_dict:
                        item_coverages_across_samples_dict[gene_callers_id] = dict([(sample_name, 0) for sample_name in self.sample_names_in_profile_db])

                    item_coverages_across_samples_dict[gene_callers_id][sample_name] = coverage

        self.progress.end()

        return item_coverages_across_samples_dict


class PopulateContigsDatabaseWithTaxonomy(TerminologyHelper):
    """Search sequences in corresponding databases, and store consensus taxonomy per seqeunce in the database.

       Depending on the context, this class currently deals with single-copy core genes or tRNA sequences. The
       key function here is the BLAST/DIAMOND search against the matching database, get back all the significant
       hits (or as many as `max_target_seqs` if hits exceed that), and try to find a LCA of all matches (hence it
       is important to not set `max_target_seqs` a small number).

       The class will recover sequences of interest from an anvi'o contigs database given the context. That said,
       it can also be used in an ad-hoc manner, where the user can search a single sequence.
    """
    def __init__(self, args):

        self.taxonomy_dict = OrderedDict()

        self.mutex = multiprocessing.Lock()

        TerminologyHelper.__init__(self)

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)
            contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
            self.contigs_db_project_name = contigs_db.meta['project_name']
            contigs_db.disconnect()
        else:
            self.contigs_db_project_name = "UNKNOWN"


    def get_sequences_dict_from_contigs_db(self):
        """Returns a dictionary of all HMM hits per SCG/anticodon of interest"""

        contigs_db = ContigsSuperclass(self.args, r=run_quiet, p=progress_quiet)
        splits_dict = {contigs_db.a_meta['project_name']: list(contigs_db.splits_basic_info.keys())}

        if self.scgs_focus:
            hmm_source = self.ctx.hmm_source_for_scg_taxonomy
            default_items_for_taxonomy = self.ctx.default_scgs_for_taxonomy
            return_amino_acid_sequences = True
        elif self.trna_focus:
            hmm_source = self.ctx.hmm_source_for_trna_genes
            default_items_for_taxonomy = self.ctx.default_anticodons_for_taxonomy
            return_amino_acid_sequences = False
        else:
            raise ConfigError("Anvi'o is lost.")

        s = hmmops.SequencesForHMMHits(self.args.contigs_db, sources=hmm_source, run=run_quiet, progress=progress_quiet)
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits(splits_dict, return_amino_acid_sequences=return_amino_acid_sequences)

        # a trick for trnas specifically since tRNA gene names contain two features and formed as `AMINOACID_ANTICODON`.
        # here we replace those gene names with `ANTICODON` only
        if self.trna_focus:
            for entry in hmm_sequences_dict:
                hmm_sequences_dict[entry]['gene_name'] = hmm_sequences_dict[entry]['gene_name'].split('_')[1]

        hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set(default_items_for_taxonomy))

        if not len(hmm_sequences_dict):
            return None

        self.progress.reset()
        self.run.info(f'Num relevant {self._ITEMS} in contigs db', '%s' % (pp(len(hmm_sequences_dict))))

        item_sequences_dict = {}
        for entry_id in hmm_sequences_dict:
            entry = hmm_sequences_dict[entry_id]

            item_name = entry['gene_name']
            if item_name in item_sequences_dict:
                item_sequences_dict[item_name][entry_id] = entry
            else:
                item_sequences_dict[item_name] = {entry_id: entry}

        return item_sequences_dict


    def populate_contigs_database(self):
        """Populates SCG/tRNA taxonomy tables in a contigs database"""

        # get an instnce for the tables for taxonomy early on.
        if self.scgs_focus:
            anvio_taxonomy_table_name = t.scg_taxonomy_table_name
            self.tables_for_taxonomy = TableForSCGTaxonomy(self.contigs_db_path, self.run, self.progress)
            database_version = self.ctx.scg_taxonomy_database_version
        elif self.trna_focus:
            anvio_taxonomy_table_name = t.trna_taxonomy_table_name
            self.tables_for_taxonomy = TableForTRNATaxonomy(self.contigs_db_path, self.run, self.progress)
            database_version = self.ctx.trna_taxonomy_database_version
        else:
            anvio_taxonomy_table_name = None
            self.tables_for_taxonomy = None
            database_version = None

        # get the dictionary that shows all hits for each self._ITEM of interest
        self.progress.new('Contigs bleep bloop')
        self.progress.update(f'Recovering the {self._ITEMS} dictionary')
        item_sequences_dict = self.get_sequences_dict_from_contigs_db()
        self.progress.end()

        if not item_sequences_dict:
            self.run.warning(f"This contigs database contains no {self._SOURCE_DATA} sequences that are used by the "
                             f"anvi'o taxonomy headquarters in Lausanne. Somewhat disappointing but totally OK.")

            # even if there are no SCGs to use for taxonomy later, we did attempt ot populate the
            # contigs database, so we shall note that in the self table to make sure the error from
            # `anvi-estimate-genome-taxonomy` is not "you seem to have not run taxonomy".
            self.tables_for_taxonomy.update_db_self_table_values(taxonomy_was_run=True, database_version=database_version)

            # return empty handed like a goose in the job market in 2020
            return None

        log_file_path = filesnpaths.get_temp_file_path()

        self.run.info('Taxonomy', self.ctx.accession_to_taxonomy_file_path)
        self.run.info('Database reference', self.ctx.search_databases_dir_path)
        self.run.info(f'Number of {self._ITEMS}', len(item_sequences_dict))

        self.run.warning('', header='Parameters for search', lc='green')
        self.run.info('Max number of target sequences', self.max_target_seqs)
        self.run.info('Max e-value to report alignments', self.evalue)
        self.run.info('Min percent identity to report alignments', self.min_pct_id)
        self.run.info('Num aligment tasks running in parallel', self.num_parallel_processes)
        self.run.info('Num CPUs per aligment task', self.num_threads)
        self.run.info('Log file path', log_file_path)

        self.tables_for_taxonomy.delete_contents_of_table(anvio_taxonomy_table_name, warning=False)
        self.tables_for_taxonomy.update_db_self_table_values(taxonomy_was_run=False, database_version=None)

        total_num_processes = len(item_sequences_dict)

        self.progress.new('Performing search', progress_total_items=total_num_processes)
        self.progress.update('Initializing %d process...' % int(self.num_parallel_processes))

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        error_queue = manager.Queue()

        search_output = []

        item_sequences_fate = {}

        for item_name in item_sequences_dict:
            item_sequences_fate[item_name] = {'in_contigs_db': len(item_sequences_dict[item_name]), 'in_search_results': 0}

            sequence = ""
            for entry in item_sequences_dict[item_name].values():
                if 'sequence' not in entry or 'gene_name' not in entry:
                    raise ConfigError("The `get_filtered_dict` function got a parameter that does not look like "
                                      "the way we expected it. This function expects a dictionary that contains "
                                      "keys `gene_name` and `sequence`.")

                sequence = sequence + ">" + str(entry['gene_callers_id']) + "\n" + entry['sequence'] + "\n"
                entry['hits'] = []

            input_queue.put([item_name, sequence])

        workers = []
        for i in range(0, int(self.num_parallel_processes)):
            worker = multiprocessing.Process(target=self.blast_search_worker, args=(input_queue, output_queue, error_queue, log_file_path))

            workers.append(worker)
            worker.start()

        num_finished_processes = 0
        while num_finished_processes < total_num_processes:
            # check error
            error_text = error_queue.get()
            if error_text:
                self.progress.reset()

                for worker in workers:
                    worker.terminate()

                if error_text and error_text is type(str) and 'incompatible' in error_text:
                    raise ConfigError(f"Your current databases are incompatible with the diamond version you have on your computer. "
                                      f"Please run the command `{self._SETUP_PROGRAM} --redo-databases` and come back.")
                else:
                    if self.scgs_focus:
                        raise ConfigError(f"Bad news. The database search operation failed somewhere :( It is very hard for anvi'o "
                                          f"to know what happened, but this is what we heard last: '{error_text}'. If this error message "
                                          f"makes sense ot you, great. Otherwise a LIKELY reason for this failure is that you have a diamond version "
                                          f"installed on your system that is incompatible with anvi'o :/ The best course of action for that "
                                          f"is to make sure running `diamond --version` on your terminal returns `0.9.14`. If not, "
                                          f"try to upgrade/downgrade your diamond to match this version. If you are in a conda environmnet "
                                          f"you can try running `conda install diamond=0.9.14`. Please feel free to contact us if the problem "
                                          f"persists. We apologize for the inconvenience.")
                    else:
                        raise ConfigError(f"Bad news. The database search operation failed somewhere :( It is very hard for anvi'o "
                                          f"to know what happened, but this is what we heard last: '{error_text}'. If this error message "
                                          f"makes sense ot you, great. Otherwise a LIKELY reason for this is that you have a BLAST version "
                                          f"installed on your system that is incompatible with anvi'o :/ If you see `2.10.1` or higher "
                                          f"version numbers when you type `blastn -version` in your terminal, it may be wortwhile to "
                                          f"get in touch with anvi'o developers. Note for developers: if the user has a newer version of "
                                          f"blast, this may be due to changing database extensions generated by `makeblastdb`.")

            try:
                search_output += output_queue.get()

                if self.write_buffer_size > 0 and len(search_output) % self.write_buffer_size == 0:
                    self.tables_for_taxonomy.add(search_output)
                    for s in search_output:
                        item_sequences_fate[s[1]]['in_search_results'] += 1

                    search_output = []

                num_finished_processes += 1

                self.progress.increment(increment_to=num_finished_processes)
                self.progress.update(f"%s of %s {self._ITEMS} are finished in %s processes with %s threads." \
                                        % (num_finished_processes, total_num_processes, int(self.num_parallel_processes), self.num_threads))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        # finally the remaining hits are written to the database, and we are done
        self.tables_for_taxonomy.add(search_output)
        for s in search_output:
            item_sequences_fate[s[1]]['in_search_results'] += 1

        # time to update the self table:
        self.tables_for_taxonomy.update_db_self_table_values(taxonomy_was_run=True, database_version=database_version)

        self.progress.end()

        # let user know about what happened to their sequences in the database
        header = [f"{self._ITEM}", "in_contigs_db", "in_search_results", "percent_annotation"]
        table = []

        for item_name in item_sequences_fate:
            table.append([item_name, item_sequences_fate[item_name]['in_contigs_db'], item_sequences_fate[item_name]['in_search_results'], round(item_sequences_fate[item_name]['in_search_results'] * 100 / item_sequences_fate[item_name]['in_contigs_db'], 1)])

        table.sort(key=lambda x: x[0])

        self.run.warning(f"Please take a careful look at the table below. It shows the number of sequences found in the contigs database "
                         f"for a given {self._ITEM}, and how many of them actually was annotated based on sequences in GTDB. The discrepancy "
                         f"between these two numbers can occur due to multiple reasons. You may have bacterial or archaeal sequences in your "
                         f"data that are way too novel to hit anything in the GTDB above the percent identity cutoff set by the parameter"
                         f"`--min-percent-identity` to report alignments (which is set to {self.min_pct_id} in the current run). Or, "
                         f"you may have a lot of eukaryotic organisms, {self._SOURCE_DATA} sequences of which may have no representation "
                         f"in the GTDB, in which case decreasing `--min-percent-identity` will only give you erronous results.")

        anvio.TABULATE(table, header, numalign="center")



    def store_hits_gene_callers_id(self, gene_callers_id, gene_sequence, item_name, hits):
        if self.trna_focus:
            header = ['db_name', 'gene_callers_id', 'anticodon', 'amino_acid', 'accession', 'percent_identity', 'bitscore'] + constants.levels_of_taxonomy
        elif self.scgs_focus:
            header = ['db_name', 'gene_callers_id', 'gene_name', 'accession', 'percent_identity', 'bitscore'] + constants.levels_of_taxonomy
        else:
            header = None

        if not os.path.exists(self.all_hits_output_file_path):
            with open(self.all_hits_output_file_path, 'w') as output:
                output.write('\t'.join(header) + '\n')

        if self.trna_focus:
            line = [str(f) for f in [self.contigs_db_project_name, gene_callers_id, item_name, constants.anticodon_to_AA[item_name]]]
        elif self.scgs_focus:
            line = [str(f) for f in [self.contigs_db_project_name, gene_callers_id, item_name]]

        with open(self.all_hits_output_file_path, 'a') as output:
            for hit in hits:
                if self.trna_focus:
                    output.write('\t'.join(line + [str(hit[k]) for k in header[4:]]) + '\n')
                elif self.scgs_focus:
                    output.write('\t'.join(line + [str(hit[k]) for k in header[3:]]) + '\n')


    def show_hits_gene_callers_id(self, gene_callers_id, gene_sequence, item_name, hits):
        self.progress.reset()
        self.run.warning(None, header=f"Hits for '{item_name}' w/gene caller id '{gene_callers_id}'", lc="green")
        self.run.info_single(gene_sequence, nl_after=1, level=0)

        if len(hits):
            header = ['%id', 'bitscore', 'accession', 'taxonomy']
            table = []

            for hit in hits:
                if hit['accession'] == 'CONSENSUS':
                    accession = terminal.c(hit['accession'], color='red')
                    taxon_text = terminal.c(' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy]), color='red')
                    table.append([str(hit['percent_identity']), str(hit['bitscore']), accession, taxon_text])
                else:
                    table.append([str(hit['percent_identity']), str(hit['bitscore']), hit['accession'], ' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy])])

            anvio.TABULATE(table, header)
        else:
            self.run.info_single(f"No hits to anything in {self.ctx.target_database_name} {self.ctx.target_database_release} at minimum percent identity of {self.min_pct_id} :/", nl_after=1, mc="red")


    def update_dict_with_taxonomy(self, d, mode=None):
        """Takes a dictionary that includes a key `accession` and populates the dictionary with taxonomy"""

        if not mode:
            if not 'accession' in d:
                raise ConfigError("`add_taxonomy_to_dict` is speaking: the dictionary sent here does not have a member "
                                  "with key `accession`.")

            if d['accession'] in self.ctx.accession_to_taxonomy_dict:
                d.update(self.ctx.accession_to_taxonomy_dict[d['accession']])
            else:
                d.update(self.ctx.accession_to_taxonomy_dict['unknown_accession'])

        elif mode == 'list_of_dicts':
            if len([entry for entry in d if 'accession' not in entry]):
                raise ConfigError("`add_taxonomy_to_dict` is speaking: you have a bad formatted data here :/")

            for entry in d:
                print(self.taxonomy_dict[entry['accession']])

        else:
            raise ConfigError("An unknown mode (%s) is set to `add_taxonomy_to_dict` :/" % (mode))

        return d


    def get_gene_estimation_output(self, item_name, fasta_formatted_sequence, log_file_path, show_all_hits=False):
        genes_estimation_output=[]

        if self.scgs_focus:
            target_database_path = self.ctx.SCGs[item_name]['db']

            diamond = Diamond(target_database_path, run=run_quiet, progress=progress_quiet)
            diamond.max_target_seqs = self.max_target_seqs
            diamond.evalue = self.evalue
            diamond.min_pct_id = self.min_pct_id
            diamond.num_threads = self.num_threads
            diamond.run.log_file_path = log_file_path

            raw_search_output = diamond.blastp_stdin_multi(fasta_formatted_sequence)
        elif self.trna_focus:
            target_database_path = self.ctx.anticodons[item_name]['db']

            blast = BLAST(None, target_database_path, run=run_quiet, progress=progress_quiet)
            blast.search_program = 'blastn'
            blast.max_target_seqs = self.max_target_seqs
            blast.evalue = self.evalue
            blast.min_pct_id = self.min_pct_id
            blast.num_threads = self.num_threads
            blast.run.log_file_path = log_file_path

            raw_search_output = blast.blast_stdin(fasta_formatted_sequence)
        else:
            raise ConfigError("You must be joking, Mr. Feynman.")

        hits_per_gene = {}
        for blastp_hit in raw_search_output.split('\n'):
            if len(blastp_hit) and not blastp_hit.startswith('Query'):
                fields = blastp_hit.split('\t')

                try:
                    gene_callers_id = int(fields[0])
                except:
                    gene_callers_id = fields[0]

                hit = dict(zip(['accession', 'percent_identity', 'bitscore'], [fields[1], float(fields[2]), float(fields[11])]))

                hit = self.update_dict_with_taxonomy(hit)

                if gene_callers_id not in hits_per_gene:
                    hits_per_gene[gene_callers_id] = {}

                if item_name not in hits_per_gene[gene_callers_id]:
                    hits_per_gene[gene_callers_id][item_name] = []

                hits_per_gene[gene_callers_id][item_name].append(hit)

        # logging and reporting block
        if (anvio.DEBUG or self.all_hits_output_file_path) and not len(hits_per_gene):
            with self.mutex:
                if anvio.DEBUG or show_all_hits:
                    self.progress.reset()
                    self.show_hits_gene_callers_id(fasta_formatted_sequence.split('\n')[0][1:], fasta_formatted_sequence.split('\n')[1], item_name, [])

                if self.all_hits_output_file_path:
                    self.store_hits_gene_callers_id(fasta_formatted_sequence.split('\n')[0][1:], fasta_formatted_sequence.split('\n')[1], item_name, [])

        for gene_callers_id, raw_hits in hits_per_gene.items():
            if len(raw_hits.keys()) > 1:
                self.run.warning("As crazy as it sounds, the gene callers id `%d` seems to have hit more than one SCG o_O Anvi'o will only use "
                                 "one of them almost absolutely randomly. Here are the SCGs the gene sequence matches: '%s'" % [s for s in raw_hits.keys()])

            item_name = list(raw_hits.keys())[0]
            raw_hits = raw_hits[item_name]

            consensus_hit = self.get_consensus_hit(raw_hits)
            consensus_hit['accession'] = 'CONSENSUS'

            if anvio.DEBUG or show_all_hits or self.all_hits_output_file_path:
                # avoid race conditions when priting this information when `--debug` is true:
                with self.mutex:
                    if anvio.DEBUG or show_all_hits:
                        self.progress.reset()
                        self.show_hits_gene_callers_id(gene_callers_id, fasta_formatted_sequence.split('\n')[1], item_name, raw_hits + [consensus_hit])

                    if self.all_hits_output_file_path:
                        self.store_hits_gene_callers_id(gene_callers_id, fasta_formatted_sequence.split('\n')[1], item_name, raw_hits + [consensus_hit])


            genes_estimation_output.append([gene_callers_id, item_name, [consensus_hit]])

        return genes_estimation_output


    def blast_search_worker(self, input_queue, output_queue, error_queue, log_file_path):
        """BLAST each SCG or tRNA sequence identified in the contigs database against the corresopinding
           target local database of GTDB seqeunces in parallel.
        """

        while True:
            item_name, fasta_formatted_sequence = input_queue.get(True)

            try:
                genes_estimation_output = self.get_gene_estimation_output(item_name, fasta_formatted_sequence, log_file_path)
                error_queue.put(None)
            except Exception as e:
                error_queue.put(e)

            output_queue.put(genes_estimation_output)


    def get_consensus_hit(self, raw_hits):
        pd.set_option('mode.chained_assignment', None)

        df = pd.DataFrame.from_records(raw_hits)

        # remove hits that are null at the phylum level if there are still hits
        # in the df that are not null:
        not_null_hits = df[df.t_phylum.notnull()]
        if len(not_null_hits):
            df = not_null_hits

        # find the max percent identity score in the df
        max_percent_identity = max(df['percent_identity'])

        # subset the data frame to those with percent identity that match to `max_percent_identity`
        df_max_identity = df.loc[df.percent_identity == max_percent_identity]

        # if some of the competing names have null species deignations, remove them from consideration
        if len(df_max_identity.t_species.unique()) > 1:
            df_max_identity = df_max_identity[df_max_identity.t_species.notnull()]

        # find the taxonomic level where the number of unique taxon names is one
        for taxonomic_level in self.ctx.levels_of_taxonomy[::-1]:
            if len(df_max_identity[taxonomic_level].unique()) == 1:
                break

        # take one of the hits from `df_max_identity`, and assign None to all taxonomic levels
        # beyond `taxonomic_level`, which, after the loop above shows the proper level of
        # assignment for this set
        final_hit = df_max_identity.head(1)
        for taxonomic_level_to_nullify in self.ctx.levels_of_taxonomy[self.ctx.levels_of_taxonomy.index(taxonomic_level) + 1:]:
            final_hit.at[0, taxonomic_level_to_nullify] = None

        # FIXME: final hit is still not what we can trust. next, we should find out whether the percent identity
        # for the level of taxonomy at `taxonomic_level` is higher than the minimum percent identity for all sequences
        # considered that are affiliated with final_hit[taxonomic_level]

        # turn it into a Python dict before returning
        final_hit_dict = final_hit.to_dict('records')[0]

        return final_hit_dict
