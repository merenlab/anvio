#!/usr/bin/env python3
# -*- coding: utf-8

import os
import re
import sys
import copy
import time
import shutil
import pickle
import tarfile
import argparse
import traceback
import multiprocessing
import pandas as pd

from copy import deepcopy
from tabulate import tabulate
from collections import Counter, OrderedDict

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as f
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.hmmopswrapper as hmmopswrapper

from anvio.dbops import ContigsSuperclass
from anvio.drivers import Aligners, driver_modules
from anvio.drivers.diamond import Diamond
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.tableops import Table
from anvio.tables.taxoestimation import TablesForTaxoestimation
from anvio.constants import levels_of_taxonomy

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Quentin Clayssen"
__email__ = "quentin.clayssen@gmail.com"
__status__ = "Development"

run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)
aligners = Aligners()


def timer(function):
    import time

    def timed_function(*args, **kwargs):
        n = 1
        start = time.time()
        for i in range(n):
            x = function(*args, **kwargs)
        end = time.time()
        print('Average time per call over {} calls for function \'{}\': {:6f} seconds'.format(
            n, function.__name__, (end - start) / n))
        return x
    return timed_function


class TaxonomyEstimation:
    def __init__(self, taxonomy_dict, dicolevel=None, reduction=False, run=run, progress=progress):
        if not dicolevel:
            self.pident_level_path = os.path.join(
                    os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY/GTDB/MIN_PCT_ID_PER_TAXONOMIC_LEVEL.pickle')

            with open(self.pident_level_path, 'rb') as handle:
                self.dicolevel = pickle.load(handle)

        if reduction:
            self.reduction=True
        else:
            self.reduction=False

        self.taxonomy_dict = taxonomy_dict

        self.taxonomic_levels_parser = dict([(l, l.split('_')[1][0] + "__") for l in levels_of_taxonomy])


    def show_table_score(self, name, selected_entrys_by_score):
        self.run.warning(None, header='%s' % (name), lc="yellow")
        header = ['Average bitscore', 'taxonomy']
        table = []

        for code, score in sorted(selected_entrys_by_score.items(), key=lambda x: (-x[1], x[0])):
            table.append([score, ' / '.join(self.taxonomy_dict[code].values())])

        print(tabulate(table, headers=header, tablefmt="fancy_grid", numalign="right"))


    def show_matrix_rank(self, name, matrix, list_position_entry, list_position_ribosomal):
        headers = list(matrix.keys())

        self.run.warning(None, header='%s' % (name), lc="blue")
        print(tabulate(matrix, headers=headers, tablefmt="fancy_grid", numalign="right"))


    def get_consensus_taxonomy(self, SCGs_hit_per_gene, name):
        consensus_taxonomy, entry = self.solo_hits(SCGs_hit_per_gene)
        SCG_hits_info = {}

        if consensus_taxonomy and entry:
            taxonomy=[{"bestSCG": name,
                       "bestident": entry['pident'],
                       "accession": entry['accession'],
                       "taxonomy": OrderedDict(consensus_taxonomy)}]
            return(consensus_taxonomy, taxonomy)
        else:
            consensus_taxonomy, taxonomy = self.rank_assignement(SCGs_hit_per_gene, name)

            return(consensus_taxonomy, taxonomy)


    def get_matching_gene(self, SCGs_hit_per_gene):
        matching_genes = [
            gene for gene in SCGs_hit_per_gene if len(SCGs_hit_per_gene[gene])]
        number_of_matching_genes = len(matching_genes)
        return matching_genes


    def fill_matrix(self, name, emptymatrix_ident, SCGs_hit_per_gene, list_position_entry, list_position_ribosomal, matchinggenes):

        for SCG in list_position_ribosomal:
            for entry in SCGs_hit_per_gene[SCG]:
                if entry['accession'] in list_position_entry:
                    emptymatrix_ident.at[entry['accession'], SCG] = float(
                        entry['pident'])

        return(emptymatrix_ident)


    def make_liste_individue(self, SCGs_hit_per_gene, matchinggenes):
        liste_individue = []
        liste_ribo = []
        for SCG in matchinggenes:
            if SCG not in liste_ribo:
                liste_ribo += [SCG]
            for entry in SCGs_hit_per_gene[SCG]:
                if entry['accession'] not in liste_individue:
                    liste_individue += [entry['accession']]
        return(liste_individue, liste_ribo)


    def make_rank_matrix(self, name, SCGs_hit_per_gene):
        matchinggenes = self.get_matching_gene(SCGs_hit_per_gene)

        list_position_entry, list_position_ribosomal = self.make_liste_individue(
            SCGs_hit_per_gene, matchinggenes)

        emptymatrix = pd.DataFrame(
            columns=list_position_ribosomal, index=list_position_entry)

        matrix_pident = self.fill_matrix(name, emptymatrix, SCGs_hit_per_gene, list_position_entry,
                                         list_position_ribosomal, matchinggenes)

        if anvio.DEBUG:
            self.show_matrix_rank(
                name, matrix_pident, list_position_entry, list_position_ribosomal)

        return(matrix_pident)


    def rank_assignement(self, SCGs_hit_per_gene, name, reduction=False):
        try:
            matrix_pident = self.make_rank_matrix(name, SCGs_hit_per_gene)

            taxonomy = self.make_list_taxonomy(matrix_pident)
            taxonomy_to_reduction = deepcopy(taxonomy)
            assignation_reduce = self.reduce_assignation_level(taxonomy_to_reduction)
            assignation = self.assign_taxonomie_solo_hit(assignation_reduce)
        except:
            traceback.print_exc()
            self.run.warning(SCGs_hit_per_gene, header='Fail matrix')
            assignation = []

        finally:
            return(assignation, taxonomy)


    def assign_taxonomie_solo_hit(self, taxonomy):
        if not taxonomy or not taxonomy[0]:
            return 0

        assignation = taxonomy[0]
        if taxonomy[1:]:
            for taxon in taxonomy[1:]:
                for level in taxon:
                    if taxon[level] not in list(assignation.values()) and taxon[level] != 'NA':
                        assignation[level] = 'NA'

        return(assignation)


    def make_list_taxonomy(self, matrix_pident):
        taxonomy = []
        matrix_rank = matrix_pident.rank(method='min', ascending=False, na_option='bottom')
        series_sum_rank = matrix_rank.sum(axis=1)
        minimum_rank = series_sum_rank.min()
        top_series = series_sum_rank.loc[series_sum_rank[:] == minimum_rank]

        for individue, row in top_series.items():
            bestSCG = pd.to_numeric(matrix_pident.loc[individue, :]).idxmax()
            bestident = matrix_pident.loc[individue, bestSCG]
            taxonomy.append({"bestSCG": bestSCG,
                             "bestident": bestident,
                             "accession": individue,
                             "taxonomy": OrderedDict(self.taxonomy_dict[individue])})

        return taxonomy


    def reduce_assignation_level(self, taxonomy):
        reduce_taxonomy = []

        for possibilitie in taxonomy:
            for level, value in reversed(possibilitie["taxonomy"].items()):
                if possibilitie["taxonomy"][level] == 'NA':
                    continue
                if possibilitie["bestident"] < float(self.dicolevel[possibilitie["bestSCG"]][self.taxonomic_levels_parser[level] + value]):
                    possibilitie["taxonomy"][level] = 'NA'
                else:
                    break

            reduce_taxonomy.append(possibilitie["taxonomy"])

        return reduce_taxonomy


    def solo_hits(self, SCGs_hit_per_gene, cutoff_solo_hit=0.90):
        consensus_taxonomy = []

        for SCG, entry in SCGs_hit_per_gene.items():
            if len(entry) == 1 and entry[0]['pident'] > cutoff_solo_hit:
                consensus_taxonomy.append(
                    self.taxonomy_dict[entry[0]['accession']])

        if consensus_taxonomy:
            assignation = self.assign_taxonomie_solo_hit(consensus_taxonomy)
            return assignation, entry[0]
        else:
            return False, False


class SCGsTaxonomy(TaxonomyEstimation):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.output_file_path = A('output_file')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.metagenome = True if A('metagenome') else False

        self.taxonomy_dict = {}
        hits_per_gene = {}

        self.initialized = False
        self.sanity_check()


    def sanity_check(self):
        if not self.contigs_db_path:
            raise ConfigError("This class needs an anvi'o contigs database to work with.")

        utils.is_contigs_db(self.contigs_db_path)

        if self.profile_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)


    def init(self, source="GTDB", number_scg=21):
        self.run.warning('', header='Taxonomy estimation for %s' % self.contigs_db_path, lc='green')
        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Source', source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('output file for taxonomy', self.output_file_path)

        self.tables_for_taxonomy = TablesForTaxoestimation(
            self.contigs_db_path, run, progress, self.profile_db_path)

        self.dictonnary_taxonomy_by_index = self.tables_for_taxonomy.get_data_for_taxonomy_estimation()

        if not (self.dictonnary_taxonomy_by_index):
            raise ConfigError("Anvi'o can't make a taxonomy estimation because aligment didn't return any\
                               match or you forgot to run 'anvi-diamond-for-taxonomy'.")


        if self.profile_db_path:
            collection_to_split = ccollections.GetSplitNamesInBins(self.args).get_dict()
            taxonomyestimation = TaxonomyEstimation.__init__(self, self.taxonomy_dict)

            self.initialized = True

            return(collection_to_split)
        else:
            self.initialized = True

            return


    def get_hits_per_bin(self,collection_to_split):

        self.tables_for_taxonomy = TablesForTaxoestimation(self.contigs_db_path, run, progress)
        self.dictonnary_taxonomy_by_index = self.tables_for_taxonomy.get_data_for_taxonomy_estimation()

        hits_per_gene={}
        split_to_gene_callers_id = dict()
        bin_to_gene_callers_id = dict()

        contigs_db = db.DB(self.contigs_db_path, anvio.__contigs__version__)

        for row in contigs_db.get_all_rows_from_table('genes_in_splits'):
            split_name, gene_callers_id = row[1], row[2]

            if split_name not in split_to_gene_callers_id:
                split_to_gene_callers_id[split_name] = set()

            split_to_gene_callers_id[split_name].add(gene_callers_id)

        for bin_name in collection_to_split:
            for split in collection_to_split[bin_name]:
                if bin_name not in bin_to_gene_callers_id:
                    bin_to_gene_callers_id[bin_name] = set()

                if split in split_to_gene_callers_id:
                    bin_to_gene_callers_id[bin_name].update(split_to_gene_callers_id[split])


        for gene_estimation in self.dictonnary_taxonomy_by_index.values():
            if gene_estimation["source"] == "Consensus":
                continue


            if gene_estimation['accession'] not in self.taxonomy_dict:
                self.taxonomy_dict[gene_estimation['accession']]= {"t_domain": gene_estimation['t_domain'],
                                                            "t_phylum": gene_estimation['t_phylum'],
                                                            "t_class": gene_estimation['t_class'],
                                                            "t_order": gene_estimation['t_order'],
                                                            "t_family": gene_estimation['t_family'],
                                                            "t_genus": gene_estimation['t_genus'],
                                                            "t_species": gene_estimation['t_species']}

            for bin_id, gene_callers_id in bin_to_gene_callers_id.items():
                if gene_estimation['gene_caller_id'] in gene_callers_id:
                    hit = [{'accession': gene_estimation['accession'], 'pident': float(gene_estimation['pourcentage_identity'])}]

                    if bin_id not in hits_per_gene:
                        hits_per_gene[bin_id] = {}

                    if gene_estimation['gene_name'] not in hits_per_gene[bin_id]:
                        hits_per_gene[bin_id][gene_estimation['gene_name']] = []

                    hits_per_gene[bin_id][gene_estimation['gene_name']] += hit

        return hits_per_gene


    def estimate_taxonomy(self):
        if self.metagenome:
            self.estimate_taxonomy_for_metagenome()

        if self.profile_db_path:
            entry_id = 0
            entries_db_profile = []
            dictionary_bin_taxonomy_estimation=dict()

            for bin_id, SCGs_hit_per_gene in hits_per_gene.items():
                consensus_taxonomy, taxonomy = self.get_consensus_taxonomy(SCGs_hit_per_gene, bin_id)
                dictionary_bin_taxonomy_estimation[bin_id]={"consensus_taxonomy": consensus_taxonomy,
                                                            "taxonomy_use_for_consensus": taxonomy}

                entries_db_profile += [(tuple([entry_id, self.collection_name, bin_id, source] + list(consensus_taxonomy.values())))]
                entry_id += 1

                self.tables_for_taxonomy.taxonomy_estimation_to_profile(
                    entries_db_profile)

            return dictionary_bin_taxonomy_estimation


    def estimate_taxonomy_for_metagenome(self, source="GTDB"):
        self.run.warning('', header='Taxonomy estimation for %s' %
        self.contigs_db_path, lc='green')
        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Source', source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('output file for taxonomy', self.output_file_path)
        self.tables_for_taxonomy = TablesForTaxoestimation(self.contigs_db_path, run, progress)
        self.dictonnary_taxonomy_by_index = self.tables_for_taxonomy.get_data_for_taxonomy_estimation()
        liste_SCG=[]
        output_genes_estimation = []
        estimate_taxonomy_presences=[]
        dictonarry_presence={}
        dictonnary_number_appear=dict()

        estimate_taxonomy_presences = [{"t_domain": "Unknow",
                                        "t_phylum": "NA",
                                        "t_class": "NA",
                                        "t_order": "NA",
                                        "t_family": "NA",
                                        "t_genus": "NA",
                                        "t_species": "NA"}]

        output_genes_estimation.append(['genes_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

        for gene_estimation in self.dictonnary_taxonomy_by_index.values():
            if gene_estimation["source"] == "GTDB" :
                continue

            if gene_estimation["gene_name"] not in liste_SCG:
                liste_SCG.append(gene_estimation["gene_name"])

            taxonomy = {"t_domain": gene_estimation['t_domain'],
                                                        "t_phylum": gene_estimation['t_phylum'],
                                                        "t_class": gene_estimation['t_class'],
                                                        "t_order": gene_estimation['t_order'],
                                                        "t_family": gene_estimation['t_family'],
                                                        "t_genus": gene_estimation['t_genus'],
                                                        "t_species": gene_estimation['t_species']}

            self.taxonomy_dict[gene_estimation['gene_caller_id']]=taxonomy

            output_genes_estimation.append([gene_estimation['gene_caller_id']] + list(taxonomy.values()))

            new=False
            for taxon_presente in taxonomy.values():
                if taxon_presente=="NA":
                    continue
                if taxon_presente not in dictonarry_presence and taxon_presente!="NA":
                    dictonarry_presence[taxon_presente]={gene_estimation['gene_name']: [gene_estimation['gene_caller_id']]}
                    new=True
                else:
                    if gene_estimation['gene_name'] not in dictonarry_presence[taxon_presente]:
                        dictonarry_presence[taxon_presente][gene_estimation['gene_name']]=[gene_estimation['gene_caller_id']]
                    else:
                        dictonarry_presence[taxon_presente][gene_estimation['gene_name']]+=[gene_estimation['gene_caller_id']]
            if new:
                for estimate_taxonomy_presence in estimate_taxonomy_presences:
                    share = {level: taxonomy[level] for level, taxon in estimate_taxonomy_presence.items() & taxonomy.items() if taxon!="Bacteria" and taxon!="Archea" and taxon!="NA"}
                    difference = { level : taxonomy[level] for level, taxon  in estimate_taxonomy_presence.items() - taxonomy.items() if taxon=="NA"}
                    if difference and share:
                        estimate_taxonomy_presence.clear()
                        estimate_taxonomy_presence.update(taxonomy)
                        new = False
                        break
                if new:
                    estimate_taxonomy_presences += [taxonomy]

        estimate_taxonomy_presences.pop(0)

        output=[['metagenome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'number of SCG\n(maximun %s)'% (len(liste_SCG))]]

        num_metagenome = 1

        for estimate_taxonomy_presence in estimate_taxonomy_presences:
            for level in estimate_taxonomy_presence.values():
                if level not in dictonnary_number_appear:
                    dictonnary_number_appear[level] = 1
                else:
                    dictonnary_number_appear[level] += 1
        if len(estimate_taxonomy_presences) > 1:
            for estimate_taxonomy_presence in estimate_taxonomy_presences:
                output+=[[self.contigs_db_path.replace(".db", "") + "_genome_" + str(num_metagenome)] + \
                         list(estimate_taxonomy_presence.values()) + \
                         [len(dictonarry_presence[list(estimate_taxonomy_presence.values())[-1]].values())]]

                num_metagenome += 1
        else:
            output += [[self.contigs_db_path.replace(".db", "")] + list(estimate_taxonomy_presence.values())]

        outpu_appear=[["taxon","number of scg"]]
        for level, appear in dictonnary_number_appear.items():
            outpu_appear += [[level], [appear]]

        for taxon, list_scgs in dictonarry_presence.items():
            if taxon=="NA":
                continue

            len_max_scg=("NA", 0)
            for SCG, list_appear_scg in list_scgs.items():
                if len(list_appear_scg) > len_max_scg[1]:
                    len_max_scg=(SCG,len(list_appear_scg))

            if dictonnary_number_appear[taxon] < len_max_scg[1]:
                # FIXME: what the fuck is going on here? What is this sentence supposed to mean?
                self.run.warning("%s is estimate %d time but it seams that the SCG %s have %d appear for this taxons, it could mean you have "\
                 % (taxon , dictonnary_number_appear[taxon], len_max_scg[0], len_max_scg[1]))

                continue

        self.show_taxonomy(output)


    def generate_output_file(self,output_data,header=False,append=False):
        if not self.output_file_path:
            return

        if filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=append):
            with open(self.output_file_path, "a") as output_file:
                output_data = ['\t'.join(line) for line in output_data]
                output_data='\n'.join(output_data)
                output_file.write(output_data)


    def show_taxonomy_estimation_bin(self):
        collection_to_split=self.init()

        possibles_taxonomy = []
        possibles_taxonomy.append(['Genome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        taxonomy_bin=self.assignation_by_bin(collection_to_split)

        hits_per_gene = self.get_hits_per_bin(collection_to_split)
        for bin_id, SCGs_hit_per_gene in hits_per_gene.items():
            consensus_taxonomy, taxonomy = self.get_consensus_taxonomy(SCGs_hit_per_gene, bin_id)
            possibles_taxonomy.append([bin_id] + list(consensus_taxonomy.values()))

        self.show_taxonomy(possibles_taxonomy)
        self.generate_output_file(possibles_taxonomy)


    def assignation_by_bin(self,collection_to_split):
        taxonomy_bin = {}

        hits_per_gene = self.get_hits_per_bin(collection_to_split)

        for bin_id, SCGs_hit_per_gene in hits_per_gene.items():
            consensus_taxonomy, taxonomy = self.get_consensus_taxonomy(SCGs_hit_per_gene, bin_id)
            taxonomy_bin[bin_id]={"taxonomy": consensus_taxonomy, "taxonomy_use": taxonomy}

        return(taxonomy_bin)


    def show_taxonomy_estimation_single_genome(self):
        self.init()
        self.SCGs_hit_per_gene={}
        for gene_estimation in self.dictonnary_taxonomy_by_index.values():
            self.taxonomy_dict[gene_estimation['accession']]= {"t_domain": gene_estimation['t_domain'],
                                                        "t_phylum": gene_estimation['t_phylum'],
                                                        "t_class": gene_estimation['t_class'],
                                                        "t_order": gene_estimation['t_order'],
                                                        "t_family": gene_estimation['t_family'],
                                                        "t_genus": gene_estimation['t_genus'],
                                                        "t_species": gene_estimation['t_species']}

            hit = [{'accession': gene_estimation['accession'], 'pident': float(gene_estimation['pourcentage_identity'])}]

            if gene_estimation['gene_name'] not in self.SCGs_hit_per_gene:
                self.SCGs_hit_per_gene[gene_estimation['gene_name']] = []

            self.SCGs_hit_per_gene[gene_estimation['gene_name']] += hit

        taxonomyestimation = TaxonomyEstimation.__init__(self, self.taxonomy_dict)

        consensus_taxonomy, taxonomy = TaxonomyEstimation.get_consensus_taxonomy(self, self.SCGs_hit_per_gene, self.contigs_db_path)

        output_full_genome=[['Genome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'], [self.contigs_db_path.replace(".db", "")] + list(consensus_taxonomy.values())]

        self.show_taxonomy(output_full_genome)
        self.generate_output_file(output_full_genome)


    def show_taxonomy(self,possibles_taxonomy):
        self.run.warning(None, header='Taxonomy estimation', lc="yellow")
        print(tabulate(possibles_taxonomy, headers="firstrow", tablefmt="fancy_grid", numalign="right"))


    def get_assignement_genes(self, source="GTDB"):
        self.run.warning('', header='Taxonomy estimation for %s' %
        self.contigs_db_path, lc='green')
        self.run.info('HMM PROFILE', "Bacteria 71")
        self.run.info('Source', source)
        self.run.info('Minimun level assigment', "species")
        self.run.info('output file for taxonomy', self.output_file_path)
        self.tables_for_taxonomy = TablesForTaxoestimation(self.contigs_db_path, run, progress)
        self.dictonnary_taxonomy_by_index = self.tables_for_taxonomy.get_data_for_taxonomy_estimation()
        possibles_taxonomy = []
        possibles_taxonomy.append(['gene_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

        for gene_estimation in self.dictonnary_taxonomy_by_index.values():
            if gene_estimation['source'] != "GTDB":
                taxonomy = [ gene_estimation['t_domain'],
                            gene_estimation['t_phylum'],
                            gene_estimation['t_class'],
                            gene_estimation['t_order'],
                            gene_estimation['t_family'],
                            gene_estimation['t_genus'],
                            gene_estimation['t_species']]

                possibles_taxonomy.append([str(gene_estimation['gene_caller_id'])] + taxonomy)

        self.show_taxonomy(possibles_taxonomy)
        self.generate_output_file(possibles_taxonomy)


    def show_taxonomy_estimation_genes(self, possibles_taxonomy):
        self.run.warning(None, header='Taxonomy estimation', lc="yellow")
        possibles_taxonomy_dataframe = pd.DataFrame(possibles_taxonomy, columns=possibles_taxonomy[0])
        possibles_taxonomy_dataframe.set_index("bin_id", inplace=True)
        possibles_taxonomy_dataframe = possibles_taxonomy_dataframe.sort_values(by=['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'], ascending=False)

        print(tabulate(possibles_taxonomy_dataframe, headers="firstrow",
                       tablefmt="fancy_grid", numalign="right"))
