#!/usr/bin/env python3
# -*- coding: utf-8
"""
Classes to setup remote SCG databases in local, use local databases to affiliate SCGs in anvi'o
contigs databases with taxon names, and estimate taxonomy for genomes and metagneomes.
"""

import os
import sys
import glob
import copy
import shutil
import hashlib
import argparse
import numpy as np
import pandas as pd
import multiprocessing
import scipy.sparse as sps

from tabulate import tabulate
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
from anvio.drivers.diamond import Diamond
from anvio.genomedescriptions import MetagenomeDescriptions
from anvio.tables.scgtaxonomy import TableForSCGTaxonomy
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.dbops import ContigsSuperclass, ContigsDatabase, ProfileSuperclass, ProfileDatabase


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print

HASH = lambda d: str(hashlib.sha224(''.join([str(d[level]) for level in constants.levels_of_taxonomy]).encode('utf-8')).hexdigest()[0:8])

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


# This is the most critical part of this entire operation. The following hard-coded dict translates
# between locally known 'HMM' names to FASTA files from GTDB. If one day you have to update this
# list, this is what you should do:
#
#   - find a FASTA file for a complete bacterial genome.
#   - generate an anvi'o contigs database, and run all default, installed SCG HMMs.
#   - export sequences for those HMMs that matches to the keys of the dictionary below (under all
#     circumstances these names must match to HMM sources in anvi'o Bacteria_71). you can do
#     something like this:
#
#             anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
#                                             -o Local_HMMs_export.fa \
#                                             --hmm-source Bacteria_71 \
#                                             --get-aa-sequences \
#                                             --return-best-hit \
#                                             --gene-names "Ribosomal_S2,Ribosomal_S3_C,Ribosomal_S6,Ribosomal_S7,Ribosomal_S8,Ribosomal_S9,Ribosomal_S11,Ribosomal_S20p,Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L6,Ribosomal_L9_C,Ribosomal_L13,Ribosomal_L16,Ribosomal_L17,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,ribosomal_L24,Ribosomal_L27A"
#             sed -i '' 's/___.*$//g' Local_HMMs_export.fa
#
#   - Then, BLAST sequences in Local_HMMs_export.fa to the entire collection of individual MSA FASTA
#     files from GTDB. For this, you could do something like this in msa_individual_genes directory
#     anvi'o generates, and carefully survey the OUTPUT.
#
#             for i in *faa; do makeblastdb -in $i -dbtype prot; done
#             for i in *faa; do echo; echo; echo $i; echo; echo; blastp -query Local_HMMs_export.fa -db $i -outfmt 6 -evalue 1e-10 -max_target_seqs 10; done > OUTPUT
#
#   - Update the list carefully based on the output.
#   - Find a FASTA file for a complete archaeal genome. Do the same :)
locally_known_HMMs_to_remote_FASTAs = {'Ribosomal_S2': ['ar122_TIGR01012.faa', 'bac120_TIGR01011.faa'],
                                       'Ribosomal_S3_C': ['ar122_TIGR01008.faa', 'bac120_TIGR01009.faa'],
                                       'Ribosomal_S6': ['bac120_TIGR00166.faa'],
                                       'Ribosomal_S7': ['ar122_TIGR01028.faa', 'bac120_TIGR01029.faa'],
                                       'Ribosomal_S8': ['ar122_PF00410.14.faa', 'bac120_PF00410.14.faa'],
                                       'Ribosomal_S9': ['ar122_TIGR03627.faa', 'bac120_PF00380.14.faa'],
                                       'Ribosomal_S11': ['ar122_TIGR03628.faa', 'bac120_TIGR03632.faa'],
                                       'Ribosomal_S20p': ['bac120_TIGR00029.faa'],
                                       'Ribosomal_L1': ['bac120_TIGR01169.faa', 'ar122_PF00687.16.faa'],
                                       'Ribosomal_L2': ['bac120_TIGR01171.faa'],
                                       'Ribosomal_L3': ['ar122_TIGR03626.faa', 'bac120_TIGR03625.faa'],
                                       'Ribosomal_L4': ['bac120_TIGR03953.faa'],
                                       'Ribosomal_L6': ['ar122_TIGR03653.faa', 'bac120_TIGR03654.faa'],
                                       'Ribosomal_L9_C': ['bac120_TIGR00158.faa'],
                                       'Ribosomal_L13': ['ar122_TIGR01077.faa', 'bac120_TIGR01066.faa'],
                                       'Ribosomal_L16': ['ar122_TIGR00279.faa', 'bac120_TIGR01164.faa'],
                                       'Ribosomal_L17': ['bac120_TIGR00059.faa'],
                                       'Ribosomal_L20': ['bac120_TIGR01032.faa'],
                                       'Ribosomal_L21p': ['bac120_TIGR00061.faa'],
                                       'Ribosomal_L22': ['ar122_TIGR01038.faa', 'bac120_TIGR01044.faa'],
                                       'ribosomal_L24': ['bac120_TIGR01079.faa', 'ar122_TIGR01080.faa'],
                                       'Ribosomal_L27A': ['bac120_TIGR01071.faa']
                                       }


class SCGTaxonomyContext(object):
    """The purpose of this base class is ot define file paths and constants for all single-copy
       core gene taxonomy operations.
    """
    def __init__(self, scgs_taxonomy_data_dir=None, scgs_taxonomy_remote_database_url=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        # hard-coded GTDB variables. poor design, but I don't think we are going do need an
        # alternative to GTDB.
        self.target_database = "GTDB"
        self.target_database_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
        self.target_database_files = ['VERSION', 'ar122_msa_individual_genes.tar.gz', 'ar122_taxonomy.tsv',
                                      'bac120_msa_individual_genes.tar.gz', 'bac120_taxonomy.tsv']

        # some variables from anvi'o constants
        self.hmm_source_for_scg_taxonomy = constants.default_hmm_source_for_scg_taxonomy
        self.default_scgs_taxonomy_data_dir = constants.default_scgs_taxonomy_data_dir
        self.default_scgs_for_taxonomy = constants.default_scgs_for_taxonomy
        self.levels_of_taxonomy = constants.levels_of_taxonomy

        # these are all the user accessible paths. defaults will serve well for all applications,
        # but these can be used for debugging.
        self.SCGs_taxonomy_data_dir = (os.path.abspath(scgs_taxonomy_data_dir) if scgs_taxonomy_data_dir else None) or (os.path.join(self.default_scgs_taxonomy_data_dir, self.target_database))
        self.msa_individual_genes_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'MSA_OF_INDIVIDUAL_SCGs')
        self.accession_to_taxonomy_file_path = os.path.join(self.SCGs_taxonomy_data_dir, 'ACCESSION_TO_TAXONOMY.txt')
        self.search_databases_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'SCG_SEARCH_DATABASES')
        self.target_database_URL = scgs_taxonomy_remote_database_url or self.target_database_URL

        # some dictionaries for convenience. we set them up here, but the proper place to sanity check
        # them may be somewhere else. for instance, when this class is inheritded by SetupLocalSCGTaxonomyData
        # the paths will not point to an actual file, but when it is inherited by PopulateContigsDatabaseWithSCGTaxonomy,
        # they better point to actual files.
        self.SCGs = dict([(SCG, {'db': os.path.join(self.search_databases_dir_path, SCG + '.dmnd'), 'fasta': os.path.join(self.search_databases_dir_path, SCG)}) for SCG in self.default_scgs_for_taxonomy])

        self.letter_to_level = dict([(l.split('_')[1][0], l) for l in self.levels_of_taxonomy])

        self.accession_to_taxonomy_dict = None
        if os.path.exists(self.accession_to_taxonomy_file_path):
            self.progress.new("Reading the accession to taxonomy file")
            self.progress.update('...')

            self.accession_to_taxonomy_dict = {}
            with open(self.accession_to_taxonomy_file_path, 'r') as taxonomy_file:
                for accession, taxonomy_text in [l.strip('\n').split('\t') for l in taxonomy_file.readlines() if not l.startswith('#') and l]:
                    # taxonomy_text kinda looks like these:
                    #
                    #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Alcaligenes;s__Alcaligenes faecalis_C
                    #    d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus_B;s__Enterococcus_B faecalis
                    #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter sp1
                    #    d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A cereus_AU
                    #    d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Tissierellales;f__Helcococcaceae;g__Finegoldia;s__Finegoldia magna_H
                    #
                    d = {}
                    for letter, taxon in [e.split('__', 1) for e in taxonomy_text.split(';')]:
                        if letter in self.letter_to_level:
                            # NOTE: This is VERY important. Here we are basically removing subclades GTDB defines for
                            # simplicity. We may have to change this behavior later. So basically, Enterococcus_B will
                            # become Enterococcus
                            if '_' in taxon:
                                if letter != 's':
                                    d[self.letter_to_level[letter]] = '_'.join(taxon.split('_')[:-1])
                                else:
                                    # special treatment for species level taxonomy string.
                                    # the genus is copied for the species level taxonomy, such as this one, 'Bacillus_A cereus', or
                                    # species itself may have a subclade, such as this one, 'Corynebacterium aurimucosum_C', so we
                                    # neeed to make sure the subclades are removed from all words in the species level
                                    # taxonomy string.
                                    d[self.letter_to_level[letter]] = ' '.join(['_'.join(word.split('_')[:-1]) if '_' in word else word for word in taxon.split(' ')])
                            else:
                                d[self.letter_to_level[letter]] = taxon
                        else:
                            self.run.warning("Some weird letter found in '%s' :(" % taxonomy_text)

                    self.accession_to_taxonomy_dict[accession] = d

            # let's add one more accession for all those missing accessions
            self.accession_to_taxonomy_dict['unknown_accession'] = dict([(taxon, None) for taxon in self.levels_of_taxonomy])

            self.progress.end()


# here we create an instance for the module. the idea is to overwrite it if
# it is necessary to overwrite some of the defaults
ctx = SCGTaxonomyContext()


class SanityCheck(object):
    def __init__(self):
        if self.skip_sanity_check:
            self.run.warning("We are skipping all sanity checks :( Dangerous stuff is happening.")
        else:
            self.sanity_check()


    def sanity_check(self):
        if sorted(list(locally_known_HMMs_to_remote_FASTAs.keys())) != sorted(self.ctx.default_scgs_for_taxonomy):
            raise ConfigError("Oh no. The SCGs designated to be used for all SCG taxonomy tasks in the constants.py "
                              "are not the same names described in locally known HMMs to remote FASTA files "
                              "conversion table definedd in SetupLocalSCGTaxonomyData module. If this makes zero "
                              "sense to you please ask a developer.")

        if not self.ctx.SCGs_taxonomy_data_dir:
            raise ConfigError("`SetupLocalSCGTaxonomyData` class is upset because it was inherited without "
                              "a directory for SCG taxonomy data to be stored :( This variable can't be None.")

        if self.user_taxonomic_level and self.user_taxonomic_level not in constants.levels_of_taxonomy:
            raise ConfigError("The taxonomic level %s is not a level anvi'o knows about. Here is the list of "
                              "taxonomic levels anvi'o recognizes: %s" % (', '.join(constants.levels_of_taxonomy)))

        # sanity checks specific to classes start below
        if self.__class__.__name__ in ['SetupLocalSCGTaxonomyData']:
            if self.reset and self.redo_databases:
                raise ConfigError("You can't ask anvi'o to both `--reset` and `--redo-databases` at the same time. Well. "
                                  "You can, but then this happens :/")

        if self.__class__.__name__ in ['SetupLocalSCGTaxonomyData', 'PopulateContigsDatabaseWithSCGTaxonomy']:
            if self.user_taxonomic_level:
                raise ConfigError("There is no need to set a taxonomic level while working with the class SetupLocalSCGTaxonomyData "
                                  "or PopulateContigsDatabaseWithSCGTaxonomy. Something fishy is going on :/")

        if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy', 'SCGTaxonomyEstimatorSingle', 'SCGTaxonomyEstimatorMulti']:
            if not os.path.exists(self.ctx.SCGs_taxonomy_data_dir):
                raise ConfigError("Anvi'o could not find the data directory for the single-copy core genes taxonomy "
                                  "setup. You may need to run `anvi-setup-scg-databases`, or provide a directory path "
                                  "where SCG databases are set up. This is the current path anvi'o is considering (which "
                                  "can be changed via the `--scgs-taxonomy-data-dir` parameter): '%s'" % (self.ctx.SCGs_taxonomy_data_dir))

            if not os.path.exists(self.ctx.accession_to_taxonomy_file_path):
                raise ConfigError("While your SCG taxonomy data dir seems to be in place, it is missing at least one critical "
                                  "file (in this case, the file to resolve accession IDs to taxon names). You may need to run "
                                  "the program `anvi-setup-scg-databases` with the `--reset` flag to set things right again.")

            ###########################################################
            # PopulateContigsDatabaseWithSCGTaxonomy
            ###########################################################
            if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy']:
                missing_SCG_databases = [SCG for SCG in self.ctx.SCGs if not os.path.exists(self.ctx.SCGs[SCG]['db'])]
                if len(missing_SCG_databases):
                    raise ConfigError("OK. It is very likley that if you run `anvi-setup-scg-databases` first you will be golden. "
                                      "Because even though anvi'o found the directory for taxonomy headquarters, "
                                      "your setup seems to be missing %d of %d databases required for everything to work "
                                      "with the current genes configuration of this class (sources say this is a record, FYI)." % \
                                                (len(missing_SCG_databases), len(self.ctx.SCGs)))

            ###########################################################
            # SCGTaxonomyEstimatorSingle
            #
            # Note: if something down below complains about a paramter
            #       because that actually belongs to the multi estimator
            #       class, you may need to set it to null in the class
            #       SCGTaxonomyArgs for single estimator
            #       initiation if clause
            ###########################################################
            if self.__class__.__name__ in ['SCGTaxonomyEstimatorSingle']:
                if self.metagenomes:
                    raise ConfigError("Taxonomy estimation classes have been initiated with a single contigs database, but your "
                            "arguments also include input for metagenomes. It is a no no. Please choose either. ")

                if self.output_file_prefix:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file path "
                                      "than an output file prefix.")

                if self.output_file_path:
                    filesnpaths.is_output_file_writable(self.output_file_path)

                if self.raw_output or self.matrix_format:
                    raise ConfigError("Haha in this mode you can't ask for the raw output or matrix format .. yet (we know that "
                                      "the parameter space of this program is like a mine field and we are very upset about it "
                                      "as well).")

                if not self.contigs_db_path:
                    raise ConfigError("For these things to work, you need to provide a contigs database for the anvi'o SCG "
                                      "taxonomy workflow :(")

                utils.is_contigs_db(self.contigs_db_path)

                scg_taxonomy_was_run = ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet).meta['scg_taxonomy_was_run']
                if not scg_taxonomy_was_run:
                    raise ConfigError("It seems the SCG taxonomy tables were not populated in this contigs database :/ Luckily it "
                                      "is easy to fix that. Please see the program `anvi-run-scg-taxonomy`.")

                if self.profile_db_path:
                    utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

                if self.collection_name and not self.profile_db_path:
                    raise ConfigError("If you are asking anvi'o to estimate taxonomy using a collection, you must also provide "
                                      "a profile database to this program.")

                if self.metagenome_mode and self.collection_name:
                    raise ConfigError("You can't ask anvi'o to treat your contigs database as a metagenome and also give it a "
                                      "collection.")

                if self.scg_name_for_metagenome_mode and not self.metagenome_mode:
                    raise ConfigError("If you are not running in `--metagenome-mode`, there is no use to define a SCG name for "
                                      "this mode :/")

                if self.scg_name_for_metagenome_mode and self.scg_name_for_metagenome_mode not in self.ctx.SCGs:
                    raise ConfigError("We understand that you wish to work with '%s' to study the taxonomic make up of your contigs "
                                      "database in metagenome mode. But then this gene is not one of those anvi'o recognizes as "
                                      "suitable SCGs to do that. Here is a list for you to choose from: '%s'." \
                                                            % (self.scg_name_for_metagenome_mode, ', '.join(self.ctx.SCGs.keys())))

                if self.compute_scg_coverages and not self.profile_db_path:
                    raise ConfigError("The flag `--compute-scg-coverages` is only good if there is a non-blank profile database around "
                                      "from which anvi'o can learn coverage statistics of genes across one or more samples :/")

                if self.profile_db_path and self.metagenome_mode and not self.compute_scg_coverages:
                    raise ConfigError("You have a profile database and you have asked anvi'o to estimate taxonomy in metagenome mode, "
                                      "but you are not asking anvi'o to compute SCG coverages which doesn't make much sense :/ Removing "
                                      "the profile database from this command or addint the flag `--compute-scg-coverages` would have "
                                      "made much more sense.")

                if self.profile_db_path and not self.metagenome_mode and not self.collection_name:
                    raise ConfigError("You have a profile database, and you are not in metagenome mode. In this case anvi'o will try to "
                                      "estimate coverages of SCGs in bins after estimating their taxonomy, but for that, you need to "
                                      "also provide a collection name. You can see what collections are available in your profile database "
                                      "you can use the program `anvi-show-collections-and-bins`, and then use the parameter "
                                      "`--collection-name` to tell anvi'o which one to use.")

                if self.update_profile_db_with_taxonomy:
                    if not self.metagenome_mode:
                        raise ConfigError("Updating the profile database with taxonomy layer data is only possible in metagenome "
                                          "mode :/ And not only that, you should also instruct anvi'o to compute single-copy core "
                                          "gene coverages.")

                    if not self.compute_scg_coverages:
                        raise ConfigError("You wish to update the profile database with taxonomy, but this will not work if anvi'o "
                                          "is computing coverages values of SCGs across samples (pro tip: you can ask anvi'o to do "
                                          "it by adding the flag `--compute-scg-coverages` to your command line).")

            ###########################################################
            # SCGTaxonomyEstimatorMulti
            ###########################################################
            if self.__class__.__name__ in ['SCGTaxonomyEstimatorMulti']:
                if self.args.contigs_db or self.args.profile_db:
                    raise ConfigError("Taxonomy estimation classes have been initiated with files for metagenomes, but your arguments "
                                      "include also a single contigs or profile database path. You make anvi'o nervous. "
                                      "Please run this program either with a metagenomes file or contigs/profile databases.")

                if self.output_file_path:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file prefix rather "
                                      "than an output file path. Anvi'o will use your prefix and will generate many files that start "
                                      "with that prefix but ends with different names for each taxonomic level.")

                if not self.output_file_prefix:
                    raise ConfigError("When using SCG taxonomy estimation in this mode, you must provide an output file prefix :/")

                if self.raw_output and self.matrix_format:
                    raise ConfigError("Please don't request anvi'o to report the output both in raw and matrix format. Anvi'o shall "
                                      "not be confused :(")

                if self.output_file_prefix:
                    filesnpaths.is_output_file_writable(self.output_file_prefix)


class SCGTaxonomyArgs(object):
    def __init__(self, args, format_args_for_single_estimator=False):
        """A base class to fill in common arguments for SCG Taxonomy classes.

        The purpose of this class is to reduce the complexity of setting member variables
        for various classes in this module so we can get away with a single multi-talented
        sanity check base class without any complaints regarding missing member variables.

        Parameters
        ==========
        format_args_for_single_estimator: bool
            This is a special case where an args instance is generated to be passed to the
            single estimator from within multi estimator. More specifically, the multi estimator
            class is nothing but one that iterates through all metagenomes
            given to it using the single estimator class. So it needs to create instances of
            single estimators, and collect results at upstream. The problem is, if a single
            estimtor is initiated with the args of a multi estimator, the sanity check will
            go haywire. This flag nullifies most common offenders.
        """

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A('output_file')
        self.output_file_prefix = A('output_file_prefix')
        self.just_do_it = A('just_do_it')
        self.simplify_taxonomy_information = A('simplify_taxonomy_information')
        self.metagenome_mode = True if A('metagenome_mode') else False
        self.scg_name_for_metagenome_mode = A('scg_name_for_metagenome_mode')
        self.compute_scg_coverages = A('compute_scg_coverages')
        self.report_scg_frequencies_path = A('report_scg_frequencies')
        self.metagenomes = A('metagenomes')
        self.user_taxonomic_level = A('taxonomic_level')
        self.matrix_format = A('matrix_format')
        self.raw_output = A('raw_output')

        if format_args_for_single_estimator:
            # so you're here to get an args instance to fool a single estimator class.
            # very cute. we shall make that happen.
            self.metagenomes = None
            self.output_file_path = None
            self.output_file_prefix = None
            self.matrix_format = None
            self.raw_output = None

        self.skip_sanity_check = A('skip_sanity_check')


class SCGTaxonomyEstimatorMulti(SCGTaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        """Iterate through metagenome descriptions using SCGTaxonomyEstimatorSingle"""

        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        SCGTaxonomyArgs.__init__(self, self.args)

        # set your context
        self.ctx = ctx

        # intiate sanity check
        SanityCheck.__init__(self)

        self.metagenomes = None
        self.profile_dbs_available = False


    def init_metagenomes(self):
        self.progress.new("Initializing contigs DBs")
        self.progress.update("...")
        g = MetagenomeDescriptions(self.args, run=run_quiet, progress=self.progress)
        g.load_metagenome_descriptions()

        # NOTE some enforced flags here.
        self.compute_scg_coverages = g.profile_dbs_available
        if not self.metagenome_mode and g.profile_dbs_available:
            self.metagenome_mode = True

        metagenomes_without_scg_taxonomy = [m for m in g.metagenomes if not g.metagenomes[m]['scg_taxonomy_was_run']]
        if metagenomes_without_scg_taxonomy:
            if len(metagenomes_without_scg_taxonomy) == len(g.metagenomes):
                self.progress.end()
                raise ConfigError("Surprise! None of the %d genomes had no SCG taxonomy information." % len(g.genomes))
            else:
                self.progress.end()
                raise ConfigError("%d of your %d genomes had no SCG taxonomy information. Here is the list: '%s'." % \
                        (len(metagenomes_without_scg_taxonomy), len(g.genomes), ', '.join(metagenomes_without_scg_taxonomy)))

        self.metagenomes = copy.deepcopy(g.metagenomes)
        self.metagenome_names = copy.deepcopy(g.metagenome_names)
        self.profile_dbs_available = g.profile_dbs_available

        self.progress.end()


    def estimate(self):
        if not self.metagenomes:
            self.init_metagenomes()
            self.run.info("Num metagenomes", len(self.metagenome_names))

        self.run.info("Taxonomic level of interest", self.user_taxonomic_level or "(None specified by the user, so 'all levels')")
        self.run.info("Output file prefix", self.output_file_prefix)
        self.run.info("Output in matrix format", self.matrix_format)
        self.run.info("Output raw data", self.raw_output)
        self.run.info("SCG coverages will be computed?", self.compute_scg_coverages)

        if self.report_scg_frequencies_path:
            self.report_scg_frequencies_as_TAB_delimited_file()
            return

        if not self.scg_name_for_metagenome_mode:
            self.scg_name_for_metagenome_mode = self.get_best_scg_name_for_metagenome_mode()

            self.run.warning("Please not that anvi'o just set the SCG for metagenome mode as '%s' since it was the most "
                             "frequent SCG occurring across all %d contigs databases involved in this analysis. But this "
                             "is nothing more than some heuristic for your convenience, and we strongly advice you to "
                             "run this program with the parameter `--report-scg-frequencies` and examine the output "
                             "to see if there is a better choice. As you can imagine, the most frequent SCG may not be "
                             "the one that is more common across all genomes or metagenomes you are interested." % \
                                                (self.scg_name_for_metagenome_mode, len(self.metagenomes)))

            self.run.info("SCG [determined by anvi'o]", self.scg_name_for_metagenome_mode, nl_after=1, mc="green")
        else:
            self.run.info("SCG [chosen by the user]", self.scg_name_for_metagenome_mode, nl_after=1, mc="green")

        scg_taxonomy_super_dict_multi = self.get_scg_taxonomy_super_dict_for_metagenomes()

        self.store_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)


    def get_print_friendly_scg_taxonomy_super_dict_multi(self, scg_taxonomy_super_dict_multi, as_data_frame=False):
        """Extract a more print-friendly data structure from `scg_taxonomy_super_dict_multi`

        Returns
        =======
        d: dict
            This is a dictionary that summarizes taxonomy and coverages per unique SCG in
            a given metagenome and looks like this:

            ----8<----8<----8<----8<----8<----8<----8<----8<----
            {
              "USA0001": {
                "Ribosomal_L16_69898": {
                  "gene_callers_id": 69898,
                  "gene_name": "Ribosomal_L16",
                  "accession": "CONSENSUS",
                  "percent_identity": "100.0",
                  "t_domain": "Bacteria",
                  "t_phylum": "Bacteroidota",
                  "t_class": "Bacteroidia",
                  "t_order": "Bacteroidales",
                  "t_family": "Rikenellaceae",
                  "t_genus": "Alistipes",
                  "t_species": "Alistipes shahii",
                  "tax_hash": "7e3cc7fc",
                  "coverages": {
                    "USA0001_01": 59.13119533527697
                  }
                },
                "Ribosomal_L16_55413": {
                  "gene_callers_id": 55413,
                  "gene_name": "Ribosomal_L16",
                  "accession": "CONSENSUS",
                  "percent_identity": "100.0",
                  "t_domain": "Bacteria",
                  "t_phylum": "Bacteroidota",
                  "t_class": "Bacteroidia",
                  "t_order": "Bacteroidales",
                  "t_family": "Tannerellaceae",
                  "t_genus": "Parabacteroides",
                  "t_species": null,
                  "tax_hash": "d16f9017",
                  "coverages": {
                    "USA0001_01": 48.82825484764543
                  }
                }

                (...)

            }
            ----8<----8<----8<----8<----8<----8<----8<----8<----

        """

        scg_taxonomy_super_dict_multi_print_friendly = {}

        for metagenome_name in scg_taxonomy_super_dict_multi:
            args = SCGTaxonomyArgs(self.args, format_args_for_single_estimator=True)
            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            if self.metagenome_mode:
                args.metagenome_mode = True
            else:
                args.metagenome_mode = False

            if self.profile_dbs_available:
                args.profile_db = self.metagenomes[metagenome_name]['profile_db_path']
                args.compute_scg_coverages = True

            d = SCGTaxonomyEstimatorSingle(args, run=run_quiet).get_print_friendly_scg_taxonomy_super_dict(scg_taxonomy_super_dict_multi[metagenome_name])

            # NOTE: what is happening down below might look stupid, because it really is. items of `d` here
            # contain a dictionary with the key 'coverages' where the coverage value of a given gene or bin
            # is kept per sample. to simplify downstream data wrangling operations, here we replace the
            # coverages dictionary with a single 'coverage': value entry to `d`. The
            # only reason we can do this here is becasue we only allow the user to use single profiles
            # with their metagenome descriptions and we are certain that there will be only a single entry in the
            # coverages dictionary.
            if self.compute_scg_coverages:
                for item in d:
                    coverages_list = list(d[item]['coverages'].items())

                    if len(coverages_list) > 1:
                        raise ConfigError("The codebase is not ready to handle this :(")

                    d[item]['coverage'] = coverages_list[0][1]
                    d[item].pop('coverages')

            scg_taxonomy_super_dict_multi_print_friendly[metagenome_name] = d

        if as_data_frame:
            return self.print_friendly_scg_taxonomy_super_dict_multi_to_data_frame(scg_taxonomy_super_dict_multi_print_friendly)
        else:
            return scg_taxonomy_super_dict_multi_print_friendly


    def print_friendly_scg_taxonomy_super_dict_multi_to_data_frame(self, scg_taxonomy_super_dict_multi_print_friendly):
        """Take a `scg_taxonomy_super_dict_multi_print_friendly`, and turn it into a neat data frame.

        Not working with dataframes from the get go sounds like a silly thing to do, but Python
        objects have been the way we passed around all the data (including to the interactive
        interface) so it feels weird to change that completely now. But obviously there are many
        advantages of working with dataframes, so the purpose of this function is to benefit from
        those in the context of multi estimator class.

        If the global `anvio.DEBUG` is True (via `--debug`), the dataframe is stored as a pickle
        object to play.

        Parameters
        ==========
        scg_taxonomy_super_dict_multi_print_friendly: dict
             This data structure is defined in `get_print_friendly_scg_taxonomy_super_dict_multi`

        Returns
        =======
        DF: pandas dataframe
            A neatly sorted dataframe. Looks like this:

            ----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<------
                metagenome_name                         taxon     coverage taxonomic_level
            0           USA0001                      Bacteria  1077.277012        t_domain
            1           USA0001               Unknown_domains     0.000000        t_domain
            2           USA0003                      Bacteria  2668.593534        t_domain
            3           USA0003               Unknown_domains     9.177677        t_domain
            4           USA0001                  Bacteroidota   932.478308        t_phylum
            5           USA0001                    Firmicutes   126.837011        t_phylum
            6           USA0001                Proteobacteria    17.961694        t_phylum
            7           USA0001                 Unknown_phyla     0.000000        t_phylum
            8           USA0001             Verrucomicrobiota     0.000000        t_phylum
            9           USA0003                    Firmicutes  1602.189160        t_phylum
            10          USA0003                  Bacteroidota  1034.471363        t_phylum
            11          USA0003             Verrucomicrobiota    22.668447        t_phylum
            12          USA0003                Proteobacteria     9.264563        t_phylum
            13          USA0003                 Unknown_phyla     9.177677        t_phylum
            14          USA0001                   Bacteroidia   932.478308         t_class
            15          USA0001                    Clostridia   116.441241         t_class
            16          USA0001           Gammaproteobacteria    17.961694         t_class
            17          USA0001                 Negativicutes    10.395770         t_class
            18          USA0001                       Bacilli     0.000000         t_class
            19          USA0001                 Lentisphaeria     0.000000         t_class
            ..              ...                           ...          ...             ...
            296         USA0003           CAG-302 sp001916775    11.435762       t_species
            297         USA0003   Paramuribaculum sp001689565    11.149206       t_species
            298         USA0003          CAG-1435 sp000433775    11.078189       t_species
            299         USA0003           CAG-269 sp001916005     9.848404       t_species
            300         USA0003            UBA737 sp002431945     7.497512       t_species
            301         USA0003           CAG-345 sp000433315     6.733333       t_species
            302         USA0003           UBA1829 sp002405835     6.684840       t_species
            303         USA0003           TF01-11 sp001916135     6.060109       t_species
            304         USA0003       Barnesiella sp002161555     0.000000       t_species
            305         USA0003           CAG-115 sp000432175     0.000000       t_species
            306         USA0003           CAG-452 sp000434035     0.000000       t_species
            307         USA0003       Duncaniella sp002494015     0.000000       t_species
            308         USA0003   Paramuribaculum sp001689535     0.000000       t_species
            309         USA0003  Succiniclasticum sp002342505     0.000000       t_species
            ----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<----8<------

        """

        self.progress.new("Data dict to dataframe")
        self.progress.update("Bleep very complex stuff bloop")

        taxonomic_levels = [self.user_taxonomic_level] if self.user_taxonomic_level else self.ctx.levels_of_taxonomy

        df = pd.DataFrame.from_dict({(i,j): scg_taxonomy_super_dict_multi_print_friendly[i][j]
                                for i in scg_taxonomy_super_dict_multi_print_friendly.keys()
                                for j in scg_taxonomy_super_dict_multi_print_friendly[i].keys()}, orient='index')

        df.reset_index(inplace=True)
        df.rename(columns={"level_0": "metagenome_name"}, inplace=True)
        df.drop(['level_1'], axis=1, inplace=True)

        if not self.compute_scg_coverages:
            # NOTE: this is a bit critical. if the user are working with external genomes where there is no
            # coverage information is availble for SCGs, the following code will explode for obvious
            # reasons (which include the fact that `scg_taxonomy_super_dict_multi` will not have an entry
            # for `coverage`). What we will do here is a little trick. First, we will add a coverage column
            # to `df`. And at the very end, we will replace it with `times_observed`.
            df["coverage"] = 1

        DFx = pd.DataFrame(columns=['metagenome_name', 'taxon', 'coverage', 'taxonomic_level'])
        for taxonomic_level in taxonomic_levels:
            x = df.fillna(constants.levels_of_taxonomy_unknown).groupby(['metagenome_name', taxonomic_level])['coverage'].sum().reset_index()
            x.columns = ['metagenome_name', 'taxon', 'coverage']
            x['taxonomic_level'] = taxonomic_level

            # adding for each metagenome the missing taxon names
            list_of_unique_taxon_names = x['taxon'].unique()
            for metagenome_name in x['metagenome_name'].unique():
                known_taxa = set(x[x['metagenome_name'] == metagenome_name]['taxon'])
                missing_taxa = [taxon for taxon in list_of_unique_taxon_names if taxon not in known_taxa]
                for taxon in missing_taxa:
                    row = {'metagenome_name': metagenome_name, 'taxon': taxon, 'coverage': 0.0, 'taxonomic_level': taxonomic_level}
                    x = x.append(row, ignore_index=True)

            DFx = DFx.append(x, ignore_index=True)

        DFx.sort_values(by=['metagenome_name', 'taxonomic_level', 'coverage'], ascending=[True, True, False], inplace=True)

        DF = pd.DataFrame(columns=['metagenome_name', 'taxon', 'coverage', 'taxonomic_level'])
        for taxonomic_level in constants.levels_of_taxonomy:
            DF = DF.append(DFx[DFx['taxonomic_level'] == taxonomic_level], ignore_index=True)

        del DFx

        if anvio.DEBUG:
            import pickle
            with open('DataFrame.pickle', 'wb') as output:
                pickle.dump(DF, output)

            self.progress.reset()
            self.run.info_single("The `--debug` flag made anvi'o to dump the data frame that was generated by the "
                                 "function `print_friendly_scg_taxonomy_super_dict_multi_to_data_frame` (lol) as "
                                 "a picke object stored in file name 'DataFrame.pickle'. If you would like to play "
                                 "with it, you can start an ipython session, and run the following: import pandas "
                                 "as pd; import pickle; df = pickle.load(open('DataFrame.pickle', 'rb'))",
                                 nl_before=1, nl_after=1)

        if not self.compute_scg_coverages:
            DF.rename(columns={"coverage": "times_observed"}, inplace = True)
            DF['times_observed'] = DF['times_observed'].astype(int)

        self.progress.end()

        return DF


    def store_scg_taxonomy_super_dict_multi(self, scg_taxonomy_super_dict_multi):
        if self.raw_output:
            self.store_scg_taxonomy_super_dict_raw(scg_taxonomy_super_dict_multi)
        else:
            if self.matrix_format:
                self.store_scg_taxonomy_super_dict_multi_matrix_format(scg_taxonomy_super_dict_multi)
            else:
                self.store_scg_taxonomy_super_dict_multi_long_format(scg_taxonomy_super_dict_multi)


    def store_scg_taxonomy_super_dict_multi_long_format(self, scg_taxonomy_super_dict_multi):
        df = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi, as_data_frame=True)

        output_file_path = self.output_file_prefix + '-LONG-FORMAT.txt'
        df.to_csv(output_file_path, index=True, index_label="entry_id", sep='\t')

        self.run.info("Long-format output", output_file_path)


    def store_scg_taxonomy_super_dict_raw(self, scg_taxonomy_super_dict_multi):
        d = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi)

        taxonomic_levels = [self.user_taxonomic_level] if self.user_taxonomic_level else self.ctx.levels_of_taxonomy

        header = ['metagenome_name', 'gene_name', 'gene_callers_id', 'percent_identity']

        if self.compute_scg_coverages:
            header += ['coverage']

        header += taxonomic_levels

        output_file_path = self.output_file_prefix + '-RAW-LONG-FORMAT.txt'
        with open(output_file_path, 'w') as output:
            output.write('\t'.join(header) + '\n')

            for metagenome_name in d:
                for gene_name in d[metagenome_name]:
                    output.write('\t'.join([metagenome_name] + [str(d[metagenome_name][gene_name][h]) for h in header[1:]]) + '\n')

        self.run.info("Raw output", output_file_path)


    def store_scg_taxonomy_super_dict_multi_matrix_format(self, scg_taxonomy_super_dict_multi):
        df = self.get_print_friendly_scg_taxonomy_super_dict_multi(scg_taxonomy_super_dict_multi, as_data_frame=True)

        taxonomic_levels = [self.user_taxonomic_level] if self.user_taxonomic_level else self.ctx.levels_of_taxonomy

        for taxonomic_level in taxonomic_levels:
            dfx = df[df['taxonomic_level'] == taxonomic_level]
            dfx.set_index(['metagenome_name', 'taxon'], inplace=True)

            if self.compute_scg_coverages:
                matrix = sps.coo_matrix((dfx.coverage, (dfx.index.codes[0], dfx.index.codes[1]))).todense().tolist()
            else:
                matrix = sps.coo_matrix((dfx.times_observed, (dfx.index.codes[0], dfx.index.codes[1]))).todense().tolist()

            rows = dfx.index.levels[0].tolist()
            cols = ['taxon'] + dfx.index.levels[1].tolist()

            output_file_path = '%s-%s-MATRIX.txt' % (self.output_file_prefix, taxonomic_level)
            temp_file_path = filesnpaths.get_temp_file_path()
            with open(temp_file_path, 'w') as output:
                output.write('\t'.join(cols) + '\n')
                for i in range(0, len(matrix)):
                    output.write('\t'.join([rows[i]] + ['%.2f' % c for c in matrix[i]]) + '\n')

            utils.transpose_tab_delimited_file(temp_file_path, output_file_path, remove_after=True)

            self.run.info('Output matrix for "%s"' % taxonomic_level, output_file_path)


    def report_scg_frequencies_as_TAB_delimited_file(self):
        scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies = self.get_scg_frequencies()

        utils.store_dict_as_TAB_delimited_file(scg_frequencies, self.report_scg_frequencies_path, headers=['genome'] + scgs_ordered_based_on_frequency, keys_order=contigs_dbs_ordered_based_on_num_scgs)

        self.run.info('SCG frequencies across contigs dbs', self.report_scg_frequencies_path)


    def get_best_scg_name_for_metagenome_mode(self):
        """Identifies the best SCG to use for taxonomy.

        Although in reality 'best' means pretty much nothing here as this
        function simply returns the most frequent SCG. In an ideal world
        there would be mulitple modes, one that takes into consideration
        the prevalence of each SCG across each `self.genome`. If we ever
        implement that, let's leave a FIXME here only to be remove by
        that hero.
        """

        scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies = self.get_scg_frequencies()

        return scgs_ordered_based_on_frequency[0]


    def get_scg_taxonomy_super_dict_for_metagenomes(self):
        """Generates the taxonomy super dict, the primary data structure for this class.

        The difference between this function and `get_scg_taxonomy_super_dict` in `SCGTaxonomyEstimatorSingle`
        is that this one aggregates multiple `scg_taxonomy_super_dict` instances into a larger dictionary.

        Returns
        =======
        scg_taxonomy_super_dict: dict
            This is a complex dictionary that looks like this:

            ----8<----8<----8<----8<----8<----8<----8<----8<----
                {
                  "USA0001": { # <- name given for this entry in internal/external genomes file
                    "taxonomy": {
                      "USA0001_01": {  # <- project name in contigs database
                        "scgs": {
                          "69898": {
                            "gene_callers_id": 69898,
                            "gene_name": "Ribosomal_L16",
                            "accession": "CONSENSUS",
                            "percent_identity": "100.0",
                            "t_domain": "Bacteria",
                            "t_phylum": "Bacteroidota",
                            "t_class": "Bacteroidia",
                            "t_order": "Bacteroidales",
                            "t_family": "Rikenellaceae",
                            "t_genus": "Alistipes",
                            "t_species": "Alistipes shahii",
                            "tax_hash": "7e3cc7fc"
                          },
                          "55413": {
                            "gene_callers_id": 55413,
                            "gene_name": "Ribosomal_L16",
                            "accession": "CONSENSUS",
                            "percent_identity": "100.0",
                            "t_domain": "Bacteria",
                            "t_phylum": "Bacteroidota",
                            "t_class": "Bacteroidia",
                            "t_order": "Bacteroidales",
                            "t_family": "Tannerellaceae",
                            "t_genus": "Parabacteroides",
                            "t_species": null,
                            "tax_hash": "d16f9017"
                          }
                        },
                        "metagenome_mode": true
                      }
                    },
                    "coverages": {
                      "69898": {
                        "USA0001_01": 18.146179401993354
                      },
                      "55413": {
                        "USA0001_01": 48.82825484764543
                      }
                    }
                  },

                  (...)

                }
            ----8<----8<----8<----8<----8<----8<----8<----8<----

            Coverages will be there only if `--compute-scg-coverages` flag was passed to the class.
        """

        scg_taxonomy_super_dict = {}

        if self.profile_dbs_available:
            self.run.info_single("Your metagenome file DOES contain profile databases, so anvi'o will turn on `--metagenome-mode`, "
                                 "set the SCG name to %s, and turn on `--compute-scg-coverages` flag. If this doesn't make sense, "
                                 "please adjust your input parameters." % (self.scg_name_for_metagenome_mode), nl_after=1)
        else:
            if self.metagenome_mode:
                self.run.info_single("Your metagenome file DOES NOT contain profile databases, but you asked anvi'o to estimate SCG "
                                     "taxonomy in metagenome mode. So be it. SCG name is set to %s." \
                                        % (self.scg_name_for_metagenome_mode), nl_after=1)
            else:
                self.run.info_single("Your metagenome file DOES NOT contain profile databases, and you haven't asked anvi'o to "
                                     "work in `--metagenome-mode`. Your contigs databases will be treated as geomes rather than "
                                     "metagenomes.", nl_after=1)

        self.progress.new("Recovering tax super dict", progress_total_items=len(self.metagenome_names))
        total_num_metagenomes = len(self.metagenome_names)
        for metagenome_name in self.metagenome_names:
            args = SCGTaxonomyArgs(self.args, format_args_for_single_estimator=True)

            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            if self.profile_dbs_available:
                args.metagenome_mode = True
                args.profile_db = self.metagenomes[metagenome_name]['profile_db_path']
                args.compute_scg_coverages = True
                args.scg_name_for_metagenome_mode = self.scg_name_for_metagenome_mode
            else:
                if self.metagenome_mode:
                    args.metagenome_mode = True
                    args.scg_name_for_metagenome_mode = self.scg_name_for_metagenome_mode
                else:
                    args.metagenome_mode = False
                    args.scg_name_for_metagenome_mode = None

            self.progress.update("[%d of %d] %s" % (self.progress.progress_current_item + 1, total_num_metagenomes, metagenome_name))
            scg_taxonomy_super_dict[metagenome_name] = SCGTaxonomyEstimatorSingle(args, progress=progress_quiet, run=run_quiet).get_scg_taxonomy_super_dict()

            self.progress.increment()

        self.progress.end()

        return scg_taxonomy_super_dict


    def get_scg_frequencies(self):
        """Calculates the SCG frequencies across `self.metagenomes`."""

        scg_frequencies = {}

        self.progress.new("Learning SCG frequencies across genomes", progress_total_items=len(self.metagenomes))
        for metagenome_name in self.metagenomes:
            self.progress.update("%s ..." % metagenome_name)
            scg_frequencies[metagenome_name] = {}

            args = SCGTaxonomyArgs(self.args, format_args_for_single_estimator=True)
            args.compute_scg_coverages = False
            args.contigs_db = self.metagenomes[metagenome_name]['contigs_db_path']

            e = SCGTaxonomyEstimatorSingle(args, progress=progress_quiet, run=run_quiet)
            for scg_name in self.ctx.default_scgs_for_taxonomy:
                scg_frequencies[metagenome_name][scg_name] = e.frequency_of_scgs_with_taxonomy[scg_name]

            self.progress.increment()

        self.progress.update("Finalizing data")

        scg_frequencies_across_contigs_dbs = [(scg_name, sum([scg_frequencies[genome_name][scg_name] for genome_name in scg_frequencies])) for scg_name in self.ctx.default_scgs_for_taxonomy]
        scgs_ordered_based_on_frequency = [frequency_tuple[0] for frequency_tuple in sorted(scg_frequencies_across_contigs_dbs, key = lambda x: x[1], reverse=True)]

        num_scgs_for_each_contigs_db = [(genome_name, sum(scg_frequencies[genome_name].values())) for genome_name in scg_frequencies]
        contigs_dbs_ordered_based_on_num_scgs = [frequency_tuple[0] for frequency_tuple in sorted(num_scgs_for_each_contigs_db, key = lambda x: x[1], reverse=True)]

        self.progress.end()

        return scgs_ordered_based_on_frequency, contigs_dbs_ordered_based_on_num_scgs, scg_frequencies


class SCGTaxonomyEstimatorSingle(SCGTaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.update_profile_db_with_taxonomy = A('update_profile_db_with_taxonomy')
        self.bin_id = A('bin_id')

        SCGTaxonomyArgs.__init__(self, self.args)

        self.ctx = ctx

        SanityCheck.__init__(self)

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Profile DB', self.profile_db_path)
        self.run.info('Metagenome mode', self.metagenome_mode)
        if self.metagenome_mode:
            self.run.info('SCG for metagenome', self.scg_name_for_metagenome_mode)

        # these dictionaries that will be initiated later
        self.contigs_db_project_name = "Unknown"
        self.scg_name_to_gene_caller_id_dict = {}
        self.frequency_of_scgs_with_taxonomy = {}
        self.gene_callers_id_to_scg_taxonomy_dict = {}
        self.split_name_to_gene_caller_ids_dict = {}
        self.gene_callers_id_to_split_name_dict = {}
        self.sample_names_in_profile_db = None

        self.initialized = False

        if not skip_init:
            self.init()


    def init(self):
        self.init_scg_data()

        if self.report_scg_frequencies_path:
            with open(self.report_scg_frequencies_path, 'w') as output:
                for scg_name, frequency in self.frequency_of_scgs_with_taxonomy.items():
                    output.write("%s\t%d\n" % (scg_name, frequency))

            self.run.info('SCG frequencies within the contigs db', self.report_scg_frequencies_path, nl_before=1)
            sys.exit()

        if self.profile_db_path:
            self.sample_names_in_profile_db = ProfileDatabase(self.profile_db_path).samples

        self.initialized = True


    def init_scg_data(self):
        """Initialize SCG taxonomy for the entire contigs database"""

        if not self.contigs_db_path:
            return None

        self.progress.new('Initializing')
        self.progress.update('SCG taxonomy dictionary')

        for scg_name in self.ctx.SCGs:
            self.scg_name_to_gene_caller_id_dict[scg_name] = set([])

        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        self.contigs_db_project_name = contigs_db.meta['project_name']
        scg_taxonomy_table = contigs_db.db.get_table_as_dict(t.scg_taxonomy_table_name)
        genes_in_splits = contigs_db.db.get_some_columns_from_table(t.genes_in_splits_table_name, "split, gene_callers_id")
        min_contig_length_in_contigs_db = contigs_db.db.get_max_value_in_column(t.contigs_info_table_name, "length", return_min_instead=True)
        contigs_db.disconnect()

        # this is important. before we begin, we need to filter out gene caller ids and splits from main dictionaries if
        # they shouldn't be there. read the warning below to see the utility of this step.
        if self.profile_db_path and self.compute_scg_coverages:
            split_names_in_profile_db = set(utils.get_all_item_names_from_the_database(self.profile_db_path))
            split_names_in_contigs_db = set([tpl[0] for tpl in genes_in_splits])
            splits_missing_in_profile_db = split_names_in_contigs_db.difference(split_names_in_profile_db)

            min_contig_length_in_profile_db = ProfileDatabase(self.profile_db_path).meta['min_contig_length']

            if len(splits_missing_in_profile_db):
                self.progress.reset()
                self.run.warning("Please note that anvi'o found %s splits in your contigs database. But only %s of them "
                                 "appeared ot be in the profile database. As a result, anvi'o will now remove the %s splits "
                                 "that occur only in the contigs db from all downstream analyses here (if you didn't use the flag "
                                 "`--compute-scg-coverages` this wouldn't have been necessary, but with the current settings "
                                 "this is really the best for everyone). Where is this difference coming from though? Well. This "
                                 "is often the case because the 'minimum contig length parameter' set during the `anvi-profile` "
                                 "step can exclude many contigs from downstream analyses (often for good reasons, too). For "
                                 "instance, in your case the minimum contig length goes as low as %s nts in your contigs database. "
                                 "Yet, the minimum contig length set in the profile databaes is %s nts. Hence the difference. Anvi'o "
                                 "hopes that this explaines some things." % (pp(len(split_names_in_contigs_db)),
                                                                             pp(len(split_names_in_profile_db)),
                                                                             pp(len(splits_missing_in_profile_db)),
                                                                             pp(min_contig_length_in_contigs_db),
                                                                             pp(min_contig_length_in_profile_db)))

                self.progress.update("Removing %s splits missing form the profile db" % pp(len(splits_missing_in_profile_db)))
                genes_in_splits = [tpl for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db]

                # so now we know the final list of split names and gene caller ids as they are stored in the updated
                # `genes_in_splits` variable. time to clean up the `scg_taxonomy_table` dictionary as well.
                final_set_of_gene_caller_ids = set([tpl[1] for tpl in genes_in_splits if tpl[0] not in splits_missing_in_profile_db])
                entry_ids_to_remove = [entry for entry in scg_taxonomy_table if scg_taxonomy_table[entry]['gene_callers_id'] not in final_set_of_gene_caller_ids]
                [scg_taxonomy_table.pop(e) for e in entry_ids_to_remove]

        # NOTE: This will modify the taxonomy strings read from the contigs database. see the
        # function header for `trim_taxonomy_dict_entry` for more information.
        if self.simplify_taxonomy_information:
            self.progress.update('SCG taxonomy dicts ... trimming main')
            for key in scg_taxonomy_table:
                scg_taxonomy_table[key] = self.trim_taxonomy_dict_entry(scg_taxonomy_table[key])

        self.progress.update('SCG taxonomy dicts ... building g->tax')
        for entry in scg_taxonomy_table.values():
            gene_callers_id = entry['gene_callers_id']
            self.gene_callers_id_to_scg_taxonomy_dict[gene_callers_id] = entry

        self.progress.update('SCG taxonomy dicts ... building tax->g')
        for entry in self.gene_callers_id_to_scg_taxonomy_dict.values():
            scg_gene_name = entry['gene_name']
            gene_callers_id = entry['gene_callers_id']
            self.scg_name_to_gene_caller_id_dict[scg_gene_name].add(gene_callers_id)

        self.progress.update('SCG taxonomy dicts ... building s->g + g->s')
        for split_name, gene_callers_id in genes_in_splits:
            if gene_callers_id not in self.gene_callers_id_to_scg_taxonomy_dict:
                continue

            if split_name not in self.split_name_to_gene_caller_ids_dict:
                self.split_name_to_gene_caller_ids_dict[split_name] = set()

            self.split_name_to_gene_caller_ids_dict[split_name].add(gene_callers_id)
            self.gene_callers_id_to_split_name_dict[gene_callers_id] = split_name

        self.progress.end()

        self.frequency_of_scgs_with_taxonomy = OrderedDict(sorted([(g, len(self.scg_name_to_gene_caller_id_dict[g])) for g in self.scg_name_to_gene_caller_id_dict], key = lambda x: x[1], reverse=True))

        if self.metagenome_mode or anvio.DEBUG:
            self.run.info_single("A total of %s single-copy core genes with taxonomic affiliations were successfully initialized "
                                 "from the contigs database 🎉 Following shows the frequency of these SCGs: %s." % \
                                            (pp(len(self.gene_callers_id_to_scg_taxonomy_dict)),
                                             ', '.join(["%s (%d)" % (g, self.frequency_of_scgs_with_taxonomy[g]) \
                                                                for g in self.frequency_of_scgs_with_taxonomy])), nl_before=1)


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


    def get_consensus_taxonomy(self, scg_taxonomy_dict):
        """Takes in a scg_taxonomy_dict, returns a final taxonomic string that summarize all"""

        if not len(scg_taxonomy_dict):
            return dict([(l, None) for l in self.ctx.levels_of_taxonomy])

        pd.set_option('mode.chained_assignment', None)

        scg_hits = list([v for v in scg_taxonomy_dict.values() if v['t_domain']])

        if not len(scg_hits):
            return self.get_blank_hit_template_dict()

        df = pd.DataFrame.from_records(scg_hits)

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


    def print_scg_taxonomy_hits_in_splits(self, hits, bin_name=None):
        self.progress.reset()
        self.run.warning(None, header='Hits for %s' % (bin_name if bin_name else "a bunch of splits"), lc="green")

        if len(hits):
            header = ['SCG', 'gene', 'pct id', 'taxonomy']
            table = []

            for hit in hits:
                taxon_text = ' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy])

                # if the hit we are working on sent here as 'consensus', we will color it up a bit so it shows up
                # more clearly in the debug output.
                if hit['gene_name'] == 'CONSENSUS':
                    taxon_text = terminal.c(taxon_text, color='red')

                    for field_name in ['gene_name', 'percent_identity', 'gene_callers_id']:
                        hit[field_name] = terminal.c(hit[field_name], color='red')

                table.append([hit['gene_name'], str(hit['gene_callers_id']), str(hit['percent_identity']), taxon_text])

            print(tabulate(table, headers=header, tablefmt="fancy_grid", numalign="right"))
        else:
            self.run.info_single("No hits :/")


    def get_scg_taxonomy_dict(self, gene_caller_ids, bin_name=None):
        scg_taxonomy_dict = {}

        improper_gene_caller_ids = [g for g in gene_caller_ids if g not in self.gene_callers_id_to_scg_taxonomy_dict]
        if improper_gene_caller_ids:
            raise ConfigError("Something weird is going on. Somehow anvi'o has a bunch of gene caller ids for which it is "
                              "supposed to estimate taxonomy. However, %d of them do not occur in a key dictionary. The code "
                              "here does not know what to suggest :( Apologies." % len(improper_gene_caller_ids))

        for gene_callers_id in gene_caller_ids:
            scg_taxonomy_dict[gene_callers_id] = self.gene_callers_id_to_scg_taxonomy_dict[gene_callers_id]
            scg_taxonomy_dict[gene_callers_id]["tax_hash"] = HASH(self.gene_callers_id_to_scg_taxonomy_dict[gene_callers_id])

        return scg_taxonomy_dict


    def estimate_for_list_of_splits(self, split_names=None, bin_name=None):
        """Estimate SCG taxonomy for a bunch of splits that belong to a single population.

           The purpose of this function is to to do critical things: identify SCGs we use for taxonomy in `split_names`,
           and generate a consensus taxonomy with the assumption that these are coming from splits that represents a
           single population.

           It will return a dictionary with multiple items, including a dictionary that contains the final consensus\
           taxonomy, another one that includes every SCG and their raw associations with taxon names (from which the\
           consensus taxonomy was computed), as well as information about how many SCGs were analyzed and supported the\
           consesnus.
        """

        if self.metagenome_mode:
            raise ConfigError("Someone is attempting to estimate taxonomy for a set of splits using a class inherited in "
                              "`metagenome mode`. If you are a programmer please note that it is best to use the member "
                              "function `estimate` directly.")

        consensus_taxonomy = None

        gene_caller_ids_of_interest = self.get_gene_caller_ids_for_splits(split_names)
        scg_taxonomy_dict = self.get_scg_taxonomy_dict(gene_caller_ids_of_interest)

        try:
            consensus_taxonomy = self.get_consensus_taxonomy(scg_taxonomy_dict)
            consensus_taxonomy['gene_name'] = 'CONSENSUS'
            consensus_taxonomy['percent_identity'] = '--'
            consensus_taxonomy['gene_callers_id'] = '--'

        except Exception as e:
            self.print_scg_taxonomy_hits_in_splits(list(scg_taxonomy_dict.values()))

            raise ConfigError("While trying to sort out the consensus taxonomy for %s anvi'o failed :( The list of SCG taxon hits that "
                              "caused the failure is printed in your terminal. But the actual error message that came from the depths "
                              "of the codebase was this: '%s'." % (('the bin "%s"' % bin_name) if bin_name else 'a bunch of splits', e))

        if anvio.DEBUG:
            self.print_scg_taxonomy_hits_in_splits(list(scg_taxonomy_dict.values()) + [consensus_taxonomy], bin_name)

        # set some useful information. `total_scgs` is the number of SCGs with taxonomy found in the collection of splits. the
        # `supporting_scgs` communicate how many of them supports the consensus taxonomy fully
        total_scgs = len(scg_taxonomy_dict)
        supporting_scgs = 0

        consensus_taxonomy_levels_occupied = [level for level in self.ctx.levels_of_taxonomy if consensus_taxonomy[level]]
        consensus_taxonomy_str = ' / '.join([consensus_taxonomy[level] for level in consensus_taxonomy_levels_occupied])

        for scg_taxonomy_hit in scg_taxonomy_dict.values():
            scg_taxonomy_hit_str = ' / '.join([str(scg_taxonomy_hit[level]) for level in consensus_taxonomy_levels_occupied])

            if scg_taxonomy_hit_str == consensus_taxonomy_str:
                scg_taxonomy_hit['supporting_consensus'] = True
                supporting_scgs += 1
            else:
                scg_taxonomy_hit['supporting_consensus'] = False

        return {'consensus_taxonomy': consensus_taxonomy,
                'scgs': scg_taxonomy_dict,
                'total_scgs': total_scgs,
                'supporting_scgs': supporting_scgs,
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

        scg_frequencies = self.frequency_of_scgs_with_taxonomy.values()
        if len([sf for sf in scg_frequencies if sf > 1]) * 100 / len(scg_frequencies) > 20:
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

        splits_in_contigs_database = self.split_name_to_gene_caller_ids_dict.keys()
        contigs_db_taxonomy_dict[self.contigs_db_project_name] = self.estimate_for_list_of_splits(split_names=splits_in_contigs_database,
                                                                                                  bin_name=self.contigs_db_project_name)
        return contigs_db_taxonomy_dict


    def estimate_for_contigs_db_for_metagenome(self):
        """Treat a given contigs database as a metagenome.

           This function deserves some attention. It relies on a single SCG to estimate the composition of a metagenome.
           For instance, its sister function, `estimate_for_contigs_db_for_genome`, works with a list of splits that are
           assumed to belong to the same genome. In which case a consensus taxonomy learned from all SCGs is most
           appropriate. In this case, however, we don't know which split will go together, hence, we can't pull together
           SCGs to learn a consensus taxonomy for independent populations in the metagenome. The best we can do is to stick
           with a single SCG with the hope that (1) it will cut through as many populations as possible and (2) will have
           reasonable power to resolve taxonomy all by itself. These independent assumptions will both work in some cases
           and both fail in others.
        """

        # we first need to decide which SCG we should use to survey taxonomy
        most_frequent_scg = next(iter(self.frequency_of_scgs_with_taxonomy))
        if self.scg_name_for_metagenome_mode:
            frequency_of_user_chosen_scg = self.frequency_of_scgs_with_taxonomy[self.scg_name_for_metagenome_mode]
            frequency_of_most_frequent_scg = self.frequency_of_scgs_with_taxonomy[most_frequent_scg]

            if frequency_of_user_chosen_scg < frequency_of_most_frequent_scg:
                additional_note = " And just so you know, there is another SCG that was observed more times (i.e., %s; %d times)\
                                   in this metagenome compared to yours (i.e., %d times). You're the boss, of course." %\
                                            (most_frequent_scg, frequency_of_most_frequent_scg, frequency_of_user_chosen_scg)
            else:
                additional_note = ""

            self.run.warning("As per your request anvi'o set '%s' to be THE single-copy core gene to survey your metagenome for its "
                             "taxonomic composition.%s" % (self.scg_name_for_metagenome_mode, additional_note))
        else:
            self.scg_name_for_metagenome_mode = most_frequent_scg

            self.run.warning("Anvi'o automatically set '%s' to be THE single-copy core gene to survey your metagenome for its "
                             "taxonomic composition. If you are not happy with that, you could change it with the parameter "
                             "`--scg-name-for-metagenome-mode`." % (self.scg_name_for_metagenome_mode))

        gene_caller_ids_of_interest = self.scg_name_to_gene_caller_id_dict[self.scg_name_for_metagenome_mode]
        scg_taxonomy_dict = self.get_scg_taxonomy_dict(gene_caller_ids=gene_caller_ids_of_interest,
                                                       bin_name=self.contigs_db_project_name)

        return {self.contigs_db_project_name: {'scgs': scg_taxonomy_dict,
                                               'metagenome_mode': True}}


    def get_scg_taxonomy_super_dict(self):
        """Function that returns the `scg_taxonomy_super_dict`.

           `scg_taxonomy_super_dict` contains a wealth of information regarding samples, SCGs,
           SCG taxonomic affiliations, consensus taxonomy, and coverages of SCGs across samples.
        """
        scg_taxonomy_super_dict = {}

        if not self.initialized:
            self.init()

        if self.profile_db_path and not self.metagenome_mode:
            scg_taxonomy_super_dict['taxonomy'] = self.estimate_for_bins_in_collection()
        elif not self.profile_db_path and not self.metagenome_mode:
            scg_taxonomy_super_dict['taxonomy'] = self.estimate_for_contigs_db_for_genome()
        elif self.metagenome_mode:
            scg_taxonomy_super_dict['taxonomy'] = self.estimate_for_contigs_db_for_metagenome()
        else:
            raise ConfigError("This class doesn't know how to deal with that yet :/")

        if self.compute_scg_coverages and self.metagenome_mode:
            scg_taxonomy_super_dict['coverages'] = self.get_scg_coverages_across_samples_dict_in_metagenome_mode(scg_taxonomy_super_dict)
        elif self.compute_scg_coverages and not self.metagenome_mode:
            scg_taxonomy_super_dict['coverages'] = self.get_scg_coverages_across_samples_dict_in_genome_mode(scg_taxonomy_super_dict)
        else:
            scg_taxonomy_super_dict['coverages'] = None

        return scg_taxonomy_super_dict


    def estimate(self):
        scg_taxonomy_super_dict = self.get_scg_taxonomy_super_dict()

        if self.update_profile_db_with_taxonomy:
            self.add_taxonomy_as_additional_layer_data(scg_taxonomy_super_dict)

        self.print_scg_taxonomy_super_dict(scg_taxonomy_super_dict)

        if self.output_file_path:
            self.store_scg_taxonomy_super_dict(scg_taxonomy_super_dict)


    def print_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        self.progress.reset()

        if self.collection_name:
            self.run.warning(None, header='Estimated taxonomy for collection "%s"' % self.collection_name, lc="green")
        elif self.metagenome_mode:
            self.run.warning(None, header='Taxa in metagenome "%s"' % self.contigs_db_project_name, lc="green")
        else:
            self.run.warning(None, header='Estimated taxonomy for "%s"' % self.contigs_db_project_name, lc="green")

        d = self.get_print_friendly_scg_taxonomy_super_dict(scg_taxonomy_super_dict)

        ordered_bin_names = sorted(list(d.keys()))

        if self.metagenome_mode:
            header = ['percent_identity', 'taxonomy']
        else:
            header = ['', 'total_scgs', 'supporting_scgs', 'taxonomy']

        # if we are in `--compute-scg-coverages` mode, and more than 5 sample names, we are in trouble since they will\
        # unlikely fit into the display while printing them. so here we will cut it to make sure things look OK.
        samples_not_shown = 0
        sample_names_to_display = None
        if self.compute_scg_coverages:
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
                sorted_bin_coverage_tuples = sorted([(bin_name, sum([(d[bin_name]['coverages'][sample_name] if d[bin_name]['supporting_scgs'] else 0) for sample_name in self.sample_names_in_profile_db])) for bin_name in d], key=lambda x: x[1], reverse=True)
            ordered_bin_names = [tpl[0] for tpl in sorted_bin_coverage_tuples]


        table = []
        for bin_name in ordered_bin_names:
            bin_data = d[bin_name]

            # set the taxonomy text depending on how much room we have. if there are sample coverages, keep it simple,
            # otherwise show the entire taxonomy text.
            if self.compute_scg_coverages:
                taxon_text_l = ['(%s) %s' % (l.split('_')[1][0], bin_data[l]) for l in self.ctx.levels_of_taxonomy[::-1] if bin_data[l]]
                taxon_text = taxon_text_l[0] if taxon_text_l else '(NA) NA'
            else:
                taxon_text = ' / '.join([bin_data[l] if bin_data[l] else '' for l in self.ctx.levels_of_taxonomy])

            # setting up the table columns here.
            if self.metagenome_mode:
                row = [bin_name, str(bin_data['percent_identity']), taxon_text]
            else:
                row = [bin_name, str(bin_data['total_scgs']), str(bin_data['supporting_scgs']), taxon_text]

            # if there are coverages, add samples to the display too
            if self.compute_scg_coverages:
                row += [d[bin_name]['coverages'][sample_name] for sample_name in sample_names_to_display]

            if samples_not_shown:
                row += ['... %d more' % len(samples_not_shown)]

            table.append(row)

        # if we are not in metagenome mode let's sort the output table based on total and
        # supporting SCGs
        if not self.metagenome_mode:
            table = sorted(table, key=lambda x: (int(x[1]), int(x[2])), reverse=True)

        print(tabulate(table, headers=header, tablefmt="fancy_grid", numalign="right"))


    def store_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        d = self.get_print_friendly_scg_taxonomy_super_dict(scg_taxonomy_super_dict)

        if self.metagenome_mode:
            headers = ['scg_name', 'percent_identity']
        else:
            headers = ['bin_name', 'total_scgs', 'supporting_scgs']

        headers += self.ctx.levels_of_taxonomy

        if self.compute_scg_coverages:
            headers_for_samples = sorted(self.sample_names_in_profile_db)
        else:
            headers_for_samples = []

        with open(self.output_file_path, 'w') as output:
            output.write('\t'.join(headers + headers_for_samples) + '\n')
            for item in d:
                line = [item] + [d[item][h] for h in headers[1:]]

                if self.compute_scg_coverages:
                    for sample_name in headers_for_samples:
                        line.append(d[item]['coverages'][sample_name])

                output.write('\t'.join([str(f) for f in line]) + '\n')

        self.run.info("Output file", self.output_file_path, nl_before=1)


    def get_print_friendly_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        d = {}

        if self.metagenome_mode:
            for scg_hit in scg_taxonomy_super_dict['taxonomy'][self.contigs_db_project_name]['scgs'].values():
                scg_hit_name = '%s_%d' % (scg_hit['gene_name'], scg_hit['gene_callers_id'])
                d[scg_hit_name] = scg_hit

                if self.compute_scg_coverages:
                    d[scg_hit_name]['coverages'] = scg_taxonomy_super_dict['coverages'][scg_hit['gene_callers_id']]
        else:
            for bin_name in scg_taxonomy_super_dict['taxonomy']:
                d[bin_name] = scg_taxonomy_super_dict['taxonomy'][bin_name]['consensus_taxonomy']
                d[bin_name]['total_scgs'] = scg_taxonomy_super_dict['taxonomy'][bin_name]['total_scgs']
                d[bin_name]['supporting_scgs'] = scg_taxonomy_super_dict['taxonomy'][bin_name]['supporting_scgs']

                if self.compute_scg_coverages:
                    d[bin_name]['coverages'] = scg_taxonomy_super_dict['coverages'][bin_name]

        return d


    def add_taxonomy_as_additional_layer_data(self, scg_taxonomy_super_dict):
        """A function that adds taxonomy to additional data tables of a given profile
           database. This will only work in metagenome mode."""

        if not self.metagenome_mode or not self.compute_scg_coverages:
            return

        self.progress.new("Adding summary taxonomy for samples")
        self.progress.update('...')

        scgs_dict = list(scg_taxonomy_super_dict['taxonomy'].values())[0]['scgs']

        # at this stage each scgs_dict entry will look like this, and most critically will
        # have the same SCG:
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
        # this will enable us to learn the which SCG has been used to calculate coverage
        # information by only looking at a single entry:
        scg_name = list(scgs_dict.values())[0]['gene_name']

        # the might for loop to go through all taxonomic levels one by one
        for level in constants.ctx.levels_of_taxonomy[::-1]:
            # setting the data group early on:
            data_group = '%s_%s' % (scg_name, level[2:])
            self.progress.update('Working on %s-level data' % level)
            data_dict = {}
            data_keys_list = set([])
            for sample_name in self.sample_names_in_profile_db:
                data_dict[sample_name] = Counter()
                for gene_callers_id in scgs_dict:
                    # starting with a tiny hack to fill in missing values. here we first find
                    # the most highly resolved level of taxonomy that is not null for this
                    # particular scg taxonomy
                    i = 0
                    for i in range(constants.ctx.levels_of_taxonomy.index(level), 0, -1):
                        if scgs_dict[gene_callers_id][constants.ctx.levels_of_taxonomy[i]]:
                            break

                    # just some abbreviations
                    l = constants.ctx.levels_of_taxonomy[i][2:]
                    m = scgs_dict[gene_callers_id][constants.ctx.levels_of_taxonomy[i]]

                    # if the best level we found in the previous step is matching to the level
                    # set by the main for loop, we're good to go with that name:
                    if level == constants.ctx.levels_of_taxonomy[i]:
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
                    data_dict[sample_name][key] += scg_taxonomy_super_dict['coverages'][gene_callers_id][sample_name]
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


    def get_split_names_for_scg_taxonomy_super_dict(self, scg_taxonomy_super_dict):
        """Returns a list of split names associated with SCGs found in a scg_taxonomy_super_dict."""

        if 'scgs' not in list(scg_taxonomy_super_dict['taxonomy'].values())[0]:
            raise ConfigError("Someone called this function with something that doesn't look like the kind "
                              "of input data it was expecting (sorry for the vagueness of the message, but "
                              "anvi'o hopes that will be able to find out why it is happening).")

        split_names = set([])

        for entry_name in scg_taxonomy_super_dict['taxonomy']:
            for gene_callers_id in scg_taxonomy_super_dict['taxonomy'][entry_name]['scgs']:
                split_names.add(self.gene_callers_id_to_split_name_dict[gene_callers_id])

        return split_names


    def get_scg_coverages_across_samples_dict_in_genome_mode(self, scg_taxonomy_super_dict):
        self.progress.reset()
        self.run.info_single("Anvi'o will now attempt to recover SCG coverages in GENOME MODE from the profile "
                             "database, which contains %d samples." % (len(self.sample_names_in_profile_db)), nl_before=1, nl_after=1)

        scg_coverages_across_samples_dict = self.get_scg_coverages_across_samples_dict(scg_taxonomy_super_dict)

        bin_avg_coverages_across_samples_dict = {}
        for bin_name in scg_taxonomy_super_dict['taxonomy']:
            bin_avg_coverages_across_samples_dict[bin_name] = dict([(sample_name, None) for sample_name in self.sample_names_in_profile_db])
            for sample_name in self.sample_names_in_profile_db:
                average_coverage_across_samples = [scg_coverages_across_samples_dict[gene_callers_id][sample_name] for gene_callers_id in scg_taxonomy_super_dict['taxonomy'][bin_name]['scgs']]
                if average_coverage_across_samples:
                    bin_avg_coverages_across_samples_dict[bin_name][sample_name] = np.mean(average_coverage_across_samples)

        self.run.warning("Anvi'o has just finished recovering SCG coverages from the profile database to estimate "
                         "the average coverage of your bins across your samples. Please note that anvi'o SCG taxonomy "
                         "framework is using only %d SCGs to estimate taxonomy. Which means, even a highly complete bin "
                         "may be missing all of them. In which case, the coverage of that bin will be `None` across all "
                         "your samples. The best way to prevent any misleading insights is take these results with a "
                         "huge grain of salt, and use the `anvi-summarize` output for critical applications.",
                         header="FRIENDLY REMINDER", lc="blue")

        return bin_avg_coverages_across_samples_dict


    def get_scg_coverages_across_samples_dict_in_metagenome_mode(self, scg_taxonomy_super_dict):
        """Get SCG coverages in metagenome mode."""

        if not self.metagenome_mode:
            raise ConfigError("You're calling the wrong function. Your class is not in metagenome mode.")

        self.progress.reset()
        self.run.info_single("Anvi'o will now attempt to recover SCG coverages from the profile database, which "
                             "contains %d samples." % (len(self.sample_names_in_profile_db)), nl_before=1, nl_after=1)

        return self.get_scg_coverages_across_samples_dict(scg_taxonomy_super_dict)


    def get_scg_coverages_across_samples_dict(self, scg_taxonomy_super_dict):
        """Get SCG coverages"""
        scg_coverages_across_samples_dict = {}

        self.progress.new('Recovering coverages')
        self.progress.update('Learning all split names affiliated with SCGs ..')
        split_names_of_interest = self.get_split_names_for_scg_taxonomy_super_dict(scg_taxonomy_super_dict)
        self.progress.end()

        # initialize split coverages for splits that have anything to do with our SCGs
        args = copy.deepcopy(self.args)
        args.split_names_of_interest = split_names_of_interest
        args.collection_name = None
        profile_db = ProfileSuperclass(args, p=self.progress, r=run_quiet)
        profile_db.init_split_coverage_values_per_nt_dict()

        # recover all gene caller ids that occur in our taxonomy estimation dictionary
        # and ge their coverage stats from the profile super
        gene_caller_ids_of_interest = set([])
        for bin_name in scg_taxonomy_super_dict['taxonomy']:
            for gene_callers_id in scg_taxonomy_super_dict['taxonomy'][bin_name]['scgs']:
                gene_caller_ids_of_interest.add(gene_callers_id)

        # at this point we have everything. splits of interest are loaded in memory in `profile_db`, and we know
        # which gene caller ids we are interested in recovering coverages for. the way to access to gene coverages
        # is a bit convoluted in the dbops for historical reasons, but it is quite straightforward. the most
        # weird part is that we need a copy of a contigs super. so we will start with that:
        self.progress.new("Recovering SCG coverages")
        self.progress.update("Initiating the contigs super class")
        contigs_db = ContigsSuperclass(self.args, r=run_quiet, p=progress_quiet)

        for split_name in split_names_of_interest:
            self.progress.update("Working with %s" % split_name)
            # note for the curious: yes, here we are sending the same gene caller ids of interest over and over to
            # the `get_gene_level_coverage_stats` for each split, but that function is smart enough to not spend any
            # time on those gene caller ids that do not occur in the split name we are interested in.
            all_scg_stats_in_split = profile_db.get_gene_level_coverage_stats(split_name, contigs_db, gene_caller_ids_of_interest=gene_caller_ids_of_interest)

            for scg_stats in all_scg_stats_in_split.values():
                for entry in scg_stats.values():
                    gene_callers_id = int(entry['gene_callers_id'])
                    sample_name = entry['sample_name']
                    coverage = entry['non_outlier_mean_coverage']

                    if gene_callers_id not in scg_coverages_across_samples_dict:
                        scg_coverages_across_samples_dict[gene_callers_id] = dict([(sample_name, 0) for sample_name in self.sample_names_in_profile_db])

                    scg_coverages_across_samples_dict[gene_callers_id][sample_name] = coverage

        self.progress.end()

        return scg_coverages_across_samples_dict


class SetupLocalSCGTaxonomyData(SCGTaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        SCGTaxonomyArgs.__init__(self, self.args)

        # user accessible variables
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.reset = A("reset") # complete start-over (downloading everything from GTDB)
        self.redo_databases = A("redo_databases") # just redo the databaes
        self.num_threads = A('num_threads')

        self.ctx = ctx

        SanityCheck.__init__(self)


    def setup(self):
        """This function downloads all GTDB files necessary to setup the SCG databases anvi'o will rely upon.

           In addition to downloading the original files, the setup will make sure everything, including the
           DIAMOND search databases are in place.
        """

        if not anvio.DEBUG:
            self.check_initial_directory_structure()

        self.run.warning("Please remember that the data anvi'o uses for SCG taxonomy is a "
                         "courtesy of The Genome Taxonomy Database (GTDB), an initiative to establish a "
                         "standardised microbial taxonomy based on genome phylogeny, primarly funded by "
                         "tax payers in Australia. Please don't forget to cite the original work, "
                         "doi:10.1038/nbt.4229 by Parks et al to explicitly mention the source of databases "
                         "anvi'o relies upon to estimate genome level taxonomy. If you are not sure how it "
                         "should look like in your methods sections, anvi'o developers will be happy to "
                         "help you if you can't find any published example to get inspiration.", lc = 'yellow')

        # let's try to figure out what do we want to do here. download files or just create/re-create
        # databases.
        a_scg_fasta, a_scg_database = list(self.ctx.SCGs.values())[0]['fasta'] + '.gz', list(self.ctx.SCGs.values())[0]['db']

        if os.path.exists(a_scg_fasta) and os.path.exists(a_scg_database) and self.redo_databases:
            self.run.warning("Anvi'o is removing all the previous databases so it can regenerate them from their "
                             "ashes.")

            db_paths = [v['db'] for v in self.ctx.SCGs.values()]
            for db_path in db_paths:
                os.remove(db_path) if os.path.exists(db_path) else None

        elif os.path.exists(a_scg_fasta) and os.path.exists(a_scg_database) and not self.reset:
            raise ConfigError("It seems you have both the FASTA files and the search databases for anvi'o SCG taxonomy "
                              "in place. If you want to start from scratch and download everything from GTDB clean, use "
                              "the flag `--reset`.")

        elif self.reset:
            self.run.info("Local directory to setup", self.ctx.SCGs_taxonomy_data_dir)
            self.run.info("Reset the directory first", self.reset, mc="red")
            self.run.info("Remote database", self.ctx.target_database, nl_before=1, mc="green")
            self.run.info("Remote URL to download files", self.ctx.target_database_URL)
            self.run.info("Remote files of interest", ', '.join(self.ctx.target_database_files))

            self.progress.new("%s setup" % self.ctx.target_database)
            self.progress.update("Reading the VERSION file...")
            content = utils.get_remote_file_content(self.ctx.target_database_URL + 'VERSION')
            version, release_date  = content.strip().split('\n')[0].strip(), content.strip().split('\n')[2].strip()
            self.progress.end()

            self.run.info("%s release found" % self.ctx.target_database, "%s (%s)" % (version, release_date), mc="green")

            self.download_and_format_files()

        elif os.path.exists(a_scg_fasta) and not os.path.exists(a_scg_database):
            self.run.warning("Anvi'o found your FASTA files in place, but not the databases. Now it will generate all "
                             "the search databases using the existing FASTA files.")

        self.create_search_databases()

        if not anvio.DEBUG:
            self.clean_up()


    def check_initial_directory_structure(self):
        if os.path.exists(self.ctx.SCGs_taxonomy_data_dir):
            if self.reset:
                shutil.rmtree(self.ctx.SCGs_taxonomy_data_dir)
                self.run.warning('The existing directory for SCG taxonomy data dir has been removed. Just so you know.')
                filesnpaths.gen_output_directory(self.ctx.SCGs_taxonomy_data_dir)
        else:
            filesnpaths.gen_output_directory(self.ctx.SCGs_taxonomy_data_dir)


    def download_and_format_files(self):
        # let's be 100% sure.
        os.remove(self.ctx.accession_to_taxonomy_file_path) if os.path.exists(self.ctx.accession_to_taxonomy_file_path) else None

        for remote_file_name in self.ctx.target_database_files:
            remote_file_url = '/'.join([self.ctx.target_database_URL, remote_file_name])
            local_file_path = os.path.join(self.ctx.SCGs_taxonomy_data_dir, remote_file_name)

            utils.download_file(remote_file_url, local_file_path, progress=self.progress, run=self.run)

            if local_file_path.endswith('individual_genes.tar.gz'):
                self.progress.new("Downloaded file patrol")
                self.progress.update("Unpacking file '%s'..." % os.path.basename(local_file_path))
                shutil.unpack_archive(local_file_path, extract_dir=self.ctx.msa_individual_genes_dir_path)
                os.remove(local_file_path)
                self.progress.end()

            if local_file_path.endswith('_taxonomy.tsv'):
                with open(self.ctx.accession_to_taxonomy_file_path, 'a') as f:
                    f.write(open(local_file_path).read())
                    os.remove(local_file_path)

        fasta_file_paths = glob.glob(self.ctx.msa_individual_genes_dir_path + '/*.faa')

        if not fasta_file_paths:
            raise ConfigError("Something weird happened while anvi'o was trying to take care of the files "
                              "it downloaded from GTDB. Please let a developer know about this unless it is "
                              "not already reported in our issue tracker at Github :(")

        # files are done, but some of the FASTA files contain alignments solely composed of
        # gap characters :/ we will have to remove them to avoid fuck-ups in downstream
        # analyses
        self.progress.new("Clean up")
        for fasta_file_path in fasta_file_paths:
            self.progress.update("Looking for only-gap sequences from '%s'..." % os.path.basename(fasta_file_path))
            total_num_sequences, num_sequences_removed = utils.remove_sequences_with_only_gaps_from_fasta(fasta_file_path, fasta_file_path + '_CLEAN.fa', inplace=True)

            if num_sequences_removed:
                self.progress.reset()
                self.run.info_single('%d of %d seq in %s were all gaps and removed.' % (num_sequences_removed, total_num_sequences, os.path.basename(fasta_file_path)))

        # Start formatting things.
        self.progress.update("Checking output directory to store files ...")
        filesnpaths.is_output_dir_writable(os.path.dirname(self.ctx.search_databases_dir_path))
        filesnpaths.gen_output_directory(self.ctx.search_databases_dir_path, delete_if_exists=True, dont_warn=True)

        # We will be working with the files downloaded in whatever directory before. The first thing is to check
        # whether whether FASTA files in the directory are suitable for the conversion
        self.progress.update("Checking the conversion dict and FASTA files ...")
        msa_individual_gene_names_required = []
        [msa_individual_gene_names_required.extend(n) for n in locally_known_HMMs_to_remote_FASTAs.values()]

        fasta_file_paths = glob.glob(self.ctx.msa_individual_genes_dir_path + '/*.faa')
        msa_individual_gene_names_downloaded = [os.path.basename(f) for f in fasta_file_paths]

        missing_msa_gene_names = [n for n in msa_individual_gene_names_required if n not in msa_individual_gene_names_downloaded]
        if missing_msa_gene_names:
            self.progress.reset()
            raise ConfigError("Big trouble :( Anvi'o uses a hard-coded dictionary to convert locally known "
                              "HMM models to FASTA files reported by GTDB project. It seems something has changed "
                              "and %d of the FASTA files expected to be in the download directory are not there. "
                              "Here is that list: '%s'. Someone needs to update the codebase by reading the "
                              "appropriate documentation. If you are a user, you can't do much at this point but "
                              "contacting the developers :( Anvi'o will keep the directory that contains all the "
                              "downloaded files to update the conversion dictionary. Here is the full path to the "
                              "output: %s" % (len(missing_msa_gene_names), ', '.join(missing_msa_gene_names), self.ctx.msa_individual_genes_dir_path))
        else:
            self.progress.reset()
            self.run.info_single("Good news! The conversion dict and the FASTA files it requires seem to be in place. "
                                 "Anvi'o is now ready to to merge %d FASTA files that correspond to %d SCGs, and "
                                 "create individual search databases for them." % \
                                        (len(msa_individual_gene_names_required), len(locally_known_HMMs_to_remote_FASTAs)), nl_before=1, nl_after=1, mc="green")

        # Merge FASTA files that should be merged. This is defined in the conversion dictionary.
        for SCG in locally_known_HMMs_to_remote_FASTAs:
            self.progress.update("Working on %s ..." % (SCG))

            files_to_concatenate = [os.path.join(self.ctx.msa_individual_genes_dir_path, f) for f in locally_known_HMMs_to_remote_FASTAs[SCG]]
            FASTA_file_for_SCG = os.path.join(self.ctx.search_databases_dir_path, SCG)

            # concatenate from the dictionary into the new destination
            utils.concatenate_files(FASTA_file_for_SCG, files_to_concatenate)

            # compress the FASTA file
            utils.gzip_compress_file(FASTA_file_for_SCG)

        self.progress.end()


    def create_search_databases(self):
        """Creates all the search databases"""

        self.progress.new("Creating search databases")
        self.progress.update("Removing any database that still exists in the output directory...")
        [os.remove(database_path) for database_path in [s['db'] for s in self.ctx.SCGs.values()] if os.path.exists(database_path)]

        # compresssing and decompressing FASTA files changes their hash and make them look like
        # modified in git. to avoid that, we will do the database generation in a temporary directory.
        temp_dir = filesnpaths.get_temp_directory_path()

        self.progress.update("Copying FASTA files to %s ..." % (temp_dir))
        # the following line basically returns a dictionary that shows the new path
        # of the FASTA file under temp_dir for a given SCG .. apologies for the
        # incomprehensible list comprehension
        new_paths = dict([(os.path.basename(fasta_path), shutil.copy((fasta_path + '.gz'), os.path.join(temp_dir, os.path.basename(fasta_path) + '.gz'))) for fasta_path in [s['fasta'] for s in self.ctx.SCGs.values()]])

        missing_FASTA_files = [SCG for SCG in self.ctx.SCGs if not os.path.exists(new_paths[SCG])]
        if len(missing_FASTA_files):
            raise ConfigError("Weird news :( Anvi'o is missing some FASTA files that were supposed to be somewhere. Since this "
                              "can't be your fault, it is not easy to advice what could be the solution to this. But you can "
                              "always try to re-run `anvi-setup-scg-databases` with `--reset` flag.")

        self.progress.update("Decompressing FASTA files in %s" % (temp_dir))
        new_paths = dict([(SCG, utils.gzip_decompress_file(new_paths[SCG], keep_original=False)) for SCG in new_paths])

        # Merge FASTA files that should be merged. This is defined in the conversion dictionary.
        for SCG in self.ctx.SCGs:
            self.progress.update("Working on %s in %d threads" % (SCG, self.num_threads))

            FASTA_file_path_for_SCG = new_paths[SCG]

            # create a diamond search database for `FASTA_file_path_for_SCG`
            diamond = Diamond(query_fasta=FASTA_file_path_for_SCG, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
            diamond.makedb(output_file_path=self.ctx.SCGs[SCG]['db'])

            if not os.path.exists(self.ctx.SCGs[SCG]['db']):
                raise ConfigError("Something went wrong and DIAMOND did not create the database file it was supposed to "
                                  "for %s :(" % SCG)

        shutil.rmtree(temp_dir)

        self.progress.end()
        self.run.info_single("Every FASTA is now turned into a fancy search database. It means you are now allowed to run "
                             "`anvi-run-scg-taxonomy` on anvi'o contigs databases. This workflow is very new, and there are "
                             "caveats to it just like every other computational approach you use to make sense of complex 'omics "
                             "data. To better understand those caveats you should read our online documentation a bit. If you see "
                             "things that concerns you, please let anvi'o developers know. They love bad news. If you get good "
                             "results from this workflow, thank to those who contributed to the GTDB.", nl_after=1, mc="green")


    def clean_up(self):
        for file_path in [os.path.join(self.ctx.SCGs_taxonomy_data_dir, 'diamond-log-file.txt')]:
            if os.path.exists(file_path):
                os.remove(file_path)

        for dir_path in [self.ctx.msa_individual_genes_dir_path]:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)


class PopulateContigsDatabaseWithSCGTaxonomy(SCGTaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        SCGTaxonomyArgs.__init__(self, self.args)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.write_buffer_size = int(A('write_buffer_size') if A('write_buffer_size') is not None else 1000)
        self.contigs_db_path = A('contigs_db')
        self.num_parallel_processes = int(A('num_parallel_processes')) if A('num_parallel_processes') else 1
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        self.ctx = ctx

        self.max_target_seqs = 20
        self.evalue = float(A('e_value')) if A('e_value') else 1e-05
        self.min_pct_id = float(A('min_percent_identity')) if A('min_percent_identity') else 90

        SanityCheck.__init__(self)

        self.taxonomy_dict = OrderedDict()

        self.mutex = multiprocessing.Lock()


    def get_SCG_sequences_dict_from_contigs_db(self):
        """Returns a dictionary of all HMM hits per SCG of interest"""

        contigs_db = ContigsSuperclass(self.args, r=run_quiet, p=progress_quiet)
        splits_dict = {contigs_db.a_meta['project_name']: list(contigs_db.splits_basic_info.keys())}

        s = hmmops.SequencesForHMMHits(self.args.contigs_db, sources=self.ctx.hmm_source_for_scg_taxonomy, run=run_quiet, progress=progress_quiet)
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits(splits_dict, return_amino_acid_sequences=True)
        hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set(self.ctx.default_scgs_for_taxonomy))

        if not len(hmm_sequences_dict):
            return None

        self.progress.reset()
        self.run.info('Num relevant SCGs in contigs db', '%s' % (pp(len(hmm_sequences_dict))))

        scg_sequences_dict = {}
        for entry_id in hmm_sequences_dict:
            entry = hmm_sequences_dict[entry_id]

            scg_name = entry['gene_name']
            if scg_name in scg_sequences_dict:
                scg_sequences_dict[scg_name][entry_id] = entry
            else:
                scg_sequences_dict[scg_name] = {entry_id: entry}

        return scg_sequences_dict


    def populate_contigs_database(self):
        """Populates SCG taxonomy tables in a contigs database"""

        # get an instnce for the tables for taxonomy early on.
        self.tables_for_taxonomy = TableForSCGTaxonomy(self.contigs_db_path, self.run, self.progress)

        # get the dictionary that shows all hits for each SCG of interest
        self.progress.new('Contigs bleep bloop')
        self.progress.update('Recovering the SCGs dictionary')
        scg_sequences_dict = self.get_SCG_sequences_dict_from_contigs_db()
        self.progress.end()

        if not scg_sequences_dict:
            self.run.warning("This contigs database contains no single-copy core genes that are used by the "
                             "anvi'o taxonomy headquarters in Lausanne. Somewhat disappointing but totally OK.")

            # even if there are no SCGs to use for taxonomy later, we did attempt ot populate the
            # contigs database, so we shall note that in the self table to make sure the error from
            # `anvi-estimate-genome-taxonomy` is not "you seem to have not run taxonomy".
            self.tables_for_taxonomy.update_self_value()

            # return empty handed like a goose in the job market in 2020
            return None

        log_file_path = filesnpaths.get_temp_file_path()

        self.run.info('Taxonomy', self.ctx.accession_to_taxonomy_file_path)
        self.run.info('Database reference', self.ctx.search_databases_dir_path)
        self.run.info('Number of SCGs', len(scg_sequences_dict))

        self.run.warning('', header='Parameters for DIAMOND blastp', lc='green')
        self.run.info('Max number of target sequences', self.max_target_seqs)
        self.run.info('Max e-value to report alignments', self.evalue)
        self.run.info('Min percent identity to report alignments', self.min_pct_id)
        self.run.info('Num aligment tasks running in parallel', self.num_parallel_processes)
        self.run.info('Num CPUs per aligment task', self.num_threads)
        self.run.info('Log file path', log_file_path)

        self.tables_for_taxonomy.delete_contents_of_table(t.scg_taxonomy_table_name)
        self.tables_for_taxonomy.update_self_value(value=False)

        total_num_processes = len(scg_sequences_dict)

        self.progress.new('Computing SCGs aligments', progress_total_items=total_num_processes)
        self.progress.update('Initializing %d process...' % int(self.num_parallel_processes))

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        error_queue = manager.Queue()

        blastp_search_output = []

        for SCG in scg_sequences_dict:
            sequence = ""
            for entry in scg_sequences_dict[SCG].values():
                if 'sequence' not in entry or 'gene_name' not in entry:
                    raise ConfigError("The `get_filtered_dict` function got a parameter that "
                                      "does not look like the way we expected it. This function "
                                      "expects a dictionary that contains keys `gene_name` and `sequence`.")

                sequence = sequence + ">" + str(entry['gene_callers_id']) + "\n" + entry['sequence'] + "\n"
                entry['hits'] = []

            input_queue.put([SCG, sequence])

        workers = []
        for i in range(0, int(self.num_parallel_processes)):
            worker = multiprocessing.Process(target=self.blast_search_scgs_worker, args=(input_queue, output_queue, error_queue, log_file_path))

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

                if 'incompatible' in error_text:
                    raise ConfigError("Your current databases are incompatible with the diamond version you have on your computer. "
                                      "Please run the command `anvi-setup-scg-databases --redo-databases` and come back.")
                else:
                    raise ConfigError("Bad news. The database search operation failed somewhere :( It is very hard for anvi'o "
                                      "to know what happened, but the MOST LIKELY reason is that you have a diamond version "
                                      "installed on your system that is incompatible with anvi'o :/ The best course of action for that "
                                      "is to make sure running `diamond --version` on your terminal returns `0.9.14`. If not, "
                                      "try to upgrade/downgrade your diamond to match this version. If you are in a conda environmnet "
                                      "you can try running `conda install diamond=0.9.14`. Please feel free to contact us if the problem "
                                      "persists. We apologize for the inconvenience.")

            try:
                blastp_search_output += output_queue.get()

                if self.write_buffer_size > 0 and len(blastp_search_output) % self.write_buffer_size == 0:
                    self.tables_for_taxonomy.add(blastp_search_output)
                    blastp_search_output = []

                num_finished_processes += 1

                self.progress.increment(increment_to=num_finished_processes)
                self.progress.update("%s of %s SCGs are finished in %s processes with %s threads." \
                                        % (num_finished_processes, total_num_processes, int(self.num_parallel_processes), self.num_threads))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        # finally the remaining hits are written to the database, and we are done
        self.tables_for_taxonomy.add(blastp_search_output)

        # time to update the self table:
        self.tables_for_taxonomy.update_self_value()

        self.progress.end()


    def show_hits_gene_callers_id(self, gene_callers_id, scg_name, hits):
        self.progress.reset()
        self.run.warning(None, header='Hits for gene caller id %s' % gene_callers_id, lc="green")

        if len(hits):
            header = ['%id', 'bitscore', 'accession', 'taxonomy']
            table = []

            self.run.info_single("For '%s'" % scg_name, nl_before=1, nl_after=1)

            for hit in hits:
                table.append([str(hit['percent_identity']), str(hit['bitscore']), hit['accession'], ' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy])])

            print(tabulate(table, headers=header, tablefmt="fancy_grid", numalign="right"))
        else:
            self.run.info_single("No hits :/")


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


    def blast_search_scgs_worker(self, input_queue, output_queue, error_queue, log_file_path):
        """BLAST each SCG identified in the contigs database against the corresopinding
           target local database of GTDB seqeunces
        """

        while True:
            scg_name, fasta_formatted_scg_sequence = input_queue.get(True)
            target_database_path = self.ctx.SCGs[scg_name]['db']

            diamond = Diamond(target_database_path, run=run_quiet, progress=progress_quiet)
            diamond.max_target_seqs = self.max_target_seqs
            diamond.evalue = self.evalue
            diamond.min_pct_id = self.min_pct_id
            diamond.num_threads = self.num_threads
            diamond.run.log_file_path = log_file_path

            blastp_search_output = diamond.blastp_stdin_multi(fasta_formatted_scg_sequence)

            hits_per_gene = {}
            genes_estimation_output=[]

            for blastp_hit in blastp_search_output.split('\n'):
                if len(blastp_hit) and not blastp_hit.startswith('Query'):
                    fields = blastp_hit.split('\t')

                    try:
                        gene_callers_id = int(fields[0])
                        error_queue.put(None)
                    except:
                        error_queue.put(blastp_search_output)

                    hit = dict(zip(['accession', 'percent_identity', 'bitscore'], [fields[1], float(fields[2]), float(fields[11])]))
                    hit = self.update_dict_with_taxonomy(hit)

                    if gene_callers_id not in hits_per_gene:
                        hits_per_gene[gene_callers_id] = {}

                    if scg_name not in hits_per_gene[gene_callers_id]:
                        hits_per_gene[gene_callers_id][scg_name] = []

                    hits_per_gene[gene_callers_id][scg_name].append(hit)
                else:
                    error_queue.put(None)

            for gene_callers_id, scg_raw_hits in hits_per_gene.items():
                if len(scg_raw_hits.keys()) > 1:
                    self.run.warning("As crazy as it sounds, the gene callers id `%d` seems to have hit more than one SCG o_O Anvi'o will only use "
                                     "one of them almost absolutely randomly. Here are the SCGs the gene sequence matches: '%s'" % [s for s in scg_raw_hits.keys()])

                scg_name = list(scg_raw_hits.keys())[0]
                scg_raw_hits = scg_raw_hits[scg_name]

                scg_consensus_hit = self.get_consensus_hit(scg_raw_hits)
                scg_consensus_hit['accession'] = 'CONSENSUS'

                if anvio.DEBUG:
                    # avoid race conditions when priting this information when `--debug` is true:
                    with self.mutex:
                        self.progress.reset()
                        self.show_hits_gene_callers_id(gene_callers_id, scg_name, scg_raw_hits + [scg_consensus_hit])

                genes_estimation_output.append([gene_callers_id, scg_name, [scg_consensus_hit]])

            output_queue.put(genes_estimation_output)


    def get_consensus_hit(self, scg_raw_hits):
        pd.set_option('mode.chained_assignment', None)

        df = pd.DataFrame.from_records(scg_raw_hits)

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
