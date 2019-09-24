#!/usr/bin/env python3
# -*- coding: utf-8
"""
Classes to setup remote SCG databases in local, use local databases to affiliate SCGs in anvi'o
contigs databases with taxon names, and estimate taxonomy for genomes and metagneomes.
"""

import os
import glob
import shutil
import pickle
import pandas as pd
import multiprocessing

from tabulate import tabulate
from collections import OrderedDict

import anvio
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.dbops import ContigsSuperclass
from anvio.drivers.diamond import Diamond
from anvio.tables.scgtaxonomy import TableForSCGTaxonomy


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)


_author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Quentin Clayssen"
__email__ = "quentin.clayssen@gmail.com"


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
    def __init__(self, args, skip_sanity_check=False):
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
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.SCGs_taxonomy_data_dir = (os.path.abspath(A("scgs_taxonomy_data_dir")) if A("scgs_taxonomy_data_dir") else None) or (os.path.join(self.default_scgs_taxonomy_data_dir, self.target_database))
        self.msa_individual_genes_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'MSA_OF_INDIVIDUAL_SCGs')
        self.accession_to_taxonomy_file_path = os.path.join(self.SCGs_taxonomy_data_dir, 'ACCESSION_TO_TAXONOMY.txt')
        self.search_databases_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'SCG_SEARCH_DATABASES')
        self.target_database_URL = A("scgs_taxonomy_remote_database_url") or self.target_database_URL

        # some dictionaries for convenience. we set them up here, but the proper place to sanity check
        # them may be somewhere else. for instance, when this class is inheritded by SetupLocalSCGTaxonomyData
        # the paths will not point to an actual file, but when it is inherited by PopulateContigsDatabaseWithSCGTaxonomy,
        # they better point to actual files.
        self.SCGs = dict([(SCG, {'db': os.path.join(self.search_databases_dir_path, SCG + '.dmnd'), 'fasta': os.path.join(self.search_databases_dir_path, SCG)}) for SCG in self.default_scgs_for_taxonomy])

        self.letter_to_level = dict([(l.split('_')[1][0], l) for l in self.levels_of_taxonomy])

        self.accession_to_taxonomy_dict = None
        if self.accession_to_taxonomy_file_path:
            self.progress.new("Reading the accession to taxonomy file")
            self.progress.update('...')

            self.accession_to_taxonomy_dict = {}
            with open(self.accession_to_taxonomy_file_path, 'r') as taxonomy_file:
                for accession, taxonomy_text in [l.strip('\n').split('\t') for l in taxonomy_file.readlines() if not l.startswith('#') and l]:
                    # taxonomy_text kinda looks like this:
                    #
                    #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter sp1
                    #
                    d = {}
                    for letter, taxon in [e.split('__', 1) for e in taxonomy_text.split(';')]:
                        if letter in self.letter_to_level:
                            # if 'species' level taxon name does not include an underscore, it usually is simply the
                            # assignment of the genus name as species name. which is not quite useful for our needs.
                            # so the following if statement removes those species level annotations.
                            if letter == 's' and '_' not in taxon:
                                d[self.letter_to_level[letter]] = None
                            else:
                                d[self.letter_to_level[letter]] = taxon
                        else:
                            self.run.warning("Some weird letter found in '%s'. " % taxonomy_text)

                    self.accession_to_taxonomy_dict[accession] = d
            
            # let's add one more accession for all those missing accessions
            self.accession_to_taxonomy_dict['unknown_accession'] = dict([(taxon, None) for taxon in self.levels_of_taxonomy])

            self.progress.end()

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        if sorted(list(locally_known_HMMs_to_remote_FASTAs.keys())) != sorted(self.default_scgs_for_taxonomy):
            raise ConfigError("Oh no. The SCGs designated to be used for all SCG taxonomy tasks in the constants.py\
                               are not the same names described in locally known HMMs to remote FASTA files\
                               conversion table definedd in SetupLocalSCGTaxonomyData module. If this makes zero\
                               sense to you please ask a developer.")

        if not self.SCGs_taxonomy_data_dir:
            raise ConfigError("`SetupLocalSCGTaxonomyData` class is upset because it was inherited without\
                               a directory for SCG taxonomy data to be stored :( This variable can't be None.")

        # sanity checks specific to each class
        if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy', 'SCGTaxonomyEstimator']:
            if not os.path.exists(self.SCGs_taxonomy_data_dir):
                raise ConfigError("Anvi'o could not find the data directory for the single-copy core genes taxonomy\
                                   setup. You may need to run `anvi-setup-scg-databases`, or provide a directory path\
                                   where SCG databases are set up. This is the current path anvi'o is considering (which\
                                   can be changed via the `--scgs-taxonomy-data-dir` parameter): '%s'" % (self.SCGs_taxonomy_data_dir))

            if not os.path.exists(self.accession_to_taxonomy_file_path):
                raise ConfigError("While your SCG taxonomy data dir seems to be in place, it is missing at least one critical\
                                   file (in this case, the file to resolve accession IDs to taxon names). You may need to run\
                                   the program `anvi-setup-scg-databases` with the `--reset` flag to set things right again.")

            if not self.contigs_db_path:
                raise ConfigError("For these things to work, you need to provide a contigs database for the anvi'o SCG\
                                   taxonomy workflow :(")

            utils.is_contigs_db(self.contigs_db_path)

            if self.__class__.__name__ in ['PopulateContigsDatabaseWithSCGTaxonomy']:
                missing_SCG_databases = [SCG for SCG in self.SCGs if not os.path.exists(self.SCGs[SCG]['db'])]
                if len(missing_SCG_databases):
                    raise ConfigError("Even though anvi'o found the directory for databases for taxonomy stuff,\
                                       your setup seems to be missing %d of %d databases required for everything to work\
                                       with the current genes configuration of this class. Here are the list of\
                                       genes for which we are missing databases: '%s'." % \
                                                (len(missing_SCG_databases), len(self.SCGs), ', '.join(missing_SCG_databases)))

            if self.__class__.__name__ in ['SCGTaxonomyEstimator']:
                scg_taxonomy_was_run = ContigsDatabase(self.contigs_db_path).meta['scg_taxonomy_was_run']
                if not scg_taxonomy_was_run:
                    raise ConfigError("It seems the SCG taxonomy tables were not populatd in this contigs database :/ Luckily it\
                                       is easy to fix that. Please see the program `anvi-run-scg-taxonomy`.")

                if self.profile_db_path:
                    utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

                if self.collection_name:
                    raise ConfigError("If you are asking anvi'o to estimate taxonomy using a collection, you must also provide\
                                       a profile database to this program.")

                if self.treat_as_metagenome and self.collection_name:
                    raise ConfigError("You can't ask anvi'o to treat your contigs database as a metagenome and also give it a\
                                       collection.")


    def update_dict_with_taxonomy(self, d, mode=None):
        """Takes a dictionary that includes a key `accession` and populates the dictionary with taxonomy"""

        if not mode:
            if not 'accession' in d:
                raise ConfigError("`add_taxonomy_to_dict` is speaking: the dictionary sent here does not have a member\
                                   with key `accession`.")

            if d['accession'] in self.accession_to_taxonomy_dict:
                d.update(self.accession_to_taxonomy_dict[d['accession']])
            else:
                d.update(self.accession_to_taxonomy_dict['unknown_accession'])

        elif mode == 'list_of_dicts':
            if len([entry for entry in d if 'accession' not in entry]):
                raise ConfigError("`add_taxonomy_to_dict` is speaking: you have a bad formatted data here :/")

            for entry in d:
                print(self.taxonomy_dict[entry['accession']])

        else:
            raise ConfigError("An unknown mode (%s) is setn to `add_taxonomy_to_dict` :/" % (mode))

        return d



        if not self.SCGs_taxonomy_data_dir:
            raise ConfigError("`SetupLocalSCGTaxonomyData` class is upset because it was inherited without\
                               a directory for SCG taxonomy data to be stored :( This variable can't be None.")


class SetupLocalSCGTaxonomyData(SCGTaxonomyContext):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # user accessible variables
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.reset = A("reset")
        self.num_threads = A('num_threads')

        SCGTaxonomyContext.__init__(self, self.args)

        self.run.info("Local directory to setup", self.SCGs_taxonomy_data_dir)
        self.run.info("Reset the directory first", self.reset, mc="red")
        self.run.info("Remote database", self.target_database, nl_before=1, mc="green")
        self.run.info("Remote URL to download files", self.target_database_URL)
        self.run.info("Remote files of interest", ', '.join(self.target_database_files))


    def setup(self):
        """This function downloads all GTDB files necessary to setup the SCG databases anvi'o will rely upon.

           In addition to downloading the original files, the setup will make sure everything, including the
           DIAMOND search databases are in place.
        """

        if not anvio.DEBUG:
            self.check_initial_directory_structure()

        self.run.warning("Please remember that the data anvi'o attempts do download on behalf of you are\
                          courtesy of The Genome Taxonomy Database (GTDB), an initiative to establish a \
                          standardised microbial taxonomy based on genome phylogeny, primarly funded by\
                          tax payers in Australia. Please don't forget to cite the original work,\
                          doi:10.1038/nbt.4229 by Parks et al to explicitly mention the source of databases\
                          anvi'o will use to estimate genome level taxonomy. If you are not sure how it\
                          should look like in your methods sections, anvi'o developers will be happy to\
                          help you if you can't find any published example to get inspiration.", lc = 'yellow')

        self.progress.new("%s setup" % self.target_database)

        self.progress.update("Reading the VERSION file...")
        content = utils.get_remote_file_content(self.target_database_URL + 'VERSION')
        version, release_date  = content.strip().split('\n')[0].strip(), content.strip().split('\n')[2].strip()
        self.progress.end()

        self.run.info("%s release found" % self.target_database, "%s (%s)" % (version, release_date), mc="green")

        # FIXME
        if not anvio.DEBUG:
            self.download_and_format_files()

        self.create_search_databases()

        if not anvio.DEBUG:
            self.clean_up()


    def check_initial_directory_structure(self):
        if os.path.exists(self.SCGs_taxonomy_data_dir):
            if self.reset:
                shutil.rmtree(self.SCGs_taxonomy_data_dir)
                self.run.warning('The existing directory for SCG taxonomy data dir has been removed. Just so you know.')
                filesnpaths.gen_output_directory(self.SCGs_taxonomy_data_dir)

            else:
                raise ConfigError("You already seem to have a directory where anvi'o intends to use for setup. If you wish to\
                                   re-run the setup, please use the flag `--reset` and BE VERY CAREFUL that this\
                                   directory does not contain anything you don't want to lose: '%s'." % self.SCGs_taxonomy_data_dir)
        else:
            filesnpaths.gen_output_directory(self.SCGs_taxonomy_data_dir)


    def download_and_format_files(self):
        # let's be 100% sure.
        os.remove(self.accession_to_taxonomy_file_path) if os.path.exists(self.accession_to_taxonomy_file_path) else None

        for remote_file_name in self.target_database_files:
            remote_file_url = '/'.join([self.target_database_URL, remote_file_name])
            local_file_path = os.path.join(self.SCGs_taxonomy_data_dir, remote_file_name)

            utils.download_file(remote_file_url, local_file_path, progress=self.progress, run=self.run)

            if local_file_path.endswith('individual_genes.tar.gz'):
                self.progress.new("Downloaded file patrol")
                self.progress.update("Unpacking file '%s'..." % os.path.basename(local_file_path))
                shutil.unpack_archive(local_file_path, extract_dir=self.msa_individual_genes_dir_path)
                os.remove(local_file_path)
                self.progress.end()

            if local_file_path.endswith('_taxonomy.tsv'):
                with open(self.accession_to_taxonomy_file_path, 'a') as f:
                    f.write(open(local_file_path).read())
                    os.remove(local_file_path)

        fasta_file_paths = glob.glob(self.msa_individual_genes_dir_path + '/*.faa')

        if not fasta_file_paths:
            raise ConfigError("Something weird happened while anvi'o was trying to take care of the files\
                               it downloaded from GTDB. Please let a developer know about this unless it is\
                               not already reported in our issue tracker at Github :(")

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

        self.progress.end()


    def create_search_databases(self):
        self.progress.new("Creating search databases")
        self.progress.update("Checking output directory to store files ...")
        filesnpaths.is_output_dir_writable(os.path.dirname(self.search_databases_dir_path))
        filesnpaths.gen_output_directory(self.search_databases_dir_path, delete_if_exists=True, dont_warn=True)

        # We will be working with the files downloaded in whatever directory before. The first thing is to check
        # whether whether FASTA files in the directory are suitable for the conversion
        self.progress.update("Checking the conversion dict and FASTA files ...")
        msa_individual_gene_names_required = []
        [msa_individual_gene_names_required.extend(n) for n in locally_known_HMMs_to_remote_FASTAs.values()]

        fasta_file_paths = glob.glob(self.msa_individual_genes_dir_path + '/*.faa')
        msa_individual_gene_names_downloaded = [os.path.basename(f) for f in fasta_file_paths]

        missing_msa_gene_names = [n for n in msa_individual_gene_names_required if n not in msa_individual_gene_names_downloaded]
        if missing_msa_gene_names:
            self.progress.reset()
            raise ConfigError("Big trouble :( Anvi'o uses a hard-coded dictionary to convert locally known\
                               HMM models to FASTA files reported by GTDB project. It seems something has changed\
                               and %d of the FASTA files expected to be in the download directory are not there.\
                               Here is that list: '%s'. Someone needs to update the codebase by reading the\
                               appropriate documentation. If you are a user, you can't do much at this point but\
                               contacting the developers :( Anvi'o will keep the directory that contains all the\
                               downloaded files to update the conversion dictionary. Here is the full path to the\
                               output: %s" % (len(missing_msa_gene_names), ', '.join(missing_msa_gene_names), self.msa_individual_genes_dir_path))
        else:
            self.progress.reset()
            self.run.info_single("Good news! Conversion dict and FASTA files it requires seem to be all in place.\
                                  Anvi'o will now proceed to merge %d FASTA files that correspond to %d SCGs, and\
                                  create individual search databases for them." % \
                                        (len(msa_individual_gene_names_required), len(locally_known_HMMs_to_remote_FASTAs)), nl_before=1, nl_after=1, mc="green")

        # Merge FASTA files that should be merged. This is defined in the conversion dictionary.
        for SCG in locally_known_HMMs_to_remote_FASTAs:
            self.progress.update("Working on %s in %d threads" % (SCG, self.num_threads))

            files_to_concatenate = [os.path.join(self.msa_individual_genes_dir_path, f) for f in locally_known_HMMs_to_remote_FASTAs[SCG]]
            destination_file_path = os.path.join(self.search_databases_dir_path, SCG)

            # concatenate from the dictionary into the new destination
            utils.concatenate_files(destination_file_path, files_to_concatenate)

            # create a diamond search database for `destination_file_path`
            diamond = Diamond(query_fasta=destination_file_path, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
            diamond.makedb(output_file_path=destination_file_path + ".dmnd")

        self.progress.reset()
        self.run.info_single("Another good news. All FASTA files that were supposed to be merged were merged and\
                              turned into fancy search databases.", nl_after=1, mc="green")


    def clean_up(self):
        for file_path in [os.path.join(self.SCGs_taxonomy_data_dir, 'diamond-log-file.txt')]:
            if os.path.exists(file_path):
                os.remove(file_path)

        for dir_path in [self.msa_individual_genes_dir_path]:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)


class PopulateContigsDatabaseWithSCGTaxonomy(SCGTaxonomyContext):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        SCGTaxonomyContext.__init__(self, self.args)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.write_buffer_size = int(A('write_buffer_size') if A('write_buffer_size') is not None else 1000)
        self.contigs_db_path = A('contigs_db')
        self.num_parallel_processes = int(A('num_parallel_processes')) if A('num_parallel_processes') else 1
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        self.max_target_seqs = 20
        self.evalue = 1e-05
        self.min_pct_id = 90

        self.taxonomy_dict = OrderedDict()

        self.mutex = multiprocessing.Lock()


    def get_SCG_sequences_dict_from_contigs_db(self):
        """Returns a dictionary of all HMM hits per SCG of interest"""

        contigs_db = ContigsSuperclass(self.args, r=self.run, p=self.progress)
        splits_dict = {contigs_db.a_meta['project_name']: list(contigs_db.splits_basic_info.keys())} 

        s = hmmops.SequencesForHMMHits(self.args.contigs_db, sources=self.hmm_source_for_scg_taxonomy)
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits(splits_dict, return_amino_acid_sequences=True)
        hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set(self.default_scgs_for_taxonomy))

        if not len(hmm_sequences_dict):
            raise ConfigError("Your selections returned an empty list of genes to work with :/")

        self.run.info('Hits', '%d hits for %d source(s)' % (len(hmm_sequences_dict), len(s.sources)))

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

        # this is the dictionary that shows all hits for each SCG of interest
        scg_sequences_dict = self.get_SCG_sequences_dict_from_contigs_db()

        self.run.info('Taxonomy', self.accession_to_taxonomy_file_path)
        self.run.info('Database reference', self.search_databases_dir_path)
        self.run.info('Number of SCGs', len(scg_sequences_dict))

        self.run.warning('', header='Parameters for DIAMOND search', lc='green')
        self.run.info('Blast type', "Blastp")
        self.run.info('Maximum number of target sequences', self.max_target_seqs)
        self.run.info('Minimum bit score to report alignments', self.min_pct_id)
        self.run.info('Number aligment running at same time', self.num_parallel_processes)
        self.run.info('Number of CPUs will be used for each aligment', self.num_threads)

        self.tables_for_taxonomy = TableForSCGTaxonomy(self.contigs_db_path, self.run, self.progress)
        self.tables_for_taxonomy.delete_contents_of_table(t.scg_taxonomy_table_name)

        total_num_processes = len(scg_sequences_dict)

        self.progress.new('Computing SCGs aligments', progress_total_items=total_num_processes)
        self.progress.update('Initializing %d process...' % int(self.num_parallel_processes))

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        blastp_search_output = []

        for SCG in scg_sequences_dict:
            sequence = ""
            for entry in scg_sequences_dict[SCG].values():
                if 'sequence' not in entry or 'gene_name' not in entry:
                    raise ConfigError("The `get_filtered_dict` function got a parameter that\
                                       does not look like the way we expected it. This function\
                                       expects a dictionary that contains keys `gene_name` and `sequence`.")

                sequence = sequence + ">" + str(entry['gene_callers_id']) + "\n" + entry['sequence'] + "\n"
                entry['hits'] = []

            input_queue.put([SCG, sequence])

        workers = []
        for i in range(0, int(self.num_parallel_processes)):
            worker = multiprocessing.Process(target=self.blast_search_scgs_worker, args=(input_queue, output_queue))

            workers.append(worker)
            worker.start()

        num_finished_processes = 0
        while num_finished_processes < total_num_processes:
            try:
                blastp_search_output += output_queue.get()

                if self.write_buffer_size > 0 and len(blastp_search_output) % self.write_buffer_size == 0:
                    self.tables_for_taxonomy.add(blastp_search_output)
                    blastp_search_output = []

                num_finished_processes += 1

                self.progress.increment(increment_to=num_finished_processes)
                self.progress.update("Processed %s of %s SGCs aligment in %s processus with %s cores." \
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
                table.append([hit['percent_identity'], hit['bitscore'], hit['accession'], ' / '.join([hit[l] if hit[l] else '' for l in self.levels_of_taxonomy])])

            print(tabulate(table, headers=header, tablefmt="fancy_grid", numalign="right"))
        else:
            self.run.info_single("No hits :/")


    def blast_search_scgs_worker(self, input_queue, output_queue):
        """BLAST each SCG identified in the contigs database against the corresopinding
           target local database of GTDB seqeunces
        """

        while True:
            scg_name, fasta_formatted_scg_sequence = input_queue.get(True)
            target_database_path = self.SCGs[scg_name]['db']

            diamond = Diamond(target_database_path, run=run_quiet, progress=progress_quiet)
            diamond.max_target_seqs = self.max_target_seqs
            diamond.evalue = self.evalue
            diamond.min_pct_id = self.min_pct_id
            diamond.num_threads = self.num_threads

            blastp_search_output = diamond.blastp_stdin_multi(fasta_formatted_scg_sequence)

            hits_per_gene = {}
            genes_estimation_output=[]

            for blastp_hit in blastp_search_output.split('\n'):
                if len(blastp_hit) and not blastp_hit.startswith('Query'):
                    fields = blastp_hit.split('\t')
                    gene_callers_id = int(fields[0])
                    hit = dict(zip(['accession', 'percent_identity', 'bitscore'], [fields[1], float(fields[2]), float(fields[11])]))
                    hit = self.update_dict_with_taxonomy(hit)

                    if gene_callers_id not in hits_per_gene:
                        hits_per_gene[gene_callers_id] = {}

                    if scg_name not in hits_per_gene[gene_callers_id]:
                        hits_per_gene[gene_callers_id][scg_name] = []

                    hits_per_gene[gene_callers_id][scg_name].append(hit)

            for gene_callers_id, scg_raw_hits in hits_per_gene.items():
                if len(scg_raw_hits.keys()) > 1:
                    self.run.warning("As crazy as it sounds, the gene callers id `%d` seems to have hit more than one SCG o_O Anvi'o will only use\
                                      one of them almost absolutely randomly. Here are the SCGs the gene sequence matches: '%s'" % [s for s in scg_raw_hits.keys()])

                scg_name = list(scg_raw_hits.keys())[0]
                scg_raw_hits = scg_raw_hits[scg_name]

                scg_consensus_hit = self.get_consensus_hit(scg_raw_hits)
                scg_consensus_hit['accession'] = 'CONSENSUS'

                scg_all_hits = scg_raw_hits + [scg_consensus_hit]

                if anvio.DEBUG:
                    # avoid race conditions when priting this information when `--debug` is true:
                    with self.mutex:
                        self.progress.reset()
                        self.show_hits_gene_callers_id(gene_callers_id, scg_name, scg_raw_hits)

                genes_estimation_output.append([gene_callers_id, scg_name, scg_all_hits])

            output_queue.put(genes_estimation_output)


    def get_consensus_hit(self, scg_raw_hits):
        df = pd.DataFrame.from_records(scg_raw_hits)

        # all this does at this point is to find the top hit based on percent identity and return it.
        # there are so many things that must be done here.
        hit_with_max_percent_identity = (df.loc[df['percent_identity'].idxmax()])

        return dict(hit_with_max_percent_identity)


class lowident():
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None

        self.scgs_directory = A('scgs_directory')
        self.path_tsv_taxonomy = A('taxonomy_files')
        self.output_file = A('output_file')
        self.num_threads = A('num_threads')
        self.num_parallel_processes = A('num_parallel_processes')

        self.classic_input_directory = os.path.join(
            os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY/GTDB')

        if not self.path_tsv_taxonomy:
            self.path_tsv_taxonomy = os.path.join(
                self.classic_input_directory, 'ACCESSION_TO_TAXONOMY.txt')

        self.num_threads = args.num_threads

        if not self.num_threads:
            self.num_threads = 1

        if not self.num_parallel_processes:
            self.num_parallel_processes = 1

        if not self.scgs_directory:
            self.scgs_directory = os.path.join(
                self.classic_input_directory, 'SCGs/')

        if not self.output_file:
            self.output_file = os.path.join(self.classic_input_directory,
             'MIN_PCT_ID_PER_TAXONOMIC_LEVEL.pickle')

        self.genes_files = [files for files in os.listdir(
            self.scgs_directory) if not files.endswith(".dmnd")]

        self.dicolevel = {}


    def process(self):
        self.output_directory = filesnpaths.get_temp_directory_path()
        self.pathpickle_dico_taxo = os.path.join(
            self.output_directory, 'dico_taxo_code_species.pickle')

        if not os.path.exists(self.pathpickle_dico_taxo):
            self.make_dicolevel(self.path_tsv_taxonomy)
            with open(self.pathpickle_dico_taxo, 'wb') as handle:
                pickle.dump(self.dicolevel, handle,
                            protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(self.pathpickle_dico_taxo, 'rb') as handle:
                self.dicolevel = pickle.load(handle)

        dico_low_ident = {}
        for genes in self.genes_files:
            dico_low_ident_genes = {}
            path_fasta = os.path.join(self.scgs_directory, genes)
            dico_low_ident_genes = self.creatsubfa_ident(
                path_fasta, genes, dico_low_ident_genes)
            pathpickle_dico_ident = os.path.join(
                self.output_directory, genes+'_dico_low_ident.pickle')
            with open(pathpickle_dico_ident, 'wb') as handle:
                pickle.dump(dico_low_ident_genes, handle,
                            protocol=pickle.HIGHEST_PROTOCOL)
            if genes not in dico_low_ident:
                dico_low_ident[genes] = dico_low_ident_genes
            else:
                dico_low_ident[genes] = dico_low_ident[genes].update(
                    dico_low_ident_genes)

        pathpickle_dico_ident = os.path.join(
            self.output_file)
        with open(self.output_file, 'wb') as handle:
            pickle.dump(dico_low_ident, handle,
                        protocol=pickle.HIGHEST_PROTOCOL)


    def make_dicolevel(self, path_tsv_taxonomy):
        with open(path_tsv_taxonomy, 'r') as taxotsv:
            linestaxo = taxotsv.readlines()
            for line in linestaxo:
                names = line.split('\t')
                code = names[0]
                leveltaxos = names[1].split(';')
                for leveltaxo in leveltaxos:
                    if leveltaxo.endswith("Bacteria"):
                        continue
                    leveltaxo = leveltaxo.rstrip()
                    if leveltaxo not in self.dicolevel:
                        self.dicolevel[leveltaxo] = [code]
                        continue
                    if code in self.dicolevel[leveltaxo]:
                        continue
                    else:
                        self.dicolevel[leveltaxo] = self.dicolevel[leveltaxo] + [code]


    def match_ident(self, fasta, codes, listindex):
        for code in codes:
            if code in fasta.ids:
                index = fasta.ids.index(code)
                if index:
                    listindex.append(index)
        return(listindex)


    def multidiamond(self, listeprocess, dico_low_ident):
        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        total_num_processes = len(listeprocess)

        self.progress.new('Aligning amino acid sequences for genes in gene clusters',
                     progress_total_items=total_num_processes)
        self.progress.update('...')

        for pathquery in listeprocess:
            input_queue.put(pathquery)

        workers = []
        for i in range(0, int(self.num_parallel_processes)):
            worker = multiprocessing.Process(target=self.diamond, args=(input_queue, output_queue))

            workers.append(worker)
            worker.start()

        num_finished_processes = 0
        while num_finished_processes < total_num_processes:
            try:
                taxo_ident_item = output_queue.get()

                if taxo_ident_item:
                    dico_low_ident[taxo_ident_item['taxonomy']] = taxo_ident_item['cutoff']

                num_finished_processes += 1
                self.progress.increment()
                self.progress.update("Processed %d of %d non-singleton GCs in 10 threads." % (num_finished_processes, total_num_processes))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        self.progress.end()

        return dico_low_ident


    def creatsubfa_ident(self, path_fasta, genes, dico_low_ident):

        fasta = u.ReadFasta(path_fasta, quiet=True)
        listeprocess = []
        for taxonomy, codes in self.dicolevel.items():
            if taxonomy.startswith("Archaea") or taxonomy.startswith("Bacteria"):
                continue
            listindex = []
            riboname = genes.replace(".faa", "")
            path_new_fasta_SCG = os.path.join(self.output_directory, taxonomy)

            if not os.path.exists(path_new_fasta_SCG):
                listindex = self.match_ident(fasta, codes, listindex)
                if len(listindex) > 1:
                    listeprocess.append(path_new_fasta_SCG)
                    with open(path_new_fasta_SCG, 'w+') as file:
                        for index in listindex:
                            file.write(">" + fasta.ids[index] +
                                       "\n" + fasta.sequences[index] + "\n")
                else:
                    if genes not in dico_low_ident:
                        dico_low_ident[riboname] = {}
                        dico_low_ident[riboname][taxonomy] = 100
                    else:
                        dico_low_ident[riboname][taxonomy] = 100

        if listeprocess:
            dico_low_ident = self.multidiamond(listeprocess, dico_low_ident)
            return(dico_low_ident)
        else:
            dico_low_ident = {}
            return(dico_low_ident)


    def diamond(self, input_queue, output_queue):
        while True:
            pathquery = input_queue.get(True)
            pathdb = pathquery+".dmnd"
            path_diamond = pathquery+'.txt'
            taxonomy = os.path.basename(pathquery)

            if not os.path.exists(path_diamond):
                self.diamonddb(pathquery, pathdb)
                ouputdiamond = self.run_diamond(pathquery, pathdb)

                os.remove(pathdb)
                os.remove(pathquery + 'log_file')

            low_ident = self.select_low_ident(ouputdiamond)
            os.remove(pathquery)
            output = {'taxonomy': taxonomy, 'cutoff': low_ident}
            output_queue.put(output)


    def diamonddb(self, pathquery, pathdb):
        diamond = Diamond(query_fasta=pathquery, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)

        diamond.query_fasta = pathquery
        diamond.run.log_file_path = pathquery + 'log_file'
        diamond.target_fasta = pathquery
        diamond.num_threads = self.num_threads

        diamond.makedb()

    def run_diamond(self, pathquery, pathdb):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)

        diamond.evalue = None
        diamond.run.log_file_path = pathquery+'log_file'
        diamond.target_fasta = pathdb
        diamond.query_fasta = pathquery
        diamond.max_target_seqs = None
        diamond.search_output_path = pathquery
        diamond.tabular_output_path = pathquery + '.txt'

        output = diamond.blastp_stdout()

        return output


    def select_low_ident(self, str_diamond_output, lowest_ident=100):
        "Select the lowest percent identity on aligment output"

        for line in str_diamond_output.split('\n'):
            if line:
                ident = line.strip().split('\t')[2]
                if float(ident) < float(lowest_ident):
                    lowest_ident = ident

        return(lowest_ident)
