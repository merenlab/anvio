#!/usr/bin/env python3
# -*- coding: utf-8
"""
    This file contains SCGsSetup, SCGsDataBase, LowestIdentity classes.

"""

import re
import os
import sys
import gzip
import glob
import shutil
import pickle
import shutil
import subprocess
import multiprocessing

from collections import Counter

import anvio
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.pfam as pfam
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.errors import FilesNPathsError

from anvio.drivers.diamond import Diamond
from anvio.constants import default_scgs_taxonomy_data_dir, default_scgs_for_taxonomy


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
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

if sorted(list(locally_known_HMMs_to_remote_FASTAs.keys())) != sorted(default_scgs_for_taxonomy):
    raise ConfigError("Oh no. The SCGs designated to be used for all SCG taxonomy tasks in the constants.py\
                       are not the same names described in locally known HMMs to remote FASTA files\
                       conversion table definedd in SetUpSCGTaxonomyDatabase module. If this makes zero\
                       sense to you please ask a developer.")


class SetUpSCGTaxonomyDatabase:
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        # hard-coded GTDB variables
        self.target_database = "GTDB"
        self.target_database_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
        self.target_database_files = ['VERSION', 'ar122_msa_individual_genes.tar.gz', 'ar122_taxonomy.tsv',
                                    'bac120_msa_individual_genes.tar.gz', 'bac120_taxonomy.tsv']

        # user accessible variables
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.SCGs_taxonomy_data_dir = (os.path.abspath(A("scgs_taxonomy_data_dir")) if A("scgs_taxonomy_data_dir") else None) or (os.path.join(default_scgs_taxonomy_data_dir, self.target_database))
        self.reset = A("reset")
        self.target_database_URL = A("scgs_taxonomy_remote_database_url") or self.target_database_URL

        self.sanity_check()

        self.run.info("Local directory to setup", self.SCGs_taxonomy_data_dir)
        self.run.info("Reset the directory first", self.reset, mc="red")
        self.run.info("Remote database", self.target_database, nl_before=1, mc="green")
        self.run.info("Remote URL to download files", self.target_database_URL)
        self.run.info("Remote files of interest", ', '.join(self.target_database_files))


    def sanity_check(self):
        if not self.SCGs_taxonomy_data_dir:
            raise ConfigError("`SetUpSCGTaxonomyDatabase` class is upset because it was inherited without\
                               a directory for SCG taxonomy data to be stored :( This variable can't be None.")



    def setup(self):
        """This function downloads all GTDB files necessary to setup the SCG databases anvi'o will rely upon.

           In addition to downloading the original files, the setup will make sure everything, including the
           DIAMOND search databases are in place.
        """

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

        self.download_and_format_files()

        self.create_search_databases()


    def download_and_format_files(self):
        msa_individual_genes_dir_path = os.path.join(self.SCGs_taxonomy_data_dir, 'msa_individual_genes')
        accession_to_taxonomy_file_path = os.path.join(self.SCGs_taxonomy_data_dir, 'ACCESSION_TO_TAXONOMY.txt')

        # let's be 100% sure.
        os.remove(accession_to_taxonomy_file_path) if os.path.exists(accession_to_taxonomy_file_path) else None

        for remote_file_name in self.target_database_files:
            remote_file_url = '/'.join([self.target_database_URL, remote_file_name])
            local_file_path = os.path.join(self.SCGs_taxonomy_data_dir, remote_file_name)

            utils.download_file(remote_file_url, local_file_path, progress=self.progress, run=self.run)

            if local_file_path.endswith('individual_genes.tar.gz'):
                self.progress.new("Downloaded file patrol")
                self.progress.update("Unpacking file '%s'..." % os.path.basename(local_file_path))
                shutil.unpack_archive(local_file_path, extract_dir=msa_individual_genes_dir_path)
                os.remove(local_file_path)
                self.progress.end()

            if local_file_path.endswith('_taxonomy.tsv'):
                with open(accession_to_taxonomy_file_path, 'a') as f:
                    f.write(open(local_file_path).read())
                    os.remove(local_file_path)

        fasta_file_paths = glob.glob(msa_individual_genes_dir_path + '/*.faa')

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

        # FINALLY, checking whether downloaded FASTA files are suitable for the conversion
        msa_individual_gene_names_required = []
        [msa_individual_gene_names_required.extend(n) for n in locally_known_HMMs_to_remote_FASTAs.values()]

        msa_individual_gene_names_downloaded = [os.path.basename(f) for f in fasta_file_paths]

        missing_msa_gene_names = [n for n in msa_individual_gene_names_required if n not in msa_individual_gene_names_downloaded]
        if missing_msa_gene_names:
            raise ConfigError("Big trouble :( Anvi'o uses a hard-coded dictionary to convert locally known\
                               HMM models to FASTA files reported by GTDB project. It seems something has changed\
                               and %d of the FASTA files expected to be in the download directory are not there.\
                               Here is that list: '%s'. Someone needs to update the codebase by reading the\
                               appropriate documentation. If you are a user, you can't do much at this point but\
                               contacting the developers :( Anvi'o will keep the directory that contains all the\
                               downloaded files to update the conversion dictionary. Here is the full path to the\
                               output: %s" % (len(missing_msa_gene_names), ', '.join(missing_msa_gene_names), msa_individual_genes_dir_path))



    def create_search_databases(self):
        # FIXME: this function is next.
        pass


class SCGsDataBase():
    """Generate diamond databases of GTDB or import's Single Copy core genes Fasta whom matching with anvi'o HMM's"""

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None

        self.genes_files_directory = A('genes_files_directory')
        self.scgs_directory = A('scgs_directory')
        self.hmms = A('hmms')
        self.path_tsv_taxonomy = A('taxonomy_files')
        self.output_directory = A('out_put_directory')
        self.num_threads = A('num_threads')
        self.taxonomy_level = A('taxonomy_level')
        self.database_url = A('taxonomy_level') or "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"

        self.classic_input_directory = os.path.join(
            os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY/GTDB')

        if not self.path_tsv_taxonomy:
            self.path_tsv_taxonomy = os.path.join(
                self.classic_input_directory, 'ACCESSION_TO_TAXONOMY.txt')

        self.num_threads = args.num_threads

        if not self.num_threads:
            self.num_threads = 2

        if not self.genes_files_directory:
            self.genes_files_directory = os.path.join(
                self.classic_input_directory, 'msa_individual_genes')

        if not self.output_directory:
            self.output_directory = self.classic_input_directory

        self.tsv_output = os.path.join(
            self.output_directory, 'ACCESSION_TO_TAXONOMY.txt')

        self.SCGs_fasta = [files for files in os.listdir(
            self.scgs_directory) if not files.endswith(".dmnd")]  # FIXME useless?


        self.sanity_check()


    def sanity_check(self):
        if not filesnpaths.is_file_exists(self.genes_files_directory, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. If you did not provided any as a parameter \
                               anvi'o looks in '%s'. You can download file by using the commande 'anvi-setup-scgs'."
                                                        % (self.genes_files_directory, self.classic_input_directory))

        self.genes_files = os.listdir(self.genes_files_directory)

        if not filesnpaths.is_file_exists(self.path_tsv_taxonomy, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. If you did not provided any as a parameter \
                               anvi'o looks in '%s'. You can download file by using the commande 'anvi-setup-scgs'."
                                                            % (self.path_tsv_taxonomy, self.classic_input_directory))

        if not filesnpaths.is_file_exists(self.hmms, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. You must declare one before\
                               continue." % self.hmms)

        if filesnpaths.is_output_dir_writable(self.output_directory):
            filesnpaths.gen_output_directory(self.output_directory)



    def make_scg_db(self):
        self.run.info('SCG genes directory', self.genes_files_directory)
        self.run.info('Taxonomy file', self.path_tsv_taxonomy)
        self.run.info('Fasta for HMM reference', self.hmms)
        self.run.info('Output directory', self.output_directory)

        self.dictionary_correspondance_SCGs = {}

        self.matrix_taxonomy = self.do_taxonomy_dictonnrary_with_tsv(
            self.path_tsv_taxonomy)

        self.pathpickle_dictionary_correspondance_SCGs = os.path.join(
            self.output_directory, 'dictionary_correspondance_SCGs.pickle')

        if not os.path.exists(self.pathpickle_dictionary_correspondance_SCGs):
            self.make_correpondance_dictonnary_SCG()

        else:
            with open(self.pathpickle_dictionary_correspondance_SCGs, 'rb') as handle:
                self.dictionary_correspondance_SCGs = pickle.load(handle)

        self.dictionnary_pickel_taxo = {}
        if not self.taxonomy_level:
            name_out_put_final_SCG_db = "SCG_diamond_db"
        else:
            name_out_put_final_SCG_db = self.taxonomy_level+"SCG_diamond_db"
        for keycorres in self.dictionary_correspondance_SCGs.keys():

            self.output_directory_SCG = os.path.join(
                self.output_directory, name_out_put_final_SCG_db)

            self.path_new_fasta_SCG = os.path.join(
                self.output_directory_SCG, keycorres)

            self.path_original_genes = os.path.join(
                self.genes_files_directory, str(keycorres))

            filesnpaths.gen_output_directory(self.output_directory_SCG)

            self.create_SCG_db(keycorres, self.taxonomy_level)
            path_dictionnary_pickel_taxo_anvi = os.path.join(
                self.output_directory, 'dictionnary_msa_taxonomy_anvio.pickle')

            with open(path_dictionnary_pickel_taxo_anvi, 'wb') as handle:
                pickle.dump(self.dictionnary_pickel_taxo, handle,
                            protocol=pickle.HIGHEST_PROTOCOL)


    def create_SCG_db(self, keycorres, dictionnary_level=None):
        self.listtaxo = []
        unmatch = []
        newfasta = ""
        with open(self.tsv_output, 'a') as tsv:
            for refence in self.dictionary_correspondance_SCGs[keycorres]:
                self.path_new_fasta_SCG = os.path.join(
                    self.output_directory_SCG, keycorres+".faa")
                fasta = u.ReadFasta(os.path.join(
                    self.genes_files_directory, str(refence)+".faa"), quiet=True)
                if dictionnary_level != None:
                    for matrix_code, matrix_taxon in self.matrix_taxonomy.items():
                        index, name = self.match(
                            fasta, matrix_code, matrix_taxon, dictionnary_level)
                        if not index:
                            unmatch.append(name)
                            continue
                        else:
                            self.listtaxo.append(name)
                            print(
                                index, fasta.ids[index], self.matrix_taxonomy[fasta.ids[index]])
                            newfasta = newfasta+">" + \
                                fasta.ids[index]+"\n" + \
                                fasta.sequences[index]+"\n"
                            tsv.write(
                                fasta.ids[index]+"\t"+';'.join(self.matrix_taxonomy[fasta.ids[index]])+"\n")
                            self.dictionnary_pickel_taxo[fasta.ids[index]
                                                         ] = self.matrix_taxonomy[fasta.ids[index]]
                else:
                    i = 0
                    for fasta_id in fasta.ids:
                        if fasta_id in self.matrix_taxonomy:
                            newfasta = newfasta+">"+fasta_id + \
                                "\n"+fasta.sequences[i]+"\n"
                            if fasta_id not in self.dictionnary_pickel_taxo:
                                self.dictionnary_pickel_taxo[fasta_id] = self.matrix_taxonomy[fasta_id]
                                tsv.write(
                                    fasta_id+"\t"+';'.join(self.matrix_taxonomy[fasta_id])+"\n")
                        else:
                            unmatch.append(fasta_id)
                        i += 1
            with open(self.path_new_fasta_SCG, 'a') as file:
                file.write(newfasta)
            self.diamonddb_stdin(newfasta, self.path_new_fasta_SCG)


    def match(self, fasta, matrix_code, matrix_taxon, dictionnary_level):
        if dictionnary_level:
            name = str(matrix_taxon[self.taxonomy_level]).rstrip()
        else:
            name = str(matrix_taxon[-1]).rstrip()
        if name not in self.listtaxo and matrix_code in fasta.ids:
            index = fasta.ids.index(matrix_code)
            return(index, name)
        else:
            index = False
            name = (name+"\t"+matrix_taxon[0]+"\n")
            return(index, name)


    def make_correpondance_dictonnary_SCG(self):
        self.diamonddb(self.hmms)
        number_genes_files = len(self.genes_files)

        self.progress.new('Searching for corresponding SCG with Anvio reference', progress_total_items=number_genes_files)
        self.progress.update('...')

        genes_files_number = 0

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        for gene_file in self.genes_files:
            genes_files_number += 1

            self.name = gene_file.replace('.faa', '')
            self.path_diamondb = os.path.join(
                self.output_directory, "diamonddb")
            self.path_refund_genes = os.path.join(
                self.path_diamondb, gene_file)
            self.path_original_genes = os.path.join(
                self.genes_files_directory, gene_file)
            filesnpaths.gen_output_directory(self.path_diamondb)
            self.sequence_to_blast = self.refund_genes()

            #self.progress.update("Reference %s : %d / %d"% self.name, genes_files_number, number_genes_files)
            self.diamond_output = self.diamondblast_stdin(
                self.sequence_to_blast, self.hmms)

            if not self.diamond_output:
                continue

            self.hmm_gene = str(self.diamond_output).split('\n')[
                0].split('\t')[1]
            if self.hmm_gene not in self.dictionary_correspondance_SCGs:
                self.dictionary_correspondance_SCGs[self.hmm_gene] = [
                    self.name]

            else:
                self.dictionary_correspondance_SCGs[self.hmm_gene] += [
                    self.name]

            os.remove(self.path_refund_genes)
            os.remove(self.pathrefpath_refund_genesundgenes+'.dmnd')

        self.progress.end()
        if not len(self.dictionary_correspondance_SCGs):
            raise ConfigError(
                "There is no match with SCGs fasta provide and anvi'o ressource")

        self.run.info("Path correspondance dictonnary", path_tsv_taxonomy)

        with open(self.pathpickle_dictionary_correspondance_SCGs, 'wb') as handle:
            pickle.dump(self.dictionary_correspondance_SCGs,
                        handle, protocol=pickle.HIGHEST_PROTOCOL)


    def do_taxonomy_dictionary_with_tsv(self, path_tsv_taxonomy, filter=False):
        """Format TSV in dictonary with code as key and list for taxonomy for value"""

        self.run.info("Correspondance taxonomy file", path_tsv_taxonomy)

        self.matrix = {}
        individues_filterd = []
        number_individues_filterd = 0
        number_individues_selected = 0
        linestaxo = [line for line in open(path_tsv_taxonomy).readlines()]

        self.run.info("Total Indivues", " %d" % len(linestaxo))
        self.progress.new('Parse tsv file')
        self.progress.update('...')
        for line in linestaxo:
            taxo = []
            for taxonomy in line.split("\t")[1].split(";"):
                if filter:
                    level_taxo = taxonomy.rstrip()
                    level_taxo = re.sub(
                        r'(\\*[a-z])\_[A-Z]', r'\1', level_taxo)
                    level_taxo = re.sub(r"sp[0-9]*", "", level_taxo)
                    level_taxo = re.sub(" ", "_", level_taxo)
                    if re.match(r'[a-z]{3,}__[A-Z][a-z]*_[a-z]{3,}', level_taxo) or re.match(r'[a-z]__[A-Z][a-z]{3,}', level_taxo):
                        taxo.append(level_taxo)
                    if not len(taxo):
                        number_individues_filterd += 1
                        continue

                else:
                    level_taxo = re.sub(r'[a-z]__', '', taxonomy).rstrip()
                    level_taxo = re.sub(r' ', '_', level_taxo).rstrip()
                    taxo.append(level_taxo)

            number_individues_selected += 1
            self.progress.update('Number of indivues selected %s' %
                                 number_individues_selected)

            self.matrix[str(line.split("\t")[0])] = taxo

        self.progress.end()

        if number_individues_selected == 0:
            raise ConfigError("There's a problem with the file '%s'. Anvi'o couldn't use it. \
                              If you are sure that the file is not empty and in the right format try the option '--no-filter'."
                              % self.path_tsv_taxonomy)

        self.run.info("Indivues Selected", " %d" %
                      int(number_individues_selected))
        return(self.matrix)


    def diamonddb(self, pathgdb):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.query_fasta = pathgdb
        diamond.makedb(pathgdb)


    def diamondblast_stdin(self, sequencetoblast, path_reference):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.target_fasta = path_reference+".dmnd"
        diamond.blastp_stdin(sequencetoblast)


    def diamonddb_stdin(self, sequencetoblast, path_reference):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.makedb_stdin(sequencetoblast, output_file_path=path_reference)


    def diamondblast_stdou(self, path_fasta_db, max_target_seqs=20, evalue=1e-05, min_pct_id=90):
        pathdb = path_fasta_db + ".dmnd"
        diamond = Diamond(pathdb, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond_output = diamond.blastp_stdout(query_fasta=self.path_refund_genes,
                                               max_target_seqs=max_target_seqs,
                                               min_pct_id=min_pct_id)
        return(diamond_output)


    def refund_genes(self, percent_gap=0.5):
        """Use only the sequence with less than percent_gap"""

        fasta = u.ReadFasta(self.path_original_genes, quiet=True)

        i = 0
        sequence_to_blast = ""
        while i < len(fasta.ids):
            if fasta.sequences[i].count('-') < len(fasta.sequences[i])*percent_gap:
                sequence_to_blast = sequence_to_blast+">" + \
                    fasta.ids[i]+"\n"+fasta.sequences[i]+"\n"
            i += 1

        return(sequence_to_blast)


class lowident():
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        def A(x): return args.__dict__[x] if x in args.__dict__ else None

        self.scgs_directory = A('scgs_directory')
        self.path_tsv_taxonomy = A('taxonomy_files')
        self.output_file = A('output_file')
        self.num_threads = A('num_threads')
        self.num_process = A('num_process')

        self.classic_input_directory = os.path.join(
            os.path.dirname(anvio.__file__), 'data/misc/SCG_TAXONOMY/GTDB')

        if not self.path_tsv_taxonomy:
            self.path_tsv_taxonomy = os.path.join(
                self.classic_input_directory, 'ACCESSION_TO_TAXONOMY.txt')

        self.num_threads = args.num_threads

        if not self.num_threads:
            self.num_threads = 1

        if not self.num_process:
            self.num_process = 1

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
                    if leveltaxo.endwith("Bacteria"):
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

        num_listeprocess = len(listeprocess)

        progress.new('Aligning amino acid sequences for genes in gene clusters',
                     progress_total_items=num_listeprocess)
        progress.update('...')

        for pathquery in listeprocess:
            input_queue.put(pathquery)

        workers = []
        for i in range(0, int(self.num_process)):
            worker = multiprocessing.Process(target=self.diamond,
                                             args=(input_queue, output_queue))

            workers.append(worker)
            worker.start()

        finish_process = 0
        while finish_process < num_listeprocess:
            try:
                taxo_ident_item = output_queue.get()

                if taxo_ident_item:
                    dico_low_ident[taxo_ident_item['taxonomy']
                                   ] = taxo_ident_item['cutoff']

                finish_process += 1
                progress.increment()
                progress.update("Processed %d of %d non-singleton GCs in 10 threads." %
                                (finish_process, num_listeprocess))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        progress.end()

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
            pathpickle_dico_ident = path_new_fasta_SCG + "_dico_low_ident.pickle"

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
