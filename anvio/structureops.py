# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

import os
import shutil

import pandas as pd
import anvio.dbops as dbops
import anvio.fastalib as u
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.drivers.MODELLER as MODELLER

from anvio.errors import ConfigError, FilesNPathsError

run = terminal.Run()
progress = terminal.Progress()

class Structure:

    def __init__(self, args):
        self.args = args

        # initialize self.arg parameters
        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.genes_of_interest_path = A('genes_of_interest', null)
        self.splits_of_interest_path = A('splits_of_interest', null)
        self.bin_id = A('bin_id', null)
        self.collection_name = A('collection_name', null)
        self.gene_caller_ids = A('gene_caller_ids', null)
        self.output_dir = A('output_dir', null)
        self.full_output = A('black_no_sugar', bool)

        # MODELLER params
        self.modeller_database = A('database', null)
        self.best = A('best', null)
        self.max_matches = A('max_number_templates', null)
        self.min_proper_pident = A('percent_identical_cutoff', null)
        self.num_models = A('num_models', null)
        self.deviation = A('deviation', null)
        self.very_fast = A('very_fast', bool)

        # check output and define absolute path
        self.output_dir = filesnpaths.check_output_directory(self.output_dir, ok_if_exists=False)

        # identify which genes user wants to model structures for
        self.get_genes_of_interest()


    def get_genes_of_interest(self):
        """
        nabs the genes of interest based on user arguments (self.args)
        """
        self.genes_of_interest = None

        # identify the gene caller ids of all genes available
        self.genes_in_database = set(dbops.ContigsSuperclass(self.args).genes_in_splits.keys())

        if not self.genes_in_database:
            raise ConfigError("This contigs database does not contain any identified genes...")

        # settling genes of interest
        if self.genes_of_interest_path and self.gene_caller_ids:
            raise ConfigError("You can't provide a gene caller id from the command line, and a list of gene caller ids\
                               as a file at the same time, obviously.")

        if self.gene_caller_ids:
            self.gene_caller_ids = set([x.strip() for x in self.gene_caller_ids.split(',')])

            self.genes_of_interest = []
            for gene in self.gene_caller_ids:
                try:
                    self.genes_of_interest.append(int(gene))
                except:
                    raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % str(gene))

            self.genes_of_interest = set(self.genes_of_interest)

        elif self.genes_of_interest_path:
            filesnpaths.is_file_tab_delimited(self.genes_of_interest_path, expected_number_of_fields=1)

            try:
                self.genes_of_interest = set([int(s.strip()) for s in open(self.genes_of_interest_path).readlines()])
            except ValueError:
                raise ConfigError("Well. Anvi'o was working on your genes of interest ... and ... those gene IDs did not\
                                   look like anvi'o gene caller ids :/ Anvi'o is now sad.")

        if not self.genes_of_interest:
            # no genes of interest are specified. Assuming all, which could be innumerable--raise warning
            self.genes_of_interest = self.genes_in_database
            run.warning("You did not specify any genes of interest, so anvi'o will assume all of them are of interest.")


        # check for genes that do not appear in the contigs database
        bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in self.genes_in_database]
        if bad_gene_caller_ids:
            raise ConfigError(("This gene caller id you provided is" if len(bad_gene_caller_ids) == 1 else \
                               "These gene caller ids you provided are") + " not known to this contigs database: {}.\
                               You have only 2 lives left. 2 more mistakes, and anvi'o will automatically uninstall \
                               itself. Yes, seriously :(".format(", ".join([str(x) for x in bad_gene_caller_ids])))

        # Finally, raise warning if number of genes is greater than 20
        if len(self.genes_of_interest) > 20:
            import time
            run.warning("Modelling protein structures is no joke. The number of genes you want protein structures for is \
                         {}, which is a lot (of time!). I'm putting you in timeout for 15 seconds, then I'm going to do \
                         what you said to do. CTRL + C to cancel.".format(len(self.genes_of_interest)))
            time.sleep(15)


    def pick_best_model(self):
        """
        MODELLER has modelled some models. Assuming not
        """
        pass


    def download_structures(self):
        """
        Downloads structure files for self.top_seq_seq_matches using Biopython
        If the 4-letter code is `wxyz`, the downloaded file is `pdbwxyz.ent`.
        """
        progress.new("Downloading homologs from PDB")

        # define directory path name to store the template PDBs (it can already exist)
        self.template_pdbs = os.path.join(self.modeller.directory, "{}_TEMPLATE_PDBS".format(self.modeller.gene_id))

        downloaded = utils.download_protein_structures([code[0] for code in self.top_seq_matches], self.template_pdbs)

        # redefine self.top_seq_matches in case not all were downloaded
        self.top_seq_matches = [(code, chain_code) for code, chain_code in self.top_seq_matches if code in downloaded]

        if not len(self.top_seq_matches):
            run.warning("No structures of the homologous proteins (templates) were downloadable. Probably something \
                         is wrong. Maybe you are not connected to the internet. Stopping here.")
            raise self.EndModeller

        progress.end()
        run.info("structures downloaded for", ", ".join([code[0] for code in self.top_seq_matches]))


    def parse_search_results(self):
        """
        Parses search results and filters for best homologs to use as structure templates.

        parameters used :: self.min_proper_pident, self.max_matches
        """
        if not self.min_proper_pident or not self.max_matches:
            raise ConfigError("parse_search_results::You initiated this class without providing values for min_proper_pident \
                               and max_matches, which is required for this function.")

        progress.new("PARSE AND FILTER HOMOLOGS")
        progress.update("Finding those with percent identicalness > {}%".format(self.min_proper_pident))

        # put names to the columns
        column_names = (      "idx"      ,  "code_and_chain"  ,       "type"     ,  "iteration_num"  ,
                            "seq_len"    , "start_pos_target" , "end_pos_target" , "start_pos_dbseq" ,
                        "send_pos_dbseq" ,   "align_length"   ,   "seq_identity" ,      "evalue"      )

        # load the table as a dataframe
        search_df = pd.read_csv(self.modeller.search_results_path, sep="\s+", comment="#", names=column_names, index_col=False)

        # matches found (-1 because target is included)
        matches_found = len(search_df) - 1

        if not matches_found:
            progress.end()
            run.warning("No proteins with homologous sequence were found for {}. No structure will be modelled".\
                        format(self.modeller.gene_id))
            raise self.modeller.EndModeller

        # add some useful columns
        search_df["proper_pident"] = search_df["seq_identity"] * search_df["align_length"] / \
                                     search_df.iloc[0, search_df.columns.get_loc("seq_len")]
        search_df["code"] = search_df["code_and_chain"].str[:-1]
        search_df["chain"] = search_df["code_and_chain"].str[-1]

        # filter results by self.min_proper_pident.
        max_pident_found = search_df["proper_pident"].max()
        search_df = search_df[search_df["proper_pident"] >= self.min_proper_pident]

        # Order them and take the first self.modeller.max_matches.
        matches_after_filter = len(search_df)
        if not matches_after_filter:
            progress.end()
            run.warning("Gene {} did not have a search result with percent identicalness above or equal \
                         to {}. The max found was {}%. No structure will be modelled.".\
                         format(self.modeller.gene_id, self.min_proper_pident, max_pident_found))
            raise self.modeller.EndModeller

        progress.update("Keeping top {} matches as the template homologs".format(self.max_matches))

        # of those filtered, get up to self.modeller.max_matches of those with the highest proper_ident scores.
        search_df = search_df.sort_values("proper_pident", ascending=False)
        search_df = search_df.iloc[:min([len(search_df), self.max_matches])]

        # Get their chain and 4-letter ids
        self.top_seq_matches = list(zip(search_df["code"], search_df["chain"]))

        progress.end()
        run.info("Max number of templates allowed", self.max_matches)
        run.info("Number of candidate templates", matches_found)
        run.info("After >{}% identical filter".format(self.min_proper_pident), matches_after_filter)
        run.info("Number accepted as templates", len(self.top_seq_matches))
        for i in range(len(self.top_seq_matches)):
            run.info("Template {}".format(i+1), 
                     "Protein ID: {}, Chain {} ({:.1f}% identical)".format(self.top_seq_matches[i][0],
                                                                           self.top_seq_matches[i][1],
                                                                           search_df["proper_pident"].iloc[i]))


    def process(self):
        """
        This is the workflow for a standard protein search for homologous templates and then
        modelling the structure of the target protein using the homologous protein structures.
        """

        """
        sqlite-migration branch has a parameter passed to dbops.export_aa_sequences_from_contigs_db
        that lets you pass genes of interest. When these branches are merged, the code will look
        like this:

        for gene in self.genes_of_interest:

            # MODELLER outputs a lot of stuff into its working directory. A temporary directory is made
            # for each instance of MODELLER (i.e. each protein), and files are moved into
            # self.output_dir afterwards. If --black-no-sugar is provided, everything is moved.
            # Otherwise, only pertinent files are moved. See move_results_to_output_dir()
            self.modeller_dir = filesnpaths.get_temp_directory_path()

            run.warning("Working directory: {}".format(self.modeller.directory),
                         header='MODELLING STRUCTURE FOR GENE ID {}'.format(self.modeller.gene_id),
                         lc="green")

            self.self.gene_fasta_path = os.path.join(self.modeller_dir, "{}.fasta".format(gene))

            dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, self.gene_fasta_path, set([gene]))

            self.run_modeller()

            self.move_results_to_output_dir()

        """
        #vvvvvvvvvvvvvvvvvvvvvvv UGLY DONT LOOK, WILL BE REPLACED WITH ABOVE CODE vvvvvvvvvvvvvvvvvvvv
        #vvvvvvvvvvvvvvvvvvvvvvv UGLY DONT LOOK, WILL BE REPLACED WITH ABOVE CODE vvvvvvvvvvvvvvvvvvvv
        #vvvvvvvvvvvvvvvvvvvvvvv UGLY DONT LOOK, WILL BE REPLACED WITH ABOVE CODE vvvvvvvvvvvvvvvvvvvv
        all_genes_fasta_path = filesnpaths.get_temp_file_path()
        dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, all_genes_fasta_path)
        fasta = u.SequenceSource(all_genes_fasta_path)
        while next(fasta):
            if int(fasta.id) not in self.genes_of_interest:
                continue

            # MODELLER outputs a lot of stuff into its working directory. A temporary directory is
            # made for each instance of MODELLER (i.e. each protein), and files are moved into
            # self.output_dir afterwards. If --black-no-sugar is provided, everything is moved.
            # Otherwise, only pertinent files are moved. See move_results_to_output_dir()
            self.modeller_dir = filesnpaths.get_temp_directory_path()
            self.gene_fasta_path = filesnpaths.get_temp_file_path()
            gene_fasta = u.FastaOutput(self.gene_fasta_path)
            gene_fasta.write_id(fasta.id)
            gene_fasta.write_seq(fasta.seq, split = False)
            gene_fasta.close()

            self.run_modeller()
        #^^^^^^^^^^^^^^^^^^^^^^^^ UGLY DONT LOOK, WILL BE REPLACED WITH ABOVE CODE ^^^^^^^^^^^^^^^^^^^^
        #^^^^^^^^^^^^^^^^^^^^^^^^ UGLY DONT LOOK, WILL BE REPLACED WITH ABOVE CODE ^^^^^^^^^^^^^^^^^^^^
        #^^^^^^^^^^^^^^^^^^^^^^^^ UGLY DONT LOOK, WILL BE REPLACED WITH ABOVE CODE ^^^^^^^^^^^^^^^^^^^^


    def move_results_to_output_dir(self):
        """
        if --black-no-sugar, all files from MODELLERs directory are recursively moved into
        output_gene_dir.  Otherwise, the list of files we care about are defined in this function
        and moved into output_gene_dir.
        """
        output_gene_dir = os.path.join(self.output_dir, self.modeller.gene_id)
        filesnpaths.check_output_directory(output_gene_dir)

        if self.full_output:
            shutil.move(self.modeller.directory, output_gene_dir)

        else:
            filesnpaths.gen_output_directory(output_gene_dir)
            best_structure_filepath = self.modeller.pick_best_model(self.best)
            list_to_keep = [best_structure_filepath,
                            self.modeller.alignment_pap_path,
                            self.modeller.search_results_path,
                            self.modeller.model_info_path]
            for filepath in list_to_keep:
                shutil.move(filepath, output_gene_dir)

        run.warning("results folder: {}".format(output_gene_dir),
                     header='FINISHED STRUCTURE FOR GENE ID {}'.format(self.modeller.gene_id),
                     lc="green")


    def run_modeller(self):
            self.modeller = MODELLER.MODELLER(self.gene_fasta_path, directory = self.modeller_dir)

            run.warning("Working directory: {}".format(self.modeller.directory),
                         header='MODELLING STRUCTURE FOR GENE ID {}'.format(self.modeller.gene_id),
                         lc="blue")

            try:
                self.modeller.run_fasta_to_pir()

                self.modeller.check_database()

                self.modeller.run_search()

                self.parse_search_results()

                self.download_structures()

                self.modeller.run_align_to_templates(self.top_seq_matches)

                self.modeller.run_get_model(self.num_models, self.deviation, self.very_fast)

                self.modeller.tidyup()

                self.move_results_to_output_dir()

            except self.modeller.EndModeller as e:
                print(e)
                self.modeller.abort()



