# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

import os
import shutil

import anvio.dbops as dbops
import anvio.fastalib as u
import anvio.filesnpaths as filesnpaths
import anvio.drivers.MODELLER as MODELLER

from anvio.errors import ConfigError, FilesNPathsError

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

        # check output and define absolute path
        self.output_dir = filesnpaths.check_output_directory(self.output_dir, ok_if_exists=False)

        # identify which genes user wants to model structures for
        self.get_genes_of_interest()


    def model_structures(self):

        # MODELLER outputs a lot of stuff into its working directory. By default, a temporary
        # directory is made for MODELLER, and pertinent files are moved into the self.output_dir
        # afterwards. If --black-no-sugar is provided, MODELLER's directory is self.output_dir
        if self.full_output:
            MODELLER_dir = self.output_dir
            os.mkdir(MODELLER_dir)
        else:
            MODELLER_dir = filesnpaths.get_temp_directory_path()

        """
        sqlite-migration branch has a parameter passed to dbops.export_aa_sequences_from_contigs_db
        that lets you pass genes of interest. When these branches are merged, the code will look
        like this:

        for gene in self.genes_of_interest:
            gene_fasta_path = os.path.join(MODELLER_DIR, "{}.fasta".format(gene))
            dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, gene_fasta_path, set([gene]))
            modeller = MODELLER.MODELLER(gene_fasta_path, directory = MODELLER_dir)
            modeller.process()
        """
        # create fasta file of all genes in database
        all_genes_fasta_path = os.path.join(MODELLER_dir, "all_genes.fasta")
        dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, all_genes_fasta_path)
        fasta = u.SequenceSource(all_genes_fasta_path)
        while next(fasta):
            if int(fasta.id) not in self.genes_of_interest:
                continue
            single_gene_fasta_path = os.path.join(MODELLER_dir, "{}.fasta".format(fasta.id))
            single_gene_fasta = u.FastaOutput(single_gene_fasta_path)
            single_gene_fasta.write_id(fasta.id)
            single_gene_fasta.write_seq(fasta.seq, split = False)
            single_gene_fasta.close()
            modeller = MODELLER.MODELLER(single_gene_fasta_path, directory = MODELLER_dir)
            modeller.process()

        # finished modelling, move pertinent contents FIXME I'm just moving everything
        if not self.full_output:
            shutil.move(MODELLER_dir, self.output_dir)


    def get_genes_of_interest(self):
        """
        nabs the genes of interest based on user arguments (self.args)
        """
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
            time.sleep(10)
            run.warning("Modelling protein structures is no joke. The number of genes you want protein structures for is \
                         {}, which is a lot (of time!). I'm putting you in timeout for 15 seconds, then I'm going to do \
                         what you said to do. CTRL + C to cancel.".format(len(self.genes_of_interest)))
            time.sleep(5)
            run.warning("YOU'RE CRAZY!!! (5 seconds left)")



