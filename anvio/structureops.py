# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

class Structure:

    def __init__(self, args):
        self.args = args

        # initialize self.arg parameters
        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.genes_of_interest_path = A('genes_of_interest', null)
        self.splits_of_interest_path = A('splits_of_interest', null)
        self.bin_id = A('bin_id', null)
        self.collection_name = A('collection_name', null)
        self.gene_caller_id = A('gene_caller_id', null)
        self.output_dir = A('output_dir', null)
        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None

        # identify which genes user wants to model structures of
        self.get_genes_of_interest()

        #



    def model_structure(gene_id):
        ######### CREATE TEMP PATH FOR MODELLER ##################

        tmp_dir = filesnpaths.get_temp_directory_path()

        ############# WRITE FASTA OF ALL GENES ###########
        """
        The ideal workflow will be to rewrite export_aa_sequences_from_contigs_db to accept
        a list of genes of interest. I will then loop through the genes of interest and 
        write single-gene FASTA files. Currently, I parse through a FASTA of ALL genes and 
        ignore those that are not in the genes of interest.
        """
        # create fasta file of all genes in database
        progress.new('GETTING GENE SEQUENCES')
        progress.update("Exporting AA gene sequences from contigs database...")
        all_genes_fasta_path = os.path.join(tmp_dir, "all_genes.fasta")
        dbops.export_aa_sequences_from_contigs_db(args.contigs_db, all_genes_fasta_path)
        progress.end()

        ######## LOOP THROUGH EACH GENE IN FASTA ########

        fasta = u.SequenceSource(all_genes_fasta_path)

        while next(fasta):

            # do nothing if not gene of interest
            if fasta.id not in genes:
                continue

            # define temp file with format <gene_callers_id>.fasta
            single_gene_fasta_path = os.path.join(tmp_dir, "{}.fasta".format(fasta.id))

            # write single gene fasta
            single_gene_fasta = u.FastaOutput(single_gene_fasta_path)
            single_gene_fasta.write_id(fasta.id)
            single_gene_fasta.write_seq(fasta.seq, split = False)

            # close fasta
            single_gene_fasta.close()

            ############### INITIALIZE MODELLER DRIVER ############

            mod = MODELLER.MODELLER(single_gene_fasta_path, directory=tmp_dir)
            mod.process()


    def get_genes_of_interest(self):
        """
        nabs the genes of interest based on user arguments (self.args)
        """
        # identify the gene caller ids of all genes available
        self.genes_in_database = set(dbops.ContigsSuperclass(args).genes_in_splits.keys())
        if not self.genes_in_database:
            raise ConfigError("This contigs database does not contain any identified genes...")

        # settling genes of interest
        if self.genes_of_interest_path and self.gene_caller_id:
            raise ConfigError("You can't provide a gene caller id from the command line, and a list of gene caller ids\
                               as a file at the same time, obviously.")

        if self.gene_caller_id:
            try:
                self.gene_caller_id = int(self.gene_caller_id)
            except:
                raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % str(self.gene_caller_id))

            self.genes_of_interest = set([self.gene_caller_id])
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
        bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in self.gene_callers_id_to_split_name_dict]
        if bad_gene_caller_ids:
            raise ConfigError("The gene caller id you provided is not known to this contigs database. You have only 2 lives\
                               left. 2 more mistakes, and anvi'o will automatically uninstall itself. Yes, seriously :(")

        # Finally, raise warning if number of genes is greater than 20
        if len(self.genes_of_interest) > 20:
            run.warning("Modelling protein structures is no joke. The number of genes you are interested in modelling \
                         protein structures for is {}, which is a lot. If even half of them turn out to be modellable, \
                         this will take some time, so you shouldn't do it (but it will be so cool, so you should totally \
                         do it).".format(len(self.genes_of_interest)))


