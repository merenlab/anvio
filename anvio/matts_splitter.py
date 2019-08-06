def export_locus(self, gene_callers_id, output_path_prefix):
        """Takes a gene callers ID, and exports a contigs database.

           Output path prefix should be unique for every export locus call. If the prefix you provide
           looks like this:

                >>> output_path_prefix = '/path/to/dir/file_name_prefix'

           the output files will be stored as this:

                >>> '/path/to/dir/file_name_prefix.fa'
                >>> '/path/to/dir/file_name_prefix.db'

           """

        if os.path.isdir(output_path_prefix):
            raise ConfigError("Output path prefix can't be a directory name...")

        filesnpaths.is_output_file_writable(output_path_prefix + '.fa')

        if not self.contigs_db:
            self.contigs_db = dbops.ContigsSuperclass(self.args, r=self.run_object)
            self.contigs_db.init_functions()

        gene_call = self.contigs_db.genes_in_contigs_dict[gene_callers_id]
        contig_name = self.contigs_db.genes_in_contigs_dict[gene_callers_id]['contig']
        genes_in_contig_sorted = sorted(list(self.contigs_db.contig_name_to_genes[contig_name]))

        D = lambda: 1 if gene_call['direction'] == 'f' else -1
        premature = False

        self.run.info("Contig name", contig_name)
        self.run.info("Contig length", self.contigs_db.contigs_basic_info[contig_name]['length'])
        self.run.info("Num genes in contig", len(genes_in_contig_sorted))
        self.run.info("Target gene call", gene_callers_id)
        self.run.info("Target gene direction", "Forward" if D() == 1 else "Reverse", mc = 'green' if D() == 1 else 'red')

        gene_1 = gene_callers_id - self.num_genes_list[0] * D()
        gene_2 = gene_callers_id + self.num_genes_list[1] * D()
        first_gene_of_the_block = min(gene_1, gene_2)
        last_gene_of_the_block = max(gene_1, gene_2)

        self.run.info("First and last gene of the locus (raw)", "%d and %d" % (first_gene_of_the_block, last_gene_of_the_block))

        # getting the ids for the first and last genes in the contig
        last_gene_in_contig = genes_in_contig_sorted[-1][0]
        first_gene_in_contig = genes_in_contig_sorted[0][0]

        if last_gene_of_the_block > last_gene_in_contig:
            last_gene_of_the_block = last_gene_in_contig
            premature = True

        if first_gene_of_the_block < first_gene_in_contig:
            first_gene_of_the_block = first_gene_in_contig
            premature = True

        if premature and self.remove_partial_hits:
            self.run.info_single("A premature locus is found .. the current configuration says 'skip'. Skipping.", mc="red", nl_before=1)
            return
        elif premature and not self.remove_partial_hits:
            self.run.info_single("A premature locus is found .. the current configuration says 'whatevs'. Anvi'o will continue.", mc="yellow", nl_before=1, nl_after=1)

        self.run.info("First and last gene of the locus (final)", "%d and %d" % (first_gene_of_the_block, last_gene_of_the_block))

        locus_start = self.contigs_db.genes_in_contigs_dict[first_gene_of_the_block]['start']
        locus_stop = self.contigs_db.genes_in_contigs_dict[last_gene_of_the_block]['stop']

        # being a performance nerd here yes
        contig_sequence = db.DB(self.input_contigs_db_path, None, ignore_version=True) \
                            .get_some_rows_from_table(t.contig_sequences_table_name,
                                                      where_clause="contig='%s'" % contig_name)[0][1]
        locus_sequence = contig_sequence[locus_start:locus_stop]

        ################################
        # Flanking gene locus extraction
        ################################

        # Parse flanking genes argument

        flanking_genes_list = self.flanking_genes.split(',')

        flank_1 = flanking_genes_list[0]
        flank_2 = flanking_genes_list[1]

        # Get the contigs_db ready
        self.run.info('Mode', 'Function search')
        contigs_db = dbops.ContigsSuperclass(self.args, r=self.run_object) # This somehow gives us the original contigs_db

        # Use flanking gene search terms method to retrieve gene_caller_Id

        ## flank_1
        contigs_db.init_functions() # This somehow allows us to search the function in the contigs_db
        self.run.info('Search term', flank_1, mc='green')
        self.run.info('Function calls being used', ', '.join(contigs_db.gene_function_call_sources))

        foo, search_report_flank_1 = contigs_db.search_for_gene_functions([flank_1], requested_sources=self.annotation_sources, verbose=True)

        gene_caller_ids_of_interest_flank_1 = search_report_flank_1[0][0]

        ## flank_2
        #contigs_db.init_functions() # dont need to repeat
        self.run.info('Search term', flank_2, mc='green')
        self.run.info('Function calls being used', ', '.join(contigs_db.gene_function_call_sources))

        foo, search_report_flank_2 = contigs_db.search_for_gene_functions([flank_2], requested_sources=self.annotation_sources, verbose=True)

        gene_caller_ids_of_interest_flank_2 = search_report_flank_2[0][0]

        # Use flanking gene gene_caller_IDs to cut out locus

        # Get positions to cut at
        locus_start_new = self.contigs_db.genes_in_contigs_dict[gene_caller_ids_of_interest_flank_1]['start']
        locus_stop_new = self.contigs_db.genes_in_contigs_dict[gene_caller_ids_of_interest_flank_2]['stop']
       
       # Get original contig sequence (probably the entire genome)
        contig_sequence = db.DB(self.input_contigs_db_path, None, ignore_version=True) \
                            .get_some_rows_from_table(t.contig_sequences_table_name,
                                                      where_clause="contig='%s'" % contig_name)[0][1]
        #... and cut!

        # Lets make a conditional so that if there are locus coordinates from either
        # flanking genes or search term the locus will be cut... sometime in the future... its working done touch!

        locus_sequence_new = contig_sequence[locus_start_new:locus_stop_new]


        # Copy of code from above but with the .....
        # here we will create a gene calls dict for genes that are specific to our locus. since we trimmed
        # the contig sequence to the locus of interest, we will have to adjust start and stop positions of
        # genes in teh gene calls dict.
        
        if first_gene_of_the_block and last_gene_of_the_block: 

            # Original
            locus_gene_calls_dict = {}
            for g in range(first_gene_of_the_block, last_gene_of_the_block + 1):
                locus_gene_calls_dict[g] = copy.deepcopy(self.contigs_db.genes_in_contigs_dict[g])
                excess = self.contigs_db.genes_in_contigs_dict[first_gene_of_the_block]['start']
                locus_gene_calls_dict[g]['start'] -= excess
                locus_gene_calls_dict[g]['stop'] -= excess

        else:
            # New with flanking
            locus_gene_calls_dict = {}
            for g in range(gene_caller_ids_of_interest_flank_1, gene_caller_ids_of_interest_flank_2 + 1):
                locus_gene_calls_dict[g] = copy.deepcopy(self.contigs_db.genes_in_contigs_dict[g])
                excess = self.contigs_db.genes_in_contigs_dict[gene_caller_ids_of_interest_flank_1]['start']
                locus_gene_calls_dict[g]['start'] -= excess
                locus_gene_calls_dict[g]['stop'] -= excess



        self.run.info("Locus gene call start/stops excess (nts)", excess)

        if D() != 1 and self.reverse_complement_if_necessary:
            reverse_complement = True
        else:
            reverse_complement = False

        self.run.info('Reverse complementing everything', reverse_complement, mc='green')

        # report a stupid FASTA file.
        if self.include_fasta_output:
            fasta_file_path = output_path_prefix + ".fa"

            self.run.info("Output FASTA file", fasta_file_path)
            with open(fasta_file_path, 'w') as f:
                locus_header = contig_name + ' ' + \
                               '|'.join(['target:%s' % ','.join(self.targets),
                                         'sources:%s' % ','.join(self.sources),
                                         'query:%s' % self.search_term or 'None',
                                         'hit_contig:%s' % contig_name,
                                         'hit_gene_callers_id:%s' % str(gene_callers_id),
                                         'project_name:%s' % self.contigs_db.a_meta['project_name'].replace(' ', '_').replace("'", '_').replace('"', '_'),
                                         'locus:%s,%s' % (str(first_gene_of_the_block), str(last_gene_of_the_block)),
                                         'nt_positions_in_contig:%s:%s' % (str(locus_start), str(locus_stop)),
                                         'premature:%s' % str(premature),
                                         'reverse_complemented:%s' % str(reverse_complement)])

                f.write('>%s\n' % locus_header)
                f.write('%s\n' % utils.rev_comp(locus_sequence) if reverse_complement else locus_sequence)

        # report a fancy anvi'o contigs database
        self.store_locus_as_contigs_db(contig_name,
                                       locus_sequence,
                                       locus_gene_calls_dict,
                                       output_path_prefix,
                                       reverse_complement)
