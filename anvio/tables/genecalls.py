
import anvio.tables as t
import anvio.dbops as dbops
import anvio.terminal as terminal

run = terminal.Run()
progress = terminal.Progress()


class TablesForGeneCalls(object):
    def __init__(self, db, run=run, progress=progress, debug=False):
        self.run = run
        self.progress = progress
        self.db = db
        self.debug = debug


    def store(self, gene_calls_dict, protein_sequences, append=False):
        if not append:
            self.db.clear_table(t.genes_in_contigs_table_name)
            self.db.clear_table(t.gene_protein_sequences_table_name)
        else:
            # so we are in the append mode. We must remove all the previous entries from genes in contigs
            # that matches to the incoming sources. otherwhise we may end up with many duplicates in the db.
            sources = set([v['source'] for v in gene_calls_dict.values()])
            for source in sources:
                contigs_db._exec('''DELETE FROM %s WHERE source = "%s"''' % (t.genes_in_contigs_table_name, source))

        self.progress.new('Processing')
        self.progress.update('Entering %d gene calls into the db ...' % (len(gene_calls_dict)))

        db_entries = [tuple([entry_id] + [gene_calls_dict[entry_id][h] for h in t.genes_in_contigs_table_structure[1:]]) for entry_id in gene_calls_dict]
        self.db.insert_many(t.genes_in_contigs_table_name, db_entries)

        db_entries = [tuple([entry_id] + [protein_sequences[entry_id]]) for entry_id in gene_calls_dict]
        self.db.insert_many(t.gene_protein_sequences_table_name, db_entries)

        self.progress.end()



    def populate_genes_in_splits_tables(self):
        #Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)
        gene_calls_dict = self.db.get_table_as_dict(t.genes_in_contigs_table_name)

        genes_in_splits = dbops.GenesInSplits()
        # build a dictionary for fast access to all proteins identified within a contig
        gene_calls_in_contigs_dict = {}
        for gene_callers_id in self.gene_calls_dict:
            contig = self.gene_calls_dict[gene_callers_id]['contig']
            if contig in gene_calls_in_contigs_dict:
                gene_calls_in_contigs_dict[contig].add(gene_callers_id)
            else:
                gene_calls_in_contigs_dict[contig] = set([gene_callers_id])

        contigs_without_any_gene_calls = list(set(self.contigs_info.keys()) - set(gene_calls_in_contigs_dict.keys()))
        run.info('Contigs with at least one gene call', '%d of %d (%.1f%%)' % (len(gene_calls_in_contigs_dict),
                                                                               len(self.contigs_info),
                                                                               len(gene_calls_in_contigs_dict) * 100.0 / len(self.contigs_info)))

        for contig in contigs_without_any_gene_calls:
            gene_calls_in_contigs_dict[contig] = set([])

        splits_dict = {}
        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                gene_start_stops = []
                # here we go through all genes in the contig and identify the all the ones that happen to be in
                # this particular split to generate summarized info for each split. BUT one important that is done
                # in the following loop is genes_in_splits.add call, which populates GenesInSplits class.
                for gene_callers_id in gene_calls_in_contigs_dict[contig]:
                    if self.gene_calls_dict[gene_callers_id]['stop'] > start and self.gene_calls_dict[gene_callers_id]['start'] < stop:
                        gene_start_stops.append((self.gene_calls_dict[gene_callers_id]['start'], self.gene_calls_dict[gene_callers_id]['stop']), )
                        genes_in_splits.add(split_name, start, stop, gene_callers_id, self.gene_calls_dict[gene_callers_id]['start'], self.gene_calls_dict[gene_callers_id]['stop'])

                # here we identify genes that are associated with a split even if one base of the gene spills into
                # the defined start or stop of a split, which means, split N, will include genes A, B and C in this
                # scenario:
                #
                # contig: (...)------[ gene A ]--------[     gene B    ]----[gene C]---------[    gene D    ]-----(...)
                #         (...)----------x---------------------------------------x--------------------------------(...)
                #                        ^ (split N start)                       ^ (split N stop)
                #                        |                                       |
                #                        |<-              split N              ->|
                #
                # however, when looking at the coding versus non-coding nucleotide ratios in a split, we have to make
                # sure that only the relevant portion of gene A and gene C is counted:
                total_coding_nts = 0
                for gene_start, gene_stop in gene_start_stops:
                    total_coding_nts += (gene_stop if gene_stop < stop else stop) - (gene_start if gene_start > start else start)

                splits_dict[split_name] = {'num_genes': len(gene_start_stops),
                                           'avg_gene_length': numpy.mean([(l[1] - l[0]) for l in gene_start_stops]) if len(gene_start_stops) else 0.0,
                                           'ratio_coding': total_coding_nts * 1.0 / (stop - start),
                                           }

        # open connection
        contigs_db = ContigsDatabase(self.db_path)
        # push raw entries for splits table
        db_entries = [tuple([split] + [splits_dict[split][h] for h in t.genes_in_splits_summary_table_structure[1:]]) for split in splits_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.genes_in_splits_summary_table_name, db_entries)

        # push entries for genes in splits table
        db_entries = [tuple([entry_id] + [genes_in_splits.splits_to_prots[entry_id][h] for h in t.genes_in_splits_table_structure[1:]]) for entry_id in genes_in_splits.splits_to_prots]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.genes_in_splits_table_name, db_entries)
        # disconnect
        contigs_db.disconnect()

