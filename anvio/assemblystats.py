# -*- coding: utf-8

import json

import anvio
import anvio.dbops as dbops
import anvio.hmmops as hmmops

class AssemblyStats():
    def __init__(self, args):
        self.contigs_db_path = args.contigs_db

    def get_summary_dict(self):
        summary = {}
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        hmm = hmmops.SequencesForHMMHits(self.contigs_db_path)

        contig_lengths = sorted(contigs_db.db.get_single_column_from_table('contigs_basic_info', 'length'), reverse=True)
        total_length = sum(contig_lengths)
        size = len(contig_lengths)

        summary['total_length'] = total_length
        summary['size'] = size
        summary['n_values'] = self.calculate_N_values(contig_lengths, total_length, N=100)
        summary['single_copy_gene_counts'] = hmm.get_single_copy_gene_counts()

        return summary


    def calculate_N_values(self, contig_lengths, total_length, N=100):
        results = []

        temp_length = 0
        contigs_index = 0
        n_index = 1

        while n_index <= N:
            if (temp_length >= ((total_length / N) * n_index)):
                results.append({
                        'num_contigs': contigs_index,
                        'length':      contig_lengths[contigs_index - 1]
                    })
                n_index += 1
            else:
                temp_length += contig_lengths[contigs_index]
                contigs_index += 1

        return results


class AssemblyInteractive():
    def __init__(self, args):
        self.mode = 'assembly'
        self.assembly_stats = AssemblyStats(args)

    def get_assembly_stats(self):
        return json.dumps(self.assembly_stats.get_summary_dict())
