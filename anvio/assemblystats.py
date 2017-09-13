# -*- coding: utf-8

import json

import anvio
import anvio.dbops as dbops

class AssemblyStats():
    def __init__(self, args):
        self.contigs_db = args.contigs_db

    def get_summary_dict(self):
        summary = {}

        contig_lengths = sorted(dbops.ContigsDatabase(self.contigs_db).db.get_single_column_from_table('contigs_basic_info', 'length'), reverse=True)
        total_length = sum(contig_lengths)
        size = len(contig_lengths)

        summary['total_length'] = total_length
        summary['size'] = size
        summary['n_values'] = []

        temp_length = 0
        n_index = 1;
        for i in range(size):
            temp_length += contig_lengths[i]
            print(temp_length, total_length, n_index)
            if (temp_length > ((total_length / 100) * n_index)):
                summary['n_values'].append({
                    'num_contigs': i + 1,
                    'length': contig_lengths[i]
                    })

                n_index += 1
                print(n_index)

        return summary


class AssemblyInteractive():
    def __init__(self, args):
        self.mode = 'assembly'
        self.assembly_stats = AssemblyStats(args)

    def get_assembly_stats(self):
        return json.dumps(self.assembly_stats.get_summary_dict())
