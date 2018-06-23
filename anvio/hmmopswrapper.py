# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    HMM related operations.
"""

import anvio.terminal as terminal

from anvio.genomedescriptions import GenomeDescriptions
from anvio.hmmops import SequencesForHMMHits


run = terminal.Run()
progress = terminal.Progress()


class SequencesForHMMHitsWrapperForMultipleContigs(SequencesForHMMHits, GenomeDescriptions):
    """A class that generates an instance of SequencesForHMMHits with multiple contigs databases.

       An instance from this class will be a fully operational SequencesForHMMHits instance."""

    def __init__(self, args, hmm_sources):
        self.args = args
        self.hmm_sources = hmm_sources

        self.splits_dict = {}

        # initialize the super
        SequencesForHMMHits.__init__(self, None, sources=hmm_sources)

        # process genome descriptions
        GenomeDescriptions.__init__(self, args)
        self.load_genomes_descriptions(skip_functions=True, init=False)
        hmm_sources_in_all_genomes = self.get_HMM_sources_common_to_all_genomes()

        # very hacky code follows. here we generate a self SequencesForHMMHits object,
        # and we will fill everything in it with slightly modified information so multiple
        # contigs databases could be processed by this talented class seamlessly.
        hmm_hits_splits_counter = 0
        for genome_name in self.genomes:
            contigs_db_path = self.genomes[genome_name]['contigs_db_path']
            contigs_db_hash = self.genomes[genome_name]['contigs_db_hash']

            current = SequencesForHMMHits(contigs_db_path, sources = hmm_sources)

            for hmm_hit_id in current.hmm_hits:
                hit = current.hmm_hits[hmm_hit_id]
                hit['gene_callers_id'] = '%s_%d' % (contigs_db_hash, hit['gene_callers_id'])
                self.hmm_hits['%s_%d' % (contigs_db_hash, hmm_hit_id)] = hit

            if not self.hmm_hits_info:
                for hmm_source in hmm_sources_in_all_genomes:
                    self.hmm_hits_info[hmm_source] = current.hmm_hits_info[hmm_source]

            for hit in current.hmm_hits_splits.values():
                hit['split'] = '%s_%s' % (contigs_db_hash, hit['split'])
                hit['hmm_hit_entry_id'] = '%s_%d'% (contigs_db_hash, hit['hmm_hit_entry_id'])
                self.hmm_hits_splits[hmm_hits_splits_counter] = hit
                hmm_hits_splits_counter += 1

            for seq in current.contig_sequences:
                self.contig_sequences['%s_%s' % (contigs_db_hash, seq)] = current.contig_sequences[seq]

            for seq in current.aa_sequences:
                self.aa_sequences['%s_%s' % (contigs_db_hash, seq)] = current.aa_sequences[seq]

            for gene_callers_id in current.genes_in_contigs:
                entry = current.genes_in_contigs[gene_callers_id]
                entry['contig'] = '%s_%s' % (contigs_db_hash, entry['contig'])
                self.genes_in_contigs['%s_%d' % (contigs_db_hash, gene_callers_id)] = entry

            self.splits_dict[genome_name] = ['%s_%s' % (contigs_db_hash, s) for s in current.splits_in_contigs]
