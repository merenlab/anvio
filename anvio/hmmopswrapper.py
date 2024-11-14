# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    HMM related operations.
"""

import argparse
import hashlib

import anvio.terminal as terminal
import anvio.ccollections as ccollections

from anvio.genomedescriptions import GenomeDescriptions
from anvio.hmmops import SequencesForHMMHits


from anvio.errors import ConfigError


class SequencesForHMMHitsWrapperForMultipleContigs(
    SequencesForHMMHits, GenomeDescriptions
):
    """A class that generates an instance of SequencesForHMMHits with multiple contigs databases.

    An instance from this class will be a fully operational SequencesForHMMHits instance.
    """

    def __init__(
        self, args, hmm_sources, run=terminal.Run(), progress=terminal.Progress()
    ):
        self.args = args
        self.run = run
        self.progress = progress
        self.hmm_sources = hmm_sources

        self.splits_dict = {}

        # process genome descriptions
        GenomeDescriptions.__init__(self, args, run=self.run, progress=self.progress)
        self.load_genomes_descriptions(skip_functions=True, init=False)

        hmm_sources_in_all_genomes = self.get_HMM_sources_common_to_all_genomes(
            sources_that_must_be_common=hmm_sources
        )

        if not len(hmm_sources_in_all_genomes):
            raise ConfigError(
                "There are no HMM sources among your external genomes that occur in every genome :/"
            )

        # initialize the super
        SequencesForHMMHits.__init__(
            self, None, sources=hmm_sources, run=self.run, progress=self.progress
        )

        num_internal_genomes = len(
            set(
                [
                    g["genome_hash"]
                    for g in self.genomes.values()
                    if "profile_db_path" in g
                ]
            )
        )
        collection_names = set(
            [g["collection_id"] for g in self.genomes.values() if "collection_id" in g]
        )

        if num_internal_genomes:
            self.run.warning(
                "SequencesForHMMHitsWrapperForMultipleContigs class is speaking (yes, the class is "
                "quite aware of its very long name thankyouverymuch). Of the total %d genome descriptions "
                "it was given, %d seem to represent internal genomes with bins in collection(s) '%s'. Anvi'o "
                "will make sure HMM hits to be used for downstream analyses are only those that match to contigs "
                "that were included in those selections."
                % (
                    len(self.genomes),
                    num_internal_genomes,
                    ", ".join(collection_names),
                ),
                lc="green",
            )

        # very hacky code follows. here we generate a self SequencesForHMMHits object,
        # and we will fill everything in it with slightly modified information so multiple
        # contigs databases could be processed by this talented class seamlessly.
        hmm_hits_splits_counter = 0
        for genome_name in self.genomes:
            g = self.genomes[genome_name]
            contigs_db_path = g["contigs_db_path"]
            contigs_db_hash = g["contigs_db_hash"]

            # this is an important variable and allows us to track origins of HMM hits for bins
            # and individual contigs databases seamlessly. if you want to understand truly what
            # the hell does this mean, look at `get_genome_hash_for_external_genome` and
            # `get_genome_hash_for_internal_genome` functions in `genomedescriptions.py`.
            genome_hash = None

            # here we check if the genome descriptions contain reference to a collection name,
            # because if it is the case, we need to focus only on hmm hits that are relevant
            # to splits in this collection:
            if "collection_id" in g:
                if ("bin_id" not in g) or ("profile_db_path" not in g):
                    raise ConfigError(
                        "There is something VERY weird going on. Your genome descriptions object contains "
                        "a collection name, yet it doesn't know anything about a bin name or profile database "
                        "path. While this is very interesting because it should never happen, anvi'o will say "
                        "goodbye and abruptly quit in confusion :("
                    )

                # setup an args object, and recover the split names of interest
                args = argparse.Namespace(
                    profile_db=g["profile_db_path"],
                    contigs_db=g["contigs_db_path"],
                    bin_id=g["bin_id"],
                    collection_name=g["collection_id"],
                )
                split_names_of_interest = ccollections.GetSplitNamesInBins(
                    args
                ).get_split_names_only()
                genome_hash = hashlib.sha224(
                    "_".join(
                        ["".join(split_names_of_interest), contigs_db_hash]
                    ).encode("utf-8")
                ).hexdigest()[0:12]

                # current hmm hits now will match to the collection
                current = SequencesForHMMHits(
                    contigs_db_path,
                    sources=hmm_sources,
                    split_names_of_interest=split_names_of_interest,
                )
            else:
                current = SequencesForHMMHits(contigs_db_path, sources=hmm_sources)
                genome_hash = contigs_db_hash

            for hmm_hit_id in current.hmm_hits:
                hit = current.hmm_hits[hmm_hit_id]
                hit["gene_callers_id"] = "%s_%d" % (
                    contigs_db_hash,
                    hit["gene_callers_id"],
                )
                hit["genome_hash"] = genome_hash
                self.hmm_hits["%s_%d" % (contigs_db_hash, hmm_hit_id)] = hit

            if not self.hmm_hits_info:
                for hmm_source in hmm_sources_in_all_genomes:
                    self.hmm_hits_info[hmm_source] = current.hmm_hits_info[hmm_source]

            for hit in current.hmm_hits_splits.values():
                hit["split"] = "%s_%s" % (contigs_db_hash, hit["split"])
                hit["hmm_hit_entry_id"] = "%s_%d" % (
                    contigs_db_hash,
                    hit["hmm_hit_entry_id"],
                )
                self.hmm_hits_splits[hmm_hits_splits_counter] = hit
                hmm_hits_splits_counter += 1

            for seq in current.contig_sequences:
                self.contig_sequences["%s_%s" % (contigs_db_hash, seq)] = (
                    current.contig_sequences[seq]
                )

            for seq in current.aa_sequences:
                self.aa_sequences["%s_%s" % (contigs_db_hash, seq)] = (
                    current.aa_sequences[seq]
                )

            for gene_callers_id in current.genes_in_contigs:
                entry = current.genes_in_contigs[gene_callers_id]
                entry["contig"] = "%s_%s" % (contigs_db_hash, entry["contig"])
                self.genes_in_contigs["%s_%d" % (contigs_db_hash, gene_callers_id)] = (
                    entry
                )

            self.splits_dict[genome_name] = [
                "%s_%s" % (contigs_db_hash, s) for s in current.splits_in_contigs
            ]
