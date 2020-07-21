# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes for basic sequence properties and manipulations."""

from hashlib import sha224

from itertools import permutations

import anvio
import anvio.constants as constants

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class Codon:
    def __init__(self):
        pass

    def get_codon_to_codon_sequence_trajectory(self, start_codon, end_codon, as_amino_acids=False):
        '''
        this returns a list of all possible sequence trajectories to get from one codon to another
        assuming the least amount of mutations necessary. if as_amino_acids, the trajectories will
        be converted into amino acid space
        '''
        indices_of_variation = []
        for i in range(3):
            if start_codon[i] != end_codon[i]:
                indices_of_variation.append(i)
        index_trajectories = list(permutations(indices_of_variation))

        all_trajectories = []
        for index_trajectory in index_trajectories:
            sequence_trajectory = [start_codon]
            mutate = list(start_codon)
            for index in index_trajectory:
                mutate[index] = end_codon[index]
                sequence_trajectory.append(''.join(mutate))
            all_trajectories.append(sequence_trajectory)

        if as_amino_acids:
            # each codon is converted to an amino acid. if two adjacent codons are the same amino
            # acid then only one is kept
            for i, trajectory in enumerate(all_trajectories):
                for j, node in enumerate(trajectory):
                    trajectory[j] = constants.codon_to_AA[node]

                # gets rid of duplicates
                all_trajectories[i] = list(dict.fromkeys(trajectory))

        return all_trajectories


    def get_codon_to_codon_dist_dictionary(self):
        """
        Returns a dictionary containing the number the nucleotide difference
        between two codons, and the number of transitions & transversions required to
        mutate from one codon to another.

        USAGE:
        =====
        >>> dist['TTC']['AAA']
        (3, 2, 1)

        Here, 3 is the total number of nucleotide differences, 2 is the number
        of transitions, and 1 is the number of transversions. Of course, the number of
        transitions added to the number of transversions is equal to the number of
        nucleotide differences.
        """

        codons = list(constants.codon_to_AA.keys())
        dist = {}

        mutation_type = {"AT": "transition",
                         "CG": "transition",
                         "AG": "transversion",
                         "AC": "transversion",
                         "CT": "transversion",
                         "GT": "transversion"}

        for start_codon in codons:
            dist[start_codon] = {}

            for end_codon in codons:

                # s = number transitions
                # v = number transversions
                s = 0; v = 0
                for nt_pos in range(3):

                    pair = ''.join(sorted(start_codon[nt_pos] + end_codon[nt_pos]))

                    if pair[0] == pair[1]:
                        continue
                    if mutation_type[pair] == "transition":
                        s += 1
                    if mutation_type[pair] == "transversion":
                        v += 1

                dist[start_codon][end_codon] = (s+v, s, v)

        return dist


class Composition:
    def __init__(self, sequence):
        self.sequence = sequence
        self.GC_content = 0.0

        self.report()


    def report(self):
        s = self.sequence
        raw_length = len(s)

        self.A = s.count('A')
        self.T = s.count('T')
        self.C = s.count('C')
        self.G = s.count('G')
        self.N = raw_length - (self.A + self.T + self.C + self.G)
        length = raw_length - self.N

        if not length:
            # sequence is composed of only N's
            self.GC_content = 0.0
        else:
            self.GC_content = (self.G + self.C) * 1.0 / length


class SequenceDereplicator:
    """
    This class takes in a list of sequences and performs full-length or prefix dereplication.
    A list of sequences and corresponding list of sequence IDs are required.
    A list of extra information for each sequence is optional.
    The extra list can take any form -- a list of tuples, strings, mixed data types, etc.
    """

    def __init__(self, ids, seqs, extras=None):
        if len(ids) != len(seqs):
            raise ConfigError(
                "SequenceDereplicator takes a list of sequence IDs and a list of sequence strings. "
                "Your input lists were not the same length. "
                "The ID list had a length of %d, and the sequence list had a length of %d."
                % (len(ids), len(seqs)))
        if extras:
            if len(extras) != len(seqs):
                raise ConfigError(
                    "SequenceDereplicator takes an optional extra list of sequence information. "
                    "Your extra list was not the same length as the ID and sequence lists. "
                    "The extra list had a length of %d, while the others had a length of %d."
                    % (len(extras), len(seqs)))

        self.ids = ids
        self.seqs = seqs
        self.hashed_seqs = [sha224(seq.encode('utf-8')).hexdigest() for seq in self.seqs]
        if extras:
            self.extras = extras
        else:
            self.extras = []


    def full_length_dereplicate(self):
        """
        Identical sequences are clustered, returning a list of tuples.

        If no extra sequence information (self.extras) is provided, the list format is:
        clusters = [(sequence A, [member sequence X ID, member sequence Y ID, ...]), (sequence B, [member IDs]), ...]
        If extra information is provided:
        clusters = [(sequence A, [member sequence X ID, ...], [member sequence X extra info, ...]), ...]
        """

        # Search for identical sequences of the same length.
        seq_lengths = [len(seq) for seq in self.seqs]
        full_dict = {seq_length: {} for seq_length in sorted(set(seq_lengths))} # store clusters here
        extra_iter = iter(self.extras)
        for seq_length, hashed_seq, seq_id, seq in zip(seq_lengths, self.hashed_seqs, self.ids, self.seqs):
            inner_dict = full_dict[seq_length]
            try:
                if self.extras:
                    value = inner_dict[hashed_seq]
                    value[1].append(seq_id)
                    value[2].append(next(extra_iter))
                else:
                    inner_dict[hashed_seq][1].append(seq_id)
            except KeyError:
                if self.extras:
                    inner_dict[hashed_seq] = (seq, [seq_id], [next(extra_iter)])
                else:
                    inner_dict[hashed_seq] = (seq, [seq_id])

        clusters = []
        for inner_dict in full_dict.values():
            clusters.extend(inner_dict.values())
        clusters.sort(key=lambda value: -len(value[1]))

        return clusters


    def get_prefix_dict(self, ids, seqs, min_seq_length=None):
        """
        The returned prefix_dict contains hashes of prefix subsequences from all of the target sequences.
        prefix_dict is keyed by sequence length, with lengths spanning from the min to max sequence length;
        min seq length can be adjusted downward using the parameter.
        prefix_dict values are dictionaries themselves.
        Inner dict keys are prefix sequence hashes,
        and values are lists of target sequence IDs containing those hashes.
        """

        if min_seq_length:
            if min_seq_length < 1:
                raise ConfigError("The `min_seq_length` argument of `get_prefix_dict` "
                                  "must be at least 1, not %d." % min_seq_length)

        seq_lengths = [len(seq) for seq in seqs]
        min_seq_length = min_seq_length if min_seq_length else min(seq_lengths)
        max_seq_length = max(seq_lengths)

        prefix_dict = {prefix_length: {} for prefix_length in range(min_seq_length, max_seq_length + 1)}
        for prefix_length, inner_dict in prefix_dict.items():
            for seq_id, seq in zip(ids, seqs):
                if prefix_length > len(seq):
                    continue
                hashed_prefix = sha224(seq[: prefix_length].encode('utf-8')).hexdigest()
                if hashed_prefix in inner_dict:
                    inner_dict[hashed_prefix].append((seq_id, len(seq)))
                else:
                    inner_dict[hashed_prefix] = [(seq_id, len(seq))]
        for inner_dict in prefix_dict.values():
            for hashed_seq, parent_seqs in inner_dict.items():
                # Sort the list of IDs containing the prefix subsequence in descending order of sequence length.
                inner_dict[hashed_seq] = [seq_id for seq_id, seq_length in sorted(parent_seqs, key=lambda t: -t[1])]

        return prefix_dict


    def get_memberships(self, seed_ids, seed_seqs):
        """
        This function is called after clustering to find which sequences are in which clusters.
        Prefix dereplication can result in the same sequence occurring in multiple clusters.
        """

        # Find the prefix hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed prefix subsequences,
        # which only requires searching one list of hashes representing prefix sequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_prefix_dict = self.get_prefix_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in self.seqs]))

        missed_seed_ids = []
        memberships = []
        extra_iter = iter(self.extras)
        for query_id, query_seq, query_hash in zip(self.ids, self.seqs, self.hashed_seqs):
            inner_dict = seed_prefix_dict[len(query_seq)]
            if query_hash in inner_dict:
                seed_ids = inner_dict[query_hash]
            else:
                # The absense of the query sequence in seed_prefix_dict
                # is explained in detail in the prefix_dereplicate method.
                inner_dict[query_hash] = [query_id]
                seed_ids = [query_id]
                missed_seed_ids.append(query_id)
            # # REMOVE try block
            # try:
            #     hit_ids = hit_prefix_dict[len(query_seq)][query_hash]
            # except KeyError:
            #     print(query_hash)
            #     print(query_id)
            #     print(query_seq)
            #     raise Exception
            if self.extras:
                memberships.append((query_id, seed_ids, next(extra_iter)))
            else:
                memberships.append((query_id, seed_ids))

        return memberships, missed_seed_ids


    def prefix_dereplicate(self):
        """
        Prefix dereplication clusters sequences that match the beginning of a longer (seed) sequence.
        Sequences can occur in multiple clusters, as they may be subsequences of distinct seed sequences:
        Cluster 1:
        ACGTACGTACGTACGT (seed, seq A)
        ACGTACGTACGT (seq X)
        ACGTACGT (seq Y)
        Cluster 2:
        ACGTACGTACGTACGG (seed, seq B)
        ACGTACGTACGT (seq X)
        ACGTACGT (seq Y)

        Two lists are returned, a `clusters` list and a `memberships` list.

        Here is the format of the `clusters` list when extra sequence information is not provided:
        clusters = [(seed seq A ID, seed sequence A, [(member seq A ID, member seq A length), (member seq X ID, member seq X length), ...]), ...]
        Here is the format when it is provided:
        clusters = [(seed seq A ID, seed sequence A, [(member seq A ID, member seq A length, member seq A extra info), ...]), ...]

        The membership list has an entry for each input sequence.
        Here is the format when extra sequence information is not provided:
        memberships = [(member sequence X, [seed seq A ID, seed seq B ID, ...]), ...]
        Here is the format when it is provided:
        memberships = [(member sequence X, [seed seq A ID, seed seq B ID, ...], member seq X extra info), ...]
        """

        prefix_dict = self.get_prefix_dict(self.ids, self.seqs)

        hit_dict = {}
        cluster_dict = {}
        extra_iter = iter(self.extras)
        for query_id, query_seq, query_hash in zip(self.ids, self.seqs, self.hashed_seqs):
            if self.extras:
                extra_item = next(extra_iter)
            hit_ids = prefix_dict[len(query_seq)][query_hash]
            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_dict[query_id] = hit_ids
            # # REMOVE
            # if query_id == 'c_000000064951':
            #     print('bob')
            #     print(hit_ids)
            if self.extras:
                member_item = (query_id, len(query_seq), extra_item)
            else:
                member_item = (query_id, len(query_seq))
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id in hit_ids:
                try:
                    cluster_dict[seed_id][1].append(member_item)
                except KeyError: # new cluster
                    # # REMOVE
                    # if query_id == 'c_000000064951':
                    #     print('al')
                    cluster_dict[seed_id] = ['', [member_item]]
                if query_id == seed_id:
                    cluster_dict[seed_id][0] = query_seq

        # Remove clusters with a seed sequence that is a member sequence of another cluster,
        # so that seed sequences only occur in one cluster.
        # This creates an issue with identical sequences that only hit each other and lack prefix subsequences.
        # Consider two identical sequences of this type, A and B.
        # They will form two clusters with seed sequences A and B, respectively,
        # both clusters containing A and B as member sequences.
        # These clusters will not be added to the `clusters` list,
        # since the seeds are member sequences of another cluster.
        # However, one of these clusters should be added to the `clusters` list.
        # A and B should also have entries in memberships, determined below.
        # The first step in the correction occurs in get_memberships,
        # which determines the membership of every query sequence
        # by searching for them as prefix subsequences of seed sequences in the `clusters` list.
        # A and B will not be found as prefix subsequences of any seeds.
        # Say A occurs first in the query list.
        # Then it will be added again as a seed sequence for the search in get_memberships,
        # such that both A and B will be returned in memberships as members of cluster A.
        # A's ID is also returned by get_memberships along with those of other such seeds,
        # so that their clusters, still found in `cluster_dict`, can be added retroactively to `clusters`.
        clusters = []
        for seed_id, value in cluster_dict.items():
            if len(hit_dict[seed_id]) == 1: # seed sequence only hits itself
                clusters.append((seed_id, value[0], sorted(value[1], key=lambda t: -t[1])))
                # # REMOVE
                # if hit_id == 'c_000000064951':
                #     print('cal')
                #     print(clusters[-1])

        # Now that the final clusters have been found,
        # determine the cluster membership of every query sequence.
        memberships, missed_seed_ids = self.get_memberships([cluster[0] for cluster in clusters],
                                                            [cluster[1] for cluster in clusters])

        for missed_seed_id in missed_seed_ids:
            value = cluster_dict[missed_seed_id]
            # The sort shouldn't be necessary,
            # since all of the members should be of the same length in the cluster.
            clusters.append((missed_seed_id, value[0], sorted(value[1], key=lambda t: -t[1])))

        return clusters, memberships