# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes for basic sequence properties and manipulations."""

from collections import OrderedDict
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


class Dereplicator:
    @staticmethod
    def _check_input(ids, seqs, extras=None):
        if len(ids) != len(seqs):
            raise ConfigError("Your input lists were not the same length. "
                              "`ids` had a length of %d, while `seqs` had a length of %d."
                              % (len(ids), len(seqs)))
        if extras:
            if len(extras) != len(seqs):
                raise ConfigError("Your input lists were not the same length. "
                                  "`extras` had a length of %d, while `ids` and `seqs` had a length of %d."
                                  % (len(extras), len(seqs)))


    @staticmethod
    def _hash_seqs(seqs):
        hashed_seqs = [sha224(seq.encode('utf-8')).hexdigest() for seq in seqs]
        return hashed_seqs


    @staticmethod
    def _full_length_dereplicate_with_extra_info(ids, seqs, extras):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        # Search for identical sequences of the same length.
        seq_lengths = [len(seq) for seq in seqs]
        full_dict = {seq_length: {} for seq_length in sorted(set(seq_lengths))} # store clusters here
        extra_iter = iter(extras)
        for seq_length, hashed_seq, seq_id, seq in zip(seq_lengths, hashed_seqs, ids, seqs):
            inner_dict = full_dict[seq_length]
            if hashed_seq in inner_dict:
                member_info = inner_dict[hashed_seq]
                member_info[1].append(seq_id)
                member_info[2].append(next(extra_iter))
            else:
                inner_dict[hashed_seq] = (seq, [seq_id], [next(extra_iter)])

        clusters = []
        for inner_dict in full_dict.values():
            clusters.extend(inner_dict.values())

        # Sort by cluster size and then ID of the first member (seed).
        clusters.sort(key=lambda member_info: (-len(member_info[1]), member_info[1][0]))

        return clusters


    @staticmethod
    def _full_length_dereplicate_without_extra_info(ids, seqs):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        # Search for identical sequences of the same length.
        seq_lengths = [len(seq) for seq in seqs]
        full_dict = {seq_length: {} for seq_length in sorted(set(seq_lengths))} # store clusters here
        for seq_length, hashed_seq, seq_id, seq in zip(seq_lengths, hashed_seqs, ids, seqs):
            inner_dict = full_dict[seq_length]
            if hashed_seq in inner_dict:
                inner_dict[hashed_seq][1].append(seq_id)
            else:
                inner_dict[hashed_seq] = (seq, [seq_id])

        clusters = []
        for inner_dict in full_dict.values():
            clusters.extend(inner_dict.values())

        # Sort by cluster size and then ID of the first member (seed).
        clusters.sort(key=lambda member_info: (-len(member_info[1]), member_info[1][0]))

        return clusters


    @staticmethod
    def full_length_dereplicate(ids, seqs, extras=None):
        """
        Identical sequences are clustered, returning a list of tuples.

        If no extra sequence information (`extras`) is provided, the list format is:
        clusters = [(<seed seq A>, [<IDs of identical seqs, starting with A>]),
                    (<seed seq B>, [<IDs of identical seqs, starting with B>]), ...]
        If extra information is provided:
        clusters = [(<seed seq A>, [<IDs of identical seqs, starting with A>], [<extra info for each identical seq, starting with A>]),
                    (<seed seq B>, [<IDs of identical seqs, starting with B>], [<extra info for each identical seq, starting with B>]), ...]
        """

        if extras:
            clusters = Dereplicator._full_length_dereplicate_with_extra_info(ids, seqs, extras)
        else:
            clusters = Dereplicator._full_length_dereplicate_without_extra_info(ids, seqs)

        return clusters


    @staticmethod
    def _get_unique_inputs(ids, seqs, extras=None):
        unique_ids = []
        unique_seqs = []
        unique_extras = []
        replicate_dict = {}
        if extras:
            for cluster in Dereplicator.full_length_dereplicate(ids, seqs, extras):
                unique_ids.append(cluster[1][0])
                unique_seqs.append(cluster[0])
                unique_extras.append(cluster[2][0])
                replicate_dict[cluster[1][0]] = list(zip(cluster[1][1: ], cluster[2][1: ]))
        else:
            for cluster in Dereplicator.full_length_dereplicate(ids, seqs):
                unique_ids.append(cluster[1][0])
                unique_seqs.append(cluster[0])
                replicate_dict[cluster[1][0]] = cluster[1][1: ]

        if len(unique_ids) == len(ids):
            unique_ids = ids
            unique_seqs = seqs
            unique_extras = extras
            replicate_dict = None

        return unique_ids, unique_seqs, unique_extras, replicate_dict


    @staticmethod
    def _get_empty_substring_dict(seqs, min_seq_length=None, max_seq_length=None):
        """
        The returned dict has the following structure:
        substring_dict = {min_seq_length: {},
                          min_seq_length + 1: {}, ...
                          max_seq_length: {}}
        By default, `min_seq_length` and `max_seq_length` are determined from `seqs`.
        This is overriden when optional values are provided.
        """
        if min_seq_length:
            if min_seq_length < 1:
                raise ConfigError("The `min_seq_length` argument of `_get_empty_substring_dict` "
                                  "must be at least 1, not %d." % min_seq_length)

        seq_lengths = [len(seq) for seq in seqs]

        if max_seq_length:
            if max_seq_length > max(seq_lengths):
                raise ConfigError("The `max_seq_length` argument of `_get_empty_substring_dict` "
                                  "cannot exceed the length of the longest sequence in the argument, `seqs`. "
                                  "You provided a value of %d, which exceeds the longest length, %d."
                                  % (max_seq_length, max(seq_lengths)))

        min_seq_length = min_seq_length if min_seq_length else min(seq_lengths)
        max_seq_length = max_seq_length if max_seq_length else max(seq_lengths)

        substring_dict = {substring_length: {} for substring_length in range(min_seq_length, max_seq_length + 1)}

        return substring_dict


    @staticmethod
    def _get_prefix_dict(ids, seqs, min_seq_length=None, max_seq_length=None):
        """
        The returned `prefix_dict` contains hashes of prefix subsequences from all of the target sequences in `seqs`.
        `prefix_dict` is keyed by sequence length, with lengths spanning from the min to max sequence length;
        The parameters `min_seq_length` and `max_seq_length` can be used to adjust the lengths considered inward.
        `prefix_dict` values are dictionaries themselves.
        Inner dict keys are prefix sequence hashes, and values are lists of target sequence IDs containing those hashes.
        """

        prefix_dict = Dereplicator._get_empty_substring_dict(
            seqs, min_seq_length=min_seq_length, max_seq_length=max_seq_length)

        for prefix_length, inner_dict in prefix_dict.items():
            for seq_id, seq in zip(ids, seqs):
                # Ignore prefixes of the same length as the sequence.
                if prefix_length >= len(seq):
                    continue
                hashed_prefix = sha224(seq[: prefix_length].encode('utf-8')).hexdigest()
                if hashed_prefix in inner_dict:
                    inner_dict[hashed_prefix].append((seq_id, len(seq)))
                else:
                    inner_dict[hashed_prefix] = [(seq_id, len(seq))]

        for inner_dict in prefix_dict.values():
            for hashed_seq, parent_seqs in inner_dict.items():
                # Sort in descending order of subject sequence length then in ascending order of ID.
                inner_dict[hashed_seq] = [seq_id for seq_id, seq_length
                                          in sorted(parent_seqs, key=lambda t: (-t[1], t[0]))]

        return prefix_dict


    @staticmethod
    def _prefix_dereplicate_with_replicate_seqs_and_with_extra_info(ids, seqs, extras, replicate_dict):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        prefix_dict = Dereplicator._get_prefix_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length, seq X extra),
        #                                 (seq Y ID, seq Y length, seq Y extra),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length, seq X extra),
        #                                 (seq Y ID, seq Y length, seq Y extra),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        extra_iter = iter(extras)
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            query_extra_item = next(extra_iter)

            inner_dict = prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hit_ids = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_dict[query_id] = hit_ids

            # Record information on the query replicates.
            member_info = [(query_id, len(query_seq), query_extra_item)]
            for replicate_id, replicate_extra_item in replicate_dict[query_id]:
                member_info.append((replicate_id, len(query_seq), replicate_extra_item))
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id in hit_ids:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].extend(member_info)
                else:
                    cluster_dict[seed_id] = member_info

        clusters = []
        memberships = []
        for query_id, query_seq, query_extra_item in zip(ids, seqs, extras):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            replicate_items = replicate_dict[query_id]
            cluster_items_for_seed = [(query_id, len(query_seq), query_extra_item)]
            # The seed seq and identical seqs are added as members of the seed's cluster.
            for replicate_id, replicate_extra_item in replicate_items:
                cluster_items_for_seed.append((replicate_id, len(query_seq), replicate_extra_item))

            if query_id in cluster_dict:
                # Other query sequences are prefixes of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq and identical seqs are added as members of the seed's cluster.
                member_info = cluster_items_for_seed + member_info
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed
                # and the query and its replicates as members of the cluster.
                clusters.append((query_id, query_seq, cluster_items_for_seed))

            # Seed and identical sequences are the only members of the seed's cluster.
            memberships.append((query_id, [query_id], query_extra_item))
            for replicate_id, replicate_extra_item in replicate_items:
                memberships.append((replicate_id, [query_id], replicate_extra_item))

        # Record memberships of query sequences that were not seed sequences.

        # Find the prefix hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed prefix subsequences,
        # which only requires searching one list of hashes representing prefix sequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_prefix_dict = Dereplicator._get_prefix_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash, query_extra_item in zip(ids, seqs, hashed_seqs, extras):
            inner_dict = seed_prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            seed_ids = inner_dict[query_hash]
            memberships.append((query_id, seed_ids, query_extra_item))
            for replicate_id in replicate_dict[query_id]:
                memberships.append((replicate_id, seed_ids, query_extra_item))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def _prefix_dereplicate_with_replicate_seqs_and_without_extra_info(ids, seqs, replicate_dict):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        prefix_dict = Dereplicator._get_prefix_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length),
        #                                 (seq Y ID, seq Y length),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length),
        #                                 (seq Y ID, seq Y length),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hit_ids = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_dict[query_id] = hit_ids

            # Record information on the query replicates.
            member_info = [(query_id, len(query_seq))]
            for replicate_id in replicate_dict[query_id]:
                member_info.append((replicate_id, len(query_seq)))
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id in hit_ids:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].extend(member_info)
                else:
                    cluster_dict[seed_id] = member_info

        clusters = []
        memberships = []
        for query_id, query_seq in zip(ids, seqs):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            replicate_ids = replicate_dict[query_id]
            cluster_items_for_seed = [(query_id, len(query_seq))]
            # The seed seq and identical seqs are added as members of the seed's cluster.
            for replicate_id in replicate_ids:
                cluster_items_for_seed.append((replicate_id, len(query_seq)))

            if query_id in cluster_dict:
                # Other query sequences are prefixes of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq and identical seqs are added as members of the seed's cluster.
                member_info = cluster_items_for_seed + member_info
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed
                # and the query and its replicates as members of the cluster.
                clusters.append((query_id, query_seq, cluster_items_for_seed))

            # Seed and identical sequences are the only members of the seed's cluster.
            memberships.append((query_id, [query_id]))
            for replicate_id in replicate_ids:
                memberships.append((replicate_id, [query_id]))

        # Record memberships of query sequences that were not seed sequences.

        # Find the prefix hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed prefix subsequences,
        # which only requires searching one list of hashes representing prefix sequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_prefix_dict = Dereplicator._get_prefix_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = seed_prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            seed_ids = inner_dict[query_hash]
            memberships.append((query_id, seed_ids))
            for replicate_id in replicate_dict[query_id]:
                memberships.append((replicate_id, seed_ids))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def _prefix_dereplicate_without_replicate_seqs_and_with_extra_info(ids, seqs, extras):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        prefix_dict = Dereplicator._get_prefix_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length, seq X extra),
        #                                 (seq Y ID, seq Y length, seq Y extra),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length, seq X extra),
        #                                 (seq Y ID, seq Y length, seq Y extra),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        extra_iter = iter(extras)
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            query_extra_item = next(extra_iter)
            inner_dict = prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hit_ids = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_dict[query_id] = hit_ids
            member_item = (query_id, len(query_seq), query_extra_item)
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id in hit_ids:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].append(member_item)
                else:
                    cluster_dict[seed_id] = [member_item]

        clusters = []
        memberships = []
        # The query sequences were unique.
        for query_id, query_seq, query_extra_item in zip(ids, seqs, extras):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            if query_id in cluster_dict:
                # Other query sequences are prefixes of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq is added as a member of its own cluster.
                member_info.insert(0, (query_id, len(query_seq), query_extra_item))
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed and the sole member of the cluster.
                member_info = [(query_id, len(query_seq), query_extra_item)]
                clusters.append((query_id, query_seq, member_info))

            # Seed sequences are only members of their own cluster.
            memberships.append((query_id, [query_id], query_extra_item))

        # Record memberships of query sequences that were not seed sequences.

        # Find the prefix hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed prefix subsequences,
        # which only requires searching one list of hashes representing prefix sequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_prefix_dict = Dereplicator._get_prefix_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash, extra_query_item in zip(ids, seqs, hashed_seqs, extras):
            inner_dict = seed_prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue
            seed_ids = inner_dict[query_hash]
            memberships.append((query_id, seed_ids, extra_query_item))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def _prefix_dereplicate_without_replicate_seqs_and_without_extra_info(ids, seqs):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        prefix_dict = Dereplicator._get_prefix_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length),
        #                                 (seq Y ID, seq Y length),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length),
        #                                 (seq Y ID, seq Y length),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hit_ids = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_dict[query_id] = hit_ids
            member_item = (query_id, len(query_seq))
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id in hit_ids:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].append(member_item)
                else: # new cluster
                    cluster_dict[seed_id] = [member_item]

        clusters = []
        memberships = []
        # The query sequences were unique.
        for query_id, query_seq in zip(ids, seqs):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            if query_id in cluster_dict:
                # Other query sequences are prefixes of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq is added as a member of its own cluster.
                member_info.insert(0, (query_id, len(query_seq)))
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed and the sole member of the cluster.
                member_info = [(query_id, len(query_seq))]
                clusters.append((query_id, query_seq, member_info))

            # Seed sequences are only members of their own cluster.
            memberships.append((query_id, [query_id]))

        # Record memberships of query sequences that were not seed sequences.

        # Find the prefix hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed prefix subsequences,
        # which only requires searching one list of hashes representing prefix sequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_prefix_dict = Dereplicator._get_prefix_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = seed_prefix_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue
            seed_ids = inner_dict[query_hash]
            memberships.append((query_id, seed_ids))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def prefix_dereplicate(ids, seqs, extras=None):
        """
        Sequences matching the beginning of a longer (seed) sequence are clustered.
        Sequences can occur in multiple clusters, as they may be subsequences of distinct seed sequences, e.g.:
        Cluster 1:
        ACGTACGTACGTACGT (seed, seq A)
        ACGTACGTACGT (seq X)
        ACGTACGT (seq Y)
        Cluster 2:
        ACGTACGTACGTACGG (seed, seq B)
        ACGTACGTACGT (seq X)
        ACGTACGT (seq Y)

        Two lists are returned, a `clusters` list and a `memberships` list.

        Here is an example of the format of the `clusters` list when extra sequence information is NOT provided:
        clusters = [(seed seq A ID, seed seq A, [(seq A ID, seq A length),
                                                 (seq X ID, seq X length),
                                                 (seq Y ID, seq Y length)]),
                    (seed seq B ID, seed seq B, [(seq B ID, seq B length),
                                                 (seq X ID, seq X length),
                                                 (seq Y ID, seq Y length)]),
                    ...]
        Here is an example of the format when extra sequence information is provided:
        clusters = [(seed seq A ID, seed seq A, [(seq A ID, seq A length, seq A extra info),
                                                 (seq X ID, seq X length, seq X extra info),
                                                 (seq Y ID, seq Y length, seq Y extra info)]),
                    (seed seq B ID, seed seq B, [(seq B ID, seq B length, seq B extra info),
                                                 (seq X ID, seq X length, seq X extra info),
                                                 (seq Y ID, seq Y length, seq Y extra info)]),
                    ...]

        The membership list has an entry for each input sequence.
        Here is an example of the format when extra sequence information is NOT provided:
        memberships = [(seq A, [seq A ID, ...]),
                       (seq X, [seq A ID, seq B ID, ...]),
                       (seq Y, [seq A ID, seq B ID, ...]),
                       (seq B, [seq B ID, ...]),
                       ...]
        Here is an example when extra sequence information is provided:
        memberships = [(seq A, [seq A ID, ...], seq A extra info),
                       (seq X, [seq A ID, seq B ID, ...], seq X extra info),
                       (seq Y, [seq A ID, seq B ID, ...], seq Y extra info),
                       (seq B, [seq B ID, ...], seq B extra info),
                       ...]
        """

        Dereplicator._check_input(ids, seqs, extras)

        unique_ids, unique_seqs, unique_extras, replicate_dict = Dereplicator._get_unique_inputs(ids, seqs, extras)

        if replicate_dict:
            if extras:
                clusters, memberships = Dereplicator._prefix_dereplicate_with_replicate_seqs_and_with_extra_info(
                    unique_ids, unique_seqs, unique_extras, replicate_dict)
            else:
                clusters, memberships = Dereplicator._prefix_dereplicate_with_replicate_seqs_and_without_extra_info(
                    unique_ids, unique_seqs, replicate_dict)
        else:
            if extras:
                clusters, memberships = Dereplicator._prefix_dereplicate_without_replicate_seqs_and_with_extra_info(
                    ids, seqs, extras)
            else:
                clusters, memberships = Dereplicator._prefix_dereplicate_without_replicate_seqs_and_without_extra_info(
                    ids, seqs)

        return clusters, memberships


    @staticmethod
    def _get_subseq_dict(ids, seqs, min_seq_length=None, max_seq_length=None):
        """
        The returned `subseq_dict` contains hashes of subseqs from all of the target sequences in `seqs`.
        `subseq_dict` is keyed by sequence length, with lengths spanning from the min to max sequence length;
        The parameters `min_seq_length` and `max_seq_length` can be used to adjust the lengths considered inward.
        `subseq_dict` values are dictionaries themselves.
        Inner dict keys are subseq hashes, and values are lists of target sequence IDs containing those hashes.
        """

        subseq_dict = Dereplicator._get_empty_substring_dict(
            seqs, min_seq_length=min_seq_length, max_seq_length=max_seq_length)

        for subseq_length, inner_dict in subseq_dict.items():
            for seq_id, seq in zip(ids, seqs):
                # Ignore subseqs of the same length as the sequence.
                if subseq_length >= len(seq):
                    continue
                hashed_subseqs = []
                prelim_subseq_items = []
                for start_index, stop_index in zip(range(0, len(seq) - subseq_length + 1),
                                                   range(subseq_length, len(seq) + 1)):
                    hashed_subseq = sha224(seq[start_index: stop_index].encode('utf-8')).hexdigest()
                    hashed_subseqs.append(hashed_subseq)
                    prelim_subseq_items.append([seq_id, [start_index], len(seq)])

                if len(set(hashed_subseqs)) < len(hashed_subseqs):
                    # Consolidate identical subseqs of the same sequence.
                    prev_hashed_subseq = hashed_subseqs[0]
                    prev_subseq_item = prelim_subseq_items[0]
                    for hashed_subseq, prelim_item in sorted(zip(hashed_subseqs, prelim_subseq_items),
                                                             key=lambda t: t[0])[1: ]:
                        if hashed_subseq == prev_hashed_subseq:
                            # Add a start index for an identical subsequence within the target sequence.
                            prev_subseq_item[1] += prelim_item[1]
                        else:
                            if prev_hashed_subseq in inner_dict:
                                inner_dict[prev_hashed_subseq].append(prev_subseq_item)
                            else:
                                inner_dict[prev_hashed_subseq] = [prev_subseq_item]
                            prev_subseq_item = prelim_item
                        prev_hashed_subseq = hashed_subseq
                    # Add the last subseq.
                    if prev_hashed_subseq in inner_dict:
                        inner_dict[prev_hashed_subseq].append(prev_subseq_item)
                    else:
                        inner_dict[prev_hashed_subseq] = [prev_subseq_item]
                else:
                    for hashed_subseq, subseq_item in zip(hashed_subseqs, prelim_subseq_items):
                        if hashed_subseq in inner_dict:
                            inner_dict[hashed_subseq].append(subseq_item)
                        else:
                            inner_dict[hashed_subseq] = [subseq_item]

        for inner_dict in subseq_dict.values():
            for hashed_seq, parent_seqs in inner_dict.items():
                # Sort in descending order of subject sequence length then in ascending order of ID.
                inner_dict[hashed_seq] = [(seq_id, seq_start_indices) for seq_id, seq_start_indices, seq_length
                                          in sorted(parent_seqs, key=lambda t: (-t[2], t[0]))]

        return subseq_dict


    @staticmethod
    def _subseq_dereplicate_with_replicate_seqs_and_with_extra_info(ids, seqs, extras, replicate_dict):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        subseq_dict = Dereplicator._get_subseq_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length, [<seq X start indices in seq A>], seq X extra),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq A>], seq Y extra),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length, [<seq X start indices in seq B>], seq X extra),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq B>], seq Y extra),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        extra_iter = iter(extras)
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            query_extra_item = next(extra_iter)
            inner_dict = subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hits = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a subsequence.
            # A query sequence may hit a target sequence at different indices, producing multiple hits.
            hit_dict[query_id] = hits
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id, seed_start_indices in hits:
                if seed_id in cluster_dict:
                    member_info = cluster_dict[seed_id]
                    member_info.append((query_id, len(query_seq), seed_start_indices, query_extra_item))
                else:
                    member_info = [(query_id, len(query_seq), seed_start_indices, query_extra_item)]
                    cluster_dict[seed_id] = member_info
                # Record information on the query replicates.
                for replicate_id, replicate_extra_item in replicate_dict[query_id]:
                    member_info.append((replicate_id, len(query_seq), seed_start_indices, replicate_extra_item))

        clusters = []
        memberships = []
        for query_id, query_seq, query_extra_item in zip(ids, seqs, extras):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            replicate_items = replicate_dict[query_id]
            member_info_for_seed = [(query_id, len(query_seq), [0], query_extra_item)]
            # The seed seq and identical seqs are added as members of the seed's cluster.
            for replicate_id, replicate_extra_item in replicate_items:
                member_info_for_seed.append((replicate_id, len(query_seq), [0], replicate_extra_item))

            if query_id in cluster_dict:
                # Other query sequences are subseqs of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq and identical seqs are added as members of the seed's cluster.
                member_info = member_info_for_seed + member_info
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed
                # and the query and its replicates as members of the cluster.
                member_info = member_info_for_seed
                clusters.append((query_id, query_seq, member_info))

            # Seed and identical sequences are the only members of the seed's cluster.
            memberships.append((query_id, [query_id], [0], query_extra_item))
            for replicate_id, replicate_extra_item in replicate_items:
                memberships.append((replicate_id, [query_id], [0], replicate_extra_item))

        # Record memberships of query sequences that were not seed sequences.

        # Find the subseq hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed subsequences,
        # which only requires searching one list of hashes representing subsequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_subseq_dict = Dereplicator._get_subseq_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash, query_extra_item in zip(ids, seqs, hashed_seqs, extras):
            inner_dict = seed_subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue
            seed_ids, seed_start_indices = zip(*inner_dict[query_hash])
            memberships.append((query_id, seed_ids, seed_start_indices, query_extra_item))
            for replicate_id in replicate_dict[query_id]:
                memberships.append((replicate_id, seed_ids, seed_start_indices, query_extra_item))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def _subseq_dereplicate_with_replicate_seqs_and_without_extra_info(ids, seqs, replicate_dict):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        subseq_dict = Dereplicator._get_subseq_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length, [<seq X start indices in seq A>]),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq A>]),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length, [<seq X start indices in seq B>]),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq B>]),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hits = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a subsequence.
            # A query sequence may hit a target sequence at different indices, producing multiple hits.
            hit_dict[query_id] = hits
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id, seed_start_indices in hits:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].append((query_id, len(query_seq), seed_start_indices))
                else:
                    cluster_dict[seed_id] = [(query_id, len(query_seq), seed_start_indices)]

        clusters = []
        memberships = []
        for query_id, query_seq in zip(ids, seqs):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            replicate_ids = replicate_dict[query_id]
            member_info_for_seed = [(query_id, len(query_seq), [0])]
            # The seed seq and identical seqs are added as members of the seed's cluster.
            for replicate_id in replicate_ids:
                member_info_for_seed.append((replicate_id, len(query_seq), [0]))

            if query_id in cluster_dict:
                # Other query sequences are subseqs of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq and identical seqs are added as members of the seed's cluster.
                member_info = member_info_for_seed + member_info
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed
                # and the query and its replicates as members of the cluster.
                member_info = member_info_for_seed
                clusters.append((query_id, query_seq, member_info))

            # Seed and identical sequences are the only members of the seed's cluster.
            memberships.append((query_id, [query_id], [0]))
            for replicate_id in replicate_ids:
                memberships.append((replicate_id, [query_id], [0]))

        # Record memberships of query sequences that were not seed sequences.

        # Find the subseq hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed subsequences,
        # which only requires searching one list of hashes representing subsequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_subseq_dict = Dereplicator._get_subseq_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = seed_subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue
            seed_ids, seed_start_indices = zip(*inner_dict[query_hash])
            memberships.append((query_id, seed_ids, seed_start_indices))
            for replicate_id in replicate_dict[query_id]:
                memberships.append((replicate_id, seed_ids, seed_start_indices))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def _subseq_dereplicate_without_replicate_seqs_and_with_extra_info(ids, seqs, extras):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        subseq_dict = Dereplicator._get_subseq_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length, [<seq X start indices in seq A>], seq X extra),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq A>], seq Y extra),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length, [<seq X start indices in seq B>], seq X extra),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq B>], seq Y extra),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        extra_iter = iter(extras)
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            query_extra_item = next(extra_iter)
            inner_dict = subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hits = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a subsequence.
            # A query sequence may hit a target sequence at different indices, producing multiple hits.
            hit_dict[query_id] = hits
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id, seed_start_indices in hits:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].append((query_id, len(query_seq), seed_start_indices, query_extra_item))
                else:
                    cluster_dict[seed_id] = [(query_id, len(query_seq), seed_start_indices, query_extra_item)]

        clusters = []
        memberships = []
        for query_id, query_seq, query_extra_item in zip(ids, seqs, extras):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            # The seed seq and identical seqs are added as members of the seed's cluster.

            if query_id in cluster_dict:
                # Other query sequences are subseqs of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq and identical seqs are added as members of the seed's cluster.
                member_info = [(query_id, len(query_seq), [0], query_extra_item)] + member_info
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed and sole member of the cluster.
                member_info = [(query_id, len(query_seq), [0], query_extra_item)]
                clusters.append((query_id, query_seq, member_info))

            # Seed and identical sequences are the only members of the seed's cluster.
            memberships.append((query_id, [query_id], [0], query_extra_item))

        # Record memberships of query sequences that were not seed sequences.

        # Find the subseq hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed subsequences,
        # which only requires searching one list of hashes representing subsequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_subseq_dict = Dereplicator._get_subseq_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash, query_extra_item in zip(ids, seqs, hashed_seqs, extras):
            inner_dict = seed_subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue
            seed_ids, seed_start_indices = zip(*inner_dict[query_hash])
            memberships.append((query_id, seed_ids, seed_start_indices, query_extra_item))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def _subseq_dereplicate_without_replicate_seqs_and_without_extra_info(ids, seqs):

        hashed_seqs = Dereplicator._hash_seqs(seqs)

        subseq_dict = Dereplicator._get_subseq_dict(ids, seqs)

        # Example `hit_dict` format:
        # hit_dict = {seq X ID: [seq A ID, seq B ID, ...],
        #             seq Y ID: [seq A ID, seq B ID, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format:
        # cluster_dict = {seed seq A ID: [(seq X ID, seq X length, [<seq X start indices in seq A>]),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq A>]),
        #                                 ...]],
        #                 seed seq B ID: [(seq X ID, seq X length, [<seq X start indices in seq B>]),
        #                                 (seq Y ID, seq Y length, [<seq Y start indices in seq B>]),
        #                                 ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}
        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue

            hits = inner_dict[query_hash]
            # Record which target sequences contain the query sequence as a subsequence.
            # A query sequence may hit a target sequence at different indices, producing multiple hits.
            hit_dict[query_id] = hits
            # Make preliminary clusters for each target sequence containing the query sequence.
            for seed_id, seed_start_indices in hits:
                if seed_id in cluster_dict:
                    cluster_dict[seed_id].append((query_id, len(query_seq), seed_start_indices))
                else:
                    cluster_dict[seed_id] = [(query_id, len(query_seq), seed_start_indices)]

        clusters = []
        memberships = []
        for query_id, query_seq in zip(ids, seqs):
            if query_id in hit_dict:
                # The query is not a seed because it is a prefix of other sequences.
                continue

            # The seed seq and identical seqs are added as members of the seed's cluster.

            if query_id in cluster_dict:
                # Other query sequences are subseqs of this query.
                member_info = cluster_dict[query_id]
                # Sort members in descending order of sequence length.
                member_info.sort(key=lambda member_item: -member_item[1])
                # The seed seq and identical seqs are added as members of the seed's cluster.
                member_info = [(query_id, len(query_seq), [0])] + member_info
                clusters.append((query_id, query_seq, member_info))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed and sole member of the cluster.
                member_info = [(query_id, len(query_seq), [0])]
                clusters.append((query_id, query_seq, member_info))

            # Seed and identical sequences are the only members of the seed's cluster.
            memberships.append((query_id, [query_id], [0]))

        # Record memberships of query sequences that were not seed sequences.

        # Find the subseq hashes from the final cluster seed sequences.
        # It is more efficient to search hashed input sequences against these hashed subsequences,
        # which only requires searching one list of hashes representing subsequences of a certain length,
        # than to search for each input sequence ID in every cluster's list of member sequence IDs.
        seed_ids, seed_seqs, _ = zip(*clusters)
        seed_subseq_dict = Dereplicator._get_subseq_dict(seed_ids, seed_seqs, min_seq_length=min([len(seq) for seq in seqs]))

        for query_id, query_seq, query_hash in zip(ids, seqs, hashed_seqs):
            inner_dict = seed_subseq_dict[len(query_seq)]
            if query_hash not in inner_dict:
                continue
            seed_ids, seed_start_indices = zip(*inner_dict[query_hash])
            memberships.append((query_id, seed_ids, seed_start_indices))

        # Sort by cluster/membership size and then by seed/member ID.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


    @staticmethod
    def subseq_dereplicate(ids, seqs, extras=None):
        """
        Sequences matching one or more subsequences of a longer (seed) sequence are clustered.
        Sequences can occur in multiple clusters, as they may be subsequences of distinct seed sequences, e.g.:
        Cluster 1:
        ACGTACGTACGTACGT (seed, seq A)
         CGTACGTACGT     (seq X, start index 1)
             CGTACGTACGT (seq X, start index 5)
        ACGTACGTACGTACG  (seq Y, start index 0)
        Cluster 2:
        ACGTACGTACGTACGG (seed, seq B)
         CGTACGTACGT     (seq X, start index 1)
        ACGTACGTACGTACG  (seq Y, start index 0)

        Two lists are returned, a `clusters` list and a `memberships` list.

        Here is an example of the format of the `clusters` list when extra sequence information is NOT provided:
        clusters = [(seed seq A ID, seed seq A, [(seq A ID, seq A length, [seq A start position in A = 0]),
                                                 (seq X ID, seq X length, [seq X start position in A = 1, seq X start position in A = 5]),
                                                 (seq Y ID, seq Y length, [seq Y start position in A = 0])]),
                    (seed seq B ID, seed seq B, [(seq B ID, seq B length, [seq B start position in B = 0]),
                                                 (seq X ID, seq X length, [seq X start position in B = 1]),
                                                 (seq Y ID, seq Y length, [seq Y start position in B = 0])]),
                    ...]
        Here is an example of the format when extra sequence information is provided:
        clusters = [(seed seq A ID, seed seq A, [(seq A ID, seq A length, [seq A start position in A = 0], seq A extra info),
                                                 (seq X ID, seq X length, [seq X start position in A = 1, seq X start position in A = 5], seq X extra info),
                                                 (seq Y ID, seq Y length, [seq Y start position in A = 0], seq Y extra info)]),
                    (seed seq B ID, seed seq B, [(seq B ID, seq B length, [seq B start position in B = 0], seq B extra info),
                                                 (seq X ID, seq X length, [seq X start position in B = 1], seq X extra info),
                                                 (seq Y ID, seq Y length, [seq Y start position in B = 0], seq Y extra info)]),
                    ...]

        The membership list has an entry for each input sequence.
        Here is an example of the format when extra sequence information is NOT provided:
        memberships = [(seq A, [seq A ID, ...], [[seq A start position in A = 0], ...]),
                       (seq X, [seq A ID, seq B ID, ...], [[seq X start position in A = 1, seq X start position in A = 5], (seq X start position in B = 1], ...]),
                       (seq Y, [seq A ID, seq B ID, ...], [[seq Y start position in A = 0], [seq Y start position in A = 0], ...]),
                       (seq B, [seq B ID, ...], [[seq B start position in B = 0], ...]),
                       ...]
        Here is an example when extra sequence information is provided:
        memberships = [(seq A, [seq A ID, ...], [[seq A start position in A = 0], ...], seq A extra info),
                       (seq X, [seq A ID, seq B ID, ...], [[seq X start position in A = 1, seq X start position in A = 5], [seq X start position in B = 1], ...], seq X extra info),
                       (seq Y, [seq A ID, seq B ID, ...], [[seq Y start position in A = 0], [seq Y start position in A = 0], ...], seq Y extra info),
                       (seq B, [seq B ID, ...], [[seq B start position in B = 0], ...], seq B extra info),
                       ...]
        """

        Dereplicator._check_input(ids, seqs, extras)

        unique_ids, unique_seqs, unique_extras, replicate_dict = Dereplicator._get_unique_inputs(ids, seqs, extras)

        if replicate_dict:
            if extras:
                clusters, memberships = Dereplicator._subseq_dereplicate_with_replicate_seqs_and_with_extra_info(
                    unique_ids, unique_seqs, unique_extras, replicate_dict)
            else:
                clusters, memberships = Dereplicator._subseq_dereplicate_with_replicate_seqs_and_without_extra_info(
                    unique_ids, unique_seqs, replicate_dict)
        else:
            if extras:
                clusters, memberships = Dereplicator._subseq_dereplicate_without_replicate_seqs_and_with_extra_info(
                    ids, seqs, extras)
            else:
                clusters, memberships = Dereplicator._subseq_dereplicate_without_replicate_seqs_and_without_extra_info(
                    ids, seqs)

        return clusters, memberships
