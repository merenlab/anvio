# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes for basic sequence properties and manipulations."""

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

    def __init__(self, id_list, seq_list, extra_list=None):
        if len(id_list) != len(seq_list):
            raise ConfigError(
                "SequenceDereplicator takes a list of sequence IDs and a list of sequence strings. "
                "Your input lists were not the same length. "
                "The ID list had a length of %d, and the sequence list had a length of %d."
                % (len(id_list), len(seq_list)))
        if extra_list:
            if len(extra_list) != len(seq_list):
                raise ConfigError(
                    "SequenceDereplicator takes an optional extra list of sequence information. "
                    "Your extra list was not the same length as the ID and sequence lists. "
                    "The extra list had a length of %d, while the others had a length of %d."
                    % (len(extra_list), len(seq_list)))

        self.id_list = id_list
        self.seq_list = seq_list
        self.hashed_seq_list = [hash(seq) for seq in self.seq_list]
        if extra_list:
            self.extra_list = extra_list
        else:
            self.extra_list = []


    def full_length_dereplicate(self):
        """
        Identical sequences are clustered, returning a list of tuples.

        If no extra sequence information (self.extra_list) is provided, the list format is:
        cluster_list = [(sequence A, [member sequence X ID, member sequence Y ID, ...]), (sequence B, [member IDs]), ...]
        If extra information is provided:
        cluster_list = [(sequence A, [member sequence X ID, ...], [member sequence X extra info, ...]), ...]
        """

        # Search for identical sequences of the same length.
        seq_length_list = [len(seq) for seq in self.seq_list]
        full_dict = {seq_length: {} for seq_length in sorted(set(seq_length_list))} # store clusters here
        extra_iter = iter(self.extra_list)
        for seq_length, hashed_seq, seq_id, seq in zip(seq_length_list, self.hashed_seq_list, self.id_list, self.seq_list):
            inner_dict = full_dict[seq_length]
            try:
                if self.extra_list:
                    value = inner_dict[hashed_seq]
                    value[1].append(seq_id)
                    value[2].append(next(extra_iter))
                else:
                    inner_dict[hashed_seq][1].append(seq_id)
            except KeyError:
                if self.extra_list:
                    inner_dict[hashed_seq] = (seq, [seq_id], [next(extra_iter)])
                else:
                    inner_dict[hashed_seq] = (seq, [seq_id])

        cluster_list = []
        for inner_dict in full_dict.values():
            cluster_list.extend(inner_dict.values())
        cluster_list.sort(key=lambda value: -len(value[1]))

        return cluster_list


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

        Two lists are returned, a cluster list and a cluster membership list.

        Here is the format of the cluster list when extra sequence information is not provided:
        cluster_list = [(seed seq A ID, seed sequence A, [(member seq A ID, member seq A length), (member seq X ID, member seq X length), ...]), ...]
        Here is the format when it is provided:
        cluster_list = [(seed seq A ID, seed sequence A, [(member seq A ID, member seq A length, member seq A extra info), ...]), ...]

        The membership list has an entry for each input sequence.
        Here is the format when extra sequence information is not provided:
        membership_list = [(member sequence X, [seed seq A ID, seed seq B ID, ...]), ...]
        Here is the format when it is provided:
        membership_list = [(member sequence X, [seed seq A ID, seed seq B ID, ...], member seq X extra info), ...]
        """

        def get_prefix_dict(id_list, seq_list, min_seq_length=None):
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
                    raise ConfigError("The `min_seq_length` argument of `get_prefix_dict` must be at least 1, not %d." % min_seq_length)

            seq_length_list = [len(seq) for seq in seq_list]
            min_seq_length = min_seq_length if min_seq_length else min(seq_length_list)
            max_seq_length = max(seq_length_list)
            del seq_length_list

            prefix_dict = {prefix_length: {} for prefix_length
                        in range(min_seq_length, max_seq_length + 1)}
            for prefix_length, inner_dict in prefix_dict.items():
                for seq_id, seq in zip(id_list, seq_list):
                    if prefix_length > len(seq):
                        continue
                    hashed_prefix = hash(seq[: prefix_length])
                    try:
                        inner_dict[hashed_prefix].append((seq_id, len(seq)))
                    except KeyError:
                        inner_dict[hashed_prefix] = [(seq_id, len(seq))]
            for inner_dict in prefix_dict.values():
                for hashed_seq, parent_seq_list in inner_dict.items():
                    # Sort the list of IDs containing the prefix subsequence in descending order of sequence length.
                    inner_dict[hashed_seq] = [seq_id for seq_id, seq_length in sorted(parent_seq_list, key=lambda t: -t[1])]

            return prefix_dict


        def get_membership_list(hit_id_list, hit_seq_list):
            """
            This function is called after clustering to find which sequences are in which clusters.
            Prefix dereplication can result in the same sequence occurring in multiple clusters.
            """

            # Find the prefix hashes from the final cluster seed sequences.
            # It is more efficient to search hashed input sequences against these hashed prefix subsequences,
            # which only requires searching one list of hashes
            # representing prefix sequences of a certain length,
            # than to search for each input sequence ID in every cluster's list of member sequence IDs.
            hit_prefix_dict = get_prefix_dict(hit_id_list, hit_seq_list, min_seq_length=min([len(seq) for seq in self.seq_list]))

            membership_list = []
            extra_iter = iter(self.extra_list)
            for query_id, query_seq, query_hash in zip(self.id_list, self.seq_list, self.hashed_seq_list):
                hit_id_list = hit_prefix_dict[len(query_seq)][query_hash]
                if self.extra_list:
                    membership_list.append((query_id, hit_id_list, next(extra_iter)))
                else:
                    membership_list.append((query_id, hit_id_list))

            return membership_list


        prefix_dict = get_prefix_dict(self.id_list, self.seq_list)

        hit_dict = {}
        cluster_dict = {}
        extra_iter = iter(self.extra_list)
        for query_id, query_seq, query_hash in zip(self.id_list, self.seq_list, self.hashed_seq_list):
            if self.extra_list:
                extra_item = next(extra_iter)
            hit_id_list = prefix_dict[len(query_seq)][query_hash]
            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_dict[query_id] = hit_id_list
            if self.extra_list:
                member_item = (query_id, len(query_seq), extra_item)
            else:
                member_item = (query_id, len(query_seq))
            # Make preliminary clusters for each target sequence containing the query sequence.
            for hit_id in hit_id_list:
                try:
                    cluster_dict[hit_id][1].append(member_item)
                except KeyError: # new cluster
                    cluster_dict[hit_id] = ['', [member_item]]
                if query_id == hit_id:
                    cluster_dict[hit_id][0] = query_seq

        # Remove clusters with seed sequences that are member sequences of another cluster.
        # Seed sequences should only hit themselves.
        cluster_list = []
        for hit_id, value in cluster_dict.items():
            if len(hit_dict[hit_id]) == 1: # seed sequence only hits itself
                cluster_list.append((hit_id, value[0], sorted(value[1], key=lambda t: -t[1])))
        del hit_dict

        # Now that the final clusters have been found,
        # determine the cluster membership of every query sequence.
        membership_list = get_membership_list([cluster[0] for cluster in cluster_list], [cluster[1] for cluster in cluster_list])

        return cluster_list, membership_list


