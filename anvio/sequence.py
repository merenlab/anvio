# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes for basic sequence properties and manipulations."""

import functools
import itertools
import multiprocessing

from collections import OrderedDict
from copy import deepcopy
from hashlib import sha224

import anvio
import anvio.constants as constants
import anvio.terminal as terminal

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
        index_trajectories = list(itertools.permutations(indices_of_variation))

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


def get_kmers(name_seq_pair, kmer_size, include_full_length=True):
    """Get information on each k-mer of a certain size in the input sequence.

    This function is located outside `Kmerizer` to allow multithreading.

    Parameters
    ==========
    name_seq_pair : tuple
        A length two tuple of the sequence name string and sequence string

    kmer_size : int
        Length of the k-mers to be extracted

    """
    name, seq = name_seq_pair

    if include_full_length:
        # Include input seqs as k-mers.
        if kmer_size > len(seq):
            return []
    elif kmer_size >= len(seq):
        # Do not include input seqs as k-mers.
        return []

    hashed_kmers = []

    # List every occurrence of every k-mer in the sequence before consolidating items for identical k-mers.
    prelim_kmer_items = []
    for start_index, stop_index in zip(range(0, len(seq) - kmer_size + 1), range(kmer_size, len(seq) + 1)):
        hashed_kmer = sha224(seq[start_index: stop_index].encode('utf-8')).hexdigest()
        hashed_kmers.append(hashed_kmer)
        prelim_kmer_items.append([hashed_kmer, name, [start_index], len(seq)])

    kmer_items = []
    if len(set(hashed_kmers)) < len(hashed_kmers):
        # Consolidate identical k-mers drawn from the same sequence.
        prelim_kmer_items.sort(key=lambda prelim_kmer_item: prelim_kmer_item[0])
        prev_kmer_item = prelim_kmer_items[0]
        prev_hashed_kmer = prev_kmer_item[0]
        for prelim_kmer_item in prelim_kmer_items[1: ]:
            hashed_kmer = prelim_kmer_item[0]
            if hashed_kmer == prev_hashed_kmer:
                # Add the start index of the identical k-mer to the list of start indices for the k-mer.
                prev_kmer_item[2] += prelim_kmer_item[2]
            else:
                prev_kmer_item[2] = tuple(prev_kmer_item[2])
                kmer_items.append(tuple(prev_kmer_item))
                prev_kmer_item = prelim_kmer_item
            prev_hashed_kmer = hashed_kmer
        # Add the last k-mer.
        prev_kmer_item[2] = tuple(prev_kmer_item[2])
        kmer_items.append(tuple(prev_kmer_item))
    else:
        # The sequence did not contain multiple instances of the k-mer.
        for kmer_item in prelim_kmer_items:
            kmer_item[2] = tuple(kmer_item[2])
            kmer_items.append(tuple(kmer_item))

    return kmer_items


class Kmerizer:
    def __init__(self, names, seqs, num_threads=1, progress=terminal.Progress()):
        """Gets hashed k-mers from input sequences.

        Parameters
        ==========
        names : list
            Name strings corresponding to sequences in `seqs`

        seqs : list
            Sequence strings

        num_threads : int, 1
            Threads available for multithreaded operations
        """

        if len(names) != len(seqs):
            raise ConfigError("Your Kmerizer input lists were not the same length. "
                              "`names` had a length of %d, while `seqs` had a length of %d."
                              % (len(names), len(seqs)))

        self.names = names
        self.seqs = seqs
        self.num_threads = num_threads
        self.progress = progress


    def get_prefix_kmer_dict1(self, kmer_size):
        kmer_dict = {}
        for name, seq in zip(self.names, self.seqs):
            if kmer_size > len(seq):
                continue

            hashed_kmer = sha224(seq[: kmer_size].encode('utf-8')).hexdigest()
            if hashed_kmer in kmer_dict:
                kmer_dict[hashed_kmer][name] = seq
            else:
                kmer_dict[hashed_kmer] = {name: seq}

        return kmer_dict


    def get_prefix_kmer_dict(self, kmer_size, include_full_length=True):
        """Get k-mers of a given size from the beginning of input sequences.

        Parameters
        ==========
        kmer_size : int
            Length of the k-mers to be extracted

        include_full_length : bool, True
            Whether to include or ignore full-length input sequences among the hashed k-mers

        Returns
        =======
        kmer_dict : dict
            Keys are unique k-mer strings identified from the input sequences.
            Values are lists of the names of input sequences containing the k-mer.
            These lists are sorted in descending order of sequence length and then alphanumerically by ID.
        """

        kmer_dict = {}
        for name, seq in zip(self.names, self.seqs):
            if include_full_length:
                # Include input seqs as k-mers.
                if kmer_size > len(seq):
                    continue
            elif kmer_size >= len(seq):
                # Do not include input seqs as k-mers.
                continue

            hashed_kmer = sha224(seq[: kmer_size].encode('utf-8')).hexdigest()
            if hashed_kmer in kmer_dict:
                kmer_dict[hashed_kmer].append((name, len(seq)))
            else:
                kmer_dict[hashed_kmer] = [(name, len(seq))]

        for hashed_kmer, parent_seqs in kmer_dict.items():
            kmer_dict[hashed_kmer] = [name for name, seq_length in sorted(parent_seqs, key=lambda t: (-t[1], t[0]))]

        return kmer_dict


    def get_kmer_dict(self, kmer_size, include_full_length=True, sort_kmer_items=False):
        """Get k-mers of a given size from the input sequences.

        Parameters
        ==========
        kmer_size : int
            Length of the k-mers to be extracted

        include_full_length : bool, True
            Whether to include or ignore full-length input sequences among the hashed k-mers

        sort_kmer_items : bool, False

        Returns
        =======
        kmer_dict : dict
            Keys are unique k-mer strings identified from the input sequences.
            Values are tuples of length two tuples.
            The first item in an inner tuple is an input sequence name.
            The second item in an inner tuple is a list of start indices of the k-mer in the input sequence,
            as a given k-mer can occur multiple times in the sequence.
        """

        # self.progress.new("Extracting k-mers of length %d" % kmer_size)

        pool = multiprocessing.Pool(self.num_threads)
        all_seq_kmer_items = []
        # total_seq_count = len(self.names)
        # num_extracted_kmers = 0
        for seq_kmer_item in pool.imap_unordered(functools.partial(get_kmers,
                                                                   kmer_size=kmer_size,
                                                                   include_full_length=include_full_length),
                                                 zip(self.names, self.seqs),
                                                 chunksize=int(len(self.names) / self.num_threads) + 1):
            all_seq_kmer_items.append(seq_kmer_item)
            # num_extracted_kmers += len(seq_kmer_item)
            # if num_processed_seqs % 50000 == 0:
            #     self.progress.update("%d k-mers have been extracted from %d/%d sequences" % (num_extracted_kmers,
            #                                                                                  num_processed_seqs,
            #                                                                                  total_seq_count))
        pool.close()
        pool.join()

        # self.progress.update("Grouping k-mers")
        kmer_dict = {}
        for hashed_kmer, kmer_items in itertools.groupby(sorted([kmer_item
                                                                 for seq_kmer_items in all_seq_kmer_items
                                                                 for kmer_item in seq_kmer_items],
                                                                key=lambda kmer_item: kmer_item[0]),
                                                         key=lambda kmer_item: kmer_item[0]):
            if sort_kmer_items:
                kmer_dict[hashed_kmer] = [(name, seq_start_indices)
                                          for hashed_kmer, name, seq_start_indices, seq_length
                                          in sorted(kmer_items,
                                                    key=lambda kmer_item: (-kmer_item[3], kmer_item[1]))]
            else:
                kmer_dict[hashed_kmer] = [(name, seq_start_indices)
                                          for hashed_kmer, name, seq_start_indices, seq_length
                                          in kmer_items]

        # self.progress.end()

        return kmer_dict


    def _check_kmer_size_bounds(self, min_kmer_size=None, max_kmer_size=None):
        """Checks that the range of k-mer sizes to be extracted makes sense in terms of input sequence length.

        Parameters
        ==========
        min_kmer_size : int
            Min k-mer length to be extracted from input sequences

        max_kmer_size : int
            Max k-mer length to be extracted from input sequences
        """
        if min_kmer_size:
            if min_kmer_size < 1:
                raise ConfigError("`min_kmer_size` in Kmerizer must be at least 1, not %d." % min_kmer_size)

        seq_lengths = [len(seq) for seq in self.seqs]

        if max_kmer_size:
            if max_kmer_size > max(seq_lengths):
                raise ConfigError("`max_kmer_size` in Kmerizer cannot exceed the length of the longest sequence in `seqs`. "
                                  "You provided a value of %d, exceeding the longest length, %d."
                                  % (max_kmer_size, max(seq_lengths)))


    def get_kmer_superdict(self, min_kmer_size=None, max_kmer_size=None, only_prefixes=False, include_full_length=True):
        """Get k-mers of a range of lengths from input sequences.

        Parameters
        ==========
        min_kmer_size : int
            Min k-mer length to be extracted from input sequences

        max_kmer_size : int
            Max k-mer length to be extracted from input sequences,
            which cannot exceed the length of the longest sequence in `seqs`

        only_prefixes : False
            Whether to only retrieve k-mers from the beginning of the input sequences

        include_full_length : bool, True
            Whether to include or ignore full-length input sequences among the hashed k-mers

        Returns
        =======
        kmer_superdict : dict
            This is a dict of dicts.
            The outer dict is keyed by k-mer size, ints spanning [min_kmer_size, max_kmer_size].
            The inner dict values are output from `get_kmer_dict`.
        """

        # By default, min_kmer_size and max_kmer_size are determined from self.seqs.
        self._check_kmer_size_bounds(min_kmer_size, max_kmer_size)

        if min_kmer_size is None or max_kmer_size is None:
            if min_kmer_size is None:
                min_kmer_size = min(map(len, self.seqs))
            if max_kmer_size is None:
                # If `include_full_length` is False, then there will be an empty dict for the max kmer size.
                max_kmer_size = max(map(len, self.seqs))

        kmer_superdict = {}
        for kmer_size in range(min_kmer_size, max_kmer_size + 1):
            if only_prefixes:
                kmer_superdict[kmer_size] = self.get_kmer_dict(kmer_size, include_full_length=include_full_length)
            else:
                kmer_superdict[kmer_size] = self.get_kmer_dict(kmer_size, include_full_length=include_full_length)

        return kmer_superdict


class Dereplicator:
    def __init__(self, names, seqs, extras=None, num_threads=1):
        """This class has methods to dereplicate input sequences in different ways.

        Parameters
        ==========
        names : list
            Name strings corresponding to sequences in `seqs`

        seqs : list
            Sequence strings

        extras : list, None
            Extra items of any type associated with each input sequence

        num_threads : int, 1
            Threads available for multithreaded operations
        """
        if len(names) != len(seqs):
            raise ConfigError("Your Dereplicator input lists were not the same length. "
                              "`names` had a length of %d, while `seqs` had a length of %d."
                              % (len(names), len(seqs)))

        if extras:
            if len(extras) != len(seqs):
                raise ConfigError("Your Dereplicator input lists were not the same length. "
                                  "`extras` had a length of %d, while `names` and `seqs` had a length of %d."
                                  % (len(extras), len(seqs)))

        self.names = names
        self.seqs = seqs
        self.extras = extras
        self.num_threads = num_threads


    def full_length_dereplicate(self):
        """ Cluster identical sequences.

        Returns
        =======
        clusters : list
            There is a tuple for each cluster of identical sequences.
            If there extra sequence information is provided in `self.extras`,
            the tuple has 3 items, otherwise it has 2 items.
            Here is an example of the format of `clusters` WITH extra information:
            clusters = [(<seed seq A>, [<Names of identical seqs, starting with A>]),
                        (<seed seq B>, [<Names of identical seqs, starting with B>]), ...]
            Here is an example of the format of `clusters` WITHOUT extra information:
            clusters = [(<seed seq A>, [<Names of identical seqs, starting with A>], [<extra items for each identical seq, starting with A>]),
                        (<seed seq B>, [<Names of identical seqs, starting with B>], [<extra items for each identical seq, starting with B>]), ...]
            The seed sequence is selected by order of occurrence in `seqs`.
        """

        clusters = []
        if self.extras:
            for seq, group in itertools.groupby(
                sorted(zip(self.seqs, zip(self.names, self.extras)), key=lambda x: x[0]), key=lambda x: x[0]):
                group = list(group)
                member_names = [member[1][0] for member in group]
                member_extras = [member[1][1] for member in group]
                clusters.append((seq, member_names, member_extras))
        else:
            for seq, group in itertools.groupby(
                sorted(zip(self.seqs, self.names), key=lambda x: x[0]), key=lambda x: x[0]):
                group = list(group)
                member_names = [member[1] for member in group]
                clusters.append((seq, member_names))

        # Sort by cluster size and then name of the first member (seed).
        clusters.sort(key=lambda member_info: (-len(member_info[1]), member_info[1][0]))

        return clusters


    def prefix_dereplicate(self):
        """Cluster input sequences matching the beginning of a longer (seed) input sequence.

        Returns
        =======
        clusters : list
        memberships : list

        Sequences can occur in multiple clusters, as they may be subsequences of distinct seed sequences, e.g.:
        Cluster 1:
        ACGTACGTACGTACGT (seed, seq A)
        ACGTACGTACGT     (seq X)
        ACGTACGT         (seq Y)
        Cluster 2:
        ACGTACGTACGTACGG (seed, seq B)
        ACGTACGTACGT     (seq X)
        ACGTACGT         (seq Y)

        Here is an example of the format of the `clusters` list when extra sequence information is NOT provided:
        clusters = [(seed seq A name, seed seq A, [(seq A name, seq A length),
                                                   (seq X name, seq X length),
                                                   (seq Y name, seq Y length)]),
                    (seed seq B name, seed seq B, [(seq B name, seq B length),
                                                   (seq X name, seq X length),
                                                   (seq Y name, seq Y length)]),
                    ...]
        Here is an example of the format when extra sequence information is provided:
        clusters = [(seed seq A name, seed seq A, [(seq A name, seq A length, seq A extra info),
                                                   (seq X name, seq X length, seq X extra info),
                                                   (seq Y name, seq Y length, seq Y extra info)]),
                    (seed seq B name, seed seq B, [(seq B name, seq B length, seq B extra info),
                                                   (seq X name, seq X length, seq X extra info),
                                                   (seq Y name, seq Y length, seq Y extra info)]),
                    ...]

        The `memberships` list has an entry for each input sequence.
        Here is an example of the format when extra sequence information is NOT provided:
        memberships = [(seq A, [seq A name, ...]),
                       (seq X, [seq A name, seq B name, ...]),
                       (seq Y, [seq A name, seq B name, ...]),
                       (seq B, [seq B name, ...]),
                       ...]
        Here is an example when extra sequence information is provided:
        memberships = [(seq A, [seq A name, ...], seq A extra info),
                       (seq X, [seq A name, seq B name, ...], seq X extra info),
                       (seq Y, [seq A name, seq B name, ...], seq Y extra info),
                       (seq B, [seq B name, ...], seq B extra info),
                       ...]
        """
        # self.progress.new("Prefix dereplicating")
        kmer_size = min(map(len, self.seqs))

        # self.progress.update("Extracting prefix k-mers")
        kmer_dict = Kmerizer(self.names, self.seqs).get_prefix_kmer_dict1(kmer_size)

        # Example `hit_dict` format:
        # hit_dict = {seq X name: [seq A name, seq B name, ...],
        #             seq Y name: [seq A name, seq B name, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format (with extra seq info):
        # cluster_dict = {seed seq A name: [(seq X name, seq X length, seq X extra),
        #                                   (seq Y name, seq Y length, seq Y extra),
        #                                   ...]],
        #                 seed seq B name: [(seq X name, seq X length, seq X extra),
        #                                   (seq Y name, seq Y length, seq Y extra),
        #                                   ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        # self.progress.update("Forming dereplicated sequence clusters")
        hashed_prefixes = [sha224(seq[: kmer_size].encode('utf-8')).hexdigest() for seq in self.seqs]

        names_of_queries_with_hits = []
        cluster_dict = {}

        for query_name, query_seq, query_extra_item, query_prefix_hash in zip(self.names, self.seqs, self.extras, hashed_prefixes):
            if query_prefix_hash not in kmer_dict:
                continue

            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_names = []
            for target_name, target_seq in kmer_dict[query_prefix_hash].items():
                if query_seq == target_seq[: len(query_seq)]:
                    if len(query_seq) != len(target_seq):
                        hit_names.append(target_name)

            if hit_names:
                names_of_queries_with_hits.append(query_name)

                member_item = (query_name, len(query_seq), query_extra_item)
                # Make preliminary clusters for each target sequence containing the query sequence.
                for name in hit_names:
                    if name in cluster_dict:
                        cluster_dict[name].append(member_item)
                    else:
                        cluster_dict[name] = [member_item]

        clusters = []
        for query_name, query_seq, query_extra_item in zip(self.names, self.seqs, self.extras):
            if query_name in names_of_queries_with_hits:
                # The query is not a cluster seed because it is a prefix of other sequences.
                continue

            if query_name in cluster_dict:
                # Other query sequences are prefixes of this query.
                member_items = cluster_dict[query_name]
                # Sort members in descending order of sequence length.
                member_items.sort(key=lambda member_item: -member_item[1])
                # The seed seq is added as a member of its own cluster.
                clusters.append((query_name, query_seq, [(query_name, len(query_seq), query_extra_item)] + member_items))
            else:
                # No other query sequences were prefixes of this query.
                # Make a cluster with the query as the seed and the sole member of the cluster.
                clusters.append((query_name, query_seq, [(query_name, len(query_seq), query_extra_item)]))

        # self.progress.update("Sorting clusters")
        # Sort by cluster size and then by seed name.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))

        # self.progress.end()

        return clusters


    def subseq_dereplicate(self):
        """Cluster input sequences matching one or more subsequences of a longer (seed) input sequence.

        Returns
        =======
        clusters : list
        memberships : list

        Sequences can occur in multiple clusters, as they may be subsequences of distinct seed sequences, e.g.:
        Cluster 1:
        ACGTACGTACGTACGT    (seed, seq A)
         CGTACGTACGT        (seq X, start index 1)
             CGTACGTACGT    (seq X, start index 5)
        ACGTACGTACGTACG     (seq Y, start index 0)
        Cluster 2:
        ACGTACGTACGTACGG    (seed, seq B)
         CGTACGTACGT        (seq X, start index 1)
        ACGTACGTACGTACG     (seq Y, start index 0)

        Here is an example of the format of the `clusters` list when extra sequence information is NOT provided:
        clusters = [(seed seq A name, seed seq A, [(seq A name, seq A length, [seq A start position in A = 0]),
                                                   (seq X name, seq X length, [seq X start position in A = 1, seq X start position in A = 5]),
                                                   (seq Y name, seq Y length, [seq Y start position in A = 0])]),
                    (seed seq B name, seed seq B, [(seq B name, seq B length, [seq B start position in B = 0]),
                                                   (seq X name, seq X length, [seq X start position in B = 1]),
                                                   (seq Y name, seq Y length, [seq Y start position in B = 0])]),
                    ...]
        Here is an example of the format when extra sequence information is provided:
        clusters = [(seed seq A name, seed seq A, [(seq A name, seq A length, [seq A start position in A = 0], seq A extra info),
                                                   (seq X name, seq X length, [seq X start position in A = 1, seq X start position in A = 5], seq X extra info),
                                                   (seq Y name, seq Y length, [seq Y start position in A = 0], seq Y extra info)]),
                    (seed seq B name, seed seq B, [(seq B name, seq B length, [seq B start position in B = 0], seq B extra info),
                                                   (seq X name, seq X length, [seq X start position in B = 1], seq X extra info),
                                                   (seq Y name, seq Y length, [seq Y start position in B = 0], seq Y extra info)]),
                    ...]

        The membership list has an entry for each input sequence.
        Here is an example of the format when extra sequence information is NOT provided:
        memberships = [(seq A, [seq A name, ...], [[seq A start position in A = 0], ...]),
                       (seq X, [seq A name, seq B name, ...], [[seq X start position in A = 1, seq X start position in A = 5], (seq X start position in B = 1], ...]),
                       (seq Y, [seq A name, seq B name, ...], [[seq Y start position in A = 0], [seq Y start position in A = 0], ...]),
                       (seq B, [seq B name, ...], [[seq B start position in B = 0], ...]),
                       ...]
        Here is an example when extra sequence information is provided:
        memberships = [(seq A, [seq A name, ...], [[seq A start position in A = 0], ...], seq A extra info),
                       (seq X, [seq A name, seq B name, ...], [[seq X start position in A = 1, seq X start position in A = 5], [seq X start position in B = 1], ...], seq X extra info),
                       (seq Y, [seq A name, seq B name, ...], [[seq Y start position in A = 0], [seq Y start position in A = 0], ...], seq Y extra info),
                       (seq B, [seq B name, ...], [[seq B start position in B = 0], ...], seq B extra info),
                       ...]
        """

        # Get unique input sequences to speed up the subsequence k-mer search.
        unique_names = []
        unique_seqs = []
        unique_extras = []
        replicate_dict = {}
        for cluster in self.full_length_dereplicate():
            unique_name = cluster[1][0]
            unique_names.append(unique_name)
            unique_seqs.append(cluster[0])
            if self.extras:
                unique_extras.append(cluster[2][0])
                replicate_dict[unique_name] = list(zip(cluster[1][1: ], cluster[2][1: ]))
            else:
                replicate_dict[unique_name] = cluster[1][1: ]

        hashed_seqs = [sha224(seq.encode('utf-8')).hexdigest() for seq in self.seqs]

        subseq_dict = Kmerizer(self.names, self.seqs).get_kmer_superdict()

        # Example `hit_dict` format:
        # hit_dict = {seq X name: [seq A name, seq B name, ...],
        #             seq Y name: [seq A name, seq B name, ...],
        #             ...}
        # Seed seqs -- A and B in the example -- do not have entries in `hit_dict`.

        # Example `cluster_dict` format (with extra seq info):
        # cluster_dict = {seed seq A name: [(seq X name, seq X length, [<seq X start indices in seq A>], seq X extra),
        #                                   (seq Y name, seq Y length, [<seq Y start indices in seq A>], seq Y extra),
        #                                   ...]],
        #                 seed seq B name: [(seq X name, seq X length, [<seq X start indices in seq B>], seq X extra),
        #                                   (seq Y name, seq Y length, [<seq Y start indices in seq B>], seq Y extra),
        #                                   ...]],
        #                 ...}
        # Seed seqs -- A and B in the example -- do not have entries in their own clusters in `cluster_dict`.

        hit_dict = {}
        cluster_dict = {}

        if unique_extras:
            extra_iter = iter(unique_extras)
            for query_name, query_seq, query_hash in zip(unique_names, unique_seqs, hashed_seqs):
                query_extra_item = next(extra_iter)

                kmer_dict = subseq_dict[len(query_seq)]
                if query_hash not in kmer_dict:
                    continue

                hits = kmer_dict[query_hash]
                # Record which target sequences contain the query sequence as a subsequence.
                # A query sequence may hit a target sequence at different indices, producing multiple hits.
                hit_dict[query_name] = hits

                # Make preliminary clusters for each target sequence containing the query sequence.
                for seed_name, seed_start_indices in hits:
                    if seed_name in cluster_dict:
                        member_items = cluster_dict[seed_name]
                        member_items.append((query_name, len(query_seq), seed_start_indices, query_extra_item))
                    else:
                        member_items = [(query_name, len(query_seq), seed_start_indices, query_extra_item)]
                        cluster_dict[seed_name] = member_items

                    if replicate_dict:
                        # Record information on the query replicates.
                        for replicate_name, replicate_extra_item in replicate_dict[query_name]:
                            member_items.append((replicate_name, len(query_seq), seed_start_indices, replicate_extra_item))
        else:
            for query_name, query_seq, query_hash in zip(unique_names, unique_seqs, hashed_seqs):
                kmer_dict = subseq_dict[len(query_seq)]
                if query_hash not in kmer_dict:
                    continue

                hits = kmer_dict[query_hash]
                # Record which target sequences contain the query sequence as a subsequence.
                # A query sequence may hit a target sequence at different indices, producing multiple hits.
                hit_dict[query_name] = hits

                # Make preliminary clusters for each target sequence containing the query sequence.
                for seed_name, seed_start_indices in hits:
                    if seed_name in cluster_dict:
                        member_items = cluster_dict[seed_name]
                        member_items.append((query_name, len(query_seq), seed_start_indices))
                    else:
                        member_items = [(query_name, len(query_seq), seed_start_indices)]
                        cluster_dict[seed_name] = member_items

                    if replicate_dict:
                        # Record information on the query replicates.
                        for replicate_name in replicate_dict[query_name]:
                            member_items.append((replicate_name, len(query_seq), seed_start_indices))

        clusters = []
        memberships = []

        if unique_extras:
            for query_name, query_seq, query_extra_item in zip(unique_names, unique_seqs, unique_extras):
                if query_name in hit_dict:
                    # The query is not a seed because it is a prefix of other sequences.
                    continue

                if replicate_dict:
                    replicate_items = replicate_dict[query_name]
                    member_items_for_seed = [(query_name, len(query_seq), [0], query_extra_item)]
                    # The seed seq and identical seqs are added as members of the seed's cluster.
                    for replicate_name, replicate_extra_item in replicate_items:
                        member_items_for_seed.append((replicate_name, len(query_seq), [0], replicate_extra_item))

                if query_name in cluster_dict:
                    # Other query sequences are subseqs of this query.
                    member_items = cluster_dict[query_name]
                    # Sort members in descending order of sequence length.
                    member_items.sort(key=lambda member_item: -member_item[1])
                    # The seed seq and identical seqs are added as members of the seed's cluster.
                    member_items = member_items_for_seed + member_items
                    clusters.append((query_name, query_seq, member_items))
                else:
                    # No other query sequences were prefixes of this query.
                    # Make a cluster with the query as the seed
                    # and the query and its replicates as members of the cluster.
                    member_items = member_items_for_seed
                    clusters.append((query_name, query_seq, member_items))

                # Seed and identical sequences are the only members of the seed's cluster.
                memberships.append((query_name, [query_name], [0], query_extra_item))

                if replicate_dict:
                    for replicate_name, replicate_extra_item in replicate_items:
                        memberships.append((replicate_name, [query_name], [0], replicate_extra_item))
        else:
            for query_name, query_seq in zip(unique_names, unique_seqs):
                if query_name in hit_dict:
                    # The query is not a seed because it is a prefix of other sequences.
                    continue

                if replicate_dict:
                    replicate_names = replicate_dict[query_name]
                    member_items_for_seed = [(query_name, len(query_seq), [0])]
                    # The seed seq and identical seqs are added as members of the seed's cluster.
                    for replicate_name in replicate_names:
                        member_items_for_seed.append((replicate_name, len(query_seq), [0]))

                if query_name in cluster_dict:
                    # Other query sequences are subseqs of this query.
                    member_items = cluster_dict[query_name]
                    # Sort members in descending order of sequence length.
                    member_items.sort(key=lambda member_item: -member_item[1])
                    # The seed seq and identical seqs are added as members of the seed's cluster.
                    member_items = member_items_for_seed + member_items
                    clusters.append((query_name, query_seq, member_items))
                else:
                    # No other query sequences were prefixes of this query.
                    # Make a cluster with the query as the seed
                    # and the query and its replicates as members of the cluster.
                    member_items = member_items_for_seed
                    clusters.append((query_name, query_seq, member_items))

                # Seed and identical sequences are the only members of the seed's cluster.
                memberships.append((query_name, [query_name], [0]))

                if replicate_dict:
                    for replicate_name in replicate_names:
                        memberships.append((replicate_name, [query_name], [0]))

        # Record memberships of query sequences that were not found to be seed sequences.
        # To determine cluster memberships, search hashed input sequences against hashed subsequences of cluster seed sequences.
        # This is more efficient than searching input sequence names against every cluster's list of member sequence names.
        seed_names, seed_seqs, _ = zip(*clusters)
        seed_subseq_dict = Kmerizer(seed_names, seed_seqs).get_kmer_superdict(min_kmer_size=min([len(seq) for seq in unique_seqs]))

        if unique_extras:
            for query_name, query_seq, query_hash, query_extra_item in zip(unique_names, unique_seqs, hashed_seqs, unique_extras):
                kmer_dict = seed_subseq_dict[len(query_seq)]
                if query_hash not in kmer_dict:
                    continue
                seed_names, seed_start_indices = zip(*kmer_dict[query_hash])
                memberships.append((query_name, seed_names, seed_start_indices, query_extra_item))

                if replicate_dict:
                    for replicate_name, replicate_extra_item in replicate_dict[query_name]:
                        memberships.append((replicate_name, seed_names, seed_start_indices, replicate_extra_item))
        else:
            for query_name, query_seq, query_hash in zip(unique_names, unique_seqs, hashed_seqs):
                kmer_dict = seed_subseq_dict[len(query_seq)]
                if query_hash not in kmer_dict:
                    continue
                seed_names, seed_start_indices = zip(*kmer_dict[query_hash])
                memberships.append((query_name, seed_names, seed_start_indices))

                if replicate_dict:
                    for replicate_name in replicate_dict[query_name]:
                        memberships.append((replicate_name, seed_names, seed_start_indices))

        # Sort by cluster/membership size and then by seed/member name.
        clusters.sort(key=lambda cluster: (-len(cluster[2]), cluster[0]))
        memberships.sort(key=lambda membership: (-len(membership[1]), membership[0]))

        return clusters, memberships


class AlignedQuery:
    def __init__(self, seq, name=None):
        self.seq = seq
        self.name = name
        self.alignments = []


    def add_alignment(self, alignment):
        self.alignments.append(alignment)


class AlignedTarget:
    def __init__(self, seq, name=None):
        self.seq = seq
        self.name = name
        self.alignments = []


    def add_alignment(self, alignment):
        self.alignments.append(alignment)


class Alignment:
    def __init__(self, query_start, target_start, cigartuples, aligned_query=None, aligned_target=None):
        self.query_start = query_start
        self.target_start = target_start
        self.cigartuples = cigartuples
        # For now, assume there are no indels in the alignment.
        self.alignment_length = sum(cigartuple[1] for cigartuple in cigartuples)
        self.query_end = self.query_start + self.alignment_length
        self.target_end = self.target_start + self.alignment_length

        self.aligned_query = aligned_query
        self.aligned_target = aligned_target


class Aligner:
    def __init__(self, query_names, query_seqs, target_names, target_seqs, num_threads=1, progress=terminal.Progress()):
        """Align query sequences to target sequences

        Parameters
        ==========
        query_names : list
            name strings corresponding to sequences in `query_seqs`

        query_seqs : list
            Query sequence strings

        target_names : list
            name strings corresponding to sequences in `target_seqs`

        target_seqs : list
            Target sequence strings

        num_threads : int
            Threads available for multithreaded operations

        progress : terminal.Progress object
        """
        if len(query_names) != len(query_seqs):
            raise ConfigError("Your Mapper query input lists were not the same length. "
                              "`query_names` had a length of %d, while `query_seqs` had a length of %d."
                              % (len(query_names), len(query_seqs)))

        if len(target_names) != len(target_seqs):
            raise ConfigError("Your Mapper target input lists were not the same length. "
                              "`target_names` had a length of %d, while `target_seqs` had a length of %d."
                              % (len(target_names), len(target_seqs)))

        self.query_names = query_names
        self.query_seqs = query_seqs
        self.target_names = target_names
        self.target_seqs = target_seqs
        self.num_threads = num_threads
        self.progress = progress


    def align(self, max_mismatch_freq=0, only_prefixes=False):
        """Perform end-to-end alignment of queries to targets.

        Parameters
        ==========
        max_mismatch_freq : float
            Max mismatch frequency of alignment, on interval [0, 1)

        Returns
        =======
        aligned_query_dict : dict
            Keys are query names and values are AlignedQuery objects.

        aligned_target_dict : dict
            Keys are target names and values are AlignedTarget objects.
        """

        if not 0 <= max_mismatch_freq < 1:
            raise ConfigError("The `max_mismatch_freq` argument of `Aligner.align` must lie in the interval, [0, 1). "
                              "You provided a value of %d." % max_mismatch_freq)

        short_seed_size = min(map(len, self.query_seqs))
        if max_mismatch_freq == 0:
            long_seed_size = None
            short_query_names = self.query_names
            short_query_seqs = self.query_seqs
            long_query_names = []
            long_query_seqs = []
        else:
            long_seed_size = int(round(0.5 / max_mismatch_freq, 0))
            short_query_threshold = 2 * long_seed_size
            short_query_names = []
            short_query_seqs = []
            long_query_names = []
            long_query_seqs = []
            for query_name, query_seq in zip(self.query_names, self.query_seqs):
                if len(query_seq) <= short_query_threshold:
                    short_query_names.append(query_name)
                    short_query_seqs.append(query_seq)
                else:
                    long_query_names.append(query_name)
                    long_query_seqs.append(query_seq)

        if short_query_seqs:
            aligned_short_query_dict, aligned_short_target_dict = self.align_without_indels(short_query_names,
                                                                                            short_query_seqs,
                                                                                            self.target_names,
                                                                                            self.target_seqs,
                                                                                            short_seed_size,
                                                                                            max_mismatch_freq=0)
        else:
            aligned_short_query_dict, aligned_short_target_dict = {}, {}

        if long_query_seqs:
            aligned_long_query_dict, aligned_long_target_dict = self.align_without_indels(long_query_names,
                                                                                          long_query_seqs,
                                                                                          self.target_names,
                                                                                          self.target_seqs,
                                                                                          long_seed_size,
                                                                                          max_mismatch_freq=max_mismatch_freq)
        else:
            aligned_long_query_dict, aligned_long_target_dict = {}, {}

        aligned_short_query_dict.update(aligned_long_query_dict)
        aligned_short_target_dict.update(aligned_long_target_dict)
        aligned_query_dict = aligned_short_query_dict
        aligned_target_dict = aligned_short_target_dict

        return aligned_query_dict, aligned_target_dict


    def match_prefixes(self, max_matches_per_query=float('inf')):
        # self.progress.new("Matching queries to target prefixes")

        kmer_size = min(map(len, self.query_seqs))

        # self.progress.update("Extracting prefix k-mers from target sequences")
        kmer_dict = Kmerizer(self.target_names, self.target_seqs).get_prefix_kmer_dict1(kmer_size)

        pool = multiprocessing.Pool(self.num_threads)
        # num_queries_matching_targets = 0
        # total_query_count = len(self.query_names)
        all_query_matches = []
        for (num_processed_queries,
             query_matches) in enumerate(pool.imap(functools.partial(find_prefix_match,
                                                                     kmer_size=kmer_size,
                                                                     kmer_dict=kmer_dict,
                                                                     max_matches_per_query=max_matches_per_query),
                                                   self.query_seqs,
                                                   chunksize=int(len(self.query_names) / self.num_threads) + 1)):
            all_query_matches.append(query_matches)
            # if query_matches:
            #     num_queries_matching_targets += 1
            # if num_processed_queries % 50000 == 0:
            #     self.progress.update("%d queries of %d/%d processed match target prefixes"
            #                          % (num_queries_matching_targets, num_processed_queries, total_query_count))
        # self.progress.end()

        return all_query_matches


    def match_subseqs(self):
        """Match input queries to input targets without mismatches or indels.

        Returns
        =======
        aligned_query_dict : dict
            Keys are query names and values are AlignedQuery objects

        aligned_target_dict : dict
            Keys are target names and values are AlignedTarget objects
        """

        target_seq_dict = dict(zip(self.target_names, self.target_seqs))

        aligned_query_dict = {}
        aligned_target_dict = {}

        queries_grouped_by_length = [list(group) for length, group in
                                     itertools.groupby(sorted(zip(self.query_names, self.query_seqs), key=lambda t: len(t[1])),
                                                       key=lambda t: len(t[1]))]
        for query_group in queries_grouped_by_length:
            kmer_size = len(query_group[1][1])
            kmer_dict = Kmerizer(self.target_names, self.target_seqs, num_threads=self.num_threads).get_kmer_dict(kmer_size)
            for query_name, query_seq in query_group:
                query_hash = sha224(query_seq.encode('utf-8')).hexdigest()

                if query_hash not in kmer_dict:
                    continue
                hits = kmer_dict[query_hash]

                aligned_query = AlignedQuery(query_seq, name=query_name)
                aligned_query_dict[query_name] = aligned_query

                for target_name, target_start_indices in hits:
                    if target_name in aligned_target_dict:
                        aligned_target = aligned_target_dict[target_name]
                    else:
                        target_seq = target_seq_dict[target_name]
                        aligned_target = AlignedTarget(target_seq, name=target_name)
                        aligned_target_dict[target_name] = aligned_target

                    for target_start_index in target_start_indices:
                        alignment = Alignment(0,
                                            target_start_index,
                                            [(7, len(query_seq))],
                                            aligned_query=aligned_query,
                                            aligned_target=aligned_target)
                        aligned_query.add_alignment(alignment)
                        aligned_target.add_alignment(alignment)

        return aligned_query_dict, aligned_target_dict


    def align_without_indels(self,
                             query_names,
                             query_seqs,
                             target_names,
                             target_seqs,
                             seed_size,
                             max_mismatch_freq=0):
        """Match input queries to input targets without indels.

        Parameters
        ==========
        query_names : list

        query_seqs : list

        target_names : list

        target_seqs : list

        seed_size : int

        max_mismatch_freq : float, 0

        Returns
        =======
        aligned_query_dict : dict
            Keys are query names and values are AlignedQuery objects

        aligned_target_dict : dict
            Keys are target names and values are AlignedTarget objects
        """

        kmer_dict = Kmerizer(target_names, target_seqs, num_threads=self.num_threads).get_kmer_dict(seed_size)

        # self.progress.new("Aligning sequences without indels")
        # self.progress.update("...")

        aligned_query_dict = {}
        aligned_target_dict = {}

        target_seq_dict = dict(zip(target_names, target_seqs))

        pool = multiprocessing.Pool(self.num_threads)
        # num_processed_queries = 0
        # total_query_count = len(query_names)
        for (query_name,
             query_seq,
             alignment_info) in pool.imap_unordered(functools.partial(find_alignment,
                                                                      kmer_dict=kmer_dict,
                                                                      seed_size=seed_size,
                                                                      target_seq_dict=target_seq_dict,
                                                                      max_mismatch_freq=max_mismatch_freq),
                                                    zip(query_names, query_seqs),
                                                    chunksize=int(len(query_names) / self.num_threads) + 1):
            if alignment_info:
                aligned_query = AlignedQuery(query_seq, name=query_name)
                aligned_query_dict[query_name] = aligned_query
                for target_name, target_seq, alignment_target_start, cigartuples in alignment_info:
                    if target_name in aligned_target_dict:
                        aligned_target = aligned_target_dict[target_name]
                    else:
                        aligned_target = AlignedTarget(target_seq, name=target_name)
                        aligned_target_dict[target_name] = aligned_target
                    alignment = Alignment(0, alignment_target_start, cigartuples, aligned_query, aligned_target)
                    aligned_query.alignments.append(alignment)
                    aligned_target.alignments.append(alignment)

            # num_processed_queries += 1
            # if num_processed_queries % 1000 == 0:
            #     self.progress.update("%d/%d query sequences have been processed" % (num_processed_queries, total_query_count))
        pool.close()
        pool.join()

        # self.progress.end()

        return aligned_query_dict, aligned_target_dict


def find_prefix_match(query_seq, kmer_size, kmer_dict, max_matches_per_query=float('inf')):
    query_hash = sha224(query_seq[: kmer_size].encode('utf-8')).hexdigest()

    matches = []
    if query_hash not in kmer_dict:
        return matches

    for target_name, target_seq in kmer_dict[query_hash].items():
        if query_seq == target_seq[: len(query_seq)]:
            matches.append(target_name)
            if len(matches) >= max_matches_per_query:
                break

    return matches


def find_alignment(name_seq_pair, kmer_dict, seed_size, target_seq_dict, max_mismatch_freq=0):
    query_name, query_seq = name_seq_pair
    mismatch_limit = int(max_mismatch_freq * len(query_seq))

    encountered_alignment_sites = []
    alignment_info = []
    for seed_query_start, seed_query_end in zip(range(0, len(query_seq) - seed_size + 1),
                                                range(seed_size, len(query_seq) + 1)):
        seed_hash = sha224(query_seq[seed_query_start: seed_query_end].encode('utf-8')).hexdigest()

        if seed_hash not in kmer_dict:
            continue

        for target_name, seed_target_starts in kmer_dict[seed_hash]:
            for seed_target_start in seed_target_starts:
                if seed_target_start < seed_query_start:
                    continue
                target_seq = target_seq_dict[target_name]
                seed_target_end = seed_target_start + seed_size
                if len(target_seq) - seed_target_end < len(query_seq) - seed_query_end:
                    continue
                alignment_target_start = seed_target_start - seed_query_start

                if (target_name, alignment_target_start) in encountered_alignment_sites:
                    continue

                cigartuples = []
                mismatch_count = 0
                for query_base, target_base in zip(query_seq[0: seed_query_start],
                                                    target_seq[alignment_target_start: seed_target_start]):
                    if query_base == target_base:
                        if cigartuples:
                            if cigartuples[-1][0] == 7:
                                cigartuples[-1][1] += 1
                            else:
                                cigartuples.append([7, 1])
                        else:
                            cigartuples.append([7, 1])
                    else:
                        mismatch_count += 1
                        if cigartuples:
                            if cigartuples[-1][0] == 8:
                                cigartuples[-1][1] += 1
                            else:
                                cigartuples.append([8, 1])
                        else:
                            cigartuples.append([8, 1])
                if mismatch_count > mismatch_limit:
                    continue

                if cigartuples:
                    if cigartuples[-1][0] == 7:
                        cigartuples[-1][1] += seed_size
                    else:
                        cigartuples.append([7, seed_size])
                else:
                    cigartuples.append([7, seed_size])

                for query_base, target_base in zip(query_seq[seed_query_end: len(query_seq)],
                                                   target_seq[seed_target_end: seed_target_end + len(query_seq) - seed_query_end]):
                    if query_base == target_base:
                        if cigartuples[-1][0] == 7:
                            cigartuples[-1][1] += 1
                        else:
                            cigartuples.append([7, 1])
                    else:
                        mismatch_count += 1
                        if cigartuples[-1][0] == 8:
                            cigartuples[-1][1] += 1
                        else:
                            cigartuples.append([8, 1])
                if mismatch_count > mismatch_limit:
                    continue

                alignment_info.append((target_name,
                                       target_seq,
                                       alignment_target_start,
                                       [tuple(cigartuple) for cigartuple in cigartuples]))
                encountered_alignment_sites.append((target_name, alignment_target_start))

    return query_name, query_seq, alignment_info