# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes for basic sequence properties and manipulations."""

import functools
import itertools
import multiprocessing
import numpy as np

from itertools import groupby
from hashlib import sha1
from operator import itemgetter

import anvio
import anvio.constants as constants
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


# Multiprocessing currently uses a version of pickle with a maximum byte string size of 4 GB.
# See https://bugs.python.org/issue17560#msg289548
# Large alignment and dereplication results may exceed this limit, so chunk the tasks.
MP_CHUNKSIZE = 10000


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


class Kmerizer:
    def __init__(self, names, seqs, num_threads=1, progress=None):
        """Gets hashed k-mers from input sequences.

        Parameters
        ==========
        names : list
            Name strings corresponding to sequences in `seqs`

        seqs : list
            Sequence strings or numpy arrays (only supported now with `get_prefix_kmer_dict`)

        num_threads : int, 1
            Available threads

        progress : anvio.terminal.Progress object, None
            Updates existing Progress if provided
        """

        if len(names) != len(seqs):
            raise ConfigError("Your Kmerizer input lists were not the same length. "
                              "`names` had a length of %d, while `seqs` had a length of %d."
                              % (len(names), len(seqs)))

        self.names = names
        self.seqs = seqs
        if isinstance(seqs[0], np.ndarray):
            self.as_array = True
        else:
            self.as_array = False
        self.num_threads = num_threads
        self.progress = progress


    def get_prefix_full_seq_dict(self, kmer_size):
        """Relates hashed prefix subsequence k-mers from input sequences to sequences containing k-mers

        Parameters
        ==========
        kmer_size : int
            Length of the prefix subsequence k-mers to be extracted

        Returns
        =======
        kmer_dict : dict
            Nested dict, with the outer dict keyed by hashed prefix subsequence k-mers,
            and each inner dict mapping names of input sequences containing the k-mer to input sequence strings
        """

        if self.progress:
            self.progress.update("Relating prefix subsequences to parent sequences")

        kmer_dict = {}
        for name, seq_string in zip(self.names, self.seqs):
            if kmer_size > len(seq_string):
                continue

            hashed_kmer = sha1(seq_string[: kmer_size].encode('utf-8')).hexdigest()
            if hashed_kmer in kmer_dict:
                kmer_dict[hashed_kmer][name] = seq_string
            else:
                kmer_dict[hashed_kmer] = {name: seq_string}

        return kmer_dict


    def get_prefix_target_dict(self, kmer_size):
        kmer_dict = {}
        targets = []
        for name, seq_string in zip(self.names, self.seqs):
            if kmer_size > len(seq_string):
                continue

            target = MappableAlignedTarget(seq_string)
            targets.append(target)

            hashed_kmer = sha1(seq_string[: kmer_size].encode('utf-8')).hexdigest()
            if hashed_kmer in kmer_dict:
                kmer_dict[hashed_kmer][name] = target
            else:
                kmer_dict[hashed_kmer] = {name: target}

        return kmer_dict, targets


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

            hashed_kmer = sha1(seq[: kmer_size].encode('utf-8')).hexdigest()
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
            Sort k-mer entries in kmer_dict
            in descending order of input sequence length then by input sequence name.

        Returns
        =======
        kmer_dict : dict
            Keys are unique k-mer strings identified from the input sequences.
            Values are tuples of length two tuples.
            The first item in an inner tuple is an input sequence name.
            The second item in an inner tuple is a list of start indices of the k-mer in the input sequence,
            as a given k-mer can occur multiple times in the sequence.
        """

        pool = multiprocessing.Pool(self.num_threads)
        all_seq_kmer_items = pool.map(functools.partial(get_kmer_worker,
                                                        kmer_size=kmer_size,
                                                        include_full_length=include_full_length,
                                                        as_array=self.as_array),
                                      zip(self.names, self.seqs),
                                      chunksize=int(len(self.names) / self.num_threads) + 1)
        pool.close()
        pool.join()

        kmer_dict = {}
        for hashed_kmer, kmer_items in groupby(sorted([kmer_item
                                                       for seq_kmer_items in all_seq_kmer_items
                                                       for kmer_item in seq_kmer_items],
                                                      key=itemgetter(0)),
                                               key=itemgetter(0)):
            if sort_kmer_items:
                kmer_dict[hashed_kmer] = tuple((name, seq_start_indices)
                                               for hashed_kmer, name, seq_start_indices, seq_length
                                               in sorted(kmer_items,
                                                         key=lambda kmer_item: (-kmer_item[3], kmer_item[1])))
            else:
                kmer_dict[hashed_kmer] = tuple((name, seq_start_indices)
                                               for hashed_kmer, name, seq_start_indices, seq_length
                                               in kmer_items)

        return kmer_dict


def get_kmer_worker(name_seq_pair, kmer_size, include_full_length=True, as_array=False):
    """Get information on each k-mer of a certain size in the input sequence.

    This function is located outside `Kmerizer` to allow multithreading.

    Parameters
    ==========
    name_seq_pair : tuple
        A length two tuple of the sequence name string and sequence string

    kmer_size : int
        Length of the k-mers to be extracted

    include_full_length : bool, True
        Include input sequences as k-mers

    as_array : bool, False
        Sequences are provided as a numpy array

    Returns
    =======
    kmer_items : list
        Contains a tuple for each unique k-mer in the input sequence
        Tuples contain four items:
            k-mer hash string,
            input sequence name,
            tuple of k-mer start indices in the input sequence (the k-mer may occur more than once),
            length of the input sequence
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
    if as_array:
        for start_index, stop_index in zip(range(0, seq.size - kmer_size + 1),
                                           range(kmer_size, seq.size + 1)):
            hashed_kmer = sha1(seq[start_index: stop_index].tobytes()).hexdigest()
            hashed_kmers.append(hashed_kmer)
            prelim_kmer_items.append([hashed_kmer, name, [start_index], seq.size])
    else:
        for start_index, stop_index in zip(range(0, len(seq) - kmer_size + 1),
                                           range(kmer_size, len(seq) + 1)):
            hashed_kmer = sha1(seq[start_index: stop_index].encode('utf-8')).hexdigest()
            hashed_kmers.append(hashed_kmer)
            prelim_kmer_items.append([hashed_kmer, name, [start_index], len(seq)])

    kmer_items = []
    if len(set(hashed_kmers)) < len(hashed_kmers):
        # Consolidate identical k-mers drawn from the same sequence.
        prelim_kmer_items.sort(key=itemgetter(0))

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


class Cluster:
    __slots__ = ['representative_name',
                 'representative_seq_string',
                 'representative_extra',
                 'member_names',
                 'member_seqs',
                 'member_extras']

    def __init__(self):
        self.representative_name = None
        self.representative_seq_string = None
        self.representative_extra = None
        self.member_names = []
        self.member_seqs = []
        self.member_extras = []


class Dereplicator:
    def __init__(self, names, seq_strings, extras=None, num_threads=1, progress=None):
        """This class has methods to dereplicate input sequences in different ways.

        Parameters
        ==========
        names : list
            Name strings corresponding to sequences in `seq_strings`

        seq_strings : list
            Sequence strings

        extras : list, None
            Extra items of any type associated with each input sequence

        num_threads : int, 1
            Threads available for multithreaded operations

        progress : anvio.terminal.Progress object, None
            Updates existing Progress if provided
        """
        if len(names) != len(seq_strings):
            raise ConfigError("Your Dereplicator input lists were not the same length. "
                              "`names` had a length of %d, while `seq_strings` had a length of %d."
                              % (len(names), len(seq_strings)))

        if extras:
            if len(extras) != len(seq_strings):
                raise ConfigError("Your Dereplicator input lists were not the same length. "
                                  "`extras` had a length of %d, while `names` and `seq_strings` had a length of %d."
                                  % (len(extras), len(seq_strings)))

        self.names = names
        self.seq_strings = seq_strings
        self.extras = extras
        self.num_threads = num_threads
        self.progress = progress


    def full_length_dereplicate(self):
        """Cluster identical sequences.

        Returns
        =======
        clusters : list
            List of Cluster objects for each dereplicated sequence
        """

        if self.progress:
            self.progress.update("Dereplicating identical sequences")

        clusters = []

        if self.extras:
            seq_info = sorted(zip(self.seq_strings, zip(self.names, self.extras)), key=lambda x: (x[0], x[1][0]))
            for _, group in groupby(seq_info, key=itemgetter(0)):
                # Ex.: list(group) = [(<seq_string 1>, (<name 1>, <extra 1>)), (<seq_string 2>, (<name 2>, <extra 2>)), ...]
                cluster = Cluster()
                for member_info in group:
                    cluster.member_names.append(member_info[1][0])
                    cluster.member_extras.append(member_info[1][1])
                cluster.representative_name = cluster.member_names[0]
                cluster.representative_seq_string = member_info[0]
                clusters.append(cluster)

        else:
            seq_info = sorted(zip(self.seq_strings, self.names), key=itemgetter(0, 1))
            for _, group in groupby(seq_info, key=itemgetter(0)):
                # Ex.: list(group) = [(<seq_string 1>, <name 1>), (<seq_string 2>, <name 2>), ...]
                cluster = Cluster()
                for member_info in group:
                    cluster.member_names.append(member_info[1])
                cluster.representative_name = cluster.member_names[0]
                cluster.representative_seq_string = member_info[0]
                clusters.append(cluster)

        clusters.sort(key=lambda cluster: (-len(cluster.member_names), cluster.representative_name))

        return clusters


    def prefix_dereplicate(self):
        """Form clusters with member sequences being prefix subsequences of seed sequences.

        Member sequences can occur in multiple clusters.

        Returns
        =======
        clusters : list
            List of Cluster objects
        """

        if self.progress:
            self.progress.update("Prefix-dereplicating sequences")

        kmer_size = min(map(len, self.seq_strings))

        kmer_dict, targets = Kmerizer(self.names, self.seq_strings).get_prefix_target_dict(kmer_size)

        hashed_prefixes = [sha1(seq_string[: kmer_size].encode('utf-8')).hexdigest() for seq_string in self.seq_strings]

        for query_name, query_seq_string, query_extra_item, prefix_hash, query_as_target in zip(self.names, self.seq_strings, self.extras, hashed_prefixes, targets):
            if prefix_hash not in kmer_dict:
                continue

            # Record which target sequences contain the query sequence as a prefix subsequence.
            hit_found = False
            for target_name, candidate_target in kmer_dict[prefix_hash].items():
                if query_seq_string == candidate_target.seq_string[: len(query_seq_string)]:
                    if len(query_seq_string) != len(candidate_target.seq_string):
                        candidate_target.alignments.append((query_name, len(query_seq_string), query_extra_item))
                        hit_found = True
            if hit_found:
                query_as_target.hit_another_target = True

        clusters = []
        for query_name, query_extra_item, query_as_target in zip(self.names, self.extras, targets):
            if query_as_target.hit_another_target:
                continue

            cluster = Cluster()
            cluster.representative_name = query_name
            cluster.representative_seq_string = query_as_target.seq_string
            cluster.member_names.append(query_name)
            cluster.member_extras.append(query_extra_item)
            if query_as_target.alignments:
                for item_for_different_query in sorted(query_as_target.alignments,
                                                       key=lambda query_item: (-query_item[1], query_item[0])):
                    cluster.member_names.append(item_for_different_query[0])
                    cluster.member_extras.append(item_for_different_query[2])
            clusters.append(cluster)

        clusters.sort(key=lambda cluster: (-len(cluster.member_names), cluster.representative_name))

        return clusters


class AlignedQuery:
    __slots__ = ['seq_string', 'name', 'alignments']

    def __init__(self, seq_string, name=None):
        self.seq_string = seq_string
        self.name = name
        self.alignments = []


class AlignedTarget:
    __slots__ = ['seq_string', 'name', 'alignments']

    def __init__(self, seq_string, name=None):
        self.seq_string = seq_string
        self.name = name
        self.alignments = []


class MappableAlignedTarget:
    __slots__ = ['seq_string', 'name', 'alignments', 'hit_another_target']

    def __init__(self, seq_string, name=None):
        self.seq_string = seq_string
        self.name = name
        self.alignments = []
        self.hit_another_target = False


class Alignment:
    __slots__ = ['query_start', 'target_start', 'cigartuples', 'alignment_length', 'aligned_query', 'aligned_target']

    def __init__(self, query_start, target_start, cigartuples, aligned_query=None, aligned_target=None):
        self.query_start = query_start
        self.target_start = target_start
        self.cigartuples = cigartuples
        # For now, assume there are no indels in the alignment.
        self.alignment_length = sum(cigartuple[1] for cigartuple in cigartuples)
        self.aligned_query = aligned_query
        self.aligned_target = aligned_target


class Aligner:
    def __init__(self, query_names, query_seq_strings, target_names, target_seq_strings, num_threads=1, progress=None):
        """Align query sequences to target sequences

        Parameters
        ==========
        query_names : list
            name strings corresponding to sequences in `query_seq_strings`

        query_seq_strings : list
            Query sequence strings

        target_names : list
            name strings corresponding to sequences in `target_seq_strings`

        target_seq_strings : list
            Target sequence strings

        num_threads : int
            Available threads

        progress : anvio.terminal.Progress object, None
            Updates existing Progress if provided
        """
        if len(query_names) != len(query_seq_strings):
            raise ConfigError("Your `Aligner` query input lists were not the same length. "
                              "`query_names` had a length of %d, while `query_seq_strings` had a length of %d."
                              % (len(query_names), len(query_seq_strings)))

        if len(target_names) != len(target_seq_strings):
            raise ConfigError("Your `Aligner` target input lists were not the same length. "
                              "`target_names` had a length of %d, while `target_seq_strings` had a length of %d."
                              % (len(target_names), len(target_seq_strings)))

        self.query_names = query_names
        self.query_seq_strings = query_seq_strings
        self.target_names = target_names
        self.target_seq_strings = target_seq_strings
        self.num_threads = num_threads
        self.progress = progress


    def align(self, max_mismatch_freq=0):
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

        short_seed_size = min(map(len, self.query_seq_strings))
        if max_mismatch_freq == 0:
            long_seed_size = None
            short_query_names = self.query_names
            short_query_seq_arrays = [np.frombuffer(s.encode('ascii'), np.uint8) for s in self.query_seq_strings]
            long_query_names = []
            long_query_seq_arrays = []
        else:
            long_seed_size = int(round(0.5 / max_mismatch_freq, 0))
            short_query_threshold = 2 * long_seed_size
            short_query_names = []
            short_query_seq_arrays = []
            long_query_names = []
            long_query_seq_arrays = []
            for query_name, query_seq_string in zip(self.query_names, self.query_seq_strings):
                if len(query_seq_string) <= short_query_threshold:
                    short_query_names.append(query_name)
                    short_query_seq_arrays.append(np.frombuffer(query_seq_string.encode('ascii'), np.uint8))
                else:
                    long_query_names.append(query_name)
                    long_query_seq_arrays.append(np.frombuffer(query_seq_string.encode('ascii'), np.uint8))

        target_seq_arrays = [np.frombuffer(s.encode('ascii'), np.uint8) for s in self.target_seq_strings]

        if short_query_seq_arrays:
            if self.progress:
                if max_mismatch_freq == 0:
                    self.progress.update("Aligning...")
                else:
                    self.progress.update("Aligning shorter queries...")
            aligned_short_query_dict, aligned_short_target_dict = self.align_without_indels(short_query_names,
                                                                                            short_query_seq_arrays,
                                                                                            self.target_names,
                                                                                            target_seq_arrays,
                                                                                            short_seed_size,
                                                                                            max_mismatch_freq=0)
        else:
            aligned_short_query_dict, aligned_short_target_dict = {}, {}

        if long_query_seq_arrays:
            if self.progress:
                self.progress.update("Aligning longer queries...")
            aligned_long_query_dict, aligned_long_target_dict = self.align_without_indels(long_query_names,
                                                                                          long_query_seq_arrays,
                                                                                          self.target_names,
                                                                                          target_seq_arrays,
                                                                                          long_seed_size,
                                                                                          max_mismatch_freq=max_mismatch_freq)
        else:
            aligned_long_query_dict, aligned_long_target_dict = {}, {}

        aligned_short_query_dict.update(aligned_long_query_dict)
        aligned_short_target_dict.update(aligned_long_target_dict)
        aligned_query_dict = aligned_short_query_dict
        aligned_target_dict = aligned_short_target_dict

        return aligned_query_dict, aligned_target_dict


    def align_without_indels(self,
                             query_names,
                             query_seq_arrays,
                             target_names,
                             target_seq_arrays,
                             seed_size,
                             max_mismatch_freq=0):
        """Match input queries to input targets without indels.

        Parameters
        ==========
        query_names : list

        query_seq_arrays : list

        target_names : list

        target_seq_arrays : list

        seed_size : int

        max_mismatch_freq : float, 0

        Returns
        =======
        aligned_query_dict : dict
            Keys are query names and values are AlignedQuery objects.

        aligned_target_dict : dict
            Keys are target names and values are AlignedTarget objects.
        """

        kmer_dict = Kmerizer(target_names, target_seq_arrays, num_threads=self.num_threads).get_kmer_dict(seed_size)

        aligned_query_dict = {}
        aligned_target_dict = {}

        target_seq_dict = dict(zip(target_names, target_seq_arrays))

        pool = multiprocessing.Pool(self.num_threads)
        num_processed_queries = 0
        total_query_count = len(query_names)
        for (query_name,
             query_seq_array,
             alignment_info) in pool.imap_unordered(functools.partial(align_without_indels_worker,
                                                                      kmer_dict=kmer_dict,
                                                                      seed_size=seed_size,
                                                                      target_seq_dict=target_seq_dict,
                                                                      max_mismatch_freq=max_mismatch_freq),
                                                    zip(query_names, query_seq_arrays),
                                                    chunksize=MP_CHUNKSIZE):
            if alignment_info:
                aligned_query = AlignedQuery(''.join(map(chr, query_seq_array)), name=query_name)
                aligned_query_dict[query_name] = aligned_query
                for target_name, target_seq_array, alignment_target_start, cigartuples in alignment_info:
                    if target_name in aligned_target_dict:
                        aligned_target = aligned_target_dict[target_name]
                    else:
                        aligned_target = AlignedTarget(''.join(map(chr, target_seq_array)), name=target_name)
                        aligned_target_dict[target_name] = aligned_target
                    alignment = Alignment(0, alignment_target_start, cigartuples, aligned_query, aligned_target)
                    aligned_query.alignments.append(alignment)
                    aligned_target.alignments.append(alignment)

            num_processed_queries += 1
            if self.progress:
                if num_processed_queries % 10000 == 0:
                    self.progress.update("%d/%d queries aligned without indels"
                                         % (num_processed_queries, total_query_count))
        pool.close()
        pool.join()
        if self.progress:
            self.progress.update("%d/%d queries aligned without indels"
                                 % (total_query_count, total_query_count))

        return aligned_query_dict, aligned_target_dict


    def match_prefixes(self, max_matches_per_query=float('inf')):
        """Match query sequences to target prefix subsequences.

        Parameters
        ==========
        max_matches_per_query : numeric, float('inf')
        Max number of targets to which query matches before satisfied with result

        Returns
        =======
        matched_target_names : list
            List of target names with prefix subsequence match to query; max length equal to `int(max_matches_per_query)`
        """

        if self.progress:
            self.progress.update("Matching prefix sequences")

        kmer_size = min(map(len, self.query_seq_strings))

        # Do not update Progress with prefix k-mer extraction.
        kmer_dict = Kmerizer(self.target_names, self.target_seq_strings).get_prefix_full_seq_dict(kmer_size)

        pool = multiprocessing.Pool(self.num_threads)
        matched_target_names = pool.map(functools.partial(prefix_match_worker,
                                                          kmer_size=kmer_size,
                                                          kmer_dict=kmer_dict,
                                                          max_matches_per_query=max_matches_per_query),
                                        self.query_seq_strings,
                                        chunksize=int(len(self.query_names) / self.num_threads) + 1)
        pool.close()
        pool.join()

        return matched_target_names


def align_without_indels_worker(name_seq_pair, kmer_dict, seed_size, target_seq_dict, max_mismatch_freq=0):
    query_name, query_seq_array = name_seq_pair
    mismatch_limit = int(max_mismatch_freq * query_seq_array.size)

    encountered_alignment_sites = []
    alignment_info = []
    for seed_query_start, seed_query_stop in zip(range(0, query_seq_array.size - seed_size + 1),
                                                 range(seed_size, query_seq_array.size + 1)):
        seed_hash = sha1(query_seq_array[seed_query_start:  ].tobytes()).hexdigest()

        if seed_hash not in kmer_dict:
            continue

        for target_name, seed_target_starts in kmer_dict[seed_hash]:
            for seed_target_start in seed_target_starts:
                if seed_target_start < seed_query_start:
                    continue

                target_seq_array = target_seq_dict[target_name]
                seed_target_stop = seed_target_start + seed_size

                if target_seq_array.size - seed_target_stop < query_seq_array.size - seed_query_stop:
                    continue

                alignment_target_start = seed_target_start - seed_query_start

                if (target_name, alignment_target_start) in encountered_alignment_sites:
                    continue

                cigartuples = []
                mismatch_count = 0
                # Process the part of the alignment preceding the seed.
                for is_match, group in groupby(query_seq_array[0: seed_query_start]
                                               == target_seq_array[alignment_target_start: seed_target_start]):
                    if is_match:
                        cigartuples.append((7, sum(1 for _ in group)))
                    else:
                        num_mismatches = sum(1 for _ in group)
                        cigartuples.append((8, num_mismatches))
                        mismatch_count += num_mismatches
                if mismatch_count > mismatch_limit:
                    continue

                # Process the matching seed in the alignment.
                if cigartuples:
                    if cigartuples[-1][0] == 7:
                        cigartuples[-1] = (7, cigartuples[-1][1] + seed_size)
                    else:
                        cigartuples.append((7, seed_size))
                else:
                    cigartuples.append((7, seed_size))

                # Process the part of the alignment following the seed.
                for is_match, group in groupby(query_seq_array[seed_query_stop: ]
                                               == target_seq_array[seed_target_stop: seed_target_stop + query_seq_array.size - seed_query_stop]):
                    if is_match:
                        cigartuples.append((7, sum(1 for _ in group)))
                    else:
                        num_mismatches = sum(1 for _ in group)
                        cigartuples.append((8, num_mismatches))
                        mismatch_count += num_mismatches
                if mismatch_count > mismatch_limit:
                    continue

                alignment_info.append((target_name,
                                       target_seq_array,
                                       alignment_target_start,
                                       cigartuples))
                encountered_alignment_sites.append((target_name, alignment_target_start))

    return query_name, query_seq_array, alignment_info


def prefix_match_worker(query_seq, kmer_size, kmer_dict, max_matches_per_query=float('inf')):
    """Module-level function intended for multiprocessing the matching of query sequences to target prefix subsequences.

    Parameters
    ==========
    query_seq : str
        Sequence to be matched in full

    kmer_size : int
        Length of prefix k-mer seeds considered from query and target sequences

    kmer_dict : dict
        Output of `Kmerizer.get_prefix_full_seq_dict`

    max_matches_per_query : numeric, float('inf')
        Max number of targets to which query matches before satisfied with result

    Returns
    =======
    matched_target_names : list
        List of target names with prefix subsequence match to query; max length equal to `int(max_matches_per_query)`
    """
    query_hash = sha1(query_seq[: kmer_size].encode('utf-8')).hexdigest()

    matched_target_names = []
    if query_hash not in kmer_dict:
        return matched_target_names

    for target_name, target_seq_string in kmer_dict[query_hash].items():
        if query_seq == target_seq_string[: len(query_seq)]:
            matched_target_names.append(target_name)
            if len(matched_target_names) >= max_matches_per_query:
                break

    return matched_target_names
