# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to deal with sequence features"""

import argparse

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2021, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


pp = terminal.pretty_print


class Palindromes:
    def __init__(self, args=argparse.Namespace(), run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.min_palindrome_length = A('min_palindrome_length') or 10
        self.max_num_mismatches = A('max_num_mismatches') or 0
        self.verbose = A('verbose') or False
        self.contigs_db_path = A('contigs_db')
        self.fasta_file_path = A('fasta_file')
        self.output_file_path = A('output_file')

        self.translate = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        self.sanity_check()

        self.run.warning(None, header="SEARCH SETTINGS", lc="green")
        self.run.info('Minimum palindrome length', self.min_palindrome_length)
        self.run.info('Number of mismatches allowed', self.max_num_mismatches)
        self.run.info('Be verbose?', 'No' if not self.verbose else 'Yes', nl_after=1)

        self.palindromes = {}


    def sanity_check(self):
        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                              "class, not both :/")

        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)
        else:
            self.verbose = True

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)

        if self.fasta_file_path:
            filesnpaths.is_file_fasta_formatted(self.fasta_file_path)

        try:
            self.min_palindrome_length = int(self.min_palindrome_length)
        except:
            raise ConfigError("Minimum palindrome length must be an integer.")

        try:
            self.max_num_mismatches = int(self.max_num_mismatches)
        except:
            raise ConfigError("Maximum number of mismatches must be an integer.")

        if self.min_palindrome_length % 2 != 0:
            raise ConfigError("OK. The minimum palindrome length parameter must be an `even` integer. You know, "
                              "because that's how palindromes work.")

        if self.max_num_mismatches > self.min_palindrome_length:
            raise ConfigError("The maximum number of mismatches you want to allow in your palindromes is larger "
                              "than the minimum palindrome length you ask for. This makes no sense :(")

        if self.min_palindrome_length < 5:
            raise ConfigError("For everyone's sake, we set the minimum value for the minimum palindrome length to "
                              "4. You have a problem with that? WELL, WELCOME TO THE CLUB, YOU'LL FIT RIGHT IN -- "
                              "WE HAVE A PROBLEM WITH LOGIC TOO.")


    def process(self):
        """Processes all sequences in a given contigs database or a FASTA file.

        What this function does depends on the configuration of the class. The member funciton `find`
        is a more appropriate one to call if there is a single sequence to process.
        """

        if self.contigs_db_path:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
            contig_sequences_dict = contigs_db.db.get_table_as_dict(anvio.tables.contig_sequences_table_name)

            self.progress.new('Searching', progress_total_items=len(contig_sequences_dict))
            for sequence_name in contig_sequences_dict:
                self.progress.update(f"{sequence_name} ({pp(len(contig_sequences_dict[sequence_name]['sequence']))} nts)", increment=True)
                self.find(contig_sequences_dict[sequence_name]['sequence'], sequence_name=sequence_name)
            self.progress.end()

        elif self.fasta_file_path:
            num_sequences = utils.get_num_sequences_in_fasta(self.fasta_file_path)
            fasta = anvio.fastalib.SequenceSource(self.fasta_file_path)
            self.progress.new('Searching', progress_total_items=num_sequences)
            while next(fasta):
                self.progress.update(f"{fasta.id} ({pp(len(fasta.seq))} nts)", increment=True)
                self.find(fasta.seq, sequence_name=fasta.id)
            self.progress.end()

        else:
            raise ConfigError("You called the `process` function of the class `Palindromes` without a FASTA "
                              "file or contigs database to process :(")

        self.report()


    def find(self, sequence, sequence_name="(a sequence does not have a name)", display_palindromes=False):
        """Find palindromes in a single sequence, and populate `self.palindromes`

        The member function `process` may be a better one to call with an `args` object. See `anvi-search-palindromes`
        for example usage.
        """

        if sequence_name in self.palindromes:
            raise ConfigError(f"The sequence '{sequence_name}' is already in `self.palindromes`.")

        sequence = sequence.upper()
        sequence_length = len(sequence)

        START, END = 0, 1

        if sequence_length < self.min_palindrome_length:
            self.progress.reset()
            self.run.warning(f"The sequence '{sequence_name}', which is only {sequence_length} nts long, is too short to find palindromes "
                             f"that are at least {self.min_palindrome_length} nts, and anvi'o will skip it :/ ")

        get_pairs = lambda: (current - distance, current + distance - 1)
        matching_nts = lambda: sequence[l] == ('N' if sequence[r] not in self.translate else self.translate[sequence[r]])

        #############################################################
        #
        # Find stretches of palindrome-like sequences
        #
        #############################################################

        stretches = []
        max_distance = int(self.min_palindrome_length / 2)
        current = max_distance
        while 1:
            # distance between pairs of nucleotides
            distance = 1
            num_mismatches = 0
            # this variable will keep track of mismatches and their location as we expand from a
            # position to form a palindrome
            mismatch_map = []

            if current + max_distance > sequence_length + 1:
                break

            while 1:
                l, r = get_pairs()

                if r == sequence_length:
                    break

                if matching_nts():
                    # we enter this if if as we expand frmo the center out, we run into a pair
                    # of sequences that match.
                    distance += 1
                    mismatch_map.append(0)

                    if current - distance < 0:
                        break

                elif num_mismatches < self.max_num_mismatches:
                    # this `elif` block enables the inclusion of palindromes with mismatches,
                    # but we never should start a palindrome with a mismatch as we expand from
                    # the center, thus we should kill it now if this is the beginning of a new
                    # palindrome (`distance == 1`) and it is a mismatch:
                    if distance == 1:
                        break

                    distance += 1
                    num_mismatches += 1
                    mismatch_map.append(1)

                    if current - distance < 0:
                        break
                else:
                    break

            if distance == 1:
                # nothing exciting here, moving on to the next nucleotide
                current += 1
                continue

            # OK! if we are here, it means we do have a palindrome-like sequence! but we are not
            # done yet: if the function was called with a `self.max_num_mismatches` parameter,
            # larger than 0, which means our candidate palindrome here may be surrounded with some
            # useless mismatches like this one here:
            #
            #     >>> TACTTAAGTC
            #     >>>  ||||||||
            #     >>> GACTTAAGTA
            #
            # we don't want that. so here we will find those that are surrounded by mismatches, and
            # trim them.
            mismatch_map.reverse()
            surrounding_mismatches = 0
            for surrounding_mismatches in range(0, len(mismatch_map)):
                if mismatch_map[surrounding_mismatches] == 0:
                    break

            # if `surrounding_mismatches` is not 0 here, it means there is
            # at least one mismatch surrounding our palindrome-like sequence.
            # we will decrease the distance by the number of mismatches, and
            # finalize things:
            if surrounding_mismatches:
                distance = distance - surrounding_mismatches

            # if after all that crazy, our distance from center still matches to minimum
            # palindrome length, this is a valid one:
            if distance - 1 >= (self.min_palindrome_length / 2):
                stretches.append((current - distance + 1, current + distance - 1),)

            current += distance - 1

            if (current + max_distance) > sequence_length + 1:
                break

        ###############################################################################
        #
        # We now have the stretches of palindrome-like sequences, but they will include
        # overlapping fragments of varying lengths that belong to the same larger
        # palindrome. They need to go away. So, as our next act, we remove those that
        # are simply fake news.
        #
        ###############################################################################

        while 1:
            indices_to_remove = set([])

            for i in range(0, len(stretches) - 1):
                curr_p = stretches[i]
                next_p = stretches[i+1]
                curr_p_len = curr_p[END] - curr_p[START]
                next_p_len = next_p[END] - next_p[START]

                if next_p[START] < curr_p[END]:
                    # these two stretches of sequences overlap. one of them has to go.
                    # let's see which one is bigger
                    if next_p_len >= curr_p_len:
                        indices_to_remove.add(i)
                    else:
                        indices_to_remove.add(i + 1)

            if not len(indices_to_remove):
                break
            else:
                stretches = [stretches[i] for i in range(0, len(stretches)) if i not in indices_to_remove]

        #################################################################################
        #
        # Of course, THERE HAS TO BE some dumb projects somewhre with some dumb sequences
        # in them that are all Ns. Let's remove any sequence in which more than half of
        # the palindrome nucleotides are Ns.
        #
        #################################################################################

        indices_to_remove = set([])
        for i in range(0, len(stretches)):
            start, end = stretches[i]
            if sequence[start:end].count('N') > ((end - start) / 2):
                indices_to_remove.add(i)

        if len(indices_to_remove):
            stretches = [stretches[i] for i in range(0, len(stretches)) if i not in indices_to_remove]

        ###############################################################################
        #
        # Time to summarize it all
        #
        ###############################################################################

        # if there are no stretches to repot, then we can return
        if not len(stretches):
            return

        # there are stretches to report. good.
        self.palindromes[sequence_name] = []

        for start, end in stretches:
            fwd = sequence[start:end]
            rev = sequence[start:end].translate(constants.complements)[::-1]
            length = end - start

            matches = ''.join(['|' if (fwd[i] == rev[i]) else 'x' for i in range(0, len(fwd))])
            num_mismatches = int(matches.count('x') / 2)

            self.palindromes[sequence_name].append({'start': start,
                                                    'end': end,
                                                    'length': end - start,
                                                    'palindrome': fwd,
                                                    'matches': matches,
                                                    'num_mismatches': num_mismatches})

            if anvio.DEBUG or display_palindromes or self.verbose:
                self.progress.reset()
                self.run.warning(None, header=f'{length} nts palindrome at "{start}:{end}"', lc='yellow')
                self.run.info('Sequence', sequence_name)
                self.run.info('Num mismatches', f"{num_mismatches}")
                self.run.info('FDW', fwd, mc='green')
                self.run.info('ALN', matches, mc='green')
                self.run.info('REV', rev, mc='green')


    def report(self):
        num_sequences = 0
        num_palindromes = 0
        longest_palindrome = 0

        for sequence_name in self.palindromes:
            num_sequences += 1
            for p in self.palindromes[sequence_name]:
                if p['length'] > longest_palindrome:
                    longest_palindrome = p['length']
                num_palindromes += 1

        self.run.warning(None, header="SEARCH RESULTS", lc="green")
        self.run.info('Total number of sequences processed', num_sequences)
        self.run.info('Total number of palindromes found', num_palindromes)
        self.run.info('Longest palindrome', longest_palindrome)

        if self.output_file_path:
            with open(self.output_file_path, 'w') as output_file:
                output_file.write("sequence_name\tstart\tend\tlength\tpalindrome\tmatches\tnum_mismatches\n")
                for sequence_name in self.palindromes:
                    for p in self.palindromes[sequence_name]:
                        output_file.write(f"{sequence_name}\t{p['start']}\t{p['end']}\t{p['length']}\t{p['palindrome']}\t{p['matches']}\t{p['num_mismatches']}\n")

            self.run.info('Output file', self.output_file_path, mc='green', nl_before=1, nl_after=1)


