# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to deal with sequence features"""

import os
import argparse
import xml.etree.ElementTree as ET

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2021, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


P = terminal.pluralize
pp = terminal.pretty_print
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)


class Palindrome:
    def __init__(self, run=terminal.Run()):
        self.run=run

        self.sequence_name = None
        self.first_start = None
        self.first_end = None
        self.first_sequence = None
        self.second_start = None
        self.second_end = None
        self.second_sequence = None
        self.num_mismatches = None
        self.length = None
        self.distance = None
        self.midline = ''


    def __str__(self):
        return f"Len: {self.length}; Dist: {self.distance}; {self.first_sequence} ({self.first_start}:{self.first_end}) :: {self.second_sequence} ({self.second_start}:{self.second_end})"


    def display(self):

        # we don't care what `verbose` variable the original instance may have. if
        # the user requests to `display` things, we will display it, and then store
        # the original state again.
        verbose = self.run.verbose
        self.run.verbose = True

        self.run.warning(None, header=f'{self.length} nts palindrome', lc='yellow')
        self.run.info('1st sequence [start:stop]', f"[{self.first_start}:{self.first_end}]", mc='green')
        self.run.info('2nd sequence [start:stop]', f"[{self.second_start}:{self.second_end}]", mc='green')
        self.run.info('Number of mismatches', f"{self.num_mismatches}", mc='red')
        self.run.info('Distance between', f"{self.distance}", mc='yellow')
        self.run.info('1st sequence', self.first_sequence, mc='green')
        self.run.info('ALN', self.midline, mc='green')
        self.run.info('2nd sequence', self.second_sequence, mc='green')

        # store the original verbose state
        self.run.verbose = verbose


class Palindromes:
    def __init__(self, args=argparse.Namespace(), run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.min_palindrome_length = A('min_palindrome_length') or 10
        self.max_num_mismatches = A('max_num_mismatches') or 0
        self.min_distance = A('min_distance') or 0
        self.verbose = A('verbose') or False
        self.contigs_db_path = A('contigs_db')
        self.fasta_file_path = A('fasta_file')
        self.output_file_path = A('output_file')

        self.num_threads = int(A('num_threads')) if A('num_threads') else 1
        self.blast_word_size = A('blast_word_size') or 10

        self.translate = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        self.sanity_check()

        self.run.warning(None, header="SEARCH SETTINGS", lc="green")
        self.run.info('Minimum palindrome length', self.min_palindrome_length)
        self.run.info('Number of mismatches allowed', self.max_num_mismatches)
        self.run.info('Minimum gap length', self.min_distance)
        self.run.info('Be verbose?', 'No' if not self.verbose else 'Yes', nl_after=1)
        self.run.info('Number of threads for BLAST', self.num_threads)
        self.run.info('BLAST word size', self.blast_word_size, nl_after=1)

        self.user_is_warned_for_potential_performance_issues = False

        self.palindromes = {}


    def sanity_check(self):
        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                              "class, not both :/")

        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=False)
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

        if self.blast_word_size < 4:
            raise ConfigError("For everyone's sake, we set the minimum value for the minimum word size for BLAST to "
                              "5. If you need this to change, please let us know (or run the same command with `--debug` "
                              "flag, find the location of this control, and hack anvi'o by replacing that 4 with something "
                              "smaller -- anvi'o doesn't mind being hacked).")

        if self.min_palindrome_length < 5:
            raise ConfigError("For everyone's sake, we set the minimum value for the minimum palindrome length to "
                              "4. You have a problem with that? WELL, WELCOME TO THE CLUB, YOU'LL FIT RIGHT IN -- "
                              "WE HAVE A PROBLEM WITH LOGIC TOO.")


    def process(self):
        """Processes all sequences in a given contigs database or a FASTA file.

        What this function does depends on the configuration of the class. Member functions `find_gapless`
        or `find_with_gaps` may be more appropriate to call if there is a single sequence to process.
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
        else:
            self.palindromes[sequence_name] = []

        sequence = sequence.upper()
        sequence_length = len(sequence)

        if sequence_length < self.min_palindrome_length * 2 + self.min_distance:
            self.progress.reset()
            self.run.warning(f"The sequence '{sequence_name}', which is only {sequence_length} nts long, is too short "
                             f"to find palindromes that are at least {self.min_palindrome_length} nts, with "
                             f"{self.min_distance} nucleoties in between :/ Anvi'o will skip it.")

        # setup BLAST job
        BLAST_search_tmp_dir = filesnpaths.get_temp_directory_path()
        fasta_file_path = os.path.join(BLAST_search_tmp_dir, 'sequence.fa')
        log_file_path = os.path.join(BLAST_search_tmp_dir, 'blast-log.txt')
        results_file_path = os.path.join(BLAST_search_tmp_dir, 'hits.xml')
        with open(fasta_file_path, 'w') as fasta_file:
            fasta_file.write(f'>sequence\n{sequence}\n')

        # run blast
        blast = BLAST(fasta_file_path, search_program='blastn', run=run_quiet, progress=progress_quiet)
        blast.evalue = 10
        blast.num_threads = self.num_threads
        blast.min_pct_id = 100 - self.max_num_mismatches
        blast.search_output_path = results_file_path
        blast.log_file_path = log_file_path
        blast.makedb(dbtype='nucl')

        if self.min_palindrome_length < 20 and len(sequence) > 10000 and not self.user_is_warned_for_potential_performance_issues:
            self.progress.reset()
            self.run.warning(f"Please note, you are searching for palindromes that are as short as {self.min_palindrome_length} "
                             f"in a sequence that is {pp(len(sequence))} nts long. If your palindrome search takes a VERY long time "
                             f"you may want to go for longer palindromes by setting a different `--min-palindrome-length` parameter "
                             f"and by increasing the BLAST word size using `--blast-word-size` parameter (please read the help menu first). "
                             f"This part of the code does not know if you have many more seqeunces to search, but anvi'o will not "
                             f"continue displaying this warning for additional seqeunces to minimize redundant informatio in your "
                             f"log files (because despite the popular belief anvi'o can actually sometimes be like nice and all).",
                             header="ONE-TIME PERFORMANCE WARNING")
            self.user_is_warned_for_potential_performance_issues = True

        blast.blast(outputfmt='5', word_size=self.blast_word_size, strand='minus')

        # parse the BLAST XML output
        root = ET.parse(blast.search_output_path).getroot()
        for query_sequence_xml in root.findall('BlastOutput_iterations/Iteration'):
            for hit_xml in query_sequence_xml.findall('Iteration_hits/Hit'):

                for hsp_xml in hit_xml.findall('Hit_hsps/Hsp'):
                    p = Palindrome(run=self.run)

                    p.sequence_name = sequence_name
                    p.first_start = int(hsp_xml.find('Hsp_query-from').text) - 1
                    p.first_end = int(hsp_xml.find('Hsp_query-to').text)
                    p.first_sequence = hsp_xml.find('Hsp_qseq').text
                    p.second_start = int(hsp_xml.find('Hsp_hit-to').text) - 1
                    p.second_end = int(hsp_xml.find('Hsp_hit-from').text)
                    p.second_sequence = hsp_xml.find('Hsp_hseq').text
                    p.distance = p.second_start - p.first_start

                    # for each hit, there will be a copy of its reverse complement.
                    # the first half of the if statement below is to control for that
                    # and make sure we keep only one of them. the other half is to
                    # remove those that do not meet the minimum distance criterion.
                    if p.distance < 0 or p.distance < self.min_distance:
                        continue

                    # before we continue, we will test for a special case: internal palindromes
                    # within larger palindromes of 0 distance. IT DOES HAPPEN I PROM.
                    if p.distance == 0:
                        internal_palindrome = False
                        for _p in self.palindromes[sequence_name]:
                            if p.first_start > _p.first_start and p.first_start < _p.first_end:
                                internal_palindrome = True
                                break

                        if internal_palindrome:
                            continue

                    p.length = int(hsp_xml.find('Hsp_align-len').text)

                    if p.length < self.min_palindrome_length:
                        # buckle your seat belt Dorothy, 'cause Kansas is going bye-bye:
                        continue

                    p.num_gaps = int(hsp_xml.find('Hsp_gaps').text)
                    p.num_mismatches = int(hsp_xml.find('Hsp_align-len').text) - int(hsp_xml.find('Hsp_identity').text)
                    p.midline = ''.join(['|' if p.first_sequence[i] == p.second_sequence[i] else 'x' for i in range(0, len(p.first_sequence))])

                    if p.num_mismatches > self.max_num_mismatches or p.num_gaps > 0:
                        # this is the crazy part: read the function docstring for `get_split_palindromes`.
                        # briefly, we conclude that there are too many mismatches in this match, we will
                        # try and see if there is anything we can salvage from it.
                        p_list = self.get_split_palindromes(p, display_palindromes=display_palindromes)
                    else:
                        # there aren't too many mismatches, and the length checks out. we will continue
                        # processing this hit as a sole palindrome
                        p_list = [p]

                    for sp in p_list:
                        if anvio.DEBUG or display_palindromes or self.verbose:
                            self.progress.reset()
                            sp.display()

                        self.palindromes[sequence_name].append(sp)

        # clean after yourself
        if anvio.DEBUG:
            self.run.info("BLAST temporary dir kept", BLAST_search_tmp_dir, nl_before=1, mc='red')
        else:
            filesnpaths.shutil.rmtree(BLAST_search_tmp_dir)


    def resolve_mismatch_map(self, s, min_palindrome_length=15, max_num_mismatches=3):
        """Takes a mismatch map, returns longest palindrome start/ends that satisfy all constraints

        If you would like to test this function, a palindrome mismatch map looks
        like this:

            >>> s = 'ooxooooooooooooooooxoooxoxxoxxoxoooooooooxoxxoxooooooooo--ooo-o-xoxoooxoooooooooooooooooo'

        The rest of the code will document itself using the mismatch map `s` as an example :)
        """

        if len(s) < min_palindrome_length:
            return []

        # The purpose here is ot find the longest stretches of 'o'
        # characters in a given string of mismatch map `s` (which
        # correspond to matching nts in a palindromic sequence)
        # without violating the maximum number of mismatches allowed
        # by the user (which correspond to `c`) characters. The gaps
        # `-` are not allowed by any means, so they are ignored.
        #
        # A single pass of this algorithm over a `s` to
        # identify sections of it with allowed number of mismathces
        # will not give an opportunity to identify longest possible
        # stretches of palindromic sequences. however, a moving
        # start position WILL consider all combinations of substrings,
        # which can later be considered to find best ones. it is
        # difficult to imagine without an example, so for the string `s`
        # shown above, `starts` will look like this:
        #
        # >>> [0, 3, 20, 24, 26, 27, 29, 30, 32, 42, 44, 45, 47, 65, 67, 71]
        starts = [0] + [pos + 1 for pos, c in enumerate(s) if c == 'x']

        # running the state machine below for a given `start` will collect
        # every matching stretch that contains an acceptable number of
        # mismatches. For instance, for a max mismatch requirement of two,
        # the iterations over `s` will identify start-end positions that
        # will yeild the following substrings:
        #
        # >>>  s: ooxooooooooooooooooxoooxoxxoxxoxoooooooooxoxxoxooooooooo--ooo-o-xoxoooxoooooooooooooooooo
        #
        # >>>  0: ooxooooooooooooooooxooo oxxo xoxooooooooo oxxo ooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>>  3: .. ooooooooooooooooxoooxo xox oxoooooooooxo xoxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 20: ................... oooxox oxxo oooooooooxox oxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 24: ....................... oxxo xoxooooooooo oxxo ooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 26: ......................... xox oxoooooooooxo xoxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 27: .......................... oxxo oooooooooxox oxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 29: ............................ xoxooooooooo oxxo ooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 30: ............................. oxoooooooooxo xoxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 32: ............................... oooooooooxox oxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 42: ......................................... oxxo ooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 44: ........................................... xoxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 45: ............................................ oxooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 47: .............................................. ooooooooo  ooo o xoxooo oooooooooooooooooo
        # >>> 65: ................................................................ oxoooxoooooooooooooooooo
        # >>> 67: .................................................................. oooxoooooooooooooooooo
        # >>> 71: ...................................................................... oooooooooooooooooo
        #
        # thus, the list `W` will contain start-end positions for every
        # possible stretch that do not include gap characters.
        W = []
        for start in starts:
            end = start
            num_mismatches = 0

            while 1:
                if s[start] == 'x':
                    start += 1
                    continue

                if s[end] == 'o':
                    end += 1
                elif s[end] == 'x':
                    num_mismatches += 1

                    if num_mismatches > max_num_mismatches:
                        W.append((start, end), )
                        num_mismatches = 0
                        start = end + 1
                        end = start
                    else:
                        end += 1
                else:
                    W.append((start, end), )
                    num_mismatches = 0
                    start = end + 1
                    end = start

                if end + 1 > len(s):
                    if end > start:
                        W.append((start, end), )

                    break

        # remove all the short ones:
        W = [(start, end) for (start, end) in W if end - start > min_palindrome_length]

        # sort all based on substring length
        W = sorted(W, key=lambda x: x[1] - x[0], reverse=True)

        # the following code will pop the longest substring from W[0], then
        # remove all the ones that overlap with it, and take the longest one
        # among those that remain, until all items in `W` are considered.
        F = []
        while 1:
            if not len(W):
                break

            _start, _end = W.pop(0)
            F.append((_start, _end), )

            W = [(start, end) for (start, end) in W if (start < _start and end < _start) or (start > _end and end > _end)]

        # now `F` contains the longest substrings with maximum number of
        # mismatches allowed, which will look like this for the `s` with
        # various minimum length (ML) and max mismatches (MM) parameters:
        #
        # >>>    input s: ooxooooooooooooooooxoooxoxxoxxoxoooooooooxoxxoxooooooooo--ooo-o-xoxoooxoooooooooooooooooo
        #
        # >>> ML 5; MM 0:    oooooooooooooooo             ooooooooo      ooooooooo               oooooooooooooooooo
        # >>> ML 5; MM 1:    ooooooooooooooooxooo         oooooooooxo  oxooooooooo           oooxoooooooooooooooooo
        # >>> ML 5; MM 2: ooxooooooooooooooooxooo       oxoooooooooxo xoxooooooooo         oxoooxoooooooooooooooooo
        # >>> ML15; MM 3: ooxooooooooooooooooxoooxo                                        oxoooxoooooooooooooooooo
        # >>> (...)

        return(sorted(F))


    def get_split_palindromes(self, p, display_palindromes=False):
        """Takes a palindrome object, and splits it into multiple.

        The goal here is to make use of BLAST matches that may include too many mismatches or
        gaps, and find long palindromic regions that still fit our criteria of maximum number
        of mismatches and minimum length.

        We go through the mismatch map, resolve regions that still could be good candidates
        to be considered as palindromes, and return a list of curated palindrome objects.
        """

        split_palindromes = []
        mismatch_map = []

        for i in range(0, len(p.first_sequence)):
            if p.first_sequence[i] == p.second_sequence[i]:
                mismatch_map.append('o')
            elif '-' in [p.first_sequence[i], p.second_sequence[i]]:
                mismatch_map.append('-')
            else:
                mismatch_map.append('x')

        mismatch_map = ''.join(mismatch_map)

        # get all the best substrings by resolving the mismatch map
        substrings = self.resolve_mismatch_map(mismatch_map,
                                               min_palindrome_length=self.min_palindrome_length,
                                               max_num_mismatches=self.max_num_mismatches)
        # if we don't get any substrings, it means it is time to go back
        if not len(substrings):
            return []

        if anvio.DEBUG or display_palindromes or self.verbose:
            self.progress.reset()
            self.run.warning(None, header='SPLITTING A HIT', lc='red')
            self.run.info('1st sequence', p.first_sequence, mc='green')
            self.run.info('ALN', p.midline, mc='green')
            self.run.info('2nd sequence', p.second_sequence, mc='green')

        # using these substrings we will generate a list of `Palindrome` objects
        # to replace the mother object.
        for start, end in substrings:
            split_p = Palindrome()
            split_p.sequence_name = p.sequence_name
            split_p.first_start = p.first_start + start
            split_p.first_end = p.first_start + end
            split_p.first_sequence = p.first_sequence[start:end]
            split_p.second_start = p.second_end - end
            split_p.second_end = p.second_end - start
            split_p.second_sequence = p.second_sequence[start:end]
            split_p.midline = p.midline[start:end]
            split_p.num_gaps = split_p.midline.count('-')
            split_p.num_mismatches = split_p.midline.count('x')
            split_p.length = end - start
            split_p.distance = p.distance

            split_palindromes.append(split_p)

            if anvio.DEBUG or display_palindromes or self.verbose:
                self.progress.reset()
                self.run.info_single(f"    Split [{start}:{end}]", nl_before=1, level=2, mc="red")
                self.run.info('    1st sequence', split_p.first_sequence, mc='green')
                self.run.info('    ALN', split_p.midline, mc='green')
                self.run.info('    2nd sequence', split_p.second_sequence, mc='green')

        return split_palindromes


    def report(self):
        num_sequences = 0
        num_palindromes = 0
        longest_palindrome = 0
        most_distant_palindrome = 0

        for sequence_name in self.palindromes:
            num_sequences += 1
            for palindrome in self.palindromes[sequence_name]:
                if palindrome.length > longest_palindrome:
                    longest_palindrome = palindrome.length
                if palindrome.distance > most_distant_palindrome:
                    most_distant_palindrome = palindrome.distance
                num_palindromes += 1

        if num_palindromes == 0:
            self.run.warning(f"Anvi'o searched {P('sequence', num_sequences)} you have provided and found no "
                             f"palindromes that satisfy your input criteria :/ No output file will be generated.")

            return

        self.run.warning(None, header="SEARCH RESULTS", lc="green")
        self.run.info('Total number of sequences processed', num_sequences)
        self.run.info('Total number of palindromes found', num_palindromes)
        self.run.info('Longest palindrome', longest_palindrome)
        self.run.info('Most distant palindrome', most_distant_palindrome)

        headers = ["sequence_name", "length", "distance", "num_mismatches", "first_start", "first_end", "first_sequence", "second_start", "second_end", "second_sequence", "midline"]
        if self.output_file_path:
            with open(self.output_file_path, 'w') as output_file:
                output_file.write('\t'.join(headers) + '\n')
                for sequence_name in self.palindromes:
                    for palindrome in self.palindromes[sequence_name]:
                        output_file.write('\t'.join([f"{getattr(palindrome, h)}" for h in headers]) + '\n')

            self.run.info('Output file', self.output_file_path, mc='green', nl_before=1, nl_after=1)


class FindPalindrome(object):
    def __init__(self, sequence, min_len=5, max_mis=0):
        self.m = min_len
        self.seq = sequence
        self.rev = utils.rev_comp(self.seq)
        self.N = max_mis

        self.run = anvio.terminal.Run()


    def find(self):
        m, N = self.m, self.N
        L = len(self.seq)

        palindromes = []

        for i in range(L-m+1):
            for j in range(L-i-m):
                n, k = 0, 0
                is_palindrome = False

                while True:
                    if self.rev[j+k] != self.seq[i+k]:
                        # mismatch
                        n += 1

                    if n > N:
                        if is_palindrome:
                            palindromes.append((i, i+k, L-j-k, L-j))
                        break

                    if k == m-1:
                        is_palindrome = True

                    k += 1

        for p in palindromes:
            self.run.info(f"{p}", f"{self.seq[p[0]:p[1]]}, {self.seq[p[2]:p[3]]}")

        return palindromes



