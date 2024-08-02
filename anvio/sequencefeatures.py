# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to deal with sequence features"""

import os
import re
import copy
import argparse
import subprocess
import xml.etree.ElementTree as ET

from numba import jit
from numba.typed import List

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
import IlluminaUtils.lib.fastqlib as u


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


P = terminal.pluralize
pp = terminal.pretty_print
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)


class PrimerSearch:
    """A class designed to search primers in (meta)genomic short-reads"""

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.samples_txt = A('samples_txt')
        self.primers_file_path = A('primers_txt')
        self.output_directory_path = A('output_dir')
        self.min_remainder_length = A('min_remainder_length') or 60
        self.stop_after = A('stop_after')
        self.only_report_primer_matches = A('only_report_primer_matches')
        self.only_report_remainders = A('only_report_remainders')

        if self.only_report_primer_matches and self.only_report_remainders:
            raise ConfigError("You can't ask anvi'o to report only primer matches AND only remainders "
                              "at the same time. Please take a look at the help menu.")

        if not (A('samples_dict') or A('samples_txt')):
            raise ConfigError("This class is being initialized incorrectly :/ The `args` object must either include "
                              "path for a `samples-txt` file through `samples_txt` parameter, or a dictionary for "
                              "a samples dictionary through `samples_dict` parameter. See online help for more.")

        if not (A('primers_txt') or A('primers_dict')):
            raise ConfigError("This class is being initialized incorrectly :/ The `args` object must either include "
                              "path for a `primers-txt` file through `primers_txt` parameter, or a dictionary for "
                              "a primers dictionary through `primers_dict` parameter. See online help for more.")

        if self.samples_txt:
            self.samples_dict = utils.get_samples_txt_file_as_dict(self.samples_txt, run=run)
        else:
            self.samples_dict = A('samples_dict')
            if not isinstance(self.samples_dict, dict):
                raise ConfigError("The `primers_dict` parameter must be a literal dictionary.")

        if self.primers_file_path:
            self.primers_dict = utils.get_primers_txt_file_as_dict(self.primers_file_path)
        else:
            self.primers_dict = A('primers_dict')
            if not isinstance(self.primers_dict, dict):
                raise ConfigError("The `primers_dict` parameter must be a literal dictionary.")

        if self.only_report_primer_matches:
            self.min_remainder_length = 0

        if self.output_directory_path:
            filesnpaths.check_output_directory(self.output_directory_path)


        self.reads_are_processed = False

        # dictionary to keep track of some global counts for summary purposes
        self.stats = {'total_read_counter': 0,
                      'total_hits': 0,
                      'sample_counter': 0,
                      'samples': {}}

        # fill in default values
        for sample_name in self.samples_dict:
            self.stats['samples'][sample_name] = {'reads': 0,
                                                  'hits': 0,
                                                  'hits_in_rc': 0,
                                                  'primers': {}}

            for primer in self.primers_dict:
                self.stats['samples'][sample_name]['primers'][primer] = {'raw_hits': 0,
                                                                         'final_hits': 0,
                                                                         'shortest_seq_length_after_match': 0,
                                                                         'longest_remainder_length': 0,
                                                                         'avg_remainder_length': 0}

        self.sanity_check()
        self.print_class_setup()


    def sanity_check(self):
        if not self.samples_dict:
            raise ConfigError("The `self.samples_dict` object is empty :/ This class is initialized incorrectly. You either must ")


    def print_class_setup(self):
        """Display what's up"""

        self.run.info("Samples found", f"({len(self.samples_dict)}) {' '.join(self.samples_dict.keys())}")
        self.run.info("Primers found", f"({len(self.primers_dict)}) {'; '.join(self.primers_dict.keys())}")
        self.run.info('Output directory set', self.output_directory_path, nl_after=1)
        self.run.info('Min remainder length', self.min_remainder_length)
        if self.stop_after:
            self.run.info('Only report primer matches', self.only_report_primer_matches)
            self.run.info('Only report remainders', self.only_report_remainders)
            self.run.info('Stop after', self.stop_after, mc='red', nl_after=1)
        else:
            self.run.info('Only report remainders', self.only_report_remainders)
            self.run.info('Only report primer matches', self.only_report_primer_matches, nl_after=1)


    def process_sample(self, sample_name):
        """Process a single sample.

        Returns
        =======
        sample_dict : dict
            A dictionary that holds `hits` for each primer.
        primers_dict : dict
            A dictionary that holds `matching_seqeunces` for each primer
        """

        if sample_name not in self.samples_dict:
            raise ConfigError(f"Someone is calling `process_sample` with a sample name ('{sample_name}') that "
                              f"does not occur in `self.samples_dict`. What kind of black magic is this?")

        self.stats['sample_counter'] += 1

        self.progress.new("Tick tock")
        self.progress.update('...')

        sample_dict = copy.deepcopy(self.samples_dict[sample_name])
        sample_dict['hits'] = {}

        sample_stats = self.stats['samples'][sample_name]

        primers_dict = copy.deepcopy(self.primers_dict)

        for primer_name in primers_dict:
            primers_dict[primer_name]['matching_sequences'] = {}
            primers_dict[primer_name]['primer_length'] = len(primers_dict[primer_name]['primer_sequence'])

            sample_dict['hits'][primer_name] = 0
            primers_dict[primer_name]['matching_sequences'] = []


        # go through
        removed_due_to_remainder_length = 0
        for pair in ['r1', 'r2']:
            input_fastq_file_path = sample_dict[pair]
            input_fastq = u.FastQSource(input_fastq_file_path, compressed=input_fastq_file_path.endswith('.gz'))

            while input_fastq.next(raw=True) and (sample_stats['hits'] < self.stop_after if self.stop_after else True):
                sample_stats['reads'] += 1

                if sample_stats['reads'] % 10000 == 0:
                    self.progress.update(f"{sample_name} ({self.stats['sample_counter']} of {len(self.samples_dict)}) / {pair} / Reads: {pp(sample_stats['reads'])} / Hits: {pp(sample_stats['hits'])} (in RC: {pp(sample_stats['hits_in_rc'])})")

                found_in_RC = False

                for primer_name in primers_dict:
                    v = primers_dict[primer_name]
                    seq = input_fastq.entry.sequence
                    primer_sequence = v['primer_sequence']
                    match = re.search(primer_sequence, seq)

                    if not match:
                        # no match here. but how about the reverse complement of it?
                        seq = utils.rev_comp(seq)

                        match = re.search(primer_sequence, seq)

                        if match:
                            # aha. the reverse complement sequence that carries our match found.
                            # will continue as if nothing happened
                            sample_stats['hits_in_rc'] += 1
                            found_in_RC = True

                    if not match:
                        continue

                    sample_stats['primers'][primer_name]['raw_hits'] += 1

                    if len(seq) - match.end() < self.min_remainder_length:
                        removed_due_to_remainder_length += 1
                        continue

                    v['matching_sequences'].append((match.start(), match.end(), seq), )

                    sample_stats['hits'] += 1
                    sample_stats['primers'][primer_name]['final_hits'] += 1
                    sample_dict['hits'][primer_name] += 1

                    if anvio.DEBUG:
                        self.progress.end()
                        print("\n%s -- %s -- %s| %s [%s] %s" % (sample_name, pair, 'RC ' if found_in_RC else '   ', seq[:match.start()], primer_sequence, seq[match.end():]))
                        self.progress.new("Tick tock")
                        self.progress.update(f"{sample_name} ({self.stats['sample_counter']} of {len(self.samples_dict)}) / {pair} / Reads: {pp(sample_stats['reads'])} / Hits: {pp(sample_stats['hits'])} (in RC: {pp(sample_stats['hits_in_rc'])})")

            self.stats['total_read_counter'] += sample_stats['reads']
            self.stats['total_hits'] += sample_stats['hits']

            # calculate and store stats for primer hits
            for primer_name in primers_dict:
                matching_sequence_hits = primers_dict[primer_name]['matching_sequences']
                seq_lengths_after_match = [len(sequence[end:]) for start, end, sequence in matching_sequence_hits]

                if len(seq_lengths_after_match):
                    sample_stats['primers'][primer_name]['shortest_seq_length_after_match'] = min(seq_lengths_after_match)
                    sample_stats['primers'][primer_name]['longest_remainder_length'] = max(seq_lengths_after_match)
                    sample_stats['primers'][primer_name]['avg_remainder_length'] = sum(seq_lengths_after_match) / len(seq_lengths_after_match)

        self.progress.end()

        return sample_dict, primers_dict


    def process(self):
        """Processes everything."""


        for sample_name in self.samples_dict:
            sample_dict, primers_dict = self.process_sample(sample_name)

            if self.output_directory_path:
                self.store_sequences(sample_name, sample_dict, primers_dict)

            # call Batman
            del sample_dict
            del primers_dict

        self.reads_are_processed = True


    def print_summary(self):
        """Prints a fancy summary of the results"""

        if not self.reads_are_processed:
            raise ConfigError("You first need to call the member function `process`.")

        self.run.warning(None, header="FINAL SUMMARY", lc='yellow')
        self.run.info_single(f"After processing {pp(self.stats['total_read_counter'])} individual reads in {P('sample', len(self.samples_dict))}, "
                             f"anvi'o found {P('hit', self.stats['total_hits'])} for your {P('sequence', len(self.primers_dict), pfs='only')} "
                             f"in the primer sequences file. What is shown below breaks these numbers down per sample because that's how "
                             f"anvi'o rolls. There are some acronyms below, and they are very creatively named. RH: number of raw hits (the actual "
                             f"number of times a primer matched to a sequence). FH: number of final hits (after testing whether the remainder length "
                             f"was longer than the user-set or default minimum value). SRL: shortest remainder length. LRL: Longest remainder length. "
                             f"ARL: Average remainder length after match).", level=0)

        for sample_name in self.samples_dict:
            self.run.info(f"{sample_name}", f"{pp(self.stats['samples'][sample_name]['reads'])} total reads", nl_before=1, lc="green", mc="green")
            for primer_sequence in self.primers_dict:
                s = self.stats['samples'][sample_name]['primers'][primer_sequence]
                self.run.info(f"    {primer_sequence}", f"RH: {pp(s['raw_hits']) if s['raw_hits'] else '--'} / "
                                                        f"FH: {pp(s['final_hits']) if s['final_hits'] else '--'} / "
                                                        f"SRL: {s['shortest_seq_length_after_match'] if s['shortest_seq_length_after_match'] else '--'} / "
                                                        f"LRL: {s['longest_remainder_length'] if s['longest_remainder_length'] else '--'} / "
                                                        f"ARL: {s['avg_remainder_length']:.2f}", mc=('yellow' if s['final_hits'] else 'red'))

        if not self.stats['total_read_counter']:
            self.run.info_single('No hits were found :/', mc='red', nl_before=1)


    def get_sequences(self, primer_name, primers_dict, target):
        """For a given list of `matching_sequences`, recover sequences of various nature from matches.

        Parameters
        ==========
        primer_name : string
            This is essentially a key that exists in `primers_dict`.
        primers_dict : dict
            The sample-specific primers dictionary (recovered from `self.process_sample`).
        target : string
            It can be one of the following: 'remainders', 'primer_match', 'trimmed', 'gapped'. Remainders are the downstream
            sequences after primer match, excluding the primer sequence. Primer matches are the primer-matching part of the
            match sequences (useful if one is working with degenerate primers and wishes to see the diversity of matching
            seqeunces). Trimmed sequences are trimmed to the shortest length (and include primer match). Gapped sequences
            are not trimmed, but shorter ones are padded with gaps.

        Returns
        =======
        sequences : list
            Where each item is the desired part of the match sequences.
        """

        valid_targets = ['remainders', 'primer_match', 'trimmed', 'gapped']

        if target not in valid_targets:
            raise ConfigError(f"You gon to the wrong neighborhood (like the way Tom Hanks did it in Berlin and lost his jacket). "
                              f"Try one of these targets: \"{', '.join(valid_targets)}\".")

        l = []

        primer_length = len(self.primers_dict[primer_name]['primer_sequence'])
        match_sequences = primers_dict[primer_name]['matching_sequences']

        if target in ['trimmed', 'gapped']:
            seq_lengths_after_match = [len(sequence[end:]) for start, end, sequence in match_sequences]
            max_seq_length = (max(seq_lengths_after_match) + primer_length) if len(seq_lengths_after_match) else 0
            min_seq_length = (min(seq_lengths_after_match) + primer_length) if len(seq_lengths_after_match) else 0
        else:
            seq_lengths_after_match = None
            max_seq_length = None
            min_seq_length = None

        if target == 'remainders':
            for start, end, sequence in match_sequences:
                if self.min_remainder_length:
                    l.append(sequence[end:end + self.min_remainder_length])
                else:
                    l.append(sequence[end:])
        elif target == 'primer_match':
            for start, end, sequence in match_sequences:
                l.append(sequence[start:end])
        elif target == 'trimmed':
            for start, end, sequence in match_sequences:
                l.append(sequence[:min_seq_length])
        elif target == 'gapped':
            for start, end, sequence in match_sequences:
                sequence = sequence[start:]
                l.append(sequence + '-' * (max_seq_length - len(sequence)))
        else:
            pass

        return l


    def store_sequences(self, sample_name, sample_dict, primers_dict):
        """Store sequence files for a given sample under `self.output_directory_path`"""

        if not self.output_directory_path:
            return

        if not os.path.exists(self.output_directory_path):
            filesnpaths.gen_output_directory(self.output_directory_path)

        self.run.info_single(f"Output files for {sample_name}:", nl_before=1, mc="green")

        if self.only_report_remainders:
            self.progress.new("Generating the remainders file")
            self.progress.update('...')
            for primer_name in primers_dict:
                remainder_sequences = self.get_sequences(primer_name, primers_dict, target='remainders')
                output_file_path = os.path.join(self.output_directory_path, '%s-%s-REMAINDERS.fa' % (sample_name, primer_name))
                with open(output_file_path, 'w') as output:
                    counter = 1
                    for sequence in remainder_sequences:
                        output.write(f'>{sample_name}_{primer_name}_{counter:05d}\n{sequence}\n')
                        counter += 1
            self.progress.end()

            self.run.info('    Remainders sequences', os.path.join(self.output_directory_path, f'{sample_name}-*-REMAINDERS.fa'))

            return

        self.progress.new("Generating the primer matches files")
        self.progress.update('...')
        for primer_name in primers_dict:
            primer_matching_sequences = self.get_sequences(primer_name, primers_dict, target='primer_match')

            output_file_path = os.path.join(self.output_directory_path, '%s-%s-PRIMER-MATCHES.fa' % (sample_name, primer_name))
            with open(output_file_path, 'w') as output:
                counter = 1
                for sequence in primer_matching_sequences:
                    output.write(f'>{sample_name}_{primer_name}_{counter:05d}\n{sequence}\n')
                    counter += 1
        self.progress.end()

        self.run.info('    Primer matches', os.path.join(self.output_directory_path, f'{sample_name}-*-PRIMER-MATCHES.fa'))

        if self.only_report_primer_matches:
            return

        self.progress.new("Generating the fancy hits files")
        self.progress.update('...')
        for primer_name in self.primers_dict:
            trimmed_output_file_path = os.path.join(self.output_directory_path, '%s-%s-HITS-TRIMMED.fa' % (sample_name, primer_name))
            sequences = self.get_sequences(primer_name, primers_dict, target='trimmed')
            with open(trimmed_output_file_path, 'w') as trimmed:
                counter = 1
                for sequence in sequences:
                    trimmed.write(f'>{sample_name}_{primer_name}_{counter:05d}\n{sequence}\n')
                    counter += 1

            gapped_output_file_path = os.path.join(self.output_directory_path, '%s-%s-HITS-WITH-GAPS.fa' % (sample_name, primer_name))
            sequences = self.get_sequences(primer_name, primers_dict, target='gapped')
            with open(gapped_output_file_path, 'w') as gapped:
                counter = 1
                for sequence in sequences:
                    gapped.write(f'>{sample_name}_{primer_name}_{counter:05d}\n{sequence}\n')
                    counter += 1

        self.progress.end()

        self.run.info('    Trimmed hits', os.path.join(self.output_directory_path, f'{sample_name}-*-HITS.fa'))
        self.run.info('    Hits with gaps', os.path.join(self.output_directory_path, f'{sample_name}-*-HITS.fa'))

        return



class Palindrome:
    def __init__(self, sequence_name=None, first_start=None, first_end=None, first_sequence=None, second_start=None,
                 second_end=None, second_sequence=None, num_mismatches=None, length=None, distance=None, num_gaps=None,
                 midline='', method=None, run=terminal.Run()):
        self.run = run
        self.sequence_name = sequence_name
        self.first_start = first_start
        self.first_end = first_end
        self.first_sequence = first_sequence
        self.second_start = second_start
        self.second_end = second_end
        self.second_sequence = second_sequence
        self.num_mismatches = num_mismatches
        self.length = length
        self.distance = distance
        self.num_gaps = num_gaps
        self.midline = midline
        self.method = method


    def __str__(self):
        return f"Len: {self.length}; Dist: {self.distance}; {self.first_sequence} ({self.first_start}:{self.first_end}) :: {self.second_sequence} ({self.second_start}:{self.second_end})"


    def display(self):

        # we don't care what `verbose` variable the original instance may have. if
        # the user requests to `display` things, we will display it, and then store
        # the original state again.
        verbose = self.run.verbose
        self.run.verbose = True

        self.run.warning(None, header=f'{self.length} nts palindrome', lc='yellow')
        self.run.info('Method', self.method, mc='red')
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
        self.palindrome_search_algorithm = A('palindrome_search_algorithm')
        self.min_palindrome_length = 10 if A('min_palindrome_length') == None else A('min_palindrome_length')
        self.max_num_mismatches = A('max_num_mismatches') or 0
        self.min_distance = A('min_distance') or 0
        self.min_mismatch_distance_to_first_base = A('min_mismatch_distance_to_first_base') or 1
        self.verbose = A('verbose') or False
        self.contigs_db_path = A('contigs_db')
        self.fasta_file_path = A('fasta_file')
        self.output_file_path = A('output_file')

        self.palindrome_search_algorithms = {
            'BLAST': self._find_BLAST,
            'numba': self._find_numba,
        }

        self.num_threads = int(A('num_threads')) if A('num_threads') else 1
        self.blast_word_size = A('blast_word_size') or 10

        self.translate = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        self.sanity_check()

        self.run.warning(None, header="SEARCH SETTINGS", lc="green")
        self.run.info('Palindrome search algorithm', self.palindrome_search_algorithm or "[will be dynamically determined based on sequence length]", mc="red")
        self.run.info('Minimum palindrome length', self.min_palindrome_length)
        self.run.info('Number of mismatches allowed', self.max_num_mismatches)
        self.run.info('Min mismatch distance to the first base', self.min_mismatch_distance_to_first_base)
        self.run.info('Minimum dist for each seq of a palindrome', self.min_distance)
        self.run.info('Be verbose?', 'No' if not self.verbose else 'Yes', nl_after=1)
        self.run.info('Number of threads for BLAST', self.num_threads)
        self.run.info('Word size for BLAST', self.blast_word_size, nl_after=1)

        self.user_is_warned_for_potential_performance_issues = False

        self.palindromes = {}


    def sanity_check(self):
        if self.contigs_db_path and self.fasta_file_path:
            raise ConfigError("You should either choose a FASTA file or a contigs db to send to this "
                              "class, not both :/")

        if self.min_mismatch_distance_to_first_base < 1:
            raise ConfigError("The minimum mismatch thistance to the first base from either of the palindrome "
                              "must be greater than 0.")

        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=False)

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

        if self.min_palindrome_length < 4:
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


    def find(self, sequence, sequence_name="N/A", display_palindromes=False, **kwargs):
        """Find palindromes in a single sequence, and populate `self.palindromes`

        This method finds palindromes by delegating to either `_find_BLAST` and `_find_numba`, unless the user
        has explicitly defined which algorithm to use for search.

        Notes
        =====
        - The method `process` may be a better one to call if you have an `args` object. See `anvi-search-palindromes`
          for example usage.
        """

        if sequence_name in self.palindromes:
            raise ConfigError(f"The sequence '{sequence_name}' is already in `self.palindromes`.")
        else:
            self.palindromes[sequence_name] = []

        sequence_length = len(sequence)
        if sequence_length < self.min_palindrome_length * 2 + self.min_distance:
            self.progress.reset()
            if sequence_name == 'N/A':
                friendly_sequence_name = "The sequence you have provided"
            else:
                friendly_sequence_name = f"The sequence '{sequence_name}'"
            self.run.warning(f"{friendly_sequence_name} is only {sequence_length} nts long, and so it is too "
                             f"short to find any palindromes in it that are at least {self.min_palindrome_length} nts with "
                             f"{self.min_distance} nucleoties in between :/ Anvi'o will most likely skip it.")

        # determine which palindrome search algorithm to use.
        if self.palindrome_search_algorithm:
            # which means the search algorithm is set by the user:
            method = self.palindrome_search_algorithms[self.palindrome_search_algorithm]
        else:
            # which means we are going to dynamically determine which algorith to use
            # as a function of the sequence length
            method = self.palindrome_search_algorithms['BLAST'] if len(sequence) >= 5000 else self.palindrome_search_algorithms['numba']

        # make sure the sequence is uppercase:
        sequence = sequence.upper()

        # get palindromes found in the sequence
        palindromes = method(sequence, **kwargs)

        for palindrome in palindromes:
            if anvio.DEBUG or display_palindromes or self.verbose:
                self.progress.reset()
                palindrome.display()

            palindrome.sequence_name = sequence_name
            self.palindromes[sequence_name].append(palindrome)


    def _find_BLAST(self, sequence):
        """Find palindromes in a single sequence using BLAST

        This method of palindrome finding is slow for short sequences, but scales very well with
        arbitrarily sized palindrome searching, e.g. entire assemblies or genomes. If the sequence
        length > 5000, this is a great choice.

        Here are the timings as a function of sequence length.

            length 100:     51.7 ms ± 2.34 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 278:     48.4 ms ± 947 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 774:     47.9 ms ± 570 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 2154:    49.8 ms ± 541 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 5994:    52 ms ± 2.53 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 16681:   51.7 ms ± 555 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 46415:   56.4 ms ± 623 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 129154:  63.8 ms ± 510 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 359381:  100 ms ± 848 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 1000000: 324 ms ± 4.95 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

        Returns
        =======
        out : list of anvio.sequencefeatures.Palindrome
            Returns a list of Palindrome objects. If no palindromes are found, empty list is returned.
        """
        palindromes = []
        sequence = sequence.upper()

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

        # setup BLAST job
        BLAST_search_tmp_dir = filesnpaths.get_temp_directory_path()
        fasta_file_path = os.path.join(BLAST_search_tmp_dir, 'sequence.fa')
        with open(fasta_file_path, 'w') as fasta_file:
            fasta_file.write(f'>sequence\n{sequence}\n')

        search_command = f"""blastn -query {fasta_file_path} \
                                    -subject {fasta_file_path} \
                                    -evalue 10 \
                                    -outfmt 5 \
                                    -num_threads {self.num_threads} \
                                    -word_size {self.blast_word_size} \
                                    -strand minus"""

        p = subprocess.Popen(search_command, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
        BLAST_output = p.communicate()[0]

        # parse the BLAST XML output
        root = ET.fromstring(BLAST_output.decode())
        for query_sequence_xml in root.findall('BlastOutput_iterations/Iteration'):
            for hit_xml in query_sequence_xml.findall('Iteration_hits/Hit'):

                for hsp_xml in hit_xml.findall('Hit_hsps/Hsp'):
                    p = Palindrome(method='BLAST', run=self.run)

                    p.first_start = int(hsp_xml.find('Hsp_query-from').text) - 1
                    p.first_end = int(hsp_xml.find('Hsp_query-to').text)
                    p.first_sequence = hsp_xml.find('Hsp_qseq').text

                    p.second_start = int(hsp_xml.find('Hsp_hit-to').text) - 1
                    p.second_end = int(hsp_xml.find('Hsp_hit-from').text)
                    p.second_sequence = hsp_xml.find('Hsp_hseq').text

                    # Calculating the 'distance' next. But it is a bit tricky. Imagine this as your genomic context for
                    # this 'in-place' palindrome:
                    #
                    #    >>> 0        1
                    #    >>> 1234567890
                    #    >>> ...TCGA...
                    #
                    # where you indeed have a proper palindrome here. the start and end of both sequences of this
                    # palindrome will be the same: TCGA (3:7) :: TCGA (3:7). In this case, we can't simply calculate
                    # 'distance' by substracting the start of the second sequence from the end of the first, OR we
                    # can't simply remove it from our consideration because p.second_start - p.first_end is a negative
                    # value.
                    #
                    # In contrast, consider this as your genomic context for this 'distance palindrome':
                    #
                    #    >>> 0        1
                    #    >>> 12345678901234567
                    #    >>> ...ATCC...GGAT...
                    #
                    # This also is a proper palindrome. But the start and the end of each sequence will be different this
                    # time in the BLAST results: ATCC (3:7) :: ATCC (10:14). And for such distant palindromes, BLAST results
                    # will ALWAYS include the same result for its reverse complement, where p.second_start - p.first_end will
                    # be negative, which we will want to remove. So the following few lines consider all these scenarios
                    # to not always remove 'in-place' palindromes.

                    if p.first_start == p.second_start:
                        # this is an in-place palindrome. which means, the distance
                        # between these sequences is 0 and we have to manually set it
                        p.distance = 0
                    else:
                        # this is a distant palindrome, so we calculate the distance
                        # from actual positions:
                        p.distance = p.second_start - p.first_end

                    # for each distant palindrome in the sequence, there will be a copy of the reverse complement of the first
                    # hit. now we have set the distance properly, we can remove those from hits to be considered:
                    if p.distance < 0:
                        continue

                    # time to check the remaining ones for minimum distance, if it is defined:
                    if p.distance < self.min_distance:
                        continue

                    # before we continue, we will test for a special case: internal palindromes
                    # within larger palindromes of 0 distance. IT DOES HAPPEN I PROM.
                    if p.distance == 0:
                        internal_palindrome = False
                        for _p in palindromes:
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
                        p_list = self.get_split_palindromes(p)
                    else:
                        # there aren't too many mismatches, and the length checks out. we will continue
                        # processing this hit as a sole palindrome
                        p_list = [p]

                    for sp in p_list:
                        palindromes.append(sp)

        # clean after yourself
        if anvio.DEBUG:
            self.run.info("BLAST temporary dir kept", BLAST_search_tmp_dir, nl_before=1, mc='red')
        else:
            filesnpaths.shutil.rmtree(BLAST_search_tmp_dir)

        return palindromes


    def _find_numba(self, sequence, coords_only=False):
        """Find palindromes in a single sequence using a numba state machine

        This method of palindrome specializes in finding palindromes for short sequences (<5000). If
        the sequence length < 5000, this is a great choice.

        Here are timings as a function of sequence size:

            length 100:   17.6 µs ± 480 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
            length 166:   42.8 µs ± 431 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)
            length 278:   139 µs ± 1.88 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)
            length 464:   429 µs ± 4.18 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
            length 774:   1.3 ms ± 10.5 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
            length 1291:  4.01 ms ± 128 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
            length 2154:  11.2 ms ± 492 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
            length 3593:  31.3 ms ± 1.05 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 5994:  88.8 ms ± 410 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
            length 10000: 241 ms ± 1.86 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

        Parameters
        ==========
        coords_only : bool, False
            If True

        Returns
        =======
        out : list of anvio.sequencefeatures.Palindrome
            By default, returns a list of Palindrome objects. But if coords_only is True, a list of
            4-length tuples are returned. Each tuple is (first_start, first_second, second_start, second_end).
            If no palindromes are found, an empty list is returned.
        """
        sequence_array = utils.nt_seq_to_nt_num_array(sequence)
        sequence_array_RC = utils.nt_seq_to_RC_nt_num_array(sequence)

        palindrome_coords = _find_palindromes(
            seq = sequence_array,
            rev = sequence_array_RC,
            m = self.min_palindrome_length,
            N = self.max_num_mismatches,
            D = self.min_distance,
            Q = self.min_mismatch_distance_to_first_base,
        )

        # sort by longest palindrome first
        palindrome_coords = sorted(palindrome_coords, key = lambda y: y[1]-y[0], reverse=True)

        if coords_only:
            return palindrome_coords

        palindromes = []
        for first_start, first_end, second_start, second_end in palindrome_coords:
            first_sequence = sequence[first_start:first_end]
            second_sequence = utils.rev_comp(sequence[second_start:second_end])
            num_mismatches = sum(1 for a, b in zip(first_sequence, second_sequence) if a != b)
            midline = ''.join('|' if a == b else 'x' for a, b in zip(first_sequence, second_sequence))

            palindrome = Palindrome(
                first_start = first_start,
                first_end = first_end,
                first_sequence = first_sequence,
                second_start = second_start,
                second_end = second_end,
                second_sequence = second_sequence,
                num_mismatches = num_mismatches,
                midline = midline,
                distance = second_start - first_end,
                length = len(first_sequence),          # length of first will always be length of second
                num_gaps = 0,                          # via the behavior of _find_palindrome
                method = 'numba',
                run = self.run,
            )

            palindromes.append(palindrome)

        return palindromes


    def resolve_mismatch_map(self, s, min_palindrome_length=15, max_num_mismatches=3, min_mismatch_distance_to_first_base=1):
        """Takes a mismatch map, returns longest palindrome start/ends that satisfy all constraints

        If you would like to test this function, a palindrome mismatch map looks
        like this:

            >>> s = 'ooxooooooooooooooooxoooxoxxoxxoxoooooooooxoxxoxooooooooo--ooo-o-xoxoooxoooooooooooooooooo'

        The rest of the code will document itself using the mismatch map `s` as an example. One could quickly test
        how a mismatch map is resolved for a given `s` by running this function this way:

            >>> from anvio.sequencefeatures import Palindromes
            >>> resolution = Palindromes().resolve_mismatch_map(s)

        Here is a more elaborate example:

            >>> from anvio.sequencefeatures import Palindromes as P
            >>> from anvio.terminal import Run
            >>> 
            >>> p = P(run=Run(verbose=False))
            >>> 
            >>> s = 'ooxooooooooooooooooxoooxoxxoxxoxoooooooooxoxxoxooooooooo--ooo-o-xoxoooxoooooooooooooooooo'
            >>> 
            >>> RESOLUTION = lambda: [s[start:end] for start, end in p.resolve_mismatch_map(s, max_num_mismatches=max_num_mismatches, min_mismatch_distance_to_first_base=min_mismatch_distance_to_first_base, min_palindrome_length=min_palindrome_length)]
            >>> 
            >>> max_num_mismatches=1
            >>> min_mismatch_distance_to_first_base=1
            >>> min_palindrome_length=5
            >>> print(RESOLUTION())
            ['ooooooooooooooooxooo', 'oxooooooooo', 'oxooooooooo', 'oooxoooooooooooooooooo']
            >>> 
            >>> max_num_mismatches=1
            >>> min_mismatch_distance_to_first_base=2
            >>> min_palindrome_length=5
            >>> print(RESOLUTION())
            ['ooooooooooooooooxooo', 'ooooooooo', 'ooooooooo', 'oooxoooooooooooooooooo']
            >>> 
            >>> max_num_mismatches=1
            >>> min_mismatch_distance_to_first_base=5
            >>> min_palindrome_length=5
            >>> print(RESOLUTION())
            ['oooooooooooooooo', 'ooooooooo', 'ooooooooo', 'oooooooooooooooooo']
            >>> 
            >>> max_num_mismatches=5
            >>> min_mismatch_distance_to_first_base=1
            >>> min_palindrome_length=5
            >>> print(RESOLUTION())
            ['ooxooooooooooooooooxoooxoxxo', 'oxoooooooooxoxxoxooooooooo', 'oxoooxoooooooooooooooooo']
            >>> 
            >>> max_num_mismatches=5
            >>> min_mismatch_distance_to_first_base=2
            >>> min_palindrome_length=5
            >>> print(RESOLUTION())
            ['ooxooooooooooooooooxooo', 'oooooooooxoxxoxooooooooo', 'oooxoooooooooooooooooo']
        """

        if len(s) < min_palindrome_length:
            return []

        # a helper function.
        ADD = lambda _: _.add((start, end), ) if end - start > min_palindrome_length else None

        # The purpose here is to find the longest stretches of 'o'
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
        # possible stretch that do not include gap characters and long
        # enough to be worth considering.
        W = set([])
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
                        ADD(W)
                        num_mismatches = 0
                        start = end + 1
                        end = start
                    else:
                        end += 1
                else:
                    ADD(W)
                    num_mismatches = 0
                    start = end + 1
                    end = start

                if end + 1 > len(s):
                    if end > start:
                        ADD(W)

                    break

        W = sorted(W)

        # make sure all stretches respect the `min_mismatch_distance_to_first_base` criterion.
        L = set([])
        for start, end in W:
            while 1:
                dirty = False
                for i in range(0, min_mismatch_distance_to_first_base):
                    if s[end-i-1] == 'x':
                        end = end - i - 1
                        dirty = True

                    if s[start + i] == 'x':
                        start = start + i + 1
                        dirty = True

                    # if no mismatches left to work with, call it a day
                    if 'x' not in s[start:end]:
                        dirty = False
                        break

                if not dirty:
                    ADD(L)
                    break
        W = L

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
                                               max_num_mismatches=self.max_num_mismatches,
                                               min_mismatch_distance_to_first_base=self.min_mismatch_distance_to_first_base)
        # if we don't get any substrings, it means it is time to go back
        if not len(substrings):
            return []

        if anvio.DEBUG or display_palindromes or self.verbose:
            self.progress.reset()
            self.run.warning(None, header=f'SPLITTING A HIT [{p.first_start}:{p.first_end} / {p.first_start}:{p.second_end}]', lc='red')
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

            if split_p.first_start > split_p.second_start:
                # this can only be the case if we are splitting an in-place palindrome
                # so we will not report the mirrored pair.
                continue

            split_p.method = p.method
            split_p.midline = p.midline[start:end]
            split_p.num_gaps = split_p.midline.count('-')
            split_p.num_mismatches = split_p.midline.count('x')
            split_p.length = end - start
            split_p.distance = split_p.second_start - split_p.first_end

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


@jit(nopython=True, cache=True)
def _find_palindromes(seq, rev, m, N, D, Q):
    L = len(seq)

    palindrome_coords = List()
    palindrome_coords.append((0,0,0,0))

    i = 0
    while i < L-m+1:
        j = 0
        while j < L-i-m:
            skip_j_amount = 0
            if rev[j] != seq[i]:
                # The (i, j) scenario doesn't even _start_ with a match. Game over.
                pass
            else:
                n, k = 0, 0
                last_match = 0
                is_palindrome = False
                while True:
                    if (i+k+1) > (L-j-k-1):
                        # Stop of left has exceeded start of right.
                        if is_palindrome:
                            candidate = (i, i+last_match+1, L-j-last_match-1, L-j)
                            candidate = trim_ends(candidate, seq, rev, Q)
                            if candidate[1] - candidate[0] >= m:
                                # After trimming, candidate still meets required length
                                if not is_internal(candidate, palindrome_coords):
                                    # candidate is internal to another palindrome.
                                    palindrome_coords.append(candidate)
                                    skip_j_amount = last_match
                        break

                    if rev[j+k] == seq[i+k]:
                        last_match = k
                    else:
                        # mismatch
                        n += 1

                    if n > N:
                        if is_palindrome:
                            candidate = (i, i+last_match+1, L-j-last_match-1, L-j)
                            candidate = trim_ends(candidate, seq, rev, Q)
                            if candidate[1] - candidate[0] >= m:
                                # After trimming, candidate still meets required length
                                if not is_internal(candidate, palindrome_coords):
                                    # candidate is internal to another palindrome.
                                    palindrome_coords.append(candidate)
                                    skip_j_amount = last_match
                        break

                    if last_match == m-1:
                        is_palindrome = True

                    if (L-j-k-1) - (i+k+1) < D:
                        break

                    k += 1
            j += skip_j_amount + 1
        i += 1

    return palindrome_coords[1:]


@jit(nopython=True, cache=True)
def trim_ends(candidate, seq, rev, Q):
    left_start, left_end, right_start, right_end = candidate
    seq_left = seq[left_start:left_end]
    seq_right = rev[::-1][right_start:right_end][::-1]
    L = len(seq_left)

    clip_left = 0
    for q in range(Q+1):
        if seq_left[q] != seq_right[q]:
            clip_left = q + 1

    clip_right = 0
    for q in range(Q+1):
        if seq_left[L-1-q] != seq_right[L-1-q]:
            clip_right = q + 1

    left_start += clip_left
    right_end -= clip_left

    left_end -= clip_right
    right_start += clip_right

    return left_start, left_end, right_start, right_end


@jit(nopython=True, cache=True)
def is_internal(candidate, palindromes):
    a, b, c, d = candidate
    for left_start, left_end, right_start, right_end in palindromes:
        if a >= left_start and b <= left_end and c >= right_start and d <= right_end:
            return True
    else:
        return False


def _find_palindromes_test(seq, rev, m, N, D, Q):
    L = len(seq)
    palindrome_coords = [(0,0,0,0)]

    i = 0
    while i < L-m+1:
        j = 0
        while j < L-i-m:
            skip_j_amount = 0
            if rev[j] != seq[i]:
                get_state(seq, i, j, 0)
                print('No match between i and j')
                # The (i, j) scenario doesn't even _start_ with a match. Game over.
                pass
            else:
                n, k = 0, 0
                last_match = 0
                is_palindrome = False
                while True:
                    get_state(seq, i, j, k)
                    if (i+k+1) > (L-j-k-1):
                        print('Stop of left has exceeded start of right.')
                        # Stop of left has exceeded start of right.
                        if is_palindrome:
                            candidate = (i, i+last_match+1, L-j-last_match-1, L-j)
                            candidate = trim_ends(candidate, seq, rev, Q)
                            if candidate[1] - candidate[0] >= m:
                                print(candidate)
                                # After trimming, candidate still meets required length
                                if not is_internal(candidate, palindrome_coords):
                                    print(candidate)
                                    palindrome_coords.append(candidate)
                                    skip_j_amount = last_match
                                    print(f'Palindrome {candidate}.')
                                else:
                                    print(f'INTERNAL Palindrome {candidate}.')
                        break

                    if rev[j+k] == seq[i+k]:
                        last_match = k
                    else:
                        # mismatch
                        n += 1

                    if n > N:
                        print('Max # mismatches met. Calling quits')
                        if is_palindrome:
                            candidate = (i, i+last_match+1, L-j-last_match-1, L-j)
                            candidate = trim_ends(candidate, seq, rev, Q)
                            if candidate[1] - candidate[0] >= m:
                                # After trimming, candidate still meets required length
                                if not is_internal(candidate, palindrome_coords):
                                    palindrome_coords.append(candidate)
                                    skip_j_amount = last_match
                                    print(f'Palindrome {candidate}.')
                                else:
                                    print(f'INTERNAL Palindrome {candidate}.')
                        break

                    if last_match == m-1:
                        print('Min palindrome length met')
                        is_palindrome = True

                    if (L-j-k-1) - (i+k+1) < D:
                        print('Min gap distance not satisfied')
                        break

                    k += 1
            print(f'skipping j ahead {skip_j_amount}')
            j += skip_j_amount + 1

        i += 1

    return palindrome_coords[1:]


def get_state(seq, i, j, k):
    print()
    print(''.join([str(s) for s in seq]))
    state = [' ']*len(seq)
    state[i] = 'i'
    state[len(seq)-1-j] = 'j'
    state[k+i] = 'k'
    print(''.join(state))


