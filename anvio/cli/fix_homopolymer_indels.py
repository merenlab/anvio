#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import xml.etree.ElementTree as ET

import anvio
import anvio.terminal as terminal
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.drivers.blast import BLAST


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['fasta']
__provides__ = ['fasta']
__description__ = ("Corrects homopolymer-region associated INDELs in a given genome based on a reference genome. The "
                   "most effective use of this script is when the input genome is a genome reconstructed by minION "
                   "long reads, and the reference genome is one that is of high-quality. Essentially, this script "
                   "will BLAST the genome you wish to correct against the reference genome you provide, identify "
                   "INDELs in the BLAST results that are exclusively associated with homopolymer regions, and will "
                   "take the reference genome as a guide to correct the input sequences, and report a new "
                   "FASTA file. You can use the output FASTA file that is fixed as the input FASTA file over and over "
                   "again to see if you can eliminate all homopolymer-associated INDELs")


# some test sequences to play with. in each tuple, the first
# sequence is the reference, the second is the query (to be
# corrected), and the third is the expected sequence after
# correction. The first two must be aligned. The third may not
# need to match the reference
test_sequences = [('ATCGATCGATCG-AAATCGATCGATCG',
                   'ATCGATCGATCGAAAATCGATCGATCG',
                   'ATCGATCGATCGAAATCGATCGATCG'),

                  ('ATCGATCGATCGAA-TCGATCGATCG',
                   'ATCGATCGATCGAAAATCGATCGATCG',
                   'ATCGATCGATCGAAATCGATCGATCG'),

                  ('ATCGATCGATCGAAAATCGATCGATCG',
                   'ATCGATCGATCGAAA-TCGATCGATCG',
                   'ATCGATCGATCGAAAATCGATCGATCG'),

                  ('ATCGATCGATCGAAAATCGATCGATCG',
                   'A-C-AT-GATCG-AAATCGATCGATCG',
                   'ACATGATCGAAAATCGATCGATCG'),

                  ('ATCGATCGATCGAA-TCGATCGATCG',
                   'ATCGATCGATCGAAATCGATCGATCG',
                   'ATCGATCGATCGAATCGATCGATCG'),

                  ]


def main():
    args = get_args()

    try:
        c = RefBasedHomopolymerINDELCorrector(args)
        c.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


class RefBasedHomopolymerINDELCorrector:
    def __init__(self, args, skip_sanity_check=False, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.progress = p
        self.run = r

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_fasta = A('input_fasta')
        self.reference_fasta = A('reference_fasta')
        self.output_fasta = A('output_fasta')
        self.num_threads = A('num_threads')
        self.verbose = A('verbose')
        self.test_run = A('test_run')
        self.min_homopolymer_length = int(A('min_homopolymer_length')) or 3

        self.resolutions_dict = {}
        self.actions_by_nucleotide = {}
        self.stats = {}

        if self.test_run:
            # don't do anything, but run our tests
            self.verbose = True
            self.do_test_run()
            sys.exit()

        self.blast_output_path = 'blast_output.xml'

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        if not self.input_fasta or not self.reference_fasta or not self.output_fasta:
            raise ConfigError("You need to provide an input, reference, and output FASTA "
                              "file names here :/")

        filesnpaths.is_file_fasta_formatted(self.input_fasta)
        filesnpaths.is_file_fasta_formatted(self.reference_fasta)
        filesnpaths.is_output_file_writable(self.output_fasta)

        try:
            self.num_threads = int(self.num_threads)
            assert(self.num_threads > 0)
        except:
            raise ConfigError("Num threads must be a positive integer. Obviously.")


    def process(self):
        """Do the thing."""

        # warn the user
        self.run.warning("You must be extremely careful with this program since it reports edited sequences. In addition, "
                         "if you need a comprehensive solution to correct frameshift errors in your long-read sequencing "
                         "data for mission critical applications, you should NOT use this script. To find out about other "
                         "solutions that may be more suitale for your needs, please see the online help for thsi script "
                         "(the URL for the online help will appear at the end of the `--help` output).")

        # let's BLAST all all input sequences against reference sequences
        self.gen_blast_output()

        # this will fill up self.resolutions_dict
        self.parse_blast_output()

        # this will report a FASTA file with corrected seqeunces in input FASTA
        self.edit_input_sequences()

        # clean after yourself unless --debug
        if not anvio.DEBUG:
            for f in ['blast-log-file.txt', 'blast_output.xml', self.reference_fasta + '.nin', self.reference_fasta + '.nhr', self.reference_fasta + '.nsq']:
                try:
                    os.remove(f)
                except:
                    pass


    def edit_sequence(self, input_seq, actions_by_nucleotide):
        """Takes in a sequence and a dict of actions by nucleotide position.

        Parameters
        ==========
        input_seq : str
            The original sequence to be edited.
        actions_by_nucleotide : dict of tuples
            That looks like this: {5: ('DEL', 'C'), 13: ('DEL', 'G'), 20: ('INS', 'G'), ...}
        """

        # we will build an edited version of the `input_seq` into this variable:
        edited_sequence = ""

        num_gaps_in_input_seq = input_seq.count('-')
        if num_gaps_in_input_seq:
            raise ConfigError(f"The input sequence to be edited seems to contain gaps (anvi'o found {num_gaps_in_input_seq}). For "
                              f"this program to work, the original input sequences should not contain any artifacts of alignment.")

        list_of_actionable_positions = sorted(list(actions_by_nucleotide.keys()))

        index = 0
        while len(list_of_actionable_positions):
            action_position = list_of_actionable_positions.pop(0)

            edited_sequence += input_seq[index:action_position]

            index = action_position

            if actions_by_nucleotide[action_position][0] == 'DEL':
                index += 1
            elif actions_by_nucleotide[action_position][0] == 'INS':
                edited_sequence += actions_by_nucleotide[action_position][1]
                index += 0
            else:
                raise ConfigError(f"Your `actions_by_nucleotide` dict includes an unknown action "
                                  f"called '{actions_by_nucleotide[action_position][0]}'. Try "
                                  f"'DEL' or 'INS'.")

        # add the remainder
        edited_sequence += input_seq[index:]

        return edited_sequence


    def populate_actions_by_nucleotide(self):
        """Populates `self.actions_by_nucleotide` from `self.resolutions_dict`."""

        if not self.resolutions_dict:
            raise ConfigError("There is nothing in resolutions dict to work with.")

        for sequence_name in self.resolutions_dict:
            self.actions_by_nucleotide[sequence_name] = {}
            self.stats[sequence_name] = {'num_homopolymers': 0, 'num_insertions': 0, 'num_deletions': 0}

            for resolution in self.resolutions_dict[sequence_name]:
                self.stats[sequence_name]['num_homopolymers'] += 1

                for nucleotide_position in resolution['positions']:
                    self.actions_by_nucleotide[sequence_name][nucleotide_position] = (resolution['action'], resolution['nt'])

                    if resolution['action'] == 'DEL':
                        self.stats[sequence_name]['num_deletions'] += 1
                    elif resolution['action'] == 'INS':
                        self.stats[sequence_name]['num_insertions'] += 1


    def edit_input_sequences(self):
        if not self.actions_by_nucleotide:
            raise ConfigError("There is nothing in the actions by nucleotide dict to work with.")

       # time to go throught he input sequences, and format them
        input_fasta = fastalib.SequenceSource(self.input_fasta)

        with open(self.output_fasta, 'w') as output_fasta:
            while next(input_fasta):
                if input_fasta.id in self.actions_by_nucleotide and len(self.actions_by_nucleotide[input_fasta.id]):
                    edited_sequence = self.edit_sequence(input_fasta.seq, self.actions_by_nucleotide[input_fasta.id])

                    output_fasta.write(f">{input_fasta.id}\n")
                    output_fasta.write(f"{edited_sequence}\n")
                else:
                    output_fasta.write(f">{input_fasta.id}\n")
                    output_fasta.write(f"{input_fasta.seq}\n")

        # lets gather some stats
        total_number_of_homopolymers = sum([e['num_homopolymers'] for e in self.stats.values()])
        total_number_of_deletions = sum([e['num_deletions'] for e in self.stats.values()])
        total_number_of_insertions = sum([e['num_insertions'] for e in self.stats.values()])
        total_number_of_actions = total_number_of_deletions + total_number_of_insertions

        self.run.warning(None, header="OVERALL & PER-SEQUENCE STATS")
        self.run.info('Num input sequences', len(self.resolutions_dict))
        self.run.info('Num homopolymers associated with INDELs', total_number_of_homopolymers)
        self.run.info('Num actions', total_number_of_actions)
        self.run.info('Num insertions', total_number_of_insertions)
        self.run.info('Num deletions', total_number_of_deletions, nl_after=1)

        for sequence_name in self.stats:
            self.run.info_single(f'{sequence_name}')
            self.run.info_single(f'Homopolymers associated with INDELs: {self.stats[sequence_name]["num_homopolymers"]}', level=2)
            self.run.info_single(f'Insertions: {self.stats[sequence_name]["num_insertions"]}', level=2)
            self.run.info_single(f'Deletions: {self.stats[sequence_name]["num_deletions"]}', level=2, nl_after=1)

        self.run.info('Corrected output FASTA', self.output_fasta, mc='green', nl_after=1)


    def gen_blast_output(self):
        blast = BLAST(self.input_fasta, target_fasta=self.reference_fasta)
        blast.search_program = 'blastn'
        blast.search_output_path = self.blast_output_path
        blast.makedb(dbtype='nucl')
        blast.num_threads = self.num_threads
        blast.max_target_seqs = 1
        blast.blast(outputfmt='5')


    def parse_blast_output(self):
        """Goes through each HSP, and solves homopolymer crimes"""

        root = ET.parse(self.blast_output_path).getroot()

        for query_sequence_xml in root.findall('BlastOutput_iterations/Iteration'):

            query_sequence_name = query_sequence_xml.find('Iteration_query-def').text

            if query_sequence_name not in self.resolutions_dict:
                self.resolutions_dict[query_sequence_name] = []

            for hit_xml in query_sequence_xml.findall('Iteration_hits/Hit'):

                hit_num =int(hit_xml.find('Hit_num').text)

                for hsp_xml in hit_xml.findall('Hit_hsps/Hsp'):
                    hsp_num = int(hsp_xml.find('Hsp_num').text)
                    num_gaps = int(hsp_xml.find('Hsp_gaps').text)

                    if not num_gaps:
                        continue

                    hsp_start = int(hsp_xml.find('Hsp_query-from').text)
                    hsp_query = hsp_xml.find('Hsp_qseq').text
                    hsp_ref = hsp_xml.find('Hsp_hseq').text

                    resolutions = self.get_resolutions_for_homopolymers(hsp_ref, hsp_query, query_sequence_name=query_sequence_name, hit_num=hit_num, hsp_num=hsp_num)

                    # here we need to update the relative start position with the absolute one
                    for resolution in resolutions:
                        resolution['positions'] = [p + hsp_start for p in resolution['positions']]

                    # now we are ready to add this in the global dict
                    self.resolutions_dict[query_sequence_name].extend(resolutions)

        # we now have a resolutions dict populated with explicit nucleotide positions. time to populate
        # actions_by_nucleotide dict
        self.populate_actions_by_nucleotide()


    def do_test_run(self):
        """Runs stuff on `test_sequences` and prints messages"""
        counter = 0

        for ref_seq, query_seq, expected_seq in test_sequences:
            seq_id = f'test_sequence_{counter}'

            self.run.info_single('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>', nl_before=1, nl_after=1, mc='red')

            self.run.info('Query sequence', query_seq)
            self.run.info('Reference sequence', ref_seq)

            self.resolutions_dict[seq_id] = self.get_resolutions_for_homopolymers(ref_seq, query_seq)

            self.populate_actions_by_nucleotide()

            edited_sequence = self.edit_sequence(query_seq.replace('-', ''), self.actions_by_nucleotide[seq_id])

            self.run.info('Edited query sequence', edited_sequence)
            self.run.info('Query sequence edited correctly?', 'Yes' if edited_sequence == expected_seq else 'No')

            counter += 1


    def get_resolutions_for_homopolymers(self, ref_seq, query_seq, query_sequence_name=None, hit_num=None, hsp_num=None):
        resolutions = []

        # variables of convenience
        alignment_length = len(ref_seq)

        # the following variable is necessary to convert nucleotide positions
        # from relative positions with aligned characters (since gaps add extra)
        # characters, the index counting actual number of characters will not
        # match to nucleotides in original, unaligned sequence later. instead of
        # setting it to zero, we will set it to the number of gaps observed between
        # 0 and min_homopolymer_length, since otherwise gaps that occur before
        # min_homopolymer_length will not be included in our position calculations
        # and will screw up everything downstream
        num_gaps_observed_in_query = query_seq[0:self.min_homopolymer_length].count('-')

        # we start counting from self.min_homopolymer_length, so we can work with downstream bases
        # at the very beginning of the sequence and upstream bases at the very end of it without
        # getting key/index out of bound errors
        index = self.min_homopolymer_length

        while True:
            gap_in_ref = ref_seq[index] == '-'
            gap_in_query = query_seq[index] == '-'

            if gap_in_ref and gap_in_query:
                raise ConfigError("Something weird happened. You have an alignment where both the ref and query "
                                  "has a gap character.")

            if not gap_in_ref and not gap_in_query:
                index += 1
                if index > alignment_length - self.min_homopolymer_length:
                    break
                else:
                    continue

            # if we are here, it means there is at least one gap character

            # this trick will allow us to work with our sequences identically regardless of
            # which sequence has the gap. we will check whether it is the ref that has a gap
            # or the query right before we finalize our resolution
            if gap_in_ref:
                sequence_with_deletions, sequence_with_insertions  = ref_seq, query_seq
            else:
                sequence_with_deletions, sequence_with_insertions  = query_seq, ref_seq

            gap_start = index
            gap_length = 1

            # find out if there are multiple gaps following each other:
            while True:
                if sequence_with_deletions[gap_start + gap_length] == '-':
                    gap_length += 1
                else:
                    break

            # so we have a gap that is 1 or more nts long in either ref or query. here we will make note of the
            # surrounding context for the conditionals below to find out what the hell is going on.
            upstream_nts_until_gap_in_gapped_seqeunce = sequence_with_deletions[gap_start - self.min_homopolymer_length:gap_start]
            upstream_nts_through_gap_in_the_other_sequence = sequence_with_insertions[gap_start - self.min_homopolymer_length:gap_start + gap_length]
            downstream_nts_after_gap_in_gapped_seqeunce = sequence_with_deletions[gap_start + gap_length:gap_start + gap_length + self.min_homopolymer_length]
            downstream_nts_through_gap_in_the_other_sequence = sequence_with_insertions[gap_start:gap_start + gap_length + self.min_homopolymer_length]

            if len(set(upstream_nts_until_gap_in_gapped_seqeunce)) == 1 and \
               len(set(upstream_nts_through_gap_in_the_other_sequence)) == 1 and \
               len(set(upstream_nts_until_gap_in_gapped_seqeunce + upstream_nts_through_gap_in_the_other_sequence)) == 1:
                """If we are here, it means we are handling either of these situations:

                     >>> one sequence : AAA--
                     >>> the other one: AAAAA
                """

                repeat_nucleotide = upstream_nts_until_gap_in_gapped_seqeunce[0]

                # here we determine the resolution for our query sequence
                if gap_in_ref:
                    resolution = {'action': 'DEL', 'positions': [p - num_gaps_observed_in_query for p in range(gap_start, gap_start + gap_length)], 'nt': repeat_nucleotide}
                else:
                    resolution = {'action': 'INS', 'positions': [p - num_gaps_observed_in_query for p in range(gap_start, gap_start + gap_length)], 'nt': repeat_nucleotide}

                resolutions.append(resolution)

                if self.verbose:
                    self.run.warning(None, header=f'GAPS AFTER REPEATS OF "{repeat_nucleotide}"', lc='yellow')
                    self.run.info('Query sequence name', query_sequence_name, lc='yellow', mc='yellow') if query_sequence_name else None
                    self.run.info('Hit num', hit_num, lc='yellow', mc='yellow') if hit_num else None
                    self.run.info('Hsp num', hsp_num, lc='yellow', mc='yellow') if hsp_num else None

                    if gap_in_ref:
                        self.run.info('Query', upstream_nts_through_gap_in_the_other_sequence, lc='yellow', mc='yellow')
                        self.run.info('Reference', upstream_nts_until_gap_in_gapped_seqeunce + '-' * gap_length, lc='yellow', mc='yellow')
                    else:
                        self.run.info('Query', upstream_nts_until_gap_in_gapped_seqeunce + '-' * gap_length, lc='yellow', mc='yellow')
                        self.run.info('Reference', upstream_nts_through_gap_in_the_other_sequence, lc='yellow', mc='yellow')

                    self.run.info('Resolution', str(resolution), nl_after=1, lc='yellow', mc='yellow')

            elif len(set(downstream_nts_after_gap_in_gapped_seqeunce)) == 1 and \
                 len(set(downstream_nts_through_gap_in_the_other_sequence)) == 1 and \
                 len(set(downstream_nts_after_gap_in_gapped_seqeunce + downstream_nts_through_gap_in_the_other_sequence)) == 1:
                """If we are here, it means we are handling either of these situations:

                     >>> one sequence : --AAA
                     >>> the other one: AAAAA
                """

                repeat_nucleotide = downstream_nts_after_gap_in_gapped_seqeunce[0]

                # that time of the code again
                if gap_in_ref:
                    resolution = {'action': 'DEL', 'positions': [p - num_gaps_observed_in_query for p in range(gap_start, gap_start + gap_length)], 'nt': repeat_nucleotide}
                else:
                    resolution = {'action': 'INS', 'positions': [p - num_gaps_observed_in_query for p in range(gap_start, gap_start + gap_length)], 'nt': repeat_nucleotide}

                resolutions.append(resolution)

                if self.verbose:
                    self.run.warning(None, header=f'GAPS BEFORE REPEATS OF "{repeat_nucleotide}"', lc='cyan')
                    self.run.info('Query sequence name', query_sequence_name, lc='yellow', mc='yellow') if query_sequence_name else None
                    self.run.info('Hit num', hit_num, lc='yellow', mc='yellow') if hit_num else None
                    self.run.info('Hsp num', hsp_num, lc='yellow', mc='yellow') if hsp_num else None

                    if gap_in_ref:
                        self.run.info('Query', downstream_nts_through_gap_in_the_other_sequence, lc='cyan', mc='cyan')
                        self.run.info('Reference', '-' * gap_length + downstream_nts_after_gap_in_gapped_seqeunce, lc='cyan', mc='cyan')
                    else:
                        self.run.info('Query', '-' * gap_length + downstream_nts_after_gap_in_gapped_seqeunce, lc='cyan', mc='cyan')
                        self.run.info('Reference', downstream_nts_through_gap_in_the_other_sequence, lc='cyan', mc='cyan')

                    self.run.info('Resolution', str(resolution), nl_after=1, lc='cyan', mc='cyan')
                else:
                    """So this is just a lonely regions of gaps without any any homopolymer region surrounding them"""
                    pass

            index += gap_length

            if gap_in_query:
                num_gaps_observed_in_query += gap_length

            if index > alignment_length - self.min_homopolymer_length:
                break

        return resolutions


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('FILES THAT MATTER', "Here you provide file paths for input sequence(s) to be corrected, "
                                       "reference sequence(s) to be used for correction, and the edited input sequences to be stored. "
                                       "UNLESS you just want to do a test run. In which case you don't need any of these but the "
                                       "parameter `--test-run`")
    groupA.add_argument('-i', '--input-fasta', metavar="FASTA", help="A FASTA file of sequences you wish to fix")
    groupA.add_argument('-r', '--reference-fasta', metavar="FASTA", help="A FASTA file for reference sequences")
    groupA.add_argument('-o', '--output-fasta', metavar="FASTA", help="Corrected FASTA file")

    groupB = parser.add_argument_group('STUFF UNDER THE HOOD', "Like how should we define homopolymers or how much information should "
                                        "we share with you as we go.")
    groupB.add_argument('--min-homopolymer-length', metavar="INT", type=int, default=3, help="Minimum number of "
                            "identical nucleotides next to each other PLUS THE GAP CHARACTER to be considered a homopolymer when these "
                            "nucleotides are aligned to a region in the other sequnce that is all composed of the same nucleotides. Confused? "
                            "Read on. The default is %(default)s. So when this value is 2, the program would consider the following to match the "
                            "definition of minimum homopolymer length to be considered for fixing: (R)eference: 'AA-' and (Q)uery: 'AAA'. "
                            "The same would be true for R: 'AA---' / Q: 'AAAAA' but not R: 'A-' / Q: 'AA' In contrast, when this value is 3, "
                            "then the minimum that would work would be this: R: 'AAA-', Q: 'AAAA'. Obviously, you shouldn't go any lower "
                            "than 2, but then why should you listen to a computer?")
    groupB.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))

    groupC = parser.add_argument_group('PERFORMANCE', "For the BLAST search")
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    groupD = parser.add_argument_group('TEST RUN', "To have an idea about what is corrected and what is not")
    groupD.add_argument('--test-run', default=False, action="store_true", help="Just do a test run and nothing more.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

