# -*- coding: utf-8
"""Parser for HMMER's various outputs"""

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.parsers.base import Parser

import os
import numpy as np
import pandas as pd


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class HMMERStandardOutput(object):
    """Parse the standard output of HMMER programs (NOTE: currently only works with hmmsearch)

    The main meat of this class is to produce the attributes:

        (1) self.seq_hits
        (2) self.dom_hits
        (3) self.ali_info

    (1) self.seq_hits is a dataframe that looks like this:

        |              query         acc target  query_len        evalue  score  bias  \
        | 0       3Beta_HSD  PF01073.18   1998        282  5.200000e-23   76.2   0.0
        | 1       3Beta_HSD  PF01073.18   1723        282  1.300000e-07   25.7   0.0
        | ...           ...         ...    ...        ...           ...    ...   ...
        | 3128  Voltage_CLC  PF00654.19    320        354  7.200000e-65  214.3  37.1
        | 3129         YkuD  PF03734.13     30        146  1.700000e-14   49.3   0.2

        |       best_dom_evalue  best_dom_score  best_dom_bias  expected_doms  num_doms
        | 0        6.600000e-22            72.6            0.0            2.0         1
        | 1        1.700000e-07            25.3            0.0            1.2         1
        | ...               ...             ...            ...            ...       ...
        | 3128     7.800000e-64           210.9           29.1            2.0         1
        | 3129     3.800000e-14            48.2            0.2            1.7         1

    (2) self.dom_hits is a frame that looks like this:

        |               query         acc target  domain qual  score  bias      c-evalue  \
        | 0       3Beta_HSD  PF01073.18   1998       1    !   72.6   0.0  2.900000e-24
        | 1       3Beta_HSD  PF01073.18   1723       1    !   25.3   0.0  7.300000e-10
        | ...           ...         ...    ...     ...  ...    ...   ...           ...
        | 2896  Voltage_CLC  PF00654.19    320       1    !  210.9  29.1  1.700000e-66
        | 2897         YkuD  PF03734.13     30       1    !   48.2   0.2  8.400000e-17
        |
        |           i-evalue  hmm_start  hmm_stop hmm_bounds  ali_start  ali_stop  \
        | 0     6.600000e-22          1       237         [.          4       243
        | 1     1.700000e-07          1        95         [.          4        92
        | ...            ...        ...       ...        ...        ...       ...
        | 2896  7.800000e-64          3       352         ..         61       390
        | 2897  3.800000e-14          2       146         .]        327       459
        |
        |      ali_bounds  env_start  env_stop env_bounds  mean_post_prob  \
        | 0            ..          4       254         ..            0.74
        | 1            ..          4       148         ..            0.72
        | ...         ...        ...       ...        ...             ...
        | 2896         ..         59       392         ..            0.94
        | 2897         ..        326       459         ..            0.78
        |
        |       match_state_align              comparison_align             sequence_align
        | 0     vvtGggGFlGrrivkeLlrl...  +v+Gg+G++G++ v +L++ ...  LVLGGAGYIGSHAVDQLISK...
        | 1     vvtGggGFlGrrivkeLlrl...  ++ Gg+GFlG++i k L+++...  IIFGGSGFLGQQIAKILVQR...
        | ...                       ...                      ...                      ...
        | 2896  gllagllvkrvapeaagsGi...  g++  +++ r+  + a  G ...  GVVFTYFYTRF-GKNASRGN...
        | 2897  kyivvdlaeqrllvlyengk...  +yi++dl++q++ +++ +gk...  NYIEIDLKDQKM-YCFIDGK...

    If you're confused about the meaning of these columns, please see starting from page 32
    of the HMMER guide http://eddylab.org/software/hmmer/Userguide.pdf. There you will be able
    to with relative ease correlate the column names in these tables to what is described
    meticulously in the tutorial. For example, `best_dom_bias` refers to the 'bias (best 1
    domain)' column.

    (3) ali_info is a nested dictionary that can be used to access on a per-hit basis which residues
        in a sequence aligned to which residues in the HMM.

    Parameters
    ==========
    hmmer_std_out : str
        Path to output of HMMER.

    context : str, None
        If provided, operations specific to a context will also be carried out. Choose from
        {'interacdome'}
    """

    def __init__(self, hmmer_std_out, context=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        self.hmmer_std_out = hmmer_std_out
        self.context = context

        self.set_names()

        self.ali_info = {}

        # This is converted to a dataframe after populating
        self.seq_hits = {
            self.query_col: [],
            self.acc_col: [],
            self.target_col: [],
            self.query_len_col: [],
            'evalue': [],
            'score': [],
            'bias': [],
            'best_dom_evalue': [],
            'best_dom_score': [],
            'best_dom_bias': [],
            'expected_doms': [],
            'num_doms': [],
        }

        self.seq_hits_dtypes = {
            self.query_col: str,
            self.acc_col: str,
            self.target_col: str,
            self.query_len_col: int,
            'evalue': float,
            'score': float,
            'bias': float,
            'best_dom_evalue': float,
            'best_dom_score': float,
            'best_dom_bias': float,
            'expected_doms': float,
            'num_doms': int,
        }

        # This is converted to a dataframe after populating
        self.dom_hits = {
            self.query_col: [],
            self.acc_col: [],
            self.target_col: [],
            'domain': [],
            'qual': [],
            'score': [],
            'bias': [],
            'c-evalue': [],
            'i-evalue': [],
            'hmm_start': [],
            'hmm_stop': [],
            'hmm_bounds': [],
            'ali_start': [],
            'ali_stop': [],
            'ali_bounds': [],
            'env_start': [],
            'env_stop': [],
            'env_bounds': [],
            'mean_post_prob': [],
            'match_state_align': [],
            'comparison_align': [],
            'sequence_align': [],
        }

        self.dom_hits_dtypes = {
            self.query_col: str,
            self.acc_col: str,
            self.target_col: str,
            'domain': int,
            'qual': str,
            'score': float,
            'bias': float,
            'c-evalue': float,
            'i-evalue': float,
            'hmm_start': int,
            'hmm_stop': int,
            'hmm_bounds': str,
            'ali_start': int,
            'ali_stop': int,
            'ali_bounds': str,
            'env_start': int,
            'env_stop': int,
            'env_bounds': str,
            'mean_post_prob': float,
            'match_state_align': str,
            'comparison_align': str,
            'sequence_align': str,
        }

        self.delim_query = '//\n'
        self.delim_seq = '>>'
        self.delim_domain = '=='

        self.load()


    def load(self):
        self.progress.new('Processing HMMER output')
        self.progress.update('Parsing %s' % self.hmmer_std_out)

        with open(self.hmmer_std_out) as f:
            for i, query in enumerate(utils.get_chunk(f, separator=self.delim_query, read_size=32768)):

                if i % 500 == 0:
                    self.progress.update('%d done' % i)
                    self.progress.increment(increment_to=i)

                self.process_query(query)

        self.seq_hits = pd.DataFrame(self.seq_hits).astype(self.seq_hits_dtypes)
        self.dom_hits = pd.DataFrame(self.dom_hits).astype(self.dom_hits_dtypes)
        self.progress.end()

        self.additional_processing()

        self.run.info('Loaded HMMER results from', self.hmmer_std_out)


    def find_line(self, condition):
        for line in self.query_lines[self.line_no:]:
            self.line_no += 1

            if line.startswith('#'):
                continue

            if condition(line):
                return line
        else:
            return False


    def read_lines_until(self, condition, include_last=False, store=True):
        lines = []
        return_value = lines if store else True

        for line in self.query_lines[self.line_no:]:
            self.line_no += 1

            if line.startswith('#'):
                continue

            if condition(line):
                if include_last and store:
                    lines.append(line)

                return lines

            if store:
                lines.append(line)
        else:
            if store:
                return lines
            else:
                return False


    def process_query(self, query):
        if self.delim_seq not in query:
            # This query had no hits
            return

        self.query_lines = query.split('\n')
        self.line_no = 0

        line = self.find_line(lambda line: line.startswith('Query:'))
        line_split = line.split()
        query_name = line_split[1]
        query_len = int(line_split[2][line_split[2].find('=')+1:-1])

        line = self.find_line(lambda line: line.startswith('Accession:'))
        acc = line.split()[1]

        line = self.find_line(lambda line: line.lstrip().startswith('E-value'))
        description_index = line.find('Desc')
        fields = line[:description_index].split() # ignore last 'Description' field

        assert len(fields) == 9, "Please report this on github with your HMMER version"

        self.read_lines_until(lambda line: line.lstrip().startswith('-------'), store=False)
        seq_score_lines = self.read_lines_until(lambda line: line == '')

        num_doms_per_seq = {}

        for seq_score_line in seq_score_lines:
            seq_scores = seq_score_line[:description_index].split()

            self.seq_hits[self.query_col].append(query_name)
            self.seq_hits[self.query_len_col].append(query_len)
            self.seq_hits[self.acc_col].append(acc)
            self.seq_hits['evalue'].append(float(seq_scores[0]))
            self.seq_hits['score'].append(float(seq_scores[1]))
            self.seq_hits['bias'].append(float(seq_scores[2]))
            self.seq_hits['best_dom_evalue'].append(float(seq_scores[3]))
            self.seq_hits['best_dom_score'].append(float(seq_scores[4]))
            self.seq_hits['best_dom_bias'].append(float(seq_scores[5]))
            self.seq_hits['expected_doms'].append(float(seq_scores[6]))
            self.seq_hits['num_doms'].append(int(seq_scores[7]))
            self.seq_hits[self.target_col].append(seq_scores[8])

            num_doms_per_seq[seq_scores[8]] = int(seq_scores[7])

        num_seq_hits = len(seq_score_lines)

        for _ in range(num_seq_hits):
            target_name = self.find_line(lambda line: line.startswith(self.delim_seq)).split()[1]

            if num_doms_per_seq[target_name] == 0:
                continue

            self.line_no += 2
            for __ in range(num_doms_per_seq[target_name]):
                dom_score_summary = self.find_line(lambda line: True).split()

                self.dom_hits[self.query_col].append(query_name)
                self.dom_hits[self.acc_col].append(acc)
                self.dom_hits[self.target_col].append(target_name)
                self.dom_hits['domain'].append(dom_score_summary[0])
                self.dom_hits['qual'].append(dom_score_summary[1])
                self.dom_hits['score'].append(dom_score_summary[2])
                self.dom_hits['bias'].append(dom_score_summary[3])
                self.dom_hits['c-evalue'].append(dom_score_summary[4])
                self.dom_hits['i-evalue'].append(dom_score_summary[5])
                self.dom_hits['hmm_start'].append(dom_score_summary[6])
                self.dom_hits['hmm_stop'].append(dom_score_summary[7])
                self.dom_hits['hmm_bounds'].append(dom_score_summary[8])
                self.dom_hits['ali_start'].append(dom_score_summary[9])
                self.dom_hits['ali_stop'].append(dom_score_summary[10])
                self.dom_hits['ali_bounds'].append(dom_score_summary[11])
                self.dom_hits['env_start'].append(dom_score_summary[12])
                self.dom_hits['env_stop'].append(dom_score_summary[13])
                self.dom_hits['env_bounds'].append(dom_score_summary[14])
                self.dom_hits['mean_post_prob'].append(dom_score_summary[15])

            for __ in range(num_doms_per_seq[target_name]):
                self.find_line(lambda line: line.lstrip().startswith(self.delim_domain))

                if __ == num_doms_per_seq[target_name] - 1:
                    if _ == num_seq_hits - 1:
                        # This is the last alignment in the summary_info. Go to end of string
                        ali_lines = self.read_lines_until(lambda line: False)
                    else:
                        # This is the last alignment in the sequence. Go to next sequence delimiter
                        ali_lines = self.read_lines_until(lambda line: line.lstrip().startswith(self.delim_seq))
                        self.line_no -= 1
                else:
                    ali_lines = self.read_lines_until(lambda line: line.lstrip().startswith(self.delim_domain))
                    self.line_no -= 1

                consensus = []
                match = []
                target = []
                line_index = 0
                while True:
                    if line_index >= len(ali_lines):
                        break

                    line = ali_lines[line_index]

                    if not line.lstrip().startswith(query_name + ' '):
                        line_index += 1
                        continue

                    cons_seq_fragment = line.split()[2]
                    frag_len = len(cons_seq_fragment)
                    ali_index = line.find(cons_seq_fragment)

                    consensus.append(cons_seq_fragment)
                    match.append(ali_lines[line_index + 1][ali_index: ali_index + frag_len])
                    target.append(ali_lines[line_index + 2][ali_index: ali_index + frag_len])

                    line_index += 2

                self.dom_hits['match_state_align'].append(''.join(consensus))
                self.dom_hits['comparison_align'].append(''.join(match))
                self.dom_hits['sequence_align'].append(''.join(target))


    def set_names(self):
        """Set the column names depending on self.context"""

        if self.context is None:
            self.query_col = 'query'
            self.acc_col = 'acc'
            self.query_len_col = 'query_len'
            self.target_col = 'target'

        elif self.context == 'interacdome':
            self.query_col = 'pfam_name'
            self.acc_col = 'pfam_id'
            self.query_len_col = 'pfam_len'
            self.target_col = 'corresponding_gene_call'


    def additional_processing(self):
        """Further process raw data"""

        if self.context is None:
            self.get_ali_info()

        elif self.context == 'interacdome':
            self.seq_hits['corresponding_gene_call'] = self.seq_hits['corresponding_gene_call'].astype(int)
            self.dom_hits['corresponding_gene_call'] = self.dom_hits['corresponding_gene_call'].astype(int)

            if self.dom_hits.empty:
                self.dom_hits['version'] = []
            else:
                self.dom_hits[['pfam_id', 'version']] = self.dom_hits['pfam_id'].str.split('.', n=1, expand=True)

            if self.seq_hits.empty:
                self.seq_hits['version'] = []
            else:
                self.seq_hits[['pfam_id', 'version']] = self.seq_hits['pfam_id'].str.split('.', n=1, expand=True)

            # For convenience this is done after pfam_id has been split
            self.get_ali_info()


    def get_ali_info(self):
        """Creates self.ali_info. See class docstring for description

        Notes
        =====
        - This function is very slow.
        - EDIT: This function is not _that_ slow
        """

        if self.dom_hits.empty:
            return

        unique_targets = self.dom_hits[self.target_col].nunique()
        self.progress.new('Processing alignment info', progress_total_items=unique_targets)

        gap_chars = {'-', '.'}

        processed = 0
        for target, subset in self.dom_hits.groupby(self.target_col):
            if processed % 50 == 0:
                self.progress.update('%d/%d done' % (processed, unique_targets))
                self.progress.increment(increment_to=processed)

            self.ali_info[target] = {}

            for acc, subsubset in subset.groupby(self.acc_col):
                for i, row in subsubset.iterrows():
                    seq_positions, seq_chars, hmm_positions, hmm_chars, comparison_chars = [], [], [], [], []

                    seq_pos, hmm_pos = row['ali_start'], row['hmm_start']
                    sequence, match_state, comparison = row['sequence_align'], row['match_state_align'], row['comparison_align']

                    assert len(sequence) == len(match_state)

                    for i in range(len(sequence)):
                        seq_char, hmm_char, comparison_char = sequence[i], match_state[i], comparison[i]
                        if (seq_char not in gap_chars) and (hmm_char not in gap_chars):
                            # there is alignment (non-gap characters)
                            seq_positions.append(seq_pos)
                            seq_chars.append(seq_char)
                            hmm_positions.append(hmm_pos)
                            hmm_chars.append(hmm_char.upper())
                            comparison_chars.append(comparison_char.upper())
                            seq_pos += 1
                            hmm_pos += 1
                        elif (seq_char in gap_chars) and (hmm_char not in gap_chars):
                            # gap in seq
                            hmm_pos += 1
                        elif (seq_char not in gap_chars) and (hmm_char in gap_chars):
                            # gap in match state
                            seq_pos += 1
                        else:
                            # this happens with 0 probability
                            pass

                    # The HMM state and sequence positions are 1-indexed. We subtract by 1 to make them zero-indexed
                    self.ali_info[target][(acc, row['domain'])] = pd.DataFrame({
                        'seq': seq_chars,
                        'hmm': hmm_chars,
                        'comparison': comparison_chars,
                        'seq_positions': np.array(seq_positions) - 1,
                        'hmm_positions': np.array(hmm_positions) - 1,
                    })

            processed += 1
        self.progress.end()



class HMMERTableOutput(Parser):
    """Parse --tblout or --domtblout output formats for hmmer programs

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    <rant>
    Parsing of HMMER tabular output needs to be redesigned. See https://github.com/merenlab/anvio/issues/1670
    for details
    </rant>

    Output specifictions of HMMER can be found in the user guide. At time of writing this,
    http://eddylab.org/software/hmmer/Userguide.pdf hosts the user guide.

    Parameters
    ==========
    hmmer_table_txt: str
        Path to the tabular output file generated by the anvio.driver.HMMER.run_hmmer() function. This is not direct output from HMMER, but
        a slightly pre-processed version (header lines removed, and output from multiple threads, if any, are concatenated)

    alphabet: str, 'AA'
        Which alphabet do the HMMs use? Pick from {'AA', 'DNA', 'RNA'}

    context: str, 'GENE'
        This tells the class how the output should be parsed. Pick from {'GENE', 'CONTIG',
        'DOMAIN'}. Before being preprocessed by anvio.driver.HMMER.hmmscan_worker (see this module's
        docstring), the header of the file should look like so, based on which context you use:

        GENE:
            #                                                               |-- full sequence ---| |-- best 1 domain ---| |-- domain number estimation ---|
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description
            #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- -----------

        DOMAIN:
            #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
            # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description
            #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- -----------

        CONTIG:
            #
            # target name        accession  query name           accession  hmmfrom hmm to  alifrom   ali to  envfrom   env to   modlen strand   E-value  score  bias  description of target
            #------------------- ---------- -------------------- ---------- ------- -------  -------  -------  -------  -------  ------- ------ --------- ------ ----- ---------------------

        `DOMAIN` is untested.

    program: str, 'hmmscan'
        Which HMMER program was used to generate the output we are parsing? Pick from {'hmmscan', 'hmmsearch'}
    """

    def __init__(self, hmmer_table_txt, alphabet='AA', context='GENE', program='hmmscan', run=terminal.Run()):
        self.alphabet = alphabet
        self.context = context
        self.program = program

        self.run = run

        if self.context == "GENE":
            if self.alphabet == 'AA':
                col_info = self.get_col_info_for_GENE_AA_context()
            else:
                # this is a relatively recent addition to accommodate for GENE context using
                # the DNA alphabet for HMM searches. for details, please see the issue
                # https://github.com/merenlab/anvio/issues/2089. special thanks go to
                # Gio Kanaan for providing a test dataset to address the problem.
                col_info = self.get_col_info_for_GENE_DNA_context()
        elif self.context == "CONTIG" and (self.alphabet == "DNA" or self.alphabet == "RNA"):
            col_info = self.get_col_info_for_CONTIG_context()
        elif self.context == "DOMAIN" and self.alphabet == "AA":
            if program != 'hmmsearch':
                raise ConfigError("HMMScan :: the 'DOMAIN' context is only available for hmmsearch.")
            col_info = self.get_col_info_for_DOMAIN_context()
        else:
            raise ConfigError("HMMScan driver is confused. Yor context and alphabet pair ('%s' and '%s') "
                              "does not seem to be implemented in the parser module. If you think this is "
                              "not a mistake on your part, please get in touch with the anvi'o developers "
                              "and watch them fix it like actual pros." % (self.context, self.alphabet))

        col_names, col_mapping = col_info

        files_structure = {
            'hits': {
                'col_names': col_names,
                'col_mapping': col_mapping,
                'indexing_field': -1,
                'no_header': True,
            },
        }

        ## Here we have some sad extra parsing that will hopefully go away once HMMER4 comes out. Due to stupid
        ## space-delimited and column-aligned HMMER output with a description field that can contain internal spaces,
        ## when anvio.driver.HMMER converts spaces into to tabs we can end up with lines of variable columns. Which
        ## means that before we can send this file to the base parser, we have to combine the split description fields
        ## into one column. Yeah, it sucks doing it this way. But we tried asking the HMMER people to just give us tab-
        ## delimited output already (https://github.com/EddyRivasLab/hmmer/issues/235) and it seems that we'll have to
        ## wait for HMMER4 to have this feature. So it has to be this way until HMMER 4 comes out.
        ##
        ## We could just parse everything with pandas right here, but the base parser is already set up nicely to do it,
        ## and it will be easy to switch back to using that parser directly once the promised tab-delimited output is
        ## implemented. So for now we just make a slight detour to fix the shitty format before going back to the base
        ## parser.
        fixed_hmmer_table_txt =  self.fix_sad_hmmer_table_output(hmmer_table_txt, col_names)
        files_expected = {'hits': fixed_hmmer_table_txt}

        Parser.__init__(self, self.program, [fixed_hmmer_table_txt], files_expected, files_structure)


    def get_col_info_for_GENE_AA_context(self):
        """Get column names and types for GENE context for AA alphabet

        See class docstring for details of the fields for AA sequence search, and DNA sequence search.
        """

        if self.program == 'hmmscan':
            #                                                               |-- full sequence ---| |-- best 1 domain ---| |-- domain number estimation ---|
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
            #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
            col_info = [
                ('gene_name', str),         # target name
                ('gene_hmm_id', str),       # accession
                ('gene_callers_id', int),   # query name
                ('gene_call_acc', str),     # accession
                ('e_value', float),         # e-value (full sequence)
                ('bit_score', float),       # score (full sequence)
                ('bias', str),              # bias (full sequence)
                ('dom_e_value', float),       # e-value (best 1 domain)
                ('dom_bit_score', float),   # score (best 1 domain)
                ('dom_bias', str),          # bias (best 1 domain)
                ('dom_exp', str),           # exp (domain number estimation)
                ('dom_reg', str),           # reg (domain number estimation)
                ('dom_clu', str),           # clu (domain number estimation)
                ('dom_ov', str),            # ov (domain number estimation)
                ('dom_env', str),           # env (domain number estimation)
                ('dom_dom', str),           # dom (domain number estimation)
                ('dom_rep', str),           # rep (domain number estimation)
                ('dom_inc', str),           # inc (domain number estimation)
                ('description', str)        # description
            ]

        elif self.program == 'hmmsearch':
            #                                                               |-- full sequence ---| |-- best 1 domain ---| |-- domain number estimation ---|
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
            #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
            col_info = [
                ('gene_callers_id', int),   # target name
                ('gene_call_acc', str),     # accession
                ('gene_name', str),         # query name
                ('gene_hmm_id', str),       # accession
                ('e_value', float),         # e-value (full sequence)
                ('bit_score', float),       # score (full sequence)
                ('bias', str),              # bias (full sequence)
                ('dom_e_value', float),       # e-value (best 1 domain)
                ('dom_bit_score', float),   # score (best 1 domain)
                ('dom_bias', str),          # bias (best 1 domain)
                ('dom_exp', str),           # exp (domain number estimation)
                ('dom_reg', str),           # reg (domain number estimation)
                ('dom_clu', str),           # clu (domain number estimation)
                ('dom_ov', str),            # ov (domain number estimation)
                ('dom_env', str),           # env (domain number estimation)
                ('dom_dom', str),           # dom (domain number estimation)
                ('dom_rep', str),           # rep (domain number estimation)
                ('dom_inc', str),           # inc (domain number estimation)
                ('description', str)        # description
            ]

        else:
            raise ConfigError("The HMMScan Parser class is not sure if you know what you are doing. You told it that you wanted to "
                                "parse HMM hits from the program %s, but this class doesn't know how to handle those." % (self.program))

        return list(zip(*col_info))


    def get_col_info_for_GENE_DNA_context(self):
        """Get column names and types for GENE context for DNA alphabet

        See class docstring for details of the fields for AA sequence search, and DNA sequence search.
        """

        # target name        accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  modlen strand   E-value  score  bias  description of target
        #------------------- ---------- -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
        col_info = [
            ('gene_name', str),         # target name
            ('gene_hmm_id', str),       # accession
            ('gene_callers_id', int),   # query name
            ('gene_call_acc', str),     # accession
            ('hmm_from', int),          # hmmfrom
            ('hmm_to', int),            # hmm to
            ('alignment_from', int),    # alifrom
            ('alignment_to', int),      # ali to
            ('envelope_from', int),     # envfrom
            ('envelope_to', int),       # env to
            ('seq_len', str),           # modlen
            ('strand', str),            # strand
            ('e_value', float),         # E-value
            ('score', str),             # score
            ('bias', str),              # bias
            ('description', str)        # description of target
        ]

        return list(zip(*col_info))


    def get_col_info_for_CONTIG_context(self):
        """Get column names and types for CONTIG context

        See class docstring for details of the fields for AA sequence search, and DNA sequence search.
        """

        # target name        accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  modlen strand   E-value  score  bias  description of target
        #------------------- ---------- -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
        col_info = [
            ('gene_name', str),         # target name
            ('gene_hmm_id', str),       # accession
            ('contig_name', str),       # query name
            ('contig_acc', str),        # accession
            ('hmm_from', int),          # hmmfrom
            ('hmm_to', int),            # hmm to
            ('alignment_from', int),    # alifrom
            ('alignment_to', int),      # ali to
            ('envelope_from', int),     # envfrom
            ('envelope_to', int),       # env to
            ('seq_len', str),           # modlen
            ('strand', str),            # strand
            ('e_value', float),         # E-value
            ('score', str),             # score
            ('bias', str),              # bias
            ('description', str)        # description of target
        ]

        return list(zip(*col_info))


    def get_col_info_for_DOMAIN_context(self):
        """Get column names and types for DOMAIN context

        See class docstring for details of the fields
        """
        #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
        # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
        #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
        col_info = [
            ('gene_callers_id', int),   # target name
            ('gene_call_acc',   str),   # accession
            ('gene_length',     int),   # tlen
            ('gene_name',       str),   # query name
            ('gene_hmm_id',     str),   # accession
            ('hmm_length',      int),   # qlen
            ('e_value',         float), # E-value (full sequence)
            ('bit_score',       float), # score (full sequence)
            ('bias',            float), # bias (full sequence)
            ('match_num',       int),   # # (this domain)
            ('num_matches',     int),   # of (this domain)
            ('dom_c_evalue',    float), # c-Evalue (this domain)
            ('dom_i_evalue',    float), # i-Evalue (this domain)
            ('dom_bit_score',   str),   # score (this domain)
            ('dom_bias',        float), # bias (this domain)
            ('hmm_start',       int),   # from (hmm coord)
            ('hmm_stop',        int),   # to (hmm coord)
            ('gene_start',      int),   # from (ali coord)
            ('gene_stop',       int),   # to (ali coord)
            ('env_start',       str),   # from (env coord)
            ('env_stop',        str),   # to (env coord)
            ('mean_post_prob',  float), # acc
            ('description',     str)    # description of target
        ]

        return list(zip(*col_info))


    def fix_sad_hmmer_table_output(self, hmmer_table_txt, column_names):
        """"Fixes the problem of inconsistent column numbers by combining multiple description-related columns into one.

        This is sadly necessary because the original description column in space-delimited HMMER output can internally
        contain spaces and will be split up into a variable number of columns across different rows. Well. It is not
        _strictly_ necessary (since one alternative is to just drop the description column entirely). But that is what
        we are doing right now.

        Parameters
        ==========
        hmmer_table_txt, str
            Path to the tabular output file generated by the anvio.driver.HMMER.run_hmmer() function.

        column_names, tuple of str
            A list of column names in the hmmer_table_txt file, which depends on the context.

        Returns
        =======
        fixed_table_path, str
            Path to a fixed version of hmmer_table_txt that can be read using utils.get_TAB_delimited_file_as_dictionary()
        """

        tmp_dir = os.path.dirname(hmmer_table_txt)
        fixed_table_path = os.path.join(tmp_dir, "hmm.table.fixed")

        hmmer_df = pd.DataFrame([line.strip().split() for line in open(hmmer_table_txt, 'r')])

        # we tack on extra columns to catch up all words split from the original description column
        num_cols = len(hmmer_df.columns)
        num_extra_cols = num_cols - len(column_names)
        extra_desc_cols = [str(x) for x in range(num_extra_cols)]
        col_names_plus_description_cols = list(column_names) + extra_desc_cols
        all_desc_cols = ['description'] + extra_desc_cols

        hmmer_df.columns = col_names_plus_description_cols
        hmmer_df.fillna('', inplace=True)

        # join the description columns into one
        hmmer_df['full_description'] = hmmer_df[all_desc_cols[0]]
        for c in all_desc_cols[1:]:
            # Concatenate two string columns
            hmmer_df['full_description'] += ' ' + hmmer_df[c]
        hmmer_df['full_description'] = hmmer_df['full_description'].str.strip()
        hmmer_df.drop(all_desc_cols, axis=1, inplace=True)
        hmmer_df.rename(columns={'full_description': 'description'}, inplace=True)

        hmmer_df.to_csv(fixed_table_path, sep="\t", index=False, header=False)

        return fixed_table_path


    def get_search_results(self, noise_cutoff_dict=None):
        """Goes through the hits provided by `hmmscan` and generates an annotation dictionary with the relevant information about each hit.

        This function makes sure only hits with a high enough bit score make it into the annotation dictionary.

        Parameters
        ==========
        noise_cutoff_dict : dictionary
            dictionary of noise cutoff terms; see setup_ko_dict in kofam.py for an example

        Returns
        =======
        annotations_dict : dictionary
            dictionary of annotations, one annotation per HMM hit. Includes full and domain-level bitscore
            for GENE context.
        """

        annotations_dict = {}

        # this is the stuff we are going to try to fill with this:
        # search_table_structure = ['entry_id', 'source', 'alphabet', 'contig', 'gene_callers_id' 'gene_name', 'gene_hmm_id', 'e_value']
        # plus we also add bitscore values for GENE context

        entry_id = 0
        num_hits_removed = 0 # a counter for the number of hits we don't add to the annotation dictionary

        for hit in list(self.dicts['hits'].values()):
            entry = None

            if self.context == 'GENE' or self.context == 'DOMAIN':
                # Here we only add the hit to the annotations_dict if the appropriate bit score is above the
                # threshold set in noise_cutoff_dict (which is indexed by profile name (aka gene_name in the hits dict)
                if noise_cutoff_dict and hit['gene_name'] in noise_cutoff_dict.keys():
                    hmm_entry_name =  hit['gene_name']
                    score_type = noise_cutoff_dict[hmm_entry_name]['score_type']
                    threshold = noise_cutoff_dict[hmm_entry_name]['threshold']
                    keep = True
                    if score_type == 'full':
                        if hit['bit_score'] < float(threshold):
                            keep = False
                    elif score_type == 'domain':
                        if hit['dom_bit_score'] < float(threshold):
                            keep = False
                    else:
                        raise ConfigError(f"Oh dear. The HMM profile {hit['gene_name']} has a strange score_type value: "
                                          f"'{score_type}'. The only accepted values for this type are 'full' or 'domain', "
                                          "so anvi'o cannot parse the hits to this profile. Please contact a developer for "
                                          "help.")

                    if keep:
                        entry = {'entry_id': entry_id,
                                 'gene_name': hit['gene_name'],
                                 'gene_hmm_id': hit['gene_hmm_id'],
                                 'gene_callers_id': hit['gene_callers_id'],
                                 'e_value': hit['e_value'],
                                 'bit_score': hit['bit_score'],
                                 'domain_bit_score': hit['dom_bit_score'],
                                 'domain_e_value': hit['dom_e_value'] if 'dom_e_value' in hit else None}
                    else:
                        num_hits_removed += 1

                elif noise_cutoff_dict and hit['gene_name'] not in noise_cutoff_dict.keys():
                    # this should never happen, in an ideal world where everything is filled with butterflies and happiness
                    raise ConfigError(f"Hmm. While parsing your HMM hits, it seems the HMM profile {hit['gene_name']} was not "
                                     "found in the noise cutoff dictionary. This should probably not ever happen, and you should "
                                     "contact a developer as soon as possible to figure out what is going on. But for now, anvi'o "
                                     "is going to fail in order to avoid adding many garbage hits to your database.")

                else:
                    if self.alphabet == 'DNA':
                        entry = {'entry_id': entry_id,
                                 'gene_name': hit['gene_name'],
                                 'gene_hmm_id': hit['gene_hmm_id'],
                                 'gene_callers_id': hit['gene_callers_id'],
                                 'start': hit['alignment_from'],
                                 'stop': hit['alignment_to'],
                                 'e_value': hit['e_value']}
                    else:
                        entry = {'entry_id': entry_id,
                                 'gene_name': hit['gene_name'],
                                 'gene_hmm_id': hit['gene_hmm_id'],
                                 'gene_callers_id': hit['gene_callers_id'],
                                 'e_value': hit['e_value'],
                                 'bit_score': hit['bit_score'],
                                 'domain_bit_score': hit['dom_bit_score'],
                                 'domain_e_value': hit['dom_e_value'] if 'dom_e_value' in hit else None}

            elif self.context == 'CONTIG' and (self.alphabet == 'DNA' or self.alphabet == 'RNA'):
                entry = {'entry_id': entry_id,
                         'gene_name': hit['gene_name'],
                         'gene_hmm_id': hit['gene_hmm_id'],
                         'contig_name': hit['contig_name'],
                         'start': hit['alignment_from'],
                         'stop': hit['alignment_to'],
                         'e_value': hit['e_value']}
            else:
                raise ConfigError("Anvi'o does not know how to parse %s:%s" % (self.alphabet, self.context))

            if entry:
                entry_id += 1
                annotations_dict[entry_id] = entry

        self.run.info("Number of weak hits removed by HMMER parser", num_hits_removed)
        self.run.info("Number of hits in annotation dict ", len(annotations_dict.keys()))


        return annotations_dict
