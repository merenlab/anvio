# -*- coding: utf-8
"""Parser for HMMER's various outputs"""

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.parsers.base import Parser

import pandas as pd


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class HMMERStandardOutput(object):
    """Parse the standard output of HMMER programs

    The main meat of this class is to produce the attributes:

        (1) self.seq_hits
        (2) self.dom_hits

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
        |               consensus_align              match_align             target_align
        | 0     vvtGggGFlGrrivkeLlrl...  +v+Gg+G++G++ v +L++ ...  LVLGGAGYIGSHAVDQLISK...
        | 1     vvtGggGFlGrrivkeLlrl...  ++ Gg+GFlG++i k L+++...  IIFGGSGFLGQQIAKILVQR...
        | ...                       ...                      ...                      ...
        | 2896  gllagllvkrvapeaagsGi...  g++  +++ r+  + a  G ...  GVVFTYFYTRF-GKNASRGN...
        | 2897  kyivvdlaeqrllvlyengk...  +yi++dl++q++ +++ +gk...  NYIEIDLKDQKM-YCFIDGK...

    If you're confused about the meaning of these columns, please see starting from page 32
    of the HMMER guide http://eddylab.org/software/hmmer/Userguide.pdf. There you will be able
    to with relative ease correlate the column names in these tables to what is described
    meticulously in the tutorial. For example, `best_dom_bias` refers to the the 'bias (best 1
    domain)' column.

    Parameters
    ==========
    hmmer_std_out : str
        Path to output of HMMER.
    """

    def __init__(self, hmmer_std_out, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        self.hmmer_std_out = hmmer_std_out

        # This is converted to a dataframe after populating
        self.seq_hits = {
            'query': [],
            'acc': [],
            'target': [],
            'query_len': [],
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
            'query': str,
            'acc': str,
            'target': str,
            'query_len': int,
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
            'query': [],
            'acc': [],
            'target': [],
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
            'consensus_align': [],
            'match_align': [],
            'target_align': [],
        }

        self.dom_hits_dtypes = {
            'query': str,
            'acc': str,
            'target': str,
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
            'consensus_align': str,
            'match_align': str,
            'target_align': str,
        }

        self.delim_query = '//\n'
        self.delim_seq = '>>'
        self.delim_domain = '=='

        self.load()


    def load(self):
        self.progress.new('Processing HMMER output')
        self.progress.update('Loading %s' % self.hmmer_std_out)

        with open(self.hmmer_std_out) as f:
            for i, query in enumerate(utils.get_chunk(f, separator=self.delim_query, read_size=32768)):

                if i % 500 == 0:
                    self.progress.update('%d done' % i)
                    self.progress.increment(increment_to=i)

                self.process_query(query)

        self.seq_hits = pd.DataFrame(self.seq_hits).astype(self.seq_hits_dtypes)
        self.dom_hits = pd.DataFrame(self.dom_hits).astype(self.dom_hits_dtypes)

        self.progress.end()
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

            self.seq_hits['query'].append(query_name)
            self.seq_hits['query_len'].append(query_len)
            self.seq_hits['acc'].append(acc)
            self.seq_hits['evalue'].append(float(seq_scores[0]))
            self.seq_hits['score'].append(float(seq_scores[1]))
            self.seq_hits['bias'].append(float(seq_scores[2]))
            self.seq_hits['best_dom_evalue'].append(float(seq_scores[3]))
            self.seq_hits['best_dom_score'].append(float(seq_scores[4]))
            self.seq_hits['best_dom_bias'].append(float(seq_scores[5]))
            self.seq_hits['expected_doms'].append(float(seq_scores[6]))
            self.seq_hits['num_doms'].append(int(seq_scores[7]))
            self.seq_hits['target'].append(seq_scores[8])

            num_doms_per_seq[seq_scores[8]] = int(seq_scores[7])

        num_seq_hits = len(seq_score_lines)

        for _ in range(num_seq_hits):
            target_name = self.find_line(lambda line: line.startswith(self.delim_seq)).split()[1]

            if num_doms_per_seq[target_name] == 0:
                continue

            self.line_no += 2
            for __ in range(num_doms_per_seq[target_name]):
                dom_score_summary = self.find_line(lambda line: True).split()

                self.dom_hits['query'].append(query_name)
                self.dom_hits['acc'].append(acc)
                self.dom_hits['target'].append(target_name)
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

                self.dom_hits['consensus_align'].append(''.join(consensus))
                self.dom_hits['match_align'].append(''.join(match))
                self.dom_hits['target_align'].append(''.join(target))


    def reformat_for_interacdome(self):
        """Reformat for use with interacdome"""

        self.dom_hits.rename(columns={
            'query': 'pfam_name',
            'acc': 'pfam_id',
            'target': 'corresponding_gene_call',
        }, inplace=True)

        self.seq_hits.rename(columns={
            'query': 'pfam_name',
            'acc': 'pfam_id',
            'query_len': 'pfam_len',
            'target': 'corresponding_gene_call',
        }, inplace=True)

        self.seq_hits['corresponding_gene_call'] = self.seq_hits['corresponding_gene_call'].astype(int)
        self.dom_hits['corresponding_gene_call'] = self.dom_hits['corresponding_gene_call'].astype(int)

        self.seq_hits[['pfam_id', 'version']] = self.seq_hits['pfam_id'].str.split('.', n=1, expand=True)
        self.dom_hits[['pfam_id', 'version']] = self.dom_hits['pfam_id'].str.split('.', n=1, expand=True)


class HMMERTableOutput(Parser):
    """Parse --tblout or --domtblout output formats for hmmer programs

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME NOTE FIXME
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    <rant>
    Parsing of HMMER tabular output needs to be redesigned. This code does not actually take output from hmmer
    and parse it. It parses the output file of anvio.driver.HMMER.hmmscan_worker which preprocesses the
    output format. The responsibility of HMMER output parsing needs to be consolidated in one spot. Biopython, a
    dependency of anvi'o, has an HMMER parser. See https://biopython.org/DIST/docs/api/Bio.SearchIO.HmmerIO-module.html.
    Perhaps this is more robust solution. This design is currently hanging on by a thread.
    </rant>

    Output specifictions of HMMER can be found in the user guide. At time of writing this,
    http://eddylab.org/software/hmmer/Userguide.pdf hosts the user guide.

    Parameters
    ==========
    hmmer_table_txt: ???
        Undocumented FIXME

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
            Undocumented FIXME
    """

    def __init__(self, hmmer_table_txt, alphabet='AA', context='GENE', program='hmmscan', run=terminal.Run()):
        self.alphabet = alphabet
        self.context = context
        self.program = program

        self.run = run

        files_expected = {'hits': hmmer_table_txt}

        if self.context == "GENE":
            col_info = self.get_col_info_for_GENE_context()
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

        Parser.__init__(self, 'HMMScan', [hmmer_table_txt], files_expected, files_structure)


    def get_col_info_for_GENE_context(self):
        """Get column names and types for GENE context

        See class docstring for details of the fields for AA sequence search, and DNA sequence search.
        """

        if self.program == 'hmmscan':
            #                                                               |-- full sequence ---| |-- best 1 domain ---| |-- domain number estimation ---|
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc
            #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- ---
            col_names = ['gene_name', 'gene_hmm_id', 'gene_callers_id', 'f', 'e_value', 'bit_score', 'f', 'f', 'dom_bit_score', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
            col_mapping = [str, str, int, str, float, float, str, str, float, str, str, str, str, str, str, str, str, str]
        elif self.program == 'hmmsearch':
            #                                                               |-- full sequence ---| |-- best 1 domain ---| |-- domain number estimation ---|
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc
            #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- ---
            col_names = ['gene_callers_id', 'f', 'gene_name', 'gene_hmm_id', 'e_value', 'bit_score', 'f', 'f', 'dom_bit_score', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
            col_mapping = [int, str, str, str, float, float, str, str, float, str, str, str, str, str, str, str, str, str]
        else:
            raise ConfigError("The HMMScan Parser class is not sure if you know what you are doing. You told it that you wanted to "
                                "parse HMM hits from the program %s, but this class doesn't know how to handle those." % (self.program))

        return col_names, col_mapping


    def get_col_info_for_CONTIG_context(self):
        """Get column names and types for GENE context

        See class docstring for details of the fields for AA sequence search, and DNA sequence search.
        """

        # 'hmm_target', 'hmm_acc', 'query_id', 'query_acc', 'hmm_from', 'hmm_to', 'alignment_from', 'alignment_to', 'envelope_from', 'envelope_to', 'seq_len', 'strand', 'e_value', 'score', 'bias',]
        col_names = ['gene_name', 'gene_hmm_id', 'contig_name', 'f', 'hmm_from', 'hmm_to', 'alignment_from', 'alignment_to', 'envelope_from', 'envelope_to', 'f', 'f', 'e_value', 'f', 'f']
        col_mapping = [str, str, str, str, str, str, int, int, int, int, str, str, float, str, str]

        return col_names, col_mapping


    def get_col_info_for_DOMAIN_context(self):
        """Get column names and types for DOMAIN context

        See class docstring for details of the fields
        """

        col_info = [
            ('gene_callers_id', str), # target name
            ('f',               str), # accession
            ('gene_length',     str), # tlen
            ('hmm_name',        str), # query name
            ('hmm_id',          str), # accession
            ('hmm_length',      str), # qlen
            ('evalue',          str), # E-value (full sequence)
            ('bitscore',        str), # score (full sequence)
            ('bias',            str), # bias (full sequence)
            ('match_num',       str), # # (this domain)
            ('num_matches',     str), # of (this domain)
            ('dom_c_evalue',    str), # c-Evalue (this domain)
            ('dom_i_evalue',    str), # i-Evalue (this domain)
            ('dom_bitscore',    str), # score (this domain)
            ('dom_bias',        str), # bias (this domain)
            ('hmm_start',       str), # from (hmm coord)
            ('hmm_stop',        str), # to (hmm coord)
            ('gene_start',      str), # from (ali coord)
            ('gene_stop',       str), # to (ali coord)
            ('f',               str), # from (env coord)
            ('f',               str), # to (env coord)
            ('mean_post_prob',  str), # acc
        ]

        return list(zip(*col_info))


    def get_search_results(self, noise_cutoff_dict = None):
        """Goes through the hits provided by `hmmscan` and generates an annotation dictionary with the relevant information about each hit.

        This function makes sure only hits with a high enough bit score make it into the annotation dictionary.

        Parameters
        ==========
        noise_cutoff_dict : dict
            dictionary of noise cutoff terms; see setup_ko_dict in kofam.py for an example

        Returns
        =======
        annotations_dict : dict
            dictionary of annotations
        """

        annotations_dict = {}

        # this is the stuff we are going to try to fill with this:
        # search_table_structure = ['entry_id', 'source', 'alphabet', 'contig', 'gene_callers_id' 'gene_name', 'gene_hmm_id', 'e_value']

        entry_id = 0
        num_hits_removed = 0 # a counter for the number of hits we don't add to the annotation dictionary

        for hit in list(self.dicts['hits'].values()):
            entry = None
            if self.context == 'GENE':
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
                        self.run.warning("Oh dear. The HMM profile %s has a strange score_type value: %s. The only accepted values "
                                         "for this type are 'full' or 'domain', so anvi'o cannot parse the hits to this profile. All hits "
                                         "will be kept regardless of bit score. You have been warned." % (hit['gene_name'], score_type))

                    if keep:
                        entry = {'entry_id': entry_id,
                                 'gene_name': hit['gene_name'],
                                 'gene_hmm_id': hit['gene_hmm_id'],
                                 'gene_callers_id': hit['gene_callers_id'],
                                 'e_value': hit['e_value']}
                    else:
                        num_hits_removed += 1

                elif noise_cutoff_dict and hit['gene_name'] not in noise_cutoff_dict.keys():
                    # this should never happen, in an ideal world where everything is filled with butterflies and happiness
                    self.run.warning("Hmm. While parsing your HMM hits, it seems the HMM profile %s was not found in the noise cutoff dictionary. "
                                     "This should probably not ever happen, and you should contact a developer as soon as possible to figure out what "
                                     "is going on. But for now, anvi'o is going to keep all hits to this profile. Consider those hits with a grain of salt, "
                                     "as not all of them may be good." % hit['gene_name'])
                    entry = {'entry_id': entry_id,
                             'gene_name': hit['gene_name'],
                             'gene_hmm_id': hit['gene_hmm_id'],
                             'gene_callers_id': hit['gene_callers_id'],
                             'e_value': hit['e_value']}

                else:
                    entry = {'entry_id': entry_id,
                             'gene_name': hit['gene_name'],
                             'gene_hmm_id': hit['gene_hmm_id'],
                             'gene_callers_id': hit['gene_callers_id'],
                             'e_value': hit['e_value']}
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

        self.run.info("Number of weak hits removed", num_hits_removed)
        self.run.info("Number of hits in annotation dict ", len(annotations_dict.keys()))

        return annotations_dict

