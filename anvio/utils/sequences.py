import os
import copy

import numpy as np

from numba import jit

import anvio
import anvio.constants as constants

from anvio.sequence import Composition
from anvio.errors import ConfigError


def rev_comp(seq):
    return seq.translate(constants.complements)[::-1]


def rev_comp_gene_calls_dict(gene_calls_dict, contig_sequence):
    contig_length = len(contig_sequence)
    gene_caller_ids = list(gene_calls_dict.keys())

    gene_caller_id_conversion_dict = dict([(gene_caller_ids[-i - 1], i) for i in range(0, len(gene_caller_ids))])
    G = lambda g: gene_caller_id_conversion_dict[g]

    reverse_complemented_gene_calls = {}
    for gene_callers_id in gene_calls_dict:
        g = copy.deepcopy(gene_calls_dict[gene_callers_id])
        g['start'], g['stop'] = contig_length - g['stop'], contig_length - g['start']
        g['direction'] = 'f' if g['direction'] == 'r' else 'r'

        reverse_complemented_gene_calls[G(gene_callers_id)] = g

    return reverse_complemented_gene_calls, gene_caller_id_conversion_dict


def get_GC_content_for_sequence(sequence):
    return Composition(sequence).GC_content



def get_synonymous_and_non_synonymous_potential(list_of_codons_in_gene, just_do_it=False):
    """
    When calculating pN/pS or dN/dS, the number of variants classified as synonymous or non
    synonymous need to be normalized by the sequence's potential for synonymous and
    non-synonymous changes. That is calculated by mutating each position to the other 3
    nucleotides and calculating whether the mutation is synonymous or non synonymous. Each
    mutation gets a score of 1/3, since there are 3 possible mutations for each site. If the
    sequence is of length L, the nonsynonymous and synonymous potentials sum to L.

    list_of_codons_in_gene is a list of the codons as they appear in the gene sequence, e.g.
    ['ATG', ..., 'TAG'], which can be generated from utils.get_list_of_codons_for_gene_call
    """
    if not any([list_of_codons_in_gene[-1] == x for x in ['TAG', 'TAA', 'TGA']]) and not just_do_it:
        raise ConfigError("The sequence `get_synonymous_and_non_synonymous_potential` received does "
                          "end with a stop codon and may be irrelevant for this analysis. If you "
                          "want to continue anyways, include the flag `--just-do-it` in your call "
                          "(if you are a programmer see the function header).")

    synonymous_potential = 0
    num_ambiguous_codons = 0 # these are codons with Ns or other characters than ATCG

    for codon in list_of_codons_in_gene:
        # first test if it is proper codon
        if not codon:
            num_ambiguous_codons += 1
            continue

        # if we are here, this is a proper codon
        for i, nt in enumerate(codon):
            for mutant_nt in [m for m in 'ACGT' if m != nt]:

                mutant_codon = list(codon)
                mutant_codon[i] = mutant_nt
                mutant_codon = ''.join(mutant_codon)

                if constants.codon_to_AA[mutant_codon] == constants.codon_to_AA[codon]:
                    synonymous_potential += 1/3

    non_synonymous_potential = 3 * (len(list_of_codons_in_gene) - num_ambiguous_codons) - synonymous_potential

    return synonymous_potential, non_synonymous_potential, num_ambiguous_codons



def get_list_of_AAs_for_gene_call(gene_call, contig_sequences_dict):

    list_of_codons = get_list_of_codons_for_gene_call(gene_call, contig_sequences_dict)
    list_of_AAs = []

    for codon in list_of_codons:

        # if concensus sequence contains shitty characters, we will not continue
        if codon not in constants.codon_to_AA:
            continue

        # genes in the reverse direction are already handled in get_list_of_codons_for_gene_call so
        # all we do is transform codons to AAs
        list_of_AAs.append(constants.codon_to_AA[codon])

    return list_of_AAs



def get_list_of_codons_for_gene_call(gene_call, contig_sequences_dict, **kwargs):
    """Get a list of the codons for a gene call

    Parameters
    ==========
    contig_sequences_dict : dict
        An object that looks like that ContigsSuperclass.contig_sequences (initialized with
        ContigsSuperclass.init_contig_sequences)
    """

    codon_order_to_nt_positions = get_codon_order_to_nt_positions_dict(gene_call, **kwargs)

    if gene_call['contig'] not in contig_sequences_dict:
        raise ConfigError("get_list_of_AAs_for_gene_call: The contig sequences dict sent to "
                           "this function does contain the contig name that appears in the gene call. "
                           "Something is wrong here...")

    try:
        contig_sequence = contig_sequences_dict[gene_call['contig']]['sequence']
    except:
        raise ConfigError("get_list_of_AAs_for_gene_call: The contig sequences dict sent to "
                           "this function does not seem to be an anvi'o contig sequences dict :/ It "
                           "doesn't have the item 'sequence' in it.")

    list_of_codons = []
    for codon_order in codon_order_to_nt_positions:
        nt_positions = codon_order_to_nt_positions[codon_order]

        # here we cut it from the contig sequence
        reference_codon_sequence = contig_sequence[nt_positions[0]:nt_positions[2] + 1]

        # NOTE: here we make sure the codon sequence is composed of unambiguous nucleotides.
        # and we will not inlcude those that contain anything other than proper
        # nucleotides in the resulting list of codons.
        if set(reference_codon_sequence).issubset(constants.unambiguous_nucleotides):
            list_of_codons.append(constants.codon_to_codon_RC[reference_codon_sequence] if gene_call['direction'] == 'r' else reference_codon_sequence)
        else:
            list_of_codons.append(None)

    return list_of_codons



def get_translated_sequence_for_gene_call(sequence, gene_callers_id, return_with_stops=False):
    try:
        translated_sequence = translate(sequence)
    except ConfigError:
        raise ConfigError("The sequence corresponding to the gene callers id '%s' has %d nucleotides, "
                          "which is indivisible by 3. This is bad because it is now ambiguous which codon "
                          "frame should be used for translation into an amino acid sequence. Here is "
                          "the culprit sequence: %s" % (gene_callers_id, len(sequence), sequence))

    if translated_sequence.endswith('*'):
        if return_with_stops:
            pass
        else:
            translated_sequence = translated_sequence[:-1]

    return translated_sequence



def translate(sequence):
    """Translate a sequence. As stupid as possible.

    Returns
    =======
    amino_acid_sequence : str
        Amino acid sequence of sequence. If translation of codon is unknown, X is used. All stop
        codons are included and represented as *.

    Notes
    =====
    - Raises error if indivisible by 3
    - Consider smarter functions: utils.get_translated_sequence_for_gene_call,
      utils.get_most_likely_translation_frame
    """

    N = len(sequence)
    sequence = sequence.upper()
    translated_sequence = []

    if N % 3:
        raise ConfigError("utils.translate :: sequence is not divisible by 3: %s" % sequence)

    for i in range(0, N, 3):
        aa = constants.AA_to_single_letter_code[constants.codon_to_AA[sequence[i:i + 3]]] or 'X'
        translated_sequence.append(aa)

    return ''.join(translated_sequence)



def get_most_likely_translation_frame(sequence, model=None, null_prob=None, stop_prob=None, log_likelihood_cutoff=2):
    """Predict the translation frame with a markov model of amino acid sequences

    Parameters
    ==========
    sequence : str
        A DNA sequence

    model : numpy array, None
        A numpy array of transition probabilities. For an example, see
        anvio/data/seq_transition_models/AA/3rd_order.npy

    null_prob : float, None
        When a markov state contains an unspecified amino acid (X), what probability transition should
        be applied to it? To be as uninformative as possible, if None, null_prob is set to the median
        transition probability of the model.

    stop_prob : float, None
        When a markov state contains a stop codon, what transition probability should
        be applied to it? Since internal stop codons are exceedingly rare, if None, stop_prob is set
        to be 1/1,000,000th the probability of the minimum transition probability of the model.

    log_likelihood_cutoff : float, 2
        If the best frame has a log likelihood with respect to the second best frame that is less
        than this value, the frame is set to None, which is to say, anvi'o is not confident enough
        to tell you the frame. The amino acid sequence is still returned. The default is 2, which
        means the probability of the first should be at least 10^2 times the probability of the
        competing frame

    Returns
    =======
    frame, amino_acid_sequence : int, str
        frame is the number of shifted nucleotides that produced the most likely frame and is either
        0, 1, or 2. amino_acid_sequence is the translated sequence. If less than log_likelihood_cutoff,
        None is returned as the frame
    """

    N = len(sequence)
    if N == 3:
         # Save ourselves the effort
        return 0, translate(sequence)
    elif N < 3:
        raise ConfigError("utils.get_most_likely_translation_frame :: sequence has a length less than 3 "
                          "so there is nothing to translate.")

    if model is None:
        default_model_path = os.path.join(os.path.dirname(anvio.__file__), 'data/seq_transition_models/AA/fourth_order.npy')
        model = np.load(default_model_path)

    order = len(model.shape)
    null_prob = null_prob if null_prob is not None else np.median(model)
    stop_prob = stop_prob if stop_prob is not None else model.min()/1e6

    aas = [constants.AA_to_single_letter_code[aa] for aa in constants.amino_acids if aa != 'STP']
    aa_to_array_index = {aa: i for i, aa in enumerate(aas)}

    # Collect all of the candidate sequences

    candidates = {}
    for frame in range(3):
        frame_seq = sequence[frame:]

        remainder = len(frame_seq) % 3
        if remainder:
            frame_seq = frame_seq[:-remainder]

        if not frame_seq:
            continue

        candidates[frame] = {
            'sequence': translate(frame_seq),
            'log_prob': 0,
        }

    # Calculate the log probability of each candidate

    smallest_seq_length = min([len(candidate['sequence']) for candidate in candidates.values()])

    for frame in candidates:
        # Some of the candidates will be one AA smaller. To not skew values, we truncate each
        # candidate to the length of the smallest candidate
        seq = candidates[frame]['sequence'][:smallest_seq_length]

        trans_probs = np.zeros(smallest_seq_length - order)

        for codon_order in range(smallest_seq_length - order):
            state = seq[codon_order:codon_order+order]

            if '*' in state:
                trans_probs[codon_order] = stop_prob
            elif 'X' in state:
                trans_probs[codon_order] = null_prob
            else:
                state_as_indices = tuple([aa_to_array_index[aa] for aa in state])
                trans_probs[codon_order] = model[state_as_indices]

        candidates[frame]['log_prob'] = np.sum(np.log10(trans_probs))

    frame_second, frame_best = sorted(candidates, key=lambda frame: candidates[frame]['log_prob'])[-2:]
    log_prob_best = candidates[frame_best]['log_prob']
    log_prob_second = candidates[frame_second]['log_prob']

    if (log_prob_best - log_prob_second) < log_likelihood_cutoff:
        # Frame is not league's better than the competing frame, which it should be if we are to
        # have any confidence in it. The sequence is returned
        return None, candidates[frame_best]['sequence']

    amino_acid_sequence = candidates[frame_best]['sequence']

    # if the best amino acid sequence ends with a stop codon, remove it.
    amino_acid_sequence = amino_acid_sequence[:-1] if amino_acid_sequence.endswith('*') else amino_acid_sequence

    return frame_best, amino_acid_sequence



def get_codon_order_to_nt_positions_dict(gene_call, subtract_by=0):
    """Returns a dictionary to translate codons in a gene to nucleotide positions

    Parameters
    ==========
    subtract_by : int, 0
        Subtract the start and stop of the gene call by this amount. This could be useful if the
        gene call start/stop are defined in terms of the contig, but you want the start/stop in
        terms of the split. Then you could supply subtract_by=split_start, where split_start is the
        start of the split
    """

    if gene_call['call_type'] != constants.gene_call_types['CODING']:
        raise ConfigError("utils.get_codon_order_to_nt_positions_dict :: this simply will not work "
                           "for noncoding gene calls, and gene caller id %d is noncoding." % gene_call['gene_callers_id'])

    start = gene_call['start'] - subtract_by
    stop = gene_call['stop'] - subtract_by

    codon_order_to_nt_positions = {}
    codon_order = 0

    if gene_call['direction'] == 'r':
        for nt_pos in range(stop - 1, start - 1, -3):
            codon_order_to_nt_positions[codon_order] = [nt_pos - 2, nt_pos - 1, nt_pos]
            codon_order += 1
    else:
        for nt_pos in range(start, stop, 3):
            codon_order_to_nt_positions[codon_order] = [nt_pos, nt_pos + 1, nt_pos + 2]
            codon_order += 1

    return codon_order_to_nt_positions



def nt_seq_to_nt_num_array(seq, is_ord=False):
    """Convert a string of sequence into an array of numbers

    Performance compared to {list comprehension with dictionary lookup} depends on sequence length.
    See Examples

    Parameters
    ==========
    seq : str
        string with A, C, T, G, N as its characters, e.g. 'AATGCN'

    is_ord : bool, False
        set True if seq is already a numpy array, where each element is the ord of the sequence. E.g.
        if `seq` is passed as array([65, 65, 67]), then it is already the ordinal representation of
        'AAC'

    Returns
    =======
    output : numpy array
        E.g. if seq = 'AATGCN', output = array([0, 0, 2, 3, 1, 4])

    Examples
    ========

    Init an environment

    >>> import anvio.constants as constants
    >>> import anvio.utils as utils
    >>> seq_short = ''.join(list(np.random.choice(constants.nucleotides, size=100)))
    >>> seq_long = ''.join(list(np.random.choice(constants.nucleotides, size=100_000_000)))
    >>> nt_to_num =  {'A': 0, 'C': 1, 'T': 2, 'G': 3, 'N': 4}

    Time short sequence:

    >>> %timeit utils.nt_seq_to_nt_num_array(seq_short)
    2.36 µs ± 20.9 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    >>> %timeit [nt_to_num[s] for s in seq_short]
    5.83 µs ± 20.7 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)

    Time long sequence:

    >>> %timeit utils.nt_seq_to_nt_num_array(seq_long)
    653 ms ± 1.02 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    >>> %timeit [nt_to_num[s] for s in seq_long]
    5.27 s ± 13.4 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    """

    return constants.nt_to_num_lookup[seq if is_ord else np.frombuffer(seq.encode('ascii'), np.uint8)]



def nt_seq_to_RC_nt_num_array(seq, is_ord=False):
    """Convert a string of sequence into an array of numbers, reverse-complemented

    Performance compared to {list comprehension with dictionary lookup} depends on sequence length.
    See Examples

    Parameters
    ==========
    seq : str
        string with A, C, T, G, N as its characters, e.g. 'AATGCN'

    is_ord : bool, False
        set True if seq is already a numpy array, where each element is the ord of the sequence. E.g.
        if `seq` is passed as array([65, 65, 67]), then it is already the ordinal representation of
        'AAC'

    Returns
    =======
    output : numpy array
        E.g. if seq = 'AATGCN', output = array([4, 2, 0, 1, 3, 3])

    Examples
    ========
    See `nt_seq_to_nt_num_array` docstring for examples
    """

    return constants.nt_to_RC_num_lookup[seq if is_ord else np.frombuffer(seq.encode('ascii'), np.uint8)][::-1]



def nt_seq_to_codon_num_array(seq, is_ord=False):
    """Convert a sequence into an array of numbers corresponding to codons

    Parameters
    ==========
    seq : str
        string with A, C, T, G as its characters, e.g. 'AATGCT'. seq must be divisible by 3

    is_ord : bool, False
        set True if seq is already a numpy array, where each element is the ord of the sequence. E.g.
        if `seq` is passed as array([65, 65, 67]), then it is already the ordinal representation of
        'AAC'

    Notes
    =====
    - Delegates to just-in-time compiled function
    """

    return _nt_seq_to_codon_num_array(
        seq if is_ord else np.frombuffer(seq.encode('ascii'), np.uint8),
        constants.codon_to_num_lookup,
    )



def nt_seq_to_RC_codon_num_array(seq, is_ord=False):
    """Convert a sequence into an array of numbers corresponding to codons, reverse-complemented

    Parameters
    ==========
    seq : str
        string with A, C, T, G as its characters, e.g. 'AATGCT'. seq must be divisible by 3

    is_ord : bool, False
        set True if seq is already a numpy array, where each element is the ord of the sequence. E.g.
        if `seq` is passed as array([65, 65, 67]), then it is already the ordinal representation of
        'AAC'

    Notes
    =====
    - Delegates to just-in-time compiled function
    """

    return _nt_seq_to_codon_num_array(
        seq if is_ord else np.frombuffer(seq.encode('ascii'), np.uint8),
        constants.codon_to_RC_num_lookup,
    )[::-1]



@jit(nopython=True)
def _nt_seq_to_codon_num_array(seq_as_ascii_ints, lookup_codon):
    """Should be called through its parent functions `nt_seq_to_codon_num_array` and `nt_seq_to_RC_codon_num_array`"""

    output = np.zeros(len(seq_as_ascii_ints)//3, dtype=np.uint8)

    for i in range(0, seq_as_ascii_ints.shape[0], 3):
        output[i//3] = lookup_codon[seq_as_ascii_ints[i], seq_as_ascii_ints[i+1], seq_as_ascii_ints[i+2]]

    return output



def is_amino_acid_functionally_conserved(amino_acid_residue_1, amino_acid_residue_2):
    """Checks if two amino acid residues are part of the same biochemical property group"""
    group = constants.amino_acid_property_group[amino_acid_residue_1]
    conserved_group = constants.conserved_amino_acid_groups[group]

    if amino_acid_residue_2 in conserved_group:
        return True

    if group == 'Polar and Nonpolar':
        #they fall in more than one group, multiple tests needed
        if amino_acid_residue_1 == 'H' and (amino_acid_residue_2 in constants.conserved_amino_acid_groups['Nonpolar'] \
                                            or amino_acid_residue_2 in constants.conserved_amino_acid_groups['Bases']):
            return True

        if amino_acid_residue_1 == 'Y' and (amino_acid_residue_2 in constants.conserved_amino_acid_groups['Aromatic']):
            return True

    return False



def convert_sequence_indexing(index, source="M0", destination="M1"):
    """
    Anvi'o zero-indexes sequences. For example, the methionine that every ORF starts with has the
    index 0 (M0). This is in contrast to the most conventions, in which the methionine is indexed by
    1 (M1). This function converts between the two.

    index : integer, numpy array, pandas series, list
        The sequence index/indices you are converting.
    source : string
        The convention you are converting from. Must be either "M0" (anvio) or
        "M1" (not anvio)
    destination : string
        The convention you are converting to. Must be either "M0" (anvio) or
        "M1" (not anvio)
    """
    convert = lambda x, a: [i + a for i in x] if type(x) == list else x + a

    if source not in ["M0", "M1"] or destination not in ["M0", "M1"]:
        raise ValueError("Must be 'M0' or 'M1'.")

    if source == "M0" and destination == "M1":
        return convert(index, 1)

    if source == "M1" and destination == "M0":
        return convert(index, -1)

    return index



def summarize_alignment(sequence):
    """Takes an alignment, and returns its summary.

        >>> alignment = '----AA---TTTT-----CC-GGGGGGGG----------------ATCG--'
        >>> sequence = alignment.replace('-')
        >>> summarize_alignment(alilgnment)
        '-|4|2|3|4|5|2|1|8|16|4|2'
        >>> summary = summarize_alignment(alignment)
        >>> restore_alignment(sequence, summary)
        '----AA---TTTT-----CC-GGGGGGGG----------------ATCG--'
    """
    alignment_summary = []

    starts_with_gap = sequence[0] == '-'
    in_gap, in_nt = (True, False) if starts_with_gap else (False, True)

    gap, nt = 0, 0
    for i in range(0, len(sequence)):
        if sequence[i] == '-':
            if in_nt:
                alignment_summary.append(nt) if nt else None
                in_gap, in_nt = True, False
                nt = 0
                gap = 1
            else:
                gap += 1
        else:
            if in_gap:
                alignment_summary.append(gap) if gap else None
                in_gap, in_nt = False, True
                gap = 0
                nt = 1
            else:
                nt += 1

    alignment_summary.append(gap or nt)

    return  '|'.join(['-' if starts_with_gap else '.'] + [str(s) for s in alignment_summary])



def restore_alignment(sequence, alignment_summary, from_aa_alignment_summary_to_dna=False):
    """Restores an alignment from its sequence and alignment summary.

       See `summarize_alignment` for the `alignment_summary` compression.
    """

    if not alignment_summary:
        return sequence

    if isinstance(sequence, bytes):
        sequence = list(sequence.decode('utf-8'))
    elif isinstance(sequence, str):
        sequence = list(sequence)
    else:
        raise ConfigError("Sequence must be of type str or bytes. What you sent is of %s :/" % type(sequence))

    in_gap = alignment_summary[0] == '-'

    alignment = ''
    for part in [(int(p) * 3) if from_aa_alignment_summary_to_dna else int(p) for p in alignment_summary.split('|')[1:]]:
        if in_gap:
            alignment += '-' * part
            in_gap = False
        else:
            for i in range(0, part):
                alignment += sequence.pop(0)
            in_gap = True

    if from_aa_alignment_summary_to_dna:
        return alignment + ''.join(sequence)
    else:
        return alignment

