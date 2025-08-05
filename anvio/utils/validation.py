import anvio.constants as constants

from anvio.errors import ConfigError


def check_sample_id(sample_id):
    if sample_id:
        if sample_id[0] in constants.digits:
            raise ConfigError("The sample name ('%s') is not a valid one. Sample names can't start with digits. "
                              "Long story. Please specify a sample name that starts with an ASCII letter (if "
                              "there are no parameters available to you to set the sample name, it may be the "
                              "case that sample name is determined automatically from the input files you have "
                              "provided to whatever anvi'o workflow you were using, in which case you may need "
                              "to change your input file names or something :/)." % sample_id)

        allowed_chars_for_samples = constants.allowed_chars.replace('-', '').replace('.', '')
        if len([c for c in sample_id if c not in allowed_chars_for_samples]):
            raise ConfigError("The sample name ('%s') contains characters anvi'o does not like. Please "
                              "limit the characters that make up the project name to ASCII letters, "
                              "digits, and the underscore character ('_')." % sample_id)



def check_collection_name(collection_name):
    try:
        check_sample_id(collection_name)
    except:
        raise ConfigError('"%s" is not a proper collection name. A proper one should be a single word and not contain '
                           'ANY characters but digits, ASCII letters and underscore character(s). There should not be '
                           'any space characters, and the collection name should not start with a digit.' % collection_name)




def is_this_name_OK_for_database(variable_name, content, stringent=True, additional_chars_allowed=''):
    if not content:
        raise ConfigError("But the %s is empty? Come on :(" % variable_name)

    if content[0] in constants.digits:
        raise ConfigError("Sorry, %s can't start with a digit. Long story. Please specify a name "
                           "that starts with an ASCII letter." % variable_name)

    if stringent:
        allowed_chars = constants.allowed_chars.replace('.', '').replace('-', '')
    else:
        allowed_chars = constants.allowed_chars.replace('.', '')

    if len(additional_chars_allowed):
        allowed_chars += additional_chars_allowed

    if len([c for c in content if c not in allowed_chars]):
        raise ConfigError("Well, the %s contains characters that anvi'o does not like :/ Please limit the characters "
                           "to ASCII letters, digits, and the underscore ('_') character." % variable_name)



def check_contig_names(contig_names, dont_raise=False):
    all_characters_in_contig_names = set(''.join(contig_names))
    characters_anvio_doesnt_like = [c for c in all_characters_in_contig_names if c not in constants.allowed_chars]
    if len(characters_anvio_doesnt_like):
        if dont_raise:
            return False

        raise ConfigError("The name of at least one contig in your BAM file %s anvio does not "
                           "like (%s). Please go back to your original files and make sure that "
                           "the characters in contig names are limited to to ASCII letters, "
                           "digits. Names can also contain underscore ('_'), dash ('-') and dot ('.') "
                           "characters. anvio knows how much work this may require for you to go back and "
                           "re-generate your BAM files and is very sorry for asking you to do that, however, "
                           "it is critical for later steps in the analysis." \
                                % ("contains multiple characters" if len(characters_anvio_doesnt_like) > 1 else "contains a character",
                                   ", ".join(['"%s"' % c for c in characters_anvio_doesnt_like])))

    return True



def check_misc_data_keys_for_format(data_keys_list):
    """Ensure user-provided misc data keys are compatible with the current version of anvi'o"""

    if not data_keys_list:
        return

    # find out whether the user data contains the older implementation of stacked
    # bar data type
    obsolete_stackedbar_keys = [k for k in data_keys_list if k.find('!') > -1 and k.find(';') > -1]

    if len(obsolete_stackedbar_keys):
        key_violates_new_rule = obsolete_stackedbar_keys[0]
        main_key, data_items = key_violates_new_rule.split('!')
        new_rule_compatible_data_keys = ['%s!%s' % (main_key, d) for d in data_items.split(';')]

        raise ConfigError("Oh no :( We recently changed the description of the stacked bar data type, and your input data "
                          "file still has the older version. Here is the list of those that are violating the new format: "
                          "%s. To avoid this issue and to turn them into the new format, you could take '%s', and present "
                          "it as %d separate TAB-delimited entries that look like this: %s. Sorry!" % \
                                            (', '.join(['"%s"' % k for k in obsolete_stackedbar_keys]),
                                             key_violates_new_rule,
                                             len(new_rule_compatible_data_keys),
                                             ', '.join(['"%s"' % k for k in new_rule_compatible_data_keys])))



def is_gene_caller_id(gene_caller_id, raise_if_fail=True):
    """Test whether a given `gene_caller_id` looks like a legitimate anvi'o gene caller id"""
    try:
        assert(int(gene_caller_id) >= 0)
    except:
        if raise_if_fail:
            raise ConfigError(f"Anvi'o gene caller ids are represented by integers between 0 and infinity. "
                              f"and what you provided ('{gene_caller_id}') doesn't look like one :/")
        else:
            return False

    return True



def is_gene_sequence_clean(seq, amino_acid=False, can_end_with_stop=False, must_start_with_met=True):
    """Returns True if gene sequence is clean (amino acid or nucleotide), otherwise raises ConfigError

    Parameters
    ==========
    seq : str
        A string of amino acid or nucleotide sequence
    amino_acid : bool, False
        If True, the sequence is assumed to be an amino acid sequence
    can_end_with_stop : bool, False
        If True, the sequence can, but does not have to, end with * if amino_acid=True, or one of
        <TAG, TGA, TAA> if amino_acid=False.
    must_start_with_met : bool, True
        If True, the sequence must start with ATG if amino_acid=False or Met if amino_acid=True

    Returns
    =======
    value : bool

    Notes
    =====
    - A 'clean gene' depends on `amino_acid`. If amino_acid=True, must contain only the 20 1-letter
      codes (case insensitive) and start with M. If amino_acid=False, must contain only A,C,T,G
      (case insenstive), start with ATG, and have length divisible by 3. If can_end_with_stop=True,
      `seq` can end with a stop. If any intermediate  and in-frame stop codons are found, the gene
      is not clean
    """
    error_msg_template = "The gene sequence is not clean. Reason: %s"
    seq = seq.upper()

    start_char = 'M' if amino_acid else 'ATG'
    end_chars = ['*'] if amino_acid else ['TAG', 'TGA', 'TAA']

    permissible_chars = (set(constants.AA_to_single_letter_code.values())
                         if amino_acid
                         else set(constants.codons)) - set(end_chars)

    if not amino_acid:
        if len(seq) % 3:
            raise ConfigError(error_msg_template % "The number of nucleotides is not divisible by 3")

        new_seq = [] # list of length-3 strings
        for i in range(0, len(seq), 3):
            new_seq.append(seq[i:i+3])

        seq = new_seq

    if not seq[0] == start_char and must_start_with_met:
        raise ConfigError(error_msg_template % "Should start with methionine but instead starts with %s" % seq[0])

    for i, element in enumerate(seq[:-1]):
        if element in end_chars:
            l, r = min([i, 3]), min([len(seq[:-1])-i, 3])
            error_msg = error_msg_template % "Premature stop codon at %dth codon position (counting from 0).\
                                              Here is the position in the context of the sequence: ...%s[%s]%s..." \
                                                % (i, ''.join(seq[:-1][i-l:i]), element, ''.join(seq[:-1][i+1:i+r+1]))
            raise ConfigError(error_msg)

        if element not in permissible_chars:
            l, r = min([i, 3]), min([len(seq[:-1])-i, 3])
            error_msg = error_msg_template % "%s at %dth codon position (counting from zero) isn't a valid sequence\
                                              element. Here is the position in the context of the sequence: ...%s[%s]%s..." \
                                                % (element, i, ''.join(seq[:-1][i-l:i]), element, ''.join(seq[:-1][i+1:i+r+1]))
            raise ConfigError(error_msg)

    if seq[-1] in end_chars:
        if not can_end_with_stop:
            raise ConfigError(error_msg_template % "Sequence should not contain an explicit stop codon")
    elif seq[-1] not in permissible_chars:
        raise ConfigError(error_msg_template % "Last codon is not a valid character: %s" % seq[-1])

    return True


def is_ascii_only(text):
    """test whether 'text' is composed of ASCII characters only"""
    return all(ord(c) < 128 for c in text)

