import os
import gzip

import anvio
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

from anvio.utils.files import get_TAB_delimited_file_as_dictionary


def anvio_hmm_target_term_to_alphabet_and_context(target):
    """Alphabet and context recovery from the target term in anvi'o HMM source directories."""

    alphabet = None
    context = None
    fields = target.split(':')

    if len(fields) == 2:
        alphabet, context = fields
    elif len(fields) == 1:
        alphabet = fields[0]
    else:
        raise ConfigError("HMM stuff is upset with you. There are unexpected number of fields in the target "
                           "file.")

    if alphabet not in ['AA', 'DNA', 'RNA']:
        raise ConfigError("The alphabet in the target file (%s) isnot one of the alphabets anvi'o knows how to "
                           "work with. Here is a list for you to choose from: 'DNA', 'RNA', or 'AA'" % alphabet)

    if context not in ['GENE', 'CONTIG', None]:
        raise ConfigError("The context you defined in the target file (%s) does not make any sense to anvi'o. "
                           "It would have, if you had chosen one of these: 'GENE', 'CONTIG'." % context)

    if alphabet == 'AA' and context == 'CONTIG':
        raise ConfigError("You can't use the AA alphabet with the CONTIGS context :/ You need to set your target "
                           "again. 'AA' or 'AA:GENE' would have worked much better.")

    if not context:
        context = 'GENE'

    return alphabet, context



def get_pruned_HMM_hits_dict(hmm_hits_dict):
    """This function will identify HMM hits that are almost identical and keep only the most significant hit.

       This is an example situation where this problem occurs:

            http://i.imgur.com/2ZxDchp.png

       And this is how that context looks like after this function does its magic:

            http://i.imgur.com/cAPKR0E.png

       The data shown in the first screenshot resolves to an input dictionary like this one:

           {
                1: {'entry_id': 0, 'gene_name': 'Bacterial_23S_rRNA','contig_name': 'c_split_00001', 'start': 3175, 'stop': 267, 'e_value': 0.0},
                2: {'entry_id': 1, 'gene_name': 'Bacterial_16S_rRNA','contig_name': 'c_split_00001', 'start': 4996, 'stop': 3439, 'e_value': 0.0},
                3: {'entry_id': 2, 'gene_name': 'Archaeal_23S_rRNA', 'contig_name': 'c_split_00001', 'start': 3162, 'stop': 275, 'e_value': 0.0},
                4: {'entry_id': 3, 'gene_name': 'Archaeal_16S_rRNA', 'contig_name': 'c_split_00001', 'start': 4988, 'stop': 3441, 'e_value': 7.7e-240}
           }

       where entry 1 and entry 2 should be removed (becuse they overlap witth 3 and 4, respectively, and they are shorter).
    """

    # first create a simpler data structure where all hits in a single contig are accessible directly.
    hits_per_contig = {}
    for entry in hmm_hits_dict:
        e = hmm_hits_dict[entry]
        contig_name = e['contig_name']
        start = e['start'] if e['start'] < e['stop'] else e['stop']
        stop = e['stop'] if e['start'] < e['stop'] else e['start']
        length = stop - start

        if contig_name not in hits_per_contig:
            hits_per_contig[contig_name] = []

        hits_per_contig[contig_name].append((length, entry, start, stop), )

    # go through hits in each contig to find overlapping hits
    entry_ids_to_remove = set([])
    for hits in hits_per_contig.values():
        indices_with_matches = set([])
        for i in range(0, len(hits)):
            if i in indices_with_matches:
                # this one is already processed and is matching
                # with something else. no need to waste time
                continue

            overlapping_hits_indices = set([])
            for j in range(i + 1, len(hits)):
                alignment_start = max(hits[i][2], hits[j][2])
                alignment_end = min(hits[i][3], hits[j][3])
                shortest_of_the_two = min(hits[i][0], hits[j][0])

                if alignment_end - alignment_start > shortest_of_the_two / 2:
                    # the overlap between these two is more than the half of the lenght of the
                    # shorter one. this is done
                    overlapping_hits_indices.add(i)
                    overlapping_hits_indices.add(j)
                    indices_with_matches.add(j)

            if overlapping_hits_indices:
                # here we have a set of overlapping indices. we will ort them based on length,
                # and add the entry id of every match except the longest one into the shitkeeping
                # variable
                [entry_ids_to_remove.add(r) for r in sorted([hits[ind][1] for ind in overlapping_hits_indices], reverse=True)[1:]]


    # time to remove all the entry ids from the actual dictionary
    for entry_id in entry_ids_to_remove:
        hmm_hits_dict.pop(entry_id)

    return hmm_hits_dict



def get_HMM_sources_dictionary(source_dirs=[]):
    """An anvi'o HMM source directory importer.

       The directory must have five files:

       - genes.hmm.gz: compressed HMM for each gene.
       - genes.txt: three column file lists all gene names appear in the genes.hmm.gz, accession numbers if there
                    are any, and HMM source for those.
       - kind.txt: the kind of genes are there in this source. i.e., 'antibiotic_genes', or 'transporters'. the
                   term 'singlecopy' is a special one, and should be used with a domain term: 'singlecopy:bacteria',
                   'singlecopy:archaea', etc. Anvi'o utilizes single-copy sources to assess the completion of MAGs
                   later.
       - reference.txt: Where is it coming from?
       - target.txt: the target term. see `anvio_hmm_target_term_to_alphabet_and_context` for details.
       - noise_cutoff_terms.txt: how the noisy hits should be dealt with? see this for details: https://github.com/merenlab/anvio/issues/498

       For an example HMM source directory, take a look at an example in the codebase:

                https://github.com/meren/anvio/tree/master/anvio/data/hmm/Campbell_et_al

    """
    if not isinstance(source_dirs, type([])):
        raise ConfigError("source_dirs parameter must be a list (get_HMM_sources_dictionary).")

    sources = {}
    allowed_chars_for_proper_sources = constants.allowed_chars.replace('.', '').replace('-', '')
    PROPER = lambda w: not len([c for c in w if c not in allowed_chars_for_proper_sources]) \
                       and len(w) >= 3 \
                       and w[0] not in '_0123456789'

    R = lambda f: open(os.path.join(source, f), 'r').readlines()[0].strip()
    for source in source_dirs:
        if source.endswith('/'):
            source = source[:-1]

        if not PROPER(os.path.basename(source)):
            raise ConfigError(f"One of the search database directories ({os.path.basename(source)}) contains characters "
                               "in its name anvio does not like. Directory names should be at least three characters long "
                               "and must not contain any characters but ASCII letters, digits and underscore")

        expected_files = ['reference.txt', 'kind.txt', 'genes.txt', 'genes.hmm.gz', 'target.txt', 'noise_cutoff_terms.txt']

        missing_files = [f for f in expected_files if not os.path.exists(os.path.join(source, f))]
        if missing_files:
            raise ConfigError(f"The HMM source '{os.path.basename(source)}' makes anvi'o unhappy. Each HMM source directory "
                              f"must contain a specific set of {len(expected_files)} files, and nothing more. See this URL "
                              f"for details: https://anvio.org/help/{anvio.anvio_version_for_help_docs}/artifacts/hmm-source/")

        empty_files = [f for f in expected_files if os.stat(os.path.join(source, f)).st_size == 0]
        if empty_files:
            raise ConfigError("One or more files for the HMM source '%s' seems to be empty. Which creates lots of "
                              "counfusion around these parts of the code. Anvi'o could set some defualts for you, "
                              "but it would be much better if you set your own defaults explicitly. You're not "
                              "sure what would make a good default for your HMM collection? Reach out to "
                              "a developer, and they will help you! Here are the files that are empty: %s." % \
                                    (os.path.basename(source), ', '.join(empty_files)))

        ref = R('reference.txt')
        kind = R('kind.txt')
        target = R('target.txt')
        noise_cutoff_terms = R('noise_cutoff_terms.txt')
        anvio_hmm_target_term_to_alphabet_and_context(target)

        domain = None
        if kind == 'singlecopy' and kind.count(':') == 0:
            raise ConfigError("This HMM profile seems to be a collection of single-copy core genes. Great. But for "
                              "this kind, you must also declare a 'domain' in your 'kind.txt' file. It is simple. "
                              "For instance, you could use 'singlecopy:bacteria', or 'singlecopy:archaea', or "
                              "'singlecopy:myspecificbranch'.")
        if kind.count(':') == 1:
            kind, domain = kind.split(':')

        if not PROPER(kind):
            raise ConfigError("'kind.txt' defines the kind of search this database offers. The kind term must be a single "
                               "word that is at least three characters long, and must not contain any characters but "
                               "ASCII letters, digits, and underscore. Here are some nice examples: 'singlecopy', "
                               "or 'pathogenicity', or 'noras_selection'. But yours is '%s'." % (kind))

        if domain and not PROPER(domain):
            raise ConfigError("That's lovely that you decided to specify a domain extension for your HMM collection in the "
                               "'kind.txt'. Although, your domain term is not a good one, as it must be a single "
                               "word that is at least three characters long, and without any characters but "
                               "ASCII letters, digits, and underscore. Confused? That's fine. Send an e-mail to the anvi'o "
                               "developers, and they will help you!")

        genes = get_TAB_delimited_file_as_dictionary(os.path.join(source, 'genes.txt'), column_names=['gene', 'accession', 'hmmsource'])

        sanity_check_hmm_model(os.path.join(source, 'genes.hmm.gz'), genes)

        sources[os.path.basename(source)] = {'ref': ref,
                                             'kind': kind,
                                             'domain': domain,
                                             'genes': list(genes.keys()),
                                             'target': target,
                                             'noise_cutoff_terms': noise_cutoff_terms,
                                             'model': os.path.join(source, 'genes.hmm.gz')}

    return sources



def get_attribute_from_hmm_file(file_path, attribute):
    """
    Retrieves the value of attribute from an HMMER/3 formatted HMM file.
    - file_path: (str) absolute file path to the .HMM file
    - attribute: (str) the attribute to get from the HMM file
    """
    filesnpaths.is_file_exists(file_path)
    value = None
    with open(file_path) as hmm:
        for line in hmm.readlines():
            if line.startswith(attribute):
                value = [f.strip() for f in line.split(attribute) if len(f)][0]
                break

    if value is None:
        raise ValueError(f"The HMM file {file_path} did not contain {attribute}.")

    return value



def sanity_check_hmm_model(model_path, genes):
    genes = set(genes)
    genes_in_model = set([])
    accession_ids_in_model = []

    with gzip.open(model_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('NAME'):
                genes_in_model.add(line.split()[1])
            if line.startswith('ACC'):
                accession_ids_in_model.append(line.split()[1])

    if len(accession_ids_in_model) != len(set(accession_ids_in_model)):
        raise ConfigError("Accession IDs in your HMM model should be unique, however, the `genes.hmm.gz` "
                          "file for `%s` seems to have the same accession ID (the line that starts with `ACC`) "
                          "more than once :(" % (os.path.abspath(model_path).split('/')[-2]))

    if len(genes.difference(genes_in_model)):
        raise ConfigError("Some gene names in genes.txt file does not seem to be appear in genes.hmm.gz. "
                          "Here is a list of missing gene names: %s" % ', '.join(list(genes.difference(genes_in_model))))

    if len(genes_in_model.difference(genes)):
        raise ConfigError("Some gene names in genes.hmm.gz file does not seem to be appear in genes.txt. "
                          "Here is a list of missing gene names: %s" % ', '.join(list(genes_in_model.difference(genes))))



def sanity_check_pfam_accessions(pfam_accession_ids):
    """This function sanity checks a list of Pfam accession IDs

    Parameters
    ==========
    pfam_accession_ids: list
        list of possible Pfam accessions
    """
    not_pfam_accession_ids = [pfam_accession_id for pfam_accession_id in pfam_accession_ids if not pfam_accession_id.startswith("PF")]

    if len(not_pfam_accession_ids):
        raise ConfigError(f"The following accessions do not appear to be from Pfam because they do not "
                          f"start with \"PF\", please double check the following: {','.join(not_pfam_accession_ids)}")

