import os
import shutil
import hashlib
import textwrap

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Run
from anvio.errors import ConfigError

from anvio.utils.sequences import get_GC_content_for_sequence
from anvio.utils.files import get_TAB_delimited_file_as_dictionary

def split_fasta(input_file_path, parts=1, file_name_prefix=None, shuffle=False, output_dir=None, return_number_of_sequences=False):
    """Splits a given FASTA file into multiple parts.

    Please note that this function will not clean after itself. You need to take care of the
    output files in context.

    Parameters
    ==========
    input_file_path : str
        FASTA-formatted flat text file to be split
    parts : int
        Number of parts the input file to be split into
    file_name_prefix : str
        Preferably a single-word prefix for the output files
    shuffle : bool
        Whether input sequences should be randomly shuffled (so the input sequences
        randomly distribute across output files)
    output_dir : str, path
        Output directory. By default, anvi'o will store things in a new directory under
        the system location for temporary files
    return_number_of_sequences : bool
        Whether to return the number of sequences in the original fasta file

    Returns
    =======
    output_file_paths : list
        Array with `parts` number of elements where each item is an output file path
    length : int, optional
        The number of sequences in the original fasta file. Only returnd if 'return_number_of_sequences' is True
    """
    if not file_name_prefix:
        file_name_prefix = os.path.basename(input_file_path)
    else:
        if '/' in file_name_prefix:
            raise ConfigError("File name prefix for split fasta can't contain slash characters. It is not "
                              "supposed to be a path after all :/")

    # check input
    filesnpaths.is_file_fasta_formatted(input_file_path)

    # check output
    if not output_dir:
        output_dir = filesnpaths.get_temp_directory_path()
    else:
        filesnpaths.gen_output_directory(output_dir)
        filesnpaths.is_output_dir_writable(output_dir)

    source = u.ReadFasta(input_file_path, quiet=True)
    length = len(source.ids)

    if length < parts:
        parts = length

    chunk_size = length // parts

    output_file_paths = []

    GET_OUTPUT_FILE_PATH = lambda p: os.path.join(output_dir, ".".join([file_name_prefix, str(p)]))

    if shuffle:
        output_file_paths = [f'{GET_OUTPUT_FILE_PATH(part_no)}' for part_no in range(parts)]
        output_fastas = [u.FastaOutput(file_name) for file_name in output_file_paths]

        # The first sequence goes to the first outfile, the second seq to the second outfile, and so on.
        for seq_idx, (seq_id, seq) in enumerate(zip(source.ids, source.sequences)):
            which = seq_idx % parts
            output_fastas[which].write_id(seq_id)
            output_fastas[which].write_seq(seq)

        for output_fasta in output_fastas:
            output_fasta.close()
    else:
        for part_no in range(parts):
            output_file = GET_OUTPUT_FILE_PATH(part_no)

            output_fasta = u.FastaOutput(output_file)

            chunk_start = chunk_size * part_no
            chunk_end   = chunk_start + chunk_size

            if (part_no + 1 == parts):
                # if this is the last chunk make sure it contains everything till end.
                chunk_end = length

            for i in range(chunk_start, chunk_end):
                output_fasta.write_id(source.ids[i])
                output_fasta.write_seq(source.sequences[i])

            output_fasta.close()
            output_file_paths.append(output_file)

    source.close()

    if return_number_of_sequences:
        return output_file_paths, length

    else:
        return output_file_paths



def get_num_sequences_in_fasta(input_file):
    fasta = u.SequenceSource(input_file)
    num_sequences = 0

    while next(fasta):
        num_sequences += 1

    return num_sequences



def get_all_ids_from_fasta(input_file):
    fasta = u.SequenceSource(input_file)
    ids = []

    while next(fasta):
        ids.append(fasta.id)

    return ids



def check_fasta_id_formatting(fasta_path):
    fasta = u.SequenceSource(fasta_path)

    while next(fasta):
        characters_anvio_doesnt_like = [
            c for c in set(fasta.id) if c not in constants.allowed_chars]

        if len(characters_anvio_doesnt_like):
            raise ConfigError(
                "At least one of the deflines in your FASTA file "
                "does not comply with the 'simple deflines' requirement of Anvi'o. "
                "You can either use the script, `anvi-script-reformat-fasta`, "
                "to take care of this issue, or read this section in the tutorial "
                "to understand the reason behind this requirement "
                "(Anvi'o is very upset for making you do this): %s"
                % "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file")

        try:
            int(fasta.id)
            is_int = True
        except:
            is_int = False
        if is_int:
            raise ConfigError(
                "At least one of the deflines in your FASTA file "
                "(well, this one to be precise: '%s') looks like a number. "
                "For reasons we can't really justify, "
                "Anvi'o does not like those numeric names, "
                "and hereby asks you to make sure every tRNA-seq name "
                "contains at least one alphanumeric character :/ "
                "Meanwhile we, the Anvi'o developers, are both surprised by and thankful for "
                "your endless patience with such eccentric requests. "
                "You the real MVP." % fasta.id)

    fasta.close()



def check_fasta_id_uniqueness(fasta_path):
    all_ids_in_FASTA = get_all_ids_from_fasta(fasta_path)
    total_num_seqs = len(all_ids_in_FASTA)
    if total_num_seqs != len(set(all_ids_in_FASTA)):
        raise ConfigError(
            "Every sequence in the input FASTA file must have a unique ID. You know...")



def get_read_lengths_from_fasta(input_file):
    contig_lengths = {}

    fasta = u.SequenceSource(input_file)
    while next(fasta):
        contig_lengths[fasta.id] = len(fasta.seq)

    fasta.close()
    return contig_lengths



def remove_sequences_with_only_gaps_from_fasta(input_file_path, output_file_path, inplace=True):
    filesnpaths.is_file_fasta_formatted(input_file_path)
    filesnpaths.is_output_file_writable(output_file_path)

    total_num_sequences = 0
    num_sequences_removed = 0
    input_fasta = u.SequenceSource(input_file_path)
    clean_fasta = u.FastaOutput(output_file_path)

    while next(input_fasta):
        total_num_sequences += 1
        if input_fasta.seq.count('-') == len(input_fasta.seq):
            num_sequences_removed += 1
        else:
            clean_fasta.store(input_fasta, split=False)

    if inplace:
        if num_sequences_removed:
            shutil.move(output_file_path, input_file_path)
        else:
            os.remove(output_file_path)

    return total_num_sequences, num_sequences_removed



def get_FASTA_file_as_dictionary(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    d = {}

    fasta = u.SequenceSource(file_path)
    while next(fasta):
        d[fasta.id] = fasta.seq

    return d



def get_GC_content_for_FASTA_entries(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    GC_content_dict = {}

    fasta = u.SequenceSource(file_path)
    while next(fasta):
        GC_content_dict[fasta.id] = get_GC_content_for_sequence(fasta.seq)

    return GC_content_dict



def unique_FASTA_file(input_file_path, output_fasta_path=None, names_file_path=None, store_frequencies_in_deflines=True):
    filesnpaths.is_file_exists(input_file_path)

    if not output_fasta_path:
        output_fasta_path = input_file_path + '.unique'

    if not names_file_path:
        names_file_path = output_fasta_path + '.names'

    if output_fasta_path == names_file_path:
        raise ConfigError("I can't unique this. Output FASTA file path can't be identical to "
                           "the names file path...")

    if output_fasta_path == input_file_path or names_file_path == input_file_path:
        raise ConfigError("Anvi'o will not unique this. Output FASTA path and names file path should "
                           "be different from the input file path...")

    filesnpaths.is_output_file_writable(output_fasta_path)
    filesnpaths.is_output_file_writable(names_file_path)

    input_fasta = u.SequenceSource(input_file_path, unique=True)
    output_fasta = u.FastaOutput(output_fasta_path)
    names_file = open(names_file_path, 'w')

    names_dict = {}
    while next(input_fasta):
        output_fasta.store(input_fasta, split=False, store_frequencies=store_frequencies_in_deflines)
        names_file.write('%s\t%s\n' % (input_fasta.id, ','.join(input_fasta.ids)))

        names_dict[input_fasta.id] = input_fasta.ids

    output_fasta.close()
    names_file.close()

    return output_fasta_path, names_file_path, names_dict



def store_dict_as_FASTA_file(d, output_file_path, wrap_from=200):
    filesnpaths.is_output_file_writable(output_file_path)
    output = open(output_file_path, 'w')

    for key in d:
        output.write('>%s\n' % key)
        if wrap_from:
            output.write('%s\n' % textwrap.fill(d[key], wrap_from, break_on_hyphens=False))
        else:
            output.write('%s\n' % (d[key]))

    output.close()
    return True



def export_sequences_from_contigs_db(contigs_db_path, output_file_path, seq_names_to_export=None, splits_mode=False, rna_alphabet=False, truncate=True, just_do_it=False, run=Run()):
    """Export sequences from a contigs database."""
    filesnpaths.is_output_file_writable(output_file_path)

    contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)
    contig_sequences_dict = contigs_db.get_table_as_dict(t.contig_sequences_table_name, string_the_key = True)
    splits_info_dict = contigs_db.get_table_as_dict(t.splits_info_table_name)
    contigs_db.disconnect()

    output_fasta = u.FastaOutput(output_file_path)

    FORMAT = lambda seq: seq.replace('T', 'U') if rna_alphabet else seq

    if not seq_names_to_export:
        if splits_mode:
            seq_names_to_export = sorted(splits_info_dict.keys())
        else:
            seq_names_to_export = sorted(contig_sequences_dict.keys())
    else:
        contig_names = [contig_name for contig_name in seq_names_to_export if contig_name in contig_sequences_dict]
        split_names = [split_name for split_name in seq_names_to_export if split_name in splits_info_dict]
        missing_names = [name for name in seq_names_to_export if name not in contig_names and name not in split_names]

        if splits_mode:
          mode = "splits"
          appropriate_seq_names = split_names

        else:
          mode = "contigs"
          appropriate_seq_names = contig_names

        if len(appropriate_seq_names) < len(seq_names_to_export):
          if just_do_it:
              run.warning("Not all the sequences you requested are %s in this CONTIGS.db. %d names are contigs, "
                          "%d are splits, and %d are neither. BUT you're in just-do-it mode and we know you're in charge, so we'll "
                          "proceed using any appropriate names." % \
                          (mode, len(contig_names), len(split_names), len(missing_names),))
              seq_names_to_export = appropriate_seq_names
          else:
              raise ConfigError("Not all the sequences you requested are %s in this CONTIGS.db. %d names are contigs, "
                                "%d are splits, and %d are neither. If you want to live on the edge and try to "
                                "proceed using any appropriate names, try out the `--just-do-it` flag." % \
                                (mode, len(contig_names), len(split_names), len(missing_names)))

    for seq_name in seq_names_to_export:
        if splits_mode:
            s = splits_info_dict[seq_name]
            sequence = FORMAT(contig_sequences_dict[s['parent']]['sequence'][s['start']:s['end']])
        else:
            sequence = FORMAT(contig_sequences_dict[seq_name]['sequence'])

        output_fasta.write_id(seq_name)
        output_fasta.write_seq(sequence, split=truncate)

    return True



def create_fasta_dir_from_sequence_sources(genome_desc, fasta_txt=None):
    """genome_desc is an instance of GenomeDescriptions"""

    from anvio.summarizer import ArgsTemplateForSummarizerClass, ProfileSummarizer, Bin

    if genome_desc is None and fasta_txt is None:
        raise ConfigError("Anvi'o was given no internal genomes, no external genomes, and no fasta "
                          "files. Although anvi'o can technically go ahead and create a temporary "
                          "FASTA directory, what's the point if there's nothing to do?")

    temp_dir = filesnpaths.get_temp_directory_path()
    hash_to_name = {}
    name_to_path = {}
    genome_names = set([])
    file_paths = set([])

    if genome_desc is not None:
        # first figure out internal genomes that are bound by the same collection name and
        # profile db path
        genome_subsets = {}
        for entry in genome_desc.genomes.values():
            if 'bin_id' in entry:
                # if we are here, this entry represents an internal genome. we will add this genome
                # to genome_subsets data structure to be processed later.
                profile_and_collection_descriptor = '_'.join([entry['profile_db_path'], entry['collection_id']])
                if profile_and_collection_descriptor in genome_subsets:
                    genome_subsets[profile_and_collection_descriptor]['genome_name_bin_name_tpl'].add((entry['name'], entry['bin_id']),)
                else:
                    genome_subsets[profile_and_collection_descriptor] = {'genome_name_bin_name_tpl': set([(entry['name'], entry['bin_id'])]),
                                                                         'profile_db_path': entry['profile_db_path'],
                                                                         'contigs_db_path': entry['contigs_db_path'],
                                                                         'collection_id': entry['collection_id']}
            else:
                # If we are here, this means this is an external genome, so we can basically take care of it here immediately.
                genome_name = entry['name']
                genome_names.add(genome_name)
                contigs_db_path = genome_desc.genomes[genome_name]['contigs_db_path']
                hash_for_output_file = hashlib.sha256(genome_name.encode('utf-8')).hexdigest()
                hash_to_name[hash_for_output_file] = genome_name

                path = os.path.join(temp_dir, hash_for_output_file + '.fa')
                file_paths.add(path)

                name_to_path[genome_name] = path

                export_sequences_from_contigs_db(contigs_db_path, path)

        # when we are here, all we have are interanl genomes as genome subsets.
        for genome_subset in genome_subsets.values():
            args = ArgsTemplateForSummarizerClass()
            args.contigs_db = genome_subset['contigs_db_path']
            args.profile_db = genome_subset['profile_db_path']
            args.collection_name = genome_subset['collection_id']
            args.output_dir = filesnpaths.get_temp_directory_path(just_the_path=True)
            args.quick_summary = True

            # note that we're initializing the summary class only for once for a given
            # genome subset
            summary = ProfileSummarizer(args, r=Run(verbose=False))
            summary.init()

            for genome_name, bin_name in genome_subset['genome_name_bin_name_tpl']:
                genome_names.add(genome_name)

                hash_for_output_file = hashlib.sha256(genome_name.encode('utf-8')).hexdigest()
                hash_to_name[hash_for_output_file] = genome_name
                path = os.path.join(temp_dir, hash_for_output_file + '.fa')
                file_paths.add(path)
                name_to_path[genome_name] = path

                bin_summary = Bin(summary, bin_name)

                with open(path, 'w') as fasta:
                    fasta.write(bin_summary.get_bin_sequence())


    if fasta_txt is not None:
        fastas = get_TAB_delimited_file_as_dictionary(fasta_txt, expected_fields=['name', 'path'], only_expected_fields=True)

        # make sure every entry in the fasta_txt has a path that exists
        genomes_missing_fasta_files = [g for g, e in fastas.items() if not os.path.exists(e['path'])]

        if len(genomes_missing_fasta_files):
            if len(genomes_missing_fasta_files) == 1:
                msg = (f"One of the genome entries in your fasta-txt file, namely '{genomes_missing_fasta_files[0]}' does "
                       f"not seem to have its FASTA file at the location it is mentioned in the file :/ ")
            else:
                msg = (f"Multiple genome entries in your fasta-txt file have a FASTA file path that don't match to an "
                       f"existing FASTA file :/ Here are the list of offenders: {', '.join(genomes_missing_fasta_files)}. ")

            msg += "Please correct your fasta-txt, and try again."

            raise ConfigError(f"{msg}")

        for name in fastas.keys():
            genome_names.add(name)
            hash_for_output_file = hashlib.sha256(name.encode('utf-8')).hexdigest()
            hash_to_name[hash_for_output_file] = name

            source = fastas[name]['path']
            path = os.path.join(temp_dir, hash_for_output_file + '.fa')
            file_paths.add(path)

            name_to_path[name] = path

            with open(path, 'w') as dest:
                with open(source, 'r') as src:
                    dest.write(src.read())

    return temp_dir, hash_to_name, genome_names, name_to_path

