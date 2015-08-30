# -*- coding: utf-8

"""Lonely, helper functions that are broadly used and don't fit anywhere"""

import os
import sys
import time
import json
import copy
import socket
import textwrap
import subprocess
import multiprocessing

import anvio
import anvio.fastalib as u
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Progress
from anvio.errors import ConfigError
from anvio.sequence import Composition
from anvio.constants import IS_ESSENTIAL_FIELD, allowed_chars, digits, complements


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


# Mock progress object that will not report anything, for general clarity.
progress = Progress()
progress.verbose = False


def rev_comp(seq):
    return seq.translate(complements)[::-1]


class Multiprocessing:
    def __init__(self, target_function, num_thread = None):
        self.cpu_count = multiprocessing.cpu_count()
        self.num_thread = num_thread or (self.cpu_count - (int(round(self.cpu_count / 10.0)) or 1))
        self.target_function = target_function
        self.processes = []
        self.manager = multiprocessing.Manager()


    def get_data_chunks(self, data_array, spiral = False):
        data_chunk_size = (len(data_array) / self.num_thread) or 1
        data_chunks = []
        
        if len(data_array) <= self.num_thread:
            return [[chunk] for chunk in data_array]

        if spiral:
            for i in range(0, self.num_thread):
                data_chunks.append([data_array[j] for j in range(i, len(data_array), self.num_thread)])
            
            return data_chunks
        else:
            for i in range(0, self.num_thread):
                if i == self.num_thread - 1:
                    data_chunks.append(data_array[i * data_chunk_size:])
                else:
                    data_chunks.append(data_array[i * data_chunk_size:i * data_chunk_size + data_chunk_size])

        return data_chunks

                
    def run(self, args, name = None):
        t = multiprocessing.Process(name = name,
                                    target = self.target_function,
                                    args = args)
        self.processes.append(t)
        t.start()


    def get_empty_shared_array(self):
        return self.manager.list()


    def get_empty_shared_dict(self):
        return self.manager.dict()

    
    def get_shared_integer(self):
        return self.manager.Value('i', 0)


    def run_processes(self, processes_to_run, progress = Progress(verbose=False)):
        tot_num_processes = len(processes_to_run)
        sent_to_run = 0
        while 1:
            NumRunningProceses = lambda: len([p for p in self.processes if p.is_alive()])
            
            if NumRunningProceses() < self.num_thread and processes_to_run:
                for i in range(0, self.num_thread - NumRunningProceses()):
                    if len(processes_to_run):
                        sent_to_run += 1
                        self.run(processes_to_run.pop())

            if not NumRunningProceses() and not processes_to_run:
                # let the blastn program finish writing all output files.
                # FIXME: this is ridiculous. find a better solution.
                time.sleep(1)
                break

            progress.update('%d of %d done in %d threads (currently running processes: %d)'\
                                                         % (sent_to_run - NumRunningProceses(),
                                                            tot_num_processes,
                                                            self.num_thread,
                                                            NumRunningProceses()))
            time.sleep(1)


def get_available_port_num(start = 8080, look_upto_next_num_ports = 100, ip='0.0.0.0'):
    """Starts from 'start' and incrementally looks for an available port
       until 'start + look_upto_next_num_ports', and returns the first
       available one."""
    for p in range(start, start + look_upto_next_num_ports):
        if not is_port_in_use(p, ip):
            return p

    return None


def is_port_in_use(port, ip='0.0.0.0'):
    in_use = False
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex((ip, port))

    if(result == 0) :
        in_use = True

    sock.close()
    return in_use


def is_program_exists(program):
    IsExe = lambda p: os.path.isfile(p) and os.access(p, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        if IsExe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.expanduser(path).strip('"')
            exe_file = os.path.join(path, program)
            if IsExe(exe_file):
                return True

    raise ConfigError, "A anvio function needs '%s' to be installed on your system, but it doesn't seem to appear\
                        in your path :/ If you are certain you have it on your system (for instance you can run it\
                        by typing '%s' in your terminal window), you may want to send a detailed bug report. Sorry!"\
                        % (program, program)


def run_command(cmdline):
    try:
        ret_val = subprocess.call(cmdline, shell = True)
        if ret_val < 0:
            raise ConfigError, "command was terminated"
        else:
            return ret_val
    except OSError, e:
        raise ConfigError, "command was failed for the following reason: '%s' ('%s')" % (e, cmdline)


def store_array_as_TAB_delimited_file(a, output_path, header, exclude_columns = []):
    filesnpaths.is_output_file_writable(output_path)

    num_fields = len(a[0])

    if len(header) != num_fields:
        raise ConfigError, "store array: header length (%d) differs from data (%d)..." % (len(header), num_fields)

    for col in exclude_columns:
        if not col in header:
            raise ConfigError, "store array: column %s is not in the header array..."

    exclude_indices = set([header.index(c) for c in exclude_columns])

    header = [header[i] for i in range(0, len(header)) if i not in exclude_indices]

    f = open(output_path, 'w')
    f.write('%s\n' % '\t'.join(header))

    for row in a:
        f.write('\t'.join([str(row[i]) for i in range(0, num_fields) if i not in exclude_indices]) + '\n')

    f.close()
    return output_path


def store_dict_as_TAB_delimited_file(d, output_path, headers, file_obj = None):
    if not file_obj:
        filesnpaths.is_output_file_writable(output_path)

    if not file_obj:
        f = open(output_path, 'w')
    else:
        f = file_obj

    f.write('%s\n' % '\t'.join(headers))

    for k in sorted(d.keys()):
        line = [str(k)]
        for header in headers[1:]:
            try:
                val = d[k][header]
            except KeyError:
                raise ConfigError, "Header ('%s') is not found in the dict :/" % (header)
            except TypeError:
                raise ConfigError, "Your dictionary is not properly formatted to be exported\
                                    as a TAB-delimited file :/ You ask for '%s', but it is not\
                                    even a key in the dictionary" % (header)

            line.append(str(val) if type(val) != type(None) else '')

        f.write('%s\n' % '\t'.join(line))

    f.close()
    return output_path


def is_all_columns_present_in_TAB_delim_file(columns, file_path):
    columns = get_columns_of_TAB_delim_file(file_path)
    return False if len([False for c in columns if c not in columns]) else True


def HTMLColorToRGB(colorstring, scaled = True):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]

    if scaled:
        return (r / 255.0, g / 255.0, b / 255.0)
    else:
        return (r, g, b)


def get_random_colors_dict(keys):
    # FIXME: someone's gotta implement this
    # keys   : set(1, 2, 3, ..)
    # returns: {1: '#ffffff', 2: '#888888', 3: '#222222', ...}
    return dict([(k, None) for k in keys])


def get_columns_of_TAB_delim_file(file_path, include_first_column=False):
    if include_first_column:
        return open(file_path).readline().strip('\n').split('\t')
    else:
        return open(file_path).readline().strip('\n').split('\t')[1:]


def get_vectors_from_TAB_delim_matrix(file_path, cols_to_return=None, rows_to_return = []):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_tab_delimited(file_path)

    rows_to_return = set(rows_to_return)
    vectors = []
    id_to_sample_dict = {}
    sample_to_id_dict = {}

    input_matrix = open(file_path)
    columns = input_matrix.readline().strip().split('\t')[1:]

    fields_of_interest = []
    if cols_to_return:
        fields_of_interest = [columns.index(col) for col in cols_to_return]
    else:
        fields_of_interest = [f for f in range(0, len(columns)) if IS_ESSENTIAL_FIELD(columns[f])]

    # update columns:
    columns = [columns[i] for i in fields_of_interest]

    if not len(columns):
        raise ConfigError, "Only a subset (%d) of fields were requested by the caller, but none of them was found\
                            in the matrix (%s) :/" % (len(cols_to_return), file_path)

    id_counter = 0
    for line in input_matrix.readlines():
        row_name = line.strip().split('\t')[0]
        if rows_to_return and row_name not in rows_to_return:
                continue
        id_to_sample_dict[id_counter] = row_name
        fields = line.strip().split('\t')[1:]

        if fields_of_interest:
            vector = [float(fields[i]) for i in fields_of_interest]
        else:
            vector = [float(f) for f in fields]

        vectors.append(vector)

        id_counter += 1

    sample_to_id_dict = dict([(v, k) for k, v in id_to_sample_dict.iteritems()])

    return id_to_sample_dict, sample_to_id_dict, columns, vectors


def get_all_ids_from_fasta(input_file):
    fasta = u.SequenceSource(input_file)
    ids = []
    
    while fasta.next():
        ids.append(fasta.id) 

    return ids


def get_read_lengths_from_fasta(input_file):
    contig_lengths = {}

    fasta = u.SequenceSource(input_file)
    while fasta.next():
        contig_lengths[fasta.id] = len(fasta.seq)

    fasta.close()
    return contig_lengths


def get_GC_content_for_FASTA_entries(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    GC_content_dict = {}

    fasta = u.SequenceSource(file_path)
    while fasta.next():
        GC_content_dict[fasta.id] = get_GC_content_for_sequence(fasta.seq)

    return GC_content_dict


def get_GC_content_for_sequence(sequence):
    return Composition(sequence).GC_content


def get_N50(contig_lengths):
    h, S = sum(contig_lengths) / 2.0, 0

    for l in sorted(contig_lengths, reverse=True):
        S += l
        if h < S:
            return l


def get_cmd_line():
    c_argv = []
    for i in sys.argv:
        if ' ' in i:
            c_argv.append('"%s"' % i)
        else:
            c_argv.append(i)
    return ' '.join(c_argv)


def get_time_to_date(local_time, fmt='%Y-%m-%d %H:%M:%S'):
    try:
        local_time = float(local_time)
    except ValueError:
        raise ConfigError, "utils::get_time_to_date is called with bad local_time."

    return time.strftime(fmt, time.localtime(local_time))


def concatenate_files(dest_file, file_list):
    if not dest_file:
        raise ConfigError, "Destination cannot be empty."
    if not len(file_list):
        raise ConfigError, "File list cannot be empty."
    for f in file_list:
        filesnpaths.is_file_exists(f)
    filesnpaths.is_output_file_writable(dest_file)

    dest_file_obj = open(dest_file, 'w')
    for chunk_path in file_list:
        for line in open(chunk_path):
            dest_file_obj.write(line)

    dest_file_obj.close()
    return dest_file


def get_split_start_stops(contig_length, split_length):
    """Returns split start stop locations for a given contig length"""
    num_chunks = contig_length / split_length

    if num_chunks < 2:
        return [(0, contig_length)]

    chunks = []
    for i in range(0, num_chunks):
        chunks.append((i * split_length, (i + 1) * split_length),)
    chunks.append(((i + 1) * split_length, contig_length),)

    if (chunks[-1][1] - chunks[-1][0]) < (split_length / 2):
        # last chunk is too small :/ merge it to the previous one.
        last_tuple = (chunks[-2][0], contig_length)
        chunks.pop()
        chunks.pop()
        chunks.append(last_tuple)

    return chunks


def get_contigs_splits_dict(split_ids, splits_basic_info):
    """
    For a given list of split ids, create a dictionary of contig names
    that represents all parents as keys, and ordered splits as items.

    split_ids is a set of split IDs, splits_basic_info comes from the contigs database:
 
     >>> contigs_db = dbops.ContigsDatabase(contigs_db_path)
     >>> splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)
     >>> znnotation_db.disconnect()
     >>> x = get_contigs_splits_dict(set([contig_A_split_00001, contig_A_split_00002, contig_A_split_00004,
                                         contig_C_split_00003, contig_C_split_00004, contig_C_split_00005]),
                                    splits_basic_info)
     >>> print x
         {
             'contig_A': {
                             0: 'contig_A_split_00001',
                             1: 'contig_A_split_00002',
                             4: 'contig_A_split_00004'
                         },
             'contig_C': {
                             3: 'contig_C_split_00003',
                             4: 'contig_C_split_00004',
                             5: 'contig_C_split_00005'
                         }
         }
    """

    contigs_splits_dict = {}

    for split_id in split_ids:
        s = splits_basic_info[split_id]
        if s['parent'] in contigs_splits_dict:
            contigs_splits_dict[s['parent']][s['order_in_parent']] = split_id
        else:
            contigs_splits_dict[s['parent']] = {s['order_in_parent']: split_id}

    return contigs_splits_dict

def get_contig_name_to_splits_dict(splits_basic_info_dict, contigs_basic_info_dict):
    """
    Returns a dict for contig name to split name conversion.

    Here are the proper source of the input params:

        contigs_basic_info_dict = database.get_table_as_dict(t.contigs_info_table_name, string_the_key = True)
        splits_basic_info_dict  = database.get_table_as_dict(t.splits_info_table_name)
    """
    contig_name_to_splits_dict = {}

    for split_name in splits_basic_info_dict:
        parent = splits_basic_info_dict[split_name]['parent']
        if contig_name_to_splits_dict.has_key(parent):
            contig_name_to_splits_dict[parent].append(split_name)
        else:
            contig_name_to_splits_dict[parent] = [split_name]

    return contig_name_to_splits_dict


def check_sample_id(sample_id):
    if sample_id:
        if sample_id[0] in digits:
            raise ConfigError, "Sample names can't start with digits. Long story. Please specify a sample name\
                                that starts with an ASCII letter (you may want to check '-s' parameter to set\
                                a sample name if your client permits (otherwise you are going to have to edit\
                                your input files))."

        allowed_chars_for_samples = allowed_chars.replace('-', '').replace('.', '')
        if len([c for c in sample_id if c not in allowed_chars_for_samples]):
            raise ConfigError, "Sample name ('%s') contains characters that anvio does not like. Please\
                                limit the characters that make up the project name to ASCII letters,\
                                digits, and the underscore character ('_')." % sample_id


def is_this_name_OK_for_database(variable_name, content, allowed_chars = allowed_chars.replace('.', '')):
    if content[0] in digits:
        raise ConfigError, "Sorry, '%s' can't start with a digit. Long story. Please specify a sample name\
                            that starts with an ASCII letter." % variable_name

    if len([c for c in content if c not in allowed_chars]):
        raise ConfigError, "Well, '%s' parameter contains characters that anvi'o does not like. Please\
                            limit the characters to ASCII letters, digits, the underscore and dash\
                            characters ('_', '-')." % variable_name



def check_contig_names(contig_names, dont_raise = False):
    all_characters_in_contig_names = set(''.join(contig_names))
    characters_anvio_doesnt_like = [c for c in all_characters_in_contig_names if c not in allowed_chars]
    if len(characters_anvio_doesnt_like):
        if dont_raise:
            return False

        raise ConfigError, "The name of at least one contig in your BAM file %s anvio does not\
                            like (%s). Please go back to your original files and make sure that\
                            the characters in contig names are limited to to ASCII letters,\
                            digits. Names can also contain underscore ('_'), dash ('-') and dot ('.')\
                            characters. anvio knows how much work this may require for you to go back and\
                            re-generate your BAM files and is very sorry for asking you to do that, however,\
                            it is critical for later steps in the analysis." \
                                % ("contains multiple characters" if len(characters_anvio_doesnt_like) > 1 else "contains a character",
                                   ", ".join(['"%s"' % c for c in characters_anvio_doesnt_like]))

    return True


def get_FASTA_file_as_dictionary(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    d = {}

    fasta = u.SequenceSource(file_path)
    while fasta.next():
        d[fasta.id] = fasta.seq

    return d


def store_dict_as_FASTA_file(d, output_file_path, wrap_from = 200):
    filesnpaths.is_output_file_writable(output_file_path)
    output = open(output_file_path, 'w')

    for key in d:
        output.write('>%s\n' % key)
        output.write('%s\n' % textwrap.fill(d[key], wrap_from, break_on_hyphens = False))

    output.close()
    return True


def is_ascii_only(text):
    """test whether 'text' is composed of ASCII characters only"""
    try:
        text.decode('ascii')
    except UnicodeDecodeError:
        return False
    return True


def get_TAB_delimited_file_as_dictionary(file_path, expected_fields = None, dict_to_append = None, column_names = None,\
                                        column_mapping = None, indexing_field = 0, assign_none_for_missing = False,\
                                        separator = '\t', no_header = False, ascii_only = False):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_tab_delimited(file_path, separator = separator)

    f = open(file_path)

    # learn the number of fields and reset the file:
    num_fields = len(f.readline().strip('\n').split(separator))
    f.seek(0)

    # if there is no file header, make up a columns list:
    if no_header and not column_names:
        column_names = ['column_%05d' % i for i in range(0, num_fields)]

    if column_names:
        columns = column_names

        if num_fields != len(columns):
            raise  ConfigError, "Number of column names declared (%d) differs from the number of columns\
                                 found (%d) in the matrix ('%s') :/" % (len(columns), num_fields, file_path)

        # now we set the column names. if the file had its header, we must discard
        # the first line. so here we go:
        if not no_header:
            f.readline()
    else:
        columns = f.readline().strip('\n').split(separator)

    if expected_fields:
        for field in expected_fields:
            if field not in columns:
                raise ConfigError, "The file '%s' does not contain the right type of header. It was expected\
                                    to have these: '%s', however it had these: '%s'" % (file_path,
                                                                                        ', '.join(expected_fields),
                                                                                        ', '.join(columns[1:]))

    d = {}
    line_counter = 0

    for line in f.readlines():
        if ascii_only:
            if not is_ascii_only(line):
                raise ConfigError, "The input file conitans non-ascii characters at line number %d. Those lines\
                                    either should be removed, or edited." % (line_counter + 2)

        line_fields = line.strip('\n').split(separator)

        if column_mapping:
            updated_line_fields = []
            for i in range(0, len(line_fields)):
                try:
                    updated_line_fields.append(column_mapping[i](line_fields[i]))
                except NameError:
                    raise ConfigError, "Mapping function '%s' did not work on value '%s'. These functions can be native\
                                        Python functions, such as 'str', 'int', or 'float', or anonymous functions\
                                        defined using lambda notation." % (column_mapping[i], line_fields[i])
                except TypeError:
                    raise ConfigError, "Mapping function '%s' does not seem to be a proper Python function :/" % column_mapping[i]
                except ValueError:
                    raise ConfigError, "Mapping funciton '%s' did not like the value '%s' in column number %d\
                                        of the matrix :/" % (column_mapping[i], line_fields[i], i + 1)
            line_fields = updated_line_fields 

        if indexing_field == -1:
            entry_name = 'line__%09d__' % line_counter
        else:
            entry_name = line_fields[indexing_field]

        d[entry_name] = {}

        for i in range(0, len(columns)):
            if i == indexing_field:
                continue
            d[entry_name][columns[i]] = line_fields[i]

        line_counter += 1

    # we have the dict, but we will not return it the way it is if its supposed to be appended to an
    # already existing dictionary.
    if dict_to_append:
        # we don't want to through keys in d each time we want to add stuff to 'dict_to_append', so we keep keys we
        # find in the first item in the dict in another variable. this is potentially very dangerous if not every
        # item in 'd' has identical set of keys.
        keys = d.values()[0].keys()

        for entry in dict_to_append:
            if entry not in d:
                # so dict to append is missing a key that is in the dict to be appended. if the user did not
                # ask us to add None for these entries via none_for_missing, we are going to make a noise,
                # otherwise we will tolerate it.
                if not assign_none_for_missing:
                    raise ConfigError, "Appending entries to the already existing dictionary from file '%s' failed\
                                        as the entry %s does not appear to be in the file." % (file_path, entry)
                else:
                    for key in keys:
                        dict_to_append[entry][key] = None
            else:
                for key in keys:
                    dict_to_append[entry][key] = d[entry][key]

        return dict_to_append

    return d


def get_filtered_dict(d, property, values):
    # removes any entry from d, where the value of the 'property' of items in d does not match
    # with desired 'values'
    d = copy.deepcopy(d)
    if type(values) != type(set([])):
        raise ConfigError, "get_filtered_dict: values must be type of set([])."

    entry_ids_to_remove = [entry_id for entry_id in d if d[entry_id][property] not in values]
    for entry_id in entry_ids_to_remove:
        d.pop(entry_id)

    return d


def get_HMM_sources_dictionary(source_dirs=[]):
    if type(source_dirs) != type([]):
        raise ConfigError, "source_dirs parameter must be a list (get_HMM_sources_dictionary)."

    sources = {}
    allowed_chars_for_proper_sources = allowed_chars.replace('.', '').replace('-', '')
    PROPER = lambda w: not len([c for c in w if c not in allowed_chars_for_proper_sources]) \
                       and len(w) >= 3 \
                       and w[0] not in '_0123456789'

    for source in source_dirs:
        if source.endswith('/'):
            source = source[:-1]

        if not PROPER(os.path.basename(source)):
            raise ConfigError, "One of the search database directories ('%s') contains characters in its name\
                                anvio does not like. Directory names should be at least three characters long\
                                and must not contain any characters but ASCII letters, digits and\
                                underscore" % os.path.basename(source)

        for f in ['reference.txt', 'kind.txt', 'genes.txt', 'genes.hmm.gz']:
            if not os.path.exists(os.path.join(source, f)):
                raise ConfigError, "Each search database directory must contain following files:\
                                    'kind.txt', 'reference.txt', 'genes.txt', and 'genes.hmm.gz'. %s does not seem\
                                    to be a proper source." % os.path.basename(source)

        ref = open(os.path.join(source, 'reference.txt')).readlines()[0].strip()
        kind = open(os.path.join(source, 'kind.txt')).readlines()[0].strip()
        if not PROPER(kind):
            raise ConfigError, "'kind.txt' defines the kind of search this database offers. This file must contain a single\
                                word that is at least three characters long, and must not contain any characters but\
                                ASCII letters, digits, and underscore. Here are some nice examples: 'singlecopy',\
                                or 'pathogenicity', or 'noras_selection'. But yours is '%s'." % (kind)

        genes = get_TAB_delimited_file_as_dictionary(os.path.join(source, 'genes.txt'), column_names = ['gene', 'accession', 'hmmsource'])

        sources[os.path.basename(source)] = {'ref': ref,
                                             'kind': kind,
                                             'genes': genes.keys(),
                                             'model': os.path.join(source, 'genes.hmm.gz')}

    return sources


def get_missing_programs_for_hmm_analysis():
    missing_programs = []
    for p in ['prodigal', 'hmmscan']:
        try:
            is_program_exists(p)
        except ConfigError:
            missing_programs.append(p)
    return missing_programs


