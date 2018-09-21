# -*- coding: utf-8
# pylint: disable=line-too-long

"""Lonely, helper functions that are broadly used and don't fit anywhere"""

import os
import sys
import gzip
import time
import socket
import shutil
import psutil
import smtplib
import hashlib
import textwrap
import webbrowser
import subprocess
import configparser
import multiprocessing
import urllib.request, urllib.error, urllib.parse

import numpy as np
import pandas as pd

from email.mime.text import MIMEText
from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Run, Progress, SuppressAllOutput, get_date
from anvio.errors import ConfigError, FilesNPathsError
from anvio.sequence import Composition

with SuppressAllOutput():
    from ete3 import Tree


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"

# for full output
pd.options.display.max_columns=100
pd.options.display.max_rows=100

# Mock progress object that will not report anything, for general clarity.
progress = Progress()
progress.verbose = False

run = Run()
run.verbose = False

def rev_comp(seq):
    return seq.translate(constants.complements)[::-1]


class Multiprocessing:
    def __init__(self, target_function, num_thread=None):
        self.cpu_count = multiprocessing.cpu_count()
        self.num_thread = num_thread or (self.cpu_count - (int(round(self.cpu_count / 10.0)) or 1))
        self.target_function = target_function
        self.processes = []
        self.manager = multiprocessing.Manager()


    def get_data_chunks(self, data_array, spiral=False):
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


    def run(self, args, name=None):
        t = multiprocessing.Process(name=name,
                                    target=self.target_function,
                                    args=args)
        self.processes.append(t)
        t.start()


    def get_empty_shared_array(self):
        return self.manager.list()


    def get_empty_shared_dict(self):
        return self.manager.dict()


    def get_shared_integer(self):
        return self.manager.Value('i', 0)

    def run_processes(self, processes_to_run, progress=Progress(verbose=False)):
        tot_num_processes = len(processes_to_run)
        sent_to_run = 0
        while True:
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


def get_total_memory_usage():
    current_process = psutil.Process(os.getpid())
    mem = current_process.memory_info().rss
    for child in current_process.children(recursive=True):
        try:
            mem += child.memory_info().rss
        except:
            pass

    return human_readable_file_size(mem)


def get_predicted_type_of_items_in_a_dict(d, key):
    """Gets a dictionary `d` and a `key` in it, and returns a type function.

    It is a bit counter intuitive. dictionary should look like this:

        d = {'x': {'key': item, (...)},
             'y': {'key': item, (...)},
             (...),
            }

    This is a shitty funciton, but there was a real need for it, so here we are :/
    """

    items = [x[key] for x in d.values()]

    if(set(items) == set([None])):
        # all items is of type None.
        return None

    not_float = False
    for item in items:
        if type(item) == float:
            return float

        if type(item) == int:
            return int

        try:
            float(item or 0)
        except ValueError:
            not_float = True
            break

    if not_float:
        return str
    else:
        for item in items:
            try:
                if int(item or 0) == float(item or 0):
                    continue
                else:
                    return float
            except ValueError:
                return float

        return int


def human_readable_file_size(nbytes):
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])


def get_port_num(port_num = 0, ip='0.0.0.0', run=run):
    """Get a port number for the `ip` address."""

    try:
        port_num = int(port_num) if port_num else 0
    except Exception as e:
        raise ConfigError("Not a happy port number :/ %s." % e)

    if not port_num:
        port_num = get_next_available_port_num(constants.default_port_number)

        if not port_num:
            raise ConfigError("Anvi'o searched a bunch of port numbers starting from %d, but failed\
                                to find an available one for you. Maybe you should specify one :/")
    else:
        if is_port_in_use(port_num):
            raise ConfigError("The port number %d seems to be in use :/" % port_num)

    if os.getuid() and port_num < 1024:
        run.warning("Using the port number %d requires superuser priviliges, which your user does not\
                     seem to have. Since anvi'o does not know anything about your system configuraiton,\
                     you are free to go for now. But be prepared for a failed attempt to use this port\
                     number to serve stuff." % port_num)

    return port_num


def get_next_available_port_num(start=constants.default_port_number, look_upto_next_num_ports=100, ip='0.0.0.0'):
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

    if result == 0:
        in_use = True

    sock.close()
    return in_use


def is_program_exists(program, dont_raise=False):
    IsExe = lambda p: os.path.isfile(p) and os.access(p, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        if IsExe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.expanduser(path).strip('"')
            exe_file = os.path.join(path, program)
            if IsExe(exe_file):
                return exe_file

    if dont_raise:
        return False

    raise ConfigError("An anvi'o function needs '%s' to be installed on your system, but it doesn't seem to appear\
                        in your path :/ If you are certain you have it on your system (for instance you can run it\
                        by typing '%s' in your terminal window), you may want to send a detailed bug report. Sorry!"\
                        % (program, program))


def format_cmdline(cmdline):
    """Takes a cmdline for `run_command` or `run_command_STDIN`, and makes it beautiful."""
    if not cmdline or (not isinstance(cmdline, str) and not isinstance(cmdline, list)):
        raise ConfigError("You made ultis::format_cmdline upset. The parameter you sent to run kinda sucks. It should be string\
                            or list type. Note that the parameter `shell` for subprocess.call in this `run_command` function\
                            is always False, therefore if you send a string type, it will be split into a list prior to being\
                            sent to subprocess.")

    if isinstance(cmdline, str):
        cmdline = [str(x) for x in cmdline.split(' ')]
    else:
        cmdline = [str(x) for x in cmdline]

    return cmdline


def gzip_compress_file(input_file_path, output_file_path=None, keep_original=False):
    filesnpaths.is_file_exists(input_file_path)

    if not output_file_path:
        output_file_path = input_file_path + '.gz'

    filesnpaths.is_output_file_writable(output_file_path)

    import gzip
    with open(input_file_path, 'rb') as f_in, gzip.open(output_file_path, 'wb') as f_out:
        f_out.writelines(f_in)

    if not keep_original:
        os.remove(input_file_path)


def gzip_decompress_file(input_file_path, output_file_path=None, keep_original=True):
    filesnpaths.is_file_exists(input_file_path)

    if not input_file_path.endswith('.gz'):
        raise ConfigError("gzip_decompress_file function is upset because your input file ('%s') does not\
                           end with a '.gz' extension :(")

    if not output_file_path:
        output_file_path = input_file_path[:-3]

    filesnpaths.is_output_file_writable(output_file_path)

    import gzip
    with gzip.open(input_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
        f_out.writelines(f_in)

    if not keep_original:
        os.remove(input_file_path)


def run_command(cmdline, log_file_path, first_line_of_log_is_cmdline=True, remove_log_file_if_exists=True):
    """Uses subprocess.call to run your `cmdline`"""
    cmdline = format_cmdline(cmdline)

    filesnpaths.is_output_file_writable(log_file_path)

    if remove_log_file_if_exists and os.path.exists(log_file_path):
        os.remove(log_file_path)

    try:
        if first_line_of_log_is_cmdline:
            with open(log_file_path, "a") as log_file: log_file.write('# DATE: %s\n# CMD LINE: %s\n' % (get_date(), ' '.join(cmdline)))

        log_file = open(log_file_path, 'a')
        ret_val = subprocess.call(cmdline, shell=False, stdout=log_file, stderr=subprocess.STDOUT)
        log_file.close()

        if ret_val < 0:
            raise ConfigError("command was terminated")
        else:
            return ret_val
    except OSError as e:
        raise ConfigError("command was failed for the following reason: '%s' ('%s')" % (e, cmdline))


def run_command_STDIN(cmdline, log_file_path, input_data, first_line_of_log_is_cmdline=True, remove_log_file_if_exists=True):
    """Uses subprocess.Popen and sends data to your `cmdline` through STDIN"""
    cmdline = format_cmdline(cmdline)

    filesnpaths.is_output_file_writable(log_file_path)

    if remove_log_file_if_exists and os.path.exists(log_file_path):
        os.remove(log_file_path)

    try:
        if first_line_of_log_is_cmdline:
            with open(log_file_path, "a") as log_file: log_file.write('# DATE: %s\n# CMD LINE: %s\n' % (get_date(), ' '.join(cmdline)))

        p = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        ret_val = p.communicate(input=input_data.encode('utf-8'))[0]
        return ret_val.decode()
    except OSError as e:
        raise ConfigError("command was failed for the following reason: '%s' ('%s')" % (e, cmdline))


def get_command_output_from_shell(cmd_line):
    ret_code = 0

    try:
        out_bytes = subprocess.check_output(cmd_line.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        out_bytes = e.output.decode("utf-8")
        ret_code = e.returncode

    return out_bytes, ret_code


def store_array_as_TAB_delimited_file(a, output_path, header, exclude_columns=[]):
    filesnpaths.is_output_file_writable(output_path)

    num_fields = len(a[0])

    if len(header) != num_fields:
        raise ConfigError("store array: header length (%d) differs from data (%d)..." % (len(header), num_fields))

    for col in exclude_columns:
        if not col in header:
            raise ConfigError("store array: column %s is not in the header array...")

    exclude_indices = set([header.index(c) for c in exclude_columns])

    header = [header[i] for i in range(0, len(header)) if i not in exclude_indices]

    f = open(output_path, 'w')
    f.write('%s\n' % '\t'.join(header))

    for row in a:
        f.write('\t'.join([str(row[i]) for i in range(0, num_fields) if i not in exclude_indices]) + '\n')

    f.close()
    return output_path


def store_dataframe_as_TAB_delimited_file(d, output_path, columns=None, include_index=False, index_label="index", naughty_characters=[-np.inf, np.inf], rep_str=""):
    """
    Stores a pandas DataFrame as a tab-delimited file.

    PARAMS
    ======
    d: pandas DataFrame
        DataFrame you want to save.
    output_path: string
        Output_path for the file. Checks if file is writable.
    columns: list, pandas.Index, tuple (default = d.columns)
        Columns in DataFrame to write. Default is all, in the order they appear.
    include_index: Boolean (default = False)
        Should the index be included as the first column? Default is no.
    index_label: String (default = "index")
        If include_index is True, this is the header for the index.
    naughty_characters: list (default = [np.inf, -np.inf])
        A list of elements that are replaced with rep_str. Note that all np.nan's (aka NaN's) are also replaced with
        rep_str.
    rep_str: String (default = "")
        The string that elements belonging to naughty_characters are replaced by.

    RETURNS
    =======
    output_path
    """

    filesnpaths.is_output_file_writable(output_path)

    if not columns:
        columns = d.columns

    d.replace(naughty_characters, np.nan, inplace=True)

    d.to_csv(output_path, sep="\t", columns=columns, index=include_index, index_label=index_label, na_rep=rep_str)
    return output_path


def store_dict_as_TAB_delimited_file(d, output_path, headers=None, file_obj=None):
    if not file_obj:
        filesnpaths.is_output_file_writable(output_path)

    if not file_obj:
        f = open(output_path, 'w')
    else:
        f = file_obj

    if not headers:
        headers = ['key'] + sorted(list(d.values())[0].keys())

    f.write('%s\n' % '\t'.join(headers))

    for k in sorted(d.keys()):
        line = [str(k)]
        for header in headers[1:]:
            try:
                val = d[k][header]
            except KeyError:
                raise ConfigError("Header ('%s') is not found in the dict :/" % (header))
            except TypeError:
                raise ConfigError("Your dictionary is not properly formatted to be exported\
                                    as a TAB-delimited file :/ You ask for '%s', but it is not\
                                    even a key in the dictionary" % (header))

            line.append(str(val) if not isinstance(val, type(None)) else '')

        f.write('%s\n' % '\t'.join(line))

    f.close()
    return output_path



def convert_numpy_array_to_binary_blob(array, compress=True):
    if compress:
        return gzip.compress(memoryview(array), compresslevel=1)
    else:
        return memoryview(array)


def convert_binary_blob_to_numpy_array(blob, dtype, decompress=True):
    if decompress:
        return np.frombuffer(gzip.decompress(blob), dtype=dtype)
    else:
        return np.frombuffer(blob, dtype=dtype)


def is_all_columns_present_in_TAB_delim_file(columns, file_path):
    columns = get_columns_of_TAB_delim_file(file_path)
    return False if len([False for c in columns if c not in columns]) else True


def HTMLColorToRGB(colorstring, scaled=True):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError("input #%s is not in #RRGGBB format" % colorstring)
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]

    if scaled:
        return (r / 255.0, g / 255.0, b / 255.0)
    else:
        return (r, g, b)


def transpose_tab_delimited_file(input_file_path, output_file_path):
    filesnpaths.is_file_exists(input_file_path)
    filesnpaths.is_file_tab_delimited(input_file_path)
    filesnpaths.is_output_file_writable(output_file_path)

    file_content = [line.strip('\n').split('\t') for line in open(input_file_path, 'rU').readlines()]

    output_file = open(output_file_path, 'w')
    for entry in zip(*file_content):
        output_file.write('\t'.join(entry) + '\n')
    output_file.close()

    return output_file_path


def get_random_colors_dict(keys):
    # FIXME: someone's gotta implement this
    # keys   : set(1, 2, 3, ..)
    # returns: {1: '#ffffff', 2: '#888888', 3: '#222222', ...}
    return dict([(k, None) for k in keys])


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


def get_column_data_from_TAB_delim_file(input_file_path, column_indices=[], expected_number_of_fields=None, separator='\t'):
    """Returns a dictionary where keys are the column indices, and items are the list of entries
    found in that that column"""
    filesnpaths.is_file_exists(input_file_path)
    filesnpaths.is_file_tab_delimited(input_file_path, expected_number_of_fields=expected_number_of_fields)

    d = {}

    for index in column_indices:
        d[index] = []

    with open(input_file_path, "rU") as input_file:
        for line in input_file.readlines():
            fields = line.strip('\n').split(separator)

            for index in column_indices:
                try:
                    d[index].append(fields[index])
                except:
                    raise ConfigError("get_column_data_from_TAB_delim_file is speaking: The file you sent\
                                        does not have data for the column index %d. Something is wrong :/" % (index))

    return d


def get_columns_of_TAB_delim_file(file_path, include_first_column=False):
    filesnpaths.is_file_exists(file_path)

    if include_first_column:
        return open(file_path, 'rU').readline().strip('\n').split('\t')
    else:
        return open(file_path, 'rU').readline().strip('\n').split('\t')[1:]


def get_names_order_from_newick_tree(newick_tree, newick_format=1, reverse=False):
    filesnpaths.is_proper_newick(newick_tree)

    tree = Tree(newick_tree, format=newick_format)

    names = [n.name for n in tree.get_leaves()]

    return list(reversed(names)) if reverse else names


def get_vectors_from_TAB_delim_matrix(file_path, cols_to_return=None, rows_to_return=[], transpose=False):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_tab_delimited(file_path)

    if transpose:
        transposed_file_path = filesnpaths.get_temp_file_path()
        transpose_tab_delimited_file(file_path, transposed_file_path)
        file_path = transposed_file_path

    rows_to_return = set(rows_to_return)
    vectors = []
    id_to_sample_dict = {}
    sample_to_id_dict = {}

    input_matrix = open(file_path, 'rU')
    columns = input_matrix.readline().strip().split('\t')[1:]

    fields_of_interest = []
    if cols_to_return:
        fields_of_interest = [columns.index(col) for col in cols_to_return]
    else:
        fields_of_interest = [f for f in range(0, len(columns)) if constants.IS_ESSENTIAL_FIELD(columns[f])]

    # update columns:
    columns = [columns[i] for i in fields_of_interest]

    if not len(columns):
        raise ConfigError("Only a subset (%d) of fields were requested by the caller, but none of them was found\
                            in the matrix (%s) :/" % (len(cols_to_return), file_path))

    id_counter = 0
    for line in input_matrix.readlines():
        row_name = line.strip().split('\t')[0]
        if rows_to_return and row_name not in rows_to_return:
                continue
        id_to_sample_dict[id_counter] = row_name
        fields = line.strip().split('\t')[1:]

        try:
            if fields_of_interest:
                vector = [float(fields[i]) for i in fields_of_interest]
            else:
                vector = [float(f) for f in fields]
        except ValueError:
            raise ConfigError("Matrix should contain only numerical values.")

        vectors.append(vector)

        id_counter += 1

    input_matrix.close()

    if transpose:
        # remove clutter
        os.remove(file_path)

    sample_to_id_dict = dict([(v, k) for k, v in id_to_sample_dict.items()])

    return id_to_sample_dict, sample_to_id_dict, columns, vectors


def apply_and_concat(df, fields, func, column_names, func_args=tuple([])):
    """ This function has been taken from https://tinyurl.com/y9ylqy4l
        and has been modified for speed considerations using this blog post:
        https://tinyurl.com/ya4e5tz3. Its utility is to append multiple columns to an existing
        dataframe row by row. This is usually a bad idea because usually operations can be
        vectorized. However when they cannot, looping through each row becomes a necessary evil.

        df: pandas DataFrame object
            An existing dataframe to loop through append columns to.
        fields: list
            A list of columns in the existing dataframe used to calculate the new columns
        func: function
            A function that takes as its first argument a row of `df` (i.e. a pd.Series
            object) and potential additional positional arguments `func_args`. It should return a
            tuple of values with with the same length as `column_names`.
        func_args: tuple
            A tuple of arguments passed to `func` besides the assumed first argument (a pd.Series
            object). For example, is `def func(row, a)`, then `func_args = (a,)`. If func_args is an
            empty tuple, `func` should take no other args.
        column_names: list
            A list of column headers for the newly appended columns
    """
    d = {column_name: [] for column_name in column_names}
    for _, row in df[fields].iterrows():
        out_values = func(row, *func_args)
        for ind, column_name in enumerate(column_names):
            d[column_name].append(out_values[ind])

    df2 = pd.DataFrame(d, index=df.index)
    return pd.concat((df, df2), axis=1, sort=True)


def get_values_of_gene_level_coverage_stats_as_dict(gene_level_coverage_stats_dict, key, genes_of_interest=None, samples_of_interest=None, as_pandas=False):
    """
        This function takes the gene_level_coverage_stats_dict and return one of the values
        as a matrix-like dict of dicts.
        THIS FUNCTION IS IN utils AND NOT IN summarizer, or dbops, because it used to be in summarizer
        and why should it be in summarizer?!? that makes no sense. And also mcg-classifier doesn't want
        to initialize summarizer, it wants to be able to just get the gene_level_coverage_stats_dict as
        input and then deal with it.

        There is also an option to as to get the data back as a pandas dataframe.
    """
    legal_keys = {'mean_coverage', 'detection', 'non_outlier_mean_coverage', 'non_outlier_coverage_std'}
    if key not in legal_keys and as_pandas:
        raise ConfigError("%s is not a valid key for creating a pandas dataframe of values of gene_level_coverage_stats_dict.\
                            Here is a list of the valid keys: %s" % (key, list(legal_keys)))

    gene_callers_ids = set(gene_level_coverage_stats_dict.keys())
    samples = set(next(iter(gene_level_coverage_stats_dict.values())).keys())

    if genes_of_interest is not None:
        missing_genes = [g for g in genes_of_interest if g not in gene_callers_ids]
        if len(missing_genes):
            raise ConfigError("The following genes are not in the gene_level_coverage_stats_dict, and yet you are asking for them: %s" % missing_genes)
    else:
        genes_of_interest = gene_callers_ids

    if samples_of_interest is not None:
        missing_samples = [s for s in samples_of_interest if s not in samples]
        if len(missing_samples):
            raise ConfigError("The following samples are not in the gene_level_coverage_stats_dict, and yet you are asking for them: %s" % missing_samples)
    else:
        samples_of_interest = samples

    d = {}

    for gene_callers_id in genes_of_interest:
        d[gene_callers_id] = {}
        for sample_name in samples_of_interest:
            d[gene_callers_id][sample_name] = gene_level_coverage_stats_dict[gene_callers_id][sample_name][key]

    if as_pandas:
        # This option is used by the mcg-classifier.
        import pandas as pd
        return pd.DataFrame.from_dict(d, orient='index')
    else:
        return d


def get_gene_caller_ids_from_args(gene_caller_ids, delimiter):
    gene_caller_ids_set = set([])
    if gene_caller_ids:
        if os.path.exists(gene_caller_ids):
            gene_caller_ids_set = set([g.strip() for g in open(gene_caller_ids, 'rU').readlines()])
        else:
            gene_caller_ids_set = set([g.strip() for g in gene_caller_ids.split(delimiter)])

    try:
        gene_caller_ids_set = set([int(g) for g in gene_caller_ids_set])
    except:
        g = gene_caller_ids_set.pop()
        raise ConfigError("The gene calls you provided do not look like gene callers anvi'o is used to working with :/ Here is\
                           one of them: '%s' (%s)." % (g, type(g)))
    return gene_caller_ids_set


def get_all_ids_from_fasta(input_file):
    fasta = u.SequenceSource(input_file)
    ids = []

    while next(fasta):
        ids.append(fasta.id)

    return ids


def get_ordinal_from_integer(num):
    """append 'st', 'nd', or 'th' to integer to make categorical. num must be integer"""
    return'%d%s' % (num, {11:'th', 12:'th', 13:'th'}.get(num%100, {1:'st', 2:'nd', 3:'rd'}.get(num%10,'th')))


def get_read_lengths_from_fasta(input_file):
    contig_lengths = {}

    fasta = u.SequenceSource(input_file)
    while next(fasta):
        contig_lengths[fasta.id] = len(fasta.seq)

    fasta.close()
    return contig_lengths


def get_GC_content_for_FASTA_entries(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    GC_content_dict = {}

    fasta = u.SequenceSource(file_path)
    while next(fasta):
        GC_content_dict[fasta.id] = get_GC_content_for_sequence(fasta.seq)

    return GC_content_dict


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
        raise ConfigError("get_synonymous_and_non_synonymous_potential :: sequence does not end \
                           with a stop codon and is therefore probably not what you want. If you \
                           want to continue anyways, use the just_do_it flag")

    synonymous_potential = 0
    for codon in list_of_codons_in_gene:
        for i, nt in enumerate(codon):
            for mutant_nt in [m for m in 'ACGT' if m != nt]:

                mutant_codon = list(codon)
                mutant_codon[i] = mutant_nt
                mutant_codon = ''.join(mutant_codon)

                if constants.codon_to_AA[mutant_codon] == constants.codon_to_AA[codon]:
                    synonymous_potential += 1/3

    non_synonymous_potential = 3 * len(list_of_codons_in_gene) - synonymous_potential
    return synonymous_potential, non_synonymous_potential


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
        raise ConfigError("utils::get_time_to_date is called with bad local_time.")

    return time.strftime(fmt, time.localtime(local_time))


def concatenate_files(dest_file, file_list):
    if not dest_file:
        raise ConfigError("Destination cannot be empty.")
    if not len(file_list):
        raise ConfigError("File list cannot be empty.")
    for f in file_list:
        filesnpaths.is_file_exists(f)
    filesnpaths.is_output_file_writable(dest_file)

    dest_file_obj = open(dest_file, 'w')
    for chunk_path in file_list:
        for line in open(chunk_path, 'rU'):
            dest_file_obj.write(line)

    dest_file_obj.close()
    return dest_file


def get_split_start_stops(contig_length, split_length, gene_start_stops=None):
    """Wrapper funciton for get_split_start_stops_with_gene_calls and get_split_start_stops_without_gene_calls"""
    if gene_start_stops:
        return get_split_start_stops_with_gene_calls(contig_length, split_length, gene_start_stops)
    else:
        return get_split_start_stops_without_gene_calls(contig_length, split_length)


def get_split_start_stops_with_gene_calls(contig_length, split_length, gene_start_stops):
    """Here we have a contig of `contig_length`, and a desired split length of `split_length`. also
       we know where genes start and stop in this contigs. we would like to split this contig into
       smaller pieces, i.e. sizes of `splits_length`, but in such a way that that splits do not
       break contigs in the middle of a gene."""

    # if the contig is too short, return it back.
    if contig_length < 2 * split_length:
        return [(0, contig_length)]

    non_coding_positions_in_contig = set(range(0, contig_length))

    # trim from the beginning and the end. we don't want to end up creating very short pieces
    non_coding_positions_in_contig = non_coding_positions_in_contig.difference(set(range(0, int(split_length / 2))))
    non_coding_positions_in_contig = non_coding_positions_in_contig.difference(set(range(contig_length - int(split_length / 2), contig_length)))

    # remove positions that code for genes
    for gene_unique_id, start, stop in gene_start_stops:
        start = start - 5
        stop = stop + 5
        non_coding_positions_in_contig = non_coding_positions_in_contig.difference(set(range(start, stop)))

    # what would be our break points in an ideal world? compute an initial list of break
    # points based on the length of the contig and desired split size:
    optimal_number_of_splits = int(contig_length / split_length)

    optimal_split_length = int(contig_length / optimal_number_of_splits)
    optimal_break_points = list(range(optimal_split_length, contig_length - optimal_split_length + 1, optimal_split_length))

    # now we will identify the very bad break points that we can't find a way to split around
    bad_break_points = set([])
    for i in range(0, len(optimal_break_points)):
        break_point = optimal_break_points[i]

        if break_point not in non_coding_positions_in_contig:
            # optimal break point hits a gene. we shall search towards both directions
            # to find a better break point:
            new_break_point = None
            for s in range(0, int(split_length / 2)):
                if break_point + s in non_coding_positions_in_contig:
                    new_break_point = break_point + s
                    break
                if break_point - s in non_coding_positions_in_contig:
                    new_break_point = break_point - s
                    break

            if not new_break_point:
                # nope. we failed. this is a bad bad break point.
                bad_break_points.add(break_point)
            else:
                # we are satisfied with the new one we found for now. it may be a shitty one,
                # but we will learn about that later. for now, let's replace the previous
                # optimal break pont with this one.
                optimal_break_points[i] = new_break_point

    # remove all the bad breakpoints from our 'optimal' break points:
    optimal_break_points = [p for p in optimal_break_points if p not in bad_break_points]

    if not len(optimal_break_points):
        # we have nothing left to work with after removal of the crappy break points. we will
        # keep this bad boy the way it is.
        return [(0, contig_length)]

    # create start/stop positions from these break points
    chunks = list(zip([0] + optimal_break_points[:-1], optimal_break_points)) + [(optimal_break_points[-1], contig_length)]

    return chunks


def get_split_start_stops_without_gene_calls(contig_length, split_length):
    """Returns split start stop locations for a given contig length."""
    num_chunks = int(contig_length / split_length)

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


def get_split_and_contig_names_of_interest(contigs_db_path, gene_caller_ids):
    """Takes a set of gene caller ids, returns all split and contig names in a
       contigs database that are affiliated with them.
    """

    if not isinstance(gene_caller_ids, set):
        raise ConfigError("`gene_caller_ids` must be of type `set`.")

    is_contigs_db(contigs_db_path)

    contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)

    where_clause_genes = "gene_callers_id in (%s)" % ', '.join(['%d' % g for g in gene_caller_ids])
    genes_in_contigs = contigs_db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause_genes)
    contig_names_of_interest = set([e['contig'] for e in genes_in_contigs.values()])

    where_clause_contigs = "parent in (%s)" % ', '.join(['"%s"' % c for c in contig_names_of_interest])
    splits_info = contigs_db.get_some_rows_from_table_as_dict(t.splits_info_table_name, where_clause=where_clause_contigs)
    split_names_of_ineterest = set(splits_info.keys())

    contigs_db.disconnect()

    return (split_names_of_ineterest, contig_names_of_interest)


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


def get_variabile_item_frequencies(e, engine='NT'):
    """
    e is a row from variable_nucleotide_positions table defined in tables.
    this function extends dictionary with consensus and departure from consensus.
    """

    items = constants.nucleotides if engine=='NT' else constants.amino_acids
    frequency_dict = Counter(dict([(item, e[item]) for item in items]))
    return frequency_dict.most_common()


def get_consensus_and_departure_data(variable_item_frequencies):
    """Make sense of `variable_item_frequencies`.

       The format of `variable_item_frequencies` follows this:

           >>> [('A', 45), ('T', 5), ('G', 0), ('N', 0), ('C', 0)]

       For a given entry of the variable_XX_frequencies table, the `variable_item_frequencies`
       tuple can be obtained via `get_variabile_item_frequencies`.
    """

    frequency_of_consensus = variable_item_frequencies[0][1]
    total_frequency_of_all_but_the_consensus = sum([tpl[1] for tpl in variable_item_frequencies[1:]])
    coverage = total_frequency_of_all_but_the_consensus + frequency_of_consensus

    n2n1ratio = variable_item_frequencies[1][1] / frequency_of_consensus if frequency_of_consensus else -1
    consensus = variable_item_frequencies[0][0]
    departure_from_consensus = total_frequency_of_all_but_the_consensus / coverage if coverage else -1

    return (n2n1ratio, consensus, departure_from_consensus)


def get_codon_order_to_nt_positions_dict(gene_call):
    """Returns a dictionary to translate codons in a gene to nucleotide positions"""

    if gene_call['partial']:
        raise ConfigError("get_codon_order_to_nt_positions_dict: this simply will not work\
                            for partial gene calls, and this on *is* a partial one.")

    start = gene_call['start']
    stop = gene_call['stop']

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


def convert_sequence_indexing(index, source="M0", destination="M1"):
    """
    Anvi'o zero-indexes sequences. For example, the methionine that every
    ORF starts with has the index 0 (M0). This is in contrast to the rest of the
    world, in which the methionine is indexed by 1 (M1). This function converts
    between the two.

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

def convert_SSM_to_single_accession(matrix_data):
    """
    The substitution scores from the SSM dictionaries created in anvio.data.SSMs are accessed via a dictionary of
    dictionaries, e.g.  data["Ala"]["Trp"]. This returns a new dictionary accessed via the concatenated sequence element
    pair, e.g. data["AlaTrp"], data["AT"], etc.  where they are ordered alphabetically.
    """
    items = matrix_data.keys()
    new_data = {}

    for row in items:
        for column in items:

            if row > column:
                continue
            new_data[''.join([row, column])] = matrix_data[row][column]
    return new_data


def get_DNA_sequence_translated(sequence, gene_callers_id, return_with_stops=False):
    sequence = sequence.upper()

    if len(sequence) % 3.0 != 0:
        raise ConfigError("The sequence corresponds to the gene callers id '%s' does not seem to\
                            have proper number of nucleotides to be translated :/ Here it is: %s" % (gene_callers_id, sequence))

    translated_sequence = ''

    for i in range(0, len(sequence), 3):
        single_letter_code = constants.AA_to_single_letter_code[constants.codon_to_AA[sequence[i:i + 3]]]

        if not single_letter_code:
            single_letter_code = 'X'

        translated_sequence += single_letter_code

    if translated_sequence.endswith('*'):
        if return_with_stops:
            pass
        else:
            translated_sequence = translated_sequence[:-1]

    return translated_sequence


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


def get_list_of_codons_for_gene_call(gene_call, contig_sequences_dict):
    codon_order_to_nt_positions = get_codon_order_to_nt_positions_dict(gene_call)

    if gene_call['contig'] not in contig_sequences_dict:
        raise ConfigError("get_list_of_AAs_for_gene_call: The contig sequences dict sent to\
                            this function does contain the contig name that appears in the gene call.\
                            Something is wrong here...")

    try:
        contig_sequence = contig_sequences_dict[gene_call['contig']]['sequence']
    except:
        raise ConfigError("get_list_of_AAs_for_gene_call: The contig sequences dict sent to\
                            this function does not seem to be an anvi'o contig sequences dict :/ It\
                            doesn't have the item 'sequence' in it.")

    list_of_codons = []
    for codon_order in codon_order_to_nt_positions:
        nt_positions = codon_order_to_nt_positions[codon_order]
        reference_codon_sequence = contig_sequence[nt_positions[0]:nt_positions[2] + 1]

        list_of_codons.append(constants.codon_to_codon_RC[reference_codon_sequence] if gene_call['direction'] == 'r' else reference_codon_sequence)

    return list_of_codons


def is_amino_acid_functionally_conserved(amino_acid_residue_1, amino_acid_residue_2):
    """Checks if two amino acid residues are part of the same biochemical property group"""
    group = constants.amino_acid_property_group[amino_acid_residue_1]
    conserved_group = constants.conserved_amino_acid_groups[group]

    if amino_acid_residue_2 in conserved_group:
        return True
    if group == 'Polar and Nonpolar': #they fall in more than one group, multiple tests needed
        if amino_acid_residue_1 == 'H' and (amino_acid_residue_2 in constants.conserved_amino_acid_groups['Nonpolar'] \
                                            or amino_acid_residue_2 in constants.conserved_amino_acid_groups['Bases']):
            return True
        if amino_acid_residue_1 == 'Y' and (amino_acid_residue_2 in constants.conserved_amino_acid_groups['Aromatic']):
            return True
    return False


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
        if parent in contig_name_to_splits_dict:
            contig_name_to_splits_dict[parent].append(split_name)
        else:
            contig_name_to_splits_dict[parent] = [split_name]

    return contig_name_to_splits_dict


def check_sample_id(sample_id):
    if sample_id:
        if sample_id[0] in constants.digits:
            raise ConfigError("Sample name ('%s') is not a valid name. Sample names can't start with digits.\
                                Long story. Please specify a sample name\
                                that starts with an ASCII letter (you may want to check '-s' parameter to set\
                                a sample name if your client permits (otherwise you are going to have to edit\
                                your input files))." % sample_id)

        allowed_chars_for_samples = constants.allowed_chars.replace('-', '').replace('.', '')
        if len([c for c in sample_id if c not in allowed_chars_for_samples]):
            raise ConfigError("Sample name ('%s') contains characters that anvio does not like. Please\
                                limit the characters that make up the project name to ASCII letters,\
                                digits, and the underscore character ('_')." % sample_id)


def is_this_name_OK_for_database(variable_name, content, stringent=True):
    if not content:
        raise ConfigError("But the %s is empty? Come on :(" % variable_name)

    if content[0] in constants.digits:
        raise ConfigError("Sorry, %s can't start with a digit. Long story. Please specify a name\
                            that starts with an ASCII letter." % variable_name)

    if stringent:
        allowed_chars = constants.allowed_chars.replace('.', '').replace('-', '')
    else:
        allowed_chars = constants.allowed_chars.replace('.', '')

    if len([c for c in content if c not in allowed_chars]):
        raise ConfigError("Well, the %s contains characters that anvi'o does not like :/ Please limit the characters\
                            to ASCII letters, digits, and the underscore ('_') character." % variable_name)


def check_contig_names(contig_names, dont_raise=False):
    all_characters_in_contig_names = set(''.join(contig_names))
    characters_anvio_doesnt_like = [c for c in all_characters_in_contig_names if c not in constants.allowed_chars]
    if len(characters_anvio_doesnt_like):
        if dont_raise:
            return False

        raise ConfigError("The name of at least one contig in your BAM file %s anvio does not\
                            like (%s). Please go back to your original files and make sure that\
                            the characters in contig names are limited to to ASCII letters,\
                            digits. Names can also contain underscore ('_'), dash ('-') and dot ('.')\
                            characters. anvio knows how much work this may require for you to go back and\
                            re-generate your BAM files and is very sorry for asking you to do that, however,\
                            it is critical for later steps in the analysis." \
                                % ("contains multiple characters" if len(characters_anvio_doesnt_like) > 1 else "contains a character",
                                   ", ".join(['"%s"' % c for c in characters_anvio_doesnt_like])))

    return True


def get_FASTA_file_as_dictionary(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    d = {}

    fasta = u.SequenceSource(file_path)
    while next(fasta):
        d[fasta.id] = fasta.seq

    return d


def unique_FASTA_file(input_file_path, output_fasta_path=None, names_file_path=None, store_frequencies_in_deflines=True):
    filesnpaths.is_file_exists(input_file_path)

    if not output_fasta_path:
        output_fasta_path = input_file_path + '.unique'

    if not names_file_path:
        names_file_path = output_fasta_path + '.names'

    if output_fasta_path == names_file_path:
        raise ConfigError("I can't unique this. Output FASTA file path can't be identical to\
                            the names file path...")

    if output_fasta_path == input_file_path or names_file_path == input_file_path:
        raise ConfigError("Anvi'o will not unique this. Output FASTA path and names file path should\
                            be different from the the input file path...")

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


def ununique_BLAST_tabular_output(tabular_output_path, names_dict):
    new_search_output_path = tabular_output_path + '.ununiqued'
    new_tabular_output = open(new_search_output_path, 'w')

    for line in open(tabular_output_path):
        fields = line.strip().split('\t')
        for query_id in names_dict[fields[0]]:
            for subject_id in names_dict[fields[1]]:
                new_tabular_output.write('%s\t%s\t%s\n' % (query_id, subject_id, '\t'.join(fields[2:])))


    new_tabular_output.close()

    shutil.move(tabular_output_path, tabular_output_path + '.unique')
    shutil.move(new_search_output_path, tabular_output_path)

    return tabular_output_path


def get_BLAST_tabular_output_as_dict(tabular_output_path, target_id_parser_func=None, query_id_parser_func=None):
    """Takes a BLAST output, returns a dict where each query appears only once!!

       If there are multiple hits for a given query, the one with lower e-value.
       remains in the dict.
    """

    results_dict = {}

    for line in open(tabular_output_path):
        fields = line.strip().split('\t')
        query_id = fields[0] if not query_id_parser_func else query_id_parser_func(fields[0])
        target_id = fields[1] if not target_id_parser_func else target_id_parser_func(fields[1])
        e_value = float(fields[10])

        if query_id in results_dict:
            if e_value > results_dict[query_id]['evalue']:
                continue

        results_dict[query_id] = {'hit': target_id, 'evalue': e_value}

    return results_dict


def store_dict_as_FASTA_file(d, output_file_path, wrap_from=200):
    filesnpaths.is_output_file_writable(output_file_path)
    output = open(output_file_path, 'w')

    for key in d:
        output.write('>%s\n' % key)
        output.write('%s\n' % textwrap.fill(d[key], wrap_from, break_on_hyphens=False))

    output.close()
    return True


def export_sequences_from_contigs_db(contigs_db_path, output_file_path, seq_names_to_export=None, splits_mode=False, rna_alphabet=False, truncate=True):
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

    for seq_name in seq_names_to_export:
        if splits_mode:
            s = splits_info_dict[seq_name]
            sequence = FORMAT(contig_sequences_dict[s['parent']]['sequence'][s['start']:s['end']])
        else:
            sequence = FORMAT(contig_sequences_dict[seq_name]['sequence'])

        output_fasta.write_id(seq_name)
        output_fasta.write_seq(sequence, split=truncate)

    return True


def gen_gexf_network_file(units, samples_dict, output_file, sample_mapping_dict=None,
                               unit_mapping_dict=None, project=None, sample_size=8, unit_size=2,
                               skip_sample_labels=False, skip_unit_labels=False):
    """A function that generates an XML network description file for Gephi.

       Two minimum required inputs are `units`, and `samples_dict`.

       Simply, `samples_dict` is a dictionary that shows the distribution of `units` and their
       frequencies across samples. Here is an example `units` variable (which is a type of `list`):

            units = ['unit_1', 'unit_2', ... 'unit_n']

       and a corresponding `samples_dict` would look like this:

            samples_dict = {'sample_1': {'unit_1': 0.5,
                                        'unit_2': 0.2,
                                         ...,
                                         'unit_n': 0.1
                                        },
                            'sample_2': { (...)
                                            },
                            (...),
                            'sample_n': { (...)
                                            }
                        }
    """

    filesnpaths.is_output_file_writable(output_file)

    output = open(output_file, 'w')

    samples = sorted(samples_dict.keys())
    sample_mapping_categories = sorted([k for k in list(sample_mapping_dict.values())[0].keys() if k != 'colors']) if sample_mapping_dict else None
    unit_mapping_categories = sorted([k for k in list(unit_mapping_dict.keys()) if k not in ['colors', 'labels']]) if unit_mapping_dict else None

    sample_mapping_category_types = []
    if sample_mapping_dict:
        for category in sample_mapping_categories:
            if RepresentsFloat(list(sample_mapping_dict.values())[0][category]):
                sample_mapping_category_types.append('double')
            else:
                sample_mapping_category_types.append('string')

    output.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
    output.write('''<gexf xmlns:viz="http:///www.gexf.net/1.1draft/viz" xmlns="http://www.gexf.net/1.2draft" version="1.2">\n''')
    output.write('''<meta lastmodifieddate="2010-01-01+23:42">\n''')
    output.write('''    <creator>Oligotyping pipeline</creator>\n''')
    if project:
        output.write('''    <creator>Network description for %s</creator>\n''' % (project))
    output.write('''</meta>\n''')
    output.write('''<graph type="static" defaultedgetype="undirected">\n\n''')

    if sample_mapping_dict:
        output.write('''<attributes class="node" type="static">\n''')
        for i in range(0, len(sample_mapping_categories)):
            category = sample_mapping_categories[i]
            category_type = sample_mapping_category_types[i]
            output.write('''    <attribute id="%d" title="%s" type="%s" />\n''' % (i, category, category_type))
        output.write('''</attributes>\n\n''')

    # FIXME: IDK what the hell is this one about:
    if unit_mapping_dict:
        output.write('''<attributes class="edge">\n''')
        for i in range(0, len(unit_mapping_categories)):
            category = unit_mapping_categories[i]
            output.write('''    <attribute id="%d" title="%s" type="string" />\n''' % (i, category))
        output.write('''</attributes>\n\n''')

    output.write('''<nodes>\n''')
    for sample in samples:
        if skip_sample_labels:
            output.write('''    <node id="%s">\n''' % (sample))
        else:
            output.write('''    <node id="%s" label="%s">\n''' % (sample, sample))

        output.write('''        <viz:size value="%d"/>\n''' % sample_size)

        if sample_mapping_dict and 'colors' in sample_mapping_dict[sample]:
            output.write('''        <viz:color r="%d" g="%d" b="%d" a="1"/>\n''' %\
                                             HTMLColorToRGB(sample_mapping_dict[sample]['colors'], scaled=False))

        if sample_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(sample_mapping_categories)):
                category = sample_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, sample_mapping_dict[sample][category]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    for unit in units:
        if skip_unit_labels:
            output.write('''    <node id="%s">\n''' % (unit))
        else:
            if unit_mapping_dict and 'labels' in unit_mapping_dict:
                output.write('''    <node id="%s" label="%s">\n''' % (unit, unit_mapping_dict['labels'][unit]))
            else:
                output.write('''    <node id="%s">\n''' % (unit))
        output.write('''        <viz:size value="%d" />\n''' % unit_size)

        if unit_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(unit_mapping_categories)):
                category = unit_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, unit_mapping_dict[category][unit]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    output.write('''</nodes>\n''')

    edge_id = 0
    output.write('''<edges>\n''')
    for sample in samples:
        for i in range(0, len(units)):
            unit = units[i]
            if samples_dict[sample][unit] > 0.0:
                if unit_mapping_dict:
                    output.write('''    <edge id="%d" source="%s" target="%s" weight="%f">\n''' % (edge_id, unit, sample, samples_dict[sample][unit]))
                    if unit_mapping_categories:
                        output.write('''        <attvalues>\n''')
                        for i in range(0, len(unit_mapping_categories)):
                            category = unit_mapping_categories[i]
                            output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, unit_mapping_dict[category][unit]))
                        output.write('''        </attvalues>\n''')
                    output.write('''    </edge>\n''')
                else:
                    output.write('''    <edge id="%d" source="%s" target="%s" weight="%f" />\n''' % (edge_id, unit, sample, samples_dict[sample][unit]))


                edge_id += 1
    output.write('''</edges>\n''')
    output.write('''</graph>\n''')
    output.write('''</gexf>\n''')

    output.close()


def is_ascii_only(text):
    """test whether 'text' is composed of ASCII characters only"""
    return all(ord(c) < 128 for c in text)


def get_TAB_delimited_file_as_dictionary(file_path, expected_fields=None, dict_to_append=None, column_names=None,\
                                        column_mapping=None, indexing_field=0, separator='\t', no_header=False,\
                                        ascii_only=False, only_expected_fields=False, assign_none_for_missing=False,\
                                        none_value=None, empty_header_columns_are_OK=False):
    """Takes a file path, returns a dictionary."""

    if expected_fields and (not isinstance(expected_fields, list) and not isinstance(expected_fields, set)):
        raise ConfigError("'expected_fields' variable must be a list (or a set).")

    if only_expected_fields and not expected_fields:
        raise ConfigError("'only_expected_fields' variable guarantees that there are no more fields present\
                            in the input file but the ones requested with 'expected_fields' variable. If you\
                            need to use this flag, you must also be explicit about what fields you expect to\
                            find in the file.")

    filesnpaths.is_file_plain_text(file_path)
    filesnpaths.is_file_tab_delimited(file_path, separator=separator)

    f = open(file_path, 'rU')

    # learn the number of fields and reset the file:
    num_fields = len(f.readline().strip('\n').split(separator))
    f.seek(0)

    # if there is no file header, make up a columns list:
    if no_header and not column_names:
        column_names = ['column_%05d' % i for i in range(0, num_fields)]

    if column_names:
        columns = column_names

        if num_fields != len(columns):
            raise  ConfigError("Number of column names declared (%d) differs from the number of columns\
                                 found (%d) in the matrix ('%s') :/" % (len(columns), num_fields, file_path))

        # now we set the column names. if the file had its header, we must discard
        # the first line. so here we go:
        if not no_header:
            f.readline()
    else:
        columns = f.readline().strip('\n').split(separator)

    if not empty_header_columns_are_OK and min(map(len, columns)) == 0:
        raise ConfigError("At least one of the column headers in your tab delimited file '%s'\
                           is empty." % file_path)

    if expected_fields:
        for field in expected_fields:
            if field not in columns:
                raise ConfigError("The file '%s' does not contain the right type of header. It was expected\
                                    to have these: '%s', however it had these: '%s'" % (file_path,
                                                                                        ', '.join(expected_fields),
                                                                                        ', '.join(columns[1:])))

    if only_expected_fields:
        for field in columns:
            if field not in expected_fields:
                raise ConfigError("There are more fields in the file '%s' than the expected fields :/\
                                    Anvi'o is telling you about this because get_TAB_delimited_file_as_dictionary\
                                    funciton is called with `only_expected_fields` flag turned on.")

    d = {}
    line_counter = 0

    for line in f.readlines():
        if ascii_only:
            if not is_ascii_only(line):
                raise ConfigError("The input file conitans non-ascii characters at line number %d. Those lines\
                                    either should be removed, or edited." % (line_counter + 2))

        line_fields = [f if f else None for f in line.strip('\n').split(separator)]

        if line_fields and line_fields[0] == None:
            raise ConfigError("The line number %d in '%s' has no data in its first column, and this doesn't\
                               seem right at all :/" % (line_counter + 1, file_path))

        if column_mapping:
            updated_line_fields = []
            for i in range(0, len(line_fields)):
                try:
                    if line_fields[i] == None and column_mapping[i] in [float, int]:
                        updated_line_fields.append(column_mapping[i](0))
                    else:
                        updated_line_fields.append(column_mapping[i](line_fields[i]))
                except NameError:
                    raise ConfigError("Mapping function '%s' did not work on value '%s'. These functions can be native\
                                        Python functions, such as 'str', 'int', or 'float', or anonymous functions\
                                        defined using lambda notation." % (column_mapping[i], line_fields[i]))
                except TypeError:
                    raise ConfigError("Mapping function '%s' does not seem to be a proper Python function :/" % column_mapping[i])
                except ValueError:
                    raise ConfigError("Mapping funciton '%s' did not like the value '%s' in column number %d\
                                        of the input matrix '%s' :/" % (column_mapping[i], line_fields[i], i + 1, file_path))
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
        keys = list(d.values())[0].keys()

        for entry in dict_to_append:
            if entry not in d:
                # so dict to append is missing a key that is in the dict to be appended. if the user did not
                # ask us to add None for these entries via none_for_missing, we are going to make a noise,
                # otherwise we will tolerate it.
                if not assign_none_for_missing:
                    raise ConfigError("Appending entries to the already existing dictionary from file '%s' failed\
                                        as the entry %s does not appear to be in the file." % (file_path, entry))
                else:
                    for key in keys:
                        dict_to_append[entry][key] = none_value
            else:
                for key in keys:
                    dict_to_append[entry][key] = d[entry][key]

        return dict_to_append

    return d


def get_filtered_dict(input_dict, item, accepted_values_set):
    # removes any entry from d, where the value of the 'item' of items in d does not match
    # with 'accepted_values'
    if not isinstance(accepted_values_set, type(set([]))):
        raise ConfigError("get_filtered_dict: values must be type of set([]).")

    filtered_dict = {}

    for entry_id in input_dict:
        if input_dict[entry_id][item] not in accepted_values_set:
            continue
        else:
            filtered_dict[entry_id] = input_dict[entry_id]

    return filtered_dict


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
        raise ConfigError("HMM stuff is upset with you. There are unexpected number of fields in the target\
                            file.")

    if alphabet not in ['AA', 'DNA', 'RNA']:
        raise ConfigError("The alphabet in the target file (%s) isnot one of the alphabets anvi'o knows how to\
                            work with. Here is a list for you to choose from: 'DNA', 'RNA', or 'AA'" % alphabet)

    if context not in ['GENE', 'CONTIG', None]:
        raise ConfigError("The context you defined in the target file (%s) does not make any sense to anvi'o.\
                            It would have, if you had chosen one of these: 'GENE', 'CONTIG'." % context)

    if alphabet == 'AA' and context == 'CONTIG':
        raise ConfigError("You can't use the AA alphabet with the CONTIGS context :/ You need to set your target\
                            again. 'AA' or 'AA:GENE' would have worked much better.")

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

    R = lambda f: open(os.path.join(source, f), 'rU').readlines()[0].strip()
    for source in source_dirs:
        if source.endswith('/'):
            source = source[:-1]

        if not PROPER(os.path.basename(source)):
            raise ConfigError("One of the search database directories ('%s') contains characters in its name\
                                anvio does not like. Directory names should be at least three characters long\
                                and must not contain any characters but ASCII letters, digits and\
                                underscore" % os.path.basename(source))

        expected_files = ['reference.txt', 'kind.txt', 'genes.txt', 'genes.hmm.gz', 'target.txt', 'noise_cutoff_terms.txt']

        missing_files = [f for f in expected_files if not os.path.exists(os.path.join(source, f))]
        if missing_files:
            raise ConfigError("Each search database directory must contain following files: %s'. Yet, the HMM source '%s' seems to\
                               be missing the follwoing one(s): %s. See this blog post to make sure you are doing it the way it\
                               should be done: http://merenlab.org/2016/05/21/archaeal-single-copy-genes/" % \
                                            (', '.join(expected_files), os.path.basename(source), ', '.join(missing_files)))

        empty_files = [f for f in expected_files if os.stat(os.path.join(source, f)).st_size == 0]
        if empty_files:
            raise ConfigError("One or more files for the HMM source '%s' seems to be empty. Which creates lots of\
                               counfusion around these parts of the code. Anvi'o could set some defualts for you,\
                               but it would be much better if you set your own defaults explicitly. You're not\
                               sure what would make a good default for your HMM collection? Reach out to\
                               a developer, and they will help you! Here are the files that are empty: %s." % \
                                    (os.path.basename(source), ', '.join(empty_files)))

        ref = R('reference.txt')
        kind = R('kind.txt')
        target = R('target.txt')
        noise_cutoff_terms = R('noise_cutoff_terms.txt')
        anvio_hmm_target_term_to_alphabet_and_context(target)

        domain = None
        if kind == 'singlecopy' and kind.count(':') == 0:
            raise ConfigError("This HMM profile seems to be a collection of single-copy core genes. Great. But for\
                               this kind, you must also declare a 'domain' in your 'kind.txt' file. It is simple.\
                               For instance, you could use 'singlecopy:bacteria', or 'singlecopy:archaea', or\
                               'singlecopy:myspecificbranch'.")
        if kind.count(':') == 1:
            kind, domain = kind.split(':')

        if not PROPER(kind):
            raise ConfigError("'kind.txt' defines the kind of search this database offers. The kind term must be a single\
                                word that is at least three characters long, and must not contain any characters but\
                                ASCII letters, digits, and underscore. Here are some nice examples: 'singlecopy',\
                                or 'pathogenicity', or 'noras_selection'. But yours is '%s'." % (kind))

        if domain and not PROPER(domain):
            raise ConfigError("That's lovely that you decided to specify a domain extension for your HMM collection in the\
                                'kind.txt'. Although, your domain term is not a good one, as it must be a single\
                                word that is at least three characters long, and without any characters but\
                                ASCII letters, digits, and underscore. Confused? That's fine. Send an e-mail to the anvi'o\
                                developers, and they will help you!")

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


def check_misc_data_keys_for_format(data_keys_list):
    """A function to make sure user-provided misc data keys are compatible
       with the current version of anvi'o. Housekeeping BS."""

    if not data_keys_list:
        return

    # findout whether the user data contains the older implementation of stacked
    # bar data type
    obsolete_stackedbar_keys = [k for k in data_keys_list if k.find('!') > -1 and k.find(';') > -1]

    if len(obsolete_stackedbar_keys):
        key_violates_new_rule = obsolete_stackedbar_keys[0]
        main_key, data_items = key_violates_new_rule.split('!')
        new_rule_compatible_data_keys = ['%s!%s' % (main_key, d) for d in data_items.split(';')]

        raise ConfigError("Oh no :( We recently changed the description of the stacked bar data type, and your input data\
                           file still has the older version. Here is the list of those that are violating the new format:\
                           %s. To avoid this issue and to turn them into the new format, you could take '%s', and present\
                           it as %d separate TAB-delimited entries that look like this: %s. Sorry!" % \
                                            (', '.join(['"%s"' % k for k in obsolete_stackedbar_keys]),
                                             key_violates_new_rule,
                                             len(new_rule_compatible_data_keys),
                                             ', '.join(['"%s"' % k for k in new_rule_compatible_data_keys])))


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
        raise ConfigError("Accession IDs in your HMM model should be unique, however, the `genes.hmm.gz`\
                           file for `%s` seems to have the same accession ID (the line that starts with `ACC`)\
                           more than once :(" % (os.path.abspath(model_path).split('/')[-2]))

    if len(genes.difference(genes_in_model)):
        raise ConfigError("Some gene names in genes.txt file does not seem to be appear in genes.hmm.gz.\
                           Here is a list of missing gene names: %s" % ', '.join(list(genes.difference(genes_in_model))))

    if len(genes_in_model.difference(genes)):
        raise ConfigError("Some gene names in genes.hmm.gz file does not seem to be appear in genes.txt.\
                           Here is a list of missing gene names: %s" % ', '.join(list(genes_in_model.difference(genes))))


def get_missing_programs_for_hmm_analysis():
    missing_programs = []
    for p in ['prodigal', 'hmmscan']:
        try:
            is_program_exists(p)
        except ConfigError:
            missing_programs.append(p)
    return missing_programs


def get_db_type(db_path):
    filesnpaths.is_file_exists(db_path)
    database = db.DB(db_path, None, ignore_version=True)

    tables = database.get_table_names()
    if 'self' not in tables:
        database.disconnect()
        raise ConfigError("'%s' does not seem to be a anvi'o database..." % db_path)

    db_type = database.get_meta_value('db_type')
    database.disconnect()

    return db_type


def get_required_version_for_db(db_path):
    db_type = get_db_type(db_path)

    if db_type not in t.versions_for_db_types:
        raise ConfigError("Anvi'o was trying to get the version of the -alleged- anvi'o database '%s', but it failed\
                            because it turns out it doesn't know anything about this '%s' type." % (db_path, db_type))

    return t.versions_for_db_types[db_type]


def get_all_sample_names_from_the_database(db_path):
    """Returns all 'sample' names from a given database. At least it tries."""

    db_type = get_db_type(db_path)
    database = db.DB(db_path, get_required_version_for_db(db_path))

    if db_type == 'profile':
        samples = []
        try:
            samples = [s.strip() for s in database.get_meta_value('samples').split(',')]
        except:
            pass

        return set(samples)

    elif db_type == 'pan':
        internal_genome_names, external_genome_names = [], []
        try:
            internal_genome_names = [g.strip() for g in database.get_meta_value('internal_genome_names').split(',')]
        except:
            pass

        try:
            external_genome_names = [g.strip() for g in database.get_meta_value('external_genome_names').split(',')]
        except:
            pass

        return set([s for s in internal_genome_names + external_genome_names if s])

    else:
        raise ConfigError("`get_all_sample_names_from_the_database` function does not know how to deal\
                            with %s databases." % db_type)


def get_all_item_names_from_the_database(db_path, run=run):
    """Return all split names or gene cluster names in a given database"""

    all_items = set([])

    database = db.DB(db_path, get_required_version_for_db(db_path))
    db_type = database.get_meta_value('db_type')

    if db_type == 'profile':
        if is_blank_profile(db_path):
            run.warning("Someone asked for the split names in a blank profile database. Sadly, anvi'o does not keep track\
                         of split names in blank profile databases. This function will return an empty set as split names\
                         to not kill your mojo, but whatever you were trying to do will not work :(")
            return set([])
        elif int(database.get_meta_value('merged')):
            all_items = set(database.get_single_column_from_table('mean_coverage_Q2Q3_splits', 'contig'))
        else:
            all_items = set(database.get_single_column_from_table('atomic_data_splits', 'contig'))
    elif db_type == 'pan':
        all_items = set(database.get_single_column_from_table(t.pan_gene_clusters_table_name, 'gene_cluster_id'))
    elif db_type == 'contigs':
        all_items = set(database.get_single_column_from_table(t.splits_info_table_name, 'split'))
    else:
        database.disconnect()
        raise ConfigError("You wanted to get all items in the database %s, but no one here knows about its type. Seriously,\
                            what is '%s' anyway?" % (db_path, db_type))

    if not len(all_items):
        database.disconnect()
        raise ConfigError("utils::get_all_item_names_from_the_database speaking. Something that should never happen happened :/\
                           There seems to be nothing in this %s database. Anvi'o is as confused as you are. Please get in touch\
                           with a developer. They will love this story." % db_path)

    database.disconnect()

    return all_items


def get_variability_table_engine_type(table_path, dont_raise=False):
    """A non-extensive test to determine if a file was generated by anvi-gen-variability-profile,
       and if it was, what engine (NT, CDN, or AA) was used.
    """
    filesnpaths.is_file_tab_delimited(table_path)
    columns_names = set(pd.read_csv(table_path, sep="\t", nrows = 0).columns)

    if set(constants.nucleotides) < columns_names:
        return "NT"

    elif set(constants.codons) < columns_names:
        return "CDN"

    elif set(constants.amino_acids) < columns_names:
        return "AA"

    else:
        if dont_raise:
            return ""
        raise ConfigError("anvi'o does not recognize %s as being a variability table generated by\
                           anvi-gen-variability-profile." % table_path)


def is_contigs_db(db_path):
    filesnpaths.is_file_exists(db_path)
    if get_db_type(db_path) != 'contigs':
        raise ConfigError("'%s' is not an anvi'o contigs database." % db_path)
    return True


def is_pan_or_profile_db(db_path):
    if get_db_type(db_path) not in ['pan', 'profile']:
        raise ConfigError("'%s' is neither a pan nor a profile database :/ Someone is in trouble." % db_path)
    return True


def is_profile_db(db_path):
    if get_db_type(db_path) != 'profile':
        raise ConfigError("'%s' is not an anvi'o profile database." % db_path)
    return True


def is_structure_db(db_path):
    if get_db_type(db_path) != 'structure':
        raise ConfigError("'%s' is not an anvi'o structure database." % db_path)
    return True


def is_blank_profile(db_path):
    if get_db_type(db_path) != 'profile':
        return False

    database = db.DB(db_path, None, ignore_version=True)
    blank = int(database.get_meta_value('blank'))
    database.disconnect()

    return True if blank == 1 else False


def get_two_sample_z_test_statistic(p1, p2, n1, n2):
    '''
        Compute a two sample z-test statistic

        If one group has no hits (e.g. p1=0) then we compute an upper bound
        for the p-value by pretending that it had one hit.

        If one group has 100% hits (e.g. p1=1) then we compute an upper bound
        for the p-value by pretending that p1=1-1/n1 hits
    '''
    import numpy
    if p1 == 0 and p2 == 0:
        return (0, 0)

    # This is done in order to estimate an upper bound
    # for the p-value
    p1 = max(p1, 1/n1) # in case p1 is zero
    p2 = max(p2, 1/n2)
    p1 = min(p1, 1 - 1/n1) # in case p1 is 1
    p2 = min(p2, 1 - 1/n2)

    p = (n1*p1 + n2*p2) / (n1 + n2)

    z = (p1 - p2) / numpy.sqrt(p*(1 - p) * (1/n1 + 1/n2))
    p_value = get_p_value_for_z_test(z)
    return (z, p_value)


def get_p_value_for_z_test(z):
    from scipy.stats import norm
    return 2*norm.cdf(-abs(z))


def is_pan_db(db_path):
    if get_db_type(db_path) != 'pan':
        raise ConfigError("'%s' is not an anvi'o pan database." % db_path)
    return True


def is_profile_db_merged(profile_db_path):
    is_profile_db(profile_db_path)

    profile_db = db.DB(profile_db_path, get_required_version_for_db(profile_db_path))
    merged = int(profile_db.get_meta_value('merged'))
    profile_db.disconnect()

    return merged


def is_profile_db_and_contigs_db_compatible(profile_db_path, contigs_db_path):
    is_profile_db(profile_db_path)
    is_contigs_db(contigs_db_path)

    profile_db = db.DB(profile_db_path, get_required_version_for_db(profile_db_path))
    contigs_db = db.DB(contigs_db_path, get_required_version_for_db(contigs_db_path))

    p_hash = profile_db.get_meta_value('contigs_db_hash')
    a_hash = contigs_db.get_meta_value('contigs_db_hash')
    merged = int(profile_db.get_meta_value('merged'))

    profile_db.disconnect()
    contigs_db.disconnect()

    if a_hash != p_hash:
        raise ConfigError('The contigs database and the profile database does not\
                           seem to be compatible. More specifically, this contigs\
                           database is not the one that was used when %s generated\
                           this profile database (%s != %s).'\
                               % ('anvi-merge' if merged else 'anvi-profile', a_hash, p_hash))

    return True


def is_structure_db_and_contigs_db_compatible(structure_db_path, contigs_db_path):
    is_structure_db(structure_db_path)
    is_contigs_db(contigs_db_path)

    structure_db = db.DB(structure_db_path, get_required_version_for_db(structure_db_path))
    contigs_db = db.DB(contigs_db_path, get_required_version_for_db(contigs_db_path))

    p_hash = structure_db.get_meta_value('contigs_db_hash')
    a_hash = contigs_db.get_meta_value('contigs_db_hash')

    structure_db.disconnect()
    contigs_db.disconnect()

    if a_hash != p_hash:
        raise ConfigError('The contigs and structure databases do not seem compatible.\
                           More specifically, the contigs database is not the one that\
                           was used when the structure database was created (%s != %s).'\
                               % (a_hash, p_hash))

    return True


def download_file(url, output_file_path, progress=progress, run=run):
    filesnpaths.is_output_file_writable(output_file_path)

    try:
        response = urllib.request.urlopen(url)
    except Exception as e:
        raise ConfigError("Something went wrong with your download attempt. Here is the\
                            problem: '%s'" % e)

    file_size = int(response.headers['Content-Length'])

    f = open(output_file_path, 'wb')

    progress.new('Downloading "%s"' % os.path.basename(output_file_path))
    progress.update('...')

    downloaded_size = 0
    while True:
        buffer = response.read(10000)

        if buffer:
            downloaded_size += len(buffer)
            f.write(buffer)
            progress.update('%.1f%%' % (downloaded_size * 100.0 / file_size))
        else:
            break

    f.close()

    progress.end()
    run.info('Downloaded succesfully', output_file_path)


def download_protein_structures(protein_code_list, output_dir):
    """
    Downloads protein structures using Biopython. protein_code_list is a list
    of 4-letter protein codes. Returns list of successful downloads
    """
    import Bio.PDB as PDB

    progress.new("Downloading proteins from PDB")

    filesnpaths.gen_output_directory(output_dir)

    pdb_list = PDB.PDBList()

    # this rule may one day change
    get_protein_path = lambda x: os.path.join(output_dir, "pdb" + x + ".ent")

    for protein_code in protein_code_list:
        progress.update("Downloading protein structure: {}".format(protein_code))

        with SuppressAllOutput(): # FIXME SuppressAllOutput gives error
            pdb_list.retrieve_pdb_file(protein_code, file_format="pdb", pdir=output_dir, overwrite=True)

        # raise warning if structure was not downloaded
        if not filesnpaths.is_file_exists(get_protein_path(protein_code), dont_raise=True):
            run.warning("The protein {} could not be downloaded. Are you connected to internet?".format(protein_code))
            protein_code_list.remove(protein_code)

    progress.end()
    return protein_code_list


def get_file_md5(file_path):
    hash_md5 = hashlib.md5()

    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    return hash_md5.hexdigest()


def run_selenium_and_export_svg(url, output_file_path, browser_path=None, run=run):
    if filesnpaths.is_file_exists(output_file_path, dont_raise=True):
        raise FilesNPathsError("The output file already exists. Anvi'o does not like overwriting stuff.")

    filesnpaths.is_output_file_writable(output_file_path)

    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        from selenium.common.exceptions import TimeoutException
    except:
        raise ConfigError("You want to export SVGs? Well, you need the Python library 'selenium' to be able to\
                           do that but you don't have it. If you are lucky, you probably can install it by\
                           typing 'pip install selenium' or something :/")

    if browser_path:
        filesnpaths.is_file_exists(browser_path)
        run.info_single('You are launching an alternative browser. Keep an eye on things!', mc='red', nl_before=1)
        driver = webdriver.Chrome(executable_path=browser_path)
    else:
        driver = webdriver.Chrome()

    driver.wait = WebDriverWait(driver, 10)
    driver.set_window_size(1920, 1080)
    driver.get(url)

    try:
        WebDriverWait(driver, 300).until(EC.text_to_be_present_in_element((By.ID, "title-panel-second-line"), "Current view"))
    except TimeoutException:
        print("Timeout occured, could not get the SVG drawing in 600 seconds.")
        driver.quit()
    time.sleep(1)

    driver.execute_script("exportSvg(true);")
    time.sleep(1)

    svg = driver.find_element_by_id('panel-center')

    svg_file = open(output_file_path, 'w')
    svg_file.write(svg.get_attribute('innerHTML'))
    svg_file.close()
    driver.quit()

    run.info_single('\'%s\' saved successfully.' % output_file_path)


def open_url_in_browser(url, browser_path=None, run=run):
    if browser_path:
        filesnpaths.is_file_exists(browser_path)
        run.info_single('You are launching an alternative browser. Keep an eye on things!', mc='red', nl_before=1)
        webbrowser.register('users_preferred_browser', None, webbrowser.BackgroundBrowser(browser_path))
        webbrowser.get('users_preferred_browser').open_new(url)
    else:
        webbrowser.open_new(url)


def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def RepresentsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class Mailer:
    def __init__(self, from_address='admin@localhost', server_address='localhost', server_port=25,
                 init_tls=False, username=None, password=None, run=Run(verbose=False),
                 progress=Progress(verbose=False)):
        self.from_address = from_address
        self.server_address = server_address
        self.server_port = server_port
        self.init_tls = init_tls
        self.username = username
        self.password = password

        self.server = None
        self.config_ini_path = None

        self.run = run
        self.progress = progress

        self.config_template = {
                'SMTP': {
                        'from_address': {'mandatory': True, 'test': lambda x: str(x)},
                        'server_address': {'mandatory': True, 'test': lambda x: str(x)},
                        'server_port': {'mandatory': True, 'test': lambda x: RepresentsInt(x) and int(x) > 0, 'required': 'an integer'},
                        'init_tls': {'mandatory': True, 'test': lambda x: x in ['True', 'False'], 'required': 'True or False'},
                        'username': {'mandatory': True, 'test': lambda x: str(x)},
                        'password': {'mandatory': True, 'test': lambda x: str(x)},
                    },
            }


    def init_from_config(self, config_ini_path):
        def get_option(self, config, section, option, cast):
            try:
                return cast(config.get(section, option).strip())
            except configparser.NoOptionError:
                return None

        filesnpaths.is_file_exists(config_ini_path)

        self.config_ini_path = config_ini_path

        config = configparser.ConfigParser()

        try:
            config.read(self.config_ini_path)
        except Exception as e:
            raise ConfigError("Well, the file '%s' does not seem to be a config file at all :/ Here\
                                is what the parser had to complain about it: %s" % (self.config_ini_path, e))

        section = 'SMTP'

        if section not in config.sections():
            raise ConfigError("The config file '%s' does not seem to have an 'SMTP' section, which\
                                is essential for Mailer class to learn server and authentication\
                                settings. Please check the documentation to create a proper config\
                                file." % self.config_ini_path)


        for option, value in config.items(section):
            if option not in list(self.config_template[section].keys()):
                raise ConfigError('Unknown option, "%s", under section "%s".' % (option, section))
            if 'test' in self.config_template[section][option] and not self.config_template[section][option]['test'](value):
                if 'required' in self.config_template[section][option]:
                    r = self.config_template[section][option]['required']
                    raise ConfigError('Unexpected value ("%s") for option "%s", under section "%s".\
                                        What is expected is %s.' % (value, option, section, r))
                else:
                    raise ConfigError('Unexpected value ("%s") for option "%s", under section "%s".' % (value, option, section))

        self.run.warning('', header="SMTP Configuration is read", lc='cyan')
        for option, value in config.items(section):
            self.run.info(option, value if option != 'password' else '*' * len(value))
            setattr(self, option, value)


    def test(self):
        self.connect()
        self.disconnect()


    def connect(self):
        if not self.server_address or not self.server_port:
            raise ConfigError("SMTP server has not been configured to send e-mails :/")

        try:
           self.server = smtplib.SMTP(self.server_address, self.server_port)

           if self.init_tls:
               self.server.ehlo()
               self.server.starttls()

           if self.username:
               self.server.login(self.username, self.password)

        except Exception as e:
            raise ConfigError("Something went wrong while connecting to the SMTP server :/ This is what we\
                                know about the problem: %s" % e)


    def disconnect(self):
        if self.server:
            self.server.quit()

        self.server = None


    def send(self, to, subject, content):
        self.progress.new('E-mail')
        self.progress.update('Establishing a connection ..')
        self.connect()

        self.progress.update('Preparing the package ..')
        msg = MIMEText(content)
        msg['To'] = to
        msg['Subject'] = subject
        msg['From'] = self.from_address
        msg['Reply-to'] = self.from_address

        try:
            self.progress.update('Sending the e-mail to "%s" ..' % to)
            self.server.sendmail(self.from_address, [to], msg.as_string())
        except Exception as e:
            self.progress.end()
            raise ConfigError("Something went wrong while trying to connet send your e-mail :(\
                                This is what we know about the problem: %s" % e)


        self.progress.update('Disconnecting ..')
        self.disconnect()
        self.progress.end()

        self.run.info('E-mail', 'Successfully sent to "%s"' % to)

