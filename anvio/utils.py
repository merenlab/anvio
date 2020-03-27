# -*- coding: utf-8
# pylint: disable=line-too-long

"""Lonely, helper functions that are broadly used and don't fit anywhere"""

import os
import sys
import gzip
import time
import copy
import socket
import shutil
import smtplib
import hashlib
import Bio.PDB as PDB
import textwrap
import linecache
import webbrowser
import subprocess
import tracemalloc
import configparser
import urllib.request, urllib.error, urllib.parse

import numpy as np
import pandas as pd

from numba import jit
from collections import Counter
from email.mime.text import MIMEText

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Run, Progress, SuppressAllOutput, get_date, TimeCode
from anvio.errors import ConfigError, FilesNPathsError
from anvio.sequence import Composition

with SuppressAllOutput():
    from ete3 import Tree

# psutil is causing lots of problems for lots of people :/
with SuppressAllOutput():
    try:
        import psutil
        PSUTIL_OK=True
    except:
        PSUTIL_OK=False

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


def get_total_memory_usage(keep_raw=False):
    """Get the total memory, including children

    Parameters
    ==========
    keep_raw : bool, False
        A human readable format is returned, e.g. "1.41 GB". If keep_raw, the raw number is
        returned, e.g. 1515601920
    """
    if not PSUTIL_OK:
        return None

    current_process = psutil.Process(os.getpid())
    mem = current_process.memory_info().rss
    for child in current_process.children(recursive=True):
        try:
            mem += child.memory_info().rss
        except:
            pass

    return mem if keep_raw else human_readable_file_size(mem)


def display_top_memory_usage(snapshot, key_type='lineno', limit=10):
    """A pretty-print for the tracemalloc memory usage module

    Modified from https://docs.python.org/3/library/tracemalloc.html

    Examples
    ========
    >>> import tracemalloc
    >>> import anvio.utils as utils
    >>> tracemalloc.start()
    >>> snap = tracemalloc.take_snapshot
    >>> utils.display_top_memory_usage(snap)
    Top 10 lines
    #1: anvio/bamops.py:160: 4671.3 KiB
        constants.cigar_consumption,
    #2: anvio/bamops.py:96: 2571.6 KiB
        self.cigartuples = np.array(read.cigartuples)
    #3: python3.6/linecache.py:137: 1100.0 KiB
        lines = fp.readlines()
    #4: <frozen importlib._bootstrap_external>:487: 961.4 KiB
    #5: typing/templates.py:627: 334.3 KiB
        return type(base)(name, (base,), dct)
    #6: typing/templates.py:923: 315.7 KiB
        class Template(cls):
    #7: python3.6/_weakrefset.py:84: 225.2 KiB
        self.data.add(ref(item, self._remove))
    #8: targets/npyimpl.py:411: 143.2 KiB
        class _KernelImpl(_Kernel):
    #9: _vendor/pyparsing.py:3349: 139.7 KiB
        self.errmsg = "Expected " + _ustr(self)
    #10: typing/context.py:456: 105.1 KiB
        def on_disposal(wr, pop=self._globals.pop):
    3212 other: 4611.9 KiB
    Total allocated size: 15179.4 KiB
    """

    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))


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


def serialize_args(args, single_dash=False, use_underscore=False, skip_keys=None, translate=None):
    cmdline = []
    for param, value in args.__dict__.items():
        if isinstance(skip_keys, list):
            if param in skip_keys:
                continue

        if translate and param in translate:
            param = translate[param]

        dash = '-' if single_dash else '--'

        if not use_underscore:
            param = param.replace('_', '-')

        if value is True:
            cmdline.append('%s%s' % (dash, param))
        elif value is not False and value is not None:
            cmdline.append('%s%s' % (dash, param))
            cmdline.append(str(value))

    return cmdline


def get_predicted_type_of_items_in_a_dict(d, key):
    """Gets a dictionary `d` and a `key` in it, and returns a type function.

    It is a bit counter intuitive. dictionary should look like this:

        d = {'x': {'key': item, (...)},
             'y': {'key': item, (...)},
             (...),
            }

    This is a shitty function, but there was a real need for it, so here we are :/
    """

    items = [x[key] for x in d.values()]

    if(set(items) == set([None])):
        # all items is of type None.
        return None

    not_float = False
    for item in items:
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
            raise ConfigError("Anvi'o searched a bunch of port numbers starting from %d, but failed "
                               "to find an available one for you. Maybe you should specify one :/")
    else:
        if is_port_in_use(port_num):
            raise ConfigError("The port number %d seems to be in use :/" % port_num)

    if os.getuid() and port_num < 1024:
        run.warning("Using the port number %d requires superuser priviliges, which your user does not "
                    "seem to have. Since anvi'o does not know anything about your system configuraiton, "
                    "you are free to go for now. But be prepared for a failed attempt to use this port "
                    "number to serve stuff." % port_num)

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

    raise ConfigError("An anvi'o function needs '%s' to be installed on your system, but it doesn't seem to appear "
                       "in your path :/ If you are certain you have it on your system (for instance you can run it "
                       "by typing '%s' in your terminal window), you may want to send a detailed bug report. Sorry!"\
                        % (program, program))


def format_cmdline(cmdline):
    """Takes a cmdline for `run_command` or `run_command_STDIN`, and makes it beautiful."""
    if not cmdline or (not isinstance(cmdline, str) and not isinstance(cmdline, list)):
        raise ConfigError("You made utils::format_cmdline upset. The parameter you sent to run kinda sucks. It should be string "
                           "or list type. Note that the parameter `shell` for subprocess.call in this `run_command` function "
                           "is always False, therefore if you send a string type, it will be split into a list prior to being "
                           "sent to subprocess.")

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
        raise ConfigError("gzip_decompress_file function is upset because your input file ('%s') does not "
                          "end with a '.gz' extension :(")

    if not output_file_path:
        output_file_path = input_file_path[:-3]

    filesnpaths.is_output_file_writable(output_file_path)

    import gzip
    with gzip.open(input_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
        f_out.writelines(f_in)

    if not keep_original:
        os.remove(input_file_path)

    return output_file_path


class RunInDirectory(object):
    """ Run any block of code in a specified directory. Return to original directory

    Parameters
    ==========
    run_dir : str or Path-like
        The directory the block of code should be run in
    """

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.cur_dir = os.getcwd()
        if not os.path.isdir(self.run_dir):
            raise ConfigError("RunInDirectory :: %s is not a directory." % str(self.run_dir))


    def __enter__(self):
        os.chdir(self.run_dir)


    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.cur_dir)


def run_command(cmdline, log_file_path, first_line_of_log_is_cmdline=True, remove_log_file_if_exists=True):
    """ Uses subprocess.call to run your `cmdline`

    Parameters
    ==========
    cmdline : str or list
        The command to be run, e.g. "echo hello" or ["echo", "hello"]
    log_file_path : str or Path-like
        All stdout from the command is sent to this filepath
    """
    cmdline = format_cmdline(cmdline)

    if anvio.DEBUG:
        Progress().reset()
        Run().info("[DEBUG] `run_command` is running", \
                   ' '.join(['%s' % (('"%s"' % str(x)) if ' ' in str(x) else ('%s' % str(x))) for x in cmdline]), \
                   nl_before=1, nl_after=1, mc='red', lc='yellow')

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
            raise ConfigError("Command failed to run. What command, you say? This: '%s'" % ' '.join(cmdline))
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
    """ Stores a pandas DataFrame as a tab-delimited file.

    Parameters
    ==========
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

    Returns
    =======
    output_path
    """

    filesnpaths.is_output_file_writable(output_path)

    if not columns:
        columns = d.columns

    d.replace(naughty_characters, np.nan, inplace=True)

    d.to_csv(output_path, sep="\t", columns=columns, index=include_index, index_label=index_label, na_rep=rep_str)
    return output_path


def store_dict_as_TAB_delimited_file(d, output_path, headers=None, file_obj=None, key_header=None, keys_order=None):
    '''
        Store a dictionary of dictionaries as a TAB-delimited file.
        input:
            d - dictionary of dictionaries
            output_path - output path
            headers - headers of the secondary dictionary to include (by default include all)
                      these are the columns that will be included in the output file (this
                      doesn't include the first column which is the keys of the major dictionary)
            file_obj - file_object to use for writing
            key_header - the header for the first column (by default: 'key')
    '''
    if not file_obj:
        filesnpaths.is_output_file_writable(output_path)

    if not file_obj:
        f = open(output_path, 'w')
    else:
        f = file_obj

    key_header = key_header if key_header else 'key'
    if not headers:
        headers = [key_header] + sorted(list(d.values())[0].keys())

    f.write('%s\n' % '\t'.join(headers))

    if not keys_order:
        keys_order = sorted(d.keys())
    else:
        missing_keys = [k for k in keys_order if k not in d]
        if len(missing_keys):
            if anvio.DEBUG:
                if len(missing_keys) > 10:
                    raise ConfigError("Some keys (n=%d) are not in your dictionary :/ Here is the first ten "
                                      " of them: %s" % (len(missing_keys), missing_keys[:10].__str__()))
                else:
                    raise ConfigError("Some keys are not in your dictionary :/ Here they are: %s" % missing_keys.__str__())
            else:
                raise ConfigError("Some keys are not in your dictionary :/ Use `--debug` to see where this "
                                  "error is coming from the codebase with a list of example keys that are "
                                  "missing.")

    for k in keys_order:
        line = [str(k)]
        for header in headers[1:]:
            try:
                val = d[k][header]
            except KeyError:
                raise ConfigError("Header ('%s') is not found in the dict :/" % (header))
            except TypeError:
                raise ConfigError("Your dictionary is not properly formatted to be exported "
                                   "as a TAB-delimited file :/ You ask for '%s', but it is not "
                                   "even a key in the dictionary" % (header))

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


@jit(nopython=True)
def add_to_2D_numeric_array(x, y, a, count=1):
    """just-in-time compiled function

    Parameters
    ==========
    x : array
        array of row indices
    y : array
        array of corresponding y indices
    count : int, 1
        How much to add to each coordinate

    Examples
    ========

    Make a 5x20000 array (a) and define 95 coordinate positions to update (i and p)

    >>> a = np.zeros((5, 20000))
    >>> i = np.random.choice(range(5), size=95, replace=True)
    >>> p = np.random.choice(range(100), size=95, replace=False) + 1000

    For comparison, define the slow method

    >>> def add_to_2D_numeric_array_slow(x, y, a, count=1):
    >>>     for idx, pos in zip(x, y):
    >>>         a[idx, pos] += count
    >>>     return a

    Compare the speeds

    >>> %timeit add_to_2D_numeric_array_slow(i, p, a)
    74.5 µs ± 4.42 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)
    >>> %timeit _add_to_2D_numeric_array(i, p, a)
    798 ns ± 12.7 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
    """
    for idx, pos in zip(x, y):
        a[idx, pos] += count

    return a


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


def transpose_tab_delimited_file(input_file_path, output_file_path, remove_after=False):
    filesnpaths.is_file_tab_delimited(input_file_path)
    filesnpaths.is_output_file_writable(output_file_path)

    file_content = [line.strip('\n').split('\t') for line in open(input_file_path, 'rU').readlines()]

    output_file = open(output_file_path, 'w')
    for entry in zip(*file_content):
        output_file.write('\t'.join(entry) + '\n')
    output_file.close()

    if remove_after:
        os.remove(input_file_path)

    return output_file_path


def split_fasta(input_file_path, parts=1, prefix=None):
    if not prefix:
        prefix = os.path.abspath(input_file_path)

    filesnpaths.is_file_exists(input_file_path)
    filesnpaths.is_file_fasta_formatted(input_file_path)

    source = u.ReadFasta(input_file_path, quiet=True)
    length = len(source.ids)

    if length < parts:
        parts = length

    chunk_size = length // parts

    output_files = []

    for part_no in range(parts):
        output_file = prefix + '.' + str(part_no)

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
        output_files.append(output_file)

    return output_files


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
                    raise ConfigError("get_column_data_from_TAB_delim_file is speaking: The file you sent "
                                       "does not have data for the column index %d. Something is wrong :/" % (index))

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
    columns = input_matrix.readline().strip('\n').split('\t')[1:]

    fields_of_interest = []
    if cols_to_return:
        fields_of_interest = [columns.index(col) for col in cols_to_return]
    else:
        fields_of_interest = [f for f in range(0, len(columns)) if constants.IS_ESSENTIAL_FIELD(columns[f])]

    # update columns:
    columns = [columns[i] for i in fields_of_interest]

    if not len(columns):
        raise ConfigError("Only a subset (%d) of fields were requested by the caller, but none of them was found "
                           "in the matrix (%s) :/" % (len(cols_to_return), file_path))

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


def get_required_packages_for_enrichment_test():
    ''' Return a dict with the packages as keys and installation instrucstions as values'''
    packages = ["tidyverse", "stringi", "magrittr", "qvalue", "optparse"]

    installation_instructions = ["conda install -c r r-tidyverse",
                                 "conda install -c r r-stringi",
                                 "conda install -c bioconda r-magrittr",
                                 "conda install -c bioconda bioconductor-qvalue",
                                 "conda install -c conda-forge r-optparse"]

    return dict(zip(packages,installation_instructions))


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
        raise ConfigError("%s is not a valid key for creating a pandas dataframe of values of gene_level_coverage_stats_dict. "
                           "Here is a list of the valid keys: %s" % (key, list(legal_keys)))

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


def get_indices_for_outlier_values(c):
    is_outlier = get_list_of_outliers(c)
    return set([p for p in range(0, c.size) if is_outlier[p]])


def get_list_of_outliers(values, threshold=None, zeros_are_outliers=False, median=None):
    """Return boolean array of whether values are outliers (True means yes)

    Modified from Joe Kington's (https://stackoverflow.com/users/325565/joe-kington)
    implementation computing absolute deviation around the median.

    Parameters
    ==========
    values : array-like
        An numobservations by numdimensions array of observations

    threshold : number, None
        The modified z-score to use as a thresholdold. Observations with
        a modified z-score (based on the median absolute deviation) greater
        than this value will be classified as outliers.

    median : array-like, None
        Pass median value of values if you already calculated it to save time

    Returns
    =======
    mask : numpy array (dtype=bool)
        A numobservations-length boolean array.

    References
    ==========
    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

    http://www.sciencedirect.com/science/article/pii/S0022103113000668
    """

    if threshold is None:
        threshold = 1.5

    if len(values.shape) == 1:
        values = values[:, None]

    if not median: median = np.median(values, axis=0)

    diff = np.sum((values - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    median_absolute_deviation = np.median(diff)

    if not median_absolute_deviation:
       if values[0] == 0:
            # A vector of all zeros is considered "all outliers"
            return np.array([True] * values.size)
       else:
            # A vector of uniform non-zero values is "all non-outliers"
            # This could be important for silly cases (like in megahit) in which there is a maximum value for coverage
            return np.array([False] * values.size)

    modified_z_score = 0.6745 * diff / median_absolute_deviation
    non_outliers = modified_z_score > threshold

    if not zeros_are_outliers:
        return non_outliers
    else:
        zero_positions = [x for x in range(len(values)) if values[x] == 0]
        for i in zero_positions:
            non_outliers[i] = True
        return non_outliers


def get_gene_caller_ids_from_args(gene_caller_ids, delimiter=','):
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
        raise ConfigError("The gene calls you provided do not look like gene callers anvi'o is used to working with :/ Here is "
                          "one of them: '%s' (%s)." % (g, type(g)))
    return gene_caller_ids_set


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
        raise ConfigError("get_synonymous_and_non_synonymous_potential :: sequence does not end "
                          "with a stop codon and is therefore probably not what you want. If you "
                          "want to continue anyways, use the just_do_it flag")

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


def compare_times(calls, as_matrix=False, iterations_per_call=1):
    """Compare times between function calls

    Parameters
    ==========
    calls : list of tuples
        Each element should be a (name, function, args, kwargs) tuples. If there are no args or
        kwargs, the element should look like (name, function, [], {})

    as_matrix : bool, False
        If True, results are output as a pandas matrix, where each element is a time difference between
        calls. Otherwise, a dictionary is returned

    iterations_per_call : int, 1
        How many times should each function call be ran? Time will be averaged

    Returns
    =======
    times : pd.DataFrame or dict
        If as_matrix, pd.DataFrame is returned, where times[i, j] is how much faster i is than j.
        Otherwise, dictionary of {name: time} is returned
    """

    call_times = np.zeros((len(calls), iterations_per_call))
    names, *_ = zip(*calls)
    for i, call in enumerate(calls):
        name, function, args, kwargs = call

        for j in range(iterations_per_call):
            try:
                with TimeCode(quiet=True) as t:
                    function(*args, **kwargs)
            except:
                raise ConfigError("compare_times :: function call with name '%s' failed." % name)

            call_times[i, j] = t.time.total_seconds()

    averaged_call_times = np.mean(call_times, axis=1)

    if not as_matrix:
        return dict(zip(names, averaged_call_times))

    matrix = []
    for i, _time in enumerate(call_times):
        row = []

        for j, _time in enumerate(call_times):
            row.append(averaged_call_times[j] - averaged_call_times[i] if i > j else 'NA')

        matrix.append(row)

    return pd.DataFrame(matrix, columns=names, index=names)


def concatenate_files(dest_file, file_list, remove_concatenated_files=False):
    if not dest_file:
        raise ConfigError("Destination cannot be empty.")

    filesnpaths.is_output_file_writable(dest_file)

    if not len(file_list):
        raise ConfigError("File list cannot be empty.")

    for f in file_list:
        filesnpaths.is_file_exists(f)

    dest_file_obj = open(dest_file, 'w')
    for chunk_path in file_list:
        for line in open(chunk_path, 'rU'):
            dest_file_obj.write(line)

    dest_file_obj.close()

    if remove_concatenated_files:
        for f in file_list:
            os.remove(f)

    return dest_file


def get_split_start_stops(contig_length, split_length, gene_start_stops=None):
    """Wrapper function for get_split_start_stops_with_gene_calls and get_split_start_stops_without_gene_calls"""
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
    split_names_of_interest = set(splits_info.keys())

    contigs_db.disconnect()

    return (split_names_of_interest, contig_names_of_interest)


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


@jit(nopython=True)
def get_constant_value_blocks(array, value):
    """Generator that returns blocks of consecutive numbers

    Parameters
    ==========
    array : array
        a numerical numpy array. If a list is passed, this function is very slow

    value : number
        The number you want to get constant blocks for.

    Examples
    ========

    >>> a = np.array([47, 47, 47, 49, 50, 47, 47, 99])
    >>> for i in get_constant_value_blocks(a, 47): print(i)
    (0, 3)
    (5, 7)
    """
    ans = []

    matching = False
    for i in range(len(array)):
        if array[i] == value:
            if not matching:
                start = i
                matching = True
        else:
            if matching:
                matching = False
                ans.append((start, i))

    if matching:
        ans.append((start, i + 1))

    return ans


@jit(nopython=True)
def find_value_index(x, val, reverse_search=False):
    """Returns first instance of indices where a value is found

    Created this because unlike np.where, this stops after the first instance is found. If you only
    want the first instance, this algorithm is therefore preferrable for array sizes < 1000 (see
    examples)

    Parameters
    ==========
    x : 1D array

    val : number
        return index of x where the value == val.

    reverse_search : bool, False
        Search in reverse order

    Examples
    ========
    >>> import numpy as np
    >>> import anvio.utils as utils
    >>> x = np.arange(1000)
    >>> %timeit utils.find_value_index(x, 999, rev=True)
    574 ns ± 15.8 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
    >>> %timeit utils.find_value_index(x, 999)
    2.21 µs ± 36.7 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    >>> %timeit np.where(x == 999)[0][0]
    2.91 µs ± 563 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    """

    for i in range(len(x)) if not reverse_search else range(len(x)-1, -1, -1):
        if x[i] == val:
            return i


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
        An object that looks like that ContigsSuperClass.contig_sequences (initialized with
        ContigsSuperClass.init_contig_sequences)
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


def get_DNA_sequence_translated(sequence, gene_callers_id, return_with_stops=False):
    sequence = sequence.upper()

    if len(sequence) % 3.0 != 0:
        raise ConfigError("The sequence corresponds to the gene callers id '%s' does not seem to "
                           "have proper number of nucleotides to be translated :/ Here it is: %s" % (gene_callers_id, sequence))

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

    if gene_call['partial']:
        raise ConfigError("get_codon_order_to_nt_positions_dict: this simply will not work "
                           "for partial gene calls, and this on *is* a partial one.")

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


def get_bin_name_from_item_name(anvio_db_path, item_name, collection_name=None):
    is_pan_or_profile_db(anvio_db_path, genes_db_is_also_accepted=True)
    database = db.DB(anvio_db_path, None, ignore_version=True)

    if t.collections_splits_table_name not in database.get_table_names():
        raise ConfigError("The database %s does not contain a collections table :/")

    if collection_name:
        where_clause = 'split = "%s" and collection_name = "%s"' % (item_name, collection_name)
    else:
        where_clause = 'split = "%s"' % (item_name)

    rows = database.get_some_rows_from_table(t.collections_splits_table_name, where_clause=where_clause)

    database.disconnect()

    return rows


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



def is_this_name_OK_for_database(variable_name, content, stringent=True):
    if not content:
        raise ConfigError("But the %s is empty? Come on :(" % variable_name)

    if content[0] in constants.digits:
        raise ConfigError("Sorry, %s can't start with a digit. Long story. Please specify a name "
                           "that starts with an ASCII letter." % variable_name)

    if stringent:
        allowed_chars = constants.allowed_chars.replace('.', '').replace('-', '')
    else:
        allowed_chars = constants.allowed_chars.replace('.', '')

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


def create_fasta_dir_from_sequence_sources(genome_desc, fasta_txt=None):
    """genome_desc is an instance of GenomeDescriptions"""

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
            args = anvio.summarizer.ArgsTemplateForSummarizerClass()
            args.contigs_db = genome_subset['contigs_db_path']
            args.profile_db = genome_subset['profile_db_path']
            args.collection_name = genome_subset['collection_id']
            args.output_dir = filesnpaths.get_temp_directory_path(just_the_path=True)
            args.quick_summary = True

            # note that we're initializing the summary class only for once for a given
            # genome subset
            summary = anvio.summarizer.ProfileSummarizer(args, r=Run(verbose=False))
            summary.init()

            for genome_name, bin_name in genome_subset['genome_name_bin_name_tpl']:
                genome_names.add(genome_name)

                hash_for_output_file = hashlib.sha256(genome_name.encode('utf-8')).hexdigest()
                hash_to_name[hash_for_output_file] = genome_name
                path = os.path.join(temp_dir, hash_for_output_file + '.fa')
                file_paths.add(path)
                name_to_path[genome_name] = path

                bin_summary = anvio.summarizer.Bin(summary, bin_name)

                with open(path, 'w') as fasta:
                    fasta.write(bin_summary.get_bin_sequence())


    if fasta_txt is not None:
        fastas = get_TAB_delimited_file_as_dictionary(fasta_txt, expected_fields=['name', 'path'], only_expected_fields=True)
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


def gen_NEXUS_format_partition_file_for_phylogenomics(partition_file_path, sequence_lengths, separator='', run=run, progress=progress):
    """ Generates a NEXUS-formatted partition file for phylogenomics. See
        https://github.com/merenlab/anvio/issues/1333 for details

    Parameters
    ==========
    partition_file_path: `str`
        File path to be generated.
    sequence_lengths: `list` of `tuples`
        A list that contins sequence names and lenghts as tuples. I.e.,
        [('seq_1', 100), ('seq_2', 42), ...]
    separator: `str`
        Characters used to separate sequences from each other in a multi-alignment
        file.
    run: `object`
        Anvi'o run object
    run: `progress`
        Anvi'o progress object

    Returns
    =======
    None
    """

    filesnpaths.is_output_file_writable(partition_file_path)

    if not isinstance(sequence_lengths, list):
        raise ConfigError("Sequence lengths must be passed as a list of tuples.")

    if not isinstance(sequence_lengths[0], tuple):
        raise ConfigError("Sequence lengths must be passed as a list of tuples.")

    with open(partition_file_path, 'w') as partition_file:
        partition_file.write("#nexus\nbegin sets;\n")
        index = 1
        for sequence_name, sequence_length in sequence_lengths:
            partition_file.write("    charset %s = %d-%d;\n" % (sequence_name, index, index + sequence_length - 1))
            index += (sequence_length + len(separator))
        partition_file.write("end;\n")

    progress.reset()
    run.info("Partition file", partition_file_path, mc='yellow')
    run.info_single("Your partition file is ready. Please do not forget to replace placeholders for model names ('[MODEL]') "
                    "in this file with appropriate model names prior to your phylogenomic analysis.", nl_before=1, nl_after=1)


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
        raise ConfigError("I can't unique this. Output FASTA file path can't be identical to "
                           "the names file path...")

    if output_fasta_path == input_file_path or names_file_path == input_file_path:
        raise ConfigError("Anvi'o will not unique this. Output FASTA path and names file path should "
                           "be different from the the input file path...")

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
        if wrap_from:
            output.write('%s\n' % textwrap.fill(d[key], wrap_from, break_on_hyphens=False))
        else:
            output.write('%s\n' % (d[key]))

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
                                        none_value=None, empty_header_columns_are_OK=False, return_failed_lines=False):
    """Takes a file path, returns a dictionary.

       - If `return_failed_lines` is True, it the function will not throw an exception, but instead
         return a list of `failed_lines` along with a dictionary of final results.
    """

    if expected_fields and (not isinstance(expected_fields, list) and not isinstance(expected_fields, set)):
        raise ConfigError("'expected_fields' variable must be a list (or a set).")

    if only_expected_fields and not expected_fields:
        raise ConfigError("'only_expected_fields' variable guarantees that there are no more fields present "
                           "in the input file but the ones requested with 'expected_fields' variable. If you "
                           "need to use this flag, you must also be explicit about what fields you expect to "
                           "find in the file.")

    filesnpaths.is_file_plain_text(file_path)
    filesnpaths.is_file_tab_delimited(file_path, separator=separator)

    failed_lines = []
    column_mapping_for_line_failed = None

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
            raise  ConfigError("Number of column names declared (%d) differs from the number of columns "
                                "found (%d) in the matrix ('%s') :/" % (len(columns), num_fields, file_path))

        # now we set the column names. if the file had its header, we must discard
        # the first line. so here we go:
        if not no_header:
            f.readline()
    else:
        columns = f.readline().strip('\n').split(separator)

    if not empty_header_columns_are_OK and min(map(len, columns)) == 0:
        raise ConfigError("At least one of the column headers in your tab delimited file '%s' "
                          "is empty." % file_path)

    if expected_fields:
        for field in expected_fields:
            if field not in columns:
                raise ConfigError("The file '%s' does not contain the right type of header. It was expected "
                                   "to have these: '%s', however it had these: '%s'" % (file_path,
                                                                                        ', '.join(expected_fields),
                                                                                        ', '.join(columns[1:])))

    if only_expected_fields:
        for field in columns:
            if field not in expected_fields:
                raise ConfigError("There are more fields in the file '%s' than the expected fields :/ "
                                  "Anvi'o is telling you about this because get_TAB_delimited_file_as_dictionary "
                                  "function is called with `only_expected_fields` flag turned on." % (file_path))

    d = {}
    line_counter = 0

    for line in f.readlines():
        if ascii_only:
            if not is_ascii_only(line):
                raise ConfigError("The input file conitans non-ascii characters at line number %d. Those lines "
                                   "either should be removed, or edited." % (line_counter + 2))

        line_fields = [f if f else None for f in line.strip('\n').split(separator)]

        if line_fields and line_fields[0] == None:
            raise ConfigError("The line number %d in '%s' has no data in its first column, and this doesn't "
                              "seem right at all :/" % (line_counter + 1, file_path))

        if column_mapping:
            column_mapping_for_line_failed = False
            updated_line_fields = []
            for i in range(0, len(line_fields)):
                try:
                    if line_fields[i] == None and column_mapping[i] in [float, int]:
                        updated_line_fields.append(column_mapping[i](0))
                    else:
                        updated_line_fields.append(column_mapping[i](line_fields[i]))
                except NameError:
                    if return_failed_lines:
                        failed_lines.append(line_counter + 1)
                        column_mapping_for_line_failed = True
                    else:
                        raise ConfigError("Mapping function '%s' did not work on value '%s'. These functions can be native "
                                          "Python functions, such as 'str', 'int', or 'float', or anonymous functions "
                                          "defined using lambda notation." % (column_mapping[i], line_fields[i]))
                except TypeError:
                    if return_failed_lines:
                        failed_lines.append(line_counter + 1)
                        column_mapping_for_line_failed = True
                    else:
                        raise ConfigError("Mapping function '%s' does not seem to be a proper Python function :/" % column_mapping[i])
                except ValueError:
                    if return_failed_lines:
                        failed_lines.append(line_counter + 1)
                        column_mapping_for_line_failed = True
                    else:
                        raise ConfigError("Mapping function '%s' did not like the value '%s' in column number %d "
                                          "of the input matrix '%s' :/" % (column_mapping[i], line_fields[i], i + 1, file_path))

            line_fields = updated_line_fields

        if column_mapping_for_line_failed:
            continue

        if indexing_field == -1:
            entry_name = 'line__%09d__' % line_counter
        else:
            entry_name = line_fields[indexing_field]

        if entry_name in d:
            raise ConfigError("The entry name %s appears more than once in the TAB-delimited file '%s'. We assume that you "
                              "did not do it that purposefully, but if you need this file in this form, then feel free to "
                              "contact us so we can try to find a solution for you. But if you have gotten this error while "
                              "working with HMMs, do not contact us since helping you in that case is beyond us (see the issue "
                              "#1206 for details))." % (entry_name, file_path))

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
                    raise ConfigError("Appending entries to the already existing dictionary from file '%s' failed "
                                       "as the entry %s does not appear to be in the file." % (file_path, entry))
                else:
                    for key in keys:
                        dict_to_append[entry][key] = none_value
            else:
                for key in keys:
                    dict_to_append[entry][key] = d[entry][key]

        return dict_to_append

    # this is here for backward compatibility.
    if return_failed_lines:
        return d, failed_lines

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

    R = lambda f: open(os.path.join(source, f), 'rU').readlines()[0].strip()
    for source in source_dirs:
        if source.endswith('/'):
            source = source[:-1]

        if not PROPER(os.path.basename(source)):
            raise ConfigError("One of the search database directories ('%s') contains characters in its name "
                               "anvio does not like. Directory names should be at least three characters long "
                               "and must not contain any characters but ASCII letters, digits and "
                               "underscore" % os.path.basename(source))

        expected_files = ['reference.txt', 'kind.txt', 'genes.txt', 'genes.hmm.gz', 'target.txt', 'noise_cutoff_terms.txt']

        missing_files = [f for f in expected_files if not os.path.exists(os.path.join(source, f))]
        if missing_files:
            raise ConfigError("Each search database directory must contain following files: %s'. Yet, the HMM source '%s' seems to "
                              "be missing the follwoing one(s): %s. See this blog post to make sure you are doing it the way it "
                              "should be done: http://merenlab.org/2016/05/21/archaeal-single-copy-genes/" % \
                                            (', '.join(expected_files), os.path.basename(source), ', '.join(missing_files)))

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

        raise ConfigError("Oh no :( We recently changed the description of the stacked bar data type, and your input data "
                          "file still has the older version. Here is the list of those that are violating the new format: "
                          "%s. To avoid this issue and to turn them into the new format, you could take '%s', and present "
                          "it as %d separate TAB-delimited entries that look like this: %s. Sorry!" % \
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
        raise ConfigError("Accession IDs in your HMM model should be unique, however, the `genes.hmm.gz` "
                          "file for `%s` seems to have the same accession ID (the line that starts with `ACC`) "
                          "more than once :(" % (os.path.abspath(model_path).split('/')[-2]))

    if len(genes.difference(genes_in_model)):
        raise ConfigError("Some gene names in genes.txt file does not seem to be appear in genes.hmm.gz. "
                          "Here is a list of missing gene names: %s" % ', '.join(list(genes.difference(genes_in_model))))

    if len(genes_in_model.difference(genes)):
        raise ConfigError("Some gene names in genes.hmm.gz file does not seem to be appear in genes.txt. "
                          "Here is a list of missing gene names: %s" % ', '.join(list(genes_in_model.difference(genes))))


def get_missing_programs_for_hmm_analysis():
    missing_programs = []
    for p in ['prodigal', 'hmmscan']:
        try:
            is_program_exists(p)
        except ConfigError:
            missing_programs.append(p)
    return missing_programs


def get_genes_database_path_for_bin(profile_db_path, collection_name, bin_name):
    if not collection_name or not bin_name:
        raise ConfigError("Genes database must be associted with a collection name and a bin name :/")

    return os.path.join(os.path.dirname(profile_db_path), 'GENES', '%s-%s.db' % (collection_name, bin_name))


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
        raise ConfigError("Anvi'o was trying to get the version of the -alleged- anvi'o database '%s', but it failed "
                           "because it turns out it doesn't know anything about this '%s' type." % (db_path, db_type))

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

    elif db_type == 'genes':
        return set([str(i) for i in database.get_single_column_from_table(t.gene_level_coverage_stats_table_name, 'sample_name')])

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
        raise ConfigError("`get_all_sample_names_from_the_database` function does not know how to deal "
                           "with %s databases." % db_type)


def get_all_item_names_from_the_database(db_path, run=run):
    """Return all split names or gene cluster names in a given database"""

    all_items = set([])

    database = db.DB(db_path, get_required_version_for_db(db_path))
    db_type = database.get_meta_value('db_type')

    if db_type == 'profile':
        if is_blank_profile(db_path):
            run.warning("Someone asked for the split names in a blank profile database. Sadly, anvi'o does not keep track "
                        "of split names in blank profile databases. This function will return an empty set as split names "
                        "to not kill your mojo, but whatever you were trying to do will not work :(")
            return set([])
        elif int(database.get_meta_value('merged')):
            all_items = set(database.get_single_column_from_table('mean_coverage_Q2Q3_splits', 'contig'))
        else:
            all_items = set(database.get_single_column_from_table('atomic_data_splits', 'contig'))
    elif db_type == 'pan':
        all_items = set(database.get_single_column_from_table(t.pan_gene_clusters_table_name, 'gene_cluster_id'))
    elif db_type == 'contigs':
        all_items = set(database.get_single_column_from_table(t.splits_info_table_name, 'split'))
    elif db_type == 'genes':
        all_items = set([str(i) for i in database.get_single_column_from_table(t.gene_level_coverage_stats_table_name, 'gene_callers_id')])
    else:
        database.disconnect()
        raise ConfigError("You wanted to get all items in the database %s, but no one here knows about its type. Seriously,\
                            what is '%s' anyway?" % (db_path, db_type))

    if not len(all_items):
        database.disconnect()
        raise ConfigError("utils::get_all_item_names_from_the_database speaking. Something that should never happen happened :/ "
                          "There seems to be nothing in this %s database. Anvi'o is as confused as you are. Please get in touch "
                          "with a developer. They will love this story." % db_path)

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
        raise ConfigError("anvi'o does not recognize %s as being a variability table generated by "
                          "anvi-gen-variability-profile." % table_path)


def is_contigs_db(db_path):
    filesnpaths.is_file_exists(db_path)
    if get_db_type(db_path) != 'contigs':
        raise ConfigError("'%s' is not an anvi'o contigs database." % db_path)
    return True


def is_pan_or_profile_db(db_path, genes_db_is_also_accepted=False):
    ok_db_types = ['pan', 'profile']

    if genes_db_is_also_accepted:
        ok_db_types += ['genes']

    db_type = get_db_type(db_path)

    if db_type not in ok_db_types:
        if genes_db_is_also_accepted:
            raise ConfigError("'%s' is not a pan, profile, or a genes database :/ Anvi'o wants what it wants and this "
                              "'%s' database is not it." % (db_path, db_type))
        else:
            raise ConfigError("'%s' is neither a pan nor a profile database :/ Someone is in trouble (*cough* 'someone' "
                               "being whoever sent this %s database as a parameter to that command *cough*)." % (db_path, db_type))

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


def get_enriched_groups(props, reps):
    '''
        Accepts a vector of proportions and number of replicates per group and
        returns a boolean vector where each group that has proportion above
        the "expected" (i.e. the overall proportion) is True and the rest are False.
    '''
    # if the function doesn't occur at all then test_statistic is zero and p-value is 1
    if not np.count_nonzero(props):
        return np.zeros(len(props))
    overall_portion = np.sum(np.multiply(props, reps)) / np.sum(reps)

    return props > overall_portion


def is_pan_db(db_path):
    if get_db_type(db_path) != 'pan':
        raise ConfigError("'%s' is not an anvi'o pan database." % db_path)
    return True


def is_genes_db(db_path):
    if get_db_type(db_path) != 'genes':
        raise ConfigError("'%s' is not an anvi'o genes database." % db_path)
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
        raise ConfigError('The contigs database and the profile database does not '
                          'seem to be compatible. More specifically, this contigs '
                          'database is not the one that was used when %s generated '
                          'this profile database (%s != %s).'\
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
        raise ConfigError('The contigs and structure databases do not seem compatible. '
                          'More specifically, the contigs database is not the one that '
                          'was used when the structure database was created (%s != %s).'\
                               % (a_hash, p_hash))

    return True


def download_file(url, output_file_path, progress=progress, run=run):
    filesnpaths.is_output_file_writable(output_file_path)

    try:
        response = urllib.request.urlopen(url)
    except Exception as e:
        raise ConfigError("Something went wrong with your download attempt. Here is the "
                           "problem: '%s'" % e)

    file_size = 0
    if 'Content-Length' in response.headers:
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

            if file_size:
                progress.update('%.1f%%' % (downloaded_size * 100.0 / file_size))
            else:
                progress.update('%s' % human_readable_file_size(downloaded_size))
        else:
            break

    f.close()

    progress.end()
    run.info('Downloaded succesfully', output_file_path)


def get_remote_file_content(url, gzipped=False):
    import requests
    from io import BytesIO

    remote_file = requests.get(url)

    if remote_file.status_code == 404:
        raise ConfigError("Bad news. The remote file at '%s' was not found :(" % url)

    if gzipped:
        buf = BytesIO(remote_file.content)
        fg = gzip.GzipFile(fileobj=buf)
        return fg.read().decode('utf-8')

    return remote_file.content.decode('utf-8')


def download_protein_structure(protein_code, output_path=None, chain=None, raise_if_fail=True):
    """Downloads protein structures using Biopython.

    Parameters
    ==========
    protein_code : str
        Each element is a 4-letter protein code

    output_path : str
        Path where structure is written to. Temporary directory is chosen if None

    chain : str, None
        If None, all chains remain in the PDB file. If specified, only the chain with the chain ID
        `chain` will be saved.

    raise_if_fail : bool, True
        If the file does not download, raise an error

    Returns
    =======
    output : output_path
        Returns the filepath of the written file. Returns None if download failed
    """

    output_dir = os.path.dirname(output_path)
    if output_dir == '': output_dir = '.'

    pdb_list = PDB.PDBList()

    try:
        with SuppressAllOutput():
            # We suppress output that looks like this:
            # >>> WARNING: The default download format has changed from PDB to PDBx/mmCif
            # >>> Downloading PDB structure '5w6y'...
            temp_output_path = pdb_list.retrieve_pdb_file(protein_code, file_format='pdb', pdir=output_dir, overwrite=True)
    except:
        pass

    if not filesnpaths.is_file_exists(temp_output_path, dont_raise=True):
        # The file wasn't downloaded
        if raise_if_fail:
            raise ConfigError("The protein %s could not be downloaded. Are you connected to internet?" % protein_code)
        else:
            return None

    if chain is not None:
        class ChainSelect(PDB.Select):
            def accept_chain(self, chain_obj):
                x = 1 if chain_obj.get_id() == chain else 0
                return 1 if chain_obj.get_id() == chain else 0

        p = PDB.PDBParser()
        try:
            structure = p.get_structure(None, temp_output_path)
        except:
            # FIXME Something very rare happened on Biopython's end. We silently return the whole
            # file instead of only the chain. Here is one such reason for failure we stumbled upon:
            # https://github.com/biopython/biopython/issues/2819
            shutil.move(temp_output_path, output_path)
            return output_path

        # Overwrite file with chain-only structure
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(temp_output_path, ChainSelect())

    shutil.move(temp_output_path, output_path)

    return output_path


def get_hash_for_list(l):
    return 'hash' + str(hashlib.sha224(''.join(sorted(list(l))).encode('utf-8')).hexdigest()[0:8])


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
        raise ConfigError("You want to export SVGs? Well, you need the Python library 'selenium' to be able to "
                          "do that but you don't have it. If you are lucky, you probably can install it by "
                          "typing 'pip install selenium' or something :/")

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


def check_h5py_module():
    """To make sure we do have the h5py module.

       The reason this function is here is becasue we removed h5py from anvi'o dependencies,
       but some migration scripts may still need it if the user has very old databases. In
       those cases the user must install it manually."""

    try:
        import h5py
    except:
        raise ConfigError("Please install the Python module `h5py` manually for this migration task to continue. "
                          "The reason why the standard anvi'o installation did not install module is complicated, "
                          "and really unimportant. If you run `pip install h5py` in your Python virtual environmnet "
                          "for anvi'o, and try running the migration program again things should be alright.")


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
            raise ConfigError("Well, the file '%s' does not seem to be a config file at all :/ Here "
                               "is what the parser had to complain about it: %s" % (self.config_ini_path, e))

        section = 'SMTP'

        if section not in config.sections():
            raise ConfigError("The config file '%s' does not seem to have an 'SMTP' section, which "
                               "is essential for Mailer class to learn server and authentication "
                               "settings. Please check the documentation to create a proper config "
                               "file." % self.config_ini_path)


        for option, value in config.items(section):
            if option not in list(self.config_template[section].keys()):
                raise ConfigError('Unknown option, "%s", under section "%s".' % (option, section))
            if 'test' in self.config_template[section][option] and not self.config_template[section][option]['test'](value):
                if 'required' in self.config_template[section][option]:
                    r = self.config_template[section][option]['required']
                    raise ConfigError('Unexpected value ("%s") for option "%s", under section "%s". '
                                       'What is expected is %s.' % (value, option, section, r))
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
            raise ConfigError("Something went wrong while connecting to the SMTP server :/ This is what we "
                               "know about the problem: %s" % e)


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
            raise ConfigError("Something went wrong while trying to connet send your e-mail :( "
                               "This is what we know about the problem: %s" % e)


        self.progress.update('Disconnecting ..')
        self.disconnect()
        self.progress.end()

        self.run.info('E-mail', 'Successfully sent to "%s"' % to)

