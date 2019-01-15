# -*- coding: utf-8
# pylint: disable=line-too-long
"""Relations with the console output, Progress and Run classes"""

import os
import re
import sys
import time
import fcntl
import struct
import termios
import textwrap

from colored import fore, back, style

import anvio.constants as constants
import anvio.dictio as dictio

from anvio.errors import TerminalError
from anvio.ttycolors import color_text as c

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


# clean garbage garbage:
ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
non_ascii_escape = re.compile(r'[^\x00-\x7F]+')
CLEAR = lambda line: ansi_escape.sub('', non_ascii_escape.sub('', line.strip()))


class SuppressAllOutput(object):
    def __enter__(self):
        sys.stderr.flush()
        self.old_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        sys.stdout.flush()
        self.old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.flush()
        sys.stderr = self.old_stderr
        sys.stdout.flush()
        sys.stdout = self.old_stdout


def remove_spaces(text):
    while True:
        if text.find("  ") > -1:
            text = text.replace("  ", " ")
        else:
            break

    return text


class Progress:
    def __init__(self, verbose=True):
        self.pid = None
        self.verbose = verbose
        self.terminal_width = None
        self.is_tty = sys.stdout.isatty()

        self.get_terminal_width()

        self.current = None

        self.progress_total_items = None
        self.progress_current_item = 0

        self.LEN = lambda s: len(s.encode('utf-16-le')) // 2


    def get_terminal_width(self):
        try:
            self.terminal_width = max(get_terminal_size()[0], 120)
        except:
            self.terminal_width = 120


    def new(self, pid, discard_previous_if_exists=False, progress_total_items=None):
        if self.pid:
            if discard_previous_if_exists:
                self.end()
            else:
                raise TerminalError("Progress.new() can't be called before ending the previous one (Existing: '%s', Competing: '%s')." % (self.pid, pid))

        if not self.verbose:
            return

        self.pid = '%s %s' % (get_date(), pid)
        self.get_terminal_width()
        self.current = None
        self.step = None
        self.progress_total_items = progress_total_items
        self.progress_current_item = 0


    def increment(self, increment_to=None):
        if increment_to:
            self.progress_current_item = increment_to
        else:
            self.progress_current_item += 1


    def write(self, c, dont_update_current=False):
        surpass = self.terminal_width - self.LEN(c)

        if surpass < 0:
            c = c[0:-(-surpass + 5)] + ' (...)'
        else:
            if not dont_update_current:
                self.current = c

            c = c + ' ' * surpass

        if self.verbose:
            if self.progress_total_items and self.is_tty:
                p_text = ' %d%% ⚙  ' % (self.progress_current_item * 100 / self.progress_total_items)
                p_length = self.LEN(p_text)

                msg_length = self.LEN(c)
                break_point = round(msg_length * self.progress_current_item / self.progress_total_items)

                # see a full list of color codes: https://gitlab.com/dslackw/colored
                if p_length >= break_point:
                    sys.stderr.write(back.CYAN + fore.BLACK + c[:break_point] + \
                                     back.GREY_30 + fore.WHITE + c[break_point:] + \
                                     style.RESET)
                else:
                    sys.stderr.write(back.CYAN + fore.BLACK + c[:break_point - p_length] + \
                                     back.SALMON_1 + fore.BLACK + p_text + \
                                     back.GREY_30 + fore.WHITE + c[break_point:] + \
                                     style.RESET)
                sys.stderr.flush()
            else:
                sys.stderr.write(back.CYAN + fore.BLACK + c + style.RESET)
                sys.stderr.flush()


    def reset(self):
        self.clear()


    def clear(self):
        if not self.verbose:
            return

        null = '\r' + ' ' * (self.terminal_width)
        sys.stderr.write(null)
        sys.stderr.write('\r')
        sys.stderr.flush()
        self.current = None
        self.step = None


    def append(self, msg):
        if not self.verbose:
            return
        self.write('%s%s' % (self.current, msg))


    def step_start(self, step, symbol="⚙ "):
        if not self.pid:
            raise TerminalError("You don't have an active progress to do it :/")

        if not self.current:
            raise TerminalError("You don't have a current progress bad :(")

        if self.step:
            raise TerminalError("You already have an unfinished step :( Here it is: '%s'." % self.step)

        if not self.verbose:
            return

        self.step = " / %s " % (step)

        self.write(self.current + self.step + symbol, dont_update_current=True)


    def step_end(self, symbol="👍"):
        if not self.step:
            raise TerminalError("You don't have an ongoing step :(")

        if not self.verbose:
            return

        self.write(self.current + self.step + symbol)

        self.step = None


    def update(self, msg):
        self.msg = msg

        if not self.verbose:
            return

        if not self.pid:
            raise TerminalError('Progress with null pid will not update for msg "%s"' % msg)

        self.clear()
        self.write('\r[%s] %s' % (self.pid, msg))


    def end(self):
        self.pid = None
        if not self.verbose:
            return
        self.clear()


class Run:
    def __init__(self, log_file_path=None, verbose=True, width=45):
        self.log_file_path = log_file_path

        self.info_dict = {}
        self.verbose = verbose
        self.width = width

        self.single_line_prefixes = {1: '* ',
                                     2: '    - ',
                                     3: '        > '}


    def log(self, line):
        if not self.log_file_path:
            self.warning("The run object got a logging request, but it was not inherited with\
                          a log file path :(")
            return

        with open(self.log_file_path, "a") as log_file: log_file.write('[%s] %s\n' % (get_date(), CLEAR(line)))


    def write(self, line, quiet=False, overwrite_verbose=False):
        if self.log_file_path:
            self.log(line)

        if (self.verbose and not quiet) or overwrite_verbose:
            try:
                sys.stderr.write(line)
            except:
                sys.stderr.write(line.encode('utf-8'))


    def info(self, key, value, quiet=False, display_only=False, nl_before=0, nl_after=0, lc='cyan', mc='yellow', progress=None):
        if not display_only:
            self.info_dict[key] = value

        if isinstance(value, bool):
            pass
        elif isinstance(value, str):
            value = remove_spaces(value)
        elif isinstance(value, int):
            value = pretty_print(value)

        label = constants.get_pretty_name(key)

        info_line = "%s%s %s: %s\n%s" % ('\n' * nl_before, c(label, lc),
                                         '.' * (self.width - len(label)),
                                         c(str(value), mc), '\n' * nl_after)

        if progress:
            progress.clear()
            self.write(info_line, quiet=quiet)
            progress.update(progress.msg)

        else:
            self.write(info_line, quiet=quiet)


    def info_single(self, message, mc='yellow', nl_before=0, nl_after=0, cut_after=80, level=1, progress=None):
        if isinstance(message, str):
            message = remove_spaces(message)

        if level not in self.single_line_prefixes:
            raise TerminalError("the `info_single` function does not know how to deal with a level of %d :/" % level)

        if cut_after:
            message_line = c("%s%s\n" % (self.single_line_prefixes[level], textwrap.fill(str(message), cut_after)), mc)
        else:
            message_line = c("%s%s\n" % (self.single_line_prefixes[level], str(message)), mc)

        message_line = ('\n' * nl_before) + message_line + ('\n' * nl_after)

        if progress:
            progress.clear()
            self.write(message_line)
            progress.update(progress.msg)

        else:
            self.write(message_line)


    def warning(self, message, header='WARNING', lc='red', raw=False, overwrite_verbose=False, nl_before=0, nl_after=0):
        if isinstance(message, str):
            message = remove_spaces(message)

        message_line = ''
        header_line = c("%s\n%s\n%s\n" % (('\n' * nl_before), header,
                                          '=' * (self.width + 2)), lc)
        if raw:
            message_line = c("%s\n\n%s" % ((message), '\n' * nl_after), lc)
        else:
            message_line = c("%s\n\n%s" % (textwrap.fill(str(message), 80), '\n' * nl_after), lc)

        self.write((header_line + message_line) if message else header_line, overwrite_verbose=overwrite_verbose)


    def store_info_dict(self, destination, strip_prefix=None):
        if strip_prefix:
            # mostly to get rid of output_dir prefix in output file names.
            # surprisingly enough, this is the best place to do it. live
            # and learn :/
            self.info_dict = dictio.strip_prefix_from_dict_values(self.info_dict, strip_prefix)

        dictio.write_serialized_object(self.info_dict, destination)


    def quit(self):
        if self.log_file_path:
            self.log('Bye.')


def pretty_print(n):
    """Pretty print function for very big integers"""
    if not isinstance(n, int):
        return n

    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)


def get_date():
    return time.strftime("%d %b %y %H:%M:%S", time.localtime())


def get_terminal_size():
    """function was taken from http://stackoverflow.com/a/566752"""
    def ioctl_GWINSZ(fd):
        try:
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])
