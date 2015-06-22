# -*- coding: utf-8
"""Relations with the console output, Progress and Run classes"""

import os
import sys
import time
import fcntl
import struct
import termios
import textwrap

import anvio.constants as constants
import anvio.dictio as dictio

from anvio.errors import TerminalError
from anvio.ttycolors import color_text as c

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class SuppressAllOutput(object):
    def __enter__(self):
        sys.stderr.flush()
        self.old_stderr = sys.stderr
        sys.stderr = open('/dev/null', 'a+', 0)
        sys.stdout.flush()
        self.old_stdout = sys.stdout
        sys.stdout = open('/dev/null', 'a+', 0)
 
    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.flush()
        sys.stderr = self.old_stderr
        sys.stdout.flush()
        sys.stdout = self.old_stdout


def remove_spaces(text):
    while 1:
        if text.find("  ") > -1:
            text = text.replace("  ", " ")
        else:
            break

    return text


class Progress:
    def __init__(self, verbose = True):
        self.pid = None
        self.verbose = verbose
        self.terminal_width = None
        self.is_tty = sys.stdout.isatty()

        self.get_terminal_width()
        self.color_prefix = '\033[0;30m\033[46m'
        self.color_postfix = '\033[0m'
        
        self.currently_shown = None


    def get_terminal_width(self):
        try:
            self.terminal_width = get_terminal_size()[0]
        except:
            self.terminal_width = 120


    def new(self, pid):
        if self.pid:
            raise TerminalError, "Progress.new() can't be called before ending the previous one (Existing: '%s', Competing: '%s')." % (self.pid, pid)

        if not self.verbose:
            return

        self.pid = '%s %s' % (get_date(), pid)
        self.get_terminal_width()
        self.currently_shown = None


    def write(self, c):
        surpass = self.terminal_width - len(c)
        
        if surpass < 0:
            c = c[0:-(-surpass + 5)] + ' (...)'
        else:
            self.currently_shown = c
            c = c + ' ' * surpass

        if self.verbose:
            sys.stderr.write(self.color_prefix + c + self.color_postfix)
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
        self.currently_shown = None


    def append(self, msg):
        if not self.verbose:
            return
        self.write('%s%s' % (self.currently_shown, msg))


    def update(self, msg):
        if not self.verbose:
            return
        self.clear()
        self.write('\r[%s] %s' % (self.pid, msg))


    def end(self):
        self.pid = None
        if not self.verbose:
            return
        self.clear()


class Run:
    def __init__(self, info_file_path = None, verbose = True, width = 45):
        if info_file_path:
            self.init_info_file_obj(info_file_path)
        else:
            self.info_file_obj = None

        self.info_dict = {}
        self.verbose = verbose
        self.width = width


    def init_info_file_obj(self, info_file_path):
            self.info_file_obj = open(info_file_path, 'w')


    def info(self, key, value, quiet = False, display_only = False, lc = 'cyan', mc = 'yellow'):
        if not display_only:
            self.info_dict[key] = value

        if quiet:
            return True

        if type(value) == str:
            value = remove_spaces(value)
        if type(value) == int:
            value = pretty_print(value)

        label = constants.get_pretty_name(key)

        info_line = "%s %s: %s\n" % (c(label, lc), '.' * (self.width - len(label)), c(str(value), mc))

        if self.info_file_obj:
            self.info_file_obj.write(info_line)

        if self.verbose:
            sys.stderr.write(info_line)


    def info_single(self, message, mc = 'yellow', nl_before = 0, nl_after = 0, cut_after = 80):
        if type(message) == str:
            message = remove_spaces(message)

        if cut_after:
            message_line = c("* %s\n" % (textwrap.fill(str(message), cut_after)), mc)
        else:
            message_line = c("* %s\n" % str(message), mc)

        if self.verbose:
            sys.stderr.write('\n' * nl_before)
            sys.stderr.write(message_line)
            sys.stderr.write('\n' * nl_after)


    def warning(self, message, header='WARNING', lc = 'red', raw = False):
        if type(message) == str:
            message = remove_spaces(message)

        header_line = c("\n%s\n%s\n" % (header, '=' * (self.width + 2)), lc)
        if raw:
            message_line = c("%s\n\n" % (message), lc)
        else:
            message_line = c("%s\n\n" % (textwrap.fill(str(message), 80)), lc)

        if self.verbose:
            sys.stderr.write(header_line)
            if message:
                sys.stderr.write(message_line)


    def store_info_dict(self, destination, strip_prefix = None):

        if strip_prefix:
            # mostly to get rid of output_dir prefix in output file names.
            # surprisingly enough, this is the best place to do it. live 
            # and learn :/
            self.info_dict = dictio.strip_prefix_from_dict_values(self.info_dict, strip_prefix)

        dictio.write_serialized_object(self.info_dict, destination)


    def quit(self):
        if self.info_file_obj:
            self.info_file_obj.close()


def pretty_print(n):
    """Pretty print function for very big integers"""
    if type(n) != int:
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
