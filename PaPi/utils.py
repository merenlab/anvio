# -*- coding: utf-8
#
# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys
import time
import fcntl
import string
import struct
import cPickle
import termios 
import itertools

from PaPi.constants import pretty_names


complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV',\
                               'tgcayrkmvhdbTGCAYRKMVHDB')


def rev_comp(seq):
    return seq.translate(complements)[::-1]


class KMers:
    def __init__(self):
        self.kmers = {}
        
        self.get_kmers(k = 4)

    def get_kmers(self, k):
        arg = ['ATCG'] * k
        kmers = set()
        
        for item in itertools.product(*arg):
            kmer = ''.join(item)
            if rev_comp(kmer) not in kmers:
                kmers.add(kmer)
        
        self.kmers[k] = kmers


    def get_kmer_frequency(self, sequence, k = 4):
        sequence = sequence.upper()

        if len(sequence) < k:
            return None

        if not self.kmers.has_key(k):
            self.get_kmers(k)
        
        kmers = self.kmers[k]
        frequencies = dict(zip(kmers, [0] * len(kmers)))
        
        for i in range(0, len(sequence) - (k - 1)):
            kmer = sequence[i:i + k]
            
            # FIXME: this can be faster/better
            if len([n for n in kmer if n not in 'ATCG']):
                continue

            if frequencies.has_key(kmer):
                frequencies[kmer] += 1
            else:
                frequencies[rev_comp(kmer)] += 1
           
        return frequencies


class Progress:
    def __init__(self):
        self.pid = None
        self.verbose = True
        self.terminal_width = None

        self.get_terminal_width()
        self.color_prefix = '\033[0;30m\033[46m'
        self.color_postfix = '\033[0m'
        
        self.currently_shown = None


    def get_terminal_width(self):
        try:
            self.terminal_width = get_terminal_size()[0]
        except:
            self.terminal_width = 80 


    def new(self, pid):
        if self.pid:
            raise LibError, "Progress.new() can't be called before ending the previous one (Existing: '%s', Competing: '%s')." % (self.pid, pid)

        if not self.verbose:
            return

        self.pid = '%s %s' % (get_date(), pid)
        self.get_terminal_width()
        self.currently_shown = None


    def write(self, c):
        surpass = self.terminal_width - len(c)
        
        if surpass < 0:
            c = c[0:-(-surpass + 4)] + ' (...)'
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
    """a class that keeps info about an oligotyping run, and deal with the console output"""
    def __init__(self, info_file_path = None, verbose = True):
        if info_file_path:
            self.init_info_file_obj(info_file_path)
        else:
            self.info_file_obj = None

        self.info_dict = {}
        self.verbose = verbose


    def init_info_file_obj(self, info_file_path):
            self.info_file_obj = open(info_file_path, 'w')


    def info(self, key, value, quiet = False):
        self.info_dict[key] = value
        
        if quiet:
            return True
        
        if type(value) == int:
            value = pretty_print(value)

        label = get_pretty_name(key)

        info_line = "%s %s: %s\n" % (label, '.' * (65 - len(label)), str(value))
        if self.info_file_obj:
            self.info_file_obj.write(info_line)

        if self.verbose:
            sys.stderr.write(info_line)


    def store_info_dict(self, destination):
        cPickle.dump(self.info_dict, open(destination, 'w'))


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


def get_pretty_name(key):
    if pretty_names.has_key(key):
        return pretty_names[key]
    else:
        return key


class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Config Error: %s' % self.e
