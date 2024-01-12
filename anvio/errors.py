# -*- coding: utf-8
# pylint: disable=line-too-long

"""Exceptions"""

import sys
import textwrap
import traceback

import anvio
from anvio.ttycolors import color_text

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


def remove_spaces(text):
    if not text:
        return ""

    while True:
        if text.find("  ") > -1:
            text = text.replace("  ", " ")
        else:
            break

    return text


class AnvioError(Exception, object):
    def __init__(self, e=None):
        Exception.__init__(self)
        return

    def __str__(self):
        max_len = max([len(l) for l in textwrap.fill(textwrap.dedent(self.e), 80).split('\n')])
        error_lines = ['%s%s' % (l, ' ' * (max_len - len(l))) for l in textwrap.fill(textwrap.dedent(self.e), 80).split('\n')]

        error_message = ['%s: %s' % (color_text(self.error_type, 'red'), error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(self.error_type) + 2), error_line))

        if anvio.DEBUG:
            exc_type, exc_value, exc_traceback = sys.exc_info()

            sep = color_text('=' * 80, 'red')

            print(color_text('\nTraceback for debugging', 'red'))
            print(sep)
            traceback.print_tb(exc_traceback, limit=100, file=sys.stdout)
            print(sep)

        return '\n\n' + '\n'.join(error_message) + '\n\n'


    def clear_text(self):
        return self.e


class CommandError(AnvioError):
    """Use this when a command (e.g., something run with utils.run_command) fails."""

    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Command Error'
        AnvioError.__init__(self)


class ConfigError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Config Error'
        AnvioError.__init__(self)


class StupidHMMError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Stupid HMM Error'
        AnvioError.__init__(self)


class GenesDBError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Genes DB Error'
        AnvioError.__init__(self)


class UpgradeError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Upgrade Error'
        AnvioError.__init__(self)


class RefineError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Refine Error'
        AnvioError.__init__(self)


class TerminalError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Terminal Error'
        AnvioError.__init__(self)


class FilesNPathsError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'File/Path Error'
        AnvioError.__init__(self)


class DictIOError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Dict IO Error'
        AnvioError.__init__(self)


class HDF5Error(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'HDF5 Error'
        AnvioError.__init__(self)

class AuxiliaryDataError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Auxiliary Data Error'
        AnvioError.__init__(self)


class AnviServerError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Anvio Server Error'
        AnvioError.__init__(self)


class ModellerError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Modeller Error'
        AnvioError.__init__(self)


class ModellerScriptError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Modeller Script Error'
        AnvioError.__init__(self)


class TRNAIdentifierError(AnvioError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'tRNA Identifier Error'
        AnvioError.__init__(self)
