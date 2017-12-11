# -*- coding: utf-8
# pylint: disable=line-too-long

"""TTY colors"""

import sys


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


tty_colors = {
    'gray'    :{'normal': '\033[1;30m%s\033[1m', 'bold': '\033[0;30m%s\033[0m'},
    'red'     :{'normal': '\033[1;31m%s\033[1m', 'bold': '\033[0;31m%s\033[0m'},
    'green'   :{'normal': '\033[1;32m%s\033[1m', 'bold': '\033[0;32m%s\033[0m'},
    'yellow'  :{'normal': '\033[1;33m%s\033[1m', 'bold': '\033[0;33m%s\033[0m'},
    'blue'    :{'normal': '\033[1;34m%s\033[1m', 'bold': '\033[0;34m%s\033[0m'},
    'magenta' :{'normal': '\033[1;35m%s\033[1m', 'bold': '\033[0;35m%s\033[0m'},
    'cyan'    :{'normal': '\033[1;36m%s\033[1m', 'bold': '\033[0;36m%s\033[0m'},
    'white'   :{'normal': '\033[1;37m%s\033[1m', 'bold': '\033[0;37m%s\033[0m'},
    'crimson' :{'normal': '\033[1;38m%s\033[1m', 'bold': '\033[0;38m%s\033[0m'}
}


def color_text(text, color="crimson", weight="bold"):
    if sys.stdout.isatty():
        return tty_colors[color][weight] % text
    else:
        return text
