# -*- coding: utf-8

import os
import sys
import pkg_resources

import anvio.tables as t
from anvio.terminal import Run


run = Run()


def set_version():
    try:
        __version__ = pkg_resources.require("anvio")[0].version
    except:
        # maybe it is not installed but being run from the codebase dir?
        try:
            __version__ = open(os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/../VERSION').read().strip()
        except:
            __version__ = 'unknown'

    return __version__, t.annotation_db_version, t.profile_db_version


def print_version():
    run.info("Anvi'o version", __version__, mc = 'green')
    run.info("Annotation DB version", __annotation__version__)
    run.info("Profile DB version", __profile__version__)


__version__, __annotation__version__, __profile__version__ = set_version()


if '-v' in sys.argv or '--version' in sys.argv:
    print_version()
    sys.exit()
