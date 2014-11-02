# -*- coding: utf-8 -*-

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

IS_ESSENTIAL_FIELD = lambda f: (not f.startswith('__')) and (f not in ["contigs", "GC_content", "length"])
IS_AUXILIARY_FIELD = lambda f: f.startswith('__')

pretty_names = {}
