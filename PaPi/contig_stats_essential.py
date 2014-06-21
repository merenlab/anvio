#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import numpy


class Essential:
    def __init__(self, reference, pileup):
        """Gets a pileup, returns simple stats"""
        self.reference = reference
        self.pileup = pileup
        self.basics = {}

    def report(self):
        coverage = [pileupcolumn.n for pileupcolumn in self.pileup]

        if not coverage:
            self.basics['coverage'] = [] 
            self.basics['min_coverage'] = 0
            self.basics['max_coverage'] = 0
            self.basics['median_coverage'] = 0
            self.basics['mean_coverage'] = 0
            self.basics['std_coverage'] = 0
        else:
            self.basics['coverage'] = coverage
            self.basics['min_coverage'] = numpy.min(coverage)
            self.basics['max_coverage'] = numpy.max(coverage)
            self.basics['median_coverage'] = numpy.median(coverage)
            self.basics['mean_coverage'] = numpy.mean(coverage)
            self.basics['std_coverage'] = numpy.std(coverage)

        return self.basics