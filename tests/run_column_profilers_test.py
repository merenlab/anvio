#!/usr/bin/env python
# -*- coding: utf-8

import time

from anvio.columnprofile import ColumnProfile as ColumnProfile_in_C
from anvio.variability import ColumnProfile as ColumnProfile_in_Python

# setup the test factory:
from anvio.variability import VariablityTestFactory
variability_test_class = VariablityTestFactory()

# pretty output
from anvio.summaryhtml import pretty
from anvio.terminal import Run
run = Run(width=55)

# available profiles dict:
ColumnProfile = {'Python': ColumnProfile_in_Python,
                 'C'     : ColumnProfile_in_C}


# test function:
def test(column, reference = 'A', quiet = False, test_class = variability_test_class):
    coverage = len(column)

    results = {}

    if not quiet:
        run.warning('', 'Profiling results for %s nts [TestFactory = %s]' % (pretty(len(column)), 'True' if test_class else 'False'))

    for method in ['C', 'Python']:
        results[method] = {}

        start = time.time()
        results[method]['profile'] = ColumnProfile[method](column, reference = reference, coverage = coverage, pos = 0, test_class = test_class).profile
        end = time.time()
    
        results[method]['delta_time'] = end - start
    
        if not quiet:
            run.info('%s profile' % method, results[method]['profile'])
            run.info('%s response time' % method, results[method]['delta_time'])

    run.info('Result', 'C is ~%.2f times faster' % (results['Python']['delta_time'] / results['C']['delta_time']), mc = 'green')


# tests:

# with and without test class:
column = ('A' * 1000) + ('T' * 1000000)
test(column)

column = ('A' * 1000) + ('T' * 1000000)
test(column, test_class = None)

column = ('A' * 1000000) + ('T' * 1000)
test(column)

column = ('A' * 1000000) + ('T' * 1000)
test(column, test_class = None)

column = ('C' * 1000) + ('T' * 1000)
test(column)

column = ('C' * 1000) + ('T' * 1000) + ('A' * 1)
test(column)
