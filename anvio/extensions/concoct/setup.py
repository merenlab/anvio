#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("vbgmm",
                             sources=["vbgmm.pyx", "c_vbgmm_fit.c"],
			     libraries=['gsl',  'gslcblas'],
                             include_dirs=[numpy.get_include(), '/opt/local/include/'])],
)
