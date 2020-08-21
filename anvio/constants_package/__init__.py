# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


DB_TYPES = ['contigs',
            'pan',
            'profile',
            'genes',
            'auxiliary data for coverages',
            'genomestorage',
            'structure',
            'keggmodules',
            'trnaseq']

VERSIONS_FOR_DB_TYPES = {'contigs': '17',
                         'profile': '34',
                         'genes': '6',
                         'pan': '14',
                         'structure': '2',
                         'genomestorage': '7',
                         'auxiliary data for coverages': '2',
                         'trnaseq': '1',
                         'config': '1',
                         'keggmodules': '2'}