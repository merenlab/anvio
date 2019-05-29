# -*- coding: utf-8
# pylint: disable=line-too-long
"""Anvi'o - CONCOCT interface for unsupervised clustering

   There are two classes in this file, `CONCOCT` and `CONCOCT_INTERFACE`.
   `CONCOCT` is pretty straightforward to use using anvi'o resources. Here
   is an example:

   >>> import anvio.concoct as concoct
   >>> class Args:
          pass
   >>> args = Args()
   >>> args.profile_db = '/path/to/PROFILE.db'
   >>> args.contigs_db = '/path/to/CONTIGS.db'
   >>> c = concoct.CONCOCT(args)
   >>> c.cluster()
   >>> print c.clusters

   The other class `CONCOCT_INTERFACE`, handles more low-level access
   to the vbgmm module.
"""

try:
    import vbgmm
    __CONCOCT_AVAILABLE__ = True
except ImportError:
    __CONCOCT_AVAILABLE__ = False

import random
import collections
import numpy as np

from sklearn.decomposition import PCA

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tables.collections import TablesForCollections


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class CONCOCT:
    arguments = {
        'seed': (
                ['--seed'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Seed for random numbers"}
                    ),
        'threads': (
                ['-T', '--threads'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Number of threads"}
                    ),
    }

    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        if not __CONCOCT_AVAILABLE__:
            raise ConfigError("CONCOCT is not installed on your system, please install CONCOCT to be \
                               able to use this driver. See https://concoct.readthedocs.io/en/latest/ \
                               for details.")

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.num_clusters_requested = A('num_clusters_requested') or 80

        utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.clusters = {}

        self.lengths = {}
        self.kmers = {}
        self.coverages = {}

        self.debug = anvio.DEBUG

        self.progress.new('Init')

        import anvio.dbops as dbops

        self.progress.update('accessing the profile database ...')
        profile_db = dbops.ProfileDatabase(args.profile_db)

        if not int(profile_db.meta['merged']):
            self.progress.end()
            raise ConfigError('CONCOCT can only be used to cluster merged runs...')

        self.coverages = profile_db.db.get_table_as_dict('mean_coverage_contigs', columns_of_interest=profile_db.samples)
        profile_db.disconnect()

        self.progress.update('accessing the profile database ...')
        contigs_db = dbops.ContigsDatabase(args.contigs_db, quiet=True)
        self.kmers = contigs_db.db.get_table_as_dict('kmer_contigs', keys_of_interest=list(self.coverages.keys()))
        splits_basic_info = contigs_db.db.get_table_as_dict('splits_basic_info', keys_of_interest=list(self.coverages.keys()))
        contigs_db.disconnect()

        self.progress.update('computing split lengths ...')
        for split_name in splits_basic_info:
            self.lengths[split_name] = splits_basic_info[split_name]['length']

        self.progress.end()


    def cluster(self):
        try:
            self.clusters = CONCOCT_INTERFACE(self.kmers, self.coverages, self.lengths, self.debug, NClusters=self.num_clusters_requested).cluster()
        except Exception as e:
            self.run.warning("CONCOCT is upset :/ There will be no CONCOCT binning results for you. Before\
                              anvi'o continues with whatever it was doing before this, here is why CONCOCT\
                              failed in case you want to go after this: '%s'" % e)
            return {}

        # be nice.
        return self.clusters

