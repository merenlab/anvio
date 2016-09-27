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

import random
import collections
import numpy as np

from sklearn.decomposition import PCA

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.vbgmm as vbgmm
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Christopher Quince"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class CONCOCT:
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.num_clusters_requested = A('num_clusters_requested') or 80

        dbops.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.clusters = {}

        self.lengths = {}
        self.kmers = {}
        self.coverages = {}

        self.debug = args.debug

        self.progress.new('Init')

        self.progress.update('accessing the profile database ...')
        profile_db = dbops.ProfileDatabase(args.profile_db)

        if not int(profile_db.meta['merged']):
            self.progress.end()
            raise ConfigError, 'CONCOCT can only be used to cluster merged runs...'

        self.coverages = profile_db.db.get_table_as_dict('mean_coverage_contigs', columns_of_interest=profile_db.samples)
        profile_db.disconnect()

        self.progress.update('accessing the profile database ...')
        contigs_db = dbops.ContigsDatabase(args.contigs_db, quiet=True)
        self.kmers = contigs_db.db.get_table_as_dict('kmer_contigs', keys_of_interest=self.coverages.keys())
        splits_basic_info = contigs_db.db.get_table_as_dict('splits_basic_info', keys_of_interest=self.coverages.keys())
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


    def store_clusters_as_TAB_delimited_text(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        self.progress.new('Storing clusters as TAB-delimited file')
        self.progress.update('creating the clusters dictionary ...')
        clusters_dict = {}
        for contig_name in self.clusters:
            clusters_dict[contig_name] = {'concoct_bin': self.clusters[contig_name]}

        self.progress.update('writing the file ...')
        utils.store_dict_as_TAB_delimited_file(clusters_dict, output_file_path, ['contig', 'concoct_bin'])
        self.progress.end()

        self.run.info('CONCOCT results in txt', output_file_path, display_only=True)


    def store_clusters_in_db(self, collection_name='CONCOCT'):
       # convert id -> bin mapping dict into a bin -> ids dict
        data = {}
        bin_info_dict = {}

        if not len(self.clusters):
            self.run.info('CONCOCT results in db', 'Nope. CONCOCT clusters are empty. Skipping!', mc='red', display_only=True)
            return

        for split_name in self.clusters:
            bin_id = self.clusters[split_name]
            if bin_id in data:
                data[bin_id].add(split_name)
            else:
                data[bin_id] = set([split_name])
                bin_info_dict[bin_id] = {'html_color': '#' + ''.join(['%02X' % random.randint(50, 230) for i in range(0, 3)]), 'source': 'CONCOCT'}
                                                            # ^
                                                            #  \
                                                            #    poor man's random color generator

        c = dbops.TablesForCollections(self.profile_db_path)
        c.append(collection_name, data, bin_info_dict)

        self.run.info('CONCOCT results in db', self.profile_db_path, display_only=True)


class CONCOCT_INTERFACE():
    def __init__(self, kmers, coverages, lengths, debug, NClusters=80, kmer_length=4, read_length=100, bNormaliseByContig=True, nc=0.90, r=run, p=progress):
        self.run = r
        self.progress = p

        if NClusters < 2 or NClusters > 2000:
            raise ConfigError, "Number of clusters requested must be between 2 and 2,000."

        self.progress.new('CONCOCT')

        self.debug = debug

        self.progress.update('Checking the number of samples ...')
        first_cov = coverages.itervalues().next()
        sample_names = first_cov.keys()
        NS = len(sample_names)

        self.progress.update('Checking the number of k-mers ...')
        first_kmer = kmers.itervalues().next()
        kmer_names = first_kmer.keys()
        NK = len(kmer_names)

        self.progress.update('Checking the number of contigs to cluster ...')
        self.contig_names = coverages.keys()
        self.NC = len(coverages.keys())

        cov_array = np.zeros((self.NC, NS))
        kmer_array = np.zeros((self.NC, NK))

        self.progress.update('Setting up kmer and coverage arrays ...')
        p = 0
        for k in self.contig_names:
            vk = kmers[k]
            valuesk = [vk[x] for x in kmer_names]
            kmer_array[p, :] = valuesk[:]

            vc = coverages[k]
            valuesc = [vc[x] for x in sample_names]
            cov_array[p, :] = valuesc[:]
            p = p + 1

        #this is not really valid since split contigs have inherited composition in Anvio
        contig_lengths = kmer_array.sum(axis=1) + kmer_length - 1
        kmer_array = kmer_array + np.ones((self.NC, NK))
        kmer_row_sums = kmer_array.sum(axis=1)
        kmer_array = np.log(kmer_array / kmer_row_sums[:, np.newaxis])

        cov_array = cov_array + (read_length / contig_lengths)[:, np.newaxis]

        #normalise by Sample maybe weight by contig length in future
        cov_col_sums = cov_array.sum(axis=0)
        cov_array = cov_array / cov_col_sums[np.newaxis, :]

        if bNormaliseByContig:
            cov_row_sums = cov_array.sum(axis=1)
            cov_array = cov_array / cov_row_sums[:, np.newaxis]
            cov_row_sums = cov_row_sums.reshape((cov_array.shape[0], 1))
            cov_array = np.append(cov_array, cov_row_sums, 1)

        cov_array = np.log(cov_array)

        self.progress.update('Joining coverage and composition arrays ...')
        self.original_data = np.concatenate((kmer_array, cov_array), axis=1)

        self.progress.update('Performing PCA ...')
        pca_object = PCA(n_components=nc).fit(self.original_data)
        self.transformed_data = pca_object.transform(self.original_data)

        self.NClusters = NClusters
        self.assign = np.zeros((self.NC), dtype=np.int32)

        self.progress.end()
        self.run.info('CONCOCT INIT', 'Complete for %d splits' % len(self.contig_names))


    def cluster(self):
        self.run.info('CONCOCT VGBMM', 'Requesting %d clusters' % (self.NClusters))

        self.progress.new('VBGMM')
        self.progress.update('Clustering ...')
        vbgmm.fit(self.transformed_data, self.assign, self.NClusters, self.debug)
        self.progress.end()

        self.run.info('CONCOCT VGBMM', 'Returning %d final clusters' % (len(set(self.assign))))

        # construct and return a results dictionary:
        bin_name_conversion_dict = {}
        bin_counter = 1
        for bin_id, num_contigs in collections.Counter(self.assign).most_common():
            bin_name_conversion_dict[bin_id] = 'Bin_%d' % bin_counter
            bin_counter += 1

        return dict(zip(self.contig_names, [bin_name_conversion_dict[g] for g in self.assign]))
