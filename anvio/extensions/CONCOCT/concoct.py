import sys, getopt
import cPickle
import numpy as np
import vbgmm

from sklearn.decomposition import PCA

class Concoct():
    def __init__(self, kmers, coverages, NClusters = 80, kmer_length = 4, read_length = 100, bNormaliseByContig = True, nc = 0.90):
        #calc number of samples
        first_cov = coverages.itervalues().next()
        sample_names = first_cov.keys()
        NS = len(sample_names)
    
        #calc number of kmers
        first_kmer = kmers.itervalues().next()
        kmer_names = first_kmer.keys()
        NK = len(kmer_names)
    
        #number of contigs to cluster
        contig_names = coverages.keys()
        self.NC = len(coverages.keys())

        cov_array = np.zeros((self.NC,NS))
        kmer_array = np.zeros((self.NC,NK))
    
        p = 0
        for k in contig_names:
            vk = kmers[k]
            valuesk = [vk[x] for x in kmer_names]
            kmer_array[p,:] = valuesk[:]
    
            vc = coverages[k]
            valuesc = [vc[x] for x in sample_names]
            cov_array[p,:] = valuesc[:]
            p = p + 1
    
        #this is not really valid since split contigs have inherited composition in Anvio
        contig_lengths = kmer_array.sum(axis=1) + kmer_length - 1
        kmer_array = kmer_array + np.ones((self.NC,NK))
        kmer_row_sums = kmer_array.sum(axis=1)
        kmer_array = np.log(kmer_array / kmer_row_sums[:, np.newaxis])

        cov_array = cov_array + (read_length/contig_lengths)[:,np.newaxis]
        
        #normalise by Sample maybe weight by contig length in future
        cov_col_sums = cov_array.sum(axis=0)
        cov_array = cov_array/cov_col_sums[np.newaxis,:]
        
        if bNormaliseByContig:
            cov_row_sums = cov_array.sum(axis=1)
            cov_array = cov_array/cov_row_sums[:,np.newaxis]
            cov_row_sums=cov_row_sums.reshape((cov_array.shape[0],1))
            cov_array = np.append(cov_array,cov_row_sums,1)
        
        cov_array = np.log(cov_array)
            
        #join log transformed composition and coverage together
        self.original_data = np.concatenate((kmer_array,cov_array),axis=1)
        
        #perform PCA
        pca_object = PCA(n_components=nc).fit(self.original_data)
        self.transformed_data = pca_object.transform(self.original_data)
    
        self.NClusters = NClusters
        self.assign = np.zeros((self.NC),dtype=np.int32)
    
    def cluster(self): 
        vbgmm.fit(self.transformed_data,self.assign,self.NClusters)
        
        return self.assign
        
def main(args):

    kmers = cPickle.load(open('kmers.cPickle'))
    coverages = cPickle.load(open('coverages.cPickle'))

    import ipdb; ipdb.set_trace()

    #import concoct
    c = Concoct(kmers, coverages)
    concoct_clustering = c.cluster()

if __name__ == "__main__":
    main(sys.argv[1:])

