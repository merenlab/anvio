# -*- coding: utf-8
# pylint: disable=line-too-long

"""
    Classes to compute completeness estimates based on the information stored in search tables in the
    contigs database.
"""

import numpy
import argparse
from collections import Counter

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.scgdomainclassifier as scgdomainclassifier

from anvio.errors import ConfigError, remove_spaces


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


class Completeness:
    def __init__(self, contigs_db_path, scg_domain_classifier_path=None, source_requested=None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.SCG_comain_predictor = scgdomainclassifier.Predict(argparse.Namespace(), run=terminal.Run(verbose=False), progress=self.progress)

        # hi db
        contigs_db = dbops.ContigsDatabase(contigs_db_path)

        # read info table to get what is available in the db
        info_table = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)

        # identify and remove non-single-copy sources of hmm search results:
        non_singlecopy_sources = set([k for k in list(info_table.keys()) if info_table[k]['search_type'] != 'singlecopy'])
        singlecopy_sources = set([k for k in list(info_table.keys()) if info_table[k]['search_type'] == 'singlecopy'])
        for non_singlecopy_source in non_singlecopy_sources:
            info_table.pop(non_singlecopy_source)

        # get the hmm hits table
        self.hmm_hits_table = contigs_db.db.get_table_as_dict(t.hmm_hits_table_name)

        # read search table (which holds hmmscan hits for splits).
        self.hmm_hits_splits_table = utils.get_filtered_dict(contigs_db.db.get_table_as_dict(t.hmm_hits_splits_table_name), 'source', singlecopy_sources)

        # an example entry in self.hmm_hits_splits_table looks loke this:
        #
        # {
        #    'percentage_in_split'   : 69.6763202725724,
        #    'source'                : u'Campbell_et_al',
        #    'split'                 : u'ANTARCTICAAQUATIC_SMPL_SITE231_3.0UMcontig18439_split_00001',
        #    'hmm_hit_entry_id'      : 1
        # }
        #

        # a little convenience for potential clients:
        self.http_refs = {}
        for source_in_db in info_table:
            self.http_refs[source_in_db] = [h for h in info_table[source_in_db]['ref'].split() if h.startswith('http')][0]

        self.genes_in_db = dict([(s, info_table[s]['genes'].split(', ')) for s in info_table])

        # we're done with the db
        contigs_db.disconnect()

        self.sources = list(info_table.keys())
        self.domains = set([info_table[source]['domain'] for source in self.sources])
        self.source_to_domain = dict([(source, info_table[source]['domain']) for source in self.sources])
        self.domain_to_sources = [(domain, [source for source in self.sources if info_table[source]['domain'] == domain]) for domain in self.domains]

        self.missing_SCG_HMMs = [s for s in self.SCG_comain_predictor.SCG_sources if s not in self.sources]
        if len(self.missing_SCG_HMMs):
            self.run.warning("Sad news :( Anvi'o completion estimations require all single-copy core gene (SCG) collections to be run on contigs\
                              databases. Yet your contigs database '%s' is lacking these default ones: '%s'. If you need completion estimates,\
                              the easiest solution is to run `anvi-run-hmms` on your contigs database. This will make sure all SCG HMM hits are \
                              stored in your contigs database. Alternatively you can run the same program with `-I XXX` parameter where XXX is one\
                              of the missing SCG HMM collection (in case you are only missing one of the SCG collections)." \
                                                                                    % (contigs_db_path, ', '.join(self.missing_SCG_HMMs)))

        if source_requested:
            if source_requested not in self.sources:
                raise ConfigError('Requested source "%s" is not one of the single-copy gene sources found in the database.' % source_requested)

            # filter out sources that are not requested
            self.sources = [source_requested]
            self.genes_in_db = {source_requested: self.genes_in_db[source_requested]}
            self.hmm_hits_splits_table = utils.get_filtered_dict(self.hmm_hits_splits_table, 'source', set([source_requested]))

        self.unique_gene_id_to_gene_name = {}
        self.splits_unique_gene_id_occurs = {}
        # these will be very useful later. trust me.
        for entry in list(self.hmm_hits_splits_table.values()):
            hmm_hit = self.hmm_hits_table[entry['hmm_hit_entry_id']]
            gene_unique_identifier = hmm_hit['gene_unique_identifier']

            if gene_unique_identifier not in self.unique_gene_id_to_gene_name:
                self.unique_gene_id_to_gene_name[gene_unique_identifier] = hmm_hit['gene_name']

            if gene_unique_identifier not in self.splits_unique_gene_id_occurs:
                self.splits_unique_gene_id_occurs[gene_unique_identifier] = [entry['split']]
            else:
                self.splits_unique_gene_id_occurs[gene_unique_identifier].append(entry['split'])


    def list_hmm_sources(self):
        self.run.warning('', 'HMM SOURCES FOUND', lc='yellow')
        for source in self.sources:
            self.run.info_single(source)


    def get_average_domain_completion_and_redundancy(self, d, domain):
        """For a given results dict `d` obtained from 'get_info_for_splits', and a domain, returns
           average percent completion and redundancy for the domain."""

        percent_completion = numpy.mean([d[domain][s]['percent_completion'] for s in d[domain]])
        percent_redundancy = numpy.mean([d[domain][s]['percent_redundancy'] for s in d[domain]])

        return percent_completion, percent_redundancy


    def get_best_matching_domain(self, d, observed_genes_per_domain):
        """Returns the best matcing domain by using a random forest classifier.

           The input dict `d` has no role in predicting the best matching domain,
           but useful to print out some helpful messages when necessary.

           It returns a tuple for best matching domain and how confident is the match
           according to the random forest classifier.
        """

        # learn everything from the random forest.
        domain_predictions, prob_mixed_domains, prob_blank_domain = self.SCG_comain_predictor.predict_from_observed_genes_per_domain(observed_genes_per_domain)
        domain_specific_estimates = []

        if anvio.DEBUG:
            self.run.warning(None, header="SCG DATA FOR C/R ESTIMTES", lc='green')
            self.run.info_single("Probability BLANK: %.2f" % (prob_blank_domain), mc='cyan')
            self.run.info_single("Probability COMPOSITE  %.2f" % (prob_mixed_domains), mc='cyan')

        domain_cr = {}
        for confidence, domain in domain_predictions:
            source = self.SCG_comain_predictor.SCG_domain_to_source[domain]
            domain_cr[domain] = {'percent_completion': d[domain][source]['percent_completion'],
                                 'percent_redundancy': d[domain][source]['percent_redundancy']}

            if anvio.DEBUG:
                self.run.info_single("Domain '%8s' (probability: %.2f) C/R: %.2f/%.2f" % (domain, confidence, domain_cr[domain]['percent_completion'], domain_cr[domain]['percent_redundancy']), mc='green')

        domain_matching_confidence, best_matching_domain = domain_predictions[0]
        best_mathcing_domain_completion, best_matching_domain_redundancy = domain_cr[domain]['percent_completion'], domain_cr[domain]['percent_redundancy']

        # it he probabiity of this domain to be blank is so high, and the best matching domain confidence
        # is so low, do not return any completion or best matching domain estimates.
        if domain_matching_confidence < 0.15 and prob_blank_domain > 0.5:
            domain_matching_confidence, best_matching_domain = 0.0, None
            best_mathcing_domain_completion, best_matching_domain_redundancy = None, None

        info_text = ''
        max_confidence = max([t[0] for t in domain_predictions])
        if prob_blank_domain > 0.5:
            info_text = "CRAP DOMAIN EST BECAUSE NO SIGNAL. So anvi'o is having hard time determining any domain\
                         for this selection either because the number of contigs are very little, or there are no\
                         SCGs among the selected ones. This may happen if you are working with genomes that are\
                         extremely low completion, or alternatively coming from parts of life that are very \
                         understudied (such as viruses or plasmids, etc). If these do not apply to you, and you are\
                         sure your selection represents a proper genome, then either anvi'o made a mistake, or you\
                         stumbled upon a graet story."
        elif max_confidence < 0.6 and prob_mixed_domains > 0.4:
            info_text = "CRAP DOMAIN EST BECAUSE STUFF IS MIXED. Please note that anvi'o determined '%s' \
                         as the best matching domain for your contigs BUT actually the probability of these \
                         contigs to be coming from mixed domains is crazy high (%.2f). Which means, neither \
                         this domain prediction, nor the completion and redundancy estimates should mean much,\
                         and you should take a look at the entire list of C/R estimates from all domain SCGs :/\
                         The good news is that more refined the input contigs (such as you make more and more\
                         precise selections or provdide more and more refined genomes), this situation will\
                         likely correct itself." % (best_matching_domain, prob_mixed_domains)
        elif max_confidence < 0.6 and prob_mixed_domains < 0.4:
            info_text = "CRAP DOMAIN EST BECAUSE WHO KNOWS WHY. Please note that anvi'o determined '%s' as the\
                         best matching domain for your contigs to predict the completion and redundancy\
                         estimates through sigle-copy core genes, however, since the confidence is as low\
                         as %.2f, you should take this estimate with a grain of salt. This low confidence\
                         is most likely due to a very small number of contigs to offer reliable estimates\
                         of domain." % (best_matching_domain, domain_matching_confidence)
        else:
            if best_matching_domain_redundancy < 10:
                info_text = "GREAT DOMAIN EST. YOU'RE GOLDEN. The very high confidence of random forest indicates that this set of contigs\
                             are almost certainly coming from a population that belongs to domain %s. IN ADDITION, the low redundancy of SCGs\
                             do not predict any serious contamination. But please remember that this information does not mean there is NO\
                             contamination in your genome bin." \
                                        % (best_matching_domain)
            elif best_matching_domain_redundancy > 10 and best_matching_domain_redundancy < 100:
                info_text = "GREAT DOMAIN EST. YOU'RE GOLDEN. The very high confidence of random forest indicates that this set of contigs\
                             are almost certainly coming from a population that belongs to domain %s. BUT you almost certainly are looking at\
                             a composite genome. You should consider refining this particular collection of contigs to lower the redundancy." \
                                        % (best_matching_domain)
            elif best_matching_domain_redundancy > 100:
                info_text = "GREAT DOMAIN EST. YOU'RE GOLDEN. The very high confidence of random forest indicates that the very large fraction\
                             of this set of contigs are almost certainly coming from populations that belong to the domain %s. HOWEVER, the\
                             extremely high amount of redundancy of SCGs in domain %s suggests that you either are working with a set of contigs\
                             that are extremely composite, or you are looking basically an entire metagenome or something." \
                                        % (best_matching_domain, best_matching_domain)
            else:
                info_text = "GREAT DOMAIN EST. YOU'RE GOLDEN. But your redundancy estimates are weird :("

        if anvio.DEBUG:
            self.run.warning(info_text)

        return (best_matching_domain, domain_matching_confidence, remove_spaces(info_text))


    def get_info_for_splits(self, split_names, min_e_value=1e-5):
        """This function takes a bunch of split names, and returns three things:

            - Average percent completion for best matching domain
            - Average redundancy for best matching domain
            - Best matching domain for this collection of splits,
            - Domain matching confidence (see get_average_domain_completion_and_redundancy for details)
            - And a comprehensive results dictionary that explains each HMM source in each domain,

        For your convenience, you can call this function this way:

        p_completion, p_redundancy, domain, domain_confidence, results_dict = get_info_for_splits(s)
        """

        hmm_hits_splits_table = utils.get_filtered_dict(self.hmm_hits_splits_table, 'split', split_names)

        # FIXME: the design here is turning into a bad case of spaghetti code. we should reimplement
        # the SCG / completion stuff around the random forest domain predictor. the previous code is
        # quite inefficient, and the late addition of random forest domain predictor is making things
        # even less clear. The following dictionary is to predict the domain:
        observed_genes_per_domain = {}
        for domain in self.SCG_comain_predictor.SCG_domains:
            observed_genes_per_domain[domain] = Counter()

        # we need to restructure 'hits' into a dictionary that gives access to sources and genes in a more direct manner
        info_dict, gene_name_to_unique_id = {}, {}
        for source in self.sources:
            info_dict[source], gene_name_to_unique_id[source] = {}, {}

        # here we go through every hit and populate 'info_dict' and 'gene_name_to_unique_id':
        for entry in list(hmm_hits_splits_table.values()):
            hmm_hit = self.hmm_hits_table[entry['hmm_hit_entry_id']]

            if hmm_hit['e_value'] > min_e_value:
                continue

            source = hmm_hit['source']
            domain = self.source_to_domain[source]
            e_value = hmm_hit['e_value']
            gene_name = hmm_hit['gene_name']
            percentage = entry['percentage_in_split']
            gene_unique_id = hmm_hit['gene_unique_identifier']

            if gene_unique_id in info_dict[source]:
                info_dict[source][gene_unique_id]['percentage'] += percentage
            else:
                info_dict[source][gene_unique_id] = {}
                info_dict[source][gene_unique_id] = {'gene_name': gene_name, 'percentage': percentage, 'e_value': e_value}

            if gene_name in gene_name_to_unique_id[source]:
                gene_name_to_unique_id[source][gene_name].add(gene_unique_id)
            else:
                gene_name_to_unique_id[source][gene_name] = set([gene_unique_id])

            observed_genes_per_domain[domain][gene_name] += 1

        # here we generate the results information
        results_dict = {}
        for domain in self.domains:
            results_dict[domain] = {}

        for source in self.sources:
            domain = self.source_to_domain[source]
            results_dict[domain][source] = {'domain': domain, 'source': source}

            genes_count = Counter([v['gene_name'] for v in list(info_dict[source].values())])

            # report num genes in the model and the num of those with hits (note that htis doesn't
            # care whether those hits are contributing to redundance or not --instad here we are
            # intrested only in the 'coverage' of the model)
            results_dict[domain][source]['num_genes_in_model'] = len(self.genes_in_db[source])
            results_dict[domain][source]['num_genes_in_model_with_hits ']= len(genes_count)
            results_dict[domain][source]['model_coverage']= len(genes_count) / len(self.genes_in_db[source])

            results_dict[domain][source]['percent_completion'] = len(genes_count) * 100.0 / len(self.genes_in_db[source])

            # report redundancy:
            genes_that_occur_multiple_times = [g for g in genes_count if genes_count[g] > 1]
            results_dict[domain][source]['percent_redundancy'] = sum([genes_count[g] - 1 for g in genes_that_occur_multiple_times]) * 100.0 / len(self.genes_in_db[source])

            # identify splits that contribute the same single_copy_gene
            redundants = {}
            for gene_name in genes_that_occur_multiple_times:
                redundants[gene_name] = [self.splits_unique_gene_id_occurs[unique_gene_id] for unique_gene_id in gene_name_to_unique_id[source][gene_name]]
            results_dict[domain][source]['redundants'] = redundants

        if not len(results_dict):
            return (None, None, None, None, results_dict)

        best_matching_domain, domain_matching_confidence, info_text = self.get_best_matching_domain(results_dict, observed_genes_per_domain)

        if best_matching_domain:
            percent_completion, percent_redundancy = self.get_average_domain_completion_and_redundancy(results_dict, best_matching_domain)
        else:
            percent_completion, percent_redundancy = 0.0, 0.0

        return (percent_completion, percent_redundancy, best_matching_domain, domain_matching_confidence, results_dict)
