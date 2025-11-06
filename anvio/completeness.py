# -*- coding: utf-8
# pylint: disable=line-too-long

"""
    Classes to compute completeness estimates based on the information stored in search tables in the
    contigs database.
"""

import argparse
from collections import Counter

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.scgdomainclassifier as scgdomainclassifier

from anvio.errors import ConfigError, remove_spaces

with terminal.SuppressAllOutput():
    import anvio.data.hmm as hmm_data


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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
        self.initialized_properly = True

        self.SCG_domain_predictor = scgdomainclassifier.Predict(argparse.Namespace(), run=terminal.Run(verbose=False), progress=self.progress)

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
        #    'source'                : u'Bacteria_74',
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

        # compatibility sanity checks 1/2: make sure domains between domain predictor and the contigs database match
        self.domains_missing_in_SCG_domain_predictor = [d for d in self.domains if d not in self.SCG_domain_predictor.SCG_domains]
        self.domains_missing_in_SCGs_run_for_contigs = [d for d in self.SCG_domain_predictor.SCG_domains if d not in self.domains]

        if len(self.domains_missing_in_SCG_domain_predictor):
            num_domains_missing = len(self.domains_missing_in_SCG_domain_predictor)
            self.progress.reset()
            self.run.warning("OK. We have a problem. You seem to have single-copy core gene collections for among your HMM hits %s that "
                             "are not included when the anvi'o domain predictor was trained :/ Here is the list of domains that are making "
                             "us upset here: \"%s\". This means either you put a new HMM single-copy core gene collection to the anvi'o HMMs "
                             "directory, or gave it as a parameter, and run `anvi-run-hmms` without updating the classifier anvi'o uses to "
                             "resolve domains for proper completion/redundancy estimates." % \
                                           ('a domain' if num_domains_missing == 1 else '%s domains' % num_domains_missing,
                                            ', '.join(self.domains_missing_in_SCG_domain_predictor)))
            self.initialized_properly = False

        if len(self.domains_missing_in_SCGs_run_for_contigs):
            num_domains_missing = len(self.domains_missing_in_SCGs_run_for_contigs)

            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning("It seems %d of the domains that are known to the classifier anvi'o uses to predict "
                                 "domains for completion estimation are missing from this contigs database. This means, the user didn't run the "
                                 "program `anvi-run-hmms` with default parameters, or removed some essential SCG domains from it later. Here is "
                                 "the list of domains that are making us upset here: \"%s\". Running `anvi-run-hmms` on this your contigs database "
                                 "will likely address this warning." % (num_domains_missing, ', '.join(self.domains_missing_in_SCG_domain_predictor)))

            # since we just established that the user did not run these domains for their contigs database,
            # we will update our self.domains variable to make sure the fucked uppery that will likely take
            # place later is to a convenient minumum:
            self.domains.discard(set(self.domains_missing_in_SCGs_run_for_contigs))

            self.initialized_properly = False

        # compatibility sanity checks 2/2: make sure sources in domain predictor to those in the contigs database
        self.sources_missing_in_SCGs_run_for_contigs = [s for s in self.SCG_domain_predictor.SCG_sources if s not in self.sources]
        self.sources_missing_in_SCG_domain_predictor = [s for s in self.sources if s not in self.SCG_domain_predictor.SCG_sources]
        if len(self.sources_missing_in_SCGs_run_for_contigs):
            num_sources_missing = len(self.sources_missing_in_SCGs_run_for_contigs)

            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning("All the SCG domains necessary to run the predictor covered in the contigs database, however, %s that are used "
                                 "during the training of the domain predictor does not seem to occur in it :/ Here is the list of HMM sources "
                                 "that are making us upset here: \"%s\". This most likely means that either a new version of anvi'o are used with "
                                 "an older set of single-copy core gene sources, or someone is exploring new single-copy core gene sources to see "
                                 "how they behave. That's all good and very exciting, but unfortunately anvi'o will not be able to predict domains "
                                 "due to this incompatibility here. Running `anvi-run-hmms` on this contigs database would've solved this problem "
                                 "but it is not an absolute necessity as anvi'o will continue running by not utilizing domain-specific HMMs for "
                                 "completion/redundancy estimates, and report all the results all at once without prioritizing a single domain." % \
                                               ('an HMM source' if num_sources_missing == 1 else '%s HMM sources' % num_sources_missing,
                                                ', '.join(self.sources_missing_in_SCGs_run_for_contigs)))
            self.initialized_properly = False


        if source_requested:
            if source_requested not in self.sources:
                raise ConfigError('Requested source "%s" is not one of the single-copy gene sources found in the database.' % source_requested)

            # filter out sources that are not requested
            self.sources = [source_requested]
            self.genes_in_db = {source_requested: self.genes_in_db[source_requested]}
            self.hmm_hits_splits_table = utils.get_filtered_dict(self.hmm_hits_splits_table, 'source', set([source_requested]))

        # these will be very useful later. trust me.
        self.unique_gene_id_to_gene_name = {}
        self.splits_unique_gene_id_occurs = {}
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


    def get_best_matching_domain(self, hmm_hits, observed_genes_per_domain, bin_name='UNKNOWN'):
        """Returns the best matching and other domain info using a random forest.

           The input dict `hmm_hits` has no role in predicting the best matching domain,
           but useful to print out some helpful messages when necessary.

           It returns a tuple for best matching domain and how confident is the match
           according to the random forest classifier.
        """

        # learn domains anvi'o hmm hits know about .. NOTE: the key problem here is that due to this line, anvi'o currently allows
        # a single SCG collection per domain. This can be changed, and we will think about that when it is a necessity. No need
        # to be rocket scientists before the need arises.
        domains_in_hmm_hits = sorted([d for d in hmm_hits if list(hmm_hits[d].values())[0]['num_genes_in_model_with_hits']])
        sources_in_hmm_hits = sorted(list(set([list(hmm_hits[s].keys())[0] for s in hmm_hits])))

        # if this class is initialized improperly, it means there are SCG domains in the contigs database anvi'o does not
        # recognize. but for a given bin, all HMM hits may be coming only from domains anvi'o recognizes. in those cases
        # we can predict the domain nicely, and move on with our lives. if hmm hits for a given bin includes hits from the
        # mysterious hmm collection the user defined, then there is not much we can do. here we will test the presence of
        # any HMM hits with a domain we don't recognize, and act accordingly.
        if not self.initialized_properly:
            hits_contain_a_domain_missing_from_SCG_domain_predictor = len(set(domains_in_hmm_hits).intersection(set(self.domains_missing_in_SCG_domain_predictor))) > 0
            hits_contain_a_source_missing_from_SCG_domain_predictor = len(set(sources_in_hmm_hits).intersection(set(self.sources_missing_in_SCG_domain_predictor))) > 0

            if hits_contain_a_domain_missing_from_SCG_domain_predictor or hits_contain_a_source_missing_from_SCG_domain_predictor:
                info_text = "NO DOMAIN ESTIMATION BECAUSE THERE IS WEIRD STUFF GOING ON. Anvi'o is having hard time determining the domain for\
                             this particular genomic bin because it includes HMM hits coming from single-copy core gene collection anvi'o did not\
                             know about when the domain predictor was trained :/ This is not a very big deal as anvi'o will continue showing you\
                             the completion/redundancy estimates for every SCG collection you have in this contigs database, but it will predict\
                             the proper domain for you."
                return ('', {}, {}, info_text)

        # learn domain predictions from anvi'o random forest
        domain_probabilities, actual_domains, control_domains = self.SCG_domain_predictor.predict_from_observed_genes_per_domain(observed_genes_per_domain)

        if anvio.DEBUG:
            self.run.warning(None, header="DOMAIN ESTIMTES FOR '%s'" % bin_name, lc='green')
            for domain in control_domains:
                self.run.info_single("Probability %s %.2f" % (domain.upper(), domain_probabilities[domain]), mc='cyan')
            for domain in actual_domains:
                source = self.SCG_domain_predictor.SCG_domain_to_source[domain]
                if domain in domains_in_hmm_hits:
                    self.run.info_single("Domain '%8s' (probability: %.2f) C/R: %.2f/%.2f" % (domain,
                                                                                              domain_probabilities[domain],
                                                                                              hmm_hits[domain][source]['percent_completion'],
                                                                                              hmm_hits[domain][source]['percent_redundancy']), mc='green')
                else:
                    self.run.info_single("Domain '%8s' (probabiity: %.2f) (HMMs were not run for this / had 0 hits)" % (domain, domain_probabilities[domain]), mc='red')

        # figure out the best matching domain and its confidence by simply sorting
        # actual domains first.
        best_matching_domain, domain_matching_confidence = sorted([d for d in domain_probabilities.items() if d[0] in actual_domains], key = lambda x: x[1], reverse=True)[0]

        # if the confidence is less than 0.2, then we are in the world of noise.
        # pick the control domain that matches best:
        if domain_matching_confidence < 0.20:
            best_matching_domain, domain_matching_confidence = sorted([d for d in domain_probabilities.items() if d[0] in control_domains], key = lambda x: x[1], reverse=True)[0]

        # figure out the completion and redundancy given the best matching domain
        # for further filtering down below.
        if best_matching_domain in domains_in_hmm_hits:
            source = self.SCG_domain_predictor.SCG_domain_to_source[best_matching_domain]
            best_mathcing_domain_completion, best_matching_domain_redundancy = hmm_hits[best_matching_domain][source]['percent_completion'], \
                                                                               hmm_hits[best_matching_domain][source]['percent_redundancy']
        else:
            best_mathcing_domain_completion, best_matching_domain_redundancy = None, None

        # figure shit out
        info_text = ''
        max_confidence = max(domain_probabilities.values())
        if best_matching_domain in control_domains:
            if best_matching_domain == "mixed":
                info_text = "ANVI'O IS GETTING MIXED DOMAIN SIGNAL. This often means that the set of contigs here probably are\
                             originating from a populations that belong to different domains of life. This happens when your\
                             genomic bin includes tremendous amount of contamination that spans through archaea to bacteria to\
                             who knows what. That's OK, but anvi'o is unable to offer a completion or redundancy estimate for this\
                             selection."
            elif best_matching_domain == "blank":
                info_text = "NO DOMAIN ESTIMATION BECAUSE THERE IS NO SIGNAL. So anvi'o is having hard time determining any domain\
                             for this set of contigs either because the number of contigs are very little, or there are no\
                             SCGs among the selected ones. This may happen if you are working with genomes that are\
                             extremely low completion, or alternatively coming from parts of life that are very \
                             understudied (such as viruses or plasmids, etc). If these do not apply to you, and you are\
                             sure your set of contigs represents a proper genome, then either anvi'o made a mistake, or you\
                             stumbled upon a graet story."
            else:
                info_text = "ANVI'O IS CONFUSED. Your best predicted domain for this set of contigs seem to be a 'control domain'\
                             yet the code does not recognize that. So this is a question for the programmers :/"
        else:
            if max_confidence < 0.5 and domain_probabilities['mixed'] >= 0.25:
                info_text = "CRAP DOMAIN EST BECAUSE STUFF IS MIXED. Please note that anvi'o determined '%s' \
                             as the best matching domain for your contigs BUT actually the probability of these \
                             contigs to be coming from mixed domains is crazy high (%.2f). Which means, neither \
                             this domain prediction, nor the completion and redundancy estimates should mean much,\
                             and you should take a look at the entire list of C/R estimates from all domain SCGs :/\
                             The good news is that more refined the input contigs (such as you make more and more\
                             precise selections or provdide more and more refined genomes), this situation will\
                             likely correct itself." % (best_matching_domain, domain_probabilities['mixed'])
            elif max_confidence < 0.5 and domain_probabilities['mixed'] < 0.25:
                info_text = "CRAP DOMAIN EST BECAUSE WHO KNOWS WHY. Please note that anvi'o determined '%s' as the\
                             best matching domain for your contigs to predict the completion and redundancy\
                             estimates through sigle-copy core genes, however, since the confidence is as low\
                             as %.2f, you should take this estimate with a grain of salt. This low confidence\
                             is most likely due to a very small number of contigs to offer reliable estimates\
                             of domain." % (best_matching_domain, domain_matching_confidence)
            else:
                if best_matching_domain_redundancy is None:
                    if best_matching_domain in hmm_data.scg_domain_to_source:
                        info_text = "HOUSTON, WE HAVE A PROBLEM. The very high confidence of domain prediction indicates that this set of contigs\
                                     contain enough signal to classify it as %(domain)s. HOWEVER, your contigs database says that there are no hits for\
                                     HMMs described in the collection that serves the domain %(domain)s. This is only possible if you haven't run\
                                     `anvi-run-hmms` for some domains. This problem should go away if you were to run that program on your contigs\
                                     database with the parameter `-I %(collection)s`. If you don't do anuthing it's OK, too. The only problem is that you will\
                                     not get any completion/redundancy estiamtes for domain %(domain)s." % {'domain': best_matching_domain,
                                                                                                            'collection': hmm_data.scg_domain_to_source[best_matching_domain]}
                    else:
                        info_text = "HOUSTON, WE HAVE A VERY BIG PROBLEM. The very high confidence of domain prediction indicates that this set of contigs\
                                     contain enough signal to classify it as %(domain)s. HOWEVER, your contigs database says that there are no hits for\
                                     HMMs described in the collection that serves the domain %(domain)s, which is fine, BUT THEN your anvi'o installatio does not\
                                     seem to have a collection that can be used to estimate the completion of this domain. This is all very very confusing\
                                     and if you let a developer know, they will be happy to investigate how did you end up here :(" % {'domain': best_matching_domain}
                elif best_matching_domain_redundancy < 10:
                    info_text = "GREAT DOMAIN EST & YOU'RE GOLDEN. The very high confidence of domain prediction indicates that this set of contigs\
                                 are almost certainly coming from a population that belongs to %s. IN ADDITION, the low redundancy of SCGs\
                                 do not predict any serious contamination. But please remember that this information does not mean there is NO\
                                 contamination in your genome bin. If you want to take a more careful look, you can try `anvi-refine`." \
                                            % (best_matching_domain)
                elif best_matching_domain_redundancy >= 10 and best_matching_domain_redundancy <= 100:
                    info_text = "GREAT DOMAIN CONFIDENCE (YAY) BUT SOME SERIOUS REDUNDANCY (BOO). The very high confidence of domain prediction indicates that this\
                                 set of contigs are almost certainly coming from a population that belongs to domain %s. BUT you almost certainly are\
                                 looking at a composite genome. You should consider refining this particular collection of contigs to lower the\
                                 redundancy." % (best_matching_domain)
                elif best_matching_domain_redundancy > 100:
                    info_text = "GREAT DOMAIN CONFIDENCE (YAY) BUT A TON OF REDUNDANCY (NOPE). The very high confidence of random forest indicates\
                                 that the very large fraction of this set of contigs are almost certainly coming from populations that belong to the\
                                 domain %s. HOWEVER, the extremely high amount of redundancy of SCGs in domain %s suggests that you either are working\
                                 with a set of contigs that are extremely composite, or you are looking basically an entire metagenome or something." \
                                            % (best_matching_domain, best_matching_domain)
                else:
                    info_text = "GREAT DOMAIN CONFIDENCE BUT ANVI'O MADE A CONFUSE. Your redundancy estimates are weird and anvi'o needs an adult :("

        if anvio.DEBUG:
            self.run.warning(info_text)

        return (best_matching_domain, domain_probabilities, control_domains, remove_spaces(info_text))


    def get_info_for_splits(self, split_names, min_e_value=1e-5, bin_name='UNKNOWN'):
        """This function takes a bunch of split names, and returns three things:

            - Average percent completion for best matching domain
            - Average redundancy for best matching domain
            - Best matching domain for this collection of splits,
            - And a comprehensive results dictionary that explains each HMM source in each domain,

        For your convenience, you can call this function this way:

            p_completion, p_redundancy, domain, domain_probabilities, info_text, hmm_hits_dict = get_info_for_splits(s)
            domain_confidence = domain_probabilities[domain] if domain else 0.0
        """

        hmm_hits_splits_table = utils.get_filtered_dict(self.hmm_hits_splits_table, 'split', split_names)

        # FIXME: the design here is turning into a bad case of spaghetti code. we should reimplement
        # the SCG / completion stuff around the random forest domain predictor. the previous code is
        # quite inefficient, and the late addition of random forest domain predictor is making things
        # even less clear. The following dictionary is to predict the domain:
        observed_genes_per_domain = {}
        for domain in self.domains:
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
        scg_hmm_hits = {}
        for domain in self.domains:
            scg_hmm_hits[domain] = {}

        for source in self.sources:
            domain = self.source_to_domain[source]
            scg_hmm_hits[domain][source] = {'domain': domain, 'source': source}

            genes_count = Counter([v['gene_name'] for v in list(info_dict[source].values())])

            # report num genes in the model and the num of those with hits (note that htis doesn't
            # care whether those hits are contributing to redundance or not --instad here we are
            # intrested only in the 'coverage' of the model)
            scg_hmm_hits[domain][source]['num_genes_in_model'] = len(self.genes_in_db[source])
            scg_hmm_hits[domain][source]['num_genes_in_model_with_hits']= len(genes_count)
            scg_hmm_hits[domain][source]['model_coverage']= len(genes_count) / len(self.genes_in_db[source])

            scg_hmm_hits[domain][source]['percent_completion'] = len(genes_count) * 100.0 / len(self.genes_in_db[source])

            # report redundancy:
            genes_that_occur_multiple_times = [g for g in genes_count if genes_count[g] > 1]
            scg_hmm_hits[domain][source]['percent_redundancy'] = sum([genes_count[g] - 1 for g in genes_that_occur_multiple_times]) * 100.0 / len(self.genes_in_db[source])

            # identify splits that contribute the same single_copy_gene
            redundants = {}
            for gene_name in genes_that_occur_multiple_times:
                redundants[gene_name] = [self.splits_unique_gene_id_occurs[unique_gene_id] for unique_gene_id in gene_name_to_unique_id[source][gene_name]]
            scg_hmm_hits[domain][source]['redundants'] = redundants

        if not len(scg_hmm_hits):
            return (None, None, None, None, "ANVI'O FOUND NO SCG HMM HITS :/", scg_hmm_hits)

        best_matching_domain, domain_probabilities, control_domains, info_text = self.get_best_matching_domain(scg_hmm_hits, observed_genes_per_domain, bin_name)

        if best_matching_domain and best_matching_domain not in control_domains:
            if best_matching_domain not in scg_hmm_hits:
                self.progress.reset()
                self.run.warning("Just so you know: Your process branched into a part of the anvi'o code that run into a weird situation. "
                                 "and wishes to tell you something. This may come accross confusing, and we apologize for that. The thing is, "
                                 "anvi'o is trying to estimate the completion and redundancy of a set of contigs here. Maybe it is doing it "
                                 "because you are running the interactive interface and clicked on something, or you are summarizing a "
                                 "collection, or doing something we never envisioned you would do. The bottom line is this: anvi'o predicts "
                                 "that this set of contigs belong to the domain %s. However, it seems the single-copy core genes for that "
                                 "domain were not run on this contigs database. Hence, you will not get any completion and redundancy "
                                 "estimates for this one :( That is fine and things will continue to run smoothly, but we thought you should "
                                 "know. Because knowledge is power .. even when you're not sure what it means." % best_matching_domain)
                percent_completion, percent_redundancy = 0.0, 0.0
            else:
                source = self.SCG_domain_predictor.SCG_domain_to_source[best_matching_domain]
                percent_completion = scg_hmm_hits[best_matching_domain][source]['percent_completion']
                percent_redundancy = scg_hmm_hits[best_matching_domain][source]['percent_redundancy']
        else:
            percent_completion, percent_redundancy = 0.0, 0.0

        return (percent_completion, percent_redundancy, best_matching_domain, domain_probabilities, info_text, scg_hmm_hits)
