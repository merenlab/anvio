# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import glob
import random
import itertools
from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.learning import RF
from anvio.errors import ConfigError

with terminal.SuppressAllOutput():
    import anvio.data.hmm as hmm_data

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class SCGDomainClassifier(object):
    """The base class for training and predicting"""
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        A = lambda x: (args.__dict__[x] if x in args.__dict__ else None) if args else None

        if self.mode == 'train':
            self.genomes_dir = os.path.abspath(A('genomes_dir'))
            self.classifier_output_path = os.path.abspath(A('output'))

            if A('classifier'):
                raise ConfigError("You should not initialize the domain training class with a input classifier path (`args.classifier`).")

            if not self.genomes_dir:
                raise ConfigError("You must provide a genomes directory. Please read the help menu if you are not sure "
                                  "how the contents of this directory should look like.")

            filesnpaths.is_output_file_writable(self.classifier_output_path)
            filesnpaths.is_file_exists(self.genomes_dir)

        elif self.mode == 'predict':
            if A('output'):
                raise ConfigError("You should not initialize the domain prediction class with an output classifier path (`args.output`).")

            default_classifier_path = 'misc/SCGDOMAINCLASSIFIER.rf'
            self.input_classifier_path = A('classifier') or os.path.join(os.path.dirname(anvio.data.__file__), default_classifier_path)

            if A('classifier'):
                filesnpaths.is_file_exists(self.input_classifier_path)
            else:
                if not filesnpaths.is_file_exists(self.input_classifier_path, dont_raise=True):
                    raise ConfigError("Somehow, this anvi'o installation dose not seem to have a SCG domain classifier. This is one of "
                                      "those anvi'o things that should never happen. If you are an anvi'o user, please feel free to panic :( "
                                      "If you are an anvi'o developer, what you need to do is to follow the instructions in "
                                      "`anvi-script-gen-scg-domain-classifier` with a reasonable set of genomes and store the resulting "
                                      "classifier at the default anvi'o path of /blah/blah/anvio/data/%s." % (default_classifier_path))

            self.rf = RF(self.input_classifier_path, r=self.run, p=self.progress)
            self.rf.initialize_classifier()

        else:
            raise ConfigError("Someone initialized the SCG domain classifier class without an explicit mode :(")

        self.SCG_sources = [d for d in hmm_data.sources if hmm_data.sources[d]['kind'] == 'singlecopy']
        self.SCG_domains = sorted([hmm_data.sources[source]['domain'] for source in self.SCG_sources])
        self.SCG_domain_to_source = dict([(hmm_data.sources[source]['domain'], source) for source in self.SCG_sources])

        if not len(self.SCG_sources):
            raise ConfigError("There is something wrong :( There is not even a single SCG source found. Usually "
                              "anvi'o comes with multiple of them :/")

        if len(self.SCG_sources) == 1:
            raise ConfigError("There is only a single SCG source in your anvi'o installation. It is OK if you are "
                              "being a hacker and playing with things, but there is no logic behind creating a "
                              "classifier with a single class.")

        if len(self.SCG_domains) != len(set(self.SCG_domains)):
            raise ConfigError("Something is wrong. For each domain, there must be a single sinlge-copy core gene "
                              "source.")

        self.data, self.labels, self.features  = [], [], []

        for domain in self.SCG_domains:
            self.features.extend(sorted(hmm_data.sources[self.SCG_domain_to_source[domain]]['genes']))

        self.run.info('SCG domain classifier mode', self.mode)
        self.run.info("SCG domains found", ', '.join(self.SCG_domains))
        self.run.info("Num features", len(self.features))


    def get_SCG_vector_for_contigs_db(self, contigs_db_path):
        SCG_features = []
        for domain in self.SCG_domains:
            where_clause = "source = '%s'" % (self.SCG_domain_to_source[domain])

            d = db.DB(contigs_db_path, None, ignore_version=True)
            hits = Counter(d.get_single_column_from_table(t.hmm_hits_table_name, 'gene_name', unique=True, where_clause=where_clause))
            d.disconnect()

            SCG_features.extend([hits[f] for f in self.features])

        return SCG_features


class Train(SCGDomainClassifier):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.mode = 'train'
        self.contigs_dbs = {}

        SCGDomainClassifier.__init__(self, args, run, progress)

        # initialize contigs databases and paths
        self.init()

        # fill in self.data and self.labels for training
        self.process()


    def init(self):
        """Initializes informaiton about the contigs databases."""

        D = lambda domain: os.path.join(self.genomes_dir, domain)
        domain_dirs = [os.path.basename(os.path.abspath(f.path)) for f in os.scandir(self.genomes_dir) if f.is_dir()]

        self.progress.new('Training')
        self.progress.update("Making sure all domain subdirectories are present")
        missing_domain_dirs = [domain for domain in self.SCG_domains if domain not in domain_dirs]
        if len(missing_domain_dirs):
            raise ConfigError("Genomes directory is missing subdirectories for these domains: '%s'." % ', '.join(missing_domain_dirs))

        unexpected_domain_dirs = [domain for domain in domain_dirs if domain not in self.SCG_domains]
        if len(unexpected_domain_dirs):
            self.progress.reset()
            self.run.warning("THIS IS VERY IMPORTANT! In the directory where you have all the domain directories to "
                             "train a new domain classifier, anvi'o found domain directories that did not match any "
                             "known domains. Here is the list of orphan domains we are talking about here: \"%s.\".\
                              This process will continue to train a classifier, but this is a serious problem as anvi'o\
                              will simply ignore all orphan domains. This is becuase the current single copy-core gene\
                              collections anvi'o knows and cares about do not include those orphan domains. Hence, the\
                              training step will not take these orphan domains into consideration :/ If you don't care\
                              about this, you should feel free to move on. If you want to include those one or more\
                              orphan domains, you should first copy the HMM directory you have for that domain into\
                              the directory (which seems to be at '%s' for your anvi'o instance) so it looks just\
                              like another HMM profile for anvi'o, and then re-run the training." % (', '.join(unexpected_domain_dirs), hmm_data.dir_path))

        self.progress.update("Learning about the number of contigs databases in each domain subdirectory")
        for domain in self.SCG_domains:
            self.contigs_dbs[domain] = glob.glob(os.path.join(D(domain), '*.db'))

            if len(self.contigs_dbs[domain]) == 0:
                self.progress.end()
                raise ConfigError("Each domain subdirectory must include at least one contigs database in it :/")

            if len(self.contigs_dbs[domain]) < 20:
                self.progress.reset()
                self.run.warning("The number of contigs databases found for the domain '%s' is %d. You should consider "
                            "increasing the number of genomes you include for this domain. A robust classifier "
                            "will require similar number of genomes for each domain that capture the diversity "
                            "of the domain they represent. Say, at least 20 genomes per domain is a good start." \
                                    % (domain, len(self.contigs_dbs[domain])))

            self.progress.update("Making sure contigs dbs are contigs dbs")
            for contigs_db_path in self.contigs_dbs[domain]:
                utils.is_contigs_db(contigs_db_path)

        self.progress.end()


    def process(self):
        """Processes contigs databases to learn features, data, and labels."""

        for domain in self.SCG_domains:
            for contigs_db_path in self.contigs_dbs[domain]:
                self.labels.append(domain)
                self.data.append(self.get_SCG_vector_for_contigs_db(contigs_db_path))

        # add noise for mixed. here we take the combination of multiple domains
        for num_domains_to_combine in range (2, len(self.SCG_domains)):
            for combined_domains in list(itertools.combinations(self.SCG_domains, num_domains_to_combine)):
                # do it 5 times for signal. totally arbitrary.
                for _ in range(0, 5):
                    mixed_contigs_dbs = []
                    for domain in combined_domains:
                        mixed_contigs_dbs.append(random.choice(self.contigs_dbs[domain]))

                    vectors = []
                    mixed_vector = [0] * len(self.SCG_domains) * len(self.features)
                    for contigs_db_path in mixed_contigs_dbs:
                        vectors.append(self.get_SCG_vector_for_contigs_db(contigs_db_path))

                    # here we have all the vectors from mixed domains, and we will do a
                    # poor man's logical or down below to generate mixed noise.
                    for v in vectors:
                        for i in range(0, len(mixed_vector)):
                            mixed_vector[i] = mixed_vector[i] or v[i]

                    self.labels.append('mixed')
                    self.data.append(mixed_vector)

        # add noise for no domain. this is to cover the case of no or extremely few SCGs.
        for i in range(0, 10):
            blank_vector = [0] * len(self.SCG_domains) * len(self.features)

            for j in random.sample(range(0, len(blank_vector)), i):
                blank_vector[j] = 1

            self.labels.append('blank')
            self.data.append(blank_vector)


    def train(self):
        rf = RF(self.classifier_output_path, r=self.run, p=self.progress)
        rf.train(self.features, self.data, self.labels)


class Predict(SCGDomainClassifier):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.mode = 'predict'

        SCGDomainClassifier.__init__(self, args, run, progress)

        if self.rf.features != self.features:
            raise ConfigError("There is something terribly wrong here :/ The SCG features anvi'o learned from your contigs database "
                              "and those that are stored in your random forest classifier at '%s' seem to differ from each other. "
                              "This can only happen if you have an older or newer classifier installed on your system compared to "
                              "your contigs database. How you got there is quite curious really. But if you want to solve it quickly "
                              "you can try to re-run `anvi-run-hmms` program on your contigs database. If that also doesn't solve your "
                              "problem, you should get in touch with an anvi'o developer :(")

        missing_classes = [c for c in self.SCG_domains if c not in self.rf.classes]
        if len(missing_classes):
            raise ConfigError("The classifier on your disk is missing some of the mandatory information in it :/ For instance, the "
                              "following SCG domains are not defined in it: '%s'. One way to fix this could be re-training the "
                              "random forest with `anvi-script-gen-scg-domain-classifier` with a reasonable set of genomes. Don't "
                              "forget to store the resulting classifier at the default anvi'o path of /blah/blah/anvio/data/%s." \
                                        % (', '.join(missing_classes), self.input_classifier_path))


    def predict_from_observed_genes_per_domain(self, observed_genes_per_domain):
        features_vector = []
        for domain in self.SCG_domains:
            for gene_name in self.features:
                if domain in observed_genes_per_domain and gene_name in observed_genes_per_domain[domain]:
                    features_vector.append(1)
                else:
                    features_vector.append(0)

        # control domains are those that are not to predict actual domains but
        # the absence of any predictable domain.
        control_domains = ['mixed', 'blank']
        actual_domains = [t for t in self.rf.classes if t not in control_domains]

        domain_probabilities = dict(zip(self.rf.classes, self.rf.classifier.predict_proba([features_vector])[0]))

        return domain_probabilities, actual_domains, control_domains
