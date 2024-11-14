# -*- coding: utf-8
# pylint: disable=line-too-long
"""A simple module with classes for learning operations"""

import pickle
import numpy as np

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

with terminal.SuppressAllOutput():
    import sklearn.ensemble

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class RF:
    def __init__(self, classifier_object_path="rf.classifier", r=run, p=progress):
        self.run = r
        self.progress = p
        self.classifier_object_path = classifier_object_path

        self.classifier_initialized = False
        self.classifier = None
        self.features = None
        self.classes = None

    def train(self, features, data, labels, n_trees=64):
        self.run.info(
            "RF Train",
            "%d observations with %d features grouped into %d classes."
            % (len(data), len(features), len(set(labels))),
        )
        filesnpaths.is_output_file_writable(self.classifier_object_path)

        self.progress.new("Training")
        self.progress.update("...")
        rf = sklearn.ensemble.RandomForestClassifier(n_estimators=n_trees)
        rf.fit(np.array(data), labels)
        self.progress.end()

        pickle.dump(
            {"features": features, "classes": rf.classes_, "classifier": rf},
            open(self.classifier_object_path, "wb"),
        )
        self.run.info("Classifier output", self.classifier_object_path)

    def predict_from_TAB_delimited_file(self, file_path):
        cols = utils.get_columns_of_TAB_delim_file(file_path)
        return self.predict(
            utils.get_TAB_delimited_file_as_dictionary(
                file_path, column_mapping=[str] + [float] * len(cols)
            )
        )

    def predict(self, data_dict):
        if not self.classifier_initialized:
            self.initialize_classifier()

        samples = list(data_dict.keys())
        self.run.info("Num samples to classify", "%d." % (len(samples)))

        data = []
        for sample in samples:
            datum = []
            for feature in self.features:
                if feature not in data_dict[sample]:
                    raise ConfigError(
                        "RF prediction run into an issue. All features described in the classifier should be present "
                        "for all observations in the data. However, that is not the case. For instance, feature "
                        "'%s' is in the classifier, but the entry '%s' in the input data does not have an observation "
                        "for it :/ Not good." % (feature, sample)
                    )
                datum.append(data_dict[sample][feature])
            data.append(datum)

        predictions = self.classifier.predict_proba(data)

        predictions_dict = {}
        for i in range(0, len(samples)):
            sample = samples[i]
            predictions_dict[sample] = {}
            for j in range(0, len(self.classes)):
                _class = self.classes[j]
                predictions_dict[sample][_class] = predictions[i][j]

        return predictions_dict

    def initialize_classifier(self):
        filesnpaths.is_file_exists(self.classifier_object_path)

        try:
            if anvio.DEBUG:
                classifier_obj = pickle.load(open(self.classifier_object_path, "rb"))
            else:
                with terminal.SuppressAllOutput():
                    classifier_obj = pickle.load(
                        open(self.classifier_object_path, "rb")
                    )
        except UnicodeDecodeError:
            raise ConfigError(
                "Your classifier object is broken. Probably because you generated is using a different verison "
                "of anvi'o. Please create a new one from scratch, and you will probably be golden."
            )

        try:
            self.features = classifier_obj["features"]
            self.classes = classifier_obj["classes"]
            self.classifier = classifier_obj["classifier"]
        except:
            raise ConfigError(
                "RF class does not like the classifier object it was sent for processing :/ Are you sure you "
                "generated it the way you were supposed to?"
            )

        self.classifier_initialized = True

        self.run.info(
            "Random Forest Classifier",
            "Initialized with %d features grouped into %d classes."
            % (len(self.features), len(self.classes)),
        )
