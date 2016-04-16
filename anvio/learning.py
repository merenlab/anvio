# -*- coding: utf-8
"""A simple module with classes for learning operations"""

import cPickle
import numpy as np
import sklearn.ensemble

import anvio
import anvio.utils as utils 
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class RF:
    def __init__(self, model_object_path = "rf.model", r = run, p = progress):
        self.run = r
        self.progress = p
        self.model_object_path = model_object_path

        self.model_initialized = False
        self.model = None
        self.features = None
        self.classes = None


    def train(self, features, data, labels, n_trees = 64):
        self.run.info('RF Train', "%d observations with %d features grouped into %d classes." % (len(data), len(features), len(set(labels))))
        filesnpaths.is_output_file_writable(self.model_object_path)

        self.progress.new('Training')
        self.progress.update('...')
        rf = sklearn.ensemble.RandomForestClassifier(n_estimators=n_trees)
        rf.fit(np.array(data), labels)
        self.progress.end()

        cPickle.dump({'features': features, 'classes': rf.classes_, 'model': rf}, open(self.model_object_path, 'w'))
        self.run.info('Model output', self.model_object_path)

    def predict_from_TAB_delimited_file(self, file_path):
        cols = utils.get_columns_of_TAB_delim_file(file_path)
        return self.predict(utils.get_TAB_delimited_file_as_dictionary(file_path, column_mapping = [str] + [float] * len(cols)))

    def predict(self, data_dict):
        if not self.model_initialized:
            self.initialize_model()

        samples = data_dict.keys()
        self.run.info('Num samples to classify', "%d." % (len(samples)))

        data = []
        for sample in samples:
            datum = []
            for feature in self.features:
                if feature not in data_dict[sample]:
                    raise ConfigError, "RF prediction run into an issue. All features described in the model should be present\
                                        for all observations in the data. However, that is not the case. For instance, feature\
                                        '%s' is in the model, but the entry '%s' in the input data does not have an observation\
                                        for it :/ Not good." % (feature, sample)
                datum.append(data_dict[sample][feature])
            data.append(datum)

        predictions = self.model.predict_proba(data)

        predictions_dict = {}
        for i in range(0, len(samples)):
            sample = samples[i]
            predictions_dict[sample] = {}
            for j in range(0, len(self.classes)):
                _class = self.classes[j]
                predictions_dict[sample][_class] = predictions[i][j]

        return predictions_dict


    def initialize_model(self):
        filesnpaths.is_file_exists(self.model_object_path)
        model_obj = cPickle.load(open(self.model_object_path))

        try:
            self.features = model_obj['features']
            self.classes = model_obj['classes']
            self.model = model_obj['model']
        except:
            raise ConfigError, "RF class does not like the model object it was sent for processing :/ Are you sure you\
                                generated it the way you were supposed to?"

        self.model_initialized = True

        self.run.info('Model', "Initialized with %d features grouped into %d classes." % (len(self.features), len(self.classes)))

