# -*- coding: utf-8
# pylint: disable=line-too-long

import numpy as np

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table
from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TableForGeneLevelCoverages(Table):
    def __init__(self, db_path, parameters, mode, split_names=None, ignore_splits_name_check=False, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.db_path = db_path
        self.parameters = parameters
        self.split_names = split_names
        self.ignore_splits_name_check = ignore_splits_name_check
        self.mode = mode

        if self.mode == 'INSEQ':
            self.table_name = t.gene_level_inseq_stats_table_name
            self.table_structure = t.gene_level_inseq_stats_table_structure
        elif self.mode == 'STANDARD':
            self.table_name = t.gene_level_coverage_stats_table_name
            self.table_structure = t.gene_level_coverage_stats_table_structure
        else:
            raise ConfigError("TableForGeneLevelCoverages class is speaking: you came here with an improper 'mode' "
                              "when you were expected to come with 'INSEQ' or 'STANDARD' modes. Your mode, '%s', "
                              "is not welcome here :(" % (self.mode))

        if not isinstance(parameters, dict):
            raise ConfigError("Parameters must be of type. These are basically the parameters such as "
                              "min_cov_for_detection, outliers_threshold, or zeros_are_outliers, that "
                              "are used to establish gene-level coverage data. Anvi'o stores them in "
                              "the gene database, so it can warn the user what they're about to read "
                              "from the database is not what they actually want to read (because the "
                              "parameters have changed at some point).")


        Table.__init__(self, self.db_path, utils.get_required_version_for_db(db_path), run=self.run, progress=self.progress)

        self.num_entries = 0

        self.collection_name = db.DB(self.db_path, None, ignore_version=True).get_meta_value('collection_name')
        self.bin_name = db.DB(self.db_path, None, ignore_version=True).get_meta_value('bin_name')


    def check_params(self):
        """Make sure params to generate gene-level stats match across the board"""
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        non_matching_parameters = []
        for parameter in self.parameters:
            try:
                parameter_in_db = database.get_meta_value(parameter)
            except:
                database.disconnect()
                raise ConfigError("Bad news of the day: You have a genes database for the collection %s and bin %s. But "
                                  "clearly the parameters you used to generate these gene-level coverage data has little "
                                  "to do with the parameters you are using now. For instance, parameter '%s' was not even "
                                  "stored in the database :/" % \
                                        (self.collection_name, self.bin_name , str(parameter)))

            parameter_user_set = self.parameters[parameter]
            try:
                parameter_in_db, parameter_user_set = float(parameter_in_db), float(parameter_user_set)
            except:
                pass

            if parameter_in_db != parameter_user_set:
                non_matching_parameters.append((parameter, parameter_in_db, parameter_user_set))

        if len(non_matching_parameters):
            e = non_matching_parameters[0]
            database.disconnect()
            raise ConfigError("OK. You have a genes database for the collection %s and bin %s. But %d "
                              "of the parameters you used to generate these gene-level coverage data is not "
                              "matching the matching parameters you are using now. For instance, the database "
                              "has %s for %s, but the same parameter is currently set to %s in your workflow. "
                              "The best solution to this is to remove this database (which is at '%s'), and let "
                              "anvi'o generate another one for you." % \
                                    (self.collection_name, self.bin_name, len(non_matching_parameters), str(e[1]),
                                    e[0], str(e[2]), self.db_path))

        database.disconnect()


    def check_split_names(self):
        """Make sure split names in the genes database match to the expected split names"""

        if not self.ignore_splits_name_check:
            if not self.split_names:
                raise ConfigError("So you want to read gene-level coverage data from this genes database "
                                  "but there is a problem. Here anvi'o is talking to the programmer: there "
                                  "are two modes reading from the genes database. You either create an instance of "
                                  "TableForGeneLevelCoverages with a list of `split_names` so anvi'o can make "
                                  "sure the splits you are looking for are certainly those the database knows "
                                  "about, OR, you set the parameter `ignore_splits_name_check` to True, so anvi'o "
                                  "doesn't care about making sure everything is in order. Well. What is going on "
                                  "here is that someone called the `read` function, but the instance of this "
                                  "class does not know any splits, and the `ignore_splits_name_check` is False.")

            splits_hash = utils.get_hash_for_list(self.split_names)

            db_hash = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path)).get_meta_value('splits_hash')

            if splits_hash != db_hash:
                raise ConfigError("Terrible news of the day: You have a genes database for the collection %s and bin %s. But "
                                  "it seems the splits your collection and bin contained when you generated this database "
                                  "has changed after its creation. Maybe you used `anvi-refine` to add or remove some? Or you "
                                  "imported other data with the same collection and bin name? We can't know. You are the one "
                                  "who is creative. But what we know is that this genes database at '%s' is not one that you "
                                  "can use anymore. The easy solution is this: remove this database, and let anvi'o generate "
                                  "another one for you. Alternatively you can run the same exact command you run right before "
                                  "you get this error. Sometimes that works too." % \
                                        (self.collection_name, self.bin_name, self.db_path))


    def read(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        if not database.get_meta_value('gene_level_coverages_stored'):
            # we don't have any gene-level coverage data stored in this database
            database.disconnect()
            return {}

        self.check_split_names()
        self.check_params()

        self.progress.new("Database bleep bloop")
        self.progress.update("Recovering %s stats from the genes database..." % self.mode)

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        raw_data = database.get_table_as_dict(self.table_name)
        data = {}

        # here we are converting the data as it is stored in the database into something that
        # the rest of anvi'o expects to see how gene-level coverage data should look like
        for entry in raw_data.values():
            gene_callers_id, sample_name = entry['gene_callers_id'], entry['sample_name']

            if gene_callers_id not in data:
                data[gene_callers_id] = {}

            if sample_name not in data[gene_callers_id]:
                data[gene_callers_id][sample_name] = entry

            g, n = data[gene_callers_id][sample_name]['gene_coverage_values_per_nt'], data[gene_callers_id][sample_name]['gene_coverage_values_per_nt']
            data[gene_callers_id][sample_name]['gene_coverage_values_per_nt'] = utils.convert_binary_blob_to_numpy_array(g, 'uint16')

            if n:
                data[gene_callers_id][sample_name]['non_outlier_positions'] = utils.convert_binary_blob_to_numpy_array(n, 'uint16')
            else:
                data[gene_callers_id][sample_name]['non_outlier_positions'] = None

        database.disconnect()
        self.progress.end()

        self.run.warning(None, header="GENE LEVEL COVERAGE STATS RECOVERED (yay)", lc="green")
        self.run.info("Mode", self.mode, mc="red")
        self.run.info("Num genes", len(data))
        self.print_info()

        return data


    def store(self, data):
        self.progress.new("Database bleep bloop")
        self.progress.update("Adding %s stats into the genes database..." % self.mode)

        db_entries = []
        for gene_callers_id in data:
            for sample_name in data[gene_callers_id]:
                entry = data[gene_callers_id][sample_name]

                d = []
                for h in self.table_structure:
                    if h in ['gene_coverage_values_per_nt', 'non_outlier_positions']:
                        d.append(utils.convert_numpy_array_to_binary_blob(np.array(entry[h]), 'uint16'))
                    else:
                        d.append(entry[h])

                db_entries.append(tuple(d), )

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % self.table_name, db_entries)

        for parameter in self.parameters:
            database.remove_meta_key_value_pair(parameter)
            database.set_meta_value(parameter, self.parameters[parameter])

        database.update_meta_value('gene_level_coverages_stored', True)
        database.disconnect()

        self.progress.end()

        self.run.warning(None, header="GENE LEVEL COVERAGE STATS STORED", lc="green")
        self.run.info("Mode", self.mode, mc="red")
        self.run.info("Num genes", len(data))
        self.run.info("Num entries", len(db_entries))
        self.print_info()


    def print_info(self):
        self.run.info("Genes database", self.db_path)
        self.run.info("Collection name", self.collection_name, mc="green")
        self.run.info("Bin name", self.bin_name, mc="green")
