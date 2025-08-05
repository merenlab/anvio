import os

import anvio.db as db
import anvio.tables as t

from anvio.dbinfo import DBInfo as dbi
from anvio.version import versions_for_db_types
from anvio.errors import ConfigError

from anvio.terminal import Run
from anvio.dbinfo import is_blank_profile


def get_db_type_and_variant(db_path, dont_raise=False):
    database = dbi(db_path, dont_raise=dont_raise)
    return (database.db_type, database.variant)



def get_db_type(db_path):
    return get_db_type_and_variant(db_path)[0]



def get_db_variant(db_path):
    return get_db_type_and_variant(db_path)[1]



def get_required_version_for_db(db_path):
    db_type = get_db_type(db_path)

    if db_type not in versions_for_db_types:
        raise ConfigError("Anvi'o was trying to get the version of the -alleged- anvi'o database '%s', but it failed "
                           "because it turns out it doesn't know anything about this '%s' type." % (db_path, db_type))

    return versions_for_db_types[db_type]



def get_all_sample_names_from_the_database(db_path):
    """Returns all 'sample' names from a given database. At least it tries."""

    db_type = get_db_type(db_path)
    database = db.DB(db_path, get_required_version_for_db(db_path))

    if db_type == 'profile':
        samples = []
        try:
            samples = [s.strip() for s in database.get_meta_value('samples').split(',')]
        except:
            pass

        return set(samples)

    elif db_type == 'genes':
        return set([str(i) for i in database.get_single_column_from_table(t.gene_level_coverage_stats_table_name, 'sample_name')])

    elif db_type == 'pan':
        internal_genome_names, external_genome_names = [], []
        try:
            internal_genome_names = [g.strip() for g in database.get_meta_value('internal_genome_names').split(',')]
        except:
            pass

        try:
            external_genome_names = [g.strip() for g in database.get_meta_value('external_genome_names').split(',')]
        except:
            pass

        return set([s for s in internal_genome_names + external_genome_names if s])

    else:
        raise ConfigError("`get_all_sample_names_from_the_database` function does not know how to deal "
                           "with %s databases." % db_type)



def get_all_item_names_from_the_database(db_path, run=Run()):
    """Return all split names or gene cluster names in a given database"""

    all_items = set([])

    database = db.DB(db_path, get_required_version_for_db(db_path))
    db_type = database.get_meta_value('db_type')

    if db_type == 'profile':
        if is_blank_profile(db_path):
            run.warning("Someone asked for the split names in a blank profile database. Sadly, anvi'o does not keep track "
                        "of split names in blank profile databases. This function will return an empty set as split names "
                        "to not kill your mojo, but whatever you were trying to do will not work :(")
            return set([])
        else:
            all_items = set(database.get_single_column_from_table('mean_coverage_Q2Q3_splits', 'item'))
    elif db_type == 'pan':
        all_items = set(database.get_single_column_from_table(t.pan_gene_clusters_table_name, 'gene_cluster_id'))
    elif db_type == 'contigs':
        all_items = set(database.get_single_column_from_table(t.splits_info_table_name, 'split'))
    elif db_type == 'genes':
        all_items = set([str(i) for i in database.get_single_column_from_table(t.gene_level_coverage_stats_table_name, 'gene_callers_id')])
    else:
        database.disconnect()
        raise ConfigError("You wanted to get all items in the database %s, but no one here knows about its type. Seriously,\
                            what is '%s' anyway?" % (db_path, db_type))

    if not len(all_items):
        database.disconnect()
        raise ConfigError("utils::get_all_item_names_from_the_database speaking. Something that should never happen happened :/ "
                          "There seems to be nothing in this %s database. Anvi'o is as confused as you are. Please get in touch "
                          "with a developer. They will love this story." % db_path)

    database.disconnect()

    return all_items



def get_genes_database_path_for_bin(profile_db_path, collection_name, bin_name):
    if not collection_name or not bin_name:
        raise ConfigError("Genes database must be associated with a collection name and a bin name :/")

    return os.path.join(os.path.dirname(profile_db_path), 'GENES', '%s-%s.db' % (collection_name, bin_name))

