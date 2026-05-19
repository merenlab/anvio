# pylint: disable=line-too-long

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


run = terminal.Run()
progress = terminal.Progress()


VALID_CLASSES = {0, 1, 2, 3, 4, 5}
CLASS_NAMES = {0: 'non-eukaryotic', 1: 'eukaryotic', 2: 'virus', 3: 'plasmid', 4: 'organelle', 5: 'unclassified'}


class TablesForContigClassification:
    def __init__(self, contigs_db_path, run=run, progress=progress):
        self.contigs_db_path = contigs_db_path
        self.run = run
        self.progress = progress

        utils.is_contigs_db(self.contigs_db_path)


    def create(self, data_dict, just_do_it=False):
        """Store contig classification entries in the contigs database.

        Parameters
        ==========
        data_dict : dict
            A dictionary as returned by utils.get_TAB_delimited_file_as_dictionary with
            indexing_field=-1. Keys are integers; values are dicts with keys matching
            t.contig_classification_table_structure.
        just_do_it : bool, False
            If True, existing rows for any conflicting source will be deleted before
            inserting. If False, a ConfigError is raised when a source already exists.
        """

        if not data_dict:
            raise ConfigError("TablesForContigClassification.create() received an empty dictionary.")

        entries = [tuple(row[col] for col in t.contig_classification_table_structure)
                   for row in data_dict.values()]

        database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))

        known_contigs = set(database.get_single_column_from_table(t.contigs_info_table_name, 'contig'))

        unknown_contigs = set(e[0] for e in entries) - known_contigs
        if unknown_contigs:
            database.disconnect()
            raise ConfigError(f"{len(unknown_contigs)} contig name(s) in your input do not appear in the contigs "
                              f"database. Every contig name in the input file must match a contig in the database. "
                              f"Here is one example that was not found: '{next(iter(unknown_contigs))}'. Please "
                              f"make sure you are using the correct contigs database.")

        sources_to_import = set(e[2] for e in entries)

        existing_sources_in_db = set(database.get_single_column_from_table(
            t.contig_classification_table_name, 'source', unique=True))
        conflicting_sources = sources_to_import & existing_sources_in_db

        if conflicting_sources:
            if not just_do_it:
                database.disconnect()
                raise ConfigError(f"The following source(s) already have classification data in this contigs "
                                  f"database: {', '.join(sorted(conflicting_sources))}. If you would like to "
                                  f"overwrite the existing data for these sources, re-run this program with "
                                  f"the --just-do-it flag.")
            else:
                for source in conflicting_sources:
                    database._exec(f'DELETE FROM {t.contig_classification_table_name} WHERE source = ?', (source,))
                self.run.warning(f"Existing classification data for the following source(s) has been deleted "
                                 f"and will be replaced: {', '.join(sorted(conflicting_sources))}.")

        database.insert_many(t.contig_classification_table_name, entries)

        existing_sources_val = database.get_meta_value('contig_classification_sources',
                                                        return_none_if_not_in_table=True)
        existing_sources = set(existing_sources_val.split(',')) if existing_sources_val else set()
        all_sources = existing_sources | sources_to_import
        database.set_meta_value('contig_classification_sources', ','.join(sorted(all_sources)))

        database.disconnect()

        self.run.warning(None, header="CONTIG CLASSIFICATION STORED", lc='green')
        for source in sorted(sources_to_import):
            source_entries = [e for e in entries if e[2] == source]
            class_counts = Counter(CLASS_NAMES[e[1]] for e in source_entries)
            self.run.info(f"Source '{source}'", f"{len(source_entries)} contigs: " +
                          ", ".join(f"{v} {k}" for k, v in sorted(class_counts.items())))


    def get(self, source=None):
        """Retrieve contig classification entries from the contigs database.

        Parameters
        ==========
        source : str, optional
            The source (tool name) to retrieve. If None, all sources are returned.

        Returns
        =======
        list of dict
            Each dict has keys matching the contig_classification table columns.
        """

        database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))

        if source:
            rows = database.get_table_as_dict(t.contig_classification_table_name,
                                               where_clause=f'source = "{source}"')
        else:
            rows = database.get_table_as_dict(t.contig_classification_table_name)

        database.disconnect()

        return list(rows.values()) if rows else []
