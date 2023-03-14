# -*- coding: utf-8
# pylint: disable=line-too-long
"""Reference databases of protein properties."""

import os
import re
import time
import tarfile
import pandas as pd
import multiprocessing as mp

from math import gcd
from glob import glob
from shutil import rmtree
from functools import reduce
from fractions import Fraction
from typing import Dict, List, Tuple
from abc import ABC, abstractmethod, abstractproperty

import anvio.proteinorthology.protein as protein

from anvio.errors import ConfigError
from anvio.utils import download_file
from anvio.terminal import Progress, Run
from anvio.filesnpaths import check_output_directory
from anvio import __file__ as ANVIO_PATH, __version__ as VERSION

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


class ProteinReferenceDatabase(ABC):
    """Protein reference database framework."""
    # By default, files for each database are stored in a subdirectory with the name of the database
    # (e.g., 'modelseed', 'kegg') of the following superdirectory.
    default_superdir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/MISC/PROTEIN_DATA')
    db_name: str
    pretty_db_name: str
    # These are the final files stored in the database subdirectory.
    files: List[str]

    @property
    def default_db_dir(self) -> str:
        return os.path.join(self.default_superdir, self.db_name)

    @abstractproperty
    def loaded(self) -> bool:
        raise NotImplementedError

    @abstractmethod
    def download(self, reset: bool = False) -> None:
        """Download database files."""
        raise NotImplementedError

    @abstractmethod
    def load(self) -> None:
        """Load database files into memory."""
        raise NotImplementedError

    def get_missing_files(self) -> List[str]:
        """Find missing files that should have been downloaded to the database directory."""
        missing = []
        for f in self.files:
            path = os.path.join(self.db_dir, f)
            if os.path.isfile(path):
                continue
            missing.append(path)
        return missing

    def raise_missing_files(self, missing: List[str]) -> None:
        """Raise an exception if there are missing database files."""
        if len(missing) == len(self.files):
            raise ConfigError(
                f"No {self.pretty_db_name} reference database files were found in the database "
                f"directory, '{self.db_dir}'. Download the reference database to a default "
                "directory with the command, 'anvi-get-metabolic-model-file --download-references "
                f"{self.pretty_db_name}."
            )
        elif 0 < len(missing) < len(self.files):
            raise ConfigError(
                f"{len(self.files) - len(missing)} of {len(self.files)} reference database files "
                f"were found in the database directory, '{self.db_dir}'. Re-download the reference "
                "database to a default directory with the command, 'anvi-get-metabolic-model-file "
                f"--download-references {self.pretty_db_name}."
            )

    @property
    def subclass_db_names(self) -> List[str]:
        return sorted(db.db_name for db in ProteinReferenceDatabase.__subclasses__())

    @property
    def subclass_pretty_db_names(self) -> List[str]:
        return sorted(db.pretty_db_name for db in ProteinReferenceDatabase.__subclasses__())

    def _set_up_db_dir(self, reset: bool) -> None:
        if os.path.exists(self.db_dir):
            if reset:
                rmtree(self.db_dir)
            else:
                raise ConfigError(
                    f"The database directory, {self.db_dir}, already exists. The 'reset' option can "
                    "be used to remove the database and set it up again."
                )
        os.mkdir(self.db_dir)

    def _check_reference_database_initialization(self) -> None:
        if not self.loaded:
            raise ConfigError(
                f"The input {self.pretty_db_name} database is not initialized. The 'load' method "
                "must be called."
            )

class ModelSEEDDatabase(ProteinReferenceDatabase):
    """
    The ModelSEED biochemistry database is designed for use in metabolic modeling of plants, fungi,
    and microbes.

    The database is set up in a default directory if a directory is not provided.
    """
    db_name = 'modelseed'
    dl_root = 'https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/'
    # These files have the same names as the downloaded files but are changed by setup.
    files = ('compounds.tsv', 'reactions.tsv')
    # Compounds are identified as cytosolic or extracellular in reactions.
    compartment_ids = {0: 'c', 1: 'e'}

    def __init__(self, db_superdir: str = None, run: Run = Run(), progress: Progress = Progress()) -> None:
        if db_superdir:
            check_output_directory(db_superdir, ok_if_exists=True)
            self.db_dir = os.path.join(db_superdir, self.db_name)
        else:
            self.db_dir = self.default_db_dir
        self.run = run
        self.progress = progress
        self.reactions_table: pd.DataFrame = None
        self.compounds_table: pd.DataFrame = None

    def download(self, reset: bool = False) -> None:
        """Download and set up biochemistry tables."""
        self._set_up_db_dir(reset=reset)
        for f in self.files:
            url = os.path.join(self.dl_root, f)
            path = os.path.join(self.db_dir, f)
            download_file(url, path, progress=self.progress)
        self._set_up_reactions_table()
        self._set_up_compounds_table()

    def load(self) -> None:
        """Load the reaction and compound tables as DataFrame attributes."""
        missing = self.get_missing_files()
        self.raise_missing_files(missing)
        self.reactions_table = self._load_reactions()
        self.compounds_table = self._load_compounds()

    @property
    def loaded(self) -> bool:
        if self.reactions_table is None or self.compounds_table is None:
            return False
        return True

    def _set_reaction_lookup_table(self, cross_reference: str) -> None:
        """
        Store a modified version of the reactions table that can be used to look up reactions in the
        cross-referenced database of interest. The new table is stored in an attribute called
        'reaction_lookup_tables', which is created if it does not already exist. The name of the
        cross-referenced database, such as 'KEGG' or 'ec_numbers', must correspond to a column of
        the reactions table.

        Parameters
        ==========
        cross_reference : str
            The cross-referenced database name found in the reactions table header.
        """
        col_names = self.reactions_table.columns.tolist()
        alias_col_names = col_names[col_names.index('ec_numbers') + 1: ]
        if cross_reference in alias_col_names:
            formatted_reactions_table = self._get_reactions_table_per_alias(
                self.reactions_table, cross_reference
            )
        elif cross_reference == 'ec_numbers':
            formatted_reactions_table = self._get_reactions_table_per_ec_number(self.reactions_table)
        else:
            raise ConfigError(
                f"The source, '{cross_reference}', is not recognized as the name of a database "
                "cross-referenced to ModelSEED."
            )
        if hasattr(self, 'reaction_lookup_tables'):
            self.reaction_lookup_tables[cross_reference] = formatted_reactions_table
        else:
            self.reaction_lookup_tables = {cross_reference: formatted_reactions_table}

    def get_reaction(self, reaction_data: Dict):
        """
        Get a reaction object from information in the ModelSEED database.

        Parameters
        ==========
        reaction_data : Dict
            A dictionary representation of a ModelSEED reactions table row.

        Returns
        =======
        anvio.proteinorthology.protein.Reaction
        """
        self._check_reference_database_initialization()
        stoichiometry: str = reaction_data['stoichiometry']
        if pd.isna(stoichiometry):
            # Ignore a reaction if it does not have a chemical equation for some reason.
            return None
        reaction = protein.Reaction()
        # Prefer to ID the reaction by BiGG ID, and if there are multiple BiGG IDs and one
        # corresponds to the ModelSEED reaction ID, by that BiGG ID.
        if pd.isna(reaction_data['select_bigg_id']):
            reaction.id: str = reaction_data['id']
        else:
            reaction.id: str = reaction_data['select_bigg_id']
        reaction.id.replace(' ', '_')
        reaction.name = '' if pd.isna(reaction_data['name']) else reaction_data['name']
        reversibility = reaction_data['reversibility']
        if reversibility == '=' or reversibility == '?':
            # Assume that reactions lacking data ('?') are reversible.
            reaction.reversibility = True
        else:
            reaction.reversibility = False
        decimal_reaction_coefficients = []
        for entry in stoichiometry.split(';'):
            decimal_reaction_coefficients.append(entry.split(':')[0])
        reaction_coefficients = self._to_lcm_denominator(decimal_reaction_coefficients)
        direction = reaction_data['direction']
        if (direction == '>' and reversibility == '<') or (direction == '<' and reversibility == '>'):
            # The way the reaction is written is the opposite of the way the reaction proceeds.
            reaction_coefficients = [-c for c in reaction_coefficients]
        for chemical_entry, int_coefficient in zip(stoichiometry.split(';'), reaction_coefficients):
            split_entry = chemical_entry.split(':')
            reaction.coefficients.append(int_coefficient)
            reaction.compartments.append(self.compartment_ids[int(split_entry[2])])
            chemical = protein.Chemical()
            compound_id = split_entry[1]
            chemical.modelseed_compound_id = compound_id
            chemical_data = self.compounds_table.loc[compound_id].to_dict()
            if chemical_data['BiGG']:
                chemical.select_bigg_id = chemical_data['select_bigg_id']
            if pd.notna(chemical_data['name']):
                chemical.name = chemical_data['name']
            if pd.notna(chemical_data['charge']):
                chemical.charge = chemical_data['charge']
            if pd.notna(chemical_data['formula']):
                chemical.formula = chemical_data['formula']
            if pd.notna(chemical_data['inchikey']):
                chemical.inchi_key = chemical_data['inchikey']
            if pd.notna(chemical_data['smiles']):
                chemical.smiles_string = chemical_data['smiles']
            reaction.chemicals.append(chemical)
        return reaction

    def _to_lcm_denominator(self, floats) -> Tuple[int]:
        def lcm(a, b):
            return a * b // gcd(a, b)
        rationals = [Fraction(f).limit_denominator() for f in floats]
        lcm_denom = reduce(lcm, [r.denominator for r in rationals])
        return tuple(int(r.numerator * lcm_denom / r.denominator) for r in rationals)

    def _set_up_reactions_table(self) -> None:
        """Reorganize the downloaded reaction table, storing in the same location."""
        reactions = self._load_reactions()
        reactions = self._expand_aliases(reactions)
        # Select a BiGG ID for each reaction, inserting the column of select BiGG IDs to the left of
        # the new alias columns.
        cols = reactions.columns.tolist()
        reactions.insert(
            cols.index('source') + 1,
            'select_bigg_id',
            self._select_bigg_ids(reactions)
        )
        path = os.path.join(self.db_dir, 'reactions.tsv')
        reactions.to_csv(path, sep='\t', index=None)
        self.run.info("Set up reactions table", path)

    def _set_up_compounds_table(self) -> None:
        """Change the stored compound table from the one downloaded."""
        path = os.path.join(self.db_dir, 'compounds.tsv')
        compounds = pd.read_csv(path, sep='\t', header=0, low_memory=False)
        compounds = self._expand_aliases(compounds)
        # Select a BiGG ID for each compound, inserting the column of select BiGG IDs to the left of
        # the new alias columns.
        cols = compounds.columns.tolist()
        compounds.insert(
            cols.index('source') + 1,
            'select_bigg_id',
            self._select_bigg_ids(compounds)
        )
        compounds.to_csv(path, sep='\t', index=None)
        self.run.info("Set up compounds table", path)

    def _load_reactions(self) -> pd.DataFrame:
        """Load the reaction table as a DataFrame."""
        path = os.path.join(self.db_dir, 'reactions.tsv')
        reactions = pd.read_csv(path, sep='\t', header=0, low_memory=False)
        return reactions

    def _load_compounds(self) -> pd.DataFrame:
        """Load the compound table as a DataFrame."""
        path = os.path.join(self.db_dir, 'compounds.tsv')
        compounds = pd.read_csv(path, sep='\t', header=0, index_col='id', low_memory=False)
        return compounds

    def _expand_aliases(self, table: pd.DataFrame) -> pd.DataFrame:
        """The downloaded reaction and compound tables each have a column of aliases: IDs
        from different databases and common names. Split these IDs into separate columns."""
        rows = []
        for aliases in table.aliases:
            aliases: str
            row = {}
            if pd.isna(aliases):
                rows.append(row)
                continue
            split_aliases = aliases.split('|')
            for alias in split_aliases:
                sep_index = alias.index(': ')
                alias_key = alias[: sep_index]
                alias_value = alias[sep_index + 2: ].lstrip()
                row[alias_key] = alias_value
            rows.append(row)
        alias_df = pd.DataFrame(rows)
        alias_df.fillna('')
        expanded_df = pd.concat([table.drop('aliases', axis=1), alias_df], axis=1)
        return expanded_df

    def _select_bigg_ids(self, table: pd.DataFrame) -> List[str]:
        """
        Select a single BiGG ID per compound or reaction.

        If there are multiple BiGG IDs, prefer one that matches the compound or reaction's
        abbreviation entry, else return the first BiGG ID.
        """
        select_ids = []
        for bigg_entry, abbreviation in zip(table.BiGG, table.abbreviation):
            bigg_entry: str
            abbreviation: str
            if pd.isna(bigg_entry):
                select_ids.append(None)
                continue
            bigg_ids = [b.strip() for b in bigg_entry.split(';')]
            if len(bigg_ids) == 1:
                select_ids.append(bigg_ids[0])
                continue
            for bigg_id in bigg_ids:
                if not abbreviation:
                    continue
                if bigg_id == abbreviation:
                    select_ids.append(bigg_id)
                    break
            else:
                select_ids.append(bigg_ids[0])
        return select_ids

    def _get_reactions_table_per_alias(
        self,
        reactions: pd.DataFrame,
        alias: str,
        sep: str = '; '
    ) -> pd.DataFrame:
        """
        Modify a ModelSEED reactions DataFrame, dropping reaction rows without IDs for the alias of
        interest and expanding rows with multiple alias IDs so there is a row per ID.

        Inspection of the reactions table indicates, and it is therefore assumed, that IDs for all
        aliases are delimited by '; '. However, this can be manually changed using 'sep'.
        """
        reactions = reactions.dropna(subset=[alias])
        expanded = []
        alias_col = []
        for ids, row in zip(reactions[alias], reactions.drop(alias, axis=1).itertuples(index=False)):
            ids: str
            for id in ids.split(sep):
                alias_col.append(id)
                expanded.append(row)
        reactions = pd.DataFrame(expanded)
        reactions[alias] = alias_col
        return reactions

    def _get_reactions_table_per_ec_number(self, reactions: pd.DataFrame) -> pd.DataFrame:
        """
        Modify a ModelSEED reactions DataFrame, dropping reaction rows without EC number references
        and expanding rows with multiple EC numbers so there is a row per EC number.

        Unlike alias IDs, EC numbers are delimited by '|'.
        """
        reactions = self._get_reactions_table_per_alias(reactions, 'ec_numbers', sep='|')
        return reactions

class KEGGDatabase(ProteinReferenceDatabase):
    db_name = 'kegg'
    dl_root = 'https://rest.kegg.jp/' # See: https://www.kegg.jp/kegg/rest/keggapi.html
    # These files are set up from downloaded files.
    files = ('ko_data.tsv', 'reaction_data.tsv')
    # Download files from the following KEGG databases.
    db_categories = {'ko': 'KO', 'reaction': 'Reaction'}

    def __init__(
        self,
        db_superdir: str = None,
        num_threads: int = 1,
        run: Run = Run(),
        progress: Progress = Progress()
    ) -> None:
        if db_superdir:
            check_output_directory(db_superdir, ok_if_exists=True)
            self.db_dir = os.path.join(db_superdir, self.db_name)
        else:
            self.db_dir = self.default_db_dir
        self.num_threads = num_threads
        self.run = run
        self.progress = progress
        self.ko_data: pd.DataFrame = None
        self.reaction_data: pd.DataFrame = None

    def download(self, reset: bool = False) -> None:
        """Download KEGG files and set up relational tables."""
        self._set_up_db_dir(reset=reset)
        if self.num_threads == 1:
            self.run.warning(
                "Only 1 thread is being used to download from KEGG. Set a higher number of threads "
                "to download faster."
            )
        for db_category in self.db_categories:
            self._download_kegg_db_txt_files(db_category)
            self._make_kegg_db_table(db_category)
            self._archive_kegg_db(db_category)

    def load(self):
        """Load the KO and reaction tables as DataFrame attributes."""
        missing = self.get_missing_files()
        self.raise_missing_files(missing)
        self.ko_data = self._load_ko_data()
        self.reaction_data = self._load_reaction_data()

    @property
    def loaded(self):
        if self.ko_data is None or self.reaction_data is None:
            return False
        return True

    def _load_ko_data(self) -> pd.DataFrame:
        """Load the KO data table and set it up as a DataFrame."""
        path = os.path.join(self.db_dir, 'ko_data.tsv')
        ko_data = pd.read_csv(path, sep='\t', header=0, index_col=0, low_memory=False)
        return ko_data

    def _load_reaction_data(self) -> pd.DataFrame:
        """Load the reaction data table and set it up as a DataFrame."""
        path = os.path.join(self.db_dir, 'reaction_data.tsv')
        reaction_data = pd.read_csv(path, sep='\t', header=0, index_col=0, low_memory=False)
        return reaction_data

    def _download_kegg_db_txt_files(self, db_category: str) -> None:
        """
        Download flat files for all entries in a KEGG database.

        Parameters
        ==========
        db_category : str
            Lowercase name of KEGG database with downloadable flat files for each entry, e.g., 'ko',
            'reaction', 'compound'.
        """
        kegg_ids = self._get_kegg_ids(db_category)
        category_dir = os.path.join(self.db_dir, db_category)
        os.mkdir(category_dir)
        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for kegg_id in kegg_ids:
            url = f'{self.dl_root}get/{kegg_id}'
            path = os.path.join(category_dir, f'{kegg_id}.txt')
            input_queue.put((url, path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()
        self.progress.new(f"Downloading KEGG {self.db_categories[db_category]} entry files")
        num_downloaded = 0
        total = len(kegg_ids)
        while num_downloaded < total:
            output_queue.get()
            num_downloaded += 1
            self.progress.update(f"{num_downloaded} / {total}")
        self.progress.end()
        for worker in workers:
            worker.terminate()

    def _get_kegg_ids(self, db_category: str) -> List[str]:
        """
        Get all KEGG entry IDs from the database of the given category.

        Parameters
        ==========
        db_category : str
            Lowercase name of KEGG database, e.g., 'ko', 'reaction', 'compound'.

        Returns
        =======
        list
            List of entry IDs. The 'KO', 'Reaction', and 'Compound' databases have IDs formatted as
            the first letter of the database followed by five digits, e.g., 'K00001', 'R00010'.
        """
        url = f'{self.dl_root}list/{db_category}'
        path = os.path.join(self.db_dir, f'{db_category}.txt')
        download_file(url, path)
        kegg_ids = []
        f = open(path)
        for line in f:
            line.split()[0]
            kegg_ids.append(line[: 6])
        f.close()
        os.remove(path)
        return kegg_ids

    def _make_kegg_db_table(self, db_category: str) -> None:
        """
        Store a tab-delimited file for a KEGG database (e.g., KO, Reaction, Compound) derived
        from downloaded text files for database entries.

        Parameters
        ==========
        db_category : str
            Lowercase name of KEGG database, e.g., 'ko', 'reaction', 'compound'.
        """
        kegg_db_dir = os.path.join(self.db_dir, db_category)
        select_data = {}
        self.progress.new(f"Processing KEGG {self.db_categories[db_category]} files")
        txt_files = glob(os.path.join(kegg_db_dir, '*'))
        total = len(txt_files)
        for num_processed, path in enumerate(txt_files):
            id = os.path.splitext(os.path.basename(path))[0]
            if db_category == 'ko':
                select_data[id] = self._get_ko_data(path)
            elif db_category == 'reaction':
                select_data[id] = self._get_reaction_data(path)
            self.progress.update(f"{num_processed} / {total}")
        self.progress.end()
        if db_category == 'ko':
            header = ['name', 'reactions', 'ec_numbers']
        elif db_category == 'reaction':
            header = ['orthology']
        columns = {h: [] for h in header}
        for data in select_data.values():
            for h, column in columns.items():
                try:
                    value = data[h]
                except KeyError:
                    value = None
                column.append(value)
        table: pd.DataFrame = pd.DataFrame.from_dict(columns)
        table.index = select_data
        table = table.sort_index()
        table_path = os.path.join(self.db_dir, f'{db_category}_data.tsv')
        table.to_csv(table_path, sep='\t')

    def _get_ko_data(self, path: str) -> Dict:
        """
        Get data from a KO database entry.

        Parameters
        ==========
        path : str
            Flat file for KO entry.

        Returns
        =======
        dict
            Data of interest extracted from the file.
        """
        data = {}
        section = None
        f = open(path)
        for line in f:
            if line[0] == ' ':
                pass
            else:
                section = line.split()[0]
            if section == 'NAME':
                # The name value follows 'NAME' at the beginning of the line.
                data['name'] = line[4: ].lstrip().rstrip()
                # EC numbers associated with the KO are recorded at the end of the name value.
                ec_string = re.search('\[EC:.*\]', line)
                if ec_string:
                    data['ec_numbers'] = ec_string[0][4: -1]
            elif section == 'DBLINKS':
                # There is a row for each linked databaes in this section. There can be a row for
                # KEGG Reaction database entries. The first line of the section starts with
                # 'DBLINKS' and is followed by a value for a linked database. Values from the linked
                # database are separated by ': ' from the name of the database, e.g., 'RN: R00001'.
                split_line = line.split()
                try:
                    rn_index = split_line.index('RN:')
                except ValueError:
                    continue
                data['reactions'] = ' '.join(split_line[rn_index + 1: ])
        f.close()
        return data

    def _get_reaction_data(self, path : str) -> Dict:
        """
        Get data from a Reaction database entry.

        Parameters
        ==========
        path : str
            Flat file for Reaction entry.

        Returns
        =======
        dict
            Data of interest extracted from the file.
        """
        data = {}
        section = None
        f = open(path)
        for line in f:
            if line[0] == ' ':
                pass
            else:
                section = line.split()[0]
            if section == 'ORTHOLOGY':
                # A reaction may or may not be associated with KOs that can be involved in its
                # catalysis. Each KO ID, formatted 'Kxxxxx', where each 'x' is a digit, is on a
                # separate line in the section. The first line starts with 'ORTHOLOGY'.
                if line[: 9] == 'ORTHOLOGY':
                    ko = line[9: ].lstrip()[: 6]
                else:
                    ko = line.lstrip()[: 6]
                try:
                    data['orthology'] += f' {ko}'
                except KeyError:
                    data['orthology'] = ko
        f.close()
        return data

    def _archive_kegg_db(self, db_category: str) -> None:
        """
        Turn the directory of downloaded KEGG database files into a tarball.

        Parameters
        ==========
        db_category : str
            Lowercase name of KEGG database, e.g., 'ko', 'reaction', 'compound'.
        """
        self.progress.new(f"Compressing downloaded KEGG {self.db_categories[db_category]} files")
        self.progress.update("...")
        tar_path = os.path.join(self.db_dir, f'{db_category}.tar.gz')
        db_path = os.path.join(self.db_dir, db_category)
        with tarfile.open(tar_path, mode='w:gz') as tar:
            tar.add(db_path, arcname='.')
        self.progress.end()
        rmtree(db_path)
        self.run.info(f"Downloaded KEGG {self.db_categories[db_category]} files", tar_path)

def _download_worker(
    input_queue: mp.Queue,
    output_queue: mp.Queue,
    max_num_tries: int = 10,
    wait_secs: float = 10.0) -> None:
    """Multiprocessing download worker."""
    while True:
        url, path = input_queue.get()
        num_tries = 0
        while True:
            try:
                download_file(url, path)
                break
            except ConfigError as e:
                num_tries += 1
                if num_tries > max_num_tries:
                    raise e
                time.sleep(wait_secs)
        output_queue.put(True)
