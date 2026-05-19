"""
Classes to make sense of GlobDB derived functional annotations.
"""

import os
import glob
import shutil

import yaml

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.diamond import Diamond
from anvio.tables.genefunctions import TableForGeneFunctions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize

# internal version to detect incompatible setups
GLOBAA_DATA_VERSION = '1'

# the location of the database
GLOBAA_DOWNLOAD_URL = 'https://fileshare.lisc.univie.ac.at/e01a9b4189ff4767-1779787805/GlobAA.tar.gz'

# DIAMOND output format used for GlobAA searches. The 'score' column (raw alignment
# score, not bitscore) is essential for LASR calculation and must remain the last field.
GLOBAA_DIAMOND_OUTFMT = '6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore score'

# BLOSUM matrix diagonal values used to compute theoretical self-alignment scores.
# Taken from the GlobAA database developer tools (parser.py).
BLOSUM_DIAGONALS = {
    'BLOSUM45': {
        'A': 5, 'R': 7, 'N': 6, 'D': 7, 'C': 12, 'Q': 6, 'E': 6, 'G': 7,
        'H': 10, 'I': 5, 'L': 5, 'K': 5, 'M': 6, 'F': 8, 'P': 9, 'S': 4,
        'T': 5, 'W': 15, 'Y': 8, 'V': 5, 'B': 5, 'J': 4, 'Z': 5, 'X': 0,
    },
    'BLOSUM62': {
        'A': 4, 'R': 5, 'N': 6, 'D': 6, 'C': 9, 'Q': 5, 'E': 5, 'G': 6,
        'H': 8, 'I': 4, 'L': 4, 'K': 5, 'M': 5, 'F': 6, 'P': 7, 'S': 4,
        'T': 5, 'W': 11, 'Y': 7, 'V': 4, 'B': 4, 'J': 3, 'Z': 4, 'X': 0,
    },
}


class Args():
    pass


def calculate_max_score(sequence, matrix='BLOSUM45'):
    """Calculate the theoretical maximum alignment score for an amino acid sequence.

    This is computed by summing the BLOSUM diagonal value for each amino acid,
    which equals the score a sequence would receive when aligned against itself.
    It serves as the denominator when calculating the Local Alognment Score Ratio (LASR).

    Parameters
    ==========
    sequence : str
        An amino acid sequence (case-insensitive).
    matrix : str
        BLOSUM matrix name. Must be one of the keys in BLOSUM_DIAGONALS.

    Returns
    =======
    int
        The summed BLOSUM diagonal score for the sequence.
    """

    if matrix not in BLOSUM_DIAGONALS:
        raise ConfigError(f"calculate_max_score :: Unknown BLOSUM matrix '{matrix}'. "
                          f"Known matrices are: {', '.join(BLOSUM_DIAGONALS.keys())}.")

    diag = BLOSUM_DIAGONALS[matrix]
    return sum(diag.get(aa, 0) for aa in sequence.upper())


def classify_globaa_hit(max_score, hit_score, lasr_cutoff, selfmin, selfmax):
    """Classify a GlobAA DIAMOND hit using gene-family-level cutoffs.

    A hit is classified based on whether the query sequence's theoretical
    self-alignment score (max_score) falls within the expected size range
    for the gene family (selfmin–selfmax), and whether the alignment score
    (hit_score) passes the LASR threshold.

    Parameters
    ==========
    max_score : float
        Theoretical self-alignment score of the query sequence, computed via
        calculate_max_score().
    hit_score : float
        Raw DIAMOND alignment score for the hit (the 'score' column, not 'bitscore').
    lasr_cutoff : float
        Minimum Local Alignment Score Ratio (LASR) threshold from the gene family YAML cutoffs.
    selfmin : int or float
        Minimum expected max_score for full-length members of this gene family.
    selfmax : int or float
        Maximum expected max_score for full-length members of this gene family.

    Returns
    =======
    str
        One of:
        - 'correct_length': hit passes and sequence is within the expected size range.
        - 'too_short': hit passes but sequence is shorter than expected for this family.
        - 'too_long': hit passes but sequence is longer than expected for this family.
        - 'below_cutoff': hit does not meet the LASR threshold and should be discarded.

    Notes
    =====
    - Only 'below_cutoff' hits should be discarded; all others represent valid
      annotations, though 'too_short' and 'too_long' flag unusual sequence lengths.
    - The halfway point between selfmin and selfmax acts as a pivot: hits for sequences
      longer than halfway are held to the LASR threshold computed at halfway rather than
      at max_score, preventing very long sequences from being unfairly penalized.
    - This algorithm is based on a draft implementation from the GlobAA database
      developers and may be revised as the database matures.
    """

    halfway = selfmin + (selfmax - selfmin) / 2
    lasr_line_value = lasr_cutoff * max_score
    halfway_line_value = lasr_cutoff * halfway

    if max_score < selfmin:
        return 'too_short' if hit_score >= lasr_line_value else 'below_cutoff'

    if max_score <= halfway:
        return 'correct_length' if hit_score >= lasr_line_value else 'below_cutoff'

    if max_score <= selfmax:
        return 'correct_length' if hit_score >= halfway_line_value else 'below_cutoff'

    # max_score > selfmax
    return 'too_long' if hit_score >= halfway_line_value else 'below_cutoff'


def validate_yaml_entry(gaa_id, entry):
    """Validate a single GlobAA gene family YAML entry for required fields and sane values.

    Used both during database setup (anvi-setup-globdb-functions) and at the start
    of annotation runs (anvi-run-globdb-functions) to catch corrupted or hand-edited
    data early.

    Parameters
    ==========
    gaa_id : str
        The GAA identifier (e.g., 'GAA00000001').
    entry : dict
        The parsed YAML content for this gene family (the value under the gaa_id key).

    Raises
    ======
    ConfigError
        If any required field is missing, is of the wrong type, or has an invalid value.
    """

    for field in ['gene_family', 'description', 'version', 'cutoffs']:
        if field not in entry:
            raise ConfigError(f"The YAML entry for '{gaa_id}' is missing the required field '{field}'.")

    for field in ['gene_family', 'description']:
        if not isinstance(entry[field], str) or not entry[field].strip():
            raise ConfigError(f"The YAML entry for '{gaa_id}' has an empty or invalid '{field}' value.")

    version = entry['version']
    if not isinstance(version, dict):
        raise ConfigError(f"The YAML entry for '{gaa_id}' has an invalid 'version' field — expected a dict.")

    for field in ['globdb', 'gene_family']:
        if field not in version:
            raise ConfigError(f"The YAML entry for '{gaa_id}' is missing 'version.{field}'.")
        if not isinstance(version[field], (int, float)) or version[field] <= 0:
            raise ConfigError(f"The YAML entry for '{gaa_id}' has an invalid 'version.{field}' "
                              f"value: {version[field]}. It must be a positive number.")

    cutoffs = entry['cutoffs']
    if not isinstance(cutoffs, dict):
        raise ConfigError(f"The YAML entry for '{gaa_id}' has an invalid 'cutoffs' field — expected a dict.")

    for field in ['lasr', 'selfmax', 'selfmin', 'matrix']:
        if field not in cutoffs:
            raise ConfigError(f"The YAML entry for '{gaa_id}' is missing 'cutoffs.{field}'.")

    lasr = cutoffs['lasr']
    selfmax = cutoffs['selfmax']
    selfmin = cutoffs['selfmin']
    matrix = cutoffs['matrix']

    if not isinstance(lasr, (int, float)) or not (0 < lasr <= 1):
        raise ConfigError(f"The YAML entry for '{gaa_id}' has an invalid 'cutoffs.lasr' value: {lasr}. "
                          f"It must be a number strictly between 0 and 1.")

    if not isinstance(selfmax, (int, float)) or selfmax <= 0:
        raise ConfigError(f"The YAML entry for '{gaa_id}' has an invalid 'cutoffs.selfmax' value: {selfmax}. "
                          f"It must be a positive number.")

    if not isinstance(selfmin, (int, float)) or selfmin <= 0:
        raise ConfigError(f"The YAML entry for '{gaa_id}' has an invalid 'cutoffs.selfmin' value: {selfmin}. "
                          f"It must be a positive number.")

    if selfmax <= selfmin:
        raise ConfigError(f"The YAML entry for '{gaa_id}' has 'cutoffs.selfmax' ({selfmax}) <= "
                          f"'cutoffs.selfmin' ({selfmin}). selfmax must be strictly greater than selfmin.")

    if not isinstance(matrix, str) or matrix.strip() not in BLOSUM_DIAGONALS:
        raise ConfigError(f"The YAML entry for '{gaa_id}' has an unknown 'cutoffs.matrix' value: "
                          f"'{matrix}'. Known matrices are: {', '.join(BLOSUM_DIAGONALS.keys())}.")


def parse_globaa_diamond_output(tabular_path):
    """Parse GlobAA DIAMOND output and return the best hit per query sequence.

    Parameters
    ==========
    tabular_path : str
        Path to the DIAMOND tabular output file produced using GLOBAA_DIAMOND_OUTFMT.

    Returns
    =======
    dict
        Keys are query sequence IDs. Values are dicts with:
        target_id, pident (float), qlen (int), score (float), evalue (float).

    Notes
    =====
    - Best hit is the one with the highest raw alignment score ('score' column).
    - Ties are broken by lower evalue.
    - Lines with fewer than 15 fields are silently skipped.
    """

    results = {}

    for line in open(tabular_path):
        fields = line.strip().split('\t')

        if len(fields) < 15:
            continue

        query_id = fields[0]
        target_id = fields[1]
        pident = float(fields[2])
        qlen = int(fields[3])
        evalue = float(fields[12])
        score = float(fields[14])

        if query_id in results:
            existing = results[query_id]
            if score < existing['score']:
                continue
            if score == existing['score'] and evalue >= existing['evalue']:
                continue

        results[query_id] = {
            'target_id': target_id,
            'pident': pident,
            'qlen': qlen,
            'score': score,
            'evalue': evalue,
        }

    return results


class GlobAASetup:
    """Download, validate, and set up the GlobAA functional annotation database like a boss."""

    def __init__(self, args=Args(), globdb_data_dir=None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.reset = A('reset')

        # many ways to find the data directory
        if globdb_data_dir:
            self.globdb_base_dir = globdb_data_dir
            self.data_source = 'The function call.'
        elif A('globdb_data_dir'):
            self.globdb_base_dir = A('globdb_data_dir')
            self.data_source = 'The command line parameter.'
        elif 'ANVIO_GLOBAA_DATA_DIR' in os.environ:
            self.globdb_base_dir = os.environ['ANVIO_GLOBAA_DATA_DIR']
            self.data_source = 'The environmental variable.'
        else:
            self.globdb_base_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/GlobAA')
            self.data_source = "The anvi'o default."

        self.globdb_base_dir = os.path.abspath(os.path.expanduser(self.globdb_base_dir))
        self.globdb_data_dir_version = os.path.join(self.globdb_base_dir, '.VERSION')

        self.run.info('GlobAA data source', self.data_source)
        self.run.info('GlobAA data directory', self.globdb_base_dir)


    def get_db_path(self):
        """Return the DIAMOND database path (without the .dmnd extension)."""
        return os.path.join(self.globdb_base_dir, 'DB_DIAMOND', 'GlobAA')


    def get_yaml_path(self):
        """Return the path to the master YAML file."""
        return os.path.join(self.globdb_base_dir, 'GlobAA.yaml')


    def is_database_exists(self):
        """Return True if the GlobAA data directory looks complete and up to date."""

        if not os.path.exists(self.globdb_base_dir):
            return False

        if not os.path.exists(self.globdb_data_dir_version):
            return False

        if open(self.globdb_data_dir_version).read().strip() != GLOBAA_DATA_VERSION:
            return False

        if not os.path.exists(self.get_yaml_path()):
            return False

        if not os.path.exists(self.get_db_path() + '.dmnd'):
            return False

        return True


    def create(self):
        """Download, validate, and set up the GlobAA database."""

        if not self.reset and self.is_database_exists():
            raise ConfigError(f"The GlobAA data is already set up in '{self.globdb_base_dir}'. "
                              f"If you want to re-download and rebuild everything from scratch, "
                              f"use the `--reset` flag.")

        filesnpaths.check_output_directory(self.globdb_base_dir, ok_if_exists=True)

        if not os.path.exists(self.globdb_base_dir):
            try:
                os.makedirs(self.globdb_base_dir)
            except Exception as e:
                raise ConfigError(f"Anvi'o could not create the GlobAA data directory at "
                                  f"'{self.globdb_base_dir}' :/ This is unexpected (since "
                                  f"anvi'o already checked if it could be generated and there "
                                  f"was nothing concerning), but here is what went wrong: {e}.")
        else:
            filesnpaths.is_output_dir_writable(self.globdb_base_dir)

        if self.reset:
            self.run.warning("The `--reset` flag is set. Anvi'o will wipe the GlobAA data directory "
                             "and re-download everything from scratch.")
            shutil.rmtree(self.globdb_base_dir)
            os.makedirs(self.globdb_base_dir)

        open(self.globdb_data_dir_version, 'w').write(GLOBAA_DATA_VERSION)

        # download and extract the data package
        raw_data_dir = self.download_and_extract()

        # validate all per-family YAML files in the raw data
        self.validate_all_yamls(raw_data_dir)

        # concatenate all FASTA files with GAA-ID-prefixed headers
        fasta_path = os.path.join(self.globdb_base_dir, 'GlobAA.faa')
        self.concatenate_fastas(raw_data_dir, fasta_path)

        # merge all per-family YAMLs into a single master YAML
        yaml_path = self.get_yaml_path()
        self.concatenate_yamls(raw_data_dir, yaml_path)

        # validate the master YAML we just wrote (sanity check on our own output)
        self.validate_master_yaml(yaml_path)

        # build the DIAMOND search database
        self.build_diamond_db(fasta_path)

        # clean up the raw extracted data to save space
        shutil.rmtree(os.path.join(self.globdb_base_dir, 'RAW_DATA'))

        self.run.info_single("GlobAA database setup is complete! You can now use `anvi-run-globdb-functions` "
                             "to annotate your contigs databases 🎉", mc='green', nl_after=1, nl_before=1)


    def download_and_extract(self):
        """Download the GlobAA data tarball and extract it.

        Returns
        =======
        str
            Path to the directory that directly contains the GAA* subdirectories.
        """

        tarball_path = os.path.join(self.globdb_base_dir, 'GlobAA_download.tar.gz')

        self.run.warning(None, header="DOWNLOADING GLOBAA DATA", lc='cyan')
        utils.download_file(GLOBAA_DOWNLOAD_URL, tarball_path, progress=self.progress, run=self.run)

        self.progress.new('Extracting GlobAA data package')
        self.progress.update('...')

        extract_dir = os.path.join(self.globdb_base_dir, 'RAW_DATA')

        if os.path.exists(extract_dir):
            shutil.rmtree(extract_dir)

        os.makedirs(extract_dir)

        utils.tar_extract_file(tarball_path, output_file_path=extract_dir, keep_original=False)

        self.progress.end()

        # here we do a bit of overengineering with the assumption that the tarball
        # may or may not wrap everything in a top-level directory; so, we will now
        # search up to 3 levels deep for the parent of the GAA* subdirectories in
        # case that happens (even though it should NEVER happen):
        gaa_parent = self._find_gaa_parent_dir(extract_dir)
        if not gaa_parent:
            raise ConfigError(f"Anvi'o could not find any GAA* gene family directories in the "
                              f"downloaded data package (searched under '{extract_dir}') :/ The "
                              f"package structure may have changed substantially. Please report "
                              f"this on anvi'o GitHub repository or the Discord channel.")

        return gaa_parent


    def _find_gaa_parent_dir(self, search_dir, depth=0):
        """Recursively search for the directory that directly contains GAA* subdirectories."""

        if depth > 3:
            return None

        gaa_dirs = [d for d in glob.glob(os.path.join(search_dir, 'GAA*')) if os.path.isdir(d)]
        if gaa_dirs:
            return search_dir

        for entry in sorted(os.listdir(search_dir)):
            entry_path = os.path.join(search_dir, entry)
            if os.path.isdir(entry_path):
                result = self._find_gaa_parent_dir(entry_path, depth + 1)
                if result:
                    return result

        return None


    def validate_all_yamls(self, data_dir):
        """Validate every per-family info.yaml in the extracted data package."""
        gaa_dirs = sorted([d for d in glob.glob(os.path.join(data_dir, 'GAA*')) if os.path.isdir(d)])

        if not gaa_dirs:
            raise ConfigError(f"No GAA* gene family directories were found in '{data_dir}'. "
                              f"The downloaded data package appears to be empty or incorrectly structured.")

        self.progress.new('Validating GlobAA YAML files', progress_total_items=len(gaa_dirs))

        seen_ids = set()

        for gaa_dir in gaa_dirs:
            gaa_id = os.path.basename(gaa_dir)
            yaml_path = os.path.join(gaa_dir, 'info.yaml')

            self.progress.increment()
            self.progress.update(f"Validating {gaa_id} ...")

            if not os.path.exists(yaml_path):
                self.progress.end()
                raise ConfigError(f"No 'info.yaml' file was found in '{gaa_dir}'. The data package "
                                  f"appears to be missing metadata for gene family '{gaa_id}'.")

            try:
                with open(yaml_path) as f:
                    entry = yaml.safe_load(f)
            except yaml.YAMLError as e:
                self.progress.end()
                raise ConfigError(f"The YAML file for gene family '{gaa_id}' could not be parsed. "
                                  f"Here is what Python said: {e}.")

            if not isinstance(entry, dict):
                self.progress.end()
                raise ConfigError(f"The YAML file in '{gaa_dir}' does not contain a dict at the "
                                  f"top level. Something is off with its structure.")

            if gaa_id not in entry:
                self.progress.end()
                raise ConfigError(f"The YAML file in '{gaa_dir}' does not contain the expected "
                                  f"top-level key '{gaa_id}'. Found: {', '.join(entry.keys())}.")

            if gaa_id in seen_ids:
                self.progress.end()
                raise ConfigError(f"The gene family identifier '{gaa_id}' appears more than once "
                                  f"in the data package. All identifiers must be unique.")

            seen_ids.add(gaa_id)
            validate_yaml_entry(gaa_id, entry[gaa_id])

        self.progress.end()
        self.run.info('Gene families with valid YAML files', len(gaa_dirs), mc='green')


    def validate_master_yaml(self, yaml_path):
        """Validate the merged master YAML file we produced during setup."""

        self.progress.new('Validating master YAML')
        self.progress.update('...')

        try:
            with open(yaml_path) as f:
                master = yaml.safe_load(f)
        except yaml.YAMLError as e:
            self.progress.end()
            raise ConfigError(f"The master GlobAA YAML at '{yaml_path}' could not be parsed. "
                              f"This is unexpected and may indicate a bug. Here is the error: {e}.")

        if not isinstance(master, dict) or not master:
            self.progress.end()
            raise ConfigError(f"The master GlobAA YAML at '{yaml_path}' is empty or malformed.")

        for gaa_id, entry in master.items():
            validate_yaml_entry(gaa_id, entry)

        self.progress.end()
        self.run.info('Gene families in master YAML', len(master))


    def concatenate_fastas(self, data_dir, output_path):
        """Concatenate all per-family FASTA files into one, prefixing headers with the GAA ID.

        Each sequence header becomes GAA<ID>___<original_header> so that DIAMOND
        hits can be mapped back to their gene family by splitting on '___'.
        """

        gaa_dirs = sorted([d for d in glob.glob(os.path.join(data_dir, 'GAA*')) if os.path.isdir(d)])

        self.progress.new('Concatenating FASTA files', progress_total_items=len(gaa_dirs))

        total_seqs = 0

        with open(output_path, 'w') as out:
            for gaa_dir in gaa_dirs:
                gaa_id = os.path.basename(gaa_dir)
                fasta_path = os.path.join(gaa_dir, 'sequences.faa')

                self.progress.increment()
                self.progress.update(f"Processing {gaa_id} ...")

                if not os.path.exists(fasta_path):
                    self.progress.end()
                    raise ConfigError(f"No 'sequences.faa' was found in '{gaa_dir}'. The data package "
                                      f"appears to be missing sequences for gene family '{gaa_id}'.")

                fasta = fastalib.SequenceSource(fasta_path)
                while next(fasta):
                    out.write(f">{gaa_id}___{fasta.id}\n{fasta.seq}\n")
                    total_seqs += 1
                fasta.close()

        self.progress.end()
        self.run.info('Total sequences in concatenated GlobAA FASTA', pp(total_seqs), mc='green')


    def concatenate_yamls(self, data_dir, output_path):
        """Merge all per-family info.yaml files into a single master YAML dict."""

        gaa_dirs = sorted([d for d in glob.glob(os.path.join(data_dir, 'GAA*'))
                           if os.path.isdir(d)])

        self.progress.new('Merging YAML files')
        self.progress.update('...')

        master = {}
        for gaa_dir in gaa_dirs:
            gaa_id = os.path.basename(gaa_dir)
            with open(os.path.join(gaa_dir, 'info.yaml')) as f:
                entry = yaml.safe_load(f)
            master[gaa_id] = entry[gaa_id]

        with open(output_path, 'w') as f:
            yaml.dump(master, f, default_flow_style=False, allow_unicode=True)

        self.progress.end()
        self.run.info('Master YAML written to', output_path)


    def build_diamond_db(self, fasta_path):
        """Build the DIAMOND search database from the concatenated GlobAA FASTA."""
        if not utils.is_program_exists('diamond', dont_raise=True):
            raise ConfigError("DIAMOND does not appear to be installed on this system. Anvi'o needs "
                              "DIAMOND to build and search the GlobAA database. Please install "
                              "DIAMOND and try again.")

        output_dir = os.path.join(self.globdb_base_dir, 'DB_DIAMOND')
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)

        db_path = self.get_db_path()
        log_file_path = os.path.join(output_dir, 'log.txt')

        self.run.info('DIAMOND log', log_file_path)

        diamond = Diamond(fasta_path, run=self.run, progress=self.progress, num_threads=self.num_threads)
        diamond.run.log_file_path = log_file_path
        diamond.makedb(db_path)


class GlobAAData:
    """Load and provide access to GlobAA gene family metadata."""

    def __init__(self, args=Args(), globdb_data_dir=None, run=run, progress=progress, panic_on_failure_to_init=False):
        self.run = run
        self.progress = progress

        self.setup = GlobAASetup(args, globdb_data_dir=globdb_data_dir, run=self.run, progress=self.progress)

        self.gene_families = None
        self.initialized = False

        if not self.setup.is_database_exists():
            if panic_on_failure_to_init:
                raise ConfigError("It seems the GlobAA database has not been set up on this system yet. "
                                  "You can fix this by running `anvi-setup-globdb-functions`. If you "
                                  "already ran it but used a custom data directory, make sure to pass "
                                  "that same path here via `--globdb-data-dir`.")
            return

        self.init()


    def init(self):
        """Load the master YAML file and validate all entries."""

        self.progress.new('Initializing GlobAA data')
        self.progress.update('Reading master YAML file ...')

        yaml_path = self.setup.get_yaml_path()

        try:
            with open(yaml_path) as f:
                self.gene_families = yaml.safe_load(f)
        except yaml.YAMLError as e:
            self.progress.end()
            raise ConfigError(f"The GlobAA master YAML at '{yaml_path}' could not be parsed. "
                              f"You may need to re-run `anvi-setup-globdb-functions --reset`. "
                              f"Here is what Python said: {e}.")

        if not isinstance(self.gene_families, dict) or not self.gene_families:
            self.progress.end()
            raise ConfigError(f"The GlobAA master YAML at '{yaml_path}' is empty or malformed. "
                              f"Please re-run `anvi-setup-globdb-functions --reset`.")

        self.progress.update('Validating entries ...')
        for gaa_id, entry in self.gene_families.items():
            validate_yaml_entry(gaa_id, entry)

        self.progress.end()
        self.initialized = True


    def get_cutoffs(self, gaa_id):
        """Return the cutoffs dict for a given GAA identifier.

        Parameters
        ==========
        gaa_id : str
            A GAA identifier (e.g., 'GAA00000001').

        Returns
        =======
        dict
            A dict with keys 'lasr', 'selfmax', 'selfmin', 'matrix'.
        """

        if not self.initialized:
            raise ConfigError("GlobAAData.get_cutoffs :: called before data was initialized.")

        if gaa_id not in self.gene_families:
            raise ConfigError(f"GlobAAData.get_cutoffs :: '{gaa_id}' not found in the master YAML. "
                              f"Anvi'o is confused and sorry.")

        return self.gene_families[gaa_id]['cutoffs']


    def get_annotation_text(self, gaa_id):
        """Return a formatted annotation string for the given GAA identifier.

        Parameters
        ==========
        gaa_id : str
            A GAA identifier (e.g., 'GAA00000001').

        Returns
        =======
        str
            A string like 'nxrA: periplasmic nitrite oxidoreductase subunit ...'.
        """

        if not self.initialized:
            raise ConfigError("GlobAAData.get_annotation_text :: called before data was initialized.")

        if gaa_id not in self.gene_families:
            raise ConfigError(f"GlobAAData.get_annotation_text :: '{gaa_id}' not found in master YAML.")

        family = self.gene_families[gaa_id]
        gene_family = family.get('gene_family', 'unknown')
        description = family.get('description', 'no description available')

        return f"{gene_family}: {description}"


class GlobAA:
    """A class to run GlobAA functional annotation against a contigs database or FASTA file."""

    def __init__(self, args=Args(), run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.contigs_db_path = A('contigs_db')
        self.temp_dir_path = A('temporary_dir_path')
        self.output_file_path = A('output_file')

        self.log_file_path = None
        self.aa_sequence_file_input = None

        if not utils.is_program_exists('diamond', dont_raise=True):
            raise ConfigError("DIAMOND does not appear to be installed on this system (wtf). Anvi'o needs "
                              "DIAMOND to search the GlobAA database, and it is a part of the standard "
                              "installation procedures. If you did something funny, please install DIAMOND "
                              "or consider reinstalling your anvi'o environment and try again.")

        if len(args.__dict__):
            self.initialize()


    def initialize(self):
        """Set up GlobAASetup and GlobAAData from args."""

        self.hits = None

        self.globaa_setup = GlobAASetup(self.args, run=self.run, progress=self.progress)
        self.globaa_data = GlobAAData(self.args, run=self.run, progress=self.progress, panic_on_failure_to_init=True)

        if not self.globaa_setup.is_database_exists():
            raise ConfigError("The GlobAA database does not seem to be fully set up on this system. "
                              "But you can fix it by simply running `anvi-setup-globdb-functions`.")


    def process(self, aa_sequences_file_path=None):
        """Annotate genes using GlobAA.

        Parameters
        ==========
        aa_sequences_file_path : str, optional
            Path to a FASTA file of amino acid sequences. If not provided, sequences
            are extracted from the contigs database set in self.contigs_db_path.
        """

        if not aa_sequences_file_path and not self.contigs_db_path:
            raise ConfigError("You must provide either an anvi'o contigs database (`--contigs-db`) "
                              "or a FASTA file with amino acid sequences (`--fasta-file`).")

        if aa_sequences_file_path and self.contigs_db_path:
            raise ConfigError("You provided both an amino acid FASTA file and a contigs database. "
                              "Please choose one input source.")

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)

        if aa_sequences_file_path and not self.output_file_path:
            raise ConfigError("When annotating a FASTA file you must also provide an output file "
                              "path via `--output-file`.")

        if aa_sequences_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)
            filesnpaths.is_file_exists(aa_sequences_file_path)
            self.aa_sequence_file_input = aa_sequences_file_path

        if not aa_sequences_file_path and self.output_file_path:
            self.run.warning("You provided an `--output-file`, but the input is a contigs database. "
                             "Annotations will be stored in the database itself, so the output file "
                             "will not be used.")

        if not self.temp_dir_path:
            self.temp_dir_path = filesnpaths.get_temp_directory_path()
            self.remove_temp_dir_path = True
        else:
            filesnpaths.is_file_exists(self.temp_dir_path)
            filesnpaths.is_output_dir_writable(self.temp_dir_path)
            self.run.warning("Because you set `--temporary-dir-path` by hand, anvi'o will not "
                             "remove its contents when done. Please clean up those files later.")
            self.remove_temp_dir_path = False

        self.run.info('GlobAA data directory', self.globaa_setup.globdb_base_dir)
        self.run.info('Temporary files directory', self.temp_dir_path)
        self.run.info('Number of threads', self.num_threads)

        if not aa_sequences_file_path:
            aa_sequences_file_path = os.path.join(self.temp_dir_path, 'aa_sequences.fa')
            dbops.ContigsSuperclass(self.args, r=terminal.Run(verbose=False)).get_sequences_for_gene_callers_ids(
                output_file_path=aa_sequences_file_path,
                report_aa_sequences=True,
                simple_headers=True)

        # compute theoretical self-alignment scores for all query sequences
        max_scores = self.calculate_max_scores(aa_sequences_file_path)

        # search against the GlobAA DIAMOND database
        search_results_tabular = self.search_with_diamond(aa_sequences_file_path)

        # parse the DIAMOND output (best hit per query)
        self.hits = parse_globaa_diamond_output(search_results_tabular)

        # apply GlobAA cutoffs and store annotations
        self.store_hits(max_scores)

        if self.remove_temp_dir_path:
            shutil.rmtree(self.temp_dir_path)


    def calculate_max_scores(self, aa_sequences_file_path):
        """Compute theoretical self-alignment scores for every sequence in the FASTA.

        Parameters
        ==========
        aa_sequences_file_path : str
            Path to the amino acid FASTA file.

        Returns
        =======
        dict
            Maps sequence ID → integer max score.
        """

        self.progress.new('Calculating maximum alignment scores')
        self.progress.update('...')

        max_scores = {}
        fasta = fastalib.SequenceSource(aa_sequences_file_path)
        while next(fasta):
            max_scores[fasta.id] = calculate_max_score(fasta.seq)
        fasta.close()

        self.progress.end()
        self.run.info('Query sequences processed', pp(len(max_scores)))

        return max_scores


    def search_with_diamond(self, aa_sequences_file_path):
        """Run DIAMOND blastp with the GlobAA-specific parameters.

        Parameters
        ==========
        aa_sequences_file_path : str
            Path to the amino acid FASTA file to search.

        Returns
        =======
        str
            Path to the DIAMOND tabular output file.
        """

        diamond = Diamond(aa_sequences_file_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads)

        diamond.target_fasta = self.globaa_setup.get_db_path()
        diamond.outfmt = GLOBAA_DIAMOND_OUTFMT
        diamond.max_target_seqs = 1
        diamond.evalue = None   # replaced by --min-score in additional_params
        diamond.additional_params_for_blastp = '--matrix blosum45 --masking 0 --comp-based-stats 0 --min-score 50'

        diamond.tabular_output_path = os.path.join(self.temp_dir_path, 'globaa-search-results.txt')

        self.run.log_file_path = self.log_file_path or os.path.join(self.temp_dir_path, 'log.txt')

        diamond.blastp()

        return diamond.tabular_output_path


    def store_hits(self, max_scores):
        """Apply GlobAA cutoffs to DIAMOND hits and store passing annotations.

        Parameters
        ==========
        max_scores : dict
            Maps query sequence ID → integer max score (from calculate_max_scores).
        """

        if not self.hits:
            self.run.warning("GlobAA found no hits for any of your genes. Returning empty-handed, "
                             "but still registering GlobAA as a functional annotation source.")
            if self.contigs_db_path:
                gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
                gene_function_calls_table.add_empty_sources_to_functional_sources({'GlobAA'})
            return

        functions_dict = {}
        self.__entry_id = 0

        def add_entry(gene_callers_id, accession, function, e_value):
            gene_caller_id = int(gene_callers_id) if self.contigs_db_path else str(gene_callers_id)
            functions_dict[self.__entry_id] = {
                'gene_callers_id': gene_caller_id,
                'source': 'GlobAA',
                'accession': accession,
                'function': function,
                'e_value': float(e_value),
            }
            self.__entry_id += 1

        hits_kept = 0
        hits_below_cutoff = 0
        hits_missing_max_score = 0
        missing_families = set()

        for gene_callers_id, hit in self.hits.items():
            target_id = hit['target_id']

            # the target sequence ID is prefixed with the GAA identifier during setup
            # (format: GAA00000001___original_seq_id), so we split on '___' to get it
            # if Daan et al changes this, we are fcked. this may require an additional
            # sanity check at some place else
            gaa_id = target_id.split('___')[0]

            if gaa_id not in self.globaa_data.gene_families:
                missing_families.add(gaa_id)
                continue

            query_max_score = max_scores.get(str(gene_callers_id), 0)
            if not query_max_score:
                hits_missing_max_score += 1
                continue

            cutoffs = self.globaa_data.get_cutoffs(gaa_id)

            classification = classify_globaa_hit(max_score=query_max_score,
                                                 hit_score=hit['score'],
                                                 lasr_cutoff=cutoffs['lasr'],
                                                 selfmin=cutoffs['selfmin'],
                                                 selfmax=cutoffs['selfmax'])

            if classification == 'below_cutoff':
                hits_below_cutoff += 1
                continue

            annotation_text = self.globaa_data.get_annotation_text(gaa_id)
            if classification in ('too_short', 'too_long'):
                annotation_text = f"{annotation_text} ({classification})"

            add_entry(gene_callers_id, gaa_id, annotation_text, hit['evalue'])
            hits_kept += 1

        self.run.info('Hits that passed GlobAA cutoffs', pp(hits_kept), mc='green')
        self.run.info('Hits discarded (below cutoff)', pp(hits_below_cutoff))

        if hits_missing_max_score:
            self.run.warning(f"{hits_missing_max_score} {P('hit', hits_missing_max_score)} could not be "
                             f"evaluated because anvi'o was unable to find the query sequence in the "
                             f"max-score dictionary. This should not happen and likely reflects a "
                             f"bug :/ Please consider reporting it.")

        if missing_families:
            self.run.warning(f"Anvi'o encountered {len(missing_families)} GAA "
                             f"{P('identifier', len(missing_families))} in the DIAMOND output that "
                             f"were not present in the master YAML. This is unexpected and likely "
                             f"means the GlobAA database is corrupted. The affected "
                             f"{P('identifier', len(missing_families))}: "
                             f"{', '.join(sorted(missing_families))}. Consider re-running "
                             f"`anvi-setup-globdb-functions --reset`.")

        if not functions_dict and self.contigs_db_path:
            gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
            gene_function_calls_table.add_empty_sources_to_functional_sources({'GlobAA'})
            return

        if self.contigs_db_path:
            gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
            gene_function_calls_table.create(functions_dict)
        elif self.aa_sequence_file_input:
            utils.store_dict_as_TAB_delimited_file(d=functions_dict,
                                                   output_path=self.output_file_path,
                                                   do_not_write_key_column=True)
