# pylint: disable=line-too-long
"""
Setup and run custom functional annotation databases (HMM or DIAMOND-based).

This module provides two classes:

  UserAnnotationDBSetup  --  prepares user-provided HMM or FASTA files as
                             annotation sources and writes them to a structured
                             output directory with a manifest.json index.

  UserAnnotationRunner   --  reads that manifest and runs the appropriate
                             search for each database:
                               * HMM sources  → hmmscan/hmmsearch, results stored
                                 in the gene_functions table with source suffix
                                 '_HMM' (e.g. 'MyModels_HMM')
                               * FASTA sources → diamond blastp, results stored
                                 in the gene_functions table with source suffix
                                 '_DIAMOND' (e.g. 'MyProteins_DIAMOND')
                             The source-name suffix lets users instantly identify
                             the search method used for any annotation.
"""

import os
import json
import gzip
import shutil
import datetime

import anvio
import anvio.db as db
import anvio.dbops as dbops
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.drivers.diamond import Diamond
from anvio.tables.genefunctions import TableForGeneFunctions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['lgallucc']


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize


MANIFEST_FILENAME = 'manifest.json'
HMM_SUBDIR = 'hmm'
DIAMOND_SUBDIR = 'diamond'

# Source name suffixes appended to the user's database name in the gene_functions table.
# These let users know immediately which search method produced a given annotation.
HMM_SOURCE_SUFFIX = '_HMM'
DIAMOND_SOURCE_SUFFIX = '_DIAMOND'


class UserAnnotationDBSetup:
    """Prepares user-provided HMM or FASTA databases for use with `anvi-run-user-annotation`.

    Parameters
    ==========
    args : argparse.Namespace
        Must contain: input_tsv, output_dir, num_threads, reset.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.input_tsv = A('input_tsv')
        self.output_dir = A('output_dir')
        self.num_threads = A('num_threads') or 1
        self.reset = A('reset') or False
        self.list_databases_only = A('list') or False
        self.remove_database = A('remove')

        self.sanity_check()


    def sanity_check(self):
        if not self.output_dir:
            raise ConfigError("You must provide an output directory with --output-dir.")

        # --list and --remove only need the output dir to exist with a manifest
        if self.list_databases_only or self.remove_database:
            if not os.path.isdir(self.output_dir):
                raise ConfigError(f"The directory '{self.output_dir}' does not exist.")
            manifest_path = os.path.join(self.output_dir, MANIFEST_FILENAME)
            if not os.path.exists(manifest_path):
                raise ConfigError(f"No manifest found in '{self.output_dir}'. "
                                  f"Is this a directory created by `anvi-setup-user-annotation-db`?")
            return

        if not self.input_tsv:
            raise ConfigError("You must provide an input TSV file with --input-tsv.")

        filesnpaths.is_file_exists(self.input_tsv)


    def list_databases(self):
        """Print a summary of all databases registered in the manifest."""
        manifest_path = os.path.join(self.output_dir, MANIFEST_FILENAME)
        with open(manifest_path) as f:
            manifest = json.load(f)

        if not manifest:
            self.run.warning(f"The manifest in '{self.output_dir}' is empty.")
            return

        self.run.warning(None, header=f"DATABASES IN {self.output_dir}", lc="cyan")

        for name, entry in manifest.items():
            db_type = entry.get('type', 'unknown').upper()
            added = entry.get('added_on', 'unknown')

            if db_type == 'HMM':
                detail = (f"{entry.get('num_models', '?')} models | "
                          f"cutoff: {entry.get('noise_cutoff_terms', '?')} | "
                          f"source: {name}{HMM_SOURCE_SUFFIX}")
                if entry.get('companion_diamond'):
                    detail += f" + {name}{DIAMOND_SOURCE_SUFFIX} (companion)"
            else:
                detail = (f"{entry.get('num_sequences', '?')} sequences | "
                          f"source: {name}{DIAMOND_SOURCE_SUFFIX}")

            self.run.info(f"[{db_type}] {name}", f"{detail} (added {added})")


    def remove_database_entry(self, name):
        """Remove a database from the manifest and delete its prepared files."""
        manifest_path = os.path.join(self.output_dir, MANIFEST_FILENAME)
        with open(manifest_path) as f:
            manifest = json.load(f)

        if name not in manifest:
            available = ', '.join(f"'{k}'" for k in manifest) or 'none'
            raise ConfigError(f"Database '{name}' not found in the manifest. "
                              f"Available: {available}.")

        entry = manifest.pop(name)
        db_type = entry.get('type')

        deleted = []

        if db_type == 'hmm':
            hmm_dir = entry.get('hmm_dir', '')
            if os.path.isdir(hmm_dir):
                shutil.rmtree(hmm_dir)
                deleted.append(hmm_dir)
            companion = entry.get('companion_diamond', {})
            if companion:
                dmnd_path = companion.get('dmnd_path', '')
                if os.path.exists(dmnd_path):
                    os.remove(dmnd_path)
                    deleted.append(dmnd_path)
        elif db_type == 'diamond':
            dmnd_path = entry.get('dmnd_path', '')
            if os.path.exists(dmnd_path):
                os.remove(dmnd_path)
                deleted.append(dmnd_path)

        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

        self.run.info(f"Removed '{name}' from manifest", '')
        for path in deleted:
            self.run.info('Deleted', path)
        self.run.info('Remaining databases', len(manifest))


    def process(self):
        """Parse the input TSV, detect file types, validate content, and prepare each database."""
        if self.list_databases_only:
            self.list_databases()
            return

        if self.remove_database:
            self.remove_database_entry(self.remove_database)
            return

        entries = self.parse_input_tsv()

        if self.reset:
            filesnpaths.gen_output_directory(self.output_dir, run=self.run, delete_if_exists=True, dont_warn=True)
        else:
            filesnpaths.gen_output_directory(self.output_dir, run=self.run)

        manifest_path = os.path.join(self.output_dir, MANIFEST_FILENAME)
        manifest = {}
        if os.path.exists(manifest_path) and not self.reset:
            with open(manifest_path) as f:
                manifest = json.load(f)

        for name, info in entries.items():
            path = info['path']
            companion_fasta = info['companion_fasta']

            self.run.warning(None, header=f"SETTING UP: {name}", lc="green")

            db_type = self.validate_and_detect_db_type(path)
            self.run.info('Database name', name)
            self.run.info('Source path', path)
            self.run.info('Type', db_type)

            if db_type == 'hmm':
                entry = self.setup_hmm_source(name, path)

                if companion_fasta:
                    companion_type = self._sniff_file_type(companion_fasta)
                    if companion_type != 'diamond':
                        raise ConfigError(f"The companion FASTA provided for '{name}' does not look like "
                                          f"a protein FASTA file (it starts with content expected from an "
                                          f"HMM profile or has unrecognised format). The companion file "
                                          f"must be a protein FASTA.")
                    self.run.warning(None, header=f"SETTING UP COMPANION DIAMOND: {name}", lc="cyan")
                    self.run.info('Companion FASTA', companion_fasta)
                    companion_entry = self.setup_diamond_source(name, companion_fasta)
                    entry['companion_diamond'] = companion_entry
            else:
                if companion_fasta:
                    self.run.warning(f"A companion FASTA was provided for '{name}', but this entry is a "
                                     f"FASTA/DIAMOND database (not an HMM profile). Companion FASTA is only "
                                     f"used with HMM databases. The companion will be ignored.")
                entry = self.setup_diamond_source(name, path)

            manifest[name] = entry

            # Write manifest after every successful entry so a failure mid-run does
            # not orphan already-prepared databases without a registry entry.
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)

        self.run.info('Manifest written to', manifest_path)
        self.run.info('Databases set up', f"{len(entries)} ({', '.join(entries.keys())})")


    def parse_input_tsv(self):
        """Parse the user-provided TSV and return a dict of {name: {path, companion_fasta}}.

        Accepted column layouts:
          2 columns: name<TAB>path
          3 columns: name<TAB>path<TAB>companion_fasta

        The optional third column (companion_fasta) is only meaningful for HMM databases.
        When provided, the setup command also builds a DIAMOND database from that FASTA file,
        and the run command searches both HMM and DIAMOND in a single call. This lets users
        cross-validate HMM hits with DIAMOND hits from the sequences that built the profiles.

        An optional header row and '#'-prefixed comment lines are accepted.
        """
        entries = {}
        header_seen = False

        with open(self.input_tsv) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 2:
                    raise ConfigError(f"Line {line_num} in '{self.input_tsv}' has fewer than 2 "
                                      f"tab-separated fields. Expected format: name<TAB>path "
                                      f"(optional third column: companion_fasta for HMM databases)")

                name = fields[0].strip()
                path = fields[1].strip()
                companion_fasta = fields[2].strip() if len(fields) >= 3 else None

                if not header_seen and name.lower() in ('name', 'db_name', 'database'):
                    header_seen = True
                    continue

                if not name:
                    raise ConfigError(f"Line {line_num} in '{self.input_tsv}' has an empty database name.")

                if not path:
                    raise ConfigError(f"Line {line_num} in '{self.input_tsv}' has an empty path for '{name}'.")

                path = os.path.abspath(path)
                filesnpaths.is_file_exists(path)

                if companion_fasta:
                    companion_fasta = os.path.abspath(companion_fasta)
                    filesnpaths.is_file_exists(companion_fasta)

                if name in entries:
                    raise ConfigError(f"Database name '{name}' appears more than once in '{self.input_tsv}'. "
                                      f"Names must be unique.")

                entries[name] = {'path': path, 'companion_fasta': companion_fasta}

        if not entries:
            raise ConfigError(f"No valid entries found in '{self.input_tsv}'. The file must contain two "
                              f"tab-separated columns (name and path). Lines starting with '#' are comments.")

        return entries


    def _sniff_file_type(self, file_path):
        """Peek at the file content and return 'hmm', 'diamond', or None.

        Parameters
        ==========
        file_path : str
            Path to the file to inspect.
        """
        is_gzipped = file_path.lower().endswith('.gz')
        opener = gzip.open if is_gzipped else open
        with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
            header = f.read(64)
        if header.startswith('HMMER3'):
            return 'hmm'
        if header.lstrip().startswith('>'):
            return 'diamond'
        return None


    def validate_and_detect_db_type(self, file_path):
        """Determine whether a file is an HMM profile or a FASTA file, and verify the content.

        The function first guesses the type from the file extension. It then sniffs the
        actual file content. If the two disagree, a ConfigError is raised so the user
        gets a clear, actionable message instead of a cryptic downstream failure.

        Parameters
        ==========
        file_path : str
            Path to the database file.

        Returns
        =======
        db_type : str
            'hmm' or 'diamond'.
        """
        lower = file_path.lower()
        basename = os.path.basename(file_path)

        if lower.endswith('.hmm') or lower.endswith('.hmm.gz'):
            extension_type = 'hmm'
        elif any(lower.endswith(ext) for ext in ['.faa', '.fasta', '.fa', '.fas', '.fna']):
            extension_type = 'diamond'
        else:
            extension_type = None

        content_type = self._sniff_file_type(file_path)

        # If extension and content both agree (or content is ambiguous), all good
        if extension_type is not None and content_type is not None and extension_type != content_type:
            type_labels = {'hmm': 'an HMM profile (HMMER3 format)',
                           'diamond': 'a protein FASTA file'}
            raise ConfigError(f"The file extension of '{basename}' suggests it is {type_labels[extension_type]}, "
                              f"but the actual file content looks like {type_labels[content_type]}. "
                              f"Please verify that the file is what you think it is, and rename or correct it "
                              f"before running setup again.")

        detected = extension_type or content_type

        if detected is None:
            raise ConfigError(f"Cannot determine the type of '{basename}'. Its extension is not recognised "
                              f"and the content does not match either a valid HMM profile (which should start "
                              f"with 'HMMER3/f') or a protein FASTA file (which should contain sequences "
                              f"starting with '>'). Supported extensions are: .hmm or .hmm.gz for HMM profiles; "
                              f".faa, .fasta, .fa, .fas, or .fna for protein FASTA files.")

        return detected


    def parse_hmm_models(self, hmm_path):
        """Parse an HMMER3 profile file and return a list of per-model metadata dicts.

        Each dict contains: name, acc, has_tc, has_ga, has_nc.

        Parameters
        ==========
        hmm_path : str
            Path to the HMM file (plain text or gzipped).
        """
        models = []
        opener = gzip.open if hmm_path.lower().endswith('.gz') else open

        current = None
        try:
            with opener(hmm_path, 'rt', encoding='utf-8', errors='replace') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.rstrip()
                    if line.startswith('HMMER3/'):
                        current = {'name': None, 'acc': None, 'has_tc': False, 'has_ga': False, 'has_nc': False}
                    elif current is None:
                        continue
                    elif line.startswith('NAME'):
                        parts = line.split()
                        if len(parts) < 2:
                            raise ConfigError(f"Malformed NAME line at line {line_num} in '{hmm_path}': '{line}'")
                        current['name'] = parts[1]
                    elif line.startswith('ACC'):
                        parts = line.split()
                        if len(parts) >= 2:
                            current['acc'] = parts[1]
                    elif line.startswith('TC'):
                        current['has_tc'] = True
                    elif line.startswith('GA'):
                        current['has_ga'] = True
                    elif line.startswith('NC'):
                        current['has_nc'] = True
                    elif line == '//':
                        if current['name'] is None:
                            raise ConfigError(f"Found an HMM model without a NAME field near line {line_num} "
                                              f"in '{hmm_path}'. Every model in a valid HMMER3 profile must "
                                              f"have a NAME line.")
                        models.append(current)
                        current = None
        except (OSError, EOFError) as e:
            raise ConfigError(f"Could not read '{hmm_path}': {e}. Is this a valid HMMER3 profile file?")

        if not models:
            raise ConfigError(f"No HMMER3 models found in '{hmm_path}'. A valid HMMER3 profile starts "
                              f"with a 'HMMER3/f' header and contains at least one model terminated by '//'.")

        # Duplicate model names would silently corrupt genes.txt (dict overwrite downstream).
        names_seen = set()
        duplicates = []
        for m in models:
            if m['name'] in names_seen:
                duplicates.append(m['name'])
            names_seen.add(m['name'])
        if duplicates:
            raise ConfigError(f"The HMM file '{hmm_path}' contains {P('model', len(duplicates))} with "
                              f"duplicate NAME fields: {', '.join(sorted(set(duplicates)))}. Every model in "
                              f"an HMM profile must have a unique NAME. Please deduplicate your profiles "
                              f"before running setup.")

        return models


    def determine_noise_cutoff(self, models, db_name):
        """Select the best HMMER noise cutoff strategy for a set of models.

        Priority: TC (trusted cutoffs) > GA (gathering thresholds) > NC (noise cutoffs) > -E 1e-5.
        A cutoff type is only used when ALL models in the file carry that annotation, because
        HMMER applies a single flag per search run.

        Parameters
        ==========
        models : list
            Output of parse_hmm_models().
        db_name : str
            Database name (for messages only).
        """
        n = len(models)
        num_tc = sum(1 for m in models if m['has_tc'])
        num_ga = sum(1 for m in models if m['has_ga'])
        num_nc = sum(1 for m in models if m['has_nc'])

        if num_tc == n:
            return '--cut_tc'
        if num_ga == n:
            return '--cut_ga'
        if num_nc == n:
            return '--cut_nc'

        partial = []
        if num_tc:
            partial.append(f"TC in {num_tc}/{n} models")
        if num_ga:
            partial.append(f"GA in {num_ga}/{n} models")
        if num_nc:
            partial.append(f"NC in {num_nc}/{n} models")

        if partial:
            self.run.warning(f"Not all models in '{db_name}' share the same cutoff annotation type "
                             f"({', '.join(partial)}). Because HMMER applies a single cutoff flag per "
                             f"run, anvi'o will fall back to -E 1e-5. To use profile-specific cutoffs, "
                             f"ensure every model in the file has the same annotation type (TC, GA, or NC). "
                             f"You can also manually edit noise_cutoff_terms.txt inside the prepared HMM "
                             f"source directory before running the annotation search.")
        else:
            self.run.warning(f"No TC, GA, or NC cutoff annotations found in any model of '{db_name}'. "
                             f"Falling back to -E 1e-5, which may yield more false positives than "
                             f"profile-specific cutoffs. Edit noise_cutoff_terms.txt in the prepared "
                             f"directory to override this before running the annotation search.")

        return '-E 1e-5'


    def setup_hmm_source(self, db_name, hmm_path):
        """Create a standard anvi'o HMM source directory from a raw HMMER3 profile file.

        Parameters
        ==========
        db_name : str
            Name for this database (also becomes the directory name).
        hmm_path : str
            Absolute path to the input HMM file.

        Returns
        =======
        entry : dict
            Manifest entry for this database.
        """
        self.progress.new(f"Setting up HMM source '{db_name}'")
        self.progress.update("Parsing models...")

        models = self.parse_hmm_models(hmm_path)
        noise_cutoff_terms = self.determine_noise_cutoff(models, db_name)

        hmm_dir = os.path.join(self.output_dir, HMM_SUBDIR, db_name)
        filesnpaths.gen_output_directory(hmm_dir, run=self.run, delete_if_exists=True, dont_warn=True)

        self.progress.update("Writing genes.txt...")
        with open(os.path.join(hmm_dir, 'genes.txt'), 'w') as f:
            f.write('gene\taccession\thmmsource\n')
            for m in models:
                acc = m['acc'] if m['acc'] else ''
                f.write(f"{m['name']}\t{acc}\tuser-provided\n")

        self.progress.update("Writing metadata files...")
        with open(os.path.join(hmm_dir, 'kind.txt'), 'w') as f:
            f.write('user_annotation')

        with open(os.path.join(hmm_dir, 'target.txt'), 'w') as f:
            f.write('AA:GENE')

        with open(os.path.join(hmm_dir, 'noise_cutoff_terms.txt'), 'w') as f:
            f.write(noise_cutoff_terms)

        with open(os.path.join(hmm_dir, 'reference.txt'), 'w') as f:
            f.write(f"User-provided HMM database. Original file: {hmm_path}")

        self.progress.update("Copying and compressing HMM file...")
        hmm_gz_path = os.path.join(hmm_dir, 'genes.hmm.gz')
        if hmm_path.lower().endswith('.gz'):
            shutil.copy2(hmm_path, hmm_gz_path)
        else:
            with open(hmm_path, 'rb') as f_in:
                with gzip.open(hmm_gz_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        self.progress.end()

        self.run.info('HMM models found', len(models))
        self.run.info('Noise cutoff strategy', noise_cutoff_terms)
        self.run.info('HMM source directory', hmm_dir)

        return {
            'type': 'hmm',
            'source_path': hmm_path,
            'hmm_dir': hmm_dir,
            'num_models': len(models),
            'noise_cutoff_terms': noise_cutoff_terms,
            'added_on': str(datetime.date.today()),
        }


    def setup_diamond_source(self, db_name, fasta_path):
        """Build a DIAMOND database from a protein FASTA file.

        Parameters
        ==========
        db_name : str
            Name for this database.
        fasta_path : str
            Absolute path to the input protein FASTA file.

        Returns
        =======
        entry : dict
            Manifest entry for this database.
        """
        utils.is_program_exists('diamond')

        diamond_dir = os.path.join(self.output_dir, DIAMOND_SUBDIR)
        filesnpaths.gen_output_directory(diamond_dir, run=self.run)

        dmnd_base = os.path.join(diamond_dir, db_name)

        diamond = Diamond(query_fasta=fasta_path, run=self.run, progress=self.progress, num_threads=self.num_threads)
        diamond.run.log_file_path = os.path.join(self.output_dir, f'{db_name}_diamond_makedb.log')
        diamond.makedb(output_file_path=dmnd_base)

        dmnd_path = dmnd_base + '.dmnd'
        if not os.path.exists(dmnd_path):
            raise ConfigError(f"DIAMOND makedb completed but '{dmnd_path}' was not created. "
                              f"Check the log at '{diamond.run.log_file_path}' for details.")

        num_seqs = 0
        try:
            with open(fasta_path) as f:
                for line in f:
                    if line.startswith('>'):
                        num_seqs += 1
        except Exception:
            pass

        self.run.info('DIAMOND database', dmnd_path)
        self.run.info('Sequences in source FASTA', num_seqs or 'unknown')

        return {
            'type': 'diamond',
            'source_path': fasta_path,
            'dmnd_path': dmnd_path,
            'dmnd_base': dmnd_base,
            'num_sequences': num_seqs,
            'added_on': str(datetime.date.today()),
        }


class UserAnnotationRunner:
    """Runs functional annotation using databases set up by `anvi-setup-user-annotation-db`.

    Results for all database types are stored in the gene_functions table.
    The annotation source name encodes the search method used:
      - HMM profile databases  → source = '{db_name}_HMM'
      - Protein FASTA databases → source = '{db_name}_DIAMOND'

    For HMM sources the gene_functions fields are:
      accession  — HMM model accession (ACC field) or the model name if ACC is absent
      function   — HMM model name (NAME field)
      e_value    — full-sequence E-value from HMMER

    For DIAMOND sources:
      accession  — target sequence ID (best BLASTP hit)
      function   — enriched string with identity %, alignment length, and bit score
      e_value    — BLASTP E-value

    Parameters
    ==========
    args : argparse.Namespace
        Must contain: contigs_db, annotation_db_dir.
        Optional: db_names (comma-separated list), num_threads, just_do_it,
                  hmmer_program, evalue.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.contigs_db_path = A('contigs_db')
        self.annotation_db_dir = A('annotation_db_dir')
        self.database = A('database')
        self.num_threads = A('num_threads') or 1
        self.just_do_it = A('just_do_it') or False
        self.hmmer_program = A('hmmer_program') or 'hmmscan'
        self.evalue = A('evalue')
        self.min_pident = A('min_pident')
        self.diamond_sensitivity = A('diamond_sensitivity')

        self.sanity_check()


    def sanity_check(self):
        if not self.annotation_db_dir:
            raise ConfigError("You must provide the annotation database directory with --annotation-db-dir.")

        if not os.path.isdir(self.annotation_db_dir):
            raise ConfigError(f"The annotation database directory does not exist: '{self.annotation_db_dir}'.")

        manifest_path = os.path.join(self.annotation_db_dir, MANIFEST_FILENAME)
        if not os.path.exists(manifest_path):
            raise ConfigError(f"No manifest file found in '{self.annotation_db_dir}'. "
                              f"Is this a directory created by `anvi-setup-user-annotation-db`? "
                              f"Expected to find a file named '{MANIFEST_FILENAME}'.")

        if not self.contigs_db_path:
            raise ConfigError("You must provide a contigs database with --contigs-db.")

        utils.is_contigs_db(self.contigs_db_path)


    def load_manifest(self):
        manifest_path = os.path.join(self.annotation_db_dir, MANIFEST_FILENAME)
        with open(manifest_path) as f:
            return json.load(f)


    def process(self):
        """Run annotation for all (or selected) databases in the manifest."""
        manifest = self.load_manifest()

        if not manifest:
            raise ConfigError(f"The manifest in '{self.annotation_db_dir}' is empty. "
                              f"Run `anvi-setup-user-annotation-db` first.")

        if self.database and self.database.lower() != 'all':
            requested = [n.strip() for n in self.database.split(',') if n.strip()]
            missing = [n for n in requested if n not in manifest]
            if missing:
                available = ', '.join(f"'{k}'" for k in manifest)
                raise ConfigError(f"The following requested database {P('name', len(missing))} "
                                  f"{P('is', len(missing), alt='are')} not in the manifest: "
                                  f"{', '.join(missing)}. Available: {available}.")
            selected = {n: manifest[n] for n in requested}
        else:
            selected = manifest

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Annotation DB dir', self.annotation_db_dir)
        self.run.info(f"Databases to run ({len(selected)})", ', '.join(selected.keys()))

        for name, entry in selected.items():
            self.run.warning(None, header=f"ANNOTATING WITH: {name}", lc="green")
            db_type = entry.get('type')

            if db_type == 'hmm':
                self.run_hmm_annotation(name, entry)
            elif db_type == 'diamond':
                self.run_diamond_annotation(name, entry)
            else:
                raise ConfigError(f"Unknown database type '{db_type}' for '{name}'. "
                                  f"Expected 'hmm' or 'diamond'. Was this manifest created by "
                                  f"`anvi-setup-user-annotation-db`?")


    def run_hmm_annotation(self, db_name, entry):
        """Search an HMM source and store hits in the gene_functions table.

        Results are written with source name '{db_name}_HMM' so users can
        identify this as an HMM-derived annotation. The HMM model name (NAME field)
        becomes the 'function', the model accession (ACC field) becomes 'accession',
        and the full-sequence E-value is stored in 'e_value'.

        Parameters
        ==========
        db_name : str
            User-defined database name.
        entry : dict
            Manifest entry for this database.
        """
        from anvio.tables.hmmhits import TablesForHMMHits

        hmm_dir = entry['hmm_dir']

        if not os.path.isdir(hmm_dir):
            raise ConfigError(f"The HMM source directory for '{db_name}' no longer exists at '{hmm_dir}'. "
                              f"Please run `anvi-setup-user-annotation-db` again.")

        sources = utils.get_HMM_sources_dictionary([hmm_dir])

        # Rename the source key to carry the _HMM suffix. This suffix will appear
        # in the gene_functions table as the annotation source, letting users distinguish
        # HMM hits from DIAMOND hits at a glance.
        internal_source_name = f'{db_name}{HMM_SOURCE_SUFFIX}'
        sources[internal_source_name] = sources.pop(db_name)

        # `TablesForHMMHits.check_sources` with `add_to_functions_table=True` always
        # raises a ConfigError when the source already exists — even with `--just-do-it`.
        # We handle that case ourselves here so reruns behave consistently with DIAMOND.
        if self.just_do_it:
            contigs_database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))
            existing_sources_raw = contigs_database.get_meta_value('gene_function_sources', return_none_if_not_in_table=True)
            existing_sources = set(existing_sources_raw.split(',') if existing_sources_raw else [])
            if internal_source_name in existing_sources:
                self.run.warning(f"Source '{internal_source_name}' already exists in the gene_functions table. "
                                 f"Removing it first because --just-do-it was passed.")
                gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path,
                                                                   run=terminal.Run(verbose=False),
                                                                   progress=terminal.Progress(verbose=False))
                gene_function_calls_table.drop_functions(contigs_database, sources_to_drop=[internal_source_name])
            contigs_database.disconnect()

        search_tables = TablesForHMMHits(self.contigs_db_path,
                                         num_threads_to_use=self.num_threads,
                                         run=self.run,
                                         progress=self.progress,
                                         just_do_it=self.just_do_it,
                                         hmm_program_to_use=self.hmmer_program,
                                         add_to_functions_table=True)
        search_tables.populate_search_tables(sources)

        # Prefix the function field with '[HMM] ' so the search method is immediately
        # visible in any tabular export alongside the source name suffix.
        contigs_database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))
        contigs_database._exec(
            f"UPDATE {t.gene_function_calls_table_name} "
            f"SET function = '[HMM] ' || function "
            f"WHERE source = ?",
            [internal_source_name]
        )
        contigs_database.disconnect()

        self.run.info('Results stored in', f"gene_functions table (source: '{internal_source_name}')")

        # If a companion DIAMOND database was built at setup time, run it now.
        companion = entry.get('companion_diamond')
        if companion:
            self.run.warning(None, header=f"COMPANION DIAMOND VALIDATION: {db_name}", lc="cyan")
            self.run_diamond_annotation(db_name, companion)


    def run_diamond_annotation(self, db_name, entry):
        """Search a DIAMOND database and store best hits in the gene_functions table.

        Results are written with source name '{db_name}_DIAMOND'. The function field
        carries an enriched string: '{target_id} [pident: {pct}%, aln_len: {len} aa,
        bitscore: {score}]'. The accession field carries the target sequence ID.

        Parameters
        ==========
        db_name : str
            User-defined database name.
        entry : dict
            Manifest entry for this database.
        """
        dmnd_path = entry.get('dmnd_path', '')
        dmnd_base = entry.get('dmnd_base', dmnd_path.removesuffix('.dmnd'))

        if not os.path.exists(dmnd_path):
            raise ConfigError(f"The DIAMOND database for '{db_name}' no longer exists at '{dmnd_path}'. "
                              f"Please run `anvi-setup-user-annotation-db` again.")

        temp_dir = filesnpaths.get_temp_directory_path()

        try:
            aa_sequences_path = os.path.join(temp_dir, 'aa_sequences.faa')

            self.progress.new(f"Exporting AA sequences from contigs DB")
            self.progress.update("...")

            contig_sc = dbops.ContigsSuperclass(self.args,
                                                r=terminal.Run(verbose=False),
                                                p=terminal.Progress(verbose=False))
            contig_sc.get_sequences_for_gene_callers_ids(output_file_path=aa_sequences_path,
                                                          simple_headers=True,
                                                          rna_alphabet=False,
                                                          report_aa_sequences=True)

            self.progress.end()

            if not os.path.exists(aa_sequences_path):
                raise ConfigError(f"Failed to export amino acid sequences from '{self.contigs_db_path}'.")

            tabular_output_path = os.path.join(temp_dir, f'{db_name}_diamond_results.txt')

            diamond = Diamond(query_fasta=aa_sequences_path,
                              run=self.run,
                              progress=self.progress,
                              num_threads=self.num_threads)
            diamond.target_fasta = dmnd_base
            diamond.tabular_output_path = tabular_output_path
            diamond.max_target_seqs = 1
            logs_dir = os.path.join(self.annotation_db_dir, 'logs')
            filesnpaths.gen_output_directory(logs_dir, run=terminal.Run(verbose=False))
            diamond_log_path = os.path.join(logs_dir, f'{db_name}_diamond_blastp.log')
            diamond.run.log_file_path = diamond_log_path

            diamond.evalue = self.evalue if self.evalue is not None else 1e-15

            if self.min_pident is not None:
                diamond.min_pct_id = self.min_pident

            if self.diamond_sensitivity:
                diamond.additional_params_for_blastp = f'--{self.diamond_sensitivity}'

            diamond.blastp()

            if not os.path.exists(tabular_output_path):
                raise ConfigError(f"DIAMOND blastp produced no output for '{db_name}'. "
                                  f"Check the log file at '{diamond_log_path}' for details.")

            functions_dict = self.parse_diamond_tabular_output(tabular_output_path, db_name)

            if not functions_dict:
                self.run.warning(f"No DIAMOND hits found for '{db_name}'. The gene_functions table will "
                                 f"not be updated for this source. Make sure your FASTA file contains "
                                 f"sequences related to the genes in the contigs database.")
                return

            source_name = f'{db_name}{DIAMOND_SOURCE_SUFFIX}'
            gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path,
                                                               run=self.run,
                                                               progress=self.progress)
            gene_function_calls_table.create(functions_dict)

            self.run.info('Results stored in', f"gene_functions table (source: '{source_name}')")

        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)


    def parse_diamond_tabular_output(self, tabular_output_path, db_name):
        """Parse DIAMOND blastp outfmt 6 output and return a gene_functions dict.

        The function field packs all search statistics so no information is lost
        even though the gene_functions schema has a fixed number of columns:
          function = "{target_id} [pident: {pct:.1f}%, aln_len: {length} aa, bitscore: {bitscore:.1f}]"

        Only the best hit per query (lowest e-value) is kept, which matches the
        diamond blastp --max-target-seqs 1 setting used during the search.

        Parameters
        ==========
        tabular_output_path : str
            Path to the diamond blastp tabular output (-outfmt 6).
        db_name : str
            User-defined database name.

        Returns
        =======
        functions_dict : dict
            Dict of {entry_id: {gene_callers_id, source, accession, function, e_value}}.
        """
        source_name = f'{db_name}{DIAMOND_SOURCE_SUFFIX}'

        # outfmt 6 columns: qseqid sseqid pident length mismatch gapopen
        #                   qstart qend sstart send evalue bitscore
        best_hits = {}

        with open(tabular_output_path) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue

                qseqid   = fields[0]
                sseqid   = fields[1]
                pident   = float(fields[2])
                aln_len  = int(fields[3])
                evalue   = float(fields[10])
                bitscore = float(fields[11])

                if qseqid in best_hits:
                    if evalue > best_hits[qseqid]['evalue']:
                        continue

                best_hits[qseqid] = {
                    'sseqid':   sseqid,
                    'pident':   pident,
                    'aln_len':  aln_len,
                    'evalue':   evalue,
                    'bitscore': bitscore,
                }

        functions_dict = {}
        entry_id = 0

        for gene_callers_id_str, hit in best_hits.items():
            try:
                gene_callers_id = int(gene_callers_id_str)
            except ValueError:
                continue

            sseqid   = hit['sseqid']
            pident   = hit['pident']
            aln_len  = hit['aln_len']
            evalue   = hit['evalue']
            bitscore = hit['bitscore']

            # The '[DMND] ' prefix mirrors the '[HMM] ' prefix added to HMM hits,
            # so any tabular export of the function field makes the search method
            # immediately visible without having to inspect the source column.
            function_str = (f"[DMND] {sseqid} "
                            f"[pident: {pident:.1f}%, "
                            f"aln_len: {aln_len} aa, "
                            f"bitscore: {bitscore:.1f}]")

            functions_dict[entry_id] = {
                'gene_callers_id': gene_callers_id,
                'source':          source_name,
                'accession':       sseqid,
                'function':        function_str,
                'e_value':         evalue,
            }
            entry_id += 1

        return functions_dict
