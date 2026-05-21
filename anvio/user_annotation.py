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
import re
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
from anvio import constants


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

HMM_SOURCE_SUFFIX = '_HMM'
DIAMOND_SOURCE_SUFFIX = '_DIAMOND'

# Inherit anvio standard evalue (matches anvio.K('min-e-value')['default'])
DEFAULT_HMM_EVALUE = 1e-15
DEFAULT_DIAMOND_EVALUE = 1e-15


def _auto_normalize_id(raw_id):
    """Auto-detect and normalize common protein database ID formats.

    Handles:
      - UniProt  sp|P12345|GENE_SPECIES  → gene name (e.g. RECA)
      - UniProt  tr|A0A000|GENE_SPECIES  → gene name
      - NCBI     ref|NP_123.1|           → NP_123  (strip version)
      - NCBI     gb|AAA12.1|             → AAA12
      - NCBI     WP_012345678.1          → WP_012345678
      - Pfam     PF00001.23              → PF00001
      - PDB      pdb|4HHB|A             → 4HHB_A
      - lcl|id                           → id
      - unknown pipe-delimited           → shortest informative field
      - plain ids                        → as-is
    """
    if not raw_id or raw_id == '-':
        return raw_id

    if '|' in raw_id:
        parts = raw_id.split('|')
        prefix = parts[0].lower()

        if prefix in ('sp', 'tr'):
            # UniProt SwissProt / TrEMBL: extract gene name from GENE_SPECIES token
            if len(parts) >= 3 and parts[2]:
                gene_species = parts[2].strip()
                return gene_species.split('_')[0] if '_' in gene_species else gene_species
            return (parts[1] or raw_id).strip()

        if prefix in ('ref', 'gb', 'emb', 'dbj', 'prf', 'pir', 'tpe', 'tpg', 'tpd'):
            acc = (parts[1] if len(parts) > 1 else raw_id).strip()
            # Strip numeric version suffix (e.g. NP_123456.1 → NP_123456)
            return re.sub(r'\.\d+$', '', acc) if acc else raw_id

        if prefix == 'pdb':
            struct = (parts[1] or '').strip()
            chain = (parts[2] if len(parts) > 2 else '').strip()
            return f"{struct}_{chain}" if chain else struct or raw_id

        if prefix == 'lcl':
            return (parts[1] if len(parts) > 1 else raw_id).strip()

        # Unknown pipe-delimited: take the first short alphanumeric-looking field
        for p in parts:
            p = p.strip()
            if p and 2 <= len(p) <= 25 and re.match(r'^[A-Za-z0-9_.:-]+$', p):
                return p
        return parts[-1].strip() or raw_id

    # Versioned NCBI-style accessions without pipes (WP_012345678.1, NP_123456.1)
    if re.match(r'^[A-Z]{1,3}[_]?\d+\.\d+$', raw_id):
        return re.sub(r'\.\d+$', '', raw_id)

    # Pfam / TIGRFAM version suffixes  PF00001.23 → PF00001
    if re.match(r'^(?:PF|TIGR)\d+\.\d+$', raw_id):
        return raw_id.rsplit('.', 1)[0]

    return raw_id


def _normalize_hmm_acc(acc):
    """Strip HMM file extensions from accession strings.

    File  GT5.hmm      → GT5
          GT5.hmm.gz   → GT5
    PF/TIGR version suffixes (PF00001.23, TIGR00001.1) are intentionally
    preserved — they carry version information that may matter for lookups.
    Other formats returned unchanged.
    """
    if not acc:
        return acc
    return re.sub(r'\.hmm(?:\.gz)?$', '', acc, flags=re.IGNORECASE)


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
        self.output_dir = A('output_dir') or constants.default_user_annotation_data_dir
        self.num_threads = A('num_threads') or 1
        self.reset = A('reset') or False
        self.list_databases_only = A('list') or False
        self.remove_database = A('remove')

        self.sanity_check()


    def sanity_check(self):
        if (self.list_databases_only or self.remove_database is not None) and self.input_tsv:
            flag = '--list' if self.list_databases_only else f'--remove {self.remove_database}'
            raise ConfigError(f"'{flag}' is a standalone operation and cannot be combined with --input-tsv. "
                              f"Use it on its own to inspect or modify an existing annotation directory "
                              f"(see `anvi-setup-user-annotation-db --help`).")

        if (self.list_databases_only or self.remove_database is not None) and self.reset:
            flag = '--list' if self.list_databases_only else '--remove'
            raise ConfigError(f"'{flag}' is a standalone operation and cannot be combined with --reset "
                              f"(see `anvi-setup-user-annotation-db --help`).")

        if self.list_databases_only or self.remove_database is not None:
            if not os.path.isdir(self.output_dir):
                raise ConfigError(f"The directory '{self.output_dir}' does not exist. "
                                  f"You need to run `anvi-setup-user-annotation-db --input-tsv <FILE>` "
                                  f"first to prepare your annotation databases.")
            manifest_path = os.path.join(self.output_dir, MANIFEST_FILENAME)
            if not os.path.exists(manifest_path):
                raise ConfigError(f"No manifest found in '{self.output_dir}'. "
                                  f"You need to run `anvi-setup-user-annotation-db --input-tsv <FILE>` "
                                  f"first to prepare your annotation databases, or point --output-dir "
                                  f"to an existing annotation directory.")
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
                cutoff_groups = entry.get('cutoff_groups', {})
                if len(cutoff_groups) > 1:
                    cutoff_str = 'mixed (' + ', '.join(f"{k}: {len(v)} models" for k, v in cutoff_groups.items()) + ')'
                else:
                    cutoff_str = entry.get('noise_cutoff_terms', '?')
                detail = (f"{entry.get('num_models', '?')} models | "
                          f"cutoff: {cutoff_str} | "
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
        manifest, manifest_path = self._load_manifest()
        failures, skipped = {}, []

        for name, info in entries.items():
            self._process_entry(name, info, manifest, manifest_path, failures, skipped)

        self._report_results(entries, failures, skipped, manifest_path)


    def _load_manifest(self):
        """Prepare output directory and load existing manifest. Returns (manifest dict, manifest path)."""
        if self.reset:
            filesnpaths.gen_output_directory(self.output_dir, run=self.run, delete_if_exists=True, dont_warn=True)
        else:
            filesnpaths.gen_output_directory(self.output_dir, run=self.run)

        manifest_path = os.path.join(self.output_dir, MANIFEST_FILENAME)
        manifest = {}
        if os.path.exists(manifest_path) and not self.reset:
            with open(manifest_path) as f:
                manifest = json.load(f)

        return manifest, manifest_path


    def _process_entry(self, name, info, manifest, manifest_path, failures, skipped):
        """Attempt setup for one TSV entry. Appends to failures or skipped on error/skip."""
        path = info['path']
        companion_fasta = info['companion_fasta']

        if name in manifest:
            self.run.warning(f"'{name}' is already in the manifest — skipping. "
                             f"Use `--remove {name}` first if you want to replace it.")
            skipped.append(name)
            return

        self.run.warning(None, header=f"SETTING UP: {name}", lc="green")

        try:
            db_type = self.validate_and_detect_db_type(path)
            self.run.info('Database name', name)
            self.run.info('Source path', path)
            self.run.info('Type', db_type)

            entry = self._setup_entry(name, path, db_type, companion_fasta)
            manifest[name] = entry

            # Write manifest after every successful entry so a failure mid-run does
            # not orphan already-prepared databases without a registry entry.
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)

        except (ConfigError, FilesNPathsError, Exception) as e:
            failures[name] = str(e)
            self.run.warning(f"Anvi'o could not set up '{name}' and will skip it. Error: {e}")


    def _setup_entry(self, name, path, db_type, companion_fasta):
        """Build and return manifest entry dict for one database (HMM or DIAMOND)."""
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
                entry['companion_diamond'] = self.setup_diamond_source(name, companion_fasta)
        else:
            if companion_fasta:
                raise ConfigError(f"A companion FASTA was provided for '{name}', but this entry is a "
                                  f"FASTA/DIAMOND database (not an HMM profile). The companion FASTA "
                                  f"column is only meaningful for HMM databases. Please remove the "
                                  f"companion FASTA path from the TSV row for '{name}' and try again.")
            entry = self.setup_diamond_source(name, path)

        return entry


    def _report_results(self, entries, failures, skipped, manifest_path):
        """Log summary of setup results."""
        self.run.info('Manifest written to', manifest_path)
        num_ok = len(entries) - len(failures) - len(skipped)
        newly_set_up = [k for k in entries if k not in failures and k not in skipped]
        self.run.info('Databases set up', f"{num_ok} ({', '.join(newly_set_up)})" if num_ok else "0")

        if skipped:
            self.run.info('Already in manifest (skipped)', ', '.join(skipped))

        if failures:
            self.run.warning(
                f"{len(failures)} database(s) failed and were skipped: "
                f"{', '.join(failures.keys())}. The manifest only contains the databases that succeeded. "
                f"Re-run setup for the failed databases after fixing the issues above."
            )


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


    def _copy_hmm_inject_acc(self, src_path, dst_gz_path, models_no_acc):
        """Copy HMM to dst_gz_path (gzipped), injecting 'ACC  {name}' after NAME for models lacking it."""
        models_no_acc_set = set(models_no_acc)
        opener_in = gzip.open if src_path.lower().endswith('.gz') else open

        with opener_in(src_path, 'rt', encoding='utf-8', errors='replace') as f_in, \
             gzip.open(dst_gz_path, 'wt', encoding='utf-8') as f_out:
            for line in f_in:
                f_out.write(line)
                if line.startswith('NAME'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] in models_no_acc_set:
                        f_out.write(f"ACC  {parts[1]}\n")


    def _sniff_file_type(self, file_path):
        """Peek at the file content and return 'hmm', 'diamond', or None."""
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
        """Determine whether a file is an HMM profile or a FASTA file, and verify the content."""
        lower = file_path.lower()
        basename = os.path.basename(file_path)

        if lower.endswith('.hmm') or lower.endswith('.hmm.gz'):
            extension_type = 'hmm'
        elif any(lower.endswith(ext) for ext in ['.faa', '.fasta', '.fa', '.fas', '.fna']):
            extension_type = 'diamond'
        else:
            extension_type = None

        content_type = self._sniff_file_type(file_path)

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

        Each dict contains: name, acc, has_tc, has_ga, has_nc, cutoff_type.
        cutoff_type is the best available cutoff annotation ('tc', 'ga', 'nc', or None).
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
                        models.append(self._finalize_hmm_model(current, line_num, hmm_path))
                        current = None
        except (OSError, EOFError) as e:
            raise ConfigError(f"Could not read '{hmm_path}': {e}. Is this a valid HMMER3 profile file?")

        if not models:
            raise ConfigError(f"No HMMER3 models found in '{hmm_path}'. A valid HMMER3 profile starts "
                              f"with a 'HMMER3/f' header and contains at least one model terminated by '//'.")

        self._check_duplicate_model_names(models, hmm_path)
        return models


    def _finalize_hmm_model(self, current, line_num, hmm_path):
        """Validate a completed HMM model record and assign its cutoff_type. Returns the model dict."""
        if current['name'] is None:
            raise ConfigError(f"Found an HMM model without a NAME field near line {line_num} "
                              f"in '{hmm_path}'. Every model in a valid HMMER3 profile must "
                              f"have a NAME line.")
        if current['has_tc']:
            current['cutoff_type'] = 'tc'
        elif current['has_ga']:
            current['cutoff_type'] = 'ga'
        elif current['has_nc']:
            current['cutoff_type'] = 'nc'
        else:
            current['cutoff_type'] = None
        return current


    def _check_duplicate_model_names(self, models, hmm_path):
        """Raise ConfigError if any model NAME appears more than once in the parsed model list."""
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


    def determine_noise_cutoff(self, models, db_name):
        """Select the best HMMER noise cutoff strategy for a set of models.

        Returns (cutoff_flag, cutoff_groups) where:
          cutoff_flag   – single flag string when all models share a cutoff type, else None
          cutoff_groups – dict mapping cutoff flag to list of model names
                          (one key when uniform, multiple keys when mixed)

        For mixed files, each model is assigned to the best available cutoff group
        rather than falling back globally to the evalue threshold.
        """
        n = len(models)
        num_tc = sum(1 for m in models if m['has_tc'])
        num_ga = sum(1 for m in models if m['has_ga'])
        num_nc = sum(1 for m in models if m['has_nc'])

        if num_tc == n:
            return '--cut_tc', {'--cut_tc': [m['name'] for m in models]}
        if num_ga == n:
            return '--cut_ga', {'--cut_ga': [m['name'] for m in models]}
        if num_nc == n:
            return '--cut_nc', {'--cut_nc': [m['name'] for m in models]}

        # Each model gets its best available cutoff (TC > GA > NC > evalue fallback)
        fallback_flag = f'-E {DEFAULT_HMM_EVALUE}'
        cutoff_groups = {}
        for m in models:
            if m['has_tc']:
                key = '--cut_tc'
            elif m['has_ga']:
                key = '--cut_ga'
            elif m['has_nc']:
                key = '--cut_nc'
            else:
                key = fallback_flag
            cutoff_groups.setdefault(key, []).append(m['name'])

        if len(cutoff_groups) == 1:
            # All models share the same (possibly fallback) cutoff — not truly mixed
            single_flag = list(cutoff_groups.keys())[0]
            if single_flag == fallback_flag:
                self.run.warning(
                    f"No TC, GA, or NC cutoff annotations found in any model of '{db_name}'. "
                    f"Using {fallback_flag} for all {n} model(s). To use profile-specific cutoffs, "
                    f"add TC/GA/NC lines to your HMM profiles before running setup."
                )
            return single_flag, cutoff_groups

        # Genuinely mixed: multiple cutoff types across models
        group_summary = ', '.join(f"{k}: {len(v)} model(s)" for k, v in cutoff_groups.items())
        self.run.info_single(
            f"Models in '{db_name}' have mixed cutoff types. Anvi'o will run each group separately "
            f"at annotation time using the optimal cutoff for each: {group_summary}.",
            cut_after=None
        )

        return None, cutoff_groups


    def setup_hmm_source(self, db_name, hmm_path):
        """Create a standard anvi'o HMM source directory from a raw HMMER3 profile file."""
        self.progress.new(f"Setting up HMM source '{db_name}'")
        self.progress.update("Parsing models...")

        models = self.parse_hmm_models(hmm_path)
        cutoff_flag, cutoff_groups = self.determine_noise_cutoff(models, db_name)

        # noise_cutoff_terms.txt in the main dir uses the primary flag; for mixed
        # files we write the fallback so the dir is still a valid HMM source.
        primary_cutoff = cutoff_flag or f'-E {DEFAULT_HMM_EVALUE}'

        hmm_dir = os.path.join(self.output_dir, HMM_SUBDIR, db_name)
        filesnpaths.gen_output_directory(hmm_dir, run=self.run, delete_if_exists=True, dont_warn=True)

        models_no_acc = [m['name'] for m in models if not m['acc']]
        if models_no_acc:
            self.run.warning(
                f"{len(models_no_acc)} model(s) in '{db_name}' have no ACC line in the HMM profile: "
                f"{', '.join(models_no_acc)}. Anvi'o will use the model name as the accession for these models."
            )

        self.progress.update("Writing genes.txt...")
        with open(os.path.join(hmm_dir, 'genes.txt'), 'w') as f:
            f.write('gene\taccession\thmmsource\n')
            for m in models:
                raw_acc = m['acc'] if m['acc'] else m['name']
                acc = _normalize_hmm_acc(raw_acc)
                f.write(f"{m['name']}\t{acc}\tuser-provided\n")

        self.progress.update("Writing metadata files...")
        with open(os.path.join(hmm_dir, 'kind.txt'), 'w') as f:
            f.write('user_annotation')

        with open(os.path.join(hmm_dir, 'target.txt'), 'w') as f:
            f.write('AA:GENE')

        with open(os.path.join(hmm_dir, 'noise_cutoff_terms.txt'), 'w') as f:
            f.write(primary_cutoff)

        with open(os.path.join(hmm_dir, 'reference.txt'), 'w') as f:
            f.write(f"User-provided HMM database. Original file: {hmm_path}")

        self.progress.update("Copying and compressing HMM file...")
        hmm_gz_path = os.path.join(hmm_dir, 'genes.hmm.gz')
        if models_no_acc:
            self._copy_hmm_inject_acc(hmm_path, hmm_gz_path, models_no_acc)
        elif hmm_path.lower().endswith('.gz'):
            shutil.copy2(hmm_path, hmm_gz_path)
        else:
            with open(hmm_path, 'rb') as f_in:
                with gzip.open(hmm_gz_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        self.progress.end()

        self.run.info('HMM models found', len(models))
        if cutoff_flag:
            self.run.info('Noise cutoff strategy', cutoff_flag)
        else:
            self.run.info('Noise cutoff strategy', f"mixed ({len(cutoff_groups)} groups — applied dynamically at run time)")
        self.run.info('HMM source directory', hmm_dir)

        return {
            'type': 'hmm',
            'source_path': hmm_path,
            'hmm_dir': hmm_dir,
            'num_models': len(models),
            'noise_cutoff_terms': primary_cutoff,
            'cutoff_groups': cutoff_groups,
            'added_on': str(datetime.date.today()),
        }


    def _build_sequence_id_mapping(self, fasta_path):
        """Parse FASTA headers and return {sseqid: {'clean_id': ..., 'description': ...}}.

        sseqid    = FASTA header up to first space (what DIAMOND reports in column 2)
        clean_id  = auto-normalized version of sseqid
        description = rest of header after first space (human-readable protein name)
        """
        mapping = {}
        opener = gzip.open if fasta_path.lower().endswith(('.gz', '.gzip')) else open
        try:
            with opener(fasta_path, 'rt', encoding='utf-8', errors='replace') as f:
                for line in f:
                    if not line.startswith('>'):
                        continue
                    header = line[1:].strip()
                    if ' ' in header:
                        sseqid, description = header.split(' ', 1)
                    else:
                        sseqid, description = header, ''
                    mapping[sseqid] = {
                        'clean_id': _auto_normalize_id(sseqid),
                        'description': description.strip(),
                    }
        except (OSError, EOFError) as e:
            self.run.warning(f"Could not parse FASTA headers from '{fasta_path}' for ID normalization: {e}. "
                             f"Raw sequence IDs will be used as accessions.")
        return mapping


    def setup_diamond_source(self, db_name, fasta_path):
        """Build a DIAMOND database from a protein FASTA file."""
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

        # Build and store ID-normalization mapping from FASTA headers
        self.progress.new(f"Building ID mapping for '{db_name}'")
        self.progress.update("Parsing FASTA headers...")
        id_mapping = self._build_sequence_id_mapping(fasta_path)
        self.progress.end()

        num_seqs = len(id_mapping)

        mapping_path = os.path.join(diamond_dir, f'{db_name}_id_mapping.json')
        with open(mapping_path, 'w') as f:
            json.dump(id_mapping, f)

        has_descriptions = sum(1 for v in id_mapping.values() if v['description'])

        self.run.info('DIAMOND database', dmnd_path)
        self.run.info('Sequences in source FASTA', num_seqs or 'unknown')
        self.run.info('IDs with descriptions', f"{has_descriptions}/{len(id_mapping)}")

        return {
            'type': 'diamond',
            'source_path': fasta_path,
            'dmnd_path': dmnd_path,
            'dmnd_base': dmnd_base,
            'num_sequences': num_seqs,
            'id_mapping_path': mapping_path,
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
        self.annotation_db_dir = A('annotation_db_dir') or constants.default_user_annotation_data_dir
        self.database = A('database')
        self.num_threads = A('num_threads') or 1
        self.just_do_it = A('just_do_it') or False
        self.force_overwrite = A('force_overwrite') or False
        self.hmmer_program = A('hmmer_program') or 'hmmscan'
        self.evalue = A('evalue') if A('evalue') is not None else DEFAULT_DIAMOND_EVALUE
        self.min_pident = A('min_pident')
        self.qcov = A('qcov')
        self.max_hsps = A('max_hsps')
        self.diamond_sensitivity = A('diamond_sensitivity')
        self.custom_tc_map = self._parse_custom_tc_tsv(A('cut_tc')) if A('cut_tc') else {}

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


    def _handle_existing_source(self, source_name):
        """Check whether source_name already exists in gene_functions.

        Returns True if the caller should proceed (source absent, or force_overwrite dropped it).
        Returns False if the caller should skip (source exists and force_overwrite not set).
        """
        contigs_database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))
        existing_raw = contigs_database.get_meta_value('gene_function_sources', return_none_if_not_in_table=True)
        existing = set(existing_raw.split(',') if existing_raw else [])

        if source_name not in existing:
            contigs_database.disconnect()
            return True

        if self.force_overwrite:
            self.run.warning(
                f"Source '{source_name}' already in the gene_functions table. "
                f"Dropping and re-annotating (--force-overwrite was passed)."
            )
            gf_table = TableForGeneFunctions(self.contigs_db_path,
                                              run=terminal.Run(verbose=False),
                                              progress=terminal.Progress(verbose=False))
            gf_table.drop_functions(contigs_database, sources_to_drop=[source_name])
            contigs_database.disconnect()
            return True

        self.run.warning(
            f"Source '{source_name}' is already in the gene_functions table. Skipping. "
            f"Pass --force-overwrite to drop the existing annotation and re-annotate."
        )
        contigs_database.disconnect()
        return False


    def _parse_custom_tc_tsv(self, path):
        """Parse a custom TC file into {db_or_None: {model_name: (seq_tc, dom_tc)}}.

        Optional first column is a database name (must not be parseable as a number).
        Entries with a db name apply only to that database; entries without apply to all.

        Formats accepted per line (auto-detected):
          model<TAB>seq_tc[<TAB>dom_tc]             — no db restriction
          db<TAB>model<TAB>seq_tc[<TAB>dom_tc]      — db-specific (tab-delimited)
          db<TAB>model: seq_tc[ dom_tc]              — db-specific with colon model expr
          model: seq_tc[ dom_tc]                     — no db restriction, colon format
          e.g. "[FeFe]: 15.9"  or  "HydDB\t[FeFe]: 15.9"

        dom_tc defaults to seq_tc when omitted. Lines starting with '#' and blank lines ignored.
        """
        result = {}
        if not os.path.isfile(path):
            raise FilesNPathsError(f"Custom TC file not found: '{path}'")

        total = 0
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parsed = self._parse_tc_line(line)
                if parsed is None:
                    continue
                db_name, name, seq_tc, dom_tc = parsed
                result.setdefault(db_name, {})[name] = (seq_tc, dom_tc)
                total += 1

        global_n = sum(len(v) for k, v in result.items() if k is None)
        db_n = total - global_n
        self.run.info(
            "Custom TC overrides loaded",
            f"{total} model(s) from '{path}' "
            f"({global_n} global, {db_n} db-specific)"
        )
        return result


    def _parse_tc_line(self, line):
        """Parse one line from a custom TC file. Returns (db_name, model_name, seq_tc, dom_tc) or None to skip."""
        db_name = None

        if '\t' in line:
            parts = line.split('\t')
            # Detect db column: col[1] not parseable as float → first col is a db name.
            has_db_col = False
            if len(parts) >= 2:
                try:
                    float(parts[1])
                except ValueError:
                    has_db_col = True

            if has_db_col:
                db_name = parts[0].strip()
                rest = '\t'.join(parts[1:])
            else:
                rest = line

            if '\t' in rest:
                rest_parts = rest.split('\t')
                name = rest_parts[0].strip()
                values = [p.strip() for p in rest_parts[1:] if p.strip()]
            else:
                colon_pos = rest.rfind(':')
                if colon_pos != -1:
                    name = rest[:colon_pos].strip()
                    values = rest[colon_pos + 1:].strip().split()
                else:
                    name = rest.strip()
                    values = []
        else:
            colon_pos = line.rfind(':')
            if colon_pos == -1:
                return None
            name = line[:colon_pos].strip()
            values = line[colon_pos + 1:].strip().split()

        if not values:
            return None
        try:
            seq_tc = float(values[0])
            dom_tc = float(values[1]) if len(values) >= 2 else seq_tc
        except ValueError:
            return None

        return db_name, name, seq_tc, dom_tc


    def _apply_custom_tc_to_groups(self, db_name, cutoff_groups):
        """Rewrite cutoff_groups to move custom-TC models into a dedicated group.

        db_name is the user-defined database name (used to look up db-specific entries).
        DB-specific entries take precedence over global (no-db) entries for the same model.

        Returns (custom_tc_models, new_cutoff_groups) where:
          custom_tc_models  – {model_name: (seq_tc, dom_tc)} for models found in this DB
          new_cutoff_groups – rebuilt cutoff_groups with a '--cut_tc_custom' key for overrides
        """
        # Merge: global entries first, then db-specific entries override
        relevant = {}
        relevant.update(self.custom_tc_map.get(None, {}))
        relevant.update(self.custom_tc_map.get(db_name, {}))

        all_models = {name for names in cutoff_groups.values() for name in names}
        custom_tc_models = {m: v for m, v in relevant.items() if m in all_models}

        # Warn about entries intended for this db that don't match any model
        unmatched = sorted(set(relevant) - set(custom_tc_models))
        if unmatched:
            self.run.warning(
                f"{len(unmatched)} name(s) from your --cut-tc file were not found among the models "
                f"of '{db_name}' and will be ignored: {', '.join(unmatched)}. "
                f"Use `grep '^NAME' your_file.hmm` or inspect genes.txt in the HMM source "
                f"directory to find the exact model names."
            )

        if not custom_tc_models:
            return {}, cutoff_groups

        custom_set = set(custom_tc_models)
        new_groups = {}
        for flag, names in cutoff_groups.items():
            remaining = [n for n in names if n not in custom_set]
            if remaining:
                new_groups[flag] = remaining
        new_groups['--cut_tc_custom'] = list(custom_tc_models)

        overridden = ', '.join(sorted(custom_tc_models))
        self.run.info_single(
            f"Custom TC overrides applied to {len(custom_tc_models)} model(s) in '{db_name}': {overridden}.",
            cut_after=None
        )
        return custom_tc_models, new_groups


    def _extract_models_from_hmm(self, hmm_gz_path, model_names, tc_overrides=None):
        """Extract a subset of models from a gzipped HMM file into a temp plain HMM file.

        tc_overrides: optional {model_name: (seq_tc, dom_tc)} — injects or replaces TC lines
        for specified models so they carry the custom values when searched with --cut_tc.

        Returns the path to the temporary (plain, uncompressed) HMM file.
        Caller is responsible for deleting it.
        """
        model_names_set = set(model_names)
        tc_overrides = tc_overrides or {}
        temp_dir = filesnpaths.get_temp_directory_path()
        temp_hmm_path = os.path.join(temp_dir, 'extracted.hmm')

        opener = gzip.open if hmm_gz_path.lower().endswith('.gz') else open

        with opener(hmm_gz_path, 'rt', encoding='utf-8', errors='replace') as f_in, \
             open(temp_hmm_path, 'w') as f_out:
            in_model = False
            keep = False
            current_name = None
            current_lines = []
            saw_tc = False

            for line in f_in:
                if line.startswith('HMMER3/'):
                    in_model = True
                    keep = False
                    current_name = None
                    current_lines = [line]
                    saw_tc = False
                elif in_model:
                    if line.startswith('NAME') and len(line.split()) >= 2:
                        current_name = line.split()[1]
                        keep = current_name in model_names_set
                        current_lines.append(line)
                    elif line.startswith('TC ') or line.startswith('TC\t'):
                        saw_tc = True
                        current_lines.append(self._get_tc_line(current_name, tc_overrides, line))
                    elif line.rstrip() == 'HMM' or line.startswith('HMM '):
                        saw_tc = self._inject_missing_tc(keep, current_name, tc_overrides, saw_tc, current_lines)
                        current_lines.append(line)
                    elif line.strip() == '//':
                        current_lines.append(line)
                        if keep:
                            f_out.writelines(current_lines)
                        in_model = False
                        keep = False
                        current_name = None
                        current_lines = []
                        saw_tc = False
                    else:
                        current_lines.append(line)

        return temp_hmm_path


    def _get_tc_line(self, current_name, tc_overrides, original_line):
        """Return TC line to write: formatted override if present, else the original line."""
        if current_name in tc_overrides:
            seq_tc, dom_tc = tc_overrides[current_name]
            return f'TC    {seq_tc:.2f}  {dom_tc:.2f};\n'
        return original_line


    def _inject_missing_tc(self, keep, current_name, tc_overrides, saw_tc, current_lines):
        """Inject TC line before HMM parameters if an override exists and TC was not yet written. Returns updated saw_tc."""
        if keep and current_name in tc_overrides and not saw_tc:
            seq_tc, dom_tc = tc_overrides[current_name]
            current_lines.append(f'TC    {seq_tc:.2f}  {dom_tc:.2f};\n')
            return True
        return saw_tc


    def _run_hmm_annotation_single_cutoff(self, internal_source_name, hmm_dir):
        """Run HMM annotation for a database with a single uniform cutoff."""
        from anvio.tables.hmmhits import TablesForHMMHits

        sources = utils.get_HMM_sources_dictionary([hmm_dir])
        original_key = list(sources.keys())[0]
        sources[internal_source_name] = sources.pop(original_key)

        search_tables = TablesForHMMHits(self.contigs_db_path,
                                         num_threads_to_use=self.num_threads,
                                         run=self.run,
                                         progress=self.progress,
                                         just_do_it=True,
                                         hmm_program_to_use=self.hmmer_program,
                                         add_to_functions_table=True)
        search_tables.populate_search_tables(sources)


    def _run_hmm_annotation_mixed_cutoffs(self, db_name, entry, cutoff_groups, internal_source_name, custom_tc_models=None):
        """Run HMM annotation for a database with mixed cutoff types (or custom TC overrides).

        Each cutoff group is searched separately with its optimal cutoff flag,
        then all results are merged under internal_source_name in the database.
        Groups keyed '--cut_tc_custom' use injected TC values and are searched with --cut_tc.
        """
        custom_tc_models = custom_tc_models or {}
        from anvio.tables.hmmhits import TablesForHMMHits

        hmm_dir = entry['hmm_dir']
        hmm_gz_path = os.path.join(hmm_dir, 'genes.hmm.gz')

        # Read genes.txt once for filtering
        genes_txt_path = os.path.join(hmm_dir, 'genes.txt')
        genes_by_name = {}
        with open(genes_txt_path) as f:
            next(f)  # skip header
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if parts:
                    genes_by_name[parts[0]] = line

        tmp_source_names = []
        temp_dirs = []

        try:
            for i, (cutoff_flag, model_names) in enumerate(cutoff_groups.items()):
                tmp_source = f'uatmp{i}x{db_name[:8]}'.replace('-', '').replace('.', '')[:20]
                # Ensure PROPER: starts with letter, only alnum+underscore
                if not tmp_source[0].isalpha():
                    tmp_source = 'ua' + tmp_source
                tmp_source_names.append(tmp_source)

                is_custom_tc = cutoff_flag == '--cut_tc_custom'
                display_flag = '--cut_tc (custom values injected)' if is_custom_tc else cutoff_flag
                self.run.info(f"Running cutoff group '{display_flag}'",
                              f"{len(model_names)} model(s)")

                # Extract models; inject custom TC values when this is a custom-TC group
                tc_overrides_for_group = custom_tc_models if is_custom_tc else None
                extracted_hmm = self._extract_models_from_hmm(hmm_gz_path, model_names, tc_overrides=tc_overrides_for_group)

                # Create temp HMM source dir named tmp_source (basename = source key)
                temp_parent = filesnpaths.get_temp_directory_path()
                temp_dirs.append(temp_parent)
                tmp_source_dir = os.path.join(temp_parent, tmp_source)
                os.makedirs(tmp_source_dir)

                # Write genes.txt filtered to this group's models
                model_names_set = set(model_names)
                with open(os.path.join(tmp_source_dir, 'genes.txt'), 'w') as f:
                    f.write('gene\taccession\thmmsource\n')
                    for gname, gline in genes_by_name.items():
                        if gname in model_names_set:
                            f.write(gline)

                # Copy static metadata files
                for fname in ['kind.txt', 'target.txt', 'reference.txt']:
                    shutil.copy2(os.path.join(hmm_dir, fname), os.path.join(tmp_source_dir, fname))

                # Custom-TC group runs with --cut_tc (TC values are now baked into the extracted HMM)
                noise_cutoff_for_file = '--cut_tc' if is_custom_tc else cutoff_flag
                with open(os.path.join(tmp_source_dir, 'noise_cutoff_terms.txt'), 'w') as f:
                    f.write(noise_cutoff_for_file)

                # Compress extracted HMM into the source dir
                with open(extracted_hmm, 'rb') as f_in:
                    with gzip.open(os.path.join(tmp_source_dir, 'genes.hmm.gz'), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(extracted_hmm)

                sources = utils.get_HMM_sources_dictionary([tmp_source_dir])

                search_tables = TablesForHMMHits(self.contigs_db_path,
                                                 num_threads_to_use=self.num_threads,
                                                 run=self.run,
                                                 progress=self.progress,
                                                 just_do_it=True,
                                                 hmm_program_to_use=self.hmmer_program,
                                                 add_to_functions_table=True)
                search_tables.populate_search_tables(sources)

            # Merge all tmp source names into internal_source_name
            contigs_database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))
            placeholders = ','.join('?' * len(tmp_source_names))

            contigs_database._exec(
                f"UPDATE {t.gene_function_calls_table_name} SET source = ? WHERE source IN ({placeholders})",
                [internal_source_name] + tmp_source_names
            )
            contigs_database._exec(
                f"UPDATE {t.hmm_hits_table_name} SET source = ? WHERE source IN ({placeholders})",
                [internal_source_name] + tmp_source_names
            )

            # Merge hmm_hits_info rows: collect unique genes, delete tmp rows, insert one merged row
            cursor = contigs_database._exec(
                f"SELECT DISTINCT gene_name FROM {t.hmm_hits_table_name} WHERE source = ?",
                [internal_source_name]
            )
            all_gene_names = [row[0] for row in cursor.fetchall()]

            cursor2 = contigs_database._exec(
                f"SELECT ref, search_type, domain FROM {t.hmm_hits_info_table_name} WHERE source IN ({placeholders})",
                tmp_source_names
            )
            info_rows = cursor2.fetchall()

            contigs_database._exec(
                f"DELETE FROM {t.hmm_hits_info_table_name} WHERE source IN ({placeholders})",
                tmp_source_names
            )

            existing_info = contigs_database._exec(
                f"SELECT source FROM {t.hmm_hits_info_table_name} WHERE source = ?",
                [internal_source_name]
            ).fetchone()

            if not existing_info and info_rows:
                ref = info_rows[0][0]
                search_type = info_rows[0][1]
                domain = info_rows[0][2]
                genes_str = ','.join(sorted(set(all_gene_names)))
                contigs_database._exec(
                    f"INSERT INTO {t.hmm_hits_info_table_name} VALUES (?,?,?,?,?)",
                    [internal_source_name, ref, search_type, domain, genes_str]
                )

            # Update gene_function_sources meta value
            existing_raw = contigs_database.get_meta_value('gene_function_sources', return_none_if_not_in_table=True)
            existing_set = set(existing_raw.split(',') if existing_raw else [])
            existing_set -= set(tmp_source_names)
            existing_set.add(internal_source_name)
            contigs_database.remove_meta_key_value_pair('gene_function_sources')
            contigs_database.set_meta_value('gene_function_sources', ','.join(existing_set))

            contigs_database.disconnect()

        finally:
            for tmp_dir in temp_dirs:
                shutil.rmtree(tmp_dir, ignore_errors=True)


    def run_hmm_annotation(self, db_name, entry):
        """Search an HMM source and store hits in the gene_functions table.

        Results are written with source name '{db_name}_HMM'. If the source already
        exists, it is dropped and replaced (same behavior as other anvio annotation programs).
        For databases with mixed cutoff types, each group is searched separately with its
        optimal cutoff flag and results are merged under the final source name.

        Parameters
        ==========
        db_name : str
            User-defined database name.
        entry : dict
            Manifest entry for this database.
        """
        hmm_dir = entry['hmm_dir']
        cutoff_groups = entry.get('cutoff_groups', {})
        internal_source_name = f'{db_name}{HMM_SOURCE_SUFFIX}'

        if not os.path.isdir(hmm_dir):
            raise ConfigError(f"The HMM source directory for '{db_name}' no longer exists at '{hmm_dir}'. "
                              f"Please run `anvi-setup-user-annotation-db` again.")

        if not self._handle_existing_source(internal_source_name):
            return

        self.run.warning("Anvi'o will use 'HMMER' by Eddy (doi:10.1371/journal.pcbi.1002195) to search your "
                         "contigs database against the HMM profiles in this database. When you publish your "
                         "findings, please do not forget to properly credit their work.", lc='green', header="CITATION")

        custom_tc_models = {}
        if self.custom_tc_map:
            custom_tc_models, cutoff_groups = self._apply_custom_tc_to_groups(db_name, cutoff_groups)

        if len(cutoff_groups) > 1 or custom_tc_models:
            self._run_hmm_annotation_mixed_cutoffs(db_name, entry, cutoff_groups, internal_source_name, custom_tc_models)
        else:
            self._run_hmm_annotation_single_cutoff(internal_source_name, hmm_dir)

        contigs_database = db.DB(self.contigs_db_path, utils.get_required_version_for_db(self.contigs_db_path))

        # Collect models with missing accession BEFORE prefixing function field
        cursor = contigs_database._exec(
            f"SELECT DISTINCT function FROM {t.gene_function_calls_table_name} "
            f"WHERE source = ? AND accession = '-'",
            [internal_source_name]
        )
        missing_acc_models = [row[0] for row in cursor.fetchall()]
        if missing_acc_models:
            self.run.warning(
                f"{len(missing_acc_models)} model(s) reported no accession ('-') from HMMER for source "
                f"'{internal_source_name}': {', '.join(missing_acc_models)}. "
                f"Anvi'o will use the model name as the accession for these hits."
            )

        # Normalize accessions and function names: strip .hmm/.hmm.gz file extensions that
        # may have come through from the HMM model's ACC or NAME field via HMMER tblout output.
        for col in ('accession', 'function'):
            cursor = contigs_database._exec(
                f"SELECT DISTINCT {col} FROM {t.gene_function_calls_table_name} WHERE source = ?",
                [internal_source_name]
            )
            for (raw_val,) in cursor.fetchall():
                normalized = _normalize_hmm_acc(raw_val)
                if normalized != raw_val:
                    contigs_database._exec(
                        f"UPDATE {t.gene_function_calls_table_name} "
                        f"SET {col} = ? WHERE source = ? AND {col} = ?",
                        [normalized, internal_source_name, raw_val]
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

        Results are written with source name '{db_name}_DIAMOND'. If the source already
        exists, TableForGeneFunctions.create() will warn and replace it automatically.

        Parameters
        ==========
        db_name : str
            User-defined database name.
        entry : dict
            Manifest entry for this database.
        """
        dmnd_path = entry.get('dmnd_path', '')
        dmnd_base = entry.get('dmnd_base', dmnd_path.removesuffix('.dmnd'))

        source_name = f'{db_name}{DIAMOND_SOURCE_SUFFIX}'
        if not self._handle_existing_source(source_name):
            return

        self.run.warning("Anvi'o will use 'DIAMOND' by Buchfink et al. (doi:10.1038/nmeth.3176) to search your "
                         "contigs database against this protein database. When you publish your findings, "
                         "please do not forget to properly credit their work.", lc='green', header="CITATION")

        if not os.path.exists(dmnd_path):
            raise ConfigError(f"The DIAMOND database for '{db_name}' no longer exists at '{dmnd_path}'. "
                              f"Please run `anvi-setup-user-annotation-db` again.")

        temp_dir = filesnpaths.get_temp_directory_path()

        try:
            aa_sequences_path = os.path.join(temp_dir, 'aa_sequences.faa')

            self.progress.new("Exporting AA sequences from contigs DB")
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
                              num_threads=self.num_threads,
                              outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen')
            diamond.target_fasta = dmnd_base
            diamond.tabular_output_path = tabular_output_path
            diamond.max_target_seqs = 1
            logs_dir = os.path.join(self.annotation_db_dir, 'logs')
            filesnpaths.gen_output_directory(logs_dir, run=terminal.Run(verbose=False))
            diamond_log_path = os.path.join(logs_dir, f'{db_name}_diamond_blastp.log')
            diamond.run.log_file_path = diamond_log_path

            diamond.evalue = self.evalue

            if self.min_pident is not None:
                diamond.min_pct_id = self.min_pident

            extra_params = []
            if self.diamond_sensitivity:
                extra_params.append(f'--{self.diamond_sensitivity}')
            if self.qcov is not None:
                extra_params.append(f'--query-cover {self.qcov}')
            if self.max_hsps is not None:
                extra_params.append(f'--max-hsps {self.max_hsps}')
            if extra_params:
                diamond.additional_params_for_blastp = ' '.join(extra_params)

            diamond.blastp()

            if not os.path.exists(tabular_output_path):
                raise ConfigError(f"DIAMOND blastp produced no output for '{db_name}'. "
                                  f"Check the log file at '{diamond_log_path}' for details.")

            # Load ID-normalization mapping if available
            id_mapping = {}
            mapping_path = entry.get('id_mapping_path')
            if mapping_path and os.path.exists(mapping_path):
                with open(mapping_path) as f:
                    id_mapping = json.load(f)

            functions_dict = self.parse_diamond_tabular_output(tabular_output_path, db_name, id_mapping)

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


    def parse_diamond_tabular_output(self, tabular_output_path, db_name, id_mapping=None):
        """Parse DIAMOND blastp outfmt 6 output and return a gene_functions dict.

        Uses id_mapping (built at setup time from FASTA headers) to:
          - normalize the accession (clean_id from auto-detected format)
          - populate the function field with the protein description from the FASTA header

        Schema stored per hit:
          accession = auto-normalized ID  (e.g. gene name, stripped accession)
          function  = '[DMND] <description or clean_id> [pident: X%, aln_len: Y aa, bitscore: Z]'

        Only the best hit per query (lowest e-value) is kept.
        """
        source_name = f'{db_name}{DIAMOND_SOURCE_SUFFIX}'
        if id_mapping is None:
            id_mapping = {}

        # outfmt 6 columns: qseqid sseqid pident length mismatch gapopen
        #                   qstart qend sstart send evalue bitscore
        best_hits = {}

        with open(tabular_output_path) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue

                qseqid = fields[0]
                sseqid = fields[1]
                pident = float(fields[2])
                aln_len = int(fields[3])
                evalue = float(fields[10])
                bitscore = float(fields[11])
                qlen = int(fields[12]) if len(fields) > 12 and fields[12] else None

                if qseqid in best_hits and evalue >= best_hits[qseqid]['evalue']:
                    continue

                best_hits[qseqid] = {
                    'sseqid': sseqid,
                    'pident': pident,
                    'aln_len': aln_len,
                    'evalue': evalue,
                    'bitscore': bitscore,
                    'qlen': qlen,
                }

        functions_dict = {}
        entry_id = 0

        for gene_callers_id_str, hit in best_hits.items():
            try:
                gene_callers_id = int(gene_callers_id_str)
            except ValueError:
                continue

            sseqid = hit['sseqid']
            pident = hit['pident']
            aln_len = hit['aln_len']
            evalue = hit['evalue']
            bitscore = hit['bitscore']

            meta = id_mapping.get(sseqid, {})
            clean_id    = meta.get('clean_id', _auto_normalize_id(sseqid))
            description = meta.get('description', '')

            label = description if description else clean_id
            qlen = hit.get('qlen')
            qcov_str = f", qcov: {aln_len / qlen * 100:.1f}%" if qlen else ""
            function_str = (f"{label}"
                            #f" [pident: {pident:.1f}%, "
                            #f"aln_len: {aln_len} aa, "
                            #f"bitscore: {bitscore:.1f}"
                            #f"{qcov_str}]"
            )

            functions_dict[entry_id] = {
                'gene_callers_id': gene_callers_id,
                'source': source_name,
                'accession': clean_id,
                'function': function_str,
                'e_value': evalue,
            }
            entry_id += 1

        return functions_dict
