# pylint: disable=line-too-long
"""
Generic loader for the anvi'o `structures-txt` artifact: a two-column TSV that
maps a gene identifier to the path of a predicted-protein-structure file.

Two consumers exist today:
- `anvi-pan-genome` (structure-informed pangenomics) where the IDs are gene
  cluster IDs (the FASTA defline of the GC representative).
- `anvi-gen-structure-database` (single-genome structure DB) where the IDs are
  gene-caller-ids in a contigs-db. That consumer subclasses StructuresTxt as
  `ExternalStructuresFile` and adds contigs-db-aware sanity checks.

The legacy column header `gene_callers_id` from the `external-structures`
artifact is accepted as an alias of `gene_id`.
"""

import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from pathlib import Path
from anvio.errors import ConfigError, FilesNPathsError


class StructuresTxt:
    REQUIRED_HEADERS = ('gene_id', 'path')
    LEGACY_ID_HEADERS = ('gene_callers_id',)
    ALLOWED_EXTENSIONS = ('.pdb', '.cif', '.mmcif', '.pdb.gz', '.cif.gz', '.mmcif.gz', '.fcz')

    def __init__(self, path, run=None, progress=None, skip_check_sanity=False):
        self.run = run or terminal.Run()
        self.progress = progress or terminal.Progress()

        self.path = str(path)

        # Populated by _parse(): gene_id (str) -> resolved absolute path (str).
        # Order is preserved (insertion order matches file order).
        self._data = {}

        # The header label actually seen for the ID column ('gene_id' or one of
        # LEGACY_ID_HEADERS). Useful for error messages.
        self._id_header_seen = None

        self._parse()

        if not skip_check_sanity:
            self.check_sanity()

    def _base_dir(self):
        return Path(self.path).resolve().parent

    @staticmethod
    def _resolve_path(p, base_dir):
        pp = Path(str(p)).expanduser()
        if not pp.is_absolute():
            pp = (base_dir / pp).resolve()
        else:
            pp = pp.resolve()
        return str(pp)

    def _parse(self):
        filesnpaths.is_file_tab_delimited(self.path)

        with open(self.path, 'r') as f:
            header_line = f.readline().rstrip('\n')
            headers = header_line.split('\t')

            if len(headers) != 2:
                raise ConfigError(f"A `structures-txt` file must have exactly two tab-separated columns, but "
                                  f"'{self.path}' has {len(headers)}. Expected headers: "
                                  f"{', '.join(self.REQUIRED_HEADERS)}.")

            id_header, path_header = headers[0], headers[1]

            if id_header == self.REQUIRED_HEADERS[0]:
                self._id_header_seen = id_header
            elif id_header in self.LEGACY_ID_HEADERS:
                self._id_header_seen = id_header
            else:
                raise ConfigError(f"The first column of '{self.path}' must be '{self.REQUIRED_HEADERS[0]}' "
                                  f"(or the legacy alias '{self.LEGACY_ID_HEADERS[0]}'). Got '{id_header}'.")

            if path_header != self.REQUIRED_HEADERS[1]:
                raise ConfigError(f"The second column of '{self.path}' must be '{self.REQUIRED_HEADERS[1]}'. "
                                  f"Got '{path_header}'.")

            base_dir = self._base_dir()
            id_label = self._id_header_seen

            for line_no, raw in enumerate(f, start=2):
                line = raw.rstrip('\n')
                if not line.strip():
                    continue

                fields = line.split('\t')
                if len(fields) != 2:
                    raise ConfigError(f"Line {line_no} of '{self.path}' does not have exactly two tab-separated "
                                      f"fields: {raw!r}.")

                gene_id, raw_path = fields[0].strip(), fields[1].strip()

                if not gene_id:
                    raise ConfigError(f"Line {line_no} of '{self.path}' has an empty {id_label}.")
                if not raw_path:
                    raise ConfigError(f"Line {line_no} of '{self.path}' has an empty path for {id_label} '{gene_id}'.")

                if gene_id in self._data:
                    raise ConfigError(f"{id_label} '{gene_id}' appears more than once in '{self.path}'. Each "
                                      f"{id_label} must map to a single structure file.")

                self._data[gene_id] = self._resolve_path(raw_path, base_dir)

    def check_sanity(self):
        """Sanity checks that don't depend on any external context (no contigs-db, no pan-db)."""

        if not self._data:
            raise ConfigError(f"The structures-txt file '{self.path}' has a header but no data rows.")

        self._check_files_exist()
        self._check_extensions()

    def _check_files_exist(self):
        missing = []
        for gene_id, p in self._data.items():
            if not filesnpaths.is_file_exists(p, dont_raise=True):
                missing.append((gene_id, p))

        if missing:
            preview = ', '.join(f"{g} -> {p}" for g, p in missing[:5])
            raise FilesNPathsError(f"{len(missing)} structure file(s) listed in '{self.path}' do not exist on "
                                   f"disk. The first {min(5, len(missing))}: {preview}.")

    def _check_extensions(self):
        bad = []
        for gene_id, p in self._data.items():
            lower = p.lower()
            if not any(lower.endswith(ext) for ext in self.ALLOWED_EXTENSIONS):
                bad.append((gene_id, p))

        if bad:
            preview = ', '.join(f"{g} -> {Path(p).name}" for g, p in bad[:5])
            raise ConfigError(f"{len(bad)} entry/entries in '{self.path}' point to files with unsupported "
                              f"extensions. Allowed: {', '.join(self.ALLOWED_EXTENSIONS)}. First offenders: "
                              f"{preview}.")

    def gene_ids(self):
        """List of gene_ids in the order they appeared in the file."""
        return list(self._data.keys())

    def paths(self):
        """List of resolved absolute structure paths in the order they appeared in the file."""
        return list(self._data.values())

    def get_path(self, gene_id):
        if gene_id not in self._data:
            raise ConfigError(f"{self._id_header_seen or 'gene_id'} '{gene_id}' not found in '{self.path}'.")
        return self._data[gene_id]

    def as_dict(self):
        """Shallow copy of the gene_id -> path mapping."""
        return dict(self._data)

    def __len__(self):
        return len(self._data)

    def __contains__(self, gene_id):
        return gene_id in self._data
