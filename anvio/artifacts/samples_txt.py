# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from pathlib import Path
from anvio.errors import ConfigError


class SamplesTxt:
    """
    Unified handler for anvi'o samples-txt files.

    Modes
    -----
    - paired_end : each row must have r1 & r2 (counts must match). Disallows lr.
    - single_end : each row must have r1 only (no r2, no lr).
    - long_read  : each row must have lr only (no r1/r2).
    - hybrid     : each row must have r1 & r2 (matched) *and* lr.
    - free       : per-row flexibility; allows single-end (r1 only), paired (r1+r2),
                   long-read (lr), or hybrid (r1+r2+lr). Minimal consistency checks.

    Internals
    ---------
    self.data = {
        <sample_name>: {
            'group': <str or None>,
            'r1': [<paths>] or [],
            'r2': [<paths>] or [],
            'lr': [<paths>] or [],
        },
        ...
    }
    """

    DEFAULT_VALID_FORMATS = {
        "paired_end": ["r1", "r2"],
        "single_end": ["r1"],
        "long_read":  ["lr"],
        "hybrid":     ["r1", "r2", "lr"],
        "free":       [],  # only first column + optional others
    }

    def __init__(self, path_or_data, expected_format="free", valid_formats=None, skip_check_sanity=False, run=terminal.Run()):
        self.run = run

        self.expected_format = expected_format
        self.valid_formats = dict(self.DEFAULT_VALID_FORMATS)
        if valid_formats:
            self.valid_formats.update(valid_formats)

        # Detect source type
        if isinstance(path_or_data, (str, Path)):
            self._source = "path"
            self.artifact_path = str(path_or_data)
        else:
            self._source = "dict"
            self.artifact_path = None

        # Placeholders used later
        self._raw_text = None
        self._columns_found = None
        self._first_col = None
        self._raw_rows = None   # mapping-like, sample_key -> row dict
        self._data = None       # normalized

        # Load/parse by source
        if self._source == "path":
            self._load_raw_text()
            self._inspect_headers_from_file()
            self._parse_rows_from_file()
        else:
            self._inspect_headers_from_data(path_or_data)
            self._parse_rows_from_data(path_or_data)

        if not skip_check_sanity:
            self.check_sanity()

    def check_sanity(self):
        """Perform all sanity checks and warnings."""
        # Validate mode key
        if self.expected_format not in self.valid_formats:
            raise ConfigError(f"SamplesTxt speaking: expected_format must be one of {list(self.valid_formats)}, got '{self.expected_format}'.")

        # File checks only when source is a path:
        if self._source == "path":
            filesnpaths.is_file_tab_delimited(self.artifact_path)

        # Header checks (first col, required columns for the chosen mode)
        if self._first_col not in ("sample", "name"):
            raise ConfigError("The first column of any samples-txt must be either `sample` or `name` :/")

        expected_columns = [self._first_col] + self.valid_formats[self.expected_format]
        if not set(expected_columns).issubset(set(self._columns_found)):
            raise ConfigError(f"A samples txt file is supposed to have at least the columns {', '.join(expected_columns)}.")

        # Warn about extras columns
        possible_columns = set([self._first_col, "r1", "r2", "lr", "group"])
        extra_columns = set(self._columns_found) - possible_columns
        if extra_columns:
            self.run.warning(f"Your samples txt file contains {utils.pluralize('extra column', len(extra_columns))}: "
                             f"{', '.join(sorted(extra_columns))} compared to what is expected of a `samples-txt` file, "
                             f"which is absolutely fine. You're reading this message because anvi'o wanted to make sure you "
                             f"know that it knows that it is the case. Classic anvi'o virtue signaling.", lc="yellow")

        # Per-row validations
        self._validate_rows_by_mode()
        self._check_file_existence_and_identical_pairs()
        self._warn_on_unconventional_fastq_suffixes()

    def as_raw(self):
        """If constructed from a file, return raw file text; otherwise a synthesized TSV."""
        if self._source == "path":
            return self._raw_text

        # synthesize TSV in canonical order: first_col, group, r1, r2, lr, extras...
        cols = [self._first_col, "group", "r1", "r2", "lr"] + list(self._extra_columns)
        # keep only those that exist anywhere
        present = [c for c in cols if any(c in row for row in self._raw_rows.values()) or c == self._first_col]
        lines = ["\t".join(present)]
        for key, row in self._raw_rows.items():
            vals = []
            for c in present:
                if c == self._first_col:
                    vals.append(str(row.get(self._first_col, key)))
                elif c in {"r1","r2","lr"}:
                    v = row.get(c, "")
                    # leave as originally provided (string or list → comma-join for display)
                    if isinstance(v, (list, tuple)):
                        vals.append(",".join(str(x) for x in v))
                    else:
                        vals.append(str(v) if v is not None else "")
                else:
                    vals.append("" if row.get(c) in (None,) else str(row.get(c)))
            lines.append("\t".join(vals))
        return "\n".join(lines)

    def as_dict(self, include_extras=True):
        """Return {sample: {'group','r1','r2','lr', [extras...]}}.

        Set include_extras=False to get only the canonical keys.
        """
        if include_extras:
            return self._data
        # only with classic keys
        return {
            s: {k: v for k, v in info.items() if k in {"group", "r1", "r2", "lr"}}
            for s, info in self._data.items()
        }

    def as_df(self, include_extras=True):
        """DataFrame view for downstream code expecting a df."""
        import pandas as pd
        rows = []
        for s, info in self._data.items():
            row = {
                "sample": s,
                "group": info.get("group"),
                "r1": ",".join(info.get("r1", [])) if info.get("r1") else "",
                "r2": ",".join(info.get("r2", [])) if info.get("r2") else "",
                "lr": ",".join(info.get("lr", [])) if info.get("lr") else "",
            }
            if include_extras:
                for col in self._extra_columns:
                    row[col] = info.get(col, "")
            rows.append(row)

        df = pd.DataFrame(rows)
        base = ["sample", "group", "r1", "r2", "lr"]
        if include_extras:
            cols = base + self._extra_columns
        else:
            cols = base
        return df[[c for c in cols if c in df.columns]]

    def extra_columns(self):
        """Return a list of user-supplied extra column names in file order."""
        return list(self._extra_columns)

    # Optional helper for CLI code that still passes args:
    @classmethod
    def from_args(cls, args, **kwargs):
        return cls(getattr(args, "samples_txt", None), **kwargs)

    @classmethod
    def from_dict(cls, data, **kwargs):
        """Build a SamplesTxt from an in-memory dict/list (no file)."""
        return cls(data, **kwargs)

    def _load_raw_text(self):
        with open(self.artifact_path, "r", encoding="utf-8") as f:
            self._raw_text = f.read()

    def _inspect_headers_from_file(self):
        self._columns_found = utils.get_columns_of_TAB_delim_file(self.artifact_path, include_first_column=True)
        self._first_col = self._columns_found[0] if self._columns_found else None

        self._extra_columns = [
            c for c in (self._columns_found or [])
            if c not in {self._first_col, "group", "r1", "r2", "lr"}
        ]

    def _parse_rows_from_file(self):
        # Do not pass expected_fields here; 'free' and partial headers should still parse.
        self._raw_rows = utils.get_TAB_delimited_file_as_dictionary(self.artifact_path)
        self._data = self._normalize_rows(self._raw_rows, self._first_col)

    def _inspect_headers_from_data(self, data):
        """
        Accept either:
          A) dict[str, dict]: {sample -> {r1/r2/lr/group/EXTRA...}}
          B) list[dict]:      [{'sample':..., 'r1':..., ...}, ...]
        Compute synthetic header list: first column, canonical fields, then extras (sorted).
        """
        if isinstance(data, dict):
            # dict-of-dicts → first column will be 'sample'
            keys = set().union(*[set((v or {}).keys()) for v in data.values()]) if data else set()
            self._first_col = "sample"
            cols = ["sample"]
            # Prefer canonical order
            for c in ["group", "r1", "r2", "lr"]:
                if c in keys:
                    cols.append(c)
            extras = sorted([c for c in keys if c not in {"group","r1","r2","lr"}])
            cols.extend(extras)
            self._columns_found = cols
            self._extra_columns = extras
        elif isinstance(data, list):
            keys = set().union(*[set((row or {}).keys()) for row in data]) if data else set()
            self._first_col = "sample" if "sample" in keys else ("name" if "name" in keys else "sample")
            cols = [self._first_col]
            for c in ["group", "r1", "r2", "lr"]:
                if c in keys:
                    cols.append(c)
            extras = sorted([c for c in keys if c not in {self._first_col,"group","r1","r2","lr"}])
            cols.extend(extras)
            self._columns_found = cols
            self._extra_columns = extras
        else:
            raise ConfigError(f"Unsupported data type for SamplesTxt: {type(data)}")

    def _parse_rows_from_data(self, data):
        """Coerce the in-memory object to the same _raw_rows shape used for file input."""
        self._raw_rows = self._coerce_data_to_raw_rows(data, self._first_col)
        self._data = self._normalize_rows(self._raw_rows, self._first_col)

    def _coerce_data_to_raw_rows(self, data, first_col):
        """
        Return a mapping: sample_key -> row dict (row includes first_col key).
        Accepts dict[str, dict] and list[dict].
        Values for r1/r2/lr may be str (comma-separated) or list; we pass them through;
        splitting is handled later by _normalize_rows.
        """
        rows = {}
        if isinstance(data, dict):
            for sample, info in data.items():
                row = dict(info or {})
                # ensure first_col exists
                row[first_col] = sample if first_col == "sample" else row.get(first_col, sample)
                rows[sample] = row
        else:  # list[dict]
            for row in data:
                if first_col not in row and "sample" in row:
                    # promote 'sample' into first_col if header chose 'name'
                    row = dict(row)
                    row[first_col] = row["sample"]
                key = str(row.get(first_col, "")).strip()
                if not key:
                    raise ConfigError("Encountered a row without a sample/name.")
                rows[key] = dict(row)
        return rows

    @staticmethod
    def _split_paths(value):
        """
        Accept str (optionally comma-separated), list/tuple/set of strings,
        or nested mixtures thereof; return a flat list of clean path strings.
        """
        if not value:
            return []

        # string: split on commas
        if isinstance(value, str):
            return [p.strip() for p in value.split(",") if p.strip()]

        # iterables: flatten; also split comma-containing elements
        if isinstance(value, (list, tuple, set)):
            out = []
            for v in value:
                if not v:
                    continue
                if isinstance(v, (list, tuple, set)):
                    out.extend(SamplesTxt._split_paths(v))
                else:
                    s = str(v).strip()
                    if not s:
                        continue
                    if "," in s:
                        out.extend([p.strip() for p in s.split(",") if p.strip()])
                    else:
                        out.append(s)
            return out

        # fallback: treat as a single scalar path
        s = str(value).strip()
        return [s] if s else []

    def _normalize_rows(self, raw, first_col):
        data = {}
        for key, row in raw.items():
            sample = str(row.get(first_col, key)).strip()
            if not sample:
                # Keep this non-fatal here; will be enforced in check_sanity()
                sample = key if key else ""
            if not sample:
                # If truly empty, raise now—there's no way to key the dict.
                raise ConfigError("Encountered an empty sample name.")

            # Validate sample id
            utils.check_sample_id(sample)

            if sample in data:
                raise ConfigError(f"Names of samples in your samples_txt file must be unique. Found duplicate: {sample}")

            info = {
                "group": (str(row.get("group")).strip() if row.get("group") not in (None, "") else None),
                "r1": self._split_paths(row.get("r1")),
                "r2": self._split_paths(row.get("r2")),
                "lr": self._split_paths(row.get("lr")),
            }

            # Flatten user extras (strings as-is; keep empty as "")
            for col in self._extra_columns:
                val = row.get(col)
                info[col] = "" if val in (None,) else str(val)

            data[sample] = info
        return data

    def _validate_rows_by_mode(self):
        mode = self.expected_format
        for sample, info in self._data.items():
            r1, r2, lr = info["r1"], info["r2"], info["lr"]
            if mode == "paired_end":
                if not r1 or not r2:
                    raise ConfigError(f"[{sample}] Paired-end expected: require both 'r1' and 'r2'.")
                if len(r1) != len(r2):
                    raise ConfigError(f"[{sample}] Paired-end expected: number of r1 files ({len(r1)}) must match r2 ({len(r2)}).")
                if lr:
                    raise ConfigError(f"[{sample}] Paired-end mode: unexpected long-read ('lr') paths provided.")
            elif mode == "single_end":
                if not r1:
                    raise ConfigError(f"[{sample}] Single-end expected: require 'r1' with at least one path.")
                if r2:
                    raise ConfigError(f"[{sample}] Single-end mode: 'r2' must be empty/absent.")
                if lr:
                    raise ConfigError(f"[{sample}] Single-end mode: 'lr' must be empty/absent.")
            elif mode == "long_read":
                if not lr:
                    raise ConfigError(f"[{sample}] Long-read expected: require 'lr' with at least one path.")
                if r1 or r2:
                    raise ConfigError(f"[{sample}] Long-read mode: 'r1'/'r2' must be empty/absent.")
            elif mode == "hybrid":
                if not r1 or not r2:
                    raise ConfigError(f"[{sample}] Hybrid expected: require both 'r1' and 'r2'.")
                if len(r1) != len(r2):
                    raise ConfigError(f"[{sample}] Hybrid expected: number of r1 files ({len(r1)}) must match r2 ({len(r2)}).")
                if not lr:
                    raise ConfigError(f"[{sample}] Hybrid expected: require 'lr' with at least one path.")
            elif mode == "free":
                if not (r1 or r2 or lr):
                    raise ConfigError(f"[{sample}] Free mode: provide at least one of 'lr' or 'r1' (with optional 'r2').")
                if r2 and not r1:
                    raise ConfigError(f"[{sample}] Free mode: 'r2' provided without 'r1'.")
                if r1 and r2 and len(r1) != len(r2):
                    raise ConfigError(f"[{sample}] Free mode: number of r1 files ({len(r1)}) must match r2 ({len(r2)}).")
            else:
                raise ConfigError(f"Unknown validation mode '{mode}'")

    def _check_file_existence_and_identical_pairs(self):
        missing = set()
        identical = set()
        for sample, info in self._data.items():
            for p in info.get("r1", []) + info.get("r2", []) + info.get("lr", []):
                if not os.path.exists(p):
                    missing.add(sample)
            r1, r2 = info.get("r1", []), info.get("r2", [])
            if r1 and r2:
                if len(r1) != len(r2):
                    # Redundant with _validate_rows_by_mode, but keep the legacy-friendly message around:
                    raise ConfigError(f"Uh oh. The sample {sample} has a different number of R1 ({len(r1)}) and R2 ({len(r2)}) paths. "
                                      f"Anvi'o expects these to be the same, so please fix this in your samples-txt file.")
                for i in range(len(r1)):
                    if r1[i] == r2[i]:
                        identical.add(sample)

        if missing:
            raise ConfigError(f"Bad news. Your samples txt contains {utils.pluralize('sample', len(missing))} "
                              f"({', '.join(sorted(missing))}) with missing files (by which we mean that the "
                              f"r1/r2/lr paths are there, but the files they point to are not).")

        if identical:
            raise ConfigError(f"Interesting. Your samples txt contains {utils.pluralize('sample', len(identical))} "
                              f"({', '.join(sorted(identical))}) where r1 and r2 file paths are identical. Not OK.")

    def _warn_on_unconventional_fastq_suffixes(self):
        # Check ALL paths (SR and LR). Accept .fastq, .fastq.gz, .fq, .fq.gz
        all_paths = []
        for info in self._data.values():
            all_paths.extend(info.get("r1", []))
            all_paths.extend(info.get("r2", []))
            all_paths.extend(info.get("lr", []))

        allowed = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
        bad = [p for p in all_paths if p and not p.endswith(allowed)]
        if bad:
            self.run.warning(f"We noticed some of your sequence files in '{self.artifact_path}' do not end with one "
                             f"of the expected extensions (which include '.fastq', '.fastq.gz', '.fq', or '.fq.gz'). "
                             f"That's okay, but anvi'o decided that it still should warn you. Here are the first 5 "
                             f"such files that have unconventional extensions: {', '.join(bad[:5])}.")

    def samples(self):
        """List of sample names."""
        return list(self._data.keys())

    def has_groups(self):
        """True if at least one row has a non-empty 'group' value."""
        return any(bool(info.get("group")) for info in self._data.values())

    def groups(self):
        """{group_value: [samples...]} for rows with 'group' set."""
        g = {}
        for s, info in self._data.items():
            grp = info.get("group")
            if grp:
                g.setdefault(grp, []).append(s)
        return g

    def group_names(self):
        """List of group names present (unique)."""
        return list(self.groups().keys())

    def group_sizes(self):
        """{group_value: size}."""
        return {g: len(members) for g, members in self.groups().items()}

    def short_reads_dict(self, include_single_end=True):
        """
        {sample: {'r1': [...], 'r2': [...]}} for samples that provide short reads.
        If include_single_end=False, only include rows where both r1 and r2 exist.
        """
        out = {}
        for s, info in self._data.items():
            r1, r2 = info.get("r1", []), info.get("r2", [])
            if r1 and r2:
                out[s] = {"r1": list(r1), "r2": list(r2)}
            elif include_single_end and r1 and not r2:
                out[s] = {"r1": list(r1), "r2": []}
        return out

    def long_reads_dict(self):
        """{sample: {'lr': [...]}} for samples that provide long reads."""
        return {s: {"lr": list(info.get("lr", []))}
                for s, info in self._data.items() if info.get("lr")}

    def get_sample(self, sample):
        """Return the normalized dict for one sample:

        {'group': ..., 'r1': [...], 'r2': [...], 'lr': [...]}
        """
        return self._data[sample]

    @staticmethod
    def _has_sr(info: dict) -> bool:
        """A row has short-reads if r1 is provided (r2 may be empty in 'free' mode)."""
        return bool(info.get("r1"))

    @staticmethod
    def _has_lr(info: dict) -> bool:
        """A row has long-reads if lr is provided."""
        return bool(info.get("lr"))

    def has_any_sr(self) -> bool:
        """Return True if at least one sample row has short-reads."""
        return any(self._has_sr(info) for info in self._data.values())

    def has_any_lr(self) -> bool:
        """Return True if at least one sample row has long-reads."""
        return any(self._has_lr(info) for info in self._data.values())

    def sample_types(self) -> dict:
        """
        Return {base_sample: {'has_sr': bool, 'has_lr': bool}} for quick inspection.
        """
        out = {}
        for s, info in self._data.items():
            out[s] = {"has_sr": self._has_sr(info), "has_lr": self._has_lr(info)}
        return out

    def iter_readsets(self, read_type_suffix: str = "auto"):
        """
        Yield canonical 'readset' records for downstream workflows.

        A readset is the concrete unit that will be mapped/profiled.
        Fields:
          - id           : identifier used in file names / Snakemake wildcards.
                           Suffix rules (default 'auto'):
                             * If a sample has BOTH SR and LR → emit two readsets:
                               '{sample}_SR' and '{sample}_LR'.
                             * Otherwise (only one type) → use the unsuffixed '{sample}'.
                           (Future-proof: read_type_suffix='force' always adds _SR/_LR.)
          - base_sample  : the unsuffixed sample name from samples.txt
          - type         : 'SR' or 'LR'
          - mixed_type   : True if the base sample has both SR and LR (useful for logic)
          - group        : the raw group value (may be None if not provided)
          - reads        : dict of lists; for SR → {'r1': [...], 'r2': [...]}, for LR → {'lr': [...]}

        Parameters
        ----------
         : {'auto', 'force'}
            'auto'  (default): only suffix when a sample has both SR and LR.
            'force': always suffix by type (_SR/_LR), even if only one type exists.

        Yields
        ------
        dict
            Readset record as described above.
        """
        if read_type_suffix not in {"auto", "force"}:
            raise ConfigError("read_type_suffix must be 'auto' or 'force'.")

        for sample, info in self._data.items():
            has_sr = self._has_sr(info)
            has_lr = self._has_lr(info)
            mixed = has_sr and has_lr

            # Nothing to emit (shouldn't happen due to sanity checks in 'free' mode).
            if not (has_sr or has_lr):
                continue

            # Helper to build one readset dict
            def _mk(id_, typ, reads_dict):
                return {"id": id_,
                        "base_sample": sample,
                        "type": typ,                # 'SR' or 'LR'
                        "mixed_type": mixed,        # True if this base sample has both SR and LR
                        "group": info.get("group"),
                        "reads": reads_dict, }      # only keys relevant to the type

            # Decide ids under chosen read-type suffix
            if mixed:
                yield _mk(f"{sample}_SR", "SR", {"r1": info.get("r1", []), "r2": info.get("r2", [])})
                yield _mk(f"{sample}_LR", "LR", {"lr":  info.get("lr",  [])})
            else:
                # Only SR
                if has_sr:
                    rs_id = f"{sample}_SR" if read_type_suffix == "force" else sample
                    yield _mk(rs_id, "SR", {"r1": info.get("r1", []), "r2": info.get("r2", [])})
                # Only LR
                if has_lr:
                    rs_id = f"{sample}_LR" if read_type_suffix == "force" else sample
                    yield _mk(rs_id, "LR", {"lr": info.get("lr", [])})

    def get_readsets(self, read_type_suffix: str = "auto"):
        """
        Convenience wrapper that returns list(self.iter_readsets(...)).
        """
        return list(self.iter_readsets(read_type_suffix=read_type_suffix))

