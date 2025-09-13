# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

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

    def __init__(self, path, expected_format="free", valid_formats=None, skip_check_sanity=False):
        self.artifact_path = path
        self.expected_format = expected_format
        self.valid_formats = dict(self.DEFAULT_VALID_FORMATS)
        if valid_formats:
            self.valid_formats.update(valid_formats)

        # Light, non-failing prep. All raising happens in check_sanity().
        self._raw_text = None
        self._columns_found = None
        self._first_col = None
        self._raw_rows = None  # as returned by get_TAB_delimited_file_as_dictionary
        self._data = None      # normalized dict we expose

        self._load_raw_text()                 # populate _raw_text
        self._inspect_headers()      # populate _columns_found / _first_col
        self._parse_rows()           # populate _raw_rows and normalized _data (no strict validation yet)

        if not skip_check_sanity:
            self.check_sanity()

    def check_sanity(self):
        """Perform all sanity checks and warnings."""
        # Validate mode key
        if self.expected_format not in self.valid_formats:
            raise ConfigError(
                f"SamplesTxt speaking: expected_format must be one of {list(self.valid_formats)}, got '{self.expected_format}'."
            )

        # File is tab-delimited (raises on failure)
        filesnpaths.is_file_tab_delimited(self.artifact_path)

        # Header checks (first col, required columns for the chosen mode)
        if self._first_col not in ("sample", "name"):
            raise ConfigError("The first column of any samples-txt must be either `sample` or `name` :/")

        expected_columns = [self._first_col] + self.valid_formats[self.expected_format]
        if not set(expected_columns).issubset(set(self._columns_found)):
            raise ConfigError(
                f"A samples txt file is supposed to have at least the columns {', '.join(expected_columns)}."
            )

        # Warn about extras columns
        possible_columns = set([self._first_col, "r1", "r2", "lr", "group"])
        extra_columns = set(self._columns_found) - possible_columns
        if extra_columns:
            run.warning(
                "Your samples txt file contains %s: %s compared to what is expected of a `samples-txt` file, "
                "which is absolutely fine. You're reading this message because anvi'o wanted to make sure you "
                "know that it knows that it is the case. Classic anvi'o virtue signaling."
                % (utils.pluralize('extra column', len(extra_columns)), ', '.join(sorted(extra_columns))),
                lc="yellow"
            )

        # Per-row validations, file existence, identical pairs
        self._validate_rows_by_mode()
        self._check_file_existence_and_identical_pairs()
        self._warn_on_unconventional_fastq_suffixes()

    def as_raw(self):
        """Return the raw TSV text of the samples-txt file."""
        return self._raw_text

    def as_dict(self, include_extras=True):
        """
        Return {sample: {'group','r1','r2','lr', [extras...]}}.
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

    def _load_raw_text(self):
        with open(self.artifact_path, "r", encoding="utf-8") as f:
            self._raw_text = f.read()

    def _inspect_headers(self):
        self._columns_found = utils.get_columns_of_TAB_delim_file(self.artifact_path, include_first_column=True)
        self._first_col = self._columns_found[0] if self._columns_found else None

        # Columns that are not part of the canonical schema (keep header order)
        self._extra_columns = [
            c for c in (self._columns_found or [])
            if c not in {self._first_col, "group", "r1", "r2", "lr"}
        ]

    def _parse_rows(self):
        # Do not pass expected_fields here; 'free' and partial headers should still parse.
        self._raw_rows = utils.get_TAB_delimited_file_as_dictionary(self.artifact_path)
        self._data = self._normalize_rows(self._raw_rows, self._first_col)

    @staticmethod
    def _split_paths(cell):
        if cell is None:
            return []
        return [p.strip() for p in str(cell).split(",") if p.strip()]

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
                    raise ConfigError(
                        f"Uh oh. The sample {sample} has a different number of R1 ({len(r1)}) and R2 ({len(r2)}) paths. "
                        f"Anvi'o expects these to be the same, so please fix this in your samples-txt file."
                    )
                for i in range(len(r1)):
                    if r1[i] == r2[i]:
                        identical.add(sample)

        if missing:
            raise ConfigError(
                "Bad news. Your samples txt contains %s (%s) with missing files (by which we mean that the "
                "r1/r2/lr paths are there, but the files they point to are not)."
                % (utils.pluralize('sample', len(missing)), ", ".join(sorted(missing)))
            )

        if identical:
            raise ConfigError(
                "Interesting. Your samples txt contains %s (%s) where r1 and r2 file paths are identical. Not OK."
                % (utils.pluralize('sample', len(identical)), ", ".join(sorted(identical)))
            )

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
            run.warning(
                "We noticed some of your sequence files in '%s' do not end with one of "
                "'.fastq', '.fastq.gz', '.fq', or '.fq.gz'. That's okay, but anvi'o decided it "
                "should warn you. Here are the first 5 such files that have unconventional "
                "file extensions: %s."
                % (self.artifact_path, ", ".join(bad[:5]))
            )

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
        """
        Return the normalized dict for one sample:
        {'group': ..., 'r1': [...], 'r2': [...], 'lr': [...]}
        """
        return self._data[sample]

