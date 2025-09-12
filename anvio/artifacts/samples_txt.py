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
        "free":       [],          # only the first column is guaranteed; others are optional
    }

    def __init__(self, args, run=terminal.Run(), expected_format="paired_end", valid_formats=None):
        """
        Parameters
        ----------
        args : argparse.Namespace with args.samples_txt
        expected_format : str in DEFAULT_VALID_FORMATS (incl. "free")
        valid_formats : optional dict[str, list[str]] to extend/override DEFAULT_VALID_FORMATS
        """
        self.samples_txt = getattr(args, "samples_txt", None)
        if not self.samples_txt:
            raise ConfigError("SamplesTxt expects args.samples_txt to point to a samples-txt file.")

        # Ensure tab-delimited and get headers (including first column)
        filesnpaths.is_file_tab_delimited(self.samples_txt)
        columns_found = utils.get_columns_of_TAB_delim_file(self.samples_txt, include_first_column=True)

        # Accept first column 'sample' or 'name'
        first_col = columns_found[0]
        if first_col not in ("sample", "name"):
            raise ConfigError("The first column of any samples-txt must be either `sample` or `name` :/")

        # Allow caller to extend supported modes
        self.valid_formats = dict(self.DEFAULT_VALID_FORMATS)
        if valid_formats:
            self.valid_formats.update(valid_formats)

        if expected_format not in self.valid_formats:
            raise ConfigError(
                f"SamplesTxt: expected_format must be one of {list(self.valid_formats)}, got '{expected_format}'."
            )
        self.expected_format = expected_format

        # Compute expected/potential columns for header-level checks
        minimal_required = [first_col] + (["group"] if "group" in columns_found else [])
        # Mode-required columns (free => none here)
        mode_required = self.valid_formats[self.expected_format]
        expected_columns = [first_col] + mode_required

        # Header sanity (keep legacy behavior: warn on extras, require mode columns)
        possible_columns = set([first_col, "r1", "r2", "lr", "group"])
        if not set(expected_columns).issubset(set(columns_found)):
            # phrase similar to legacy error tone
            raise ConfigError(f"A samples txt file is supposed to have at least the columns {', '.join(expected_columns)}.")

        extra_columns = set(columns_found) - possible_columns
        if extra_columns:
            run.warning(
                "Your samples txt file contains %s: %s compared to what is expected of a `samples-txt` file, "
                "which is absolutely fine. You're reading this message because anvi'o wanted to make sure you "
                "know that it knows that it is the case. Classic anvi'o virtue signaling."
                % (utils.pluralize('extra column', len(extra_columns)), ', '.join(sorted(extra_columns))),
                lc="yellow"
            )

        # Read as dict keyed by first column via existing utils function
        # NOTE: we do not pass expected_fields here because 'free' should not force r1/r2/lr presence.
        raw = utils.get_TAB_delimited_file_as_dictionary(self.samples_txt)

        # Normalize into a new dict keyed by canonical 'sample' (even if file used 'name')
        self.data = self._normalize_and_validate(raw, first_col)

        # Soft warning about unconventional fastq suffixes for short-reads (keeps workflow vibes)
        self._warn_on_unconventional_fastq_suffixes()

    def to_dict(self):
        """Return the full normalized mapping."""
        return self.data

    def as_dataframe(self):
        """Return a pandas.DataFrame resembling the legacy workflow expectation."""
        import pandas as pd
        rows = []
        for s, info in self.data.items():
            row = {
                "sample": s,
                "group": info.get("group"),
                "r1": ",".join(info.get("r1", [])) if info.get("r1") else "",
                "r2": ",".join(info.get("r2", [])) if info.get("r2") else "",
                "lr": ",".join(info.get("lr", [])) if info.get("lr") else "",
            }
            rows.append(row)
        df = pd.DataFrame(rows)
        # Preserve legacy columns order if present
        cols = ["sample", "group", "r1", "r2", "lr"]
        return df[[c for c in cols if c in df.columns]]

    def samples(self):
        return list(self.data.keys())

    def groups(self):
        """Return {group_value: [samples...]} for rows with a non-empty group."""
        g = {}
        for s, info in self.data.items():
            grp = info.get("group")
            if grp:
                g.setdefault(grp, []).append(s)
        return g

    def has_groups(self):
        return any(bool(info.get("group")) for info in self.data.values())

    def short_reads_dict(self, include_single_end=True):
        """
        {sample: {'r1': [...], 'r2': [...]}} for samples with short reads.
        If include_single_end=False, only include rows with both r1 and r2.
        """
        out = {}
        for s, info in self.data.items():
            r1, r2 = info.get("r1", []), info.get("r2", [])
            if r1 and r2:
                out[s] = {"r1": list(r1), "r2": list(r2)}
            elif include_single_end and r1 and not r2:
                out[s] = {"r1": list(r1), "r2": []}
        return out

    def long_reads_dict(self):
        """{sample: {'lr': [...]}} for samples that provide long reads."""
        return {s: {"lr": list(info.get("lr", []))}
                for s, info in self.data.items() if info.get("lr")}

    @staticmethod
    def _split_paths(cell):
        if cell is None:
            return []
        # utils.get_TAB_delimited_file_as_dictionary returns strings
        return [p.strip() for p in str(cell).split(",") if p.strip()]

    def _normalize_and_validate(self, raw, first_col):
        """
        - Normalize keys to 'sample'.
        - Normalize columns to {group, r1[], r2[], lr[]}.
        - Validate per expected_format, check file existence, matching counts, identical r1/r2.
        """
        data = {}
        missing_files = set()
        identical_r1_r2 = set()

        # Build normalized map
        for key, row in raw.items():
            # row is a dict of column->value (strings)
            # Canonical sample id from first column (sample/name)
            sample = str(row.get(first_col, key)).strip()
            if not sample:
                raise ConfigError("Encountered an empty sample name.")
            # Validate sample id via existing util
            utils.check_sample_id(sample)
            if sample in data:
                raise ConfigError(f"Names of samples in your samples_txt file must be unique. Found duplicate: {sample}")

            r1 = self._split_paths(row.get("r1"))
            r2 = self._split_paths(row.get("r2"))
            lr = self._split_paths(row.get("lr"))
            group = row.get("group", None)
            group = (str(group).strip() if group is not None else None) or None

            info = {"group": group, "r1": r1, "r2": r2, "lr": lr}
            data[sample] = info

        # Mode-specific validation
        for sample, info in data.items():
            self._validate_row_by_mode(sample, info, self.expected_format)

        # file presence + identical r1/r2 checks (extended to lr)
        for sample, info in data.items():
            # Paired consistency checks already done; now existence:
            for path in info.get("r1", []) + info.get("r2", []) + info.get("lr", []):
                if not os.path.exists(path):
                    missing_files.add(sample)

            # Identical r1/r2 paths (pairwise)
            r1, r2 = info.get("r1", []), info.get("r2", [])
            if r1 and r2:
                if len(r1) != len(r2):
                    # Should have been caught earlier, but keep a helpful error here too.
                    raise ConfigError(
                        f"Uh oh. The sample {sample} has a different number of R1 ({len(r1)}) and R2 ({len(r2)}) paths."
                        " Anvi'o expects these to be the same, so please fix this in your samples-txt file."
                    )
                for i in range(len(r1)):
                    if r1[i] == r2[i]:
                        identical_r1_r2.add(sample)

        if missing_files:
            raise ConfigError(
                "Bad news. Your samples txt contains %s (%s) with missing files (by which we mean that the "
                "r1/r2/lr paths are there, but the files they point to are not)." %
                (utils.pluralize('sample', len(missing_files)), ', '.join(sorted(missing_files)))
            )

        if identical_r1_r2:
            raise ConfigError(
                "Interesting. Your samples txt contains %s (%s) where r1 and r2 file paths are identical. Not OK." %
                (utils.pluralize('sample', len(identical_r1_r2)), ', '.join(sorted(identical_r1_r2)))
            )

        return data

    def _validate_row_by_mode(self, sample, info, mode):
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
            # must have paired AND lr
            if not r1 or not r2:
                raise ConfigError(f"[{sample}] Hybrid expected: require both 'r1' and 'r2'.")
            if len(r1) != len(r2):
                raise ConfigError(f"[{sample}] Hybrid expected: number of r1 files ({len(r1)}) must match r2 ({len(r2)}).")
            if not lr:
                raise ConfigError(f"[{sample}] Hybrid expected: require 'lr' with at least one path.")

        elif mode == "free":
            # must provide at least one modality
            if not (r1 or r2 or lr):
                raise ConfigError(f"[{sample}] Free mode: provide at least one of 'lr' or 'r1' (with optional 'r2').")
            # if r2 is present, r1 must be present and counts must match
            if r2 and not r1:
                raise ConfigError(f"[{sample}] Free mode: 'r2' provided without 'r1'.")
            if r1 and r2 and len(r1) != len(r2):
                raise ConfigError(
                    f"[{sample}] Free mode: number of r1 files ({len(r1)}) must match r2 ({len(r2)})."
                )
        else:
            raise ConfigError(f"Unknown validation mode '{mode}'")

    def _warn_on_unconventional_fastq_suffixes(self):
        # Check ALL paths (short-reads and long-reads)
        all_paths = []
        for info in self.data.values():
            all_paths.extend(info.get("r1", []))
            all_paths.extend(info.get("r2", []))
            all_paths.extend(info.get("lr", []))

        allowed_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
        bad = [p for p in all_paths if p and not p.endswith(allowed_suffixes)]

        if bad:
            run.warning(
                "We noticed some of your sequence files in '%s' do not end with one of "
                "'.fastq', '.fastq.gz', '.fq', or '.fq.gz'. That's okay, but anvi'o decided it "
                "should warn you. Here are the first 5 such files that have unconventional "
                "file extensions: %s."
                % (self.samples_txt, ", ".join(bad[:5]))
            )

