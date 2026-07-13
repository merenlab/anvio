"""
    QCModule — a mixin providing quality-control steps for short-read and long-read workflows.
"""

import os
import gzip
import anvio
import anvio.utils as u

from anvio.workflows import WorkflowSuperClass
from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


class QCModule(WorkflowSuperClass):
    """Mixin providing quality-control steps reusable across workflows.

    Provides:
      - illumina-utils SR filtering (iu-gen-configs + iu-filter-quality-minoche)
      - QC report aggregation (gen_qc_report)
      - Optional gzip of QC'd SR reads
      - Optional FastQC on SR reads
      - Optional Filtlong filtering of LR reads
      - Optional NanoPlot quality assessment on LR reads (on the raw and/or filtered reads)
      - Optional FastQC quality assessment on SR reads (on the raw and/or filtered reads)
      - Optional MultiQC aggregation of all QC outputs

    Subclasses must set before any QC methods are called:
      - self.readsets: list[dict] from SamplesTxt.iter_readsets()
      - self.dirs_dict: includes "QC_DIR"
      - self.run_qc: bool (SR QC enabled)
      - self.run_filtlong: bool
      - self.run_nanoplot: bool (LR QC enabled)
      - self.run_multiqc: bool
    """

    def __init__(self):
        # Register the QC rule names this mixin provides. Their accepted parameters, types,
        # and default values are declared in the consuming workflow's params.json (e.g.
        # anvio/workflows/metagenomics/params.json), which is the single source of truth —
        # WorkflowSuperClass.init() populates rule_acceptable_params_dict and the default
        # config from that schema, so we intentionally do not duplicate them here.
        self.rules.extend([
            'iu_gen_configs',
            'iu_filter_quality_minoche',
            'gen_qc_report',
            'gzip_fastqs',
            'fastqc_sr',
            'filtlong',
            'nanoplot',
            'multiqc',
        ])

    def filtered_lr_path(self, readset):
        """Single source of truth for the Filtlong output path of a long-read readset.

        Used by get_fastq(), get_nanoplot_input_files(), get_qc_target_files(), and the filtlong
        rule (lr.smk) so the filtered-reads filename can never drift between the rule that writes it
        and the code that consumes it.
        """
        return os.path.join(self.dirs_dict["QC_DIR"], f"{readset}-FILTERED_LR.fastq.gz")

    def get_fastq(self, readset, pre_ref_removal=False):
        """Return FASTQ paths for a readset.

        For LR: returns filtlong output (if run_filtlong), otherwise raw LR.
        For SR: delegates to _resolve_sr_path() — override there for ref-removal logic.
        """
        rs = self.readsets_by_id.get(readset)
        if rs['type'] == 'LR':
            if getattr(self, 'run_filtlong', False):
                return {'lr': [self.filtered_lr_path(readset)]}
            return {'lr': self.get_lr_files_for_readset(readset)}
        elif rs['type'] == 'SR':
            return self._resolve_sr_path(readset, pre_ref_removal=pre_ref_removal)

    def _resolve_sr_path(self, readset, pre_ref_removal=False):
        """Return SR FASTQ paths; uses QC'd outputs when run_qc is set."""
        if getattr(self, 'run_qc', False):
            zipped = self.get_param_value_from_config(['gzip_fastqs', 'run']) == True
            ext = ".fastq.gz" if zipped else ".fastq"
            r1 = os.path.join(self.dirs_dict["QC_DIR"], f"{readset}-QUALITY_PASSED_R1{ext}")
            r2 = os.path.join(self.dirs_dict["QC_DIR"], f"{readset}-QUALITY_PASSED_R2{ext}")
            return {'r1': [r1], 'r2': [r2]}
        return self.get_sr_files_for_readset(readset)

    def _conda_options_set(self, tool):
        """Return the conda-source options set for a tool, as a list of their config-key names.

        The three (mutually exclusive) ways a rule can be given its program via conda:
        'conda_yaml' (a user YAML path), 'conda_env' (an existing env name), and
        'use_anvio_conda_yaml' (the env file anvi'o ships). Single source of truth for
        _tool_provided_by_conda(), MetagenomicsWorkflow.ensure_tool_in_path_or_conda(), and the
        conda mutual-exclusivity check in MetagenomicsWorkflow.init().
        """
        y = self.get_param_value_from_config([tool, 'conda_yaml'])
        n = self.get_param_value_from_config([tool, 'conda_env'])
        a = self.get_param_value_from_config([tool, 'use_anvio_conda_yaml']) == True
        options = []
        if y and y.strip():
            options.append('conda_yaml')
        if n and n.strip():
            options.append('conda_env')
        if a:
            options.append('use_anvio_conda_yaml')
        return options

    def _tool_provided_by_conda(self, tool):
        """Whether a tool is supplied via a conda env rather than $PATH.

        True when any conda-source option is set (see _conda_options_set). When it is, we must not
        check for the executable (or its Python modules) on the current $PATH / interpreter —
        Snakemake will run the rule inside the configured environment.
        """
        return bool(self._conda_options_set(tool))

    def check_qc_program_dependencies(self):
        """Raise ConfigError for any program required by an enabled QC tool that is missing.

        Only checks tools that are actually enabled in the config, avoiding the
        false-positive noise from the generic check_workflow_program_dependencies. Tools provided
        via conda (conda_yaml/conda_env) are not checked on the current $PATH/interpreter.
        """
        missing = []

        if self.get_param_value_from_config(['iu_filter_quality_minoche', 'run']) == True:
            for prog in ['iu-gen-configs', 'iu-filter-quality-minoche']:
                if not u.is_program_exists(prog, dont_raise=True):
                    missing.append((prog, 'iu_filter_quality_minoche'))

        if self.get_param_value_from_config(['filtlong', 'run']) == True:
            if not self._tool_provided_by_conda('filtlong') and not u.is_program_exists('filtlong', dont_raise=True):
                missing.append(('filtlong', 'filtlong'))

        if self.get_param_value_from_config(['nanoplot', 'run']) == True:
            if not self._tool_provided_by_conda('nanoplot') and not u.is_program_exists('NanoPlot', dont_raise=True):
                missing.append(('NanoPlot', 'nanoplot'))

        if self.get_param_value_from_config(['fastqc_sr', 'run']) == True:
            if not self._tool_provided_by_conda('fastqc_sr') and not u.is_program_exists('fastqc', dont_raise=True):
                missing.append(('fastqc', 'fastqc_sr'))

        if self.get_param_value_from_config(['multiqc', 'run']) == True:
            if not self._tool_provided_by_conda('multiqc') and not u.is_program_exists('multiqc', dont_raise=True):
                missing.append(('multiqc', 'multiqc'))

        if missing:
            details = '\n  '.join(f"{prog}  (required by: {rule})" for prog, rule in missing)
            raise ConfigError(
                f"Anvi'o found {len(missing)} missing program(s) required by enabled QC steps. "
                f"Please install them before running the workflow:\n\n  {details}"
            )

    def check_lr_readset_no_duplicate_names(self, readset):
        """Raise ConfigError if a long-read readset's FASTQ(s) contain duplicate read names.

        Filtlong aborts mid-run on the first duplicate it finds, giving a cryptic error, so this
        validates a readset up front with a clear, actionable message. It scans the full input
        file(s), which is expensive — so it is invoked as its own Snakemake rule (see the
        check_lr_read_names rule in lr.smk) that gates filtlong, rather than at parse time. Doing
        it at parse time would re-read every long-read file on every dry run and DAG rebuild.
        """
        problematic = []
        # A readset can have several input files, and filtlong receives all of them together — so a
        # read name is a duplicate if it repeats ANYWHERE in the readset, not just within one file.
        # Track `seen` across every file of the readset (not reset per file).
        seen = set()
        for path in self.get_lr_files_for_readset(readset):
            if not os.path.exists(path):
                continue
            opener = gzip.open if path.endswith('.gz') else open
            with opener(path, 'rt') as f:
                for i, line in enumerate(f):
                    if i % 4 != 0:
                        continue
                    line = line.strip()
                    if not line:
                        # tolerate blank lines (e.g. a trailing newline at end of file)
                        continue
                    tokens = line[1:].split()
                    if not tokens:
                        # a header line with no read name → malformed FASTQ; fail with a clear
                        # message rather than an opaque IndexError from indexing an empty split
                        raise ConfigError(
                            f"Anvi'o ran into a malformed FASTQ record while checking read names in "
                            f"'{path}' (readset '{readset}'): the header on line {i + 1} has no read "
                            f"name. Please make sure this file is valid FASTQ before retrying."
                        )
                    name = tokens[0]
                    if name in seen:
                        problematic.append((path, name))
                        break
                    seen.add(name)

        if problematic:
            details = '\n  '.join(f"{p}: first duplicate: '{n}'" for p, n in problematic)
            raise ConfigError(
                f"Anvi'o found duplicate read names in the long-read input for readset '{readset}'. "
                f"Filtlong will crash mid-run when it encounters them, so anvi'o refuses to continue. "
                f"Fix the input file(s) with 'seqkit rename' before retrying:\n\n"
                f"  seqkit rename <input.fastq.gz> -o <fixed.fastq.gz>\n\n"
                f"Affected file(s):\n  {details}"
            )

    def _qc_stages_for(self, tool):
        """Return the list of QC stages ('raw' and/or 'filtered') enabled for a stats tool.

        A stats tool (nanoplot / fastqc_sr) can be asked to run on the raw input reads
        ('run_on_raw'), on the post-filter reads ('run_on_filtered'), or both. The returned
        list drives both the target-file generation and the per-stage rule wildcards.
        Assumes sanity_check_qc_stage_flags() has already validated the combination.
        """
        stages = []
        if self.get_param_value_from_config([tool, 'run_on_raw']) == True:
            stages.append('raw')
        if self.get_param_value_from_config([tool, 'run_on_filtered']) == True:
            stages.append('filtered')
        return stages

    def get_nanoplot_input_files(self, readset, stage):
        """Return the long-read FASTQ files NanoPlot should run on for a readset and stage.

        stage='raw'      → the readset's original long-read files.
        stage='filtered' → the Filtlong output for the readset (requires filtlong to be enabled;
                           depending on this path also creates the DAG edge that forces filtlong
                           to finish first).
        """
        if stage == 'raw':
            return self.get_lr_files_for_readset(readset)
        elif stage == 'filtered':
            return [self.filtered_lr_path(readset)]
        else:
            raise ConfigError(f"get_nanoplot_input_files :: unknown stage '{stage}' (expected 'raw' or 'filtered').")

    def get_fastqc_sr_input_files(self, readset, stage):
        """Return the short-read FASTQ files FastQC should run on for a readset and stage.

        stage='raw'      → the readset's original r1/r2 files.
        stage='filtered' → the quality-controlled QUALITY_PASSED_R{1,2} files (extension depends on
                           whether gzip_fastqs is on; requires iu_filter_quality_minoche to be
                           enabled). Depending on these paths also creates the DAG edge that forces
                           SR QC to finish first.
        """
        if stage == 'raw':
            d = self.get_sr_files_for_readset(readset)
            return list(d.get('r1', [])) + list(d.get('r2', []))
        elif stage == 'filtered':
            ext = '.fastq.gz' if self.get_param_value_from_config(['gzip_fastqs', 'run']) == True else '.fastq'
            qc_dir = self.dirs_dict["QC_DIR"]
            return [os.path.join(qc_dir, f"{readset}-QUALITY_PASSED_R1{ext}"),
                    os.path.join(qc_dir, f"{readset}-QUALITY_PASSED_R2{ext}")]
        else:
            raise ConfigError(f"get_fastqc_sr_input_files :: unknown stage '{stage}' (expected 'raw' or 'filtered').")

    def _check_stage_flags(self, tool, tool_label, filter_rule, filter_label):
        """Validate the run_on_raw / run_on_filtered combination for one stats tool.

        Only enforced when the stats tool itself is enabled. Raises ConfigError when the tool is
        asked to do nothing, or asked to run on filtered reads that will never be produced because
        its matching filter is disabled.
        """
        if self.get_param_value_from_config([tool, 'run']) != True:
            return

        raw = self.get_param_value_from_config([tool, 'run_on_raw']) == True
        filtered = self.get_param_value_from_config([tool, 'run_on_filtered']) == True
        filter_on = self.get_param_value_from_config([filter_rule, 'run']) == True

        if not raw and not filtered:
            raise ConfigError(
                f"{tool_label} is enabled ('{tool}' → run: true), but both 'run_on_raw' and "
                f"'run_on_filtered' are set to false — so it has nothing to run on. Set at least one "
                f"of them to true, or disable {tool_label} altogether."
            )

        if filtered and not filter_on:
            raise ConfigError(
                f"{tool_label} is set to run on filtered reads ('{tool}' → run_on_filtered: true), but "
                f"{filter_label} ('{filter_rule}' → run) is not enabled, so there will be no filtered reads "
                f"for {tool_label} to look at. Either enable {filter_label}, or set '{tool}' → "
                f"run_on_filtered: false and run_on_raw: true to run {tool_label} on the raw reads instead."
            )

    def sanity_check_qc_stage_flags(self):
        """Validate the raw/filtered stage flags for every stats tool that supports them."""
        self._check_stage_flags('nanoplot', 'NanoPlot', 'filtlong', 'Filtlong long-read filtering')
        self._check_stage_flags('fastqc_sr', 'FastQC', 'iu_filter_quality_minoche', 'illumina-utils short-read quality filtering')

    def sanity_check_filtlong_has_filtering(self):
        """Raise ConfigError if filtlong is enabled but no filtering criteria are given.

        Filtlong is a filter: with no criteria it has nothing to do (it errors out, or at best
        copies the reads through unchanged), so 'run: true' without any parameters is always a
        mistake. Accept any of the explicit length/bases params, or a non-empty additional_params
        (where filtlong's other filters such as --keep_percent / --min_mean_q are passed).
        """
        if self.get_param_value_from_config(['filtlong', 'run']) != True:
            return

        # a param counts as "set" when it has an actual value — including 0. Only None (unset in
        # params.json) or "" (blank) mean "not given", so we must not use plain truthiness here.
        explicit = [self.get_param_value_from_config(['filtlong', k])
                    for k in ['--min-length', '--max-length', '--target-bases']]
        any_explicit = any(v not in (None, "") for v in explicit)
        additional = self.get_param_value_from_config(['filtlong', 'additional_params'])

        if not any_explicit and not (additional and str(additional).strip()):
            raise ConfigError(
                "You set 'filtlong' → run: true, but you didn't give it anything to filter on. "
                "Filtlong is a filter, so with no criteria it has nothing to do (it will error out "
                "or just copy your reads through unchanged). Please set at least one of "
                "'--min-length', '--max-length', or '--target-bases' in the 'filtlong' section of "
                "your config — or, for filtlong's other options (e.g. --keep_percent, --min_mean_q), "
                "put them in 'additional_params'. If you don't actually want to filter your long "
                "reads, set 'filtlong' → run: false instead."
            )

    def qc_producers(self):
        """Return the QC stats producers whose output MultiQC should aggregate.

        Each entry is (parent_dir, readset_ids, stages) for a stats tool (fastqc_sr / nanoplot)
        that will ACTUALLY create output: it is enabled, has matching readsets, and has at least
        one selected stage. A tool that is enabled but produces nothing (no matching readsets, or
        no selected stage) is omitted — so MultiQC is never pointed at a directory that no rule
        creates. Single source of truth for get_qc_target_files() (the per-stage report targets
        and the MultiQC gate) and the multiqc rule (multiqc.smk); adding a QC tool means editing
        this one place.
        """
        producers = []
        if self.get_param_value_from_config(['fastqc_sr', 'run']) == True:
            sr = self.get_sr_readset_ids()
            stages = self._qc_stages_for('fastqc_sr')
            if sr and stages:
                producers.append((os.path.join(self.dirs_dict["QC_DIR"], "fastqc"), sr, stages))
        if self.get_param_value_from_config(['nanoplot', 'run']) == True:
            lr = self.get_lr_readset_ids()
            stages = self._qc_stages_for('nanoplot')
            if lr and stages:
                producers.append((os.path.join(self.dirs_dict["QC_DIR"], "nanoplot"), lr, stages))
        return producers

    def get_qc_target_files(self):
        """Return the list of all QC target files based on enabled options."""
        self.check_qc_program_dependencies()
        self.sanity_check_qc_stage_flags()
        self.sanity_check_filtlong_has_filtering()
        targets = []

        if getattr(self, 'run_qc', False):
            targets.append(os.path.join(self.dirs_dict["QC_DIR"], "qc-report.txt"))

        # A stats tool enabled with no matching readsets produces nothing (and is omitted from
        # qc_producers below); warn so the skip is visible rather than silent.
        if self.get_param_value_from_config(['fastqc_sr', 'run']) == True and not self.get_sr_readset_ids():
            self.run.warning(
                "'fastqc_sr' is enabled, but there are no short-read samples in your samples-txt "
                "for FastQC to run on — it will be skipped."
            )
        if self.get_param_value_from_config(['nanoplot', 'run']) == True and not self.get_lr_readset_ids():
            self.run.warning(
                "'nanoplot' is enabled, but there are no long-read samples in your samples-txt "
                "for NanoPlot to run on — it will be skipped."
            )

        if getattr(self, 'run_filtlong', False):
            # NB: the duplicate-read-name validation is NOT done here — it runs as the
            # check_lr_read_names rule (which gates filtlong), so full long-read files are not
            # re-scanned at parse time on every dry run / DAG rebuild.
            for rs_id in self.get_lr_readset_ids():
                targets.append(self.filtered_lr_path(rs_id))

        # FastQC / NanoPlot write into a per-readset, per-stage directory (see their rules). Derive
        # the target dirs from the single producer list so they can't drift from the multiqc rule.
        producers = self.qc_producers()
        for parent_dir, readset_ids, stages in producers:
            for stage in stages:
                for rs_id in readset_ids:
                    targets.append(os.path.join(parent_dir, rs_id, stage))

        if getattr(self, 'run_multiqc', False):
            # Gate on whether any producer will actually create output (not merely on the 'run'
            # flags): MultiQC scheduled against a never-created directory would fail the run.
            if producers:
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], "multiqc", "multiqc_report.html"))
            else:
                self.run.warning(
                    "MultiQC is enabled but there are no FastQC or NanoPlot outputs for it to "
                    "aggregate — MultiQC will be skipped."
                )

        return targets
