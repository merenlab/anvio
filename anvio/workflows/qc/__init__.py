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

    def get_fastq(self, readset, pre_ref_removal=False):
        """Return FASTQ paths for a readset.

        For LR: returns filtlong output (if run_filtlong), otherwise raw LR.
        For SR: delegates to _resolve_sr_path() — override there for ref-removal logic.
        """
        rs = self.readsets_by_id.get(readset)
        if rs['type'] == 'LR':
            if getattr(self, 'run_filtlong', False):
                path = os.path.join(self.dirs_dict["QC_DIR"], f"{readset}-FILTERED_LR.fastq.gz")
                return {'lr': [path]}
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

    def _tool_provided_by_conda(self, tool):
        """Whether a tool is supplied via a conda env (conda_yaml or conda_env) rather than $PATH.

        When it is, we must not check for the executable (or its Python modules) on the current
        $PATH / interpreter — Snakemake will run the rule inside the configured environment. Mirrors
        the logic in MetagenomicsWorkflow.ensure_tool_in_path_or_conda().
        """
        y = self.get_param_value_from_config([tool, 'conda_yaml'])
        n = self.get_param_value_from_config([tool, 'conda_env'])
        return bool((y and y.strip()) or (n and n.strip()))

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
            if not u.is_program_exists('fastqc', dont_raise=True):
                missing.append(('fastqc', 'fastqc_sr'))

        if self.get_param_value_from_config(['multiqc', 'run']) == True:
            if not u.is_program_exists('multiqc', dont_raise=True):
                missing.append(('multiqc', 'multiqc'))

        if missing:
            details = '\n  '.join(f"{prog}  (required by: {rule})" for prog, rule in missing)
            raise ConfigError(
                f"Anvi'o found {len(missing)} missing program(s) required by enabled QC steps. "
                f"Please install them before running the workflow:\n\n  {details}"
            )

    def check_lr_readsets_no_duplicate_names(self):
        """Raise ConfigError if any LR readset FASTQ has duplicate read names.

        Filtlong aborts mid-run on the first duplicate it finds, giving a cryptic
        error. This check catches the problem upfront with a clear, actionable message.
        """
        problematic = []
        for rs_id in self.get_lr_readset_ids():
            for path in self.get_lr_files_for_readset(rs_id):
                if not os.path.exists(path):
                    continue
                seen = set()
                opener = gzip.open if path.endswith('.gz') else open
                with opener(path, 'rt') as f:
                    for i, line in enumerate(f):
                        if i % 4 == 0:
                            name = line[1:].split()[0]
                            if name in seen:
                                problematic.append((rs_id, path, name))
                                break
                            seen.add(name)

        if problematic:
            details = '\n  '.join(
                f"readset '{rs}' ({p}): first duplicate: '{n}'"
                for rs, p, n in problematic
            )
            raise ConfigError(
                f"Anvi'o found duplicate read names in {len(problematic)} long-read input "
                f"file(s). Filtlong will crash mid-run when it encounters them, so anvi'o "
                f"refuses to start. Fix the input files with 'seqkit rename' before retrying:\n\n"
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
            return [os.path.join(self.dirs_dict["QC_DIR"], f"{readset}-FILTERED_LR.fastq.gz")]
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

    def get_qc_target_files(self):
        """Return the list of all QC target files based on enabled options."""
        self.check_qc_program_dependencies()
        self.sanity_check_qc_stage_flags()
        targets = []

        if getattr(self, 'run_qc', False):
            targets.append(os.path.join(self.dirs_dict["QC_DIR"], "qc-report.txt"))

        run_fastqc_sr = self.get_param_value_from_config(['fastqc_sr', 'run']) == True
        if run_fastqc_sr:
            sr_readset_ids = self.get_sr_readset_ids()
            if not sr_readset_ids:
                self.run.warning(
                    "'fastqc_sr' is enabled, but there are no short-read samples in your samples-txt "
                    "for FastQC to run on — it will be skipped."
                )
            fastqc_dir = os.path.join(self.dirs_dict["QC_DIR"], "fastqc")
            # FastQC writes into a per-stage, per-readset directory (see the fastqc_sr rule); target the dir.
            for stage in self._qc_stages_for('fastqc_sr'):
                for rs_id in sr_readset_ids:
                    targets.append(os.path.join(fastqc_dir, stage, rs_id))

        if getattr(self, 'run_filtlong', False):
            self.check_lr_readsets_no_duplicate_names()
            for rs_id in self.get_lr_readset_ids():
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], f"{rs_id}-FILTERED_LR.fastq.gz"))

        run_nanoplot = self.get_param_value_from_config(['nanoplot', 'run']) == True
        if run_nanoplot:
            lr_readset_ids = self.get_lr_readset_ids()
            if not lr_readset_ids:
                self.run.warning(
                    "'nanoplot' is enabled, but there are no long-read samples in your samples-txt "
                    "for NanoPlot to run on — it will be skipped."
                )
            nanoplot_dir = os.path.join(self.dirs_dict["QC_DIR"], "nanoplot")
            # NanoPlot writes into a per-stage, per-readset directory (see the nanoplot rule); target the dir.
            for stage in self._qc_stages_for('nanoplot'):
                for rs_id in lr_readset_ids:
                    targets.append(os.path.join(nanoplot_dir, stage, rs_id))

        if getattr(self, 'run_multiqc', False):
            if run_fastqc_sr or run_nanoplot:
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], "multiqc", "multiqc_report.html"))
            else:
                self.run.warning(
                    "MultiQC is enabled but neither 'fastqc_sr' nor 'nanoplot' is — MultiQC has no "
                    "compatible inputs to aggregate and will be skipped."
                )

        return targets
