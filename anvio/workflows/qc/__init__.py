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
      - Optional MultiQC aggregation of all QC outputs

    Subclasses must set before any QC methods are called:
      - self.readsets: list[dict] from SamplesTxt.iter_readsets()
      - self.dirs_dict: includes "QC_DIR"
      - self.run_qc: bool (SR QC enabled)
      - self.run_filtlong: bool
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

    def get_fastqc_sr_input_files(self, readset):
        """Return the short-read FASTQ files FastQC should run on for a readset.

        When short-read QC is enabled, these are the quality-controlled QUALITY_PASSED_R{1,2} files
        (with the extension depending on whether gzip_fastqs is on) — depending on them also creates
        the DAG edge that forces SR QC to finish first. When SR QC is disabled, they are the
        readset's raw r1/r2 files, so FastQC still has something to report on rather than the
        workflow failing on missing QUALITY_PASSED inputs.
        """
        if self.get_param_value_from_config(['iu_filter_quality_minoche', 'run']) == True:
            ext = '.fastq.gz' if self.get_param_value_from_config(['gzip_fastqs', 'run']) == True else '.fastq'
            qc_dir = self.dirs_dict["QC_DIR"]
            return [os.path.join(qc_dir, f"{readset}-QUALITY_PASSED_R1{ext}"),
                    os.path.join(qc_dir, f"{readset}-QUALITY_PASSED_R2{ext}")]

        # SR QC disabled: run FastQC on the raw short reads instead
        d = self.get_sr_files_for_readset(readset)
        return list(d.get('r1', [])) + list(d.get('r2', []))

    def get_qc_target_files(self):
        """Return the list of all QC target files based on enabled options."""
        self.check_qc_program_dependencies()
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
            # FastQC writes into a per-readset directory (see the fastqc_sr rule); target the dir.
            for rs_id in sr_readset_ids:
                targets.append(os.path.join(fastqc_dir, rs_id))

        if getattr(self, 'run_filtlong', False):
            self.check_lr_readsets_no_duplicate_names()
            for rs_id in self.get_lr_readset_ids():
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], f"{rs_id}-FILTERED_LR.fastq.gz"))

        if getattr(self, 'run_multiqc', False):
            if run_fastqc_sr:
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], "multiqc", "multiqc_report.html"))
            else:
                self.run.warning(
                    "MultiQC is enabled but 'fastqc_sr' is not — MultiQC has no compatible inputs "
                    "to aggregate and will be skipped."
                )

        return targets
