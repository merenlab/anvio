"""
    QCModule — a mixin providing short-read quality-control steps for workflows.
"""

import os
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

    Subclasses must set before any QC methods are called:
      - self.readsets: list[dict] from SamplesTxt.iter_readsets()
      - self.dirs_dict: includes "QC_DIR"
      - self.run_qc: bool (SR QC enabled)
    """

    def __init__(self):
        self.rules.extend([
            'iu_gen_configs',
            'iu_filter_quality_minoche',
            'gen_qc_report',
            'gzip_fastqs',
        ])

        rule_acceptable_params_dict = {}
        rule_acceptable_params_dict['iu_gen_configs'] = ["--r1-prefix", "--r2-prefix"]
        rule_acceptable_params_dict['iu_filter_quality_minoche'] = [
            'run', '--visualize-quality-curves', '--ignore-deflines',
            '--limit-num-pairs', '--print-qual-scores', '--store-read-fate',
        ]
        rule_acceptable_params_dict['gen_qc_report'] = []
        rule_acceptable_params_dict['gzip_fastqs'] = ["run"]

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.default_config.update({
            'iu_filter_quality_minoche': {"run": True, "--ignore-deflines": True},
            'gzip_fastqs': {"run": True},
        })

    def get_fastq(self, readset, pre_ref_removal=False):
        """Return FASTQ paths for a readset.

        For LR: returns the raw long-read files.
        For SR: delegates to _resolve_sr_path() — override there for ref-removal logic.
        """
        rs = self.readsets_by_id.get(readset)
        if rs['type'] == 'LR':
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

    def check_qc_program_dependencies(self):
        """Raise ConfigError for any program required by an enabled QC tool that is missing.

        Only checks tools that are actually enabled in the config, avoiding the
        false-positive noise from the generic check_workflow_program_dependencies.
        """
        missing = []

        if self.get_param_value_from_config(['iu_filter_quality_minoche', 'run']) == True:
            for prog in ['iu-gen-configs', 'iu-filter-quality-minoche']:
                if not u.is_program_exists(prog, dont_raise=True):
                    missing.append((prog, 'iu_filter_quality_minoche'))

        if missing:
            details = '\n  '.join(f"{prog}  (required by: {rule})" for prog, rule in missing)
            raise ConfigError(
                f"Anvi'o found {len(missing)} missing program(s) required by enabled QC steps. "
                f"Please install them before running the workflow:\n\n  {details}"
            )

    def get_qc_target_files(self):
        """Return the list of QC target files based on enabled options."""
        self.check_qc_program_dependencies()
        targets = []

        if getattr(self, 'run_qc', False):
            targets.append(os.path.join(self.dirs_dict["QC_DIR"], "qc-report.txt"))

        return targets
