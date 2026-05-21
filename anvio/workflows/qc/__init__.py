"""
    QCModule — a mixin providing quality-control steps for short-read and long-read workflows.
"""

import os
import gzip
import importlib.util
import anvio
import anvio.utils as u

from anvio.workflows import WorkflowSuperClass, LR_TECHNOLOGY_MAP
from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


class QCModule(WorkflowSuperClass):
    """Mixin providing quality-control steps reusable across workflows.

    Provides:
      - illumina-utils SR filtering (iu-gen-configs + iu-filter-quality-minoche)
      - QC report aggregation (metagenomics_gen_qc_report)
      - Optional gzip of QC'd SR reads
      - Optional FastQC on SR reads
      - Optional LongQC quality assessment on LR reads
      - Optional Filtlong filtering of LR reads

    Subclasses must set before any QC methods are called:
      - self.readsets: list[dict] from SamplesTxt.iter_readsets()
      - self.dirs_dict: includes "QC_DIR"
      - self.run_qc: bool (SR QC enabled)
      - self.run_lr_qc: bool (LongQC enabled)
      - self.run_filtlong: bool
    """

    def __init__(self):
        self.rules.extend([
            'iu_gen_configs',
            'iu_filter_quality_minoche',
            'metagenomics_gen_qc_report',
            'gzip_fastqs',
            'fastqc_sr',
            'longqc',
            'filtlong',
        ])

        rule_acceptable_params_dict = {}
        rule_acceptable_params_dict['iu_gen_configs'] = ["--r1-prefix", "--r2-prefix"]
        rule_acceptable_params_dict['iu_filter_quality_minoche'] = [
            'run', '--visualize-quality-curves', '--ignore-deflines',
            '--limit-num-pairs', '--print-qual-scores', '--store-read-fate',
        ]
        rule_acceptable_params_dict['metagenomics_gen_qc_report'] = []
        rule_acceptable_params_dict['gzip_fastqs'] = ["run"]
        rule_acceptable_params_dict['fastqc_sr'] = ["run", "additional_params"]
        rule_acceptable_params_dict['longqc'] = [
            "run", "conda_yaml", "conda_env", "additional_params",
        ]
        rule_acceptable_params_dict['filtlong'] = [
            "run", "conda_yaml", "conda_env",
            "--min-length", "--max-length", "--target-bases", "additional_params",
        ]

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.default_config.update({
            'iu_filter_quality_minoche': {"run": True, "--ignore-deflines": True},
            'gzip_fastqs': {"run": True},
            'fastqc_sr': {"run": False, "additional_params": ""},
            'longqc': {"run": False, "conda_yaml": "", "conda_env": "", "additional_params": ""},
            'filtlong': {
                "run": False, "conda_yaml": "", "conda_env": "",
                "--min-length": "", "--max-length": "", "--target-bases": "",
                "additional_params": "",
            },
        })

    def get_longqc_platform(self, readset_id):
        """Return the LongQC --sample_type (-x) value for this readset.

        Supported presets: ont-ligation, ont-rapid, ont-1dsq, pb-rs2, pb-sequel.

        pb-hifi is intentionally not supported: LongQC's pb-hifi preset always triggers
        spike-in control filtering against an instrument reference, but PacBio HiFi library
        preparation does not include spike-in controls. Biological HiFi data produces an
        empty coverage file which causes LongQC to crash (EmptyDataError). pb-hifi readsets
        are blocked earlier in check_qc_program_dependencies() with an actionable error.
        """
        rs = self.readsets_by_id.get(readset_id)
        if rs is None:
            raise ConfigError(
                f"Anvi'o was asked to look up the LongQC platform for a readset called '{readset_id}', "
                f"but that readset id does not exist in the workflow. This is most likely a bug — please "
                f"let a developer know."
            )
        tech = rs.get('lr_technology')
        if not tech:
            raise ConfigError(
                f"Anvi'o tried to determine the LongQC platform for the long-read readset '{readset_id}', "
                f"but no 'lr_technology' value was found for that sample. LongQC requires knowing the "
                f"sequencing platform to properly assess read quality. Please make sure your samples-txt "
                f"file has a 'lr_technology' column and that every long-read sample has a valid entry. "
                f"Valid values are: {', '.join(sorted(LR_TECHNOLOGY_MAP.keys()))}."
            )
        if tech not in LR_TECHNOLOGY_MAP:
            raise ConfigError(
                f"Anvi'o doesn't know what to do with the lr_technology value '{tech}' for readset "
                f"'{readset_id}'. This value needs to map to a LongQC platform and a minimap2 preset, "
                f"and '{tech}' is not in anvi'o's technology map. The valid values are: "
                f"{', '.join(sorted(LR_TECHNOLOGY_MAP.keys()))}. If you believe a technology is missing "
                f"from this list, please open an issue on the anvi'o GitHub."
            )
        return LR_TECHNOLOGY_MAP[tech]['longqc']

    def get_minimap2_preset(self, readset_id):
        """Return the minimap2 preset for this readset.

        Uses lr_technology column if available, otherwise falls back to config.
        """
        rs = self.readsets_by_id.get(readset_id)
        if rs is not None:
            tech = rs.get('lr_technology')
            if tech and tech in LR_TECHNOLOGY_MAP:
                return LR_TECHNOLOGY_MAP[tech]['minimap2']
        return self.get_param_value_from_config(['minimap2', 'preset'])

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

        if self.get_param_value_from_config(['filtlong', 'run']) == True:
            if not u.is_program_exists('filtlong', dont_raise=True):
                missing.append(('filtlong', 'filtlong'))

        if self.get_param_value_from_config(['fastqc_sr', 'run']) == True:
            if not u.is_program_exists('fastqc', dont_raise=True):
                missing.append(('fastqc', 'fastqc_sr'))

        if self.get_param_value_from_config(['longqc', 'run']) == True:
            hifi_readsets = [
                rs_id for rs_id in self.get_lr_readset_ids()
                if self.readsets_by_id.get(rs_id, {}).get('lr_technology', '') == 'pb-hifi'
            ]
            if hifi_readsets:
                raise ConfigError(
                    f"LongQC is enabled but {len(hifi_readsets)} readset(s) use the 'pb-hifi' "
                    f"technology ({', '.join(hifi_readsets)}). LongQC's pb-hifi preset always "
                    f"attempts to filter spike-in control reads, but PacBio HiFi library "
                    f"preparation does not include instrument spike-in controls, so no control "
                    f"reads are present in biological HiFi data. This causes LongQC to produce "
                    f"an empty coverage file and crash. LongQC does work with PacBio RS II "
                    f"(pb-rs2) and Sequel/Sequel II (pb-sequel) data, which do carry spike-in "
                    f"controls. Please set lr_technology to 'pb-rs2' or 'pb-sequel' if "
                    f"appropriate, or disable LongQC ('longqc': {{'run': false}}) for HiFi "
                    f"samples."
                )

            if not u.is_program_exists('LongQC.py', dont_raise=True):
                missing.append(('LongQC.py', 'longqc'))
            else:
                for mod in ['mixem', 'pysam', 'numpy', 'scipy', 'matplotlib']:
                    if importlib.util.find_spec(mod) is None:
                        missing.append((f"Python module '{mod}'", 'longqc'))
            longqc_threads = self.get_param_value_from_config(['longqc', 'threads'])
            if longqc_threads and int(longqc_threads) < 4:
                raise ConfigError(
                    f"LongQC requires at least 4 threads (-p/--ncpu >= 4) but 'longqc.threads' "
                    f"is set to {longqc_threads}. Please set it to 4 or higher in your config."
                )

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

    def get_qc_target_files(self):
        """Return the list of all QC target files based on enabled options."""
        self.check_qc_program_dependencies()
        targets = []

        if getattr(self, 'run_qc', False):
            targets.append(os.path.join(self.dirs_dict["QC_DIR"], "qc-report.txt"))

        run_fastqc_sr = self.get_param_value_from_config(['fastqc_sr', 'run']) == True
        if run_fastqc_sr:
            fastqc_dir = os.path.join(self.dirs_dict["QC_DIR"], "fastqc")
            for rs_id in self.get_sr_readset_ids():
                targets.append(os.path.join(fastqc_dir, f"{rs_id}-QUALITY_PASSED_R1_fastqc.html"))
                targets.append(os.path.join(fastqc_dir, f"{rs_id}-QUALITY_PASSED_R2_fastqc.html"))

        if getattr(self, 'run_lr_qc', False):
            for rs_id in self.get_lr_readset_ids():
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], "longqc", rs_id, "log.txt"))

        if getattr(self, 'run_filtlong', False):
            self.check_lr_readsets_no_duplicate_names()
            for rs_id in self.get_lr_readset_ids():
                targets.append(os.path.join(self.dirs_dict["QC_DIR"], f"{rs_id}-FILTERED_LR.fastq.gz"))

        return targets
