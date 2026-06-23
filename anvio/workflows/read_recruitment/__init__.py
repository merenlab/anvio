"""
    Classes for the read recruitment module — a mixin that provides
    indexing, mapping, BAM processing, profiling, and merging rules
    reusable across workflows.
"""

import os
import anvio
from anvio.workflows import WorkflowSuperClass


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


class ReadRecruitmentModule(WorkflowSuperClass):
    """Mixin providing read recruitment: indexing -> mapping -> BAM -> profile -> merge.

    Subclasses must set up (before calling this __init__):
      - self.readsets: list[dict] with keys id, type (SR/LR), reads, base_sample
      - self.group_names: list[str]
      - self.group_sizes: dict[str, int]
      - self.references_mode: bool
      - self.references_for_removal: dict (can be empty)
      - self.remove_short_reads_based_on_references: bool
      - self.run_qc: bool
      - self.fasta_information: dict[str, dict]
    """

    def __init__(self):
        self.rules.extend([
            'bowtie_build', 'bowtie',
            'minimap2_index', 'minimap2',
            'samtools_view', 'anvi_init_bam',
            'anvi_profile', 'anvi_merge',
            'import_percent_of_reads_mapped',
        ])

        rule_acceptable_params_dict = {}

        rule_acceptable_params_dict['bowtie'] = ["conda_yaml", "conda_env", "additional_params"]
        rule_acceptable_params_dict['bowtie_build'] = ["additional_params"]
        rule_acceptable_params_dict['minimap2_index'] = ["additional_params"]
        rule_acceptable_params_dict['minimap2'] = ["preset","conda_yaml", "conda_env", "additional_params", "threads"]
        rule_acceptable_params_dict['samtools_view'] = ["additional_params"]
        rule_acceptable_params_dict['anvi_profile'] = ["--overwrite-output-destinations", "--report-variability-full",
                                                        "--skip-SNV-profiling", "--profile-SCVs", "--description",
                                                        "--skip-hierarchical-clustering", "--distance", "--linkage", "--min-contig-length",
                                                        "--min-mean-coverage", "--min-coverage-for-variability", "--cluster-contigs",
                                                        "--contigs-of-interest", "--queue-size", "--write-buffer-size", "--write-buffer-size-per-thread",
                                                        "--fetch-filter", "--min-percent-identity", "--max-contig-length"]
        rule_acceptable_params_dict['anvi_merge'] = ["--sample-name", "--description", "--skip-hierarchical-clustering",
                                                     "--enforce-hierarchical-clustering", "--distance", "--linkage",
                                                     "--overwrite-output-destinations"]
        rule_acceptable_params_dict['import_percent_of_reads_mapped'] = ["run"]

        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.dirs_dict.update({"QC_DIR": "01_QC",
                                "MAPPING_DIR": "04_MAPPING",
                                "PROFILE_DIR": "05_ANVIO_PROFILE",
                                "MERGE_DIR": "06_MERGED"})

        self.default_config.update({
            "bowtie": {"additional_params": "--no-unal", "threads": 3},
            "minimap2_index": {"additional_params": ""},
            "minimap2": {"threads": 3, "preset": "map-hifi", "additional_params": "--secondary-seq"},
            "samtools_view": {"additional_params": "-F 4"},
            "anvi_profile": {"threads": 3, "--overwrite-output-destinations": True},
            "anvi_merge": {"--sample-name": "{group}", "--overwrite-output-destinations": True},
            "import_percent_of_reads_mapped": {"run": True},
        })

    def get_minimap2_preset(self, readset_id):
        """Return the minimap2 preset from the config (same for every readset)."""
        return self.get_param_value_from_config(['minimap2', 'preset'])

    def get_sr_readset_ids(self):
        return [rs['id'] for rs in self.readsets if rs['type'] == 'SR']

    def get_lr_readset_ids(self):
        return [rs['id'] for rs in self.readsets if rs['type'] == 'LR']

    def get_readset_ids(self):
        return [rs['id'] for rs in self.readsets]

    def get_sr_files_for_readset(self, readset_id):
        rs = self.readsets_by_id.get(readset_id)
        return {
            "r1": list(rs['reads'].get('r1', [])),
            "r2": list(rs['reads'].get('r2', [])),
        }

    def get_lr_files_for_readset(self, readset_id):
        rs = self.readsets_by_id.get(readset_id)
        return list(rs['reads'].get('lr', []))

    def get_fastq(self, readset, pre_ref_removal=False):
        """Return FASTQ paths for a readset.

        Default implementation returns raw input paths.
        Subclasses should override to apply QC/filtered path logic.
        """
        rs = self.readsets_by_id.get(readset)
        if rs['type'] == 'SR':
            return self.get_sr_files_for_readset(readset)
        elif rs['type'] == 'LR':
            return {'lr': self.get_lr_files_for_readset(readset)}

    def get_input_fasta_path(self, wildcards, remove_gz_suffix=True):
        """Return the input FASTA path for a group. Subclasses may override."""
        group = wildcards.group
        if group in self.fasta_information:
            path = self.fasta_information[group]['path']
            if remove_gz_suffix and path.endswith('.gz'):
                return path[:-3]
            return path
        return os.path.join(self.dirs_dict.get("FASTA_DIR", "01_FASTA"), group, "final.contigs.fa")

    def get_fasta(self, wildcards):
        """Return the contigs FASTA path for a group. Subclasses may override."""
        return self.get_input_fasta_path(wildcards)

    def get_contigs_db_path(self):
        return os.path.join(self.dirs_dict["CONTIGS_DIR"], "{group}.db")

    def get_readsets_for_mapping_to_group(self, group_id):
        """Return list of readset IDs to map/profile for a group.

        Subclasses should override with specific grouping logic
        (e.g., all_against_all, group-based).
        """
        return self.get_readset_ids()

    def get_merge_optional_inputs(self):
        """Override to inject additional named inputs into merge/README rules.

        Returns dict of {input_name: (flag_file_name, run_flag)}.
        If run_flag is True the flag file is required as a real dependency;
        otherwise a placeholder (contigs DB) is substituted.
        """
        return {}

    def input_for_anvi_merge(self, wildcards):
        """Return the list of per-readset PROFILE.db paths to merge for {group}."""
        if self.get_param_value_from_config(['all_against_all']):
            member_readsets = self.get_readset_ids()
        else:
            member_readsets = self.get_readsets_for_mapping_to_group(wildcards.group)
        return [os.path.join(self.dirs_dict["PROFILE_DIR"], f"{wildcards.group}/{rs}", "PROFILE.db")
                for rs in member_readsets]

    @property
    def readsets_by_id(self):
        if not hasattr(self, '_readsets_by_id_cache'):
            self._readsets_by_id_cache = {rs['id']: rs for rs in self.readsets}
        return self._readsets_by_id_cache
