# QC module — MultiQC aggregation rule
#
# Expects the following in the including Snakefile scope:
#   M          — a QCModule-enabled workflow instance
#   dirs_dict  — M.dirs_dict (or equivalent)
#   rule_log() — canonical log-path helper
#   SR_READSETS     — list of SR readset ids
#   LR_READSETS     — list of LR readset ids
#   run_fastqc_sr   — bool: whether fastqc_sr was run
#   run_lr_qc       — bool: whether longqc was run
#
# Tool:
#   MultiQC — Ewels et al. 2016, PMID 27312411 (https://multiqc.info)

import os


def _multiqc_inputs(wildcards):
    """Collect all QC output files that MultiQC should aggregate."""
    inputs = []

    if run_fastqc_sr:
        fastqc_dir = os.path.join(dirs_dict["QC_DIR"], "fastqc")
        for rs in SR_READSETS:
            inputs.append(os.path.join(fastqc_dir, f"{rs}-QUALITY_PASSED_R1_fastqc.zip"))
            inputs.append(os.path.join(fastqc_dir, f"{rs}-QUALITY_PASSED_R2_fastqc.zip"))

    if run_lr_qc:
        for rs in LR_READSETS:
            inputs.append(os.path.join(dirs_dict["QC_DIR"], "longqc", rs, "log.txt"))

    return inputs


rule multiqc:
    """Aggregate QC outputs from FastQC and LongQC into a single MultiQC report.

    Citation: Ewels et al. 2016, PMID 27312411
    """
    input:
        _multiqc_inputs,
    output:
        report=os.path.join(dirs_dict["QC_DIR"], "multiqc", "multiqc_report.html"),
    log:
        rule_log("multiqc", "multiqc"),
    threads: M.T("multiqc")
    resources:
        nodes=M.T("multiqc"),
    params:
        outdir=os.path.join(dirs_dict["QC_DIR"], "multiqc"),
        indir=dirs_dict["QC_DIR"],
        additional_params=M.get_param_value_from_config(["multiqc", "additional_params"]),
    shell:
        r"""
        multiqc {params.indir} -o {params.outdir} --force {params.additional_params} >> {log} 2>&1
        """
