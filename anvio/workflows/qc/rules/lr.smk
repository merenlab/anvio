# QC module — long-read QC rules (Filtlong)
#
# Expects the following in the including Snakefile scope:
#   M           — a QCModule-enabled workflow instance
#   dirs_dict   — M.dirs_dict (or equivalent)
#   rule_log()  — canonical log-path helper
#   w           — anvio.workflows (imported in parent Snakefile)
#   LR_RS_RE    — wildcard constraint regex for LR readsets
#
# Tools:
#   Filtlong — Wick R, GitHub only (https://github.com/rrwick/Filtlong)


rule filtlong:
    """
    Filter long reads using Filtlong to remove low-quality or short reads.

    Output is written to {QC_DIR}/{readset}-FILTERED_LR.fastq.gz.
    When filtlong is enabled, M.get_fastq() returns these filtered reads
    for downstream mapping and assembly.
    """
    input:
        reads=lambda wildcards: M.get_lr_files_for_readset(wildcards.readset),
    output:
        filtered=os.path.join(dirs_dict["QC_DIR"], "{readset}-FILTERED_LR.fastq.gz"),
    log:
        rule_log("filtlong", "{readset}-filtlong"),
    wildcard_constraints:
        readset=LR_RS_RE,
    conda:
        w.get_conda_yaml_path(M, "filtlong")
    threads: M.T("filtlong")
    resources:
        nodes=M.T("filtlong"),
    params:
        reads=lambda wildcards, input: " ".join(input.reads),
        env_prefix=w.get_conda_env_prefix(M, "filtlong"),
        # Config keys use anvi'o dashes (--min-length) but Filtlong CLI requires underscores (--min_length).
        min_length=(
            '--min_length ' + str(M.get_param_value_from_config(["filtlong", "--min-length"]))
            if M.get_param_value_from_config(["filtlong", "--min-length"]) else ''
        ),
        max_length=(
            '--max_length ' + str(M.get_param_value_from_config(["filtlong", "--max-length"]))
            if M.get_param_value_from_config(["filtlong", "--max-length"]) else ''
        ),
        target_bases=(
            '--target_bases ' + str(M.get_param_value_from_config(["filtlong", "--target-bases"]))
            if M.get_param_value_from_config(["filtlong", "--target-bases"]) else ''
        ),
        additional_params=M.get_param_value_from_config(["filtlong", "additional_params"]),
    shell:
        r"""
        {params.env_prefix} filtlong {params.min_length} {params.max_length} \
            {params.target_bases} {params.additional_params} \
            {params.reads} 2>> {log} | gzip > {output.filtered}
        """
