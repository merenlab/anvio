# QC module — long-read QC rules (LongQC + Filtlong)
#
# Expects the following in the including Snakefile scope:
#   M           — a QCModule-enabled workflow instance
#   dirs_dict   — M.dirs_dict (or equivalent)
#   rule_log()  — canonical log-path helper
#   w           — anvio.workflows (imported in parent Snakefile)
#   LR_RS_RE    — wildcard constraint regex for LR readsets
#
# Tools:
#   LongQC   — Fukasawa et al. 2020, PMID 32663312 (https://github.com/yfukasawa/LongQC)
#   Filtlong — Wick R, GitHub only (https://github.com/rrwick/Filtlong)


longqc_output_dir = os.path.join(dirs_dict["QC_DIR"], "longqc")


rule longqc:
    """
    Run LongQC quality assessment on long reads.

    Platform is resolved per-readset from the lr_technology column in samples.txt
    via M.get_longqc_platform(). LongQC must be installed and accessible via
    conda_yaml/conda_env or directly on $PATH.
    """
    input:
        reads=lambda wildcards: M.get_lr_files_for_readset(wildcards.readset),
    output:
        log_file=os.path.join(longqc_output_dir, "{readset}", "log.txt"),
    log:
        rule_log("longqc", "{readset}-longqc"),
    wildcard_constraints:
        readset=LR_RS_RE,
    conda:
        w.get_conda_yaml_path(M, "longqc")
    threads: M.T("longqc")
    resources:
        nodes=M.T("longqc"),
    params:
        platform=lambda wildcards: M.get_longqc_platform(wildcards.readset),
        outdir=lambda wildcards: os.path.join(longqc_output_dir, wildcards.readset),
        reads=lambda wildcards, input: " ".join(input.reads),
        env_prefix=w.get_conda_env_prefix(M, "longqc"),
        additional_params=M.get_param_value_from_config(["longqc", "additional_params"]),
    shell:
        r"""
        # LongQC refuses to run if its output directory already exists.
        rm -rf {params.outdir}
        # conda installs LongQC.py without execute permission, so we locate it with `which` and run
        # it via python. The `bash -c` wrapper makes the command substitution run INSIDE the conda
        # environment ({params.env_prefix}); without it, `$(which LongQC.py)` would resolve in the
        # outer shell and break the `conda_env` case. When env_prefix is empty (PATH or conda_yaml,
        # where Snakemake already activated the env), it simply runs in the current environment.
        {params.env_prefix} bash -c 'python "$(which LongQC.py)" sampleqc \
            -x {params.platform} \
            -o {params.outdir} \
            -p {threads} \
            {params.additional_params} \
            {params.reads}' >> {log} 2>&1
        """


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
