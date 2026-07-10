# QC module — long-read QC rules (Filtlong + NanoPlot)
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
#   NanoPlot — De Coster & Rademakers 2023, PMID 37171891 (https://github.com/wdecoster/NanoPlot)


nanoplot_output_dir = os.path.join(dirs_dict["QC_DIR"], "nanoplot")


rule nanoplot:
    """
    Long-read quality assessment with NanoPlot (one report per readset per stage).

    The {stage} wildcard selects which reads to assess: 'raw' for the readset's original long
    reads, or 'filtered' for the Filtlong output (see M.get_nanoplot_input_files()). Which stages
    actually run is controlled by the 'run_on_raw' / 'run_on_filtered' flags in the nanoplot config
    (validated in QCModule.sanity_check_qc_stage_flags). The output is a per-stage, per-readset
    directory: NanoPlot writes its report, plots and NanoStats there, and MultiQC (if enabled)
    aggregates the NanoStats by scanning the parent nanoplot directory. NanoPlot needs no
    sequencing-technology preset, and must be available on $PATH or via conda_yaml/conda_env.
    """
    input:
        reads=lambda wildcards: M.get_nanoplot_input_files(wildcards.readset, wildcards.stage),
    output:
        report_dir=directory(os.path.join(nanoplot_output_dir, "{stage}", "{readset}")),
    log:
        rule_log("nanoplot", "{stage}-{readset}-nanoplot"),
    wildcard_constraints:
        readset=LR_RS_RE,
        stage="raw|filtered",
    conda:
        w.get_conda_yaml_path(M, "nanoplot")
    threads: M.T("nanoplot")
    resources:
        nodes=M.T("nanoplot"),
    params:
        env_prefix=w.get_conda_env_prefix(M, "nanoplot"),
        reads=lambda wildcards, input: " ".join(input.reads),
        additional_params=M.get_param_value_from_config(["nanoplot", "additional_params"]),
    shell:
        r"""
        {params.env_prefix} NanoPlot -t {threads} --fastq {params.reads} \
            -o {output.report_dir} --tsv_stats {params.additional_params} >> {log} 2>&1
        """


rule check_lr_read_names:
    """
    Pre-flight validation: fail fast if a long-read readset has duplicate read names.

    Filtlong aborts mid-run on the first duplicate with a cryptic error, so this scans the raw
    long reads and raises a clear anvi'o ConfigError (with a seqkit fix) if any are found. It runs
    as its own rule — NOT at parse time — so full long-read files are not re-scanned on every dry
    run / DAG rebuild; on a real run it executes once per readset and gates the filtlong rule. On
    success it writes a small sentinel that filtlong depends on.
    """
    input:
        reads=lambda wildcards: M.get_lr_files_for_readset(wildcards.readset),
    output:
        names_ok=os.path.join(dirs_dict["QC_DIR"], "{readset}-LR_NAMES_OK.flag"),
    log:
        rule_log("check_lr_read_names", "{readset}-check_lr_read_names"),
    wildcard_constraints:
        readset=LR_RS_RE,
    run:
        # raises ConfigError (failing this job) if the readset has duplicate read names
        M.check_lr_readset_no_duplicate_names(wildcards.readset)
        with open(output.names_ok, "w") as f:
            f.write("no duplicate read names found\n")


rule filtlong:
    """
    Filter long reads using Filtlong to remove low-quality or short reads.

    Output is written to {QC_DIR}/{readset}-FILTERED_LR.fastq.gz.
    When filtlong is enabled, M.get_fastq() returns these filtered reads
    for downstream mapping and assembly. Depends on the check_lr_read_names sentinel so the
    duplicate-read-name validation runs (once) before filtlong touches the reads.
    """
    input:
        reads=lambda wildcards: M.get_lr_files_for_readset(wildcards.readset),
        names_ok=os.path.join(dirs_dict["QC_DIR"], "{readset}-LR_NAMES_OK.flag"),
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
