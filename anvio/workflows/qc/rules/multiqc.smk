# QC module — MultiQC aggregation rule
#
# Expects the following in the including Snakefile scope:
#   M             — a QCModule-enabled workflow instance
#   dirs_dict     — M.dirs_dict (or equivalent)
#   rule_log()    — canonical log-path helper
#   w             — anvio.workflows (imported in parent Snakefile)
#   SR_READSETS   — list of SR readset ids
#   LR_READSETS   — list of LR readset ids
#   run_fastqc_sr — bool: whether fastqc_sr was run
#   run_nanoplot  — bool: whether nanoplot was run
#
# Tool:
#   MultiQC — Ewels et al. 2016, PMID 27312411 (https://multiqc.info)

def get_multiqc_inputs(wildcards):
    """Collect all QC output files that MultiQC should aggregate."""
    inputs = []

    if run_fastqc_sr:
        fastqc_dir = os.path.join(dirs_dict["QC_DIR"], "fastqc")
        # fastqc_sr writes reports into a per-readset, per-stage directory; depend on those directories.
        for stage in M._qc_stages_for("fastqc_sr"):
            for rs in SR_READSETS:
                inputs.append(os.path.join(fastqc_dir, rs, stage))

    if run_nanoplot:
        nanoplot_dir = os.path.join(dirs_dict["QC_DIR"], "nanoplot")
        # nanoplot writes reports (incl. NanoStats) into a per-readset, per-stage directory; depend on them.
        for stage in M._qc_stages_for("nanoplot"):
            for rs in LR_READSETS:
                inputs.append(os.path.join(nanoplot_dir, rs, stage))

    return inputs


# get_multiqc_inputs (above) lists specific FILES for Snakemake DAG tracking.
# _multiqc_input_dirs lists the same sources as DIRECTORIES for the MultiQC CLI.
# Both must be updated together whenever a new QC tool is added. Passing the top-level
# parent (fastqc/ or nanoplot/) lets MultiQC recurse into every <readset>/<stage> report so a
# before/after comparison lands in a single report.
_multiqc_input_dirs = []
if run_fastqc_sr:
    _multiqc_input_dirs.append(os.path.join(dirs_dict["QC_DIR"], "fastqc"))
if run_nanoplot:
    _multiqc_input_dirs.append(os.path.join(dirs_dict["QC_DIR"], "nanoplot"))


rule multiqc:
    """Aggregate QC outputs (FastQC and/or NanoPlot) into a single MultiQC report."""
    input:
        get_multiqc_inputs,
    output:
        report=os.path.join(dirs_dict["QC_DIR"], "multiqc", "multiqc_report.html"),
    log:
        rule_log("multiqc", "multiqc"),
    conda:
        w.get_conda_yaml_path(M, "multiqc")
    threads: M.T("multiqc")
    resources:
        nodes=M.T("multiqc"),
    params:
        env_prefix=w.get_conda_env_prefix(M, "multiqc"),
        outdir=os.path.join(dirs_dict["QC_DIR"], "multiqc"),
        indirs=" ".join(_multiqc_input_dirs),
        additional_params=M.get_param_value_from_config(["multiqc", "additional_params"]),
    shell:
        # --dirs/--dirs-depth 2 prefixes each sample name with its <readset>/<stage> directories,
        # yielding '<readset> | <stage>' (e.g. 'S1 | raw', 'S1 | filtered'). Two reasons this matters:
        #   1) NanoPlot writes a fixed 'NanoStats.txt' in every report dir, so without directory
        #      prefixing MultiQC would collapse a readset's raw and filtered reports into one
        #      colliding sample.
        #   2) With the readset first, alphabetical sorting keeps a sample's raw and filtered rows
        #      next to each other (they'd be split apart if the stage led the name).
        # Within a sample the rows sort alphabetically, so 'filtered' lands just above 'raw'.
        # Users can still override naming via multiqc additional_params.
        r"""
        {params.env_prefix} multiqc {params.indirs} --dirs --dirs-depth 2 -o {params.outdir} --force {params.additional_params} >> {log} 2>&1
        """
