# QC module — MultiQC aggregation rule
#
# Expects the following in the including Snakefile scope:
#   M             — a QCModule-enabled workflow instance
#   dirs_dict     — M.dirs_dict (or equivalent)
#   rule_log()    — canonical log-path helper
#   w             — anvio.workflows (imported in parent Snakefile)
#
# Tool:
#   MultiQC — Ewels et al. 2016, PMID 27312411 (https://multiqc.info)

# The enabled QC producers whose output MultiQC should aggregate: (parent_dir, readset_ids,
# stages), one entry per tool that will actually create output. Computed by QCModule.qc_producers()
# — the SAME source get_qc_target_files() uses to build its targets and MultiQC gate, so the two
# can't drift. Both the per-file DAG dependencies (get_multiqc_inputs) and the parent dirs for the
# MultiQC CLI (_multiqc_input_dirs) are derived from it.
_qc_producers = M.qc_producers()


def get_multiqc_inputs(wildcards):
    """Per-readset, per-stage report directories MultiQC depends on (for Snakemake DAG tracking)."""
    inputs = []
    for parent_dir, readsets, stages in _qc_producers:
        for stage in stages:
            for rs in readsets:
                inputs.append(os.path.join(parent_dir, rs, stage))
    return inputs


# Parent directories passed to the MultiQC CLI; --dirs-depth 2 makes it recurse into every
# <readset>/<stage> report so a before/after comparison lands in a single report.
_multiqc_input_dirs = [parent_dir for parent_dir, _, _ in _qc_producers]


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
