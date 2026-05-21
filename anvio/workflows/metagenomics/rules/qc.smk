rule iu_gen_configs:
    """
Generate QC .ini files for SHORT-READ readsets.
We run iu-gen-configs once (it writes {QC_DIR}/{base_sample}.ini),
then copy/symlink to {QC_DIR}/{readset}.ini so suffixed SR readsets exist.
"""

    # one-time rule -> no wildcard constraints
    input:
        source=ancient(M.get_param_value_from_config(["samples_txt"])),
    output:
        # ONLY SR readsets (not LR)
        files=expand(
            "{DIR}/{readset}.ini",
            DIR=dirs_dict["QC_DIR"],
            readset=M.get_sr_readset_ids(),
        ),
    log:
        rule_log("iu_gen_configs", "iu_gen_configs"),
    threads: M.T("iu_gen_configs")
    resources:
        nodes=M.T("iu_gen_configs"),
    params:
        dir=dirs_dict["QC_DIR"],
        r1_prefix=M.get_rule_param("iu_gen_configs", "--r1-prefix"),
        r2_prefix=M.get_rule_param("iu_gen_configs", "--r2-prefix"),
    run:
        import os, shutil

        # 1) run iu-gen-configs once; it writes {QC_DIR}/{base_sample}.ini
        shell(
            "iu-gen-configs {input} -o {params.dir} {params.r2_prefix} {params.r1_prefix} >> {log} 2>&1"
        )
        # 2) ensure every SR readset has its own .ini path
        #    (copy base_sample.ini -> readset.ini when readset != base_sample)
        #    This keeps {readset}.ini consistent with the rest of the pipeline.
        for rs in M.get_sr_readset_ids():
            base = M.readsets_by_id[rs]["base_sample"]  # e.g., 'sample_01'
            src = os.path.join(dirs_dict["QC_DIR"], f"{base}.ini")
            dst = os.path.join(dirs_dict["QC_DIR"], f"{rs}.ini")
            if os.path.abspath(src) == os.path.abspath(dst):
                # unsuffixed SR readset (no copy needed)
                if not os.path.exists(dst):
                    # iu-gen-configs should have produced it already; just sanity
                    raise WorkflowError(
                        f"iu-gen-configs did not produce expected ini: {dst}"
                    )
            else:
                # suffixed SR readset; make a copy (or symlink if you prefer)
                if not os.path.exists(src):
                    raise WorkflowError(
                        f"iu-gen-configs did not produce expected base ini: {src}"
                    )
                shutil.copyfile(src, dst)


def input_for_qc(readset):
    """Return dict of inputs for QC for a given SR readset id (handles suffixed names)."""
    d = {"ini": ancient(os.path.join(dirs_dict["QC_DIR"], f"{readset}.ini"))}
    sr = M.get_sr_files_for_readset(readset)  # {'r1': [...], 'r2': [...]}
    d.update(sr)
    return d


qc_output_r1 = os.path.join(dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_R1.fastq")
qc_output_r2 = os.path.join(dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_R2.fastq")


rule iu_filter_quality_minoche:
    """Run QC using iu-filter-quality-minoche"""
    input:
        # ini is read by iu-filter-quality-minoche; r1/r2 here are just for dependency tracking
        unpack(lambda wildcards: input_for_qc(wildcards.readset)),
    output:
        r1=(
            temp(qc_output_r1)
            if M.remove_short_reads_based_on_references
            else qc_output_r1
        ),
        r2=(
            temp(qc_output_r2)
            if M.remove_short_reads_based_on_references
            else qc_output_r2
        ),
        stats=dirs_dict["QC_DIR"] + "/{readset}-STATS.txt",
    log:
        rule_log("iu_filter_quality_minoche", "{readset}-iu_filter_quality_minoche"),
    wildcard_constraints:
        readset=SR_RS_RE,
    threads: M.T("iu_filter_quality_minoche")
    resources:
        nodes=M.T("iu_filter_quality_minoche"),
    params:
        ignore_deflines=M.get_rule_param(
            "iu_filter_quality_minoche", "--ignore-deflines"
        ),
        visualize_quality_curves=M.get_rule_param(
            "iu_filter_quality_minoche", "--visualize-quality-curves"
        ),
        limit_num_pairs=M.get_rule_param(
            "iu_filter_quality_minoche", "--limit-num-pairs"
        ),
        print_qual_scores=M.get_rule_param(
            "iu_filter_quality_minoche", "--print-qual-scores"
        ),
        store_read_fate=M.get_rule_param(
            "iu_filter_quality_minoche", "--store-read-fate"
        ),
        qcdir=dirs_dict["QC_DIR"],
        base_sample=lambda w: M.readsets_by_id[w.readset]["base_sample"],
    shell:
        r"""
        # Run QC using the .ini generated for this readset (contents may still use base sample name)
        iu-filter-quality-minoche {input.ini} {params.store_read_fate} {params.print_qual_scores} \
                                  {params.limit_num_pairs} {params.visualize_quality_curves} \
                                  {params.ignore_deflines} >> {log} 2>&1

        # If the tool wrote files with the *base sample* prefix, copy them to the expected *readset* names.
        # This is a no-op for unsuffixed SR readsets (base==readset).
        if [ ! -f "{output.r1}" ] && [ -f "{params.qcdir}/{params.base_sample}-QUALITY_PASSED_R1.fastq" ]; then
            cp "{params.qcdir}/{params.base_sample}-QUALITY_PASSED_R1.fastq" "{output.r1}"
        fi
        if [ ! -f "{output.r2}" ] && [ -f "{params.qcdir}/{params.base_sample}-QUALITY_PASSED_R2.fastq" ]; then
            cp "{params.qcdir}/{params.base_sample}-QUALITY_PASSED_R2.fastq" "{output.r2}"
        fi
        if [ ! -f "{output.stats}" ] && [ -f "{params.qcdir}/{params.base_sample}-STATS.txt" ]; then
            cp "{params.qcdir}/{params.base_sample}-STATS.txt" "{output.stats}"
        fi
        """


rule gen_qc_report:
    """Aggregate quality-control statistics across readsets or samples."""
    input:
        targets=expand(
            dirs_dict["QC_DIR"] + "/{readset}-STATS.txt", readset=SR_READSETS
        ),
    output:
        report=dirs_dict["QC_DIR"] + "/qc-report.txt",
    log:
        rule_log("gen_qc_report", "gen_qc_report"),
    threads: M.T("gen_qc_report")
    resources:
        nodes=M.T("gen_qc_report"),
    run:
        report_dict = {}
        report_column_headers = [
            "number of pairs analyzed",
            "total pairs passed",
            "total pairs passed (percent of all pairs)",
            "total pair_1 trimmed",
            "total pair_1 trimmed (percent of all passed pairs)",
            "total pair_2 trimmed",
            "total pair_2 trimmed (percent of all passed pairs)",
            "total pairs failed",
            "total pairs failed (percent of all pairs)",
            "pairs failed due to pair_1",
            "pairs failed due to pair_1 (percent of all failed pairs)",
            "pairs failed due to pair_2",
            "pairs failed due to pair_2 (percent of all failed pairs)",
            "pairs failed due to both",
            "pairs failed due to both (percent of all failed pairs)",
            "FAILED_REASON_P",
            "FAILED_REASON_P (percent of all failed pairs)",
            "FAILED_REASON_N",
            "FAILED_REASON_N (percent of all failed pairs)",
            "FAILED_REASON_C33",
            "FAILED_REASON_C33 (percent of all failed pairs)",
        ]
        for filename in input:
            sample = os.path.basename(filename).split("-STATS.txt")[0]
            report_dict[sample] = dict.fromkeys(report_column_headers, 0)
            with open(filename, "r") as f:
                firstline = True
                for line in f.readlines():
                    s1 = line.split(":")
                    numeric_header = s1[0].strip()
                    s2 = s1[1].split("(")
                    numeric = s2[0].strip()
                    report_dict[sample][numeric_header] = numeric
                    if not firstline:
                        s3 = s2[1].split(" ")
                        percent = s3[0].strip("%")
                        percent_header = (
                            numeric_header + " (percent " + " ".join(s3[1:])
                        )
                        percent_header = percent_header.strip()
                        report_dict[sample][percent_header] = percent
                    else:
                        firstline = False
        u.store_dict_as_TAB_delimited_file(
            report_dict, output[0], headers=["sample"] + report_column_headers
        )


gzip_fastq_output = os.path.join(
    dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_{R}.fastq.gz"
)


rule gzip_fastqs:
    """Compressing the quality controlled fastq files"""
    input:
        fastq=os.path.join(dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_{R}.fastq"),
    output:
        target=(
            temp(gzip_fastq_output)
            if M.remove_short_reads_based_on_references
            else gzip_fastq_output
        ),
    log:
        rule_log("gzip_fastqs", "{readset}-{R}-gzip"),
    wildcard_constraints:
        readset=SR_RS_RE,
    threads: M.T("gzip_fastqs")
    resources:
        nodes=M.T("gzip_fastqs"),
    shell:
        "gzip {input.fastq} >> {log} 2>&1"
