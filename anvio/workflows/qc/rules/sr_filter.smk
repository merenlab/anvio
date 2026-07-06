# QC module — short-read filtering rules
#
# Expects the following in the including Snakefile scope:
#   M               — a QCModule-enabled workflow instance
#   dirs_dict       — M.dirs_dict (or equivalent)
#   rule_log()      — canonical log-path helper
#   u               — anvio.utils (imported in parent Snakefile)
#   SR_READSETS     — list of SR readset ids
#   SR_RS_RE        — wildcard constraint regex for SR readsets
#   run_gzip_fastqs — bool: whether gzip_fastqs is enabled
#   run_fastqc_sr   — bool: whether fastqc_sr is enabled


rule iu_gen_configs:
    """
    Generate QC .ini files for SHORT-READ readsets.

    Runs iu-gen-configs once (writes {QC_DIR}/{base_sample}.ini), then copies
    to {QC_DIR}/{readset}.ini so suffixed SR readsets also have their own file.
    """
    input:
        source=ancient(M.get_param_value_from_config(["samples_txt"])),
    output:
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
        import shutil

        shell(
            "iu-gen-configs {input.source} -o {params.dir} {params.r2_prefix} {params.r1_prefix} >> {log} 2>&1"
        )
        for rs in M.get_sr_readset_ids():
            base = M.readsets_by_id[rs]["base_sample"]
            src = os.path.join(dirs_dict["QC_DIR"], f"{base}.ini")
            dst = os.path.join(dirs_dict["QC_DIR"], f"{rs}.ini")
            if os.path.abspath(src) == os.path.abspath(dst):
                if not os.path.exists(dst):
                    raise ConfigError(
                        f"Anvi'o ran iu-gen-configs and expected it to produce a .ini file at '{dst}' "
                        f"for the unsuffixed SR readset '{rs}', but that file is nowhere to be found. "
                        f"This is almost certainly a bug in iu-gen-configs or an unexpected change in "
                        f"how it names its output files. Please check the log at {log} for clues."
                    )
            else:
                if not os.path.exists(src):
                    raise ConfigError(
                        f"Anvi'o ran iu-gen-configs and expected it to produce a base .ini file at "
                        f"'{src}' for the base sample '{base}' (parent of readset '{rs}'), but that "
                        f"file is missing. This usually means iu-gen-configs failed or wrote its output "
                        f"somewhere unexpected. Please check the log at {log} for the full story."
                    )
                shutil.copyfile(src, dst)


def input_for_qc(readset):
    """Return dict of inputs for iu-filter-quality-minoche for a given SR readset id."""
    input_dict = {"ini": ancient(os.path.join(dirs_dict["QC_DIR"], f"{readset}.ini"))}
    fastq_files = M.get_sr_files_for_readset(readset)
    input_dict.update(fastq_files)
    return input_dict


qc_output_r1 = os.path.join(dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_R1.fastq")
qc_output_r2 = os.path.join(dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_R2.fastq")


rule iu_filter_quality_minoche:
    """Run QC using iu-filter-quality-minoche."""
    input:
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
        base_sample=lambda wildcards: M.readsets_by_id[wildcards.readset]["base_sample"],
    shell:
        r"""
        iu-filter-quality-minoche {input.ini} {params.store_read_fate} {params.print_qual_scores} \
                                  {params.limit_num_pairs} {params.visualize_quality_curves} \
                                  {params.ignore_deflines} >> {log} 2>&1

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
    """Aggregate iu-filter-quality-minoche statistics across SR readsets."""
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
        for filename in input.targets:
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
            report_dict, output.report, headers=["sample"] + report_column_headers
        )


gzip_fastq_output = os.path.join(
    dirs_dict["QC_DIR"], "{readset}-QUALITY_PASSED_{R}.fastq.gz"
)


rule gzip_fastqs:
    """Compress quality-controlled short-read FASTQ files."""
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
    threads: M.T("gzip_fastqs")
    resources:
        nodes=M.T("gzip_fastqs"),
    shell:
        "gzip {input.fastq} >> {log} 2>&1"


fastqc_sr_output_dir = os.path.join(dirs_dict["QC_DIR"], "fastqc")

if run_fastqc_sr:
    rule fastqc_sr:
        """Run FastQC on short reads for MultiQC input.

        The input is resolved by M.get_fastqc_sr_input_files(): the quality-controlled
        QUALITY_PASSED reads when SR QC is enabled (which also creates the DAG edge that forces
        iu_filter_quality_minoche / gzip_fastqs to finish first), or the raw short reads when SR QC
        is disabled. The output is a per-readset directory rather than named files, because FastQC
        derives report filenames from the input basenames (which vary for raw / multi-file readsets);
        we let it write whatever it produces into {readset}/ and MultiQC aggregates by scanning the
        parent fastqc directory.
        """
        input:
            reads=lambda wildcards: M.get_fastqc_sr_input_files(wildcards.readset),
        output:
            report_dir=directory(os.path.join(fastqc_sr_output_dir, "{readset}")),
        log:
            rule_log("fastqc_sr", "{readset}-fastqc_sr"),
        wildcard_constraints:
            readset=SR_RS_RE,
        threads: M.T("fastqc_sr")
        resources:
            nodes=M.T("fastqc_sr"),
        params:
            additional_params=M.get_param_value_from_config(["fastqc_sr", "additional_params"]),
        shell:
            r"""
            mkdir -p {output.report_dir}
            fastqc -o {output.report_dir} -t {threads} {params.additional_params} {input.reads} >> {log} 2>&1
            """
