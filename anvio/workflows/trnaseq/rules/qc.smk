rule make_iu_input:
    """Create an Illumina-utils samples file for each sequence library."""
    input:
        source=ancient(M.get_param_value_from_config(["samples_txt"])),
    output:
        target=os.path.join(
            os.path.join(dirs_dict["QC_DIR"], "{sample_name}"), "iu_samples_input.txt"
        ),
    log:
        rule_log("make_iu_input", "make_iu_input"),
    run:
        library_num = M.sample_names.index(wildcards.sample_name)
        r1_path = M.r1_paths[library_num]
        r2_path = M.r2_paths[library_num]
        iu_samples_input_df = pd.DataFrame(
            [[wildcards.sample_name, r1_path, r2_path]],
            columns=["sample", "r1", "r2"],
        )
        iu_samples_input_df.to_csv(output[0], sep="\t", index=False)


rule iu_gen_configs:
    """Create an Illumina-utils config file from each Illumina-utils samples file representing each sequence library."""
    input:
        samples=rules.make_iu_input.output.target,
    output:
        target=os.path.join(
            os.path.join(dirs_dict["QC_DIR"], "{sample_name}"), "{sample_name}.ini"
        ),
    log:
        rule_log("iu_gen_configs", "iu_gen_configs"),
    run:
        out_dir = os.path.join(dirs_dict["QC_DIR"], wildcards.sample_name)
        library_num = M.sample_names.index(wildcards.sample_name)
        r1_prefix = M.r1_prefixes[library_num] if M.r1_prefixes else None
        r2_prefix = M.r2_prefixes[library_num] if M.r2_prefixes else None
        shell(
            "iu-gen-configs {input} -o {out_dir} %s %s >> {log} 2>&1"
            % (
                f"{'--r1-prefix ' + r1_prefix if r1_prefix else ''}",
                f"{'--r2-prefix ' + r2_prefix if r2_prefix else ''}",
            )
        )


rule iu_merge_pairs:
    """Merge paired-end reads using Illumina-utils."""
    input:
        source=ancient(
            os.path.join(
                os.path.join(dirs_dict["QC_DIR"], "{sample_name}"), "{sample_name}.ini"
            )
        ),
    output:
        done=touch(
            os.path.join(
                os.path.join(dirs_dict["QC_DIR"], "{sample_name}"), "MERGE.done"
            )
        ),
    log:
        rule_log("iu_merge_pairs", "iu_merge_pairs"),
    threads: M.T("iu_merge_pairs")
    params:
        marker_gene_stringent=M.get_rule_param(
            "iu_merge_pairs", "--marker-gene-stringent"
        ),
        # allows for both full and partial overlap of inserts, by default trimming trailing adapters following fully overlapping inserts
        max_num_mismatches=M.get_rule_param("iu_merge_pairs", "--max-num-mismatches"),
        report_r1_prefix=M.get_rule_param("iu_merge_pairs", "--report-r1-prefix"),  # flag for reporting the actual sequence of an adapter tag at the beginning of the forward read, which should be specified in iu_gen_configs
        report_r2_prefix=M.get_rule_param("iu_merge_pairs", "--report-r2-prefix"),  # flag for reporting the actual sequence of an adapter tag at the beginning of the reverse read
    run:
        shell(
            "iu-merge-pairs {input} {params.marker_gene_stringent} {params.max_num_mismatches} {params.report_r1_prefix} {params.report_r2_prefix} --num-threads {threads} >> {log} 2>&1"
        )
        if M.gzip_iu_merge_pairs_output:
            merged_fasta = os.path.join(
                os.path.join(dirs_dict["QC_DIR"], wildcards.sample_name),
                wildcards.sample_name + "_MERGED",
            )
            out_dir = os.path.join(dirs_dict["QC_DIR"], wildcards.sample_name)
            if M.run_anvi_reformat_fasta:
                shell("gzip -c %s > %s.gz 2>>{log}" % (merged_fasta, merged_fasta))
            else:
                shell("gzip -f %s >> {log} 2>&1" % merged_fasta)
            shell(
                "gzip -f %s >> {log} 2>&1"
                % os.path.join(out_dir, wildcards.sample_name + "_FAILED")
            )
            shell(
                "gzip -f %s >> {log} 2>&1"
                % os.path.join(out_dir, wildcards.sample_name + "_FAILED_WITH_Ns")
            )
            if params.report_r1_prefix:
                shell(
                    "gzip -f %s >> {log} 2>&1"
                    % os.path.join(
                        out_dir, wildcards.sample_name + "_MERGED_R1_PREFIX"
                    )
                )
            if params.report_r2_prefix:
                shell(
                    "gzip -f %s >> {log} 2>&1"
                    % os.path.join(
                        out_dir, wildcards.sample_name + "_MERGED_R2_PREFIX"
                    )
                )


rule gen_qc_report:
    """Report all quality control statistics (run once)."""
    input:
        targets=expand(
            os.path.join(
                os.path.join(dirs_dict["QC_DIR"], "{sample_name}"), "MERGE.done"
            ),
            sample_name=M.sample_names,
        ),
    output:
        report=os.path.join(dirs_dict["QC_DIR"], "qc_report.txt"),  # optional target file that triggers Illumina-utils QC steps
    log:
        rule_log("gen_qc_report", "gen_qc_report"),
    run:
        report_dict = {}
        headers = []
        for i, input_filepath in enumerate(input):
            input_dirname = os.path.dirname(input_filepath)
            sample_name = os.path.basename(input_dirname)
            stats_filepath = os.path.join(input_dirname, sample_name + "_STATS")
            report_dict[sample_name] = {}
            file_headers = []
            with open(stats_filepath) as f:
                firstline = True
                for line in f:
                    if line == "\n":
                        break
                    line_frags = line.rstrip().split(" ...")
                    header = line_frags[0]
                    file_headers.append(header)
                    number = line_frags[1].split("\t")[1]
                    report_dict[sample_name][header] = number
            if i == 0:
                headers = file_headers
            else:
                if file_headers != headers:
                    raise ConfigError(
                        "The difference in output headers between STATS files "
                        "indicates an inconsistency in how files were processed by 'iu-merge-pairs'. "
                        "These files, for example, have a difference between their headers: "
                        "%s and %s" % (input[i], input[i - 1])
                    )
        u.store_dict_as_TAB_delimited_file(
            report_dict, output.report, headers=["sample"] + headers
        )
