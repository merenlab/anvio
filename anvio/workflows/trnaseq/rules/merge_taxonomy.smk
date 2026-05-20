rule anvi_merge_trnaseq:
    """Merge per-sample tRNA-seq databases into a project database."""
    input:
        targets=expand(
            os.path.join(
                os.path.join(dirs_dict["IDENT_DIR"], "{sample_name}"), "IDENT.done"
            ),
            sample_name=M.sample_names,
        ),
    output:
        done=touch(os.path.join(dirs_dict["CONVERT_DIR"], "CONVERT.done")),
    log:
        rule_log("anvi_merge_trnaseq", "anvi_merge_trnaseq"),
    threads: M.T("anvi_merge_trnaseq")
    params:
        trnaseq_dbs=expand(
            os.path.join(
                os.path.join(dirs_dict["IDENT_DIR"], "{sample_name}"),
                "{sample_name}-TRNASEQ.db",
            ),
            sample_name=M.sample_names,
        ),
        project_name=M.get_rule_param("anvi_merge_trnaseq", "--project-name"),
        output_dir=os.path.join(
            dirs_dict["CONVERT_DIR"],
            M.get_param_value_from_config(["anvi_merge_trnaseq", "--project-name"]),
        ),
        max_reported_trna_seeds=M.get_rule_param(
            "anvi_merge_trnaseq", "--max-reported-trna-seeds"
        ),
        overwrite_output_destinations=M.get_rule_param(
            "anvi_merge_trnaseq", "--overwrite-output-destinations"
        ),
        description=M.get_rule_param("anvi_merge_trnaseq", "--description"),
        feature_threshold=M.get_rule_param("anvi_merge_trnaseq", "--feature-threshold"),
        preferred_treatment=M.get_rule_param(
            "anvi_merge_trnaseq", "--preferred-treatment"
        ),
        nonspecific_output=M.get_rule_param(
            "anvi_merge_trnaseq", "--nonspecific-output"
        ),
        min_variation=M.get_rule_param("anvi_merge_trnaseq", "--min-variation"),
        min_third_fourth_nt=M.get_rule_param(
            "anvi_merge_trnaseq", "--min-third-fourth-nt"
        ),
        min_indel_fraction=M.get_rule_param(
            "anvi_merge_trnaseq", "--min-indel-fraction"
        ),
        distance=M.get_rule_param("anvi_merge_trnaseq", "--distance"),
        linkage=M.get_rule_param("anvi_merge_trnaseq", "--linkage"),
    shell:
        "anvi-merge-trnaseq {params.trnaseq_dbs} --output-dir {params.output_dir} {params.project_name} {params.max_reported_trna_seeds} {params.overwrite_output_destinations} {params.description} {params.feature_threshold} {params.preferred_treatment} {params.nonspecific_output} {params.min_variation} {params.min_third_fourth_nt} {params.min_indel_fraction} {params.distance} {params.linkage} -T {threads} >> {log} 2>&1"


rule anvi_run_trna_taxonomy:
    """Assign taxonomy to tRNA sequences in the merged tRNA-seq project."""
    input:
        done=rules.anvi_merge_trnaseq.output.done,
    output:
        done=touch(os.path.join(dirs_dict["CONVERT_DIR"], "TAXONOMY.done")),
    log:
        rule_log("anvi_run_trna_taxonomy", "anvi_run_trna_taxonomy"),
    threads: M.T("anvi_run_trna_taxonomy")
    params:
        contigs_db=os.path.join(
            os.path.join(
                dirs_dict["CONVERT_DIR"],
                M.get_param_value_from_config(
                    ["anvi_merge_trnaseq", "--project-name"]
                ),
            ),
            "CONTIGS.db",
        ),
        min_percent_identity=M.get_rule_param(
            "anvi_run_trna_taxonomy", "--min-percent-identity"
        ),
        max_num_target_sequences=M.get_rule_param(
            "anvi_run_trna_taxonomy", "--max-num-target-sequences"
        ),
        write_buffer_size=M.get_rule_param(
            "anvi_run_trna_taxonomy", "--write-buffer-size"
        ),
        all_hits_output_file=os.path.join(
            os.path.join(
                dirs_dict["CONVERT_DIR"],
                M.get_param_value_from_config(
                    ["anvi_merge_trnaseq", "--project-name"]
                ),
            ),
            "TAXONOMY-HITS.txt",
        ),
    shell:
        "anvi-run-trna-taxonomy -c {params.contigs_db} {params.min_percent_identity} {params.max_num_target_sequences} {params.write_buffer_size} --all-hits-output-file {params.all_hits_output_file} -T {threads} >> {log} 2>&1"


rule anvi_tabulate_trnaseq:
    """Export tabular tRNA-seq results from the merged project databases."""
    input:
        done=rules.anvi_run_trna_taxonomy.output.done,
    output:
        done=touch(os.path.join(dirs_dict["CONVERT_DIR"], "TABULATE.done")),
    log:
        rule_log("anvi_tabulate_trnaseq", "anvi_tabulate_trnaseq"),
    params:
        contigs_db=os.path.join(
            os.path.join(
                dirs_dict["CONVERT_DIR"],
                M.get_param_value_from_config(
                    ["anvi_merge_trnaseq", "--project-name"]
                ),
            ),
            "CONTIGS.db",
        ),
        specific_profile_db=os.path.join(
            os.path.join(
                os.path.join(
                    dirs_dict["CONVERT_DIR"],
                    M.get_param_value_from_config(
                        ["anvi_merge_trnaseq", "--project-name"]
                    ),
                ),
                "SPECIFIC_COVERAGE",
            ),
            "PROFILE.db",
        ),
        nonspecific_profile_db=os.path.join(
            os.path.join(
                os.path.join(
                    dirs_dict["CONVERT_DIR"],
                    M.get_param_value_from_config(
                        ["anvi_merge_trnaseq", "--project-name"]
                    ),
                ),
                "NONSPECIFIC_COVERAGE",
            ),
            "PROFILE.db",
        ),
        output_dir=os.path.join(
            dirs_dict["CONVERT_DIR"],
            M.get_param_value_from_config(["anvi_merge_trnaseq", "--project-name"]),
        ),
        overwrite_output_destinations=M.get_rule_param(
            "anvi_tabulate_trnaseq", "--overwrite-output-destinations"
        ),
    shell:
        "anvi-tabulate-trnaseq -c {params.contigs_db} -s {params.specific_profile_db} -n {params.nonspecific_profile_db} --output-dir {params.output_dir} {params.overwrite_output_destinations} >> {log} 2>&1"
