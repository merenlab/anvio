rule anvi_import_collection:
    input:
        profile=lambda wildcards: M.profile_databases[wildcards.group],
        contigs=ancient(M.get_contigs_db_path()),
        collection_file=lambda wildcards: M.collections[wildcards.group][
            "collection_file"
        ],
    output:
        imported_collection_flag=touch(
            os.path.join(dirs_dict["MERGE_DIR"], "{group}", "collection-import.done")
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_import_collection.log",
    threads: 1
    resources:
        nodes=1,
    params:
        bins_info=lambda wildcards: M.collections[wildcards.group]["bins_info"],
        collection_name=lambda wildcards: M.collections[wildcards.group][
            "collection_name"
        ],
        contigs_mode=lambda wildcards: M.collections[wildcards.group]["contigs_mode"],
    shell:
        """
           anvi-import-collection -p {input.profile} \
                                  -c {input.contigs} \
                                  -C {params.collection_name} \
                                  {params.bins_info} \
                                  {params.contigs_mode} \
                                  {input.collection_file} >> {log} 2>&1
        """


rule anvi_import_default_collection:
    input:
        profile=lambda wildcards: M.profile_databases[wildcards.group],
        contigs=ancient(M.get_contigs_db_path()),
    output:
        imported_collection_flag=touch(
            os.path.join(
                dirs_dict["MERGE_DIR"], "{group}", "default-collection-import.done"
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_import_default_collection.log",
    threads: 1
    resources:
        nodes=1,
    shell:
        """
           anvi-script-add-default-collection -p {input.profile} \
                                              -b {wildcards.group} \
                                              -c {input.contigs} >> {log} 2>&1
        """


rule anvi_summarize:
    input:
        profile=lambda wildcards: M.profile_databases[wildcards.group],
        contigs=ancient(M.get_contigs_db_path()),
        imported_collection_flag=lambda wildcards: M.get_collection_import_flag(
            str(wildcards.group)
        ),
    output:
        summary_done=temp(
            touch(os.path.join(dirs_dict["SUMMARY_DIR"], "{group}-summary.touch"))
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{group}-anvi_summarize.log"),
    threads: M.T("anvi_summarize")
    resources:
        nodes=M.T("anvi_summarize"),
    params:
        collection_name=lambda wildcards: M.collections[wildcards.group][
            "collection_name"
        ],
        summary_temp_dir=os.path.join(dirs_dict["SUMMARY_DIR"], "{group}-TEMP"),
        additional_params=M.get_param_value_from_config(
            ["anvi_summarize", "additional_params"]
        ),
    shell:
        """
            anvi-summarize -p {input.profile} \
                           -c {input.contigs} \
                           -C {params.collection_name} \
                           {params.additional_params} \
                           -o {params.summary_temp_dir} >> {log} 2>&1
        """


rule anvi_split:
    input:
        profile=lambda wildcards: M.profile_databases[wildcards.group],
        contigs=ancient(M.get_contigs_db_path()),
        imported_collection_flag=os.path.join(
            dirs_dict["MERGE_DIR"], "{group}", "collection-import.done"
        ),
    output:
        split_done=touch(
            os.path.join(dirs_dict["SPLIT_PROFILES_DIR"], "{group}-split.done")
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{group}-anvi_split.log"),
    threads: M.T("anvi_split")
    resources:
        nodes=M.T("anvi_split"),
    params:
        collection_name=lambda wildcards: M.collections[wildcards.group][
            "collection_name"
        ],
        split_dir=os.path.join(dirs_dict["SPLIT_PROFILES_DIR"], "{group}"),
        additional_params=M.get_param_value_from_config(
            ["anvi_summarize", "additional_params"]
        ),
    shell:
        """
            anvi-split -p {input.profile} \
                       -c {input.contigs} \
                       -C {params.collection_name} \
                       -o {params.split_dir} >> {log} 2>&1
        """


rule dummy_rule_for_anvi_summarize:
    input:
        summary_done=os.path.join(dirs_dict["SUMMARY_DIR"], "{group}-summary.touch"),
    output:
        summary=directory(os.path.join(dirs_dict["SUMMARY_DIR"], "{group}-SUMMARY")),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-dummy_rule_for_anvi_summarize.log",
    threads: 1
    resources:
        nodes=1,
    params:
        summary_temp_dir=os.path.join(dirs_dict["SUMMARY_DIR"], "{group}-TEMP"),
    shell:
        "mv {params.summary_temp_dir} {output.summary} >> {log} 2>&1"
