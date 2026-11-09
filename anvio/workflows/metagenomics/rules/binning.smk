rule anvi_cluster_contigs:
    """Cluster contigs into bins with the configured anvi-cluster-contigs driver."""
    input:
        contigs=ancient(M.get_contigs_db_path()),
        profile=lambda wildcards: M.profile_databases[wildcards.group],
    output:
        done=touch(
            os.path.join(
                dirs_dict["MERGE_DIR"], "{group}" + "-" + "{driver}" + ".done"
            )
        ),
    log:
        rule_log("anvi_cluster_contigs", "{group}-{driver}-anvi_cluster_contigs"),
    threads: M.T("anvi_cluster_contigs")
    resources:
        nodes=M.T("anvi_cluster_contigs"),
    params:
        driver="{driver}",
        collection_name=M.get_rule_param("anvi_cluster_contigs", "--collection-name"),
        just_do_it=M.get_rule_param("anvi_cluster_contigs", "--just-do-it"),
        additional_params=lambda wildcards: M.get_param_value_from_config(
            [
                "anvi_cluster_contigs",
                M.get_param_name_for_binning_driver(wildcards.driver),
            ]
        ),
    shell:
        w.r("""anvi-cluster-contigs -c {input.contigs} \
                                       -p {input.profile} \
                                       --driver {params.driver} \
                                       {params.collection_name} \
                                       {params.just_do_it} \
                                       {params.additional_params} >> {log} 2>&1""")
