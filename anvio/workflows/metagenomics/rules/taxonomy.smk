rule krakenuniq:
    """Classify reads taxonomically with KrakenUniq."""
    input:
        unpack(lambda wildcards: M.get_fastq(wildcards.readset)),
    output:
        raw=os.path.join(dirs_dict["TAXONOMY_DIR"], "{readset}.krakenuniq"),
        tsv=os.path.join(dirs_dict["TAXONOMY_DIR"], "{readset}-krakenuniq-report.tsv"),
    log:
        rule_log("krakenuniq", "{readset}-krakenuniq"),
    wildcard_constraints:
        readset=SR_RS_RE,
    threads: M.T("krakenuniq")
    resources:
        nodes=M.T("krakenuniq"),
    params:
        db=M.get_rule_param("krakenuniq", "--db"),
        gzip_compressed=M.get_rule_param("krakenuniq", "--gzip-compressed"),
        additional_params=M.get_param_value_from_config(
            ["krakenuniq", "additional_params"]
        ),
    shell:
        """
            krakenuniq \
                     {params.db} \
                     --threads {threads} \
                     --fastq-input \
                     {params.gzip_compressed} \
                     --paired \
                     --output {output.raw} \
                     --report-file {output.tsv} \
                     {input.r1} \
                     {input.r2} \
                     {params.additional_params} >> {log} 2>&1
        """


rule krakenuniq_mpa_report:
    """configure input for import_percent_of_reads_mapped."""
    input:
        source=os.path.join(dirs_dict["TAXONOMY_DIR"], "{readset}.krakenuniq"),
    output:
        txt=os.path.join(dirs_dict["TAXONOMY_DIR"], "{readset}-krakenuniq.tsv"),
    log:
        rule_log("krakenuniq_mpa_report", "{readset}-krakenuniq_mpa_report"),
    threads: M.T("krakenuniq_mpa_report")
    resources:
        nodes=M.T("krakenuniq_mpa_report"),
    # FIXME: I need to see if more params should be available for krakenuniq_mpa_report
    params:
        db=M.get_rule_param("krakenuniq", "--db"),
    shell:
        """
            krakenuniq-mpa-report {params.db} \
                                  --header-line \
                                  {input} \
                                  > {output} 2>> {log}
        """


def input_for_import_krakenuniq_taxonomy(wildcards):
    """Select KrakenUniq annotations and optional mapping statistics for taxonomy import."""
    D = {}
    if M.kraken_annotation_dict:
        # if the user supplied a file with kraken annotations then use it
        D["tsv"] = M.kraken_annotation_dict[wildcards.readset]["path"]
    else:
        # if the user didn't provide a file, then annotation needs to run
        D["tsv"] = os.path.join(
            dirs_dict["TAXONOMY_DIR"], wildcards.readset + "-krakenuniq.tsv"
        )

    if run_import_percent_of_reads_mapped:
        # we include percent_of_reads_mapped_imported_flag as input here so that kraken and import_percent_of_reads_mapped will never run in parallel.
        # If they run in parallel it could cause trouble as they will try to write to the profile database at the same time. very bad.
        D["percent_of_reads_mapped_imported_flag"] = (
            os.path.join(
                dirs_dict["PROFILE_DIR"],
                wildcards.group,
                wildcards.readset,
                "layers-additional-data.txt",
            ),
        )
    return D


rule import_krakenuniq_taxonomy:
    """Import KrakenUniq taxonomy assignments into anvi-o profile layers."""
    input:
        unpack(input_for_import_krakenuniq_taxonomy),
        profile=ancient(
            os.path.join(
                dirs_dict["PROFILE_DIR"], "{group}", "{readset}", "PROFILE.db"
            )
        ),
    output:
        done=touch(
            os.path.join(
                dirs_dict["PROFILE_DIR"],
                "{group}/{readset}/import_krakenuniq_taxonomy.done",
            )
        ),
    log:
        rule_log(
            "import_krakenuniq_taxonomy",
            "{group}-{readset}-import_krakenuniq_taxonomy",
        ),
    threads: M.T("import_krakenuniq_taxonomy")
    resources:
        nodes=M.T("import_krakenuniq_taxonomy"),
    params:
        min_abundance=M.get_rule_param("import_krakenuniq_taxonomy", "--min-abundance"),
    shell:
        """
            anvi-import-taxonomy-for-layers -p {input.profile} \
                                            -i {input.tsv} \
                                            --parser krakenuniq \
                                            {params.min_abundance} >> {log} 2>&1
        """
