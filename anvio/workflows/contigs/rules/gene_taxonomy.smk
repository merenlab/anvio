rule export_gene_calls_for_centrifuge:
    """Export gene calls and use for centrifuge"""

    # marking the input as ancient in order to ignore timestamps.
    input:
        contigs=ancient(M.get_contigs_db_path()),
    # output is temporary. No need to keep this file.
    output:
        fasta=temp(dirs_dict["CONTIGS_DIR"] + "/{group}-gene-calls.fa"),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-export_gene_calls_for_centrifuge.log",
    threads: M.T("export_gene_calls_for_centrifuge")
    resources:
        nodes=M.T("export_gene_calls_for_centrifuge"),
    shell:
        "anvi-get-sequences-for-gene-calls -c {input} -o {output} >> {log} 2>&1"


rule centrifuge:
    """Run centrifuge on the exported gene calls of the contigs.db"""
    input:
        fasta=rules.export_gene_calls_for_centrifuge.output.fasta,
    output:
        hits=dirs_dict["CONTIGS_DIR"] + "/{group}-centrifuge_hits.tsv",
        report=dirs_dict["CONTIGS_DIR"] + "/{group}-centrifuge_report.tsv",
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-centrifuge.log",
    threads: M.T("centrifuge")
    resources:
        nodes=M.T("centrifuge"),
    params:
        db=M.get_param_value_from_config(["centrifuge", "db"]),
    shell:
        w.r("centrifuge -f \
                           -x {params.db} \
                           {input} \
                           -S {output.hits} \
                           --report-file {output.report} \
                           --threads {threads} >> {log} 2>&1")


rule anvi_import_taxonomy_for_genes:
    """Run anvi-import-taxonomy-for-genes"""
    input:
        hits=rules.centrifuge.output.hits,
        report=rules.centrifuge.output.report,
        # marking the contigs.db as ancient in order to ignore timestamps.
        contigs=ancient(M.get_contigs_db_path()),
    # using a flag file because no file is created by this rule.
    # for more information see:
    # http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#flag-files
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"],
                "{group}-steps",
                "anvi_anvi_import_taxonomy_for_genes.done",
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_import_taxonomy_for_genes.log",
    threads: M.T("anvi_import_taxonomy_for_genes")
    resources:
        nodes=M.T("anvi_import_taxonomy_for_genes"),
    params:
        parser="centrifuge",
    shell:
        w.r("anvi-import-taxonomy-for-genes -c {input.contigs} \
                                               -i {input.report} {input.hits} \
                                               -p {params.parser} >> {log} 2>&1")


rule anvi_run_hmms:
    """Run anvi-run-hmms"""

    # marking the input as ancient in order to ignore timestamps.
    input:
        contigs=ancient(M.get_contigs_db_path()),
    # using a snakemake flag file as an output since no file is generated
    # by the rule.
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_hmms.done"
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_run_hmms.log",
    threads: M.T("anvi_run_hmms")
    resources:
        nodes=M.T("anvi_run_hmms"),
    params:
        installed_hmm_profile=M.get_rule_param(
            "anvi_run_hmms", "--installed-hmm-profile"
        ),
        hmm_profile_dir=M.get_rule_param("anvi_run_hmms", "--hmm-profile-dir"),
        also_scan_trnas=M.get_rule_param("anvi_run_hmms", "--also-scan-trnas"),
        add_to_functions_tbl=M.get_rule_param(
            "anvi_run_hmms", "--add-to-functions-table"
        ),
        just_do_it=M.get_rule_param("anvi_run_hmms", "--just-do-it"),
    shell:
        w.r("anvi-run-hmms -c {input} \
                              -T {threads} \
                              {params.hmm_profile_dir} \
                              {params.installed_hmm_profile} \
                              {params.also_scan_trnas} \
                              {params.add_to_functions_tbl} \
                              {params.just_do_it} >> {log} 2>&1")


rule anvi_run_pfams:
    """Annotate contigs with Pfam protein family hits."""
    input:
        contigs=ancient(M.get_contigs_db_path()),
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_pfams.done"
            )
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{group}-anvi_run_pfams.log"),
    threads: M.T("anvi_run_pfams")
    resources:
        nodes=M.T("anvi_run_pfams"),
    params:
        pfam_data_dir=M.get_rule_param("anvi_run_pfams", "--pfam-data-dir"),
    shell:
        " anvi-run-pfams -c {input} {params.pfam_data_dir} -T {threads} >> {log} 2>&1"


rule anvi_run_kegg_kofams:
    """Annotate contigs with KEGG KOfam functional hits."""
    input:
        contigs=ancient(M.get_contigs_db_path()),
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_kegg_kofams.done"
            )
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{group}-anvi_run_kegg_kofams.log"),
    threads: M.T("anvi_run_kegg_kofams")
    resources:
        nodes=M.T("anvi_run_kegg_kofams"),
    params:
        kegg_data_dir=M.get_rule_param("anvi_run_kegg_kofams", "--kegg-data-dir"),
        hmmer_program=M.get_rule_param("anvi_run_kegg_kofams", "--hmmer-program"),
        keep_all_hits=M.get_rule_param("anvi_run_kegg_kofams", "--keep-all-hits"),
        log_bitscores=M.get_rule_param("anvi_run_kegg_kofams", "--log-bitscores"),
        just_do_it=M.get_rule_param("anvi_run_kegg_kofams", "--just-do-it"),
    shell:
        "anvi-run-kegg-kofams -c {input} {params.kegg_data_dir} {params.hmmer_program} {params.keep_all_hits} {params.log_bitscores} -T {threads} {params.just_do_it} >> {log} 2>&1"


rule anvi_run_ncbi_cogs:
    """Annotate contigs with NCBI COG functional hits."""
    input:
        contigs=ancient(M.get_contigs_db_path()),
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_ncbi_cogs.done"
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_run_ncbi_cogs.log",
    threads: M.T("anvi_run_ncbi_cogs")
    resources:
        nodes=M.T("anvi_run_ncbi_cogs"),
    params:
        # anvi-run-ncbi-cogs params. See anvi-run-ncbi-cogs help menu for more info.
        cog_data_dir=M.get_rule_param("anvi_run_ncbi_cogs", "--cog-data-dir"),
        temporary_dir_path=M.get_rule_param(
            "anvi_run_ncbi_cogs", "--temporary-dir-path"
        ),
        search_with=M.get_rule_param("anvi_run_ncbi_cogs", "--search-with"),
    shell:
        w.r("""anvi-run-ncbi-cogs -c {input.contigs} \
                                     -T {threads} \
                                     {params.cog_data_dir} \
                                     {params.temporary_dir_path} \
                                     {params.search_with} >> {log} 2>&1""")


rule anvi_run_scg_taxonomy:
    """Estimate taxonomy from single-copy core genes in the contigs database."""
    input:
        # make sure HMMs were done before running this rule
        hmms_done=ancient(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_hmms.done"
            )
        ),
        contigs=ancient(M.get_contigs_db_path()),
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_scg_taxonomy.done"
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_run_scg_taxonomy.log",
    threads: M.T("anvi_run_scg_taxonomy")
    resources:
        nodes=M.T("anvi_run_scg_taxonomy"),
    params:
        scg_taxonomy_data_dir=M.get_rule_param(
            "anvi_run_scg_taxonomy", "--scgs-taxonomy-data-dir"
        ),
    shell:
        w.r(
            """anvi-run-scg-taxonomy -c {input.contigs} \
                                        -T {threads} \
                                        {params.scg_taxonomy_data_dir} >> {log} 2>&1"""
        )


rule anvi_run_trna_scan:
    """Scan contigs for tRNA genes."""
    input:
        contigs=ancient(M.get_contigs_db_path()),
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"], "{group}-steps", "anvi_run_trna_scan.done"
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_run_trna_scan.log",
    threads: M.T("anvi_run_trna_scan")
    resources:
        nodes=M.T("anvi_run_trna_scan"),
    params:
        # FIXME not sure how to handle optional --trna-hits-file parameter
        #trna_scan_hits_file = M.get_rule_param('anvi_run_trna_scan', '--trna-hits-file'),
        trna_cutoff_score=M.get_rule_param("anvi_run_trna_scan", "--trna-cutoff-score"),
        trna_model=M.get_rule_param("anvi_run_trna_scan", "--trna-model"),
    shell:
        w.r("""anvi-scan-trnas -c {input.contigs} \
                                  -T {threads} \
                                  {params.trna_cutoff_score} \
                                  {params.trna_model} >> {log} 2>&1""")


rule anvi_get_sequences_for_gene_calls:
    """Export amino acid sequences for gene calls before external annotation."""
    input:
        contigs=ancient(M.get_contigs_db_path()),
    output:
        fasta=temp(dirs_dict["CONTIGS_DIR"] + "/{group}-contigs-aa-sequences.fa"),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_get_sequences_for_gene_calls.log",
    threads: M.T("anvi_get_sequences_for_gene_calls")
    resources:
        nodes=M.T("anvi_get_sequences_for_gene_calls"),
    shell:
        "anvi-get-sequences-for-gene-calls -c {input} -o {output} --get-aa-sequences >> {log} 2>&1"


rule emapper:
    """Run eggNOG-mapper on exported gene-call amino acid sequences."""
    input:
        fasta=rules.anvi_get_sequences_for_gene_calls.output.fasta,
    output:
        target=dirs_dict["CONTIGS_DIR"] + "/{group}-contigs.emapper.annotations",
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-emapper.log",
    threads: M.T("emapper")
    resources:
        nodes=M.T("emapper"),
    # TODO: add other emapper params
    params:
        contigs_path_without_extension=dirs_dict["CONTIGS_DIR"] + "/{group}-contigs",
        path_to_emapper_dir=M.get_param_value_from_config(
            ["emapper", "path_to_emapper_dir"]
        ),
        database=M.get_rule_param("emapper", "--database"),
        usemem=M.get_rule_param("emapper", "--usemem"),
        override=M.get_rule_param("emapper", "--override"),
    run:
        # running emapper
        shell(
            "python {params.path_to_emapper_dir}/emapper.py -i {input} --output {params.contigs_path_without_extension} "
            + "--cpu {threads} {params.database} {params.usemem} {params.override} >> {log} 2>&1"
        )


rule anvi_script_run_eggnog_mapper:
    """Import eggNOG-mapper annotations into the contigs database."""
    input:
        eggnog_output=rules.emapper.output.target,
        contigs=dirs_dict["CONTIGS_DIR"] + "/{group}.db",
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"],
                "{group}-steps",
                "anvi_script_run_eggnog_mapper.done",
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_script_run_eggnog_mapper.log",
    threads: M.T("anvi_script_run_eggnog_mapper")
    resources:
        nodes=M.T("anvi_script_run_eggnog_mapper"),
    params:
        use_version=M.get_rule_param("anvi_script_run_eggnog_mapper", "--use-version"),
    run:
        # Adding a 'g' prefix before every gene id (the anvi'o emapper driver requires this)
        shell(
            "sed 's/^[0-9]/g&/' {input.eggnog_output} > {input.eggnog_output}.temp 2>{log}"
        )
        shell(
            "anvi-script-run-eggnog-mapper -c {input.contigs} --annotation {input.eggnog_output}.temp {params.use_version}  >> {log} 2>&1"
        )
        shell("rm {input.eggnog_output}.temp 2>{log}")


# generate external genomes storage using `project_description` as name for contigs.
# this rule is only used when this workflow is inherited (e.g. by pangenomics or phylogenomics workflows)
