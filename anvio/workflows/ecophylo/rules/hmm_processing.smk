def get_hmm_threads(wildcards):
    """This function conditionally selects threads for anvi-run-hmms based on
    if a contigs-db is a metagenome or not.

    Parameters
    ==========
    wildcards: snakemake object
        allows you to programmatically access snakemake wildcards

    Returns
    =======
    threads : int
        number of threads to be used
    """
    threads = M.T("anvi_run_hmms_hmmsearch")
    max_threads = M.get_param_value_from_config("max_threads")
    if not max_threads:
        max_threads = float("Inf")

    if M.metagenomes:
        if wildcards.sample_name in M.metagenomes_name_list:
            threads = M.get_param_value_from_config(
                ["anvi_run_hmms_hmmsearch", "threads_metagenomes"]
            )
    else:
        threads = M.get_param_value_from_config(
            ["anvi_run_hmms_hmmsearch", "threads_genomes"]
        )

    if threads:
        try:
            if int(threads) > float(max_threads):
                return int(max_threads)
            else:
                return int(threads)
        except:
            raise ConfigError(
                f'"threads" must be an integer number. In your config file you provided "{threads}" for '
                f"the number of threads for rule anvi_run_hmms_hmmsearch"
            )
    else:
        return 1


rule anvi_run_hmms_hmmsearch:
    """Run hmmsearch with input hmms to get domtblout"""
    output:
        done=os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "contigs-hmmsearch.done",
        ),
        hmm_hits=os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "hmm_hits.txt",
        ),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            "anvi_run_hmms_hmmsearch-{sample_name}-{hmm_source}.log",
        ),
    threads: get_hmm_threads
    params:
        additional_params=M.get_param_value_from_config(
            ["anvi_run_hmms_hmmsearch", "additional_params"]
        ),
    run:
        contigs_db_path = os.path.join(
            M.contigs_db_name_path_dict[wildcards.sample_name]
        )
        hmm_source = wildcards.hmm_source
        hmm_dir = os.path.join(M.unique_hmm_source[hmm_source])
        hmmer_output_dir = os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            f"{wildcards.sample_name}",
            f"{hmm_source}-dom-hmmsearch",
        )
        if M.metagenomes:
            if wildcards.sample_name in M.metagenomes_name_list:
                threads = (
                    M.get_param_value_from_config(
                        ["anvi_run_hmms_hmmsearch", "threads_metagenomes"]
                    ),
                )
                threads = int(threads[0])
        else:
            threads = (
                M.get_param_value_from_config(
                    ["anvi_run_hmms_hmmsearch", "threads_genomes"]
                ),
            )
        domtblout = os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            f"{wildcards.sample_name}-{hmm_source}",
            "dom-hmmsearch",
            "hmm.domtable",
        )
        # Run different hmm search depending on whether a hmm is internal or external because anvio
        if hmm_source in M.internal_hmm_sources:
            if not os.path.exists(domtblout):
                print(f"Running internal hmm dataset: {hmm_source}")
                shell("anvi-run-hmms -c {contigs_db_path} \
                                     --hmmer-program hmmsearch \
                                     --hmmer-output-dir {hmmer_output_dir} \
                                     --installed-hmm-profile {hmm_source} \
                                     --domain-hits-table \
                                     --just-do-it \
                                     -T {threads} >> {log} 2>&1")
        else:
            if not os.path.exists(domtblout):
                print(f"Running external hmm dataset: {hmm_source}")
                shell("anvi-run-hmms -c {contigs_db_path} \
                                     --hmmer-program hmmsearch \
                                     --hmm-profile-dir {hmm_dir} \
                                     --hmmer-output-dir {hmmer_output_dir} \
                                     --domain-hits-table \
                                     --just-do-it \
                                     -T {threads} >> {log} 2>&1")
        # Get hmm_hits.txt
        get_hmm_hits_txt(contigs_db_path, output.hmm_hits)
        shell("touch {output.done}")


def get_hmm_hits_txt(contigs_db, out_file):
    """This function extracts the hmm_hits table from a contigs-db

    Parameters
    ==========
    contigs_db

    out_file : str
        Output file path for hmm_hits.tsv

    Returns
    =======
    hmm_hits.tsv : tsv
    """

    database = db.DB(contigs_db, None, ignore_version=True)
    tables_in_database = database.get_table_names()
    args_table = "hmm_hits"

    if args_table not in tables_in_database:
        # Make empty hmm_hits.txt so that the rule the next rule has something to work with
        column_names = [
            "entry_id",
            "source",
            "gene_unique_identifier",
            "gene_callers_id",
            "gene_name",
            "gene_hmm_id",
            "e_value",
        ]
        df = pd.DataFrame(columns=column_names)
        df.to_csv(out_file, sep="\t", index=False, header=True)
    else:
        table_columns = database.get_table_structure(args_table)
        table_content = database.get_table_as_dataframe(
            args_table, columns_of_interest=table_columns, error_if_no_data=False
        )
        u.store_dataframe_as_TAB_delimited_file(table_content, out_file)


rule filter_hmm_hits_by_model_coverage:
    """Remove weak hmm hits using model alignment coverage filter"""
    input:
        #done = ancient(rules.anvi_run_hmms_hmmsearch.output.done),
        hmm_hits=ancient(rules.anvi_run_hmms_hmmsearch.output.hmm_hits),
    output:
        done=os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "{sample_name}-{hmm_source}-DB_filtered.done",
        ),
        hmm_hits_filtered=os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "hmm_hits_filtered.txt",
        ),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            "filter_hmm_hits_by_model_coverage-{sample_name}-{hmm_source}.log",
        ),
    threads: M.T("filter_hmm_hits_by_model_coverage")
    params:
        model_coverage=M.get_param_value_from_config(
            ["filter_hmm_hits_by_model_coverage", "--min-model-coverage"]
        ),
        partial_ORFs=M.get_rule_param(
            "filter_hmm_hits_by_model_coverage", "--filter-out-partial-gene-calls"
        ),
        additional_params=M.get_param_value_from_config(
            ["filter_hmm_hits_by_model_coverage", "additional_params"]
        ),
    run:
        contigs_db = os.path.join(
            M.contigs_db_name_path_dict[wildcards.sample_name]
        )
        hmm_source = wildcards.hmm_source
        hmm_dir = os.path.join(M.unique_hmm_source[hmm_source])
        domtblout = os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            f"{wildcards.sample_name}",
            f"{hmm_source}-dom-hmmsearch",
            "hmm.domtable",
        )
        # import hmm_hits to find out how many hmm_hits we got from the hmm model
        # if we don't have any then we can skippp all dis
        df = pd.read_csv(input.hmm_hits, sep="\t")
        df = df[df.source == hmm_source]
        gene_name_list = df["gene_name"].tolist()
        # There is no point in running the following if we don't have the target hmm in the contigs_db!
        # TODO: better catch missing hmm hits
        if gene_name_list:
            print("we have a hit!")
            if hmm_source in M.internal_hmm_sources:
                shell(
                    "anvi-script-filter-hmm-hits-table -c {contigs_db} \
                                                        --domain-hits-table {domtblout} \
                                                        --hmm-source {hmm_source} \
                                                        --min-model-coverage {params.model_coverage} \
                                                        {params.partial_ORFs} \
                                                        {params.additional_params} >> {log} 2>&1"
                )
            else:
                shell(
                    "anvi-script-filter-hmm-hits-table -c {contigs_db} \
                                                        --domain-hits-table {domtblout} \
                                                        --hmm-profile-dir {hmm_dir} \
                                                        --hmm-source {hmm_source} \
                                                        --min-model-coverage {params.model_coverage} \
                                                        {params.partial_ORFs} \
                                                        {params.additional_params} >> {log} 2>&1"
                )
        else:
            no_hmm_hits_string = f"The hmm source {hmm_source} was not found in the hmm_hits of the contigs_db: {contigs_db}"
            print(no_hmm_hits_string)
            with open(log, "w") as logfile:
                logfile.write(no_hmm_hits_string)
        # Get hmm_hits.txt
        get_hmm_hits_txt(contigs_db, output.hmm_hits_filtered)
        shell("touch {output.done}")


rule process_hmm_hits:
    """Extract AA/NT fastas and external-gene-calls from contigs-db then rename external-gene-calls
names with reformated names.
"""
    input:
        done=ancient(rules.filter_hmm_hits_by_model_coverage.output.done),
    output:
        done=os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            "{sample_name}",
            "{hmm_source}",
            "{hmm_name}",
            "{sample_name}-{hmm_name}-processed.done",
        ),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            "process_hmm_hits-{sample_name}-{hmm_source}-{hmm_name}.log",
        ),
    threads: M.T("process_hmm_hits")
    params:
        hmm_hits=os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "hmm_hits_filtered.txt",
        ),
    run:
        contigs_db = os.path.join(
            M.contigs_db_name_path_dict[wildcards.sample_name]
        )
        hmm_source = wildcards.hmm_source
        fasta_output_dir = os.path.join(
            dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
            f"{wildcards.sample_name}",
            f"{wildcards.hmm_source}",
            f"{wildcards.hmm_name}",
        )
        if not os.path.exists(fasta_output_dir):
            os.mkdir(fasta_output_dir)
        faa = os.path.join(
            fasta_output_dir,
            f"{wildcards.sample_name}-{wildcards.hmm_name}-hmm_hits.faa",
        )
        fna = os.path.join(
            fasta_output_dir,
            f"{wildcards.sample_name}-{wildcards.hmm_name}-hmm_hits.fna",
        )
        # anvi-script-reformat-fasta output files
        fasta_NT = (
            os.path.join(
                fasta_output_dir,
                f"{wildcards.sample_name}-{wildcards.hmm_name}-hmm_hits_renamed.fna",
            ),
        )
        fasta_AA = (
            os.path.join(
                fasta_output_dir,
                f"{wildcards.sample_name}-{wildcards.hmm_name}-hmm_hits_renamed.faa",
            ),
        )
        report_file_NT = (
            os.path.join(
                fasta_output_dir,
                f"{wildcards.sample_name}-{wildcards.hmm_name}-reformat_report_nt.txt",
            ),
        )
        report_file_AA = os.path.join(
            fasta_output_dir,
            f"{wildcards.sample_name}-{wildcards.hmm_name}-reformat_report_AA.txt",
        )
        # anvi-get-sequences-for-gene-calls input files
        external_gene_calls = os.path.join(
            fasta_output_dir,
            f"{wildcards.sample_name}-{wildcards.hmm_name}-external_gene_calls.tsv",
        )
        fasta = os.path.join(
            fasta_output_dir,
            f"{wildcards.sample_name}-{wildcards.hmm_name}-orfs.fna",
        )
        # rename external gene calls
        external_gene_calls_renamed = os.path.join(
            fasta_output_dir,
            f"{wildcards.sample_name}-{wildcards.hmm_name}-external_gene_calls_renamed.tsv",
        )
        # import hmm_hits to find out how many hmm_hits we got from the hmm model
        df = pd.read_csv(params.hmm_hits, sep="\t")
        df = df[df.source == hmm_source]
        gene_name_list = df["gene_name"].tolist()
        # There is no point in running the following if we don't have any hmm_hits to begin with!
        if wildcards.hmm_name in gene_name_list:
            if hmm_source in M.internal_hmm_sources:
                # Get AA and NT fasta files
                shell(
                    "echo -e 'anvi-get-sequences-for-hmm-hits AA stdout:\n' >> {log}"
                )
                shell("echo -e '' >> {log}")
                shell(
                    "anvi-get-sequences-for-hmm-hits -c {contigs_db} \
                                                       --hmm-sources  {hmm_source} \
                                                       --gene-names {wildcards.hmm_name} \
                                                       --get-aa-sequences \
                                                       -o {faa} --just-do-it >> {log} 2>&1"
                )
                shell("echo -e '' >> {log}")
                shell(
                    "echo -e 'anvi-get-sequences-for-hmm-hits NT stdout:\n' >> {log}"
                )
                shell("echo -e '' >> {log}")
                shell(
                    "anvi-get-sequences-for-hmm-hits -c {contigs_db} \
                                                       --hmm-sources  {hmm_source} \
                                                       --gene-names {wildcards.hmm_name} \
                                                       -o {fna} --just-do-it >> {log} 2>&1"
                )
                shell("echo -e '' >> {log}")
            else:
                # Get AA and NT fasta files
                shell(
                    "echo -e 'anvi-get-sequences-for-hmm-hits AA stdout:\n' >> {log}"
                )
                shell(
                    "anvi-get-sequences-for-hmm-hits -c {contigs_db} \
                                                        --hmm-sources  {hmm_source} \
                                                        --get-aa-sequences \
                                                        -o {faa} --just-do-it >> {log} 2>&1"
                )
                shell("echo -e '' >> {log}")
                shell(
                    "echo -e 'anvi-get-sequences-for-hmm-hits NT stdout:\n' >> {log}"
                )
                shell(
                    "anvi-get-sequences-for-hmm-hits -c {contigs_db} \
                                                        --hmm-sources  {hmm_source} \
                                                        -o {fna} --just-do-it >> {log} 2>&1"
                )
            # Clean up fasta headers for tree calculation
            prefix_list = [
                wildcards.sample_name,
                wildcards.hmm_source,
                wildcards.hmm_name,
            ]
            prefix = "_".join(prefix_list)
            shell("echo -e 'anvi-script-reformat-fasta AA stdout:\n' >> {log}")
            shell("anvi-script-reformat-fasta {faa} \
                                    --simplify-names \
                                    --prefix {prefix} \
                                    --report-file {report_file_AA} \
                                    -o {fasta_AA} >> {log} 2>&1")
            shell("echo -e '' >> {log}")
            shell("echo -e 'anvi-script-reformat-fasta NT stdout:\n' >> {log}")
            shell("anvi-script-reformat-fasta {fna} \
                                    --simplify-names \
                                    --prefix {prefix} \
                                    --report-file {report_file_NT} \
                                    -o {fasta_NT} >> {log} 2>&1")
            shell("echo -e 'anvi-get-sequences-for-gene-calls stdout:\n' >> {log}")
            shell(
                "anvi-get-sequences-for-gene-calls -c {contigs_db} \
                                                     --external-gene-calls {external_gene_calls} \
                                                     -o {fasta} >> {log} 2>&1"
            )
            # Import tables
            external_gene_calls = pd.read_csv(
                external_gene_calls, delim_whitespace=True, index_col=False
            )
            reformat_file = os.path.join(
                fasta_output_dir,
                f"{wildcards.sample_name}-{wildcards.hmm_name}-reformat_report_AA.txt",
            )
            reformat_report = pd.read_csv(
                reformat_file,
                sep="\t",
                index_col=False,
                names=["new_header", "header"],
            )
            # Parse input files
            # -----------------
            # Parse reformat_report to get gene-callers-ids
            reformat_report["gene_callers_id"] = (
                reformat_report["header"]
                .str.split("gene_callers_id:|\|start:", expand=True)[1]
                .astype(str)
            )
            # Parse external-gene-calls contig column to get gene-callers-ids
            external_gene_calls[["name", "contig_number", "gene_callers_id"]] = (
                external_gene_calls["contig"].str.rsplit("_", 2, expand=True)
            )
            # Join external-gene-calls with reformat-report on gene-callers-id
            # -----------------
            external_gene_calls = external_gene_calls.merge(
                reformat_report, on="gene_callers_id", how="inner"
            )
            # Replace contig name from anvi-export-gene-calls with the simplified, new header from anvi-script-reformat-fasta
            external_gene_calls = external_gene_calls[
                [
                    "gene_callers_id",
                    "new_header",
                    "start",
                    "stop",
                    "direction",
                    "partial",
                    "call_type",
                    "source",
                    "version",
                    "aa_sequence",
                ]
            ]
            external_gene_calls = external_gene_calls.rename(
                columns={"new_header": "contig"}
            )
            # Write file
            external_gene_calls.to_csv(
                external_gene_calls_renamed, sep="\t", index=False, header=True
            )
            shell("touch {output.done}")
        else:
            shell("touch {output.done}")
