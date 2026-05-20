def extract_misc_data(mmseqs_cluster_rep_index, final_sequences_headers, output):
    """This function creates a tsv of metadata for mmseqs cluster tsv outfile (*_cluster.tsv)

    Parameters
    ==========
    mmseqs_cluster_rep_index: str
        path to tsv containing mmseqs cluster representatives and member index

    final_sequences_headers: str
        path to tsv with fasta file headers

    output:
        target=str
        path to output tsv, metadata values include ['split_name', 'contigs_db_type', 'genomic_seq_in_cluster', 'cluster_size']
        - 'split_name': ID name in interactive interface
        - 'contigs_db_type': cluster representative origin: MAG, SAG, isolate genome, metagenome.
        - 'genomic_seq_in_cluster': detects if a cluster member came from an external-genome
        - 'cluster_size': number of sequences in cluster
    """

    # Import data
    # ------------
    cluster_rep_index = pd.read_csv(
        mmseqs_cluster_rep_index,
        sep="\t",
        index_col=False,
        names=["representative", "cluster_members"],
    )

    final_sequences_headers = pd.read_csv(
        final_sequences_headers, sep="\t", index_col=False, names=["identifier"]
    )

    # Clean data
    # ------------

    # Detect if there is a genomic reference protein in cluster
    cluster_rep_index_dict = (
        cluster_rep_index.groupby("representative")["cluster_members"]
        .apply(list)
        .to_dict()
    )

    cluster_reps_with_genomic_references_list = []
    for seq in final_sequences_headers.iloc[:, 0].tolist():
        cluster_members_list = cluster_rep_index_dict[seq]
        for external_genome in M.external_genomes_names_list:
            check = any(external_genome in s for s in cluster_members_list)
            if check is True:
                cluster_reps_with_genomic_references_list.append(seq)

    # Count size of clusters
    df = cluster_rep_index.groupby("representative").apply(count_cluster_size)[
        ["cluster_members", "cluster_size"]
    ]

    # Make split names for anvi-interactive
    df["split_name"] = df["cluster_members"].astype(str) + "_split_00001"

    # subset misc data to final set of proteins
    df = pd.merge(
        df,
        final_sequences_headers,
        left_on="cluster_members",
        right_on="identifier",
        how="inner",
    )

    df["genomic_seq_in_cluster"] = np.where(
        df["cluster_members"].isin(cluster_reps_with_genomic_references_list),
        "yes",
        "no",
    )

    # Determine contigs_db type: metagenome or genomes
    # FIXME: This will need to be changed in the future to accomidate SAGs, MAGs, and other genomic sources
    # If metagenome_name is in external_genomes_names_list then it's a genome
    contigs_db_type_dict = {}
    for name in list(df.cluster_members):
        if any(x in name for x in M.external_genomes_names_list):
            contigs_db_type_dict[name] = "genome"
        else:
            contigs_db_type_dict[name] = "metagenome"

    df["contigs_db_type"] = df.cluster_members.apply(
        lambda name: contigs_db_type_dict[name]
    )

    # grab the final columns
    df = df[["split_name", "contigs_db_type", "genomic_seq_in_cluster", "cluster_size"]]

    # Export
    # -------
    df.to_csv(output, sep="\t", index=None, na_rep="NA")


def count_cluster_size(group):
    c = group["cluster_members"].count()
    group["cluster_size"] = c

    return group


rule make_misc_data:
    """Make misc data file for clustered sequences"""
    input:
        final_list_of_sequences_for_mapping_headers=os.path.join(
            dirs_dict["MSA"], "{group}", "{group}_headers.tmp"
        ),
    output:
        misc_data_final=os.path.join(
            dirs_dict["MISC_DATA"], "{group}", "{group}_misc.tsv"
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "add_contigs_db_type_{group}.log"),
    threads: M.T("add_misc_data_to_taxonomy")
    params:
        mmseqs_cluster_rep_index=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-mmseqs_NR_cluster.tsv",
        ),
        coverage_cluster_rep_index=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-coverage_cluster.tsv",
        ),
    run:
        """Here we determine the origin of each SCG (which kind of contigs_db): metagenome, isolate genome, etc."""
        if M.cluster_representative_method == "mmseqs":
            extract_misc_data(
                mmseqs_cluster_rep_index=params.mmseqs_cluster_rep_index,
                final_sequences_headers=input.final_list_of_sequences_for_mapping_headers,
                output=output.misc_data_final,
            )
        if M.cluster_representative_method == "cluster_rep_with_coverages":
            extract_misc_data(
                mmseqs_cluster_rep_index=params.coverage_cluster_rep_index,
                final_sequences_headers=input.final_list_of_sequences_for_mapping_headers,
                output=output.misc_data_final,
            )


rule anvi_run_scg_taxonomy:
    """Run anvi-run-scg-taxonomy"""
    input:
        targets=expand(
            os.path.join(
                dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                "{{sample_name}}",
                "{hmm_source}-dom-hmmsearch",
                "contigs-hmmsearch.done",
            ),
            hmm_source=M.unique_hmm_source.keys(),
        ),
    output:
        done=touch(
            os.path.join(
                dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                "{sample_name}",
                "{sample_name}_scg_taxonomy.done",
            )
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "anvi_run_scg_taxonomy-{sample_name}.log"),
    threads: M.T("anvi_run_scg_taxonomy")
    params:
        scg_taxonomy_version=M.get_param_value_from_config(
            ["scg_taxonomy_database_version"]
        ),
        additional_params=M.get_param_value_from_config(
            ["anvi_run_scg_taxonomy", "additional_params"]
        ),
    run:
        # check if anvi-run-scg was already run
        contigs_db = M.contigs_db_name_path_dict[wildcards.sample_name]
        print(contigs_db)
        with db.DB(contigs_db, None, ignore_version=True) as database:
            if (
                database.get_meta_value("scg_taxonomy_database_version")
                != params.scg_taxonomy_version
            ):
                print(f"Running anvi-run-scg-taxonomy on {contigs_db}")
                shell(
                    "anvi-run-scg-taxonomy -c {contigs_db} --num-threads {threads} {params.additional_params} > {log} 2>&1"
                )


rule anvi_estimate_scg_taxonomy:
    """Run anvi-estimate-SCG-taxonomy and import the resulting taxonomy misc data to profileDB for internal hmms only"""
    input:
        final_list_of_sequences_for_mapping_headers=os.path.join(
            dirs_dict["MSA"], "{group}", "{group}_headers.tmp"
        ),
    output:
        done=touch(
            os.path.join(
                dirs_dict["MISC_DATA"],
                "{group}",
                "anvi_estimate_scg_taxonomy_for_SCGs.done",
            )
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "anvi_scg_taxonomy_{group}.log"),
    threads: M.T("anvi_estimate_scg_taxonomy")
    params:
        reformat_file=rules.combine_sequence_data.output.reformat_report_all,
        misc_data_dir=os.path.join(dirs_dict["MISC_DATA"], "{group}"),
        taxonomy=os.path.join(
            dirs_dict["MISC_DATA"], "{group}", "{group}_estimate_scg_taxonomy_results"
        ),
        taxonomy_long=os.path.join(
            dirs_dict["MISC_DATA"],
            "{group}",
            "{group}_estimate_scg_taxonomy_results-RAW-LONG-FORMAT.txt",
        ),
        tax_data_final=os.path.join(
            dirs_dict["MISC_DATA"], "{group}", "{group}_scg_taxonomy_data.tsv"
        ),
    run:
        # get hmm name, if more than one, exit for now
        hmm_name, hmm_source = "", ""
        for hmm, value in M.hmm_dict.items():
            if value["group"] == wildcards.group:
                hmm_name = value["name"]
                hmm_source = value["source"]
        # Concatenate metagenomes.txt and external-genomes.txt
        contigs_db_name_path_list = list(M.contigs_db_name_path_dict.items())
        contigs_db_name_path_df = pd.DataFrame(
            contigs_db_name_path_list, columns=["name", "contigs_db_path"]
        )
        combined_genomes_df_path = f"{wildcards.group}_combined_genomes.txt"
        contigs_db_name_path_df.to_csv(
            combined_genomes_df_path, sep="\t", index=False, header=True
        )
        shell("mkdir -p {params.misc_data_dir}; \
        anvi-estimate-scg-taxonomy -M {combined_genomes_df_path} \
                                   --metagenome-mode \
                                   --scg-name-for-metagenome-mode {hmm_name} \
                                   -T {threads} \
                                   --raw-output \
                                   -O {params.taxonomy} &> {log}")
        # Import data
        # ------------
        scg_taxonomy = pd.read_csv(params.taxonomy_long, sep="\t", index_col=False)
        reformat_report = pd.read_csv(
            params.reformat_file,
            sep="\t",
            index_col=False,
            names=["new_header", "header"],
        )
        final_sequences_headers = pd.read_csv(
            input.final_list_of_sequences_for_mapping_headers,
            sep="\t",
            index_col=False,
            names=["identifier"],
        )
        # Clean Data
        # -----------
        reformat_report = pd.merge(
            reformat_report,
            final_sequences_headers,
            left_on="new_header",
            right_on="identifier",
            how="inner",
        )
        reformat_report["gene_callers_id"] = (
            reformat_report["header"]
            .str.split("gene_callers_id:|\|start:", expand=True)[1]
            .astype(str)
        )
        reformat_report["new_header_tmp"] = (
            reformat_report["new_header"].str.rsplit("_", 1).str[0]
        )
        reformat_report["identifier"] = (
            reformat_report["new_header_tmp"]
            + "_"
            + reformat_report["gene_callers_id"]
        )
        scg_taxonomy["identifier"] = (
            scg_taxonomy["metagenome_name"]
            + "_"
            + hmm_source
            + "_"
            + scg_taxonomy["gene_name"]
            + "_"
            + scg_taxonomy["gene_callers_id"].astype(str)
        )
        scg_taxonomy = scg_taxonomy.merge(
            reformat_report, on="identifier", how="inner"
        )
        scg_taxonomy["split_name"] = (
            scg_taxonomy["new_header"].astype(str) + "_split_00001"
        )
        scg_taxonomy = scg_taxonomy[
            [
                "split_name",
                "identifier",
                "percent_identity",
                "t_domain",
                "t_phylum",
                "t_class",
                "t_order",
                "t_family",
                "t_genus",
                "t_species",
            ]
        ]
        # Export
        # -------
        scg_taxonomy.to_csv(
            params.tax_data_final, sep="\t", index=None, na_rep="NA"
        )
