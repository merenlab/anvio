rule combine_sequence_data:
    """Cat all NT/AA fastas, reformat files, and renamed external-gene-calls files
Narrow end of the workflow funnel!

In this rule output files from metagenomes, genomes, SAGs, and MAGs all come together here!

cat all neccessary files:
    - AA sequences
    - NT sequences
    - header reformat files
    - external-gene-calls sequences and reformat_files from seperate metagenomes, genomes, SAGs, or MAGs into one fasta

"""
    input:
        files=lambda wildcards: M.get_input_files_combine_sequence_data(wildcards.group),
    output:
        NT_all=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-all.fna"
        ),
        AA_all=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-all.faa"
        ),
        reformat_report_all=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-reformat-report-all.txt",
        ),
        external_gene_calls_all=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-external_gene_calls_all.tsv",
        ),
        done=touch(
            os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-combine_sequence_data.done",
            )
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "combine_sequence_data_{group}.log"),
    threads: M.T("combine_sequence_data")
    run:
        # this code checks for existing hits per hmm per sample.
        # TODO: do the check at the previous step and provide list of sample for this rule
        # get list of unique hmm sources and hmm names from the group
        hmm_sources_name = []
        unique_source = []
        for hmm, value in M.hmm_dict.items():
            if value["group"] == wildcards.group:
                hmm_sources_name.append((value["source"], value["name"]))
        # list of hmm_hits file
        hmm_hits = [
            os.path.join(
                dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                sample_name,
                f"{hmm_source}-dom-hmmsearch",
                "hmm_hits_filtered.txt",
            )
            for sample_name in M.names_list
            for hmm_source in unique_source
        ]
        # list for various paths, to merge later
        NT_list = []
        AA_list = []
        AA_reformat_list = []
        external_gene_calls_reformat_list = []
        # check for hmm hit per sample
        for hmm_source, hmm_name in hmm_sources_name:
            contigs_db_with_hmm_hits = []
            contigs_db_with_hmm_NO_hits = []
            for sample_name in M.names_list:
                hmm_hit = os.path.join(
                    dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                    sample_name,
                    f"{hmm_source}-dom-hmmsearch",
                    "hmm_hits_filtered.txt",
                )
                df = pd.read_csv(hmm_hit, sep="\t")
                df = df[df.source == hmm_source]
                gene_name_list = df["gene_name"].tolist()
                if hmm_name in gene_name_list:
                    contigs_db_with_hmm_hits.append(sample_name)
                    # add path to merge
                    working_dir = os.path.join(
                        dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                        f"{sample_name}",
                        f"{hmm_source}",
                        f"{hmm_name}",
                    )
                    NT_path = os.path.join(
                        working_dir,
                        f"{sample_name}-{hmm_name}-hmm_hits_renamed.fna",
                    )
                    AA_path = os.path.join(
                        working_dir,
                        f"{sample_name}-{hmm_name}-hmm_hits_renamed.faa",
                    )
                    NT_reformat_path = os.path.join(
                        working_dir,
                        f"{sample_name}-{hmm_name}-reformat_report_NT.txt",
                    )
                    AA_reformat_path = os.path.join(
                        working_dir,
                        f"{sample_name}-{hmm_name}-reformat_report_AA.txt",
                    )
                    external_gene_calls_reformat_reformat_path = os.path.join(
                        working_dir,
                        f"{sample_name}-{hmm_name}-external_gene_calls_renamed.tsv",
                    )
                    NT_list.append(NT_path)
                    AA_list.append(AA_path)
                    AA_reformat_list.append(AA_reformat_path)
                    external_gene_calls_reformat_list.append(
                        external_gene_calls_reformat_reformat_path
                    )
                else:
                    contigs_db_with_hmm_NO_hits.append(sample_name)
            outfile = os.path.join(
                dirs_dict["LOGS_DIR"],
                f"contigDBs_with_no_hmm_hit_{hmm_source}-{hmm_name}.log",
            )
            with open(outfile, "a") as outfile:
                for element in contigs_db_with_hmm_NO_hits:
                    outfile.write(element + "\n")
        # time to merge all these files
        with open(output.NT_all, "a") as NT_output:
            for f in NT_list:
                NT_output.write(open(f).read())
        with open(output.AA_all, "a") as AA_all:
            for f in AA_list:
                AA_all.write(open(f).read())
        with open(output.reformat_report_all, "a") as reformat_report_all:
            for f in AA_reformat_list:
                reformat_report_all.write(open(f).read())
        col_names = [
            "gene_callers_id",
            "contig",
            "start",
            "stop",
            "direction",
            "partial",
            "call_type",
            "source",
            "version",
            "aa_sequence",
        ]
        with open(output.external_gene_calls_all, "a") as external_gene_calls_all:
            external_gene_calls_all.write("\t".join(col_names) + "\n")
            for f in external_gene_calls_reformat_list:
                file = open(f).read().splitlines(True)
                external_gene_calls_all.writelines(file[1:])
