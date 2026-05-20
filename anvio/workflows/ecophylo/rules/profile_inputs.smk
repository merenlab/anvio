rule subset_DNA_reps_with_QCd_AA_reps_for_mapping:
    """Extract NT sequences based on final list of AA sequences for read recruitment"""
    input:
        reps=rules.remove_sequences_with_X_percent_gaps.output.fasta,
    output:
        NT_for_mapping=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-references_for_mapping_NT.fa",
        ),
        headers=os.path.join(dirs_dict["MSA"], "{group}", "{group}_headers.tmp"),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            "subset_DNA_reps_with_QCd_AA_reps_for_mapping_{group}.log",
        ),
    threads: M.T("subset_DNA_reps_with_QCd_AA_reps_for_mapping")
    params:
        fasta=rules.combine_sequence_data.output.NT_all,
    shell:
        """
        grep '>' {input.reps} | sed 's/>//g' > {output.headers}
        anvi-script-reformat-fasta {params.fasta} -I {output.headers} -o {output.NT_for_mapping} >> {log} 2>&1
        """


rule subset_external_gene_calls_file_all:
    """Subset the concatenated external_gene_calls.txt for the final set of NT sequences for profiling.
ALSO, once all external-gene-calls are concatenated, gene-callers-ids will be re-indexed so that
there are NO repeating gene-callers-id values.
"""
    input:
        headers=os.path.join(dirs_dict["MSA"], "{group}", "{group}_headers.tmp"),
    output:
        external_gene_calls_subset=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-external_gene_calls_subset.tsv",
        ),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"], "subset_external_gene_calls_file_all_{group}.log"
        ),
    threads: M.T("subset_external_gene_calls_file_all")
    params:
        external_gene_calls_all=rules.combine_sequence_data.output.external_gene_calls_all,
    script:
        "scripts/subset_external_gene_calls_file.py"


rule make_fasta_txt:
    """Format a fasta.txt with the filtered NT sequences for profiling in metagenomics workflow"""
    input:
        fasta=M.get_input_files_fasta_txt(),
    output:
        fasta_txt=os.path.join(dirs_dict["HOME"], "METAGENOMICS_WORKFLOW", "fasta.txt"),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "make_fasta_txt.log"),
    threads: M.T("make_fasta_txt")
    run:
        with open(output.fasta_txt, "a") as fasta_txt:
            # fasta-txt header
            fasta_txt.write("name\tpath\texternal_gene_calls\n")
            # per unique group name
            unique_group = []
            for hmm, value in M.hmm_dict.items():
                group = value["group"]
                if group not in unique_group:
                    # terrible hack with two .. in the relative path
                    fasta = os.path.join(
                        "..",
                        "..",
                        dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                        f"{group}",
                        f"{group}-references_for_mapping_NT.fa",
                    )
                    external_gene_calls = os.path.join(
                        "..",
                        "..",
                        dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                        f"{group}",
                        f"{group}-external_gene_calls_subset.tsv",
                    )
                    line = [group, fasta, external_gene_calls]
                    fasta_txt.write("\t".join(line) + "\n")
                    unique_group.append(group)
