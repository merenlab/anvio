rule anvi_get_sequences_for_hmm_hits:
    """Export HMM-hit sequences for phylogenomic tree building."""
    input:
        unpack(lambda wildcards: M.input_for_anvi_get_sequences_for_hmm_hits),  # The lambda function here is just a trick. from some reason without it, snakemake can't unpack the dict
    output:
        fasta=os.path.join(dirs_dict["PHYLO_DIR"], M.project_name + "-proteins.fa"),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            M.project_name + "_anvi_get_sequences_for_hmm_hits.log",
        ),
    threads: M.T("anvi_get_sequences_for_hmm_hits")
    resources:
        nodes=M.T("anvi_get_sequences_for_hmm_hits"),
    params:
        internal_genomes_argument=lambda wildcards: (
            "--internal-genomes " + M.internal_genomes_file
            if M.internal_genomes_file
            else ""
        ),
        external_genomes_argument=lambda wildcards: (
            "--external-genomes " + M.external_genomes_file
            if M.external_genomes_file
            else ""
        ),
        return_best_hit=M.get_rule_param(
            "anvi_get_sequences_for_hmm_hits", "--return-best-hit"
        ),
        separator=M.get_rule_param("anvi_get_sequences_for_hmm_hits", "--separator"),
        align_with=M.get_rule_param("anvi_get_sequences_for_hmm_hits", "--align-with"),
        min_num_bins_gene_occurs=M.get_rule_param(
            "anvi_get_sequences_for_hmm_hits", "--min-num-bins-gene-occurs"
        ),
        max_num_genes_missing_from_bin=M.get_rule_param(
            "anvi_get_sequences_for_hmm_hits", "--max-num-genes-missing-from-bin"
        ),
        concatenate_genes=M.get_rule_param(
            "anvi_get_sequences_for_hmm_hits", "--concatenate-genes"
        ),
        get_aa_sequences=M.get_rule_param(
            "anvi_get_sequences_for_hmm_hits", "--get-aa-sequences"
        ),
        gene_names=M.get_rule_param("anvi_get_sequences_for_hmm_hits", "--gene-names"),
        hmm_sources=M.get_rule_param("anvi_get_sequences_for_hmm_hits", "--hmm-sources"),
    shell:
        """
            anvi-get-sequences-for-hmm-hits {params.internal_genomes_argument} {params.external_genomes_argument} \
                                            {params.return_best_hit} {params.separator} {params.align_with} \
                                            {params.min_num_bins_gene_occurs} {params.max_num_genes_missing_from_bin} \
                                            {params.concatenate_genes} {params.get_aa_sequences} {params.gene_names} \
                                            {params.hmm_sources} -o {output} >> {log} 2>&1
        """


rule trimal:
    """Trim poorly aligned columns from the phylogenomic alignment."""
    input:
        source=M.phylogenomics_sequence_file,
    output:
        target=os.path.join(
            dirs_dict["PHYLO_DIR"], M.project_name + "-proteins_GAPS_REMOVED.fa"
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], M.project_name + "_trimal.log"),
    threads: M.T("trimal")
    resources:
        nodes=M.T("trimal"),
    params:
        gt=M.get_rule_param("trimal", "-gt"),
        additional_params=M.get_param_value_from_config(["trimal", "additional_params"]),
    shell:
        """
            trimal -in {input} \
                   -out {output} \
                   {params.gt} \
                   {params.additional_params} >> {log} 2>&1
        """


rule iqtree:
    """Infer a phylogenomic tree with IQ-TREE."""
    input:
        source=os.path.join(
            dirs_dict["PHYLO_DIR"], M.project_name + "-proteins_GAPS_REMOVED.fa"
        ),
    output:
        target=os.path.join(
            dirs_dict["PHYLO_DIR"],
            M.project_name + "-proteins_GAPS_REMOVED.fa" + ".contree",
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], M.project_name + "_iqtree.log"),
    threads: M.T("iqtree")
    resources:
        nodes=M.T("iqtree"),
    params:
        m=M.get_rule_param("iqtree", "-m"),
        bb=M.get_rule_param("iqtree", "-bb"),
        additional_params=M.get_param_value_from_config(["iqtree", "additional_params"]),
    shell:
        """
            iqtree -s {input} \
                       -nt {threads} \
                       {params.m} \
                       {params.bb} \
                       {params.additional_params} >> {log} 2>&1
        """
