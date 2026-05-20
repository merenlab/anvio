rule anvi_get_sequences_for_gene_clusters:
    """Export gene-cluster sequences for phylogenetic analysis."""
    input:
        genomes_storage=os.path.join(
            dirs_dict["PAN_DIR"], M.project_name + "-GENOMES.db"
        ),
        pan_db=M.pan_db_path,
    output:
        fasta=os.path.join(dirs_dict["PHYLO_DIR"], M.project_name + "-GC-sequences.fa"),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            M.project_name + "-anvi_get_sequences_for_gene_clusters.log",
        ),
    threads: M.T("anvi_get_sequences_for_gene_clusters")
    resources:
        nodes=M.T("anvi_get_sequences_for_gene_clusters"),
    params:
        gene_cluster_id=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--gene-cluster-id"
        ),
        gene_cluster_ids_file=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--gene-cluster-ids-file"
        ),
        collection_name=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--collection-name"
        ),
        bin_id=M.get_rule_param("anvi_get_sequences_for_gene_clusters", "--bin-id"),
        min_num_genomes_gene_cluster_occurs=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters",
            "--min-num-genomes-gene-cluster-occurs",
        ),
        max_num_genomes_gene_cluster_occurs=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters",
            "--max-num-genomes-gene-cluster-occurs",
        ),
        min_num_genes_from_each_genome=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--min-num-genes-from-each-genome"
        ),
        max_num_genes_from_each_genome=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--max-num-genes-from-each-genome"
        ),
        max_num_gene_clusters_missing_from_genome=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters",
            "--max-num-gene-clusters-missing-from-genome",
        ),
        min_functional_homogeneity_index=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters",
            "--min-functional-homogeneity-index",
        ),
        max_functional_homogeneity_index=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters",
            "--max-functional-homogeneity-index",
        ),
        min_geometric_homogeneity_index=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--min-geometric-homogeneity-index"
        ),
        max_geometric_homogeneity_index=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--max-geometric-homogeneity-index"
        ),
        add_into_items_additional_data_table=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters",
            "--add-into-items-additional-data-table",
        ),
        concatenate_gene_clusters=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--concatenate-gene-clusters"
        ),
        separator=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--separator"
        ),
        align_with=M.get_rule_param(
            "anvi_get_sequences_for_gene_clusters", "--align-with"
        ),
    shell:
        """
           anvi-get-sequences-for-gene-clusters -p {input.pan_db} -g {input.genomes_storage} \
           {params.gene_cluster_id} {params.gene_cluster_ids_file} {params.collection_name} \
           {params.bin_id} {params.min_num_genomes_gene_cluster_occurs} \
           {params.max_num_genomes_gene_cluster_occurs} {params.min_num_genes_from_each_genome} \
           {params.max_num_genes_from_each_genome} {params.max_num_gene_clusters_missing_from_genome} \
           {params.min_functional_homogeneity_index} {params.max_functional_homogeneity_index} \
           {params.min_geometric_homogeneity_index} {params.max_geometric_homogeneity_index} \
           {params.add_into_items_additional_data_table} {params.concatenate_gene_clusters} \
           {params.separator} {params.align_with} \
           -o {output} --just-do-it >> {log} 2>&1
        """


rule import_phylogenetic_tree_to_pangenome:
    """Import a computed phylogenetic tree into the pangenome database."""
    input:
        pan_db=M.pan_db_path,
        newick=os.path.join(
            M.dirs_dict["PHYLO_DIR"],
            M.project_name + "-proteins_GAPS_REMOVED.fa" + ".contree",
        ),
    output:
        layers_orders_file=temp(
            os.path.join(M.dirs_dict["PAN_DIR"], M.project_name + "-layers-orders.txt")
        ),
        phylogeny_imported=touch(M.get_phylogeny_imported_flag()),
    log:
        os.path.join(
            M.dirs_dict["LOGS_DIR"],
            M.pan_project_name + "-import_phylogenetic_tree_to_pangenome.log",
        ),
    threads: M.T("import_phylogenetic_tree_to_pangenome")
    resources:
        nodes=M.T("import_phylogenetic_tree_to_pangenome"),
    params:
        tree_name=M.tree_name,
        just_do_it=M.get_rule_param(
            "import_phylogenetic_tree_to_pangenome", "--just-do-it"
        ),
    shell:
        """
        echo -e "item_name\tdata_type\tdata_value" > {output.layers_orders_file}
        echo -e "{params.tree_name}\tnewick\t`cat {input.newick}`" >> {output.layers_orders_file}
        anvi-import-misc-data -p {input.pan_db} -t layer_orders {params.just_do_it} {output.layers_orders_file} >> {log} 2>&1
        """


rule anvi_compute_genome_similarity:
    """Compute genome similarity metrics for the pangenome inputs."""
    input:
        unpack(lambda wildcards: M.input_for_anvi_compute_genome_similarity),
    output:
        anvi_compute_genome_similarity_flag=touch(M.anvi_compute_genome_similarity_flag),
        output_dir=directory(M.anvi_compute_genome_similarity_output_dir),
    log:
        os.path.join(
            M.dirs_dict["LOGS_DIR"],
            M.pan_project_name + "-anvi_compute_genome_similarity.log",
        ),
    threads: M.T("anvi_compute_genome_similarity")
    resources:
        nodes=M.T("anvi_compute_genome_similarity"),
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
        additional_params=M.get_param_value_from_config(
            ["anvi_compute_genome_similarity", "additional_params"]
        ),
    shell:
        " anvi-compute-genome-similarity {params.internal_genomes_argument} \
                                            {params.external_genomes_argument} \
                                            -T {threads} \
                                            -o {output.output_dir} \
                                            -p {input.pan_db} \
                                            {params.additional_params} >> {log} 2>&1"
