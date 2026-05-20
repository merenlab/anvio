# run pangenome
rule anvi_pan_genome:
    """Compute a pangenome from the genomes storage database."""
    input:
        pan_db=dirs_dict["PAN_DIR"] + "/" + M.project_name + "-GENOMES.db",
    output:
        pan_db=os.path.join(dirs_dict["PAN_DIR"], M.pan_project_name + "-PAN.db"),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], M.project_name + "-anvi_pan_genome.log"),
    threads: M.T("anvi_pan_genome")
    resources:
        nodes=M.T("anvi_pan_genome"),
    params:
        output_dir=dirs_dict["PAN_DIR"],
        genome_names=M.get_rule_param("anvi_pan_genome", "--genome-names"),
        project_name=M.pan_project_name,
        skip_alignments=M.get_rule_param("anvi_pan_genome", "--skip-alignments"),
        align_with=M.get_rule_param("anvi_pan_genome", "--align-with"),
        exclude_partial_gene_calls=M.get_rule_param(
            "anvi_pan_genome", "--exclude-partial-gene-calls"
        ),
        use_ncbi_blast=M.get_rule_param("anvi_pan_genome", "--use-ncbi-blast"),
        minbit=M.get_rule_param("anvi_pan_genome", "--minbit"),
        mcl_inflation=M.get_rule_param("anvi_pan_genome", "--mcl-inflation"),
        min_occurrence=M.get_rule_param("anvi_pan_genome", "--min-occurrence"),
        min_percent_identity=M.get_rule_param(
            "anvi_pan_genome", "--min-percent-identity"
        ),
        description=M.get_rule_param("anvi_pan_genome", "--description"),
        overwrite_output_destinations=M.get_rule_param(
            "anvi_pan_genome", "--overwrite-output-destinations"
        ),
        skip_hierarchical_clustering=M.get_rule_param(
            "anvi_pan_genome", "--skip-hierarchical-clustering"
        ),
        enforce_hierarchical_clustering=M.get_rule_param(
            "anvi_pan_genome", "--enforce-hierarchical-clustering"
        ),
        distance=M.get_rule_param("anvi_pan_genome", "--distance"),
        linkage=M.get_rule_param("anvi_pan_genome", "--linkage"),
        i_know_this_is_not_a_good_idea=M.get_rule_param(
            "anvi_pan_genome", "--I-know-this-is-not-a-good-idea"
        ),
    shell:
        """
            anvi-pan-genome -g {input} --num-threads {threads} -o {params.output_dir} --project-name {params.project_name} {params.genome_names}\
            {params.skip_alignments} {params.align_with} {params.exclude_partial_gene_calls}\
            {params.use_ncbi_blast} {params.minbit} {params.mcl_inflation}\
            {params.min_occurrence} {params.min_percent_identity} \
            {params.description} {params.overwrite_output_destinations}\
            {params.skip_hierarchical_clustering} {params.enforce_hierarchical_clustering}\
            {params.distance} {params.linkage} {params.i_know_this_is_not_a_good_idea} >> {log} 2>&1
        """


# generate anvi'o genomes storage
rule anvi_gen_genomes_storage:
    """Build the genomes storage database from internal and external genome tables."""
    input:
        unpack(lambda wildcards: M.input_for_anvi_gen_genomes_storage),
    output:
        pan_db=dirs_dict["PAN_DIR"] + "/" + M.project_name + "-GENOMES.db",
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"], M.project_name + "-anvi_gen_genomes_storage.log"
        ),
    threads: M.T("anvi_gen_genomes_storage")
    resources:
        nodes=M.T("anvi_gen_genomes_storage"),
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
        gene_caller=M.get_rule_param("anvi_gen_genomes_storage", "--gene-caller"),
    shell:
        """
            anvi-gen-genomes-storage -o {output}\
                                     {params.internal_genomes_argument}\
                                     {params.external_genomes_argument}\
                                     {params.gene_caller} >> {log} 2>&1
        """
