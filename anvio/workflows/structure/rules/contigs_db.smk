# Structure database rules for contigs databases: one structure database per `{group}`.
#
# Expects the following in scope:
#   M         — a StructureModule (or compatible) instance, already initialized
#   dirs_dict — M.dirs_dict (or equivalent)
#   rule_log  — the canonical log path helper of the including Snakefile

if not M.split_msa_and_predict:

    rule anvi_gen_structure_database:
        """Predict the structures of the genes of a contigs database in a single step."""
        input:
            unpack(M.get_input_for_structure_rules),
        output:
            structure_db=M.get_structure_db_path(),
        log:
            rule_log("anvi_gen_structure_database", "{group}-anvi_gen_structure_database"),
        threads: M.T("anvi_gen_structure_database")
        resources:
            nodes=M.T("anvi_gen_structure_database"),
            gpu=M.get_structure_gpu_resource("anvi_gen_structure_database"),
        params:
            genes_of_interest=M.get_structure_genes_of_interest_param,
            engine_params=M.get_structure_engine_params(),
        shell:
            """
                anvi-gen-structure-database -c {input.contigs_db} \
                                            -o {output.structure_db} \
                                            {params.genes_of_interest} \
                                            {params.engine_params} \
                                            --num-threads {threads} >> {log} 2>&1
            """

else:

    rule colabfold_msa:
        """Generate the ColabFold MSAs of the genes of a contigs database (CPU step)."""
        input:
            unpack(M.get_input_for_structure_rules),
        output:
            # a directory() output, because Snakemake does not create it before the job runs, and
            # --only-msa refuses to write into a --dump-dir that already exists. It is temp() so that
            # it goes away once the structures are predicted, without triggering a new MSA step later
            msa_dir=temp(directory(M.get_colabfold_msa_dir())),
        log:
            rule_log("colabfold_msa", "{group}-colabfold_msa"),
        threads: M.T("colabfold_msa")
        resources:
            nodes=M.T("colabfold_msa"),
            gpu=M.get_structure_gpu_resource("colabfold_msa"),
        params:
            genes_of_interest=M.get_structure_genes_of_interest_param,
            engine_params=M.get_structure_engine_params(),
        shell:
            """
                anvi-gen-structure-database -c {input.contigs_db} \
                                            --only-msa \
                                            --dump-dir {output.msa_dir} \
                                            {params.genes_of_interest} \
                                            {params.engine_params} \
                                            --num-threads {threads} >> {log} 2>&1
            """


    rule colabfold_predict:
        """Predict the structures from the ColabFold MSAs and build the structure database (GPU step)."""
        input:
            unpack(M.get_input_for_structure_rules),
            msa_dir=M.get_colabfold_msa_dir(),
        output:
            structure_db=M.get_structure_db_path(),
        log:
            rule_log("colabfold_predict", "{group}-colabfold_predict"),
        threads: M.T("colabfold_predict")
        resources:
            nodes=M.T("colabfold_predict"),
            gpu=M.get_structure_gpu_resource("colabfold_predict"),
        params:
            genes_of_interest=M.get_structure_genes_of_interest_param,
            engine_params=M.get_structure_engine_params(msa_source=False),
        shell:
            """
                anvi-gen-structure-database -c {input.contigs_db} \
                                            -o {output.structure_db} \
                                            --only-predict \
                                            --dump-dir {input.msa_dir} \
                                            {params.genes_of_interest} \
                                            {params.engine_params} \
                                            --num-threads {threads} >> {log} 2>&1
            """
