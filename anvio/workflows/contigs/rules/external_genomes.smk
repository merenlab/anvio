rule gen_external_genome_file:
    input:
        annotation_done=expand(
            os.path.join(
                dirs_dict["CONTIGS_DIR"],
                "{sample}-steps",
                "annotate_contigs_database.done",
            ),
            sample=M.group_names,
        ),
        contigs_db_paths=expand(
            dirs_dict["CONTIGS_DIR"] + "/{sample}.db", sample=M.group_names
        ),
    output:
        target=M.external_genomes_file,
    log:
        dirs_dict["LOGS_DIR"] + "" + "/gen_external_genome_file.log",
    threads: M.T("gen_external_genome_file")
    resources:
        nodes=M.T("gen_external_genome_file"),
    run:
        contigs_db_project_names = [
            ContigsDatabase(contigs_db_path).meta["project_name"]
            for contigs_db_path in input.contigs_db_paths
        ]
        contigs_db_name_path_tuples = zip(
            contigs_db_project_names, input.contigs_db_paths
        )
        with open(output[0], "w") as f:
            f.write("name\tcontigs_db_path\n")
            for name, path in contigs_db_name_path_tuples:
                name = name.replace(".", "_").replace("-", "_").replace(" ", "_")
                f.write("%s\t%s\n" % (name, path))
