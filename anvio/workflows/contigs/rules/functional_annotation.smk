rule reformat_external_functions:
    """Remove any gene that is not included in the contigs database due to fasta reformt with --min-len"""
    input:
        unpack(M.get_input_for_reformat_external_functions),
    output:
        target=os.path.join(
            dirs_dict["FASTA_DIR"],
            "{group}",
            "{group}" + "-gene-functional-annotation.txt",
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{group}-reformat_external_functions.log"),
    threads: M.T("reformat_external_functions")
    resources:
        nodes=M.T("reformat_external_functions"),
    run:
        import pandas as pd

        gene_functional_annotation_file = input.gene_functional_annotation
        external_gene_calls_file = input.external_gene_calls
        gene_functional_annotation = pd.read_csv(
            gene_functional_annotation_file, sep="\t", index_col=0
        )
        external_gene_calls = pd.read_csv(
            external_gene_calls_file, sep="\t", index_col=0
        )
        genes_in_both = [
            g
            for g in gene_functional_annotation.index
            if g in external_gene_calls.index
        ]
        genes_only_in_functional_annotation_file = [
            g
            for g in gene_functional_annotation.index
            if g not in external_gene_calls.index
        ]
        if genes_only_in_functional_annotation_file:
            warning = (
                "When configuring gene functional annotations for %s\n"
                % wildcards.group
                + "we noticed that the following gene caller ids are not in your contigs database: \n'%s'.\n"
                % ", ".join(genes_only_in_functional_annotation_file)
                + "These are probably missing because you used the --min-len parameter when you ran\n"
                + "anvi-script-reformat-fasta. And hence, there is probably no reason for concern.\n"
            )
            shell("echo -e '%s' > %s" % (warning, log))
        new_gene_functional_annotation = gene_functional_annotation.loc[
            genes_in_both,
        ]
        new_gene_functional_annotation.to_csv(output[0], sep="\t")


rule import_external_functions:
    """Import externally supplied functional annotations into the contigs database."""
    input:
        unpack(M.get_input_for_import_external_functions),
    output:
        functions_imported=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"],
                "{group}-steps",
                "external_functions_imported.done",
            )
        ),
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{group}-import_external_functions.log"),
    threads: M.T("import_external_functions")
    resources:
        nodes=M.T("import_external_functions"),
    shell:
        "anvi-import-functions -c {input.contigs} -i {input.gene_functional_annotation} >> {log} 2>&1"
