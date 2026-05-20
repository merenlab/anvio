rule annotate_contigs_database:
    """This is a dummy rule to mark that annotation was completed for each database

This allows workflows that inherit the contigs workflow to just track the output
of this rule, instead of tracking the individual
"""
    input:
        files=lambda wildcards: M.get_contigs_target_files(),
    output:
        done=touch(
            os.path.join(
                dirs_dict["CONTIGS_DIR"],
                "{group}-steps",
                "annotate_contigs_database.done",
            )
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-annotate_contigs_database.log",


rule anvi_script_reformat_fasta_prefix_only:
    """
Reformating the headers of the contigs fasta files.

This is required to make sure taht the headers don't contain
any charachters that anvi'o doesn't like.It give contigs
meaningful names; so that if the group name is 'MYSAMPLE01', the
contigs would look like this:
> MYSAMPLE01_000000000001
> MYSAMPLE01_000000000002

We do it in two steps in this workflow with the rules:
anvi_script_reformat_fasta_prefix_only and anvi_script_reformat_fasta_min_length
The first rule anvi_script_reformat_fasta_prefix_only does all the formating
of headers without considering the min_len parameter. The second rule takes
the output of the first rule as input and considers the min_len param.
The reasoning is that this way we have a fasta file that includes all the contigs
from the raw fasta file, but with headers that match the fasta final fasta file.
Hence, if later we decide we want to include shorter contigs, we can do it using
the output of anvi_script_reformat_fasta_prefix_only and the headers would be consistant.

"""
    input:
        contigs=M.get_input_fasta_path,
    output:
        # write protecting the contigs fasta file using protected() because
        # runnig the assembly is probably the most time consuming step and
        # we don't want anyone accidentaly deleting or changing this file.
        contigs=protected(
            dirs_dict["FASTA_DIR"]
            + "/{group}/{group}-contigs-prefix-formatted-only.fa"
        ),
        report=dirs_dict["FASTA_DIR"] + "/{group}/{group}-reformat-report.txt",
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_script_reformat_fasta_prefix_only.log",
    threads: M.T("anvi_script_reformat_fasta")
    resources:
        nodes=M.T("anvi_script_reformat_fasta"),
    params:
        prefix=M.get_rule_param("anvi_script_reformat_fasta", "--prefix"),
        keep_ids=M.get_rule_param("anvi_script_reformat_fasta", "--keep-ids"),
        exclude_ids=M.get_rule_param("anvi_script_reformat_fasta", "--exclude-ids"),
        simplify_names=M.get_rule_param(
            "anvi_script_reformat_fasta", "--simplify-names"
        ),
        seq_type=M.get_rule_param("anvi_script_reformat_fasta", "--seq-type"),
    shell:
        w.r("""anvi-script-reformat-fasta {input} \
                                             -o {output.contigs} \
                                             -r {output.report} \
                                             {params.prefix} \
                                             {params.exclude_ids} \
                                             {params.keep_ids} \
                                             {params.simplify_names} \
                                             {params.seq_type} >> {log} 2>&1""")


rule anvi_script_reformat_fasta:
    """
See docummentation for anvi_script_reformat_fasta_prefix_only
"""
    input:
        contigs=dirs_dict["FASTA_DIR"]
        + "/{group}/{group}-contigs-prefix-formatted-only.fa",
    output:
        contigs=temp(
            os.path.join(dirs_dict["FASTA_DIR"], "{group}", "{group}-contigs.fa")
        ),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_script_reformat_fasta.log",
    threads: M.T("anvi_script_reformat_fasta")
    resources:
        nodes=M.T("anvi_script_reformat_fasta"),
    params:
        prefix="{group}",
        min_len=M.get_rule_param("anvi_script_reformat_fasta", "--min-len"),
        max_len=M.get_rule_param("anvi_script_reformat_fasta", "--max-len"),
        seq_type=M.get_rule_param("anvi_script_reformat_fasta", "--seq-type"),
    shell:
        w.r("""anvi-script-reformat-fasta {input} \
                                             -o {output.contigs} \
                                             {params.seq_type} \
                                             {params.min_len} \
                                             {params.max_len} >> {log} 2>&1""")


rule reformat_external_gene_calls_table:
    """Normalize external gene-call tables before contigs database creation."""
    input:
        unpack(M.get_input_for_reformat_external_gene_calls_table),
    output:
        target=os.path.join(
            dirs_dict["FASTA_DIR"], "{group}", "{group}-external-gene-calls.txt"
        ),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"], "{group}-reformat_external_gene_calls_table.log"
        ),
    threads: M.T("reformat_external_gene_calls_table")
    resources:
        nodes=M.T("reformat_external_gene_calls_table"),
    run:
        import pandas as pd
        import anvio.utils as u
        import anvio.fastalib as f

        fa = f.ReadFasta(input.contigs[0])
        contigs = fa.ids
        reformat_report = pd.read_csv(
            input.reformat_report[0], sep="\t", index_col=1, header=None
        )
        contig_dict = next(iter(reformat_report.to_dict().values()))
        external_gene_calls = u.get_TAB_delimited_file_as_dictionary(
            input.external_gene_calls
        )
        contigs_in_external_gene_calls_not_in_contigs_database = [
            c for c in contig_dict.values() if c not in contigs
        ]
        if contigs_in_external_gene_calls_not_in_contigs_database:
            warning = (
                "The following contigs are in your external gene calls "
                + "but not in your reformatted fasta file (%s): %s. This "
                % (
                    input.contigs[0],
                    ", ".join(
                        contigs_in_external_gene_calls_not_in_contigs_database
                    ),
                )
                + "is no cause for concern since they have been removed due "
                + "to the --min-len criteria of anvi-script-reformat-fasta, "
                + "but we still thought you should know. Any gene that was called "
                + "within these contigs will naturally not be included in your contigs database."
            )
            shell("echo '%s' > %s" % (warning, log))
        new_external_gene_calls = {}
        for g in external_gene_calls:
            contig_name = external_gene_calls[g]["contig"]
            if contig_name not in contig_dict:
                raise ConfigError(
                    "When proccessing external gene calls for %s "
                    "we bumped into the following problem: The "
                    'following contig name: "%s" appears in your external '
                    "gene calls file, but not in your fasta file."
                    % (wildcards.group, contig_name)
                )
            new_contig_name = contig_dict[contig_name]
            if new_contig_name in contigs:
                new_external_gene_calls[g] = external_gene_calls[g].copy()
                new_external_gene_calls[g]["contig"] = new_contig_name
        #            u.store_dict_as_TAB_delimited_file(new_external_gene_calls, output[0], key_header='gene_callers_id')
        D = pd.DataFrame.from_dict(new_external_gene_calls, orient="index")
        D.to_csv(output[0], index_label="gene_callers_id", sep="\t")


rule anvi_gen_contigs_database:
    """Generates a contigs database using anvi-gen-contigs-database"""

    # depending on whether human contamination using centrifuge was done
    # or not, the input to this rule will be the raw assembly or the
    # filtered.
    input:
        unpack(M.get_input_for_anvi_gen_contigs_database),
    output:
        db=os.path.join(dirs_dict["CONTIGS_DIR"], "{group}.db"),
    # Setting the version to the same as that of the contigs__version in anvi'o
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-anvi_gen_contigs_database.log",
    threads: M.T("anvi_gen_contigs_database")
    resources:
        nodes=M.T("anvi_gen_contigs_database"),
    params:
        description=M.get_rule_param("anvi_gen_contigs_database", "--description"),
        skip_gene_calling=M.get_rule_param(
            "anvi_gen_contigs_database", "--skip-gene-calling"
        ),
        external_gene_calls=lambda wildcards: M.get_external_gene_calls_param(wildcards),
        ignore_internal_stop_codons=M.get_rule_param(
            "anvi_gen_contigs_database", "--ignore-internal-stop-codons"
        ),
        skip_predict_frame=M.get_rule_param(
            "anvi_gen_contigs_database", "--skip-predict-frame"
        ),
        skip_mindful_splitting=M.get_rule_param(
            "anvi_gen_contigs_database", "--skip-mindful-splitting"
        ),
        contigs_fasta=M.get_rule_param("anvi_gen_contigs_database", "--contigs-fasta"),
        project_name=M.get_rule_param("anvi_gen_contigs_database", "--project-name"),
        split_length=M.get_rule_param("anvi_gen_contigs_database", "--split-length"),
        kmer_size=M.get_rule_param("anvi_gen_contigs_database", "--kmer-size"),
        prodigal_translation_table=lambda wildcards: M.get_prodigal_translation_table_param(
            wildcards
        ),
    shell:
        w.r("anvi-gen-contigs-database -f {input.fasta} \
                                          -o {output.db} \
                                          {params.kmer_size} \
                                          {params.split_length} \
                                          {params.project_name} \
                                          {params.contigs_fasta} \
                                          {params.skip_mindful_splitting} \
                                          {params.ignore_internal_stop_codons} \
                                          {params.external_gene_calls} \
                                          {params.skip_predict_frame} \
                                          {params.skip_gene_calling} \
                                          {params.prodigal_translation_table} \
                                          {params.description} \
                                          -T {threads} >> {log} 2>&1")
