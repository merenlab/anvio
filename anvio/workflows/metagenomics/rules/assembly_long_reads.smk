if M.get_param_value_from_config(["flye", "run"]) and M.has_lr:

    rule flye:
        """Assemble long-read metagenomes with metaFlye/Flye."""
        input:
            reads=lambda wildcards: M.get_lr_files_for_group(wildcards.group),
        output:
            raw_fasta=M.dirs_dict["FASTA_DIR"] + "/{group}/final.contigs.fa",
        log:
            rule_log("flye", "{group}-flye"),
        wildcard_constraints:
            group=LR_GRP_RE,
        conda:
            w.get_conda_yaml_path(M, "flye")
        threads: M.T("flye")
        resources:
            nodes=M.T("flye"),
        params:
            env_prefix=lambda wildcards: w.get_conda_env_prefix(M, "flye"),
            meta=M.get_rule_param("flye", "--meta"),
            read_type=lambda wildcards: M.get_flye_flag_for_group(wildcards.group),
            additional_params=M.get_param_value_from_config(
                ["flye", "additional_params"]
            ),
            outdir=M.dirs_dict["FASTA_DIR"] + "/{group}",
            genome_size=M.get_rule_param("flye", "--genome-size"),
            iterations=M.get_rule_param("flye", "--iterations"),
            min_overlap=M.get_rule_param("flye", "--min-overlap"),
            read_error=M.get_rule_param("flye", "--read-error"),
            keep_haplotypes=M.get_rule_param("flye", "--keep-haplotypes"),
            no_alt_contigs=M.get_rule_param("flye", "--no-alt-contigs"),
            scaffold=M.get_rule_param("flye", "--scaffold"),
            polish_target=M.get_rule_param("flye", "--polish-target"),
        shell:
            r"""
            mkdir -p {params.outdir}

            {params.env_prefix} flye {params.meta} {params.read_type} {input.reads} \
                 -o {params.outdir} -t {threads} \
                 {params.genome_size} {params.iterations} {params.min_overlap} {params.read_error} \
                 {params.keep_haplotypes} {params.no_alt_contigs} {params.scaffold} {params.polish_target} \
                 {params.additional_params} \
                 >> {log} 2>&1

            if [ -f "{params.outdir}/assembly.fasta" ]; then
                cp "{params.outdir}/assembly.fasta" {output.raw_fasta} && rm "{params.outdir}/assembly.fasta"
            elif [ -f "{params.outdir}/assembly.fa" ]; then
                cp "{params.outdir}/assembly.fa" {output.raw_fasta} && rm "{params.outdir}/assembly.fa"
            else
                echo "Could not find (meta)Flye assembly fasta in {params.outdir}" >> {log}
                exit 1
            fi
            """
