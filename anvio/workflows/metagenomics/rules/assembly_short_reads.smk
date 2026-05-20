rule merge_fastas_for_co_assembly:
    """Merge paired short reads into FASTA input for IDBA-UD co-assembly."""
    input:
        unpack(lambda wildcards: M.get_sr_fastqs_for_group(wildcards.group)),
    output:
        fasta=temp(dirs_dict["QC_DIR"] + "/{group}-merged.fa"),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-merge_fastas_for_co_assembly.log",
    threads: M.T("merge_fastas_for_co_assembly")
    resources:
        nodes=M.T("merge_fastas_for_co_assembly"),
    run:
        concatenate = False
        if len(input.r1) > 1:
            concatenate = True
        r1_list, r2_list, did_unzip = uncompress_fastqs_if_needed(
            input.r1, input.r2, log
        )
        for r1, r2 in zip(r1_list, r2_list):
            # run fq2fa
            prefix1 = os.path.basename(r1).split(".fastq")[0]
            prefix2 = os.path.basename(r2).split(".fastq")[0]
            temp_merged_name = prefix1 + prefix2 + ".temp.fa"
            fq2fa_output = os.path.join(dirs_dict["QC_DIR"], temp_merged_name)
            shell("fq2fa --merge %s %s %s >> %s 2>&1" % (r1, r2, fq2fa_output, log))
            if did_unzip:
                # if we had to unzip the fastq files then now we delete the unzipped file
                shell("rm %s %s 2>>%s" % (r1, r2, log))
            if concatenate:
                # concatenate
                shell("cat %s >> %s 2>>%s" % (fq2fa_output, output, log))
                # delete individual fasta file
                shell("rm %s 2>>%s" % (fq2fa_output, log))
            else:
                shell("mv %s %s 2>>%s" % (fq2fa_output, output, log))


rule merge_fastqs_for_co_assembly:
    """Concatenate short-read FASTQs for co-assembly with assemblers that accept FASTQ."""
    input:
        unpack(lambda wildcards: M.get_sr_fastqs_for_group(wildcards.group)),
    output:
        r1=temp(os.path.join(dirs_dict["QC_DIR"], "{group}-merged_R1.fastq")),
        r2=temp(os.path.join(dirs_dict["QC_DIR"], "{group}-merged_R2.fastq")),
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-merge_fastqs_for_co_assembly.log",
    threads: M.T("merge_fastqs_for_co_assembly")
    resources:
        nodes=M.T("merge_fastqs_for_co_assembly"),
    run:
        r1_list, r2_list, did_unzip = uncompress_fastqs_if_needed(
            input.r1, input.r2, log
        )
        for r1, r2 in zip(r1_list, r2_list):
            # concatenate
            shell("cat %s >> %s 2>>%s" % (r1, output.r1, log))
            shell("cat %s >> %s 2>>%s" % (r2, output.r2, log))
            if did_unzip:
                # if we had to unzip the fastq files then now we delete the unzipped file
                shell("rm %s %s 2>>%s" % (r1, r2, log))


if M.get_param_value_from_config(["idba_ud", "run"]) and M.has_sr:

    rule idba_ud:
        """
        The purpose of this functino is to figure out whether co-assembly is done.
        if co-assembly is done then we need to first merge the fastq files
        """
        input:
            fasta=dirs_dict["QC_DIR"] + "/{group}-merged.fa",
        output:
            contigs=temp(
                os.path.join(dirs_dict["FASTA_DIR"], "{group}", "final.contigs.fa")
            ),
            # idba_ud generates both contigs.fa and scaffolds.fa
            # the user can choose which one will be used for the rest of the analysis
            # the other (unused) file will still be kept in it's raw form just in case the user wants it.
            raw_contigs_or_scaffold=(
                os.path.join(dirs_dict["FASTA_DIR"], "{group}", "contigs.fa")
                if M.use_scaffold_from_idba_ud
                else os.path.join(dirs_dict["FASTA_DIR"], "{group}", "scaffolds.fa")
            ),
        log:
            dirs_dict["LOGS_DIR"] + "/{group}-idba_ud.log",
        wildcard_constraints:
            group=SR_GRP_RE,
        conda:
            w.get_conda_yaml_path(M, "idba_ud")
        threads: M.T("idba_ud")
        resources:
            nodes=M.T("idba_ud"),
        params:
            env_prefix=lambda wildcards: w.get_conda_env_prefix(M, "idba_ud"),
            temp_dir=temp(dirs_dict["FASTA_DIR"] + "/{group}_TEMP"),
            mink=M.get_rule_param("idba_ud", "--mink"),
            maxk=M.get_rule_param("idba_ud", "--maxk"),
            step=M.get_rule_param("idba_ud", "--step"),
            inner_mink=M.get_rule_param("idba_ud", "--inner_mink"),
            inner_step=M.get_rule_param("idba_ud", "--inner_step"),
            prefix=M.get_rule_param("idba_ud", "--prefix"),
            min_count=M.get_rule_param("idba_ud", "--min_count"),
            min_support=M.get_rule_param("idba_ud", "--min_support"),
            seed_kmer=M.get_rule_param("idba_ud", "--seed_kmer"),
            min_contig=M.get_rule_param("idba_ud", "--min_contig"),
            similar=M.get_rule_param("idba_ud", "--similar"),
            max_mismatch=M.get_rule_param("idba_ud", "--max_mismatch"),
            min_pairs=M.get_rule_param("idba_ud", "--min_pairs"),
            no_bubble=M.get_rule_param("idba_ud", "--no_bubble"),
            no_local=M.get_rule_param("idba_ud", "--no_local"),
            no_coverage=M.get_rule_param("idba_ud", "--no_coverage"),
            no_correct=M.get_rule_param("idba_ud", "--no_correct"),
            pre_correction=M.get_rule_param("idba_ud", "--pre_correction"),
            use="scaffold.fa" if M.use_scaffold_from_idba_ud else "contig.fa",
            keep="contig.fa" if M.use_scaffold_from_idba_ud else "scaffold.fa",
        shell:
            """
            echo Running the following command: >> {log} 2>&1
            echo {params.env_prefix} idba_ud -o {params.temp_dir} --read {input.fasta} --num_threads {threads} {params.mink} {params.maxk} {params.step} {params.inner_mink} {params.inner_step} {params.prefix} {params.min_count} {params.min_support} {params.seed_kmer} {params.min_contig} {params.similar} {params.max_mismatch} {params.min_pairs} {params.no_bubble} {params.no_local} {params.no_coverage} {params.no_correct} {params.pre_correction} >> {log} 2>&1
            {params.env_prefix} idba_ud -o {params.temp_dir} --read {input.fasta} --num_threads {threads} {params.mink} {params.maxk} {params.step} {params.inner_mink} {params.inner_step} {params.prefix} {params.min_count} {params.min_support} {params.seed_kmer} {params.min_contig} {params.similar} {params.max_mismatch} {params.min_pairs} {params.no_bubble} {params.no_local} {params.no_coverage} {params.no_correct} {params.pre_correction} >> {log} 2>&1
            mv {params.temp_dir}/{params.use} {output.contigs}
            mv {params.temp_dir}/{params.keep} {output.raw_contigs_or_scaffold}
            rm -rf {params.temp_dir}
            """


if M.get_param_value_from_config(["metaspades", "run"]) and M.has_sr:

    def input_for_metaspades(wildcards):
        """Return either merged or per-readset FASTQ inputs for metaSPAdes."""
        input_file_dict = {}
        d = M.get_sr_fastqs_for_group(wildcards.group)
        if len(d["r1"]) > 1:
            input_file_dict = {
                "r1": temp(
                    os.path.join(
                        dirs_dict["QC_DIR"], wildcards.group + "-merged_R1.fastq"
                    )
                ),
                "r2": temp(
                    os.path.join(
                        dirs_dict["QC_DIR"], wildcards.group + "-merged_R2.fastq"
                    )
                ),
            }
        else:
            input_file_dict = {"r1": d["r1"], "r2": d["r2"]}
        return input_file_dict

    rule metaspades:
        """Assemble short-read metagenomes with metaSPAdes."""
        input:
            unpack(input_for_metaspades),
        output:
            contigs=temp(
                os.path.join(dirs_dict["FASTA_DIR"], "{group}", "final.contigs.fa")
            ),
            # metaspades generates both contigs.fasta and scaffolds.fasta
            # the user can choose which one will be used for the rest of the analysis
            # the other (unused) file will still be kept in it's raw form just in case the user wants it.
            raw_contigs_or_scaffold=(
                os.path.join(dirs_dict["FASTA_DIR"], "{group}", "contigs.fasta")
                if M.use_scaffold_from_metaspades
                else os.path.join(dirs_dict["FASTA_DIR"], "{group}", "scaffolds.fasta")
            ),
        log:
            dirs_dict["LOGS_DIR"] + "/{group}-metaspades.log",
        wildcard_constraints:
            group=SR_GRP_RE,
        conda:
            w.get_conda_yaml_path(M, "metaspades")
        threads: M.T("metaspades")
        resources:
            nodes=M.T("metaspades"),
        params:
            env_prefix=lambda wildcards: w.get_conda_env_prefix(M, "metaspades"),
            temp_dir=dirs_dict["FASTA_DIR"] + "/{group}_TEMP",
            additional_params=M.get_param_value_from_config(
                ["metaspades", "additional_params"]
            ),
            use="scaffolds.fasta" if M.use_scaffold_from_metaspades else "contigs.fasta",
            keep=(
                "contigs.fasta"
                if M.use_scaffold_from_metaspades
                else "scaffolds.fasta"
            ),
        shell:
            """
            {params.env_prefix} metaspades.py -1 {input.r1} -2 {input.r2} -t {threads} {params.additional_params} -o {params.temp_dir} >> {log} 2>&1
            mv {params.temp_dir}/{params.use} {output.contigs}
            mv {params.temp_dir}/{params.keep} {output.raw_contigs_or_scaffold}
            rm -rf {params.temp_dir}
            """


if M.get_param_value_from_config(["megahit", "run"]) and M.has_sr:

    rule megahit:
        """
        Assembling fastq files using megahit.

        All files created by megahit are stored in a temporary folder,
        and only the fasta file is kept for later analysis.
        """

        # Making this rule a shadow rule so all extra files created by megahit
        # are not retaineded (it is not enough to define the directory as temporary
        # because when failing in the middle of a run, snakemake doesn't delete directories)
        input:
            r1=lambda wildcards: M.get_sr_fastqs_for_group(wildcards.group)["r1"],
            r2=lambda wildcards: M.get_sr_fastqs_for_group(wildcards.group)["r2"],
        # Notice that megahit requires a directory to be specified as
        # output. If the directory already exists then megahit will not
        # run. To avoid this, the for megahit is a temporary directory,
        # once megahit is done running then the contigs database is moved
        # to the final location.
        output:
            contigs=temp(dirs_dict["FASTA_DIR"] + "/{group}/final.contigs.fa"),
        log:
            dirs_dict["LOGS_DIR"] + "/{group}-megahit.log",
        wildcard_constraints:
            group=SR_GRP_RE,
        conda:
            w.get_conda_yaml_path(M, "megahit")
        threads: M.T("megahit")
        resources:
            nodes=M.T("megahit"),
        params:
            env_prefix=lambda wildcards: w.get_conda_env_prefix(M, "megahit"),
            temp_dir=dirs_dict["FASTA_DIR"] + "/{group}_TEMP",
            # the minimum length for contigs (smaller contigs will be discarded)
            min_contig_len=M.get_rule_param("megahit", "--min-contig-len"),
            min_count=M.get_rule_param("megahit", "--min-count"),
            k_min=M.get_rule_param("megahit", "--k-min"),
            k_max=M.get_rule_param("megahit", "--k-max"),
            k_step=M.get_rule_param("megahit", "--k-step"),
            k_list=M.get_rule_param("megahit", "--k-list"),
            no_mercy=M.get_rule_param("megahit", "--no-mercy"),
            no_bubble=M.get_rule_param("megahit", "--no-bubble"),
            merge_level=M.get_rule_param("megahit", "--merge-level"),
            prune_level=M.get_rule_param("megahit", "--prune-level"),
            prune_depth=M.get_rule_param("megahit", "--prune-depth"),
            low_local_ratio=M.get_rule_param("megahit", "--low-local-ratio"),
            max_tip_len=M.get_rule_param("megahit", "--max-tip-len"),
            no_local=M.get_rule_param("megahit", "--no-local"),
            kmin_1pass=M.get_rule_param("megahit", "--kmin-1pass"),
            presets=M.get_rule_param("megahit", "--presets"),
            memory=M.get_rule_param("megahit", "--memory"),
            mem_flag=M.get_rule_param("megahit", "--mem-flag"),
            use_gpu=M.get_rule_param("megahit", "--use-gpu"),
            gpu_mem=M.get_rule_param("megahit", "--gpu-mem"),
            keep_tmp_files=M.get_rule_param("megahit", "--keep-tmp-files"),
            tmp_dir=M.get_rule_param("megahit", "--tmp-dir"),
            _continue=M.get_rule_param("megahit", "--continue"),
            verbose=M.get_rule_param("megahit", "--verbose"),
            r1=lambda wildcards, input: ",".join(input.r1),
            r2=lambda wildcards, input: ",".join(input.r2),
        shell:
            """
            {params.env_prefix} megahit -1 {params.r1} -2 {params.r2} -o {params.temp_dir} -t {threads} {params.min_contig_len} {params.min_count} {params.k_min} {params.k_max} {params.k_step} {params.k_list} {params.no_mercy} {params.no_bubble} {params.merge_level} {params.prune_level} {params.prune_depth} {params.low_local_ratio} {params.max_tip_len} {params.no_local} {params.kmin_1pass} {params.presets} {params.memory} {params.mem_flag} {params.use_gpu} {params.gpu_mem} {params.keep_tmp_files} {params.tmp_dir} {params._continue} {params.verbose} >> {log} 2>&1
            mv {params.temp_dir}/final.contigs.fa {output.contigs} >> {log} 2>&1
            rm -rf {params.temp_dir}
            """
