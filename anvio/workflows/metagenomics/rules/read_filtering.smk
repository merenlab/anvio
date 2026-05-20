def need_to_uncompress_fastqs(r1_list, r2_list):
    # Checking if any of the input files are compressed
    # this should only happen if QC is not performed in this snakemake session
    # because if QC is performed in this session then the output of iu-filter-quality-minoche
    # is not compressed and snakemake will schedule fq2fa rule before gzip rule
    pair_flags = []
    for i, (r1, r2) in enumerate(zip(r1_list, r2_list)):
        gz1, gz2 = r1.endswith(".gz"), r2.endswith(".gz")
        if gz1 != gz2:
            raise ConfigError(
                f"Something seems very bad: one of the pair fastq files "
                f"is compressed and the other one is not: {r1} vs {r2}"
            )
        pair_flags.append(gz1)
    any_gz = any(pair_flags)
    all_gz = all(pair_flags)
    if any_gz and not all_gz:
        raise ConfigError(
            "Mixed of compressed/uncompressed set of r1 and r2 fastq files. "
            "Make sure they are uniform."
        )
    return all_gz  # True > need to gunzip all


def uncompress_fastqs_if_needed(r1_list, r2_list, log):
    """Return uncompressed FASTQ paths when compressed inputs need temporary expansion."""
    did_unzip = False
    out1, out2 = [], []
    if need_to_uncompress_fastqs(r1_list, r2_list):
        did_unzip = True
        for r1, r2 in zip(r1_list, r2_list):
            r1_unzipped = r1.replace(".fastq.gz", "_temp.fastq")
            r2_unzipped = r2.replace(".fastq.gz", "_temp.fastq")
            w.gunzip_file(r1, log, shell, output_path=r1_unzipped)
            w.gunzip_file(r2, log, shell, output_path=r2_unzipped)
            out1.append(r1_unzipped)
            out2.append(r2_unzipped)
    else:
        out1, out2 = r1_list, r2_list
    return out1, out2, did_unzip


def concatenate_fastqs(inputs, dest, log):
    """Concatenate multiple FASTQ files when co-assembly needs one input stream."""
    if len(inputs) == 1:
        return inputs[0]  # nothing to do
    shell("cat {files} > {dest} 2>> {log}".format(files=" ".join(inputs), dest=dest))
    return dest


def input_for_remove_short_reads_based_on_references(wildcards):
    """Build inputs for removing reads recruited to configured reference contigs."""
    d = {}
    d["bam"] = expand(
        os.path.join(dirs_dict["MAPPING_DIR"], "{group}", "{readset}.bam"),
        group=list(M.references_for_removal.keys()),
        readset=wildcards.readset,
    )
    d.update(M.get_fastq(wildcards.readset, pre_ref_removal=True))
    return d


gzip_suffix = ".gz" if run_gzip_fastqs else ""


rule remove_short_reads_based_on_references:
    """Remove or report short reads that map to configured reference contigs."""
    input:
        unpack(input_for_remove_short_reads_based_on_references),
    output:
        ids_to_remove=os.path.join(dirs_dict["QC_DIR"], "{readset}-ids-to-remove.txt"),
        # For each one of the four following outputs we have an alternative
        # file name for the case that the user chose dont_remove_just_map
        # These mock files end with '-mock-output'
        keep1=(
            os.path.join(
                dirs_dict["QC_DIR"], "{readset}-FILTERED_R1.fastq" + gzip_suffix
            )
            if M.remove_short_reads_based_on_references
            else temp(os.path.join(dirs_dict["QC_DIR"], "{readset}-mock-output1"))
        ),
        remove1=(
            temp(os.path.join(dirs_dict["QC_DIR"], "{readset}.R1.fastq.removed"))
            if M.remove_short_reads_based_on_references
            else temp(os.path.join(dirs_dict["QC_DIR"], "{readset}-mock-output2"))
        ),
        keep2=(
            os.path.join(
                dirs_dict["QC_DIR"], "{readset}-FILTERED_R2.fastq" + gzip_suffix
            )
            if M.remove_short_reads_based_on_references
            else temp(os.path.join(dirs_dict["QC_DIR"], "{readset}-mock-output3"))
        ),
        remove2=(
            temp(os.path.join(dirs_dict["QC_DIR"], "{readset}.R2.fastq.removed"))
            if M.remove_short_reads_based_on_references
            else temp(os.path.join(dirs_dict["QC_DIR"], "{readset}-mock-output4"))
        ),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            "{readset}-remove_short_reads_based_on_references.log",
        ),
    threads: M.T("remove_short_reads_based_on_references")
    resources:
        nodes=M.T("remove_short_reads_based_on_references"),
    params:
        delimiter=M.get_param_value_from_config(
            [
                "remove_short_reads_based_on_references",
                "delimiter-for-iu-remove-ids-from-fastq",
            ]
        ),
    run:
        references = M.references_for_removal.keys()
        for bam in input.bam:
            # merge the ids from matching any reference
            shell(
                "samtools view %s | cut -f 1 >> %s 2>>{log}"
                % (bam, output.ids_to_remove)
            )
        ids_to_remove = set(open(output.ids_to_remove).read().splitlines())
        with open(output.ids_to_remove, "w") as f:
            for item in ids_to_remove:
                f.write("%s\n" % item)
        if M.remove_short_reads_based_on_references:
            if ids_to_remove:
                r1_list, r2_list, did_unzip = uncompress_fastqs_if_needed(
                    input.r1, input.r2, log
                )
                # concatenate the r1 and r2 files:
                r1 = concatenate_fastqs(
                    r1_list,
                    os.path.join(
                        dirs_dict["QC_DIR"],
                        f"{wildcards.readset}_R1_concat.temp.fastq",
                    ),
                    log,
                )
                r2 = concatenate_fastqs(
                    r2_list,
                    os.path.join(
                        dirs_dict["QC_DIR"],
                        f"{wildcards.readset}_R2_concat.temp.fastq",
                    ),
                    log,
                )
                # remove ids
                intermediate_fastq = {"1": r1, "2": r2}
                output_fastq_files = {
                    "remove": {"1": output.remove1, "2": output.remove2},
                    "keep": {"1": output.keep1, "2": output.keep2},
                }
                output_prefix = dirs_dict["QC_DIR"]
                for i in intermediate_fastq:
                    # first make sure that the expected output files don't already exist
                    # if they exist, then we delete them, since they are probably left from previous failed run
                    shell(
                        "rm -rf {fq}.removed {fq}.survived".format(
                            fq=intermediate_fastq[i]
                        )
                    )
                    shell(
                        'iu-remove-ids-from-fastq -i %s -l {output.ids_to_remove} -d "{params.delimiter}"'
                        % intermediate_fastq[i]
                    )
                    shell(
                        "mv %s.removed %s"
                        % (intermediate_fastq[i], output_fastq_files["remove"][i])
                    )
                    if run_gzip_fastqs:
                        shell(
                            "gzip < %s.survived > %s"
                            % (intermediate_fastq[i], output_fastq_files["keep"][i])
                        )
                        shell("rm %s.survived" % intermediate_fastq[i])
                    else:
                        shell(
                            "mv %s.survived %s"
                            % (intermediate_fastq[i], output_fastq_files["keep"][i])
                        )
                if did_unzip or r1.endswith(".temp.fastq"):
                    shell("rm -f {r1}")
                if did_unzip or r2.endswith(".temp.fastq"):
                    shell("rm -f {r2}")
            else:
                # No reads from this sample were mapped to references so we will just:
                # change the name of the input fastq to the expected output fastq
                # and then create mock files for the "removed" fastq files (just so snakemake is happy and sees all expected output)
                shell("mv {input.r1} {output.keep1}")
                shell("mv {input.r2} {output.keep2}")
                shell("touch {output.remove1} {output.remove2}")
        else:
            # we only need the list of ids that matched the references
            # we are not removing enything so we just create mock output files to make snakemake happy
            shell(
                "touch {output.keep1} {output.keep2} {output.remove1} {output.remove2}"
            )


rule gen_report_for_mapping_to_references_for_removal:
    """Summarize read removal caused by reference mapping."""
    input:
        ids_to_remove=expand(
            os.path.join(dirs_dict["QC_DIR"], "{readset}-ids-to-remove.txt"),
            readset=SR_READSETS,
        ),
    output:
        report_file=os.path.join(dirs_dict["QC_DIR"], "short-read-removal-report.txt"),
    log:
        os.path.join(
            dirs_dict["LOGS_DIR"],
            "gen_report_for_mapping_to_references_for_removal.log",
        ),
    threads: 1
    resources:
        nodes=1,
    run:
        M.gen_report_with_references_for_removal_info(input, output[0])
