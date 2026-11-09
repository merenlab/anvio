def fastq_candidates(accession, directory):
    """Return the exact FASTQ file names `fasterq-dump --split-3` can produce for an accession.

    Matching these names exactly rather than globbing on the accession prefix is what keeps
    accessions whose names are prefixes of one another (such as SRR1195143 and SRR11951439)
    from claiming each other's FASTQ files.
    """
    return [
        os.path.join(directory, f"{accession}{suffix}.fastq")
        for suffix in ["", "_1", "_2", "_3"]
    ]


rule prefetch:
    """Prefetch data from the Sequence Read Archive (SRA).

Inputs:
    None

Outputs:
    prefetch_DONE: the empty marker file {accession}.prefetch.done

    SRA: File in the SRA_prefetch directory with the name {accession}.sra or {accession}.sralite. This file
    is not used to control the execution of subsequent rules

Params:
    SRA_OUTPUT_DIR: Output directory for the prefetched data
    MAX_SIZE: Maximum size of the SRA file to download

Threads:
    The number of threads to use is specified by the prefetch variable

NOTES:
- This is the first rule of the workflow
"""
    output:
        prefetch_DONE=touch(
            os.path.join(
                dirs_dict["SRA_prefetch"], "{accession}", "{accession}.prefetch.done"
            )
        ),
    log:
        rule_log("prefetch", "{accession}_prefetch"),
    wildcard_constraints:
        accession=ACCESSION_RE,
    threads: M.T("prefetch")
    resources:
        nodes=M.T("prefetch"),
    params:
        SRA_OUTPUT_DIR=os.path.join(dirs_dict["SRA_prefetch"]),
        MAX_SIZE=M.get_param_value_from_config(["prefetch", "--max-size"]),
    run:
        shell(
            "prefetch {wildcards.accession} --output-directory {params.SRA_OUTPUT_DIR} --max-size {params.MAX_SIZE} >> {log} 2>&1"
        )
        sra_file = os.path.join(
            dirs_dict["SRA_prefetch"],
            f"{wildcards.accession}",
            f"{wildcards.accession}.sra",
        )
        sralite_file = os.path.join(
            dirs_dict["SRA_prefetch"],
            f"{wildcards.accession}",
            f"{wildcards.accession}.sralite",
        )
        if not os.path.exists(sra_file) and not os.path.exists(sralite_file):
            raise ConfigError(
                f"Following the execution of `prefetch` for accession {wildcards.accession}, we could not find either of "
                f"the expected output files {sra_file} or {sralite_file}."
            )


rule check_md5sum:
    """Curl the md5sum file from the SRA FTP site"""
    input:
        prefetch_DONE=os.path.join(
            dirs_dict["SRA_prefetch"], "{accession}", "{accession}.prefetch.done"
        ),
    output:
        md5sum_DONE=touch(
            os.path.join(
                dirs_dict["SRA_prefetch"], "{accession}", "{accession}.md5sum.done"
            )
        ),
    log:
        rule_log("check_md5sum", "{accession}_check_md5sum"),
    wildcard_constraints:
        accession=ACCESSION_RE,
    threads: M.T("check_md5sum")
    resources:
        nodes=M.T("check_md5sum"),
    params:
        md5sum=os.path.join(
            dirs_dict["SRA_prefetch"], "{accession}", "{accession}.json"
        ),
        sra_file=os.path.join(
            dirs_dict["SRA_prefetch"], "{accession}", "{accession}.sra"
        ),
    run:
        # Identify the correct SRA file from prefetch. `prefetch` gives us either a full
        # `.sra` or a reduced-representation `.sralite`, and NCBI reports a separate md5
        # for each, so the file type we downloaded decides which entry we compare against.
        accession = wildcards.accession
        sra_file = params.sra_file
        file_type = "sra"
        sralite_file = os.path.join(
            dirs_dict["SRA_prefetch"], accession, f"{accession}.sralite"
        )
        if not os.path.exists(sra_file) and os.path.exists(sralite_file):
            sra_file = sralite_file
            file_type = "sralite"
        log_path = str(log)
        with open(log_path, "w") as log_file:
            log_file.write(f"Checking MD5sum for file: {sra_file}\n")
        # Get the md5sum from the SRA FTP site
        shell(
            "curl 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&location-type=forced&location=s3.us-east-1&accept-charges=aws&acc={wildcards.accession}' --output {params.md5sum} >> {log} 2>&1"
        )
        # Get the expected md5sum hash. The response is a JSON document that looks like
        # {"result": [{"bundle": ACC, "status": 200, "files": [{"type": "sra", "md5": ...}]}]},
        # but NCBI also returns non-200 bundles for suppressed or restricted runs, multiple
        # file entries for some accessions, and entries that carry no md5 at all.
        try:
            with open(params.md5sum) as sra_metadata_file:
                sra_metadata_dict = json.load(sra_metadata_file)
            bundles = sra_metadata_dict["result"]
        except (OSError, ValueError, KeyError, TypeError) as e:
            raise ConfigError(
                f"Anvi'o asked NCBI for the md5 checksum of the accession {accession}, but could not make "
                f"sense of the answer it got back. This usually means the request failed or NCBI returned "
                f"something unexpected rather than its usual JSON document. The raw response is at "
                f"'{params.md5sum}' and the log at '{log_path}' if you wish to take a look. This is what "
                f"we know: {e}."
            )

        bundle = next((b for b in bundles if b.get("bundle") == accession), None)
        if bundle is None and len(bundles) == 1:
            bundle = bundles[0]
        if bundle is None:
            raise ConfigError(
                f"NCBI did not return any information for the accession {accession} when anvi'o asked for "
                f"its md5 checksum, even though `prefetch` was able to download it. You can see the raw "
                f"response at '{params.md5sum}'."
            )
        if bundle.get("status") != 200:
            raise ConfigError(
                f"NCBI refused to give anvi'o the md5 checksum for the accession {accession}, and said this "
                f"about it: \"{bundle.get('msg', 'no message')}\" (status code: {bundle.get('status')}). "
                f"Accessions that are suppressed by NCBI or that require special access permissions often "
                f"end up here."
            )

        files = bundle.get("files") or []
        entry = next((f for f in files if f.get("type") == file_type), None)
        expected_md5 = entry.get("md5") if entry else None

        if not expected_md5:
            run.warning(
                f"NCBI did not report an md5 checksum for the '{file_type}' file of the accession "
                f"{accession}, so anvi'o has no way to verify that this download is intact and will "
                f"continue without checking it. Everything will likely be fine, but if you see strange "
                f"things downstream for this accession, a truncated download is a suspect worth ruling out."
            )
            with open(log_path, "a") as log_file:
                log_file.write("No md5 reported by NCBI for this file; verification skipped.\n")
        else:
            calculated_md5 = M.calculate_md5(file_path=sra_file)
            if calculated_md5 != expected_md5:
                raise ConfigError(
                    f"The md5 checksum of the file anvi'o downloaded for the accession {accession} does not "
                    f"match the one NCBI reports for it. NCBI expects '{expected_md5}' but the file at "
                    f"'{sra_file}' hashes to '{calculated_md5}'. This means the download is incomplete or "
                    f"corrupt. Removing that file and running the workflow again will download it afresh."
                )
            with open(log_path, "a") as log_file:
                log_file.write(f"Checksums match: {calculated_md5}\n")


rule fasterq_dump:
    """Use fasterq-dump to extract FASTQ file(s) from an SRA prefetch *.sra
Using the flag --split-3. For more information about fasterq-dump:
https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

Inputs:
- SRA file from the output of the prefetch rule

Outputs:
- Three possible FASTQ files: single reads, and/or R1 and R2.

Params:
- SRA_INPUT_DIR: directory containing the SRA file
- OUTPUT_DIR: directory to write the FASTQ file(s)

Threads:
- Number of threads specified in the M object
"""
    input:
        done=ancient(rules.check_md5sum.output.md5sum_DONE),
    output:
        FASTERQDUMP_DONE=touch(
            os.path.join(dirs_dict["FASTQ_DIR"], "{accession}-fasterq-dump.done")
        ),
        FASTERQDUMP_TEMP=temp(directory("FASTERQDUMP_TEMP/{accession}")),
    log:
        rule_log("fasterq_dump", "{accession}_fasterq_dump"),
    wildcard_constraints:
        accession=ACCESSION_RE,
    threads: M.T("fasterq_dump")
    resources:
        nodes=M.T("fasterq_dump"),
    params:
        SRA_INPUT_DIR=os.path.join(dirs_dict["SRA_prefetch"], "{accession}"),
        OUTPUT_DIR=dirs_dict["FASTQ_DIR"],
        REMOVE_SRA_FILES=remove_tmp,
    run:
        shell(
            "fasterq-dump {params.SRA_INPUT_DIR} -t {output.FASTERQDUMP_TEMP} --outdir {params.OUTPUT_DIR} --split-3 --verbose --progress --threads {threads} >> {log} 2>&1"
        )
        # Check if fasterq-dump encountered a Disk quota exceeded error
        error_message = "Disk quota exceeded"
        log_path = str(log)
        with open(log_path, "r") as log_file:
            log_contents = log_file.read()
            if error_message in log_contents:
                raise ConfigError(
                    f"fasterq-dump ran out of disk space while processing the accession "
                    f"{wildcards.accession}. You can see the details in the log file at '{log_path}'."
                )

        produced = [
            f
            for f in fastq_candidates(wildcards.accession, params.OUTPUT_DIR)
            if os.path.exists(f)
        ]
        if not produced:
            raise ConfigError(
                f"fasterq-dump finished working on the accession {wildcards.accession} without complaining, "
                f"yet anvi'o cannot find a single FASTQ file for it in '{params.OUTPUT_DIR}'. The log file "
                f"at '{log_path}' should have more to say about what happened here."
            )

        # The prefetched SRA archive has served its purpose once the reads are extracted from it.
        # Only the archive itself goes, so the marker files snakemake tracks in the same directory
        # survive and the workflow does not download this accession again on a subsequent run.
        if params.REMOVE_SRA_FILES:
            for extension in [".sra", ".sralite"]:
                sra_file = os.path.join(
                    dirs_dict["SRA_prefetch"],
                    wildcards.accession,
                    f"{wildcards.accession}{extension}",
                )
                if os.path.exists(sra_file):
                    os.remove(sra_file)


rule pigz:
    """Compress FASTQ file(s) using pigz in parallel!

Inputs:
- FASTQ file(s) from the output of the fasterq_dump rule

Outputs:
- gzipped FASTQ file(s) in the FASTQ directory

Params:
- READS: the exact FASTQ file names this accession can have in the FASTQ directory

Threads:
- Number of threads specified in the M object

example:
    pigz --processes 8 --verbose 02_FASTQ/ERR6450080_1.fastq >> 00_LOGS/ERR6450080_pigz.log 2>&1

"""
    input:
        done=ancient(rules.fasterq_dump.output.FASTERQDUMP_DONE),
    output:
        done=touch(os.path.join(dirs_dict["FASTQ_DIR"], "{accession}-pigz.done")),
    log:
        rule_log("pigz", "{accession}_pigz"),
    wildcard_constraints:
        accession=ACCESSION_RE,
    threads: M.T("pigz")
    resources:
        nodes=M.T("pigz"),
    params:
        READS=lambda wildcards: " ".join(
            fastq_candidates(wildcards.accession, dirs_dict["FASTQ_DIR"])
        ),
    shell:
        """
        for f in $(ls -1 {params.READS} 2>/dev/null || true) ; do
            pigz --processes {threads} --verbose "$f" >> {log} 2>&1 || exit 1
        done
        """


rule generate_samples_txt:
    """Create samples-txt files for paired-end reads and single reads samples

Inputs:
- gziped FASTQ file(s) from the output of the pigz rule

Outputs:
- samples.txt and samples_single_reads.txt in the FASTQ directory

Params:
- ACCESSION: list of accessions
- OUTPUT_DIR: directory to write the samples.txt file

Threads:
- Number of threads specified in the M object
"""
    input:
        targets=expand(
            os.path.join(dirs_dict["FASTQ_DIR"], "{accession}-pigz.done"),
            accession=M.accessions_list,
        ),
    output:
        done=touch(os.path.join(dirs_dict["FASTQ_DIR"], "generate_samples_txt.done")),
    log:
        rule_log("generate_samples_txt", "generate_samples_txt"),
    threads: 1
    resources:
        nodes=1,
    params:
        ACCESSION=M.accessions_list,
        OUTPUT_DIR=dirs_dict["FASTQ_DIR"],
    run:
        paired_reads = []
        single_reads = []
        for sample in params.ACCESSION:
            if os.path.exists(
                os.path.join(params.OUTPUT_DIR, "".join([sample, "_1.fastq.gz"]))
            ):
                paired_reads.append(sample)
            elif os.path.exists(
                os.path.join(params.OUTPUT_DIR, "".join([sample, ".fastq.gz"]))
            ):
                single_reads.append(sample)
            else:
                raise ConfigError(
                    f"Looks like sample {sample} doesn't have the expected output(s) format."
                )
        if paired_reads:
            with open("samples.txt", "w") as f:
                f.write("sample\tr1\tr2\n")
                for sample in paired_reads:
                    r1 = os.path.join(
                        os.getcwd(),
                        params.OUTPUT_DIR,
                        "".join([sample, "_1.fastq.gz"]),
                    )
                    r2 = os.path.join(
                        os.getcwd(),
                        params.OUTPUT_DIR,
                        "".join([sample, "_2.fastq.gz"]),
                    )
                    f.write("%s\t%s\t%s\n" % (sample, r1, r2))
        if single_reads:
            with open("samples_single_reads.txt", "w") as f:
                f.write("sample\tr1\n")
                for sample in single_reads:
                    r1 = os.path.join(
                        os.getcwd(),
                        params.OUTPUT_DIR,
                        "".join([sample, ".fastq.gz"]),
                    )
                    f.write("%s\t%s\n" % (sample, r1))
