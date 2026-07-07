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
        # Identify the correct SRA file from prefetch
        sralite_file = os.path.join(
            dirs_dict["SRA_prefetch"],
            f"{wildcards.accession}",
            f"{wildcards.accession}.sralite",
        )
        if not os.path.exists(params.sra_file) and os.path.exists(sralite_file):
            params.sra_file = sralite_file
        log_path = str(log)
        with open(log_path, "w") as log_file:
            log_file.write(f"Checking MD5sum for file: {params.sra_file}\n")
        # Get the md5sum from the SRA FTP site
        shell(
            "curl 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&location-type=forced&location=s3.us-east-1&accept-charges=aws&acc={wildcards.accession}' --output {params.md5sum} >> {log} 2>&1"
        )
        # Get the expected md5sum hash
        with open(params.md5sum) as sra_metadata_file:
            sra_metadata_dict = json.load(sra_metadata_file)
            expected_md5 = sra_metadata_dict["result"][0]["files"][0]["md5"]
        calculated_md5 = M.calculate_md5(file_path=params.sra_file)
        # Compare expected and calculated md5 and log the result
        with open(log_path, "a") as log_file:
            if calculated_md5 == expected_md5:
                log_file.write(f"Checksums match: {calculated_md5}\n")
            else:
                raise ValueError(
                    f"[{wildcards.accession}] MD5 checksum failed: "
                    f"expected {expected_md5}, got {calculated_md5} for file {params.sra_file}"
                )


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
            os.path.join(dirs_dict["FASTAS"], "{accession}-fasterq-dump.done")
        ),
        FASTERQDUMP_TEMP=temp(directory("FASTERQDUMP_TEMP/{accession}")),
    log:
        rule_log("fasterq_dump", "{accession}_fasterq_dump"),
    threads: M.T("fasterq_dump")
    resources:
        nodes=M.T("fasterq_dump"),
    params:
        SRA_INPUT_DIR=os.path.join(dirs_dict["SRA_prefetch"], "{accession}"),
        OUTPUT_DIR=dirs_dict["FASTAS"],
        FASTERQDUMP_TEMP=temp(os.path.join(os.getcwd(), "FASTERQDUMP_TEMP")),
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
                raise Exception(
                    "fasterq-dump encountered a Disk quota exceeded when processing {wildcards.accession}"
                )


rule pigz:
    """Compress FASTQ file(s) using pigz in parallel!

Inputs:
- FASTQ file(s) from the output of the fasterq_dump rule

Outputs:
- gzipped FASTQ file(s) in the FASTAS directory

Params:
- READS: prefix of the FASTQ file(s) in the FASTAS directory

Threads:
- Number of threads specified in the M object

example:
    pigz --processes 8 --verbose 02_FASTA/ERR6450080*.fastq >> 00_LOGS/ERR6450080_pigz.log 2>&1

"""
    input:
        done=ancient(rules.fasterq_dump.output.FASTERQDUMP_DONE),
    output:
        done=touch(os.path.join(dirs_dict["FASTAS"], "{accession}-pigz.done")),
    log:
        rule_log("pigz", "{accession}_pigz"),
    threads: M.T("pigz")
    resources:
        nodes=M.T("pigz"),
    params:
        READS=os.path.join(dirs_dict["FASTAS"], "{accession}*.fastq"),
    shell:
        "pigz --processes {threads} --verbose {params.READS} >> {log} 2>&1"


rule generate_samples_txt:
    """Create samples-txt files for paired-end reads and single reads samples

Inputs:
- gziped FASTQ file(s) from the output of the pigz rule

Outputs:
- samples.txt and samples_single_reads.txt in the FASTAS directory

Params:
- ACCESSION: list of accessions
- OUTPUT_DIR: directory to write the samples.txt file

Threads:
- Number of threads specified in the M object
"""
    input:
        targets=expand(
            os.path.join(dirs_dict["FASTAS"], "{accession}-pigz.done"),
            accession=M.accessions_list,
        ),
    output:
        done=touch(os.path.join(dirs_dict["FASTAS"], "generate_samples_txt.done")),
    log:
        rule_log("generate_samples_txt", "generate_samples_txt"),
    threads: 1
    resources:
        nodes=1,
    params:
        ACCESSION=M.accessions_list,
        OUTPUT_DIR=dirs_dict["FASTAS"],
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
