# Rules that download sequencing reads from the SRA and hand them to the rest of the workflow.
#
# Expects the following in the including Snakefile scope:
#   M           — a SRAReadsModule-enabled workflow instance
#   dirs_dict   — M.dirs_dict
#   rule_log()  — canonical log-path helper
#   w           — anvio.workflows
#
# The reads produced here are wrapped in temp() unless the user asked to keep them, so snakemake
# deletes them as soon as the last rule that reads them is done. What keeps a thousand downloads
# from happening at once is not a resource limit (a resource is released when a job ends, not
# when its output is finally consumed) but the `gate` input below: it makes a download wait for
# an earlier sample's reads to have been released.

import glob
import shutil


# Paired-end and long-read runs are dumped by two different rules, because a rule's outputs have
# to be known from its wildcards alone and these two produce a different number of files. Each
# rule only ever matches the accessions of its own kind.
SRA_PAIRED_RE = w.regex_from_ids(sorted(M.get_accessions_of_read_type("PE-SR")))
SRA_LONG_READ_RE = w.regex_from_ids(sorted(M.get_accessions_of_read_type("LR")))
SRA_ACCESSION_RE = w.regex_from_ids(sorted(M.sra_metadata.keys()))
SRA_SAMPLE_RE = w.regex_from_ids(sorted(M.sra_accessions_by_sample.keys()))
SRA_UNIT_RE = w.regex_from_ids(sorted(M.download_units))


localrules: sra_reads_released


rule sra_prefetch:
    """Download the SRA archive of a single run.

    The `gate` input is what enforces the disk budget: when it is not empty, this download waits
    for an earlier release unit to be finished with its reads. The archive is declared as a temp
    directory so that it goes away as soon as the reads have been extracted from it."""
    input:
        gate=lambda wildcards: M.get_gate_flag_for_accession(wildcards.accession),
    output:
        prefetch_dir=temp(directory(os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}", "prefetch"))),
    log:
        rule_log("sra_prefetch", "{accession}"),
    wildcard_constraints:
        accession=SRA_ACCESSION_RE,
    threads: M.T("sra_prefetch")
    resources:
        nodes=M.T("sra_prefetch"),
    params:
        max_size=M.get_rule_param("sra_prefetch", "--max-size"),
    run:
        shell("prefetch {wildcards.accession} --output-directory {output.prefetch_dir} "
              "--verify yes {params.max_size} >> {log} 2>&1")

        if not M.find_sra_archive(output.prefetch_dir):
            raise ConfigError(f"`prefetch` finished working on the accession {wildcards.accession} without "
                              f"complaining, and yet there is no SRA archive to be found in "
                              f"'{output.prefetch_dir}'. The log file at '{log}' should say more about what "
                              f"happened.")


rule sra_fasterq_dump_paired:
    """Turn the SRA archive of a paired-end run into a pair of gzipped FASTQ files."""
    input:
        prefetch_dir=os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}", "prefetch"),
    output:
        r1=temp(os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}", "{accession}_1.fastq.gz")),
        r2=temp(os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}", "{accession}_2.fastq.gz")),
    log:
        rule_log("sra_fasterq_dump", "{accession}"),
    wildcard_constraints:
        accession=SRA_PAIRED_RE,
    threads: M.T("sra_fasterq_dump")
    resources:
        nodes=M.T("sra_fasterq_dump"),
    params:
        output_dir=os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}"),
    run:
        M.run_fasterq_dump(shell, wildcards.accession, input.prefetch_dir, params.output_dir,
                           threads, str(log), "PE-SR")


rule sra_fasterq_dump_long_reads:
    """Turn the SRA archive of a long-read run into a single gzipped FASTQ file."""
    input:
        prefetch_dir=os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}", "prefetch"),
    output:
        reads=temp(os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}", "{accession}.fastq.gz")),
    log:
        rule_log("sra_fasterq_dump", "{accession}"),
    wildcard_constraints:
        accession=SRA_LONG_READ_RE,
    threads: M.T("sra_fasterq_dump")
    resources:
        nodes=M.T("sra_fasterq_dump"),
    params:
        output_dir=os.path.join(dirs_dict["SRA_DIR"], "runs", "{accession}"),
    run:
        M.run_fasterq_dump(shell, wildcards.accession, input.prefetch_dir, params.output_dir,
                           threads, str(log), "LR")


rule sra_gather_reads:
    """Collect every run belonging to a sample into the one file the workflow expects.

    A sample sequenced across several runs is described by several accessions, and this is where
    they become one set of reads. Concatenating gzipped FASTQ files is safe: the result is a
    valid gzip stream that every downstream tool reads without noticing."""
    input:
        fastqs=lambda wildcards: M.get_gather_inputs(wildcards.sample, wildcards.kind),
    output:
        reads=M.temp_unless_reads_are_kept(
            os.path.join(dirs_dict["SRA_DIR"], "reads", "{sample}_{kind}.fastq.gz"), "raw"),
    log:
        rule_log("sra_gather_reads", "{sample}_{kind}"),
    wildcard_constraints:
        sample=SRA_SAMPLE_RE,
        kind="R1|R2|LR",
    threads: M.T("sra_gather_reads")
    resources:
        nodes=M.T("sra_gather_reads"),
    shell:
        "cat {input.fastqs} > {output.reads} 2>> {log}"


rule sra_reads_released:
    """Note that everything which needed a release unit's reads is done with them.

    This rule deliberately does not take the reads as input. By the time it runs snakemake has
    already deleted them, which is exactly the point: the marker it leaves behind is what allows
    the next download to start."""
    input:
        sinks=lambda wildcards: M.get_read_consumer_targets(wildcards.unit),
    output:
        flag=os.path.join(dirs_dict["SRA_DIR"], "released", "{unit}.released"),
    wildcard_constraints:
        unit=SRA_UNIT_RE,
    shell:
        "touch {output.flag}"
