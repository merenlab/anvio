# Read recruitment rules: indexing -> mapping -> BAM processing -> profiling -> merging
#
# Expects the following in scope:
#   M   — a ReadRecruitmentModule (or compatible) instance, already initialized
#   dirs_dict — M.dirs_dict (or equivalent)


import os
import argparse
import anvio
import anvio.utils as u
import anvio.workflows as w
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tables.miscdata import TableForLayerAdditionalData


run = terminal.Run()
progress = terminal.Progress()


# convenience sets:
SR_READSETS = M.get_sr_readset_ids()
LR_READSETS = M.get_lr_readset_ids()


def regex_from_ids(ids):
    return r"(?:{})".format("|".join(re.escape(x) for x in ids)) if ids else r"DO_NOT_MATCH"


SR_RS_RE = regex_from_ids(SR_READSETS)
LR_RS_RE = regex_from_ids(LR_READSETS)


# conda helper functions:
def conda_yaml(tool):
    p = M.get_param_value_from_config([tool, 'conda_yaml'])
    if not p: return None
    return os.path.abspath(os.path.expanduser(p))


def env_prefix(tool):
    name = M.get_param_value_from_config([tool, 'conda_env'])
    return f"conda run -n {name}" if name else ""


# bowtie2 index prefix
BT2_EXT = "bt2l" if "--large-index" in M.get_param_value_from_config(["bowtie_build", "additional_params"]) else "bt2"
BT2_PREFIX = os.path.join(dirs_dict["MAPPING_DIR"], "{group}", "{group}-contigs")


def use_post_ref_filter(wildcards):
    return M.remove_short_reads_based_on_references and (wildcards.group in M.references_for_removal)


############################## INDEXING ##############################


rule bowtie_build:
    """Build bowtie2 index for contigs fasta."""
    version: 1.0
    wildcard_constraints:
        readset = SR_RS_RE
    conda: conda_yaml('bowtie')
    log: dirs_dict["LOGS_DIR"] + "/{group}-bowtie_build.log"
    input:
        contigs = M.get_fasta
    output:
        idx = [ BT2_PREFIX + f".{i}.{BT2_EXT}" for i in (1,2,3,4) ] + \
              [ BT2_PREFIX + f".rev.1.{BT2_EXT}", BT2_PREFIX + f".rev.2.{BT2_EXT}" ]
    params:
        env_prefix = lambda wildcards: env_prefix('bowtie'),
        prefix = BT2_PREFIX,
        additional_params = M.get_param_value_from_config(["bowtie_build", "additional_params"])
    threads: M.T('bowtie_build')
    resources: nodes = M.T('bowtie_build')
    shell:
        r"""
        mkdir -p "$(dirname {params.prefix})"
        {params.env_prefix} bowtie2-build --threads {threads} {params.additional_params} \
            {input.contigs} {params.prefix} >> {log} 2>&1
        """


rule minimap2_index:
    version: 1.0
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-minimap2_index.log"
    wildcard_constraints:
        readset = LR_RS_RE
    conda: conda_yaml('minimap2')
    input:
        contigs = lambda wildcards: M.get_fasta(wildcards)
    output:
        idx = M.dirs_dict["MAPPING_DIR"] + "/{group}/{group}.mmi"
    params:
        env_prefix = lambda wildcards: env_prefix('minimap2'),
        additional_params = M.get_param_value_from_config(["minimap2_index", "additional_params"]),
        outdir = M.dirs_dict["MAPPING_DIR"] + "/{group}"
    threads: M.T('minimap2_index')
    resources: nodes = M.T('minimap2_index')
    shell:
        r"""
        mkdir -p {params.outdir}
        {params.env_prefix} minimap2 -d {output.idx} {params.additional_params} {input.contigs} >> {log} 2>&1
        """


############################## MAPPING ##############################


rule bowtie:
    """Map SR readsets with bowtie2 (extension handled via get_fastq)."""
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-{readset}-bowtie.log"
    wildcard_constraints:
        readset = SR_RS_RE
    conda: conda_yaml('bowtie')
    input:
        idx = rules.bowtie_build.output,
        r1  = lambda wildcards: M.get_fastq(wildcards.readset, pre_ref_removal=use_post_ref_filter(wildcards))['r1'],
        r2  = lambda wildcards: M.get_fastq(wildcards.readset, pre_ref_removal=use_post_ref_filter(wildcards))['r2'],
    output:
        sam = temp(dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.sam")
    params:
        r1 = lambda wildcards, input: ','.join(input.r1),
        r2 = lambda wildcards, input: ','.join(input.r2),
        env_prefix = lambda wildcards: env_prefix('bowtie'),
        index_prefix = BT2_PREFIX,
        additional_params = M.get_param_value_from_config(["bowtie", "additional_params"])
    threads: M.T('bowtie')
    resources: nodes = M.T('bowtie')
    shell:
        r"""
        {params.env_prefix} bowtie2 --threads {threads} \
                -x {params.index_prefix} \
                -1 {params.r1} -2 {params.r2} \
                {params.additional_params} \
                -S {output.sam} >> {log} 2>&1
        """


rule minimap2:
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-{readset}-minimap2.log"
    wildcard_constraints:
        readset = LR_RS_RE
    conda: conda_yaml('minimap2')
    input:
        idx = M.dirs_dict["MAPPING_DIR"] + "/{group}/{group}.mmi",
        reads = lambda wildcards: M.get_fastq(wildcards.readset)['lr']
    output:
        sam = temp(M.dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.sam")
    params:
        env_prefix = lambda wildcards: env_prefix('minimap2'),
        preset = M.get_param_value_from_config(["minimap2", "preset"]),
        additional_params = M.get_param_value_from_config(["minimap2", "additional_params"])
    threads: M.T("minimap2")
    resources: nodes = M.T('minimap2')
    shell:
        r"""
        {params.env_prefix} minimap2 -x {params.preset} -t {threads} -a {params.additional_params} {input.idx} {input.reads} -o {output.sam} 2>> {log}
        """


############################## BAM PROCESSING ##############################


rule samtools_view:
    """ sort sam file with samtools and create a RAW.bam file"""
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-{readset}-samtools_view.log"
    input:  M.dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.sam"
    params: additional_params = M.get_param_value_from_config(["samtools_view", "additional_params"])
    output: temp(dirs_dict["MAPPING_DIR"] + "/{group}/{readset}-RAW.bam")
    threads: M.T('samtools_view')
    resources: nodes = M.T('samtools_view'),
    shell: "samtools view -bS {input} -o {output} {params.additional_params} >> {log} 2>&1"


rule anvi_init_bam:
    """ run anvi-init-bam on RAW bam file to create a bam file ready for anvi-profile"""
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-{readset}-anvi_init_bam.log"
    input:  M.dirs_dict["MAPPING_DIR"] + "/{group}/{readset}-RAW.bam"
    output:
        bam = dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.bam",
        bai = dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.bam.bai"
    threads: M.T('anvi_init_bam')
    resources: nodes = M.T('anvi_init_bam'),
    shell: "anvi-init-bam {input} -o {output.bam} -T {threads} >> {log} 2>&1"


############################## PROFILING ##############################


def get_cluster_contigs_param(wildcards):
    if M.get_param_value_from_config(['anvi_profile', '--cluster-contigs']) != "":
        cluster_contigs = M.get_rule_param('anvi_profile', '--cluster-contigs')
    else:
        cluster_contigs = '--cluster-contigs' if M.group_sizes[wildcards.group] == 1 else ''
    return cluster_contigs


rule anvi_profile:
    """ run anvi-profile on the bam file"""
    version: anvio.__profile__version__
    log: dirs_dict["LOGS_DIR"] + "/{group}-{readset}-anvi_profile.log"
    input:
        bam = dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.bam",
        contigs = ancient(M.get_contigs_db_path())
    output:
        profile = dirs_dict["PROFILE_DIR"] + "/{group}/{readset}/PROFILE.db",
        AUXILIARY_DATA = dirs_dict["PROFILE_DIR"] + "/{group}/{readset}/AUXILIARY-DATA.db",
        runlog = dirs_dict["PROFILE_DIR"] + "/{group}/{readset}/RUNLOG.txt"
    params:
        output_dir = dirs_dict["PROFILE_DIR"] + "/{group}/{readset}",
        cluster_contigs = get_cluster_contigs_param,
        sample_name = lambda wildcards: f'--sample-name "{wildcards.readset}"',
        overwrite_output_destinations = "--overwrite-output-destinations",
        report_variability_full = M.get_rule_param("anvi_profile", "--report-variability-full"),
        skip_SNV_profiling = M.get_rule_param("anvi_profile", "--skip-SNV-profiling"),
        profile_SCVs = M.get_rule_param("anvi_profile", "--profile-SCVs"),
        description = M.get_rule_param("anvi_profile", "--description"),
        skip_hierarchical_clustering = M.get_rule_param("anvi_profile", "--skip-hierarchical-clustering"),
        distance = M.get_rule_param("anvi_profile", "--distance"),
        linkage = M.get_rule_param("anvi_profile", "--linkage"),
        min_contig_length = M.get_rule_param("anvi_profile", "--min-contig-length"),
        min_mean_coverage = M.get_rule_param("anvi_profile", "--min-mean-coverage"),
        min_coverage_for_variability = M.get_rule_param("anvi_profile", "--min-coverage-for-variability"),
        contigs_of_interest = M.get_rule_param("anvi_profile", "--contigs-of-interest"),
        queue_size = M.get_rule_param("anvi_profile", "--queue-size"),
        write_buffer_size = M.get_rule_param("anvi_profile", "--write-buffer-size"),
        write_buffer_size_per_thread = M.get_rule_param("anvi_profile", "--write-buffer-size-per-thread"),
        min_percent_identity = M.get_rule_param("anvi_profile", "--min-percent-identity"),
        fetch_filter = M.get_rule_param("anvi_profile", "--fetch-filter"),
        max_contig_length = M.get_rule_param("anvi_profile", "--max-contig-length"),
    threads: M.T('anvi_profile')
    resources: nodes = M.T('anvi_profile'),
    shell:
        """
            anvi-profile -i {input.bam} -c {input.contigs} -o {params.output_dir} \
                             {params.cluster_contigs} {params.min_contig_length} \
                             {params.sample_name} -T {threads} {params.overwrite_output_destinations} \
                             {params.profile_SCVs} {params.report_variability_full} \
                             {params.skip_SNV_profiling} {params.description} \
                             {params.skip_hierarchical_clustering} {params.distance} \
                             {params.linkage} {params.min_mean_coverage} \
                             {params.min_coverage_for_variability} {params.contigs_of_interest} \
                             {params.queue_size} {params.write_buffer_size} {params.write_buffer_size_per_thread} \
                             {params.min_percent_identity} {params.fetch_filter} \
                             {params.max_contig_length} >> {log} 2>&1
        """


############################## MERGING ##############################


def input_for_anvi_merge(wildcards):
    if M.get_param_value_from_config(['all_against_all']):
        member_readsets = M.get_readset_ids()
    else:
        member_readsets = M.get_readsets_for_mapping_to_group(wildcards.group)

    return [
        dirs_dict['PROFILE_DIR'] + f"/{wildcards.group}/{rs}/PROFILE.db"
        for rs in member_readsets
    ]


def configure_flag_files_for_anvi_merge_optional_inputs(wildcards, flag_file_name, run_flag):
    if M.get_param_value_from_config(['all_against_all']):
        member_readsets = M.get_readset_ids()
    else:
        member_readsets = M.get_readsets_for_mapping_to_group(wildcards.group)

    flag_files = [
        os.path.join(dirs_dict["PROFILE_DIR"], f"{wildcards.group}/{rs}", flag_file_name)
        for rs in member_readsets
    ]

    if not run_flag:
        flag_files = ancient(os.path.join(dirs_dict["CONTIGS_DIR"], f"{wildcards.group}.db"))

    return flag_files


def get_merge_optional_inputs(wildcards):
    d = {}
    for input_name, (flag_fn, run_flag) in M.get_merge_optional_inputs().items():
        d[input_name] = configure_flag_files_for_anvi_merge_optional_inputs(wildcards, flag_fn, run_flag)
    return d


def create_fake_output_files(_message, output):
    for o in output:
        with open(o, 'w') as f:
            f.write(_message + '\n')


def remove_empty_profile_databases(profiles, group):
    empty_profiles = []
    progress.new("Checking for empty profile databases")
    for p in profiles:
        keys, data = TableForLayerAdditionalData(argparse.Namespace(profile_db=p)).get(['total_reads_mapped'])
        if not next(iter(data.values()))['total_reads_mapped']:
            empty_profiles.append(p)
    profiles = list(set(profiles) - set(empty_profiles))
    progress.end()

    if not profiles:
        run.warning('It seems that all your profiles are empty for the '
                    'contigs database: %s.db. And so cannot be merged.' % group)

    run.info('Number of non-empty profile databases', len(profiles))
    run.info('Number of empty profile databases', len(empty_profiles))
    if len(empty_profiles) > 0:
        run.info('The following databases are empty: ', empty_profiles)
    return profiles


rule gen_readme_file_for_unmerged_groups:
    version: 1.0
    log: os.path.join(dirs_dict["LOGS_DIR"], "{group}-gen_readme_file_for_unmerged_groups.log")
    input:
        unpack(get_merge_optional_inputs),
        contigs = ancient(M.get_contigs_db_path()),
        profiles = input_for_anvi_merge,
    output: os.path.join(M.dirs_dict["MERGE_DIR"], "{group}", "README.txt")
    threads: M.T('gen_readme_file_for_unmerged_groups')
    resources: nodes = M.T('gen_readme_file_for_unmerged_groups')
    shell:
        """
            echo -e 'The group {wildcards.group} has only one sample. Hence, there is nothing to merge, but you can find\n\
                     the profile database here: {input.profiles}. Also, just so you know, profile was done using\n\
                     --cluster-contigs so you can visualize this profile database using anvi-interactive.' > {output} 2>>{log}
        """


rule anvi_merge:
    '''
        Run create a merged profile database.

        If there are multiple profiles mapped to the same contigs database,
        then merges these profiles. For individual profile, creates a symlink
        to the profile database. The purpose is to have one folder in
        which for every contigs database there is a profile database (or
        a symlink to a profile database) that could be used together for
        anvi-interactive.
    '''
    version: anvio.__profile__version__
    log: dirs_dict["LOGS_DIR"] + "/{group}-anvi_merge.log"
    input:
        unpack(get_merge_optional_inputs),
        contigs = ancient(M.get_contigs_db_path()),
        profiles = input_for_anvi_merge,
    output:
        profile = dirs_dict["MERGE_DIR"] + "/{group}/PROFILE.db",
        AUXILIARY_DATA = dirs_dict["MERGE_DIR"] + "/{group}/AUXILIARY-DATA.db",
        runlog = dirs_dict["MERGE_DIR"] + "/{group}/RUNLOG.txt"
    threads: M.T('anvi_merge')
    resources: nodes = M.T('anvi_merge'),
    params:
        output_dir = dirs_dict["MERGE_DIR"] + "/{group}",
        name = "{group}",
        profile_dir = dirs_dict["PROFILE_DIR"] + "/{group}",
        sample_name = M.get_rule_param("anvi_merge", "--sample-name"),
        description = M.get_rule_param("anvi_merge", "--description"),
        skip_hierarchical_clustering = M.get_rule_param("anvi_merge", "--skip-hierarchical-clustering"),
        enforce_hierarchical_clustering = M.get_rule_param("anvi_merge", "--enforce-hierarchical-clustering"),
        distance = M.get_rule_param("anvi_merge", "--distance"),
        linkage = M.get_rule_param("anvi_merge", "--linkage"),
        overwrite_output_destinations = "--overwrite-output-destinations",
    run:
        input.profiles = remove_empty_profile_databases(input.profiles, wildcards.group)

        if not input.profiles:
            _message = "Nothing to merge for %s. This should " \
                       "only happen if all profiles were empty " \
                       "(you can check the log file: {log} to see " \
                       "if that is indeed the case). " \
                       "This file was created just so that your workflow " \
                       "would continue with no error (snakemake expects " \
                       "to find these output files and if we don't create " \
                       "them, then it will be upset). As we see it, " \
                       "there is no reason to throw an error here, since " \
                       "you mapped your metagenome to some fasta files " \
                       "and you got your answer: whatever you have in " \
                       "your fasta file is not represented in your  " \
                       "metagenomes. Feel free to contact us if you think " \
                       "that this is our fault. sincerely, Meren Lab" \
                       % wildcards.group
            create_fake_output_files(_message, output)
            shell("echo %s >> {log}" % _message)

        elif M.group_sizes[wildcards.group] == 1:
            _message = "Only one file was profiled with %s so there \
                       is nothing to merge. But dont worry, you can \
                       still use anvi-interacite with the single profile \
                       database that is here: %s" \
                       % (wildcards.group, input.profiles[0])
            create_fake_output_files(_message, output)
            shell("echo %s >> {log}" % _message)

        elif len(input.profiles) == 1:
            _message = "Only one sample had reads recruited to %s " \
                       "and hence merging could not occur." \
                       % wildcards.group
            create_fake_output_files(_message, output)
            shell("echo %s >> {log}" % _message)

        else:
            shell("anvi-merge {input.profiles} -o {params.output_dir} -c {input.contigs} \
                   {params.sample_name} \
                   {params.overwrite_output_destinations} {params.description} \
                   {params.skip_hierarchical_clustering} {params.enforce_hierarchical_clustering} \
                   {params.distance} {params.linkage} >> {log} 2>&1")


############################## STATS ##############################


rule count_reads_in_fastq:
    """Count the reads in fastq(s) per readset."""
    version: 1.0
    log:
        dirs_dict["LOGS_DIR"] + "/{readset}-count_reads_in_fastq.log"
    input:
        unpack(lambda wildcards: M.get_fastq(wildcards.readset)),
    output:
        txt = dirs_dict["QC_DIR"] + "/{readset}-total_num_reads.txt"
    threads: 1
    resources: nodes = 1
    run:
        import subprocess, shlex, os

        os.makedirs(dirs_dict["QC_DIR"], exist_ok=True)

        fastq_files = []
        for k in ("lr", "r1", "r2"):
            if hasattr(input, k):
                fastq_files.extend(input[k])

        reads_total = 0
        for fastq in fastq_files:
            if fastq.endswith(".gz"):
                cmd = f"gunzip -c {shlex.quote(fastq)} | awk 'END{{print NR/4}}'"
            else:
                cmd = f"awk 'END{{print NR/4}}' {shlex.quote(fastq)}"
            out = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.PIPE)
            reads_total += int(float(out.strip() or 0))

        with open(output.txt, "w") as f:
            f.write(f"{reads_total}\n")


run_import_percent_of_reads_mapped = M.get_param_value_from_config(['import_percent_of_reads_mapped', 'run']) == True


rule import_percent_of_reads_mapped:
    """Mapper-agnostic import of mapping stats into layers."""
    version: 1.0
    log:
        dirs_dict["LOGS_DIR"] + "/{group}-{readset}-import_percent_of_reads_mapped.log"
    input:
        total_reads = dirs_dict["QC_DIR"] + "/{readset}-total_num_reads.txt",
        profiledb = dirs_dict["PROFILE_DIR"] + "/{group}/{readset}/PROFILE.db"
    output:
        layers_txt = dirs_dict["PROFILE_DIR"] + "/{group}/{readset}/layers-additional-data.txt",
    params:
        bam = dirs_dict["MAPPING_DIR"] + "/{group}/{readset}.bam",
    threads: 1
    resources: nodes = 1
    run:
        import subprocess, shlex

        with open(input.total_reads) as f:
            reads_total = int(f.read().strip() or 0)

        cmd_mapped = f"samtools view -c -F 2308 {shlex.quote(params.bam)}"
        mapped_out = subprocess.check_output(cmd_mapped, shell=True, text=True, stderr=subprocess.PIPE)
        reads_mapped = int(mapped_out.strip() or 0)

        pct = (100 * reads_mapped / reads_total) if reads_total > 0 else 0

        with open(output.layers_txt, "w") as f:
            f.write("layers\ttotal_num_reads\ttotal_unique_reads_mapped\tpercent_mapped\n")
            f.write(f"{wildcards.readset}\t{reads_total}\t{reads_mapped}\t{pct:.2f}\n")

        shell(
            "anvi-import-misc-data {output.layers_txt} "
            "-p {input.profiledb} "
            "--target-data-table layers >> {log} 2>&1"
        )
