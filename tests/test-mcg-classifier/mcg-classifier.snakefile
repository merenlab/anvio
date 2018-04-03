'''
    A test file to test mcg-classifier

    To run it:
        snakemake -s mcg-classifier.snakefile -p
'''

configfile: "config.json"

def A(_list, d, default_value = ""):
    '''
        A helper function to make sense of config details.
        string_list is a list of strings (or a single string)
        d is a dictionary

        this function checks if the strings in x are nested values in y.
        For example if x = ['a','b','c'] then this function checkes if the
        value y['a']['b']['c'] exists, if it does then it is returned
    '''
    if type(_list) is not list:
        # converting to list for the cases of only one item
        _list = [_list]
    while _list:
        a = _list.pop(0)
        if a in d:
            d = d[a]
        else:
            return default_value
    return d


output_dir = A("output_dir", config, "test-output")
files_dir = A("files_dir", config, "../sandbox/mock_files_for_alons_classifier")
samples = ["hmp0041", "hmp0062", "hmp0074", "hmp0075", "hmp0079", "hmp0094"]

rule all:
    input: output_dir + "/mcg.finished"


rule gen_contigs_db:
    """ Generates a contigs database using anvi-gen-contigs-database """
    log: output_dir + "/TEST-gen_contigs_db.log"
    # depending on whether human contamination using centrifuge was done
    # or not, the input to this rule will be the raw assembly or the
    # filtered.
    input: files_dir + "/TEST.fa"
    output:
        db = output_dir + "/TEST.db",
    shell: "anvi-gen-contigs-database -f {input} -o {output.db} -n TEST >> {log} 2>&1"


rule anvi_init_bam:
    log: output_dir + "/TEST-{sample}-anvi_init_bam.log"
    input: files_dir + "/{sample}-RAW.bam"
    output:
        bam = output_dir + "/TEST/{sample}.bam",
        bai = output_dir + "/TEST/{sample}.bam.bai"
    shell: "anvi-init-bam {input} -o {output.bam} >> {log} 2>&1"


rule profile:
    log: output_dir + "/TEST-{sample}-anvi_profile.log"
    input:
        bam = output_dir + "/TEST/{sample}.bam",
        contigs = output_dir + "/TEST.db"
    output:
        profile = output_dir + "/TEST/{sample}/PROFILE.db",
        aux = output_dir + "/TEST/{sample}/AUXILIARY-DATA.db",
        runlog = output_dir + "/TEST/{sample}/RUNLOG.txt"
    params:
        name = "{sample}",
        output_dir = output_dir + "/TEST/{sample}"
    shell: "anvi-profile -c {input.contigs} -i {input.bam} -o {params.output_dir} -S {params.name} --overwrite-output-destinations"

rule merge:
    log: output_dir + "/TEST-anvi_merge.log"
    input:
        profiles = expand(output_dir + "/TEST/{sample}/PROFILE.db", sample=samples),
        contigs = rules.gen_contigs_db.output.db
    output:
        profile = output_dir + "/TEST/MERGED-SAMPLES/PROFILE.db",
        aux = output_dir + "/TEST/MERGED-SAMPLES/AUXILIARY-DATA.db",
        runlog = output_dir + "/TEST/MERGED-SAMPLES/RUNLOG.txt"
    params:
        output_dir = output_dir + "/TEST/MERGED-SAMPLES",
        name = "TEST"
    shell: "anvi-merge {input.profiles} -o {params.output_dir} -c {input.contigs} -S {params.name} --overwrite-output-destinations --skip-concoct-binning >> {log} 2>&1"


rule import_collection:
    log: output_dir + "/TEST-import_collection.log"
    input:
        profile = output_dir + "/TEST/MERGED-SAMPLES/PROFILE.db",
        contigs = output_dir + "/TEST.db"
    output: touch(output_dir + "/import_collection.done")
    params: collection = files_dir + "/TEST-COLLECTION.txt"
    shell: "anvi-import-collection -c {input.contigs} -p {input.profile} {params.collection} -C TEST >> {log} 2>&1"


rule run_mcg_classifier:
    log: output_dir + "/TEST-run_mcg_classifier.log"
    input:
        profile = output_dir + "/TEST/MERGED-SAMPLES/PROFILE.db",
        contigs = output_dir + "/TEST.db"
    output: touch(output_dir + "/mcg.finished")
    params:
        output_prefix= output_dir + "/TEST"
    shell: "anvi-mcg-classifier -p {input.profile} -c {input.contigs} -O {params.output_prefix} --outliers-threshold 1.5 --alpha 0.05 --gen-figures --zeros-are-outliers -C TEST -b Bin_1 -W --debug >> {log} 2>&1"
