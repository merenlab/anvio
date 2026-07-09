# Metagenomics workflow QC rules — sourced from the shared QCModule.
#
# This file re-exports the QC rules so the metagenomics Snakefile include
# statement stays unchanged while rules live in the reusable qc/ module.

# Expose QCModule flags to the multiqc.smk scope
run_fastqc_sr = M.get_param_value_from_config(["fastqc_sr", "run"]) == True
run_filtlong  = M.get_param_value_from_config(["filtlong",  "run"]) == True


include: w.get_workflow_rule_file_path("qc", "sr_filter.smk")

if M.has_lr:
    include: w.get_workflow_rule_file_path("qc", "lr.smk")

if M.get_param_value_from_config(["multiqc", "run"]) == True:
    include: w.get_workflow_rule_file_path("qc", "multiqc.smk")
