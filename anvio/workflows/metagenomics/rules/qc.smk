# Metagenomics workflow QC rules — sourced from the shared QCModule.
#
# This file re-exports the QC rules so the metagenomics Snakefile include
# statement stays unchanged while rules live in the reusable qc/ module.

# Expose QCModule flags to the multiqc.smk scope
run_fastqc = M.get_param_value_from_config(["fastqc", "run"]) == True
run_filtlong  = M.get_param_value_from_config(["filtlong",  "run"]) == True
run_nanoplot  = M.get_param_value_from_config(["nanoplot",  "run"]) == True


include: w.get_workflow_rule_file_path("qc", "sr_filter.smk")

if M.has_lr:
    include: w.get_workflow_rule_file_path("qc", "lr.smk")

if M.get_param_value_from_config(["multiqc", "run"]) == True:
    include: w.get_workflow_rule_file_path("qc", "multiqc.smk")
