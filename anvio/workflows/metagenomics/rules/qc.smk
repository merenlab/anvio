# Metagenomics workflow QC rules — sourced from the shared QCModule.
#
# This file re-exports the QC rules so the metagenomics Snakefile include
# statement stays unchanged while rules live in the reusable qc/ module.


include: w.get_workflow_rule_file_path("qc", "sr_filter.smk")
