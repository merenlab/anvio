#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

cmd="-pn"
# if you want the test to actually run through the pipeline
# then call it like this: bash run_merenlab_metagenomics_pipeline_test.sh sandbox/test-output full
if [ $# -ge 2 ]; then
    if [ $2 == "full" ]; then
        cmd=$3
    fi
fi
# copy latest script here
cp ../merenlab-metagenomics-pipeline.snakefile $output_dir
cp -R ../wrappers/ $output_dir/wrappers/
cp ../mock_files_for_merenlab_metagenomics_pipeline/*json $output_dir
cp -R ../mock_files_for_merenlab_metagenomics_pipeline/three_samples_example/ $output_dir/three_samples_example/
cp ../mock_files_for_merenlab_metagenomics_pipeline/samples*.txt $output_dir
cp ../mock_files_for_merenlab_metagenomics_pipeline/references.txt $output_dir

# we have to go into the test directory because snakemake requires you run the command from the directory where the snakemake is
cd $output_dir


INFO "Call snakefile with megahit" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd


INFO "Call snakefile with idba_ud" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config output_dirs='{"MERGE_DIR": "06_MERGED_idba_ud"}' \
          idba_ud='{"run": True}' \
          megahit='{"run": False}'



INFO "Call snakefile with all against all"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_ALL_AGAINST_ALL"}'


INFO "Call snakefile with all against all with no qc" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_ALL_AGAINST_ALL_USING_RAW_INPUTS"}' \
          qc='{"run": False}'


INFO "make a copy of all fastq.gz files"
for f in `ls $output_dir/three_samples_example/*.fastq.gz`; do
    s=$output_dir/three_samples_example/`echo $f | rev | cut -f 1 -d \/ | rev`
    cp $f ${s%.fastq.gz}-for-idba_ud.fastq.gz
done
INFO "uncompress all copied fastq.gz files"
gzip -d $output_dir/three_samples_example/*-for-idba_ud*gz

INFO "Call snakefile with idba_ud with no qc" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_ALL_AGAINST_ALL_USING_RAW_INPUTS"}' \
          samples_txt='samples-for-idba_ud.txt' \
          qc='{"run": False}' \
          idba_ud='{"run": True}' \
          megahit='{"run": False}'


INFO "decompress mock reference files"
gzip -d three_samples_example/*.fa.gz 


INFO "Call snakefile with group list" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt' \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE"}' \
          samples_txt='samples-no-groups.txt'


INFO "Call snakefile with group list with all against all" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt'\
          all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE_all_against_all"}'


INFO "Call snakefile with no group list in reference mode"
INFO "If you run this test in full mode then this one shouldn't do anything and just say 'Nothing to be done.'" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt' \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE_all_against_all"}' \
          samples_txt='samples-no-groups.txt'


INFO "Call snakefile with group list with all against all with no qc and no reformat_fasta" $2
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt'\
          all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE_all_against_all_USING_RAW_INPUTS"}' \
          qc='{"run": False}' \
          reformat_fasta='{"run": False}'

INFO "All tests finished succesfully"
INFO "Here is a list of all the rules in the workflow:"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile -l

# go back to the directory where we started
cd -

INFO "clear all files"
rm -rf $output_dir
