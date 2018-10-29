#!/bin/bash
source 00.sh
set -e

if [ $# -eq 0  ]
then
      echo "\

        No arguments supplied. You must call this script with one argument: 'new',
        or 'continue':

        'new'     : generate short reads, map them to the contig, gen contigs.db, run profiling.
        'continue': re-run profiling.
        "
        exit -1
fi

if [ $# -gt 1  ]
then
      echo "
        This scripts expect only one argument ('new', or 'continue').
        "
        exit -1
fi

if [ $1 = "new"  ]
then
    INFO "Creating the output directory"
    # change directory and clean the old mess if it exists
    cd sandbox
    rm -rf test-output
    mkdir test-output
    
    INFO "Anvo'o version"
    anvi-profile --version
    
    INFO "Generating a Bowtie2 ref"
    bowtie2-build single_contig.fa test-output/single_contig

    INFO "Generating anvi'o contigs database"
    anvi-gen-contigs-database -f single_contig.fa -o test-output/single_contig.db -L 60 -n "Single contig"

    for sample in 01 02 03 04 05
    do
        INFO "Generating short reads for sample $sample"
        anvi-script-gen-short-reads config_files_for_variability_testing/$sample.ini --output-file-path test-output/$sample.fa

        INFO "Mapping short reads to the ref"
        ../misc/bowtie_batch_single_fasta.sh test-output/$sample.fa test-output/$sample test-output/single_contig

        INFO "Profiling the BAM file"
        anvi-profile -i test-output/$sample.bam -c test-output/single_contig.db -o test-output/$sample-PROFILE -M 0 --profile-SCVs
    done
    

    INFO "Merging all"
    anvi-merge test-output/*/PROFILE.db -c test-output/single_contig.db -o test-output/SAMPLES-MERGED --skip-concoct-binning

elif [ $1 = "continue"  ]
then
    cd sandbox

    if [ ! -f "test-output/single_contig.db"  ]
    then
      echo "
        You asked to continue with the previously generated files,
        but the contigs database is not there. Maybe you should start
        from scratch by re-running this script with the parameter 'new'.
        "
        exit -1
    else
        INFO "Attempting to continue with the previously generated files"

        rm -rf test-output/*-PROFILE
        rm -rf test-output/*-MERGED

        for sample in 01 02 03 04 05
        do
            INFO "Re-profiling the BAM file"
            anvi-profile -i test-output/$sample.bam -c test-output/single_contig.db -o test-output/$sample-PROFILE -M 0 --profile-SCVs
        done

        INFO "Merging all"
        anvi-merge test-output/*/PROFILE.db -c test-output/single_contig.db -o test-output/SAMPLES-MERGED --skip-concoct-binning
    fi
else
      echo "
        Unknown parameter $1 :/ Try 'new', or 'continue'.
        "
        exit -1
fi

INFO "Importing a default state"
anvi-import-state -s config_files_for_variability_testing/default.json -n default -p test-output/SAMPLES-MERGED/PROFILE.db

# splits of interest:
cat << EOF > test-output/splits_of_interest.txt
recA_split_00001
recA_split_00002
recA_split_00003
EOF

INFO "anvi-gen-variability for NT"
anvi-gen-variability-profile -p test-output/SAMPLES-MERGED/PROFILE.db \
                             -c test-output/single_contig.db \
                             -o test-output/variability_NT.txt \
                             --splits-of-interest test-output/splits_of_interest.txt \
                             --quince-mode \
                             --engine NT

column -t test-output/variability_NT.txt | head

INFO "anvi-gen-variability for AA"
anvi-gen-variability-profile -p test-output/SAMPLES-MERGED/PROFILE.db \
                             -c test-output/single_contig.db \
                             -o test-output/variability_AA.txt \
                             --splits-of-interest test-output/splits_of_interest.txt \
                             --quince-mode \
                             --engine AA

column -t test-output/variability_AA.txt | head

INFO "anvi-gen-variability for CDN"
anvi-gen-variability-profile -p test-output/SAMPLES-MERGED/PROFILE.db \
                             -c test-output/single_contig.db \
                             -o test-output/variability_CDN.txt \
                             --gene-caller-ids 0 \
                             --quince-mode \
                             --engine CDN

column -t test-output/variability_CDN.txt | head


INFO "Do you want thhe interactive interface? Run the following:"

echo "anvi-interactive -p `pwd`/test-output/SAMPLES-MERGED/PROFILE.db -c `pwd`/test-output/single_contig.db"
echo
echo
