# change directory and clean the old mess if it exists
cd mini_test
rm -rf test-output
mkdir test-output


# init raw bam files.
for f in 6M 7M 9M
do
    papi-init-bam 204_3contigs_"$f".bam -o test-output/204-$f
done


# generate an annotation db using files obtained from myrast_gui using contigs.fa (contigs.fa
# is the original file all samples were mapped to) using split size 1000 (the default split
# size is better for most projects, small split size here is for testing purposes)
papi-gen-annotation --contigs contigs.fa -p myrast_gui myrast_gui/* -o test-output/ -L 1000 


# for each sample, run the profiling using the same split size used for the annotation database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 6M 7M 9M
do
    papi-profile -i test-output/204-$f.bam -o test-output/204-$f -a test-output/ANNOTATION.db -L 1000
done


# merge samples
papi-merge-multiple-runs test-output/204*/RUNINFO.cp -o test-output/204-MERGED


# generate gene and function networks for the merge
papi-gen-network test-output/204-MERGED/RUNINFO.mcp test-output/ANNOTATION.db


# fire up the browser to show how does the merged samples look like.
papi-interactive-binning -r test-output/204-MERGED/RUNINFO.mcp -a test-output/ANNOTATION.db
