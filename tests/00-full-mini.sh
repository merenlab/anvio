cd 00-full-mini
rm -rf tmp
mkdir tmp

for f in 6M 7M 9M
do
    papi-init-bam 204_3contigs_"$f".bam -o tmp/204-$f
    papi-profile -i tmp/204-$f.bam -o tmp/204-$f
done

papi-merge-multiple-runs tmp/204*/RUNINFO.cp -o tmp/204-MERGED
papi-interactive-binning -r tmp/204-MERGED/RUNINFO-mean_coverage.cp
rm -rf tmp
cd ..
