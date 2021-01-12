MAKE SURE YOU RUN EVERYTHING DIRECTLY IN THIS DIRECTORY THIS WAY:

    ./00_RUN.sh

If everything goes well, running `00_RUN.sh` will generate BAM files and a
`contigs.fa` for anvi'o mini test using all the contigs that are under
`contigs/` directory and all samples that are declared through their .ini
files under `samples` directory.

Please note that you will need `bowtie2` and the program `reads-for-assembly`
for this to run. If you don't have reads-for-assembly, you can get it this way:

    mkdir ~/github
    cd ~/github
    git clone https://github.com/merenlab/reads-for-assembly.git
    cd -

Only contigs in `_orig.fa` files are going to be used for read recruitment. Which
means, you can copy a contig into a new file that does not end with `_orig.fa`,
create a variant of it by modifying it with nucleotide-level or structural changes,
and set the desired coverage of this new variant in a given sample by mentioning
it in an `.ini` file.

Once it is done running, you can copy `output/contigs.fa` and `output/*-RAW.bam` to
`anvio/tests/sandbox` directory, and run the following to see the effect:

    cd
    anvi-self-test --suite mini -o NEW_MINI_TEST_RESULTS

If you add new samples, you will have change all test sh files accordingly.

Unless you manually move `output/contigs.fa` and `output/*-RAW.bam` into the
`anvio/tests/sandbox` directory, these files will have no effect on future
tests.
