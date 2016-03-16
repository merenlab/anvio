rm -rf test-output
mkdir test-output

echo "name	path" > test-output/contig_dbs.txt
for g in Escherichia_coli_BL21 Escherichia_coli_B_REL606 Escherichia_coli_O111
do
    gzip -d $g.gz --to-stdout > test-output/$g.fa
    anvi-gen-contigs-database -f test-output/$g.fa -o test-output/$g.db -L -1 
    echo "$g	$g.db" >> test-output/contig_dbs.txt
done
