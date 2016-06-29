rm -rf test-output
mkdir test-output

echo "name	contigs_db_path" > test-output/contig_dbs.txt
for g in E_coli_BL21 E_coli_B_REL606 E_coli_O111
do
    gzip -d $g.fa.gz --to-stdout > test-output/$g.fa
    anvi-gen-contigs-database -f test-output/$g.fa -o test-output/$g.db -L -1 
    anvi-run-hmms -c test-output/$g.db
    echo "$g	$g.db" >> test-output/contig_dbs.txt
done
