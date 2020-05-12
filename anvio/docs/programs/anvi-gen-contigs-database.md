### Create a contigs database from a FASTA file

```
anvi-gen-contigs-database -f %(contigs-fasta)s -o %(contigs-db)s
```

A FASTA file that contains one or more sequences. These sequences may belong
to a single genome, or could be contigs obtained from an assembly.

### Create a contigs database with external gene calls

```
anvi-gen-contigs-database -f %(contigs-fasta)s -o %(contigs-db)s -e %(external-gene-calls)s
```

