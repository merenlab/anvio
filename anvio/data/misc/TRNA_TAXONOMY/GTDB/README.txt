We used the genome collection of GlobDB v226 (https://globdb.org), which itself contains all genome representatives of GTDB v226.
GlobDB already provides contigs-databases with annotated tRNA, through the command anvi-scan-trnas.

We made an external-genomes-txt file with all the GTDB genomes present in GlobDB and used anvi-get-sequences-for-hmm-hits to extract the sequences for all tRNAs:
```
anvi-get-sequences-for-hmm-hits -e external-genomes.txt --hmm-sources Transfer_RNAs -o all_tRNAs.fa
```

The headers of that fasta file contains the information about the anticodon and the genome's id:
```
>Ile_GAT___Transfer_RNAs___4068027d3b59b533bad5e550e79d768a08c97ce5fbceb337f7a6171e bin_id:GCA_000008085|source:Transfer_RNAs|e_value:92.8|contig:hash112713c5_GCA_000008085_000000000001|gene_callers_id:hash112713c5_612|start:7000|stop:7074|length:74
```

Using that information, we can create as many fasta file as anticodon and only keep the genome's name as the defline:
```
mkdir -p ANTICODON_SEARCH_DATABASES

awk -v outdir="ANTICODON_SEARCH_DATABASES" '
    /^>/ {
        # parse anticodon and genome id from header
        match($0, />[^_]+_([A-Z]{3})___.*bin_id:([^|]+)/, m)
        anticodon = m[1]
        genome = m[2]
        # print sequence under genome ID
        if (anticodon != "" && genome != "") {
            fname = outdir "/" anticodon
            print ">" genome >> fname
        }
        next
    }
    {
        print >> fname
    }
' all_tRNAs.fa

# gzip all outputs
gzip ANTICODON_SEARCH_DATABASES/*
```

Finally, we can use the taxonomy file provided in GlobDB to construct the table mapping the genome's name/accession to a taxonomy:
```
grep "^>" all_tRNAs.fa | sed -E 's/.*bin_id:([^|]+).*/\1/' | sort -u > genome_ids.txt
grep -Ff genome_ids.txt globdb_r226_taxonomy.tsv > ACCESSION_TO_TAXONOMY.txt
gzip ACCESSION_TO_TAXONOMY.txt && rm genome_ids.txt
```
