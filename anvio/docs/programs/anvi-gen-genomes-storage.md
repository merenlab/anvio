### From external genomes

{{ codestart }}
anvi-gen-genomes-storage -e %(external-genomes)s \
                         -o %(genomes-storage-db)s
{{ codestop }}

### From internal genomes

{{ codestart }}
anvi-gen-genomes-storage -i %(internal-genomes)s \
                         -o %(genomes-storage-db)s
{{ codestop }}

### From internal and external genomes

{{ codestart }}
anvi-gen-genomes-storage -i %(internal-genomes)s \
                         -e %(external-genomes)s \
                         -o %(genomes-storage-db)s
{{ codestop }}

See also %(anvi-pan-genome)s, which computes a pangenome from a %(genomes-storage-db)s.
