This program returns the nucleotide-level coverage data for a specific set of the splits or gene in your %(profile-db)s. 

If you want to get the coverage data for all splits in your %(profile-db)s, run %(anvi-export-splits-and-coverages)s with the flag `--splits-mode`. 

Simply provide a %(profile-db)s and %(contigs-db)s pair and specify which splits, or gene, you want to look at. You have three ways to do this: 

1.  Provide a single split name. (You can list all splits available with `--list-splits`)

{{ codestart }}
anvi-get-split-coverages -p %(profile-db)s \
                         -c %(contigs-db)s \
                         -o %(coverages-txt)s \ 
                         --split-name Day17a_QCcontig9_split_00003
{{ codestop }}


2. Provide both the name of a %(bin)s and the %(collection)s it is contained in. 

{{ codestart }}
anvi-get-split-coverages -p %(profile-db)s \
                         -c %(contigs-db)s \
                         -o %(coverages-txt)s \ 
                         -b %(bin)s \
                         -C %(collection)s
{{ codestop }}

You can list all collections available with `--list-collections` or all bins in a collection with `--list-bins`. Alternatively, you could run %(anvi-show-collections-and-bins)s on your %(profile-db)s to get a more comprehensive overview. 

3. Provide a gene caller id and a flanking size (bp).

{{ codestart }}
anvi-get-split-coverages -p %(profile-db)s \
                         -c %(contigs-db)s \
                         -o %(coverages-txt)s \ 
                         --gene-caller-id 25961 \
                         --flank-length 500
{{ codestop }}
