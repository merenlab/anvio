This program gives you the **coverage and detection data** for all of the genes found in your %(contigs-db)s, using the short reads data in your %(profile-db)s. 

{{ codestart }}
anvi-export-gene-coverage-and-detection -c %(contigs-db)s \
                                        -p %(profile-db)s \
                                        -O MY_DATA
{{ codestop }}

This will give you a %(coverages-txt)s and a %(detection-txt)s whose file names will begin with `MY_DATA`
