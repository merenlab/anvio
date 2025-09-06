This program **provides you with the coverage information in your %(profile-db)s as external files**. Essentially, if you want to extract that information from your %(profile-db)s for use outside of anvi'o, this program is for you. 

Once you input your %(profile-db)s and the %(contigs-db)s you used to generate it, it will create a %(contigs-fasta)s that lists your contigs for you, as well as a %(coverages-txt)s, which describes your coverage information. 

{{ codestart }}
anvi-export-splits-and-coverages -p %(profile-db)s \
                                 -c %(contigs-db)s
{{ codestop }}

If your coverages are skewed by outlier positions, consider using Q2Q3-coverages instead.

{{ codestart }}
anvi-export-splits-and-coverages -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 --use-Q2Q3-coverages
{{ codestop }}

### Contigs or splits?

*Wondering what the difference is? Check out [our vocab page](http://merenlab.org/vocabulary/#split).*

By default, this program will give you the sequences of your splits, but will examine coverage data in terms of the parent contig. If you want to get coverage information for your splits, use `--splits-mode`. Alternatively, you can ask the program to `--report-contigs` to look at contig sequences instead.
