This program takes the variability data stored within a %(profile-db)s and compiles it from across samples into a single matrix that comprehensively describes your SNVs, SCVs or SAAVs (a %(variability-profile-txt)s).

This program is described on [this blog post](http://merenlab.org/2015/07/20/analyzing-variability/#the-anvio-way), so take a look at that for more details. 

## Let's talk parameters 

Here is a basic run with no bells or whisles: 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s
{{ codestop }}
