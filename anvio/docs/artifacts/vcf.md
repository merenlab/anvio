This represents a file in the Variant Call Format, which is a standard format for storing sequence variations like SNVs. 

You can convert the information in a %(variability-profile-txt)s to %(vcf)s with the program %(anvi-script-variability-to-vcf)s. 

### What's in this file? 

For more details, you can check out the [VCF wikipedia page](https://en.wikipedia.org/wiki/Variant_Call_Format). 

#### Header

Briefly, this file's header (marked by `##` at the beginning of each line) contains various metadata. This includes the date, link to the reference file, contig information, etc. It also contains 
- what information will be reported (denoted by `INFO`). 
- what additional filters will be run on each SNV (denoted by `FILTER`). For example, marking which variants are below a certain quality threshold. 
- what format to display additional data in (denoted by `FORMAT`). 

#### Body

The body of the file contains identifying information for the variation (the chromosome, position and ID), the identity of the position in the reference and alternative alleles present in your data. Following this is a quality score for your data and the additional information specified by the header. 
