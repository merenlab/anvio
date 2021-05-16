This script **converts a %(variability-profile-txt)s into %(vcf)s (Variant Call Format).** 

It is very easy to run: just provide the input and output paths as so:

{{ codestart }}
anvi-script-variability-to-vcf -i %(variability-profile-txt)s \ 
                               -o %(vcf)s 
{{ codestop }}

Note that to run this, you'll need to have run %(anvi-gen-variability-profile)s with the default nucleotide engine. 
