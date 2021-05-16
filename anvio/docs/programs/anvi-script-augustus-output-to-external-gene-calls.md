This program converts a gene call file from [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) (as an %(augustus-gene-calls)s artifact) to an anvi'o %(external-gene-calls)s artifact. 

This essentially just reformats the data in the %(augustus-gene-calls)s artifact (for example, removing the UTR information) so that it can be read by other anvi'o programs. 

A run of this program will look something like this:

{{ codestart }}
anvi-script-augustus-output-to-external-gene-calls -i %(augustus-gene-calls)s
                                                   -o %(external-gene-calls)s
{{ codestop }}

Here is an example of what the resulting %(external-gene-calls)s file will look like (from the gff file used as an example on the %(augustus-gene-calls)s page):  

    gene_callers_id    contig       start    stop    direction    partial    call_type    source      version    aa_sequence
    0                  unnamed-1    56       1252    f            0          1            AUGUSTUS    v3.3.3     MSEGNAAGEPSTPGGPRPLLTGARGLIGRRPAPPLTPGRLPSIRSRDLTLGGVKKKTFTPNIISRKIKEEPKEEVTVKKEKRERDRDRQREGHGRGRGRPEVIQSHSIFEQGPAEMMKKKGNWDKTVDVSDMGPSHIINIKKEKRETDEETKQILRMLEKDDFLDDPGLRNDTRNMPVQLPLAHSGWLFKEENDEPDVKPWLAGPKEEDMEVDIPAVKVKEEPRDEEEEAKMKAPPKAARKTPGLPKDVSVAELLRELSLTKEEELLFLQLPDTLPGQPPTQDIKPIKTEVQGEDGQVVLIKQEKDREAKLAENACTLADLTEGQVGKLLIRKSGRVQLLLGKVTLDVTMGTACSFLQELVSVGLGDSRTGEMTVLGHVKHKLVCSPDFESLLDHKHR

