This program **takes a %(functions)s artifact to create a %(functions-txt)s.** Essentially, if you want to extract the information in your %(functions)s artifact from anvi'o or share it with a fellow anvi'o user (for them to [import it](http://merenlab.org/software/anvio/help/programs/anvi-import-functions/) into their own project), you can obtain that information using this command. 

Simply provide the %(contigs-db)s that has been annotated with %(functions)s: 

{{ codestart }}
anvi-export-functions -c %(contigs-db)s 
{{ codestop }}

You can also get annotations for only specific sources. For example:

{{ codestart }}
anvi-export-functions -c %(contigs-db)s \
                      --annotation-sources source_1,source_2,source_3
{{ codestop }}

