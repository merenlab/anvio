This program **takes in a %(functions)s artifact to create a %(functions-txt)s.** Basically, if you want to take the information in your %(functions)s artifact out of anvi'o or give it to a fellow anvi'o user (for them to [import it](http://merenlab.org/software/anvio/help/programs/anvi-import-functions/) into their own project), you get that information using this command. 

Simply provide the %(contigs-db)s that has been annotated with %(functions)s: 

{{ codestart }}
anvi-export-functions -c %(contigs-db)s 
{{ codestop }}

You can also get annotations for only a specific list of sources. For example:

{{ codestart }}
anvi-export-functions -c %(contigs-db)s \
                      --annotation-sources source_1,source_2,source_3
{{ codestop }}

