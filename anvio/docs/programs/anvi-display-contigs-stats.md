This program **lets you look at various stats in your %(contigs-db)s(s)**. 

When you run with all default parameters, as so,

{{ codestart }}
anvi-display-contigs-stats %(contigs-db)s
{{ codestop }}

It will open the interactive interface. From there, you'll be able to see a breakdown of the genes in your contigs database. This will look something like this: 

![An example of the anvi'o interface for contigs stats](../../images/contigs-stats-interface-example.png)

Let's walk through the stats displayed here: 

At the top of the page are two graphs: 
* The bars in the top graph represent every integer N and L statistic from 1 to 100. The y-axis is the respective N length and the x-axis is the L length (the exact values can be seen by hovering over each bar). For a definition of the L and N statistics, see below. In other words, if you had sorted your contigs by length (from longest to shortest), and walked through each one, every time you had seen another 1 percent of your total dataset, you would add a bar to the graph showing the number of contigs that you had seen (the L statistic) and the length of the one you were looking at at the moment (the N statistic). 
* The lower part of the graph tells you about which HMM hits your contigs database has. Each column is a gene in the %(hmm-source)s, and the graph tells you how many hits each gene has in your data. (Hover your mouse over the graph to see the specifics of each gene.) The sidebar shows you how many of the genes in this graph were seen exactly that many times. For example, in the graph above, for the Bacteria_71 %(hmm-source)s, a lot of genes were detected 9-11 times, so those bars are longer. This helps you estimate about how many of these genomes there are in your contigs database (so here, there is likely around 9-11 bacteria genomes in this contigs database). 

Below this are the **contigs stats**:
- Total length in nucleotides
- Number of contigs in your database
- The number of contigs that are of varying lengths. (for example "Num Contigs > 2.5 kb" gives you the numberof contigs that are longer than 2500 base pairs)
- The length of the longest and shortest contig in your database in nucleotides
- The number of genes in your contigs (as predicted by [Prodigal](https://github.com/hyattpd/Prodigal))
- L50, L75, L90: If you ordered the contigs in your database from longest to shortest, these stats describe the *number of contigs* you would need to go through before you had looked at a certain percent of a genome. For example, L50 describes the number of contigs you would have to go through before you reached 50 percent of the entire dataset. 
- N50, N75, N90:  If you ordered the contigs in your database from longest to shortest, these stats describe the *length of the contig* you would be looking when you had looked at a certain percent of a genome. For example, N50 describes the length of contig you would be on when you reached 50 percent of the entire genome length. 
- The number of HMM hits in your contigs. This goes through every %(hmm-source)s and gives the number of hits its genes had in all of your contigs. Basically, this is the number of hits that is given in the lower graph at the top of the page. 
- The number of genomes that we predict are in your sample, based on how many hits the single copy core genes got from the various %(hmm-source)ss. See the description of the lower graph above, or [this blog post](http://merenlab.org/2015/12/07/predicting-number-of-genomes/) for more information. 

This interface is espeically useful if you want to compare multiple databases, since you can view all of their stats stimultaneously.

You can also change various server configuration settings when you run this command. 

### The interactive interface is great, but I need text

Then you're still in the right place. Just add the tag `--report-as-text` and you'll get a lovely tab-deliminated output in the path provided:

{{ codestart }}
anvi-display-contigs-stats %(contigs-db)s \ 
                    --report-as-text \
                    -o path/to/my_cool_text_file.txt
{{ codestop }}
