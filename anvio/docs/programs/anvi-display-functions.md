For a given annotation source for %(functions)s, this program will display distribution patterns of unique function names (or accession numbers) across genomes stored in anvi'o databases.

It is a powerful way to analyze differentially occurring functions for any source of annotation that is shared across all genomes.

Currently, %(anvi-display-functions)s can work with any combination of genomes from %(external-genomes)s, %(internal-genomes)s, and %(genomes-storage-db)s.

### Basic Run 

The simplest way to run this program is as follows:

{{ codestart }}
anvi-display-functions -e %(external-genomes)s \
                       --annotation-source KOfam \
                       --profile-db KOFAM-PROFILE.db
{{ codestop }}


You can replace the annotation source based on what is available across your genomes. You can use the program %(anvi-db-info)s to see all available function annotation sources in a given %(contigs-db)s or %(genomes-storage-db)s. Please see %(functions)s for more information on functions and how to obtain them.


{:.notice}
Please note that a %(profile-db)s will be automatically generated for you. Once it is generated, the same profile database can be visualized over and over again using %(anvi-interactive)s in manual mode, without having to retain any other files.

### Aggregating functions using accession IDs

Once it is run, this program essentially aggregates all function names that occur in one or more genomes in the set of genomes found in input sources. The user can ask the program to use accession IDs to aggregate functions rather than function names:

{{ codestart }}
anvi-display-functions -e %(external-genomes)s \
                       --annotation-source KOfam \
                       --profile-db KOFAM-PROFILE.db \
                       --aggregate-based-on-accession
{{ codestop }}

While the default setting which is to use function names will be appropriate for most applications, using accession IDs instead of function names may be important for specific applications. There may be an actual difference between using functions or accession to aggregate data since multiple accession IDs in various databases may correspond to the same function. This may lead to misleading enrichment analyses downstream as identical function annotations may be over-split into multiple groups. Thus, the default aggregation method uses function names.

### Aggregating functions using accession IDs

In some cases a gene may be annotated with multiple functions. This is a decision often made at the function annotation tool level. For instance %(anvi-run-ncbi-cogs)s may yield two COG annotations for a single gene because the significance score for both hits may exceed the default cutoff. While this can be useful in %(anvi-summarize)s output where things should be most comprehensive, having some genes annotated with multiple functions and others with one function may over-split them (since in this scenario a gene with COGXXX and COGXXX;COGYYY would end up in different bins). Thus, %(anvi-display-functions)s will will use the best hit for any gene that has multiple hits. But this behavior can be turned off the following way:

{{ codestart }}
anvi-display-functions -e %(external-genomes)s \
                       --annotation-source KOfam \
                       --profile-db KOFAM-PROFILE.db \
                       --aggregate-using-all-hits
{{ codestop }}

---

The user also can keep only functions that occur in more than a minimum number of genomes:

{{ codestart }}
anvi-display-functions -e %(external-genomes)s \
                       --annotation-source KOfam \
                       --profile-db KOFAM-PROFILE.db \
                       --min-occurrence 5
{{ codestop }}

### Combining genomes from multiple sources

Alternatively, you can run the program by combining genomes from multiple sources:

{{ codestart }}
anvi-display-functions -e %(external-genomes)s \
                       -i %(internal-genomes)s \
                       -g %(genomes-storage-db)s \
                       --annotation-source KOfam \
                       --profile-db KOFAM-PROFILE.db

{{ codestop }}

### A real-world example

Assume we have a list of %(external-genomes)s that include three different species of *Bifidobacterium*. Running the following command,

{{ codestart }}
anvi-display-functions --external-genomes Bifidobacterium.txt \
                       --annotation-source COG20_FUNCTION \
                       --profile-db COG20-PROFILE.db \
                       --min-occurrence 3
{{ codestop }}

Would produce the following display by default, where each layer is one of the genomes described in the %(external-genomes)s file, and each item is a unique function name that occur in `COG20_FUNCTION` (which was obtained by running %(anvi-run-ncbi-cogs)s on each %(contigs-db)s in the external genomes file) that were found in more than three genomes:

[![Example output](../../images/anvi-display-functions-01.png){:.center-img .width-50}](../../images/anvi-display-functions-01.png)

The outermost layer shows the function names:

[![Example output](../../images/anvi-display-functions-02.png){:.center-img .width-50}](../../images/anvi-display-functions-02.png)

After a quick prettification through the %(interactive)s interface, leads to a cleaner display of three distinct species in this group, and functions that are uniquely enriched in either of them: 

[![Example output](../../images/anvi-display-functions-03.png){:.center-img .width-80}](../../images/anvi-display-functions-03.png)

Now the resulting %(profile-db)s can be used by %(anvi-interactive)s to re-visualize these data, or can be shared with the community without sharing the underlying contigs databases.
