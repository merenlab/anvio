A directory of **user-defined metabolism data**, created by the user for estimating metabolism on custom metabolic pathways. The program %(anvi-setup-user-modules)s takes this directory as input and creates a %(modules-db)s out of the data within, for use by %(anvi-estimate-metabolism)s.

Instructions for creating this data directory and using it to estimate completeness of custom (ie, non-KEGG) metabolic pathways can be found below.

## A step-by-step guide to creating your own metabolic modules for anvi-estimate-metabolism

If you want to define your own metabolic pathway so that you can estimate its completeness in genomes, MAGs, and metagenomes, follow the steps below!

### 1. Find the enzymes

What you need first is a list of enzyme accession numbers. For each reaction in your metabolic pathway, figure out what enzyme(s) or enzyme complexes (if any) are required to catalyze the reaction. Then, for each of these enzymes and/or components of enzyme complexes, figure out if they are present in common databases like [NCBI COG](https://www.ncbi.nlm.nih.gov/research/cog), [KEGG KOfam](https://www.genome.jp/tools/kofamkoala/), or [PFAM](http://pfam.xfam.org/). If so, mark down their accession numbers in those databases. If not, you may need to create your own HMM profile for the enzyme (and create an accession number for it).

Also, think about how you will annotate each enzyme, because for each one you will need to write down its functional annotation source in the module file in step 3. Here is a short guide to common annotation sources:

Enzyme comes from... | annotation program | ANNOTATION_SOURCE
|:---|:---|:---|
KEGG KOfam | %(anvi-run-kegg-kofams)s | Kofam
NCBI COG (2020) | %(anvi-run-ncbi-cogs)s with `--cog-version COG20`| COG20_FUNCTION
NCBI COG (2014) | %(anvi-run-ncbi-cogs)s with `--cog-version COG14`| COG14_FUNCTION
archaeal COGs (2014) | %(anvi-run-ncbi-cogs)s with `--cog-version arCOG14`| arCOG14_FUNCTION
PFAM | %(anvi-run-pfams)s | Pfam
CAZy | %(anvi-run-cazymes)s | CAZyme
custom HMMs | %(anvi-run-hmms)s with `--hmm-source` and `--add-to-functions-table` parameters | name of directory given to `--hmm-source`
other annotation strategy | %(anvi-import-functions)s | source defined in input file

### 2. Define the module

You need to write a DEFINITION string for the module. This string should be in the style of KEGG MODULE definitions, which are described [here](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/#what-data-is-used-for-estimation). Briefly, you will put the enzyme accessions in order of their corresponding reactions in the metabolic pathway. Different steps (reactions) in the pathway should be separated by spaces, and alternative enzymes that can catalyze the same reaction should be separated by commas. You can use parentheses to distinguish alternatives with multiple steps. For enzyme complexes, all components should be in one string, with essential components separated by '+' signs and non-essential components separated by '-' signs.

### 3. Write a module file

Put all the information about your metabolic pathway into a text file. The file format and types of information you need to include are discussed [here](https://merenlab.org/software/anvio/help/main/programs/anvi-setup-user-modules/#how-do-i-format-the-module-files). At minimum, you need to pick an identifier (ENTRY) and NAME for the module, include your DEFINITION string from step 2, write an ORTHOLOGY line and an ANNOTATION_SOURCE line for each enzyme and/or enzyme component, and write a CLASS string to categorize your module into its class/category/subcategory. The module file should be given the same name as the identifier in the ENTRY line, and this identifier should not be the same as any module in the KEGG database.

{:.notice}
Check out %(anvi-script-gen-user-module-file)s for a way to automatically format your user module files.

### 4. Set up the USER_MODULES.db

Once you have created a module file for each metabolic pathway you are interested in, you should put these files within a folder called `modules`, within a parent directory (that can have any name you choose), as described [here](https://merenlab.org/software/anvio/help/main/programs/anvi-setup-user-modules/#input-directory-format). This parent directory is the %(user-modules-data)s directory. Then you should run the program %(anvi-setup-user-modules)s and provide this directory to the `--user-modules` parameter. If all goes well, you will end up with a database called `USER_MODULES.db` in this folder.

### 5. Annotate your contigs database(s)

Before you can estimate metabolism, you will need to annotate your contigs database(s) with each annotation source that you used to define your modules. This will require running one or more annotation programs, as described in the table given for step 1 above. If you want to quickly remind yourself of which annotation sources are required for your metabolic modules, you can run %(anvi-db-info)s on the `USER_MODULES.db`. But don't worry - if you forget one, you will get a helpful error message telling you what you missed when you try to run %(anvi-estimate-metabolism)s.

Since estimation will be run on KEGG data, too, you will have to make sure you also run %(anvi-run-kegg-kofams)s on your database(s), if you haven't already, _UNLESS_ you are choosing to skip KEGG estimation by using the `--only-user-modules` parameter for %(anvi-estimate-metabolism)s.

### 6. Estimate the completeness of your pathways

The last step is to run %(anvi-estimate-metabolism)s and provide this directory to the `--user-modules` parameter. This program will estimate the completeness of the metabolic modules defined in the `USER_MODULES.db` (by default, this will be in addition to the KEGG modules from the %(kegg-data)s directory. But, as mentioned above, you can specify `--only-user-modules` to only estimate on your own data).

## A toy example

For this example, we will be creating a completely FAKE, biologically-nonsensical Frankenstein of a metabolic pathway. This is not anything you should be putting in your own `USER_MODULES.db`; it only exists for demonstrating the steps above, and particularly so you have a reference for how to handle the different annotation sources mentioned in step 1.

First, let's select 5 different enzymes. Typically at this step you would use your biological knowledge of a real metabolic pathway to determine the specific set of enzymes catalyzing the reactions of the pathway - but for this toy example, we're going to use random enzymes that come from a variety of annotation sources, not enzymes that actually work together biologically in a real cell. So we'll go with a couple of KOfams, a COG, a PFAM, and a TIGRFAM (to demonstrate the 'other' annotation strategy). Here is the list of accessions: K01657, K01658, COG1362, PF06603.14, and TIGR01709.2.

It doesn't matter what they are or what they do. What matters is that we will learn how to annotate each one. So let's talk about their annotation sources. K01657 and K01658 will both come from `KOfam`, and COG1362 will come from the 2020 distribution of the COGs database, so its source will be `COG20_FUNCTION`. PF06603.14 is a PFAM, so it _could_ come from the `Pfam` source. But let's suppose we don't want to waste our precious computational resources on running %(anvi-run-pfams)s when we are only interested in one enzyme from this database. Instead, we'll make a custom HMM profile for this particular enzyme by following the directions on [creating HMM sources from ad hoc PFAM accessions](https://merenlab.org/software/anvio/help/main/artifacts/hmm-source/#creating-anvio-hmm-sources-from-ad-hoc-pfam-accessions), and then we will annotate it using %(anvi-run-hmms)s. In this case, the annotation source will be the name of the directory we make using %(anvi-script-pfam-accessions-to-hmms-directory)s, so we need to pick a name for it - let's call it `METABOLISM_HMM`. Last but not least, what about the TIGRFAM enzyme TIGR01709.2? Anvi'o doesn't have a program for annotating TIGRFAMs, but we can annotate our gene sequences with TIGRFAM using [Interproscan](https://www.ebi.ac.uk/interpro/search/sequence/), compile the results into a %(functions-txt)s, and import those annotations into our contigs database using %(anvi-import-functions)s. We'll put the source `TIGRFAM` in the %(functions-txt)s file.

Great, so now that we have our enzymes and we know how we will annotate them, it's time for step 2 - creating the module DEFINITION string. Again, this is not going to be a biologically-realistic metabolic pathway, but an example to demonstrate the different ways of representing steps in a pathway.

Let's say the first reaction in our pathway is catalyzed by an enzyme complex made up of two essential components, K01657 and K01658. We represent this step by the string "K01657+K01658" (no spaces between the components). Suppose the next part of the pathway can _either_ be one reaction catalyzed by the enzyme PF06603.14, _or_ it can be a two-step reaction in which the first reaction is catalyzed by COG1362 and the second is catalyzed by TIGR01709.2. We use a comma to separate the alternatives, and since the second option requires two different steps (that will be separated by a space), we surround the second option with parentheses to make sure both steps are considered as the alternative. It looks like this: "PF06603.14,(COG1362 TIGR01709.2)".

So our full module DEFINITION string is "K01657+K01658 PF06603.14,(COG1362 TIGR01709.2)".

Now we need to put this information into a module file. We'll give the module the identifier `UD0042` (UD for 'user-defined', and 42 because 42 is the answer to life, the universe, and everything), and this will also be the name of the file.

Here is the module file. Any information that we didn't discuss above has been filled in to demonstrate the formatting requirements:
```
ENTRY       UD0042
NAME        Frankenstein pathway for demo purposes
DEFINITION  K01657+K01658 PF06603.14,(COG1362 TIGR01709.2)
ORTHOLOGY   K01657  anthranilate synthase component I [EC:4.1.3.27]
            K01658  anthranilate synthase component II [EC:4.1.3.27]
            PF06603.14  UpxZ
            COG1362  Aspartyl aminopeptidase
            TIGR01709.2  type II secretion system protein GspL
CLASS       User modules; Demo set; Frankenstein metabolism
ANNOTATION_SOURCE  K01657  KOfam
                    K01658  KOfam
                    PF06603.14  METABOLISM_HMM
                    COG1362  COG20_FUNCTION
                    TIGR01709.2  TIGRFAM
///
```

If you were actually going to use this pathway, this is how you could create the %(user-modules-data)s directory and then the `USER_MODULES.db`:

{{ codestart }}
mkdir USER_METABOLISM
mkdir USER_METABOLISM/modules
vi USER_METABOLISM/modules/UD0042
\#copy the above into this file, save and quit
anvi-setup-user-modules --user-modules USER_METABOLISM/
{{ codestop }}

You would see the following output after %(anvi-setup-user-modules)s completed:
```
Modules database .............................: A new database, USER_METABOLISM/USER_MODULES.db, has been created.
Number of modules ............................: 1
Number of entries ............................: 14
Number of parsing errors (corrected) .........: 0
Number of parsing errors (uncorrected) .......: 0
Annotation sources required for estimation ...: COG20_FUNCTION, METABOLISM_HMM, KOfam, TIGRFAM
```
As expected, if we want to use this modules database for estimating completeness of our Frankenstein pathway, we would need to annotate our %(contigs-db)s of interest with the four annotation sources we discussed above. And that, in fact, is the next step.

The first two annotation sources we discussed are easy, because we don't need to do anything besides run the designated anvi'o program for KEGG KOfams and NCBI COGs, respectively:

{{ codestart }}
%(anvi-run-kegg-kofams)s -c CONTIGS.db \
                      --num-threads 4
%(anvi-run-ncbi-cogs)s -c CONTIGS.db \
                    --num-threads 4
{{ codestop }}

Annotating PF06603.14 requires an extra step, because we first need to create a custom HMM for this enzyme. Luckily, there is another anvi'o program to do that. We give the enzyme accession to %(anvi-script-pfam-accessions-to-hmms-directory)s, and we make sure to set the output directory name to be the same as the annotation source string that we put in the module file:

{{ codestart }}
%(anvi-script-pfam-accessions-to-hmms-directory)s --pfam-accessions-list PF06603.14 \
                                               -O METABOLISM_HMM
%(anvi-run-hmms)s -c CONTIGS.db \
               -H METABOLISM_HMM \
               --add-to-functions-table \
               --num-threads 4
{{ codestop }}

Please note that you _must_ use the `--add-to-functions-table` parameter when you use %(anvi-run-hmms)s, otherwise the annotations for PF06603.14 will not be stored in the proper database table and %(anvi-estimate-metabolism)s will not be able to find them later. Also, if you use the %(anvi-script-pfam-accessions-to-hmms-directory)s program to create your custom HMM profiles, you should make sure that the accessions in the resulting `genes.txt` file are matching to the corresponding enzyme accessions in the module file, because those are the accessions that will be put into your contigs database.

Finally, to annotate TIGR01709.2 we need to take our (hypothetical) Interproscan results and convert them into a %(functions-txt)s file. You can visit that page for a lengthier discussion of the file format, but let's say the TIGR01709.2 annotations in that file looked like this:

|gene_callers_id|source|accession|function|e_value|
|:--|:--:|:--:|:--|:--:|
|7|TIGRFAM|TIGR01709.2|type II secretion system protein GspL|1.5e-75|
|23|TIGRFAM|TIGR01709.2|type II secretion system protein GspL|3.4e-20|

The things that are especially critical here is that the `accession` matches to the accession in the module DEFINITION, ORTHOLOGY, and ANNOTATION_SOURCE lines, and that the `source` matches to the source string in the module ANNOTATION_SOURCE line(s).

Suppose the file is called `TIGRFAM_annotations.txt`. Then you can import those annotations, like so:

{{ codestart }}
%(anvi-import-functions)s -c CONTIGS.db \
                       -i TIGRFAM_annotations.txt
{{ codestop }}

Once this is done, you are ready to estimate the pathway's completeness! Here is the command:

{{ codestart }}
%(anvi-estimate-metabolism)s -c CONTIGS.db \
                          --user-modules USER_METABOLISM/ \
                          -O frankenstein
{{ codestop }}

If you did this, the results for module UD0042 would appear at the end of the resulting 'modules' output file, after the estimation results for KEGG modules.
