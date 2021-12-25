A directory of **user-defined metabolism data**, created by the user for estimating metabolism on custom metabolic pathways. The program %(anvi-setup-user-modules)s takes this directory as input and creates a %(modules-db)s out of the data within, for use by %(anvi-estimate-metabolism)s.

Instructions for creating this data directory and using it to estimate completeness of custom (ie, non-KEGG) metabolic pathways can be found below.

## A step-by-step guide to creating your own metabolic modules for anvi-estimate-metabolism

If you want to define your own metabolic pathway so that you can estimate its completeness in genomes, MAGs, and metagenomes, follow the steps below!

1. Find the enzymes

What you need first is a list of enzyme accession numbers. For each reaction in your metabolic pathway, figure out what enzyme(s) or enzyme complexes (if any) are required to catalyze the reaction. Then, for each of these enzymes and/or components of enzyme complexes, figure out if they are present in common databases like [NCBI COG](https://www.ncbi.nlm.nih.gov/research/cog), [KEGG KOfam](https://www.genome.jp/tools/kofamkoala/), or [Pfam](http://pfam.xfam.org/). If so, mark down their accession numbers in those databases. If not, you may need to create your own HMM profile for the enzyme (and create an accession number for it).

Also, think about how you will annotate each enzyme, because for each one you will need to write down its functional annotation source in the module file in step 3. Here is a short guide to common annotation sources:

Enzyme comes from... | annotation program | ANNOTATION_SOURCE
|:---|:---|:---|
KEGG KOfam | %(anvi-run-kegg-kofams)s | Kofam
NCBI COGs (2020) | %(anvi-run-ncbi-cogs)s | COG20_FUNCTION
NCBI COGs (2014) | %(anvi-run-ncbi-cogs)s | COG14_FUNCTION
Pfams | %(anvi-run-pfam)s | Pfam
custom HMMs | %(anvi-run-hmm)s with `--hmm-source` and `--add-to-functions-table` parameters | name of directory given to `--hmm-source`
other annotation strategy | %(anvi-import-functions)s | source defined in input file

2. Define the module

You need to write a DEFINITION string for the module. This string should be in the style of KEGG MODULE definitions, which are described [here](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/#what-data-is-used-for-estimation). Briefly, you will put the enzyme accessions in order of their corresponding reactions in the metabolic pathway. Different steps (reactions) in the pathway should be separated by spaces, and alternative enzymes that can catalyze the same reaction should be separated by commas. You can use parentheses to distinguish alternatives with multiple steps. For enzyme complexes, all components should be in one string, with essential components separated by '+' signs and non-essential components separated by '-' signs.

3. Write a module file

Put all the information about your metabolic pathway into a text file. The file format and types of information you need to include are discussed [here](https://merenlab.org/software/anvio/help/main/programs/anvi-setup-user-modules/#how-do-i-format-the-module-files). At minimum, you need to pick an identifier (ENTRY) and NAME for the module, include your DEFINITION string from step 2, write an ORTHOLOGY line and an ANNOTATION_SOURCE line for each enzyme and/or enzyme component, and write a CLASS string to categorize your module into its class/category/subcategory. The module file should be given the same name as the identifier in the ENTRY line, and this identifier should not be the same as any module in the KEGG database.

4. Set up the USER_MODULES.db

Once you have created a module file for each metabolic pathway you are interested in, you should put these files within a folder called `modules`, within a parent directory (that can have any name you choose), as described [here](https://merenlab.org/software/anvio/help/main/programs/anvi-setup-user-modules/#input-directory-format). This parent directory is the %(user-modules-data)s directory. Then you should run the program %(anvi-setup-user-modules)s and provide this directory to the `--input-dir` parameter. If all goes well, you will end up with a database called `USER_MODULES.db` in this folder.

5. Annotate your contigs database(s)

Before you can estimate metabolism, you will need to annotate your contigs database(s) with each annotation source that you used to define your modules. This will require running one or more annotation programs, as described in the table given for step 1 above. If you want to quickly remind yourself of which annotation sources are required for your metabolic modules, you can run %(anvi-db-info)s on the `USER_MODULES.db`. But don't worry - if you forget one, you will get a helpful error message telling you what you missed when you try to run %(anvi-estimate-metabolism)s.

Since estimation will always be run on KEGG data, too, you will have to make sure you also run %{anvi-run-kegg-kofams}s on your database(s), if you haven't already.

6. Estimate the completeness of your pathways

The last step is to run %(anvi-estimate-metabolism)s and provide this directory to the `--input-dir` parameter. This program will estimate the completeness of the metabolic modules defined in the `USER_MODULES.db`, in addition to the KEGG modules from the %(kegg-data)s directory.
