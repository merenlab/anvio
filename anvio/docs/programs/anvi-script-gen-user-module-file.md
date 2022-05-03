Given an %(enzymes-list-for-module)s file, this script will produce a properly-formatted module file for use in %(anvi-setup-user-modules)s.

## Basics
You should provide this script with an accession ID for your module (which will become the module file name) (`-I`); a name for the module (`-n`); a categorization which includes module class, category, and sub-category separated by semicolons (`-c`); an %(enzymes-list-for-module)s file listing the component enzymes, their annotation sources, and their functions (`-e`); and a definition which puts these enzymes in the proper order (`-d`):

{{ codestart }}
anvi-script-gen-user-module-file -I "UD0023" \
                  -n "Frankenstein pathway for demo purposes" \
                  -c "User modules; Demo set; Frankenstein metabolism" \
                  -e %(enzymes-list-for-module)s \
                  -d "K01657+K01658 PF06603.14,(COG1362 TIGR01709.2)"
{{ codestop }}

Then this script will produce a properly-formatted module file, which in this case would be called `UD0023` and would look like this (see %(enzymes-list-for-module)s for the example file containing these enzymes):

```
ENTRY       UD0023
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
\\\
```

## Automatically generating the module definition

The module definition parameter is not required, and if you do not provide one, the definition will be generated with each enzyme in the input file as a different 'step' of the module. This option may be especially appropriate for generating what KEGG calls "signature modules", in which each enzyme is not technically part of a metabolic pathway, but instead the module represents a functionally-related set of enzymes (like all tRNA modification enzymes, for instance).

Here is an example command without the definition parameter:

{{ codestart }}
anvi-script-gen-user-module-file -I "UD0023" \
                  -n "Frankenstein pathway for demo purposes" \
                  -c "User modules; Demo set; Frankenstein metabolism" \
                  -e %(enzymes-list-for-module)s \
{{ codestop }}

And here is the module file that this would produce:

```
ENTRY       UD0023
NAME        Frankenstein pathway for demo purposes
DEFINITION  K01657 K01658 PF06603.14 COG1362 TIGR01709.2
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
\\\
```
