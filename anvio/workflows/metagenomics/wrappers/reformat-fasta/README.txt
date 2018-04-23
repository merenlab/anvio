# Wrapper for anvi-script-reformat-fasta.

## Example:

```
rule reformat_fasta:
    '''
        Reformating the headers of the contigs fasta files in order to
        give contigs meaningful names; so that if the group name is
        'MYSAMPLE01', the contigs would look like this:
        > MYSAMPLE01_000000000001
        > MYSAMPLE01_000000000002
    '''
    version: 1.0
    log: LOGS_DIR + "/{group}-reformat_fasta.log"
    input:
        dir = ASSEMBLY_DIR + "/{group}_TEMP"
        contigs = ASSEMBLY_DIR + "/{group}_TEMP/final.contigs.fa"
    output:
        contig = protected(ASSEMBLY_DIR + "/{group}/{group}-contigs.fa"),
        report = ASSEMBLY_DIR + "/{group}/{group}-reformat-report.txt"
    params: prefix = "{group}"
    wrapper:
        "wrappers/reformat-fasta"
```
