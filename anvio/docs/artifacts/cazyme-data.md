This stores a local copy of the data from the [dbCAN2 CAZyme HMM database](https://bcb.unl.edu/dbCAN2/download/Databases/) for functional annotation.

It is required to run %(anvi-run-cazymes)s and is set up on your computer by the program %(anvi-setup-cazymes)s. 

## Database version options

The following versions of dbCAN CAZyme HMMs are currently supported by anvi'o: V9–V13.

You can find them at: `https://bcb.unl.edu/dbCAN2/download/Databases/V*/dbCAN-HMMdb-V*.txt`

Please see the help page for %(anvi-setup-cazymes)s if you want to learn how to choose which version to set up on your computer.

## Notes for developers

[dbCAN3](https://bcb.unl.edu/dbCAN2/) has multiple other files available for download `https://bcb.unl.edu/dbCAN2/download/Databases/` that could be useful for future development of cazyme annotation in anvi'o including: 
- Substrate prediction based on dbCAN-sub majority voting: `dbCAN_sub.hmm`
- CAZyme family information: `FamInfo.txt.08022020.xls`