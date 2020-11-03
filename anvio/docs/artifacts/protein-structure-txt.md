This is a Protein Data Bank (`X.pdb`) file that describes the structure of a protein as stored in your %(structure-db)s. This is the output of running %(anvi-export-structures)s.

This file format has its own [Wikipedia page](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)), as well as pages on PDB-101 [for beginners](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner's-guide-to-pdb-structures-and-the-pdbx-mmcif-format) and for [coordinates specifically](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/dealing-with-coordinates), but is also briefly explained here. 

The header describes the title (if one exists), the type of data (denoted by `EXPDTA`), and any free-form annotations (denoted by `REMARK` ). In Anvi'o, these are primarily MODELLER information calculated when you ran %(anvi-gen-structure-database)s. 

Most of the data will describe the position of individual atoms (denoted by `ATOM`) in your protein, where columns 6, 7, and 8 describe the three dimensional coordinate of the atom. The rest of the columns describe information like what position that atom is in the amno acid. 

`TER` statements separate independent chains from each other.

Here is an example: 

    EXPDTA    THEORETICAL MODEL, MODELLER 9.22 2020/10/13 14:38:54
    REMARK   6 MODELLER OBJECTIVE FUNCTION:       255.0071
    REMARK   6 MODELLER BEST TEMPLATE percent SEQ ID:  73.077
    ATOM      1  N   MET A   1       4.009  -3.600  -0.411  1.00 59.26           N
    ATOM      2  CA  MET A   1       5.250  -3.864  -1.173  1.00 59.26           C
    ATOM      3  CB  MET A   1       6.409  -3.005  -0.631  1.00 59.26           C
    ATOM      4  CG  MET A   1       6.204  -1.504  -0.854  1.00 59.26           C
    ATOM      5  SD  MET A   1       7.545  -0.444  -0.229  1.00 59.26           S
    ATOM      6  CE  MET A   1       6.982  -0.449   1.495  1.00 59.26           C
    ATOM      7  C   MET A   1       5.617  -5.306  -1.058  1.00 59.26           C
    ATOM      8  O   MET A   1       4.900  -6.175  -1.552  1.00 59.26           O
    ATOM      9  N   SER A   2       6.753  -5.603  -0.397  1.00 49.70           N
    ATOM     10  CA  SER A   2       7.165  -6.971  -0.290  1.00 49.70           C
    ATOM     11  CB  SER A   2       8.547  -7.150   0.362  1.00 49.70           C
    ATOM     12  OG  SER A   2       9.546  -6.534  -0.437  1.00 49.70           O
    ATOM     13  C   SER A   2       6.184  -7.694   0.556  1.00 49.70           C
    ATOM     14  O   SER A   2       5.954  -7.346   1.714  1.00 49.70           O
    ATOM     15  N   GLU A   3       5.553  -8.718  -0.037  1.00103.21           N
    ATOM     16  CA  GLU A   3       4.632  -9.540   0.676  1.00103.21           C
    ATOM     17  CB  GLU A   3       3.856 -10.490  -0.249  1.00103.21           C
    ATOM     18  CG  GLU A   3       4.774 -11.467  -0.988  1.00103.21           C
    ATOM     19  CD  GLU A   3       3.918 -12.407  -1.826  1.00103.21           C
    ATOM     20  OE1 GLU A   3       2.672 -12.402  -1.638  1.00103.21           O
    ATOM     21  OE2 GLU A   3       4.502 -13.146  -2.663  1.00103.21           O
    ATOM     22  C   GLU A   3       5.402 -10.410   1.614  1.00103.21           C
    ...
    ATOM    594  C   ASN A  79      19.969 -15.504   4.267  1.00 61.57           C
    ATOM    595  O   ASN A  79      21.042 -14.862   4.423  1.00 61.57           O
    ATOM    596  OXT ASN A  79      19.857 -16.742   4.474  1.00 61.57           O
    TER     597      ASN A  79
    END
