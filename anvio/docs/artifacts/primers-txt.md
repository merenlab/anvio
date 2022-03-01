A **TAB-delimited** file to describe primer sequences. A primer sequence can be exact (such as `ATCG`), or fuzzy (such as `AT.G`, which would match any combination of `ATAG`, `ATTG`, `ATCG`, or `ATGG`). Fuzzy primers are defined by regular expressions, properties of which are explained best in [the Python documentation](https://docs.python.org/3/library/re.html), or cheatsheets like [this one](https://www.debuggex.com/cheatsheet/regex/python).

This file type includes two required and any number of user-defined optional columns.

The following three columns are **required** for this file type:

* `name`: a single-word name for a given primer,
* `primer_sequence`: the primer sequence.

### An example primers-txt file

Here is an example file with three primers:

|name|primer_sequence|
|:--|:--|
|PR01|AA.A..G..G..G.CCG.C.A.C|
|PR02|AACACCGCAGTCCATGAGA|
|PR03|A[TC]A[CG]T[ATC]TCGAGC|
