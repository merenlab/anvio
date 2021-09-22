This program finds [palindromes](https://en.wikipedia.org/wiki/Palindromic_sequence) in any DNA sequence. It will search for palindromes that exceeds a minimum length defined by the user (default is 10) and that do not have more internal mismatches than the number of those allowed by the user (default is 0). It will optionally report the search results as a %(palindromes-txt)s.

Please note that this program will only report perfect palindromes, where the identity and order of nucleotides on one strand match to those on the complementary strand. The current version is **not** designed to search for special cases of palindromes that form [hairpins](https://en.wikipedia.org/wiki/Stem-loop) (because we didn't need them (but if you do, let us know and we will do something)).


{:.notice}
The algorithm processes about 1,000,000 nts in every 2 seconds, and about 10,000,000 nts in every 9 seconds on a laptop computer.

### Sequence input sources

%(anvi-search-palindromes)s can use multiple different sequence sources.

#### Contigs database

In this mode %(anvi-search-palindromes)s will go through every contig sequence in a given %(contigs-db)s.

{{ codestart }}
anvi-search-palindromes -c %(contigs-db)s \
                        --output-file %(palindromes-txt)s
{{ codestop }}

#### FASTA file

Alternatively, you can use a %(fasta)s file as input.

{{ codestart }}
anvi-search-palindromes --fasta-file %(fasta)s \
                        --output-file %(palindromes-txt)s
{{ codestop }}

#### DNA sequence

Those who are lazy can also pass a DNA sequence for quick searches:

{{ codestart }}
anvi-search-palindromes --dna-sequence (.. A DNA SEQUENCE OF ANY LENGTH ..)
{{ codestop }}


### Verbose output

If you provide an `--output-file` parameter, your results will be stored into a %(palindromes-txt)s file for downstream analyses. If you do not provide an output file, or explicitly asked for a verbose output with the flag `--verbose`, you will see all your palindromes listed on your screen.

Here is an example with a single sequence and no output file path:

{{ codestart }}
%(anvi-search-palindromes)s --dna-sequence CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG
{{ codestop }}

```
SEARCH SETTINGS
===============================================
Minimum palindrome length ....................: 10
Number of mismatches allowed .................: 0
Be verbose? ..................................: Yes


58 nts palindrome at "0:58"
===============================================
Sequence .....................................: (a sequence does not have a name)
Num mismatches ...............................: 0
FDW ..........................................: CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG
ALN ..........................................: ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
REV ..........................................: CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG

SEARCH RESULTS
===============================================
Total number of sequences processed ..........: 1
Total number of palindromes found ............: 1
Longest palindrome ...........................: 58
```

Here is another example with a %(contigs-db)s, an output file path, and the `--verbose` flag:

{{ codestart }}
%(anvi-search-palindromes)s -c CONTIGS.db \
                         --min-palindrome-length 50 \
                         --max-num-mismatches 1 \
                         --output-file palindromes.txt \
                         --verbose
{{ codestop }}

```
SEARCH SETTINGS
===============================================
Minimum palindrome length ....................: 56
Number of mismatches allowed .................: 2
Be verbose? ..................................: Yes


56 nts palindrome at "1067670:1067726"
===============================================
Sequence .....................................: Day17a_QCcontig1
Num mismatches ...............................: 2
FDW ..........................................: TAAGAAAATACCCTTGGTTTATTAAGACGACTTAATAAACCAAGGGTATTTTTTTA
ALN ..........................................: |||x||||||||||||||||||||||x||x||||||||||||||||||||||x|||
REV ..........................................: TAAAAAAATACCCTTGGTTTATTAAGTCGTCTTAATAAACCAAGGGTATTTTCTTA

58 nts palindrome at "187419:187477"
===============================================
Sequence .....................................: Day17a_QCcontig6
Num mismatches ...............................: 2
FDW ..........................................: TCAACACTAAAAAAGCAAACAAAACGCGTAAACGTTTTGTTTGCTTTTTTAGTGTTGA
ALN ..........................................: ||||||||||||||||||||||||||xx||xx||||||||||||||||||||||||||
REV ..........................................: TCAACACTAAAAAAGCAAACAAAACGTTTACGCGTTTTGTTTGCTTTTTTAGTGTTGA

58 nts palindrome at "323612:323670"
===============================================
Sequence .....................................: Day17a_QCcontig74
Num mismatches ...............................: 0
FDW ..........................................: CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG
ALN ..........................................: ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
REV ..........................................: CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG

SEARCH RESULTS
===============================================
Total number of sequences processed ..........: 4,189
Total number of palindromes found ............: 3
Longest palindrome ...........................: 58

Output file ..................................: palindromes.txt

```


### Programmer access

Just like everything else in anvi'o, you can access the functionality the program `anvi-search-palindromes` offers without using the program itself by inheriting an instance from the `Palindromes` class and use it in your own Python scripts.

Here is an example:

``` python
# import argparse to pass arguments to the class
import argparse

# `Palindromes` is the class we need
from anvio.sequencefeatures import Palindromes

# we also import `Progress` and `Run` helper classes from the terminal
# module to ask the class to print no output messages to our workspace
# (this is obviously optional)
from anvio.terminal import Progress, Run

# get an instance for the case of a contigs database, and process everything in it.
# this example is with an anvi'o contigs db, but you can also pass a FASTA file
# via `fasta_file='FILE.fa'` instead of `contigs_db='CONTIGS.db'`:
p = Palindromes(argparse.Namespace(contigs_db='CONTIGS.db', min_palindrome_length=50), run=Run(verbose=False), progress=Progress(verbose=False))
p.process()

# once the processing is done, the palindromes are stored in a member dictionary:
print(p.palindromes)

>>> {'Day17a_QCcontig23': [{'start': 51287, 'end': 51337, 'length': 50, 'palindrome': 'ATAAATAAACAGAGGCCTTAGAAATATTTCTAAGGCCTCTGTTTATTTAT', 'matches': '||||||||||||||||||||||||||||||||||||||||||||||||||', 'num_mismatches': 0}], 
     'Day17a_QCcontig33': [{'start': 309701, 'end': 309757, 'length': 56, 'palindrome': 'TAAATAAGTTACAATAATAATTGTTATCGATAACAATTATTATTGTAACTTATTTA', 'matches': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'num_mismatches': 0}], 
     'Day17a_QCcontig60': [{'start': 213563, 'end': 213613, 'length': 50, 'palindrome': 'CCTGACATGGCAAAACCCTCTACCNNGGTAGAGGGTTTTGCCATGTCAGG', 'matches': '||||||||||||||||||||||||||||||||||||||||||||||||||', 'num_mismatches': 0}], 
     'Day17a_QCcontig74': [{'start': 323612, 'end': 323670, 'length': 58, 'palindrome': 'CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG', 'matches': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'num_mismatches': 0}], 
     'Day17a_QCcontig90': [{'start': 227448, 'end': 227504, 'length': 56, 'palindrome': 'CGAGACATGATTGAGCGCCGTGACGGTCGACCGTCACGGCGCTCAATCATGTCTCG', 'matches': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'num_mismatches': 0}]}


# alternatively you can get an instance without any input files,
p = Palindromes()

# and perhaps set some values,
p.min_palindrome_length = 14
p.max_num_mismatches = 1

# to go through some sequences of your liking:
some_sequences = {'a_sequence': 'CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG',
                  'antoher_sequence': 'AAATCGGCCGATTT',
                  'sequence_with_no_palindrome': 'AAAAAAAAAAAAAA'}

# in this case (where there are no input files) you can call the function `find`,
# rather than `process`, to populate the `p.palindromes` dictionary:
for sequence_name in some_sequences:
    p.find(some_sequences[sequence_name], sequence_name=sequence_name)

# tadaaa:
print(p.palindromes)

>>> {'a_sequence': [{'start': 0, 'end': 58, 'length': 58, 'palindrome': 'CATTGACGTTGACGGCGTCCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG', 'matches': '|||||||||||||||||x||||||||||||||||||||||x|||||||||||||||||', 'num_mismatches': 1}], 
     'antoher_sequence': [{'start': 0, 'end': 14, 'length': 14, 'palindrome': 'AAATCGGCCGATTT', 'matches': '||||||||||||||', 'num_mismatches': 0}]}
```
