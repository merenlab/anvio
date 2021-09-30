This program finds [palindromes](https://en.wikipedia.org/wiki/Palindromic_sequence) in any DNA sequence. It will search for palindromes that mathes criteria listed by the user (i.e., minimum lenght of the palindromic sequences, maximum number of mismatches, and minimum distance between the two palindromic regions). The program will print out its findings (and tribulations) and will optionally report the search results as a %(palindromes-txt)s.

Please note that this program can find both perfect palindromes (i.e., the identity and order of nucleotides on one strand match to those on the complementary strand) and special cases of palindromes that form [hairpins](https://en.wikipedia.org/wiki/Stem-loop). You can use the minimum distance parameter to target any group of palindromes (i.e., minimum distance of 0 will report only perfect palindromes).

{:.notice}
The speed of the algorithm will depend on the minimum palindrome length parameter. The shorter the palindrome length, the longer the processing time. Searching for palindromes longer than 50 nts in a 10,000,000 nts long sequence takes about 4 seconds on a laptop.

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
Minimum gap length ...........................: 0
Be verbose? ..................................: Yes


58 nts palindrome"
===============================================
1st sequence [start:stop] ....................: [0:58]
2nd sequence [start:stop] ....................: [0:58]
Number of mismatches .........................: 0
Distance between .............................: 0
1st sequence .................................: CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG
ALN ..........................................: ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2nd sequence .................................: CATTGACGTTGACGGCGACCGGTCGGTGATCACCGACCGGTCGCCGTCAACGTCAATG

SEARCH RESULTS
===============================================
Total number of sequences processed ..........: 1
Total number of palindromes found ............: 1
Longest palindrome ...........................: 58
Most distant palindrome ......................: 0
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
Minimum palindrome length ....................: 50
Number of mismatches allowed .................: 1
Minimum gap length ...........................: 0
Be verbose? ..................................: Yes

147 nts palindrome"
===============================================
1st sequence [start:stop] ....................: [268872:269019]
2nd sequence [start:stop] ....................: [269631:269778]
Number of mismatches .........................: 1
Distance between .............................: 759
1st sequence .................................: TTTCGTAATACTTTTTTGCAGTAGGCATCAAATTGGTGTTGTATAGATTTCTCATTATAATTTTGTTGCATGATAATATGCTCCTTTTTCCCCTTTCCACTAATACAACAATCAGAGAGCCCCTTTTTTTCGAAAAAGCTAGAAAAA
ALN ..........................................: |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||x|||||||||
2nd sequence .................................: TTTCGTAATACTTTTTTGCAGTAGGCATCAAATTGGTGTTGTATAGATTTCTCATTATAATTTTGTTGCATGATAATATGCTCCTTTTTCCCCTTTCCACTAATACAACAATCAGAGAGCCCCTTTTTTTCGAAAAAACTAGAAAAA

SEARCH RESULTS
===============================================
Total number of sequences processed ..........: 11
Total number of palindromes found ............: 1
Longest palindrome ...........................: 147
Most distant palindrome ......................: 759

Output file ..................................: palindromes.txt
```


### Programmer access

Just like everything else in anvi'o, you can access the functionality the program `anvi-search-palindromes` offers without using the program itself by inheriting an instance from the `Palindromes` class and use it in your own Python scripts.

Here is an example, first with an input file and then an ad hoc sequence. Starting with the file (i.e., an anvi'o %(contigs-db)s):

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
```

Once the processing is done, the palindromes are stored in a member dictionary, which contains a key for each sequence:

``` python
print(p.palindromes)

>>> {'Day17a_QCcontig1' : [],
     'Day17a_QCcontig2' : [],
     'Day17a_QCcontig4' : [<anvio.sequencefeatures.Palindrome object at 0x7f8d6072f278>],
     'Day17a_QCcontig6' : [],
     'Day17a_QCcontig10': [], 
     'Day17a_QCcontig16': [],
     'Day17a_QCcontig23': [],
     'Day17a_QCcontig24': [],
     'Day17a_QCcontig45': [],
     'Day17a_QCcontig54': [],
     'Day17a_QCcontig97': []}

```

Non-empty arrays are the proper palindromes found in a given sequence, described with an instance of the class `Palindrome` which is defined as the following:

``` python
class Palindrome:
    def __init__(self, run=terminal.Run()):
        self.run=run
        self.first_start = None
        self.fisrt_end = None
        self.first_sequence = None
        self.second_start = None
        self.second_end = None
        self.second_sequence = None
        self.num_mismatches = None
        self.length = None
        self.distance = None
        self.midline = ''
```

Not only you can access to each member variable to deal with them, you can easily display the contents of one using the `display()` function:

``` python
palindrome = p.palindromes['Day17a_QCcontig4'][0]
print(palindrome)

>>> TTTCGTAATACTTTTTTGCAGTAGGCATCAAATTGGTGTTGTATAGATTTCTCATTATAATTTTGTTGCATGATAATATGCTCCTTTTTCCCCTTTCCACTAATACAACAATCAGAGAGCCCCTTTTTTTCGAAAAA (268872:269009) :: TTTCGTAATACTTTTTTGCAGTAGGCATCAAATTGGTGTTGTATAGATTTCTCATTATAATTTTGTTGCATGATAATATGCTCCTTTTTCCCCTTTCCACTAATACAACAATCAGAGAGCCCCTTTTTTTCGAAAAA (269631:269768)

palindrome.display()

>>> 137 nts palindrome"
>>> ===============================================
>>> 1st sequence [start:stop] ....................: [268872:269009]
>>> 2nd sequence [start:stop] ....................: [269631:269768]
>>> Number of mismatches .........................: 0
>>> Distance between .............................: 759
>>> 1st sequence .................................: TTTCGTAATACTTTTTTGCAGTAGGCATCAAATTGGTGTTGTATAGATTTCTCATTATAATTTTGTTGCATGATAATATGCTCCTTTTTCCCCTTTCCACTAATACAACAATCAGAGAGCCCCTTTTTTTCGAAAAA
>>> ALN ..........................................: |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
>>> 2nd sequence .................................: TTTCGTAATACTTTTTTGCAGTAGGCATCAAATTGGTGTTGTATAGATTTCTCATTATAATTTTGTTGCATGATAATATGCTCCTTTTTCCCCTTTCCACTAATACAACAATCAGAGAGCCCCTTTTTTTCGAAAAA
```

Alternatively you can process an ad hoc sequence without any input files,

``` python
p = Palindromes()

# let's set some values for fun,
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

>>> {'a_sequence': [<anvio.sequencefeatures.Palindrome object at 0x7fce807ddb00>],
     'antoher_sequence': [<anvio.sequencefeatures.Palindrome object at 0x7fce807ddc88>],
     'sequence_with_no_palindrome': []}
```

If you are a programmer and need more from this module, please let us know.