This program takes an INI file that looks like this:

```ini
[general]
short_read_length = 10
error_rate = 0.05
coverage = 100
contig = CTGTGGTTACGCCACCTTGAGAGATATTAGTCGCGTATTGCATCCGTGCCGACAAATTGCCCAACGCATCGTTCCTTCTCCTAAGTAATTTAACATGCGT
```

When this is used this way:

{{ codestart }}
anvi-script-gen-short-reads -i %(configuration-ini)s \
                            -o %(short-reads-fasta)s
{{ codestop }}

The resulting FASTA file with short reads will cover the `contig` with short reads that are 10 nts long at 100X.
