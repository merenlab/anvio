This program **associates genes in your %(contigs-db)s with viral functions from the [PHROGs database](https://phrogs.lmge.uca.fr/) by running PHROGs HMMs.**

Before you run this program, set up PHROGs data on your computer with %(anvi-setup-phrogs)s.

To run, provide a %(contigs-db)s. If you stored %(phrogs-data)s in a custom location, provide that path with `--phrogs-data-dir`. The output is a %(functions)s artifact.

Here is a default run:

{{ codestart }}
anvi-run-phrogs -c %(contigs-db)s \
                --phrogs-data-dir %(phrogs-data)s
{{ codestop }}

By default, this uses `hmmsearch` and HMMER's trusted GA thresholds (`--cut_ga`). You can choose `hmmscan` instead:

{{ codestart }}
anvi-run-phrogs -c %(contigs-db)s \
                --phrogs-data-dir %(phrogs-data)s \
                --hmmer-program hmmscan
{{ codestop }}
