This program **downloads and organizes a local copy of PHROGs HMM profiles and annotations for viral function annotation.** It generates a %(phrogs-data)s artifact, which is required to run %(anvi-run-phrogs)s.

### Set up PHROGs data
{{ codestart }}
anvi-setup-phrogs
{{ codestop }}

By default, this data is stored at `anvio/data/misc/PHROGs`. To set up this data in a non-default location, run
{{ codestart }}
anvi-setup-phrogs --phrogs-data-dir path/to/location
{{ codestop }}

If you already have a %(phrogs-data)s artifact and want to redownload this data, run

{{ codestart }}
anvi-setup-phrogs --reset
{{ codestop }}
