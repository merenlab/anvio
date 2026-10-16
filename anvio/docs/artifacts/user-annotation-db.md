A directory created by %(anvi-setup-user-annotation-db)s that contains one or more
custom functional annotation databases prepared from user-provided HMM profiles or
protein FASTA files.

## Structure

```
<annotation_db_dir>/
├── manifest.json          # database registry consumed by anvi-run-user-annotation
├── logs/                  # DIAMOND blastp log files written during anvi-run-user-annotation
├── hmm/
│   └── <db_name>/         # standard anvi'o HMM source directory per HMM database
│       ├── genes.hmm.gz
│       ├── genes.txt
│       ├── kind.txt
│       ├── noise_cutoff_terms.txt
│       ├── reference.txt
│       └── target.txt
└── diamond/
    └── <db_name>.dmnd     # binary DIAMOND database per FASTA-based database
                           # (also present when a companion_fasta was set for an HMM entry)
```

## Contents

`manifest.json` is a JSON file that maps each database name to its metadata:

```json
{
  "MyHMMs": {
    "type": "hmm",
    "source_path": "/original/path/to/models.hmm",
    "hmm_dir": "/annotation_db_dir/hmm/MyHMMs",
    "num_models": 42,
    "noise_cutoff_terms": "--cut_tc",
    "added_on": "2026-05-20",
    "companion_diamond": {
      "type": "diamond",
      "source_path": "/original/path/to/sequences.faa",
      "dmnd_path": "/annotation_db_dir/diamond/MyHMMs.dmnd",
      "num_sequences": 300,
      "added_on": "2026-05-20"
    }
  },
  "MyProteins": {
    "type": "diamond",
    "source_path": "/original/path/to/proteins.faa",
    "dmnd_path": "/annotation_db_dir/diamond/MyProteins.dmnd",
    "dmnd_base": "/annotation_db_dir/diamond/MyProteins",
    "num_sequences": 1000,
    "added_on": "2026-05-20"
  }
}
```

The `companion_diamond` sub-entry is only present when a `companion_fasta` column was supplied
in the input TSV for an HMM database. When present, %(anvi-run-user-annotation)s runs both
an HMM search and a DIAMOND search, producing two sources in the `gene_functions` table:
`MyHMMs_HMM` and `MyHMMs_DIAMOND`.

## Produced by

%(anvi-setup-user-annotation-db)s

## Consumed by

%(anvi-run-user-annotation)s
