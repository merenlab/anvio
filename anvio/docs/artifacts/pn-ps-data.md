This describes the output of %(anvi-script-calculate-pn-ps-ratio)s, which calculates the pN/pS ratio for each gene in a %(contigs-db)s. 

{:.notice}
See the page for %(anvi-script-calculate-pn-ps-ratio)s for an explanation of the pN/pS ratio 

This describes a directory that contains the following four files: 

`pNpS.txt`: a long-format table of the pN/pS values, along with the groupby variables:

|   | corresponding_gene_call | sample_id   | pNpS_reference       |
| - | ----------------------- | ----------- | -------------------- |
| 0 | 1744                    | ANE_004_05M | 0.043503524536208836 |
| 1 | 1744                    | ANE_004_40M | 0.043628712253629943 |
| 2 | 1744                    | ANE_150_05M | 0.03810623760551494  |
| 3 | 1744                    | ANE_150_40M | 0.040815421982026576 |

`pN.txt`: a long-format table of the pN values, along with the groupby variables:

|   | corresponding_gene_call | sample_id   | pN_reference       |
| - | ----------------------- | ----------- | ------------------ |
| 0 | 1744                    | ANE_004_05M | 11.827627600424583 |
| 1 | 1744                    | ANE_004_40M | 11.106801744995472 |
| 2 | 1744                    | ANE_150_05M | 9.62355553228605   |
| 3 | 1744                    | ANE_150_40M | 10.067364489809782 |

`pS.txt`: a long-format table of the pS values, along with the groupby variables:

|   | corresponding_gene_call | sample_id   | pS_reference       |
| - | ----------------------- | ----------- | ------------------ |
| 0 | 1744                    | ANE_004_05M | 271.87745651689016 |
| 1 | 1744                    | ANE_004_40M | 254.57551165909962 |
| 2 | 1744                    | ANE_150_05M | 252.54541348089631 |
| 3 | 1744                    | ANE_150_40M | 246.6558962502711  |

`num_SCVs.txt`: a long-format table of the number of SCVs belonging to each group:

|   | corresponding_gene_call | sample_id   | num_SCVs |
| - | ----------------------- | ----------- | -------- |
| 0 | 1744                    | ANE_004_05M | 180      |
| 1 | 1744                    | ANE_004_40M | 166      |
| 2 | 1744                    | ANE_150_05M | 162      |
| 3 | 1744                    | ANE_150_40M | 160      |

