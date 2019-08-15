
# Before you start

### Setup your environment

#### For Barhall

```conda activate /project2/meren/shared/virtual-envs/anvio-quentin```

```git pull```

#### For every other envrinoment

Go to taxonomy branch

```git checkout anvio-taxonomy```

decompress data

```cd /anvio/data/misc```

```tar -zxvf SCG.tar.gz```

### Upgrade your database version.

``` anvi-migrate-db CONTIGS.db```

``` anvi-migrate-db PROFILE.db```

##### if it wasn't done before

``` anvi-run-hmms -c contigs.db -I Bacteria_71```


# taxonomy assignment



The taxonomic estimation is a two-step process.



### Alignement SCGs

the first steps consists in aligning the sequences of the Single copy core genes detected by HMMer.


```anvi-diamond-for-taxonomy -c contigs.db```


Note that the program will use only two CPU and threats by default, especially if you have multiple of them available to you, you should use the --num-threads and --num-core-by-threads parameter. It significantly improves the runtime. --num-core-by-threads will be the number of core allow for diamond's alignment. --num-threads will be the number of alignment running at the same time.

For example, if you have 20 cores available, you can do :

```anvi-diamond-for-taxonomy -c contigs.db -T 4 --num-process 5 ```


### Taxonomy estimation

for this second and last step the taxnomy estimation will use the result of blast sore in the contigs database.

```anvi-estimate-taxonomy -c contigs.db -o output.tsv```

By giving just your contigs.db or using --metagenome the taxonomy assigment will be on the gene level


```anvi-estimate-taxonomy -c contigs.db \
                          -p PROFILE.db \
                          -C default \
                          -o output.tsv```
