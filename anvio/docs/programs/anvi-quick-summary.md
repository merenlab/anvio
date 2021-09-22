This program is a quicker, but less comprehensive, alternative to %(anvi-summarize)s. It is used to summarize basic read recruitment statistics (like detection and coverage) from many single profiles that are all associated with the same %(contigs-db)s.

Given a list of samples (single profiles) and a collection, `anvi-quick-summary` will compute the per-sample weighted average of each statistic for each bin in the collection. This is an average of the statistic value over each split in the bin, _weighted by the split length_.

The output will be a text file, and you can find details about its format by clicking on %(quick-summary)s.

### Basic usage

In addition to your list of %(single-profile-db)ss, you must provide this program with their corresponding contigs database and a collection name.
{{ codestart }}
anvi-quick-summary -c %(contigs-db)s -C %(collection)s PROFILE_1.db PROFILE_2.db PROFILE_3.db [...]
{{ codestop }}

The program will summarize the same collection across all of your profile databases. However, it will use only the first profile database in the argument list to learn about what is in the collection, so it is not exactly necessary to have this collection defined for all of the other profile databases (though one could argue that it is a good idea to do this regardless...). The collection name you provide to this program must be a collection that is present in at least the first profile database in the argument list. In the example above, only `PROFILE_1.db` is strictly required to include the collection you wish to summarize (though all other profiles must contain the same splits as this first profile, which should not be a problem if you generated them all in the same way).

### Choosing a different output prefix

By default, the output file will be prefixed with the collection name that you provided. If you wish to set a different prefix, you can use the `--output-file-prefix`, or `-O`, parameter:
{{ codestart }}
anvi-quick-summary -c %(contigs-db)s -C %(collection)s -O new_prefix PROFILE_1.db PROFILE_2.db PROFILE_3.db [...]
{{ codestop }}

No matter what, the output will end in `*-quick_summary.txt`. There is no option to change this. Sorry (not sorry).

### Choosing which statistics to summarize

The default statistics that will be summarized are detection and something called 'mean_coverage_Q2Q3' (which is [this](https://merenlab.org/2017/05/08/anvio-views/#mean-overage-q2q3)). You can choose which statistics to summarize by providing them as a comma-separated list (no spaces in the list) to the `--stats-to-summarize`, or `-S`, parameter:
{{ codestart }}
anvi-quick-summary -c %(contigs-db)s -C %(collection)s -S  PROFILE_1.db PROFILE_2.db PROFILE_3.db [...]
{{ codestop }}

Each statistic will get its own column in the output file.

If you are not sure which statistics are available to choose from, just provide some ridiculous, arbitrary string (that cannot possibly be a name of a statistic) to this flag, and you will get an error message that includes a list of the available statistics. Or, you can just look at this example error message (but no guarantees that the list in this example will be the same as whatever you would get by doing it yourself. Just sayin'.)
```
Config Error: The statistic you requested, cattywampus, does not exist. Here are the options
              to choose from: std_coverage, mean_coverage, mean_coverage_Q2Q3, detection,
              abundance, variability
```

If you are curious about the statistics in the list, many of them have definitions in [this blog post](https://merenlab.org/2017/05/08/anvio-views).

## Common errors

### Existing file error

If the output file already exists, you will encounter the following error:
```
File/Path Error: AppendableFile class is refusing to open your file at test-quick_summary.txt
                 because it already exists. If you are a user, you should probably give Anvi'o a
                 different file name to work with. If you are a programmer and you don't want
                 this behavior, init this class with `fail_if_file_exists=False` instead.
```
You can either provide a different file prefix using the `-O` parameter, as the error message suggests, or you can simply delete the existing file and re-run your command.

### Missing table error

If you get an error that looks like this:
```
Config Error: The database at [PROFILE.db] does not seem to have a table named
              `detection_splits` :/ Here is a list of table names this database knows:
              [...]
```

That means your profile databases are not the correct version. The tables we are accessing in this program were introduced in profile database version 36. So the solution to this error is to update your databases to at least that version, using %(anvi-migrate)s. :)
