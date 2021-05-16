This program **merges two or more %(bin)ss together** into a single %(bin)s.

To run this program, the bins that you want to merge must be contained within a single %(collection)s. Just provide the collection name, the %(pan-db)s or %(profile-db)s you're working with, the bins that you want to merge, and the name of the output bin. 

To check what collections and bins are contained in a database, you can either run this program with the flag `--list-collections` or `--list-bins`, or you can run %(anvi-show-collections-and-bins)s.

For example, if you wanted to merge the bins `first_third`, `middle_third`, and `last_third` in a pan-db into a single bin called `complete_bin`, just run 

{{ codestart }}
anvi-merge-bins -p %(pan-db)s \
                -C %(collection)s \
                -b first_third,middle_third,last_third \
                -B complete_bin
{{ codestop }}

Now your collection will contain the bin `complete_bin` and the original bins will be gone forever (unless you had run%(anvi-summarize)s, %(anvi-export-collection)s, or a similar program beforehand)
