This programs takes a %(variability-profile-txt)s and generates the information necessary to visualize its contents with %(anvi-interactive)s. 

Specifically, this program outputs a directory that contains a %(profile-db)s, a %(view-data)s artifact, and a %(dendrogram)s. For example, if you ran this program like so: 

{{ codestart }}
anvi-script-snvs-to-interactive -o OUTPUT_DIR \
                                %(variability-profile)s
{{ codestop }}

Then, you can open the interactive interface by running 

{{ codestart }}
anvi-interactive --manual-mode \
                 -p OUTPUT_DIR/profile.db \
                 --tree OUTPUT_DIR/tree.txt \
                 --view-data OUTPUT_DIR/view.txt
{{ codestop }}

## Other parameters 

### Using Only a Subset of the Input

By default, all variability positions in your variability profile are considered. However, if the input is too large (i.e. more than 25,000 variability positions), the runtime on this program will be very long and the results won't display well. So, there are several ways to remove variability positions from  the input to get under this threshold: 

1. Ignore positions with with certain departures from the consensus sequence (with `--min-departure-from-consensus` and `--max-departure-from-consensus`)
2. Ignore positions with with certain departures from the reference sequence (with `--min-departure-from-reference` and `--max-departure-from-reference`)
3. Ignore positions in all non-coding regions with the flag `--only-in-genes`. 

If you still have more positions than you can tell the program to pick a random subset of the input with the parameter `--random` followed by a seed integer.

### Modifying the Output

By the default, the output data will use the departure from consensus values. If instead you want to look at the departure from the reference, just add the falg `--display-dep-from-reference`
